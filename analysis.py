import nibabel
import numpy
import os
import sys
import csv
import pprint
import scipy

import time
import pickle
import config
import inspect
import pdb

from scipy.stats import pearsonr
from scipy.stats import spearmanr
import statsmodels.sandbox.stats.multicomp
import numpy as np
import pandas

from ontology import *

print "Base script folder:" + config.scriptLocation

R_VALUE_FUNCTION = spearmanr

def load_nifti_data(brainID):
  '''

  Takes the name of a brainID containing two NIFTI files, 'T1.nii' and 'T2.nii'.

  Returns a 3-item list of numpy arrays containing the T1,T2 and T1/T2 measurements.

  '''
  os.chdir(config.figShareFolder)
  NIFTI_FILES = ['m' + brainID + '_T1.nii.gz', 'm' + brainID + '_T2.nii.gz']
  #NIFTI_FILES = [brainID + '_T1.nii.gz', brainID + '_T2.nii.gz'] #use images that have not been bias corrected

  imgs = []
  data = []
  for f in NIFTI_FILES:
    img = nibabel.load(f)
    imgs.append(img)
    #print img.header
    data.append(img.get_data())

  data.append(numpy.divide(data[0]*1.0,data[1]*1.0))
  return data

def get_gene_exp_data_file(directory):
  '''
  Takes the name of a directory containing a txt file with the MRI data.

  Returns a file handle to the open file

  '''
  os.chdir(directory)
  return open('9861.matrix.MRI(xyz).29131 x 946.txt')

def flatten_mri_data(mri_data,coords):
  flat_mri_data = []
  for coord in coords:
    flat_mri_data.append(mri_data[coord])
  return flat_mri_data

def get_coords_and_region_ids_from_gene_exp_data(gene_exp_fh):
  '''
  Returns a list of coordinates (x,y,z tuples) and a dictionary mapping regionIDs to coords
  '''
  header_line = gene_exp_fh.readline().strip()
  coord_strings = header_line.split('\t')
  evaluate_strings = lambda x : eval(x.strip('"').split('|')[1].split(':')[1])
  coord_to_region_map = {}
  for string in coord_strings:
      components = string.strip('"').split('|')
      coord_to_region_map[components[1].split(':')[1]] = eval(components[0].split(':')[1])

  coords = map(evaluate_strings,coord_strings)

  return coords,coord_to_region_map

def get_coords_from_region_id(regionID,coords,coord_to_region_map,ontology,to_exclude_list=None):
  '''
  Returns the list of all coordinates corresponding to a regionID (including all its child regions)

  to_exclude is an optional parameter containing a regionID to exlude from the list
  '''
  regionIDs = ontology.get_all_regionIDs(regionID)
  coords_list = []
  for regionID_iter in regionIDs:
    for k,v in coord_to_region_map.iteritems():
      if v == regionID_iter:
        coords_list.append(k)
  if to_exclude_list:
    excluded_coords = []
    for to_exclude in to_exclude_list:
      to_exclude_coords = get_coords_from_region_id(to_exclude,coords,coord_to_region_map,ontology)
      excluded_coords.append(to_exclude_coords)
    excluded_coords = [item for sublist in excluded_coords for item in sublist]
    return [x for x in coords_list if x not in excluded_coords]

  return coords_list


def rank_regions_by_intensity(coords,coord_to_region_map,ontology,flat_mri_data):
  '''
  Returns a list of all brain regions from the ontology, ranked by the average T2/T1 MRI intensity at each
  point
  '''
  rank = []
  total = 0
  nonzero = 0
  all_intensities = []
  for ID,name in ontology.names.iteritems():

    flag = 0
    if ID == 9299:
      flag = 1

    intensity = 0
    coords_list = get_coords_from_region_id(ID,coords,coord_to_region_map,ontology)

    if coords_list != []:

      if flag:
        print len(coords_list)

      for coord in coords_list:

        index = coords.index(eval(coord))
        intensity += float(flat_mri_data[index])

        if flag:
          all_intensities.append(flat_mri_data[index])
          print "i:", intensity

      intensity = intensity / float(len(coords_list))
      if flag:
          all_intensities.append(intensity)
          print "i:", intensity
          print name
    rank.append((name, intensity))

  print all_intensities
  return sorted(rank,key=lambda x: x[1],reverse=True)

def get_flat_coords_from_region_id(ID,coords,coord_to_region_map,ontology,to_exclude_list=None):
  flat_coords_list = []
  coords_list = get_coords_from_region_id(ID,coords,coord_to_region_map,ontology,to_exclude_list=to_exclude_list)
  for coord in coords_list:
    flat_coords_list.append(coords.index(eval(coord)))
  return flat_coords_list

def correlate_MRI_and_gene_exp_data(flat_mri_data,gene_exp_filename,indices=None):
  '''
  Takes the name of a directory containing a txt file with the gene_expression data and two NIFTI files.

  Reads all three files.

  Correlates the expression data for each gene with the one measure of MRI intensity.

  Returns a list of genes

  Using the file handle and processing one line at a time allows us to avoid loading the entire gene expression text file into memory.

  If the indices parameter is defined, then only the entries at those flattened indices are considered.
  This allows gene_lists to be generated for specific regions.

  '''

  t1 = time.clock()

  expression_values = pandas.read_csv(gene_exp_filename,delimiter="\t")
  IDs = expression_values.iloc[:,0].apply(lambda x: '"' + x + '"')
  adjusted_indices = [x+1 for x in indices]
  expression_values = expression_values.iloc[:,adjusted_indices]
  flat_mri_data_of_interest = [flat_mri_data[i] for i in indices]
  correlate = lambda x: R_VALUE_FUNCTION(x,flat_mri_data_of_interest)
  correlations = expression_values.apply(correlate,axis=1)
  correlations = pandas.concat([IDs,correlations],axis=1).sort_values(0,ascending=False)
  result=zip(list(correlations.iloc[:,0].values),list(correlations.iloc[:,1].values))
  #import pdb; pdb.set_trace()
  t2 = time.clock()

  print "Correlation execution time was", t2-t1

  return result


def analysis(o,files,brain_ids,region_sets,MRI_data_labels,MRI_of_interest):
  #get number of genes
  expression_values = pandas.read_csv(os.path.join(config.processedOutputLocation,files[0]),delimiter="\t")
  gene_count = expression_values.shape[0]
  print "Gene count:" + str(gene_count) 
  
  # dimensions are: MRI modality, regions of interest, genes, statistics (for each 6 correlation, raw p, corrected p (not used)) 
  data_array =  np.zeros((len(MRI_data_labels),len(region_sets),gene_count,3*6+1))  

  print "Data array shape:" + str(data_array.shape)
  # i is the brain
  for i in range(len(brain_ids)):

    MRI_data = load_nifti_data(brain_ids[i])
    gene_exp_fh = open(os.path.join(config.processedOutputLocation,files[i]))
    coords,coord_to_region_map = get_coords_and_region_ids_from_gene_exp_data(gene_exp_fh)
    gene_exp_fh.close()

    for j in range(len(MRI_data)):

      if MRI_data_labels[j] not in MRI_of_interest:
          print "SKIPPING {} FOR BRAIN {}".format(MRI_data_labels[j],brain_ids[i])
          continue

      measure = MRI_data[j]
      flat_mri_data = flatten_mri_data(measure,coords)

      for k in range(len(region_sets)):

        label = brain_ids[i] + "." + MRI_data_labels[j] + "." + region_sets[k][0]
        print 'Correlating MRI and gene expression data for ' + label

        indices = get_flat_coords_from_region_id(region_sets[k][1],coords,coord_to_region_map,o,region_sets[k][2])
        if region_sets[k][0] == "hippocampus" and brain_ids[i] == 10021:
            print len (indices)

        gene_exp_filename = os.path.join(config.processedOutputLocation,files[i])
        correlated_data = correlate_MRI_and_gene_exp_data(flat_mri_data,gene_exp_filename,indices=indices)

        # Get all columns in order (sorted alphabetically by gene)

        data_in_order = map(lambda y: (y[0],y[1]), sorted(correlated_data, key=lambda x: x[0] ))
        gene_names_in_order = map(lambda y: y[0], data_in_order)
        p_values_in_order = map(lambda y: y[1].pvalue, data_in_order)
        correlations_in_order = map(lambda y: y[1].correlation, data_in_order)
        adj_p_values_in_order = statsmodels.sandbox.stats.multicomp.multipletests(p_values_in_order,method="fdr_bh")[1]

        for ind in range(len(gene_names_in_order)):
          data_array[j][k][ind][0] = ind
          data_array[j][k][ind][i+1+len(files)] = str(p_values_in_order[ind])
          data_array[j][k][ind][i+1] = str(correlations_in_order[ind])
          data_array[j][k][ind][i+1+2*len(files)] = str(adj_p_values_in_order[ind])

  print("Writing Files")
  for j in range(3):
    if MRI_data_labels[j] not in MRI_of_interest:
      continue
    for k in range(len(region_sets)):
      label = MRI_data_labels[j] + "." + region_sets[k][0]

      with open(os.path.join(config.resultFolder,  label + ".gene_list.csv"),'w') as f:
        data = data_array[j][k]
        f.write("ID,Correlation."  + (",Correlation.").join(brain_ids) + ",PValue."+ (",PValue.").join(brain_ids) +",PValueAdjusted."+ (",PValueAdjusted.").join(brain_ids) + "\n")

        ind=0
        for gene_entry in data:
          gene_name = gene_names_in_order[int(gene_entry[0])]
          f.write(gene_name + "," + ",".join(map(str,gene_entry[1:])) + "\n")
          ind +=1

def main():
  correlated_data = None
  t1 = time.clock()

  o = Ontology(os.path.join(config.scriptLocation, "data",  "Ontology.csv"))

  print('Getting coordinates from gene expression file.')

  files = config.expression_filenames

  brain_ids = [f.split(".")[0] for f in files]


  #### regions of interest are defined as a 3-tuple (name,ID, ID of excluded subregion)
  #regions_sets_of_interest = [('cortex',4008,None),('cortex_excluding_limbic_lobe',4008,4219),('full_brain',4005,None),("hippocampus",4005,None)]
  #regions_sets_of_interest = [('cortex',4008,None)]
  regions_sets_of_interest = [('cortex_excluding_piriform_hippocampus',4008, [4249, 10142] ),('full_brain',4005,None)]

  MRI_data_labels = ["T1","T2","T1T2Ratio"]
  MRI_of_interest = ["T1","T2","T1T2Ratio"]
  #MRI_of_interest = ["T1"]

  analysis(o,files,brain_ids,regions_sets_of_interest,MRI_data_labels,MRI_of_interest)

  t2 = time.clock()

  print "Total execution time was",t2-t1

if __name__ == '__main__':
  main()
