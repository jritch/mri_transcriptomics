import nibabel
import numpy
import os
import sys
import csv
import pprint
import scipy
#from subprocess import call

import time
import pickle
import config
import inspect, os

from scipy.stats import pearsonr
from scipy.stats import spearmanr
import numpy as np

#import debug

print inspect.getfile(inspect.currentframe()) # script filename (usually with path)
baseProjectFolder = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) + "/" # script directory


print "Base Allen folder:" + config.baseAllenFolder
print "Base script folder:" + baseProjectFolder

R_VALUE_FUNCTION = spearmanr

class ClassName(object):
  """docstring for ClassName"""
  def __init__(self, arg):
    super(ClassName, self).__init__()
    self.arg = arg
    
class Ontology(object):
  """An Ontology has two member variables:

  Attributes:
    names (dict): 
      An (ID,name) dictionary 
    hierarchy (dict):
      a (parentID,[childID]) dictionary 
  """
  def __init__(self,ontology_file_name):
      names = dict()
      hierarchy = dict()
      with open(ontology_file_name) as f:
        rdr = csv.reader(f)
        header = next(rdr)
        for row in rdr:
          ID = int(row[0])
          names[ID] = row[2]
          if row[3]:
            parentID = int(row[3])
            if hierarchy.get(parentID):
                hierarchy[parentID].append(ID)
            else:
              hierarchy[parentID] = [ID]

      self.names = names
      self.hierarchy = hierarchy

  def get_all_regionIDs(self,regionID):
    IDs = [regionID]
    all_IDs = []
    while IDs:
      # print 'IDs is',IDs
      ID = IDs.pop()
      all_IDs.append(ID)
      child_IDs = self.hierarchy.get(ID)
      if child_IDs:
        IDs += child_IDs
    return all_IDs

def load_nifti_data(directory):
  '''

  Takes the name of a directory containing two NIFTI files, 'T1.nii' and 'T2.nii'.

  Returns a 3-item list of numpy arrays containing the T1,T2 and T1/T2 measurements.

  '''
  os.chdir(directory)
  NIFTI_FILES = ['T1.nii','T2.nii']

  imgs = [] 
  data = []
  for f in NIFTI_FILES:
    img = nibabel.load(f)
    imgs.append(img)
    data.append(img.get_data())

  data.append(numpy.divide(data[0],data[1]))

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

def get_coords_from_region_id(regionID,coords,coord_to_region_map,ontology):
  '''
  Returns the list of all coordinates corresponding to a regionID (including all its child regions)
  '''
  regionIDs = ontology.get_all_regionIDs(regionID)
  coords_list = []
  for regionID_iter in regionIDs:
    for k,v in coord_to_region_map.iteritems():
      if v == regionID_iter:
        coords_list.append(k)
  return coords_list


def rank_regions_by_intensity(coords,coord_to_region_map,ontology,flat_mri_data):
  '''
  Returns a list of all brain regions from the ontology, ranked by the average T2/T1 MRI intensity at each
  point 
  '''
  rank = []
  total = 0
  #total = len(ontology.names)
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

  #print len(all_intensities)
  print all_intensities
  #print sum(all_intensities)
  return sorted(rank,key=lambda x: x[1],reverse=True)

def get_flat_coords_from_region_id(ID,coords,coord_to_region_map,ontology):
  flat_coords_list = []
  coords_list = get_coords_from_region_id(ID,coords,coord_to_region_map,ontology)
  for coord in coords_list:
    flat_coords_list.append(coords.index(eval(coord)))
  return flat_coords_list

def correlate_MRI_and_gene_exp_data(flat_mri_data,gene_exp_fh,indices=None):
  '''
  Takes the name of a directory containing a txt file with the gene_expression data and two NIFTI files.

  Reads all three files.

  Correlates(using the pearson correlation coefficient) the expression data for each gene with the one measure of MRI intensity.

  Returns a list of the top ten genes by (genes with p value > 0.05 are not considered)

  Using the file handle and processing one line at a time allows us to avoid loading the entire gene expression text file into memory.

  If the indices parameter is defined, then only the entries at those flattened indices are considered.
  This allows gene_lists to be generated for specific regions.

  '''
  t1 = time.clock()

  f = gene_exp_fh

  # Throw away header line
  gene_exp_fh.readline()

  line = gene_exp_fh.readline()
  IDs = list()
  correlations = list()

  while line:
    entries = line.strip().split('\t') 
    numerical_entries = map(float,entries[1:])
    ID = entries[0]
    IDs.append(ID)
    if not indices:
      correlation = R_VALUE_FUNCTION(numerical_entries,flat_mri_data)
    else:
      correlation = R_VALUE_FUNCTION([numerical_entries[i] for i in indices],
                                      [flat_mri_data[i] for i in indices])

    # look at genes with p-value less than 0.05
    #if correlation[1] <= 0.05:
      #print ID, correlation
    correlations.append([ID,correlation])

    line = gene_exp_fh.readline()

  #gene_exp_fh.close();
  #sort the significantly correlated genes based on correlation and return
  result = sorted(correlations, key=lambda entry: abs(entry[1][0]),reverse=True)#[0:9]
  t2 = time.clock()

  print "Correlation execution time was", t2-t1

  return result

def fisher_p(p_vector):
  return scipy.stats.combine_pvalues(p_vector)[1]

def visualize():
  pass

if __name__ == '__main__':
  correlated_data = None
  t1 = time.clock()

  o = Ontology(config.ontologyFolder + "Ontology.csv")

  print 'Getting coordinates from gene expression file.'

  files = ["10021.matrix.regionID.MRI(xyz).29131 x 893.txt",
           "12876.matrix.regionID.MRI(xyz).29131 x 363.txt",
           "14380.matrix.regionID.MRI(xyz).29131 x 529.txt",
           "15496.matrix.regionID.MRI(xyz).29131 x 470.txt",
           "15697.matrix.regionID.MRI(xyz).29131 x 501.txt",
           "9861.matrix.regionID.MRI(xyz).29131 x 946.txt"]
  
  brain_ids = [f.split(".")[0] for f in files]
  
  #regions_of_interest = [('cortex',4008),('full_brain',4005)]
  regions_of_interest = [('subcortex',4275),('cerebellum',4696)]
  #regions_of_interest = [('cortex',4008)]
  MRI_data_labels = ["T1","T2","T1T2Ratio"]

  #files = ["10021.matrix.regionID.MRI(xyz).29131 x 893.txt"]
  #brain_ids = [f.split(".")[0] for f in files]

  #data_array =  np.array([ [ [ [ [0] * 2 * len(brain_ids)] * 29131] * 1] * 3])
  data_array =  np.zeros((3,4,29131,3*6+1))
  print data_array.shape
  # i is the brain
  for i in range(len(brain_ids)):
    MRI_data = load_nifti_data(config.baseAllenFolder + "normalized_microarray_donor" + brain_ids[i])
    #MRI_data = MRI_data[0]
    gene_exp_fh = open(os.path.join(config.MRIFolder,files[i]))
    coords,coord_to_region_map = get_coords_and_region_ids_from_gene_exp_data(gene_exp_fh)
    gene_exp_fh.close()
    #print len(sorted(get_flat_coords_from_region_id(4696,coords,coord_to_region_map,o)))
    flat_t1t2_ratio_data = flatten_mri_data(MRI_data[2],coords)  
    for j in range(len(MRI_data)):
  #  for measure in MRI_data:
      measure = MRI_data[j]

      print 'Flattening MRI data.'
      flat_mri_data = flatten_mri_data(measure,coords)

      for k in range(len(regions_of_interest)):
        label = brain_ids[i] + "." + MRI_data_labels[j] + "." + regions_of_interest[k][0]
        print 'Correlating MRI and gene expression data for ' + label

        indices = get_flat_coords_from_region_id(regions_of_interest[k][1],coords,coord_to_region_map,o)
        gene_exp_fh = open(os.path.join(config.MRIFolder,files[i]))
        #if not correlated_data:
        correlated_data = correlate_MRI_and_gene_exp_data(flat_mri_data,gene_exp_fh,indices=indices) 
        gene_exp_fh.close()

        #with open('correlated_data.pickle', 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        #  pickle.dump(correlated_data, f, pickle.HIGHEST_PROTOCOL)

        #correlated_data = pickle.load('correlated_data.pickle')

        #ranked_list_of_gene_names = map(lambda y: (y[0]), sorted(correlated_data, key=lambda x: np.sign(x[1].correlation) * x[1].pvalue))

        ranked_list_of_gene_names = map(lambda y: (y[0]), sorted(correlated_data, key=lambda x: x[1].pvalue))
        #sort gene_list_based on name
        data_in_order = map(lambda y: (y[0],y[1]), sorted(correlated_data, key=lambda x: x[0] ))
        
        gene_names_in_order = map(lambda y: y[0], data_in_order)
        p_values_in_order = map(lambda y: y[1].pvalue, data_in_order)
        correlations_in_order = map(lambda y: y[1].correlation, data_in_order)

        for ind in range(len(gene_names_in_order)):
          #gene_names_in_order[ind]
          #print ind
          #data_array[j][k][ind][0]
          data_array[j][k][ind][0] = ind #gene_names_in_order[ind]
          data_array[j][k][ind][i+1+len(files)] = str(p_values_in_order[ind])
          #FDR-adjusted p-value
          #print p_values_in_order
          #print ranked_list_of_gene_names.index(gene_names_in_order[ind])
          data_array[j][k][ind][i+1] = str(correlations_in_order[ind])
          data_array[j][k][ind][i+1+2*len(files)] = str( p_values_in_order[ind] * len(gene_names_in_order) / (ranked_list_of_gene_names.index(gene_names_in_order[ind]) + 1))        

  print("Writing Files")
  for j in range(3):
    for k in range(len(regions_of_interest)):
      label = MRI_data_labels[j] + "." + regions_of_interest[k][0]
      with open(config.outputCSVFolder + label + ".gene_list.csv",'w') as f:
        data = data_array[j][k]
        f.write(",correlation,,,,,,raw,,,,,,adjusted,,,,,,raw_meta_p,adjusted_meta_p\n")
        f.write("ID," + (",").join(brain_ids) + "," + ",".join(brain_ids)+ "\n")
        for gene_entry in data:
          gene_name = gene_names_in_order[int(gene_entry[0])]
          f.write(gene_name + "," + ",".join(map(str,gene_entry[1:])) + "," + str(fisher_p(gene_entry[7:12])) + "," + str(fisher_p(gene_entry[13:18])) + "\n")

  t2 = time.clock()

  print "Total execution time was",t2-t1