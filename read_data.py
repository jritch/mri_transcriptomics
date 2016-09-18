import nibabel
import numpy
import os
import csv
#from subprocess import call

import time

from scipy.stats import pearsonr
from scipy.stats import spearmanr

R_VALUE_FUNCTION = spearmanr

def load_ontology():
  '''
  Reads the ontology file and returns an (ID,name) dictionary 
  and a (parentID,[childID]) dictionary that defines the ontology / hierarchy.
  '''
  names = dict()
  hierarchy = dict()
  with open('C:/Users/Jacob/large_thesis_files/Ontology.csv') as f:
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

#  print names
  print hierarchy
  return names, hierarchy

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
  header_line = gene_exp_fh.readline().strip() 
  coord_strings = header_line.split('\t')
  evaluate_strings = lambda x : eval(x.strip('"').split('|')[1].split(':')[1])
  coord_to_region_map = {}
  for string in coord_strings:
      components = string.strip('"').split('|')
      coord_to_region_map[components[1].split(':')[1]] = eval(components[0].split(':')[1])

  coords = map(evaluate_strings,coord_strings)

  return coords,coord_to_region_map

def correlate_MRI_and_gene_exp_data(flat_mri_data,gene_exp_fh):
  '''
  Takes the name of a directory containing a txt file with the gene_expression data and two NIFTI files.

  Reads all three files.

  Correlates(using the pearson correlation coefficient) the expression data for each gene with the one measure of MRI intensity.

  Returns a list of the top ten genes by (genes with p value > 0.05 are not considered)

  Using the file handle and processing one line at a time allows us to avoid loading the entire gene expression text file into memory.

  '''
  t1 = time.clock()

  f = gene_exp_fh

  # Throw away header line
  f.readline()

  line = f.readline()
  IDs = list()
  correlations = list()

  while line:
    entries = line.strip().split('\t') 
    numerical_entries = map(float,entries[1:])
    ID = entries[0]
    IDs.append(ID)
    correlation = R_VALUE_FUNCTION(numerical_entries,flat_mri_data)

    # look at genes with p-value less than 0.05
    if correlation[1] <= 0.05:
      #print ID, correlation
      correlations.append([ID,correlation])

    line = f.readline()

  f.close();
  #sort the significantly correlated genes based on correlation and return
  result = sorted(correlations, key=lambda entry: abs(entry[1][0]),reverse=True)[0:9]
  t2 = time.clock()

  print "Correlation execution time was", t2-t1

  return result

def visualize():
  pass

if __name__ == '__main__':
  
  #load_ontology()

  gene_exp_fh = open("C:/Users/Jacob/large_thesis_files/AllenHBAProcessedExpressionWithBrainID/9861.matrix.regionID.MRI(xyz).29131 x 946.txt")
  coords,dictionary = get_coords_and_region_ids_from_gene_exp_data(gene_exp_fh)
  print coords
  print dictionary
  exit()

  t1 = time.clock()

  MRI_data = load_nifti_data("C:/Users/Jacob/Downloads/AllenHBAProcessedExpressionAndMRIs/normalized_microarray_donor9861")
  
  print 'Getting coordinates from gene expression file.'
  
  gene_exp_fh = open("C:/Users/Jacob/Downloads/AllenHBAProcessedExpressionAndMRIs/normalized_microarray_donor9861/9861.matrix.MRI(xyz).29131 x 946.txt")
  coords.dictionary = get_coords_and_region_ids_from_gene_exp_data(gene_exp_fh)
  gene_exp_fh.close()

  for measure in MRI_data:
    gene_exp_fh = open("C:/Users/Jacob/Downloads/AllenHBAProcessedExpressionAndMRIs/normalized_microarray_donor9861/9861.matrix.MRI(xyz).29131 x 946.txt")

    print 'Flattening MRI data.'
    flat_mri_data = flatten_mri_data(measure,coords)
    
    print 'Correlating MRI and gene expression data.'
    print correlate_MRI_and_gene_exp_data(flat_mri_data,gene_exp_fh)
    
    gene_exp_fh.close()

  t2 = time.clock()

  print "Total execution time was",t2-t1