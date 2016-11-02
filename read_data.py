import nibabel
import numpy
import os
import csv
#from subprocess import call

import time
import config
import inspect, os

from scipy.stats import pearsonr
from scipy.stats import spearmanr


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
  for k,v in coord_to_region_map.iteritems():
    if v == regionID:
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

  for ID,name in ontology.names.iteritems():
    flag = 0
    if ID == 9222:
      flag = 1
    intensity = 0
    coords_list = get_coords_from_region_id(ID,coords,coord_to_region_map,ontology)
    if coords_list != []:
      if flag:
        print coords_list
      
      print name,ID,coords_list
      for coord in coords_list:
        index = coords_list.index(coord)
        intensity += flat_mri_data[index]
        total +=1
        if not flat_mri_data[index] in [3,4,5,6] :
          nonzero += 1
      intensity = intensity / float(len(coords_list))
    rank.append((name, intensity))

  print nonzero * 1. / total * 100.
  return sorted(rank,key=lambda x: x[1],reverse=True)

def convert_coords_to_flat_indices(coords_list):
  pass  

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
  

  t1 = time.clock()

  MRI_data = load_nifti_data(config.baseAllenFolder + "normalized_microarray_donor9861")

  o = Ontology(baseProjectFolder + "/data/Ontology.csv")

  print 'Getting coordinates from gene expression file.'
  
  gene_exp_fh = open(config.baseAllenFolder +  "normalized_microarray_donor9861/9861.matrix.regionID.MRI(xyz).29131 x 946.txt")
  coords,coord_to_region_map = get_coords_and_region_ids_from_gene_exp_data(gene_exp_fh)
  gene_exp_fh.close()

  flat_t1t2_ratio_data = flatten_mri_data(MRI_data[2],coords)

  import pprint
  ratioListFilename = baseProjectFolder + 'results/regions_ranked_by_t1overt2_ratio.txt'
  if (os.path.isfile(ratioListFilename)): os.unlink(ratioListFilename)
    
  with open(ratioListFilename,'w') as f:
    f.write(pprint.pformat(rank_regions_by_intensity(coords,coord_to_region_map,o,flat_t1t2_ratio_data)))



  for measure in MRI_data:
    gene_exp_fh = open(config.baseAllenFolder + "normalized_microarray_donor9861/9861.matrix.MRI(xyz).29131 x 946.txt")

    print 'Flattening MRI data.'
    flat_mri_data = flatten_mri_data(measure,coords)
    
    print 'Correlating MRI and gene expression data.'
    print correlate_MRI_and_gene_exp_data(flat_mri_data,gene_exp_fh)
    
    gene_exp_fh.close()

  t2 = time.clock()

  print "Total execution time was",t2-t1