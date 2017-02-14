import config
import analysis
import os, csv, sys, inspect
import numpy as np

def get_single_gene_data(gene_exp_fh,gene_name,indices=None):
  '''
  Takes a file hangle and an optional list of indices.

  '''

  gene_name = "\"" + gene_name + "\""

  gene_exp_fh.readline()

  line = gene_exp_fh.readline()
  
  flag = 0
  while line:
    entries = line.strip().split('\t') 
    numerical_entries = map(float,entries[1:])
    ID = entries[0]
    #print ID

    if ID == gene_name:
       flag = 1
       break

    line = gene_exp_fh.readline()

  # if flag is still 0, then that gene symbol was not in the list
  if flag == 0:
    return None

  if not indices:
    return numerical_entries 
  else:
    return [numerical_entries[i] for i in indices]

  return result



def main():
  gene_name =  "ZNF362"
  MRI_dimension = 2 # 0: T1, 1: T2, 2: ratio
  regionID = 4008 #cortex
  region_name = "cortex"

  files = ["10021.matrix.regionID.MRI(xyz).29131 x 893.txt",
           "12876.matrix.regionID.MRI(xyz).29131 x 363.txt",
           "14380.matrix.regionID.MRI(xyz).29131 x 529.txt",
           "15496.matrix.regionID.MRI(xyz).29131 x 470.txt",
           "15697.matrix.regionID.MRI(xyz).29131 x 501.txt",
           "9861.matrix.regionID.MRI(xyz).29131 x 946.txt"]
  
  brain_ids = [f.split(".")[0] for f in files]

  #header = "\"(x,y,z)\"," + ",".join([ '"' + f + "_MRI\",\"" + f + "_" + gene_name + '"'  for f in brain_ids])
  header = "\"(x,y,z)\"," + ",".join([ '"' + f + "_MRI\",\"" + f + "_" + gene_name + '"'  for f in brain_ids])
  baseProjectFolder = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) + "/" 
  filenames = [baseProjectFolder + "single_gene_data/" + f + "." + gene_name + "." + region_name + ".MRI(xyz).expression.csv" for f in brain_ids]
  print filenames
  

  o = analysis.Ontology(config.ontologyFolder + "Ontology.csv")

  num_brains = len(brain_ids)
  for i in range(num_brains):

        gene_exp_fh = open(os.path.join(config.MRIFolder,files[i]))
        coords,coord_to_region_map = analysis.get_coords_and_region_ids_from_gene_exp_data(gene_exp_fh)
        indices = analysis.get_flat_coords_from_region_id(regionID,coords,coord_to_region_map,o)
        single_gene_data = np.array(get_single_gene_data(gene_exp_fh,gene_name,indices))
        gene_exp_fh.close()


        mri_data = analysis.load_nifti_data(config.baseAllenFolder + "normalized_microarray_donor" + brain_ids[i])[MRI_dimension]
        flat_mri_data = np.array(analysis.flatten_mri_data(mri_data,coords))
        region_specific_flat_mri_data = [flat_mri_data[j] for j in indices]

        results = np.zeros((3,len(indices)))

        results[1,:] = region_specific_flat_mri_data
        results[2,:] = single_gene_data

        print os.getcwd()

        coord_subset = [coords[j] for j in indices]

        with open(filenames[i], "w") as f:
          f.write("\"(x,y,z)\",MRI_Intensity," + gene_name  + "\"")
          for j in range(len(indices)):

            str_coords = [str(coord) for coord in coord_subset[j]]
            coord_string = "(" + ",".join(str_coords) + ")" 
            f.write("\"" + coord_string + "\","  + str(results[1,j]) + "," + str(results[2,j]) + "\n")


if __name__ == '__main__':
    main()
