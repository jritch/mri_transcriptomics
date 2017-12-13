import config, analysis, os, csv, sys, inspect
import numpy as np

from ontology import *

def get_single_gene_data(gene_exp_fh,gene_name,indices=None):
  '''
  Takes a file handle and an optional list of indices.

  '''
  gene_name = "\"" + gene_name + "\""
  gene_exp_fh.readline()
  line = gene_exp_fh.readline()
  flag = 0
  while line:
    entries = line.strip().split('\t') 
    numerical_entries = map(float,entries[1:])
    ID = entries[0]
    if ID == gene_name or ID == gene_name.strip("\""):
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

def main():
  if len(sys.argv) > 1:
    gene_name = sys.argv[1]
  else:
    gene_name =  "TREML1"
    
  use_HCP_and_MNI = False

  if not use_HCP_and_MNI:
    MRI_data_labels = ["T1","T2","T1T2Ratio"]
    MRI_dimension = 2 # 0: T1, 1: T2, 2: ratio
  else:
    MRI_data_labels = ["HCP"]
    MRI_dimension = 0

  MRI_data_label = MRI_data_labels[MRI_dimension]
  print("Using: " + MRI_data_label)
  print("Gene:" + gene_name)
  '''
  regionIDs = [4219] 
  region_name = "limbic_lobe"
  to_exclude_list = None
  '''

  regionIDs = [4008] #cortex
  region_name = "cortex_excluding_piriform_hippocampus"
  to_exclude_list = [4249, 10142]
  '''
  regionIDs = [4005] #full brain
  region_name = "whole_brain"
  to_exclude_list = None
  '''


  files = config.get_expression_filenames(use_MNI=use_HCP_and_MNI)

  brain_ids = [f.split(".")[0] for f in files]

  outputFolder = os.path.join(config.scriptLocation, "results", "single_gene_data")
  if not os.path.exists(outputFolder):
      os.mkdir(outputFolder)

  if not use_HCP_and_MNI:
    filenames = [os.path.join(outputFolder, f + "." + gene_name + "." + region_name + ".MRI(xyz).expression."+ MRI_data_label +".csv") for f in brain_ids]
  else:
    filenames = [os.path.join(outputFolder, f + "." + gene_name + "." + region_name + ".MNI(xyz).expression."+ MRI_data_label +".csv") for f in brain_ids]
  
  o = Ontology(os.path.join(config.scriptLocation, "data",  "Ontology.csv"))
  
  cortex_divisions = [str(x) for x in o.hierarchy.get(4008)]
  
  num_brains = len(brain_ids)
  for i in range(num_brains):

        gene_exp_fh = open(os.path.join(config.processedOutputLocation,files[i]))
        print("Loading: " + os.path.join(config.processedOutputLocation,files[i]))
        coords,coord_to_region_map = analysis.get_coords_and_region_ids_from_gene_exp_data(gene_exp_fh)
        indices = []
        for j in range(len(regionIDs)):
          indices += analysis.get_flat_coords_from_region_id(regionIDs[j],coords,coord_to_region_map,o,to_exclude_list=to_exclude_list)
        assert(len(set(indices)) == len(indices))
        single_gene_data = np.array(get_single_gene_data(gene_exp_fh,gene_name,indices))
        gene_exp_fh.close()

        mri_data = analysis.load_nifti_data(brain_ids[i], use_HCP_and_MNI)[MRI_dimension]
        use_voxel_avg = False #True
        if use_voxel_avg:
            mri_data = voxel_averaged_mri.voxel_average(mri_data)

        flat_mri_data = np.array(analysis.flatten_mri_data(mri_data,coords, use_HCP_and_MNI))
        region_specific_flat_mri_data = [flat_mri_data[j] for j in indices]
        results = np.zeros((5,len(indices)))
        results[1,:] = region_specific_flat_mri_data
        results[2,:] = single_gene_data
        coord_subset = [coords[j] for j in indices]

        with open(filenames[i], "w") as f:
          f.write("\"(x,y,z)\","  + "MRI_Intensity" +"," + gene_name  + ",regionID,region_name,cortical_division\n")
          for j in range(len(indices)):
            str_coords = [str(coord) for coord in coord_subset[j]]
            coord_string = "(" + ",".join(str_coords) + ")" 
            results[3,j] = coord_to_region_map[coord_string]
            current_region_ID = int(results[3,j])
            #get the main cortical lobe
            enclosing_regions = [str(x) for x in o.get_enclosing_regions(current_region_ID)]
            cortex_subdivision = set(enclosing_regions).intersection(cortex_divisions)
            if (len(cortex_subdivision) != 0):
                cortex_subdivision = next(iter(cortex_subdivision))
                cortex_subdivision = o.names[int(cortex_subdivision)]
            else:
                cortex_subdivision = ""
            f.write("\"" + coord_string + "\","  + str(results[1,j]) + "," + str(results[2,j]) + "," + str(current_region_ID) + ",\""+ o.names[current_region_ID] + "\", " + cortex_subdivision + "\n")


if __name__ == '__main__':
    main()
    print("Done")
