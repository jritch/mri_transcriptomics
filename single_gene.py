import config, analysis, os, csv, sys, inspect
import numpy as np

import voxel_averaged_mri

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
  gene_name =  "FLJ23867"
  MRI_dimension = 2 # 0: T1, 1: T2, 2: ratio

  regionID = 4008 #cortex
  region_name = "cortex"
  to_exclude = 4219

  '''
  regionID = 4008 #cortex
  region_name = "cortex_excluding_limbic_lobe"
	to_exclude = 4219

	'''

  files = ["10021.matrix.regionID.MRI(xyz).29131 x 893.txt",
           "12876.matrix.regionID.MRI(xyz).29131 x 363.txt",
           "14380.matrix.regionID.MRI(xyz).29131 x 529.txt",
           "15496.matrix.regionID.MRI(xyz).29131 x 470.txt",
           "15697.matrix.regionID.MRI(xyz).29131 x 501.txt",
           "9861.matrix.regionID.MRI(xyz).29131 x 946.txt"]
  
  brain_ids = [f.split(".")[0] for f in files]

  header = "\"(x,y,z)\"," + ",".join([ '"' + f + "_MRI\",\"" + f + "_" + gene_name + '"'  for f in brain_ids])
  baseProjectFolder = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) + "/" 
  if not os.path.exists(baseProjectFolder + "single_gene_data_avg/"):
      os.mkdir(baseProjectFolder + "single_gene_data_avg/")

  filenames = [baseProjectFolder + "single_gene_data_avg/" + f + "." + gene_name + "." + region_name + ".MRI(xyz).expression.csv" for f in brain_ids]
  
  o = analysis.Ontology(config.ontologyFolder + "Ontology.csv")
  
  cortex_divisions = [str(x) for x in o.hierarchy.get(4008)]
  
  num_brains = len(brain_ids)
  for i in range(num_brains):

        gene_exp_fh = open(os.path.join(config.expressionFolder,files[i]))
        coords,coord_to_region_map = analysis.get_coords_and_region_ids_from_gene_exp_data(gene_exp_fh)
        indices = analysis.get_flat_coords_from_region_id(regionID,coords,coord_to_region_map,o,to_exclude=to_exclude)
        single_gene_data = np.array(get_single_gene_data(gene_exp_fh,gene_name,indices))
        gene_exp_fh.close()

        mri_data = analysis.load_nifti_data(config.baseAllenFolder + "normalized_microarray_donor" + brain_ids[i])[MRI_dimension]

        use_voxel_avg = True
        if use_voxel_avg:
            mri_data = voxel_averaged_mri.voxel_average(mri_data)

        flat_mri_data = np.array(analysis.flatten_mri_data(mri_data,coords))
        region_specific_flat_mri_data = [flat_mri_data[j] for j in indices]
        results = np.zeros((5,len(indices)))
        results[1,:] = region_specific_flat_mri_data
        results[2,:] = single_gene_data
        coord_subset = [coords[j] for j in indices]

        with open(filenames[i], "w") as f:
          f.write("\"(x,y,z)\",MRI_Intensity," + gene_name  + ",regionID,region_name,cortical_division\n")
          for j in range(len(indices)):
            str_coords = [str(coord) for coord in coord_subset[j]]
            coord_string = "(" + ",".join(str_coords) + ")" 
            results[3,j] = coord_to_region_map[coord_string]
            current_region_ID = int(results[3,j])
            #get the main cortical lobe
            enclosing_regions = [str(x) for x in o.get_enclosing_regions(current_region_ID)]
            cortex_subdivision = set(enclosing_regions).intersection(cortex_divisions)
            cortex_subdivision= next(iter(cortex_subdivision))
            cortex_subdivision = o.names[int(cortex_subdivision)]
            f.write("\"" + coord_string + "\","  + str(results[1,j]) + "," + str(results[2,j]) + "," + str(current_region_ID) + ",\""+ o.names[current_region_ID] + "\", " + cortex_subdivision + "\n")



if __name__ == '__main__':
    main()
