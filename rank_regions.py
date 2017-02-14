import config
import analysis
import os, csv, sys, inspect
import numpy as np

def main():

  MRI_dimension = 2 # 0: T1, 1: T2, 2: ratio
  
  files = ["10021.matrix.regionID.MRI(xyz).29131 x 893.txt",
           "12876.matrix.regionID.MRI(xyz).29131 x 363.txt",
           "14380.matrix.regionID.MRI(xyz).29131 x 529.txt",
           "15496.matrix.regionID.MRI(xyz).29131 x 470.txt",
           "15697.matrix.regionID.MRI(xyz).29131 x 501.txt",
           "9861.matrix.regionID.MRI(xyz).29131 x 946.txt"]
  
  brain_ids = [f.split(".")[0] for f in files]

  baseProjectFolder = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) + "/" 
  filenames = [baseProjectFolder + "region_rankings/" + f + ".region.avg_MRI" + ".csv" for f in brain_ids]
  
  o = analysis.Ontology(config.ontologyFolder + "Ontology.csv")

  print len(o.get_all_regionIDs(4005))
  sys.exit()

  num_brains = len(brain_ids)
  for i in range(num_brains):

        gene_exp_fh = open(os.path.join(config.MRIFolder,files[i]))
        coords,coord_to_region_map = analysis.get_coords_and_region_ids_from_gene_exp_data(gene_exp_fh)
        gene_exp_fh.close()


        mri_data = analysis.load_nifti_data(config.baseAllenFolder + "normalized_microarray_donor" + brain_ids[i])[MRI_dimension]
        flat_mri_data = analysis.flatten_mri_data(mri_data,coords)
        
        ranking = analysis.rank_regions_by_intensity(coords,coord_to_region_map,o,flat_mri_data)

        with open(filenames[i], "w") as f:
          f.write("\"region\",\"MRI Intensity\"\n")
          for item in ranking:
              f.write("\"" + str(item[0]) + "\","  + str(item[1])  + "\n")


if __name__ == '__main__':
    main()