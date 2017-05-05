import pandas
import numpy as np
import pdb
import time

num_runs = 0

def average_group(x):
  global num_runs
  num_runs += 1

  if num_runs % 1000 == 0:
    print "GROUP NUMBER:", num_runs

  #[x.gene_symbol.iloc[0]] + 
  return np.mean(x[range(1,894)])

def average_probes(probes,exp_data):
  merged_table = pandas.merge(probes,exp_data,left_on="probe_id",right_on=0)

  gb = merged_table.groupby("gene_symbol")
  print ("THERE ARE " + str(len(gb)) + " GROUPS")

  start_time = time.time()
    
  averaged_table = gb.apply(average_group)

  print "TOOK", time.time() - start_time, "SECONDS"

  return averaged_table


def get_header(x):
  """ This processes the SampleAnnot.csv file to get the header with RegionID and Location"""
  LOC_tuple = lambda x : "(" + str(x.mri_voxel_x)+ "," + str(x.mri_voxel_y)+ "," + str(x.mri_voxel_z) + ")"
  return pandas.Series("RegionID:"+ str(x.structure_id) + "|LOC:" +str(LOC_tuple(x)))
  
def main():
  probes = pandas.read_csv("~/Documents/microarray_data/normalized_microarray_donor10021/Probes.csv")
  annotations = pandas.read_csv("~/Documents/microarray_data/normalized_microarray_donor10021/SampleAnnot.csv")
  exp_data = pandas.read_csv("~/Documents/microarray_data/normalized_microarray_donor10021/MicroarrayExpression.csv", header=None)
  data  = average_probes(probes,exp_data)
  header = annotations.apply(get_header,axis=1)

  col_map = {}
  
  for i in range (0,893):
    col_map[i] = header.loc[i][0]
  col_map["gene_symbol"] = ""

  data.rename(columns=col_map, inplace=True)

  data.to_csv("~/munged.csv",sep="\t")

if __name__ == '__main__':
  main()



