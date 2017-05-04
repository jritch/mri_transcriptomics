import pandas
import numpy as np
import pdb
import time



num_runs = 0

def average_group(x):
  global num_runs
  num_runs += 1
  pdb.set_trace()
  if num_runs % 1000 == 0:
    print "GROUP NUMBER:", num_runs
    print time.time() - start_time, "SECONDS SO FAR"
  y = np.zeros([893])
  for i in range(1,893):
    y[i-1]=np.mean(x[i]);
  columns = ["probe_name","gene_id","gene_symbol","gene_name","entrez_id","chromosome",0]
  defaults_list = [x[col] for col in columns]
  y_val_list = [item for item in y]
  return pandas.Series(defaults_list +  y_val_list)

def average_probes(probes,exp_data):
  merged_table = pandas.merge(probes,exp_data,left_on="probe_id",right_on=0)
  #d.to_csv("~/merged.csv")

  # to make these easier to deal with I will cut the merged table down to 100 probes
  # before 

  merged_table_truncated = merged_table.loc

  gb = merged_table_truncated.groupby("gene_name")
  print "GROUPED"
  print ("THERE ARE " + str(len(gb)) + " GROUPS")

  start_time = time.time()
    
  averaged_table = gb.apply(average_group)
  print "APPLIED"

  print "TOOK", time.time() - start_time, "SECONDS"

  f.to_csv("~/munged.csv")
  print "DONE"


def get_header2(x):
  LOC_tuple = lambda x : "(" + str(x.mri_voxel_x)+ "," + str(x.mri_voxel_y)+ "," + str(x.mri_voxel_z) + ")"
  return pandas.Series("regionID="+ str(x.structure_id) + "|LOC=" +str(LOC_tuple(x)))


def get_header(annotations):
  """ This processes the SampleAnnot.csv file to get the header with RegionID and Location"""
  #pdb.set_trace()
  a = annotations.apply(get_header2,axis=1)
  print a 
  print a[0:2]
  
  return
  
def main():
 # probes = pandas.read_csv("~/Documents/microarray_data/normalized_microarray_donor10021/Probes.csv")
  annotations = pandas.read_csv("~/Documents/microarray_data/normalized_microarray_donor10021/SampleAnnot.csv")
 # pdb.set_trace()
  #exp_data = pandas.read_csv("~/Documents/microarray_data/normalized_microarray_donor10021/MicroarrayExpression.csv", header=None)
  #data  = average_probes(probes,exp_data)
  header = get_header(annotations)
  # CONCATENATE HEADER AND DATA

if __name__ == '__main__':
  main()



