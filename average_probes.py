import pandas
import numpy as np
import pdb
import time
import pickle

num_runs = 0

def sort_columns(dataframe):
  """
  Helper function to sort the columns of a dataframe
  """
  return dataframe.reindex_axis(sorted(dataframe.columns), axis=1)

def average_group(x):
  global num_runs
  num_runs += 1

  if num_runs % 1000 == 0:
    print "GROUP NUMBER:", num_runs
  return np.mean(x[range(1,894)])

def average_probes(probes,exp_data):
  merged_table = pandas.merge(probes,exp_data,left_on="probe_id",right_on=0)#.iloc[0:1]

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
  


def check_probes():
  """
  This function will check the probes where the probe symbol is a gene_symbol, to see what happens to them
  """
  # limit probes to the agilent microarray gene_symbols
  probes = pandas.read_csv("~/Documents/microarray_data/normalized_microarray_donor10021/Probes.csv")
  of_interest = probes.gene_symbol.str.match(r"A_.*_.*")
  exp_data = pandas.read_csv("~/Documents/microarray_data/normalized_microarray_donor10021/MicroarrayExpression.csv", header=None)
  avg =  average_probes(probes,exp_data)
  avg = avg.reset_index()
  of_interest = avg.gene_symbol.str.match(r"A_.*_.*")
  print avg[of_interest == True]

  #pdb.set_trace()

def main():
  probes = pandas.read_csv("~/Documents/microarray_data/normalized_microarray_donor10021/Probes.csv")
  annotations = pandas.read_csv("~/Documents/microarray_data/normalized_microarray_donor10021/SampleAnnot.csv")
  exp_data = pandas.read_csv("~/Documents/microarray_data/normalized_microarray_donor10021/MicroarrayExpression.csv", header=None)
  data  = average_probes(probes,exp_data)
  header = annotations.apply(get_header,axis=1) 

  col_map = {}
  
  for i in range (0,893):
    col_map[i+1] = header.loc[i][0]
  col_map["gene_symbol"] = ""

  data.rename(columns=col_map, inplace=True)
  data = data.reset_index() 

  # Now test this against Leon's earlier solution

  validation_data = pandas.read_csv("~/Downloads/AllenHBAProcessedExpressionWithBrainID/10021.matrix.regionID.MRI(xyz).29131 x 893.txt",sep="\t")

  data.sort_values(by="gene_symbol",inplace=True)
  
  validation_data.sort_values(by="Unnamed: 0",inplace=True)

  # check that all the columns are the same
  # the one difference is the "gene_symbol" column which is
  # unnamed in the validation_data dataframe
  a = [col in validation_data.columns for col in data.columns]
  print a.count(True)
  print a.count(False)

  # check that all the gene_symbols are the same
  a = [sym in validation_data['Unnamed: 0'].values for sym in data.gene_symbol]
  print a.count(True)
  print a.count(False)

  # Get rid of gene_symbol column name so that columns will match
  data.columns.values[0] = None
  
  #print small csv for test purposes
  data.loc[0:10].to_csv("~/test.csv",sep="\t",index=False)

  #print full csv (takes a long time)
  #data.to_csv("~/test.csv",sep="\t",index=False)

if __name__ == '__main__':
  main()
  #check_probes()

