import pandas
import numpy as np
import pdb
import time
import pickle
import config
import os 
import sys
from numpy.testing import assert_approx_equal

num_runs = 0

def sort_columns(dataframe):
  """
  Helper function to sort the columns of a dataframe
  """
  return dataframe.reindex_axis(sorted(dataframe.columns), axis=1)

def average_group(x):
  global num_runs
  num_runs += 1
  #import pdb; pdb.set_trace()
  #sys.exit()

  if num_runs % 5000 == 0:
    print "GROUP NUMBER:", num_runs
  return np.mean(x[range(1,x.columns[-1]+1)])

def average_probes(probes,exp_data):
  merged_table = pandas.merge(probes,exp_data,left_on="probe_id",right_on=0)#.iloc[0:1]

  gb = merged_table.groupby("gene_symbol")
  print ("THERE ARE " + str(len(gb)) + " GROUPS")
  start_time = time.time()
  #import pdb; pdb.set_trace()
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
  return 
  #pdb.set_trace()

def validate(data):
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
  
  return 


def spot_check_values(filename, output_filename):
  my_df = pandas.read_csv(filename, sep="\t")
  validation_df = pandas.read_csv(output_filename, sep="\t")
  my_df_value = my_df[my_df['Unnamed: 0'] =="SLC17A7"]["RegionID:4284|LOC:(80,85,62)"]  
  validation_df_value =  validation_df[validation_df['Unnamed: 0'] =="SLC17A7"]["RegionID:4284|LOC:(80,85,62)"]  
  assert(my_df_value == validation_df_value)
  sys.exit()


def check_all_values(filename, output_filename):
  my_df = pandas.read_csv(filename, sep="\t")
  my_df.sort_values(by="Unnamed: 0",inplace="True")
  my_df = my_df.reindex_axis(sorted(my_df.columns), axis=1)

  validation_df = pandas.read_csv(output_filename, sep="\t")
  validation_df.sort_values(by="Unnamed: 0",inplace="True")
  validation_df = validation_df.reindex_axis(sorted(validation_df.columns), axis=1)

  for i in range(my_df.shape[0]):
   for j in range(my_df.shape[1]-1):
     assert_approx_equal(validation_df.iloc[i,j], my_df.iloc[i,j],significant=7)
   if i % 100:
    print(i)


def process_file(probe_filename,annotation_filename,expression_data_filename,output_filename):
  """
  Gathers data 
  """
  
  global num_runs
  num_runs = 0

  probes = pandas.read_csv(probe_filename)
  annotations = pandas.read_csv(annotation_filename)
  exp_data = pandas.read_csv(expression_data_filename, header=None)

  data  = average_probes(probes,exp_data)
  header = annotations.apply(get_header,axis=1) 
  #pdb.set_trace()
  col_map = {}
  
  for i in range (0,len(header)):
    col_map[i+1] = header.loc[i][0]
  col_map["gene_symbol"] = ""

  data.rename(columns=col_map, inplace=True)
  data = data.reset_index() 

  #validate(data)

  # Get rid of gene_symbol column name so that columns will match
  data.columns.values[0] = None
  
  #print small csv for test purposes
  data.to_csv(output_filename,sep="\t",index=False)

def main():
  
  OUTPUT_FOLDER = "python_processed_expression_data"

  brain_numbers = [x.split(".")[0] for x in config.expression_filenames]

  for i in range(len(brain_numbers)):

    print ("PROCESSING BRAIN {} of 6, BRAIN ID = {}".format(i+1, brain_numbers[i]))

    probe_filename = os.path.join(config.microarrayFolder, "normalized_microarray_donor" + str(brain_numbers[i]), "Probes.csv")
    annotation_filename = os.path.join(config.microarrayFolder, "normalized_microarray_donor" + str(brain_numbers[i]),"SampleAnnot.csv")
    expression_data_filename = os.path.join(config.microarrayFolder, "normalized_microarray_donor" + str(brain_numbers[i]),"MicroarrayExpression.csv")
    output_filename = os.path.join(OUTPUT_FOLDER, config.expression_filenames[i])

    process_file(probe_filename,annotation_filename,expression_data_filename,output_filename)


if __name__ == '__main__':
  #spot_check_values("~/test_full.csv", "~/Downloads/AllenHBAProcessedExpressionWithBrainID/9861.matrix.regionID.MRI(xyz).29131 x 946.txt")
  check_all_values("~/test_full.csv", "~/Downloads/AllenHBAProcessedExpressionWithBrainID/9861.matrix.regionID.MRI(xyz).29131 x 946.txt")
  #main()
  #check_probes()

