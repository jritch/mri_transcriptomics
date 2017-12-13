import pandas
import numpy as np
import pdb
import time
import pickle
import config
import os
from openpyxl import load_workbook
import sys
import xlrd
from numpy.testing import assert_approx_equal
from os import listdir
from os.path import isfile, join

num_runs = 0
probe_processed_filename = os.path.join(config.microarrayFolder, "Probes.enhanced.csv")
sums_processed_filename = os.path.join(config.processedOutputLocation , "PA_sums.csv")


def sort_columns(dataframe):
  """
  Helper function to sort the columns of a dataframe
  """
  return dataframe.reindex_axis(sorted(dataframe.columns), axis=1)

def average_group(x, do_mean):
  global num_runs
  num_runs += 1

  if num_runs % 5000 == 0:
    print "GROUP NUMBER:", num_runs
  if do_mean:
    return np.mean(x[range(1, x.columns[-1] + 1)])
  else:
    return np.sum(x[range(1, x.columns[-1] + 1)])

def average_probes(probes, exp_data, probe_strategy):
  merged_table = pandas.merge(probes, exp_data, left_on="probe_id", right_on=0)  # .iloc[0:1]
  merged_table = merged_table.loc[merged_table['gene_symbol'] != 'na']

  print("Probes before filtering: " + str(merged_table.shape[0]))
  if probe_strategy is "all":
    None
  elif probe_strategy is "scale":
    #filter for nonNA
    merged_table = merged_table[merged_table.is_qc_pass == True] 

    #multiply and add b
    merged_table.m = merged_table.m.astype(float)
    merged_table.b = merged_table.b.astype(float)
    merged_table = merged_table[np.logical_not(pandas.isnull(merged_table['m']))]

    g=merged_table.columns.to_series().groupby(merged_table.dtypes).groups

    #remove m and b columns
    g[np.dtype('float64')] = g[np.dtype('float64')].difference(('m', 'b'))
    #only do it on select columns
    merged_table[g[np.dtype('float64')]] = merged_table[g[np.dtype('float64')]].apply(lambda x: merged_table.m * x + merged_table.b)
    #drop m and b columns
    merged_table = merged_table.drop(['m', 'b'], axis=1)
  elif probe_strategy is "highSTD":
    merged_table = merged_table[merged_table.is_most_variable == True]
  elif probe_strategy is "passQC":
    merged_table = merged_table[merged_table.is_qc_pass != False] #count NA values as passing
  print("Probes after filtering: " + str(merged_table.shape[0]))

  gb = merged_table.groupby("gene_symbol")

  print ("THERE ARE " + str(len(gb)) + " GROUPS")
  start_time = time.time()
  averaged_table = gb.apply(average_group, do_mean = probe_strategy is not "scale")

  print "TOOK", time.time() - start_time, "SECONDS"

  return averaged_table


def get_header(x):
  """ This processes the SampleAnnot.csv file to get the header with RegionID and Location in native MRI coordinates"""
  LOC_tuple = lambda x : "(" + str(x.mri_voxel_x) + "," + str(x.mri_voxel_y) + "," + str(x.mri_voxel_z) + ")"
  return pandas.Series("RegionID:" + str(x.structure_id) + "|LOC:" + str(LOC_tuple(x)))

def get_MNI_header(x):
  """ This processes the SampleAnnot.csv file to get the header with RegionID and Location in MNI coordinates"""
  LOC_tuple = lambda x : "(" + str(x.mni_x) + "," + str(x.mni_y) + "," + str(x.mni_z) + ")"
  return pandas.Series("RegionID:" + str(x.structure_id) + "|LOC:" + str(LOC_tuple(x)))

def check_probes():
  """
  This function will check the probes where the probe symbol is a gene_symbol, to see what happens to them
  """
  # limit probes to the agilent microarray gene_symbols
  probes = pandas.read_csv("~/Documents/microarray_data/normalized_microarray_donor10021/Probes.csv")
  of_interest = probes.gene_symbol.str.match(r"A_.*_.*")
  exp_data = pandas.read_csv("~/Documents/microarray_data/normalized_microarray_donor10021/MicroarrayExpression.csv", header=None)
  avg = average_probes(probes, exp_data)
  avg = avg.reset_index()
  of_interest = avg.gene_symbol.str.match(r"A_.*_.*")
  print avg[of_interest == True]
  return
  # pdb.set_trace()

def validate(data):
    # Now test this against Leon's earlier solution

  validation_data = pandas.read_csv("~/Downloads/AllenHBAProcessedExpressionWithBrainID/10021.matrix.regionID.MRI(xyz).29131 x 893.tsv", sep="\t")

  data.sort_values(by="gene_symbol", inplace=True)

  validation_data.sort_values(by="Unnamed: 0", inplace=True)

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
  my_df_value = my_df[my_df['Unnamed: 0'] == "SLC17A7"]["RegionID:4284|LOC:(80,85,62)"]
  validation_df_value = validation_df[validation_df['Unnamed: 0'] == "SLC17A7"]["RegionID:4284|LOC:(80,85,62)"]
  assert(my_df_value == validation_df_value)


def check_all_values(filename, output_filename):
  my_df = pandas.read_csv(filename, sep="\t")
  my_df.sort_values(by="Unnamed: 0", inplace="True")
  my_df = my_df.reindex_axis(sorted(my_df.columns), axis=1)

  validation_df = pandas.read_csv(output_filename, sep="\t")
  validation_df.sort_values(by="Unnamed: 0", inplace="True")
  validation_df = validation_df.reindex_axis(sorted(validation_df.columns), axis=1)

  for i in range(my_df.shape[0]):
   for j in range(my_df.shape[1] - 1):
     assert_approx_equal(validation_df.iloc[i, j], my_df.iloc[i, j], significant=7)
   if i % 100:
    print(i)


def process_file(probe_filename, annotation_filename, expression_data_filename, output_filename,use_mni_coordinates=False, probe_strategy="all"):
  """
  Gathers data
  """

  global num_runs
  num_runs = 0

  probes = pandas.read_csv(probe_filename)
  annotations = pandas.read_csv(annotation_filename)
  exp_data = pandas.read_csv(expression_data_filename, header=None)#, nrows=1000)

  data = average_probes(probes, exp_data, probe_strategy)
  if use_mni_coordinates:
      header = annotations.apply(get_MNI_header,axis=1)
  else:
      header = annotations.apply(get_header, axis=1)
  # pdb.set_trace()
  col_map = {}

  for i in range (0, len(header)):
    col_map[i + 1] = header.loc[i][0]
  col_map["gene_symbol"] = ""

  data.rename(columns=col_map, inplace=True)
  data = data.reset_index()

  # validate(data)

  # Get rid of gene_symbol column name so that columns will match
  data.columns.values[0] = None

  data.to_csv(output_filename, sep="\t", index=False)
  

  
def write_probe_information():
  donorFolders = [f for f in listdir(config.microarrayFolder) if 'normalized_microarray_donor' in f]
  
  result_probe_info = None
  
  brain_sd_col_list = []
  for donorFolder in donorFolders:
    brain_number = donorFolder.split("_donor")[1]
    brain_number = "BrainSD_" + str(brain_number)
    brain_sd_col_list.append(brain_number)

    print ("PROCESSING SD for " + donorFolder)

    probe_filename = os.path.join(config.microarrayFolder, donorFolder, "Probes.csv")
    expression_data_filename = os.path.join(config.microarrayFolder, donorFolder, "MicroarrayExpression.csv")
    
    probes = pandas.read_csv(probe_filename)
    if result_probe_info is None: 
        result_probe_info = probes.copy()
    
    exp_data = pandas.read_csv(expression_data_filename, header=None) #, nrows=200)
    merged_table = pandas.merge(probes[['probe_id', 'probe_name']], exp_data, left_on="probe_id", right_on=0) 

    merged_table = merged_table.drop(0, axis=1)
    merged_table = merged_table.drop('probe_id', axis=1)

    sdTable = merged_table.std(axis=1)

    sdTable = pandas.concat([merged_table[['probe_name']], sdTable], axis=1)
    sdTable.rename(columns={0: brain_number }, inplace=True)
    result_probe_info = pandas.merge(result_probe_info, sdTable, on='probe_name') 
  
  #compute average SD
  result_probe_info = result_probe_info.assign(mean_sd = result_probe_info[brain_sd_col_list].mean(axis=1))
  
  #find probe with highest SD for a gene
  indexes_with_max_sd = result_probe_info.groupby(['gene_symbol'])['mean_sd'].transform(max) == result_probe_info['mean_sd']
  result_probe_info['is_most_variable'] = indexes_with_max_sd 

  #merge in xlsx data
  

  print("Loading Miller et al. probe information (xlsx)")
  
  qc_table = pandas.read_excel(os.path.join(config.scriptLocation, "data",  "Miller et al. doi.org_10.1186_1471-2164-15-154 12864_2013_7016_MOESM8_ESM.xlsx"))
  qc_table = qc_table.iloc[:, [0, 1, 2, 3, 7]].drop(0)  # .rename(columns=[col_names])
  col_names = ['probe_name', 'gene_symbol', 'm', 'b', 'is_qc_pass']
  qc_table.columns = col_names
  qc_table = qc_table.drop('gene_symbol', axis=1)
  print(qc_table.head())
    
  #miller_probe_info = load_workbook(os.path.join(config.scriptLocation, "data",  "Miller et al. doi.org_10.1186_1471-2164-15-154 12864_2013_7016_MOESM8_ESM.xlsx"))

  #sheet = miller_probe_info.get_sheet_by_name("Table")
  #probe_col = sheet['A']
  #qc_col = sheet['H']

  #probe_col = probe_col[2:]
  #qc_col = qc_col[2:]
  
  #probe_col = [cell.value for cell in probe_col]
  #qc_col = [cell.value for cell in qc_col]
  
  #qc_table = pandas.DataFrame(data={'probe_name' : probe_col, 'is_qc_pass' : qc_col})
  
  result_probe_info = pandas.merge(result_probe_info, qc_table, on = 'probe_name')
  
  #calculate statistics for probes
  probes_per_gene = result_probe_info.groupby(['gene_symbol']).size().reset_index(name='n')

  genes_with_probes = probes_per_gene.query('n > 1')['gene_symbol']
  print("Number of genes that have more than one probe:" + str(len(genes_with_probes)))
  
  #probes that are in genes that have more than one probe
  filtered_probes = result_probe_info[result_probe_info['gene_symbol'].isin(genes_with_probes)]
  print("Number of probes in genes that have more than one probe:" + str(len(filtered_probes.index)))
  
  #get statistics
  print(filtered_probes.groupby(['is_most_variable', 'is_qc_pass']).size())
   
  print("Writing enhanced probe information to " + probe_processed_filename)
  result_probe_info.to_csv(probe_processed_filename, index=False)
    
def write_PA_call_sums():
  #calculate information on present/absent calls

  donorFolders = [f for f in listdir(config.microarrayFolder) if 'normalized_microarray_donor' in f]
  combined_sums = pandas.DataFrame()
  brain_sd_col_list = []
  for donorFolder in donorFolders:
    brain_number = donorFolder.split("_donor")[1]
    brain_number = "Donor " + str(brain_number)
    brain_sd_col_list.append(brain_number)

    print ("PROCESSING present/absent calls for " + donorFolder)
    annotation_filename = os.path.join(config.microarrayFolder, donorFolder, "SampleAnnot.csv")
    pa_data_filename = os.path.join(config.microarrayFolder, donorFolder, "PACall.csv")
    
    annotations = pandas.read_csv(annotation_filename)
    pa_data = pandas.read_csv(pa_data_filename, header=None)
    pa_data = pa_data.sum(axis=0)
    pa_data = pa_data.drop(0).reset_index()

    #get annotation names
    header = annotations.apply(get_header, axis=1)
    merged = pandas.concat([header, pa_data], axis=1, ignore_index=True)
    merged.drop(1, axis=1, inplace=True)
    merged.columns = ['label', 'present_calls']
    merged['donor'] = brain_number
    #merge
    combined_sums = pandas.concat([combined_sums, merged], axis=0)
    print("Result shape so far:" + str(combined_sums.shape))
    

  #write out as single file
  print("Writing PA sums to " + sums_processed_filename)
  combined_sums.to_csv(sums_processed_filename, index=False)


def main():
    
  OUTPUT_FOLDER = config.processedOutputLocation

  if not os.path.isfile(probe_processed_filename): 
    write_probe_information()

  donorFolders = [f for f in listdir(config.microarrayFolder) if 'normalized_microarray_donor' in f]

  use_mni_coordinates = True
  
  #probe_strategy = "all"
  #probe_strategy = "highSTD"
  probe_strategy = "passQC"
  #probe_strategy = "scale" #scale method is not fully setup - the range of values varies 

  print("Probe filter strategy: " + str(probe_strategy))
  print("Using MNI coordinates? " + str(use_mni_coordinates))

  for donorFolder in donorFolders:
    brain_number = donorFolder.split("_donor")[1]

    print ("PROCESSING " + donorFolder)
    brain_number = donorFolder.split("_donor")[1]

    probe_filename = probe_processed_filename
    annotation_filename = os.path.join(config.microarrayFolder, donorFolder, "SampleAnnot.csv")
    expression_data_filename = os.path.join(config.microarrayFolder, donorFolder, "MicroarrayExpression.csv")

    if use_mni_coordinates:
        output_filename = os.path.join(OUTPUT_FOLDER, str(brain_number) + ".matrix.regionID.MNI(xyz).tsv")
    else:
        output_filename = os.path.join(OUTPUT_FOLDER, str(brain_number) + ".matrix.regionID.MRI(xyz).tsv")

    process_file(probe_filename, annotation_filename, expression_data_filename, output_filename, use_mni_coordinates, probe_strategy)
    print("Matrix written to:" + output_filename)


if __name__ == '__main__':
  main()
