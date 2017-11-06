import config, analysis, os, csv, sys, inspect, pandas, itertools
import matplotlib.pyplot as plt
import numpy as np

from ontology import *
from single_gene import get_single_gene_data

def generate_tables():

    regionIDs = [4005] #full brain
    region_name = "whole_brain"
    to_exclude_list = None

    gene_names = [
                  "TAS2R43",
                  "TAS2R30",
                  "TAS2R31",
                  "TAS2R14",
                  "TAS2R50"
               ]

    files = config.expression_filenames

    brain_ids = [f.split(".")[0] for f in files]

    outputFolder = os.path.join(config.scriptLocation, "results", "correlated_genes_data")
    if not os.path.exists(outputFolder):
      os.mkdir(outputFolder)

    filenames = [os.path.join(outputFolder, f + "." + ".".join(gene_names) + "." + region_name + ".raw.csv") for f in brain_ids]

    o = Ontology(os.path.join(config.scriptLocation, "data",  "Ontology.csv"))

    cortex_divisions = [str(x) for x in o.hierarchy.get(4008)]

    num_brains = len(brain_ids)

    for i in range(num_brains):
        df = pandas.DataFrame()
        for gene_name in gene_names:
            gene_exp_fh = open(os.path.join(config.processedOutputLocation,files[i]))
            coords,coord_to_region_map = analysis.get_coords_and_region_ids_from_gene_exp_data(gene_exp_fh)
            indices = []
            for j in range(len(regionIDs)):
              indices += analysis.get_flat_coords_from_region_id(regionIDs[j],coords,coord_to_region_map,o,to_exclude_list=to_exclude_list)
            assert(len(set(indices)) == len(indices))

            single_gene_data = np.array(get_single_gene_data(gene_exp_fh,gene_name,indices))
            df[gene_name] = single_gene_data
            gene_exp_fh.close()

        df.to_csv(filenames[i],index=False)
    return filenames

def get_correlations(df_filename,genes=None):
    df = pandas.read_csv(df_filename)
    if not genes:
        genes = df.columns

    table_dict = {}
    for i in genes:
        table_dict[i] = {}

    for i in itertools.combinations(genes,2):
        corr = analysis.R_VALUE_FUNCTION(df[i[0]],df[i[1]])
        print "Correlation between",i[0],"and",i[1],":",corr
        table_dict[i[0]][i[1]] = tuple(corr)
        table_dict[i[1]][i[0]] = tuple(corr)

    filename = df_filename.replace(".raw.csv",".correlation.tsv")
    pandas.DataFrame.from_dict(table_dict).to_csv(filename,sep="\t")

    def plot_two_genes(df_filename,genes):
        df = pandas.read_csv(df_filename)
        plt.scatter(df[genes[0]],df[genes[1]])
        plt.title("Correlation between "+genes[0] +" and " + genes[1])
        plt.show()

if __name__ == '__main__':
    filenames = generate_tables()
    get_correlations(filenames[0])
    #plot_two_genes(filenames[0],("TAS2R43","TAS2R30"))
    df = pandas.read_csv(filenames[0].replace(".raw.csv",".correlation.tsv"),delimiter="\t")
