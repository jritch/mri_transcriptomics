from analysis import fisher_p
import statsmodels.sandbox.stats.multicomp
import pandas
import numpy as np
import config
from copy import deepcopy

def one_sided(x,up=True):
    y = deepcopy(x)

    for i in range(len(x)):
        if up:
            if x[i][0] < 0:
              y[i][1] = (1 - y[i][1])
            else:
               y[i][1] = (y[i][1])
        else:
            if x[i][0] > 0:
              y[i][1] = (1 - y[i][1])
            else:
               y[i][1] = (y[i][1])
    return y

def get_pairs(x):
    return [[x[i],x[i+6]]for i in range(1,7)]

def get_one_sided_lists(data):
    """
    Uses the output of analysis.py to get 2 lists of sets of 6 one-sided (but otherwise unadjusted) correlation / p-value pairs
    """

    pairs = data[1:].apply(get_pairs, axis=1)

    one_sided_pairs = (pairs.apply(one_sided,args=(True,)),pairs.apply(one_sided,args=(False,)))

    return one_sided_pairs

def get_single_two_sided_list(data):
    pairs = data[1:].apply(get_pairs, axis=1)
    return pairs

def adjust_one_sided_lists(lists):
    """
    Takes the two lists of sets returned by get_one_sided_lists() and return two lists 
    """
    adjusted_lists = []
    for l in lists:
        meta_list = []
        #for both the lists (up & down)
        for gene in l:
            p_values = [x[1] for x in gene]
            meta_p = fisher_p(p_values)
            meta_list.append(meta_p)
        adjusted_meta_list = statsmodels.sandbox.stats.multicomp.multipletests(meta_list,method="fdr_bh")[1]
        adjusted_lists.append(adjusted_meta_list)
    return adjusted_lists

def adjust_single_list(l):
    meta_list = []
    #for both the lists (up & down)
    for gene in l:
        p_values = [x[1] for x in gene]
        meta_p = fisher_p(p_values)
        meta_list.append(meta_p)
    adjusted_meta_list = statsmodels.sandbox.stats.multicomp.multipletests(meta_list,method="fdr_bh")[1]
    return adjusted_meta_list

def main():
    # TODO: change to using config.py
    # For nowm, this is a file that has T1/T2 ratio in the cortex
    filename = "/Users/jritchie/data/garbage2/T1T2Ratio.cortex.gene_list.csv"
    #filename = config.outputCSVFolder + "/T1T2Ratio.cortex.gene_list.csv"
    data = pandas.read_csv(filename, sep=",")
    lists = get_one_sided_lists(data)
    l = get_single_two_sided_list(data)
    adj_single_list = adjust_single_list(l)
    adj_lists = adjust_one_sided_lists(lists)
    adj=pandas.DataFrame(data= {"gene_symbol": data.iloc[1:,0],"up":adj_lists[0],"down":adj_lists[1],"two_sided":adj_single_list}).reindex_axis(["gene_symbol","up","down","two_sided"],axis=1)
    return adj

if __name__ == '__main__':
    adj = main()