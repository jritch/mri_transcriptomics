import os 

import config
import analysis

import numpy as np
import matplotlib.pyplot as plt

files = config.expression_filenames
brain_ids = [f.split(".")[0] for f in files]

measures = ["T1-w","T2-w","T1-w T2-w Ratio"]
for m in range(len(measures)):
    percentages = []
    for i in range(len(brain_ids)):
        MRI_data = analysis.load_nifti_data(config.basePathMRI + brain_ids[i])
        gene_exp_fh = open(os.path.join(config.expressionFolder,files[i]))
        coords, a  = analysis.get_coords_and_region_ids_from_gene_exp_data(gene_exp_fh)
        gene_exp_fh.close()

        # get the flat coordinates for the T1T2Ratio MRI
        flat_mri_data = analysis.flatten_mri_data(MRI_data[m],coords)
        unique_values = len(set(flat_mri_data))
        total_values = len(flat_mri_data)
        print "In Brain {} there are {} unique values / {} total values".format(brain_ids[i],unique_values,total_values)
        percentages.append(100 - float(unique_values)/total_values * 100)
    print percentages

    import seaborn as sns

    bar_plot = sns.barplot(x=brain_ids,y=percentages,
                            palette="muted",
                           )

    plt.xlabel("Brain ID")
    plt.ylim([0,100])
    plt.ylabel("Ties (%)")
    plt.title("Number of Ties for {} MR Image".format(measures[m]))
    plt.savefig("/Users/jritchie/git-repos/mri_transcriptomics/" +measures[m] + ".png")
    plt.clf()