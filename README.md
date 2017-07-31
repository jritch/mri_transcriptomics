# Code to replicate analysis from "Magnetic Resonance Imaging From The Transcriptomic Perspective"

This project is compatible with Python v2.7 and R v3.4.0

## Inputs 

To run this analysis, the T1-w and T2-w MRI scans (http://human.brain-map.org/mri_viewers/data) and normalized microarray data (http://human.brain-map.org/static/download) from the Allen Human Brain Atlas project must first be downloaded.

If you wish to change the regions of interest considered in the analysis, you should consult the Ontology.csv file which is provided with the normalized microarray data.

You will also need to install the following python libraries (installation using Anaconda is recommended) : 

* numpy
* scipy
* pandas
* nibabel

And the following R packages :

* AnnotationDbi
* GO.db
* annotate
* cowplot
* dplyr
* ggplot2
* homologene
* magrittr
* metap
* optparse
* org.Hs.eg.db
* readr
* tidyr
* tmod
* xlsx

## Configuration 

You should update config.py with the absolute paths to the locations of your input files and to appropriately named directories for your intermediate and output files

## How to run

### Main analysis

1. Run average_probes.py to obtain average expression values for each gene from the AHBA normalized expression data.
2. Run analysis.py on the output files of the previous step to generate .csv files summarizing the correlation between gene expression and MRI intensity in each donor for each image and region of interest.
3. Run R Code/RunSingleGO.AUROC.Analysis.R to perform gene set enrichment analysis on the .csv files from the previous step.

### Additional analyses

(All of these analyses depend on the output of average_probes.py, so you should run the main analysis first)

1. Run single_gene.py followed by R Code/RunSingleGO.AUROC.Analysis.R to obtain correlation scatterplots for individual genes (edit the files to change the R.O.Is and gene names).
2. Run one_sided_lists.py to obtain a summary of the number of significant (after FDR-correction) correlations between T1-w/T2-w ratio image intensity and gene expression.
3. Run write_t1t2_img.py to dump out an example T1-w/T2-w ratio image in NIFTI1 format.
