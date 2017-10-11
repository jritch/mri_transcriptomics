# Code to replicate analysis from "Transcriptomic  characterization  of  MRI  contrast focused  on  the T1-w/T2-w  ratio  in  the  cerebral  cortex"

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

Download the three files from the figshare dataset (https://figshare.com/s/9a5fc70f3d0b7a6ab7fb)
* AllenMRImages.zip - requires unziping, with expanded files kept in the same folder as the below files. m prefix on files denotes bias corrected images.
* AllenHBA_DK_ExpressionMatrix.tsv - expression matrix from a Freesurfer view of the cortical transcriptome
* AllPhenocartaAnnotations.downloadedOct28.2016.tsv - annotations from Phenocarta

Download the Allen expression data at http://human.brain-map.org/static/download  (six files)

Update config.py with the the paths for these two sets of downloaded files. (microarrayFolder and figShareFolder)

## How to run

### Main analysis

1. Run average_probes.py to obtain average expression values for each gene from the AHBA normalized expression data. This file is created in the python_processed_expression_data folder that is created in the parent folder of the microarrayFolder. 
2. Run analysis.py on the output files of the previous step to generate .csv files summarizing the correlation between gene expression and MRI intensity in each donor for each image and region of interest. Regions of interest and MRI measures can be specified in the main method. 
3. Run R Code/RunSingleGO.AUROC.Analysis.R to perform gene set enrichment analysis on the .csv files from the previous step. Requires setting of figshare_data_folder (to location of figshare files) and setting of the working directory (start of the script).

### Additional analyses

(All of these analyses depend on the output of average_probes.py, so you should run the main analysis first)

1. Run single_gene.py followed by R Code/PlotSingleGene.R to obtain correlation scatterplots for individual genes (edit the main function variables to change the ROIs and gene names).
2. Run write_t1t2_img.py to write out an example T1-w/T2-w ratio image in NIFTI1 format.


## References for Data sources used

Mancarci, B. O., Toker, L., Tripathy, S. J., Li, B., and Rocco, B. (2017). Cross-laboratory analysis of brain cell type transcriptomes with applications to interpretation of bulk tissue data. bioRxiv. Available at: http://www.biorxiv.org/content/early/2017/05/18/089219.abstract. Gene lists taken from https://github.com/oganm/neuroExpressoAnalysis/tree/master/analysis/01.SelectGenes/Markers_Final/CellTypes/Cortex .

Hawrylycz, M. J., Lein, E. S., Guillozet-Bongaarts, A. L., Shen, E. H., Ng, L., Miller, J. A., et al. (2012). An anatomically comprehensive atlas of the adult human brain transcriptome. Nature 489, 391–399.

French, L., and Paus, T. (2015). A FreeSurfer view of the cortical transcriptome generated from the Allen Human Brain Atlas. Front. Neurosci. 9, 323.

Darmanis, S., Sloan, S. A., Zhang, Y., Enge, M., Caneda, C., Shuer, L. M., et al. (2015). A survey of human brain transcriptome diversity at the single cell level. Proc. Natl. Acad. Sci. U. S. A. 112, 7285–7290.

Zeisel, A., Muñoz-Manchado, A. B., Codeluppi, S., Lönnerberg, P., La Manno, G., Juréus, A., et al. (2015). Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq. Science 347, 1138–1142.

Zeng, H., Shen, E. H., Hohmann, J. G., Oh, S. W., Bernard, A., Royall, J. J., et al. (2012). Large-scale cellular-resolution gene profiling in human neocortex reveals species-specific molecular signatures. Cell 149, 483–496.

Portales-Casamar, E., Ch’ng, C., Lui, F., St-Georges, N., Zoubarev, A., Lai, A. Y., et al. (2013). Neurocarta: aggregating and sharing disease-gene relations for the neurosciences. BMC Genomics 14, 129.
