import platform
import os
from os import listdir
from os.path import isfile, join



# Jacob's Personal Laptop
if 'Windows-10' in platform.platform():
    baseAllenFolder = "C:/Users/Jacob/large_thesis_files/AllenHBAProcessedExpressionAndMRIs/"
    ontologyFolder = baseAllenFolder + "../"
    expressionFolder = "C:/Users/Jacob/large_thesis_files/AllenHBAProcessedExpressionWithBrainID/"
    outputCSVFolder = "C:/Users/Jacob/Google Drive/4th Year/Thesis/gene_list_csvs_float/"
    basePathMRI = "C:/Users/Jacob/large_thesis_files/AllenHBAProcessedExpressionAndMRIs/normalized_microarray_donor"

# Jacob's CAMH laptop
elif platform.node() =='RES-CO2T6CS3GTFL.local':


    baseAllenFolder = "/Users/jritchie/data/AllenHBAProcessedExpressionAndMRIs/"
    #figShareFolder = 
    #expressionFolder = "/Users/jritchie/data/AllenHBAProcessedExpressionWithBrainID/"

    expressionFolder= "/Users/jritchie/git-repos/mri_transcriptomics/python_processed_expression_data"

    microarrayFolder = "/Users/jritchie/Documents/microarray_data"
    #outputCSVFolder = "/Users/jritchie/Google Drive/4th year/Thesis/gene_list_csvs_float/"

    if USE_BIAS_CORRECTED_IMAGES:
        outputCSVFolder = "/Users/jritchie/data/final_with_bias_correction/"
    else:
        outputCSVFolder = "/Users/jritchie/data/final/"

    if USE_BIAS_CORRECTED_IMAGES:
        basePathMRI = "/Users/jritchie/data/allen_bias_corrected/"
    
    else:
        basePathMRI = "/Users/jritchie/data/AllenHBAProcessedExpressionAndMRIs/normalized_microarray_donor"

    '''
    expressionFolder
    basePathMRI =

    '''
# Leon's CAMH Laptop
elif platform.node() =='RES-C02RF0T2.local':
    microarrayFolder = "/Users/lfrench/Desktop/results/mri_transcriptomics/data/allen expression/"
    figShareFolder = "/Users/lfrench/Desktop/results/mri_transcriptomics/data/figshare data/"
else:
    print("Please setup config.py")
    print(platform.node())


scriptLocation = os.path.dirname(os.path.realpath(__file__))
processedOutputLocation = os.path.join(scriptLocation, "data", "python_processed_expression_data")
if not os.path.exists(processedOutputLocation):
  os.makedirs(processedOutputLocation)

expression_filenames = [f for f in listdir(processedOutputLocation) if '.matrix.regionID.MRI(xyz).tsv' in f]

resultFolder = os.path.join(scriptLocation, "results")
if not os.path.exists(resultFolder):
  os.makedirs(resultFolder)


