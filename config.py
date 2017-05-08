import platform

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
    ontologyFolder = baseAllenFolder + "../"
    expressionFolder = "/Users/jritchie/data/AllenHBAProcessedExpressionWithBrainID/"

    outputCSVFolder = "/Users/jritchie/Google Drive/4th year/Thesis/gene_list_csvs_float/"
    basePathMRI = "/Users/jritchie/data/AllenHBAProcessedExpressionAndMRIs/normalized_microarray_donor"

# Leon's CAMH Laptop
else:
    baseAllenFolder = "/Users/lfrench/Desktop/data/Allen/HBA/"
    ontologyFolder = "/Users/lfrench/Desktop/results/mri_transcriptomics/data/"
    expressionFolder = "/Users/lfrench/Desktop/data/Allen/HBA/ignore.regoinIDDataForJacob/"
    outputCSVFolder = "/Users/lfrench/Desktop/results/mri_transcriptomics/results/"
    basePathMRI = "/Users/lfrench/Desktop/data/Allen/HBA/normalized_microarray_donor"