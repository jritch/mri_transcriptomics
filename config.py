import platform
if 'Windows-10' in platform.platform():
    baseAllenFolder = "C:/Users/Jacob/large_thesis_files/AllenHBAProcessedExpressionAndMRIs/"
    ontologyFolder = baseAllenFolder + "../"
    MRIFolder = "C:/Users/Jacob/large_thesis_files/AllenHBAProcessedExpressionWithBrainID/"
    outputCSVFolder = "C:/Users/Jacob/Google Drive/4th Year/Thesis/gene_list_csvs/"
else:
    baseAllenFolder = "/Users/lfrench/Desktop/data/Allen/HBA/"
    #ontologyFolder = ...
    #MRIFolder = ...
    #outputCSVFolder = ...