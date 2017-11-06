import platform
import os
from os import listdir
from os.path import isfile, join

# Leon's CAMH Laptop
if platform.node() =='RES-C02RF0T2.local':
    microarrayFolder = "/Users/lfrench/Desktop/results/mri_transcriptomics/data/allen expression/"
    figShareFolder = "/Users/lfrench/Desktop/results/mri_transcriptomics/data/figshare data/"
elif platform.node() == "Kurosawa":
    microarrayFolder = "C:/Users/Jacob/mri_transcriptomics/data/allen expression/"
    figShareFolder = "C:/Users/Jacob/mri_transcriptomics/data/figshare data/"
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
