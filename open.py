import nibabel
import numpy
import os
from subprocess import call

os.chdir("C:/Users/Jacob/Downloads/AllenHBAProcessedExpressionAndMRIs/normalized_microarray_donor10021")

call(["ls", "-l"])

NIFTI_FILES = ['T1.nii','T2.nii']

imgs = [] 
data = []
for f in NIFTI_FILES:
	img = nibabel.load(f)
	imgs.append(img)
	data.append(img.get_data())

print numpy.divide(data[0],data[1])

#print data




