import numpy
import numpy as np
import analysis
import nibabel
import config
import os

from nilearn import plotting
import matplotlib.pyplot as plt
from scipy import stats
import scipy

f, ax = plt.subplots(3);

files = get_expression_filenames(use_MNI=False)
brain_ids = [f.split(".")[0] for f in files]

brain_to_use = brain_ids[0]

# How we do it in analysis
data = analysis.load_nifti_data(brain_to_use)

#get the structure from a single file
print(config.figShareFolder + 'm'+ brain_to_use + "_T1.nii.gz")
img = nibabel.load(config.figShareFolder + 'm'+ brain_to_use + "_T1.nii.gz")
hdr = img.header
#print hdr
#print img.affine

#plotting.plot_glass_brain("T1.nii", axes=ax[0],alpha=0,threshold=0)
#plotting.plot_glass_brain("T2.nii", axes=ax[1], alpha=0,threshold=0);

#changed_data = numpy.nan_to_num(data[2]) #* np.mean( data[1][np.nonzero(data[1])] )

changed_data = data[2]

changed_data = numpy.nan_to_num(data[2])

#  pre-processing and z-scoring

changed_data[ data[1] == 0] = 2
changed_data[ data[0] == 0] = 0

# Enable this to smooth the image

#changed_data = scipy.ndimage.filters.gaussian_filter(changed_data,0.3)

#changed_data[ data[1] == 0] = 2
#changed_data[ data[0] == 0] = 0

changed_data = (0 + stats.zscore(numpy.nan_to_num(changed_data)))

#changed_data[ data[1] == 0] = 0
changed_data[ data[0] == 0] = 0

img = nibabel.Nifti1Image(data[2] , hdr.get_sform())
img.to_filename(config.figShareFolder + '/m' + brain_to_use +"_raw_ratio.nii")

img = nibabel.Nifti1Image(changed_data , hdr.get_sform())
img.to_filename(config.figShareFolder + '/m' + brain_to_use +"_z_scored_ratio.nii")

plotting.plot_glass_brain(config.figShareFolder + '/m' + brain_to_use + "_z_scored_ratio.nii",axes=ax[2], alpha=0,threshold=0);
#plotting.plot_glass_brain("ratio.nii", axes=ax[1], alpha=0);
#plotting.plot_glass_brain("ratio.nii", axes=ax[2], alpha=0);

# uncomment to display plot
#plt.show()
#plt.savefig("3_images.png")
