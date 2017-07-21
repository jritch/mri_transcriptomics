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

files = config.expression_filenames
brain_ids = [f.split(".")[0] for f in files]

example_dir = config.basePathMRI + brain_ids[0]

# How we do it in analysis
# Note: this calls " os.chdir(example_dir) "

data = analysis.load_nifti_data(example_dir)

img = nibabel.load("T1.nii")
hdr = img.header
#print hdr
#print img.affine

plotting.plot_glass_brain("T1.nii", axes=ax[0],alpha=0,threshold=0)
plotting.plot_glass_brain("T2.nii", axes=ax[1], alpha=0,threshold=0);

os.chdir(analysis.baseProjectFolder)

#changed_data = numpy.nan_to_num(data[2]) #* np.mean( data[1][np.nonzero(data[1])] )

changed_data = data[2]

changed_data = numpy.nan_to_num(data[2])

changed_data[ data[1] == 0] = 2
changed_data[ data[0] == 0] = 0

# Enable this to smooth the image

#changed_data = scipy.ndimage.filters.gaussian_filter(changed_data,0.3)

#changed_data[ data[1] == 0] = 2
#changed_data[ data[0] == 0] = 0

changed_data = (0 + stats.zscore(numpy.nan_to_num(changed_data)))

#changed_data[ data[1] == 0] = 0
changed_data[ data[0] == 0] = 0

changed_data = numpy.nan_to_num(changed_data)

img = nibabel.Nifti2Image(changed_data, hdr.get_sform())

img.to_filename("z_scored_ratio.nii")

# See if it's aligned
#plotting.plot_glass_brain("ratio.nii")
#plt.show()

plotting.plot_glass_brain("z_scored_ratio.nii", axes=ax[2], alpha=0,threshold=0);
#plotting.plot_glass_brain("ratio.nii", axes=ax[1], alpha=0);
#plotting.plot_glass_brain("ratio.nii", axes=ax[2], alpha=0);

plt.show()
#plt.savefig("3_images.png")