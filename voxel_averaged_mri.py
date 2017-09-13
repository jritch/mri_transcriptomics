import config, analysis, nibabel, os, csv, sys
import numpy as np
import scipy.spatial, scipy.ndimage.filters

def voxel_average(arr,voxel_size=1):
    kernel = np.ones((3,3,3))
    kernel /= 27.0
    return scipy.ndimage.filters.convolve(arr,kernel,mode='constant')

#deprecated
def write_avg_image(dirname="C:\Users\Jacob\large_thesis_files\AllenHBAProcessedExpressionAndMRIs\\normalized_microarray_donor9861"):
    imgs = analysis.load_nifti_data(dirname)
    ratio = nibabel.Nifti1Image(imgs[0], np.eye(4,4))
    ratio.to_filename(dirname+"\\"+"our_T1.nii")
    ratio = nibabel.Nifti1Image(imgs[2], np.eye(4,4))
    ratio.to_filename(dirname+"\\"+"ratio.nii")
    avg_ratio = nibabel.Nifti1Image(voxel_average(imgs[2]), np.eye(4,4)) 
    avg_ratio.to_filename(dirname+"\\"+"ratio_averaged.nii")

def main():
    write_avg_image()
    return 

def test():
  test_array = np.array([[[1.,2.,3.,4.,5.],
                          [1.,2.,3.,4.,5.],
                          [1.,2.,3.,4.,5.]],
                          [[1.,2.,3.,4.,5.],
                          [1.,2.,3.,4.,5.],
                          [1.,2.,3.,4.,5.]],
                          [[1.,2.,3.,4.,5.],
                          [1.,2.,3.,4.,5.],
                          [1.,2.,3.,4.,5.]]])
  voxel_average(test_array)


if __name__ == '__main__':
    main()
