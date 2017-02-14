import config
import analysis
import os, csv, sys
import numpy as np
import scipy.spatial
import scipy.ndimage.filters

def voxel_average(arr,voxel_size=1):
    """
    #print arr[30,120]
    idx =  np.vstack(np.where(np.ones(arr.size))).T
    #print idx
    #print idx.size
    
    tree = scipy.spatial.cKDTree(idx)

    neighbours = tree.query_ball_point(idx, r=1, p=np.inf)
    new_vals = np.hstack(np.mean(arr[n]) for n in neighbours)

    arr_new = arr.astype(np.double, copy=True)
    arr_new[idx] = new_vals
    """
    kernel = np.ones((3,3,3))
    kernel /= 27.0
    print scipy.ndimage.filters.convolve(arr,kernel,mode='constant')
    return None

def main():
  
  test_array = np.array([[[1.,2.,3.,4.,5.],[1.,2.,3.,4.,5.],[1.,2.,3.,4.,5.]],[[1.,2.,3.,4.,5.],[1.,2.,3.,4.,5.],[1.,2.,3.,4.,5.]],[[1.,2.,3.,4.,5.],[1.,2.,3.,4.,5.],[1.,2.,3.,4.,5.]]])
  print test_array
  voxel_average(test_array)

  sys.exit()

  MRI_dimension = 2 # 0: T1, 1: T2, 2: ratio

  files = ["10021.matrix.regionID.MRI(xyz).29131 x 893.txt",
           "12876.matrix.regionID.MRI(xyz).29131 x 363.txt",
           "14380.matrix.regionID.MRI(xyz).29131 x 529.txt",
           "15496.matrix.regionID.MRI(xyz).29131 x 470.txt",
           "15697.matrix.regionID.MRI(xyz).29131 x 501.txt",
           "9861.matrix.regionID.MRI(xyz).29131 x 946.txt"]
  
  brain_ids = [f.split(".")[0] for f in files]  
  
  num_brains = len(brain_ids)
  for i in range(num_brains):
    mri_data = analysis.load_nifti_data(config.baseAllenFolder + "normalized_microarray_donor" + brain_ids[i])[MRI_dimension]
    print voxel_average(mri_data)

if __name__ == '__main__':
    main()
