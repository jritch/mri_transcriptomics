#!python
import analysis
import sys
import itertools
import nibabel 
from PIL import Image
import numpy as np

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def check_slices(data,coord):
	
	data = np.fliplr(np.transpose(data,(2, 0, 1)))

	# take a z-slice
	my_slice = data[:,:,coord]

	# open a plot
	plt.imshow(np.invert(my_slice), cmap='Greys')
	
	# Title the plot & axes
	plt.xlabel('x')
	plt.ylabel('y')
	plt.title('Subject 9861, z='+str(coord))
	plt.show()

	my_slice = data[:,coord,:]

	# open a plot
	plt.imshow(np.invert(my_slice), cmap='Greys')
	
	# Title the plot & axes
	plt.xlabel('x')
	plt.ylabel('z')
	plt.title('Subject 9861, y='+str(coord))
	plt.show()

	my_slice = data[coord,:,:]

	# open a plot
	plt.imshow(np.invert(my_slice), cmap='Greys')
	
	# Title the plot & axes
	plt.xlabel('y')
	plt.ylabel('z')
	plt.title('Subject 9861, x='+str(coord))
	
	plt.show()


def check_slice(data,z):

	my_slice = data[:,:,z]

	# open a plot
	plt.imshow(np.invert(my_slice), cmap='Greys')
	
	# Title the plot & axes
	plt.xlabel('x')
	plt.ylabel('y')
	plt.title('Subject 9861, Z=96 slice')
	
	plt.show()

	raw_input("Now try to fix misalignment.\nPress any key to continue.")
	
	# Now try to fix misalignment 

	# transpose array axes to compensate for header misinterpretation
	transformed = np.fliplr(np.transpose(data,(2, 0, 1)))

	# take a z-slice
	my_slice = transformed[:,:,z]

	# open a plot
	plt.imshow(np.invert(my_slice), cmap='Greys')
	plt.xlabel('x')
	plt.ylabel('y')
	plt.title('Subject 9861, Z=96 slice')
	
	plt.show()

def check_brain_3d(data):
	z,x,y = (data > 375).nonzero()
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.scatter(x, y, z, zdir='z', c= 'red')
	plt.show()

def check_header(filename):
		img = nibabel.load(filename)
		print img.header()
		print img.data.shape

def print_intensities(imgs,coords):
	x,y,z = coords
	intensities = [img[x,y,z] for img in imgs]
	print "T1= ",intensities[0], "\tT2= ",intensities[1], "\tratio= ",intensities[2] 

def main():
	"""
			Args are dirname x y z.
			Prints the intensity at coord XYZ
	"""

	if not len(sys.argv) == 5:
		print "Usage: check_pixel dirname x y z\nDirname should contain files T1.nii, T2.nii."
		sys.exit(1)

	print sys.argv
	script,dirname,x,y,z = sys.argv
	print dirname	
	imgs = analysis.load_nifti_data(dirname)
	intensities = [img[x,y,z] for img in imgs]
	print "T1= ",intensities[0], "\tT2= ",intensities[1], "\tratio= ",intensities[2] 

	# Try out other permutations....
	potential_coords = set(itertools.permutations([x,y,z]))
	for coords in potential_coords:
		print coords
		print_intensities(imgs,coords)

if __name__ == '__main__':
	dirname = "C:\Users\Jacob\large_thesis_files\AllenHBAProcessedExpressionAndMRIs\\normalized_microarray_donor9861"
	imgs = analysis.load_nifti_data(dirname)
	#check_slices(imgs[0],96)
	check_slice(imgs[0],96)
	#check_brain_3d(imgs[0])