#!python
import analysis
import sys
import itertools

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
	main()