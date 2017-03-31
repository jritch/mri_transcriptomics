import sys
import math

def read_single_gene_csv(filename):
	with open(filename,"r") as f:
		f.readline()
		coord_to_intensity_map = {}
		for line in f.readlines():
			entries = line.split(",")
			coords = eval(entries[0] + "," + entries[1] + "," +  entries[2])
			#print entries, coords, entries[-1]
			coord_to_intensity_map[coords] = entries[3]
#			print entries[3]
	return coord_to_intensity_map

if __name__ == '__main__':
	#print "C:\Users\Jacob\Google Drive\4th year\Thesis\single_gene_data\9861.CAPN6.cortex.MRI(xyz).expression.csv"
	dict1 = read_single_gene_csv("C:\Users\Jacob\Google Drive\\4th year\Thesis\single_gene_data\9861.CAPN6.cortex.MRI(xyz).expression.csv")
	dict2 = read_single_gene_csv("C:\Users\Jacob\Google Drive\\4th year\Thesis\single_gene_data_new\9861.CAPN6.cortex.MRI(xyz).expression.csv")
	dict3 = read_single_gene_csv("C:\Users\Jacob\Google Drive\\4th year\Thesis\single_gene_data_avg\9861.CAPN6.cortex.MRI(xyz).expression.csv")
	
	dists = []
	diffs =[]
	for key in dict1:
		coord_tup = eval(key)
		#dist = math.sqrt((coord_tup[0] - 90) ** 2 + (coord_tup[1] - 90) ** 2 + (coord_tup[2] - 90) ** 2)
		dist = (coord_tup[1] - 90) ** 2

		dists.append(dist)
		diffs.append( float(dict1[key]) - float(dict2[key]) )


	import matplotlib.pyplot as plt

	plt.plot(dists,diffs, 'ro')
	plt.show()

