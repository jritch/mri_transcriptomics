import os
import sys

SOURCE_DIR = "./gene_list_csvs"
DEST_DIR = "./gene_list_csvs_stripped"
#os.mkdir(DEST_DIR)

csv_files = os.listdir(SOURCE_DIR)

csv_files  = filter(lambda x: "csv" in x, csv_files)

print csv_files

#sys.exit()

probes = []
unprobes = []

for filename in  csv_files:
	flag = 1
	flag2 = 0
	with open(SOURCE_DIR + "/" + filename) as f1:
		with open(DEST_DIR + "/" + filename, "w") as f2:
				for line in f1:
					if flag:
							f2.write(line)
					if flag2:
						flag2 = 0
						line = line.strip("\n") + ",10021,12876,14380,15496,15697,9861\n"
					if flag:
						flag = 0
						flag2 = 1
					if "_" in line and not "RP" in line:
							pass
							#probes.append(line.split(",")[0])
					else:
#						print line
						f2.write(line)
						#unprobes.append(line.split(",")[0]
						pass

print len(set(probes))
print len(set(unprobes))