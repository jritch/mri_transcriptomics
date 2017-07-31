import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

a = pd.read_csv("/Users/jritchie/git-repos/mri_transcriptomics/single_gene_data_avg/12876.SLC17A7.cortex.MRI(xyz).expression.csv")
print(a.columns)

a = a[a["MRI_Intensity"].notnull()]["MRI_Intensity"].values
x = np.empty([10,a.shape[0]])
x[:,:] = a
plt.contourf(x)
plt.show()

a = pd.read_csv("subset_of_expression_values.csv",sep="\t")
third = len(a)/3

f, (ax1,ax2,ax3) = plt.subplots(3)

subplots = (ax1,ax2,ax3)

for i in range(len(subplots)):
	vals = a.iloc[:,i].values.flatten()
	vals[third:2*third] = 0
	sns.barplot(range(len(a)), vals, ax=subplots[i])

plt.show()
