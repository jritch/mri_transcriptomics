library(plotly)
library(tidyr)
library(readr)
singleGene <- read_csv("/Users/lfrench/Desktop/results/mri_transcriptomics/single_gene_data/9861.EPCAM.full_brain.MRI(xyz).expression.csv")
head(singleGene$`(x,y,z)`)

(singleGene <- separate(singleGene, `(x,y,z)`, c("x","y","z"), ",|[)]"))
singleGene$x <- gsub("[(]","", singleGene$x)

plot_ly(singleGene, x = ~x, y =~y, z = ~z)