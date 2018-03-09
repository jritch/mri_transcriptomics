library(broom)
library(magrittr)
library(ggplot2)
library(readr)
library(dplyr)

regionOfInterest <- "cortex_excluding_piriform_hippocampus"
#top and bottom 5

modalityOfInterest <- "T1T2Ratio"
xLabel <- "T1-/T2-w Ratio"

if (Sys.info()['nodename'] == "RES-C02RF0T2.local") {
  #set working directory
  setwd("/Users/lfrench/Desktop/results/mri_transcriptomics/")
}
single_gene_folder <- "./results/single_gene_data/"

source("get_expression_single_gene_function.R")


genesOfInterest <- c("VAMP1","SYT2","ESRRG","AGPAT9","RHBDL3","PCDH20","FRMPD2","KIAA1644","RBP4","SCARA5")


all_genes <- NULL
for (gene in genesOfInterest) {
  single_table <- get_expression_single_gene(gene, regionOfInterest, modalityOfInterest, single_gene_folder)
  if (!is.null(all_genes)) {
    single_table %<>% dplyr::select(`(x,y,z)`, gene)
    all_genes <- inner_join(all_genes, single_table, by="(x,y,z)") 
  } else {
    all_genes <- single_table
  }
}

#check if we got all the genes
setdiff(genesOfInterest, colnames(all_genes))

cor(all_genes[, genesOfInterest])
mean(abs(cor(all_genes[, genesOfInterest])))

all_genes %>% group_by(donor) %>% do(glance(summary(lm(MRI_Intensity ~ VAMP1+SYT2+ESRRG+AGPAT9+RHBDL3+PCDH20+FRMPD2+KIAA1644+RBP4+SCARA5, data=.))))

range((all_genes %>% group_by(donor) %>% do(glance(summary(lm(MRI_Intensity ~ VAMP1+SYT2+ESRRG+AGPAT9+RHBDL3+PCDH20+FRMPD2+KIAA1644+RBP4+SCARA5, data=.)))) )$adj.r.squared )
mean((all_genes %>% group_by(donor) %>% do(glance(summary(lm(MRI_Intensity ~ VAMP1+SYT2+ESRRG+AGPAT9+RHBDL3+PCDH20+FRMPD2+KIAA1644+RBP4+SCARA5, data=.)))) )$adj.r.squared )
median((all_genes %>% group_by(donor) %>% do(glance(summary(lm(MRI_Intensity ~ VAMP1+SYT2+ESRRG+AGPAT9+RHBDL3+PCDH20+FRMPD2+KIAA1644+RBP4+SCARA5, data=.)))) )$adj.r.squared )



########################################
########################################
# Bitter taste receptors
########################################

genesOfInterest <- c("TAS2R43","TAS2R30","TAS2R31","TAS2R14","TAS2R50","TAS2R5","TAS2R10","TAS2R4","TAS2R3","TAS2R1")
all_genes <- NULL
for (gene in genesOfInterest) {
  single_table <- get_expression_single_gene(gene, regionOfInterest, modalityOfInterest, single_gene_folder)
  if (!is.null(all_genes)) {
    single_table %<>% dplyr::select(`(x,y,z)`, gene)
    all_genes <- inner_join(all_genes, single_table, by="(x,y,z)") 
  } else {
    all_genes <- single_table
  }
}
#check if we got all the genes
setdiff(genesOfInterest, colnames(all_genes))

cor(all_genes[, genesOfInterest])
mean(abs(cor(all_genes[, genesOfInterest])))

corByDonor <- all_genes %>% 
  group_by(donor) %>%
  do(data.frame(Cor=cor(.[, genesOfInterest])))
corByDonor$gene <- genesOfInterest
dplyr::filter(corByDonor, gene=="TAS2R43") %>% dplyr::select(Cor.TAS2R30)

one <- all_genes %>% filter(donor == "Donor 10021")
one 
summary(lm(MRI_Intensity ~ TAS2R43+TAS2R30+TAS2R31+TAS2R14+TAS2R50+TAS2R5+TAS2R10+TAS2R4+TAS2R3+TAS2R1, data=one))


all_genes %>% group_by(donor) %>% do(glance(summary(lm(MRI_Intensity ~ TAS2R43+TAS2R30+TAS2R31+TAS2R14+TAS2R50+TAS2R5+TAS2R10+TAS2R4+TAS2R3+TAS2R1, data=.))))
range((all_genes %>% group_by(donor) %>% do(glance(summary(lm(MRI_Intensity ~ TAS2R43+TAS2R30+TAS2R31+TAS2R14+TAS2R50+TAS2R5+TAS2R10+TAS2R4+TAS2R3+TAS2R1, data=.)))) )$adj.r.squared )
mean((all_genes %>% group_by(donor) %>% do(glance(summary(lm(MRI_Intensity ~ TAS2R43+TAS2R30+TAS2R31+TAS2R14+TAS2R50+TAS2R5+TAS2R10+TAS2R4+TAS2R3+TAS2R1, data=.)))) )$adj.r.squared )

mean((all_genes %>% group_by(donor) %>% do(glance(summary(lm(MRI_Intensity ~ TAS2R43+TAS2R30+TAS2R31+TAS2R14+TAS2R50, data=.)))) )$adj.r.squared )

########################################
########################################
# Mito genes
########################################

genesOfInterest <- c("CHCHD6","SNCA","BCKDK","TIMM10","TOMM5","NDUFA13","C19orf70","NDUFAB1","TOMM6","ATP5J2")
all_genes <- NULL
for (gene in genesOfInterest) {
  single_table <- get_expression_single_gene(gene, regionOfInterest, modalityOfInterest, single_gene_folder)
  if (!is.null(all_genes)) {
    single_table %<>% dplyr::select(`(x,y,z)`, gene)
    all_genes <- inner_join(all_genes, single_table, by="(x,y,z)") 
  } else {
    all_genes <- single_table
  }
}
#check if we got all the genes
setdiff(genesOfInterest, colnames(all_genes))

cor(all_genes[, genesOfInterest])
mean(abs(cor(all_genes[, genesOfInterest])))

corByDonor <- all_genes %>% 
  group_by(donor) %>%
  do(data.frame(Cor=cor(.[, genesOfInterest])))
corByDonor$gene <- genesOfInterest

one <- all_genes %>% filter(donor == "Donor 10021")
one 
summary(lm(MRI_Intensity ~ CHCHD6+SNCA+BCKDK+TIMM10+TOMM5+NDUFA13+C19orf70+NDUFAB1+TOMM6+ATP5J2, data=one))

all_genes %>% group_by(donor) %>% do(glance(summary(lm(MRI_Intensity ~ CHCHD6+SNCA+BCKDK+TIMM10+TOMM5+NDUFA13+C19orf70+NDUFAB1+TOMM6+ATP5J2, data=.))))
range((all_genes %>% group_by(donor) %>% do(glance(summary(lm(MRI_Intensity ~ CHCHD6+SNCA+BCKDK+TIMM10+TOMM5+NDUFA13+C19orf70+NDUFAB1+TOMM6+ATP5J2, data=.)))) )$adj.r.squared )
mean((all_genes %>% group_by(donor) %>% do(glance(summary(lm(MRI_Intensity ~ CHCHD6+SNCA+BCKDK+TIMM10+TOMM5+NDUFA13+C19orf70+NDUFAB1+TOMM6+ATP5J2, data=.)))) )$adj.r.squared )


mean((all_genes %>% group_by(donor) %>% do(glance(summary(lm(MRI_Intensity ~ CHCHD6+SNCA+BCKDK+TIMM10+TOMM5, data=.)))) )$adj.r.squared )

############################
# histone

genesOfInterest <- c("HR","KDM4C","ARID5B","KDM5A","KDM5D","JMJD1C","PHF2","KDM2B","PHF8","KDM2A")



all_genes <- NULL
for (gene in genesOfInterest) {
  single_table <- get_expression_single_gene(gene, regionOfInterest, modalityOfInterest, single_gene_folder)
  if (!is.null(all_genes)) {
    single_table %<>% dplyr::select(`(x,y,z)`, gene)
    all_genes <- inner_join(all_genes, single_table, by="(x,y,z)") 
  } else {
    all_genes <- single_table
  }
}
#check if we got all the genes
setdiff(genesOfInterest, colnames(all_genes))

cor(all_genes[, genesOfInterest])
mean(abs(cor(all_genes[, genesOfInterest])))

corByDonor <- all_genes %>% 
  group_by(donor) %>%
  do(data.frame(Cor=cor(.[, genesOfInterest])))
corByDonor$gene <- genesOfInterest

one <- all_genes %>% filter(donor == "Donor 10021")
one 
summary(lm(MRI_Intensity ~ HR+KDM4C+ARID5B+KDM5A+KDM5D+JMJD1C+PHF2+KDM2B+PHF8+KDM2A, data=one))

all_genes %>% group_by(donor) %>% do(glance(summary(lm(MRI_Intensity ~ HR+KDM4C+ARID5B+KDM5A+KDM5D+JMJD1C+PHF2+KDM2B+PHF8+KDM2A, data=.))))
range((all_genes %>% group_by(donor) %>% do(glance(summary(lm(MRI_Intensity ~ HR+KDM4C+ARID5B+KDM5A+KDM5D+JMJD1C+PHF2+KDM2B+PHF8+KDM2A, data=.)))) )$adj.r.squared )
mean((all_genes %>% group_by(donor) %>% do(glance(summary(lm(MRI_Intensity ~ HR+KDM4C+ARID5B+KDM5A+KDM5D+JMJD1C+PHF2+KDM2B+PHF8+KDM2A, data=.)))) )$adj.r.squared )

mean((all_genes %>% group_by(donor) %>% do(glance(summary(lm(MRI_Intensity ~ HR+KDM4C+ARID5B+KDM5A+KDM5D, data=.)))) )$adj.r.squared )
