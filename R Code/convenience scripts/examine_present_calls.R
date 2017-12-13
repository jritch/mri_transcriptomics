pa_calls <- read_csv("/Users/lfrench/Desktop/results/mri_transcriptomics/data/python_processed_expression_data/PA_sums.csv")

regionOfInterest <- "cortex_excluding_piriform_hippocampus"
geneOfInterest <- "TAS2R4"
#regionOfInterest <- "whole_brain"

modalityOfInterest <- "T1T2Ratio"
xLabel <- "T1-/T2-w Ratio"


single_gene_folder <- "./results/single_gene_data/"
source("get_expression_single_gene_function.R")
allDonors <- get_expression_single_gene(geneOfInterest, regionOfInterest, modalityOfInterest, single_gene_folder)

allDonors %<>% mutate(label = paste0("RegionID:", regionID, "|LOC:", `(x,y,z)`))

allDonors <- left_join(allDonors, pa_calls)
allDonors %>% dplyr::select(present_calls, everything())
allDonors %>% arrange(desc(present_calls)) %>% dplyr::select(-cortical_division)

(result <- allDonors %>% group_by(donor) %>% summarize(n = n(), cor = cor(MRI_Intensity , present_calls, m='s'), p=cor.test(MRI_Intensity , present_calls, m='s')$p.value))
median(result$cor)
result <- allDonors %>% group_by(donor) %>% summarize(n = n(), cor = cor(Expression , present_calls, m='s'), p=cor.test(Expression , present_calls, m='s')$p.value)
median(result$cor)
