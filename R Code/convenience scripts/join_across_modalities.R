library(readr)
T1 <- read_csv("/Users/lfrench/Desktop/results/mri_transcriptomics/results/T1.cortex_excluding_piriform_hippocampus.gene_list.GO.myelin.results.csv")
T2 <- read_csv("/Users/lfrench/Desktop/results/mri_transcriptomics/results/T2.cortex_excluding_piriform_hippocampus.gene_list.GO.myelin.results.csv")
T1T2Ratio <- read_csv("/Users/lfrench/Desktop/results/mri_transcriptomics/results/T1T2Ratio.cortex_excluding_piriform_hippocampus.gene_list.GO.myelin.results.csv")

joined <- inner_join(T1, T2, by='Name', suffix = c(".T1", ".T2")) %>% dplyr::select(Name, rank.T1, rank.T2, `T1 AUROC`= AUROC.T1, `T2 AUROC` = AUROC.T2)
joined <- inner_join(T1T2Ratio, joined, by='Name')
joined %<>% dplyr::select(Name, `Gene Count`, p, pFDR, AUROC, rank, rank.T1, rank.T2) #, `T1 AUROC`, `T2 AUROC`)
joined
write_csv(joined, "/Users/lfrench/Desktop/results/mri_transcriptomics/results/T1T2Ratio.cortex_excluding_piriform_hippocampus.gene_list.GO.myelin.results.addedT1T2.csv")


### AUC for top GO hit table
T1 <- read_csv("/Users/lfrench/Desktop/results/mri_transcriptomics/results/T1.cortex_excluding_piriform_hippocampus.gene_list.GO.results.csv")
T2 <- read_csv("/Users/lfrench/Desktop/results/mri_transcriptomics/results/T2.cortex_excluding_piriform_hippocampus.gene_list.GO.results.csv")
T1T2Ratio <- read_csv("/Users/lfrench/Desktop/results/mri_transcriptomics/results/T1T2Ratio.cortex_excluding_piriform_hippocampus.gene_list.GO.up10.csv")

joined <- inner_join(T1, T2, by='MainTitle', suffix = c(".T1", ".T2")) %>% dplyr::select(Name = MainTitle, `T1 AUROC`= AUC.T1, `T2 AUROC` = AUC.T2)
joined <- inner_join(T1T2Ratio, joined, by='Name')
write_csv(joined, "/Users/lfrench/Desktop/results/mri_transcriptomics/results/T1T2Ratio.cortex_excluding_piriform_hippocampus.gene_list.GO.up10.addedT1T2.csv")

### AUC for top hits
T1 <- read_csv("/Users/lfrench/Desktop/results/mri_transcriptomics/results/T1.cortex_excluding_piriform_hippocampus.gene_list.GO.results.csv")
T2 <- read_csv("/Users/lfrench/Desktop/results/mri_transcriptomics/results/T2.cortex_excluding_piriform_hippocampus.gene_list.GO.results.csv")
T1T2Ratio <- read_csv("/Users/lfrench/Desktop/results/mri_transcriptomics/results/T1T2Ratio.cortex_excluding_piriform_hippocampus.gene_list.GO.down15.csv") %>% dplyr::select(-p)

joined <- inner_join(T1, T2, by='MainTitle', suffix = c(".T1", ".T2")) %>% dplyr::select(Name = MainTitle, AUC.T1, AUC.T2)
joined <- inner_join(T1T2Ratio, joined, by='Name')
write_csv(joined, "/Users/lfrench/Desktop/results/mri_transcriptomics/results/T1T2Ratio.cortex_excluding_piriform_hippocampus.gene_list.GO.down15.addedT1T2.csv")



