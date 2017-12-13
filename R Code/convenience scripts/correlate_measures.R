if (Sys.info()['nodename'] == "RES-C02RF0T2.local") {
  #set working directory
  setwd("/Users/lfrench/Desktop/results/mri_transcriptomics/")
}

source("get_expression_single_gene_function.R")

single_gene_folder <- "./results/single_gene_data/"
regionOfInterest <- "cortex_excluding_piriform_hippocampus"
geneOfInterest <- "SCARA5"

ratio <- get_expression_single_gene(geneOfInterest, regionOfInterest, "T1T2Ratio", single_gene_folder) 
t1 <- get_expression_single_gene(geneOfInterest, regionOfInterest, "T1", single_gene_folder) 
t2 <- get_expression_single_gene(geneOfInterest, regionOfInterest, "T2", single_gene_folder) 
hcp <- get_expression_single_gene(geneOfInterest, regionOfInterest, "HCP", single_gene_folder) 

allDonors %>% group_by(donor) %>% summarize(n(), corHCP=cor(MRI_Intensity,SCARA5,m='s',use="pairwise.complete.obs"))

t1 <- dplyr::select(t1, `(x,y,z)`, donor, t1= MRI_Intensity)
t2 <- dplyr::select(t2, `(x,y,z)`, donor, t2= MRI_Intensity)
ratio <- dplyr::select(ratio, `(x,y,z)`, donor, ratio= MRI_Intensity, Expression)
hcp <- dplyr::select(hcp, `(x,y,z)`, donor, hcp_ratio= MRI_Intensity, Expression)
joined <- inner_join(ratio,hcp, by=c("donor", "Expression")) #probably should join in a safer way

joined %>% group_by(donor) %>% summarize(n(), corHCP=cor(ratio,hcp_ratio,m='s',use="pairwise.complete.obs"))
mean((joined %>% group_by(donor) %>% summarize(n(), corHCP=cor(ratio,hcp_ratio,m='s',use="pairwise.complete.obs")))$corHCP)
#18.53 for cortex correlation mean
ggplot(data= joined, aes(x=ratio, y=hcp_ratio)) + geom_point() + facet_wrap(~donor, scales = "free")

joined <- inner_join(t1, t2)
joined <- inner_join(ratio, joined)
joined %>% group_by(donor) %>% summarize(n(), cort1t2=cor(t1,t2,m='s'),cort1ratio=cor(t1,ratio,m='s'),corratiot2=cor(ratio,t2,m='s'))
