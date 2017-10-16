library(magrittr)
library(ggplot2)
library(readr)
library(dplyr)

allDonors <- NULL
#geneOfInterest <- "RBP4"
#regionOfInterest <- "cortex_excluding_piriform_hippocampus"
geneOfInterest <- "HLA-E"
regionOfInterest <- "whole_brain"

modalityOfInterest <- "T1T2Ratio"
xLabel <- "T1-/T2-w Ratio"

if (Sys.info()['nodename'] == "RES-C02RF0T2.local") {
  #set working directory
  setwd("/Users/lfrench/Desktop/results/mri_transcriptomics/")
}

single_gene_folder <- "./results/single_gene_data/"

for (file in list.files(single_gene_folder, pattern=paste0("[.]", geneOfInterest,"[.]",regionOfInterest, ".*", modalityOfInterest, "[.]"), full.names = T)) {
  print(file)
  oneDonor <- as.data.frame(read_csv(file))
  oneDonor$Expression <- oneDonor[,geneOfInterest]
  oneDonor$donor <- paste("Donor", gsub("[.].*","", basename(file)))
  allDonors <- bind_rows(allDonors, oneDonor)
}
allDonors <- tbl_df(allDonors)
#print out regions with bad values
filter(allDonors, is.nan(MRI_Intensity) | is.infinite(MRI_Intensity))
#remove those points
allDonors %<>% filter(!(is.nan(MRI_Intensity)| is.infinite(MRI_Intensity)))

allDonors %>% group_by(donor) %>% summarize(n())

correlationSummary <- allDonors %>% group_by(donor) %>% summarize(n = n(), cor = cor(MRI_Intensity , Expression, m='s'), p=cor.test(MRI_Intensity , Expression, m='s')$p.value)
correlationSummary$p <- signif(correlationSummary$p, digits=2)
correlationSummary$cor <- round(correlationSummary$cor, digits=2)
correlationSummary$label <- paste0("\n  n = ", correlationSummary$n,"\n  ρ = ", correlationSummary$cor, "  \n  p-value = ", correlationSummary$p,"\n")

#get rid of p=0
correlationSummary %<>% dplyr::mutate(label = gsub("p-value = 0\n", paste0("p-value < ", min(setdiff(p,0)), "\n") , label ))


#linear model for exploration
summary(lm(data=allDonors, MRI_Intensity ~ Expression + cortical_division + donor ))
summary(lm(data=allDonors, MRI_Intensity ~ Expression  + donor ))

#find the best spot for the correlatoin text
if ( median(correlationSummary$cor) < 0) {
  yLegend = -Inf
  vjustLegend=0
} else {
  yLegend = +Inf
  vjustLegend=1
}

if(geneOfInterest == "CAT") {
  yLegend = -Inf
  vjustLegend=0
}

old_id <- c(14380,15496,15697,9861,10021,12876)
old_id <- paste("Donor",as.character(old_id))
new_id <- c("H0351.1012","H0351.1015","H0351.1016","H0351.2001","H0351.2002","H0351.1009")
donor_id_mapping <- data.frame(old_id,new_id,stringsAsFactors = F)

allDonors %<>% left_join(donor_id_mapping,by=c("donor" = "old_id"))
correlationSummary %<>% left_join(donor_id_mapping,by=c("donor" = "old_id"))

plot(1) #to clear plot area
ggplot(allDonors, aes(x=MRI_Intensity, y = Expression)) + geom_point(alpha=0.6, aes(color = cortical_division)) + geom_smooth(method = 'loess')  +
  ylab(paste(geneOfInterest, "Expression")) + xlab(xLabel) + labs(color="Cortical Division") + 
  geom_text(data = correlationSummary, aes(label=label), x=-Inf, y=yLegend, hjust=0, vjust=vjustLegend, size = 3.5) +
  facet_wrap(~ new_id, scales="free")+ theme_bw()
#11x6 inch pdf

#for single donor/figures
unique(allDonors$donor)
singleDonor <- allDonors %>% filter(donor=="Donor 9861" )
singleCorrelationSummary <- correlationSummary%>% filter(donor=="Donor 9861" )

#### If p is too small, look it up from RunSingleGO.AUROC.Analysis.R and hard-code ###

if (singleCorrelationSummary$p == 0) {
  hard_coded_p <- 4.79 * 10 ^-17
  singleCorrelationSummary$p <- hard_coded_p
  singleCorrelationSummary$label <- paste0("\n  ρ = ", singleCorrelationSummary$cor, "  \n  p-value = ", singleCorrelationSummary$p,"\n")
  
}

#plot single
ggplot(singleDonor, aes(x=MRI_Intensity, y = Expression)) + geom_point(alpha=0.6, aes(color = cortical_division)) + geom_smooth(method = 'loess')  +
  ylab(paste(geneOfInterest, "Expression")) + xlab(xLabel) + labs(color="") + 
  geom_text(data = singleCorrelationSummary, aes(label=label), x=-Inf, y=yLegend, hjust=0, vjust=vjustLegend, size = 3.5) + theme_bw() +
  theme(legend.position="bottom") +guides(color=guide_legend(nrow=1,byrow=TRUE)) + guides(color=FALSE) + theme(aspect.ratio=1) 
#use 3.5 inches pdf for figure

#for figure/poster raster bar
donor9861 <- allDonors %>% filter(donor == "Donor 9861")
ggplot(donor9861, aes(donor, factor(regionID))) +
  geom_tile(aes(fill = MRI_Intensity), color="black") + theme_void() +
  scale_fill_gradientn(name="", colors=c("#000000", "#FFFFFF"), guide=F)  + theme(plot.margin = unit(c(1,1,1,1), "cm"))
#save as 3x100 PDF in rstudio   
