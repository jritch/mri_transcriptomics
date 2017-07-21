library(magrittr)
library(ggplot2)
library(readr)
library(dplyr)

allDonors <- NULL
geneOfInterest <- "FLJ23867" #

if(Sys.info()['nodename'] == "RES-C02RF0T2.local") {
  single_gene_folder <- "/Users/lfrench/Desktop/results/mri_transcriptomics/single_gene_data_avg/"
} else {
  single_gene_folder <- "/Users/jritchie/git-repos/mri_transcriptomics/single_gene_data_avg/"
}

for (file in list.files(single_gene_folder, pattern=paste0("[.]", geneOfInterest,"[.]"), full.names = T)) {
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

correlationSummary <- allDonors %>% group_by(donor) %>% summarize(cor = cor(MRI_Intensity , Expression, m='s'), p=cor.test(MRI_Intensity , Expression, m='s')$p.value)
correlationSummary$p <- signif(correlationSummary$p, digits=2)
correlationSummary$cor <- round(correlationSummary$cor, digits=2)
correlationSummary$label <- paste0("\n  ρ = ", correlationSummary$cor, "  \n  p-value = ", correlationSummary$p,"\n")

#linear model for exploration
summary(lm(data=allDonors, MRI_Intensity ~ Expression + cortical_division + donor ))

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

ggplot(allDonors, aes(x=MRI_Intensity, y = Expression)) + geom_point(alpha=0.6, aes(color = cortical_division)) + geom_smooth(method = 'loess')  +
  ylab(paste(geneOfInterest, "Expression")) + xlab("T1-/T2-w Ratio") + labs(color="Cortical Division") + 
  geom_text(data = correlationSummary, aes(label=label), x=-Inf, y=yLegend, hjust=0, vjust=vjustLegend, size = 3.5) +
  facet_wrap(~ donor, scales="free")+ theme_bw()

unique(allDonors$donor)
singleDonor <- allDonors %>% filter(donor=="Donor 9861" )
singleCorrelationSummary <- correlationSummary%>% filter(donor=="Donor 9861" )

ggplot(singleDonor, aes(x=MRI_Intensity, y = Expression)) + geom_point(alpha=0.6, aes(color = cortical_division)) + geom_smooth(method = 'loess')  +
  ylab(paste(geneOfInterest, "Expression")) + xlab("T1-/T2-w Ratio") + labs(color="") + 
  geom_text(data = singleCorrelationSummary, aes(label=label), x=-Inf, y=yLegend, hjust=0, vjust=vjustLegend, size = 3.5) + theme_bw() +
  theme(legend.position="bottom") +guides(color=guide_legend(nrow=1,byrow=TRUE)) #guides(color=FALSE) + theme(aspect.ratio=1) + 




#for poster
donor9861 <- allDonors %>% filter(donor == "Donor 9861")
ggplot(donor9861, aes(donor, factor(regionID))) +
  geom_tile(aes(fill = MRI_Intensity), color="black") + theme_void() +
  scale_fill_gradientn(name="", colors=c("#000000", "#FFFFFF"), guide=F)  + theme(plot.margin = unit(c(1,1,1,1), "cm"))
#save as 3x100 PDF in rstudio   
