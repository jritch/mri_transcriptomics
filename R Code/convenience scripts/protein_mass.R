#basic script to correlate protein mass with gene ranking - data from biomart and uniprot. Biomart data not used in manuscript. pH(I) was also used (from ExPASy calculator).
library(dplyr)
require(stringr)
library(readr)
library(magrittr)


masses<-read_tsv("/Users/lfrench/Downloads/uniprot-yourlist%3AM20170927B8D4AF2DBA97F39AD84F8C0D29E95A9297D353L.tab.tsv")
pI<-read_tsv("/Users/lfrench/Downloads/pitool.4555.txt", col_names = F) %>% dplyr::rename(pI = X3, massPI=X4, Entry=X2, label = X1)
pI$pI <- as.numeric(pI$pI)
gene_length <- read_tsv("/Users/lfrench/Downloads/mart_export.txt") %>% dplyr::select(geneSymbol = `HGNC symbol`,`Transcript length (including UTRs and CDS)`)

masses$geneSymbol <- gsub("_HUMAN", "", masses$`Entry name`)


masses <- inner_join(masses, pI)

masses %<>% dplyr::select(geneSymbol, proteinLength = Length, proteinMass = Mass, pI)

gene_length %<>% group_by(geneSymbol) %>% summarize(geneLengthMax = max(`Transcript length (including UTRs and CDS)`),geneLengthMean = mean(`Transcript length (including UTRs and CDS)`))

joined <- inner_join(gene_length,masses)
cor.test(joined$geneLengthMax, joined$proteinLength,m='s')
cor.test(joined$geneLengthMean, joined$proteinLength,m='s')
cor.test(joined$geneLengthMean, joined$proteinMass,m='s')
cor.test(joined$geneLengthMax, joined$proteinMass,m='s')

joined <- inner_join(geneStatistics, joined)

cor.test(joined$pValueWithDirection, joined$proteinMass,m='s')

joined %>% filter(isSig==T) %>% group_by(pValueWithDirection>0) %>% summarize(medianMass = median(proteinMass), meanMass = mean(proteinMass), variance = sd(proteinMass))

t.test(filter(joined, isSig, medianCorrelation > 0)$proteinMass, filter(joined, isSig, medianCorrelation < 0)$proteinMass)

wilcox.test(filter(joined, isSig, medianCorrelation > 0)$proteinMass, filter(joined, isSig, medianCorrelation < 0)$proteinMass)
boxplot(filter(joined, isSig, medianCorrelation > 0)$proteinMass, filter(joined, isSig, medianCorrelation < 0)$proteinMass)

cor.test(joined$pI, joined$proteinMass,m='s')

cor.test(joined$pValueWithDirection, joined$pI, m='s') #p < 0.05, wilcox not significant


#protein mass test:
cor.test(joined$pValueWithDirection, joined$proteinMass,m='s')
