#Leon is using GOSOURCEDATE: 2017-Mar29 (type GO.db to find out)

Sys.info()["nodename"]

#set figshare data folder
figshare_data_folder = "./data/figshare data/"

if(Sys.info()['nodename'] == "RES-C02RF0T2.local") {
  #set working directory
  setwd("/Users/lfrench/Desktop/results/mri_transcriptomics/")
  #filename <-  "/Users/lfrench/Desktop/results/mri_transcriptomics/results/HCP.full_brain.gene_list.csv"
  #filename <-  "/Users/lfrench/Desktop/results/mri_transcriptomics/results/HCP.cortex_excluding_piriform_hippocampus.gene_list.csv"
  #filename <-  "/Users/lfrench/Desktop/results/mri_transcriptomics/results/T1T2Ratio.full_brain.gene_list.csv"
  filename <-  "/Users/lfrench/Desktop/results/mri_transcriptomics/results/T1T2Ratio.cortex_excluding_piriform_hippocampus.gene_list.csv"
} else {
  setwd("C:/Users/Jacob/mri_transcriptomics")
  filename <- "./results/T1T2Ratio.cortex_excluding_limbic_lobe.gene_list.csv"
}

cat(paste("Working directory:", getwd()))
cat(paste("Using input file:",filename))
baseFilename <- gsub(".csv", "", filename)

otherGeneListsFolder <- "./other gene lists/"

library(ggsignif)
library(cowplot)
library(metap)
library(readr)
library(ggplot2)
library(tidyr)
library(magrittr)
library(homologene) #install via install_github('oganm/homologene')
library(org.Hs.eg.db)
library(AnnotationDbi)
library(annotate)
library(GO.db)
library(tmod)
library(metap)
library(reshape2)
library(dplyr)
library(readxl)

geneStatistics <- read_csv(filename, guess_max = 25000) 
cor(geneStatistics[,3:7])
range(geneStatistics[,3:7])
geneStatistics %>% group_by(ID) %>% summarize(ids = length(unique(entrez_id)))
length(unique(geneStatistics$entrez_id))
nrow(geneStatistics %>% filter(is.na(entrez_id)))

#update gene symbols
symbol_to_entrez <- geneStatistics %>% dplyr::select(ID, entrez_id) %>% filter(!is.nan(entrez_id)) %>% filter(!is.na(entrez_id)) %>% distinct()
symbol_to_entrez %<>% mutate(new_symbol = getSYMBOL(as.character(entrez_id), data='org.Hs.eg')) 

#melt into gene, brain
correlationLongForm <- gather(dplyr::select(geneStatistics, ID, contains("Correlation.")), brain, correlation, -ID)
(correlationLongForm %<>% mutate(brain = gsub("Correlation.", "", brain)))
pvaluesLongForm <- gather(dplyr::select(geneStatistics, ID,  contains("PValue.")), brain, pvalue, -ID)
(pvaluesLongForm %<>% mutate(brain = gsub("PValue.", "", brain)))
longFormGeneStats <- inner_join(pvaluesLongForm, correlationLongForm)

#this handling of p-values changed to reduce noise - thresholds changed below
longFormGeneStats$pvalue.neg <- two2one(longFormGeneStats$pvalue, invert =  1 == sign(longFormGeneStats$correlation)) 
longFormGeneStats$pvalue.pos <- two2one(longFormGeneStats$pvalue, invert = 1 == sign(-1*longFormGeneStats$correlation)) 


(geneStatistics <- longFormGeneStats %>% group_by(ID) %>% 
    summarise(medianCorrelation = median(correlation), 
              directionSum = (sum(sign(correlation))/2 +3)/6, 
              metaP.pos = sumlog(pvalue.pos)$p, 
              metaP.neg = sumlog(pvalue.neg)$p, 
              metaP.anydirection = sumlog(pvalue)$p))

geneStatistics %<>% dplyr::rename(geneSymbol = ID)
#filter custom and unmapped probes
geneStatistics <- geneStatistics %>% filter(!grepl("A_", geneSymbol)) %>% filter(!grepl("CUST_", geneSymbol)) 

#adjust p-values
geneStatistics %<>% mutate(metaP.pos.adj = p.adjust(metaP.pos, method="holm"))
geneStatistics %<>% mutate(metaP.neg.adj = p.adjust(metaP.neg, method="holm"))

(geneStatistics %<>% arrange(metaP.pos))
(geneStatistics %<>% arrange(metaP.neg))

range(geneStatistics$medianCorrelation)

sum(geneStatistics$metaP.pos.adj < 0.025)
sum(geneStatistics$metaP.neg.adj < 0.025)
sum(geneStatistics$metaP.neg.adj < 0.025 | geneStatistics$metaP.pos.adj < 0.025)

geneStatistics %<>% dplyr::mutate(isSig = metaP.neg.adj < 0.025 | metaP.pos.adj < 0.025)

geneStatistics %<>% mutate(pValueWithDirection = if_else(metaP.pos < metaP.neg, nrow(geneStatistics) - rank(metaP.pos), -1* nrow(geneStatistics) + rank(metaP.neg)))

freeSurferData <- read_tsv(paste0(figshare_data_folder, "AllenHBA_DK_ExpressionMatrix.tsv"))
freeSurferData <- tbl_df(melt(freeSurferData))
#add in median donor correlation
geneStatistics %<>% inner_join(filter(freeSurferData, variable == 'Average donor correlation to median') %>% dplyr::select(geneSymbol = X1, DonorCorrelation = value))

#add in average expression level (left hemisphere)
averageExpression <- filter(freeSurferData, grepl("ctx-lh-", variable)) %>% group_by(X1) %>% summarize(averageExpression = mean(value)) %>% dplyr::rename(geneSymbol = X1)
geneStatistics %<>% inner_join(averageExpression)

#significant in both directions - should be none
dplyr::filter(geneStatistics, metaP.neg.adj < 0.05 & metaP.pos.adj < 0.05)



paste("Genes with negative correlations:", dplyr::filter(geneStatistics, metaP.neg.adj < 0.05) %>% summarize(n = n()))
paste("Genes with positive correlations:", dplyr::filter(geneStatistics, metaP.pos.adj < 0.05) %>% summarize(n = n()))


#average expression for neg and positive
geneStatistics %>% filter(isSig) %>% group_by(medianCorrelation > 0 ,isSig ) %>% summarize(n=n(), averageExpressionMean = mean(averageExpression,na.rm=TRUE), averageExpressionSD = sd(averageExpression,na.rm=TRUE))
wilcox.test(filter(geneStatistics, isSig, medianCorrelation > 0)$averageExpression, filter(geneStatistics, isSig, medianCorrelation < 0)$averageExpression)

geneStatistics %>% group_by( metaP.neg.adj > metaP.pos.adj) %>% filter(isSig) %>% summarize(n=n())
#percent agree
geneStatistics %>% filter(isSig) %>% summarize(n=n(), agreement = sum(directionSum %% 1 == 0 ), proportion = agreement/n)
geneStatistics %>% group_by(metaP.neg.adj > metaP.pos.adj) %>% filter(isSig) %>% summarize(n=n(), agreement = sum(directionSum %% 1 == 0 ), proportion = agreement/n)

#agree in 5 of 6
geneStatistics %>% filter(isSig) %>% summarize(n=n(), agreement = sum(abs(directionSum -.5) > 0.33), proportion = agreement/n)

#average agreement - used to judge probe strategy
geneStatistics %>% summarize(n=n(), meanAgreement = sum(abs(directionSum -.5)), proportion = meanAgreement/n, donorCount = proportion*6+3)

#average agreement - used to judge probe strategy
geneStatistics %>% filter(isSig)%>% summarize(n=n(), meanAgreement = sum(abs(directionSum -.5)), proportion = meanAgreement/n, donorCount = proportion*6+3)

#update symbols (needs to be done after merging freesufer info)
geneStatistics <- left_join(geneStatistics, dplyr::rename(symbol_to_entrez, geneSymbol=ID))
geneStatistics %<>% mutate(legacySymbol = geneSymbol) #save old symbol
geneStatistics %<>% mutate(geneSymbol = if_else(is.na(new_symbol), geneSymbol, new_symbol)) %>% dplyr::select(geneSymbol, entrez_id,  everything(), -new_symbol)

#add gene names
geneStatistics %<>% mutate(name = unlist(lookUp(as.character(entrez_id), "org.Hs.eg", "GENENAME"))) %>% dplyr::select(geneSymbol, name, everything())

#sort 
geneStatistics <- arrange(geneStatistics, desc(pValueWithDirection))

write_csv(geneStatistics, paste0(baseFilename, ".addedStats.csv"))

sortedGenes <- geneStatistics$geneSymbol

######################################################
# AUC via tmod
######################################################

if (exists("geneSetsGO") && length(geneSetsGO$MODULES2GENES) > 1000 ) { #assume it's already loaded - needs a fix to see if the variable is declared
} else {
  go_object <- as.list(org.Hs.egGO2ALLEGS)
  
  symbolsInGO <- getSYMBOL(unique(unlist(go_object)), data='org.Hs.eg')
  
  #build GO sets for tmod -slow
  tmodNames <- data.frame()
  modules2genes <- list()
  goGroupName <- names(go_object)[1]
  showMethods(Term)
  
  goCount <- length(go_object)
  count <- 1
  for(goGroupName in names(go_object)) {
    if (count %% 1000 == 0) print(paste(count, "of", goCount))
    count <- count + 1
    
    goGroup <- go_object[goGroupName]
    geneIDs <- unique(unlist(goGroup, use.names=F))  #discard evidence codes
    genesymbols <- unique(getSYMBOL(geneIDs, data='org.Hs.eg'))
    
    genesymbols <- intersect(genesymbols, sortedGenes)
    if (!(length(genesymbols) >= 10 & length(genesymbols) <= 200)) next();
    
    modules2genes[goGroupName] <- list(genesymbols)
    
    tmodNames <- rbind(tmodNames, data.frame(ID=goGroupName, Title = Term(goGroupName)))
  }
  geneSetsGO <- makeTmod(modules = tmodNames, modules2genes = modules2genes)
}

result <- tbl_df(tmodUtest(c(sortedGenes), mset=geneSetsGO, qval = 1, filter = T))
result %<>% rowwise() %>% mutate(P.Value = P.Value * 2) %>% ungroup() %>% dplyr::select(-adj.P.Val) 
result %<>% rowwise() %>% mutate(aspect = Ontology(ID))
(result %<>% group_by(U, N1, AUC, P.Value) %>% summarize(MainTitle = first(Title),  ID=paste(ID, collapse=","), aspect= first(aspect), allNames = if_else(n() > 1, paste(Title[2:length(Title)], collapse=","), "")))
result %<>% ungroup() %>% mutate(adj.P.Value=p.adjust(P.Value, method = "holm")) %>% arrange(P.Value) #adjust again
result$rank <- 1:nrow(result)
result %<>% dplyr::select(MainTitle, geneCount = N1, AUC, P.Value, adj.P.Value, everything(), -U) 

sum(result$adj.P.Value < 0.05)

head(filter(result, AUC > 0.5),n=15)
head(filter(result, AUC < 0.5),n=15) %>% select(MainTitle)

result$adj.P.Value <- signif(result$adj.P.Value, digits=3)
result$AUC <- signif(result$AUC, digits=3)
result$P.Value <- signif(result$P.Value, digits=3)

# Output tables for top ten positively and negatively enriched GO groups for paper
(result.up <- head(filter(result, AUC > 0.5) %>% dplyr::select(-ID), n=10))
(result.down <- head(filter(result, AUC < 0.5) %>% dplyr::select(-ID),n=15))

#check one p-value to test TMOD
setGO <- filter(geneStatistics, geneSymbol %in% geneSetsGO$MODULES2GENES$`GO:0006521`)
setAll <- filter(geneStatistics, geneSymbol %in% unlist(geneSetsGO$MODULES2GENES)) %>% filter(!geneSymbol %in% setGO$geneSymbol)
wilcox.test( setGO$pValueWithDirection, setAll$pValueWithDirection)

write_csv( result,  paste(baseFilename,".GO.results.csv",sep=""))
write_csv( dplyr::select(result.up, Name = MainTitle,`Gene Count` = geneCount, AUROC = AUC,  `p` = P.Value, `pFWER` = adj.P.Value, aspect),  paste0(baseFilename,".GO.up10.csv"))
write_csv( dplyr::select(result.down, Name = MainTitle,`Gene Count` = geneCount, AUROC = AUC, `p` = P.Value, `pFWER` = adj.P.Value, aspect),  paste0(baseFilename,".GO.down15.csv"))


##### Look at what proportion of the results match certain categories

cat(paste("Total number of GO groups tested (after removing dupes):",nrow(result)))

cat(paste("Number of significant GO groups",as.character(lengths(result %>% filter (adj.P.Value < 0.05))[1])))

cat(paste("Number of significant positively enriched GO groups",as.character(lengths(result %>% filter (adj.P.Value < 0.05, AUC > 0.5))[1])))

cat(paste("Number of significant negatively enriched GO groups",as.character(lengths(result %>% filter (adj.P.Value < 0.05, AUC < 0.5))[1])))

cat(paste("Number of significant negatively enriched BP GO groups related to synapse",
          as.character(lengths(result %>% filter (adj.P.Value < 0.05, AUC < 0.5,aspect=="BP", grepl("synap",MainTitle)))[1])))
cat(paste("synapse groups tested: ", as.character(lengths(result %>% filter (aspect=="BP", grepl("synap",MainTitle)))[1])))

source("./R Code/ROCPlots.R")

dplyr::filter(tbl_df(geneSetsGO$MODULES), Title == "detection of chemical stimulus involved in sensory perception")$ID
dplyr::filter(tbl_df(geneSetsGO$MODULES), Title == "transcriptional repressor activity, RNA polymerase II transcription regulatory region sequence-specific binding")$ID

#plots <- createPlots(sortedGenes, c("GO:0098798", "GO:1905368", "GO:0005882","GO:0031490","GO:0005814","GO:0033038","GO:0032452", "GO:0001227", "GO:0045178", "GO:0001047"), geneSetsGO)

plots <- createPlots(sortedGenes, c("GO:0050907", "GO:0098798", "GO:1905368", "GO:0033038","GO:0006521"), geneSetsGO)
(bothPlots <- plot_grid(plots$AUCPlot, plots$rasterPlot, nrow = 2, align = "v", rel_heights=c(1,0.6),scale = 0.95,labels = c("A", "B"))) #add labels = c("A", "B"), for manuscript
#save as 11x11 PDF

#glycolysis
filter(result,grepl("glycolysis",allNames))
filter(result,grepl("pyruvate",allNames) | grepl("pyruvate",MainTitle))
filter(result,grepl("lactate",allNames) | grepl("lactate",MainTitle))


#filter(result,grepl("myelin",MainTitle))
myelinResult <- filter(result, grepl("myelin|ensheathment",MainTitle),!grepl("peripheral|sphingomyelin",MainTitle))
myelinResult$adj.P.Value <- p.adjust(myelinResult$P.Value, method="fdr")
myelinResult$adj.P.Value <- signif(myelinResult$adj.P.Value, digits=3)

myelinResult
length(intersect(geneSetsGO$MODULES2GENES$`GO:0007272`, geneSetsGO$MODULES2GENES$`GO:0042552`))
write_csv( dplyr::select(myelinResult, Name = MainTitle,`Gene Count` = geneCount, AUROC = AUC,  p=P.Value, `pFDR` = adj.P.Value, aspect, synonyms = allNames, rank),  
           paste0(baseFilename,".GO.myelin.results.csv"))

plots <- createPlots(sortedGenes, c("GO:0008366","GO:0043209", "GO:0043217"), geneSetsGO)
(bothPlots <- plot_grid(plots$AUCPlot, plots$rasterPlot, nrow = 2, align = "v", rel_heights=c(1,0.6),scale = 0.95, labels = c("A", "B"))) #add labels = c("A", "B"), for manuscript
#save as 10x7



#################################################################
# cell type gene lists and other lists
#################################################################
tmodNames <- data.frame()
modules2genes <- list()


for(geneListFilename in list.files(otherGeneListsFolder, pattern = ".*txt", full.names = T)) {
  print(geneListFilename)
  genesOfInterest <- read.csv(geneListFilename,header=F,stringsAsFactors = F)
  shortName <- gsub(".txt","",gsub(paste0(".*/"),"", geneListFilename))
  
  genesOfInterest$term <- shortName
  
  #already a human gene list
  if (grepl(pattern = "Darmanis.", geneListFilename  ) | grepl(pattern = "Mistry", geneListFilename) | grepl(pattern = "HouseKeeping", geneListFilename  ) | grepl(pattern = "human", geneListFilename  )) {
    modules2genes[shortName] <- list(genesOfInterest$V1)
  } else { #needs conversion from mouse
    print(" converting from mouse to human")
    modules2genes[shortName] <- list(mouse2human(genesOfInterest$V1)$humanGene)
  }

  tmodNames <- rbind(tmodNames, data.frame(ID=shortName, Title = shortName))
}
geneSets <- makeTmod(modules = tmodNames, modules2genes = modules2genes)
geneSetsCellType <- geneSets #for later reuse

result <- tmodUtest(sortedGenes, mset=geneSets, qval = 1, filter = F)
result <- tbl_df(result) %>% dplyr::select(Title, geneCount =N1,AUC,  P.Value, adj.P.Val, -ID)
result %<>% rowwise() %>% mutate(P.Value = P.Value * 2) %>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value, method="holm")) #tmod runs one-sided tests

result$adj.P.Val <- signif(result$adj.P.Val, digits=3)
result$AUC <- signif(result$AUC, digits=3)

writeTableWrapper <- function(prefixFilter, result) {
  subsetResult <- filter(result, grepl(prefixFilter, Title))
  subsetResult$oldTitle <- subsetResult$Title
  subsetResult$Title <- gsub(paste0(prefixFilter, "."), "", subsetResult$Title)
  subsetResult$Title <- gsub("[.]", " ", subsetResult$Title)
  subsetResult$Title <- gsub("[_]", " ", subsetResult$Title)
  subsetResult$Title <- gsub("Neuron interneuron", "Interneuron", subsetResult$Title)
  subsetResult$Title <- gsub("ligo", "ligodendrocyte", subsetResult$Title)
  subsetResult$adj.P.Val <- p.adjust(subsetResult$P.Value, method="holm")
  subsetResult$adj.P.Val <- signif(subsetResult$adj.P.Val, digits=3)
  subsetResult$P.Value <- signif(subsetResult$P.Value, digits=3)
  write_csv(dplyr::select(subsetResult, `Cell-type or class` = Title,`Gene Count` = geneCount, AUROC = AUC, p = P.Value, `pFWER` = adj.P.Val), paste0(baseFilename,".", prefixFilter,".csv")) 
  subsetResult
}

(zeiselResult <- writeTableWrapper("Zeisel", result))
(darmResult <- writeTableWrapper("Darmanis", result))
writeTableWrapper("NeuroExpresso.Cortex", result)
(mistryResult <- writeTableWrapper("Mistry", result))

writeTableWrapper("Thakurela", result)
writeTableWrapper("Hametner", result)
writeTableWrapper("TEMP", result)

intersect(head(sortedGenes,n=10),geneSetsCellType$MODULES2GENES$Thakurela.Myelin.cortexCompare)
intersect(head(sortedGenes,n=10),geneSetsCellType$MODULES2GENES$Thakurela.Myelin.proteome)

intersect(tail(sortedGenes,n=10),geneSetsCellType$MODULES2GENES$Thakurela.Myelin.cortexCompare)
intersect(tail(sortedGenes,n=10),geneSetsCellType$MODULES2GENES$Thakurela.Myelin.proteome)


plots <- createPlots(sortedGenes, as.character(zeiselResult$oldTitle), geneSets, customNames=zeiselResult$Title, filter=F)
(bothPlots <- plot_grid(plots$AUCPlot, plots$rasterPlot, nrow = 2, align = "v", rel_heights=c(1,0.9),scale = 0.95,labels = c("A", "B"))) #add labels = c("A", "B"), for manuscript
#save as 11x11


plots <- createPlots(sortedGenes, as.character(darmResult$oldTitle), geneSets, customNames=darmResult$Title, filter=F)
(bothPlots <- plot_grid(plots$AUCPlot, plots$rasterPlot, nrow = 2, align = "v", rel_heights=c(1,0.8), labels = c("A", "B"))) #add labels = c("A", "B"), for manuscript
plots$rasterPlot #save as 10x4 pdf
#filter(geneStatistics, geneSymbol %in% geneSets["Darmanis.Oligo"]$GENES$ID)

#plot Mistry genes - for slides
plots <- createPlots(sortedGenes, as.character(mistryResult$oldTitle), geneSets, customNames=mistryResult$Title, filter=F)
(bothPlots <- plot_grid(plots$AUCPlot, plots$rasterPlot, nrow = 2, align = "v", rel_heights=c(1,0.8),scale = 0.95)) #add labels = c("A", "B"), for manuscript


##################
# Zeng et al. 
################

zengPath <- "./other gene lists/ZengEtAl/1-s2.0-S0092867412003480-mmc2.xlsx"

zengTable <- read_xlsx(path = zengPath, sheet = "Final1000New", skip=1)

zengTable %<>% select(`Gene symbol`, `Cortical marker (human)`) %>% rename(Gene.symbol = `Gene symbol`, Cortical.marker..human.=`Cortical marker (human)`)

backGroundGenes <- zengTable$Gene.symbol

expandedZengTable <- zengTable %>% 
  mutate(Cortical.marker..human. = strsplit(as.character(Cortical.marker..human.), "[/+]|( or )")) %>% 
  unnest(Cortical.marker..human.)

as.data.frame(expandedZengTable)

expandedZengTable %<>% mutate(Cortical.marker..human. = gsub("layer( )?","", Cortical.marker..human.))
expandedZengTable %<>% mutate(Cortical.marker..human. = gsub("[?]","", Cortical.marker..human.))
expandedZengTable %<>% mutate(Cortical.marker..human. = gsub("4c","4", Cortical.marker..human.))
expandedZengTable %<>% mutate(Cortical.marker..human. = gsub("5a","5", Cortical.marker..human.))
expandedZengTable %<>% mutate(Cortical.marker..human. = gsub("6b","6", Cortical.marker..human.))
expandedZengTable %<>% mutate(Cortical.marker..human. = gsub("([0-6])","layer \\1", Cortical.marker..human.))
expandedZengTable %<>% mutate(Cortical.marker..human. = gsub("VEC","vascular endothelial cell", Cortical.marker..human.))

expandedZengTable %<>% filter(Cortical.marker..human. != 'others' & Cortical.marker..human. != "laminar")
expandedZengTable %>% group_by(Cortical.marker..human.) %>% summarise(n())

tmodNames <- data.frame()
modules2genes <- list()

for(geneSetName in unique(expandedZengTable$Cortical.marker..human.)) {
  print(geneSetName)
  
  genesOfInterest <- filter(expandedZengTable, Cortical.marker..human. == geneSetName) %>% .$Gene.symbol
  shortName <- geneSetName
  
  modules2genes[shortName] <- list(genesOfInterest)

  tmodNames <- rbind(tmodNames, data.frame(ID=shortName, Title = shortName))
}

geneSets <- makeTmod(modules = tmodNames, modules2genes = modules2genes)
sortedGenesInBackground <- sortedGenes[sortedGenes %in% backGroundGenes]

#for the whole Zeng list, is it biased?
mean(filter(geneStatistics, geneSymbol %in% sortedGenesInBackground)$medianCorrelation)

zeng_result <- tmodUtest(sortedGenesInBackground, mset=geneSets, qval = 1, filter = F)
zeng_result <- tbl_df(zeng_result) %>% dplyr::select(Title, geneCount = N1, AUC, P.Value, adj.P.Val, ID)
zeng_result %<>% rowwise() %>% mutate(P.Value = P.Value * 2) %>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value, method="holm")) #tmod runs one-sided tests
zeng_result %<>% arrange(as.character(Title))

filter(geneStatistics, geneSymbol %in% expandedZengTable$Gene.symbol)

write_csv(dplyr::select(zeng_result, `Type` = Title,`Gene Count` = geneCount, AUROC = AUC,  `pFWER` = adj.P.Val), paste0(baseFilename,".ZengEtAl.csv")) 

zeng_result$Title <- factor(zeng_result$Title, levels=rev(sort(as.character(unique(zeng_result$Title)))))
(barplot <- ggplot(zeng_result, aes(y=AUC-0.5, x=Title, fill=AUC > 0.5)) + geom_col() + 
  coord_flip() + xlab("") + ylab("AUROC") + guides(fill=FALSE) +
  scale_y_continuous(breaks=seq(-0.9, .9, by= .1), labels=seq(-0.9, .9, by= .1) + .5))
barplot <- barplot + geom_text(data = filter(zeng_result, adj.P.Val  < 0.05 & adj.P.Val > 0.005), aes(hjust = -.5*sign(AUC-0.5)+.5), label = "*")
barplot <- barplot + geom_text(data = filter(zeng_result, adj.P.Val < 0.005), aes(hjust = -.5*sign(AUC-0.5)+.5), label = "**")
barplot
ggsave(plot= barplot,paste0(baseFilename,".ZengEtAl.pdf" ), height=5, width=5)
    

################################################
#### He Z et al. cortical layer markers (not currently used in manuscript)
#    https://www.ncbi.nlm.nih.gov/pubmed/28414332
#    http://www.picb.ac.cn/Comparative/data_methods/data_layer_2017.html
################################################
HeZPath <- "./other gene lists/HeEtAl/table_s2.xlsx"

HeZTable <- read_xlsx(path = HeZPath, sheet = "Sheet1")
HeZTable %<>% rename(Gene.symbol = `Gene symbol`, Layer.marker.in.human=`Layer marker in human`) %>% dplyr::select(Gene.symbol, Layer.marker.in.human)

backGroundGenes <- unique(HeZTable$Gene.symbol)

tmodNames <- data.frame()
modules2genes <- list()

for(geneSetName in unique(HeZTable$Layer.marker.in.human)) {
  if(geneSetName == "NA") next
  print(geneSetName)
  
  genesOfInterest <- filter(HeZTable, Layer.marker.in.human == geneSetName) %>% .$Gene.symbol
  shortName <- geneSetName
  shortName <- gsub("L", "Layer ", shortName)
  
  modules2genes[shortName] <- list(genesOfInterest)
  
  tmodNames <- rbind(tmodNames, data.frame(ID=shortName, Title = shortName))
}
geneSets <- makeTmod(modules = tmodNames, modules2genes = modules2genes)
sortedGenesInBackground <- sortedGenes[sortedGenes %in% backGroundGenes]

HeZ_result <- tmodUtest(sortedGenesInBackground, mset=geneSets, qval = 1, filter = F)
HeZ_result <- tbl_df(HeZ_result) %>% dplyr::select(Title, geneCount = N1, AUC, P.Value, adj.P.Val, ID)
HeZ_result %<>% rowwise() %>% mutate(P.Value = P.Value * 2) %>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value, method="holm")) #tmod runs one-sided tests
HeZ_result %<>% arrange(as.character(Title))

HeZ_result %>% filter(adj.P.Val < 0.05, AUC>0.5)
HeZ_result %>% filter(adj.P.Val < 0.05, AUC<0.5)


HeZ_result$Title <- factor(HeZ_result$Title, levels=rev(sort(as.character(unique(HeZ_result$Title)))))
(barplot <- ggplot(HeZ_result, aes(y=AUC-0.5, x=Title, fill=AUC > 0.5)) + geom_col() + 
    coord_flip() + xlab("") + ylab("AUROC") + guides(fill=FALSE) +
    scale_y_continuous(breaks=seq(-0.9, .9, by= .1), labels=seq(-0.9, .9, by= .1) + .5))
barplot <- barplot + geom_text(data = filter(HeZ_result, adj.P.Val  < 0.05 & adj.P.Val > 0.005), aes(hjust = -.5*sign(AUC-0.5)+.5), label = "*")
barplot <- barplot + geom_text(data = filter(HeZ_result, adj.P.Val < 0.005), aes(hjust = -.5*sign(AUC-0.5)+.5), label = "**")
barplot
ggsave(plot= barplot,paste0(baseFilename,".HeEtAl.pdf" ), height=5, width=5)

plots <- createPlots(sortedGenesInBackground, HeZ_result$Title, geneSets, filter=F, useGroupGivenOrder = T)
plots$AUCPlot
plots$rasterPlot

HeZ_result$adj.P.Val <- signif(HeZ_result$adj.P.Val, digits=3)
HeZ_result$AUC <- signif(HeZ_result$AUC, digits=3)

write_csv(dplyr::select(HeZ_result, `Marker Set` = Title,`Gene Count` = geneCount, AUROC = AUC,  `pFWER` = adj.P.Val), paste0(baseFilename,".HeZEtAl.csv")) 

#################################################################
#################################################################

loadPhenocarta <- function(taxon, geneBackground) {
  #guess column types from the whole dataset, basically
  phenocarta <- read_tsv(paste0(figshare_data_folder, "AllPhenocartaAnnotations.downloadedOct28.2016.tsv"), skip = 4, guess_max = 130000)
  phenocarta$ID <- gsub("http://purl.obolibrary.org/obo/", "", phenocarta$`Phenotype URIs`)
  phenocarta <- dplyr::filter(phenocarta, Taxon == taxon) %>% dplyr::select(symbol = `Gene Symbol`, name = `Phenotype Names`, ID) %>% filter(symbol %in% geneBackground) %>% distinct()
  geneLists <- group_by(phenocarta, ID) %>% dplyr::summarise(name = paste(unique(name), collapse = ","), genes = unique(list(symbol)), size = n()) %>% filter(size > 5 & size < 200) 
  distinct(geneLists)
  namedLists <- geneLists$genes
  names(namedLists) <- geneLists$ID
  idToName <- data.frame(ID = geneLists$ID, Title = geneLists$name)
  geneSets <- makeTmod(modules = idToName, modules2genes = namedLists)
  geneSets
}

geneSetsPhenoCarta <- loadPhenocarta("human", sortedGenes)

result <- tmodUtest(c(sortedGenes), mset=geneSetsPhenoCarta, qval = 1, filter = F)
result %<>% rowwise() %>% mutate(P.Value = P.Value * 2) %>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value, method="holm")) #tmod runs one-sided tests
result <- tbl_df(result) %>% dplyr::select(ID, Title, geneCount =N1,AUC,  P.Value, adj.P.Val)

result$adj.P.Val <- signif(result$adj.P.Val, digits=3)
result$AUC <- signif(result$AUC, digits=3)
result
write_csv( dplyr::select(result, Disease = Title,`Gene Count` = geneCount, AUROC = AUC, pvalue = P.Value, `pFWER` = adj.P.Val),  
           paste0(baseFilename,".phenocarta.results.csv"))

#plots <- createPlots(sortedGenes, c("DOID_9008", "DOID_3213", "DOID_4233", "DOID_9975"), geneSets)
#(bothPlots <- plot_grid(plots$AUCPlot, plots$rasterPlot, nrow = 2, align = "v", rel_heights=c(1,0.8),scale = 0.95)) #add labels = c("A", "B"), for manuscript

#looking at one result - epilepsy
#filter(geneStatistics, geneSymbol %in% geneSetsPhenoCarta["DOID_1932"]$GENES$ID)

################
#DisGeNET - inspired by Burt et al. usage (not currently used in manuscript)
###############

disgenet <- read_tsv("./other gene lists/DisGeNET/curated_gene_disease_associations.tsv.gz")
disgenet %<>% dplyr::select(symbol = geneSymbol, name = diseaseName, ID = diseaseId)
  
geneLists <- group_by(disgenet, ID) %>% dplyr::summarise(name = paste(unique(name), collapse = ","), genes = unique(list(symbol)), size = n()) %>% filter(size > 5 & size < 200) 
distinct(geneLists)
namedLists <- geneLists$genes
names(namedLists) <- geneLists$ID
idToName <- data.frame(ID = geneLists$ID, Title = geneLists$name)
geneSets <- makeTmod(modules = idToName, modules2genes = namedLists)

result <- tmodUtest(sortedGenes, mset=geneSets, qval = 1, filter = F)
result <- tbl_df(result) %>% dplyr::select(Title, geneCount =N1,AUC,  P.Value, adj.P.Val, -ID)
result %<>% rowwise() %>% mutate(P.Value = P.Value * 2) %>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value, method="holm")) #tmod runs one-sided tests


result %>% filter(adj.P.Val < 0.05, AUC>0.5)
result %>% filter(adj.P.Val < 0.05, AUC<0.5)
result
write_csv( dplyr::select(result, Disease = Title,`Gene Count` = geneCount, AUROC = AUC,  `pFWER` = adj.P.Val),  
           paste0(baseFilename,".disgenet.results.csv"))


################
# Spaethling et al. lists (not currently used in manuscript, agrees with Darmanis)
# Primary Cell Culture of Live Neurosurgically Resected Aged Adult Human Brain Cells and Single Cell Transcriptomics
###############
SpaethlingTable <- read.xlsx("./other gene lists/Spaethling et al./1-s2.0-S2211124716317739-mmc3.xlsx", sheetName = "Sheet1", stringsAsFactors = F)
(SpaethlingTable <- tbl_df(SpaethlingTable) %>% dplyr::rename(symbol = Table.S4..Related.to.Figure.3B..Genes.enriched.in.each.cell.type, cellType = NA.) ) 

geneLists <- group_by(SpaethlingTable, cellType) %>% dplyr::summarise(genes = unique(list(symbol)), size = n()) 
namedLists <- geneLists$genes
names(namedLists) <- geneLists$cellType
idToName <- data.frame(ID = geneLists$cellType, Title = geneLists$cellType)
geneSets <- makeTmod(modules = idToName, modules2genes = namedLists)

Spaethling_result <- tmodUtest(sortedGenesInBackground, mset=geneSets, qval = 1, filter = F)
Spaethling_result <- tbl_df(Spaethling_result) %>% dplyr::select(Title, geneCount = N1, AUC, P.Value, adj.P.Val, ID)
Spaethling_result %<>% rowwise() %>% mutate(P.Value = P.Value * 2) %>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value, method="holm")) #tmod runs one-sided tests

Spaethling_result$adj.P.Val <- signif(Spaethling_result$adj.P.Val, digits=3)
Spaethling_result$AUC <- signif(Spaethling_result$AUC, digits=3)
Spaethling_result
write_csv( dplyr::select(Spaethling_result, CellType = Title,`Gene Count` = geneCount, AUROC = AUC,  `pFWER` = adj.P.Val),  
           paste0(baseFilename,".Spaethling.results.csv"))
