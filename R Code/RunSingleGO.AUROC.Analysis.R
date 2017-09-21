#Leon is using GOSOURCEDATE: 2017-Mar29 (type GO.db to find out)

Sys.info()["nodename"]

#set figshare data folder
figshare_data_folder = "./data/figshare data/"

if(Sys.info()['nodename'] == "RES-C02RF0T2.local") {
  #set working directory
  setwd("/Users/lfrench/Desktop/results/mri_transcriptomics/")
  filename <-  "/Users/lfrench/Desktop/results/mri_transcriptomics/results/T1T2Ratio.cortex_excluding_piriform_hippocampus.gene_list.csv"
} else {
  filename <- "/Users/jritchie/data/final/T1T2Ratio.full_brain.gene_list.csv"
}

cat(paste("Working directory:", getwd()))
cat(paste("Using input file:",filename))
baseFilename <- gsub(".csv", "", filename)

otherGeneListsFolder <- "./other gene lists/"

library(ggsignif)
library(cowplot)
library(readr)
library(dplyr)
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
library(xlsx)

geneStatistics <- read_csv(filename) 

#convert to one sided
#melt into gene, brain
correlationLongForm <- gather(dplyr::select(geneStatistics, ID, contains("Correlation.")), brain, correlation, -ID)
(correlationLongForm %<>% mutate(brain = gsub("Correlation.", "", brain)))
pvaluesLongForm <- gather(dplyr::select(geneStatistics, ID,  contains("PValue.")), brain, pvalue, -ID)
(pvaluesLongForm %<>% mutate(brain = gsub("PValue.", "", brain)))
longFormGeneStats <- inner_join(pvaluesLongForm, correlationLongForm)
longFormGeneStats %<>% mutate(pvalue.pos = if_else(sign(correlation) > 0, pvalue, 1 - pvalue ))
longFormGeneStats %<>% mutate(pvalue.neg = if_else(sign(correlation) < 0, pvalue, 1 - pvalue ))
(geneStatistics <- longFormGeneStats %>% group_by(ID ) %>% 
    summarise(medianCorrelation = median(correlation), 
              directionSum = (sum(sign(correlation))/2 +3)/6, 
              metaP.pos = sumlog(pvalue.pos)$p, 
              metaP.neg = sumlog(pvalue.neg)$p, 
              metaP.anydirection = sumlog(pvalue)$p))

geneStatistics %<>% dplyr::rename(geneSymbol = ID)
#filter custom and unmapped probes
geneStatistics <- geneStatistics %>% filter(!grepl("A_", geneSymbol)) %>% filter(!grepl("CUST_", geneSymbol)) 

#adjust p-values
geneStatistics %<>% mutate(metaP.pos.adj = p.adjust(metaP.pos))
geneStatistics %<>% mutate(metaP.neg.adj = p.adjust(metaP.neg))

(geneStatistics %<>% arrange(metaP.pos))
(geneStatistics %<>% arrange(metaP.neg))

sum(geneStatistics$metaP.pos.adj < 0.05)
sum(geneStatistics$metaP.neg.adj < 0.05)
geneStatistics %<>% dplyr::mutate(isSig = metaP.neg.adj < 0.05 | metaP.pos.adj < 0.05)

geneStatistics %<>% mutate(pValueWithDirection = if_else(metaP.pos < metaP.neg, nrow(geneStatistics) - rank(metaP.pos), -1* nrow(geneStatistics) + rank(metaP.neg)))

plot(geneStatistics$pValueWithDirection, geneStatistics$medianCorrelation)
ggplot(geneStatistics, aes(x=pValueWithDirection, y = medianCorrelation, color = isSig)) + geom_point()


freeSurferData <- read_tsv(paste0(figshare_data_folder, "AllenHBA_DK_ExpressionMatrix.tsv"))
freeSurferData <- tbl_df(melt(freeSurferData))
#add in median donor correlation
geneStatistics %<>% inner_join(filter(freeSurferData, variable == 'Average donor correlation to median') %>% dplyr::select(geneSymbol = X1, DonorCorrelation = value))
cor.test(geneStatistics$DonorCorrelation, abs(geneStatistics$pValueWithDirection), m='s')

#add in average expression level (left hemisphere)
averageExpression <- filter(freeSurferData, grepl("ctx-lh-", variable)) %>% group_by(X1) %>% summarize(averageExpression = mean(value)) %>% dplyr::rename(geneSymbol = X1)
geneStatistics %<>% inner_join(averageExpression)

cor.test(geneStatistics$averageExpression, geneStatistics$metaP.neg, m='s')
cor.test(geneStatistics$averageExpression, geneStatistics$metaP.pos, m='s')
plot(geneStatistics$averageExpression, geneStatistics$medianCorrelation, m='s')
plot(geneStatistics$DonorCorrelation, geneStatistics$pValueWithDirection, m='s')

ggplot(geneStatistics, aes(x=medianCorrelation, y = log(metaP.pos.adj), color = isSig)) + geom_point()

#significant in both directions - should be none
dplyr::filter(geneStatistics, metaP.neg.adj < 0.05 & metaP.pos.adj < 0.05)


paste("Genes with negative correlations:", dplyr::filter(geneStatistics, metaP.neg.adj < 0.05) %>% summarize(n = n()))
paste("Genes with positive correlations:", dplyr::filter(geneStatistics, metaP.pos.adj < 0.05) %>% summarize(n = n()))

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
result %<>% rowwise() %>% mutate(aspect = Ontology(ID))
(result %<>% group_by(U, N1, AUC, P.Value,adj.P.Val) %>% summarize(MainTitle = first(Title),  ID=paste(ID, collapse=","), aspect= first(aspect), allNames = if_else(n() > 1, paste(Title[2:length(Title)], collapse=","), "")))
result %<>% arrange(adj.P.Val)
result$rank <- 1:nrow(result)
result %<>% ungroup() %>% dplyr::select(MainTitle, geneCount = N1, AUC, P.Value, adj.P.Value = adj.P.Val, everything(), -U) %>% arrange(adj.P.Value)

result$adj.P.Value <- signif(result$adj.P.Value, digits=3)
result$AUC <- signif(result$AUC, digits=3)

# Output tables for top ten positively and negatively enriched GO groups for paper
(result.up <- head(filter(result, AUC > 0.5) %>% dplyr::select(-ID), n=10))
(result.down <- head(filter(result, AUC < 0.5) %>% dplyr::select(-ID), n=10))


write_csv( result,  paste(baseFilename,".GO.results.csv",sep=""))
write_csv( dplyr::select(result.up, Name = MainTitle,`Gene Count` = geneCount, AUROC = AUC,  `Adjusted PValue` = adj.P.Value, aspect, synonyms = allNames, rank),  paste(baseFilename,".GO.up10.csv",sep=""))
write_csv( dplyr::select(result.down, Name = MainTitle,`Gene Count` = geneCount, AUROC = AUC, `Adjusted PValue` = adj.P.Value, aspect, synonyms = allNames, rank),  paste(baseFilename,".GO.down10.csv",sep=""))


##### Look at what proportion of the results match certain categories

cat(paste("Total number of GO groups tested",as.character(lengths(result)[1])))

cat(paste("Number of significant GO groups",as.character(lengths(result %>% filter (adj.P.Value < 0.05))[1])))

cat(paste("Number of significant positively enriched GO groups",as.character(lengths(result %>% filter (adj.P.Value < 0.05, AUC > 0.5))[1])))

cat(paste("Number of significant negatively enriched GO groups",as.character(lengths(result %>% filter (adj.P.Value < 0.05, AUC < 0.5))[1])))

cat(paste("Number of significant negatively enriched CC GO groups",as.character(lengths(result %>% filter (adj.P.Value < 0.05, AUC < 0.5,aspect=="CC"))[1])))

cat(paste("Number of significant negatively enriched CC GO groups related to mitochondria",
          as.character(lengths(result %>% filter (adj.P.Value < 0.05, AUC < 0.5,aspect=="CC", grepl("mitoc",MainTitle)))[1])))
cat(paste("Mito groups tested: ", as.character(lengths(result %>% filter (aspect=="CC", grepl("mitoc",MainTitle)))[1])))

cat(paste("Number of significant negatively enriched CC GO groups related to ribosome",
          as.character(lengths(result %>% filter (adj.P.Value < 0.05, AUC < 0.5,aspect=="CC", grepl("ribos",MainTitle)))[1])))
cat(paste("ribos groups tested: ", as.character(lengths(result %>% filter (aspect=="CC", grepl("ribos",MainTitle)))[1])))

cat(paste("Number of significant negatively enriched BP GO groups related to mitochondria",
          as.character(lengths(result %>% filter (adj.P.Value < 0.05, AUC < 0.5,aspect=="BP", grepl("mitoc",MainTitle)))[1])))
cat(paste("ribos groups tested: ", as.character(lengths(result %>% filter (aspect=="BP", grepl("mitoc",MainTitle)))[1])))

cat(paste("Number of significant negatively enriched BP GO groups related to synapse",
          as.character(lengths(result %>% filter (adj.P.Value < 0.05, AUC < 0.5,aspect=="BP", grepl("synap",MainTitle)))[1])))
cat(paste("synapse groups tested: ", as.character(lengths(result %>% filter (aspect=="BP", grepl("synap",MainTitle)))[1])))

source("./R Code/ROCPlots.R")

#plots <- createPlots(sortedGenes, c("GO:0005882", "GO:0032543", "GO:0000502", "GO:0060337", "GO:0060076", "GO:0044309","GO:0007422", "GO:0042552"), geneSetsGO,customNames = as.character(1:8))
#plot(plots$rasterPlot)

#(bothPlots <- plot_grid(plots$AUCPlot, plots$rasterPlot,  nrow = 2, align = "v", rel_heights=c(1,0.9),scale = 0.95)) #add labels = c("A", "B"), for manuscript
#plot(plots$rasterPlot)

filter(result,grepl("myelin",MainTitle))
myelinResult <- filter(result, grepl("myelin|ensheathment",MainTitle),!grepl("peripheral|sphingomyelin",MainTitle))
myelinResult$adj.P.Value <- p.adjust(myelinResult$P.Value, method="BH")
myelinResult$adj.P.Value <- signif(myelinResult$adj.P.Value, digits=3)

myelinResult

write_csv( dplyr::select(myelinResult, Name = MainTitle,`Gene Count` = geneCount, AUROC = AUC,  `Adjusted PValue` = adj.P.Value, aspect, synonyms = allNames, rank),  
           paste0(baseFilename,".GO.myelin.results.csv"))

#plots <- createPlots(sortedGenes, c("GO:0043217", "GO:0043218", "GO:0022010", "GO:0008366","GO:0042552"), geneSetsGO)
#(bothPlots <- plot_grid(plots$AUCPlot, plots$rasterPlot, nrow = 2, align = "v", rel_heights=c(1,0.8),scale = 0.95)) #add labels = c("A", "B"), for manuscript

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
result <- tbl_df(result) %>% dplyr::select(ID, Title, geneCount =N1,AUC,  P.Value, adj.P.Val)

result$adj.P.Val <- signif(result$adj.P.Val, digits=3)
result$AUC <- signif(result$AUC, digits=3)

write_csv( dplyr::select(result, Disease = Title,`Gene Count` = geneCount, AUROC = AUC,  `Adjusted PValue` = adj.P.Val),  
           paste0(baseFilename,".phenocarta.results.csv"))

#plots <- createPlots(sortedGenes, c("DOID_9008", "DOID_3213", "DOID_4233", "DOID_9975"), geneSets)
#(bothPlots <- plot_grid(plots$AUCPlot, plots$rasterPlot, nrow = 2, align = "v", rel_heights=c(1,0.8),scale = 0.95)) #add labels = c("A", "B"), for manuscript

#looking at one result - epilepsy
#filter(geneStatistics, geneSymbol %in% geneSets["DOID_1932"]$GENES$ID)

#################################################################
#################################################################
tmodNames <- data.frame()
modules2genes <- list()


for(geneListFilename in list.files(otherGeneListsFolder, pattern = ".*txt", full.names = T)) {
  print(geneListFilename)
  genesOfInterest <- read.csv(geneListFilename,header=F,stringsAsFactors = F)
  shortName <- gsub(".txt","",gsub(paste0(".*/"),"", geneListFilename))
  
  genesOfInterest$term <- shortName
  
  #already a human gene list
  if (grepl(pattern = "Darmanis.", geneListFilename  ) | grepl(pattern = "PH", geneListFilename) | grepl(pattern = "HouseKeeping", geneListFilename  ) | grepl(pattern = "human", geneListFilename  )) {
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
(result <- tbl_df(result) %>% dplyr::select(Title, geneCount =N1,AUC,  P.Value, adj.P.Val, -ID))

result$adj.P.Val <- signif(result$adj.P.Val, digits=3)
result$AUC <- signif(result$AUC, digits=3)

writeTableWrapper <- function(prefixFilter, result) {
  subsetResult <- filter(result, grepl(prefixFilter, Title))
  subsetResult$Title <- gsub(paste0(prefixFilter, "."), "", subsetResult$Title)
  subsetResult$Title <- gsub("[.]", " ", subsetResult$Title)
  subsetResult$Title <- gsub("[_]", " ", subsetResult$Title)
  subsetResult$Title <- gsub("ligo", "ligodendrocyte", subsetResult$Title)
  subsetResult$adj.P.Val <- p.adjust(subsetResult$P.Value)
  write_csv(dplyr::select(subsetResult, `Cell-type or class` = Title,`Gene Count` = geneCount, AUROC = AUC,  `Adjusted PValue` = adj.P.Val), paste0(baseFilename,".", prefixFilter,".csv")) 
  subsetResult
}

writeTableWrapper("Zeisel", result)
writeTableWrapper("Darmanis", result)
writeTableWrapper("NeuroExpresso.Cortex", result)
writeTableWrapper("Mistry", result)

plots <- createPlots(sortedGenes, c("Zeisel.Oligo", "Zeisel.Neuron.CA1.pryamidal", "Zeisel.Endothelial", "Zeisel.Neuron.interneuron"), geneSets, customNames=NULL)#c("Oligodendrocytes", "CA1 Pyramidal Neurons", "Endothelial Cells", "Interneurons"))

#plots <- createPlots(sortedGenes, , geneSets)=c())
(bothPlots <- plot_grid(plots$AUCPlot, plots$rasterPlot, nrow = 2, align = "v", rel_heights=c(1,0.8),scale = 0.95)) #add labels = c("A", "B"), for manuscript
  

darm <- filter(result, grepl("Darmanis", Title))
plots <- createPlots(sortedGenes, darm$Title, geneSets)
(bothPlots <- plot_grid(plots$AUCPlot, plots$rasterPlot, nrow = 2, align = "v", rel_heights=c(1,0.8))) #add labels = c("A", "B"), for manuscript
plots$rasterPlot
#filter(geneStatistics, geneSymbol %in% geneSets["Darmanis.Oligo"]$GENES$ID)

##################
# Zeng et al. 
################

zengPath <- "./other gene lists/ZengEtAl/1-s2.0-S0092867412003480-mmc2.xlsx"


zengTable <- read.xlsx(zengPath, sheetName = "Final1000New", startRow = 2, stringsAsFactors = F)

(zengTable <- tbl_df(zengTable))

zengTable %<>% select(Gene.symbol, Cortical.marker..human.)
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
zeng_result %<>% arrange(as.character(Title))

filter(geneStatistics, geneSymbol %in% expandedZengTable$Gene.symbol)

write_csv(dplyr::select(zeng_result, `Type` = Title,`Gene Count` = geneCount, AUROC = AUC,  `Adjusted PValue` = adj.P.Val), paste0(baseFilename,".ZengEtAl.csv")) 

zeng_result$Title <- factor(zeng_result$Title, levels=rev(sort(as.character(unique(zeng_result$Title)))))
(barplot <- ggplot(zeng_result, aes(y=AUC-0.5, x=Title, fill=AUC > 0.5)) + geom_col() + 
  coord_flip() + xlab("") + ylab("AUROC") + guides(fill=FALSE) +
  scale_y_continuous(breaks=seq(-0.9, .9, by= .1), labels=seq(-0.9, .9, by= .1) + .5))
barplot <- barplot + geom_text(data = filter(zeng_result, adj.P.Val  < 0.05 & adj.P.Val > 0.005), aes(hjust = -.5*sign(AUC-0.5)+.5), label = "*")
barplot <- barplot + geom_text(data = filter(zeng_result, adj.P.Val < 0.005), aes(hjust = -.5*sign(AUC-0.5)+.5), label = "**")
barplot
ggsave(plot= barplot,paste0(baseFilename,".ZengEtAl.pdf" ), height=5, width=5)
    