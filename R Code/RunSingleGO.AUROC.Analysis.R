library(optparse)
#todo - Leon is using GOSOURCEDATE: 20160305 (type GO.db to find out)

#example call: Rscript RunSingleGO.AUROC.Analysis.R -f "/Users/lfrench/Google Drive/gene_list_csvs/T1T2Ratio.cortex.gene_list.csv"
option_list = list(
  make_option(c("-f", "--filename"), type="character", default=NULL, help="filename of per gene statistics", metavar="character")
  #more options here...
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Sys.info()["nodename"]

if (interactive()) { #set the variables manually if in Rstudio, for testing
  filename <- "/Users/lfrench/Google Drive/gene_list_csvs/T1T2Ratio.cortex.gene_list.csv" #should be passed as an argument so python can call it
  filename <- "/Users/lfrench/Google Drive/gene_list_csvs/T1.cortex.gene_list.csv" #should be passed as an argument so python can call it
  filename <- "/Users/lfrench/Google Drive/gene_list_csvs/T2.cortex.gene_list.csv" #should be passed as an argument so python can call it #looks like the ratio
  
  if(Sys.info()['nodename'] == "RES-C02RF0T2.local") {
    filename <- "/Users/lfrench/Desktop/results/mri_transcriptomics/results/T1T2Ratio.cortex.gene_list.csv"
  } else {
    filename <- "C://Users/Jacob/Google Drive/4th Year/Thesis/gene_list_csvs/T1T2Ratio.cortex.gene_list.csv"
    filename <- "C://Users/Jacob/Google Drive/4th Year/Thesis/gene_list_csvs/T2.full_brain.gene_list.csv"
    filename <- "C://Users/Jacob/Google Drive/4th Year/Thesis/gene_list_csvs/T1T2Ratio.full_brain.gene_list.csv"
  }
  
} else if (!is.null(opt$filename)) {
  filename <- opt$filename
} else {
  print_help(opt_parser)
  stop()
}

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


#needs direction statistic
geneStatistics <- read_csv(filename) 
geneStatistics <- geneStatistics %>% filter(X1 != "ID")

geneStatistics <- rowwise(geneStatistics) %>% mutate(medianCorrelation = median(c(correlation,X3,X4,X5,X6,X7)))

geneStatistics$adjustAgain <- p.adjust(geneStatistics$raw_meta_p, method="fdr")
geneStatistics$adjustAgain <- p.adjust(geneStatistics$raw_meta_p, method="BH")

geneStatistics <- rowwise(geneStatistics) %>% mutate(adj_p = sumlog(c(adjusted,X15,X16,X17,X18,X19))$p)


cor.test(geneStatistics$adjustAgain, geneStatistics$adjusted_meta_p) #todo - fix, the adjusted p-values are not lining up
plot(geneStatistics$adjustAgain, geneStatistics$adjusted_meta_p)

#comparison
plot(geneStatistics$raw_meta_p, geneStatistics$adjusted_meta_p)

plot(geneStatistics$raw_meta_p, geneStatistics$adjustAgain)

plot(geneStatistics$raw_meta_p, geneStatistics$adj_p)

cor.test(geneStatistics$raw_meta_p, geneStatistics$adjusted_meta_p,m='s')
cor.test(geneStatistics$raw_meta_p, geneStatistics$adjustAgain,m='s')


geneStatistics <- geneStatistics %>% dplyr::select(geneSymbol = X1, raw_meta_p, medianCorrelation) 
#filter custom and unmapped probes
geneStatistics <- geneStatistics %>% filter(!grepl("A_", geneSymbol)) %>% filter(!grepl("CUST_", geneSymbol)) 

#sort by median correlation
geneStatistics <- arrange(geneStatistics, desc(medianCorrelation))
#geneStatistics <- arrange(geneStatistics, raw_meta_p)

#todo - check if this is in the right order
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
    if (!(length(genesymbols) > 10 & length(genesymbols) < 200)) next();
    
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
result %<>% dplyr::select(MainTitle, N1, AUC, P.Value, adj.P.Val, everything()) %>% arrange(adj.P.Val)

result1 <- head(filter(result, AUC > 0.5) %>% dplyr::select(-ID), n=20)
resutl2 <- head(filter(result, AUC < 0.5) %>% dplyr::select(-ID), n=20)

head(filter(result, AUC < 0.5, aspect=="BP", adj.P.Val < 0.05) %>% dplyr::select(-ID), n=20)
head(filter(result, AUC < 0.5, aspect=="CC", adj.P.Val < 0.05) %>% dplyr::select(-ID), n=20)
head(filter(result, AUC < 0.5, aspect=="MF", adj.P.Val < 0.05) %>% dplyr::select(-ID), n=20)

source("./R Code/ROCPlots.R")
plots <- createPlots(sortedGenes, c("GO:0005882", "GO:0032543", "GO:0000502", "GO:0060337", "GO:0060076", "GO:0044309","GO:0007422", "GO:0042552"), geneSetsGO)
(bothPlots <- plot_grid(plots$AUCPlot, plots$rasterPlot,  nrow = 2, align = "v", rel_heights=c(1,0.9))) #add labels = c("A", "B"), for manuscript
plot(plots$rasterPlot)

myelinResult <- filter(result, grepl("myelin|ensheathment",MainTitle),!grepl("peripheral|sphingomyelin",MainTitle))
#remove synonyms/duplicate results - use the first name of a match
myelinResult %<>% group_by(U, N1, AUC, P.Value) %>% dplyr::select(rank, MainTitle, ID, everything()) %>% arrange(P.Value)
myelinResult$adj.P.Val <- p.adjust(myelinResult$P.Value, method="fdr")
myelinResult

plots <- createPlots(sortedGenes, c("GO:0042552", "GO:0043218", "GO:0022010", "GO:0043209"), geneSetsGO)
(bothPlots <- plot_grid(plots$AUCPlot, plots$rasterPlot, nrow = 2, align = "v", rel_heights=c(1,0.8))) #add labels = c("A", "B"), for manuscript


#write.csv(result, file=paste0(filename, ".enrichment.GO.csv"))

#################################################################
#################################################################

loadPhenocarta <- function(taxon, geneBackground) {
  #guess column types from the whole dataset, basically
  phenocarta <- read_tsv("/Users/lfrench/Desktop/results/mri_transcriptomics/other gene lists/phenoCarta/AllPhenocartaAnnotations.downloadedOct28.2016.tsv", skip = 4, guess_max = 130000)
  #phenocarta <- read_tsv("C://Users/Jacob/Google Drive/4th Year/Thesis/other gene lists/phenoCarta/AllPhenocartaAnnotations.downloadedOct28.2016.tsv", skip = 4, guess_max = 130000)
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

geneSets <- loadPhenocarta("human", sortedGenes)

result <- tmodUtest(c(sortedGenes), mset=geneSets, qval = 1, filter = F)
result <- tbl_df(result) %>% dplyr::select(ID, Title, geneCount =N1,AUC,  P.Value, adj.P.Val)
head(result, n=20)

#write.csv(result, file=paste0(filename, ".enrichment.PhenoCarta.csv"))

plots <- createPlots(sortedGenes, c("DOID_9008", "DOID_3213", "DOID_4233", "DOID_9975"), geneSets)
(bothPlots <- plot_grid(plots$AUCPlot, plots$rasterPlot, nrow = 2, align = "v", rel_heights=c(1,0.8))) #add labels = c("A", "B"), for manuscript

#looking at one result - epilepsy
#filter(geneStatistics, geneSymbol %in% geneSets["DOID_1932"]$GENES$ID)


#################################################################
#################################################################
tmodNames <- data.frame()
modules2genes <- list()

#need to set folder here
if(Sys.info()['sysname'] == "Darwin") {
  otherGeneListsFolder <- "/Users/lfrench/Desktop/results/mri_transcriptomics/other gene lists/"
} else {
  otherGeneListsFolder <- "C://Users/Jacob/Google Drive/4th Year/Thesis/other gene lists/"
}

for(geneListFilename in list.files(otherGeneListsFolder, pattern = ".*txt", full.names = T)) {
  print(geneListFilename)
  
  for (keyword in c("Cajigas","Axon")) {
    if (grepl(keyword,geneListFilename)) {
      print("Skipping")
      next()
    }
  }

  genesOfInterest <- read.csv(geneListFilename,header=F,stringsAsFactors = F)
  shortName <- gsub(".txt","",gsub(paste0(".*/"),"", geneListFilename))
  genesOfInterest$term <- shortName
  
  #already a human gene list
  if (grepl(pattern = "Darmanis.", geneListFilename  ) | grepl(pattern = "HouseKeeping", geneListFilename  ) | grepl(pattern = "human", geneListFilename  )) {
    modules2genes[shortName] <- list(genesOfInterest$V1)
  } else { #needs conversion from mouse
    print(" converting to mouse")
    modules2genes[shortName] <- list(mouse2human(genesOfInterest$V1)$humanGene)
  }

  tmodNames <- rbind(tmodNames, data.frame(ID=shortName, Title = shortName))
}
geneSets <- makeTmod(modules = tmodNames, modules2genes = modules2genes)


result <- tmodUtest(sortedGenes, mset=geneSets, qval = 1, filter = F)
result <- tbl_df(result) %>% dplyr::select(Title, geneCount =N1,AUC,  P.Value, adj.P.Val, ID)

(result1 <- subset(result, AUC > 0.5))
(result2 <- subset(result, AUC < 0.5))
darm <- filter(result, grepl("Darmanis", Title))
darm$adj.P.Val <- p.adjust(darm$P.Value)

zeisel <- filter(result, grepl("Zeisel", Title))
zeisel$adj.P.Val <- p.adjust(zeisel$P.Value)
zeisel


#print out the genes for the oligo list
filter(geneStatistics, geneSymbol %in% modules2genes$Darmanis.Oligo)

plots <- createPlots(sortedGenes, c("Zeisel.oligo", "Zeisel.Neuron.CA1.pryamidal", "Zeisel.Endothelial", "Zeisel.Neuron.interneuron"), geneSets)
(bothPlots <- plot_grid(plots$AUCPlot, plots$rasterPlot, nrow = 2, align = "v", rel_heights=c(1,0.8))) #add labels = c("A", "B"), for manuscript

plots <- createPlots(sortedGenes, darm$ID, geneSets)
plots$rasterPlot

#filter(geneStatistics, geneSymbol %in% geneSets["Darmanis.Oligo"]$GENES$ID)

#write.csv(result, file=paste0(filename, ".enrichment.CellTypes.csv"))




