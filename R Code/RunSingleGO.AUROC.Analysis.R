library(optparse)
#todo - Leon is using GOSOURCEDATE: 20160305 (type GO.db to find out)

#example call: Rscript RunSingleGO.AUROC.Analysis.R -f "/Users/lfrench/Google Drive/gene_list_csvs/T1T2Ratio.cortex.gene_list.csv"
option_list = list(
  make_option(c("-f", "--filename"), type="character", default=NULL, help="filename of per gene statistics", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, help="output directory", metavar="character")
  
  #more options here...
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Sys.info()["nodename"]

if (interactive()) { #set the variables manually if in Rstudio, for testing
  
  filename <- "/Users/lfrench/Google Drive/gene_list_csvs/T1T2Ratio.cortex.gene_list.csv" #should be passed as an argument so python can call it
  filename <- "/Users/lfrench/Google Drive/gene_list_csvs/T1.cortex.gene_list.csv" #should be passed as an argument so python can call it
  filename <- "/Users/lfrench/Google Drive/gene_list_csvs/T2.cortex.gene_list.csv" #should be passed as an argument so python can call it #looks like the ratio
  
  #filename <- "/Users/jritchie/data/garbage2/T1T2Ratio.cortex.gene_list.csv"
  #filename <- "/Users/jritchie/Google Drive/4th year/Thesis/gene_list_csvs_float/T1.cortex.gene_list.csv"
  #filename <- "/Users/jritchie/Google Drive/4th year/Thesis/gene_list_csvs_float/T2.cortex.gene_list.csv"
  #filename <- "/Users/jritchie/Google Drive/4th year/Thesis/gene_list_csvs_float/T1T2Ratio.cortex.gene_list.csv"
  #filename <- "/Users/jritchie/Google Drive/4th year/Thesis/gene_list_csvs_stripped/T1T2Ratio.cortex.gene_list.csv"
  #filename <- "/Users/jritchie/Google Drive/4th year/Thesis/gene_list_csvs/T1T2Ratio.cortex.gene_list.csv"
  #filename <- "/Users/jritchie/Google Drive/4th year/Thesis/gene_list_csvs/T2.full_brain.gene_list.csv"
  #filename <- "/Users/jritchie/Google Drive/4th year/Thesis/gene_list_csvs/T1.full_brain.gene_list.csv"
  #filename <- "/Users/jritchie/data/garbage4/T1T2Ratio.cortex_excluding_limbic_lobe.gene_list.csv"
  #filename <- "/Users/jritchie/data/garbage3/T1T2Ratio.full_brain.gene_list.csv"
  #filename <- "/Users/jritchie/data/garbage3/T1T2Ratio.cortex_excluding_limbic_lobe.gene_list.csv"
  #filename <- "/Users/jritchie/data/garbage3/T1T2Ratio.cortex.gene_list.csv"
  
  filename <- "/Users/jritchie/data/final_with_bias_correction/T1T2Ratio.cortex.gene_list.csv"
  filename <- "/Users/jritchie/data/final_with_bias_correction/T1T2Ratio.cortex_excluding_limbic_lobe.gene_list.csv"
  filename <- "/Users/jritchie/data/final_with_bias_correction/T1T2Ratio.full_brain.gene_list.csv"
  filename <- "/Users/jritchie/data/final_with_bias_correction/T1.hippocampus.gene_list.csv"  
  filename <- "/Users/jritchie/data/final_with_bias_correction/T2.hippocampus.gene_list.csv"
  filename <- "/Users/jritchie/data/final_with_bias_correction/T1T2Ratio.hippocampus.gene_list.csv"
  
  filename <- "/Users/jritchie/data/final/T1T2Ratio.cortex.gene_list.csv"
  filename <- "/Users/jritchie/data/final/T1T2Ratio.cortex_excluding_limbic_lobe.gene_list.csv"
  filename <- "/Users/jritchie/data/final/T1T2Ratio.full_brain.gene_list.csv"
  
  output_dir <- "/Users/jritchie/data/bias_corrected_results"
  output_dir <- "/Users/jritchie/data/results"
  
  
  if(Sys.info()['nodename'] == "RES-C02RF0T2.local") {
    filename <- "/Users/lfrench/Desktop/results/mri_transcriptomics/results/T1T2Ratio.cortex.gene_list.csv"
  } else {
    #filename <- "C://Users/Jacob/Google Drive/4th Year/Thesis/gene_list_csvs/T1T2Ratio.cortex.gene_list.csv"
    #filename <- "C://Users/Jacob/Google Drive/4th Year/Thesis/gene_list_csvs/T2.full_brain.gene_list.csv"
    #filename <- "C://Users/Jacob/Google Drive/4th Year/Thesis/gene_list_csvs/T1T2Ratio.full_brain.gene_list.csv"
  }
  
} else if (!is.null(opt$filename)) {
  filename <- opt$filename
  output_dir <- opt$output_dir
  
} else {
  print_help(opt_parser)
  stop()
}

cat(paste("USING INPUT FILE:",filename))
cat(paste("WRITING OUTPUT FILES TO:",output_dir))

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

source("./R Code/Utils.R")

#needs direction statistic
geneStatistics <- read_csv(filename) 
geneStatistics <- geneStatistics %>% filter(X1 != "ID")

geneStatistics <- rowwise(geneStatistics) %>% mutate(medianCorrelation = median(c(correlation,X3,X4,X5,X6,X7)))

#geneStatistics$adjustAgain <- p.adjust(geneStatistics$raw_meta_p, method="fdr")
geneStatistics$adjustAgain <- p.adjust(geneStatistics$raw_meta_p, method="BH")

geneStatistics <- rowwise(geneStatistics) %>% mutate(adj_p = sumlog(c(adjusted,X15,X16,X17,X18,X19))$p)

#cor.test(geneStatistics$adjustAgain, geneStatistics$adjusted_meta_p) #todo - fix, the adjusted p-values are not lining up
#cor.test(geneStatistics$raw_meta_p, geneStatistics$adjusted_meta_p,m='s')
#cor.test(geneStatistics$raw_meta_p, geneStatistics$adjustAgain,m='s')
#plot(geneStatistics$adjustAgain, geneStatistics$adjusted_meta_p)
#plot(geneStatistics$raw_meta_p, geneStatistics$adjusted_meta_p)
#plot(geneStatistics$raw_meta_p, geneStatistics$adjustAgain)

geneStatistics <- geneStatistics %>% dplyr::select(geneSymbol = X1, raw_meta_p, adjusted_meta_p, medianCorrelation) 
#filter custom and unmapped probes
geneStatistics <- geneStatistics %>% filter(!grepl("A_", geneSymbol)) %>% filter(!grepl("CUST_", geneSymbol)) 

print(paste(sum(geneStatistics$adjusted_meta_p < 0.05),"of",length(row.names(geneStatistics)),"genes survive BH multiple-test correction"))
# Re-adjust after getting rid of other probes to see if it affects analysis.
survive <- sum(p.adjust(geneStatistics$raw_meta_p, method="BH") < 0.05)
print(paste(survive,"of",length(row.names(geneStatistics)),"genes survive BH multiple-test correction after getting rid of custom probes"))
# Doesn't have much affect

significantCorrelations <- geneStatistics %>% filter(adjusted_meta_p < 0.05) %>% dplyr::select(medianCorrelation)
nonSignificantCorrelations <- geneStatistics %>% filter(adjusted_meta_p > 0.05) %>% dplyr::select(medianCorrelation)

#plot_grid(qplot(significantCorrelations, geom="histogram",bins=80),
#          qplot(nonSignificantCorrelations, geom="histogram",bins=80),
#          qplot(geneStatistics$medianCorrelation, geom="histogram",bins=80),
#          qplot(geneStatistics$adjusted_meta_p, geom="histogram",bins=100) + geom_vline(xintercept=0.05, color="red"))

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

result$adj.P.Val <- signif(result$adj.P.Val, digits=2)
result$AUC <- signif(result$AUC, digits=3)

# Output tables for top ten positively and negatively enriched GO groups

result1 <- head(filter(result, AUC > 0.5) %>% dplyr::select(-ID), n=10)
result2 <- head(filter(result, AUC < 0.5) %>% dplyr::select(-ID), n=10)

split_filename = strsplit(filename,"/")
analysis_id = gsub(".gene_list.csv","",unlist(split_filename)[length(unlist(split_filename))])

output_filename <- paste(analysis_id,".GO.enrichment.positive.results.csv",sep="")
writeTableWrapper(result1,output_dir,output_filename,col_indices = c(1:3,5,7),col_names=c("GO group",    "gene count",    "AUROC",    "adjusted p. value", "aspect"))

output_filename <- paste(analysis_id,".GO.enrichment.negative.results.csv",sep="")
writeTableWrapper(result2,output_dir,output_filename,col_indices = c(1:3,5,7),col_names=c("GO group",    "gene count",    "AUROC",    "adjusted p. value", "aspect"))

# Examine the 3 aspects in the Gene Ontology separately

head(filter(result, AUC < 0.5, aspect=="BP", adj.P.Val < 0.05) %>% dplyr::select(-ID), n=20)
head(filter(result, AUC < 0.5, aspect=="CC", adj.P.Val < 0.05) %>% dplyr::select(-ID), n=20)
head(filter(result, AUC < 0.5, aspect=="MF", adj.P.Val < 0.05) %>% dplyr::select(-ID), n=20)


##### Look at what proportion of the results match certain categories

result %>% filter (AUC < 0.5, adj.P.Val < 0.05)

cat(paste("Total number of GO groups",as.character(lengths(result)[1])))

cat(paste("Number of significant GO groups",as.character(lengths(result %>% filter (adj.P.Val < 0.05))[1])))

cat(paste("Number of significant positively enriched GO groups",as.character(lengths(result %>% filter (adj.P.Val < 0.05, AUC > 0.5))[1])))

cat(paste("Number of significant negatively enriched GO groups",as.character(lengths(result %>% filter (adj.P.Val < 0.05, AUC < 0.5))[1])))

cat(paste("Number of significant negatively enriched CC GO groups",as.character(lengths(result %>% filter (adj.P.Val < 0.05, AUC < 0.5,aspect=="CC"))[1])))

cat(paste("Number of significant negatively enriched CC GO groups related to mitochondria",
          as.character(lengths(result %>% filter (adj.P.Val < 0.05, AUC < 0.5,aspect=="CC", grepl("mitoc",MainTitle)))[1])))

cat(paste("Number of significant negatively enriched CC GO groups related to ribosome",
          as.character(lengths(result %>% filter (adj.P.Val < 0.05, AUC < 0.5,aspect=="CC", grepl("ribos",MainTitle)))[1])))

cat(paste("Number of significant negatively enriched BP GO groups related to mitochondria",
          as.character(lengths(result %>% filter (adj.P.Val < 0.05, AUC < 0.5,aspect=="BP", grepl("mitoc",MainTitle)))[1])))

cat(paste("Number of significant negatively enriched BP GO groups related to synapse",
          as.character(lengths(result %>% filter (adj.P.Val < 0.05, AUC < 0.5,aspect=="BP", grepl("synap",MainTitle)))[1])))

filter(result,grepl("myelin",MainTitle))

source("./R Code/ROCPlots.R")

#plots <- createPlots(sortedGenes, c("GO:0005882", "GO:0032543", "GO:0000502", "GO:0060337", "GO:0060076", "GO:0044309","GO:0007422", "GO:0042552"), geneSetsGO,customNames = as.character(1:8))
#plot(plots$rasterPlot)

#(bothPlots <- plot_grid(plots$AUCPlot, plots$rasterPlot,  nrow = 2, align = "v", rel_heights=c(1,0.9),scale = 0.95)) #add labels = c("A", "B"), for manuscript
#plot(plots$rasterPlot)

myelinResult <- filter(result, grepl("myelin|ensheathment",MainTitle),!grepl("peripheral|sphingomyelin",MainTitle))
#remove synonyms/duplicate results - use the first name of a match
myelinResult %<>% group_by(U, N1, AUC, P.Value) %>% dplyr::select(rank, MainTitle, ID, everything()) %>% arrange(P.Value)
myelinResult$adj.P.Val <- p.adjust(myelinResult$P.Value, method="BH")

myelinResult$adj.P.Val <- signif(myelinResult$adj.P.Val, digits=2)
myelinResult$AUC <- signif(myelinResult$AUC, digits=3)

myelinResult

output_filename <- paste(analysis_id,".myelin.enrichment.results.csv",sep="")
writeTableWrapper(myelinResult,output_dir,output_filename,col_indices = c(2,1,4,5,7,9),col_names=c("GO group",  "rank",  "gene count",    "AUROC",    "adjusted p. value", "aspect"))

#plots <- createPlots(sortedGenes, c("GO:0043217", "GO:0043218", "GO:0022010", "GO:0008366","GO:0042552"), geneSetsGO)
#(bothPlots <- plot_grid(plots$AUCPlot, plots$rasterPlot, nrow = 2, align = "v", rel_heights=c(1,0.8),scale = 0.95)) #add labels = c("A", "B"), for manuscript

#write.csv(result, file=paste0(filename, ".enrichment.GO.csv"))

#################################################################
#################################################################

loadPhenocarta <- function(taxon, geneBackground) {
  #guess column types from the whole dataset, basically
  #phenocarta <- read_tsv("/Users/lfrench/Desktop/results/mri_transcriptomics/other gene lists/phenoCarta/AllPhenocartaAnnotations.downloadedOct28.2016.tsv", skip = 4, guess_max = 130000)
  phenocarta <- read_tsv("/Users/jritchie/Google Drive/4th Year/Thesis/other gene lists/phenoCarta/AllPhenocartaAnnotations.downloadedOct28.2016.tsv", skip = 4, guess_max = 130000)
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

result$adj.P.Val <- signif(result$adj.P.Val, digits=2)
result$AUC <- signif(result$AUC, digits=3)

result <- head(result, n=10)

output_filename <- paste(analysis_id,".phenocarta.enrichment.results.csv",sep="")
writeTableWrapper(result,output_dir,output_filename,col_indices = c(2:4,6),col_names=c("Disease group",    "gene count",    "AUROC",    "adjusted p. value"))

#write.csv(result, file=paste0(filename, ".enrichment.PhenoCarta.csv"))

#plots <- createPlots(sortedGenes, c("DOID_9008", "DOID_3213", "DOID_4233", "DOID_9975"), geneSets)
#(bothPlots <- plot_grid(plots$AUCPlot, plots$rasterPlot, nrow = 2, align = "v", rel_heights=c(1,0.8),scale = 0.95)) #add labels = c("A", "B"), for manuscript

#looking at one result - epilepsy
#filter(geneStatistics, geneSymbol %in% geneSets["DOID_1932"]$GENES$ID)

#################################################################
#################################################################
tmodNames <- data.frame()
modules2genes <- list()

#need to set folder here
if(Sys.info()['user'] == "lfrench") {
  otherGeneListsFolder <- "/Users/lfrench/Desktop/results/mri_transcriptomics/other gene lists/"
} else {
  otherGeneListsFolder <- "/Users/jritchie/git-repos/mri_transcriptomics/other gene lists/"
}

for(geneListFilename in list.files(otherGeneListsFolder, pattern = ".*txt", full.names = T)) {
  print(geneListFilename)
  genesOfInterest <- read.csv(geneListFilename,header=F,stringsAsFactors = F)
  shortName <- gsub(".txt","",gsub(paste0(".*/"),"", geneListFilename))
  
  genesOfInterest$term <- shortName
  
  #already a human gene list
  if (grepl(pattern = "Darmanis.", geneListFilename  ) | grepl(pattern = "SCZ", geneListFilename) | grepl(pattern = "HouseKeeping", geneListFilename  ) | grepl(pattern = "human", geneListFilename  )) {
    modules2genes[shortName] <- list(genesOfInterest$V1)
  } else { #needs conversion from mouse
    print(" converting from mouse to human")
    modules2genes[shortName] <- list(mouse2human(genesOfInterest$V1)$humanGene)
  }

  tmodNames <- rbind(tmodNames, data.frame(ID=shortName, Title = shortName))
}
geneSets <- makeTmod(modules = tmodNames, modules2genes = modules2genes)

result <- tmodUtest(sortedGenes, mset=geneSets, qval = 1, filter = F)
result <- tbl_df(result) %>% dplyr::select(Title, geneCount =N1,AUC,  P.Value, adj.P.Val, ID)

result$adj.P.Val <- signif(result$adj.P.Val, digits=2)
result$AUC <- signif(result$AUC, digits=3)


(result1 <- subset(result, AUC > 0.5))
#(result1 <- subset(result1, adj.P.Val < 0.05))

output_filename <- paste(analysis_id,".cell.type.positive.enrichment.results.csv",sep="")
writeTableWrapper(result1,output_dir,output_filename,col_indices = c(1:3,5),col_names=c("Cell type",    "gene count",    "AUROC",    "adjusted p. value"))

(result2 <- subset(result, AUC < 0.5))

output_filename <- paste(analysis_id,".cell.type.negative.enrichment.results.csv",sep="")
writeTableWrapper(result2,output_dir,output_filename,col_indices = c(1:3,5),col_names=c("Cell type",    "gene count",    "AUROC",    "adjusted p. value"))

darm <- filter(result, grepl("Darmanis", Title))
darm$adj.P.Val <- p.adjust(darm$P.Value)
#darm %>% filter(AUC < 0.5)
darm %>% filter(AUC > 0.5)

zeisel <- filter(result, grepl("Zeisel", Title))
zeisel$adj.P.Val <- p.adjust(zeisel$P.Value)

#zeisel %>% filter(AUC < 0.5)
zeisel %>% filter(AUC > 0.5)

neuro_expresso <- filter(result, grepl("Expresso", Title))
neuro_expresso$adj.P.Val <- p.adjust(neuro_expresso$P.Value)
#neuro_expresso %>% filter(AUC < 0.5)
neuro_expresso %>% filter(AUC > 0.5)

#print out the genes for the oligo list
#colnames(geneSets)

filter(geneStatistics, geneSymbol %in% modules2genes$Darmanis.Oligo)

plots <- createPlots(sortedGenes, c("Zeisel.oligo", "Zeisel.Neuron.CA1.pryamidal", "Zeisel.Endothelial", "Zeisel.Neuron.interneuron"), geneSets, customNames=NULL)#c("Oligodendrocytes", "CA1 Pyramidal Neurons", "Endothelial Cells", "Interneurons"))

#plots <- createPlots(sortedGenes, , geneSets)=c())
(bothPlots <- plot_grid(plots$AUCPlot, plots$rasterPlot, nrow = 2, align = "v", rel_heights=c(1,0.8),scale = 0.95)) #add labels = c("A", "B"), for manuscript
  
plots <- createPlots(sortedGenes, darm$ID, geneSets)
(bothPlots <- plot_grid(plots$AUCPlot, plots$rasterPlot, nrow = 2, align = "v", rel_heights=c(1,0.8))) #add labels = c("A", "B"), for manuscript

#filter(geneStatistics, geneSymbol %in% geneSets["Darmanis.Oligo"]$GENES$ID)

#write.csv(result, file=paste0(filename, ".enrichment.CellTypes.csv"))

##################
# Zeng et al. 
################

library(xlsx)
library(dplyr)
library(magrittr)
library(tidyr)

if(Sys.info()['user'] == "lfrench") {
  zengPath <- "/Users/lfrench/Desktop/results/mri_transcriptomics/other gene lists/ZengEtAl/1-s2.0-S0092867412003480-mmc2.xlsx"
} else {
  zengPath <- "/Users/jritchie/git-repos/mri_transcriptomics/other gene lists/ZengEtAl/1-s2.0-S0092867412003480-mmc2.xlsx"
}


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

zeng_result <- tmodUtest(sortedGenesInBackground, mset=geneSets, qval = 1, filter = F)
zeng_result <- tbl_df(result) %>% dplyr::select(Title, geneCount =N1,AUC,  P.Value, adj.P.Val, ID)
