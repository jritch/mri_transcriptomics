library(optparse)

#example call: Rscript RunSingleGO.AUROC.Analysis.R -f "/Users/lfrench/Google Drive/gene_list_csvs/T1T2Ratio.cortex.gene_list.csv"
option_list = list(
  make_option(c("-f", "--filename"), type="character", default=NULL, help="filename of per gene statistics", metavar="character")
  #more options here...
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (interactive()) { #set the variables manually if in Rstudio, for testing
  filename <- "/Users/lfrench/Google Drive/gene_list_csvs/T1T2Ratio.cortex.gene_list.csv" #should be passed as an argument so python can call it
} else if (!is.null(opt$filename)) {
  filename <- opt$filename
} else {
  print_help(opt_parser)
  stop()
}

library(readr)
library(dplyr)
library(ggplot2)

library(homologene) #install via install_github('oganm/homologene')
library(org.Hs.eg.db)
library(AnnotationDbi)
library(annotate)
library(GO.db)
library(tmod)

#needs direction statistic
geneStatistics <- read_csv(filename) 
geneStatistics <- geneStatistics %>% filter(X1 != "ID") %>% dplyr::select(geneSymbol = X1, raw_meta_p) 
#filter custom and unmapped probes
geneStatistics <- geneStatistics %>% filter(!grepl("A_", geneSymbol)) %>% filter(!grepl("CUST_", geneSymbol)) 
#sort by p-value
geneStatistics <- arrange(geneStatistics, raw_meta_p)
  
######################################################
# AUC via tmod
######################################################

go_object <- as.list(org.Hs.egGO2ALLEGS)

symbolsInGO <- getSYMBOL(unique(unlist(go_object)), data='org.Hs.eg')

#todo - check if this is in the right order
sortedGenes <- geneStatistics$geneSymbol

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
geneSets <- makeTmod(modules = tmodNames, modules2genes = modules2genes)


result <- tmodUtest(c(sortedGenes), mset=geneSets, qval = 1)
head(result, n=20)
result <- tmodUtest(c(sortedGenes), mset=geneSets, qval = 1, filter = T)
head(result, n=20)
head(subset(result, AUC > 0.5), n=20)
head(subset(result, AUC < 0.5), n=20)

write.csv(result, file=paste0(filename, ".enrichment.GO.csv"))

#make AUC plots to visualize - needs a name of the gene list
#evidencePlot(sortedGenes, mset=geneSets, m="GO:0000910")

#TODO - filter out duplicate sets?
#################################################################
#################################################################

loadPhenocarta <- function(taxon, geneBackground) {
  #guess column types from the whole dataset, basically
  phenocarta <- read_tsv("/Users/lfrench/Desktop/results/mri_transcriptomics/other gene lists/phenoCarta/AllPhenocartaAnnotations.downloadedOct28.2016.tsv", skip = 4, guess_max = 130000)
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

result <- tmodUtest(c(sortedGenes), mset=geneSets, qval = 1)
head(result, n=20)
result <- tmodUtest(c(sortedGenes), mset=geneSets, qval = 1, filter = T)
head(result, n=20)

write.csv(result, file=paste0(filename, ".enrichment.PhenoCarta.csv"))

#looking at one result - epilepsy
#filter(geneStatistics, geneSymbol %in% geneSets["DOID_1932"]$GENES$ID)


#################################################################
#################################################################
tmodNames <- data.frame()
modules2genes <- list()

#need to set folder here
for(geneListFilename in list.files("/Users/lfrench/Desktop/results/mri_transcriptomics/other gene lists/", pattern = ".*txt", full.names = T)) {
  print(geneListFilename)
  genesOfInterest <- read.csv(geneListFilename,header=F,stringsAsFactors = F)
  shortName <- gsub(".txt","",gsub(paste0(".*/"),"", geneListFilename))
  genesOfInterest$term <- shortName
  
  #already a human gene list
  if (grepl(pattern = "Darmanis.", geneListFilename  ) | grepl(pattern = "House keeping", geneListFilename  )) {
    modules2genes[shortName] <- list(genesOfInterest$V1)
  } else { #needs conversion from mouse
    modules2genes[shortName] <- list(mouse2human(genesOfInterest$V1)$humanGene)
  }

  tmodNames <- rbind(tmodNames, data.frame(ID=shortName, Title = shortName))
}
geneSets <- makeTmod(modules = tmodNames, modules2genes = modules2genes)

result <- tmodUtest(sortedGenes, mset=geneSets, qval = 1, filter = T)
result
subset(result, AUC > 0.5)
subset(result, AUC < 0.5)

evidencePlot(sortedGenes, mset=geneSets, m="Zeisel.Microglia")

write.csv(result, file=paste0(filename, ".enrichment.CellTypes.csv"))

