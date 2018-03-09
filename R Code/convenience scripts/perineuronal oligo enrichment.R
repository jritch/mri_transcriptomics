library(homologene)
library(dplyr)
library(magrittr)
library(readr)
#####################################
#####################################
#perineuronal oligo marker gene list
#####################################
#####################################

tmodNames <- data.frame()
modules2genes <- list()

geneListFilename <- "./other gene lists/RatGeneLists/RatBrainPerineuronalOligoGenes.SzuchetEtAl.txt"
genesOfInterest <- read.csv(geneListFilename,header=F,stringsAsFactors = F)
shortName <- gsub(".txt","",gsub(paste0(".*/"),"", geneListFilename))

genesOfInterest$term <- shortName
print("Converting from rat to human")
modules2genes[shortName] <- list(homologene(genesOfInterest$V1, 10116, 9606)$`9606`)
tmodNames <- rbind(tmodNames, data.frame(ID=shortName, Title = shortName))

#read supplement tables
SzuchetTableS1 <- read_tsv("./other gene lists/RatGeneLists/EJN_7922_sm_DataS1.txt")
SzuchetTableS2 <- read_tsv("./other gene lists/RatGeneLists/EJN_7922_sm_DataS2.txt")
nrow(SzuchetTableS1 %<>% filter(Fold > 0))
nrow(SzuchetTableS2 %<>% filter(Fold > 0))

myGeneResult <- tbl_df(mygene::queryMany(SzuchetTableS1$`GenBank gi`,scopes="reporter"))
myGeneResult %<>% filter(taxid == 10116) %>% select(symbol, `GenBank gi` = query)
SzuchetTableS1 <- inner_join(SzuchetTableS1, myGeneResult)

myGeneResult <- tbl_df(mygene::queryMany(SzuchetTableS2$`GenBank gi`,scopes="reporter"))
myGeneResult %<>% filter(taxid == 10116) %>% select(symbol, `GenBank gi` = query)
SzuchetTableS2 <- inner_join(SzuchetTableS2, myGeneResult)

shortName <- "SzuchetD1andD2"
periGenes <- intersect(SzuchetTableS2$symbol, SzuchetTableS1$symbol)
modules2genes[shortName] <- list(homologene(periGenes, 10116, 9606)$`9606`)
tmodNames <- rbind(tmodNames, data.frame(ID=shortName, Title = shortName))

geneSets <- makeTmod(modules = tmodNames, modules2genes = modules2genes)

result <- tmodUtest(sortedGenes, mset=geneSets, qval = 1.1, filter = F)
result <- tbl_df(result) %>% dplyr::select(Title, geneCount =N1,AUC,  P.Value, adj.P.Val, -ID)
result %<>% rowwise() %>% mutate(P.Value = P.Value * 2) %>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value, method="holm")) #tmod runs one-sided tests
result

geneStatistics %>% filter(geneSymbol %in% geneSets$MODULES2GENES$RatBrainPerineuronalOligoGenes.SzuchetEtAl)
geneStatistics %>% filter(geneSymbol %in% geneSets$MODULES2GENES$SzuchetD2)
