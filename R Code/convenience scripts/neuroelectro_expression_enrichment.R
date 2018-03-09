#this needs emperical p-values, based on shuffling or something as there is a bias for these genes

electro <- read_csv("/Users/lfrench/Downloads/pcbi.1005814.s009.csv")
electroDesc <- read_csv("/Users/lfrench/Downloads/pcbi.1005814.s007.csv") %>% dplyr::select(EphysProp, EphysDescription)
electro <- inner_join(electro, electroDesc)

electro %>% group_by(sign(DiscCorr), EphysProp) %>% summarize(n())
#electro$EphysProp <- sample(electro$EphysProp) #for testing
electro$direction = sign(DiscCorr)

tmodNames <- data.frame()
modules2genes <- list()

for(geneSetName in unique(electro$EphysProp)) {
  for(sign in c(-1,1)) {
    
    if(sign <0 ) fullName <- paste0(geneSetName, ".neg")
    if(sign >0 ) fullName <- paste0(geneSetName, ".pos")
    print(fullName)
    genesOfInterest <- filter(electro, EphysProp == geneSetName, sign(DiscCorr) == sign) %>% .$GeneSymbol
    modules2genes[fullName] <- list(mouse2human(genesOfInterest)$humanGene)
    tmodNames <- rbind(tmodNames, data.frame(ID=fullName, Title = fullName))
  }
}  

geneSets <- makeTmod(modules = tmodNames, modules2genes = modules2genes)
sortedGenesInBackground <- sortedGenes[sortedGenes %in% backGroundGenes]

result <- tmodUtest(sortedGenes, mset=geneSets, qval = 1, filter = T)
result <- tbl_df(result) %>% dplyr::select(Title, geneCount =N1,AUC,  P.Value, adj.P.Val, -ID)
result %<>% rowwise() %>% mutate(P.Value = P.Value * 2) %>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value)) #tmod runs one-sided tests
result

#side analysis
IDgenes <- read_csv("/Users/lfrench/Desktop/results/IDGenes/SecondListFull.txt", col_names = F)$X1
IDgenes <- human2mouse(IDgenes) %>% dplyr::rename(GeneSymbol = mouseGene)
inner_join(IDgenes, electro)
resultTable <- inner_join(IDgenes, electro) %>% dplyr::select(humanGene, GeneSymbol, EphysProp, EphysDescription, everything())
write_csv(resultTable, "/Users/lfrench/Downloads/electro/Intersect.csv")
