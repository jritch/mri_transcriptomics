library(mygene)

#run RunSingleGO.AUROC.Analysis.r first


targetGO <- "MHC protein complex"
targetGO <- "type I interferon signaling pathway"


#write out all positive
for (targetGO in result.up$MainTitle) {
  targetGOID <- dplyr::filter(tbl_df(geneSetsGO$MODULES), Title == targetGO)$ID
  targetGenes <- filter(geneStatistics, geneSymbol %in% unlist(geneSetsGO$MODULES2GENES[targetGOID]))
  nameFrame <- tbl_df(queryMany(targetGenes$geneSymbol, scopes = 'symbol',species="human"))
  (targetGenes <- left_join(targetGenes, nameFrame %>% dplyr::select(geneSymbol = symbol, name)) %>% dplyr::select(geneSymbol, name, everything()))
  targetGenes %<>% distinct()
  write_csv(targetGenes, paste0(baseFilename, ".", targetGO, ".addedStats.csv"))
}

#all neg
for (targetGO in result.down$MainTitle) {
  targetGOID <- dplyr::filter(tbl_df(geneSetsGO$MODULES), Title == targetGO)$ID
  targetGenes <- filter(geneStatistics, geneSymbol %in% unlist(geneSetsGO$MODULES2GENES[targetGOID]))
  nameFrame <- tbl_df(queryMany(targetGenes$geneSymbol, scopes = 'symbol',species="human"))
  (targetGenes <- left_join(targetGenes, nameFrame %>% dplyr::select(geneSymbol = symbol, name)) %>% dplyr::select(geneSymbol, name, everything())) %>% arrange(pValueWithDirection)
  write_csv(targetGenes, paste0(baseFilename, ".", targetGO, ".addedStats.csv"))
}

targetGO <- 'Darmanis.Oligo'
targetGO <- 'Darmanis.Astrocytes'
targetGO <- 'Darmanis.Microglia'


targetGOID <- dplyr::filter(tbl_df(geneSetsCellType$MODULES), Title == targetGO)$ID
targetGenes <- filter(geneStatistics, geneSymbol %in% unlist(geneSetsCellType$MODULES2GENES[targetGOID]))
nameFrame <- tbl_df(queryMany(targetGenes$geneSymbol, scopes = 'symbol',species="human"))
(targetGenes <- left_join(targetGenes, nameFrame %>% dplyr::select(geneSymbol = symbol, name)) %>% dplyr::select(geneSymbol, name, everything())) %>% arrange(pValueWithDirection)

write_csv(targetGenes, paste0(baseFilename, ".", targetGO, ".addedStats.csv"))





targetGO <- 'Human immunodeficiency virus infectious disease'
targetGO <- "Behcet's disease"
targetGO <- "Stillbirth"
targetGO <- "Angelman syndrome"


targetGOID <- dplyr::filter(tbl_df(geneSetsPhenoCarta$MODULES), Title == targetGO)$ID
targetGenes <- filter(geneStatistics, geneSymbol %in% unlist(geneSetsPhenoCarta$MODULES2GENES[targetGOID]))
nameFrame <- tbl_df(queryMany(targetGenes$geneSymbol, scopes = 'symbol',species="human"))
(targetGenes <- left_join(targetGenes, nameFrame %>% dplyr::select(geneSymbol = symbol, name)) %>% dplyr::select(geneSymbol, name, everything())) %>% arrange(pValueWithDirection) %>% distinct()
targetGenes %<>% distinct()
write_csv(targetGenes, paste0(baseFilename, ".", targetGO, ".addedStats.csv"))
