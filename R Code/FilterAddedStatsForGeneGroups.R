#run RunSingleGO.AUROC.Analysis.r first
#don't run as a script - use manually

targetGO <- "MHC protein complex"
targetGO <- "type I interferon signaling pathway"
targetGO <- "ensheathment of neurons"
targetGO <- "cellular response to zinc ion"
targetGO <- "detection of chemical stimulus involved in sensory perception"
targetGO <- "keratinization"
targetGO <- "histone demethylase activity"




#write out all positive
for (targetGO in result.up$MainTitle) {
  targetGOID <- dplyr::filter(tbl_df(geneSetsGO$MODULES), Title == targetGO)$ID
  targetGenes <- filter(geneStatistics, geneSymbol %in% unlist(geneSetsGO$MODULES2GENES[targetGOID]))
  targetGenes %<>% dplyr::arrange(desc(pValueWithDirection))
  write_csv(targetGenes, paste0(baseFilename, ".", targetGO, ".addedStats.csv"))
}

#all neg
for (targetGO in result.down$MainTitle) {
  targetGOID <- dplyr::filter(tbl_df(geneSetsGO$MODULES), Title == targetGO)$ID
  targetGenes <- filter(geneStatistics, geneSymbol %in% unlist(geneSetsGO$MODULES2GENES[targetGOID]))
  targetGenes %<>% dplyr::arrange(pValueWithDirection)
  targetGO <- gsub('/', "", targetGO)
  write_csv(targetGenes, paste0(baseFilename, ".", targetGO, ".addedStats.csv"))
}


#pick one
targetGO <- 'Darmanis.Oligo'
targetGO <- 'Darmanis.Astrocytes'
targetGO <- 'Darmanis.Microglia'
targetGO <- 'Thakurela.Myelin.cortexCompare'
targetGO <- 'Hametner.Iron.associated.human'
targetGO <- 'Tripathy.electro'

targetGOID <- dplyr::filter(tbl_df(geneSetsCellType$MODULES), Title == targetGO)$ID
targetGenes <- filter(geneStatistics, geneSymbol %in% unlist(geneSetsCellType$MODULES2GENES[targetGO]))
write_csv(targetGenes, paste0(baseFilename, ".", targetGO, ".addedStats.csv"))


#disgenet
targetGO <- "CSF lactate increased"
targetGOID <- dplyr::filter(tbl_df(geneSets$MODULES), Title == targetGO)$ID
targetGenes <- filter(geneStatistics, geneSymbol %in% unlist(geneSets$MODULES2GENES[targetGOID]))
tail(targetGenes)


targetGO <- 'Human immunodeficiency virus infectious disease'
targetGO <- "Behcet's disease"
targetGO <- "Stillbirth"
targetGO <- "Angelman syndrome"

targetGOID <- dplyr::filter(tbl_df(geneSetsPhenoCarta$MODULES), Title == targetGO)$ID
targetGenes <- filter(geneStatistics, geneSymbol %in% unlist(geneSetsPhenoCarta$MODULES2GENES[targetGOID]))
write_csv(targetGenes, paste0(baseFilename, ".", targetGO, ".addedStats.csv"))


targetGO <- "L6"
targetGOID <- dplyr::filter(tbl_df(geneSets$MODULES), Title == targetGO)$ID
targetGenes <- filter(geneStatistics, geneSymbol %in% unlist(geneSets$MODULES2GENES[targetGOID]))
write_csv(targetGenes, paste0(baseFilename, ".", targetGO, ".addedStats.csv"))
