symbol_to_entrez <- geneStatistics %>% dplyr::select(ID, entrez_id) %>% filter(!is.nan(entrez_id)) %>% filter(!is.na(entrez_id)) %>% distinct()
symbol_to_entrez %<>% mutate(new_symbol = getSYMBOL(as.character(entrez_id), data='org.Hs.eg')) 

probeMap <- read_csv("/Users/lfrench/Desktop/data/Allen/HBA/normalized_microarray_donor10021/Probes.csv") %>% dplyr::select(ID = gene_symbol, entrez_id) %>% distinct()
genes <- read_tsv("/Users/lfrench/Desktop/results/FreeSurferMapDataReport/files/AllenHBA_DK_ExpessionMatrix.15496.tsv")$X1
probeMap %<>% filter(gene_symbol %in% genes)
symbol_to_entrez <- probeMap %>% filter(!is.nan(entrez_id)) %>% filter(!is.na(entrez_id)) %>% distinct()

symbol_to_entrez %<>% mutate(new_symbol = getSYMBOL(as.character(entrez_id), data='org.Hs.eg')) %>% dplyr::rename(oldGeneSymbol=ID)
write_csv(symbol_to_entrez, "/Users/lfrench/Downloads/SymbolUpdateMap.csv")
