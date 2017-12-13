#testing Mistry lists with p fdr and hypergeometric

#requires geneStatistics table
geneSetsCellType

posCorrel <- filter(geneStatistics, isSig & pValueWithDirection > 0) %>% .$geneSymbol
negCorrel <- filter(geneStatistics, isSig & pValueWithDirection < 0) %>% .$geneSymbol

#GO group tests
result <- tbl_df(tmodHGtest(fg = posCorrel, bg = geneStatistics$geneSymbol, mset = geneSetsGO, filter = T, qval = 1.1))
result <- tbl_df(tmodHGtest(fg = negCorrel, bg = geneStatistics$geneSymbol, mset = geneSetsGO, filter = T, qval = 1.1))

result <- tbl_df(tmodHGtest(fg = geneSetsCellType$MODULES2GENES$`Mistry.pH downregulated`, bg = geneStatistics$geneSymbol, mset = geneSetsGO, filter = T, qval = 1.1))
myelinResult <- filter(result, grepl("myelin|ensheathment",Title),!grepl("peripheral|sphingomyelin",Title))
myelinResult$adj.P.Value <- p.adjust(myelinResult$P.Value, method="BH")
myelinResult$adj.P.Value <- signif(myelinResult$adj.P.Value, digits=3)

result

result <- tbl_df(tmodHGtest(fg = posCorrel, bg = geneStatistics$geneSymbol, mset = geneSetsCellType, filter = F, qval = 1.1))
result %>% filter(grepl("Mistry", Title)) %>% mutate(adj.P.Val = p.adjust(P.Value))

result <- tbl_df(tmodHGtest(fg = negCorrel, bg = geneStatistics$geneSymbol, mset = geneSetsCellType, filter = F, qval = 1.1))
result %>% filter(grepl("Mistry", Title)) %>% mutate(adj.P.Val = p.adjust(P.Value))
