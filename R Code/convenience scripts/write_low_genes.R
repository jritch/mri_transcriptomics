
write.table(head(geneStatistics %>% arrange(averageExpression), n=138)$geneSymbol, "/Users/lfrench/Downloads/lowgenes138.human.txt", col.names = F, row.names = F, quote=F)
write.table(head(geneStatistics %>% arrange(averageExpression), n=10)$geneSymbol, "/Users/lfrench/Downloads/lowgenes10.human.txt", col.names = F, row.names = F, quote=F)
