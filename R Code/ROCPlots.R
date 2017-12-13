library(dplyr)
library(ggplot2)
library(tidyr)
library(plotROC)


#input 
#sorted genes, gene ID names, tmod object for getting genes <-> group
createPlots <- function(sortedGenes, groupIDs, tmodSets, customNames=NULL, filter = T) {
  
  ranking <- tbl_df(data.frame(gene_symbol = rev(sortedGenes), rank = 1: length(sortedGenes), stringsAsFactors = F))
  geneToClass <- NULL
  for(groupID in groupIDs) {
    geneToClass <- bind_rows(geneToClass, tbl_df(data.frame(gene_symbol = unlist(tmodSets$MODULES2GENES[groupID]), group = tmodSets$MODULES[groupID,]$Title, stringsAsFactors = F)))  
  }
  geneToClass %<>% distinct()
  geneToClassAUC <- left_join(ranking, geneToClass, by = "gene_symbol") %>% spread(key=group, value=group)
  geneToClassAUC %<>% gather(group, present, -gene_symbol, -rank) %>% filter(group != "<NA>")
  geneToClassAUC <- geneToClassAUC %>% mutate(present = if_else(is.na(present), 0, 1))
  
  #geneToClassAUC$dummy <- "True positive fraction" #to create a fake facet for lining things up - no longer needed
  if (length(customNames) > 0){
    names_df = data.frame(ID=groupIDs,customNames)
  }
  
  #order groups by direction of AUC
  forOrder <- tmodUtest(c(sortedGenes), mset=tmodSets, qval = 1, filter = filter)
  forOrder <- tbl_df(forOrder)
  forOrder %<>% filter(ID %in% groupIDs ) %>% arrange(desc(AUC))
  #forOrder %<>% filter(ID %in% groupIDs ) %>% arrange(desc(sign(AUC-0.5)* P.Value))
  geneToClassAUC$group <- factor(geneToClassAUC$group, levels= as.character(forOrder$Title))
  
  throwawayAUC <-geneToClassAUC 
  
  if (length(customNames) > 0){
    forOrder <- forOrder %>% left_join(names_df)
    levels(throwawayAUC$group) <- forOrder$customNames
  }
  
  (AUCPlot <- ggplot(throwawayAUC, aes(d = present, m = rank, color=group)) + ylab("") + 
    geom_roc(n.cuts=0) + 
    style_roc() + coord_cartesian(expand=F) +
    theme(legend.position = c(1,0), legend.justification = c(1, 0), legend.background= element_rect(fill = "transparent", colour = "transparent"), plot.margin=unit(c(.5,.5,.5,.5),"cm")) + 
    labs(color='Gene Group')  + 
    #facet_grid(dummy ~ ., switch="y")  + ylab("") +
    theme(strip.background = element_blank(), strip.placement = "inside", strip.text = element_blank()) )
  
  #geneToClassAUC$group <- save

  forOrder$labelWithAUC <- paste0(tmodSets$MODULES[groupID,]$Title, " (AUROC=", signif(forOrder[groupID, "AUC"],digits=2), ")")

  
  forOrder %<>% mutate(labelWithAUC = paste0(Title, " (AUROC=", signif(AUC,digits=2), ")")) %>% dplyr::select(ID,AUC,group = Title, labelWithAUC)
  
  if (length(customNames > 0)) {
    forOrder <- forOrder %>% left_join(names_df)
    forOrder <- forOrder %>% mutate(labelWithAUC = paste0(customNames, " (AUROC=", signif(AUC,digits=2), ")")) 
  }
  
  forOrder$group <- as.character(forOrder$group)
  geneToClassAUC$group<- as.character(geneToClassAUC$group)
  geneToClassAUC <- inner_join(geneToClassAUC, forOrder, by="group") %>% dplyr::select(-group) %>% dplyr::rename(group = labelWithAUC)

  geneToClassAUC$group <- factor(geneToClassAUC$group, levels= as.character(forOrder$labelWithAUC))

  geneToClassAUC$rank <- -1*geneToClassAUC$rank
  
  #print(geneToClassAUC$group)
 
  throwawayAUC <-geneToClassAUC 
 
  (rasterPlot <- ggplot(throwawayAUC, aes(x = rank, y = present, color= group)) + 
    geom_blank() + 
    geom_vline(data = filter(throwawayAUC, present == 1), aes(xintercept=rank, color=group)) + #,color="black") + #, size=0.07) + 
    theme_bw()+coord_cartesian(expand=F) +
    ylab("Transcriptomic cell type") + 
    facet_wrap(~group, strip.position="top",ncol=1) + #, switch = "both"
    theme(strip.background = element_blank(), strip.placement = "inside") + #, strip.text.y = element_text(angle = 180)) +
    theme(axis.title.y = element_blank(),  axis.text.y=element_blank(), axis.ticks.y=element_blank(),axis.ticks.x=element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + guides(color=FALSE) +
    scale_x_continuous(name = paste0("T1-/T2-w association gene ranking (",length(unique(throwawayAUC$gene_symbol))," genes)"), breaks= c(min(throwawayAUC$rank)+1000, max(throwawayAUC$rank)-1000), labels = c("Positive correlation", "Negative correlation")))
  returnPlots = list()
  returnPlots[["AUCPlot"]] <- AUCPlot
  returnPlots[["rasterPlot"]] <- rasterPlot
  returnPlots
}

#plots <- createPlots(sortedGenes, c("GO:0005882", "GO:0032543", "GO:0000502", "GO:0060337", "GO:0060076", "GO:0044309","GO:0007422", "GO:0042552"), geneSetsGO,customNames = as.character(1:8))

#plot_grid(plots$AUCPlot,plots$rasterPlot)

