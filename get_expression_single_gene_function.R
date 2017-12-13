get_expression_single_gene <- function(geneOfInterest, regionOfInterest, modalityOfInterest, single_gene_folder) {
  allDonors <- NULL
  for (file in list.files(single_gene_folder, pattern=paste0("[.]", geneOfInterest,"[.]",regionOfInterest, ".*", modalityOfInterest, "[.]"), full.names = T)) {
    print(file)
    oneDonor <- as.data.frame(read_csv(file))
    oneDonor$Expression <- oneDonor[,geneOfInterest]
    oneDonor$donor <- paste("Donor", gsub("[.].*","", basename(file)))
    allDonors <- bind_rows(allDonors, oneDonor)
  }
  allDonors <- tbl_df(allDonors)
  allDonors
}
