#/Users/lfrench/Desktop/results/mri_transcriptomics/single_gene_data/9861.ZNF362.cortex.MRI(xyz).expression.csv

allDonors <- NULL
for (file in list.files("/Users/lfrench/Desktop/results/mri_transcriptomics/single_gene_data/", full.names = T)) {
  print(file)
  oneDonor <- read_csv(file)
  oneDonor$donor <- file
  allDonors <- bind_rows(allDonors, oneDonor)
}
allDonors

ggplot(allDonors, aes(x=MRI_Intensity, y = EPCAM)) + geom_point() + geom_smooth(method=lm)  + facet_wrap(~ donor, scales="free")+ theme_bw()
+ stat_density2d(aes(fill = ..level..), geom="polygon") +
ggplot(averageAndAll, aes(x=meanCor, y=averageExpression)) + stat_density2d(aes(fill = ..level..), geom="polygon") + theme_bw()
