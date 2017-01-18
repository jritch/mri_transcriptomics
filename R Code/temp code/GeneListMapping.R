library(homologene)

genesOfInterest <- read.csv("/Users/lfrench/Downloads/test.txt",header=F,stringsAsFactors = F)[,"V1"]

z <- homologene(genes = genesOfInterest, inTax = 10116, outTax = 10090)
write.table(unique(z$`10090`), "/Users/lfrench/Desktop/data/Gene Lists/Mouse/Somata.filtered.Cajigas.txt",quote = F, row.names=F,col.names = F)
z <- homologene(genes = genesOfInterest, inTax = 10116, outTax = 9606)
write.table(unique(z$`9606`), "/Users/lfrench/Desktop/data/Gene Lists/Human/Somata.filtered.Cajigas.txt",quote = F, row.names=F,col.names = F)


library(mygene)
genesOfInterest <- read.csv("/Users/lfrench/Downloads/taylor IDs.txt",header=F,stringsAsFactors = F)[,"V1"]
head(genesOfInterest)
symbols <- queryMany(genesOfInterest,  species="rat", scopes="reporter", fields="symbol")
symbols <- subset(symbols, is.na(notfound))
symbols$symbol
symbols <- symbols[symbols$symbol != toupper(symbols$symbol),]
dim(symbols)
genesOfInterest <- unique(symbols$symbol)
z <- homologene(genes = genesOfInterest, inTax = 10116, outTax = 10090)
write.table(z$`10090`, "/Users/lfrench/Desktop/data/Gene Lists/Mouse/AxonalmRNA.Taylor.txt",quote = F, row.names=F,col.names = F)


genesOfInterest <- read.csv("/Users/lfrench/Desktop/results/CellTypesAging/data/GeneLists/MouseCustom.AdultAxons.Gumy.txt",header=F,stringsAsFactors = F)[,"V1"]
genesOfInterest <- mouse2human(genesOfInterest)
write.table(genesOfInterest$humanGene, "/Users/lfrench/Desktop/results/CellTypesAging/data/GeneLists/Custom.AdultAxons.Gumy.txt",quote = F, row.names=F,col.names = F)
