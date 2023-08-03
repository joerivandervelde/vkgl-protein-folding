library(plyr)
library(dplyr)
library(reshape2)

fullStrGenes <- c("ATM","RYR2","SCN5A","NF1","ABCA1","ABCA4","ABCG8","DYNC1H1","FAT4","VCL","SCN10A","PCSK9","MSH2","TNNT2","CFTR","SCN9A","SCN1A","ABCC8","KCNH2","IFT140","HUWE1","NLRP3","MVK","CACNA1D","GPD1","SCN8A","TPM1","ERCC6","POLE","SCN2A","SRCAP","TRPM4","RANBP2","ABCG5","USP9X","KCNQ1","SCN4A","TNNI3","ERCC2","TG","SCN3A","F8","ARID1A","POLD1","INSR","ATR","PTEN","PTCH1","ABCD1","CACNA1B","TRRAP","SNRNP200","DEPDC5","NPC1","WDR19","TAF1","PRKDC","MYL3","GSN","TERT","SCN7A","GHSR","CDK5RAP2","TUBGCP6","NALCN","LRRK2","CACNA1G","BBS9","ATP1A3","TNNC1","TMEM67","RHO","POLR3A","ACTG1","UNC80","TSC1","SMARCA4","SLC12A6","SLC12A3","SCN1B","TCAP","SLC2A1","PRPF8","MTOR","KCNQ4","KCNJ11","INTS1","DNMT3B","CTNNA1","BRAF","ABCB4","KARS1","GBA1","CTC1","BEST1","TPP1","SKIC3","PIK3R1","IDUA","GJB2","CACNA1E","SLC5A2","PMM2","PIK3CA","MET","HARS1","DNM1","CP","CNGA3","CLCN1","CHD4","WDR35","TWNK","SAR1B","RAF1","POLR3B","PEPD","MC1R","KCNJ2","CBS","ATP8B1","SLC3A1","SGSH","PKD2","IMPDH1","IFT122","HK1","FUS","EGFR","C5","TTR","SRRM2","SMC1A","SMARCC2","PC","MMACHC","KIFBP","KCND3","FIG4","DHX38","CUL7","CRYAB","CFB","ATP13A2","ANK1","XDH","WNK3","TAP2","SKIC2","LSS","FANCA","DICER1","CPSF1","CPS1","ARID2","ABCB6","ABCB11","ABCA3")
geneDF <- data.frame(list(Gene=fullStrGenes))

vkglProtLoc <- "/Users/joeri/VKGL/VKGL-prot/VKGL_apr2023_protForFolding.tsv"
vkgl <- read.table(file=vkglProtLoc, sep = '\t', header = TRUE)
geneDF$hasVKGLVEPprotannot <- geneDF$Gene %in% vkgl$Gene

geneMappingLoc <- "/Applications/AlphaFold2/hgnc-uniprot-mapping.txt"
geneMapping <- read.table(file=geneMappingLoc, sep = '\t',header = TRUE)

vkgl <- merge(vkgl, geneMapping, by.x="Gene", by.y="HGNC.symbol")
geneDF$hasUniProtID<- geneDF$Gene %in% vkgl$Gene

protWithOneFragLoc <- "/Applications/AlphaFold2/misc/uniprotIds_singleFragment_only.txt"
protWithOneFrag <- read.table(file=protWithOneFragLoc)
vkgl <- vkgl[vkgl$UniProtKB.Swiss.Prot.ID %in% protWithOneFrag$V1,]
geneDF$hasOneProtFragment <- geneDF$Gene %in% vkgl$Gene

geneCountsTable <- table(vkgl$Classification, vkgl$Gene)
geneCountsDF <- geneCountsTable %>% as.data.frame()
geneCountsCast <- dcast(geneCountsDF, Var2~Var1)

threshold <- 5
genesThr <- subset(geneCountsCast, LB >= threshold)
genesThr <- subset(genesThr, LP >= threshold)
geneDF$hasFiveOrMoreLBBAndLPP <- geneDF$Gene %in% genesThr$Var2

threshold <- 10
genesThr <- subset(geneCountsCast, LB >= threshold)
genesThr <- subset(genesThr, LP >= threshold)
geneDF$hasTenOrMoreLBBAndLPP <- geneDF$Gene %in% genesThr$Var2

threshold <- 20
genesThr <- subset(geneCountsCast, LB >= threshold)
genesThr <- subset(genesThr, LP >= threshold)
geneDF$hasTwentyOrMoreLBBAndLPP <- geneDF$Gene %in% genesThr$Var2

rootDir <- "/Users/joeri/git/vkgl-protein-folding"
setwd(rootDir)
write.table(geneDF, sep="\t",file="geneselect.txt", quote=FALSE, row.names =FALSE)
