library(plyr)
library(dplyr)
library(reshape2)

# Read VKGL variant classification data (public consensus)
vkglProtLoc <- "/Users/joeri/VKGL/VKGL-prot/VKGL_apr2023_protForFolding.tsv"
vkgl <- read.table(file=vkglProtLoc, sep = '\t', header = TRUE)
dim(vkgl)

# Gene mapping
geneMappingLoc <- "/Applications/AlphaFold2/hgnc-uniprot-mapping.txt"
geneMapping <- read.table(file=geneMappingLoc, sep = '\t',header = TRUE)

# Merge with UniProt IDs
vkgl <- merge(vkgl, geneMapping, by.x="Gene", by.y="HGNC.symbol")

# Genes with 1 protein fragment
protWithOneFragLoc <- "/Applications/AlphaFold2/misc/uniprotIds_singleFragment_only.txt"
protWithOneFrag <- read.table(file=protWithOneFragLoc)

# Leave only genes with 1 fragment
vkgl <- vkgl[!(vkgl$Gene %in% protWithOneFrag$V1),]
dim(vkgl)

# Count classifications per gene and apply thresholds
threshold <- 100
geneCountsTable <- table(vkgl$Classification, vkgl$Gene)
geneCounts <- geneCountsTable %>% as.data.frame()
geneCountsCast <- dcast(geneCounts, Var2~Var1)
sub <- subset(geneCountsCast, LB > threshold)
sub2 <- subset(geneCountsCast, LP > threshold)

# Proportion of VUS, prioritize genes
sub2$propVUS <- sub2$VUS/(sub2$LB+sub2$LP)
sorted <- sub2[order(-sub2$propVUS),]
sorted
