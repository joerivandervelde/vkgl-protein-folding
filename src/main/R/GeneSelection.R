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
vkgl <- vkgl[vkgl$UniProtKB.Swiss.Prot.ID %in% protWithOneFrag$V1,]
dim(vkgl)

# Count classifications per gene and apply thresholds
geneCountsTable <- table(vkgl$Classification, vkgl$Gene)
geneCountsDF <- geneCountsTable %>% as.data.frame()
geneCountsCast <- dcast(geneCountsDF, Var2~Var1)
threshold <- 50
genesThr <- subset(geneCountsCast, LB >= threshold)
genesThr <- subset(genesThr, LP >= threshold)

# Proportion of VUS, prioritize genes
genesThr$propVUS <- genesThr$VUS/(genesThr$LB+genesThr$LP)
sorted <- genesThr[order(-genesThr$propVUS),]
sorted

#    Gene  LB  LP VUS   propVUS
#   SCN5A 100 147 495 2.0040486
#   FGFR3  61  56 155 1.3247863
# CACNA1A 145  51 159 0.8112245
#   FGFR1  61 115 133 0.7556818
#    GNAS  84  79 121 0.7423313
#    MLH1 152 176 221 0.6737805
#   SCN1A  87 281 209 0.5679348
#   ATP7B  84 151 127 0.5404255
#    SOS1  53  56  55 0.5045872
#   MUTYH 143 120 113 0.4296578
#   MECP2  82  82  67 0.4085366
#    TSC2 649 537 480 0.4047218
#   BRCA1 596 170 241 0.3146214
#   FGFR2  74 199  75 0.2747253
#    LDLR 114 768 225 0.2551020
#    MAPT  72  69  18 0.1276596 
