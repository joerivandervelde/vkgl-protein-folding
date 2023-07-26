# Retrieve gene and exon coordinates
geneCoordsB37 <- subset(ENSGENES_37, gene_symbol==geneName)
geneChrB37 <- gsub("chr","", geneCoordsB37$chrom)
exonCoordsB37 <- subset(ENSEXONS_37, gene_symbol==geneName)
exonStartCoordsB37 <- str_split(exonCoordsB37$exon_chromstart, ",")
exonEndCoordsB37 <- str_split(exonCoordsB37$exon_chromend, ",")
exonsB37 <- data.frame(exonStart = as.numeric(unlist(exonStartCoordsB37)), exonEnd = as.numeric(unlist(exonEndCoordsB37)))
cat(paste("Build 37 gene location is ",  geneChrB37, ":", geneCoordsB37$gene_start, "-", geneCoordsB37$gene_end, " with ", dim(exonsB37)[[1]], " exons\n",sep=""))


geneCoordsB38 <- subset(ENSGENES, gene_symbol==geneName)
geneChrB38 <- gsub("chr","", geneCoordsB38$chrom)
exonCoordsB38 <- subset(ENSEXONS, gene_symbol==geneName)
exonStartCoordsB38 <- str_split(exonCoordsB38$exon_chromstart, ",")
exonEndCoordsB38 <- str_split(exonCoordsB38$exon_chromend, ",")
exonsB38 <- data.frame(exonStart = as.numeric(unlist(exonStartCoordsB38)), exonEnd = as.numeric(unlist(exonEndCoordsB38)))
cat(paste("Build 38 gene location is ",  geneChrB38, ":", geneCoordsB38$gene_start, "-", geneCoordsB38$gene_end, " with ", dim(exonsB38)[[1]], " exons\n",sep=""))

geneTabixB38 <- paste(geneChrB38, paste(geneCoordsB38$gene_start, geneCoordsB38$gene_end, sep="-"), sep=":")
