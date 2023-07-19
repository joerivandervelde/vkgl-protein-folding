################################
# Install packages (only once) #
################################
#install.packages('R.utils')
#install.packages('ggplot2')
#install.packages('reshape2') 
#install.packages('toprdata')
#install.packages('stringr')
#install.packages('scales')
#install.packages('cutpointr')
#install.packages('plyr')
#install.packages('dplyr')

#################
# Load packages #
#################
library(R.utils)
library(ggplot2)
library(reshape2)
library(toprdata)
library(stringr)
library(scales)
library(cutpointr)
library(plyr)
library(dplyr)

##############################
# Set gene name for all following steps
##############################
geneName <- "FGFR3"

##########################################
# Create gene dirs and link to resources #
##########################################
dataDir <- "/Users/joeri/git/vkgl-protein-folding/data"
geneWorkingDir <- paste(dataDir, geneName, sep="/")
mkdirs(geneWorkingDir)
tmpDir <- paste(geneWorkingDir, "tmp", sep="/")
mkdirs(tmpDir)
geneMappingLoc <- "/Applications/AlphaFold2/hgnc-uniprot-mapping.txt"
alphaFoldLoc <- "/Applications/AlphaFold2/UP000005640_9606_HUMAN_v4.tar"
vkglProtLoc <- "/Users/joeri/VKGL/VKGL-prot/VKGL_apr2023_protForFolding.tsv"
foldx <- "/Applications/FoldX/foldx5MacStd/foldx_20231231" # seems about 2.5x faster than the C11 version

#################################################
# Retrieve mapping of HGNC symbol to UniProt ID #
#################################################
geneMapping <- read.table(file=geneMappingLoc, sep = '\t',header = TRUE)
uniProtID <- geneMapping$UniProtKB.Swiss.Prot.ID[geneMapping$HGNC.symbol==geneName]
uniProtID # if multiple, which one to pick?

###
# VKGL
####
vkglAll <- read.table(file=vkglProtLoc, sep = '\t', header = TRUE)
vkgl <- subset(vkglAll, Gene == geneName)
vkgl$Classification <- revalue(vkgl$Classification, c("LB"="LB/B", "VUS"="VUS", "LP"="LP/P"))
dim(vkgl)

##############
# Find in TAR file, extract into working dir, gunzip and repair (may take a while)
################
setwd(geneWorkingDir)
if(length(list.files(pattern="*_Repair.pdb")) == 0){
  alphaFoldAll <- untar(alphaFoldLoc, list = TRUE)
  alphaFoldPDBs <- grep(".pdb.gz",alphaFoldAll,value=TRUE)
  PDBForGeneGz <- grep(uniProtID, alphaFoldPDBs,value=TRUE)
  PDBForGeneGz # if more than 1, must select correct fragment (todo)
  untar(alphaFoldLoc, files = PDBForGeneGz)
  PDBForGene <- gunzip(PDBForGeneGz, overwrite=TRUE)[[1]]
  system(paste(foldx, " --command=RepairPDB --pdb=",PDBForGene,sep=""), intern = TRUE)
}
repPDB <- list.files(pattern="*_Repair.pdb")
repPDBAbsLoc <- paste(geneWorkingDir, repPDB, sep="/")
repPDBAbsLoc

####
# Generate individual_list.txt and fold
###
for(i in 1:nrow(vkgl))
{
  # i=1  # Debug purposes
  cat(paste("Working on", vkgl[i, 7], "(", i, "of", dim(vkgl)[1],")\n", sep=" "))
  protChangeDir <- paste(tmpDir, vkgl[i, 7], sep="/")
  mkdirs(protChangeDir)
  if(length(list.files(protChangeDir, pattern="*.fxout")) > 0){
    cat("  already folded, skipping...\n")
    next
  }
  if(length(list.files(protChangeDir, pattern="exception.txt")) > 0){
    cat("  already tried before but failed, skipping...\n")
    next
  }
  cat("  folding...\n")
  file.copy(from = repPDBAbsLoc, to = protChangeDir)
  setwd(protChangeDir)
  write(paste(vkgl[i,7], ";", sep=""), file = "individual_list.txt")
  state <- system(paste(foldx, " --command=BuildModel --mutant-file=individual_list.txt --pdb=", repPDB, sep=""), intern = TRUE)
  if(any(grepl("Specified residue not found", state)))
  {
    write(state, file = "exception.txt")
  }
  file.remove(repPDB)
  cat("...done!\n")
}

###
# Cleanup rotabase files generated for each mutation in whole data dir
###
setwd(dataDir)
rotabaseFiles <- list.files(pattern="rotabase.txt", recursive=TRUE)
file.remove(rotabaseFiles)

#####
# Gather results from tmp dir
######
setwd(tmpDir)
results <- data.frame()
for(i in 1:nrow(vkgl))
{
  cat(paste("Collecting results from ", vkgl[i, 7], "(", i, "of", dim(vkgl)[1],")\n", sep=" "))
  protChangeDir <- paste(tmpDir, vkgl[i, 7], sep="/")
  if(length(list.files(protChangeDir, pattern="*.fxout")) == 0){
    cat("  no results...\n")
    next
  }
  avgDiff <- list.files(protChangeDir, pattern="Average")
  result <- read.table(file = paste(protChangeDir, avgDiff, sep="/"), header = TRUE, skip = 8, sep="\t")
  result$assembly <- vkgl[i, 1]
  result$chrom <- vkgl[i, 2]
  result$pos <- vkgl[i, 3]
  result$ref <- vkgl[i, 4]
  result$alt <- vkgl[i, 5]
  result$gene <- vkgl[i, 6]
  result$protChange <- vkgl[i, 7]
  result$classification <- vkgl[i, 8]
  results <- rbind(results, result)
}

# drop columns with only 0 and 'Pdb'
results <- results[, colSums(results != 0) > 0]
dropCols <- c("Pdb")
results<- results[ , !(names(results) %in% dropCols)]

# melt dataframe for ggplot
mResults <- melt(results, id = c("assembly", "chrom", "pos", "ref", "alt", "gene", "protChange","classification")) 
mResultsNoVUS <- mResults[mResults$classification != "VUS",]

####
# Plots
###

setwd(geneWorkingDir)

# Plot 1: overview of all FoldX terms for LB/B and LP/P variants
plotdata <- mResultsNoVUS %>%
  group_by(classification, variable) %>%
  dplyr::summarize(n = n(),
            mean = mean(value),
            median = median(value),
            sd = sd(value),
            se = sd / sqrt(n))
plotdata$variable <- gsub("\\.", "\n", plotdata$variable)
ggplot(plotdata, 
       aes(x = classification, 
           y = mean,
           color = classification)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, 
                    ymax = mean + se),
                width = .1) +
  facet_grid(. ~ variable) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.text.x = element_text(size = 6),
        axis.text.x=element_text(size=7)) +
  labs(x="", 
       y="", 
       title=paste("FoldX terms for ", geneName, " based on VKGL variant classifications", sep=""),
       subtitle = "(Means and standard errors, based on VKGL April 2023 public consensus, FoldX 5.0, and AlphaFold2 human proteome v4)") +
  scale_colour_manual(name = "Classification", values = c("LB/B" = "#28A014","VUS" = "darkgray","LP/P" = "#E41A1C"))
ggsave(paste(geneName,".png",sep=""), width=9, height=5)


######
# Plot 2: gene plot per attribute
######

# Retrieve gene and exon coordinates
geneCoords <- subset(ENSGENES_37, gene_symbol==geneName)
geneChr <- gsub("chr","", geneCoords$chrom)
exonCoords <- subset(ENSEXONS_37, gene_symbol==geneName)
exonStartCoords <- str_split(exonCoords$exon_chromstart, ",")
exonEndCoords <- str_split(exonCoords$exon_chromend, ",")
exons <- data.frame(exonStart = as.numeric(unlist(exonStartCoords)), exonEnd = as.numeric(unlist(exonEndCoords)))

# Loop over variables and create one plot for each
for(termName in unique(mResults$variable))
{

# select the data for this term
# termName <- "total.energy" # Debug purposes
selectVar <- mResults[mResults$variable==termName,]

# arrange for plot, LP/P on top
selectVar <- selectVar %>% arrange(factor(classification, levels = c("VUS","LB/B","LP/P")))

# Determine optimal threshold using Youden's Index #
cutpointDF <- subset(selectVar, classification != "VUS")
opt_cut <- cutpointr(cutpointDF, value, classification, direction = ">=", pos_class = "LP/P", neg_class = "LB/B", method = maximize_metric, metric = youden)
youdenIndex <- opt_cut$optimal_cutpoint
tp <- sum(cutpointDF[cutpointDF$classification=="LP/P",'value'] >= youdenIndex)
fp <- sum(cutpointDF[cutpointDF$classification=="LB/B",'value'] >= youdenIndex)
ppv <- 100 *tp/(tp+fp)
sens <- opt_cut$sensitivity*100

# Determine plot window
xmin <- min(exons$exonStart)
xmax <- max(exons$exonEnd)
ymin <- min(selectVar$value)
ymax <- max(selectVar$value)

# Create and save plot
ggplot() +
  theme_bw() + theme(panel.grid = element_blank(), axis.title.x=element_text(size=10)) +
  geom_rect(data = exons, aes(xmin = exonStart, xmax = exonEnd, ymin = ymin, ymax = ymax), linetype = 0, fill="lightgray", alpha = 1) +
  geom_point(data = selectVar, aes(x=pos, y=value, colour=classification), alpha=1.0, size = 1, stroke = 1) +
  geom_text(data = selectVar, aes(x=pos, y=value, label=protChange), nudge_y=((ymax-ymin)/50), check_overlap = TRUE, alpha=1.0, size = 2) +
  geom_hline(yintercept = youdenIndex) +
  scale_colour_manual(name = "Classification", values = c("LB/B" = "#28A014","VUS" = "#505050","LP/P" = "#E41A1C")) +
  scale_x_continuous(limits = c(xmin,xmax), labels = comma) +
  xlab(paste("",geneName," at GRCh37 chr",geneChr,":", xmin, "-",xmax,", lightgray: exons", sep="")) +
  ylab(termName) +
  ggtitle(paste("FoldX results for ",geneName,". At a threshold of ",round(youdenIndex, 2), " the PPV is ",round(ppv),"% and the sensitivity is ",round(sens),"%.",sep=""))
ggsave(paste(geneName,"_",termName,".png",sep=""), width=9, height=5)

}

