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
library(MASS)
library(tidyr)
library(seqminer)


##############################
# Set gene name for all following steps
##############################
geneName <- "CFTR"


##########################################
# Create gene dirs and link to resources #
##########################################
rootDir <- "/Users/joeri/git/vkgl-protein-folding/"
scriptDir <- paste(rootDir, "src/main/R", sep="/")
dataDir <- paste(rootDir, "data", sep="/")
geneWorkingDir <- paste(dataDir, geneName, sep="/")
mkdirs(geneWorkingDir)
tmpDir <- paste(geneWorkingDir, "tmp", sep="/")
mkdirs(tmpDir)
geneMappingLoc <- "/Applications/AlphaFold2/hgnc-uniprot-mapping.txt"
alphaFoldLoc <- "/Applications/AlphaFold2/UP000005640_9606_HUMAN_v4.tar"
vkglProtLoc <- "/Users/joeri/VKGL/VKGL-prot/VKGL_apr2023_protForFolding.tsv"
foldx <- "/Applications/FoldX/foldx5MacStd/foldx_20231231" # seems about 2.5x faster than the C11 version
clinVarLoc <- "/Applications/ClinVar/clinvar_20230702.vcf.gz"


#################################################
# Retrieve mapping of HGNC symbol to UniProt ID #
#################################################
geneMapping <- read.table(file=geneMappingLoc, sep = '\t',header = TRUE)
uniProtID <- geneMapping$UniProtKB.Swiss.Prot.ID[geneMapping$HGNC.symbol==geneName]
uniProtID # if multiple, which one to pick?


####################################################################################
# Find in TAR file, extract into working dir, gunzip and repair (may take a while) #
####################################################################################
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


##################################################
# Load genome coordinates and clinical variation #
##################################################
setwd(scriptDir)
source("load/LoadGenomicCoordinates.R") # Coordinates of genes and exons for build 37 and 38
setwd(scriptDir)
source("load/LoadVKGL.R") # Load VKGL data (B37)
setwd(scriptDir)
source("load/LoadClinVar.R") # Load ClinVar data (B38)


#########################################
# Generate individual_list.txt and fold #
#########################################
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

######################################################
# Cleanup all rotabase files and gene-specific plots #
######################################################
setwd(dataDir)
rotabaseFiles <- list.files(pattern="rotabase.txt", recursive=TRUE)
file.remove(rotabaseFiles)
setwd(geneWorkingDir)
pngFiles <- list.files(pattern="*.png", recursive=TRUE)
file.remove(pngFiles)
pdfFiles <- list.files(pattern="*.pdf", recursive=TRUE)
file.remove(pdfFiles)

############################################
# Gather folding results from gene tmp dir #
############################################
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

########################################################################################
# Postprocess results: drop 0/'Pdb' columns, melt for ggplot, and make VUSless version #
########################################################################################
results <- results[, colSums(results != 0) > 0]
dropCols <- c("Pdb")
results<- results[ , !(names(results) %in% dropCols)]
mResults <- melt(results, id = c("assembly", "chrom", "pos", "ref", "alt", "gene", "protChange","classification"))
mResultsNoVUS <- mResults[mResults$classification != "VUS",]


########################
# Create various plots #
########################
setwd(scriptDir)
source("plot/PlotOverview.R")

setwd(scriptDir)
source("plot/PlotProtein.R")

setwd(scriptDir)
source("plot/PlotGene.R")

