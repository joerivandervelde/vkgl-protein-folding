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
#install.packages('seqminer')
#install.packages('Cairo')


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
library(Cairo)


##############################
# Set gene name(s) for all following steps
##############################
genes <- c("MEFV", "CFTR", "MECP2", "TERT", "CACNA1A",
           "MLH1", "SOS1", "FGFR3", "ATP7B", "SCN5A",
           "MUTYH", "MAPT", "FGFR1", "LDLR", "SCN1A",
           "BRCA1", "FGFR2", "GNAS", "PLCG2", "ABCG8",
           "SCN10A", "KCNH2", "ABCA1", "SCN8A", "PTEN",
           "MET", "RAF1", "F8", "TSC2", "BEST1")
# Keep track of results per gene
columns = c("gene","nbenign","npatho","threshold","ppv","npv","sens","spec","foldingSuccessRate") 
geneResults = data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(geneResults) = columns


##############################
# Project dirs and resources #
##############################
rootDir <- "/Users/joeri/git/vkgl-protein-folding"
resourcesDir <- paste(rootDir,"resources", sep="/")
scriptDir <- paste(rootDir, "src/main/R", sep="/")
dataDir <- paste(rootDir, "data", sep="/")
outputsDir <- paste(rootDir, "out", sep="/")
mkdirs(outputsDir)

geneMappingLoc <- paste(resourcesDir, "hgnc-uniprot-mapping.txt", sep="/")
vkglProtLoc <- paste(resourcesDir, "VKGL_apr2023_protForFolding.tsv", sep="/")
alphaFoldLoc <- "/Applications/AlphaFold2/UP000005640_9606_HUMAN_v4.tar"
foldx <- "/Applications/FoldX/foldx5MacStd/foldx_20231231" # seems about 2.5x faster than the C11 version
#clinVarLoc <- "/Users/joeri/ClinVar/clinvar_20230722_protForFolding.tsv"

for (geneName in genes)
{
# geneName <- "CFTR" # To try out new genes


############################
# Create gene and tmp dirs #
############################
geneWorkingDir <- paste(dataDir, geneName, sep="/")
mkdirs(geneWorkingDir)
tmpDir <- paste(geneWorkingDir, "tmp", sep="/")
mkdirs(tmpDir)


#################################################
# Retrieve mapping of HGNC symbol to UniProt ID #
#################################################
geneMapping <- read.table(file=geneMappingLoc, sep = '\t',header = TRUE)
uniProtID <- geneMapping$UniProtKB.Swiss.Prot.ID[geneMapping$HGNC.symbol==geneName]
uniProtID
if(identical(uniProtID, character(0))){
  stop(paste("No UniProt ID found for gene", geneName))
}


####################################################################################
# Find in TAR file, extract into working dir, gunzip and repair (may take a while) #
####################################################################################
setwd(geneWorkingDir)
if(length(list.files(pattern="*_Repair.pdb")) == 0){
  alphaFoldAll <- untar(alphaFoldLoc, list = TRUE)
  alphaFoldPDBs <- grep(".pdb.gz",alphaFoldAll, value=TRUE)
  PDBForGeneGz <- grep(paste(uniProtID,collapse="|"), alphaFoldPDBs, value=TRUE)
  if(identical(PDBForGeneGz, character(0))){
    stop(paste("No PDB found for uniprot "), paste(uniProtID,collapse = " "))
  }
  if(length(PDBForGeneGz) > 1){
    stop(paste("Multiple PDB and/or fragments found for uniprot "), paste(uniProtID, collapse = " "), ": ", paste(PDBForGeneGz, collapse=" "))
  }
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
#setwd(scriptDir)
#source("load/LoadGenomicCoordinates.R") # Coordinates of genes and exons for build 37 and 38
setwd(scriptDir)
source("load/LoadVKGL.R") # Load VKGL data (B37)
#setwd(scriptDir)
#source("load/LoadClinVar.R") # Load ClinVar data (B38)
variants <- vkgl
# Optionally merge with ClinVar for more variants
#variants <- merge(x = vkgl, y = clinVar, all = TRUE, by = "ProtChange")
#variants$rowSource <- apply(variants[c("source.x", "source.y")], 1, function(x) paste(na.omit(x), collapse = ""))


#########################################
# Generate individual_list.txt and fold #
#########################################
for(i in 1:nrow(variants))
{
  #i=1  # Debug purposes
  cat(paste("Working on", variants[i, "ProtChange"], "(", i, "of", dim(variants)[1],")\n", sep=" "))
  protChangeDir <- paste(tmpDir, variants[i, "ProtChange"], sep="/")
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
  write(paste(variants[i, "ProtChange"], ";", sep=""), file = "individual_list.txt")
  state <- system(paste(foldx, " --command=BuildModel --mutant-file=individual_list.txt --pdb=", repPDB, sep=""), intern = TRUE)
  if(any(grepl("Specified residue not found", state)))
  {
    write(state, file = "exception.txt")
  }
  file.remove(repPDB)
  cat("...done!\n")
}

#############################################################
# Cleanup plots and files that can be removed after folding #
#############################################################
setwd(dataDir)
rotabaseFiles <- list.files(pattern="rotabase.txt", recursive=TRUE)
file.remove(rotabaseFiles)

setwd(geneWorkingDir)
pdfFiles <- list.files(pattern="*.pdf", recursive=TRUE)
file.remove(pdfFiles)
moleculesDir <- list.files(pattern="molecules", recursive=TRUE, include.dirs=TRUE)
file.remove(moleculesDir)

setwd(tmpDir)
pdbFiles <- list.files(pattern="*.pdb", recursive=TRUE)
file.remove(pdbFiles)
indvFiles <- list.files(pattern="*individual_list.txt", recursive=TRUE)
file.remove(indvFiles)


############################################
# Gather folding results from gene tmp dir #
############################################
setwd(tmpDir)
results <- data.frame()
for(i in 1:nrow(variants))
{
  #i=1  # Debug purposes
  cat(paste("Collecting results from ", variants[i, "ProtChange"], "(", i, "of", dim(variants)[1],")\n", sep=" "))
  protChangeDir <- paste(tmpDir, variants[i, "ProtChange"], sep="/")
  if(length(list.files(protChangeDir, pattern="*.fxout")) == 0){
    cat("  no results...\n")
    next
  }
  avgDiff <- list.files(protChangeDir, pattern="Average")
  result <- read.table(file = paste(protChangeDir, avgDiff, sep="/"), header = TRUE, skip = 8, sep="\t")
  result$protChange <- variants[i, "ProtChange"]
  result$classificationVKGL <- variants[i, "Classification"]
  # if merged with ClinVar, use below instead
  #result$classificationVKGL <- variants[i, "Classification.x"]
  #result$classificationClinVar <- variants[i, "Classification.y"]
  #result$source <- variants[i, "rowSource"]
  results <- rbind(results, result)
}

#################################################
# Check folding success rate and stop if needed #
#################################################
foldingSuccessRate = nrow(results)/nrow(variants)
if(foldingSuccessRate < 0.5)
{
  #stop(paste("Folding success rate for ",geneName," too low (",round(successRate, 2),", ",nrow(results)," out of ",nrow(variants),")",sep=""))
}


#########################################
# Merge classifications into one column #
#########################################
# Only applies if rows have a source (after ClinVar merge)
if(!is.null(results$source)){
  results$classification <- "NA"
  results["source"][results["source"] == "VKGLClinVar"] <- "Both"
  for (i in 1:nrow(results)) {
    if(is.na(results[i,]$classificationVKGL)){
      results[i,]$classification <- results[i,]$classificationClinVar
    }else if(is.na(results[i,]$classificationClinVar)){
      results[i,]$classification <- results[i,]$classificationVKGL
    }else if(results[i,]$classificationVKGL == results[i,]$classificationClinVar){
      results[i,]$classification <- results[i,]$classificationVKGL
    }else{
      results[i,]$classification <- "Conflicting"
    }
  }
}


#######################################################################################
# Postprocess results: drop 0-sum and 'Pdb' columns, melt for ggplot, add AA position #
#######################################################################################
results <- results[, colSums(results != 0, na.rm = TRUE) > 0]
results <- results[ , !(names(results) %in% c("Pdb"))]
results$aaLoc <- gsub("[A-Z]", "", results$protChange)
results$aaLoc <- as.numeric(results$aaLoc)
mResults <- melt(results, id = c("aaLoc","protChange","classificationVKGL"))
# If merged with ClinVar, use below instead
# mResults <- melt(results, id = c("aaLoc","protChange","source","classificationVKGL","classificationClinVar","classification"))
# Finally, if not merged, assign 'classification' directly from VKGL
if(is.null(results$source)){
  mResults$classification <- mResults$classificationVKGL
}

########################
# Create various plots #
########################
# Often not so interesting, skip
#setwd(scriptDir)
#source("plot/PlotOverview.R")

setwd(scriptDir)
source("plot/PlotProtein.R")

# Ability to build a 2D predictor using KDE, overfits obviously
# setwd(scriptDir)
# source("plot/PlotPredictProtein.R")

# Obsolete, difficult anyway because of combined b37/b38 data using ClinVar
#setwd(scriptDir)
#source("plot/PlotGene.R")

# export any data
#write.table(results, sep="\t",file="cftr_vkgl_clinvar_foldx_af2_results.txt", quote=FALSE, row.names =FALSE)

}

geneResults

setwd(outputsDir)
write.table(geneResults, sep="\t",file="ddg_vkgl_gene_results.txt", quote=FALSE, row.names =FALSE)
