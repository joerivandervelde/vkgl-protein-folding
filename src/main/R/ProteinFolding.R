################################
# Install packages (only once) #
################################
#install.packages('R.utils')
#install.packages('ggplot2')
#install.packages('reshape2') 
#install.packages('dplyr') 

#################
# Load packages #
#################
library(R.utils)
library(ggplot2)
library(reshape2)
library(dplyr)

##############################
# Set gene name for all following steps
##############################
geneName <- "MEFV"

######################################
# Set gene working dir and resource paths #
######################################
geneWorkingDir <- paste("/Users/joeri/git/vkgl-protein-folding/data/", geneName, sep="")
mkdirs(geneWorkingDir)
setwd(geneWorkingDir)
tmpDir <- paste(getwd(), "tmp", sep="/")
mkdirs(tmpDir)
geneMappingLoc <- "/Applications/AlphaFold2/hgnc-uniprot-mapping.txt"
alphaFoldLoc <- "/Applications/AlphaFold2/UP000005640_9606_HUMAN_v4.tar"
vkglProtLoc <- "/Users/joeri/VKGL/VKGL-prot/VKGL_prot.tsv"
foldx <- "/Applications/FoldX/foldx5MacStd/foldx_20231231" # seems about 2.5x faster than the C11 version

#################################################
# Retrieve mapping of HGNC symbol to UniProt ID #
#################################################
geneMapping <- read.table(file=geneMappingLoc, sep = '\t',header = TRUE)
uniProtID <- geneMapping$UniProtKB.Swiss.Prot.ID[geneMapping$HGNC.symbol==geneName]
uniProtID

###
# VKGL
####
vkglAll <- read.table(file=vkglProtLoc, sep = '\t', header = TRUE)
vkgl <- subset(vkglAll, Gene == geneName)
dim(vkgl)

###
# Find in TAR file, extract into working dir, and unzip
#####
alphaFoldAll <- untar(alphaFoldLoc, list = TRUE)
alphaFoldPDBs <- grep(".pdb.gz",alphaFoldAll,value=TRUE)
PDBForGeneGz <- grep(uniProtID, alphaFoldPDBs,value=TRUE)
PDBForGeneGz
untar(alphaFoldLoc, files = PDBForGeneGz)
PDBForGene <- gunzip(PDBForGeneGz, overwrite=TRUE)[[1]]

#######
# Fix the structure (may take a while)
########
repPDB <- gsub(".pdb", "_Repair.pdb", PDBForGene)
if(length(list.files(pattern=repPDB)) == 0){
  system(paste(foldx, " --command=RepairPDB --pdb=",PDBForGene,sep=""), intern = TRUE)
}
repPDBAbsLoc <- paste(getwd(), repPDB, sep="/")
repPDBAbsLoc

####
# Generate individual_list.txt
###
for(i in 1:nrow(vkgl))
{
  cat(paste("Working on", vkgl[i, 2], "(", i, "of", dim(vkgl)[1],")\n", sep=" "))
  protChangeDir <- paste(tmpDir, vkgl[i, 2], sep="/")
  mkdirs(protChangeDir)
  if(length(list.files(protChangeDir, pattern="*.fxout")) > 0){
    cat("  already folded, skipping...\n")
    next
  }
  cat("  folding...\n")
  file.copy(from = repPDBAbsLoc, to = protChangeDir)
  setwd(protChangeDir)
  write(paste(vkgl[i,2], ";", sep=""), file = "individual_list.txt")
  system(paste(foldx, " --command=BuildModel --mutant-file=individual_list.txt --pdb=", repPDB, sep=""), intern = TRUE)
  file.remove(repPDB)
  cat("...done!\n")
}

#####
# Gather results
######
setwd(tmpDir)
results <- data.frame()
for(i in 1:nrow(vkgl))
{
  cat(paste("Collecting results from ", vkgl[i, 2], "(", i, "of", dim(vkgl)[1],")\n", sep=" "))
  protChangeDir <- paste(tmpDir, vkgl[i, 2], sep="/")
  if(length(list.files(protChangeDir, pattern="*.fxout")) == 0){
    cat("  no results...\n")
    next
  }
  avgDiff <- list.files(protChangeDir, pattern="Average")
  result <- read.table(file = paste(protChangeDir, avgDiff, sep="/"), header = TRUE, skip = 8, sep="\t")
  result$protChange <- vkgl[i, 2]
  result$classification <- vkgl[i, 3]
  results <- rbind(results, result)
}

# drop columns with only 0 and 'Pdb'
results <- results[, colSums(results != 0) > 0]
dropCols <- c("Pdb")
results<- results[ , !(names(results) %in% dropCols)]

# melt dataframe for ggplot
melted <- melt(results, id = c("protChange","classification")) 


####
# plots
###

setwd(geneWorkingDir)

plotdata <- melted %>%
  group_by(classification, variable) %>%
  summarize(n = n(),
            mean = mean(value),
            sd = sd(value),
            se = sd / sqrt(n))

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
        strip.text.x = element_text(size = 6)) +
  labs(x="", 
       y="", 
       title=paste("FoldX terms for ", geneName, " based on VKGL variant classifications", sep=""),
       subtitle = "(Means and standard errors, based on VKGL April 2023 public consensus, FoldX 5.0, and AlphaFold2 human proteome v4)") +
  scale_colour_manual(name = "Classification", values = c("LB" = "#28A014","LP" = "#E41A1C"))
ggsave(paste(geneName,".png",sep=""), width=9, height=5)

