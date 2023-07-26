
setwd(dataDir)
geneClinVarData <- tabix.read(clinVarLoc, geneTabixB38)
geneClinVarDF <- read.table(text=geneClinVarData, sep = '\t', quote="", fill=FALSE, header = FALSE, fileEncoding = "UTF-16LE", col.names = c("chr","pos","id","ref","alt","qual","filter","info"))
# Only missense
clinvar <- subset(geneClinVarDF, grepl("missense_variant", geneClinVarDF$info, fixed=TRUE))
# Remove any conflicts and records without a classification
clinvar <- subset(clinvar, !grepl("Conflicting_interpretations_of_pathogenicity", geneClinVarDF$info, fixed=TRUE))
clinvar <- subset(clinvar, !grepl("no_assertion_provided", clinvar$info, fixed=TRUE))
clinvar <- subset(clinvar, !grepl("drug_response", clinvar$info, fixed=TRUE))
clinvar <- subset(clinvar, !grepl("no_assertion_criteria_provided", clinvar$info, fixed=TRUE))
clinvar <- subset(clinvar, !grepl("no_interpretation_for_the_single_variant", clinvar$info, fixed=TRUE))
# Relabel classifications into 3 groups
clinvar$classification[grepl("benign", clinvar$info, ignore.case = TRUE)] <- "LB/B"
clinvar$classification[grepl("uncertain_significance", clinvar$info, ignore.case = TRUE)] <- "VUS"
clinvar$classification[grepl("pathogenic", clinvar$info, ignore.case = TRUE)] <- "LP/P"
clinvar$Source <- "ClinVar"
clinvar$uniProtVar <- str_extract(clinvar$info, "VAR_[:digit:]+")
# there should be 0 records left without a classification, this is a way to check:
# subset(clinvar, clinvar$classification==FALSE)
cat(paste("ClinVar data has", dim(clinvar)[[1]], "rows\n", sep=" "))

