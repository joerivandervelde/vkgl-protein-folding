
clinVar <- read.table(file=clinVarLoc, sep = '\t', header = TRUE)
clinVar$Classification <- revalue(clinVar$Classification, c("LB"="LB/B", "VUS"="VUS", "LP"="LP/P", "CF"="Conflicting"))
clinVar <- subset(clinVar, Gene == geneName)
clinVar$source <- "ClinVar"
cat(paste("ClinVar data has", dim(clinVar)[[1]], "rows\n", sep=" "))
