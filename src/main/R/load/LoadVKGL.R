###
# VKGL variant classification data
####

vkglAll <- read.table(file=vkglProtLoc, sep = '\t', header = TRUE)
vkgl <- subset(vkglAll, Gene == geneName)
vkgl$Classification <- revalue(vkgl$Classification, c("LB"="LB/B", "VUS"="VUS", "LP"="LP/P", "CF"="Conflicting"))
vkgl$source <- "VKGL"
cat(paste("VKGL data has", dim(vkgl)[[1]], "rows\n", sep=" "))
