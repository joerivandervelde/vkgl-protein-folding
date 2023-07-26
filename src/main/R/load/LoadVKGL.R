###
# VKGL variant classification data
####

setwd(dataDir)
vkglAll <- read.table(file=vkglProtLoc, sep = '\t', header = TRUE)
vkgl <- subset(vkglAll, Gene == geneName)
vkgl$Classification <- revalue(vkgl$Classification, c("LB"="LB/B", "VUS"="VUS", "LP"="LP/P"))
cat(paste("VKGL data has", dim(vkgl)[[1]], "rows\n", sep=" "))
