####
# Prep protein based plot, 1 per attribute
####

setwd(geneWorkingDir)

xmin <- min(mResults$aaLoc)
xmax <- max(mResults$aaLoc)

# select term data and arrange for plot, LP/P on top
termName <- "total.energy" # Select term manually
cat(paste("Plotting term", termName,"\n", sep=" "))
selectVar <- mResults[mResults$variable==termName,]
selectVar <- selectVar %>% arrange(factor(classification, levels = c("VUS","CF","LB/B","LP/P")))

# Determine optimal threshold using Youden's Index #
cutpointDF <- subset(selectVar, classification == "LB/B" | classification == "LP/P")
opt_cut <- cutpointr(cutpointDF, value, classification, direction = ">=", pos_class = "LP/P", neg_class = "LB/B", method = maximize_metric, metric = youden)
youdenIndex <- opt_cut$optimal_cutpoint
tp <- sum(cutpointDF[cutpointDF$classification=="LP/P",'value'] >= youdenIndex)
fp <- sum(cutpointDF[cutpointDF$classification=="LB/B",'value'] >= youdenIndex)
tn <- sum(cutpointDF[cutpointDF$classification=="LB/B",'value'] < youdenIndex)
fn <- sum(cutpointDF[cutpointDF$classification=="LP/P",'value'] < youdenIndex)
ppv <- 100 *tp/(tp+fp)
npv <- 100 *tn/(tn+fn)
sens <- opt_cut$sensitivity*100
spec <- opt_cut$specificity*100

ymin <- min(selectVar$value)
ymax <- max(selectVar$value)

geneResults <- rbind(geneResults, list(gene=geneName, nbenign=sum(selectVar[,"classification"]=="LB/B"), npatho=sum(selectVar[,"classification"]=="LP/P"), threshold=youdenIndex, ppv=ppv, npv=npv, sens=sens, spec=spec, foldingSuccessRate=foldingSuccessRate))
