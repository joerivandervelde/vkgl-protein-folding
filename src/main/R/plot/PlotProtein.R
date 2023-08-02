####
# Protein based plot, 1 per attribute
####

setwd(geneWorkingDir)

xmin <- min(mResults$aaLoc)
xmax <- max(mResults$aaLoc)

# Optionally loop over all variables and create one plot for each
#for(termName in unique(mResults$variable))
#{
    
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
  
  # Create and save plot
  ggplot() +
    theme_bw() + theme(panel.grid = element_blank()) +
    geom_point(data = selectVar, aes(x=aaLoc, y=value, colour=classification), alpha=1.0, size = 1, stroke = 1) + #, shape=source
    scale_colour_manual(name = "VKGL\nclassification", values = c("LB/B" = "#28A014","VUS" = "#505050","LP/P" = "#E41A1C", "Conflicting" = "orange")) +
    #scale_shape_manual(name = "Source", values = c("ClinVar" = 6, "VKGL" = 2, "Both" = 1)) +
    geom_text(data = selectVar, aes(x=aaLoc, y=value, label=protChange), nudge_y=((ymax-ymin)/50), check_overlap = TRUE, alpha=1.0, size = 2) +
    geom_hline(yintercept = youdenIndex) +
    scale_x_continuous(limits = c(xmin,xmax), labels = comma) +
    xlab(paste("Amino acid position within the protein alpha chain of", geneName, sep=" ")) +
    ylab(paste("Difference in Gibbs free energy change (ΔΔG)\nbetween GRCh37 reference and variant protein", sep=" ")) +
    ggtitle(label = paste("Protein folding results for VKGL missense variants in the", geneName, "gene", sep=" "),
          subtitle = paste("Based on VKGL public consensus release April 2023, FoldX 5.0, and AlphaFold2 human proteome v4\nAt a ΔΔG threshold of ",round(youdenIndex, 2), " the PPV is ",round(ppv),"%, the sensitivity is ",round(sens),"%, the NPV is ",round(npv), "% and the specificity is ", round(spec), "%", sep=""))
  ggsave(paste("ddg_vkgl_",geneName,".pdf",sep=""), width=9, height=5, device=cairo_pdf)
  
#}
