####
# Protein based plot, 1 per attribute
####

setwd(geneWorkingDir)

xmin <- min(mResults$aaLoc)
xmax <- max(mResults$aaLoc)

# Loop over variables and create one plot for each
for(termName in unique(mResults$variable))
{
    
  # select term data and arrange for plot, LP/P on top
  # termName <- "total.energy" # Debug purposes
  cat(paste("Plotting term", termName,"\n", sep=" "))
  selectVar <- mResults[mResults$variable==termName,]
  selectVar <- selectVar %>% arrange(factor(classification, levels = c("VUS","CF","LB/B","LP/P")))
  
  ymin <- min(selectVar$value)
  ymax <- max(selectVar$value)
  
  # Create and save plot
  ggplot() +
    theme_bw() + theme(panel.grid = element_blank(), axis.title.x=element_text(size=10)) +
    geom_point(data = selectVar, aes(x=aaLoc, y=value, colour=classification, shape=source), alpha=1.0, size = 1, stroke = 1) +
    scale_colour_manual(name = "Classification", values = c("LB/B" = "#28A014","VUS" = "#505050","LP/P" = "#E41A1C", "Conflicting" = "orange")) +
    scale_shape_manual(name = "Source", values = c("ClinVar" = 6, "VKGL" = 2, "Both" = 1)) +
    geom_text(data = selectVar, aes(x=aaLoc, y=value, label=protChange), nudge_y=((ymax-ymin)/50), check_overlap = TRUE, alpha=1.0, size = 2) +
    scale_x_continuous(limits = c(xmin,xmax), labels = comma) +
    xlab(paste("Amino acid position within protein alpha chain of", geneName, sep=" ")) +
    ylab(paste("Wild type vs mutant difference in",termName, sep=" ")) +
    ggtitle(label = paste("FoldX", termName, "results for VKGL and ClinVar variants in", geneName, "", sep=" "),
          subtitle = "(based on VKGL public release April 2023, ClinVar 20230702, FoldX 5.0, and AlphaFold2 human proteome v4)")
  ggsave(paste("aapos_",geneName,"_",termName,".pdf",sep=""), width=9, height=5)
  
}
