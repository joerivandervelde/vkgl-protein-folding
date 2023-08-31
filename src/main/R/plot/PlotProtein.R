####
# Create and save protein based plot
####

setwd(outputsDir)
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
