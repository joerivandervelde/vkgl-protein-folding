####
# Protein based plot, 1 per attribute
####

setwd(geneWorkingDir)

# Variants ONLY present in ClinVar
mResultsClinVarExclusive <-subset(mResults, source == "ClinVar")

# Variants seen in VKGL, may also have been present in ClinVar
mResultsVKGLInclusive <- subset(mResults, source == "VKGL" | source == "Both")

xmin <- min(mResultsVKGLInclusive$aaLoc)
xmax <- max(mResultsVKGLInclusive$aaLoc)

# Loop over variables and create one plot for each
for(termName in unique(mResults$variable))
{
  
  # termName <- "total.energy" # Debug purposes
  cat(paste("Plotting term", termName,"\n", sep=" "))
  
  # select VKGL term data and arrange for plot, LP/P on top
  selectVarVKGL <- mResultsVKGLInclusive[mResultsVKGLInclusive$variable==termName,]
  selectVarVKGL <- selectVarVKGL %>% arrange(factor(classification, levels = c("VUS","CF","LB/B","LP/P")))
  
  selectVarClinVar <- mResultsClinVarExclusive[mResultsClinVarExclusive$variable==termName,]
  selectVarClinVarLP <- subset(selectVarClinVar, classification == "LP/P")
  selectVarClinVarLB <- subset(selectVarClinVar, classification == "LB/B")
  
  xRes = 200
  yRes = 200
  
  everythingIsFine <- TRUE
  tryCatch({
    smoothLP <- MASS::kde2d(selectVarClinVarLP$aaLoc, selectVarClinVarLP$value, n=c(xRes, yRes),lims=c(range(selectVarClinVar$aaLoc), range(selectVarClinVar$value))) 
  }, error = function(e) {
    cat(paste("There was a problem for computing LP densities for term",termName,"\n",sep=" "))
    everythingIsFine <- FALSE
  })
  
  tryCatch({
    smoothLB <- MASS::kde2d(selectVarClinVarLB$aaLoc, selectVarClinVarLB$value, n=c(xRes, yRes),lims=c(range(selectVarClinVar$aaLoc), range(selectVarClinVar$value))) 
  }, error = function(e) {
    cat(paste("There was a problem for computing LB densities for term",termName,"\n",sep=" "))
    everythingIsFine <- FALSE
  })
  
  if (!everythingIsFine) next
  
  smoothDiff <- list(x = smoothLP$x, y = smoothLP$y, z = smoothLP$z-smoothLB$z)
  
  df_smoothedDiff <- smoothDiff$z %>% 
    as_tibble() %>%
    pivot_longer(cols = everything(), names_to = "col", values_to = "val") %>% 
    mutate(aaLoc = rep(smoothDiff$x, each=yRes),
           value = rep(smoothDiff$y, xRes))
  
  ymin <- min(selectVarVKGL$value)
  ymax <- max(selectVarVKGL$value)
  
  # Create and save plot
  ggplot() +
    theme_bw() + theme(panel.grid = element_blank(), axis.title.x=element_text(size=10)) +
    geom_tile(data = df_smoothedDiff, aes(x=aaLoc, y=value, fill = val), na.rm = TRUE) +
    geom_point(data = selectVarVKGL, aes(x=aaLoc, y=value, colour=classification), alpha=1.0, size = 1, stroke = 1) +
    scale_fill_gradient2(name=paste("ClinVar\ndifference\nin 2D kernel\ndensity esti-\nmation of\n",gsub("\\.", "\n", termName),"\n(LP/P-LB/B)",sep=""), low = "#28A014", mid = "white", high = "#E41A1C") +
    scale_colour_manual(name = "VKGL\nclassification", values = c("LB/B" = "#28A014","VUS" = "#505050","LP/P" = "#E41A1C", "Conflicting" = "orange")) +
    geom_text(data = selectVarVKGL, aes(x=aaLoc, y=value, label=protChange), nudge_y=((ymax-ymin)/50), check_overlap = TRUE, alpha=1.0, size = 2) +
    scale_x_continuous(limits = c(xmin,xmax), labels = comma) +
    scale_y_continuous(limits = c(ymin,ymax), labels = comma) +
    xlab(paste("Amino acid position within protein alpha chain of", geneName, sep=" ")) +
    ylab(paste("Wild type vs mutant difference in",termName, sep=" ")) +
    ggtitle(label = paste("FoldX", termName, "results for VKGL variants in", geneName, "", sep=" "),
    subtitle = "(based on VKGL public release April 2023, ClinVar 20230702, FoldX 5.0, and AlphaFold2 human proteome v4)")
  ggsave(paste("dens_",geneName,"_",termName,".pdf",sep=""), width=9, height=5)
}
