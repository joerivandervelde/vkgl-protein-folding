
  
  ####
  # Protein based plot
  ####
  
  setwd(geneWorkingDir)
  

  
  # Loop over variables and create one plot for each
  for(termName in unique(mResults$variable))
  {
    
    # select term data and arrange for plot, LP/P on top
    # termName <- "Backbone.Hbond" # Debug purposes
    cat(paste("Plotting term", termName,"\n", sep=" "))
    selectVar <- mResults[mResults$variable==termName,]
    selectVar <- selectVar %>% arrange(factor(classification, levels = c("VUS","LB/B","LP/P")))
  
  # Determine amino acid position based plot window
  selectVar$aaLoc <- gsub("[A-Z]", "", selectVar$protChange)
  selectVar$aaLoc <- as.numeric(selectVar$aaLoc)
  xmin <- min(selectVar$aaLoc)
  xmax <- max(selectVar$aaLoc)
  ymin <- min(selectVar$value)
  ymax <- max(selectVar$value)
  
  selectVarLP <- subset(selectVar, classification == "LP/P")
  selectVarLB <- subset(selectVar, classification == "LB/B")
  
  xRes = 200
  yRes = 200
  
  everythingIsFine <- TRUE
  tryCatch({
    smoothLP <- MASS::kde2d(selectVarLP$aaLoc, selectVarLP$value, n=c(xRes, yRes),lims=c(range(selectVar$aaLoc), range(selectVar$value))) 
  }, error = function(e) {
    cat(paste("There was a problem for computing LP densities for term",termName,"\n",sep=" "))
    everythingIsFine <- FALSE
  })
  
  tryCatch({
    smoothLB <- MASS::kde2d(selectVarLB$aaLoc, selectVarLB$value, n=c(xRes, yRes),lims=c(range(selectVar$aaLoc), range(selectVar$value))) 
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
  
  # Create and save plot
  ggplot() +
    theme_bw() + theme(panel.grid = element_blank(), axis.title.x=element_text(size=10)) +
    geom_tile(data = df_smoothedDiff, aes(x=aaLoc, y=value, fill = val), na.rm = TRUE) +
    geom_point(data = selectVar, aes(x=aaLoc, y=value, colour=classification), alpha=1.0, size = 1, stroke = 1) +
    scale_fill_gradient2(name=paste("Difference\nin 2D kernel\ndensity esti-\nmation of\n",gsub("\\.", "\n", termName),"\n(LP/P-LB/B)",sep=""), low = "#28A014", mid = "white", high = "#E41A1C") +
    scale_colour_manual(name = "Classification", values = c("LB/B" = "#28A014","VUS" = "#505050","LP/P" = "#E41A1C")) +
    geom_text(data = selectVar, aes(x=aaLoc, y=value, label=protChange), nudge_y=((ymax-ymin)/50), check_overlap = TRUE, alpha=1.0, size = 2) +
    scale_x_continuous(limits = c(xmin,xmax), labels = comma) +
    xlab(paste("Amino acid position within protein alpha chain of", geneName, sep=" ")) +
    ylab(paste("Wild type vs mutant difference in",termName, sep=" ")) +
    ggtitle(paste("FoldX", termName, "results for clinically relevant variants in", geneName, "", sep=" "))
  ggsave(paste("aapos_",geneName,"_",termName,".pdf",sep=""), width=9, height=5)
  
  }
  
