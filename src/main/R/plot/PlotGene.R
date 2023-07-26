
  ######
  # Plot 2: gene plot per attribute
  ######
  
  setwd(geneWorkingDir)
  

  
  # Loop over variables and create one plot for each
  for(termName in unique(mResults$variable))
  {
    
    # select term data and arrange for plot, LP/P on top
    # termName <- "total.energy" # Debug purposes
    selectVar <- mResults[mResults$variable==termName,]
    selectVar <- selectVar %>% arrange(factor(classification, levels = c("VUS","LB/B","LP/P")))
    
    # Determine optimal threshold using Youden's Index #
    cutpointDF <- subset(selectVar, classification != "VUS")
    opt_cut <- cutpointr(cutpointDF, value, classification, direction = ">=", pos_class = "LP/P", neg_class = "LB/B", method = maximize_metric, metric = youden)
    youdenIndex <- opt_cut$optimal_cutpoint
    tp <- sum(cutpointDF[cutpointDF$classification=="LP/P",'value'] >= youdenIndex)
    fp <- sum(cutpointDF[cutpointDF$classification=="LB/B",'value'] >= youdenIndex)
    ppv <- 100 *tp/(tp+fp)
    sens <- opt_cut$sensitivity*100
    
    # Determine plot window
    xmin <- min(exonsB37$exonStart)
    xmax <- max(exonsB37$exonEnd)
    ymin <- min(selectVar$value)
    ymax <- max(selectVar$value)
    
    # Create and save plot
    ggplot() +
      theme_bw() + theme(panel.grid = element_blank(), axis.title.x=element_text(size=10)) +
      geom_rect(data = exonsB37, aes(xmin = exonStart, xmax = exonEnd, ymin = ymin, ymax = ymax), linetype = 0, fill="lightgray", alpha = 1) +
      geom_point(data = selectVar, aes(x=pos, y=value, colour=classification), alpha=1.0, size = 1, stroke = 1) +
      geom_text(data = selectVar, aes(x=pos, y=value, label=protChange), nudge_y=((ymax-ymin)/50), check_overlap = TRUE, alpha=1.0, size = 2) +
      geom_hline(yintercept = youdenIndex) +
      scale_colour_manual(name = "Classification", values = c("LB/B" = "#28A014","VUS" = "#505050","LP/P" = "#E41A1C")) +
      scale_x_continuous(limits = c(xmin,xmax), labels = comma) +
      xlab(paste("",geneName," at GRCh37 chr",geneChrB37,":", xmin, "-",xmax,", lightgray: exons", sep="")) +
      ylab(termName) +
      ggtitle(paste("FoldX results for ",geneName,". At a threshold of ",round(youdenIndex, 2), " the PPV is ",round(ppv),"% and the sensitivity is ",round(sens),"%.",sep=""))
    ggsave(paste("genomicpos_",geneName,"_",termName,".pdf",sep=""), width=9, height=5)
    
  }
  
