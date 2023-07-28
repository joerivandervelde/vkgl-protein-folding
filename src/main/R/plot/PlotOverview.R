
setwd(geneWorkingDir)

# Plot 1: overview of all FoldX terms for LB/B and LP/P variants
mResultsOnlyLBLP <-subset(mResults, classification == "LB/B" | classification == "LP/P")

plotdata <- mResultsOnlyLBLP %>%
  group_by(classification, variable) %>%
  dplyr::summarize(n = n(),
                   mean = mean(value),
                   median = median(value),
                   sd = sd(value),
                   se = sd / sqrt(n))
plotdata$variable <- gsub("\\.", "\n", plotdata$variable)
ggplot(plotdata, 
       aes(x = classification, 
           y = mean,
           color = classification)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, 
                    ymax = mean + se),
                width = .1) +
  facet_grid(. ~ variable) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.text.x = element_text(size = 6),
        axis.text.x=element_text(size=7)) +
  labs(x="", 
       y="", 
       title=paste("FoldX terms for ", geneName, " based on VKGL and ClinVar variant classifications", sep=""),
       subtitle = "(Means and standard errors, based on VKGL April 2023, ClinVar 20230702, FoldX 5.0, and AlphaFold2 human proteome v4)") +
  scale_colour_manual(name = "Classification", values = c("LB/B" = "#28A014","LP/P" = "#E41A1C"))
ggsave(paste("overview_",geneName,".pdf",sep=""), width=9, height=5)