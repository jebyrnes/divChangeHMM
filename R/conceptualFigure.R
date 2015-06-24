#########
#Load libraries and data
#########
source("./fishPrep.R")

#resulting in
#fish - the  data
#fishSamples - summary data
#fishPlot - fish richness over time at the plot level
#fishAll - fish richness over time with all plots aggregated

fishAll$SampleID <- "All"

##Create an intermediate data set
fishMed <- getSubData(fish, nplots, 6, 
                      sampleframe=fishSamples,
                      uniquePerms=T)

#qplot(Year, Aggregated_Richness, data=fishAll, geom="line") + 
#  xlab("Year") +ylab("Fish species Richness")

#qplot(Year, Aggregated_Richness, 
#      color=SampleID, data=fishMed, geom="line")

#bind them together
fishForFig <- plyr::rbind.fill(fishPlot, fishAll, fishMed)

fishForFig <- join(fishForFig, 
                   data.frame(Scale=c(1,6,nplots),
                              Aggregation=c("Small", "Medium", "Large")))

#make sure facets are in the right order
fishForFig$Aggregation <- factor(fishForFig$Aggregation, levels=c("Small", "Medium", "Large"))

#top panel
topPanel <- qplot(Year, Aggregated_Richness, 
      color=SampleID, data=fishForFig, geom="line",
      facets=~Aggregation) + 
  scale_color_discrete(guide="none") +
  ggtitle("Scale of Aggregation\n") +
  ylab("Richness")



#Analysis of faked data

simpleMod <- lmer(Aggregated_Richness ~ Year*log(Scale) +
                              (1+log(Scale)|SampleID),
                            data=fishForFig)

se.fixef <- function(obj){summary(obj)$coefficients[,2]}
simpleSE <- se.fixef(simpleMod)

coefDataFrame <- data.frame(Aggregation=c("Small", "Medium", "Large"),
                            Trend = c(fixef(simpleMod)[2]+fixef(simpleMod)[4],
                                      fixef(simpleMod)[2]+fixef(simpleMod)[4]*log(6),
                                      fixef(simpleMod)[2]+fixef(simpleMod)[4]*log(nplots)
                            ),
                            se=simpleSE[2]+simpleSE[4]
)
coefDataFrame$Aggregation <- factor(coefDataFrame$Aggregation, levels=c("Small", "Medium", "Large"))

bottomPanel <- ggplot(coefDataFrame, mapping=aes(x=Aggregation, y=Trend, 
                                  ymin=Trend-2*se, ymax=Trend+2*se)) +
  geom_point(size=5) +
  geom_linerange(size=1.2) +
  ylab("Trend in Diversity over Time\n(Change in # of species/year)\n")


library(gridExtra)
grid.arrange(topPanel, bottomPanel, ncol=1)
