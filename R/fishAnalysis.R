#################################################################################
# Analysis of SBC Fish Data for Change in Species Richness through time
#
# Author: Jarrett Byrnes
#################################################################################

#########
#Load libraries and data
#########

#for data aggregation
library(dplyr)

#for plotting
library(ggplot2)
theme_set(theme_bw(base_size=17))

#For analysis
library(lme4)
library(lmerTest)

source("./dataGenerationFunctions.R")

#Load the fish data
#fish <- read.csv("../../output_data/sbc_fish_species_rawdata.csv") #jarrett's data for checking
fish <- read.csv("../data/rawData_SBCfish.csv")
fish <- subset(fish, fish$Year!=2000) #too much variation in this year

#########
# So, there are differences in # of years sampled from site to site
# and we want the core set that was sampled in all years
#filter to sites that are in the whole dataset
#########

fishSamples <- fish %>% group_by(Latitude, Longitude, Year) %>%
  summarise(n=n()) %>% ungroup() %>% 
  group_by(Latitude, Longitude) %>%
  summarise(nYearsSampled=n()) %>% 
  ungroup() %>%
  filter(nYearsSampled==max(nYearsSampled))

fishSamples$SampleID <- 1:nrow(fishSamples)
  
fish <- inner_join(fish, fishSamples)

nplots <- length(unique(paste(fish$Latitude, fish$Longitude)))


#########
#What fish showed up after the 2005 bump in diversity?
#########

firstDate2005 <- fish %>% group_by(Genus, Species) %>%
  summarise(firstDate=min(Year)) %>%
  filter(firstDate>2004)


##########
#Look at small and large scales only
##########
fishPlot <- fish %>% group_by(Latitude, Longitude, Year, SampleID) %>%
  summarise(Aggregated_Richness = length(unique(paste(Genus, Species))),
            Scale=1,
            Bounded_region=1)
qplot(Year, Aggregated_Richness, color=paste(Latitude, Longitude), data=fishPlot, geom="line")


fishAll <- fish %>% group_by(Year) %>%
  summarise(Aggregated_Richness = length(unique(paste(Genus, Species))),
            Scale=length(unique(paste(Latitude, Longitude))),
            Bounded_region=getBoundingRegion(data.frame(Latitude=Latitude, Longitude=Longitude)))

qplot(Year, Aggregated_Richness, data=fishAll, geom="line") + xlab("Year") +ylab("Fish species Richness")


#########
# Create a simulated data set where for every
# scale, there are the same number of observations (total number of transects)
# And plot it relative to different predictors on the X-axis
#########

simData <- lapply(2:17, function(m) getSubData(fish, nplots, m, sampleframe=fishSamples))
simData <- plyr::ldply(simData)
simData <- plyr::rbind.fill(simData, fishPlot)

qplot(Year, Aggregated_Richness, data=simData, group=SampleID, 
      geom="line", color=Bounded_region, facets=~Scale)


qplot(Scale, Aggregated_Richness, data=simData, group=SampleID, 
      geom="point", color=Bounded_region, facets=~Year)


qplot(Scale, Aggregated_Richness, data=simData, group=Year, 
      geom="point",size=I(0), color=Year) +
  stat_smooth(fill=NA)

qplot(Bounded_region, Aggregated_Richness, data=simData, group=SampleID, 
      geom="point", color=Scale, facets=~Year)


########
# Some analyses of simulated data
##########
library(lme4)
library(lmerTest)
library(MCMCglmm)

#Model it!

#Full model - will not coverge
mod.full <- lmer(Aggregated_Richness ~ Year*Scale*Bounded_region +
                   (1+ Year*Scale*Bounded_region|sampleID),
                 data=simData)

simData$BR <- simData$Bounded_region/sd(simData$Bounded_region)

#What about only the additive slopes varying?
#will not coverge
mod.simple.ranef <- lmer(Aggregated_Richness ~ Year*log(Scale)*BR +
                           (1+ Year+log(Scale)+BR|SampleID),
                         data=simData)

#Variable intercept only - converges!
#even with log(scale) varying randmoly!
mod.novarslope.ranef<- lmer(Aggregated_Richness ~ Year*log(Scale)*Bounded_region +
                                   (1+log(Scale)|SampleID),
                                 data=simData)

summary(mod.novarslope.ranef)

#Let's look at autocorrelation structure
library(nlme)
mod.lme.full <- lme(Aggregated_Richness ~ Year*Scale*Bounded_region,
                    random = ~1+Year*Scale*Bounded_region|SampleID,
                    correlation= corCAR1 (form = ~Year | SampleID),
                    data=simData)


mod.lme.varInt <- lme(Aggregated_Richness ~ Year*Scale*Bounded_region,
                      random = ~1|sampleID,
                      correlation= corCAR1 (form = ~Year | sampleID),
                      data=simData)
summary(mod)

mod.mcmc <- MCMCglmm(Aggregated_Richness ~ Year*Scale*BR,
                     random=~us(1+Scale*BR):sampleID,
                     data=simData)

summary(mod.mcmc)

