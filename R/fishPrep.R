#################################################################################
# Prep for analysis of SBC Fish Data for Change in Species Richness through time
#
# Author: Jarrett Byrnes
#################################################################################

#########
#Load libraries and data
#########

#for plotting
library(plyr)
library(ggplot2)
theme_set(theme_bw(base_size=17))

#for data aggregation
library(dplyr)

#For analysis
library(lme4)
library(lmerTest)

#Load methods for dealing with data
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
  summarise(n=length(Latitude)) %>% ungroup() %>%  #reduce the dataset down to one line per site*year
  group_by(Latitude, Longitude) %>%
  summarise(nYearsSampled=length(Latitude)) %>% 
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
#qplot(Year, Aggregated_Richness, color=paste(Latitude, Longitude), data=fishPlot, geom="line")


fishAll <- fish %>% group_by(Year) %>%
  summarise(Aggregated_Richness = length(unique(paste(Genus, Species))),
            Scale=length(unique(paste(Latitude, Longitude))),
            Bounded_region=getBoundingRegion(data.frame(Latitude=Latitude, Longitude=Longitude)))

#qplot(Year, Aggregated_Richness, data=fishAll, geom="line") + xlab("Year") +ylab("Fish species Richness")
