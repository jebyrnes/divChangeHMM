################
# Attribution analysis
# with three predictors - temperature, wave height, help abundance
################

library(lmer)
library(lmerTest)
#library(MCMCglmm)

######
#Load and do some relabeling on the different datasets
######

#use fishPlot for analysis
source("fishPrep.R")

temp_summary_annual <- read.csv("../data/sbc_temp_summary_annual.csv")
temp_summary_annual <- temp_summary_annual %>% rename(Site=site)

sites <- read.csv("../data/LTER_Sites_latlong.csv")

######
#Combine datasets
######
#join the lat/long info with temps
tempSites <- join(sites, temp_summary_annual)

#use lat/long to join the temp info with fish
tempSites <- tempSites %>% 
  rename(Latitude = Lat, Longitude=Long, Year=interval_year)


fishWithTemp <- join(fishPlot, tempSites) %>% filter(!is.na(mean_temp_c))


######
# Modeling attribution
######

#Model the effect of temperature

tempLmer <- lmer(Aggregated_Richness ~ Year + mean_temp_c + 
                   (1|Year) + (1+mean_temp_c|Site),
                 data=fishWithTemp)
summary(tempLmer)
