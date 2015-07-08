################
# Attribution analysis
# with three predictors - temperature, wave height, help abundance
################

library(lme4)
library(lmerTest)
#library(MCMCglmm)

######
#Load and do some relabeling on the different datasets
######

#use fishPlot for analysis
source("fishPrep.R")

temp_summary_annual <- read.csv("../data/sbc_temp_summary_annual.csv")
temp_summary_annual <- temp_summary_annual %>% rename(Site=site)

waves <- read.csv("../data/sbc_wave_summary_annual.csv")
waves <- waves %>% rename(Nearest.MOP=mop)

sites <- read.csv("../data/LTER_Sites_latlong.csv")

######
#Combine datasets
######
#join the lat/long info with temps
tempSites <- join(sites, temp_summary_annual)

tempWavesSites <- join(tempSites, waves)


#use lat/long to join the temp info with fish
tempWavesSites <- tempWavesSites %>% 
  rename(Latitude = Lat, Longitude=Long, Year=interval_year)


fishWithPredictors <- join(fishPlot, tempWavesSites) %>% 
  filter(!is.na(mean_temp_c))%>% 
  filter(!is.na(mean_waveheight))


######
# Modeling attribution
######

#Model the effect of temperature

tempLmer <- lmer(Aggregated_Richness ~ Year + 
                   mean_temp_c + max_waveheight +
                   (1+max_waveheight+mean_temp_c|Site),
                 data=fishWithPredictors)
summary(tempLmer)
