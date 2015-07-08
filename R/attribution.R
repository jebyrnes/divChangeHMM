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

temp <- read.csv("../data/sbc_temp_summary_MayToJuly.csv")
temp <- temp %>% rename(Site=site)

waves <- read.csv("../data/sbc_wave_summary_MayToJuly.csv")
waves <- waves %>% rename(Nearest.MOP=mop, Year=interval_year)

kelp <- read.csv("../data/sbc_macrocystis_annual.csv")

sites <- read.csv("../data/LTER_Sites_latlong.csv")

######
#Combine datasets
######
#join the lat/long info with temps
tempSites <- join(temp, sites)

tempWavesSites <- join(tempSites, waves)


#use lat/long to join the temp info with fish
tempWavesSites <- tempWavesSites %>% 
  rename(Latitude = Lat, Longitude=Long)

tempWavesKelpSites <- join(tempWavesSites, kelp)

fishWithPredictors <- join(fishPlot, tempWavesKelpSites) %>% 
  filter(!is.na(mean_temp_c))%>% 
  filter(!is.na(mean_waveheight))


######
# Modeling attribution
######

# Model the effect of temperature, waves, and kelp
# with a temp*wave interaction as abiotic influences might both affect one 
# another

tempLmer <- lmer(Aggregated_Richness ~ scale(Year, scale=F) + 
                   mean_temp_c * mean_waveheight +
                   stipe_density +
                   (1 |Site/Transect) ,
                 data=fishWithPredictors)
summary(tempLmer)

qplot(mean_temp_c, Aggregated_Richness, group=paste(Site, Transect),
      data=fishWithPredictors, size=mean_waveheight, alpha=I(0.7), 
      shape=factor(Transect)) +
  facet_wrap(~Site, scale="free_x") +
  stat_smooth(method="lm", fill=NA, color="red", size=0.5)


qplot(mean_waveheight, Aggregated_Richness, group=paste(Site, Transect),
      data=fishWithPredictors, size=mean_temp_c) +
  facet_wrap(~Site) +
  stat_smooth(method="lm", fill=NA, color="red", size=1)


########
# Let's take a gander at kelp
#
# Although we have no expectations, as it's *winter*
# conditions that drive kelp abundance
########

kelpLmer <- lmer(stipe_density ~ scale(Year, scale=F) + 
                   mean_temp_c * mean_waveheight +
                   (1 |Site/Transect) ,
                 data=fishWithPredictors)
summary(kelpLmer)

with(fishWithPredictors, cor(cbind(mean_temp_c, mean_waveheight, stipe_density)))
