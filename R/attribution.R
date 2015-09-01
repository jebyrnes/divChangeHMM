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
kelp$log_stipe_density <- log(kelp$stipe_density+1)

hard_substrate <- read.csv("../data/sbc_hard_substrate.csv")

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

fishWithPredictors <- join(fishWithPredictors, hard_substrate)
######
# Modeling attribution
######

# Model the effect of temperature, waves, and kelp
# with a temp*wave interaction as abiotic influences might both affect one 
# another

fishLmer <- lmer(Aggregated_Richness ~ scale(Year, scale=F) + 
                   mean_temp_c* mean_waveheight +
                   log_stipe_density + Hard_Substrate_Percent +
                   (1 |Site/Transect),
                 data=fishWithPredictors)
summary(fishLmer)

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
# wave conditions that drive kelp abundance
########

kelpLmer <- lmer(log_stipe_density ~ scale(Year, scale=F) + 
                    mean_waveheight*
                   I(mean_temp_c^2)  + 
                   Hard_Substrate_Percent+
                   (1 |Site/Transect) ,
                 data=fishWithPredictors)
summary(kelpLmer)


########
# Let's take a gander at temperature
#
# Although we have no expectations, as it's *winter*
# wave conditions that drive kelp abundance
########

tempLmer <- lmer(mean_temp_c ~ scale(Year, scale=F) + 
                   (1 |Site/Transect) ,
                 data=fishWithPredictors)
summary(tempLmer)




########
# Let's take a gander at waves
#
# Although we have no expectations, as it's *winter*
# wave conditions that drive kelp abundance
########

waveLmer <- lmer(mean_waveheight ~ scale(Year, scale=F) + 
                   (1 |Site/Transect) ,
                 data=fishWithPredictors)
summary(waveLmer)


