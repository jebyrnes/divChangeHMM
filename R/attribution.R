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
temp <- temp %>% dplyr::rename(Site=site)

waves <- read.csv("../data/sbc_wave_summary_MayToJuly.csv")
waves <- waves %>% dplyr::rename(Nearest.MOP=mop, Year=interval_year)

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
  dplyr::rename(Latitude = Lat, Longitude=Long)

tempWavesKelpSites <- join(tempWavesSites, kelp)

allSitePredictors <- join(tempWavesKelpSites, hard_substrate)

fishWithPredictors <- join(fishPlot, allSitePredictors %>% select(-X, -Month)) %>% 
  filter(!is.na(mean_temp_c))%>% 
  filter(!is.na(mean_waveheight))

################################
####Make permuted data frame
#todo: add a 'envt=T' option to getSubdata averaging envt vars
#       run model as in fishAnalysis but with envt predictors
#       look at how envt predictors change over time at different scales
################################
fishWithEnvt <- join(fish, allSitePredictors %>% select(-X, -Month))

simDataEnvt <- lapply(2:13, function(m) getSubData(dataset=fishWithEnvt, 
                                                   nplots=nplots, nsamps=m, 
                                                   sampleframe=fishSamples, envt=T))
simDataEnvt <- rbind_all(simDataEnvt)
simDataEnvt <- plyr::rbind.fill(simDataEnvt, fishWithPredictors)

######
# Modeling attribution for 1 plot
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
# Modeling attribution for sims
######

# Model the effect of temperature, waves, and kelp
# with a temp*wave interaction as abiotic influences might both affect one 
# another

fishLmerSims <- lmer(Aggregated_Richness ~Bounded_region* log(Scale)*(
                     scale(Year, scale=F) + 
                   scale(mean_temp_c, scale=F)* scale(mean_waveheight, scale=F) +
                   log_stipe_density + Hard_Substrate_Percent) +
                     (1+log(Scale)+ Bounded_region|SampleID),
                 data=simDataEnvt)
summary(fishLmerSims)

########
# Let's take a gander at kelp
#
# Although we have no expectations, as it's *winter*
# wave conditions that drive kelp abundance
########

kelpLmer <- lmer(log_stipe_density ~ scale(Year, scale=F) + 
                    mean_waveheight*
                   mean_temp_c  + 
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



