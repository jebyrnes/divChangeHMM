################
# Attribution analysis
# with three predictors - temperature, wave height, help abundance
################

library(lme4)
library(lmerTest)
library(broom)
library(nlme)

#library(MCMCglmm)

######
#Load and do some relabeling on the different datasets
######

#use fishPlot for analysis
source("fishPrep.R")
fishPlot %>% group_by(SampleID) %>% summarise(nYears = n_distinct(Year))

temp <- read.csv("../data/sbc_temp_summary_MayToJuly.csv")
temp <- temp %>% dplyr::rename(Site=site)
head(temp)
temp %>% group_by(Site) %>% summarise(nYears = n_distinct(Year))

waves <- read.csv("../data/sbc_wave_summary_MayToJuly.csv")
waves <- waves %>% dplyr::rename(Nearest.MOP=mop, Year=interval_year)
head(waves)

kelp <- read.csv("../data/sbc_macrocystis_annual.csv")
kelp$log_stipe_density <- log(kelp$stipe_density+1)
head(kelp)

hard_substrate <- read.csv("../data/sbc_hard_substrate.csv")

sites <- read.csv("../data/LTER_Sites_latlong.csv")
head(sites)

######
#Combine datasets
######
#join the lat/long info with temps
tempSites <- join(temp, sites)
tempSites %>% group_by(Transect) %>% summarise(nYears = n_distinct(Year))

tempWavesSites <- join(tempSites, waves)

#use lat/long to join the temp info with fish
tempWavesSites <- tempWavesSites %>% 
  dplyr::rename(Latitude = Lat, Longitude=Long)

tempWavesKelpSites <- join(tempWavesSites, kelp)

allSitePredictors <- join(tempWavesKelpSites, hard_substrate)

allSitePredictors %>% group_by(Transect) %>% summarise(nYears = n_distinct(Year))

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
                   mean_temp_c * mean_waveheight +
                   log_stipe_density + Hard_Substrate_Percent +
                   (1 |Site/Transect),
                 data=fishWithPredictors)
summary(fishLmer)

### CHECK THE VARIANCES
# How do the standardized inputs compare?
fishWithPredictors %>% summarise(rich_SD = sd(Aggregated_Richness), 
                                 Year_SD = sd(scale(Year, scale = T)), 
                                 tempC_SD = sd(scale(mean_temp_c)), 
                                 wave_SD = sd(scale(mean_waveheight)), 
                                 stipe_SD = sd(scale(log_stipe_density)), 
                                 hard_SD = sd(scale(Hard_Substrate_Percent, scale = T))
                                 )
summary(fishWithPredictors)


# Make sure data are ordered by year
fishWithPredictors <- fishWithPredictors %>% group_by(SampleID) %>% 
  arrange(Year)

fishWithPredictors$YearZ <- scale(fishWithPredictors$Year)
head(fishWithPredictors)

df1 <- fishWithPredictors %>%
  group_by(Site, Transect) %>%
  summarise(nYears = n_distinct(Year))

View(df1)

### Why is there a variable number of years per transect?  
### They should be equal (e.g., 14 years per transect)
### Ok, I see - driver data is not available for every year

fishLME <- lme(fixed = Aggregated_Richness ~ 1 + 
                 YearZ * (scale(mean_temp_c)) + 
                 YearZ * (scale(mean_waveheight)) + 
                 YearZ * scale(log_stipe_density) + 
                 YearZ * scale(Hard_Substrate_Percent), 
            data = fishWithPredictors, method = "REML", 
            random =  ~ 1 | Site/Transect, 
            correlation = corAR1())

summary(fishLME)
summary(fishLME)$tTable
plot(fishLME)

write.csv(round(summary(fishLME)$tTable, 3), 
          "../output/fishLME_tTable.csv")

qplot(Year, Aggregated_Richness, group=paste(Site, Transect),
      data=fishWithPredictors, size=log_stipe_density, alpha=I(0.7), 
      shape=factor(Transect)) +
  facet_wrap(~Site, scale="free_x") +
  stat_smooth(method="lm", fill=NA, color="red", size=0.5)

qplot(log_stipe_density, Aggregated_Richness, group=paste(Site, Transect),
      data=fishWithPredictors, alpha=I(0.7), 
      shape=factor(Transect)) +
  facet_wrap(~Site, scale="free_x") +
  stat_smooth(method="lm", fill=NA, color="red", size=0.5)


qplot(mean_temp_c, Aggregated_Richness, group=paste(Site, Transect),
      data=fishWithPredictors, size=mean_waveheight, alpha=I(0.7), 
      shape=factor(Transect)) +
  facet_wrap(~Site, scale="free_x") +
  stat_smooth(method="lm", fill=NA, color="red", size=0.5)

qplot(mean_waveheight, Aggregated_Richness, group=paste(Site, Transect),
      data=fishWithPredictors, size=mean_temp_c) +
  facet_wrap(~Site) +
  stat_smooth(method="lm", fill=NA, color="red", size=1)

# So: fish richness is predicted directly by kelp abundance, and the rate of change in fish richness is predicted by kelp abundance

# Next question: what predicts kelp abundance?
kelpLME <- lme(fixed = log_stipe_density ~ 1 + 
                 YearZ * scale(mean_temp_c) * 
                 YearZ * scale(mean_waveheight) + 
                 YearZ * scale(Hard_Substrate_Percent), 
               data = fishWithPredictors, method = "REML", 
               random =  ~ 1 | Site/Transect, 
               correlation = corAR1())

summary(kelpLME)
summary(kelpLME)$tTable
plot(kelpLME)

names(fishWithPredictors)
qplot(mean_waveheight, log_stipe_density, group=paste(Site, Transect),
      data=fishWithPredictors, size = mean_temp_c, alpha=I(0.7), 
      shape=factor(Transect)) +
  facet_wrap(~Site, scale="free_x") +
  stat_smooth(method="lm", fill=NA, color="red", size=0.5)

qplot(YearZ, log_stipe_density, group=paste(Site, Transect),
      data=fishWithPredictors, size = mean_waveheight, alpha=I(0.7), 
      shape=factor(Transect)) +
  facet_wrap(~Site, scale="free_x") +
  stat_smooth(method="lm", fill=NA, color="red", size=0.5)



########
# Modeling attribution for sims
######

# Model the effect of temperature, waves, and kelp
# with a temp*wave interaction as abiotic influences might both affect one 
# another

simDataEnvt %>% filter(Scale == "13") %>% summarise(nSims = n_distinct(SampleID))

summary(simDataEnvt)
fishLmerSims <- lmer(Aggregated_Richness ~ scale(Year, scale = F) + 
                       mean_temp_c * mean_waveheight +
                       log(Bounded_region) + Scale + 
                       log_stipe_density + Hard_Substrate_Percent +
                     (1 | SampleID),
                 data=simDataEnvt)
summary(fishLmerSims)

### RE's models

names(allSitePredictors)

# Summarize predictors for entire SBC region
allPredictors <- allSitePredictors %>% group_by(Year) %>% 
  summarise(temp_mean = mean(mean_temp_c, na.rm = TRUE), 
            temp_CV = (sd(mean_temp_c, na.rm = TRUE))/(mean(mean_temp_c, na.rm = TRUE)),
            stipe_mean = mean(log_stipe_density, na.rm = TRUE), 
            stipe_CV = (sd(log_stipe_density, na.rm = TRUE))/
              (mean(log_stipe_density, na.rm = TRUE)),
            hard_mean = mean(Hard_Substrate_Percent, na.rm = TRUE), 
            hard_CV = (sd(Hard_Substrate_Percent, na.rm = TRUE))/
              (mean(Hard_Substrate_Percent, na.rm = TRUE)),
            waves_mean = mean(mean_waveheight, na.rm = TRUE), 
            waves_CV = (sd(mean_waveheight, na.rm = TRUE))/
              (mean(mean_waveheight, na.rm = TRUE)))
            
            
allPredictors

# Join fishAll with predictors

fishAll_predictors <- join(fishAll, allPredictors) %>% 
  filter(!is.na(temp_mean)) %>% 
  filter(!is.na(temp_CV)) %>% 
  filter(!is.na(stipe_mean)) %>% 
  filter(!is.na(stipe_CV)) %>% 
  filter(!is.na(hard_mean)) %>% 
  filter(!is.na(hard_CV)) %>% 
  filter(!is.na(waves_mean)) %>% 
  filter(!is.na(waves_CV))

head(fishAll_predictors)  
names(fishAll_predictors)

fishAll_predictors$YearZ <- scale(fishAll_predictors$Year, scale = FALSE)

fishAll_gls1 <- gls(Aggregated_Richness ~ 1 + 
                     YearZ * scale(temp_mean) + 
                     YearZ * scale(stipe_mean) + 
                     YearZ * scale(hard_mean) + 
                     YearZ * scale(waves_mean), 
                   data = fishAll_predictors, 
                   correlation = corAR1())
summary(fishAll_gls1)$tTable

fishAll_gls2 <- gls(Aggregated_Richness ~ 1 + 
                     YearZ * scale(temp_CV) + 
                     YearZ * scale(stipe_CV) +
                     YearZ * scale(hard_CV) +
                     YearZ * scale(waves_CV), 
                   data = fishAll_predictors, 
                   correlation = corAR1())
summary(fishAll_gls2)$tTable

AIC(fishAll_gls1, fishAll_gls2) # CV has lower AIC


fishAll_gls3 <- gls(Aggregated_Richness ~ YearZ + 
                      YearZ : scale(temp_mean) + 
                      YearZ : scale(stipe_mean) + 
                      YearZ : scale(hard_mean) + 
                      YearZ : scale(waves_mean) + 
                      YearZ : scale(temp_CV) + 
                      YearZ : scale(stipe_CV) +
                      YearZ : scale(hard_CV) +
                      YearZ : scale(waves_CV), 
                    data = fishAll_predictors, 
                    correlation = corAR1())
summary(fishAll_gls3)
summary(fishAll_gls3)$tTable
write.csv(round(summary(fishAll_gls3)$tTable, 3), 
          "../output/fishALL_gls_tTable.csv")
plot(fishAll_gls2)


##### LOOK AT TRENDS IN CV #####
library(tidyr)
head(fishAll_predictors)

predLong <- fishAll_predictors %>% select(Year, temp_mean:waves_CV) %>%
  gather(key = driver, value = value, temp_mean:waves_CV)

ggplot(data = predLong, aes(Year, value)) + geom_point() + 
  geom_smooth(method = "lm") + facet_wrap(~ driver, scales = "free")



fishAll_lme <- lme(fixed = Aggregated_Richness ~ 1 + 
                     YearZ * scale(mean_temp_c) * 
                     YearZ * scale(mean_waveheight) + 
                     YearZ * scale(Hard_Substrate_Percent), 
                   data = fishWithPredictors, method = "REML", 
                   random =  ~ 1 | Site/Transect, 
                   correlation = corAR1())

summary(kelpLME)
summary(kelpLME)$tTable
plot(kelpLME)

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



