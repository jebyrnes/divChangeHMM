################
# Attribution analysis, at the regional scale
# Robin Elahi, Jarrett Byrnes
# 7 March 2016
################

rm(list=ls(all=TRUE)) # removes all previous material from R's memory

# We did not detect a temporal trend at the local scale, so we cannot attribute anything
# But, there is a significant positive trend at the regional scale

library(nlme)

##### JOIN FISH AND DRIVER DATA #####

#use fishPlot for analysis
setwd("./R/")
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

##### GLS MODELS, ACCOUNTING FOR TEMPORAL AUTOCORRELATION #####
fishAll_predictors <- droplevels(fishAll_predictors)
fishAll_predictors$Year0 <- fishAll_predictors$Year - min(fishAll_predictors$Year - 1)
head(fishAll_predictors)

summary(gls(Aggregated_Richness ~ Year0, 
            data = fishAll_predictors, 
            correlation = corAR1()))

# But with fish_all?  
summary(gls(Aggregated_Richness ~ Year, 
            data = fishAll, 
            correlation = corAR1()))


#### NO SIGNIFICANT TEMPORAL TREND ---- STOP HERE, NO ATTRIBUTION? #####
# 
qplot(Year, Aggregated_Richness, data = fishAll)
qplot(Year, Aggregated_Richness, data = fishAll_predictors)

with(fishAll_predictors, sd(temp_mean))
with(fishAll_predictors, sd(stipe_mean))
with(fishAll_predictors, sd(waves_mean))

fishAll_gls1 <- gls(Aggregated_Richness ~ Year0 * stipe_mean, 
                    data = fishAll_predictors, 
                    correlation = corAR1())
summary(fishAll_gls1)
summary(fishAll_gls1)$tTable
plot(fishAll_gls1)

write.csv(round(summary(fishAll_gls1)$tTable, 3), 
          "../output/fishALL_gls_tTable.csv")
