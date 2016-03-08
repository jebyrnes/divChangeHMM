######################################################
### Analyze biodiversity trends at different spatial 
### scales using SBC LTER reef fish data
### Author: Robin Elahi
### Date: 151202

### Modified
### 160308 - added rarefied S
######################################################

##### DETAILS #####
# The original fish data were pre-processed
# See sbcFish_process_raw_data.R for details

##### LOAD PACKAGES, DATA #####
library(ggplot2)
library(dplyr)
options(dplyr.print_max = 1e9)
library(vegan)
library(nlme)

rm(list=ls(all=TRUE)) 

# setwd
setwd("./robine/")

# Source functions
source("./R/divMetricF.R")

# plotting functions
source("./R/multiplotF.R")

# Load pre-processed sbc fish species-abundance matrix
datW <- read.csv("./output/sbcFish_wide.csv", 
                header=TRUE, na.strings= c("NA", "NULL"))
names(datW)

# get div metrics
datW <- divMetF(datW, 5, 59)

##### MASSAGE RAW DATA FOR FISH ANALYSIS #####

# Include transects with data from 2001-2014
fishSamples <- datW %>% group_by(SITE, TRANSECT, YEAR) %>%
  summarise(n = n()) %>% ungroup() %>%  #reduce the dataset down to one line per site*year
  group_by(SITE, TRANSECT) %>%
  summarise(nYearsSampled = n()) %>% 
  ungroup() %>%
  filter(nYearsSampled==max(nYearsSampled))

fishSamples
head(fishSamples)
fishSamples$SampleID <- 1:nrow(fishSamples)

names(datW)
fish <- inner_join(datW, fishSamples)
dim(fish)
nplots <- length(unique(paste(fish$SITE, fish$TRANSECT)))
nplots

qplot(YEAR, rich, data = fish, color = as.factor(TRANSECT), 
      facets = ~SITE, geom = "line")

# remove sites with only one transect
fish <- filter(fish, SITE != "AHND" & 
                  SITE != "GOLB" &
                  SITE != "BULL")

qplot(YEAR, rich, data = fish, color = as.factor(TRANSECT), 
      facets = ~SITE, geom = "line")

qplot(YEAR, chao, data = fish, color = as.factor(TRANSECT), 
      facets = ~SITE, geom = "line")

head(fish)

ind_df <- fish %>% group_by(YEAR) %>% summarise(minInd = min(abund), maxInd = max(abund))

fish2 <- left_join(fish, ind_df, by = "YEAR")

##### GET SITE (GAMMA) AND REGIONAL ESTIMATES OF BIODIVERSITY #####

alphaDF <- fish2 %>% rename(site = SITE, year = YEAR, transect = TRANSECT)
alphaDF <- droplevels(alphaDF)

### Calculate sums per site, to get gamma scale metrics
names(alphaDF)
sdF <- function(x) return(as.factor(paste(x$site, x$year, sep="_")))
alphaDF$sd <- sdF(alphaDF)
names(alphaDF)

# Create a uniform data frame
temp0 <- alphaDF # MODIFY
temp1 <- temp0[, 5:59] # MODIFY!
temp2 <- select(temp0, site, year, sd)
temp3 <- cbind(temp1, temp2)
head(temp3)
dim(temp3)

# Gamma - rich, div, even, abund
gammaDF <- temp3 %>% group_by(year, site, sd) %>% 
  summarise_each(funs(sum)) %>% ungroup()
names(gammaDF)
gammaDF <- divMetF(gammaDF, 4, dim(gammaDF)[2])
gammaDF$abund

# Regional - rich, div, even, abund
str(temp3)
regDF <- temp3 %>% select(-site, -sd) %>% group_by(year) %>% 
  summarise_each(funs(sum)) %>% ungroup()
names(regDF)
regDF <- divMetF(regDF, 2, dim(regDF)[2])

names(alphaDF)
names(gammaDF)

##### TRANSECT SCALE ANALYSIS #####
lmeDat <- fish
# center year on 2007
unique(lmeDat$YEAR)
lmeDat$YEARZ <- lmeDat$YEAR - 2007
unique(lmeDat$YEARZ)
# using the centered year doesn't improve the 
# correlation between slope and int
str(lmeDat)

rand1 <- ~ 1 | SITE
rand2 <- ~ YEARZ | SITE
rand3 <- ~ YEARZ | TRANSECT

fm1 <- lme(fixed = chao ~ YEARZ, 
           data = lmeDat, method = "ML", 
           random =  list(rand3), 
           correlation = corAR1())

transect_YEARZ <- summary(fm1)$tTable[2,]

##### GAMMA SCALE ANALYSIS #####
summary(gammaDF$year)
gammaDF$yearZ <- gammaDF$year - 2007
summary(gammaDF)

fm2 <- lme(fixed = chao ~ yearZ, 
           data = gammaDF, method = "ML", 
           random =  list(~ 1 | site), 
           correlation = corAR1())
summary(fm2)$tTable
site_YEARZ <- summary(fm2)$tTable[2,]

##### REGIONAL SCALE ANALYSIS #####
summary(regDF$year)
regDF$yearZ <- regDF$year - 2007
summary(regDF)

# need a dummy variable 
regDF$dummy <- rep(1, dim(regDF)[1])

fm3 <- lme(fixed = chao ~ yearZ, 
           data = regDF, method = "ML", 
           random = ~ 1|dummy, 
           correlation = corAR1())
summary(fm3)$tTable[2,]
region_YEARZ <- summary(fm3)$tTable[2,]

##### PLOT MODEL COEFFICIENTS #####
# bind the slope and se's from the linear models
scale_yearZ <- as.data.frame(rbind(transect_YEARZ, site_YEARZ, region_YEARZ))
scale_yearZ
scale_yearZ$scale <- c("Transect", "Site", "Regional")

# relevel
scale_yearZ$scale <- factor(scale_yearZ$scale, 
                            levels = c("Transect", "Site", "Regional"))

scale_trend <- ggplot(scale_yearZ, 
                      mapping=aes(x = scale, y = Value, 
                                  ymin = Value - 2*Std.Error, 
                                  ymax = Value + 2*Std.Error)) +
  geom_point(size = 4) +
  geom_linerange(size = 1) +
  ylab("Temporal trend in richness\n(change in # of species/year)\n") + 
  xlab("Scale") +
  geom_hline(yintercept = 0, linetype = 'dashed')

scale_trend

##### PLOT TEMPORAL TRENDS FOR EACH SCALE #####
no_legend <- theme(legend.position = "none")
theme_set(theme_bw(base_size=16))

pTransects <- ggplot(data = alphaDF, 
                     aes(year, chao, group = as.factor(transect))) +
  geom_line(aes(color = site)) + facet_wrap(~ site) + 
  no_legend + ggtitle("Transect scale") + 
  xlab("Year") + ylab("S-chao") 
pTransects

pSites <- qplot(year, chao, data = gammaDF, color = site, 
                geom = "line", xlab = "Year", ylab = "S-chao") + 
  no_legend + ggtitle("Site scale")
pSitesrobine

pRegion <- qplot(year, chao, data = regDF, 
                 xlab = "Year", ylab = "S-chao", geom = "line") + 
  no_legend + ggtitle("Regional scale")
pRegion

### save as pdf
pdf("./figs/sbcFish_chaoTrend.pdf", width = 7, height = 10)

multiplot(pTransects, pSites, pRegion, scale_trend,
          layout = matrix(c(1,1,2,3,4,4), nrow = 3, byrow = TRUE))

dev.off()	

##### PLOT TEMPORAL TRENDS FOR EACH SCALE - ABUNDANCE #####
no_legend <- theme(legend.position = "none")
theme_set(theme_bw(base_size=16))

pTransects <- ggplot(data = alphaDF, 
                     aes(year, abund, group = as.factor(transect))) +
  geom_line(aes(color = site)) + facet_wrap(~ site) + 
  no_legend + ggtitle("Transect scale") + 
  xlab("Year") + ylab("Abundance")
pTransects

pSites <- qplot(year, abund, data = gammaDF, color = site, 
                geom = "line", xlab = "Year", ylab = "Abundance") + 
  no_legend + ggtitle("Site scale")
pSites

pRegion <- qplot(year, abund, data = regDF, 
                 xlab = "Year", ylab = "Abundance", geom = "line") + 
  no_legend + ggtitle("Regional scale")
pRegion

### save as pdf
pdf("./figs/sbcFish_abundTrend.pdf", width = 7, height = 7)

multiplot(pTransects, pSites, pRegion,
          layout = matrix(c(1,1,2,3), nrow = 2, byrow = TRUE))

dev.off()	


########
########
########
########
summary(alphaDF$div)
# DIVERSITY

pSites <- qplot(year, div, data = gammaDF, color = site, 
                geom = "line", xlab = "Year", ylab = "H'") + 
  no_legend + ggtitle("Site scale")
pSites

pTransects <- qplot(year, div, data = alphaDF, group = as.factor(transect), 
                    facets = ~site, xlab = "Year", ylab = "H'", geom = "point") + 
  geom_line(aes(color = site)) + 
  no_legend + ggtitle("Transect scale")

pTransects <- ggplot(data = alphaDF, 
                     aes(year, div, group = as.factor(transect))) +
  geom_line(aes(color = site)) + facet_wrap(~ site) + 
  no_legend + ggtitle("Transect scale") + 
  xlab("Year") + ylab("H'") 

pTransects

pRegion <- qplot(year, div, data = regDF, 
                 xlab = "Year", ylab = "H'", geom = "line") + 
  no_legend + ggtitle("Regional scale")
pRegion


multiplot(pTransects, pSites, pRegion, layout = matrix(c(1,1,2,3), 
                                                       nrow = 2, byrow = TRUE))

# save multiplot as 8x8 figure

########
summary(alphaDF$even)
# Evenness

pSites <- qplot(year, even, data = gammaDF, color = site, 
                geom = "line", xlab = "Year", ylab = "J") + 
  no_legend + ggtitle("Site scale")
pSites

pTransects <- qplot(year, even, data = alphaDF, group = as.factor(transect), 
                    facets = ~site, xlab = "Year", ylab = "J", geom = "point") + 
  geom_line(aes(color = site)) + 
  no_legend + ggtitle("Transect scale")

pTransects <- ggplot(data = alphaDF, 
                     aes(year, even, group = as.factor(transect))) +
  geom_line(aes(color = site)) + facet_wrap(~ site) + 
  no_legend + ggtitle("Transect scale") + 
  xlab("Year") + ylab("J") 

pTransects

pRegion <- qplot(year, even, data = regDF, 
                 xlab = "Year", ylab = "J", geom = "line") + 
  no_legend + ggtitle("Regional scale")
pRegion


multiplot(pTransects, pSites, pRegion, layout = matrix(c(1,1,2,3), 
                                                       nrow = 2, byrow = TRUE))

# save multiplot as 8x8 figure






#### OLD CODE ####
# which site year combinations need to be removed 
# so that I can calculate gamma div?
with(datW, table(st, YEAR))
head(datW)

unique(datW$st)
datW2 <- droplevels(filter(datW, 
                           st != "AQUE_6" &
                             st != "BULL_2" &
                             st != "BULL_4" &
                             st != "BULL_5" &
                             st != "BULL_7" &
                             st != "BULL_8" &
                             st != "CARP_2" &
                             st != "CARP_4" &
                             st != "CARP_8" &
                             st != "IVEE_3" &
                             st != "IVEE_4" &
                             st != "IVEE_5" &
                             st != "IVEE_6" &
                             st != "IVEE_7" &
                             st != "IVEE_8" &
                             st != "SCDI_1"))

with(datW2, table(st, YEAR))

# select good ahnd, golb, bull
ahnd <- droplevels(filter(datW2, SITE == "AHND" &
                            YEAR > 2006))
golb <- droplevels(filter(datW2, SITE == "GOLB" &
                            YEAR > 2005))
bull <- droplevels(filter(datW2, SITE == "BULL" &
                            YEAR < 2011))

# remove ahdn, golb, bull
datW3 <- filter(datW2, SITE != "AHND" & 
                  SITE != "GOLB" &
                  SITE != "BULL")

datW4 <- rbind(datW3, ahnd, golb, bull)
dim(with(datW4, table(st, YEAR)))

with(datW4, table(st, YEAR))

# quick plots
names(datW4)
qplot(YEAR, rich, data = datW4, color = as.factor(TRANSECT), 
      facets = ~SITE, geom = "line")

# LME model at the transect scale
lmeDat <- filter(datW4, SITE != "AHND" & SITE != "GOLB")
lmeDat <- datW4

# Include transects with data from 2001-2014
fishSamples <- datW4 %>% group_by(SITE, TRANSECT, YEAR) %>%
  summarise(n = n()) %>% ungroup() %>%  #reduce the dataset down to one line per site*year
  group_by(SITE, TRANSECT) %>%
  summarise(nYearsSampled = n()) %>% 
  ungroup() %>%
  filter(nYearsSampled==max(nYearsSampled))

head(fishSamples)

fishSamples$SampleID <- 1:nrow(fishSamples)

names(datW4)
fish <- inner_join(datW4, fishSamples)
dim(fish)
nplots <- length(unique(paste(fish$SITE, fish$TRANSECT)))
nplots

qplot(YEAR, rich, data = fish, color = as.factor(TRANSECT), 
      facets = ~SITE, geom = "line")

