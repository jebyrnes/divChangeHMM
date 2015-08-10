#################################################
# Prepare raw data for biodiversity change database
# Compare alpha vs gamma scales
# Prepare database for CIEE database and papers
# Author: Robin Elahi
# Date: 150625
#################################################
# Dataset: Santa Barbara Coastal LTER
# www.sbc.lternet.edu
# Ongoing time-series of kelp forest community structure

# Abundance and size of reef fish
# Aggregate cryptic fish and fish on each transect
# Only use autumn data
# Lump species they had lumped in 2005

library(reshape2)
library(ggplot2)
theme_set(theme_bw(base_size=12))
library(dplyr)
options(dplyr.print_max = 1e9)
library(vegan)
library(nlme)
rm(list=ls(all=TRUE))

# Source functions
source("./R/divMetricF.R")
# plotting functions
source("./R/multiplotF.R")

# Load SBC data from lter website
dat <- read.csv("./data/all_fish_all_years_20140903.csv", 
                header=TRUE, na.strings= c("NA", "NULL"))
dim(dat)
# REMOVE OTRI - SP_CODE
dat <- filter(dat, SP_CODE != "OTRI")
# Select autumn data only
dat <- filter(dat, SURVEY_TIMING == "A")

# Dates	
dat$dateR <- as.Date(strftime(dat$date, format = "%Y-%m-%d"))
dat <- dat[with(dat, order(SP_CODE, site, transect, dateR)), ]

# write.csv(dat, 'dat.csv')

# remove data with count =  -99999, 5 observations
dat <- filter(dat, COUNT != -99999)

################################
# CREATE RAW DATA FILE 
################################
# The raw data have sizes, but I want counts
# sum abundances by site and date
# sum across quads and swaths, this removes the visibility column
# b/c visibility is only measured for fish (not cryptic fish)

dat2 <- aggregate(COUNT ~ year + month + date + dateR + 
                    SURVEY_TIMING  + 
                    site + transect + 
                    + TAXON_GENUS + TAXON_SPECIES + SP_CODE, 
                  data = dat, sum)

dat2 <- droplevels(dat2[with(dat2, order(SP_CODE, site, transect, dateR)), ])
summary(dat2)
# write.csv(dat2, 'dat2.csv')
################################
# Remove row IF genus is no fish,
dat3 <- filter(dat2, TAXON_GENUS != "No fish")
dim(dat3)

# create temp data frame to extract sites with no fish
tempDat <- filter(dat2, TAXON_GENUS == "No fish")
dim(tempDat) # only one transect had no fish 	
tempDat <- filter(tempDat, COUNT == 1)

################################
# If TAXON_GENUS = -99999, then remove row - funky data
# Must do this before adding the one observation of 'no fish'
dat4 <- dat3[dat3$TAXON_GENUS != "-99999", ]
dim(dat4)

# Add no fish back in
dat5 <- rbind(dat4, tempDat)			

styF <- function(x) return(as.factor(paste(x$site, x$transect, x$year, sep="_")))
dat5$sty <- styF(dat5)

# Check mohawk1 in 2009
mohk1 <- filter(dat5, site == "MOHK" & transect == "1" & year == "2009")
dim(mohk1)

mohk1 <- filter(dat5, sty == "MOHK_1_2009")
View(mohk1)

glimpse(mohk1)
mohk2 <- filter(mohk1, date == "2009-07-23")
mohk2 <- slice(mohk1, c(1:8, 12))
glimpse(mohk2)

mohk2 # good data for mohawk1 in 2009

# remove original mohawk 1 2009, then add back in the revised version
dat6 <- filter(dat5, sty != "MOHK_1_2009")
dim(dat5); dim(dat6)
dat7 <- rbind(dat6, mohk2)
dim(dat7)

################################
### Abundance
dat7$Abundance <- dat7$COUNT
dat7 <- dat7[with(dat7, order(site, transect, dateR)), ]

siteTranF <- function(x) return(as.factor(paste(x$site, x$transect, sep="_")))
dat7$siteTran <- siteTranF(dat7)
unique(dat7$siteTran) 

################################
# lump species that were lumped < 2005
genSpF <- function(x) return(as.factor(paste(x$TAXON_GENUS, x$TAXON_SPECIES, sep="_")))

dat7$genSp <- genSpF(dat7)

dat04 <- droplevels(filter(dat7, year < 2005))
dat05 <- droplevels(filter(dat7, year > 2004))
dim(dat05)

spp04 <- unique(dat04$genSp)
spp05 <- unique(dat05$genSp)

spp04 <- paste(c(spp04[order(spp04)], rep('null', 20)))
spp05 <- paste(spp05[order(spp05)])

tempDat <- data.frame(spp04, spp05)

# lump Sebastes flavidus, miniatus, paucispinis, rastrelliger, 
# and serriceps to spp
dat8 <- dat7
dat8$genSp <- gsub('Sebastes_flavidus', 'Sebastes_spp.', dat8$genSp)
dat8$genSp <- gsub('Sebastes_miniatus', 'Sebastes_spp.', dat8$genSp)
dat8$genSp <- gsub('Sebastes_paucispinis','Sebastes_spp.', dat8$genSp)
dat8$genSp <- gsub('Sebastes_rastrelliger', 'Sebastes_spp.', dat8$genSp)
dat8$genSp <- gsub('Sebastes_serriceps', 'Sebastes_spp.', dat8$genSp)
unique(dat8$genSp)
unique(dat7$genSp)

# now create new species column from new genSp column
# This is a bit of a hack
species2 <- unlist(strsplit(dat8$genSp, "_"))
length(species2)
# here i take the even values only
species3 <- species2[seq(2, length(species2), by = 2)] 

dat9 <- dat8
dat9$TAXON_SPECIES <- species3
unique(dat9$genSp)

head(dat9)

# Script ready for github

################################
# RAW DATASET IS READY, PREP FOR FISH ANALYSIS
################################
# go from long to wide
datL <- dat9

# rename columns to match
datL <- rename(datL, YEAR = year, SITE = site, TRANSECT = transect)
names(datL)

datW <- dcast(datL, YEAR + SITE + TRANSECT ~ genSp, 
              value.var = "Abundance")

names(datW)
summary(datW)
datW <- datW[with(datW, order(SITE, TRANSECT, YEAR)), ]

# remove NO FISH column
datW <- datW[, -33]
names(datW)

# get rich, div, even, abund
datW <- divMetF(datW, 4, 58)

# check for duplicates
styF <- function(x) return(as.factor(paste(x$SITE, x$TRANSECT, x$YEAR, sep="_")))
datW$sty <- styF(datW)
which(duplicated(datW$sty))

# now look for replication within sites across year
stF <- function(x) return(as.factor(paste(x$SITE, x$TRANSECT, sep="_")))
datW$st <- stF(datW)

# remove year 2000
datW <- filter(datW, YEAR > 2000)

fishDatW <- datW
write.csv(fishDatW, "fishDatW.csv")

# which site year combinations need to be removed so that I can calculate gamma div?
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

######################################################
# TRANSECT SCALE ANALYSIS
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

fm1 <- lme(fixed = rich ~ YEARZ, 
              data = lmeDat, method = "ML", 
              random =  list(rand3), 
              correlation = corAR1())

transect_YEARZ <- summary(fm1)$tTable[2,]

no_legend <- theme(legend.position = "none")

pTransects <- qplot(YEAR, rich, data = lmeDat, color = as.factor(st), 
                    geom = "line", xlab = "Year", ylab = "S") + 
  no_legend + ggtitle("Transect scale")

######################################################
# SITE ANALYSIS
alphaDF <- fish %>% rename(site = SITE, year = YEAR, transect = TRANSECT)

### Calculate sums per site, to get gamma scale metrics
names(alphaDF)
sdF <- function(x) return(as.factor(paste(x$site, x$year, sep="_")))
alphaDF$sd <- sdF(alphaDF)

# Create a uniform data frame
temp0 <- alphaDF # MODIFY
temp1 <- temp0[, 4:58] # MODIFY!
temp2 <- select(temp0, site, year, sd)
temp3 <- cbind(temp1, temp2)
head(temp3)
dim(temp3)

# Gamma - rich, div, even, abund
gammaDF <- temp3 %>% group_by(year, site, sd) %>% 
  summarise_each(funs(sum)) %>% ungroup()
head(gammaDF)
gammaDF <- divMetF(gammaDF, 4, dim(gammaDF)[2])
gammaDF$rich

# Regional - rich, div, even, abund
regDF <- temp3 %>% group_by(year) %>% 
  summarise_each(funs(sum)) %>% ungroup()
head(regDF)
regDF <- divMetF(regDF, 4, dim(regDF)[2])


names(alphaDF)
names(gammaDF)

##########################################
# Mixed model for gammaDF
summary(gammaDF$year)
gammaDF$yearZ <- gammaDF$year - 2007
summary(gammaDF)

fm2 <- lme(fixed = rich ~ yearZ, 
           data = gammaDF, method = "ML", 
           random =  list(~ 1 | site), 
           correlation = corAR1())
summary(fm2)$tTable

site_YEARZ <- summary(fm2)$tTable[2,]
##########################################

##########################################
# Linear model for regDF
summary(regDF$year)
regDF$yearZ <- regDF$year - 2007
summary(regDF)

# need a dummy variable 
regDF$dummy <- rep(1, dim(regDF)[1])

fm3 <- lme(fixed = rich ~ yearZ, 
           data = regDF, method = "ML", 
           random = ~ 1|dummy, 
           correlation = corAR1())
summary(fm3)$tTable[2,]

region_YEARZ <- summary(fm3)$tTable[2,]
##########################################
# bind the slope and se's 
scale_yearZ <- as.data.frame(rbind(transect_YEARZ, site_YEARZ, region_YEARZ))
scale_yearZ
scale_yearZ$scale <- c("Transect", "Site", "Regional")

# relevel
scale_yearZ$scale <- factor(scale_yearZ$scale, levels = c("Transect", "Site", "Regional"))

scale_trend <- ggplot(scale_yearZ, 
                      mapping=aes(x = scale, y = Value, 
                                  ymin = Value - 2*Std.Error, ymax = Value + 2*Std.Error)) +
  geom_point(size = 4) +
  geom_linerange(size = 1) +
  ylab("Temporal trend in richness\n(change in # of species/year)\n") + 
  xlab("Scale") +
  geom_hline(a = 0, linetype = 'dashed')

scale_trend
ggsave("./figs/scale_trend.pdf", width = 5, height = 5)

##########################################
pSites <- qplot(year, rich, data = gammaDF, color = site, 
                    geom = "line", xlab = "Year", ylab = "S") + 
  no_legend + ggtitle("Site scale")
pSites

pTransects <- qplot(year, rich, data = alphaDF, group = as.factor(transect), 
                    facets = ~site, xlab = "Year", ylab = "S", geom = "point") + 
  geom_line(aes(color = site)) + 
  no_legend + ggtitle("Transect scale")

pTransects <- ggplot(data = alphaDF, 
                     aes(year, rich, group = as.factor(transect))) +
  geom_line(aes(color = site)) + facet_wrap(~ site) + 
  no_legend + ggtitle("Transect scale") + 
  xlab("Year") + ylab("S") 

pTransects

pRegion <- qplot(year, rich, data = regDF, 
                    xlab = "Year", ylab = "S", geom = "line") + 
  no_legend + ggtitle("Regional scale")
pRegion


multiplot(pTransects, pSites, pRegion, layout = matrix(c(1,1,2,3), 
                                               nrow = 2, byrow = TRUE))

# save multiplot as 8x8 figure

multiplot(pTransects, pSites, pRegion, scale_trend,
          layout = matrix(c(1,1,2,3,4,4), nrow = 3, byrow = TRUE))


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


