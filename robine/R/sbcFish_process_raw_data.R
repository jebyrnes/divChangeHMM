######################################################
### Process kelp forest reef fish data 
### Author: Robin Elahi
### Date: 151202
######################################################


##### DETAILS #####
# Dataset: Santa Barbara Coastal LTER
# www.sbc.lternet.edu
# Ongoing time-series of kelp forest community structure
# Abundance and size of reef fish
# Aggregate cryptic fish and fish on each transect
# Only use autumn data
# Lump species they had lumped in 2005

##### LOAD PACKAGES, PREPARE RAW SIZE DATA #####
library(reshape2)
library(dplyr)
options(dplyr.print_max = 1e9)

rm(list=ls(all=TRUE))

# setwd
setwd("./robine/")

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

# remove data with count =  -99999, 5 observations
dat <- filter(dat, COUNT != -99999)

##### CREATE ABUNDANCE DATA FILE - LONG #####
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

####
# Remove row IF genus is no fish,
dat3 <- filter(dat2, TAXON_GENUS != "No fish")
dim(dat3)

# create temp data frame to extract sites with no fish
tempDat <- filter(dat2, TAXON_GENUS == "No fish")
dim(tempDat) # only one transect had no fish 	
tempDat <- filter(tempDat, COUNT == 1)

####
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

mohk2 <- filter(mohk1, date == "2009-07-23")
mohk2 <- slice(mohk1, c(1:8, 12))

mohk2 # good data for mohawk1 in 2009

# remove original mohawk 1 2009, then add back in the revised version
dat6 <- filter(dat5, sty != "MOHK_1_2009")
dim(dat5); dim(dat6)
dat7 <- rbind(dat6, mohk2)
dim(dat7)

### Abundance
dat7$Abundance <- dat7$COUNT
dat7 <- dat7[with(dat7, order(site, transect, dateR)), ]

siteTranF <- function(x) return(as.factor(paste(x$site, x$transect, sep="_")))
dat7$siteTran <- siteTranF(dat7)
unique(dat7$siteTran) 

#### lump species that were lumped < 2005
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

##### CHANGE FROM LONG TO WIDE FORMAT #####
### i.e., get a species abundance matrix

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
write.csv(fishDatW, "output/sbcFish_wide.csv")
