#################################################################################
# Functions for analysis of SBC Fish Data for Change in Species Richness through time
# including data aggregation and determining region sizes
#
# Author: Jarrett Byrnes
#################################################################################

library(plyr)
library(dplyr)
library(sp)
library(rgeos)
library(rgdal)

getBoundingRegion <- function(dataset, 
                              inputProj=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"), 
                              proj=CRS("+init=ESRI:54009"), #projected lat/long
                              na.rm=T,
                              plotsize=0.00001){
  
  #first, filter down to unique latlongs
  uniqueLatLongs <- {dataset %>% group_by(Longitude, Latitude) %>% dplyr::summarise(n=length(Latitude))}[,1:2]
  
  #get rid of bad rows  
  if(na.rm) uniqueLatLongs <- uniqueLatLongs %>% filter(!is.na(Latitude) & !is.na(Longitude))

  #if we have <3 plots
  if(nrow(uniqueLatLongs)==1) return(1)
  if(nrow(uniqueLatLongs)==2) uniqueLatLongs <- rbind(uniqueLatLongs, uniqueLatLongs+c(plotsize,plotsize))
  
  #build a spatialPoints object
  spPoints <- SpatialPoints(coords=as.matrix(uniqueLatLongs[,1:2]), 
                            proj4string= inputProj)
  
  #reproject onto equal area projection
  spPoints <- spTransform(spPoints, proj)

  #get the convex hull
  spHull <- as(gConvexHull(spPoints), "SpatialPolygons")
  
  gArea(spHull)
}



#get unique permutations of a vector
#sampling from that vector without replacement
#to generate each permutation to preserve independence
unique_perm_sample <- function(v, l){
  m <- matrix(nrow=l)
  while(length(v)>l){
    s <- sample(v,l, replace=FALSE)
    v <- v[-which(v %in% s)]
    #print(s)
    m <- cbind(m,s)
  }
  m[,-1]
}

#
getSubData <- function(dataset, 
                       nplots, nsamps, 
                       sampleframe=NA, noSpecies="-99999",
                       uniquePerms=T){
  #if no sample frame with unique info about samples has been provided, make one
  if(is.na(sampleframe[1,1])){
    sampleframe <- dataset %>% 
      group_by(Latitude, Longitude) %>%
      summarise(nYearsSampled=n()) %>% 
      ungroup() 
    
    sampleframe$SampleID <- 1:nrow(sampleframe)
    
    dataset <- inner_join(dataset, sampleframe)
    
  }
  
  #what are the combinations of samples 
  #############consider permutations
  samps <- replicate(nplots, sample.int(nplots,nsamps))
  if(uniquePerms) samps <- unique_perm_sample(1:nplots,nsamps)
  if(nsamps==1) samps <- 1:nrow(dataset)
  if (is.null(dim(samps))) dim(samps) <- c(1, length(samps))
  
  newdata <- apply(samps, 2, function(x) makDatasetFromSampleIds(x, dataset, noSpecies=noSpecies))
  
  #turn the list back into a data frame
  newdata <- rbind_all(newdata)
  return(newdata)
}



makDatasetFromSampleIds <- function(x, dataset, noSpecies="-99999" , envt=F){
  
  #subset the data to one set of samples
  subdata <- dataset %>% filter(SampleID %in% x)
  
  #create a set of data for this group of plots 
  subdataOut <- subdata %>% group_by(Year) %>%
    summarise(Aggregated_Richness = length(unique(paste(Genus, Species))),
              Scale=length(unique(paste(Latitude, Longitude))),
              Bounded_region=getBoundingRegion(data.frame(Latitude=Latitude, Longitude=Longitude)),
              SampleID=paste(x, collapse="-"),
              noSp = as.numeric(noSpecies %in% unique(Species)))
  
  #deal with 0 species, so we don't have false 1 species samples
  noSp <- as.numeric(noSpecies %in% unique(subdata$Species))
  
  subdataOut$Aggregated_Richness <- subdataOut$Aggregated_Richness - noSp
  
  #return data frame
  subdataOut
  
}

