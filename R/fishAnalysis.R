#################################################################################
# Analysis of SBC Fish Data for Change in Species Richness through time
#
# Author: Jarrett Byrnes
#################################################################################

#########
#Load libraries and data
#########
source("./fishPrep.R")

#resulting in
#fish - the raw data
#fishSamples - the processed ready to go data
#fishPlot - fish richness over time at the plot level
#fishAll - fish richness over time with all plots aggregated

#########
# Create a simulated data set where for every
# scale, there are the same number of observations (total number of transects)
# And plot it relative to different predictors on the X-axis
#########

# We need to do half of 27 rounded down (so 13 samples). 14 plots cannot be 
# sampled without replacement (see the unique_perm_sample() function)
simData <- lapply(2:13, function(m) getSubData(dataset=fish, nplots=nplots, nsamps=m, sampleframe=fishSamples))
simData <- rbind_all(simData)
simData <- plyr::rbind.fill(simData, fishPlot)

qplot(Year, Aggregated_Richness, data=simData, group=SampleID, 
      geom="line", color=Bounded_region, facets=~Scale)


qplot(Scale, Aggregated_Richness, data=simData, group=SampleID, 
      geom="point", color=Bounded_region, facets=~Year)


qplot(Scale, Aggregated_Richness, data=simData, group=Year, 
      geom="point",size=I(0), color=Year) +
  stat_smooth(fill=NA)

qplot(Bounded_region, Aggregated_Richness, data=simData, group=SampleID, 
      geom="point", color=Scale, facets=~Year)


########
# Some analyses of simulated data
##########
library(lme4)
library(lmerTest)
library(MCMCglmm)

#Model it!

#Full model - will not coverge
mod.full <- lmer(Aggregated_Richness ~ Year*Scale*Bounded_region +
                   (1+ Year*Scale*Bounded_region|SampleID),
                 data=simData)

simData$BR <- simData$Bounded_region/sd(simData$Bounded_region)

#What about only the additive slopes varying?
#will not coverge
mod.simple.ranef <- lmer(Aggregated_Richness ~ Year*log(Scale)*BR +
                           (1+ Year+log(Scale)+BR|SampleID),
                         data=simData)

#Variable intercept only - converges!
#even with log(scale) varying randmoly!
mod.novarslope.ranef<- lmer(Aggregated_Richness ~ Year*log(Scale)*Bounded_region +
                                   (1 + log(Scale)|SampleID),
                                 data=simData)

summary(mod.novarslope.ranef)

#Let's look at autocorrelation structure
library(nlme)
mod.lme.full <- lme(Aggregated_Richness ~ Year*Scale*Bounded_region,
                    random = ~1+Year*Scale*Bounded_region|SampleID,
                    correlation= corCAR1(form = ~Year | SampleID),
                    data=simData)


mod.lme.varInt <- lme(Aggregated_Richness ~ Year*Scale*Bounded_region,
                      random = ~1|SampleID,
                      correlation= corCAR1(form = ~Year | SampleID),
                      data=simData) %>% select(distinct(SampleID, Year)))
summary(mod)

mod.mcmc <- MCMCglmm(Aggregated_Richness ~ Year*Scale*BR,
                     random=~us(1+Scale*BR):SampleID,
                     data=simData)

summary(mod.mcmc)

