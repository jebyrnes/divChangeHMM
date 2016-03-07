#################################################################################
# Analysis of SBC Fish Data for Change in Species Richness through time
#
# Author: Jarrett Byrnes
# Contributor: Robin Elahi
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


##### FINAL ANALYSES USING LME, TO ACCOUNT FOR TEMPORAL AUTOCORRELATION ####
library(nlme)

# Center Year on minimum year
range(simData$Year)
sd(simData$Year)

simData$YearCentered <- scale(simData$Year, scale = FALSE)
sd(simData$YearCentered)

simData$Year0 <- simData$Year - min(simData$Year + 1)


# Scale
range(simData$Scale)
sd(simData$Scale)
simData$ScaleLog <- log(simData$Scale)
sd(simData$ScaleLog)

qplot(Scale, Aggregated_Richness, data=simData, group=Year, 
      geom="point",size=I(0), color=Year) +
  stat_smooth(fill=NA)

qplot(ScaleLog, Aggregated_Richness, data=simData, group=Year, 
      geom="point",size=I(0), color=Year) +
  stat_smooth(fill=NA)


# Bounded region
range(simData$Bounded_region)
sd(simData$Bounded_region)

# Log transform bounded region, and change to km2
simData$BrLog <- log(simData$Bounded_region)
range(simData$BrLog)
sd(simData$BrLog)

qplot(Bounded_region, Aggregated_Richness, data=simData, group=Year, 
      geom="point",size=I(0), color=Year) +
  stat_smooth(fill=NA)

qplot(BrLog, Aggregated_Richness, data=simData, group=Year, 
      geom="point",size=I(0), color=Year) +
  stat_smooth(fill=NA)

# Scale bounded_region
simData$BrScaled <- scale(simData$Bounded_region)
range(simData$BrScaled)
sd(simData$BrScaled)

# Final model
mod1 <- lme(Aggregated_Richness ~ Year0 * Scale + 
              Year0 * BrLog, 
            random = ~ Year0 | SampleID,
            correlation = corAR1(), 
            data = simData, 
            control = lmeControl(opt = "optim"))

summary(mod1)
summary(mod1)$tTable
plot(mod1)

write.csv(round(summary(mod1)$tTable, 3), 
          "../output/simData_tTable.csv")

##### FINAL PLOTS #####

theme_set(theme_bw(base_size = 12))


# All scales
ggplot(data = simData, aes(Year, Aggregated_Richness, group = SampleID, color = BrLog)) + 
  geom_line() + facet_wrap(~ Scale) +
  ylab("Richness") + 
  scale_color_gradient(name = "Bounded\narea\nlog(m2)")

ggsave("../figs/simPlot_all.pdf", height = 5, width = 5)

# Twelve transects only (leave out 13 for plotting)
ggplot(data = subset(simData, Scale != 13), aes(Year, Aggregated_Richness, 
                                                group = SampleID, color = BrLog)) + 
  geom_line() + facet_wrap(~ Scale) +
  ylab("Richness") + 
  scale_color_gradient(name = "Bounded\narea\nlog(m2)")

ggsave("../figs/simPlot_12.pdf", height = 5, width = 7)






##### OLD ANALYSES #####

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
library(broom)
mod.lme.full <- lme(Aggregated_Richness ~ Year*Scale*Bounded_region,
                    random = ~1+Year*Scale*Bounded_region|SampleID,
                    correlation= corCAR1(form = ~Year | SampleID),
                    data=simData)
plot(mod.lme.full)
tidy(mod.lme.full)

mod.lme.varInt <- lme(Aggregated_Richness ~ Year*Scale*Bounded_region,
                      random = ~1|SampleID,
                      correlation= corCAR1(form = ~Year | SampleID),
                      data=simData) 

summary(mod.lme.varInt)
plot(mod.lme.varInt)

mod.mcmc <- MCMCglmm(Aggregated_Richness ~ Year*Scale*BR,
                     random=~us(1+Scale*BR):SampleID,
                     data=simData)

summary(mod.mcmc)

