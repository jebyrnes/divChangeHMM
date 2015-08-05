#################################################
# Analyze SBC lter temperature data
# Author: Robin Elahi
# Date: 150626

# Date: 150805
# Include six sites for which at least two transects were sampled
# from 2002-2015
#################################################

library(ggplot2)
theme_set(theme_bw(base_size = 12))
library(dplyr)
options(dplyr.print_max = 1e9)

# Load daily temperature data
dat<- read.csv("../data/sbc_temp_summary_daily.csv")
names(dat)

# get dates in R format
dat$dateR <- as.Date(with(dat, paste(Year, Month, Day, sep = "-")))
glimpse(dat)

# quick plot
qplot(dateR, mean_temp_c, data = dat, facets = ~site, geom = "line", 
      xlab = "Year", ylab = "Temperature (C)")
ggsave("../figs/dailyTemp_bySite.pdf")

# now subset 6 sites for which we have good diversity data
datSub <- dat %>% filter(site == "ABUR" |
                           site == "AQUE" |
                           site == "CARP" |
                           site == "IVEE" |
                           site == "MOHK" |
                           site == "NAPL")
qplot(dateR, mean_temp_c, data = datSub, facets = ~site, geom = "line", 
      xlab = "Year", ylab = "Temperature (C)", color = site)
ggsave("../figs/dailyTemp_bySiteSub.pdf")

# use monthly data for better visualisation?
dat2 <- read.csv("../data/sbc_temp_summary_monthly.csv")
glimpse(dat2)

# get dates in R format
dat2$dateR <- as.Date(with(dat2, paste(Year, Month, 15, sep = "-")))
glimpse(dat2)
dat2Sub <- dat2 %>% filter(site == "ABUR" |
                           site == "AQUE" |
                           site == "CARP" |
                           site == "IVEE" |
                           site == "MOHK" |
                           site == "NAPL")
names(dat2Sub)

qplot(dateR, mean_temp_c, data = dat2Sub, facets = ~site, geom = "line", 
      xlab = "Year", ylab = "Mean temperature (C)", color = site) + 
  theme(legend.position = "none")
ggsave("../figs/monthlyTemp_bySiteSub.pdf")

qplot(dateR, max_temp_c, data = dat2Sub, facets = ~site, geom = "line", 
      xlab = "Year", ylab = "Max temperature (C)", color = site) + 
  theme(legend.position = "none")
ggsave("../figs/monthlyMaxTemp_bySiteSub.pdf")

qplot(dateR, min_temp_c, data = dat2Sub, facets = ~site, geom = "line", 
      xlab = "Year", ylab = "Min temperature (C)", color = site) + 
  theme(legend.position = "none")
ggsave("../figs/monthlyMinTemp_bySiteSub.pdf")
