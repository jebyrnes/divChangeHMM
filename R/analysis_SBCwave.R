#################################################
# Analyze SBC lter wave data
# Author: Robin Elahi
# Date: 150626
#################################################

library(ggplot2)
theme_set(theme_bw(base_size = 16))
library(dplyr)
options(dplyr.print_max = 1e9)

# Load daily wave data
dat<- read.csv("../data/cdip_mop_hs_tp_daily_20130528.csv")
names(dat)
head(dat)

# get dates in R format
dat$dateR <- as.Date(dat$iso_date)
glimpse(dat)

# quick plot
qplot(dateR, max_hs_daily, data = dat, facets = ~site, geom = "line")



