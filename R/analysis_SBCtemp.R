#################################################
# Analyze SBC lter temperature data
# Author: Robin Elahi
# Date: 150626
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

# Shoot, the drop in fish richness/diversity is from 2002 to 2003
# But the temperature data start in 2002, so we dont' have 
# a good sense of before-after

# Need to find temperature data that starts in 2001, at least.  





