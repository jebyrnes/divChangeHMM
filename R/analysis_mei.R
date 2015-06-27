#################################################
# Analyze ENSO data from NOAA
# Author: Robin Elahi
# Date: 150626
#################################################

# Multivariate Enso Index from:
#http://www.esrl.noaa.gov/psd/enso/mei/
  
# load packages
library(reshape2)
library(ggplot2)
theme_set(theme_bw(base_size = 16))


# load mei data
dat <- read.csv("../data/mei_150610.csv")
head(dat)

# go from wide to long format
dat2 <- melt(data = dat, id.vars = "year")
head(dat2)

# create month column
dat2$month <- as.numeric(substr(dat2$variable, start = 2, stop = 2))
glimpse(dat2)

dat2$dateR <- ISOdate(year = dat2$year, month = dat2$month, day = 1)

head(dat2)

# plot
qplot(dateR, value, data = dat2[dat2$year > 1990, ], geom = "line", 
      xlab = "Year", ylab = "Multivariate ENSO index")
ggsave("../figs/meiPlot.pdf")

# Did the 1998 El Nino cause lagged shifts in the fish community between 2002-2003?

