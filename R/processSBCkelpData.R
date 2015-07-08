##############
# Boil the SBC kelp data down to meaningful annual values
#
#
# Author: Jarrett Byrnes
#
# Note -kelp highly processed reduced file from JEKB
# on another process. For original post-processing code
# email jarrett.byrnes@umb.edu
# File is to large for my github account
##############

library(readr)
library(lubridate)
library(tidyr)
library(dplyr)


kelp_data <- read.csv("../data/temporal_data_sbc_lter_longterm_community.csv")


kelp_data_annual<- kelp_data %>%
  filter(Taxon=="Macrocystis pyrifera") %>%
  rename(Year = Sample.Year) %>%
  separate(Sample.ID, c("Transect", "Quad", "Side", "End"), sep="\\.") %>%
  group_by(Year, Site, Transect) %>%
  summarize(stipe_density = mean(Stipe.Density.num.per.sq.m, na.rm=T))
  
qplot(Year, stipe_density, group=paste(Site, Transect), 
      data=kelp_data_annual) + 
  facet_wrap(~Site) +
  geom_line()

write_csv(kelp_data_annual, "../data/sbc_macrocystis_annual.csv")
