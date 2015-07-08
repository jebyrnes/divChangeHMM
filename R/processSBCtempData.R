##############
# Boil the SBC temperature data down to meaningful annual values
#
# Boil! Boil! I kill me!
#
# Author: Jarrett Byrnes
#
# Note - bottom temp downloaded from SBC internal data website
# File is to large for my github account
##############

library(readr)
library(lubridate)
library(dplyr)


temp_data <- read_csv("../data/bottom_temp_all_years.csv")

#Make monthly summaries
temp_summary_monthly <- temp_data %>%
  mutate(Month = month(Date), Year = year(Date)) %>%
  group_by(Year, Month, site) %>%
  summarise(mean_temp_c = mean(temp_c, na.rm=T),
            max_temp_c = max(temp_c, na.rm=T),
            min_temp_c = min(temp_c, na.rm=T)) %>% 
  ungroup()

write_csv(temp_summary_monthly, "../data/sbc_temp_summary_monthly.csv")

#Sampling is in July, so, we want to summarise from July y1 to June y2

temp_summary_annual <- temp_summary_monthly %>% 
  mutate(interval_year = ifelse(Month<7, Year, Year+1)) %>%
  group_by(interval_year, site) %>% 
  summarise(mean_temp_c = mean(mean_temp_c, na.rm=T),
            max_mean_temp_c = max(mean_temp_c, na.rm=T),
            min_mean_temp_c = min(mean_temp_c, na.rm=T))

write_csv(temp_summary_annual, "../data/sbc_temp_summary_annual.csv")

#Just the three previous months

temp_summary_3prev <- temp_summary_monthly %>% 
  filter(Month %in% c(5,6,7)) %>%
  group_by(Year, site) %>% 
  summarise(mean_temp_c = mean(mean_temp_c, na.rm=T),
            max_mean_temp_c = max(mean_temp_c, na.rm=T),
            min_mean_temp_c = min(mean_temp_c, na.rm=T))

write_csv(temp_summary_3prev, "../data/sbc_temp_summary_MayToJuly.csv")
