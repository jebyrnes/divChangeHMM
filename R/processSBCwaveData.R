##############
# Process the SBC wave data down to meaningful annual values
#
#
# Author: Jarrett Byrnes
##############

library(readr)
library(lubridate)
library(dplyr)

waves <- read_csv("../data/cdip_mop_hs_tp_daily_20130528.csv")


wave_summary_monthly <- waves %>%
  mutate(Month = month(iso_date), Year = year(iso_date),mop=site) %>%
  group_by(Year, Month, mop) %>%
  summarise(mean_waveheight = mean(avg_hs_daily, na.rm=T),
            max_waveheight = max(max_hs_daily, na.rm=T),
            mean_max_waveheight=mean(max_hs_daily, na.rm=T)) %>% 
  ungroup()

write_csv(wave_summary_monthly, "../data/sbc_wave_summary_monthly.csv")


wave_summary_annual <- wave_summary_monthly %>%
  mutate(interval_year = ifelse(Month<7, Year, Year+1)) %>%
  group_by(interval_year, mop) %>% 
  summarise(mean_waveheight = mean(mean_waveheight, na.rm=T),
            max_waveheight = max(max_waveheight, na.rm=T),
            mean_max_waveheight=mean(mean_max_waveheight, na.rm=T)) %>% 
  ungroup()

write_csv(wave_summary_annual, "../data/sbc_wave_summary_annual.csv")


wave_summary_3prev <- wave_summary_monthly %>%
  mutate(interval_year = ifelse(Month<7, Year, Year+1)) %>%
  filter(Month %in% c(5,6,7)) %>%
  group_by(interval_year, mop) %>% 
  summarise(mean_waveheight = mean(mean_waveheight, na.rm=T),
            max_waveheight = max(max_waveheight, na.rm=T),
            mean_max_waveheight=mean(mean_max_waveheight, na.rm=T)) %>% 
  ungroup()

write_csv(wave_summary_3prev, "../data/sbc_wave_summary_MayToJuly.csv")
