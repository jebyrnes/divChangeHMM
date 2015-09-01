##############
# Process the SBC wave data down to meaningful annual values
#
#
# Author: Jarrett Byrnes
##############

library(readr)
library(lubridate)
library(dplyr)

cover <- read_csv("../data/cover_all_years_20140902.csv")
cover <- cover %>% filter(GROUP =="SUBSTRATE")

cover_summary <- cover %>%
  group_by(YEAR, MONTH, TRANSECT, SITE) %>%
  filter(TAXON_GENUS != "SAND") %>%
  summarise(Hard_Substrate_Percent = sum(PERCENT_COVER))

names(cover_summary) <- c("Year", "Month", "Transect", "Site", "Hard_Substrate_Percent")
write_csv(cover_summary, "../data/sbc_hard_substrate.csv")


