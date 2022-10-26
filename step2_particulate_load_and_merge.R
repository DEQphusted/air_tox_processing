# This script prepares particulate data for odeq's air toxics report

rm(list = ls()) # remove everything

source("step0_toxics_tools.R")

root_dir <- getwd()

# load data and codes
cross_table <- load_cross_tab("./input/supporting_cross_tables")

or_state_particulate <- load_tox_dat(root_dir, "/input/particulate_for_air_tox/or_state_pm_tox/or_state_pm_tox-0.txt")
portland_particulate <- load_tox_dat(root_dir, "/input/particulate_for_air_tox/portland_state_pm_tox/portland_state_pm_tox-0.txt")

# build the dataset and link the codes
particulates <- bind_rows(or_state_particulate, portland_particulate)
names(particulates)[names(particulates) == "epa_id"] <- "particulate_link_epa_id"
particulates <- left_join(particulates, cross_table$methods_all, by = c("parameter_code", "method_code"))
particulates <- left_join(particulates, cross_table$air_toxics_sites, by = "particulate_link_epa_id")

particulates$sample_value <- as.double(particulates$sample_value)

particulate_daily <- particulates %>% group_by(epa_id, date, poc) %>%
  summarize(pm_site_epa_id = unique(particulate_link_epa_id),
            pm25 = mean(sample_value, na.rm = TRUE),
            pm_nobs = sum(!is.na(sample_value)),
            .groups = "drop")

# too few observations
bad <- particulate_daily$pm_nobs < 18 & particulate_daily$pm25 < 35
particulate_daily$pm25[bad] <- NA
bad <- particulate_daily$pm_nobs < 12
particulate_daily$pm25[bad] <- NA

particulate_daily$year <- year(particulate_daily$date)
particulate_daily$month <- month(particulate_daily$date)
particulate_daily$day <- day(particulate_daily$date)

particulate_daily$wildfire_impact <- "none"
wildfire_period <- c(particulate_daily$month == 6 | 
                     (particulate_daily$month == 7 & particulate_daily$day != 4) |
                      particulate_daily$month == 8 |
                      particulate_daily$month == 9 |
                     (particulate_daily$month == 10 & particulate_daily$day <= 12))
high_pm25 <- particulate_daily$pm25 >= 25
particulate_daily$wildfire_impact[wildfire_period&high_pm25] <- "wildfire"

# trace and status
particulate_daily$proc_date <- Sys.Date() %>% gsub("-", "", .)
particulate_daily$database <- "aqs"
particulate_daily$status <- "final"

# save 
proc <- Sys.Date() %>% gsub("-", "", .)
write.csv(particulate_daily, file = paste0(root_dir, "/output/tables/particulates_daily_", proc, ".csv"), row.names = FALSE)

