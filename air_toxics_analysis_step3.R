# compiles statistics for the odeq air toxics report
# run air_toxics_load_and_merge.R first

#####################################

# library and function

#####################################

rm(list = ls()) # remove everything

library(NADA) # for kaplan-meier non-detect means
library(dplyr)
library(tidyr)
library(lubridate)
library(readr)

do_km_mean <- function(df){
  dat <- cenfit(df$sample_value, df$below_mdl)
  dat <- mean(dat)[[1]]
  return(dat)
} 

#####################################

# load data and processing

#####################################
# aqs data
air_toxics <- read_csv("./output/tables/air_toxics_compiled_20220427.csv", col_types = cols(.default = "c"))
out_ext_dat <- ""
# repository data
#air_toxics <- read_csv("./output/tables/air_toxics_repository_20220427.csv", col_types = cols(.default = "c"))
#out_ext_dat <- "_repository_"
air_toxics$date <- as.Date(air_toxics$date)
air_toxics$sample_value <- as.double(air_toxics$sample_value)
air_toxics$mol_weight_g_mol <- as.double(air_toxics$mol_weight_g_mol)
air_toxics$alternate_method_detectable_limit <- as.double(air_toxics$alternate_method_detectable_limit)

# find non-detects and sub in 1/2 DL
air_toxics$best_val <- air_toxics$sample_value

air_toxics$below_mdl <- air_toxics$best_val < air_toxics$alternate_method_detectable_limit
void <- is.na(air_toxics$sample_value)
air_toxics$best_val[air_toxics$below_mdl & !void] <- air_toxics$alternate_method_detectable_limit[air_toxics$below_mdl & !void] * 0.5

#####################################

# do statistics

#####################################

## do quarterly average & filter count
air_toxics_quarter <- air_toxics %>% group_by(year, epa_id, analyte_group, analyte_name_deq, poc, quarter) %>%
  summarize(n_obs = sum(!is.na(best_val)),
            mean_val = mean(best_val, na.rm = TRUE),
            .groups = "drop")

air_toxics_quarter_avg <- air_toxics_quarter %>% group_by(year, epa_id, analyte_group, analyte_name_deq, poc) %>%
  summarize(quarter_mean_val = mean(mean_val, na.rm = TRUE),
            .groups = "drop")

air_toxics_quarter_tally <- air_toxics_quarter %>% select(year, epa_id, analyte_group, analyte_name_deq, poc, quarter, n_obs) %>% 
  pivot_wider(names_from = quarter, values_from = n_obs)

names(air_toxics_quarter_tally)[names(air_toxics_quarter_tally) == "1"] <- "nobs_q1"
names(air_toxics_quarter_tally)[names(air_toxics_quarter_tally) == "2"] <- "nobs_q2"
names(air_toxics_quarter_tally)[names(air_toxics_quarter_tally) == "3"] <- "nobs_q3"
names(air_toxics_quarter_tally)[names(air_toxics_quarter_tally) == "4"] <- "nobs_q4"

air_toxics_quarter_avg <- left_join(air_toxics_quarter_avg, air_toxics_quarter_tally, 
                                    by = c("year", "epa_id", "analyte_group", "analyte_name_deq", "poc"))

## summarize over the entire year
air_toxics_summary <- air_toxics %>% group_by(year, epa_id, analyte_group, analyte_name_deq, poc) %>%
  summarize(site_name = unique(station_description),
            start_date = min(date, na.rm = TRUE),
            end_date = max(date, na.rm = TRUE),
            parameter_epa = unique(parameter),
            cas_num = unique(cas_number),
            n_submitted = sum(!is.na(date)),
            n_void = sum(!is.na(null_data_code)),
            n_obs = sum(!is.na(best_val)), 
            n_below_mdl = sum(below_mdl, na.rm = TRUE),
            min_val = min(best_val, na.rm = TRUE),
            mean_val = mean(best_val, na.rm = TRUE),
            max_val = max(best_val, na.rm = TRUE),
            sd_val = sd(best_val, na.rm = TRUE),
            q10 = quantile(best_val, c(.10), na.rm = TRUE),
            q50 = median(best_val, na.rm = TRUE),
            q90 = quantile(best_val, c(.90), na.rm = TRUE),
            q95 = quantile(best_val, c(.95), na.rm = TRUE),
            unts = unique(units.unit),
            mol_wght_g_mol = unique(mol_weight_g_mol),
            pah_tox_eq_factor = unique(pah_toxicity_equivalent_factor),
            oregon_abc = unique(oregon_abc_ug_m3),
            hr24_rbc_cao = unique(hr24_rbc_cao_ug_m3),
            mdl_min = min(alternate_method_detectable_limit, na.rm = TRUE),
            mdl_mean = mean(alternate_method_detectable_limit, na.rm = TRUE),
            mdl_max = max(alternate_method_detectable_limit, na.rm = TRUE),
            .groups = "drop")

air_toxics_km_mean <- air_toxics %>% filter(!is.na(sample_value) & !is.na(below_mdl)) %>%
  group_by(year, epa_id, analyte_group, analyte_name_deq, poc) %>%
  summarize(km_mean = do_km_mean(.data),
            .groups = "drop") 
            
air_toxics_summary <- left_join(air_toxics_summary, air_toxics_quarter_avg,
                                by = c("year", "epa_id", "analyte_group", "analyte_name_deq", "poc"))

air_toxics_summary <- left_join(air_toxics_summary, air_toxics_km_mean,
                                by = c("year", "epa_id", "analyte_group", "analyte_name_deq", "poc"))

## unit conversions on annual mean and max by year
air_toxics_summary$mean_ug_m3 <- as.double(NA)
air_toxics_summary$max_ug_m3 <- as.double(NA)

good_units <- air_toxics_summary$unts == "Micrograms/cubic meter (25 C)" | air_toxics_summary$unts == "Micrograms/cubic meter (LC)"
air_toxics_summary$mean_ug_m3[good_units] <- air_toxics_summary$mean_val[good_units]
air_toxics_summary$max_ug_m3[good_units] <- air_toxics_summary$max_val[good_units]

ng_units <- air_toxics_summary$unts == "Nanograms/cubic meter (25 C)" | air_toxics_summary$unts == "Nanograms/cubic meter (LC)"
air_toxics_summary$mean_ug_m3[ng_units] <- air_toxics_summary$mean_val[ng_units] / 1000
air_toxics_summary$max_ug_m3[ng_units] <- air_toxics_summary$max_val[ng_units] / 1000

vocs <- air_toxics_summary$analyte_group == "VOCs" & air_toxics_summary$unts == "Parts per billion"
air_toxics_summary$mean_ug_m3[vocs] <- (air_toxics_summary$mean_val[vocs] * air_toxics_summary$mol_wght_g_mol[vocs] * 12.187)/298.15
air_toxics_summary$max_ug_m3[vocs] <- (air_toxics_summary$max_val[vocs] * air_toxics_summary$mol_wght_g_mol[vocs] * 12.187)/298.15

## write the data out
proc <- Sys.Date() %>% gsub("-", "", .)
write.csv(air_toxics_summary, file = paste0("./output/tables/air_toxics_summary_", out_ext_dat, proc, ".csv"), row.names = FALSE)

#############################################

## now do the wildfire analysis

#############################################

pm25 <- read_csv("./output/tables/particulates_daily_20220427.csv", col_types = cols(.default = "c"))
pm25$date <- as.Date(pm25$date)
pm25$pm25 <- as.double(pm25$pm25)

pm25 <- pm25 %>% select(epa_id, date, pm_site_epa_id, pm25, pm_nobs, wildfire_impact)

# there are data from pch, but these are sensor data --> use sel neph
good <- !is.na(pm25$epa_id)
pm25 <- pm25[good,]
# NOTE: sel neph is missing some observations

# join pm25 wildfire
air_toxics_wildfire <- left_join(air_toxics, pm25, by = c("epa_id", "date"))


## do quarterly average & filter count
air_toxics_quarter_wildfire <- air_toxics_wildfire %>% group_by(year, epa_id, analyte_group, analyte_name_deq, poc, quarter, wildfire_impact) %>%
  summarize(n_obs = sum(!is.na(best_val)),
            mean_val = mean(best_val, na.rm = TRUE),
            .groups = "drop")

air_toxics_quarter_wildfire_avg <- air_toxics_quarter_wildfire %>% group_by(year, epa_id, analyte_group, analyte_name_deq, poc, wildfire_impact) %>%
  summarize(quarter_mean_val = mean(mean_val, na.rm = TRUE),
            .groups = "drop")

air_toxics_quarter_wildfire_tally <- air_toxics_quarter_wildfire %>% select(year, epa_id, analyte_group, analyte_name_deq, poc, quarter, wildfire_impact, n_obs) %>% 
  pivot_wider(names_from = quarter, values_from = n_obs)

names(air_toxics_quarter_wildfire_tally)[names(air_toxics_quarter_wildfire_tally) == "1"] <- "nobs_q1"
names(air_toxics_quarter_wildfire_tally)[names(air_toxics_quarter_wildfire_tally) == "2"] <- "nobs_q2"
names(air_toxics_quarter_wildfire_tally)[names(air_toxics_quarter_wildfire_tally) == "3"] <- "nobs_q3"
names(air_toxics_quarter_wildfire_tally)[names(air_toxics_quarter_wildfire_tally) == "4"] <- "nobs_q4"

air_toxics_quarter_wildfire_avg <- left_join(air_toxics_quarter_wildfire_avg, air_toxics_quarter_wildfire_tally, 
                                    by = c("year", "epa_id", "analyte_group", "analyte_name_deq", "poc", "wildfire_impact"))


air_toxics_wildfire_summary <- air_toxics_wildfire %>% group_by(year, epa_id, analyte_group, analyte_name_deq, poc, wildfire_impact) %>%
  summarize(site_name = unique(station_description),
            start_date = min(date, na.rm = TRUE),
            end_date = max(date, na.rm = TRUE),
            parameter_epa = unique(parameter),
            cas_num = unique(cas_number),
            n_submitted = sum(!is.na(date)),
            n_void = sum(!is.na(null_data_code)),
            n_obs = sum(!is.na(best_val)), 
            n_below_mdl = sum(below_mdl, na.rm = TRUE),
            min_val = min(best_val, na.rm = TRUE),
            mean_val = mean(best_val, na.rm = TRUE),
            max_val = max(best_val, na.rm = TRUE),
            sd_val = sd(best_val, na.rm = TRUE),
            q10 = quantile(best_val, c(.10), na.rm = TRUE),
            q50 = median(best_val, na.rm = TRUE),
            q90 = quantile(best_val, c(.90), na.rm = TRUE),
            q95 = quantile(best_val, c(.95), na.rm = TRUE),
            unts = unique(units.unit),
            mol_wght_g_mol = unique(mol_weight_g_mol),
            pah_tox_eq_factor = unique(pah_toxicity_equivalent_factor),
            oregon_abc = unique(oregon_abc_ug_m3),
            hr24_rbc_cao = unique(hr24_rbc_cao_ug_m3),
            mean_pm25 = mean(pm25, na.rm = TRUE),
            mdl_min = min(alternate_method_detectable_limit, na.rm = TRUE),
            mdl_mean = mean(alternate_method_detectable_limit, na.rm = TRUE),
            mdl_max = max(alternate_method_detectable_limit, na.rm = TRUE),
            .groups = "drop")

air_toxics_wildfire_km_mean <- air_toxics_wildfire %>% filter(!is.na(sample_value) & !is.na(below_mdl)) %>%
  group_by(year, epa_id, analyte_group, analyte_name_deq, poc, wildfire_impact) %>%
  summarize(km_mean = do_km_mean(.data),
            .groups = "drop")

air_toxics_wildfire_summary <- left_join(air_toxics_wildfire_summary, air_toxics_quarter_wildfire_avg,
                                by = c("year", "epa_id", "analyte_group", "analyte_name_deq", "poc", "wildfire_impact"))

air_toxics_wildfire_summary <- left_join(air_toxics_wildfire_summary, air_toxics_wildfire_km_mean,
                                by = c("year", "epa_id", "analyte_group", "analyte_name_deq", "poc", "wildfire_impact"))

# unit conversions on annual mean and max by year
air_toxics_wildfire_summary$mean_ug_m3 <- as.double(NA)
air_toxics_wildfire_summary$max_ug_m3 <- as.double(NA)

good_units <- air_toxics_wildfire_summary$unts == "Micrograms/cubic meter (25 C)" | air_toxics_wildfire_summary$unts == "Micrograms/cubic meter (LC)"
air_toxics_wildfire_summary$mean_ug_m3[good_units] <- air_toxics_wildfire_summary$mean_val[good_units]
air_toxics_wildfire_summary$max_ug_m3[good_units] <- air_toxics_wildfire_summary$max_val[good_units]

ng_units <- air_toxics_wildfire_summary$unts == "Nanograms/cubic meter (25 C)" | air_toxics_wildfire_summary$unts == "Nanograms/cubic meter (LC)"
air_toxics_wildfire_summary$mean_ug_m3[ng_units] <- air_toxics_wildfire_summary$mean_val[ng_units] / 1000
air_toxics_wildfire_summary$max_ug_m3[ng_units] <- air_toxics_wildfire_summary$max_val[ng_units] / 1000

vocs <- air_toxics_wildfire_summary$analyte_group == "VOCs" & air_toxics_wildfire_summary$unts == "Parts per billion"
air_toxics_wildfire_summary$mean_ug_m3[vocs] <- (air_toxics_wildfire_summary$mean_val[vocs] * air_toxics_wildfire_summary$mol_wght_g_mol[vocs] * 12.187)/298.15
air_toxics_wildfire_summary$max_ug_m3[vocs] <- (air_toxics_wildfire_summary$max_val[vocs] * air_toxics_wildfire_summary$mol_wght_g_mol[vocs] * 12.187)/298.15

# write the data out
proc <- Sys.Date() %>% gsub("-", "", .)
write.csv(air_toxics_wildfire_summary, file = paste0("./output/tables/air_toxics_wildfire_summary_", out_ext_dat, proc, ".csv"), row.names = FALSE)

