# This script prepares data for odeq's air toxics report
# This assumes you downloaded data from the epa aqs web app. using an amp501 report
# aqs and supporting links can be found at: https://www.epa.gov/aqs
# details on the query are included with aqs output and those details are 
# in the .pdf located with the data in the input directory 
# see toxics_tools.R for technical details

##########################################
# Part 0: libraries, functions, and setup
##########################################
library(readr)
source("toxics_tools.R")

# this is top directory for processing
root_dir <- getwd()

##########################################
# Part 1: AQS processing
##########################################
# load data and codes
# cross tables are from: https://www.epa.gov/aqs/aqs-code-list
cross_table <- load_cross_tab("./input/supporting_cross_tables")

# load the EPA AQS data from the AMP501 Raw data report
or_state_air_toxics <- load_tox_dat(root_dir, "/input/air_toxics/or_state_air_toxics_a501/or_state_air_toxics_a501-0.txt")
portland_air_toxics <- load_tox_dat(root_dir, "/input/air_toxics/portland_air_toxics_a501/portland_air_toxics_a501-0.txt")

# build the dataset and link codes
air_toxics <- bind_rows(or_state_air_toxics, portland_air_toxics)
air_toxics <- left_join(air_toxics, cross_table$air_toxics_sites, by = "epa_id", suffix = c(".data", ".method"))
air_toxics <- left_join(air_toxics, cross_table$methods_all, by = c("parameter_code", "method_code"), suffix = c(".data", ".method"))
air_toxics <- left_join(air_toxics, cross_table$duration, by = c("sample_duration" = "duration_code"), suffix = c(".all", ".duration"))
air_toxics <- left_join(air_toxics, cross_table$units, by = c("unit" = "unit_code"), suffix = c(".all", ".unit"))
air_toxics <- left_join(air_toxics, cross_table$parameters, by = "parameter_code", suffix = c("", ".parameter"))
air_toxics <- left_join(air_toxics, cross_table$Analyte_Master_LookupTable, by = "cas_number", suffix = c("", ".deq"))

# filter down to air toxic's analytes
air_toxics <- air_toxics %>% filter(!is.na(analyte_name_deq)) 

# link qualifiers 
qualifiers <- cross_table$qualifiers %>% select(qualifier_code, qualifier_description, qualifier_type)

air_toxics <- left_join(air_toxics, qualifiers, by = c("null_data_code" = "qualifier_code"), suffix = c("", ".null"))
air_toxics <- left_join(air_toxics, qualifiers, by = c("qualifier___1" = "qualifier_code"), suffix = c("", ".1"))
air_toxics <- left_join(air_toxics, qualifiers, by = c("qualifier___2" = "qualifier_code"), suffix = c("", ".2"))
air_toxics <- left_join(air_toxics, qualifiers, by = c("qualifier___3" = "qualifier_code"), suffix = c("", ".3"))
air_toxics <- left_join(air_toxics, qualifiers, by = c("qualifier___4" = "qualifier_code"), suffix = c("", ".4"))
air_toxics <- left_join(air_toxics, qualifiers, by = c("qualifier___5" = "qualifier_code"), suffix = c("", ".5"))
air_toxics <- left_join(air_toxics, qualifiers, by = c("qualifier___6" = "qualifier_code"), suffix = c("", ".6"))
air_toxics <- left_join(air_toxics, qualifiers, by = c("qualifier___7" = "qualifier_code"), suffix = c("", ".7"))
air_toxics <- left_join(air_toxics, qualifiers, by = c("qualifier___8" = "qualifier_code"), suffix = c("", ".8"))
air_toxics <- left_join(air_toxics, qualifiers, by = c("qualifier___9" = "qualifier_code"), suffix = c("", ".9"))
air_toxics <- left_join(air_toxics, qualifiers, by = c("qualifier___10" = "qualifier_code"), suffix = c("", ".10"))

# trace and status
air_toxics$proc_date <- Sys.Date() %>% gsub("-", "", .)
air_toxics$database <- "aqs"
air_toxics$status <- "final"

# save 
proc <- Sys.Date() %>% gsub("-", "", .)
write.csv(air_toxics, file = paste0(root_dir, "/output/tables/air_toxics_compiled_", proc, ".csv"), row.names = FALSE)

##########################################
# Part 2: Repository data processing
##########################################
# Data was exported from the Repository database (DEQs) using Bruces tool 
# There was some light pre-processing using Bruce's tool (from the query) and in excel before import here.  Pre-processing included:
# (1) did not export invalidated data --> only exported data w/ data quality = A, B, or F
# (2) updated unit format to remove mu and superscripts, etc in excel
# (3) exported primary monitor only - poc = 7 for air toxics
repository_dat <- read_csv("./input/air_toxics/Repository_Data_2019-2021.csv", 
                           col_types = cols(Sampled_Date = col_date(format = "%d-%b-%y"), 
                           Sampled = col_character()))

###################
# 2021 fixes only!
# there are a few issues with the CASNum in the repository
# fix that before using the cas to merge
###################
missing_cas <- is.na(repository_dat$CASNum)
foc_analyte <- repository_dat$Analyte == "Hexavalent Chromium [Cr(VI)]"
repository_dat$CASNum[missing_cas&foc_analyte] <- "18540-29-9"
foc_analyte <- repository_dat$Analyte == "1,4-Dimethylbenzene + 1,3-Dimethylbenzene"
repository_dat$CASNum[missing_cas&foc_analyte] <- "108-38-3"
foc_analyte <- repository_dat$CASNum == "1985-5-0"
repository_dat$CASNum[foc_analyte] <- "198-55-0"
foc_analyte <- repository_dat$Analyte == "Benzo(e)pyrene"
repository_dat$CASNum[foc_analyte] <- "192-97-2"

# add analyte group
repository_dat$analyte_group <- as.character('')
foc_group <- repository_dat$Analysis == "Air Toxics Volatiles by GCMS TO-15"
repository_dat$analyte_group[foc_group] <- "VOCs"

foc_group <- repository_dat$Analysis == "Air Toxics Carbonyls by HPLC TO-11A"
repository_dat$analyte_group[foc_group] <- "Carbonyls"
foc_group <- repository_dat$Analysis == "Air Toxics Carbonyls by HPLC TO-11A (2019)"
repository_dat$analyte_group[foc_group] <- "Carbonyls"

foc_group <- repository_dat$Analysis == "Air Toxics PAH by GCMS D6209 13"
repository_dat$analyte_group[foc_group] <- "PAHs"
foc_group <- repository_dat$Analysis == "Air Toxics PAH by GCMS D6209 13 (LRAPA)"
repository_dat$analyte_group[foc_group] <- "PAHs"
foc_group <- repository_dat$Analysis == "Air Toxics PAH by GCMS D6209 13 (2020)"
repository_dat$analyte_group[foc_group] <- "PAHs"

foc_group <- repository_dat$Analysis == "Air Toxics Hexavalent Chromium by IC"
repository_dat$analyte_group[foc_group] <- "Metals"
foc_group <- repository_dat$Analysis == "Air Toxics Metals (HighVol) by ICPMS"
repository_dat$analyte_group[foc_group] <- "Metals"
foc_group <- repository_dat$Analysis == "Air Toxics Metals (LoVol) by ICPMS"
repository_dat$analyte_group[foc_group] <- "Metals"

repository_dat <- left_join(repository_dat, cross_table$Analyte_Master_LookupTable, 
                            by = c("analyte_group" = "analyte_group", "CASNum" = "cas_number"))

tox_sites <- cross_table$air_toxics_sites %>% select(project, epa_id)
tox_sites$epa_id <- as.numeric(tox_sites$epa_id)

# added 13JUN2022 to remove sites that we do not want to process
sites_to_process <- repository_dat$Project %in% tox_sites$project
repository_dat <- repository_dat[sites_to_process,]

repository_dat <- left_join(repository_dat, tox_sites, by = c("Project" = "project"))

repository_dat$station_description <- repository_dat$Station_Description
repository_dat$date <- repository_dat$Sampled_Date

repository_dat$year <- year(repository_dat$date)
repository_dat$month <- month(repository_dat$date)
repository_dat$quarter <- as.numeric(NA)
foc <- repository_dat$month == 1 | repository_dat$month == 2 | repository_dat$month == 3
repository_dat$quarter[foc] <- 1
foc <- repository_dat$month == 4 | repository_dat$month == 5 | repository_dat$month == 6
repository_dat$quarter[foc] <- 2
foc <- repository_dat$month == 7 | repository_dat$month == 8 | repository_dat$month == 9
repository_dat$quarter[foc] <- 3
foc <- repository_dat$month == 10 | repository_dat$month == 11 | repository_dat$month == 12
repository_dat$quarter[foc] <- 4

# non-detects
repository_dat$sample_value <- as.numeric(NA)
non_detect <- repository_dat$Result == "ND"
repository_dat$sample_value[non_detect] <- 0 # insert a 0 for non-detects
repository_dat$sample_value[!non_detect] <- as.numeric(repository_dat$Result[!non_detect])

repository_dat$alternate_method_detectable_limit <- as.numeric(repository_dat$MDL)
repository_dat$poc <- 7 # air toxics primary monitors are poc 7

repository_dat$parameter <- repository_dat$analyte_name_deq
repository_dat$cas_number <- repository_dat$CASNum
repository_dat$null_data_code <- NA

repository_dat$units.unit <- as.character('')
foc <- repository_dat$Units == "ug_m3 STP"
repository_dat$units.unit[foc] <- "Micrograms/cubic meter (25 C)"
foc <- repository_dat$Units == "ug_m3 LTP"
repository_dat$units.unit[foc] <- "Micrograms/cubic meter (LC)"
foc <- repository_dat$Units == "ng_m3 STP"
repository_dat$units.unit[foc] <- "Nanograms/cubic meter (25 C)"
foc <- repository_dat$Units == "ng_m3 LTP"
repository_dat$units.unit[foc] <- "Nanograms/cubic meter (LC)"
foc <- repository_dat$Units == "ppbv"
repository_dat$units.unit[foc] <- "Parts per billion"

# trace and status
repository_dat$proc_date <- Sys.Date() %>% gsub("-", "", .)
repository_dat$database <- "repository"
repository_dat$status <- "final"

# save 
proc <- Sys.Date() %>% gsub("-", "", .)
write.csv(repository_dat, file = paste0(root_dir, "/output/tables/air_toxics_repository_", proc, ".csv"), row.names = FALSE)

