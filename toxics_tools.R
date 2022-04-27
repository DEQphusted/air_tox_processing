# This script contains functions to support the odeq air toxics report

# The data is from the epa's aqs database 
# https://www.epa.gov/aqs
# The codes are from:
# https://www.epa.gov/aqs/aqs-code-list

# To get data from aqs:
# - go to https://www.epa.gov/aqs, launch web app, and sign in
# - set Criteria Set
# -- AMP501 report
# -- Send via Email
# -- Report Output = WORKFILE
# -- name file
# - set Data selection tab 
# -- site
# -- Pollutant Type = ALL
# -- Data Criteria - set date range

##########################################

# (1) library & functions

##########################################

library(readr)
library(dplyr)
library(lubridate)

load_tox_dat <- function(root2get, file2get){
  data <- read_delim(paste0(root2get, file2get), "|", escape_double = FALSE, 
                     col_types = cols(Parameter = col_double(), 
                                      POC = col_double(), 
                                      `Sample Duration` = col_character(), 
                                      Unit = col_double(), 
                                      Method = col_double(), 
                                      Date = col_character(), 
                                      `Start Time` = col_character(), 
                                      `Alternate Method Detectable Limit` = col_double(), 
                                      Uncertainty = col_double()), trim_ws = TRUE)
  
  # trim extra rows
  st_row <- 2
  ed_row <- dim(data)[1] - 1
  data <- data[st_row:ed_row,] # trim off first and last row
  
  # some changes to the header make this easier to work with
  names(data) <- gsub("-", " ", names(data)) 
  names(data) <- tolower(gsub(" ", "_", names(data))) 
  
  # parameter is the parameter code in data
  parameter_code_id <- grepl('\\<parameter\\>', names(data))
  names(data)[parameter_code_id] <- 'parameter_code'
  
  # method is the method code in data
  method_code_id <- grepl('\\<method\\>', names(data))
  names(data)[method_code_id] <- 'method_code'
  
  data$date <- as.Date(data$date, format = "%Y%m%d")
  data$year <- year(data$date)
  data$month <- month(data$date)
  data$quarter <- as.numeric(NA)
  data$quarter[data$month == 1 | data$month == 2 | data$month == 3] <- 1
  data$quarter[data$month == 4 | data$month == 5 | data$month == 6] <- 2
  data$quarter[data$month == 7 | data$month == 8 | data$month == 9] <- 3
  data$quarter[data$month == 10 | data$month == 11 | data$month == 12] <- 4
  
  data$epa_id <- as.numeric(data$state_code) * 10000000 +  
    as.numeric(data$county_code) * 10000 + 
    as.numeric(data$site_id)
  
  return(data)
}

load_cross_tab <- function(root2get){
  air_toxics_sites <- read_csv(paste0(root2get, "/air_toxics_sites.csv"), 
                               col_types = cols(site = col_character(),
                                                epa_id = col_double(), 
                                                station_description = col_character(), 
                                                street_address = col_character(),
                                                city = col_character(),
                                                county = col_character(),
                                                latitude = col_double(), 
                                                longitude = col_double()))
  
  names(air_toxics_sites) <- tolower(gsub(" ", "_", names(air_toxics_sites)))
  
  methods_all <- read_csv(paste0(root2get, "/methods_all.csv"), 
                          col_types = cols(`Parameter Code` = col_double(), 
                                           `Method Code` = col_double(), 
                                           `Method Type` = col_character(), 
                                           `Reference Method ID` = col_character(), 
                                           `Equivalent Method` = col_character()))
  
  names(methods_all) <- tolower(gsub(" ", "_", names(methods_all)))
  
  qualifiers <- read_csv(paste0(root2get, "/qualifiers.csv"), 
                         col_types = cols(`Qualifier Code` = col_character(), 
                                          `Qualifier Description` = col_character(), 
                                          `Qualifier Type` = col_character(), 
                                          `Qaulifier Type Code` = col_character(), 
                                          `Still Active` = col_character(),
                                          `Legacy Code` = col_double()))
  
  names(qualifiers) <- tolower(gsub(" ", "_", names(qualifiers)))
  
  duration <- read_csv(paste0(root2get, "/durations.csv"), 
                       col_types = cols(`Duration Code` = col_character(), 
                                        `Duration Description` = col_character(), 
                                        `Observed or Calculated` = col_character()))
  
  names(duration) <- tolower(gsub(" ", "_", names(duration)))
  
  units <- read_csv(paste0(root2get, "/units.csv"), 
                    col_types = cols(`Unit Code` = col_double(), 
                                     Units = col_character()))
  
  names(units) <- tolower(gsub(" ", "_", names(units)))
  
  parameters <- read_csv(paste0(root2get, "/parameters.csv"))
  names(parameters) <- tolower(gsub(" ", "_", names(parameters)))
  
  Analyte_Master_LookupTable <- read_csv(paste0(root2get, "/Analyte_Master_LookupTable.csv"), 
                                         col_types = cols(pah_toxicity_equivalent_factor = col_double(), 
                                                          oregon_abc_ug_m3 = col_double()))
  
  cross_table <- list(air_toxics_sites, methods_all, qualifiers, duration, units, parameters, Analyte_Master_LookupTable)
  names(cross_table) <- c("air_toxics_sites", "methods_all", "qualifiers", "duration", "units", "parameters", "Analyte_Master_LookupTable")
  
  return(cross_table)
}