
library(readr)
library(tidyr)

air_toxics <- read_csv("./output/tables/air_toxics_summary__repository_20220616.csv")

air_toxics <- air_toxics %>% select(year, epa_id, Project, analyte_group, analyte_name_deq, poc, conc_benchmark_mean)

air_toxics_pivot <- pivot_wider(data = air_toxics,
                               id_cols = c(epa_id, Project, analyte_group, analyte_name_deq, poc),
                               names_from = year,
                               values_from = conc_benchmark_mean)


# write the data out
proc <- Sys.Date() %>% gsub("-", "", .)
write.csv(air_toxics_pivot, file = paste0("./output/tables/air_toxics_pivot_", proc, ".csv"), row.names = FALSE)
