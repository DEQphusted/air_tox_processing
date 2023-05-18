# make time series figures 

rm(list = ls()) # remove everything

library(readr)
library(dplyr)
library(ggplot2)
library(viridis)
library(tidyr)

air_toxics <- read_csv("./output/tables/air_toxics_repository_20230104.csv", col_types = cols(.default = "c"))
air_toxics$date <- as.Date(air_toxics$date)
air_toxics$sample_value <- as.double(air_toxics$sample_value)
air_toxics$mol_weight_g_mol <- as.double(air_toxics$mol_weight_g_mol)
air_toxics$alternate_method_detectable_limit <- as.double(air_toxics$alternate_method_detectable_limit)

all_analytes <- unique(air_toxics$Analyte)
all_projects <- unique(air_toxics$Project)

for (i_site in all_projects){
  
  for (i_analyte in all_analytes){
    
    dat <- air_toxics %>% filter(Analyte == i_analyte & Project == i_site)
    dat$doy <- yday(dat$date)
    dat$year <- as.factor(dat$year)
    
    title_name <- paste0("ODEQ provisional air toxics data from ", i_site)
    y_name <- paste0(tolower(i_analyte), " (", tolower(unique(dat$units.unit)), ")")
  
    fig <- ggplot(data=dat, aes(x=doy, y=sample_value, 
                            group = year,
                            color = year)) +
           theme_bw() +
           scale_color_viridis(discrete=TRUE) +
           ggtitle(title_name) + xlab("day of year") + ylab(y_name) +
           theme(plot.title = element_text(hjust = 0.5),
                 text = element_text(size = 24),
                 axis.title.y = element_text(size=24),
                 axis.title.x = element_text(size=24),
                 legend.text = element_text(size = 22),
                 axis.text.x = element_text(size = 22),
                 axis.text.y = element_text(size = 22),
                 legend.position=c(0.9,0.9)) +
           geom_line() +
           geom_point() 
  
    fle2save <- paste0(getwd(), "/output/figures/air_tox_", i_site, "_", i_analyte, ".png")
    ggsave(fle2save, plot = fig, width = 50, height = 30, units = "cm")
    
  } # end analyte
} # end i_site

