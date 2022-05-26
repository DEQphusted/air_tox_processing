# air_tox_processing

This project provides R scripts for processing the air toxics data collected by the Air Quality Monitoring section at the ODEQ lab.

These scripts require the user to first download the data from the EPA AQS database and the Repository database located at the Hillsboro Lab.

Details on obtaining data:
downloaded data from the epa aqs web app. 
use an amp501 raw data report
aqs and supporting linksusage info can be found at: https://www.epa.gov/aqs
specific details on the amp501 query are included with aqs output file and those details are 
in the .pdf located with the data in the input directory 

download data from the Repository database (DEQs) using Bruces tool 
There was some light processing using Bruce's tool (via the query request) and in excel before import into R.  
Pre-processing included:
(1) did not export invalidated data --> only exported data w/ data quality = A, B, or F
(2) updated unit format to remove mu and superscripts, etc in excel
(3) exported primary monitor only - poc = 7 for air toxics

Put the downloaded data into the correct directories located in the input folder
Make an output folder - name it "output"
Make subfolders in the output folder for figures ["figures"] and tables ["tables"] 

Run the R scripts in the following order
(1) air_toxics_load_ad_merge_step1.R - this preps air toxics data to a consistent format for analysis
(2) particulate_load_and_merge_step2.R - this preps neph pm25est data that is used to identify wildfires for the report
(3) air_toxics_analysis_step3.R - this does the annual statistics
(4) figure_maker_step4.R - this makes time series figures for all the site X pollutant combinations
