## Main analysis script ----- 

library(cancereffectsizeR)


# download the data ----- 

## Only run once
# source("R/download_data.R")


# run cancereffectsizeR analysis on the data ---- 

## Only run if changes to data above, or changes to cesR pipeline
# source("R/cesR_analyses.R")


# load in cesR analyses ----
source("R/load_cesR_analyses.R")


# Run attributions analyses ----- 
source("R/attributions_calculations.R")









