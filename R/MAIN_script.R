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


# color choices for figures ----
source("R/color_choice.R")


# Generate figure 1 ---- 

# +script for the function to make figure 1 -----
source("output_data/plots/stacked_bar_plot.R")

# +script that makes and saves figure 1 ----- 
source("figures/fig1_attr_sample.R")


# Generate figure 2 ----- 

# +script for the function to make figure 2 ----- 
source("figures/most_likely_source.R")

# +script that makes and saves figure 2 -----
source("figures/fig2_linking_variant_effect_and_signature.R")










