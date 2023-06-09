# -------------------------------------
# Author: Jake Diamond
# Purpose: To calculate uncertainty in young water fraction for NZ rivers
# Date: 2023-03-01
# -------------------------------------

# ################## How to Use this file ##################################
# Step 1: Install the "targets" package, a library to ensure reproducibility
# Step 2: Load the "targets" library
# To run: type "tar_make()" in the console and hit enter
# You can also type "tar_visnetwork()" and run it to see how all the functions
# and data are linked to eachother

# Created by use_targets().
# Follow the manual to check and run the pipeline:
# https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(here)

# Load packages and set up  -----------------------------------------------
# Set target options:
tar_option_set(
  # packages that your targets need to run
  packages = c("plotly", "lubridate", "htmltools",
               "here", "tidyr", "skimr", "tidyverse",
               "tidytable"),
  # default storage format
  format = "rds"
)

# Run the R scripts in the R/ folder with your custom functions:
source("src/R/functions.R") # Source other scripts as needed.


# List of targets to run--------------------------------------------------
# Do not change unless you know what you're doing!!!
# Replace the target list below with your own:
list(
  # Load raw isotope and discharge data for rivers
  tar_target(data_raw, read_csv(here("data", "raw", "WQ.iso.flow.data.csv"))),
  # Clean it to be in easy-to-use format
  tar_target(data_iso, clean_raw_fun(data_raw)),
  # load results of previous analysis, inlcuding sine fits to precip data
  tar_target(data_amp, load_Rdata(here("data",
                                       "river.site.est.Amplitude.Rdata"))),
  # Clean that data to be in easy-to-use format
  tar_target(amp_old, clean_fun_results(data_amp)),
  # Run the curve fitting and extract fit parameters and uncertainty
  tar_target(fit_riv, riv_analysis(data_iso)),
  # Calculate the young water fraction (Fwy) using the amplitude ratios
  # and estimate uncertainty with gaussian error propagation
  tar_target(Fwy, amp_ratio_fun(fit_riv, amp_old)),
  # Estimate uncertainty with monte carlo approach, 10000 runs
  tar_target(Fwy_MC, monte_carlo_fun(Fwy))
)

# To run: type "tar_make()" in the console and hit enter
