# example script for fitting the pCap and abundance models (now in separate STAN files)
# 2-10-2025

# re-install latest version of SRJPEmodel and SRJPEdata
#remotes::install_github("SRJPE/SRJPEdata")
#remotes::install_github("SRJPE/SRJPEmodel")
library(rstan)
library(tidyverse)
library(SRJPEmodel)
library(SRJPEdata)

# run pCap model --------------------------------------------------

# prepare pCap inputs
# set mainstem = TRUE if you want to fit for mainstem sites
pCap_inputs <- prepare_pCap_inputs(mainstem = FALSE)

# fit pCap model
pCap <- fit_pCap_model(pCap_inputs$inputs)

# save fit
pCap_filepath <- "~/Downloads/pCap_model.rds"
# saveRDS(pCap, pCap_filepath)

# if you've already fit pCap for your data, just read it in here
# pCap <- readRDS(pCap_filepath)


# run abundance model -----------------------------------------------------

# prepare abundance inputs
site <- "ubc"
run_year <- 2018
abundance_inputs <- prepare_abundance_inputs(site = site,
                                             run_year = run_year,
                                             effort_adjust = T)

# generate lt_pCap_Us from pCap model
lt_pCap_Us <- generate_lt_pCap_Us(abundance_inputs, pCap)

# fit abundance BUGS model
bugs_abundance_filepath <- "C:/Users/Liz/Documents/SRJPEmodel/model_files/abundance_model.bug"
bugs_directory <- "C:/Users/Liz/Documents/SRJPEmodel/data-raw/WinBUGS14"

abundance <- fit_abundance_model_BUGS(abundance_inputs, lt_pCap_Us,
                                      bugs_abundance_filepath,
                                      bugs_directory)

# helpful formatting functions
parameter_table <- extract_abundance_estimates(site, run_year,
                                               abundance_inputs,
                                               abundance)

# helpful plotting functions
plot_juv_data("ubc", 2018) # raw data
generate_diagnostic_plot_juv(site, run_year, parameter_table)


