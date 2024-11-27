# write out script for Josh to show how I have been fitting the split models
library(rstan)
library(tidyverse)
library(SRJPEmodel)
library(SRJPEdata)

# send him also some model fit objects

# he can do some diagnostics and tell us what he needs for the report
# automate priors

# send him the STAN model code and the input data
# maybe just draft a function that produces the input data list/inits
# schedule meeting for next week

# example for lower clear creek, 2018

# prepare data for pCap model
# see ?prepare_inputs_pCap_abundance_STAN
inputs <- prepare_inputs_pCap_abundance_STAN(SRJPEdata::weekly_juvenile_abundance_catch_data,
                                             SRJPEdata::weekly_juvenile_abundance_efficiency_data,
                                             site = "lcc", run_year = 2018, effort_adjust = T)

# look at what we have for pCap
inputs$pCap_inputs

# call pCap model
pCap <- fit_pCap_model(inputs$pCap_inputs)

# get lt_pCap_U estimates and add to the data statement for abundance model
abundance <- fit_abundance_model(inputs$abundance_inputs, pCap)

diagnostic_plots_split("lcc", 2018, abundance)


