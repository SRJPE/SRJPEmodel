# example script for fitting the pCap and abundance models (now in separate STAN files)
# 11-27-2024

# uncomment these if you need to
#remotes::install_github("SRJPE/SRJPEdata")
#remotes::install_github("SRJPE/SRJPEmodel")
library(rstan)
library(tidyverse)
library(SRJPEmodel)
library(SRJPEdata)

# example for lower clear creek, 2018

# prepare data for pCap model
# see ?prepare_inputs_pCap_abundance_STAN
inputs <- prepare_inputs_pCap_abundance_STAN(SRJPEdata::weekly_juvenile_abundance_catch_data,
                                             SRJPEdata::weekly_juvenile_abundance_efficiency_data,
                                             site = "lcc", run_year = 2018, effort_adjust = T)

# look at what we have for pCap
inputs$pCap_inputs

# look at what we have for abundance (except lt_pCap_U)
inputs$abundance_inputs

# call pCap model. this function just takes in the inputs and calls the
# appropriate version of the pCap model
pCap <- fit_pCap_model(inputs$pCap_inputs)

# check outputs for convergence
rstan::summary(pCap, pars = c("lt_pCap_U"))$summary |>
  data.frame() |>
  filter(Rhat > 1.05)

# call abundance model - we need to pass in pCap to get the lt_pCap_U estimates
# and pass those in as data
abundance <- fit_abundance_model(inputs$abundance_inputs, pCap)

# check for convergence
rstan::summary(abundance, pars = c("lt_pCap_U"))$summary |>
  data.frame() |>
  filter(Rhat > 1.05)

rstan::summary(abundance, pars = c("N"))$summary |>
  data.frame() |>
  filter(Rhat > 1.05)

# this is an automated plot function I created. If it's not working, you can
# get data from SRJPEdata and estimates from the stanfit objects to create a plot
# that is helpful for you!
diagnostic_plots_split("lcc", 2018, abundance)


