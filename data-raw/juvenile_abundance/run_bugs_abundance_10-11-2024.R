library(rstan)
library(tidyverse)
library(SRJPEmodel)
library(SRJPEdata)
library(R2WinBUGS)

# pCap model
pCap_inputs <- prepare_pCap_inputs(mainstem = FALSE)
pCap <- fit_pCap_model(pCap_inputs$inputs)

# save, fit 10-11-2024
# saveRDS(pCap, "~/Downloads/pCap_model.rds")

# read in
# pCap <- readRDS("~/Downloads/pCap_model.rds")

abundance_inputs <- prepare_abundance_inputs("ubc", 2018, effort_adjust = T)
lt_pCap_Us <- generate_lt_pCap_Us(abundance_inputs, pCap)

# BUGS --------------------------------------------------------------------
abundance <- fit_abundance_model_BUGS(abundance_inputs, lt_pCap_Us,
                                      here::here("model_files", "abundance_model.bug"))

# STAN - not working for right now ----------------------------------------
abundance <- fit_abundance_model(abundance_inputs$inputs, pCap, lt_pCap_Us)
