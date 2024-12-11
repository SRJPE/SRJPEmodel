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
parameters <- c("lt_pCap_U", "pCap_U", "N", "Ntot", "sd.N", "sd.Ne", "lg_CumN")
Nmcmc = 2000
Nburnin = 500
Nthin = 2
Nchains = 3
abundance <- bugs(abundance_inputs$inputs$data,
                  abundance_inputs$inputs$inits,
                  parameters,
                  here::here("model_files", "abundance_model.bug"), debug = F,
                  n.chains = Nchains, n.burnin = Nburnin, n.thin = Nthin, n.iter = Nmcmc,
                  codaPkg = F, DIC = T, clearWD = T,
                  bugs.directory = here::here("data-raw", "WinBUGS14"))

# STAN - not working for right now ----------------------------------------
abundance <- fit_abundance_model(abundance_inputs$inputs, pCap, lt_pCap_Us)
