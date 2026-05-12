# test_storage_workflow.R
# Minimal script to test the upload and download functions for a pCap all-sites
# fit and a single abundance model fit.

library(tidyverse)
library(SRJPEdata)
library(SRJPEmodel)
library(rstan)
library(R2WinBUGS)

# stan options
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# file path for bugs
bugs_directory  <- "C:/Users/Liz/Documents/SRJPEmodel/data-raw/WinBUGS14"

# connect to db
cfg <- config::get()
con <- DBI::dbConnect(
  RPostgres::Postgres(),
  dbname   = cfg$db_name,
  host     = cfg$db_host,
  port     = cfg$db_port,
  user     = cfg$db_user,
  password = cfg$db_password
)
on.exit(DBI::dbDisconnect(con), add = TRUE)


# test fit, upload, download ----------------------------------------------

test_site <- "ubc"
test_run_year <- 2020

# pcap
pCap_inputs <- prepare_pCap_inputs(model_type = "all_sites")
pCap_fit <- fit_pCap_model(pCap_inputs)

store_model_fit(
  con = con,
  model_fit_object = pCap_fit,
  model_inputs = pCap_inputs,
  description = paste("TEST — pCap all-sites", Sys.Date())
)

# abundance model
abundance_inputs <- prepare_abundance_inputs(
  site = test_site,
  run_year = test_run_year,
  pCap_model_type = "all_sites",
  pCap_model_object = pCap_fit,
  min_pCap_mult = 0.5
)

abundance_fit <- fit_abundance_model_BUGS(
  abundance_inputs,
  bugs_directory
)

store_model_fit(
  con              = con,
  model_fit_object = abundance_fit,
  model_inputs     = abundance_inputs,
  description      = paste("TEST —", test_site, test_run_year, "abundance", Sys.Date())
)


# download ----------------------------------------------------------------

# pCap models (all sites, one site tisdale, one site knights landing)
pCap_all_sites <- get_model_fit("pcap_all_sites")
pCap_tisdale <- get_model_fit("pcap_one_site", con = con,
                              site_selection = "tisdale")
pCap_kdl <- get_model_fit("pcap_one_site", con = con,
                          site_selection = "knights landing")

# abundance - one fit
abund_single <- get_model_fit("abundance", con = con,
                              site = test_site, run_year = test_run_year)

# abundance - all fits for each site/run year
all_abundance_fits <- get_many_model_fits(con, model_name = "abundance")

