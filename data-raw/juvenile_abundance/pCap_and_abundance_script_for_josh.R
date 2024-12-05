# example script for fitting the pCap and abundance models (now in separate STAN files)
# 11-27-2024

# uncomment these if you need to
#remotes::install_github("SRJPE/SRJPEdata")
#remotes::install_github("SRJPE/SRJPEmodel")
library(rstan)
library(tidyverse)
library(SRJPEmodel)
library(SRJPEdata)



# run pCap general model --------------------------------------------------

# prepare data for pCap model
# see ?prepare_inputs_pCap_abundance_STAN
inputs_general <- prepare_inputs_pCap_abundance_STAN(SRJPEdata::weekly_juvenile_abundance_catch_data,
                                                     SRJPEdata::weekly_juvenile_abundance_efficiency_data,
                                                     site = "lcc", run_year = 2020, effort_adjust = T)

# look at what we have for pCap
inputs_general$pCap_inputs

# call pCap model. this function just takes in the inputs and calls the
# appropriate version of the pCap model.
# this is now modified to not produce any lt_pCap_Us
# we only need to call once
pCap <- fit_pCap_model(inputs_general$pCap_inputs)

# update this with your local filepath
local_filepath <- "~/Downloads/pCap_model.rds"
#saveRDS(pCap, file = local_filepath)


# run site-specific model -------------------------------------------------

# read in general pCap model fit object
pCap <- readRDS(local_filepath)

# start workflow that is specific to site/run
# ubc 2005 has catch for all weeks - fits
# lcc 2018 has catch for all weeks - fits
# lcc 2020 has catch for some weeks (tail end missing) - fits
# ubc 2018 has catch for some weeks (beginning missing) - DOES NOT FIT, INDEX ISSUES
# ubc 2006 has catch for some weeks (missing week 3) - DOES NOT FIT, INDEX ISSUES, lots of 0s
inputs_site <- prepare_inputs_pCap_abundance_STAN(SRJPEdata::weekly_juvenile_abundance_catch_data,
                                                  SRJPEdata::weekly_juvenile_abundance_efficiency_data,
                                                  site = "ubc", run_year = 2006,
                                                  effort_adjust = T,
                                                  # uncomment the below if you want to set a specific denominator for the prior
                                                  # default_lgN_prior_denominator = 0.0001
                                                  )

lt_pCap_Us <- generate_lt_pCap_Us(inputs_site$pCap_inputs, pCap)

# see outputs
lt_pCap_Us$lt_pCap_mu
lt_pCap_Us$lt_pCap_sd

# call abundance model - we need to pass in pCap to get the lt_pCap_U estimates
# and pass those in as data
abundance <- fit_abundance_model(inputs_site$abundance_inputs, pCap, lt_pCap_Us)

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
diagnostic_plots_split("ubc", 2018, abundance)


# table of sites to run for lcc and ubc -----------------------------------
trials_to_fit <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
  filter(life_stage %in% c("fry", "smolt")) |>
  mutate(filter_out = ifelse(is.na(life_stage) & count > 0, TRUE, FALSE)) |> # we do not want to keep NA lifestage associated with counts > 0
  filter(!filter_out,
         site %in% c("lcc", "ubc"),
         week %in% c(seq(45, 53), seq(1, 22))) |>
  mutate(count = round(count, 0),
         catch_standardized_by_hours_fished = round(catch_standardized_by_hours_fished, 0)) |>
  group_by(site, run_year) |>
  tally() |>
  arrange(desc(site), desc(run_year)) |>
  select(-n)


# functionalize and run ---------------------------------------------------

run_site_year <- function(site_arg, year_arg) {
  cli::cli_bullets(paste0("Fitting for ", site_arg, " and run year ", year_arg))

  inputs <- prepare_inputs_pCap_abundance_STAN(SRJPEdata::weekly_juvenile_abundance_catch_data,
                                               SRJPEdata::weekly_juvenile_abundance_efficiency_data,
                                               site = site_arg, run_year = year_arg, effort_adjust = T)

  pCap <- fit_pCap_model(inputs$pCap_inputs)

  # check outputs for convergence
  check_pCap <- rstan::summary(pCap, pars = c("lt_pCap_U"))$summary |>
    data.frame() |>
    filter(Rhat > 1.05)

  if(nrow(check_pCap) > 0) {
    return(check_pCap)
  }

  # call abundance model - we need to pass in pCap to get the lt_pCap_U estimates
  # and pass those in as data
  abundance <- fit_abundance_model(inputs$abundance_inputs, pCap)

  check_abundance <- bind_rows(rstan::summary(abundance, pars = c("lt_pCap_U"))$summary |>
                                 data.frame() |>
                                 filter(Rhat > 1.05),
                               rstan::summary(abundance, pars = c("N"))$summary |>
                                 data.frame() |>
                                 filter(Rhat > 1.05))

  if(nrow(check_abundance) > 0) {
    return(check_abundance)
  }

  diagnostic_plots_split(site_arg, year_arg, abundance)

  return(abundance)
}

trials_to_fit_lcc <- trials_to_fit |>
  filter(site == "lcc")

lcc <- purrr::pmap(list(trials_to_fit_lcc$site,
                        trials_to_fit_lcc$run_year),
                   run_site_year,
                   .progress = T)



