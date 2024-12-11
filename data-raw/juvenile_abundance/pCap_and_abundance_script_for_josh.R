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
                                                     site = "deer creek", run_year = 2002, effort_adjust = T)

# look at what we have for pCap
inputs_general$pCap_inputs

# call pCap model. this function just takes in the inputs and calls the
# appropriate version of the pCap model.
# this is now modified to not produce any lt_pCap_Us
# we only need to call once
pCap <- fit_pCap_model(inputs_general$pCap_inputs)

# update this with your local filepath
local_filepath <- "~/Downloads/pCap_model.rds"

# save to your local filepath
#saveRDS(pCap, file = local_filepath)

# run site-specific model -------------------------------------------------

# read in general pCap model fit object
pCap <- readRDS(local_filepath)

# start workflow that is specific to site/run
inputs_site <- prepare_inputs_pCap_abundance_STAN(SRJPEdata::weekly_juvenile_abundance_catch_data,
                                                  SRJPEdata::weekly_juvenile_abundance_efficiency_data,
                                                  site = "ubc", run_year = 2018,
                                                  effort_adjust = T,
                                                  #default_lgN_prior_denominator = 0.0001 # uncomment the below if you want to set a specific denominator for the prior
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
diagnostic_plots_split("lcc", 2018, abundance)


# table of sites to run for lcc and ubc  - LIZ CODE -----------------------------------
trials_to_fit <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
  filter(life_stage %in% c("fry", "smolt")) |>
  mutate(filter_out = ifelse(is.na(life_stage) & count > 0, TRUE, FALSE)) |> # we do not want to keep NA lifestage associated with counts > 0
  filter(!filter_out,
         #site %in% c("lcc", "ubc"),
         week %in% c(seq(45, 53), seq(1, 22))) |>
  mutate(count = round(count, 0),
         catch_standardized_by_hours_fished = round(catch_standardized_by_hours_fished, 0)) |>
  group_by(stream, site, run_year) |>
  tally() |>
  arrange(desc(site), desc(run_year)) |>
  select(-n)


# functionalize and run for lcc - LIZ CODE ---------------------------------------------------

run_site_year <- function(site_arg, year_arg, local_filepath) {
  cli::cli_bullets(paste0("Fitting for ", site_arg, " and run year ", year_arg))

  pCap <- readRDS(local_filepath)

  # start workflow that is specific to site/run
  inputs <- prepare_inputs_pCap_abundance_STAN(SRJPEdata::weekly_juvenile_abundance_catch_data,
                                                    SRJPEdata::weekly_juvenile_abundance_efficiency_data,
                                                    site = site_arg, run_year = year_arg,
                                                    effort_adjust = T,
                                                    #default_lgN_prior_denominator = default_denom # uncomment the below if you want to set a specific denominator for the prior
  )

  lt_pCap_Us <- generate_lt_pCap_Us(inputs$pCap_inputs, pCap)

  # call abundance model - we need to pass in pCap to get the lt_pCap_U estimates
  # and pass those in as data
  abundance <- fit_abundance_model(inputs$abundance_inputs, pCap, lt_pCap_Us)

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
  filter(site == "lcc",
         run_year != 2010)

# update this with your local filepath
local_filepath <- "~/Downloads/pCap_model.rds"
pCap <- readRDS(local_filepath)

lcc <- purrr::pmap(list(trials_to_fit_lcc$site,
                        trials_to_fit_lcc$run_year,
                        rep(local_filepath, nrow(trials_to_fit_lcc))),
                   run_site_year,
                   .progress = T)

# 2010 lcc
fit_index <- lapply(lcc, is.data.frame) |>
  unlist()
fit_index <- !fit_index

model_id <- trials_to_fit_lcc |>
  mutate(id = paste(site, run_year, sep = "_")) |>
  pull(id)

models_that_fit <- lcc[fit_index]
names(models_that_fit) <- model_id[fit_index]

saveRDS(models_that_fit, file = "~/Downloads/lcc_successful_fits.rds")

fits <- readRDS("~/Downloads/lcc_successful_fits.rds")

# models that converged ---------------------------------------------------

model_id[!fit_index]

# lcc 2021
diagnostic_plots_split("lcc", 2021, models_that_fit$lcc_2021)
# lcc 2018
diagnostic_plots_split("lcc", 2018, models_that_fit$lcc_2018)
# lcc 2017
diagnostic_plots_split("lcc", 2017, models_that_fit$lcc_2017)
# lcc 2016
diagnostic_plots_split("lcc", 2016, models_that_fit$lcc_2016)
# lcc 2015
diagnostic_plots_split("lcc", 2015, models_that_fit$lcc_2015)
# lcc 2011
diagnostic_plots_split("lcc", 2011, models_that_fit$lcc_2011)
# lcc 2006
diagnostic_plots_split("lcc", 2006, models_that_fit$lcc_2006)
# lcc 2005
diagnostic_plots_split("lcc", 2005, models_that_fit$lcc_2005)

# all of these seem to have rarely missing catch


# models that did not converge --------------------------------------------

model_id[!fit_index]

# lcc 2021
inputs_2021 <- prepare_inputs_pCap_abundance_STAN(SRJPEdata::weekly_juvenile_abundance_catch_data,
                                             SRJPEdata::weekly_juvenile_abundance_efficiency_data,
                                             site = "lcc", run_year = 2021,
                                             effort_adjust = T,
                                             #default_lgN_prior_denominator = 0.005 # uncomment the below if you want to set a specific denominator for the prior
)

lt_pCap_Us_2021 <- generate_lt_pCap_Us(inputs_2021$pCap_inputs, pCap)
abundance_2021 <- fit_abundance_model(inputs_2021$abundance_inputs, pCap, lt_pCap_Us_2021)

# lcc 2024
plot_juv_data("lcc", 2024)

inputs_2024 <- prepare_inputs_pCap_abundance_STAN(SRJPEdata::weekly_juvenile_abundance_catch_data,
                                                  SRJPEdata::weekly_juvenile_abundance_efficiency_data,
                                                  site = "lcc", run_year = 2024,
                                                  effort_adjust = T,
                                                  #default_lgN_prior_denominator = 0.005 # uncomment the below if you want to set a specific denominator for the prior
)

lt_pCap_Us_2024 <- generate_lt_pCap_Us(inputs_2024$pCap_inputs, pCap)
abundance_2024 <- fit_abundance_model(inputs_2024$abundance_inputs, pCap, lt_pCap_Us_2024)
rstan::summary(abundance_2024, pars = c("lt_pCap_U"))$summary |>
  data.frame() |>
  filter(Rhat > 1.05)

rstan::summary(abundance_2024, pars = c("N"))$summary |>
  data.frame() |>
  filter(Rhat > 1.05)

# lcc 2023
plot_juv_data("lcc", 2024)

# lcc 2022 - converges
plot_juv_data("lcc", 2022)
inputs_2022 <- prepare_inputs_pCap_abundance_STAN(SRJPEdata::weekly_juvenile_abundance_catch_data,
                                                  SRJPEdata::weekly_juvenile_abundance_efficiency_data,
                                                  site = "lcc", run_year = 2022,
                                                  effort_adjust = T,
                                                  #default_lgN_prior_denominator = 0.005 # uncomment the below if you want to set a specific denominator for the prior
)

lt_pCap_Us_2022 <- generate_lt_pCap_Us(inputs_2022$pCap_inputs, pCap)
abundance_2022 <- fit_abundance_model(inputs_2022$abundance_inputs, pCap, lt_pCap_Us_2022)
rstan::summary(abundance_2022, pars = "lt_pCap_U")$summary
rstan::traceplot(abundance_2022, pars = "lt_pCap_U")

# lcc 2020
plot_juv_data("lcc", 2020)

# lcc 2019
plot_juv_data("lcc", 2024)

# lcc 2014
plot_juv_data("lcc", 2014)
inputs_2014 <- prepare_inputs_pCap_abundance_STAN(SRJPEdata::weekly_juvenile_abundance_catch_data,
                                                  SRJPEdata::weekly_juvenile_abundance_efficiency_data,
                                                  site = "lcc", run_year = 2014,
                                                  effort_adjust = T,
                                                  #default_lgN_prior_denominator = 0.005 # uncomment the below if you want to set a specific denominator for the prior
)

lt_pCap_Us_2014 <- generate_lt_pCap_Us(inputs_2014$pCap_inputs, pCap)

abundance_2014 <- fit_abundance_model(inputs_2014$abundance_inputs, pCap, lt_pCap_Us_2014)
rstan::summary(abundance_2014, pars = "lt_pCap_U")$summary
# lcc 2024
plot_juv_data("lcc", 2024)

# lcc 2024
plot_juv_data("lcc", 2024)

# lcc 2024
plot_juv_data("lcc", 2024)

# lcc 2024
plot_juv_data("lcc", 2024)

# lcc 2024
plot_juv_data("lcc", 2024)

# lcc 2024
plot_juv_data("lcc", 2024)

# lcc 2024
plot_juv_data("lcc", 2024)

# lcc 2024
plot_juv_data("lcc", 2024)

# lcc 2024
plot_juv_data("lcc", 2024)



