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
                                             site = "lcc", run_year = 2002, effort_adjust = T)

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
diagnostic_plots_split("ubc", 2008, abundance)


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


