# single stream dauphin model
# libraries ---------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(googleCloudStorageR)
library(rstan)
library(rstanarm)
library(bayesplot)
library(GGally) # pairs plot
library(waterYearType)
library(car) # vif
library(glmulti)
library(tidybayes)


# pull in data prep -------------------------------------------------------
load(here::here("data-raw", "adult_model", "adult_data_objects.Rdata"))
load(here::here("data-raw", "adult_model", "best_adult_models.Rdata"))


# code stan model for a single stream -------------------------------------

single_stream_dauphin <- "
  data {
    int N;
    int upstream_count[N];
    int spawner_count[N];
    real prespawn_survival[N];
    real environmental_index[N];
  }
  parameters {
    real <lower = 0> mu_k;
    real <lower = 0> sigma_k;
    real mu_a;
    real <lower = 0> sigma_a;
  }
  model {
    // priors
    mu_k ~ gamma(3,1);
    sigma_k ~ gamma(1,2);
    mu_a ~ uniform(0, 1000);
    sigma_a ~ gamma(0.001,0.001);

    real alpha[N];
    real beta;

    beta = mu_k * sigma_k;

    vector[N] lambda;

    // calibration between upstream adults and redd count
    for(i in 1:N) {
      alpha[i] = mu_k * beta * environmental_index[N];
      upstream_count[i] ~ lognormal(mu_a, sigma_a);
      prespawn_survival[i] ~ gamma(alpha[i], beta);
      lambda[i] = upstream_count[i] * prespawn_survival[i];
      spawner_count[i] ~ poisson(lambda[i]);
    }

  }"

# function to prepare data in list format to pass into stan model
prep_data_dauphin <- function(raw_data, stream_name, predictor_vars) {
  new_dat <- raw_data |>
    filter(stream == stream_name) |>
    drop_na(all_of(predictor_vars)) # can't have any NAs in predictor variables

  if(stream_name == "deer creek"){
    return(list(N = length(new_dat$year),
                spawner_count = new_dat$holding_count, # deer uses holding count
                upstream_count = new_dat$upstream_count,
                prespawn_survival = new_dat$prespawn_survival))
  } else {
    return(list(N = length(new_dat$year),
                spawner_count = new_dat$redd_count,
                upstream_count = new_dat$upstream_count,
                prespawn_survival = new_dat$prespawn_survival))
  }
}


# battle ------------------------------------------------------------------

best_adult_models$best_battle_lm
battle_coef <- unname(coef(best_adult_models$best_battle_lm))

# prep data list
battle_data_list <- prep_data_dauphin(adult_data_objects$survival_model_data_raw, "battle creek", "water_year_type")

# incorporate environmental index. Best fitting battle model only incorporates additive effect
# of a categorical variable (water year type), so we just shift the intercept of the mean
# based on water year type == dry or wet
battle_data_list$environmental_index <- adult_data_objects$survival_model_data_raw |>
  filter(stream == "battle creek") |>
  mutate(EI = ifelse(water_year_type == "dry", battle_coef[1],
                     (battle_coef[1] + battle_coef[2]))) |>
  pull(EI)

# run stan model
battle_dauphin <- stan(model_code = single_stream_dauphin,
                       data = battle_data_list,
                       chains = 4, iter = 5000*2, seed = 84735)

# diagnostics
mcmc_trace(battle_dauphin, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a"))
mcmc_areas(battle_dauphin, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a"))
mcmc_dens_overlay(battle_dauphin, pars = c("mu_k", "sigma_k",  "mu_a", "sigma_a")) # should be indistinguishable
neff_ratio(battle_dauphin, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should be >0.1
mcmc_acf(battle_dauphin, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should drop to be low
rhat(battle_dauphin, c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should be close to 1


# clear -------------------------------------------------------------------
best_adult_models$best_clear_lm
clear_coef <- unname(coef(best_adult_models$best_clear_lm))

# prep data list
clear_data_list <- prep_data_dauphin(adult_data_objects$survival_model_data_raw, "clear creek",
                                     c("gdd_total"))

# incorporate environmental index. Best fitting claer model incorporates additive effect
# of gdd_total
clear_data_list$environmental_index <- adult_data_objects$survival_model_data_raw |>
  filter(stream == "clear creek") |>
  mutate(EI = clear_coef[1] + (clear_coef[2] * gdd_total)) |>
  drop_na(EI) |>
  pull(EI)

# run stan model
clear_dauphin <- stan(model_code = single_stream_dauphin,
                      data = clear_data_list,
                      chains = 4, iter = 5000*2, seed = 84735)

# diagnostics
mcmc_trace(clear_dauphin, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a"))
mcmc_areas(clear_dauphin, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a"))
mcmc_dens_overlay(clear_dauphin, pars = c("mu_k", "sigma_k",  "mu_a", "sigma_a")) # should be indistinguishable
neff_ratio(clear_dauphin, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should be >0.1
mcmc_acf(clear_dauphin, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should drop to be low
rhat(clear_dauphin, c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should be close to 1


# mill --------------------------------------------------------------------
best_adult_models$best_mill_lm
mill_coef <- unname(coef(best_adult_models$best_mill_lm))

# prep data list
mill_dauphin_list <- prep_data_dauphin(adult_data_objects$survival_model_data_raw, "mill creek",
                                       c("max_flow", "gdd_total"))

# incorporate environmental index. Best fitting model incorporates additive effect of
# water year type and interactive effects of water year type * max_flow and max_flow * gdd_total
mill_dauphin_list$environmental_index <- adult_data_objects$survival_model_data_raw |>
  filter(stream == "mill creek") |>
  mutate(EI = mill_coef[1] + (mill_coef[2] * max_flow) + (mill_coef[3] * max_flow * gdd_total)) |>
  drop_na(EI) |>
  pull(EI)

# run stan model
mill_dauphin <- stan(model_code = single_stream_dauphin,
                     data = mill_dauphin_list,
                     chains = 4, iter = 5000*2, seed = 84735)

# diagnostics
mcmc_trace(mill_dauphin, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a"))
mcmc_areas(mill_dauphin, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a"))
mcmc_dens_overlay(mill_dauphin, pars = c("mu_k", "sigma_k",  "mu_a", "sigma_a")) # should be indistinguishable
neff_ratio(mill_dauphin, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should be >0.1
mcmc_acf(mill_dauphin, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should drop to be low
rhat(mill_dauphin, c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should be close to 1


# deer --------------------------------------------------------------------
best_adult_models$best_deer_lm
deer_coef <- unname(coef(best_adult_models$best_deer_lm))

# prep data list
deer_dauphin_list <- prep_data_dauphin(adult_data_objects$survival_model_data_raw, "deer creek",
                                       c("gdd_total", "max_flow"))

# incorporate environmental index. Best fitting model incorporates additive effect of water year type and
# interactive effects of gdd_total, max flow, and water year type
deer_dauphin_list$environmental_index <- adult_data_objects$survival_model_data_raw |>
  filter(stream == "deer creek") |>
  mutate(EI = deer_coef[1] + (deer_coef[2] * gdd_total * max_flow)) |>
  drop_na(EI) |>
  pull(EI)

deer_dauphin <- stan(model_code = single_stream_dauphin,
                     data = deer_dauphin_list,
                     chains = 4, iter = 5000*2, seed = 84735)

mcmc_trace(deer_dauphin, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a"))
mcmc_areas(deer_dauphin, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a"))
mcmc_dens_overlay(deer_dauphin, pars = c("mu_k", "sigma_k",  "mu_a", "sigma_a")) # should be indistinguishable
neff_ratio(deer_dauphin, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should be >0.1
mcmc_acf(deer_dauphin, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should drop to be low
rhat(deer_dauphin, c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should be close to 1


# save data objects -------------------------------------------------------
single_stream_dauphin_fits <- list("battle_dauphin" = battle_dauphin,
                                   "clear_dauphin" = clear_dauphin,
                                   "mill_dauphin" = mill_dauphin,
                                   "deer_dauphin" = deer_dauphin)
save(single_stream_dauphin_fits, file = here::here("data-raw", "adult_model", "single_stream_dauphin_fits.Rdata"))

# pull together into data frame -------------------------------------------
# TODO decide on best way to summarize
summarized_results <- tibble("stream" = c("battle creek", "clear creek", "mill creek", "deer creek"),
                             "mu_k" = c(median(extract(battle_dauphin, permuted = TRUE)$mu_k),
                                        median(extract(clear_dauphin, permuted = TRUE)$mu_k),
                                        median(extract(mill_dauphin, permuted = TRUE)$mu_k),
                                        median(extract(deer_dauphin, permuted = TRUE)$mu_k)),
                             "sigma_k" = c(median(extract(battle_dauphin, permuted = TRUE)$sigma_k),
                                           median(extract(clear_dauphin, permuted = TRUE)$sigma_k),
                                           median(extract(mill_dauphin, permuted = TRUE)$sigma_k),
                                           median(extract(deer_dauphin, permuted = TRUE)$sigma_k))) |>
  glimpse()
