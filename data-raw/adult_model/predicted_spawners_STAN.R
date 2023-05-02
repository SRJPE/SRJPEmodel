# adult model - predicted redds

# libraries ---------------------------------------------------------------
library(tidyverse)
library(rstan)
library(rstanarm)
library(bayesplot)
library(tidybayes)


# pull in data prep -------------------------------------------------------
load(here::here("data-raw", "adult_model", "adult_data_objects.Rdata"))
load(here::here("data-raw", "adult_model", "best_adult_models.Rdata"))


# prep covariates ---------------------------------------------------------
full_data_for_input <- adult_data_objects$survival_model_data_raw |>
  mutate(wy_type_std = ifelse(water_year_type == "dry", 0, 1),
         max_flow_std = as.vector(scale(max_flow)),
         gdd_std = as.vector(scale(gdd_total)),
         min_passage_timing_std = as.vector(scale(min_passage_timing))) |>
  select(year, stream, upstream_count, redd_count, holding_count,
         wy_type_std, max_flow_std, gdd_std, min_passage_timing_std) |>
  glimpse()

# stan model --------------------------------------------------------------

predicted_redds <- "
  data {
    int N;
    int observed_passage[N];
    int observed_spawners[N]; // this is redd count but also holding count
    real environmental_covar[N];
  }
  parameters {
    real <lower = 0> mean_redds_per_spawner;
    real <lower = 0> tau_redds_per_spawner;
    real logit_redds_per_spawner[N];
    //real <lower = 0> redds_per_spawner[N];
    //real <lower = 0> predicted_redds[N];
    real b1_survival;

  } model {

    // priors
    mean_redds_per_spawner ~ normal(-0.69, 0.001);
    tau_redds_per_spawner ~ gamma(0.01, 0.01);
    b1_survival ~ normal(0, 0.001);

    // transform parameter
    real sigma_redds_per_spawner = pow(tau_redds_per_spawner, -0.5);

    // calculated values
    real survival_rate[N];
    real redds_per_spawner[N];
    real predicted_spawners[N];

    // calibration between upstream adults and redd count
    for(i in 1:N) {

      // predicted survival rate
      logit_redds_per_spawner[i] ~ normal(mean_redds_per_spawner, tau_redds_per_spawner);
      survival_rate[i] = inv_logit(logit_redds_per_spawner[i] + b1_survival * environmental_covar[i]);

      redds_per_spawner[i] = inv_logit(logit_redds_per_spawner[i]);

      // predicted redds is product of observed passage * survival rate
      predicted_spawners[i] = observed_passage[i] * survival_rate[i];

      // likelihood of predicted redds | observed redds
      observed_spawners[i] ~ poisson(predicted_spawners[i]);

    }

  } generated quantities {
    real survival_rate[N];
    real predicted_spawners[N];
    real redds_per_spawner[N];

    for(i in 1:N) {
      survival_rate[i] = inv_logit(logit_redds_per_spawner[i] + b1_survival * environmental_covar[i]);
      predicted_spawners[i] = observed_passage[i] * survival_rate[i];
      redds_per_spawner[i] = inv_logit(logit_redds_per_spawner[i]);
    }
}"

full_battle <- full_data_for_input |>
  filter(stream == "battle creek") |>
  drop_na(wy_type_std) |>
  glimpse()

battle_data_list <- list("N" = length(unique(full_battle$year)),
                         "observed_passage" = full_battle$upstream_count,
                         "observed_spawners" = full_battle$redd_count,
                         "environmental_covar" = full_battle$wy_type_std)

battle_pred_redd <- stan(model_code = predicted_redds,
                       data = battle_data_list,
                       chains = 3, iter = 20000*2, seed = 84735)

mcmc_trace(battle_pred_redd, pars = c("b1_survival"))
mcmc_areas(battle_pred_redd, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a"))
mcmc_dens_overlay(battle_pred_redd, pars = c("mu_k", "sigma_k",  "mu_a", "sigma_a")) # should be indistinguishable
neff_ratio(battle_pred_redd, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should be >0.1
mcmc_acf(battle_pred_redd, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should drop to be low
rhat(battle_pred_redd, c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should be close to 1

