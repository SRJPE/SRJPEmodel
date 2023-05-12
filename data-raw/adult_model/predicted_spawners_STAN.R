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
         wy_type_sd = as.vector(scale(wy_type_std)), # TODO double check this. was producing negative b1_surv estimates for wet year types
         max_flow_std = as.vector(scale(max_flow)),
         gdd_std = -as.vector(scale(gdd_total)),
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
    real log_redds_per_spawner[N];
    //real <lower = 0> redds_per_spawner[N];
    //real <lower = 0> predicted_redds[N];
    real b1_survival;

  } model {

    // priors
    //mean_redds_per_spawner ~ normal(-0.69, 0.001);
    //tau_redds_per_spawner ~ gamma(0.01, 0.01);
    //b1_survival ~ normal(0, 0.001);

    // transform parameter
    real sigma_redds_per_spawner = pow(tau_redds_per_spawner, -0.5);

    // calculated values
    real survival_rate[N];
    real redds_per_spawner[N];
    real predicted_spawners[N];

    // calibration between upstream adults and redd count
    for(i in 1:N) {

      // predicted survival rate
      log_redds_per_spawner[i] ~ normal(mean_redds_per_spawner, tau_redds_per_spawner);
      survival_rate[i] = inv_logit(log_redds_per_spawner[i] + b1_survival * environmental_covar[i]);

      redds_per_spawner[i] = exp(log_redds_per_spawner[i]);

      // predicted redds is product of observed passage * survival rate
      predicted_spawners[i] = observed_passage[i] * survival_rate[i];

      // likelihood of observed redds | predicted redds
      observed_spawners[i] ~ poisson(predicted_spawners[i]);

    }

  } generated quantities {
    real survival_rate[N];
    real predicted_spawners[N];
    real redds_per_spawner[N];

    for(i in 1:N) {
      survival_rate[i] = inv_logit(log_redds_per_spawner[i] + b1_survival * environmental_covar[i]);
      predicted_spawners[i] = observed_passage[i] * survival_rate[i];
      redds_per_spawner[i] = exp(log_redds_per_spawner[i]);
    }
}"


# battle ------------------------------------------------------------------
full_battle <- full_data_for_input |>
  filter(stream == "battle creek") |>
  drop_na(wy_type_std) |>
  glimpse()

battle_data_list <- list("N" = length(unique(full_battle$year)),
                         "observed_passage" = full_battle$upstream_count,
                         "observed_spawners" = full_battle$redd_count,
                         "environmental_covar" = as.vector(scale(full_battle$wy_type_std)))

battle_pred_redd <- stan(model_code = predicted_redds,
                       data = battle_data_list,
                       chains = 3, iter = 20000*2, seed = 84735)


# clear -------------------------------------------------------------------
full_clear <- full_data_for_input |>
  filter(stream == "clear creek") |>
  drop_na(max_flow_std) |>
  glimpse()

clear_data_list <- list("N" = length(unique(full_clear$year)),
                         "observed_passage" = full_clear$upstream_count,
                         "observed_spawners" = full_clear$redd_count,
                         "environmental_covar" = as.vector(scale(full_clear$max_flow_std)))

clear_pred_redd <- stan(model_code = predicted_redds,
                         data = clear_data_list,
                         chains = 3, iter = 20000*2, seed = 84735)


# mill --------------------------------------------------------------------
full_mill <- full_data_for_input |>
  filter(stream == "mill creek") |>
  drop_na(gdd_std) |>
  glimpse()

mill_data_list <- list("N" = length(unique(full_mill$year)),
                         "observed_passage" = full_mill$upstream_count,
                         "observed_spawners" = full_mill$redd_count,
                         "environmental_covar" = as.vector(scale(full_mill$gdd_std)))

mill_pred_redd <- stan(model_code = predicted_redds,
                         data = mill_data_list,
                         chains = 3, iter = 20000*2, seed = 84735)


# deer --------------------------------------------------------------------
full_deer <- full_data_for_input |>
  filter(stream == "deer creek") |>
  drop_na(wy_type_std) |>
  glimpse()

deer_data_list <- list("N" = length(unique(full_deer$year)),
                         "observed_passage" = full_deer$upstream_count,
                         "observed_spawners" = full_deer$holding_count,
                         "environmental_covar" = as.vector(scale(full_deer$wy_type_std)))

deer_pred_redd <- stan(model_code = predicted_redds,
                         data = deer_data_list,
                         chains = 3, iter = 20000*2, seed = 84735)



# get results -------------------------------------------------------------
# define function to pull out predicted spawners
get_predicted_spawners <- function(model_fit) {
  par_results <- summary(model_fit)$summary
  results_tibble <- as.data.frame(par_results) |>
    rownames_to_column("par_names") |>
    filter(str_detect(par_names, "predicted_spawners"))

  return(results_tibble)
}

# save as object ----------------------------------------------------------
model_fit_summaries <- list("battle_pred" = get_predicted_spawners(battle_pred_redd),
                            "clear_pred" = get_predicted_spawners(clear_pred_redd),
                            "mill_pred" = get_predicted_spawners(mill_pred_redd),
                            "deer_pred" = get_predicted_spawners(deer_pred_redd))

save(model_fit_summaries, file = here::here("data-raw", "adult_model",
                                       "model_fit_summaries.Rdata"))

pred_redd_fits <- list("battle_pred_redd" = battle_pred_redd,
                       "clear_pred_redd" = clear_pred_redd,
                       "mill_pred_redd" = mill_pred_redd,
                       "deer_pred_redd" = deer_pred_redd)

save(pred_redd_fits, file = here::here("data-raw", "adult_model",
                                       "pred_redd_fits.Rdata"))


# function to get par name ------------------------------------------------
get_par_names <- function(data_list) {
  N <- data_list$N

  names <- c()
  for(i in 1:N) {
    names[i] <- paste0("predicted_spawners[", i, "]")
  }
  #names[N+1] <- "b1_survival"
  return(names)
}

get_par_names(battle_data_list)

# diagnostics -------------------------------------------------------------

summary(battle_pred_redd)$summary[, 1:4]

summary(battle_pred_redd)$summary[str_detect(rownames(summary(battle_pred_redd)$summary), "predicted_spawners"), 1:4] |> round(4)
summary(battle_pred_redd)$summary[rownames(summary(battle_pred_redd)$summary) == "b1_survival", 1:10] |> round(4)

mcmc_trace(battle_pred_redd, pars = get_par_names(battle_data_list))
mcmc_dens_overlay(battle_pred_redd, pars = get_par_names(battle_data_list))
mcmc_dens_overlay(battle_pred_redd, pars = "b1_survival")
mcmc_areas(battle_pred_redd, pars = get_par_names(battle_data_list))
neff_ratio(battle_pred_redd, pars = get_par_names(battle_data_list)) # should be >0.1
mcmc_acf(battle_pred_redd, pars = get_par_names(battle_data_list)) # should drop to be low
rhat(battle_pred_redd, get_par_names(battle_data_list)) # should be close to 1

