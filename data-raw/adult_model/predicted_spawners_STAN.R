# adult model - predicted redds

# libraries ---------------------------------------------------------------
library(tidyverse)
library(rstan)
library(rstanarm)
library(bayesplot)
library(tidybayes)
library(googleCloudStorageR)

# pull in data from google cloud ------------------------------------------
gcs_auth(json_file = Sys.getenv("GCS_AUTH_FILE"))
# Set global bucket
gcs_global_bucket(bucket = Sys.getenv("GCS_DEFAULT_BUCKET"))
# save data as csv
survival_model_data_raw <- gcs_get_object(object_name = "jpe-model-data/adult-model/survival_model_data_raw.csv",
                                          bucket = gcs_get_global_bucket(),
                                          saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
                                                                  "survival_model_data_raw.csv"),
                                          overwrite = TRUE)
# read csvs
survival_model_data_raw <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                               "survival_model_data_raw.csv"))


# prep covariates ---------------------------------------------------------
full_data_for_input <- survival_model_data_raw |>
  mutate(wy_type_std = ifelse(water_year_type == "dry", 0, 1),
         wy_type_std = as.vector(scale(wy_type_std)), # TODO double check this. was producing negative b1_surv estimates for wet year types
         max_flow_std = as.vector(scale(max_flow)),
         gdd_std = as.vector(scale(gdd_total)),
         min_passage_timing_std = as.vector(scale(min_passage_timing))) |>
  select(year, stream, upstream_count, redd_count, holding_count,
         wy_type_std, max_flow_std, gdd_std, min_passage_timing_std) |>
  glimpse()

calculate_ss_tot <- function(data) {
  ss_tot <- 0
  data <- data |>
    mutate(observed_spawners = ifelse(stream == "deer creek", holding_count, redd_count))

  N <- length(data$year)
  for(i in 1:N) {
    ss_tot <- ss_tot + (data$observed_spawners[i] - mean(data$observed_spawners))^2
  }
  return(ss_tot)
}

# stan model --------------------------------------------------------------

predicted_spawners <- "
  data {
    int N;
    int observed_passage[N];
    int observed_spawners[N]; // this is redd count but also holding count
    real environmental_covar[N];
    real SStot;
  }

  parameters {
    real log_mean_redds_per_spawner; //no lower constraint of 0 since in log space
    real <lower = 0> sigma_redds_per_spawner;
    real b1_survival;
    real log_redds_per_spawner[N];
  }

  transformed parameters {
    real survival_rate[N];
    real redds_per_spawner[N];
    real predicted_spawners[N];

    real mean_redds_per_spawner=exp(log_mean_redds_per_spawner);

    real devre[N];// for R2_fixed computation (amount of variation in survival explained by fixed effect)

    for(i in 1:N) {

      survival_rate[i] = exp(log_redds_per_spawner[i] + b1_survival * environmental_covar[i]);
      redds_per_spawner[i] = exp(log_redds_per_spawner[i]);

      // predicted redds is product of observed passage * survival rate
      predicted_spawners[i] = observed_passage[i] * survival_rate[i];

      devre[i] = exp(log_redds_per_spawner[i]);
    }
  }

  model {
    log_redds_per_spawner ~ normal(log_mean_redds_per_spawner, sigma_redds_per_spawner);
    observed_spawners ~ poisson(predicted_spawners);
  }

  generated quantities {

    // model diagnostics
    real musurv = mean(survival_rate[]);
    real mudevre = mean(devre[]);

    real R2_data;
    real R2_fixed;
    real nom=0;
    real denom = 0;
    real SSres = 0.0;

    for(i in 1:N){
       nom = nom + (devre[i] - mudevre)^2; //sums of squares on variation in survival rate explained by redds_per_spawner
       denom = denom + (survival_rate[i] - musurv)^2; //sums of squares on total variation in survival (due to variation redds_per_spawner and due to covariate effectd)
       SSres = SSres + (observed_spawners[i] - predicted_spawners[i])^2; //residual variation for pearson r2 calculation (R2_data) to compare with R2_fixed
    }

    R2_fixed = 1.0 - nom / denom;
    R2_data = 1.0 - SSres / SStot;
  }"


# battle ------------------------------------------------------------------
full_battle <- full_data_for_input |>
  filter(stream == "battle creek") |>
  drop_na(wy_type_std) |>
  glimpse()

battle_data_list <- list("N" = length(unique(full_battle$year)),
                         "observed_passage" = full_battle$upstream_count,
                         "observed_spawners" = full_battle$redd_count,
                         "environmental_covar" = as.vector(scale(full_battle$wy_type_std)),
                         "SStot" = calculate_ss_tot(full_battle))

battle_pred_spawners <- stan(model_code = predicted_spawners,
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
                        "environmental_covar" = as.vector(scale(full_clear$max_flow_std)),
                        "SStot" = calculate_ss_tot(full_clear))

clear_pred_spawners <- stan(model_code = predicted_spawners,
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
                       "environmental_covar" = as.vector(scale(full_mill$gdd_std)),
                       "SStot" = calculate_ss_tot(full_mill))

mill_pred_spawners <- stan(model_code = predicted_spawners,
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
                       "environmental_covar" = as.vector(scale(full_deer$wy_type_std)),
                       "SStot" = calculate_ss_tot(full_deer))

deer_pred_spawners <- stan(model_code = predicted_spawners,
                         data = deer_data_list,
                         chains = 3, iter = 20000*2, seed = 84735)



# get results -------------------------------------------------------------
# define function to pull out predicted spawners
get_report_pars <- function(model_fit, stream_name) {
  par_results <- summary(model_fit)$summary
  results_tibble <- as.data.frame(par_results) |>
    rownames_to_column("par_names") |>
    filter(str_detect(par_names, "predicted_spawners"))

  b1_survival_val <- as.data.frame(par_results) |>
    rownames_to_column("par_names") |>
    filter(par_names == "b1_survival")

  results_tibble <- bind_rows(results_tibble, b1_survival_val) |>
    mutate(stream = stream_name)

  return(results_tibble)
}

# define function to pull out diagnostics pars
get_diagnostic_pars <- function(model_fit, stream_name) {
  par_results <- summary(model_fit)$summary
  results_tibble <- as.data.frame(par_results) |>
    rownames_to_column("par_names") |>
    filter(!str_detect(par_names, "predicted_spawners")) |>
    mutate(stream = stream_name)

  return(results_tibble)
}


# write model summaries ---------------------------------------------------
model_fit_summaries <- bind_rows(get_report_pars(battle_pred_spawners, "battle creek"),
                                 get_report_pars(clear_pred_spawners, "clear creek"),
                                 get_report_pars(mill_pred_spawners, "mill creek"),
                                 get_report_pars(deer_pred_spawners, "deer creek")) |>
  glimpse()

model_fit_diagnostics <- bind_rows(get_diagnostic_pars(battle_pred_spawners, "battle creek"),
                                   get_diagnostic_pars(clear_pred_spawners, "clear creek"),
                                   get_diagnostic_pars(mill_pred_spawners, "mill creek"),
                                   get_diagnostic_pars(deer_pred_spawners, "deer creek")) |>
  glimpse()

# save to google cloud ----------------------------------------------------
f <- function(input, output) write_csv(input, file = output)

gcs_upload(model_fit_summaries,
           object_function = f,
           type = "csv",
           name = "jpe-model-data/adult-model/model_fit_summaries.csv")

gcs_upload(model_fit_diagnostics,
           object_function = f,
           type = "csv",
           name = "jpe-model-data/adult-model/model_fit_diagnostic_pars.csv")

# # save as object ----------------------------------------------------------
# model_fit_summaries <- list("battle_pred" = get_pars_of_interest(battle_pred_redd, "battle creek"),
#                             "clear_pred" = get_pars_of_interest(clear_pred_redd),
#                             "mill_pred" = get_pars_of_interest(mill_pred_redd),
#                             "deer_pred" = get_pars_of_interest(deer_pred_redd))
#
# save(model_fit_summaries, file = here::here("data-raw", "adult_model",
#                                        "model_fit_summaries.Rdata"))



# scratch ------------------------------------------------
# pred_redd_fits <- list("battle_pred_redd" = battle_pred_redd,
#                        "clear_pred_redd" = clear_pred_redd,
#                        "mill_pred_redd" = mill_pred_redd,
#                        "deer_pred_redd" = deer_pred_redd)
#
# save(pred_redd_fits, file = here::here("data-raw", "adult_model",
#                                        "pred_redd_fits.Rdata"))

# get_par_names <- function(data_list) {
#   N <- data_list$N
#
#   names <- c()
#   for(i in 1:N) {
#     names[i] <- paste0("predicted_spawners[", i, "]")
#   }
#   #names[N+1] <- "b1_survival"
#   return(names)
# }
#
# get_par_names(battle_data_list)

# diagnostics -------------------------------------------------------------

# summary(battle_pred_redd)$summary[, 1:4]
#
# summary(battle_pred_redd)$summary[str_detect(rownames(summary(battle_pred_redd)$summary), "predicted_spawners"), 1:4] |> round(4)
# summary(battle_pred_redd)$summary[rownames(summary(battle_pred_redd)$summary) == "b1_survival", 1:10] |> round(4)
#
# mcmc_trace(battle_pred_redd, pars = get_par_names(battle_data_list))
# mcmc_dens_overlay(battle_pred_redd, pars = get_par_names(battle_data_list))
# mcmc_dens_overlay(battle_pred_redd, pars = "b1_survival")
# mcmc_areas(battle_pred_redd, pars = get_par_names(battle_data_list))
# neff_ratio(battle_pred_redd, pars = get_par_names(battle_data_list)) # should be >0.1
# mcmc_acf(battle_pred_redd, pars = get_par_names(battle_data_list)) # should drop to be low
# rhat(battle_pred_redd, get_par_names(battle_data_list)) # should be close to 1

