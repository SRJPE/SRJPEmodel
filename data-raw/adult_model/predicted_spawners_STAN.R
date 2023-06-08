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

# download input data
gcs_get_object(object_name = "jpe-model-data/adult-model/adult_data_input_raw.csv",
               bucket = gcs_get_global_bucket(),
               saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
                                       "adult_data_input_raw.csv"),
               overwrite = TRUE)
# download covariate data
gcs_get_object(object_name = "jpe-model-data/adult-model/adult_model_covariates_standard.csv",
               bucket = gcs_get_global_bucket(),
               saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
                                       "adult_model_covariates_standard.csv"),
               overwrite = TRUE)

# read in datasets
adult_data_input_raw <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                            "adult_data_input_raw.csv"))

adult_model_covariates <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                              "adult_model_covariates_standard.csv"))

full_data_for_input <- full_join(adult_data_input_raw,
                                 adult_model_covariates,
                                 by = c("year", "stream")) |>
  filter(!is.na(data_type)) |>
  pivot_wider(id_cols = c(year, stream, wy_type, max_flow_std, gdd_std,
                          median_passage_timing_std),
              names_from = data_type,
              values_from = count) |>
  glimpse()

# function for prepping data for model

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

# function for extracting parameter estimates
get_all_pars <- function(model_fit, stream_name) {
  par_results <- summary(model_fit)$summary

  results_tibble <- as.data.frame(par_results) |>
    rownames_to_column("par_names") |>
    mutate(stream = stream_name)

  return(results_tibble)
}

# function for running model
run_passage_to_spawner_model <- function(data, stream_name, selected_covar, seed) {

  stream_data <- data |>
    mutate(observed_spawners = ifelse(stream == "deer creek", holding_count, redd_count)) |>
    filter(stream == stream_name) |>
    drop_na(all_of(selected_covar), observed_spawners, upstream_estimate)

  covar <- stream_data |>
    pull(all_of(selected_covar))

  stream_data_list <- list("N" = length(unique(stream_data$year)),
                           "observed_passage" = stream_data$upstream_estimate,
                           "observed_spawners" = stream_data$observed_spawners,
                           "environmental_covar" = covar,
                           "SStot" = calculate_ss_tot(stream_data))

  stream_model_fit <- stan(model_code = predicted_spawners,
                               data = stream_data_list,
                               chains = 3, iter = 20000*2, seed = seed)

  return(get_all_pars(stream_model_fit, stream_name))

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


# run STAN model ----------------------------------------------------------
battle_results <- run_passage_to_spawner_model(full_data_for_input,
                                               "battle creek",
                                               "wy_type",
                                               84735)

clear_results <- run_passage_to_spawner_model(full_data_for_input,
                                               "clear creek",
                                               "max_flow_std",
                                               84735)

mill_results <- run_passage_to_spawner_model(full_data_for_input,
                                               "mill creek",
                                               "gdd_std",
                                               84735)

deer_results <- run_passage_to_spawner_model(full_data_for_input,
                                               "deer creek",
                                               "wy_type",
                                               84735)


# write model summaries ---------------------------------------------------
model_fit_summaries <- bind_rows(battle_results, clear_results,
                                 mill_results, deer_results) |>
                                   filter(str_detect(par_names, "predicted_spawners") |
                                            par_names == "b1_survival") |>
  glimpse()

model_fit_diagnostics <- bind_rows(battle_results, clear_results,
                                   mill_results, deer_results) |>
  filter(!str_detect(par_names, "predicted_spawners")) |>
  filter(par_names != "b1_survival") |>
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
