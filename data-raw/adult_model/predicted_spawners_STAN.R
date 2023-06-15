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

# stan model --------------------------------------------------------------

predicted_spawners <- "
  data {
    int N;
    int observed_passage[N];
    int observed_spawners[N]; // this is redd count but also holding count
    real environmental_covar[N];
    real ss_total;
    real average_upstream_passage;
  }

  parameters {
    real log_mean_redds_per_spawner; //no lower constraint of 0 since in log space
    real <lower = 0> sigma_redds_per_spawner;
    real b1_survival;
    real log_redds_per_spawner[N];
  }

  transformed parameters{

    real conversion_rate[N];
    real redds_per_spawner[N];
    real predicted_spawners[N];
    real mean_redds_per_spawner = exp(log_mean_redds_per_spawner);

    for(i in 1:N) {
      conversion_rate[i] = exp(log_redds_per_spawner[i] + b1_survival * environmental_covar[i]);

      // predicted redds is product of observed passage * conversion rate
      predicted_spawners[i] = observed_passage[i] * conversion_rate[i];
      redds_per_spawner[i] = exp(log_redds_per_spawner[i]); //for ouptut only
    }
  }

  model {
    log_redds_per_spawner ~ normal(log_mean_redds_per_spawner, sigma_redds_per_spawner);
    observed_spawners ~ poisson(predicted_spawners);
  }

  generated quantities {

    // calculate R2 fixed effects
    real y_pred[N];
    real RE_draws[N]; // vector for storing random draws for RE
    real var_fit = 0;
    real var_res = 0.0;
    real ss_res = 0.0;

    for(i in 1:N) {
      RE_draws[i] = normal_rng(log_mean_redds_per_spawner, sigma_redds_per_spawner);
      y_pred[i] = observed_passage[i] * exp(RE_draws[i]);
    }

    real y_pred_mu = mean(y_pred[]);
    real R2_data;
    real R2_fixed;

    // calculate R2_data and R2_fixed
    for(i in 1:N) {
      var_fit = var_fit + (y_pred[i] - y_pred_mu)^2;
      // TODO var_res is actually based off y_pred with full env effects?
      // TODO and then subtract y_pred from that?
      var_res = var_res + (observed_spawners[i] - y_pred[i])^2;
      ss_res = ss_res + (observed_spawners[i] - predicted_spawners[i])^2;
    }

    R2_data = 1.0 - (ss_res / ss_total);
    R2_fixed = 1.0 - (var_fit / (var_fit + var_res));

    // forecast error in spawner abundance at average upstream passage abundance
    real spawner_abundance_forecast[2];
    real rps_forecast = normal_rng(log_mean_redds_per_spawner, sigma_redds_per_spawner);// random draw from hyper

    spawner_abundance_forecast[1] = average_upstream_passage * exp(rps_forecast); //for dummy variable = 0 condition, or at average covariate conditions
    spawner_abundance_forecast[2] = average_upstream_passage * exp(rps_forecast + b1_survival); //for dummy variable=1 condition, or at 1 SD from mean covariate condition
  }
"

# test covariates with bayesian model -------------------------------------
compare_covars <- function(data, stream_name, seed, truncate_data) {
  final_results <- tibble("par_names" = "DELETE",
                          "stream" = "DELETE",
                          "mean" = 0.0,
                          "sd" = 0.0,
                          "covar_considered" = "DELETE")

  covars_to_consider = c("wy_type", "max_flow_std", "gdd_std",
                         "null_covar") # no more "median_passage_timing_std" bc of sample size

  # loop through covariates to test
  for(i in 1:length(covars_to_consider)) {

    selected_covar <- covars_to_consider[i]

    if(truncate_data == TRUE) {
      # drop NAs for all covariates being considered
      # makes sure you are using the same dataset for each cycle
      stream_data <- data |>
        mutate(observed_spawners = ifelse(stream == "deer creek", holding_count, redd_count)) |>
        filter(upstream_estimate > 0) |>
        filter(stream == stream_name) |>
        drop_na(wy_type, max_flow_std, gdd_std, observed_spawners, upstream_estimate) |>
        mutate(null_covar = 0)
    } else if(truncate_data == FALSE) {
      stream_data <- data |>
        mutate(observed_spawners = ifelse(stream == "deer creek", holding_count, redd_count)) |>
        filter(upstream_estimate > 0) |>
        filter(stream == stream_name) |>
        mutate(null_covar = 0) |>
        drop_na(observed_spawners, upstream_estimate, all_of(selected_covar))
    }

    print(selected_covar)
    print(unique(stream_data$stream))

    covar <- stream_data |>
      pull(all_of(selected_covar))

    stream_data_list <- list("N" = length(unique(stream_data$year)),
                             "observed_passage" = stream_data$upstream_estimate,
                             "observed_spawners" = stream_data$observed_spawners,
                             "environmental_covar" = covar,
                             "ss_total" = calculate_ss_tot(stream_data),
                             "average_upstream_passage" = mean(stream_data$upstream_estimate,
                                                               na.rm = TRUE))

    print(stream_data_list)

    stream_model_fit <- stan(model_code = predicted_spawners,
                             data = stream_data_list,
                             chains = 3, iter = 20000*2, seed = seed)

    stream_results <- get_all_pars(stream_model_fit, stream_name) |>
      filter(par_names %in% c("R2_data", "R2_fixed", "mean_redds_per_spawner",
                              "b1_survival") |
               str_detect(par_names, "spawner_abundance_forecast")) |> # TODO add pspawn after it's added to STAN code
      select(par_names, stream, mean, sd) |>
      mutate(covar_considered = selected_covar)

    final_results <- bind_rows(final_results, stream_results)
  }

  final_results <- final_results |>
    filter(par_names != "DELETE")

  return(final_results)
}

all_streams_truncated <- tibble(par_names = "DELETE",
                                stream = "DELETE",
                                mean = 0.0,
                                sd = 0.0,
                                covar_considered = "DELETE")

all_streams_full <- tibble(par_names = "DELETE",
                           stream = "DELETE",
                           mean = 0.0,
                           sd = 0.0,
                           covar_considered = "DELETE")


#for(i in c("battle creek", "clear creek", "deer creek", "mill creek")) {
for(i in c("battle creek", "deer creek", "mill creek")) {
  truncated_stream_results <- compare_covars(full_data_for_input, i, 84735, TRUE)
  #full_stream_results <- compare_covars(full_data_for_input, i, 84735, FALSE)
  all_streams_truncated <- bind_rows(all_streams_truncated, truncated_stream_results)
  #all_streams_full <- bind_rows(all_streams_full, full_stream_results)
}

# identify best covariate for each trib -----------------------------------
# The best covariate will have the
# highest R2_fixed, the
# largest absolute value of b1_surv, the
# most precise b1_surv, and the
# most precise forecast error.

get_rankings <- function(results, stream_name) {
  results <- results |>
    filter(stream == stream_name) |>
    # b1_survival and spawner_forecast[2] not valid for null model
    mutate(mean = ifelse(covar_considered == "null_covar" &
                           par_names %in% c("b1_survival", "spawner_abundance_forecast[2]"),
                         NA, mean),
           sd = ifelse(covar_considered == "null_covar" &
                           par_names %in% c("b1_survival", "spawner_abundance_forecast[2]"),
                         NA, sd)) |>
    mutate(mean = case_when(covar_considered == "null_covar" &
                              par_names %in% c("b1_survival", "spawner_abundance_forecast[2]") ~ NA,
                            covar_considered != "null_covar" &
                              par_names == "spawner_abundance_forecast[1]" ~ NA,
                            TRUE ~ mean),
           sd = case_when(covar_considered == "null_covar" &
                              par_names %in% c("b1_survival", "spawner_abundance_forecast[2]") ~ NA,
                            covar_considered != "null_covar" &
                              par_names == "spawner_abundance_forecast[1]" ~ NA,
                          TRUE ~ sd)) |>
    drop_na(mean, sd) |>
    mutate(par_names = case_when(par_names == "spawner_abundance_forecast[1]" ~ "spawner_abundance_forecast",
                                 par_names == "spawner_abundance_forecast[2]" ~ "spawner_abundance_forecast",
                                 TRUE ~ par_names))

  # assign rankings
  max_rankings <- results |>
    filter(par_names %in% c("b1_survival")) |>
    #filter(par_names %in% c("R2_fixed", "b1_survival")) |>
    mutate(mean = ifelse(par_names == "R2_fixed", mean, abs(mean))) |>
    group_by(par_names) |>
    mutate(rank = order(order(mean, decreasing = TRUE))) |>
    arrange(rank)

  min_rankings <- results |>
    filter(par_names %in% c("b1_survival", "spawner_abundance_forecast")) |>
    filter(mean != Inf) |>
    group_by(par_names) |>
    mutate(rank = order(order(sd, decreasing = FALSE))) |>
    arrange(rank)

  results_with_rankings <- bind_rows(max_rankings, min_rankings)
  # results_with_rankings_wide <- results_with_rankings |>
  #   pivot_wider(id_cols = covar_considered,
  #               names_from = par_names,
  #               values_from = rank)
  return(results_with_rankings)

  best_covar <- results_with_rankings |>
    group_by(covar_considered) |>
    summarise(sum_rank = sum(rank),
              mean_rank = mean(rank)) |>
    #slice(which.min(sum_rank)) |>
    slice(which.min(mean_rank)) |>
    pull(covar_considered)

  return(list("best_model" = best_covar,
              "best_model_values" = results_with_rankings |>
                filter(covar_considered == best_covar) |>
                select(-c(stream, covar_considered))))

  # return(list("highest_R2_fixed" = results |>
  #               filter(par_names == "R2_fixed") |>
  #               slice(which.max(mean)) |>
  #               select(covar_considered, mean),
  #             "highest_absv_b1_surv" = results |>
  #               filter(par_names == "b1_survival") |>
  #               mutate(mean = abs(mean)) |>
  #                      slice(which.max(mean)) |>
  #               select(covar_considered, mean),
  #             "b1_surv_precise" = results |>
  #               filter(par_names == "b1_survival") |>
  #               slice(which.min(sd)) |>
  #               select(covar_considered, sd),
  #             "forecast_error_precise_1" = results |>
  #               filter(par_names == "spawner_abundance_forecast[1]") |>
  #               slice(which.min(sd)) |>
  #               select(covar_considered, sd),
  #            "forecast_error_precise_2" = results |>
  #              filter(par_names == "spawner_abundance_forecast[2]") |>
  #              slice(which.min(sd)) |>
  #              select(covar_considered, sd)))
}

all_streams |>
  filter(stream == "clear creek") |>
  arrange(par_names, mean) |> print(n=Inf)

all_streams_truncated |> get_rankings("battle creek")# water year type
all_streams_truncated |> get_rankings("clear creek") # null model
all_streams_truncated |> get_rankings("deer creek") # water year type
all_streams_truncated |> get_rankings("mill creek") # null model

# TODO assign each covariate considered a ranking from 1 (best)-5 based on how they
# perform on the diagnostics
# then sum rankings, pull lowest
# compare to "1-D" analysis



# run STAN model with covars identified by previous analysis -------------------------
# function for running model
run_passage_to_spawner_model <- function(data, stream_name, selected_covar, seed) {

  stream_data <- data |>
    filter(upstream_estimate > 0) |>
    mutate(observed_spawners = ifelse(stream == "deer creek", holding_count, redd_count)) |>
    filter(stream == stream_name) |>
    # null covariate
    mutate(null_covar = 0) |>
    drop_na(all_of(selected_covar), observed_spawners, upstream_estimate)

  covar <- stream_data |>
    pull(all_of(selected_covar))

  stream_data_list <- list("N" = length(unique(stream_data$year)),
                           "observed_passage" = stream_data$upstream_estimate,
                           "observed_spawners" = stream_data$observed_spawners,
                           "environmental_covar" = covar,
                           "ss_total" = calculate_ss_tot(stream_data),
                           "average_upstream_passage" = mean(stream_data$upstream_estimate, na.rm = TRUE))

  stream_model_fit <- stan(model_code = predicted_spawners,
                           data = stream_data_list,
                           chains = 3, iter = 20000*2, seed = seed)

  return(get_all_pars(stream_model_fit, stream_name))

}

battle_results <- run_passage_to_spawner_model(full_data_for_input,
                                               "battle creek",
                                               "wy_type",
                                               84735)

clear_results <- run_passage_to_spawner_model(full_data_for_input,
                                               "clear creek",
                                               "wy_type",
                                               84735)

deer_results <- run_passage_to_spawner_model(full_data_for_input,
                                             "deer creek",
                                             "wy_type",
                                             84735)

mill_results <- run_passage_to_spawner_model(full_data_for_input,
                                               "mill creek",
                                               "null_covar",
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



# scratch to fix R2_fixed -------------------------------------------------

test_1 <- run_passage_to_spawner_model(full_data_for_input,
                                              "clear creek",
                                              "wy_type",
                                              84735)
test_2 <- run_passage_to_spawner_model(full_data_for_input,
                                       "clear creek",
                                       "max_flow_std",
                                       84735)
test_3 <- run_passage_to_spawner_model(full_data_for_input,
                                       "clear creek",
                                       "gdd_std",
                                       84735)
test_4 <- run_passage_to_spawner_model(full_data_for_input,
                                       "clear creek",
                                       "null_covar",
                                       84735)

test <- bind_rows(test_1 |>
                    mutate(covar_considered = "wy_type") |>
                    select(par_names, mean, sd, covar_considered),
                  test_2 |>
                    mutate(covar_considered = "max_flow_std") |>
                    select(par_names, mean, sd, covar_considered),
                  test_3 |>
                    mutate(covar_considered = "gdd_std") |>
                    select(par_names, mean, sd, covar_considered),
                  test_4 |>
                    mutate(covar_considered = "null_covar") |>
                    select(par_names, mean, sd, covar_considered))


