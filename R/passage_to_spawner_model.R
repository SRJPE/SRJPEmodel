# adult model - predicted redds

# libraries ---------------------------------------------------------------
library(tidyverse)
library(rstan)
library(rstanarm)
library(bayesplot)
library(tidybayes)
library(SRJPEdata)

# helper functions for P2S ------------------------------------------------
#' @title Passage to Spawner Sum of Squares
#' @description This function calculates the total sum of squares, which is required
#' as input to the `data` call for `run_passage_to_spawner_model()`.
#' @export
#' @md
calculate_ss_tot <- function(data) {
  ss_tot <- 0
  data <- data |>
    mutate(observed_spawners = case_when(stream == "deer creek" ~ holding_count,
                                         stream == "butte creek" ~ carcass_estimate,
                                         TRUE ~ redd_count))


  N <- length(data$year)
  for(i in 1:N) {
    ss_tot <- ss_tot + (data$observed_spawners[i] - mean(data$observed_spawners))^2
  }
  return(ss_tot)
}

#' @title Extract Passage to Spawner Estimates
#' @description This function extracts parameter estimates and summary statistics for
#' the Passage to Spawner model (see `?run_passage_to_spawner_model`). The function
#' calls `summary()` on the `stanfit` object produced by running the model. For details,
#' see \html{https://mc-stan.org/rstan/reference/stanfit-method-summary.html}.
#' @export
#' @md
get_all_pars <- function(model_fit, stream_name) {
  par_results <- rstan::summary(model_fit)$summary

  results_tibble <- as.data.frame(par_results) |>
    rownames_to_column("par_names") |>
    mutate(stream = stream_name)

  return(results_tibble)
}


#' @title Run Passage to Spawner (P2S) STAN Model
#' @description This function takes in `SRJPEdata::observed_adult_input` and `SRJPEdata::adult_model_covariates_standard`
#' and runs the Passage to Spawner (P2S) model for a selected stream and environmental covariate. The function
#' will return a formatted list of parameter estimates and the full fitted `STAN` object.
#' TODO link to a vignette describing model and model comparisons?
#' @param observed_adult_input: The exported data object `SRJPEdata::observed_adult_input`. See `?SRJPEdata::observed_adult_input` for
#' details on data structure.
#' @param adult_model_covariates_standard: The exported data object `SRJPEdata::adult_model_covariates_standard`. See `?SRJPEdata::adult_model_covariates_standard` for
#' details on data structure.
#' @param stream_name: The name of the stream for which you would like to run the P2S model. Can be
#' `battle creek`, `clear creek`, `deer creek`, `mill creek`, or (draft form) `butte creek`.
#' @param selected_covariate: The environmental covariate you'd like to run the model for. Can be
#' either `wy_type` (water year type), `max_flow_std` (maximum flow), `gdd_std` (growing degree days),
#' `passage_index` (total upstream passage), `median_passage_timing_std` (median passage timing),
#' or `null_covar` (no environmental covariate).
#' See \code{vignette("prep_environmental_covariates.Rmd", package = "SRJPEdata")} for more details.
#' @returns A list containing `full_object`, the full fitted `stanfit` object (see [details](https://mc-stan.org/rstan/reference/stanfit-class.html)),
#' and `formatted_pars`, a formatted data table containing the following variables:
#' * **par_names** Parameter name
#' * **mean** Mean of the posterior distribution for a parameter
#' * **se_mean** Monte Carlo standard error for summary of all chains merged (see [details](https://mc-stan.org/rstan/reference/stanfit-method-summary.html))
#' * **sd** Mean of the posterior distribution for a parameter
#' * **`2.5%`** 2.5% quantile of posterior distribution for a parameter.
#' * **`25%`** 25% quantile of posterior distribution for a parameter.
#' * **`50%`** 50% quantile of posterior distribution for a parameter.
#' * **`75%`** 75% quantile of posterior distribution for a parameter.
#' * **`97.5%`** 97.5% quantile of posterior distribution for a parameter.
#' * **n_eff** Effective sample size for a parameter
#' * **Rhat** Split Rhats for a parameter
#' * **stream** Stream name for the parameter estimate
#' @export
#' @family passage_to_spawner
#' @md
run_passage_to_spawner_model <- function(observed_adult_input, adult_model_covariates,
                                         stream_name, selected_covariate) {

  stream_name <- tolower(stream_name)
  if(!stream_name %in% c("battle creek", "clear creek", "mill creek", "butte creek")){
    cli::cli_abort("Incorrect stream name. Please pass an approved stream.")
  }
  if(!selected_covariate %in% c("wy_type", "max_flow_std", "gdd_std", "passage_index", "median_passage_timing_std", "null_covar")){
    cli::cli_abort("Incorrect/Unavailable environmental covariate. Please pass an approved covariate name.")
  }

  # combine observed counts and covariates and pivot wider
  data <- full_join(observed_adult_input,
                                   adult_model_covariates,
                                   by = c("year", "stream")) |>
    filter(!is.na(data_type)) |>
    pivot_wider(id_cols = c(year, stream, wy_type, max_flow_std, gdd_std,
                            median_passage_timing_std, passage_index),
                names_from = data_type,
                values_from = count)

  # set percent female variable to 1 for carcass and holding surveys, 0.5 for redd
  percent_female <- case_when(stream_name %in% c("battle creek", "clear creek", "mill creek") ~ 0.5,
                              stream_name %in% c("yuba river", "feather river", "butte creek", "deer creek") ~ 1)

  # use holding count for deer, carcass for butte, and redd for other streams
  stream_data <- data |>
    mutate(observed_spawners = case_when(stream == "deer creek" ~ holding_count,
                                         stream == "butte creek" ~ carcass_estimate,
                                         TRUE ~ redd_count)) |>
    filter(upstream_estimate > 0,
           stream == stream_name) |>
    # null covariate
    mutate(null_covar = 0) |>
    drop_na(all_of(selected_covariate), observed_spawners, upstream_estimate)

  covar <- stream_data |>
    pull(all_of(selected_covariate))

  stream_data_list <- list("N" = length(unique(stream_data$year)),
                           "input_years" = unique(stream_data$year),
                           "observed_passage" = stream_data$upstream_estimate,
                           "observed_spawners" = stream_data$observed_spawners,
                           "percent_female" = percent_female,
                           "environmental_covar" = covar,
                           "ss_total" = calculate_ss_tot(stream_data),
                           "average_upstream_passage" = mean(stream_data$upstream_estimate, na.rm = TRUE))

  passage_to_spawner_STAN_code <- read_file("model-files/passage_to_spawner.txt")

  cli::cli_process_start("Fitting P2S STAN model")
  stream_model_fit <- rstan::stan(model_name = paste("passage_to_spawner", stream_name, selected_covariate, sep = "_"),
                                  model_code = passage_to_spawner_STAN_code,
                                  data = stream_data_list,
                                  chains = 3, iter = 20000*2, seed = 84735)
  cli::cli_process_done("P2S STAN model fitting complete")

  # check rhat
  rhat_diagnostic <- bayesplot::rhat(stream_model_fit) |> unname()
  rhat_diagnostic <- rhat_diagnostic[rhat_diagnostic != "NaN"]
  if(sum(rhat_diagnostic > 1.05) > 0){
    convergence_message <- paste0("not converged: rhat for ", length(rhat_diagnostic > 1.05), " parameter(s) over 1.05")
  } else {
    convergence_message <- "success - convergence by Rhat metric"
  }

  cli::cli_alert(convergence_message)
  return(list("full_object" = stream_model_fit,
              "formatted_pars" = get_all_pars(stream_model_fit, stream_name)))
}

