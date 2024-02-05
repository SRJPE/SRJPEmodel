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

#' @title Extract Years Passage to Spawner Estimates
#' @description This function extracts the years the data were collected from the results of
#' the Passage to Spawner model (see `?run_passage_to_spawner_model`).
#' @export
#' @md
get_years_from_P2S_model_fits <- function(P2S_model_fits) {
  years <- P2S_model_fits |>
    filter(stringr::str_detect(par_names, "years")) |>
    mutate(year_index = readr::parse_number(par_names),
           mean = round(mean)) |>
    select(year_index, year = mean, stream)

  years
}

#' @title Extract Passage to Spawner Estimates
#' @description This function extracts parameter estimates and summary statistics for
#' the Passage to Spawner model (see `?run_passage_to_spawner_model`). The function
#' calls `summary()` on the `stanfit` object produced by running the model. For details,
#' see [details](https://mc-stan.org/rstan/reference/stanfit-method-summary.html).
#' @export
#' @md
get_all_pars <- function(model_fit, stream_name) {
  par_results <- rstan::summary(model_fit)$summary

  results_tibble <- as.data.frame(par_results) |>
    rownames_to_column("par_names") |>
    mutate(stream = stream_name)

  years_modeled <- suppressWarnings(get_years_from_P2S_model_fits(results_tibble))

  results_tibble_with_years <- results_tibble |>
    mutate(year_index = suppressWarnings(readr::parse_number(par_names))) |>
    left_join(years_modeled, by = c("year_index", "stream"))

  return(results_tibble_with_years)
}

#' @title Extract Predicted Spawners from P2S Model
#' @description This function is called within `run_passage_to_spawner_model()` and pulls only the predicted spawner
#' counts, 2.5% and 97.5% confidence intervals from the `stanfit` object.
#' @export
#' @md
get_predicted_spawners_from_P2S <- function(results_tibble_with_years) {
  predicted_spawners <- results_tibble_with_years |>
    filter(stringr::str_detect(par_names, "predicted_spawners")) |>
    select(stream, year, median_predicted_spawners = `50%`,
           lcl = `2.5%`, ucl = `97.5%`) |>
    mutate(data_type = "P2S_predicted_spawners")

  predicted_spawners
}

#' @title Run Passage to Spawner (P2S) STAN Model
#' @description This function takes in `SRJPEdata::observed_adult_input` and `SRJPEdata::adult_model_covariates_standard`
#' and runs the Passage to Spawner (P2S) model for a selected stream and environmental covariate. The function
#' will return a formatted list of parameter estimates and the full fitted `stanfit` object.
#' See \code{vignette("passage_to_spawner_submodel.Rmd", package = "SRJPEmodel")} for more details.
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
#' @param extract_predicted_spawners: If `TRUE`, will only produce a tibble with median predicted spawners with upper (97.5%)
#' and lower (2.5%) confidence intervals. If `FALSE`, returns a list with all parameter estimates and the full `stanfit` object.
#' See \code{vignette("prep_environmental_covariates.Rmd", package = "SRJPEdata")} for more details.
#' @returns If `extracted_predicted_spawners == TRUE`, Returns a tibble containing the following variables:
#' * **stream** Stream name
#' * **year** Year observed data came from
#' * **median_predicted_spawners** Median predicted spawner value for a given stream and year
#' * **lcl** 2.5% confidence interval for predicted spawner value
#' * **ucl** 97.5% confidence interval for predicted spawner value
#' * **data_type** Data type: `P2S_predicted_spawners`
#'
#' If `extracted_predicted_spawners == FALSE`, returns a list containing
#' `full_object`, the full fitted `stanfit` object (see [details](https://mc-stan.org/rstan/reference/stanfit-class.html)),
#' and `formatted_pars`, a formatted data table containing the following variables:
#' * **par_names** Parameter name
#' * **mean** Mean of the posterior distribution for a parameter
#' * **se_mean** Monte Carlo standard error for summary of all chains merged (see [details](https://mc-stan.org/rstan/reference/stanfit-method-summary.html))
#' * **sd** Standard deviation of the posterior distribution for a parameter
#' * **`2.5%`** 2.5% quantile of posterior distribution for a parameter.
#' * **`25%`** 25% quantile of posterior distribution for a parameter.
#' * **`50%`** 50% quantile of posterior distribution for a parameter.
#' * **`75%`** 75% quantile of posterior distribution for a parameter.
#' * **`97.5%`** 97.5% quantile of posterior distribution for a parameter.
#' * **n_eff** Effective sample size for a parameter
#' * **Rhat** Split Rhats for a parameter
#' * **stream** Stream name for the parameter estimate
#' * **year** Year observed data came from
#' @export
#' @family passage_to_spawner
#' @md
run_passage_to_spawner_model <- function(observed_adult_input, adult_model_covariates,
                                         stream_name, selected_covariate, extract_predicted_spawners = c(FALSE, TRUE)) {

  stream_name <- tolower(stream_name)

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

  if(!stream_name %in% c("battle creek", "clear creek", "mill creek", "deer creek", "butte creek")){
    cli::cli_abort("Incorrect stream name. Please pass an approved stream.")
    if(stream_name == "butte creek") {
      cli::cli_bullet("Butte Creek Passage to Spawner model is still in development and these results have not been tested or verified and should not be used for further analysis.")
    }
  }
  if(!selected_covariate %in% c("wy_type", "max_flow_std", "gdd_std", "passage_index", "median_passage_timing_std", "null_covar")){
    cli::cli_abort("Incorrect/Unavailable environmental covariate. Please pass an approved covariate name.")
  }

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

  formatted_pars = get_all_pars(stream_model_fit, stream_name)

  if(extract_predicted_spawners) {
    predicted_spawners = suppressWarnings(get_predicted_spawners_from_P2S(formatted_pars))
    return(predicted_spawners)
  } else {
    return(list("full_object" = stream_model_fit,
                "formatted_pars" = get_all_pars(stream_model_fit, stream_name)))
  }
}


#' @title Compare Environmental Covariates in Passage to Spawner (P2S) Model
#' @description This function iterates through adult input data for every stream and
#' environmental covariate, fiting the P2S model and pulling out diagnostic statistics. The
#' data is truncated to only include years for which all environmental covariates are available.
#' @returns A table containing the following variables:
#' * **par_names** Parameter name
#' * **stream** Mean of the posterior distribution for a parameter
#' * **mean** Monte Carlo standard error for summary of all chains merged (see [details](https://mc-stan.org/rstan/reference/stanfit-method-summary.html))
#' * **median** 50% quantile of posterior distribution for a parameter.
#' * **sd** standard deviation
#' * **covar_considered** the covariate considered
#' @export
#' @family passage_to_spawner
#' @md
compare_P2S_model_covariates <- function(observed_adult_input, adult_model_covariates) {

  cli::cli_alert("Fitting Passage to Spawner (P2S) model for all stream and covariate combinations")

  # combine observed counts and covariates and pivot wider
  data <- full_join(observed_adult_input,
                    adult_model_covariates,
                    by = c("year", "stream")) |>
    filter(!is.na(data_type)) |>
    pivot_wider(id_cols = c(year, stream, wy_type, max_flow_std, gdd_std,
                            median_passage_timing_std, passage_index),
                names_from = data_type,
                values_from = count)

  all_streams <- tibble(par_names = "DELETE",
                        stream = "DELETE",
                        mean = 0.0,
                        median = 0.0,
                        sd = 0.0,
                        lcl = 0.0,
                        ucl = 0.0,
                        covar_considered = "DELETE",
                        convergence_metric = NA)

  for(i in c("battle creek", "clear creek", "deer creek", "mill creek")) {
    stream_results <- compare_P2S_covariates_within_stream(data, i)
    all_streams <- bind_rows(all_streams, stream_results)
  }

  comparison_results <- all_streams |>
    filter(par_names != "DELETE")

  return(list("covariate_comparison_results" = comparison_results))

}

#' @title Compare Environmental Covariates in Passage to Spawner (P2S) Model Within a Stream
#' @description This function is called within `compare_P2S_model_covariates()` and iterates through all
#' environmental covariates for a stream, producing diagnostic statistics for that stream.
#' @returns A table containing the following variables:
#' * **par_names** Parameter name
#' * **stream** Mean of the posterior distribution for a parameter
#' * **mean** Monte Carlo standard error for summary of all chains merged. Seee [details](https://mc-stan.org/rstan/reference/stanfit-method-summary.html)
#' * **median** 50% quantile of posterior distribution for a parameter.
#' * **sd** standard deviation
#' * **covar_considered** the covariate considered
#' @export
#' @family passage_to_spawner
#' @md
compare_P2S_covariates_within_stream <- function(data, stream_name) {
  #TODO continue to debug

  # set percent female variable to 1 for carcass and holding surveys, 0.5 for redd
  percent_female <- ifelse(stream_name %in% c("battle creek", "clear creek", "mill creek"), 0.5, 1)

  # empty
  final_results <- tibble("par_names" = "DELETE",
                          "stream" = "DELETE",
                          "mean" = 0.0,
                          "median" = 0.0,
                          "sd" = 0.0,
                          "covar_considered" = "DELETE")

  covars_to_consider = c("wy_type", "max_flow_std", "gdd_std",
                         "null_covar", "passage_index") # no more "median_passage_timing_std" bc of sample size

  # loop through covariates to test
  for(i in 1:length(covars_to_consider)) {

    selected_covariate <- covars_to_consider[i]

    # drop NAs for all covariates being considered
    # makes sure you are using the same dataset for each cycle
    stream_data <- data |>
      mutate(observed_spawners = ifelse(stream == "deer creek", holding_count, redd_count)) |>
      filter(upstream_estimate > 0) |>
      filter(stream == stream_name) |>
      drop_na(wy_type, max_flow_std, gdd_std, observed_spawners, upstream_estimate) |>
      mutate(null_covar = 0)

    cli::cli_alert(paste0("Running model for ", stream_name, " and ", selected_covariate))

    covar <- stream_data |>
      pull(all_of(selected_covariate))

    stream_data_list <- list("N" = length(unique(stream_data$year)),
                             "input_years" = unique(stream_data$year),
                             "observed_passage" = stream_data$upstream_estimate,
                             "observed_spawners" = stream_data$observed_spawners,
                             "percent_female" = percent_female,
                             "environmental_covar" = covar,
                             "ss_total" = calculate_ss_tot(stream_data),
                             "average_upstream_passage" = mean(stream_data$upstream_estimate,
                                                               na.rm = TRUE))

    passage_to_spawner_STAN_code <- read_file("model-files/passage_to_spawner.txt")

    stream_model_fit <- rstan::stan(model_name = paste("passage_to_spawner", stream_name, selected_covariate, sep = "_"),
                                    model_code = passage_to_spawner_STAN_code,
                                    data = stream_data_list,
                                    chains = 3, iter = 20000*2, seed = 84735)

    # check rhat
    rhat_diagnostic <- bayesplot::rhat(stream_model_fit) |> unname()
    rhat_diagnostic <- rhat_diagnostic[rhat_diagnostic != "NaN"]
    rhats_over_threshold <- sum(rhat_diagnostic > 1.05) > 0

    if(rhats_over_threshold){
      convergence_rhat <- FALSE
    } else {
      convergence_rhat <- TRUE
    }

    stream_results <- get_all_pars(stream_model_fit, stream_name) |>
      filter(par_names %in% c("R2_data", "R2_fixed", "mean_redds_per_spawner",
                              "b1_survival", "sigma_redds_per_spawner") |
               str_detect(par_names, "spawner_abundance_forecast")) |> # TODO add pspawn after it's added to STAN code
      select(par_names, stream, mean, median = `50%`, sd, lcl = `2.5%`, ucl = `97.5%`) |>
      mutate(covar_considered = selected_covariate,
             convergence_metric = convergence_rhat)

    final_results <- bind_rows(final_results, stream_results)
  }

  final_results <- final_results |>
    filter(par_names != "DELETE")

  return(final_results)
}

