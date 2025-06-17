# adult model - predicted redds

#' @title Passage to Spawner Sum of Squares
#' @description This function calculates the total sum of squares, which is required
#' as input to the `data` call for `fit_passage_to_spawner_model()`.
#' @keywords internal
#' @export
#' @md
calculate_ss_tot <- function(data) {
  ss_tot <- 0

  N <- length(data$year)
  for(i in 1:N) {
    ss_tot <- ss_tot + (data$observed_spawners[i] - mean(data$observed_spawners))^2
  }
  return(ss_tot)
}

#' @title Prepare data for P2S
#' @description This function prepares data for running P2S on a single stream and covariate. It is called in `fit_passage_to_spawner_model()`.
#' @param stream_name: The name of the stream for which you would like to run the P2S model. Can be
#' `battle creek`, `clear creek`, `deer creek`, `mill creek`, or (draft form) `butte creek`.
#' @param selected_covariate: The environmental covariate you'd like to run the model for. Can be
#' either `wy_type` (water year type), `max_flow_std` (maximum flow), `gdd_std` (growing degree days),
#' or `null_covar` (no environmental covariate).
#' @param truncate_dataset a TRUE/FALSE value. If TRUE, will truncate the dataset so that only years where all covariate values are
#' present are used.
#' @returns a list containing data elements required to run the STAN passage to spawner model
#' @export
#' @md
prepare_P2S_inputs <- function(stream, selected_covariate, truncate_dataset = FALSE) {
  stream <- tolower(stream)

  # combine observed counts and covariates and pivot wider
  data <- SRJPEdata::annual_adult |>
    group_by(year, stream, data_type) |>
    summarise(count = sum(count, na.rm = T)) |> # count adipose clipped, run together
    ungroup() |>
    full_join(SRJPEdata::p2s_model_covariates_standard,
              by = c("year", "stream")) |>
    filter(!is.na(data_type)) |>
    pivot_wider(id_cols = c(year, stream, wy_type, max_flow_std, gdd_std),
                names_from = data_type,
                values_from = count) |>
    mutate(observed_spawners = case_when(stream == "deer creek" ~ holding,
                                         stream == "butte creek" ~ carcass_estimate,
                                         stream == "feather river" ~ broodstock_tag,
                                         TRUE ~ redd)) |>
    filter(upstream_estimate > 0,
           observed_spawners > 0,
           stream == !!stream) |>
    # null covariate
    mutate(null_covar = 0)

  # percent female variable is used in model_files/passage_to_spawner.stan to modify the parameter redds_per_spawner
  # as of January 2025 this model is not run on any streams where we are using non-redd adult data
  # for redd data, we assume that each redd represents two spawners (assuming a 50/50 sex ratio).
  # for carcass or holding data, a value of 1 indicates 1 spawner so we set `percent female` to 1 here.
  percent_female <- case_when(stream %in% c("battle creek", "clear creek", "mill creek") ~ 0.5,
                              stream %in% c("yuba river", "feather river", "butte creek", "deer creek") ~ 1)

  if(!stream %in% c("battle creek", "clear creek", "mill creek", "deer creek", "butte creek")){
    cli::cli_abort("Incorrect stream name. Please pass an approved stream.")
    if(stream == "butte creek") {
      cli::cli_bullet("Butte Creek Passage to Spawner model is still in development and these results have not been tested or verified and should not be used for further analysis.")
    }
  }
  if(!selected_covariate %in% c("wy_type", "max_flow_std", "gdd_std", "null_covar")){
    cli::cli_abort("Incorrect/Unavailable environmental covariate. Please pass an approved covariate name.")
  }

  # use holding count for deer, carcass for butte, and redd for other streams
  stream_data <- data

  if(truncate_dataset) {
    truncated_data <- stream_data |>
      # truncate dataset for all
      drop_na(wy_type, max_flow_std, gdd_std, observed_spawners, upstream_estimate)
  } else {
    # just drop the NAs for the selected covariate
    truncated_data <- stream_data  |>
      drop_na(all_of(selected_covariate), observed_spawners, upstream_estimate)
  }

  covar <- truncated_data |>
    dplyr::pull(all_of(selected_covariate))

  data <- list("N" = length(unique(truncated_data$year)),
               "input_years" = unique(truncated_data$year),
               "observed_passage" = truncated_data$upstream_estimate,
               "observed_spawners" = truncated_data$observed_spawners,
               "percent_female" = percent_female,
               "environmental_covar" = covar,
               "ss_total" = calculate_ss_tot(truncated_data),
               "average_upstream_passage" = mean(truncated_data$upstream_estimate, na.rm = TRUE))

  stream_data_list <- list("inputs" = list("data" = data),
                           "stream" = stream,
                           "model_name" = "p2s")

  return(stream_data_list)
}


#' @title Fit Passage to Spawner (P2S) STAN Model
#' @description This functions take in a list of data required to run the Passage to Spawner STAN model. The function
#' will return the full fitted `stanfit` object.
#' See \code{vignette("passage_to_spawner_submodel.Rmd", package = "SRJPEmodel")} for more details.
#' @param data_inputs A list object containing the following variables:
#' * **N** number of years
#' * **observed_passage** upstream passage for the associated years
#' * **observed_spawners** observed spawner counts for the associated years. Redd counts for Battle and Clear Creeks.
#' * **percent_female** Assumed sex ratio for the population and used to convert spawner counts to total population.
#' `0.5` for redds, `1` for other spawner survey types.
#' * **environmental_covar** the environmental covariate you want to run the model with.
#' * **ss_total** total sum of squares.
#' * **average_upstream_passage** average upstream passage across all years, used to forecast.
#' @returns a STANfit model object (see [details](https://mc-stan.org/rstan/reference/stanfit-class.html))
#' @export
#' @family passage_to_spawner
#' @md
fit_passage_to_spawner_model <- function(data_inputs) {

  p2s_model <- eval(parse(text = "SRJPEmodel::p2s_model_code"))

  cli::cli_process_start("Fitting P2S STAN model")
  p2s_model <- rstan::stan(model_code = p2s_model,
                           data = data_inputs$inputs$data,
                           chains = 3, iter = 2000, seed = 84735)

  cli::cli_process_done("P2S STAN model fitting complete")

  return(p2s_model)
}

#' @title Extract parameter estimates from P2S
#' @description This function takes in a STAN fit object produced by running `fit_passage_to_spawner_model()` and produces a formatted
#' table of parameter estimates from the Passage to Spawner model.
#' @param passage_to_spawner_model_object
#' @returns a table with the following variables:
#' * **parameter** Parameter name
#' * **mean** Mean of the posterior distribution for a parameter
#' * **se_mean** Monte Carlo standard error for summary of all chains merged (see [details](https://mc-stan.org/rstan/reference/stanfit-method-summary.html))
#' * **sd** Standard deviation of the posterior distribution for a parameter
#' * **`2.5`** 2.5% quantile of posterior distribution for a parameter.
#' * **`25`** 25% quantile of posterior distribution for a parameter.
#' * **`50`** 50% quantile of posterior distribution for a parameter.
#' * **`7%`** 75% quantile of posterior distribution for a parameter.
#' * **`97.5`** 97.5% quantile of posterior distribution for a parameter.
#' * **n_eff** Effective sample size for a parameter
#' * **Rhat** Split Rhats for a parameter
#' * **year** Year observed data came from
#' @export
#' @family passage_to_spawner
#' @md
extract_P2S_estimates <- function(p2s_inputs,
                                  passage_to_spawner_model_object){

  # get parameter estimates
  summary_table <- rstan::summary(passage_to_spawner_model_object)$summary |>
    data.frame() |>
    tibble::rownames_to_column("parameter") |>
    rename(`2.5` = X2.5.,
           `25` = X25.,
           `50` = X50.,
           `75` = X75.,
           `97.5` = X97.5.)

  # get years
  year_lookup <- summary_table |>
    filter(str_detect(parameter, "years")) |>
    mutate(index = readr::parse_number(parameter),
           year = round(mean)) |>
    select(index, year)

  # bind in years
  table_with_years <- summary_table |>
    filter(!str_detect(parameter, "years")) |>
    mutate(index = suppressWarnings(readr::parse_number(parameter)),
           parameter = str_remove_all(parameter, "\\[|\\]|[0-9]+")) |>
    left_join(year_lookup, by = "index") |>
    select(-index)

  if(any(table_with_years$Rhat > 1.05)) {
    cli::cli_alert_warning("Some estimates have an Rhat statistic over 1.05.")
  }

  table_with_years_long <- table_with_years |>
    pivot_longer(mean:Rhat, names_to = "statistic",
                 values_to = "value") |>
    mutate(model_name = "p2s",
           site = NA,
           week_fit = NA,
           stream = p2s_inputs$stream,
           srjpedata_version = as.character(packageVersion("SRJPEdata"))) |>
    select(model_name, site, stream, year, week_fit,
           parameter, statistic, value, srjpedata_version)

  return(table_with_years_long)

}


#' @title Compare Environmental Covariates in Passage to Spawner (P2S) Model
#' @description This function fits the Passage to Spawner (`fit_passage_to_spawner_model()`) for a given stream to all
#' environmental covariates. For a stream, the dataset is truncated to only include years for which all environmental
#' covariates are available. Parameters that are helpful for selecting a covariate are: `R2_data`, `R2_fixed`, `mean_redds_per_spawner`,
#' `b1_survival`, `sigma_redds_per_spawner`, and `spawner_abundance_forecast`.
#' @returns A table containing the following variables:
##' * **parameter** Parameter name
#' * **mean** Mean of the posterior distribution for a parameter
#' * **se_mean** Monte Carlo standard error for summary of all chains merged (see [details](https://mc-stan.org/rstan/reference/stanfit-method-summary.html))
#' * **sd** Standard deviation of the posterior distribution for a parameter
#' * **`2.5`** 2.5% quantile of posterior distribution for a parameter.
#' * **`25`** 25% quantile of posterior distribution for a parameter.
#' * **`50`** 50% quantile of posterior distribution for a parameter.
#' * **`7%`** 75% quantile of posterior distribution for a parameter.
#' * **`97.5`** 97.5% quantile of posterior distribution for a parameter.
#' * **n_eff** Effective sample size for a parameter
#' * **Rhat** Split Rhats for a parameter
#' * **year** Year observed data came from
#' @export
#' @family passage_to_spawner
#' @md
compare_P2S_model_covariates <- function(stream) {

  cli::cli_alert(paste0("Fitting Passage to Spawner (P2S) model to all covariates for", stream))

  approved_covariates <- c("wy_type", "max_flow_std", "gdd_std")

  inputs <- purrr::pmap(list(stream,
                             approved_covariates,
                             TRUE),
                        prepare_P2S_inputs)
  names(inputs) <- approved_covariates

  fits <- lapply(inputs, fit_passage_to_spawner_model)

  estimates <- lapply(fits, extract_P2S_estimates)

  estimate_table_with_covariate_names <- map2(estimates, names(estimates), function(df, name) {
                                                                            new_df <- df |>
                                                                              mutate(covariate = name)
                                                                            return(new_df)
                                                                          }) |>
    bind_rows()

  return(estimate_table_with_covariate_names)

}

#' P2S results plots
#' @details This function produces a plot with predicted spawners with uncertainty.
#' for simple comparison.
#' @param p2s_inputs inputs for the fit model, created by running `prepare_P2S_inputs()`
#' @param con a connection to the database.
#' @returns A plot showing observed and predicted spawners by year.
#' @export
#' @md
generate_results_plot_p2s <- function(p2s_inputs, con) {

  dark_JPE <- c("#F5CAC2", "#6E9881", "#9A8723", "#2D4755", "#869AA0")

  params <- get_most_recent_model_results(con) |>
    filter(model_name == "p2s",
           stream == p2s_inputs$stream) |>
    select(-c(id, week_fit, model_run_id, site, model_name))

  obsv_spawners <- tibble("year" = p2s_inputs$inputs$data$input_years,
                          "obsv_spawners" = p2s_inputs$inputs$data$observed_spawners)

  plot <- params |>
    pivot_wider(names_from = "statistic", values_from = "value") |>
    filter(parameter == "predicted_spawners") |>
    left_join(obsv_spawners, by = "year") |>
    mutate(year = factor(year)) |>
    ggplot() +
    geom_ribbon(aes(x = year, ymin = `2.5`, ymax = `97.5`, group = 1), alpha = 0.4, fill = dark_JPE[5]) +
    geom_line(aes(x = year, y = `50`, group = 1)) +
    geom_point(aes(x = year, y = obsv_spawners), size = 2, color = dark_JPE[2]) +
    theme_minimal() +
    labs(x = "Year",
         y = "Predicted spawners",
         title = paste0(str_to_title(p2s_inputs$stream)))

  return(plot)
}
