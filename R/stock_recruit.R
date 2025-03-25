#' Prepare inputs for stock-recruit model
#' @details This function prepares data for input into a stock recruit model.
#' @param con A valid connection to the model run database.
#' @param stream The stream for which you want to fit the model.
#' @param adult_data_type the type of survey the adult data came from. Either `passage`,
#' `redd`, `holding`, or `carcass`.
#' @param covariate The covariate you would like to use to fit the model. One of `rearing_max_flow`, `rearing_mean_flow`,
#' `rearing_median_flow`, `rearing_min_flow`, `spawning_max_flow`, `spawning_mean_flow`, `spawning_median_flow`, `spawning_min_flow`,
#' `spawning_above_13_temp_day`, `spawning_above_13_temp_week`, `spawning_gdd_spawn`, `spawning_weekly_max_temp_max`,
#' `spawning_weekly_max_temp_mean`, `spawning_weekly_max_temp_median`, or `null`.
#' @param truncate_dataset either `TRUE` or `FALSE`. If `TRUE`, will filter the dataset to only include years where
#' there are data for ALL covariates. If `FALSE`, will filter the dataset to include years where there are data
#' for YOUR SELECTED covariate.
#' @returns a list of tables:
#' * **stock_recruit_table** juvenile abundance and adult abundance by brood year.
#' * **full_covariate_table** full table of covariates associated with the site, stream, and brood years in the stock recruit table
#' * **truncated_covariate_table** table of covariates limited to only those years where all covariates are available (for covariate comparison analyses)
#' @export
#' @md
prepare_stock_recruit_inputs <- function(con, stream, adult_data_type,
                                         covariate, truncate_dataset = FALSE) {

  if(!DBI::dbIsValid(con)) {
    cli::cli_abort("Connection argument does not have a valid connection to the database. Please try reconnecting to the database using DBI::dbConnect")
  }

  if(!adult_data_type %in% c("redd", "holding", "passage", "carcass")) {
    cli::cli_abort("Please supply an approved adult data type value: either redd, carcass, holding, or passage.")
  }

  if(!covariate %in% c("rearing_max_flow", "rearing_mean_flow",
                       "rearing_median_flow", "rearing_min_flow", "spawning_max_flow",
                       "spawning_mean_flow", "spawning_median_flow", "spawning_min_flow",
                       "spawning_above_13_temp_day", "spawning_above_13_temp_week",
                       "spawning_gdd_spawn", "spawning_weekly_max_temp_max", "spawning_weekly_max_temp_mean",
                       "spawning_weekly_max_temp_median", "null")) {
    cli::cli_abort("Please supply an approved covariate value. See ?prepare_stock_recruit_inputs() for approved values.")
  }

  years_filter <- SRJPEdata::stock_recruit_year_lookup |>
    mutate(have_both_data = ifelse(adult & rst, TRUE, FALSE)) |>
    filter(have_both_data) |>
    distinct(brood_year, run_year, stream, site) |>
    filter(stream == !!stream)

  adult_data_and_covariates <- SRJPEdata::stock_recruit_model_inputs |>
    ungroup() |>
    rename(brood_year = year) |>
    right_join(years_filter, by = c("brood_year", "stream"))

  # get juvenile results from database
  # currently pulling most recent model run for the stream
  juv_results_from_db <- tbl(con, "model_parameters") |>
    left_join(tbl(con, "trap_location") |>
                select(location_id = id, site, stream),
              by = "location_id") |>
    left_join(tbl(con, "parameter") |>
                select(parameter_id = id, parameter = definition),
              by = "parameter_id") |>
    left_join(tbl(con, "statistic") |>
                select(statistic_id = id, statistic = definition),
              by = "statistic_id") |>
    filter(stream == !!stream,
           parameter == "Ntot",
           statistic %in% c("mean", "sd")) |>
    collect() |>
    group_by(updated_at, year) |>
    mutate(recent = cur_group_id()) |>
    ungroup() |>
    group_by(year) |>
    mutate(most_recent_by_year = max(recent)) |>
    ungroup() |>
    # filter to only keep the most recent run for each year
    filter(recent == most_recent_by_year) |>
    select(-c(id, location_id, parameter_id, statistic_id,
              updated_at, model_run_id, location_fit_id,
              recent, most_recent_by_year, week_fit)) |>
    distinct_all() # TODO this is an error in the model parameter upload. we are joining on location_id but there are multiple matches for knights landing, so we are getting duplicates of all the parameter esitimates.

  if(nrow(juv_results_from_db) == 0) {
    cli::cli_abort(paste("There are no model results in the database for", stream))
  }

  juv_results <- juv_results_from_db |>
    mutate(brood_year = year - 1) |>
    pivot_wider(names_from = "statistic",
                values_from = "value") |>
    mutate(mean_juvenile_abundance = mean,
           cv_juvenile_abundance = sd/mean) |>
    select(brood_year, site, mean_juvenile_abundance, cv_juvenile_abundance) |>
    right_join(years_filter, by = c("brood_year", "site"))

  data_with_adult_juv_and_covars <- full_join(adult_data_and_covariates, juv_results,
                                              by = c("brood_year", "stream", "run_year", "site")) |>
    mutate(null_covar = 0,
           adult_abundance = case_when(adult_data_type == "redd" ~ redd,
                                       adult_data_type == "holding" ~ holding,
                                       adult_data_type == "carcass" ~ carcass,
                                       adult_data_type == "passage" ~ passage)) |>
    filter(!is.na(adult_abundance))

  if(missing(covariate)) {
    cli::cli_bullets("Running a no-covariate model")
    covariate <- "null"
  }

  data_with_selected_covar_adult_data <- data_with_adult_juv_and_covars |>
    mutate(null = 0) |>
    select(stream, brood_year, site, run_year, adult_abundance, mean_juvenile_abundance,
           cv_juvenile_abundance, all_of(covariate)) |>
    filter(stream == !!stream) |>
    drop_na(adult_abundance, mean_juvenile_abundance, cv_juvenile_abundance,
            all_of(covariate))

  # get year lookup for matching with parameter estimates
  year_lookup <- data_with_selected_covar_adult_data |>
    mutate(year_index = row_number()) |>
    select(brood_year, year_index)

  if(nrow(data_with_selected_covar_adult_data) == 0) {
    cli::cli_abort(paste0("There is not enough data to run the model for ", stream,
                          " and covariate type ", covariate))
  }

  n_years <- nrow(data_with_selected_covar_adult_data)
  spawners <- data_with_selected_covar_adult_data$adult_abundance
  mean_log_recruits <- log(data_with_selected_covar_adult_data$mean_juvenile_abundance)
  sd_log_recruits <- sqrt(log(data_with_selected_covar_adult_data$cv_juvenile_abundance^2 + 1))
  mean_obsv_log_recruits_per_spawner <- vector(length = n_years)
  sd_obsv_log_recruits_per_spawner <- vector(length = n_years)

  # covariate
  covariate_data <- data_with_selected_covar_adult_data |>
    pull(all_of(covariate))

  # simulate variation in log recruits and divide by spawner stock to get mean and sd of observed log recruits-per-spawner
  for(i in 1:n_years) {
    simulated_recruits <- exp(rnorm(n = 1000, mean = mean_log_recruits[i], sd = sd_log_recruits[i]))
    observed_log_recruits_per_spawner <- log(simulated_recruits / spawners[i])
    mean_obsv_log_recruits_per_spawner[i] <- mean(observed_log_recruits_per_spawner)
    sd_obsv_log_recruits_per_spawner[i] <- sd(observed_log_recruits_per_spawner)
  }

  # data list for passing to STAN
  data_list <- list(Nyrs = n_years,
                    SP = spawners,
                    mu_obslgRS = mean_obsv_log_recruits_per_spawner,
                    sd_obslgRS = sd_obsv_log_recruits_per_spawner,
                    X = covariate_data)

  # inits
  if(covariate == "null") {
    reg <- lm(mean_obsv_log_recruits_per_spawner ~ spawners)
    ini_gamma <- 0
  } else {
    reg <- lm(mean_obsv_log_recruits_per_spawner ~ spawners + covariate_data)
    ini_gamma <- as.double(reg$coefficients[3])
  }

  ini_alpha <- as.double(reg$coefficients[1])
  ini_beta <- as.double(reg$coefficients[2])

  inits1 <- list(alpha = ini_alpha,
                 beta = ini_beta,
                 gamma = ini_gamma,
                 sd_pro = 1)

  inits2 <- list(alpha = ini_alpha * 1.5,
                 beta = ini_beta * 0.75,
                 gamma = ini_gamma * 0.8,
                 sd_pro = 1.25)

  inits3 <- list(alpha = ini_alpha * 0.75,
                 beta = ini_beta * 1.25,
                 gamma = ini_gamma * 1.2,
                 sd_pro = 0.75)

  inits <- list(inits1, inits2, inits3)

  return(list("data" = data_list,
              "inits" = inits,
              "year_lookup" = year_lookup))

}

#' Fit stock-recruit model in STAN
#' @details This function runs the STAN stock-recruit model.
#' @param stock_recruit_inputs the object produced by `prepare_stock_recruit_inputs()`
#' @returns a STAN object.
#' @export
#' @md
fit_stock_recruit_model <- function(stock_recruit_inputs) {

  cli::cli_process_start("Fitting stock-recruit model")
  stock_recruit_fit <- rstan::stan(model_code = SRJPEmodel::stock_recruit_model_code,# eval(parse(text = "SRJPEmodel::stock_recruit_model_code"))
                                   data = stock_recruit_inputs$data,
                                   init = stock_recruit_inputs$inits,
                                   chains = 3,
                                   iter = 1000,
                                   seed = 84735,  # TODO confirm setting seed
                                   include = T)

  cli::cli_process_done("stock-recruit fitting complete")

  return(stock_recruit_fit)
}

#' Extract stock-recruit parameter estimates
#' @details This function extracts parameter estimates from the stock recruit model object.
#' @description This function takes in a STAN fit object produced by running `fit_stock_recruit_model()` and produces a formatted
#' table of parameter estimates.
#' @param stock_recruit_inputs a named list produced by running `prepare_stock_recruit_inputs()`
#' @param stock_recruit_model_object a STANfit object produced by running `fit_stock_recruit_model()`
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
#' @md
extract_stock_recruit_estimates <- function(stock_recruit_inputs,
                                            stock_recruit_model_object) {

  summary_table <- rstan::summary(stock_recruit_model_object)$summary |>
    data.frame() |>
    tibble::rownames_to_column("parameter") |>
    mutate(year_index = suppressWarnings(readr::parse_number(parameter)),
           parameter = gsub("[0-9]+|\\[|\\]", "", parameter)) |>
    left_join(stock_recruit_inputs$year_lookup, by = "year_index") |>
    rename(`2.5%` = X2.5.,
           `25%` = X25.,
           `50%` = X50.,
           `75%` = X75.,
           `97.5%` = X97.5.) |>
    select(-year_index)

  if(any(summary_table$Rhat > 1.05)) {
    cli::cli_alert_warning("One or more parameters has an Rhat value over 1.05 :(")
  } else {
    cli::cli_bullets("No parameters have an Rhat value exceeding 1.05 :)")
  }

  summary_table_long <- summary_table |>
    pivot_longer(mean:Rhat, names_to = "statistic",
                 values_to = "value")

  return(summary_table_long)
}



