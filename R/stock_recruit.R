#' Prepare inputs for stock-recruit model
#' @details This function prepares data for input into a stock recruit model.
#' @param con A valid connection to the model run database.
#' @param adult_data_type the type of survey the adult data came from. Either `upstream_estimate`,
#' `redd_count`, `holding_count`, or `carcass_estimate`.
#' @site The stream for which you are fitting the model.
#' @model_run_id The model run from which you want to pull your juvenile abundance estimates.
#' @covariate The covariate you would like to use to fit the model. Either `gdd_spawn`,
#' `above_13_temp_day`, `above_13_temp_week`, `weekly_max_temp_max`, `weekly_max_temp_mean`, `weekly_max_temp_median`,
#' `mean_flow`, `max_flow`, `min_flow`, `median_flow`, or `null`
#' @truncate_dataset either `TRUE` or `FALSE`. If `TRUE`, will filter the dataset to only include years where
#' there are data for ALL covariates. If `FALSE`, will filter the dataset to include years where there are data
#' for YOUR SELECTED covariate.
#' @returns a list of tables:
#' * **stock_recruit_table** juvenile abundance and adult abundance by brood year.
#' * **full_covariate_table** full table of covariates associated with the site, stream, and brood years in the stock recruit table
#' * **truncated_covariate_table** table of covariates limited to only those years where all covariates are available (for covariate comparison analyses)
#' @export
#' @md
prepare_stock_recruit_inputs <- function(con, year, stream, adult_data_type,
                                         covariate, truncate_dataset = FALSE) {

  years_filter <- SRJPEdata::stock_recruit_year_lookup |>
    mutate(have_both_data = ifelse(adult & rst, TRUE, FALSE)) |>
    filter(have_both_data) |>
    distinct(brood_year, run_year, stream, site) |>
    filter(stream == !!stream)

  adult_data_and_covariates <- SRJPEdata::stock_recruit_model_inputs |>
    ungroup() |>
    rename(brood_year = year) |>
    right_join(years_filter, by = c("brood_year", "stream")) |>
    glimpse()

  # get juvenile results from database
  # currently pulling most recent model run for the stream
  # TODO we are going to have to also make sure we're pulling all years? so do we want to filter here?
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

  if(nrow(juv_results_from_db == 0)) {
    cli::cli_abort(paste("There are no model results in the database for", stream))
  }

  juv_results <- juv_results_from_db |>
    mutate(brood_year = year - 1) |>
    pivot_wider(names_from = "statistic",
                values_from = "value") |>
    mutate(mean_juvenile_abundance = mean,
           cv_juvenile_abundance = sd/mean) |>
    select(brood_year, site, mean_juvenile_abundance, cv_juvenile_abundance) |>
    right_join(years_filter, by = c("brood_year", "site")) |>
    glimpse()

  data_with_adult_juv_and_covars <- full_join(adult_data_and_covariates, juv) |>
    mutate(null_covar = 0) |>
    # TODO filter by adult data type
    glimpse()

  if(missing(covariate)) {
    cli::cli_bullets("Running a no-covariate model")
    covariate <- "null"
  }

  data_with_selected_covar <- data |>
    mutate(null = 0) |>
    select(stream, brood_year, site, juv_year, adult_abundance, mean_juvenile_abundance,
           cv_juvenile_abundance, all_of(covariate)) |>
    filter(stream == !!stream) |>
    drop_na(adult_abundance, mean_juvenile_abundance, cv_juvenile_abundance,
            all_of(covariate))

  # get year lookup for matching with parameter estimates
  year_lookup <- data_with_selected_covar |>
    mutate(year_index = row_number()) |>
    select(brood_year, year_index)

  if(nrow(data_with_selected_covar) == 0) {
    cli::cli_abort(paste0("There is not enough data to run the model for ", stream,
                          " and covariate type ", covariate))
  }

  n_years <- nrow(data_with_selected_covar)
  spawners <- data_with_selected_covar$adult_abundance
  mean_log_recruits <- log(data_with_selected_covar$mean_juvenile_abundance)
  sd_log_recruits <- sqrt(log(data_with_selected_covar$cv_juvenile_abundance^2 + 1))
  mean_obsv_log_recruits_per_spawner <- vector(length = n_years)
  sd_obsv_log_recruits_per_spawner <- vector(length = n_years)

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
                    X = data_with_selected_covar$covariate)

  # inits
  if(covariate == "null") {
    reg <- lm(mean_obsv_log_recruits_per_spawner ~ spawners)
    ini_gamma <- 0
  } else {
    reg <- lm(mean_obsv_log_recruits_per_spawner ~ spawners + standardized_covariate)
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



