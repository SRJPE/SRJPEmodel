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
    # this filters to the appropriate stream via the right join
    right_join(years_filter, by = c("brood_year", "stream"))

  # get juvenile results from database
  juv_results_from_db <- SRJPEmodel::get_most_recent_model_output(con) |>
    filter(model_name == "bt_spas_x",
           stream == !!stream,
           parameter == "Ntot",
           statistic %in% c("mean", "sd"))

  if(nrow(juv_results_from_db) == 0) {
    cli::cli_abort(paste("There are no model results in the database for", stream))
  }

  juv_results <- juv_results_from_db |>
    mutate(brood_year = year - 1) |>
    pivot_wider(names_from = "statistic",
                values_from = "value",
                id_cols = c(brood_year, stream, site)) |>
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
    filter(stream == !!stream,
           # filter out selected adult data type values of 0
           adult_abundance > 0) |>
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
  observed_recruits <- data_with_selected_covar_adult_data$mean_juvenile_abundance
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
                    R = observed_recruits,
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

  return(list(inputs = list("data" = data_list,
                            "inits" = inits),
              "year_lookup" = year_lookup,
              "stream" = stream,
              "data_type" = adult_data_type,
              "truncate_dataset" = truncate_dataset,
              "covariate_name" = covariate,
              "model_name" = "stock_recruit"))

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
                                   data = stock_recruit_inputs$inputs$data,
                                   init = stock_recruit_inputs$inputs$inits,
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
    rename(`2.5` = X2.5.,
           `25` = X25.,
           `50` = X50.,
           `75` = X75.,
           `97.5` = X97.5.) |>
    select(-year_index)

  if(any(summary_table$Rhat > 1.05)) {
    cli::cli_alert_warning("One or more parameters has an Rhat value over 1.05 :(")
  } else {
    cli::cli_bullets("No parameters have an Rhat value exceeding 1.05 :)")
  }

  summary_table_long <- summary_table |>
    pivot_longer(mean:Rhat, names_to = "statistic",
                 values_to = "value") |>
    rename(year = brood_year) |>
    mutate(model_name = "stock_recruit",
           site = NA,
           week_fit = NA,
           location_fit = stock_recruit_inputs$stream,
           srjpedata_version = as.character(packageVersion("SRJPEdata")))

  return(summary_table_long)
}

#' Stock-recruit diagnostic plots
#' @details This function produces a plot with observed spawner and estimated recruitment data,
#' alongside predicted recruitment across a range of spawning stock sizes produced using
#' posterior estimates from the model fit. It is meant to compare predicted and observed
#' recruitment as a function of spawning stock.
#' @param sr_inputs inputs for the fit model, created by running `prepare_stock_recruit_inputs()`
#' @param sr_fit the model fit object, created by running `fit_stock_recruit_model()`
#' @returns A plot.
#' @export
#' @md
generate_diagnostic_plot_sr <- function(sr_inputs, sr_fit) {
  # Load posteriors and calculate predicted recruitment across a range of
  # spawning stock sizes (rather than just observed ones as model does)
  spawner_vector <- seq(0, max(sr_inputs$data$SP) * 1.05, length.out = 50)
  posteriors <- as.data.frame(sr_fit, pars = c("alpha", "beta", "gamma", "sd_pro"))
  pred_recruitment_matrix <- matrix(nrow = 50, ncol = 3)

  for(i in 1:50) {
    pred_recruitment <- spawner_vector[i] * exp(posteriors$alpha + posteriors$beta * spawner_vector[i]) * 0.001 # scale
    pred_recruitment_matrix[i, ] <- quantile(pred_recruitment, probs = c(0.025, 0.5, 0.975))
    pred_recruitment_matrix[i, 2] <- mean(pred_recruitment) # replace median with mean
  }

  # create observed data frame that also includes effect of covariate size
  # estimated using the posteriors
  obsv_data <- tibble("obsv_spawners" = sr_inputs$data$SP,
                      "obsv_recruits" = sr_inputs$data$R * 0.001,
                      "obsv_covar" = sr_inputs$data$X,
                      "brood_year" = paste0("\'", substr(sr_inputs$year_lookup$brood_year, 3, 4)))

  # calculate covariate effects
  # TODO can we add this code into the mutate?
  pred_covar <- c()
  pred_covar_w_gamma <- c()
  for(i in 1:sr_inputs$data$Nyrs) {
    pred_covar[i] <- mean(obsv_data$obsv_spawners[i] * exp(posteriors$alpha + posteriors$beta * obsv_data$obsv_spawners[i])) * 0.001
    pred_covar_w_gamma[i] <- mean(obsv_data$obsv_spawners[i] * exp(posteriors$alpha + posteriors$beta * obsv_data$obsv_spawners[i] +
                                                                     posteriors$gamma * obsv_data$obsv_covar[i])) * 0.001
  }

  obsv_data <- obsv_data |>
    mutate(pred_covar = pred_covar,
           pred_covar_w_gamma = pred_covar_w_gamma,
           covar_effect_color = ifelse(pred_covar_w_gamma > pred_covar & (obsv_recruits > pred_covar) | pred_covar_w_gamma < pred_covar & (obsv_recruits < pred_covar),
                                       "covariate effect (consistent)",
                                       "covariate effect (inconsistent)"))

  # create data frame with a range of spawning stock sizes and
  # associated predicted recruitment across that range
  pred_data <- tibble("pred_spawners" = spawner_vector,
                      "pred_recruit_025" = pred_recruitment_matrix[, 1],
                      "pred_recruit_mean" = pred_recruitment_matrix[, 2],
                      "pred_recruit_975" = pred_recruitment_matrix[, 3]) |>
    mutate(color = "@ average cov value")


  plot <- ggplot() +
    # range of stock sizes
    geom_ribbon(data = pred_data,
                aes(spawner_vector, ymin = pred_recruit_025, ymax = pred_recruit_975),
                alpha = 0.6, fill = "grey") +
    geom_line(data = pred_data,
              aes(x = spawner_vector, y = pred_recruit_mean,
                  color = color)) +
    # observed data
    geom_errorbar(data = obsv_data,
                  aes(x = obsv_spawners,
                      ymin = pred_covar,
                      ymax = pred_covar_w_gamma,
                      color = covar_effect_color),
                  width = 0) +
    geom_point(data = obsv_data,
               aes(x = obsv_spawners, y = obsv_recruits)) +
    geom_text(data = obsv_data,
              aes(x = obsv_spawners + 5, y = obsv_recruits + 3,
                  label = brood_year),
              size = 3) +
    labs(x = "Spawner abundance",
         y = "Outmigrant abundance ('000s)") +
    theme_minimal() +
    scale_color_manual(name = "",
                       breaks = c("@ average cov value", "covariate effect (consistent)", "covariate effect (inconsistent)"),
                       values = c("@ average cov value" = "black",
                                  "covariate effect (consistent)" = "blue",
                                  "covariate effect (inconsistent)" = "red")) +
    theme(legend.position = "bottom")

  return(plot)
}


#' Stock-recruit diagnostic plots
#' @details This function produces a plot with observed spawner and estimated recruitment data,
#' alongside predicted recruitment across a range of covariate values produced using
#' posterior estimates from the model fit. It is meant to compare predicted and observed
#' recruitment as a function of spawner effect on covariate.
#' @param sr_inputs inputs for the fit model, created by running `prepare_stock_recruit_inputs()`
#' @param sr_fit the model fit object, created by running `fit_stock_recruit_model()`
#' @returns A plot.
#' @export
#' @md
generate_diagnostic_plot_sr_covar <- function(sr_inputs, sr_fit) {
  # Load posteriors and calculate predicted recruitment across a range of
  # covariate sizes (rather than just observed ones as model does)
  # TODO confirm we will still have these data, being updated?
  raw_covar <- SRJPEdata::stock_recruit_covariates |>
    ungroup() |>
    mutate(covar_nm = ifelse(lifestage == "spawning and incubation",
                             paste0("spawning_", covariate_structure),
                             paste0(lifestage, "_", covariate_structure))) |>
    filter(covar_nm == sr_inputs$covariate_name,
           stream == sr_inputs$stream,
           year %in% sr_inputs$year_lookup$brood_year) |>
    pull(value)

  # create covariate vector across range of values
  covar_vector <- seq(min(raw_covar), max(raw_covar), length.out = 50)
  std_covar_vector <- (covar_vector - mean(covar_vector)) / sd(covar_vector)
  posteriors <- as.data.frame(sr_fit, pars = c("alpha", "beta", "gamma", "sd_pro"))
  pred_recruitment_matrix <- matrix(nrow = 50, ncol = 3)
  mean_spawners <- mean(sr_inputs$data$SP)

  for(i in 1:50) {
    pred_recruitment <- mean_spawners * exp(posteriors$alpha + posteriors$beta * mean_spawners +
                                              posteriors$gamma * std_covar_vector[i]) * 0.001 # scale
    pred_recruitment_matrix[i, ] <- quantile(pred_recruitment, probs = c(0.025, 0.5, 0.975))
    pred_recruitment_matrix[i, 2] <- mean(pred_recruitment) # replace median with mean
  }

  # create observed data frame that also includes effect of covariate size
  # estimated using the posteriors
  obsv_data <- tibble("obsv_spawners" = sr_inputs$data$SP,
                      "obsv_recruits" = sr_inputs$data$R * 0.001,
                      "obsv_raw_covar" = raw_covar,
                      "obsv_covar" = sr_inputs$data$X,
                      "brood_year" = paste0("\'", substr(sr_inputs$year_lookup$brood_year, 3, 4)))

  # calculate covariate effects
  pred_covar_1 <- c()
  pred_covar_2 <- c()
  for(i in 1:sr_inputs$data$Nyrs) {
    pred_covar_1[i] <- mean(mean_spawners * exp(posteriors$alpha + posteriors$beta *
                                                  mean_spawners + posteriors$gamma * obsv_data$obsv_covar[i])) * 0.001
    pred_covar_2[i] <- mean(obsv_data$obsv_spawners[i] * exp(posteriors$alpha + posteriors$beta * obsv_data$obsv_spawners[i] +
                                                               posteriors$gamma * obsv_data$obsv_covar[i])) * 0.001
  }

  obsv_data <- obsv_data |>
    mutate(pred_covar_1 = pred_covar_1,
           pred_covar_2 = pred_covar_2)



  # create data frame with a range of spawning stock sizes and
  # associated predicted recruitment across that range
  pred_data <- tibble("covar_vector" = covar_vector,
                      "pred_recruit_025" = pred_recruitment_matrix[, 1],
                      "pred_recruit_mean" = pred_recruitment_matrix[, 2],
                      "pred_recruit_975" = pred_recruitment_matrix[, 3]) |>
    mutate(color = "@ average cov value")

  plot <- ggplot() +
    # range of covar values
    geom_ribbon(data = pred_data,
                aes(x = covar_vector, ymin = pred_recruit_025, ymax = pred_recruit_975),
                alpha = 0.6, fill = "grey") +
    geom_line(data = pred_data,
              aes(x = covar_vector, y = pred_recruit_mean,
                  color = "@ average cov value")) +
    # observed data
    geom_errorbar(data = obsv_data,
                  aes(x = obsv_raw_covar,
                      ymin = pred_covar_1,
                      ymax = pred_covar_2,
                      color = "with spawner effect"),
                  width = 0) +
    geom_point(data = obsv_data,
               aes(x = obsv_raw_covar, y = obsv_recruits))  +
    geom_text(data = obsv_data,
              aes(x = obsv_raw_covar + 0.5, y = obsv_recruits + 0.5, # TODO this "jigger" will be different for different covariates
                  label = brood_year),
              size = 3) +
    labs(x = sr_inputs$covariate_name,
         y = "Outmigrant abundance ('000s)") +
    theme_minimal() +
    scale_color_manual(name = "",
                       breaks = c("@ average cov value", "with spawner effect"),
                       values = c("@ average cov value" = "black",
                                  "with spawner effect" = "red")) +
    theme(legend.position = "bottom")

  return(plot)
}

#' Stock-recruit results plots
#' @details This function produces a plot with observed and predicted recruits
#' for simple comparison.
#' @param sr_inputs inputs for the fit model, created by running `prepare_stock_recruit_inputs()`
#' @param con a connection to the database.
#' @returns A plot.
#' @export
#' @md
generate_results_plot_sr <- function(sr_inputs, con) {

  dark_JPE <- c("#F5CAC2", "#6E9881", "#9A8723", "#2D4755", "#869AA0")

  obsv_R <- sr_inputs$year_lookup |>
    mutate(obsv_recruits = sr_inputs$data$R)

  params <- get_most_recent_model_output(con) |>
    filter(model_name == "stock_recruit",
           stream == sr_inputs$stream) |>
    rename(brood_year = year)

  errors <- params |>
    filter(parameter == "pred_R",
           statistic %in% c("2.5", "97.5")) |>
    select(-parameter) |>
    pivot_wider(names_from = "statistic",
                values_from = "value") |>
    mutate(type = "predicted")

  plot_data <- params |>
    filter(parameter == "pred_R",
           statistic == "50") |>
    left_join(obsv_R, by = "brood_year") |>
    select(-c(statistic, year_index, parameter)) |>
    pivot_longer(c(value, obsv_recruits),
                 names_to = "type",
                 values_to = "recruits") |>
    mutate(type = ifelse(type == "value", "predicted", "observed")) |>
    left_join(errors, by = c("brood_year", "type")) |>
    mutate(brood_year = factor(brood_year),
           `2.5`  = `2.5` * 0.001,
           `97.5` = `97.5` * 0.001,
           recruits = recruits * 0.001)

  plot <- plot_data |>
    ggplot(aes(x = brood_year, y = recruits, fill = type)) +
    geom_col(position = "dodge") +
    geom_errorbar(aes(x = brood_year, ymin = `2.5`, ymax = `97.5`),
                  position = position_dodge(1),
                  width = 0.2) +
    scale_fill_manual(values = dark_JPE) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(x = "Brood Year",
         y = "Recruits ('000s)")

  return(plot)

}


