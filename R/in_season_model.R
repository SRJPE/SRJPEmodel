#' Prepare inputs for in season model
#' @details This function prepares data for input into a in season model.
#' @param con A valid connection to the model run database.
#' @param stream The stream for which you want to fit the model.
#' @param site The site for which you want to fit the model.
#' @returns a list of inputs TBD
#' @export
#' @md
prepare_inseason_inputs <- function(con, stream, site, covariate_effect, autocorrelation) {

  if(!DBI::dbIsValid(con)) {
    cli::cli_abort("Connection argument does not have a valid connection to the database. Please try reconnecting to the database using DBI::dbConnect")
  }

  # TODO our dates are slightly off (3 days?) from josh's weeks
  # set up weeks processing
  number_weeks <- nrow(SRJPEmodel::julian_week_to_date_lookup)
  first_week <- 36
  last_week <- 35 #All year Sep-01 - Aug-31
  weeks_ordered <- c(first_week:53, 1:last_week)
  theta <- (1:number_weeks) / number_weeks # set proportional week (proportion of total year through week iwk)
  theta[number_weeks] <- 0.999 # set last value to 0.999 instead of 1

  # get juvenile results from database
  juv_results_from_db <- SRJPEmodel::get_most_recent_model_results(con) |>
    filter(model_name == "bt_spas_x",
           stream == !!stream,
           site == !!site,
           parameter %in% c("Ntot", "N"),
           statistic %in% c("mean", "sd")) |>
    select(-c(id, model_name))

  model_run_ids <- juv_results_from_db |>
    distinct(model_run_id, year)

  inputs <- lapply(model_run_ids$model_run_id, get_model_object,
                   con = con,
                   model_component = "model_input",
                   keyword = NULL,
                   access_key = Sys.getenv("AZ_CONTAINER_ACCESS_KEY"))
  names(inputs) <- model_run_ids$year

  objects <- lapply(model_run_ids$model_run_id, get_model_object,
                    con = con,
                    model_component = "model_fit",
                    keyword = NULL,
                    access_key = Sys.getenv("AZ_CONTAINER_ACCESS_KEY"))
  names(objects) <- model_run_ids$year

  number_years <- length(unique(juv_results_from_db$year)) # number of years fit for that stream
  number_estimated_weeks <- vector(length = number_years)
  years <- vector(length = number_years)
  mean_total_abundance <- vector(length = number_years)
  sd_total_abundance <- vector(length = number_years)
  cv_total_abundance <- vector(length = number_years)

  estimated_weeks_matrix <- matrix(nrow = number_weeks, ncol = number_years, 0) # number of weeks fit within a year
  proportional_weeks_matrix <- matrix(nrow = number_weeks, ncol = number_years, 0) # proportional week for each week-year combination. Varies by year based on when they sampled at RST
  mean_abundance_through_week <- matrix(nrow = number_weeks, ncol = number_years, 0)
  sd_abundance_through_week <- matrix(nrow = number_weeks, ncol = number_years, 0)
  cv_abundance_through_week <- matrix(nrow = number_weeks, ncol = number_years, 0)
  covariate_matrix <- matrix(nrow = number_years, ncol = 2) # annual covariates affecting phi [, 1] and lambda [,2]
  covariate_matrix_standardized <- matrix(nrow = number_years, ncol = 2)

  annual_abundance <- juv_results_from_db |>
    filter(parameter == "Ntot") |>
    pivot_wider(names_from = "statistic",
                values_from = "value") |>
    mutate(cv = sd / mean,
           log_N = log(mean))

  weekly_abundance <- juv_results_from_db |>
    filter(parameter == "N") |>
    pivot_wider(names_from = "statistic",
                values_from = "value") |>
    mutate(cv = sd / mean) |>
    arrange(year, week_fit)

  for(year in model_run_ids$year) {
    # iterating through model run ID is like iterating through years (it is a 1-1 relationship)
    input_list <- inputs[[as.character(year)]]
    weeks_fit <- input_list$weeks_fit
    catch_flow <- input_list$catch_flow_raw
    #max_flow <- max(input_list$catch_flow_raw)
    nstrata <- input_list$inputs$data$Nstrata
    catch <- input_list$inputs$data$u
    year_index <- which(model_run_ids$year == year)
    model_object <- objects[[as.character(year)]]

    #strip out some of initial and last weeks with 0 catch
    trim_run <- FALSE # TODO do we need to keep this?

    if(trim_run == TRUE){
      week_buffer <- 3
      first_week_trim <- min(which(catch > 0)) - week_buffer; if(first_week_trim <= 0) first_week_trim <- 1
      last_week_trim <- max(which(catch > 0)) + week_buffer; if(last_week_trim > nstrata) last_week_trim <- nstrata
      number_estimated_weeks[year_index] <- last_week_trim = first_week_trim + 1
    } else { #Use all weeks where abundance was estimated
      first_week_trim <- 1
      number_estimated_weeks[year_index] <- nstrata
    }

    for(week in 1:number_estimated_weeks[year_index]) {
      estimated_weeks_matrix[week, year_index] <- which(weeks_ordered == weeks_fit[week + first_week_trim - 1])
      proportional_weeks_matrix[week, year_index] <- theta[estimated_weeks_matrix[week, year_index]]

      column_index <- first_week_trim:(week + first_week_trim - 1)
      N <- model_object$sims.list$N[ , column_index]

      if(week == 1) {
        cumulative_abundance <- N
      } else {
        cumulative_abundance <- rowSums(N)
      }

      mean_abundance <- mean(cumulative_abundance)
      sd_abundance <- sd(cumulative_abundance)
      if(mean_abundance == 0) mean_abundance <- 0.01
      if(sd_abundance == 0) sd_abundance <- 0.01

      mean_abundance_through_week[week, year_index] <- log(mean_abundance) # total abundance through this week in log space
      cv_abundance_through_week[week, year_index] <- sd_abundance / mean_abundance # sd across posterior samples of cumulative abundance
      sd_abundance_through_week[week, year_index] <- sqrt(log(cv_abundance_through_week[week, year_index] ^ 2 + 1)) # convert from cv in untransformed space to sd in log space
    }

    flow_index <- which(weeks_fit <= 53 | (weeks_fit >= 1 & weeks_fit <= 5)) # flows prior to February
    covariate_matrix[year_index, 1] <- max(catch_flow[flow_index])
    covariate_matrix[year_index, 2] <- max(catch_flow[flow_index])# use same covariate for prediction of phi and lambda but this structure allows one to use different covariates
  }

  # standardize covariate
  for(j in 1:2) {
    mean_cov <- mean(covariate_matrix[ , j])
    sd_cov <- sd(covariate_matrix[, j])
    covariate_matrix_standardized[, j] <- (covariate_matrix[, j] - mean_cov) / sd_cov
  }

  cv_abundance_through_week_transposed <- t(cv_abundance_through_week)
  colnames(cv_abundance_through_week_transposed) <- paste0("lg_Nx_cv_", 1:number_weeks)

  if(covariate_effect) {
    UseCovX <- c(1, 0)
  } else {
    UseCovX <- c(0, 0)
  }

  data <- list(Nyrs = number_years,
               Nwks = number_weeks,
               Nestwks = number_estimated_weeks,
               theta = theta,
               ewk = estimated_weeks_matrix,
               Ntot_mu = annual_abundance$mean,
               Ntot_sd = annual_abundance$sd,
               Nx_mu = mean_abundance_through_week,
               Nx_sd = sd_abundance_through_week,
               CovX = covariate_matrix_standardized,
               UseCovX = UseCovX,
               pr_sd_pro = c(20, 10))# shape and rate pars for gamma prior on sd_pro, median of 2, 95 %CI's 1.2-3 which covers range for sites where model converges without a gamma prior

  # inits
  iphi <- 0.4
  ilam <- 50
  mu_phi <- log(iphi / (1 - iphi))
  mu_lambda <- log(ilam)
  irt <- cbind(rep(mu_phi, number_years), rep(mu_lambda, number_years))

  init_list <- list(muRt = c(mu_phi, mu_lambda),
                    RTpars = irt,
                    sigmaRT = c(1, 1),
                    bCov = c(0, 0),
                    sd_pro = 1)

  inits <- list(init_list, init_list, init_list)

  # model name
  model_name <- ifelse(autocorrelation, "BetaDevHBMRT_lag1", "BetaDevHBMRT")

  return(list(inputs = list(data = data,
                            inits = inits),
              covariate_effect = covariate_effect,
              autocorrelation = autocorrelation,
              model_name = model_name,
              stream = stream,
              site = site,
              year_lookup = annual_abundance$year)) # TODO update model name table in database
}

#' @title Fit In-Season STAN model
#' @description This functions take in a list of data required to run the in-season STAN model. The function
#' will return the full fitted `stanfit` object.
#' See \code{vignette("passage_to_spawner_submodel.Rmd", package = "SRJPEmodel")} for more details.
#' @param inputs a list generated by running `prepare_in_season_inputs()`
#' @returns a STANfit model object (see [details](https://mc-stan.org/rstan/reference/stanfit-class.html))
#' @export
#' @md
fit_inseason_model <- function(inputs) {

  parameters <- c("phi", "lambda", "cp", "sd_pro", "RTpars", "muRT", "sigmaRT",
                  "rho", "mu_phi", "mu_lambda", "bCov", "For_cp")

  if(inputs$autocorrelation) {
    parameters <- c(parameters, "rho_pro")
    stan_model <- eval(parse(text = "SRJPEmodel::in_season_model_code$autocorrelation"))
  } else {
    stan_model <- eval(parse(text = "SRJPEmodel::in_season_model_code$no_autocorrelation"))
  }

  fit <- rstan::stan(model_code = stan_model,
                     data = inputs$inputs$data,
                     init = inputs$inputs$inits,
                     chains = 3,
                     iter = 2000,
                     include = T,
                     pars = parameters)

  return(fit)
}

#' @title Extract parameter estimates from in-season model
#' @description This produces a formatted data frame of parameter estimates for the inseason STAN model.
#' @param inputs a list produced by running `prepare_inseason_inputs()`
#' @param fit a stanfit object produced by running `fit_inseason_model()`
#' @returns a table
#' @export
#' @md
extract_inseason_estimates <- function(inputs,
                                       fit) {

  # lookups for params
  # year-indexed variables are phi, lambda
  year_lookup <- tibble("year" = inputs$year_lookup) |>
    mutate(year_index = row_number())

  # week, year matrices is cp
  # nweeks is For_cp
  week_lookup <- tibble("week" = seq(1:inputs$inputs$data$Nwks)) |>
    mutate(week_index = row_number())

  # year, log lambda matrix is RTpars
  # logit phi, log lambda vector is muRT, sigmaRT, bCov
  par_lookup <- tibble("par_lookup" = c("logit_phi", "log_lambda"),
                       "par_index" = c(1, 2))

  summary_table <- rstan::summary(fit)$summary |>
    as.data.frame() |>
    tibble::rownames_to_column("parameter") |>
    # get indices to match to year, reach, etc. for parameters
    separate_wider_delim(parameter, delim = "[", names = c("parameter", "indices"),
                         too_few = "align_start") |>
    mutate(indices = str_remove_all(indices, "\\]")) |>
    separate_wider_delim(indices, ",", names = c("index_1", "index_2"),
                         too_few = "align_start") |>
    mutate(across(index_1:index_2, as.numeric)) |>
    # now match to lookups
    mutate(year_index = case_when(parameter %in% c("phi", "lambda", "RT_pars") ~ index_1,
                                  parameter == "cp" ~ index_2,
                                  TRUE ~ NA),
           week_index = ifelse(parameter %in% c("cp", "For_cp"), index_1, NA),
           par_index = case_when(parameter == "RTpars" ~ index_2,
                                 parameter %in% c("muRT", "sigmaRT", "bCov") ~ index_1,
                                 TRUE ~ NA)) |>
    # joins
    left_join(year_lookup, by = "year_index") |>
    left_join(week_lookup, by = "week_index") |>
    left_join(par_lookup, by = "par_index") |>
    mutate(parameter = ifelse(parameter %in% c("RTpars", "muRT", "sigmaRT", "bCov"), paste0(parameter, "_", par_lookup), parameter)) |>
    # remove indices
    select(-c(year_index, week_index, par_lookup, par_index, index_1, index_2)) |>
    rename(`2.5` = `2.5%`,
           `25` = `25%`,
           `50` = `50%`,
           `75` = `75%`,
           `97.5` = `97.5%`)

  if(any(summary_table$Rhat > 1.05)) {
    cli::cli_alert_warning("One or more parameters has an Rhat value over 1.05 :(")
  } else {
    cli::cli_bullets("No parameters have an Rhat value exceeding 1.05 :)")
  }

  summary_table_long <- summary_table |>
    pivot_longer(mean:Rhat, names_to = "statistic",
                 values_to = "value") |>
    mutate(model_name = inputs$model_name,
           site = inputs$site,
           week_fit = NA,
           stream = inputs$stream,
           srjpedata_version = as.character(packageVersion("SRJPEdata")))

  return(summary_table_long)

}
