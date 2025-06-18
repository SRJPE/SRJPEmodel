#' Prepare inputs for within season model
#' @details This function prepares data for input into a within season model.
#' @param con A valid connection to the model run database.
#' @param stream The stream for which you want to fit the model.
#' @param site The site for which you want to fit the model.
#' @returns a list of inputs TBD
#' @export
#' @md
prepare_within_season_inputs <- function(con, stream, site) {

  if(!DBI::dbIsValid(con)) {
    cli::cli_abort("Connection argument does not have a valid connection to the database. Please try reconnecting to the database using DBI::dbConnect")
  }

  # get connection
  cfg <- config::get()
  con <- DBI::dbConnect(RPostgres::Postgres(),
                        dbname = cfg$db_name,
                        host = cfg$db_host,
                        port = cfg$db_port,
                        user = cfg$db_user,
                        password = cfg$db_password)


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

  data <- list(year = model_run_ids$year,
               number_estimated_weeks = number_estimated_weeks,
               annual_abundance_cv = annual_abundance$cv,
               abundance_through_week_matrix = cv_abundance_through_week_transposed)

  return(list(inputs = list(data = data)))
}
