#' Prepare inputs for stock-recruit model
#' @details This function prepares data for input into a stock recruit model.
#' @param adult_data_type the type of survey the adult data came from. Either `upstream_estimate`,
#' `redd_count`, `holding_count`, or `carcass_estimate`.
#' @site The RST site from which you want to use juvenile abundance estimates.
#' @covariate The covariate you would like to use to fit the model. Either `gdd_spawn`,
#' `above_13_temp_day`, `above_13_temp_week`, `weekly_max_temp_max`, `weekly_max_temp_mean`, `weekly_max_temp_median`,
#' `mean_flow`, `max_flow`, `min_flow`, `median_flow`, or `null`
#' @covariate_lifestage The corresponding lifestage for the covariate: either `spawning and incubation` or `rearing`.
#' @returns a list of tables:
#' * **stock_recruit_table** juvenile abundance and adult abundance by brood year.
#' * **full_covariate_table** full table of covariates associated with the site, stream, and brood years in the stock recruit table
#' * **truncated_covariate_table** table of covariates limited to only those years where all covariates are available (for covariate comparison analyses)
#' @export
#' @md
prepare_stock_recruit_inputs <- function(adult_data_type, site, covariate,
                                         covariate_lifestage) {

  adult_stream <- SRJPEdata::site_lookup |>
    filter(site == !!site) |>
    distinct(stream) |>
    pull()

  minimum_brood_year <- 2001 # TODO confirm purpose of this with josh
  all_covariate_options <- c("si_gdd_spawn", "si_above_13_temp_day", "si_above_13_temp_week",
                             "si_weekly_max_temp_max", "si_weekly_max_temp_mean", "si_weekly_max_temp_median",
                             "si_mean_flow", "si_max_flow", "si_min_flow", "si_median_flow",
                             "rr_mean_flow", "rr_max_flow", "rr_min_flow", "rr_median_flow")

  # adult data
  adult_data <- SRJPEdata::observed_adult_input |>
    filter(stream == adult_stream,
           data_type == adult_data_type) |>
    rename(adult_year = year,
           adult_count = count,
           adult_data_type = data_type)

  # juvenile abundance estimates
  # TODO build out workflow for pulling total abundance fits
  # TODO this should be after PLAD application
  juvenile_data <- readRDS("~/Downloads/all_JPE_sites_clean.rds") |>
    filter(parameter == "Ntot",
           is.na(error)) |>
    mutate(parameter = "total_abundance") |>
    pivot_wider(names_from = "statistic",
                values_from = "value") |>
    filter(site == site) |>
    select(-c(model_name, week_fit, location_fit, srjpedata_version, error,
              Rhat)) |>
    mutate(adult_year = run_year - 1) |>
    left_join(SRJPEdata::site_lookup |>
                select(site, stream) |>
                distinct(site, stream))

  combined_data <- full_join(juvenile_data, adult_data,
                             by = c("stream", "adult_year")) |>
    mutate(both_have_data = ifelse(!is.na(`50`) & !is.na(adult_count), TRUE, FALSE),
           # TODO check calculation of mean and cv of juvenile abundance based on PLAD output
           cv_juvenile_abundance = sd/mean) |>
    filter(both_have_data) |>
    select(-both_have_data) |>
    rename(brood_year = adult_year) |>
    select(brood_year, adult_count, mean_juvenile_abundance = mean,
           cv_juvenile_abundance)

  # TODO Confirm that covariates are summarized at stream level? if not, we will have to
  # get a site lookup??
  # TODO stop here - causing lists for non-distinct values. need to figure out how to match
  # covariates with stock recruit table
  covariates <- SRJPEdata::stock_recruit_covariates |>
    ungroup() |>
    filter(stream == adult_stream, #TODO do we filter to site too, esp for feather?
           # TODO confirm year in covariate table is a brood year
           year %in% combined_data$brood_year,
           year >= minimum_brood_year,
           # TODO fix -Inf values in table
           is.finite(value)) |>
    suppressWarnings() |>
    mutate(variable = ifelse(lifestage == "spawning and incubation", "si", "rr"),
           variable = paste(variable, covariate_structure, sep = "_")) |>
    select(year, variable, value) |>
    pivot_wider(names_from = variable,
                values_from = value) |>
    mutate(null = "0") |>
    pivot_longer(-year,
                 values_to = "value",
                 names_to = "variable")

  truncated_covariates <- covariates |>
    distinct(year, variable, value) |>
    pivot_wider(names_from = variable,
                values_from = value) |>
    # truncate
    drop_na(all_of(all_covariate_options)) |>
    pivot_longer(-year,
                 names_to = "variable",
                 values_to = "value")

  # data tables

  return(list("stock_recruit_table" = combined_data,
              "full_covariate_table" = covariates,
              "truncated_covariate_table" = truncated_covariates))

  # now prep data inputs based on what model you want to run
  if(covariate == "null") {
    covariate_filter <- "null"
  } else {
    ls <- ifelse(covariate_lifestage == "spawning and incubation", "si", "rr")
    covariate_filter <- paste(ls, covariate, sep = "_")
  }

}


