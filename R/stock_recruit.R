#' Prepare inputs for stock-recruit model
#' @details This function prepares data for input into a stock recruit model.
#' @param adult_data_type the type of survey the adult data came from. Either `upstream_estimate`,
#' `redd_count`, `holding_count`, or `carcass_estimate`.
#' @param adult_stream stream from which the adult data comes from.
#' @param minimum_year to filter covariates # TODO confirm with Josh.
#' @returns a list of tables:
#' * **stock_recruit_table** juvenile abundance and adult abundance by brood year.
#' * **full_covariate_table** full table of covariates associated with the site, stream, and brood years in the stock recruit table
#' * **truncated_covariate_table** table of covariates limited to only those years where all covariates are available (for covariate comparison analyses)
#' @export
#' @md
prep_stock_recruit_inputs <- function(adult_data_type, site, adult_stream,
                                      minimum_year) {
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
           cv_juvenile_abundance = sd/mean) |>
    filter(both_have_data) |>
    select(-both_have_data) |>

  # TODO Confirm site lookup is correct
  covariates <- SRJPEdata::stock_recruit_covariates |>
    ungroup() |>
    mutate(site_lookup = tolower(gage_number)) |>
    # get site lookup
    # left_join(SRJPEdata::site_lookup |>
    #             distinct(site, site_group),
    #           by = "site_group") |>
    filter(stream == adult_stream,
           # TODO decide on how to join/filter sites, is that correct?
           site_lookup == site,
           # TODO confirm year in covariate table is a brood year
           year %in% combined_data$brood_year,
           year > minimum_year_arg,
           # TODO fix -Inf values in table
           is.finite(value)) |>
    suppressWarnings() |>
    mutate(variable = ifelse(lifestage == "spawning and incubation", "si", "rr"),
           variable = paste(variable, covariate_structure, sep = "_")) |>
    select(year, variable, value)

  truncated_covariates <- covariates |>
    distinct(year, variable, value) |>
    pivot_wider(names_from = variable,
                values_from = value) |>
    # truncate
    drop_na(si_gdd_spawn, si_above_13_temp_day, si_above_13_temp_week,
            si_weekly_max_temp_max, si_weekly_max_temp_mean,
            si_weekly_max_temp_median) |>
    pivot_longer(c(si_gdd_spawn, si_above_13_temp_day, si_above_13_temp_week,
                   si_weekly_max_temp_max, si_weekly_max_temp_mean,
                   si_weekly_max_temp_median),
                 names_to = "variable",
                 values_to = "value")

  return(list("stock_recruit_table" = combined_data,
              "full_covariate_table" = covariates,
              "truncated_covariate_table" = truncated_covariates))
}


