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
prepare_stock_recruit_inputs <- function(year, stream, adult_data_type,
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

  # TODO pull from model storage
  juv <- readRDS("~/Downloads/all_JPE_sites_clean.rds") |>
    mutate(brood_year = run_year - 1) |>
    pivot_wider(names_from = "statistic",
                values_from = "value") |>
    filter(parameter == "Ntot") |>
    mutate(mean_juvenile_abundance = mean,
           cv_juvenile_abundance = sd/mean) |>
    select(brood_year, site, mean_juvenile_abundance, cv_juvenile_abundance) |>
    right_join(years_filter, by = c("brood_year", "site")) |>
    glimpse()

  data_with_adult_juv_and_covars <- full_join(adult_data_and_covariates, juv) |>
    # TODO filter by adult data type
    glimpse()

  # TODO now filter by covariate args according to josh's code
  # TODO produce data and inits

}


