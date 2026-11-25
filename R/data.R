#' BT-SPAS-X MCMC settings
#'
#' @details MCMC settings used to fit the BT-SPAS-X models.
#' @format A named list:
#' * **number_mcmc** Numeric. Total iterations per chain for the BUGS abundance model.
#' * **number_burnin** Numeric. Iterations discarded as burn-in.
#' * **number_thin** Numeric. Thinning interval.
#' * **number_chains** Numeric. Number of chains, used by both the pCap and abundance models.
#' @source Built in `data-raw/cache-params.R`.
"bt_spas_x_bayes_params"

#' BT-SPAS-X model code
#'
#' @details Stan code for the pCap models and BUGS code for the abundance
#'   model, read from `model_files/`. `fit_pCap_model()` picks the pCap model
#'   by name.
#' @format A named list of character strings:
#' * **pCap_all_sites** Stan code for the hierarchical tributary pCap model.
#' * **pCap_one_site** Stan code for the single-site (mainstem) pCap model.
#' * **pCap_one_site_skew** As `pCap_one_site`, with skew-normal process error.
#' * **abundance_BUGS** WinBUGS code for the abundance model.
#' @source Built in `data-raw/cache-models.R`.
"bt_spas_x_model_code"

#' Julian week to date lookup
#'
#' @format A tibble with 53 rows:
#' * **Jwk** Numeric. Julian week, 1-53.
#' * **date** Character. First day of the week, as `"Mon-DD"` (e.g. `"Jan-01"`).
#' @source `data-raw/juvenile_abundance/Jwk_Dates.txt`, built in `data-raw/cache-models.R`.
"julian_week_to_date_lookup"

#' Site order, north to south
#'
#' @details Used to sort sites from north to south in model inputs and plots.
#' @format A tibble with 15 rows:
#' * **ns_order** Integer. Position from north (1) to south.
#' * **site** Character. Site name, as used in SRJPEdata.
#' @source Built in `data-raw/cache-params.R`.
"site_order_north_south"
