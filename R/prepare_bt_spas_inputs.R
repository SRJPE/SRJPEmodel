#' Prepare inputs for pCap STAN model
#' @details This function prepares data for input into a pCap STAN model.
#' @param mainstem Whether or not you want to evaluate efficiency trials for a mainstem site
#' (`knights landing`, `tisdale`, and `red bluff diversion dam`) or for a tributary site. If `FALSE`,
#' the mark recapture dataset will be filtered to exclude those mainstem sites.
#' @returns a list:
#' * **pCap_inputs** a list of data and inits for input into the pCap model.
#' * **sites_fit** a list of site names associated with `ind_trib`.
#' @export
#' @md
prepare_pCap_inputs <- function(input_catch_data = NULL, input_efficiency_data = NULL,
                                mainstem = c(FALSE, TRUE)) {

  if(is.null(input_catch_data)) {
    input_catch_data <- SRJPEdata::weekly_juvenile_abundance_catch_data
  }

  if(is.null(input_efficiency_data)) {
    input_efficiency_data <- SRJPEdata::weekly_juvenile_abundance_efficiency_data
  }

  cli::cli_bullets("Producing data for pCap model from efficiency dataset")

  # analyze efficiency trials for all relevant sites (do not filter to site)
  # set up filter - if it's a tributary-based model, we cannot use efficiencies from KDL, TIS, RBDD
  if(!mainstem) {
    remove_sites <- c("knights landing", "tisdale", "red bluff diversion dam")
  } else {
    remove_sites <- input_efficiency_data |>
      filter(!site %in% c("knights landing", "tisdale", "red bluff diversion dam")) |>
      distinct(site) |>
      pull(site)
  }

  # prepare "mark recapture" dataset - all mark-recap trials in the system
  mark_recapture_data <- input_efficiency_data |>
    # grab standardized_flow
    left_join(input_catch_data |>
                select(year, week, stream, site, run_year, standardized_flow),
              by = c("year", "week", "run_year", "stream", "site")) |>
    # or do we want to filter just no number released?
    dplyr::filter(!site %in% remove_sites &
                    !is.na(standardized_flow),
                  !is.na(number_released) &
                    !is.na(number_recaptured)) |>
    # right now there's lifestage in the dataset, so we have to do distinct()
    distinct(site, run_year, week, number_released, number_recaptured, .keep_all = TRUE)


  # get use_trib, sites_fit, and ind_trib indexing
  # first assign 1:Ntribs to the unique sites in the dataset
  site_lookup <- mark_recapture_data |>
    distinct(site) |>
    mutate(ID = row_number())

  # pull the order of those sites
  sites_fit <- site_lookup |>
    pull(site)

  # assign the IDs to the sites in the mark-recapture dataset
  indices_site_mark_recapture <- mark_recapture_data |>
    left_join(site_lookup, by = "site") |>
    pull(ID)

  # get indexing for "mark recap" dataset (pCap model)
  Ntribs <- length(sites_fit) # number of sites (for pCap calculations)
  number_efficiency_experiments <- unique(mark_recapture_data[c("site", "run_year", "week")]) |>
    nrow() # number of efficiency experiments completed, nrow(mark_recapture_data) this depends on whether you have lifestage or not

  # build data list with ALL elements
  data <- list("Nmr" = number_efficiency_experiments,
               "Ntribs" = Ntribs,
               "ind_trib" = indices_site_mark_recapture,
               "Releases" = mark_recapture_data$number_released,
               "Recaptures" = mark_recapture_data$number_recaptured,
               "mr_flow" = mark_recapture_data$standardized_flow)


  # check data list for NaNs and Infs
  cli::cli_process_start("Checking data inputs",
                         msg_done = "Data checked",
                         msg_failed = "Data check failed")
  invisible(lapply(names(data), function(x) {
    if(any(is.nan(data[[x]])) | any(is.infinite(data[[x]]))) {
      cli::cli_abort(paste0("NaNs detected in ", x, ". Please check your input data."))
    }
  }))
  cli::cli_process_done()

  pCap_parameters <-  c("trib_mu.P", "trib_sd.P", "flow_mu.P", "flow_sd.P", "pro_sd.P",
                        "b0_pCap", "b_flow")

  # initial parameter values
  ini_b0_pCap <- rep(NA, Ntribs)
  for(i in 1:Ntribs) {
    irows = which(indices_site_mark_recapture == i)
    ini_b0_pCap[i] = qlogis(sum(mark_recapture_data$number_recaptured[irows]) /
                              sum(mark_recapture_data$number_released[irows]))
    if(is.nan(ini_b0_pCap[i]) | is.infinite(ini_b0_pCap[i])) {
      # -Inf happens when number recaptured == 0, logit of 0 is -Inf
      ini_b0_pCap[i] <- -5
    }
  }

  pCap_mu_prior <- qlogis(sum(mark_recapture_data$number_recaptured) /
                            sum(mark_recapture_data$number_released))

  init_list <- list(trib_mu.P = pCap_mu_prior,
                    b0_pCap = ini_b0_pCap,
                    flow_mu.P = 0,
                    b_flow = rep(0, Ntribs),
                    trib_tau.P = 1,
                    flow_tau.P = 1,
                    pro_tau.P = 1)


  cli::cli_process_start("Checking init inputs for pCap model",
                         msg_done = "Inits checked",
                         msg_failed = "Init check failed")
  invisible(lapply(names(init_list), function(x) {
    if(any(is.nan(init_list[[x]])) | any(is.infinite(init_list[[x]]))) {
      cli::cli_abort(paste0("NaNs detected in ", x, ". Please check your input data."))
    }
  }))
  cli::cli_process_done()

  inits <- list(init_list, init_list, init_list)

  # inputs for pCap model
  inputs <- list(data = data,
                 inits = inits,
                 parameters = pCap_parameters)

  return(list("inputs" = inputs,
              "sites_fit" = sites_fit))

}

#' Prepare inputs for abundance model
#' @details This function prepares data for input into an abundance model.
#' @param site site for which you want to fit the model
#' @param run_year run year for which you want to fit the model
#' @param effort_adjust whether or not you want to use catch adjusted by effort.
#' @returns a list:
#' * **abundance_inputs** a list of data and inits for input into the abundance model.
#' * **model_name** The version of the pCap logic to use for generating lt_pCap_U values. Either
#' `all_mark_recap`, `missing_mark_recap`, `no_mark_recap`, or `no_mark_recap_no_trib`.
#' * **weeks_fit** The weeks fit for the abundance model
#' * **weeks_date** Associated dates for the weeks fit for the abundance model.
#' @export
#' @md
prepare_abundance_inputs <- function(input_catch_data, input_efficiency_data,
                                     site, run_year,
                                     effort_adjust = c(T, F)) {

  if(is.null(input_catch_data)) {
    input_catch_data <- SRJPEdata::weekly_juvenile_abundance_catch_data
  }

  if(is.null(input_efficiency_data)) {
    input_efficiency_data <- SRJPEdata::weekly_juvenile_abundance_efficiency_data
  }

  cli::cli_bullets("Filtering catch data to site, run year, and weeks")
  catch_data <- input_catch_data |>
    filter(life_stage != "yearling") |>
    filter(run_year == !!run_year,
           site == !!site,
           week %in% c(seq(45, 53), seq(1, 22))) |>
    group_by(year, week, stream, site, run_year) |>
    # keep NAs in count columns
    summarise(count = if(all(is.na(count))) NA_real_ else sum(count, na.rm = TRUE),
              mean_fork_length = mean(mean_fork_length, na.rm = T),
              hours_fished = mean(hours_fished, na.rm = T),
              flow_cfs = mean(flow_cfs, na.rm = T),
              average_stream_hours_fished = mean(average_stream_hours_fished, na.rm = T),
              standardized_flow = mean(standardized_flow, na.rm = T),
              catch_standardized_by_hours_fished = if(all(is.na(catch_standardized_by_hours_fished))) NA_real_ else sum(catch_standardized_by_hours_fished, na.rm = TRUE),
              lgN_prior = mean(lgN_prior, na.rm = T)) |>
    ungroup() |>
    mutate(count = round(count, 0),
           catch_standardized_by_hours_fished = round(catch_standardized_by_hours_fished, 0),
           # change all NaNs to NAs
           across(mean_fork_length:lgN_prior, ~ifelse(is.nan(.x), NA, .x)))

  # flag if no data are available
  if(nrow(catch_data) == 0) {
    cli::cli_abort(paste0("There is no catch data for site ", site,
                          " and run year ", run_year, ". Please try with a different combination of site and year."))
  }

  # get numbers for looping in BUGs code - abundance model
  number_weeks_catch <- nrow(catch_data) # for looping through the catch dataset
  indices_with_catch <- which(!is.na(catch_data$count)) # indices of weeks with catch data
  number_weeks_with_catch <- length(indices_with_catch) # how many weeks actually have catch

  # analyze efficiency trials for all relevant sites (do not filter to site)
  # set up filter - if it's a tributary-based model, we cannot use efficiencies from KDL, TIS, RBDD
  if(!site %in% c("knights landing", "tisdale", "red bluff diversion dam")) {
    remove_sites <- c("knights landing", "tisdale", "red bluff diversion dam")
  } else {
    remove_sites <- input_efficiency_data |>
      filter(!site %in% c("knights landing", "tisdale", "red bluff diversion dam")) |>
      distinct(site) |>
      pull(site)
  }

  # prepare "mark recapture" dataset - all mark-recap trials in the system
  mark_recapture_data <- input_efficiency_data |>
    # grab standardized_flow
    left_join(input_catch_data |>
                select(year, week, stream, site, run_year, standardized_flow),
              by = c("year", "week", "run_year", "stream", "site")) |>
    # or do we want to filter just no number released?
    dplyr::filter(!site %in% remove_sites &
                    !is.na(standardized_flow),
                  !is.na(number_released) &
                    !is.na(number_recaptured)) |>
    # right now there's lifestage in the dataset, so we have to do distinct()
    distinct(site, run_year, week, number_released, number_recaptured, .keep_all = TRUE)



  # bring together efficiency and catch data so that we can get the indices of
  # catch data (hence left join) that correspond to certain efficiency trial
  # information.
  all_data_for_indexing <- left_join(catch_data, mark_recapture_data)

  # get use_trib, sites_fit, and ind_trib indexing
  # first assign 1:Ntribs to the unique sites in the dataset
  site_lookup <- mark_recapture_data |>
    distinct(site) |>
    mutate(ID = row_number())

  # pull the order of those sites
  sites_fit <- site_lookup |>
    pull(site)

  # assign the IDs to the sites in the mark-recapture dataset
  indices_site_mark_recapture <- mark_recapture_data |>
    left_join(site_lookup, by = "site") |>
    pull(ID)

  # indices for generating lt_pCap_Us
  indices_with_mark_recapture <- which(!is.na(all_data_for_indexing$number_released) &
                                         !is.na(all_data_for_indexing$standardized_efficiency_flow)) # indices of efficiency experiments in catch data
  indices_without_mark_recapture <- which(is.na(all_data_for_indexing$number_released) |
                                            is.na(all_data_for_indexing$standardized_efficiency_flow))   # indices (in catch data) where no mark recap were performed
  weeks_with_mark_recapture <- all_data_for_indexing$week[indices_with_mark_recapture] # weeks (in catch data) where mark recapture were performed
  number_weeks_with_mark_recapture <- length(indices_with_mark_recapture) # number of weeks (in mark-recap data) where effiency experiments were performed
  number_weeks_without_mark_recapture <- length(indices_without_mark_recapture)   # number of weeks (in mark-recap data) where effiency experiments were not performed
  indices_pCap <- which(mark_recapture_data$site == site &
                          mark_recapture_data$run_year == run_year &
                          mark_recapture_data$week %in% weeks_with_mark_recapture)   # indices (in mark-recap data) for the selected site and run year, filtered to weeks where mark-recap were performed (in catch data)
  indices_sites_pCap <- which(sites_fit == site) # indices of those sites where efficiency trials were performed, can be length = 0

  # set up b-spline basis matrix
  # this corresponds to line 148-153 in josh_original_model_code.R
  if(number_weeks_catch < 4) {
    cli::cli_abort(paste0("There are fewer than 4 weeks with catch data for ",
                          site, " and run year ", run_year, ". Spline parameters cannot function with fewer than 4 data points."))
  }

  # spline parameter calculation
  spline_data <- SRJPEmodel::build_spline_data(number_weeks_catch, k_int = 4) # rule of thumb is 1 knot for every 4 data points for a cubic spline (which has 4 parameters)

  if(effort_adjust) {
    weekly_catch_data <- catch_data$catch_standardized_by_hours_fished[indices_with_catch]
  } else {
    weekly_catch_data <- catch_data$count[indices_with_catch]
  }

  # set lgN priors, using josh's code from 12-11-2024
  # data input
  lgN_max = rep(log(0.001 * (mean(weekly_catch_data, na.rm=T) + 1) / 0.01), number_weeks_catch)

  for(j in 1:number_weeks_with_catch){
    if(is.na(weekly_catch_data[j]) == F) lgN_max[indices_with_catch[j]] = log(0.001 * (weekly_catch_data[j] + 1) / 0.01)
  }

  # initial value for lgN input
  ini_lgN = rep(log(0.001 * (min(weekly_catch_data) + 1) / 0.025), number_weeks_catch)
  for(i in 1:number_weeks_with_catch) ini_lgN[indices_with_catch[i]] = log(0.001 * (weekly_catch_data[i] + 1) / 0.025)

  # build data list
  data <- list("Nstrata" = number_weeks_catch,
               "Nstrata_wc" = number_weeks_with_catch,
               "u" = weekly_catch_data,
               "K" = spline_data$K,
               "ZP" = spline_data$b_spline_matrix,
               "Uwc_ind" = indices_with_catch,
               "lgN_max" = lgN_max)

  # data needed for generating lt_pCap_Us
  lt_pCap_U_data <- list("catch_flow" = catch_data$standardized_flow,
                         "use_trib" = indices_sites_pCap,
                         "Nwmr" = number_weeks_with_mark_recapture,
                         "Nwomr" = number_weeks_without_mark_recapture,
                         "Uind_wMR" = indices_with_mark_recapture,
                         "Uind_woMR" = indices_without_mark_recapture,
                         "ind_pCap" = indices_pCap)

  # use number of experiments at site to determine which model to call
  number_experiments_at_site <- mark_recapture_data |>
    distinct(site, run_year, week) |>
    filter(site == !!site) |>
    nrow()

  if(number_experiments_at_site == 1) {
    # TODO issue will be the inter-trial variance, and maybe we need a better
    # TODO informative prior
    cli::cli_abort("There is only one experiment at the selected site, and no model
                   is available yet for that case")
  }
  # if efficiency trials occurred in the site
  if(number_experiments_at_site > 1) {
    if(number_weeks_without_mark_recapture == 0) {
      # all weeks have efficiency trials
      model_name <- "all_mark_recap"
    } else {
      # some or all strata don't have efficiency trials
      if(number_weeks_with_mark_recapture > 0) {
        # some weeks have efficiency trials
        model_name <- "missing_mark_recap"
      } else if(number_weeks_with_mark_recapture == 0) {
        # no weeks have efficiency trials
        model_name <- "no_mark_recap"
      }
    }
  } else if(number_experiments_at_site == 0) { # no efficiency trials were performed at that site
    model_name <- "no_mark_recap_no_trib"
  }

  # check data list for NaNs and Infs
  cli::cli_process_start("Checking data inputs",
                         msg_done = "Data checked",
                         msg_failed = "Data check failed")
  invisible(lapply(names(data), function(x) {
    if(any(is.nan(data[[x]])) | any(is.infinite(data[[x]]))) {
      cli::cli_abort(paste0("NaNs detected in ", x, ". Please check your input data."))
    }
  }))
  cli::cli_process_done()

  parameters <- c("tau_N", "tau_Ne", "b_sp", "lg_N", "lt_pCap_U",
                  "N", "Ntot", "lg_CumN")

  # inits
  init_list <- list(b_sp = rep(1, spline_data$K),
                    lg_N = ini_lgN)


  cli::cli_process_start("Checking init inputs for abundance model",
                         msg_done = "Inits checked",
                         msg_failed = "Init check failed")
  invisible(lapply(names(init_list), function(x) {
    if(any(is.nan(init_list[[x]])) | any(is.infinite(init_list[[x]]))) {
      cli::cli_abort(paste0("NaNs detected in ", x, ". Please check your input data."))
    }
  }))
  cli::cli_process_done()

  inits <- list(init_list, init_list, init_list)

  inputs_for_abundance <- list(data = data,
                               inits = inits,
                               parameters = parameters)

  weeks_fit <- tibble("Jwk" = catch_data$week) |>
    left_join(SRJPEmodel::julian_week_to_date_lookup, by = "Jwk")

  return(list("inputs" = inputs_for_abundance,
              "lt_pCap_U_data" = lt_pCap_U_data,
              "model_name" = model_name,
              "catch_flow_raw" = catch_data$flow_cfs,
              "mr_flow_raw" = mark_recapture_data$flow_cfs,
              "weeks_fit" = weeks_fit$Jwk,
              "week_date" = weeks_fit$date,
              "sites_fit" = sites_fit))
}
