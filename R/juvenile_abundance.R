#' Call BT-SPAS-X model
#' @details TODO
#' @param bt_spas_x_bayes_params: a list containing `number_mcmc`, `number_burnin`, `number_thin`,
#' and `number_chains`.
#' @param bt_spas_x_input_data
#' @param site
#' @param run_year
#' @param effort_adjust
#' @param multi_run_mode
#' @param mainstem_version
#' @param special_priors_data
#' @param bugs_directory
#' @returns TODO
#' @export
#' @md
run_bt_spas_x <- function(bt_spas_x_bayes_params,
                          bt_spas_x_input_data, site, run_year, effort_adjust, multi_run_mode,
                          mainstem_version, special_priors_data, bugs_directory) {
  if(!multi_run_mode) {
    run_single_bt_spas_x(bt_spas_x_bayes_params,
                         bt_spas_x_input_data, site, run_year, effort_adjust, mainstem_version,
                         special_priors_data, bugs_directory)
  } else if (multi_run_mode) {

    all_results <- list()

    for(i in unique(bt_spas_x_input_data$site)) {

      # TODO clean this up
      available_years <- bt_spas_x_input_data |>
        filter(site == i) |>
        distinct(run_year) |>
        arrange(run_year) |>
        pull(run_year)

      for(j in available_years)
        cli::cli_process_start(paste0("Running BT-SPAS-X for site ", i,
                                      " and run year ", j))

        all_results[[i]] <- run_single_bt_spas_x(bt_spas_x_bayes_params,
                                                 bt_spas_x_input_data,
                                                 site = i,
                                                 run_year = j,
                                                 effort_adjust, mainstem_version,
                                                 special_priors_data, bugs_directory)
    }
  }
}

#' Call BT-SPAS-X on a single site/run year combination
#' @details This function is called within `run_bt_spas_x()`
#' @param bt_spas_x_bayes_params: a list containing `number_mcmc`, `number_burnin`, `number_thin`,
#' and `number_chains`
#' @param bt_spas_x_input_data
#' @param site
#' @param run_year
#' @param effort_adjust
#' @param mainstem_version
#' @param special_priors_data
#' @param bugs_directory
#' @returns TODO
#' @export
#' @md
run_single_bt_spas_x <- function(bt_spas_x_bayes_params,
                                 bt_spas_x_input_data, site, run_year,
                                 effort_adjust = c(T, F), mainstem_version = c(F, T), special_priors_data,
                                 bugs_directory) {

  # filter datasets to match site, run_year, and weeks
  # catch_flow is average for julian week, standardized_efficiency_flow is average over recapture days (< 1 week)
  input_data <- bt_spas_x_input_data |>
    dplyr::filter(run_year == !!run_year,
                  site == !!site)

  if(nrow(input_data) == 0) {
    cli::cli_alert_warning(paste0("There is no catch data for site ", site,
                                  " and run year ", run_year, ". Please try with a different combination of site and year."))
    return(NULL)
  }

  # get numbers for looping in BUGs code - abundance model
  number_weeks_catch <- nrow(input_data) # for looping through the catch dataset
  # number_weeks_catch <- length(unique(input_data$week)) # number of weeks in the catch dataset
  indices_with_catch <- which(!is.na(input_data$count)) # indices of weeks with catch data
  number_weeks_with_catch <- length(indices_with_catch) # how many weeks actually have catch

  # analyze efficiency trials for all relevant sites (do not filter to site)

  # set up filter - if it's a tributary-based model, we cannot use efficiencies from KDL, TIS, RBDD
  if(!mainstem_version) {
    remove_sites <- c("knights landing", "tisdale", "red bluff diversion dam")
  } else if(mainstem_version) {
    remove_sites <- c("deer creek", "eye riffle", "live oak",
                      "okie dam", "mill creek", "yuba river", "herringer riffle", "ubc",
                      "lcc", "ucc", "hallwood", "steep riffle", "sunset pumps", "shawn's beach",
                      "gateway riffle", "lower feather river")
  }

  mark_recapture_data <- bt_spas_x_input_data |>
    dplyr::filter(!site %in% remove_sites,
           !is.na(standardized_efficiency_flow),
           !is.na(number_released),
           !is.na(number_recaptured))

  # get numbers for looping in BUGs code - pCap model
  # number of efficiency experiments completed
  number_efficiency_experiments <- nrow(mark_recapture_data) # unique(mark_recapture_data[c("site", "run_year", "week")]) |> nrow()
  # years where efficiency experiments were done
  years_with_efficiency_experiments <- unique(mark_recapture_data$run_year)
  # number of sites (for pCap calculations)
  number_sites_pCap <- length(unique(mark_recapture_data$site))
  # indices of those sites where efficiency trials were performed
  indices_sites_pCap <- which(unique(mark_recapture_data$site) == site)
  # indices of efficiency experiments in catch data
  indices_with_mark_recapture <- which(!is.na(input_data$number_released) &
                                         !is.na(input_data$standardized_efficiency_flow))
  # weeks (in catch data) where mark recapture were performed
  weeks_with_mark_recapture <- input_data$week[indices_with_mark_recapture]
  # indices (in mark-recap data) for the selected site and run year, filtered to weeks where mark-recap were performed (in catch data)
  indices_pCap <- which(mark_recapture_data$site == site &
                          mark_recapture_data$run_year == run_year &
                          mark_recapture_data$week %in% weeks_with_mark_recapture)
  # indices (in catch data) where no mark recap were performed
  indices_without_mark_recapture <- which(is.na(input_data$number_released) |
                                            is.na(input_data$standardized_efficiency_flow))
  # indices (in mark-recap data) for each site
  indices_site_mark_recapture <- mark_recapture_data |>
    group_by(site) |>
    mutate(ID = cur_group_id()) |>
    pull(ID)
  # number of weeks (in mark-recap data) where effiency experiments were performed
  number_weeks_with_mark_recapture <- length(indices_with_mark_recapture)
  # number of weeks (in mark-recap data) where effiency experiments were not performed
  number_weeks_without_mark_recapture <- length(indices_without_mark_recapture)

  # so BUGS doesn't bomb
  if(number_weeks_without_mark_recapture == 1) {
    indices_without_mark_recapture <- c(indices_without_mark_recapture, -99)
  }
  if(number_weeks_with_mark_recapture == 1) {
    indices_with_mark_recapture <- c(indices_with_mark_recapture, -99)
  }

  # set up b-spline basis matrix
  # this corresponds to line 148-153 in josh_original_model_code.R
  if(number_weeks_catch < 4) {
    cli::cli_alert_warning(paste0("There are fewer than 4 weeks with catch data for ",
                                  site, " and run year ", run_year, ". Spline parameters cannot function with fewer than 4 data points."))
    return(NULL)
  }

  # spline parameter calculation
  spline_data <- build_spline_data(number_weeks_catch, k_int = 4) # rule of thumb is 1 knot for every 4 data points for a cubic spline (which has 4 parameters)

  if(effort_adjust) {
    weekly_catch_data <- input_data$catch_standardized_by_hours_fished
  } else {
    weekly_catch_data <- input_data$catch_standardized_by_hours_fished
  }

  lgN_max <- log((input_data$catch_standardized_by_hours_fished/1000 + 1) / 0.025)
  lgN_max <- ifelse(is.na(lgN_max), mean(lgN_max, na.rm = T), lgN_max)

  # full data list
  full_data_list <- list("Nmr" = number_efficiency_experiments,
                         "Ntribs" = number_sites_pCap,
                         "ind_trib" = indices_site_mark_recapture,
                         "Releases" = mark_recapture_data$number_released,
                         "Recaptures" = mark_recapture_data$number_recaptured,
                         "Nstrata" = number_weeks_catch,
                         "u" = weekly_catch_data,
                         "K" = spline_data$K,
                         "ZP" = spline_data$b_spline_matrix,
                         "ind_pCap" = indices_pCap,
                         "Nwmr" = number_weeks_with_mark_recapture,
                         "Nwomr" = number_weeks_without_mark_recapture,
                         "Uind_wMR" = indices_with_mark_recapture,
                         "Uind_woMR" = indices_without_mark_recapture,
                         "use_trib" = indices_sites_pCap,
                         "Nstrata_wc" = number_weeks_with_catch,
                         "Uwc_ind" = indices_with_catch,
                         "mr_flow" = mark_recapture_data$standardized_efficiency_flow,
                         "catch_flow" = input_data$flow_cfs,
                         "lgN_max" = lgN_max)


  # set up models and data to pass to bugs
  number_experiments_at_site <- mark_recapture_data |>
    distinct(site, run_year, week) |>
    filter(site == !!site) |>
    nrow()

  # if efficiency trials occurred in the site
  if(number_experiments_at_site > 1) {
    if(number_weeks_without_mark_recapture == 0) {
      # all weeks have efficiency trials
      model_name <- "all_mark_recap.bug"
    } else {
      # some or all strata don't have efficiency trials
      if(number_weeks_with_mark_recapture > 0) {
        # some weeks have efficiency trials
        model_name <- "missing_mark_recap.bug"
      } else if(number_weeks_with_mark_recapture == 0) {
        # no weeks have efficiency trials
        model_name <- "no_mark_recap.bug"
      }
    }
  } else if(number_experiments_at_site == 0) { # no efficiency trials were performed at that site
    model_name <- "no_mark_recap_no_trib.bug"
  }

  data <- get_bt_spas_x_data_list(model_name, full_data_list)

  parameters <- c("trib_mu.P", "trib_sd.P", "flow_mu.P", "flow_sd.P", "pro_sd.P",
                   "b0_pCap", "b_flow", "pCap_U", "N", "Ntot", "sd.N", "sd.Ne")

  # initial parameter values
  ini_b0_pCap <- mark_recapture_data |>
    group_by(site) |>
    summarise(ini_b0_pCap = stats::qlogis(sum(number_recaptured) / sum(number_released))) |>
    pull(ini_b0_pCap)

  ini_lgN <- input_data |>
    mutate(ini_lgN = log(catch_standardized_by_hours_fished / 1000 + 2),
           ini_lgN = ifelse(is.na(ini_lgN), log(2 / 1000), ini_lgN)) |>
    pull(ini_lgN)

  pCap_mu_prior <- gtools::logit(sum(mark_recapture_data$number_recaptured) /
                                   sum(mark_recapture_data$number_released))

  init_list <- list(trib_mu.P = pCap_mu_prior,
                    b0_pCap = ini_b0_pCap,
                    flow_mu.P = 0,
                    b_flow = rep(0, number_sites_pCap),
                    trib_tau.P = 1,
                    flow_tau.P = 1,
                    pro_tau.P = 1,
                    b_sp = rep(1, spline_data$K),
                    lg_N = ini_lgN)

  #inits <- list(inits1 = init_list, inits2 = init_list, inits3 = init_list)
  inits <- list(init_list, init_list, init_list)

  # run the bugs model
  results <- bt_spas_x_bugs(data, inits, parameters, model_name, bt_spas_x_bayes_params, bugs_directory = paste0(bugs_directory))
  return(results)
}


#' Execute BT-SPAS-X in WinBUGs
#' @details This function is called within `run_single_bt_spas_x()` and calls the WinBUGS code
#' @param data a list containing the following elements:
#' * **Nmr** number of unique mark-recapture experiments performed across tributaries, years and weeks
#' * **Ntribs** number of tributaries to use for the pCap component of the model
#' * **Nstrata** number of weeks with catch data
#' * **Uwc_ind** number of weeks in catch data with associated pCap and flow data
#' * **Nwomr** number of weeks in catch data without associated pCap and flow data
#' * **Nstrata_wc** number of weeks in catch data with catch data (RST fished)
#' * **indices_site_pCap** tributary index (of all possible tributaries) to use to predict efficiency for missing strata. Note length=0 if selected tributary is not part of trib set that has mark-recap data. In this case model that samples from trib hyper will be called.
#' * **ind_trib** indices (1:Ntribs) assigned to the mark-recapture experiment table for use in BUGs
#' * **ind_pCap** indices of weeks in mark-recapture table for U strata being estimated
#' * **Uind_woMR** indices of weeks in catch data without associated pCap and flow data
#' * **Uind_wMR** indices of weeks in catch data with associated pCap and flow data
#' * **Uwc_ind** indices of weeks in catch data with catch data
#' * **Releases** number of fish released for each mark-recapture experiment
#' * **Recaptures** number of fish recaptured for each mark-recapture experiment
#' * **u** weekly abundance
#' * **mr_flow** standardized flow, averaged over recapture days (< 1 week)
#' * **catch_flow** standardized flow, averaged by week
#' * **K** number of columns in the bspline basis matrix
#' * **ZP** The b_spline_matrix (bspline basis matrix. One row for each data point (1:Nstrata), and one column for each term in the cubic polynomial function (4) + number of knots)
#' * **lgN_max** maxmimum possible value for log N across strata
#' @param inits
#' @param parameters
#' @param model_name
#' @param bt_spas_x_bayes_params: a list containing `number_mcmc`, `number_burnin`, `number_thin`,
#' and `number_chains`.
#' @param bugs_directory
#' @returns either a list of the required inputs for a WinBUGS model (if running on a Mac), or the results
#' of the model run.
#' @export
#' @md
bt_spas_x_bugs <- function(data, inits, parameters, model_name, bt_spas_x_bayes_params,
                           number_mcmc, bugs_directory) {

  # set model name directory
  model_name_full <- here::here("model_files", model_name)
 # model_name_full <- paste0("model_files/", model_name)

  # get operating system - bugs can't run on a mac without serious set-up
  operating_system <- ifelse(grepl("Mac", Sys.info()['nodename']) | grepl("MBP", Sys.info()['nodename']), "mac", "pc")
  if(operating_system == "mac") {
    # TODO update this message to be a [Y/n] and/or link to wine emulator
    cli::cli_alert_warning("This model is currently coded in WinBUGS, which cannot easily be run on a Mac.
                           All the information required to run the model will be returned, but the model will not be run.")

    return(list(data = data, inits = inits, parameters = parameters, model_name = model_name_full,
                n.chains = bt_spas_x_bayes_params$number_chains,
                n.burnin = bt_spas_x_bayes_params$number_burnin,
                n.thin = bt_spas_x_bayes_params$number_thin,
                n.iter = bt_spas_x_bayes_params$number_mcmc, debug = FALSE, codaPkg = FALSE, DIC = TRUE,
                clearWD = TRUE, bugs_directory = bugs_directory))
  } else {

    cli::cli_process_start("WinBUGS model running")
    # run bugs model
    model_results <- R2WinBUGS::bugs(data, inits, parameters, model.file = model_name_full,
                                     n.chains = bt_spas_x_bayes_params$number_chains,
                                     n.burnin = bt_spas_x_bayes_params$number_burnin,
                                     n.thin = bt_spas_x_bayes_params$number_thin,
                                     n.iter = bt_spas_x_bayes_params$number_mcmc,
                                     debug = TRUE, codaPkg = FALSE, DIC = TRUE, clearWD = TRUE,
                                     bugs.directory = bugs_directory)

    posterior_output <- model_results$sims.list
    summary_output <- round(model_results$summary, 3)
    dic_output <- c(model_results$pD, model_results$DIC)
    knots_output <- knot_positions

    return(list("posterior_output" = posterior_output,
                "summary_output" = summary_output,
                "dic_output" = dic_output,
                "knots_output" = knots_output))
  }
}

#' Prepare Data Object for BT-SPAS-X WinBUGs call
#' @details This function is called within `run_single_bt_spas_x()` and prepares the data lists
#' based on which model will be run.
#' @param model_name
#' @param full_data_list
#' @returns a named list with the required elements for that model run.
#' @export
#' @md
get_bt_spas_x_data_list <- function(model_name, full_data_list) {
  if(model_name == "all_mark_recap.bug") {
    data_needed <- c("Nmr", "Ntribs" ,"ind_trib", "Releases", "Recaptures", "Nstrata",
                     "u", "K", "ZP", "ind_pCap", "Nstrata_wc", "Uwc_ind", "mr_flow", "lgN_max")
  } else if(model_name == "missing_mark_recap.bug") {
    data_needed <- c("Nmr", "Ntribs", "ind_trib", "Releases", "Recaptures", "Nstrata",
                     "u", "K", "ZP", "ind_pCap", "Nwmr", "Nwomr",
                     "Uind_wMR", "Uind_woMR", "use_trib", "Nstrata_wc", "Uwc_ind",
                     "mr_flow", "catch_flow", "lgN_max")
  } else if(model_name == "no_mark_recap.bug") {
    data_needed <- c("Nmr", "Ntribs", "ind_trib", "Releases", "Recaptures", "Nstrata",
                     "u", "K", "ZP", "Nwomr", "Uind_woMR", "use_trib",
                     "Nstrata_wc", "Uwc_ind", "mr_flow", "catch_flow", "lgN_max")
  } else if(model_name == "no_mark_recap_no_trib.bug") {
    data_needed <- c("Nmr", "Ntribs", "ind_trib", "Releases", "Recaptures", "Nstrata",
                     "u", "K", "ZP", "Nwomr", "Uind_woMR", "Nstrata_wc",
                     "Uwc_ind", "mr_flow", "catch_flow", "lgN_max")
  }
  new_data_list <- full_data_list[sapply(names(full_data_list), function(x) x %in% data_needed)]
  return(new_data_list)
}


#' Prepare spline parameters object for BT-SPAS-X WinBUGs call
#' @details This function is called within `run_single_bt_spas_x()` and prepares spline parameters
#' to pass to the data list
#' @param number_weeks_catch
#' @param k_int
#' @returns a named list containing *K* and *b_spline_matrix* for passing to WinBUGS
#' @export
#' @md
build_spline_data <- function(number_weeks_catch, k_int) {

  number_knots <- round(number_weeks_catch / k_int, 0)
  first_knot_position <- 2
  final_knot_position <- number_weeks_catch - 1 # keep first and/or last knot positions away from tails if there are intervals with no sampling on the tails
  knot_positions <- seq(first_knot_position, final_knot_position, length.out = number_knots) # define position of b-spline knots using even interval if no missing data
  b_spline_matrix <- splines2::bSpline(x = 1:number_weeks_catch, knots = knot_positions, deg = 3, intercept = T) # bspline basis matrix. One row for each data point (1:number_weeks_catch), and one column for each term in the cubic polynomial function (4) + number of knots
  K <- ncol(b_spline_matrix)

  return(list("K" = K,
              "b_spline_matrix" = b_spline_matrix))
}
