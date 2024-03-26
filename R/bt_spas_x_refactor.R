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
        arrange() |>
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
    # filter(run_year == !!run_year,
    #        site == site_selection)
    dplyr::filter(run_year == !!run_year,
                  site == !!site)

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
           !is.na(standardized_efficiency_flow), # TODO should stay at flow_cfs?
           !is.na(number_released),
           !is.na(number_recaptured))

  # get numbers for looping in BUGs code - pCap model
  # TODO add stop for if there are no mark recapture experiments in the stream
  number_efficiency_experiments <- unique(mark_recapture_data[c("site", "run_year", "week")]) |>
    nrow()
  years_with_efficiency_experiments <- unique(mark_recapture_data$run_year)
  number_sites_pCap <- length(unique(mark_recapture_data$site))
  indices_sites_pCap <- which(unique(mark_recapture_data$site) == site) # TODO check this
  indices_pCap <- which(mark_recapture_data$site == site &
                          mark_recapture_data$run_year == run_year)
  indices_with_mark_recapture <- which(!is.na(input_data$number_released) &
                                         !is.na(input_data$standardized_efficiency_flow))
  indices_without_mark_recapture <- which(is.na(input_data$number_released) |
                                            is.na(input_data$standardized_efficiency_flow))
  # TODO check this
  indices_site_mark_recapture <- mark_recapture_data |>
    distinct(site, run_year, week) |>
    group_by(site) |>
    mutate(ID = cur_group_id()) |>
    pull(ID)
  number_weeks_with_mark_recapture <- length(indices_with_mark_recapture)
  number_weeks_without_mark_recapture <- length(indices_without_mark_recapture)

  # subset mark recapture data to be able to pull the flow, number released, and
  # number recaptured for every unique mark-recapture experiment for data inputs
  mark_recapture_data_unique_experiments <- mark_recapture_data |>
    distinct(site, run_year, week, number_released, number_recaptured,
             standardized_efficiency_flow)

  # TODO josh has something for BUGS-specific code here
  # if(nrow(weeks_without_recap) == 1 | nrow(weeks_with_recap) == 1)

  # set up b-spline basis matrix
  # this corresponds to line 148-153 in josh_original_model_code.R
  k_int <- 4 # rule of thumb is 1 knot for every 4 data points for a cubic spline (which has 4 parameters)
  # TODO create "spline parameters" object
  number_knots <- round(nrow(input_data) / k_int, 0)
  first_knot_position <- 2
  final_knot_position <- nrow(input_data) - 1 # keep first and/or last knot positions away from tails if there are intervals with no sampling on the tails
  knot_positions <- seq(first_knot_position, final_knot_position, length.out = number_knots) # define position of b-spline knots using even interval if no missing data
  b_spline_matrix <- splines2::bSpline(x = 1:nrow(input_data), knots = knot_positions, deg = 3, intercept = T) # bspline basis matrix. One row for each data point (1:number_weeks_catch), and one column for each term in the cubic polynomial function (4) + number of knots
  K <- ncol(b_spline_matrix)

  if(effort_adjust) {
    weekly_catch_data <- input_data$catch_standardized_by_effort
  } else {
    weekly_catch_data <- input_data$catch_standardized_by_effort
  }

  # full data list
  full_data_list <- list("number_efficiency_experiments" = number_efficiency_experiments,
                         "number_sites_pCap" = number_sites_pCap,
                         "indices_site_mark_recapture" = indices_site_mark_recapture,
                         "number_released" = mark_recapture_data_unique_experiments$number_released,
                         "number_recaptured" = mark_recapture_data_unique_experiments$number_recaptured,
                         "number_weeks_catch" = number_weeks_catch,
                         "weekly_catch" = weekly_catch_data,
                         "K" = K,
                         "b_spline_matrix" = b_spline_matrix,
                         "indices_pCap" = indices_pCap,
                         "number_weeks_with_mark_recapture" = number_weeks_with_mark_recapture,
                         "number_weeks_without_mark_recapture" = number_weeks_without_mark_recapture,
                         "indices_with_mark_recapture" = indices_with_mark_recapture,
                         "indices_without_mark_recapture" = indices_without_mark_recapture,
                         "indices_sites_pCap" = indices_sites_pCap,
                         "number_weeks_with_catch" = number_weeks_with_catch,
                         "indices_with_catch" = indices_with_catch,
                         "standardized_efficiency_flow" = mark_recapture_data_unique_experiments$standardized_efficiency_flow,
                         "catch_flow" = input_data$flow_cfs,
                         "lgN_max" = input_data$lgN_prior) # TODO rename?

  # set up models and data to pass to bugs
  number_experiments_at_site <- mark_recapture_data |>
    distinct(site, run_year, week) |>
    group_by(site) |>
    tally() |>
    dplyr::filter(site == !!site) |>
    pull(n)

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

  parameters <- c("pCap_mu_prior", "pCap_sd_prior", "flow_mu_prior", "flow_sd_prior", "process_error_sd_prior", "b0_pCap", "b_flow",
                  "pCap_U", "N", "N_total", "sd_N", "sd_Ne")

  # initial parameter values
  ini_b0_pCap <- mark_recapture_data |>
    dplyr::filter(site == !!site,
           run_year == !!run_year) |>
    mutate(ini_b0_pCap = stats::qlogis(sum(number_recaptured) / sum(number_released))) |>
    pull(ini_b0_pCap)

  ini_lgN <- input_data |>
    mutate(ini_lgN = log(catch_standardized_by_effort / 1000 + 2),
           ini_lgN = ifelse(is.na(ini_lgN), log(2 / 1000), ini_lgN)) |>  # TODO double check these
    pull(ini_lgN)

  pCap_mu_prior <- mark_recapture_data |>
    mutate(pCap_mu_prior = gtools::logit(sum(number_recaptured) / sum(number_released))) |>
    pull(pCap_mu_prior)

  init_list <- list("pCap_mu_prior" = pCap_mu_prior,
                    "b0_pCap" = ini_b0_pCap,
                    "flow_mu_prior" = 0,
                    "b_flow" = rep(0, number_sites_pCap),
                    "pCap_tau_prior" = 1,
                    "flow_tau_prior" = 1,
                    "process_error_tau_prior" = 1,
                    "b_sp" = rep(1, K),
                    "lg_N" = ini_lgN)

  inits <- list(inits1 = init_list, inits2 = init_list, inits3 = init_list)

  # run the bugs model
  results <- bt_spas_x_bugs(data, inits, parameters, model_name, bt_spas_x_bayes_params, bugs_directory = paste0(bugs_directory))
  return(results)
}


#' Execute BT-SPAS-X in WinBUGs
#' @details This function is called within `run_single_bt_spas_x()` and calls the WinBUGS code
#' @param data a list containing the following elements:
#' * **number_efficiency_experiments** number of unique mark-recapture experiments performed across tributaries, years and weeks
#' * **number_sites_pCap** number of tributaries to use for the pCap component of the model
#' * **number_weeks_catch** number of weeks with catch data
#' * **number_weeks_with_mark_recapture** number of weeks in catch data with associated pCap and flow data
#' * **number_weeks_without_mark_recapture** number of weeks in catch data without associated pCap and flow data
#' * **number_weeks_with_catch** number of weeks in catch data with catch data (RST fished)
#' * **indices_site_pCap** tributary index (of all possible tributaries) to use to predict efficiency for missing strata. Note length=0 if selected tributary is not part of trib set that has mark-recap data. In this case model that samples from trib hyper will be called.
#' * **indices_site_mark_recapture** indices (1:number_sites_pCap) assigned to the mark-recapture experiment table for use in BUGs
#' * **indices_pCap** indices of weeks in mark-recapture table for U strata being estimated
#' * **indices_without_mark_recapture** indices of weeks in catch data without associated pCap and flow data
#' * **indices_with_mark_recapture** indices of weeks in catch data with associated pCap and flow data
#' * **indices_with_catch** indices of weeks in catch data with catch data
#' * **number_released** number of fish released for each mark-recapture experiment
#' * **number_recaptured** number of fish recaptured for each mark-recapture experiment
#' * **weekly_catch** weekly abundance
#' * **standardized_efficiency_flow** standardized flow, averaged over recapture days (< 1 week)
#' * **catch_flow** standardized flow, averaged by week
#' * **K** number of columns in the bspline basis matrix
#' * **b_spline_matrix** The b_spline_matrix (bspline basis matrix. One row for each data point (1:number_weeks_catch), and one column for each term in the cubic polynomial function (4) + number of knots)
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

  # get operating system - bugs can't run on a mac without serious set-up
  operating_system <- ifelse(grepl("Mac", Sys.info()['nodename']) | grepl("MBP", Sys.info()['nodename']), "mac", "pc")
  if(operating_system == "mac") {
    # TODO update this message to be a [Y/n] and/or link to wine emulator
    cli::cli_alert_warning("This model is currently coded in WinBUGS, which cannot easily be run on a Mac.
                           All the information required to run the model will be returned, but the model will not be run.")

    return(list(data = data, inits = inits, parameters = parameters, model_name = model_name,
                n.chains = bt_spas_x_bayes_params$number_chains,
                n.burnin = bt_spas_x_bayes_params$number_burnin,
                n.thin = bt_spas_x_bayes_params$number_thin,
                n.iter = bt_spas_x_bayes_params$number_mcmc, debug = FALSE, codaPkg = FALSE, DIC = TRUE,
                clearWD = TRUE, bugs_directory = bugs_directory))
  } else {

    cli::cli_process_start("WinBUGS model running")
    # run bugs model
    model_results <- bugs(data, inits, parameters, model_name,
                          n.chains = bt_spas_x_bayes_params$number_chains,
                          n.burnin = bt_spas_x_bayes_params$number_burnin,
                          n.thin = bt_spas_x_bayes_params$number_thin,
                          n.iter = bt_spas_x_bayes_params$number_mcmc,
                          debug = FALSE, codaPkg = FALSE, DIC = TRUE, clearWD = TRUE,
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
    data_needed <- c("number_efficiency_experiments", "number_sites_pCap" ,"indices_site_mark_recapture", "number_released", "number_recaptured", "number_weeks_catch",
                     "weekly_catch", "K", "b_spline_matrix", "indices_pCap", "number_weeks_with_catch", "indices_with_catch", "standardized_efficiency_flow", "lgN_max")
  } else if(model_name == "missing_mark_recap.bug") {
    data_needed <- c("number_efficiency_experiments", "number_sites_pCap", "indices_site_mark_recapture", "number_released", "number_recaptured", "number_weeks_catch",
                     "weekly_catch", "K", "b_spline_matrix", "indices_pCap", "number_weeks_with_mark_recapture", "number_weeks_without_mark_recapture",
                     "indices_with_mark_recapture", "indices_without_mark_recapture", "indices_sites_pCap", "number_weeks_with_catch", "indices_with_catch",
                     "standardized_efficiency_flow", "catch_flow", "lgN_max")
  } else if(model_name == "no_mark_recap.bug") {
    data_needed <- c("number_efficiency_experiments", "number_sites_pCap", "indices_site_mark_recapture", "number_released", "number_recaptured", "number_weeks_catch",
                     "weekly_catch", "K", "b_spline_matrix", "number_weeks_without_mark_recapture", "indices_without_mark_recapture", "indices_sites_pCap",
                     "number_weeks_with_catch", "indices_with_catch", "standardized_efficiency_flow", "catch_flow", "lgN_max")
  } else if(model_name == "no_mark_recap_no_trib.bug") {
    data_needed <- c("number_efficiency_experiments", "number_sites_pCap", "indices_site_mark_recapture", "number_released", "number_recaptured", "number_weeks_catch",
                     "weekly_catch", "K", "b_spline_matrix", "number_weeks_without_mark_recapture", "indices_without_mark_recapture", "number_weeks_with_catch",
                     "indices_with_catch", "standardized_efficiency_flow", "catch_flow", "lgN_max")
  }
  new_data_list <- full_data_list[sapply(names(full_data_list), function(x) x %in% data_needed)]
  return(new_data_list)
}
