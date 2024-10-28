#' Call BT-SPAS-X model on all site/run years
#' @details This calls `run_single_bt_spas_x()` on all site and run year combinations
#' present in the input data.
#' @param bt_spas_x_bayes_params: a list containing `number_mcmc`, `number_burnin`, `number_thin`,
#' and `number_chains`. Can use `SRJPEmodel::bt_spas_x_bayes_params`.
#' @param weekly_juvenile_abundance_catch_data data frame containing weekly RST catch data. See
#' `?SRJPEdata::weekly_juvenile_abundance_catch_data`
#' @param weekly_juvenile_abundance_efficiency_data data frame containing weekly RST catch data. See
#' `?SRJPEdata::weekly_juvenile_abundance_efficiency_data`
#' @param sites_to_run a subset of site/run year/lifestage combinations in tibble format to run
#' @param effort_adjust whether or not you want to use catch adjusted by effort
#' @param bugs_directory where the `WinBUGS.exe` file can be found. Needs to end in `/WinBUGS`
#' @param debug_mode whether you want to run `bugs` in debug mode.
#' @returns a list containing the results for every unique `site` and `run year` combination
#' in the input data.
#' @export
#' @md
run_multiple_bt_spas_x <- function(bt_spas_x_bayes_params,
                                   weekly_juvenile_abundance_catch_data,
                                   weekly_juvenile_abundance_efficiency_data,
                                   sites_to_run = NULL,
                                   effort_adjust,
                                   bugs_directory, debug_mode,
                                   no_cut = c(F, T)) {

  if(!is.null(sites_to_run)) {
    site_run_year_combinations <- sites_to_run
  } else {
    site_run_year_combinations <- weekly_juvenile_abundance_catch_data |>
      distinct(site, run_year, life_stage)
  }

  all_results <- list()
  # TODO convert to purrr::map2
  for(i in 1:nrow(site_run_year_combinations)) {

    cli::cli_bullets(paste0("running bt-spas-x on ", site_run_year_combinations$site[i],
                            " run year ", site_run_year_combinations$run_year[i],
                            " and lifestage ", site_run_year_combinations$life_stage[i]))

      all_results[[i]] <- tryCatch({run_single_bt_spas_x(bt_spas_x_bayes_params,
                                                         weekly_juvenile_abundance_catch_data,
                                                         weekly_juvenile_abundance_efficiency_data,
                                                          site_run_year_combinations$site[i],
                                                          site_run_year_combinations$run_year[i],
                                                          site_run_year_combinations$life_stage[i],
                                                          effort_adjust,
                                                          bugs_directory, debug_mode,
                                                          no_cut)
        },
        error = function(e) return(1e12)
      )
      if(!is.list(all_results[[i]])) {
        all_results[[i]] <- tibble("site" = site_run_year_combinations$site[i],
                                   "run_year" = site_run_year_combinations$run_year[i],
                                   "life_stage" = site_run_year_combinations$life_stage[i],
                                   "statistic" = "error")
      } else {
        all_results[[i]] <- all_results[[i]]$final_results
      }
  }

  all_results_df <- bind_rows(all_results)

  errors <- all_results_df %>%
    filter(statistic == "error") %>%
    nrow()

  non_errors <- all_results_df %>%
    filter(statistic != "error") %>%
    distinct(life_stage, run_year, site) %>%
    nrow()

  cli::cli_bullets(paste0("Of ", nrow(site_run_year_combinations), " site, run year, and life stage
                         combinations fit, ", (errors/(errors + non_errors) * 100), "% errored out."))

  return(all_results_df)
}


#' Call BT-SPAS-X on a single site/run year combination
#' @details This function is called within `run_bt_spas_x()` or can be run by itself
#' on a single site/run year combination.
#' @param bt_spas_x_bayes_params: a list containing `number_mcmc`, `number_burnin`, `number_thin`,
#' and `number_chains`. Can use `SRJPEmodel::bt_spas_x_bayes_params`.
#' @param weekly_juvenile_abundance_catch_data data frame containing weekly RST catch data. See
#' `?SRJPEdata::weekly_juvenile_abundance_catch_data`
#' @param weekly_juvenile_abundance_efficiency_data data frame containing weekly RST catch data. See
#' `?SRJPEdata::weekly_juvenile_abundance_efficiency_data`
#' @param site site for which you want to fit the model
#' @param run_year run year for which you want to fit the model
#' @param lifetage the lifestage for which you want tor run the model. One of `yearling`, `fry`, and `smolt`.
#' @param effort_adjust whether or not you want to use catch adjusted by effort
#' @param bugs_directory where the `WinBUGS.exe` file can be found. Needs to end in `/WinBUGS`
#' @param debug_mode whether you want to run `bugs` in debug mode.
#' @returns a list:
#' * **results** model results - see `?bt_spas_x_bugs()` for details.
#' * **site** the site used to fit the model
#' * **run_year** the run year used to fit the model
#' * **weeks_fit** the weeks with catch used to fit the model (for analysis and plotting)
#' * **knots_output** knot positions for the weekly abundance spline curve
#' @export
#' @md
run_single_bt_spas_x <- function(bt_spas_x_bayes_params,
                                 weekly_juvenile_abundance_catch_data,
                                 weekly_juvenile_abundance_efficiency_data,
                                 site, run_year, lifestage,
                                 effort_adjust = c(T, F),
                                 bugs_directory, debug_mode,
                                 no_cut = c(F, T)) {

  # prepare "catch" dataset - filtered to weeks, site, run_year, and lifestage selected
  # catch_flow is average for julian week, standardized_efficiency_flow is average over recapture days (< 1 week)
  catch_data <- weekly_juvenile_abundance_catch_data |>
    mutate(filter_out = ifelse(is.na(life_stage) & count > 0, TRUE, FALSE)) |> # we do not want to keep NA lifestage associated with counts > 0
    filter(!filter_out,
           run_year == !!run_year,
           site == !!site,
           week %in% c(seq(45, 53), seq(1, 22)),
           life_stage %in% c(lifestage, NA)) |>
    mutate(count = round(count, 0),
           catch_standardized_by_hours_fished = round(catch_standardized_by_hours_fished, 0))

  if(nrow(catch_data) == 0) {
    cli::cli_alert_warning(paste0("There is no catch data for site ", site,
                                  ", run year ", run_year, ", and lifestage ", lifestage, ". Please try with a different combination of site and year."))
    return(-99)
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
    remove_sites <- weekly_juvenile_abundance_efficiency_data |>
      filter(!site %in% c("knights landing", "tisdale", "red bluff diversion dam")) |>
      distinct(site) |>
      pull(site)
  }

  # prepare "mark recapture" dataset - all mark-recap trials in the system
  mark_recapture_data <- weekly_juvenile_abundance_efficiency_data |>
    # grab standardized_flow
    left_join(weekly_juvenile_abundance_catch_data |>
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

  # get indexing for "mark recap" dataset (pCap model)
  Ntribs <- length(unique(mark_recapture_data$site)) # number of sites (for pCap calculations)
  number_efficiency_experiments <- unique(mark_recapture_data[c("site", "run_year", "week")]) |>
    nrow() # number of efficiency experiments completed, nrow(mark_recapture_data) this depends on whether you have lifestage or not
  years_with_efficiency_experiments <- unique(mark_recapture_data$run_year) # years where efficiency experiments were done
  indices_sites_pCap <- which(unique(mark_recapture_data$site) == site) # indices of those sites where efficiency trials were performed, can be length = 0
  indices_with_mark_recapture <- which(!is.na(all_data_for_indexing$number_released) &
                                         !is.na(all_data_for_indexing$standardized_flow)) # indices of efficiency experiments in catch data
  weeks_with_mark_recapture <- all_data_for_indexing$week[indices_with_mark_recapture] # weeks (in catch data) where mark recapture were performed
  indices_pCap <- which(mark_recapture_data$site == site &
                          mark_recapture_data$run_year == run_year &
                          mark_recapture_data$week %in% weeks_with_mark_recapture)   # indices (in mark-recap data) for the selected site and run year, filtered to weeks where mark-recap were performed (in catch data)
  indices_without_mark_recapture <- which(is.na(all_data_for_indexing$number_released) |
                                            is.na(all_data_for_indexing$standardized_flow))   # indices (in catch data) where no mark recap were performed
  indices_site_mark_recapture <- mark_recapture_data |>
    group_by(site) |>
    mutate(ID = cur_group_id()) |>
    pull(ID)   # indices (in mark-recap data) for each site
  number_weeks_with_mark_recapture <- length(indices_with_mark_recapture) # number of weeks (in mark-recap data) where effiency experiments were performed
  number_weeks_without_mark_recapture <- length(indices_without_mark_recapture)   # number of weeks (in mark-recap data) where effiency experiments were not performed

  # TODO keep this? so BUGS doesn't bomb
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
    return(-99)
  }

  # spline parameter calculation
  spline_data <- SRJPEmodel::build_spline_data(number_weeks_catch, k_int = 4) # rule of thumb is 1 knot for every 4 data points for a cubic spline (which has 4 parameters)

  if(effort_adjust) {
    weekly_catch_data <- catch_data$catch_standardized_by_hours_fished
  } else {
    weekly_catch_data <- catch_data$count
  }

  # build data list with ALL elements
  full_data_list <- list("Nmr" = number_efficiency_experiments,
                         "Ntribs" = Ntribs,
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
                         "mr_flow" = mark_recapture_data$standardized_flow,
                         "catch_flow" = catch_data$standardized_flow,
                         "lgN_max" = catch_data$lgN_prior)

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

  data <- get_bt_spas_x_data_list(model_name, full_data_list)

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

  parameters <- c("trib_mu.P", "trib_sd.P", "flow_mu.P", "flow_sd.P", "pro_sd.P",
                   "b0_pCap", "b_flow", "pCap_U", "N", "Ntot", "sd.N", "sd.Ne")

  # initial parameter values
  ini_b0_pCap <- rep(NA, Ntribs)
  for(i in 1:Ntribs) {
    irows = which(indices_site_mark_recapture == i)
    ini_b0_pCap[i] = gtools::logit(sum(mark_recapture_data$number_recaptured[irows]) /
                                     sum(mark_recapture_data$number_released[irows]))
    if(is.nan(ini_b0_pCap[i]) | is.infinite(ini_b0_pCap[i])) {
      # -Inf happens when number recaptured == 0, logit of 0 is -Inf
      ini_b0_pCap[i] <- -5
    }
  }

  ini_lgN <- catch_data |>
    mutate(ini_lgN = log(catch_standardized_by_hours_fished / 1000 + 2),
           ini_lgN = ifelse(is.na(ini_lgN), log(2 / 1000), ini_lgN)) |>
    pull(ini_lgN)

  pCap_mu_prior <- gtools::logit(sum(mark_recapture_data$number_recaptured) /
                                   sum(mark_recapture_data$number_released))

  init_list <- list(trib_mu.P = pCap_mu_prior,
                    b0_pCap = ini_b0_pCap,
                    flow_mu.P = 0,
                    b_flow = rep(0, Ntribs),
                    trib_tau.P = 1,
                    flow_tau.P = 1,
                    pro_tau.P = 1,
                    b_sp = rep(1, spline_data$K),
                    lg_N = ini_lgN)


  cli::cli_process_start("Checking init inputs",
                         msg_done = "Inits checked",
                         msg_failed = "Init check failed")
  invisible(lapply(names(init_list), function(x) {
    if(any(is.nan(init_list[[x]])) | any(is.infinite(init_list[[x]]))) {
      cli::cli_abort(paste0("NaNs detected in ", x, ". Please check your input data."))
    }
  }))
  cli::cli_process_done()

  inits <- list(init_list, init_list, init_list)

  # run the bugs model
  results <- bt_spas_x_bugs(data, inits, parameters, model_name, bt_spas_x_bayes_params,
                            bugs_directory = paste0(bugs_directory), debug_mode = debug_mode,
                            no_cut = no_cut)

  # get operating system - bugs can't run on a mac without serious set-up
  operating_system <- ifelse(grepl("Mac", Sys.info()['nodename']) | grepl("MBP", Sys.info()['nodename']), "mac", "pc")
  if(operating_system == "mac") {
    cli::cli_bullets("Run on a mac, not able to run full model")
    return(results)
  } else {
    final_results <- get_summary_table(results, site, run_year, lifestage,
                                       weeks_fit = catch_data$week[data$Uwc_ind],
                                       sites_fit = unique(mark_recapture_data$site),
                                       model_name = model_name)

    # TODO build workflow to upload this to the cloud
    # TODO split these out and document
    diagnostic_results <- list("input_data" = data, # TODO may not need if we store SRJPEdata version
                               "input_initial_values" = inits[[1]], # TODO may not need if storing SRJPEdata version
                               "bayes_params" = bt_spas_x_bayes_params,
                               "knots_output" = spline_data$knot_positions,
                               "full_model_object" = results)

    # return(final_results)

    # TODO delete this, just testing
    return(list("final_results" = final_results,
                "full_object" = results))

    # TODO check for convergence (rhat > 1.05)

  }
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
#' @param inits a list containing initial values for the following parameters:
#' * **trib_mu.P** mean of hyper-distribution for site-effect
#' * **b0_pCap** site effect on trap efficiency for each site
#' * **flow_mu.P** mean of hyper-distribution for flow effect
#' * **b_flow** flow effect on trap efficiency for each site
#' * **trib_tau.P** used to estimate standard deviation of hyper-distribution for site effect
#' * **flow_tau.P** used to estimate standard deviation of hyper-distribution for flow effect
#' * **pro_tau.P** used to estimate standard deviation of zero-centered normal distribution for unexplained error
#' * **b_sp** basis function of each spline node
#' * **lg_N** predicted weekly abundance
#' @param parameters a list of parameters to be estimated in the model:
#' * **trib_mu.P** mean of hyper-distribution for site-effect
#' * **trib_sd.P** standard deviation of hyper-distribution for site effect
#' * **flow_mu.P** mean of hyper-distribution for flow effect
#' * **flow_sd.P** standard deviation of hyper-distribution for flow effect
#' * **pro_sd.P** standard deviation of zero-centered normal distribution for unexplained error
#' * **b0_pCap** site effect on trap efficiency for each site
#' * **b_flow** flow effect on trap efficiency for each site
#' * **pCap_U** weekly trap efficiency (capture probability)
#' * **N** weekly juvenile abundance
#' * **Ntot** total juvenile abundance for the year
#' * **sd.N** standard deviation controlling flexibility of spline weekly abundance curve
#' * **sd.Ne** standard deviation controlling extent of non-spline variation in weekly abundance
#' @param model_name model to be called based on number of efficiency trials available for the selected site. Either
#' * **all_mark_recap.bug** all weeks with catch have corresponding efficiency trials
#' * **missing_mark_recap.bug** some weeks with catch have corresponding efficiency trials
#' * **no_mark_recap_no_trib.bug** the selected site has no efficiency data at all
#' * **no_mark_recap.bug** no weeks with catch have corresponding efficiency trials
#' @param bt_spas_x_bayes_params: a list containing `number_mcmc`, `number_burnin`, `number_thin`,
#' and `number_chains`.
#' @param bugs_directory a filepath indicating where to find the `WinBUGS14/` file. This needs
#' to be in a character format ending with `/WinBUGS14`.
#' @returns if running on a Mac operating system, returns a list of all inputs formatted
#' to pass to `bugs`. If running on a PC or another operating system capable of running WinBUGS,
#' returns a nested list containing the following elements:
#' * **model_results** the `BUGs` object from fitting the model
#' * **model_called** the model called on the data
#' * **data_inputs** the data passed to the model
#' * **init_inputs** the initial values passed to the model
#' @export
#' @md
bt_spas_x_bugs <- function(data, inits, parameters, model_name, bt_spas_x_bayes_params,
                           number_mcmc, bugs_directory, debug_mode,
                           no_cut) {

  # TODO remove - this is for testing no_cut models
  if(no_cut) {
    # set model name directory
    model_name_full <- eval(parse(text = paste0("SRJPEmodel::bt_spas_x_model_code$winbugs$no_cut$", model_name)))
  } else {
    model_name_full <- eval(parse(text = paste0("SRJPEmodel::bt_spas_x_model_code$winbugs$cut$", model_name)))
  }

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
    # TODO wrap this in a tryCatch
    # run bugs model
    model_results <- R2WinBUGS::bugs(data, inits, parameters, model.file = model_name_full,
                                     n.chains = bt_spas_x_bayes_params$number_chains,
                                     n.burnin = bt_spas_x_bayes_params$number_burnin,
                                     n.thin = bt_spas_x_bayes_params$number_thin,
                                     n.iter = bt_spas_x_bayes_params$number_mcmc,
                                     debug = debug_mode, codaPkg = FALSE, DIC = TRUE, clearWD = TRUE,
                                     bugs.directory = bugs_directory)

    return(model_results)

    # return(list("model_results" = model_results,
    #             "model_called" = sub(".*/", "", model_name_full),
    #             "data_inputs" = data,
    #             "init_inputs" = inits))
  }
}


#' Prepare Data Object for BT-SPAS-X WinBUGs call
#' @details This function is called within `run_single_bt_spas_x()` and prepares the data lists
#' based on which model will be run.
#' @param model_name which model to call on the data. Either `all_mark_recap.bug`,
#' `missing_mark_recap.bug`, `no_mark_recap.bug`, or `no_mark_recap_no_trib.bug`
#' @param full_data_list a list containing all possible data objects to use in the model.
#' @returns a named list with the required elements for that model run.
#' @keywords internal
#' @export
#' @md
get_bt_spas_x_data_list <- function(model_name, full_data_list) {
  if(str_detect(model_name, "all_mark_recap")) {
    data_needed <- c("Nmr", "Ntribs" ,"ind_trib", "Releases", "Recaptures", "Nstrata",
                     "u", "K", "ZP", "ind_pCap", "Nstrata_wc", "Uwc_ind", "mr_flow", "lgN_max")
  } else if(str_detect(model_name, "missing_mark_recap")) {
    data_needed <- c("Nmr", "Ntribs", "ind_trib", "Releases", "Recaptures", "Nstrata",
                     "u", "K", "ZP", "ind_pCap", "Nwmr", "Nwomr",
                     "Uind_wMR", "Uind_woMR", "use_trib", "Nstrata_wc", "Uwc_ind",
                     "mr_flow", "catch_flow", "lgN_max")
  } else if(str_detect(model_name, "no_mark_recap")) {
    data_needed <- c("Nmr", "Ntribs", "ind_trib", "Releases", "Recaptures", "Nstrata",
                     "u", "K", "ZP", "Nwomr", "Uind_woMR", "use_trib",
                     "Nstrata_wc", "Uwc_ind", "mr_flow", "catch_flow", "lgN_max")
  } else if(str_detect(model_name, "no_mark_recap_no_trib")) {
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
#' @param number_weeks_catch number of weeks with catch data.
#' @param k_int number of knots to use in building spline for abundance data.
#' Rule of thumb is 1 knot for every 4 data points for a cubic spline (which has 4 parameters)
#' @returns a named list containing *K* and *b_spline_matrix* for passing to WinBUGS
#' @keywords internal
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
              "b_spline_matrix" = b_spline_matrix,
              "knot_positions" = knot_positions))
}


#' @title Extract summary in table
#' @description TODO
#' @keywords internal
#' @export
#' @md
get_summary_table <- function(model_fit_object, site, run_year,
                              lifestage, weeks_fit, sites_fit,
                              model_name) {

  summary_table <- model_fit_object$summary |>
    as.data.frame() |>
    cbind(par_names = rownames(model_fit_object$summary)) |>
    janitor::clean_names() |>
    # don't extract index for parameters estimated by stream, not week
    mutate(week_index = ifelse(str_detect(par_names, "b0_pCap|b_flow"), NA,
                               suppressWarnings(readr::parse_number(par_names))),
           site_index = ifelse(str_detect(par_names, "b0_pCap|b_flow"),
                               suppressWarnings(readr::parse_number(substr(par_names, 3, length(par_names)))),
                               NA),
           site = site,
           run_year = run_year,
           life_stage = lifestage,
           model_called = model_name)

  rownames(summary_table) = NULL

  weeks_fit_lookup <- tibble(week_fit = weeks_fit) |>
    mutate(week_index = row_number())

  sites_fit_lookup <- tibble(site_fit_hierarchical = sites_fit) |>
    mutate(site_index = row_number())

  summary_table_final <- summary_table |>
    left_join(weeks_fit_lookup, by = "week_index") |>
    left_join(sites_fit_lookup, by = "site_index") |>
    select(site, run_year, life_stage, week_fit, site_fit_hierarchical,
           parameter = par_names,
           mean, sd, `2.5` = x2_5_percent,
           `25` = x25_percent, `50` = x50_percent,
           `75` = x75_percent, `97.5` = x97_5_percent,
           rhat, n_eff, model_name = model_called) |>
    mutate(srjpedata_version = as.character(packageVersion("SRJPEdata")),
           converged = ifelse(rhat <= 1.05, TRUE, FALSE)) |>
    pivot_longer(mean:n_eff, names_to = "statistic", values_to = "value")

  summary_table_final
}

#' Extract results from BT-SPAS-X bugs object
#' @details This function is called within `run_single_bt_spas_x()` and prepares spline parameters
#' to pass to the data list. These are helpful for model diagnostics for a given model run.
#' @param model_fit_object the object produced by running `run_single_bt_spas_x()`
#' for a given site and run year.
#' @returns a named list with the following elements:
#' * **posterior_output** the posterior simulations produced by the `BUGs` object
#' * **summary_output** the summary produced by the `BUGs` object
#' * **dic_output** DIC for the `BUGs` object
#' * **knots_output** the knot positions produced for the weekly abundance spline curve
#' @export
#' @md
extract_bt_spas_x_results <- function(model_fit_object) {
  posterior_output <- model_fit_object$results$model_results$sims.list
  summary_output <- round(model_fit_object$results$model_results$summary, 3)
  dic_output <- c(model_fit_object$results$model_results$pD,
                  model_fit_object$results$model_results$DIC)
  knots_output <- model_fit_object$knots_output

  return(list("posterior_output" = posterior_output,
              "summary_output" = summary_output,
              "dic_output" = dic_output,
              "knots_output" = knots_output))
}


# STAN version ------------------------------------------------------------

#' Call BT-SPAS-X using STAN
#' @details This is a draft version of the function that calls STAN instead of WinBUGS.
#' @param bt_spas_x_bayes_params: a list containing `number_mcmc`, `number_burnin`, `number_thin`,
#' and `number_chains`. Can use `SRJPEmodel::bt_spas_x_bayes_params`.
#' @param weekly_juvenile_abundance_catch_data data frame containing weekly RST catch data. See
#' `?SRJPEdata::weekly_juvenile_abundance_catch_data`
#' @param weekly_juvenile_abundance_efficiency_data data frame containing weekly RST catch data. See
#' `?SRJPEdata::weekly_juvenile_abundance_efficiency_data`
#' @param site site for which you want to fit the model
#' @param run_year run year for which you want to fit the model
#' @param lifetage the lifestage for which you want tor run the model. One of `yearling`, `fry`, and `smolt`.
#' @param effort_adjust whether or not you want to use catch adjusted by effort
#' @export
#' @md
run_single_bt_spas_x_stan <- function(bt_spas_x_bayes_params,
                                       weekly_juvenile_abundance_catch_data,
                                       weekly_juvenile_abundance_efficiency_data,
                                       site, run_year, lifestage,
                                       effort_adjust = c(T, F)) {

  # some streams do not have lifestage data (i.e. Mokelumne, American) so we summarize all lifestages together
  if(is.na(lifestage)) {
    cli::cli_bullets("No lifestage passed, so grouping all lifestages for analysis")
    catch_data <- weekly_juvenile_abundance_catch_data |>
      select(-life_stage) |>
      filter(run_year == !!run_year,
             site == !!site,
             week %in% c(seq(45, 53), seq(1, 22))) |>
      group_by(year, week, stream, site, run_year) |>
      summarise(count = sum(count, na.rm = T),
                mean_fork_length = mean(mean_fork_length, na.rm = T),
                hours_fished = mean(hours_fished, na.rm = T),
                flow_cfs = mean(flow_cfs, na.rm = T),
                average_stream_hours_fished = mean(average_stream_hours_fished, na.rm = T),
                standardized_flow = mean(standardized_flow, na.rm = T),
                catch_standardized_by_hours_fished = sum(catch_standardized_by_hours_fished, na.rm = T),
                lgN_prior = mean(lgN_prior, na.rm = T)) |>
      ungroup() |>
      mutate(count = round(count, 0),
             catch_standardized_by_hours_fished = round(catch_standardized_by_hours_fished, 0))
  } else {
    catch_data <- weekly_juvenile_abundance_catch_data |>
      mutate(filter_out = ifelse(is.na(life_stage) & count > 0, TRUE, FALSE)) |> # we do not want to keep NA lifestage associated with counts > 0
      filter(!filter_out,
             run_year == !!run_year,
             site == !!site,
             week %in% c(seq(45, 53), seq(1, 22)),
             life_stage %in% c(lifestage, NA)) |>
      mutate(count = round(count, 0),
             catch_standardized_by_hours_fished = round(catch_standardized_by_hours_fished, 0))
  }

  if(nrow(catch_data) == 0) {
    cli::cli_alert_warning(paste0("There is no catch data for site ", site,
                                  ", run year ", run_year, ", and lifestage ", lifestage, ". Please try with a different combination of site and year."))
    return(-99)
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
    remove_sites <- weekly_juvenile_abundance_efficiency_data |>
      filter(!site %in% c("knights landing", "tisdale", "red bluff diversion dam")) |>
      distinct(site) |>
      pull(site)
  }

  # prepare "mark recapture" dataset - all mark-recap trials in the system
  mark_recapture_data <- weekly_juvenile_abundance_efficiency_data |>
    # grab standardized_flow
    left_join(weekly_juvenile_abundance_catch_data |>
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

  # get indexing for "mark recap" dataset (pCap model)
  Ntribs <- length(unique(mark_recapture_data$site)) # number of sites (for pCap calculations)
  number_efficiency_experiments <- unique(mark_recapture_data[c("site", "run_year", "week")]) |>
    nrow() # number of efficiency experiments completed, nrow(mark_recapture_data) this depends on whether you have lifestage or not
  years_with_efficiency_experiments <- unique(mark_recapture_data$run_year) # years where efficiency experiments were done
  indices_sites_pCap <- which(unique(mark_recapture_data$site) == site) # indices of those sites where efficiency trials were performed, can be length = 0
  indices_with_mark_recapture <- which(!is.na(all_data_for_indexing$number_released) &
                                         !is.na(all_data_for_indexing$standardized_flow)) # indices of efficiency experiments in catch data
  weeks_with_mark_recapture <- all_data_for_indexing$week[indices_with_mark_recapture] # weeks (in catch data) where mark recapture were performed
  indices_pCap <- which(mark_recapture_data$site == site &
                          mark_recapture_data$run_year == run_year &
                          mark_recapture_data$week %in% weeks_with_mark_recapture)   # indices (in mark-recap data) for the selected site and run year, filtered to weeks where mark-recap were performed (in catch data)
  indices_without_mark_recapture <- which(is.na(all_data_for_indexing$number_released) |
                                            is.na(all_data_for_indexing$standardized_flow))   # indices (in catch data) where no mark recap were performed
  indices_site_mark_recapture <- mark_recapture_data |>
    group_by(site) |>
    mutate(ID = cur_group_id()) |>
    pull(ID)   # indices (in mark-recap data) for each site
  number_weeks_with_mark_recapture <- length(indices_with_mark_recapture) # number of weeks (in mark-recap data) where effiency experiments were performed
  number_weeks_without_mark_recapture <- length(indices_without_mark_recapture)   # number of weeks (in mark-recap data) where effiency experiments were not performed

  # TODO keep this? so BUGS doesn't bomb
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
    return(-99)
  }

  # spline parameter calculation
  spline_data <- SRJPEmodel::build_spline_data(number_weeks_catch, k_int = 4) # rule of thumb is 1 knot for every 4 data points for a cubic spline (which has 4 parameters)

  if(effort_adjust) {
    weekly_catch_data <- catch_data$catch_standardized_by_hours_fished
  } else {
    weekly_catch_data <- catch_data$count
  }

  # build data list with ALL elements
  full_data_list <- list("Nmr" = number_efficiency_experiments,
                         "Ntribs" = Ntribs,
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
                         "mr_flow" = mark_recapture_data$standardized_flow,
                         "catch_flow" = catch_data$standardized_flow,
                         "lgN_max" = catch_data$lgN_prior)

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
      model_name <- "all_mark_recap" # .stan
    } else {
      # some or all strata don't have efficiency trials
      if(number_weeks_with_mark_recapture > 0) {
        # some weeks have efficiency trials
        model_name <- "missing_mark_recap" # .stan
      } else if(number_weeks_with_mark_recapture == 0) {
        # no weeks have efficiency trials
        model_name <- "no_mark_recap" # .stan
      }
    }
  } else if(number_experiments_at_site == 0) { # no efficiency trials were performed at that site
    model_name <- "no_mark_recap_no_trib" # .stan
  }

  data <- get_bt_spas_x_data_list(model_name, full_data_list)

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

  parameters <- c("trib_mu.P", "trib_sd.P", "flow_mu.P", "flow_sd.P", "pro_sd.P",
                  "b0_pCap", "b_flow", "pCap_U", "N", "Ntot", "sd.N", "sd.Ne")

  # initial parameter values
  ini_b0_pCap <- rep(NA, Ntribs)
  for(i in 1:Ntribs) {
    irows = which(indices_site_mark_recapture == i)
    ini_b0_pCap[i] = gtools::logit(sum(mark_recapture_data$number_recaptured[irows]) /
                                     sum(mark_recapture_data$number_released[irows]))
    if(is.nan(ini_b0_pCap[i]) | is.infinite(ini_b0_pCap[i])) {
      # -Inf happens when number recaptured == 0, logit of 0 is -Inf
      ini_b0_pCap[i] <- -5
    }
  }

  ini_lgN <- catch_data |>
    mutate(ini_lgN = log(catch_standardized_by_hours_fished / 1000 + 2),
           ini_lgN = ifelse(ini_lgN %in% c(NA, Inf), log(2 / 1000), ini_lgN)) |>
           #ini_lgN = ifelse(is.na(ini_lgN), log(2 / 1000), ini_lgN)) |>
    pull(ini_lgN)

  pCap_mu_prior <- gtools::logit(sum(mark_recapture_data$number_recaptured) /
                                   sum(mark_recapture_data$number_released))

  init_list <- list(trib_mu.P = pCap_mu_prior,
                    b0_pCap = ini_b0_pCap,
                    flow_mu.P = 0,
                    b_flow = rep(0, Ntribs),
                    trib_tau.P = 1,
                    flow_tau.P = 1,
                    pro_tau.P = 1,
                    b_sp = rep(1, spline_data$K),
                    lg_N = ini_lgN)


  cli::cli_process_start("Checking init inputs",
                         msg_done = "Inits checked",
                         msg_failed = "Init check failed")
  invisible(lapply(names(init_list), function(x) {
    if(any(is.nan(init_list[[x]])) | any(is.infinite(init_list[[x]]))) {
      cli::cli_abort(paste0("NaNs detected in ", x, ". Please check your input data."))
    }
  }))
  cli::cli_process_done()

  inits <- list(init_list, init_list, init_list)

  # run the bugs model
  results <- bt_spas_x_stan(data, inits, parameters, model_name, bt_spas_x_bayes_params)

  return(results)
}

#' Execute BT-SPAS-X in STAN
#' @details This function is called within `run_single_bt_spas_x_stan()` and calls the STAN code
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
#' @param inits a list containing initial values for the following parameters:
#' * **trib_mu.P** mean of hyper-distribution for site-effect
#' * **b0_pCap** site effect on trap efficiency for each site
#' * **flow_mu.P** mean of hyper-distribution for flow effect
#' * **b_flow** flow effect on trap efficiency for each site
#' * **trib_tau.P** used to estimate standard deviation of hyper-distribution for site effect
#' * **flow_tau.P** used to estimate standard deviation of hyper-distribution for flow effect
#' * **pro_tau.P** used to estimate standard deviation of zero-centered normal distribution for unexplained error
#' * **b_sp** basis function of each spline node
#' * **lg_N** predicted weekly abundance
#' @param parameters a list of parameters to be estimated in the model:
#' * **trib_mu.P** mean of hyper-distribution for site-effect
#' * **trib_sd.P** standard deviation of hyper-distribution for site effect
#' * **flow_mu.P** mean of hyper-distribution for flow effect
#' * **flow_sd.P** standard deviation of hyper-distribution for flow effect
#' * **pro_sd.P** standard deviation of zero-centered normal distribution for unexplained error
#' * **b0_pCap** site effect on trap efficiency for each site
#' * **b_flow** flow effect on trap efficiency for each site
#' * **pCap_U** weekly trap efficiency (capture probability)
#' * **N** weekly juvenile abundance
#' * **Ntot** total juvenile abundance for the year
#' * **sd.N** standard deviation controlling flexibility of spline weekly abundance curve
#' * **sd.Ne** standard deviation controlling extent of non-spline variation in weekly abundance
#' @param model_name model to be called based on number of efficiency trials available for the selected site. Either
#' * **all_mark_recap.bug** all weeks with catch have corresponding efficiency trials
#' * **missing_mark_recap.bug** some weeks with catch have corresponding efficiency trials
#' * **no_mark_recap_no_trib.bug** the selected site has no efficiency data at all
#' * **no_mark_recap.bug** no weeks with catch have corresponding efficiency trials
#' @param bt_spas_x_bayes_params: a list containing `number_mcmc`, `number_burnin`, `number_thin`,
#' and `number_chains`.
#' @returns a stanfit object
#' @export
#' @md
bt_spas_x_stan <- function(data, inits, parameters, model_name,
                           bt_spas_x_bayes_params) {

  cli::cli_process_start("STAN model running")
  stan_model <- eval(parse(text = paste0("SRJPEmodel::bt_spas_x_model_code$stan$", model_name)))

  # options(mc.cores=parallel::detectCores())
  # rstan_options(auto_write=TRUE)

  model_results <- rstan::stan(model_code = stan_model,
                               data = data,
                               init = inits,
                               chains = bt_spas_x_bayes_params$number_chains,
                               thin = bt_spas_x_bayes_params$number_thin,
                               iter = bt_spas_x_bayes_params$number_mcmc,
                               warmup = bt_spas_x_bayes_params$number_burnin,
                               seed = 84735)

    return(model_results)
}
