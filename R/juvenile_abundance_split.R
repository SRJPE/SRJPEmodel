#' Call BT-SPAS-X (WinBUGS) on a single site/run year combination
#' @details This function prepares data for input into a pCap STAN model.
#' @param weekly_juvenile_abundance_catch_data data frame containing weekly RST catch data. See
#' `?SRJPEdata::weekly_juvenile_abundance_catch_data`
#' @param weekly_juvenile_abundance_efficiency_data data frame containing weekly RST catch data. See
#' `?SRJPEdata::weekly_juvenile_abundance_efficiency_data`
#' @param site site for which you want to fit the model
#' @param run_year run year for which you want to fit the model
#' @param effort_adjust whether or not you want to use catch adjusted by effort.
#' @returns a list:
#' * **pCap_inputs** a list of data and inits for input into the pCap model.
#' * **abundance_inputs** a list of data and inits for input into the abundance model.
#' @export
#' @md
prepare_inputs_pCap_abundance_STAN <- function(weekly_juvenile_abundance_catch_data,
                                               weekly_juvenile_abundance_efficiency_data,
                                               site, run_year,
                                               effort_adjust = c(T, F),
                                               default_lgN_prior_denominator = NULL) {

  # group and summarize for all lifestages
  # filter to site and run_year
  cli::cli_bullets("Grouping all lifestages for analysis")
  catch_data <- weekly_juvenile_abundance_catch_data |>
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
  years_with_efficiency_experiments <- unique(mark_recapture_data$run_year) # years where efficiency experiments were done
  indices_sites_pCap <- which(sites_fit == site) # indices of those sites where efficiency trials were performed, can be length = 0
  indices_with_mark_recapture <- which(!is.na(all_data_for_indexing$number_released) &
                                         !is.na(all_data_for_indexing$standardized_efficiency_flow)) # indices of efficiency experiments in catch data
  weeks_with_mark_recapture <- all_data_for_indexing$week[indices_with_mark_recapture] # weeks (in catch data) where mark recapture were performed
  indices_pCap <- which(mark_recapture_data$site == site &
                          mark_recapture_data$run_year == run_year &
                          mark_recapture_data$week %in% weeks_with_mark_recapture)   # indices (in mark-recap data) for the selected site and run year, filtered to weeks where mark-recap were performed (in catch data)
  indices_without_mark_recapture <- which(is.na(all_data_for_indexing$number_released) |
                                            is.na(all_data_for_indexing$standardized_efficiency_flow))   # indices (in catch data) where no mark recap were performed
  number_weeks_with_mark_recapture <- length(indices_with_mark_recapture) # number of weeks (in mark-recap data) where effiency experiments were performed
  number_weeks_without_mark_recapture <- length(indices_without_mark_recapture)   # number of weeks (in mark-recap data) where effiency experiments were not performed

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

  # added 12-2-2024
  # set default lgN prior denominator to the 5th percentile of
  # lincoln peterson abundance
  if(is.null(default_lgN_prior_denominator)) {
    default_lgN_denom <- mark_recapture_data |>
      filter(site == !!site,
             week %in% c(seq(45, 53), seq(1, 22))) |>
      mutate(lp = number_recaptured/number_released) |>
      pull(lp) |>
      quantile(0.05) |>
      unname()

    # set to 0.001 if it is 0 or null
    default_lgN_denom = ifelse(is.null(default_lgN_denom) | is.na(default_lgN_denom) |
                                 is.nan(default_lgN_denom) | default_lgN_denom == 0, 0.001, default_lgN_denom)
  } else {
    default_lgN_denom = default_lgN_prior_denominator
  }

  # using calculation to set both lgN_max data and lgN_max priors (inits)
  ini_lgN <- catch_data |>
    mutate(ini_lgN = log((catch_data$catch_standardized_by_hours_fished / 1000 + 2)/default_lgN_denom),
           ini_lgN = ifelse(is.na(ini_lgN) | is.nan(ini_lgN), log((2 / 1000)/default_lgN_denom), ini_lgN)) |>
    pull(ini_lgN)

  # replace values for NA weeks with default, not low
  weeks_with_no_catch <- which(is.na(catch_data$count))
  ini_lgN[weeks_with_no_catch] <- log(0001 * (min(weekly_catch_data) + 1) / 0.025)

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
                         "lgN_max" = ini_lgN)
                         #"lgN_max" = catch_data$lgN_prior + 13)

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

  cli::cli_bullets(paste0("Calling pCap model ", model_name))

  data <- get_pCap_data_list(model_name, full_data_list)

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
  inputs_for_pCap <- list(data = data,
                          inits = inits,
                          parameters = pCap_parameters,
                          model_name = model_name)

  # inputs for abundance model
  abundance_data <- get_abundance_data_list(full_data_list)

  abundance_parameters <- c("tau_N", "tau_Ne", "b_sp", "lg_N", "lt_pCap_U",
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

  abundance_inits <- list(init_list, init_list, init_list)

  inputs_for_abundance <- list(data = abundance_data,
                               inits = abundance_inits,
                               parameters = abundance_parameters,
                               model_name = model_name)

  weeks_fit <- tibble("Jwk" = catch_data$week) |>
    left_join(SRJPEmodel::julian_week_to_date_lookup, by = "Jwk")

  return(list("pCap_inputs" = inputs_for_pCap,
              "abundance_inputs" = inputs_for_abundance,
              "catch_flow_raw" = catch_data$flow_cfs,
              "mr_flow_raw" = mark_recapture_data$flow_cfs,
              "weeks_fit" = weeks_fit$Jwk,
              "week_date" = weeks_fit$date,
              "sites_fit" = sites_fit))


}


#' Prepare Data Object for pCap STAN call
#' @details This function prepares the data list for a pCap STAN model based on what the model name is.
#' @param model_name which model to call on the data. Either `all_mark_recap.bug`,
#' `missing_mark_recap.bug`, `no_mark_recap.bug`, or `no_mark_recap_no_trib.bug`
#' @param full_data_list a list containing all possible data objects to use in the model.
#' @returns a named list with the required elements for that model run.
#' @keywords internal
#' @export
#' @md
get_pCap_data_list <- function(model_name, full_data_list) {
  if(str_detect(model_name, "all_mark_recap")) {
      data_needed <- c("Ntribs", "Nmr", "Nwmr", "Nstrata", "Nstrata_wc",
                       "Releases", "Recaptures", "mr_flow",
                       "ind_trib", "ind_pCap", "Uwc_ind")
  } else if(str_detect(model_name, "missing_mark_recap")) {
    data_needed <- c("Ntribs", "Nmr", "Nwmr", "Nwomr", "Nstrata", "Nstrata_wc",
                     "use_trib", "ind_trib", "ind_pCap", "mr_flow",
                     "Recaptures", "Releases", "catch_flow", "Uind_wMR", "Uind_woMR")
  } else if(str_detect(model_name, "no_mark_recap")) {
    data_needed <- c("Ntribs", "Nmr", "Nwomr", "Nstrata", "Nstrata_wc",
                     "Releases", "Recaptures", "mr_flow", "catch_flow",
                     "ind_trib", "use_trib", "Uind_woMR")
  } else if(str_detect(model_name, "no_mark_recap_no_trib")) {
    data_needed <- c("Ntribs", "Nmr", "Nwomr", "Nstrata", "Nstrata_wc",
                     "ind_trib", "Releases", "Recaptures", "mr_flow", "catch_flow", "Uind_woMR")
  }
  new_data_list <- full_data_list[sapply(names(full_data_list), function(x) x %in% data_needed)]
  return(new_data_list)
}

#' Prepare Data Object for abundance STAN call
#' @details This function prepares the data list for the abundance STAN model based on what the model name is.
#' @param full_data_list a list containing all possible data objects to use in the model.
#' @returns a named list with the required elements for that model run.
#' @keywords internal
#' @export
#' @md
get_abundance_data_list <- function(full_data_list) {

  data_needed <- c("Nstrata", "Nstrata_wc", "lt_pCap_mu", "lt_pCap_sd",
                   "u", "Uwc_ind", "K", "ZP", "lgN_max")

  new_data_list <- full_data_list[sapply(names(full_data_list), function(x) x %in% data_needed)]
  return(new_data_list)
}

#' Fit pCap model
#' @details This function prepares the data list for the abundance STAN model based on what the model name is.
#' @param input A list containing the inputs for the pCap model.
#' @returns a named list with the required elements for that model run.
#' @keywords internal
#' @export
#' @md
fit_pCap_model <- function(input) {

  # model_name <- paste0("pCap_", input$model_name)
  # print(model_name)
  #
  # stan_model <- eval(parse(text = paste0("SRJPEmodel::bt_spas_x_model_code$", model_name)))
  stan_model <- eval(parse(text = "SRJPEmodel::bt_spas_x_model_code$pCap_all"))

  options(mc.cores=parallel::detectCores())
  cli::cli_alert("running pCap model")
  pcap <- rstan::stan(model_code = stan_model,
                      data = input$data,
                      init = input$inits,
                      # do not save logit_pCap or pro_dev_P (way too big)
                      pars = c("logit_pCap", "b0_pCap", "b_flow", "pro_sd_P", "trib_mu_P", "trib_sd_P",
                               "flow_mu_P", "flow_sd_P"), # for new efficient model
                      # pars = c("lt_pCap_U", "sim_pro_dev", "b0_pCap", "b_flow", "trib_mu_P", "trib_tau_P",
                      #          "flow_mu_P", "flow_tau_P", "pro_tau_P"),
                      chains = SRJPEmodel::bt_spas_x_bayes_params$number_chains,
                      iter = SRJPEmodel::bt_spas_x_bayes_params$number_mcmc,
                      seed = 84735)

  return(pcap)
}

#' Generate lt_pCap_U values from the pCap model object
#' @details This function prepares the data list for the abundance STAN model based on what the model name is.
#' @param pCap_inputs A list containing the inputs for the pCap model including `model_name`, `Nstrata`,
#' `catch_flow`, `use_trib`, `Ntribs`, `Nmr`, `Nwmr`, `Nwomr`, `Uind_wMR`, `Uind_woMR`,
#' `ind_pCap`.
#' @param pCap_model_object A STANfit object with the output of the pCap hierarchical model.
#' @returns a named list with the required elements for that model run.
#' @export
#' @md
generate_lt_pCap_Us <- function(abundance_inputs, pCap_model_object){

  # set up objects
  ModelName <- abundance_inputs$model_name
  Nstrata <- abundance_inputs$inputs$data$Nstrata
  catch_flow <- abundance_inputs$lt_pCap_U_data$catch_flow
  use_trib <- abundance_inputs$lt_pCap_U_data$use_trib
  Nwmr <- abundance_inputs$lt_pCap_U_data$Nwmr
  Nwomr <- abundance_inputs$lt_pCap_U_data$Nwomr
  Uind_wMR <- abundance_inputs$lt_pCap_U_data$Uind_wMR
  Uind_woMR <- abundance_inputs$lt_pCap_U_data$Uind_woMR
  ind_pCap <- abundance_inputs$lt_pCap_U_data$ind_pCap

  # extract from pCap model fit object

  samples <- rstan::extract(pCap_model_object, pars = c("logit_pCap", "b0_pCap", "b_flow", "pro_sd_P", "trib_mu_P", "trib_sd_P",
                                           "flow_mu_P", "flow_sd_P"),
                            permuted = TRUE)
  Ntrials <- dim(samples$logit_pCap)[1] # of saved posterior samples from pCap model in stan
  logit_pCap <- samples$logit_pCap # logit_pCap[1:Ntrials,1:Nmr] # The estimated logit pCap posterior for each efficiency trial
  b0_pCap <- samples$b0_pCap # b0_pCap[1:Ntrials,1:Ntribs] #mean logit pCap for each site (at mean discharge)
  b_flow <- samples$b_flow # b_flow[1:Ntrials,1:Ntribs] #flow effect for each site
  pro_sd_P <- samples$pro_sd_P # pro_sd_P[1:Ntrials]        #process error (sd)
  trib_mu_P <- samples$trib_mu_P # trib_mu_P[1:Ntrials] #hyper mean for b0_pCap
  trib_sd_P <- samples$trib_sd_P # trib_sd_P[1:Ntrials] #hyper sd for b0_pCap
  flow_mu_P <- samples$flow_mu_P # flow_mu_P[1:Ntrials] #hyper mean for b_flow
  flow_sd_P <- samples$flow_sd_P # flow_sd_P[1:Ntrials] #hyper sd for b_flow

  # calculations
  lt_pCap_U=matrix(nrow=Ntrials,ncol=Nstrata)
  sim_pro_dev=vector(length=Ntrials)
  lt_pCap_mu=matrix(nrow=Nstrata,ncol=Ntrials) #function needs to return this
  lt_pCap_sd=lt_pCap_mu                        #function needs to return this

  if(ModelName=="all_mark_recap" ){

    for(i in 1:Nstrata){
      lt_pCap_U[,i] = logit_pCap[,ind_pCap[i]];
    }

  } else if (ModelName=="missing_mark_recap"){

    for(i in 1:Nwmr){
      #Assign estimated pCaps for strata with efficiency data
      lt_pCap_U[,Uind_wMR[i]] = logit_pCap[,ind_pCap[i]];
    }
    for (i in 1:Nwomr) {
      #for weeks without efficiency trials
      for(itrial in 1:Ntrials) sim_pro_dev[itrial] = rnorm(n=1, mean=0,sd=pro_sd_P[itrial]);
      lt_pCap_U[,Uind_woMR[i]] = b0_pCap[,use_trib] + b_flow[,use_trib] * catch_flow[Uind_woMR[i]] + sim_pro_dev[1:Ntrials]
    }

  } else if (ModelName=="no_mark_recap"){

    for (i in 1:Nwomr) {
      for(itrial in 1:Ntrials) sim_pro_dev[itrial] = rnorm(n=1, mean=0, sd=pro_sd_P[itrial]);
      lt_pCap_U[,Uind_woMR[i]] = b0_pCap[,use_trib] + b_flow[,use_trib] * catch_flow[Uind_woMR[i]] + sim_pro_dev[1:Ntrials]
    }

  } else if (ModelName=="no_mark_recap_no_trib"){

    logit_b0=vector(length=Ntrials); logit_bflow=logit_b0
    for(itrial in 1:Ntrials){
      logit_b0[itrial] = rnorm(n=1, mean=trib_mu_P[itrial], sd=trib_sd_P[itrial])
      logit_bflow[itrial] = rnorm(n=1, mean=flow_mu_P[itrial], sd=flow_sd_P[itrial])
    }

    for (i in 1:Nwomr) {
      for(itrial in 1:Ntrials) sim_pro_dev[itrial] = rnorm(n=1, mean=0, sd=pro_sd_P[itrial]);
      lt_pCap_U[,Uind_woMR[i]] = logit_b0 + logit_bflow * catch_flow[Uind_woMR[i]] + sim_pro_dev[1:Ntrials];
    }

  }#end if on ModelName


  #Calculate mean and sd for each lt_pCap_U and return these from function
  for(i in 1:Nstrata){
    lt_pCap_mu[i,]=mean(lt_pCap_U[,i])
    lt_pCap_sd[i,]=sd(lt_pCap_U[,i])
    # lt_pCap_mu[,i]=mean(lt_pCap_U[i])
    # lt_pCap_sd[,i]=sd(lt_pCap_U[i])
  }

  return(list("lt_pCap_mu" = lt_pCap_mu |> rowMeans(),
              "lt_pCap_sd" = lt_pCap_sd |> rowMeans()))
}


#' Run BUGS abundance model.
#' @details This function runs the BUGS abundance model.
#' @param abundance_inputs the object produced by `prepare_abundance_inputs()`
#' @param lt_pCap_Us the object produced by `generate_lt_pCap_Us()`
#' @param bugs_model_file the filepath pointing to where your BUGS abundance model file is
#' @param bugs_directory the filepath pointing to where your WinBUGS14/ directory is
#' @returns a BUGS object.
#' @export
#' @md
fit_abundance_model_BUGS <- function(abundance_inputs, lt_pCap_Us,
                                     bugs_model_file,
                                     bugs_directory) {

  parameters <- c("lt_pCap_U", "pCap_U", "N", "Ntot", "sd.N", "sd.Ne", "lg_CumN")
  Nmcmc = 2000
  Nburnin = 500
  Nthin = 2
  Nchains = 3

  # clean up later
  data <- abundance_inputs$inputs$data
  data$lt_pCap_mu <- lt_pCap_Us$lt_pCap_mu
  data$lt_pCap_tau <- 1/lt_pCap_Us$lt_pCap_sd^2

  inits_with_lt_pCap_U <- abundance_inputs$inputs$inits[[1]]

  inits_with_lt_pCap_U$lt_pCap_U <- data$lt_pCap_mu
  inits_with_lt_pCap_U$tau.N <- 1
  inits_with_lt_pCap_U$tau.Ne <- 1

  new_inits <- list(inits_with_lt_pCap_U,
                    inits_with_lt_pCap_U,
                    inits_with_lt_pCap_U)

  abundance <- bugs(data,
                    new_inits,
                    parameters,
                    bugs_model_file, debug = F,
                    n.chains = Nchains, n.burnin = Nburnin, n.thin = Nthin, n.iter = Nmcmc,
                    codaPkg = F, DIC = T, clearWD = T,
                    bugs.directory = bugs_directory)
                    #bugs.directory = here::here("data-raw", "WinBUGS14"))

  return(abundance)
}



#' Fit abundance model
#' @details This function prepares the data list for the abundance STAN model based on what the model name is.
#' @param input A list containing the inputs for the abundance model.
#' @param pCap_fit A STANfit object resulting from running `fit_pCap_model()`
#' @returns a named list with the required elements for that model run.
#' @keywords internal
#' @export
#' @md
fit_abundance_model <- function(input, pCap_fit, lt_pCap_Us) {

  # get lt_pCap_Us for data from pcap fit
  # logit_pCaps <- rstan::summary(pCap_fit, pars = c("lt_pCap_U"))$summary |>
  #   data.frame()

  input$data$lt_pCap_mu <- lt_pCap_Us$lt_pCap_mu
  input$data$lt_pCap_sd <- lt_pCap_Us$lt_pCap_sd

  # generate inits for lt_pCap_U
  inits_with_lt_pCap_U <- input$inits[[1]]
  inits_with_lt_pCap_U$lt_pCap_U <- input$data$lt_pCap_mu
  new_inits <- list(inits = inits_with_lt_pCap_U,
                    inits = inits_with_lt_pCap_U,
                    inits = inits_with_lt_pCap_U)

  stan_model <- eval(parse(text = "SRJPEmodel::bt_spas_x_model_code$abundance"))

  options(mc.cores=parallel::detectCores())

  cli::cli_alert("running abundance model")
  abundance <- rstan::stan(model_code = stan_model,
                           data = input$data,
                           # algorithm = "HMC",
                           init = new_inits,
                           chains = SRJPEmodel::bt_spas_x_bayes_params$number_chains,
                           iter = SRJPEmodel::bt_spas_x_bayes_params$number_mcmc,
                           seed = 84735)

  return(abundance)
}

#' BT SPAS X diagnostic plots
#' @details This function produces a plot with data and results of fitting the pCap and abundance models.
#' @param site_arg The site being fit
#' @param run_year_arg The run year being fit
#' @param abundance_model The STANfit model object produced by running `fit_abundance_model()`
#' @returns A plot
#' @export
#' @md
diagnostic_plots_split <- function(site_arg, run_year_arg,
                                   abundance_model) {


  julian_week_to_date_lookup <- SRJPEmodel::julian_week_to_date_lookup

  data <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
    filter(life_stage != "yearling") |>
    filter(run_year == run_year_arg,
           site == site_arg,
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
    left_join(SRJPEdata::weekly_juvenile_abundance_efficiency_data,
              by = c("year", "run_year", "week", "stream", "site")) |>
    mutate(count = round(count, 0),
           catch_standardized_by_hours_fished = round(catch_standardized_by_hours_fished, 0),
           # change all NaNs to NAs
           across(mean_fork_length:lgN_prior, ~ifelse(is.nan(.x), NA, .x)),
           # plot things
           lincoln_peterson_abundance = count * (number_released / number_recaptured),
           lincoln_peterson_efficiency = number_recaptured / number_released) |>
    left_join(julian_week_to_date_lookup, by = c("week" = "Jwk")) |>
    mutate(year = ifelse(week >= 43, run_year - 1, run_year),
           fake_date = ymd(paste0(year, "-01-01")),
           final_date = fake_date + weeks(week - 1),
           date = format(final_date, "%b-%d"),
           week_index = row_number())

  pCap_estimates <- rstan::summary(abundance_model,pars=c("lt_pCap_U"))$summary |>
    data.frame() |>
    tibble::rownames_to_column("parameter") |>
    mutate(week_index = readr::parse_number(parameter)) |>
    select(week_index, mean, sd, `50` = X50., `2.5` = X2.5., `97.5` = X97.5.) |>
    mutate(parameter = "lt_pCap_U")

  N_estimates <- rstan::summary(abundance_model,pars=c("N"))$summary |>
    data.frame() |>
    tibble::rownames_to_column("parameter") |>
    mutate(week_index = readr::parse_number(parameter)) |>
    select(week_index, mean, sd, `50` = X50., `2.5` = X2.5., `97.5` = X97.5.) |>
    mutate(parameter = "N")


  data |>
    ggplot(aes(x = final_date, y = catch_standardized_by_hours_fished)) +
    geom_bar(stat = "identity", fill = "grey", width = 5) +
    scale_x_date(date_breaks = "1 week", date_labels = "%b %d") +
    #scale_x_date(date_breaks = "1 week", date_labels = "%V") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(x = "Date",
         y = "Raw catch",
         title = paste0("Catch at ", site_arg, " for run year ",
                        run_year_arg))

  # abundance plot
  abundance_plot <- data |>
    left_join(N_estimates) |>
    #mutate(across(c(mean, `50`, `2.5`, `97.5`), plogis)) |>
    ggplot(aes(x = final_date, y = `50`)) +
    geom_bar(stat = "identity", fill = "grey", width = 5) +
    geom_errorbar(aes(x = final_date, ymin = `2.5`, ymax = `97.5`), width = 0.2) +
    geom_point(aes(x = final_date, y = lincoln_peterson_abundance),
               shape = 1, color = "blue") +
    # geom_point(aes(x = fake_date, y = Inf, color = sampled),
    #            size = 3) +
    geom_text(aes(x = final_date, y = Inf,
                  label = paste(count),
                  angle = 90),
              hjust = 1,
              size = 3) +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "red")) +
    theme_minimal() +
    labs(x = "",
         #x = "Date",
         y = "Abundance",
         title = paste(site_arg, run_year_arg)) +
    #theme(axis.text.x=element_blank()) +
    scale_x_date(date_breaks = "1 week", date_labels = "%b %d") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  # efficiency
  efficiency_plot <- data |>
    left_join(pCap_estimates) |>
    mutate(across(c(mean, `50`, `2.5`, `97.5`), plogis),
           number_released_label = ifelse(is.na(number_released), "", number_released),
           number_recaptured_label = ifelse(is.na(number_recaptured), "", number_recaptured)) |>
    ggplot(aes(x = final_date, y = `50`)) +
    geom_bar(stat = "identity", fill = "grey", width = 4) +
    geom_errorbar(aes(x = final_date, ymin = `2.5`, ymax = `97.5`), width = 0.2) +
    geom_point(aes(x = final_date, y = lincoln_peterson_efficiency),
               shape = 1, color = "blue") +
    geom_text(aes(x = final_date, y = Inf,
                  label = paste(number_released_label, number_recaptured_label),
                  angle = 90),
              hjust = 1,
              size = 3) +
    # geom_point(aes(x = fake_date, y = Inf, color = sampled),
    #            size = 3) +
    theme_minimal() +
    labs(x = "Date", y = "Weekly Efficiency") +
    scale_x_date(date_breaks = "1 week", date_labels = "%b %d") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  # arrange together
  gridExtra::grid.arrange(abundance_plot, efficiency_plot)
}


#' BT SPAS X raw data plots
#' @details This function produces a plot with data used to fit BT SPAS X.
#' @param site_arg The site being fit
#' @param run_year_arg The run year being fit
#' @returns A plot
#' @export
#' @md
plot_juv_data <- function(site, run_year) {
  data <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
    filter(life_stage != "yearling") |>
    filter(run_year == !!run_year,
           site == !!site,
           week %in% c(seq(45, 53), seq(1, 22))) |>
    group_by(year, week, stream, site, run_year) |>
    # keep NAs in count columns
    summarise(count = if(all(is.na(count))) NA_real_ else sum(count, na.rm = TRUE),
              mean_fork_length = mean(mean_fork_length, na.rm = T),
              hours_fished = mean(hours_fished, na.rm = T),
              catch_flow_cfs = mean(flow_cfs, na.rm = T),
              average_stream_hours_fished = mean(average_stream_hours_fished, na.rm = T),
              standardized_flow = mean(standardized_flow, na.rm = T),
              catch_standardized_by_hours_fished = if(all(is.na(catch_standardized_by_hours_fished))) NA_real_ else sum(catch_standardized_by_hours_fished, na.rm = TRUE),
              lgN_prior = mean(lgN_prior, na.rm = T)) |>
    ungroup() |>
    left_join(SRJPEdata::weekly_juvenile_abundance_efficiency_data,
              by = c("year", "run_year", "week", "stream", "site")) |>
    mutate(count = round(count, 0),
           catch_standardized_by_hours_fished = round(catch_standardized_by_hours_fished, 0),
           # change all NaNs to NAs
           across(mean_fork_length:lgN_prior, ~ifelse(is.nan(.x), NA, .x)),
           # plot things
           lincoln_peterson_abundance = count * (number_released / number_recaptured),
           lincoln_peterson_efficiency = number_recaptured / number_released) |>
    left_join(SRJPEmodel::julian_week_to_date_lookup, by = c("week" = "Jwk")) |>
    mutate(year = ifelse(week >= 43, run_year - 1, run_year),
           fake_date = ymd(paste0(year, "-01-01")),
           final_date = fake_date + weeks(week - 1),
           date = format(final_date, "%b-%d"),
           week_index = row_number())

  abundance_plot <- data |>
    ggplot(aes(x = final_date, y = count)) +
    geom_bar(stat = "identity", fill = "grey", width = 5) +
    # geom_point(aes(x = final_date, y = lincoln_peterson_abundance),
    #            shape = 1, color = "blue") +
    geom_text(aes(x = final_date, y = Inf,
                  label = paste(count),
                  angle = 90),
              hjust = 1,
              size = 3) +
    theme_minimal() +
    labs(x = "",
         #x = "Date",
         y = "Abundance",
         title = paste(site, run_year)) +
    #theme(axis.text.x=element_blank()) +
    scale_x_date(date_breaks = "1 week", date_labels = "%b %d") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  efficiency_plot <- data |>
    mutate(number_released_label = ifelse(is.na(number_released), "", number_released),
           number_recaptured_label = ifelse(is.na(number_recaptured), "", number_recaptured)) |>
    ggplot(aes(x = final_date, y = lincoln_peterson_efficiency)) +
    geom_point(shape = 1, color = "blue") +
    # geom_bar(stat = "identity", fill = "grey", width = 4) +
    geom_text(aes(x = final_date, y = Inf,
                  label = paste(number_released_label, number_recaptured_label),
                  angle = 90),
              hjust = 1,
              size = 3) +
    theme_minimal() +
    labs(x = "Date", y = "Weekly Efficiency") +
    scale_x_date(date_breaks = "1 week", date_labels = "%b %d") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  # arrange together
  gridExtra::grid.arrange(abundance_plot, efficiency_plot)

}

