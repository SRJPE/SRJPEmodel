#' Prepare inputs for pCap STAN model
#' @details This function prepares data for input into a pCap STAN model.
#' @param mainstem Whether or not you want to evaluate efficiency trials for a mainstem site
#' (`knights landing`, `tisdale`, and `red bluff diversion dam`) or for a tributary site. If `FALSE`,
#' the mark recapture dataset will be filtered to exclude those mainstem sites.
#' @param input_catch_data Optional argument for weekly catch data.
#' Defaults to `SRJPEdata::weekly_juvenile_abundance_catch_data`. If
#' passed in, structure of data frame must match that of the default.
#' @param input_efficiency_data Optional argument for weekly efficiency data.
#' Defaults to `SRJPEdata::SRJPEdata::weekly_juvenile_abundance_efficiency_data`. If
#' passed in, structure of data frame must match that of the default.
#' @returns a list:
#' * **pCap_inputs** a list of data and inits for input into the pCap model.
#' * **sites_fit** a list of site names associated with `ind_trib`.
#' @export
#' @md
prepare_pCap_inputs <- function(mainstem = c(FALSE, TRUE),
                                input_catch_data = NULL,
                                input_efficiency_data = NULL) {

  if(missing(input_catch_data)) {
    input_catch_data <- SRJPEdata::weekly_juvenile_abundance_catch_data
  }

  if(missing(input_efficiency_data)) {
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
#' @param input_catch_data Optional argument for weekly catch data.
#' Defaults to `SRJPEdata::weekly_juvenile_abundance_catch_data`. If
#' passed in, structure of data frame must match that of the default.
#' @param input_efficiency_data Optional argument for weekly efficiency data.
#' Defaults to `SRJPEdata::SRJPEdata::weekly_juvenile_abundance_efficiency_data`. If
#' passed in, structure of data frame must match that of the default.
#' @returns a list:
#' * **abundance_inputs** a list of data and inits for input into the abundance model.
#' * **model_name** The version of the pCap logic to use for generating lt_pCap_U values. Either
#' `all_mark_recap`, `missing_mark_recap`, `no_mark_recap`, or `no_mark_recap_no_trib`.
#' * **weeks_fit** The weeks fit for the abundance model
#' * **weeks_date** Associated dates for the weeks fit for the abundance model.
#' @export
#' @md
prepare_abundance_inputs <- function(site, run_year,
                                     effort_adjust = c(T, F),
                                     input_catch_data = NULL,
                                     input_efficiency_data = NULL) {

  if(missing(input_catch_data)) {
    input_catch_data <- SRJPEdata::weekly_juvenile_abundance_catch_data
  }

  if(missing(input_efficiency_data)) {
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
  # TODO iterative improvement: set the denominator to be the lower quantile just for lgN_max as an argument
  # data input
  lgN_max = rep(log(0.001 * (mean(weekly_catch_data, na.rm=T) + 1) / 0.005), number_weeks_catch)

  for(j in 1:number_weeks_with_catch){
    if(is.na(weekly_catch_data[j]) == F) lgN_max[indices_with_catch[j]] = log(0.001 * (weekly_catch_data[j] + 1) / 0.005)
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

#' Prepare spline parameters object for BT-SPAS-X model
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


#' Fit pCap model
#' @details This function prepares the data list for the abundance STAN model based on what the model name is.
#' @param input A list containing the inputs for the pCap model.
#' @returns a STANfit object with the pCap model fit.
#' @export
#' @md
fit_pCap_model <- function(input) {

  stan_model <- eval(parse(text = "SRJPEmodel::bt_spas_x_model_code$pCap_all"))

  options(mc.cores=parallel::detectCores())
  cli::cli_alert("running pCap model")
  pcap <- rstan::stan(model_code = stan_model,
                      data = input$data,
                      init = input$inits,
                      # do not save logit_pCap or pro_dev_P (way too big)
                      pars = c("logit_pCap", "b0_pCap", "b_flow", "pro_sd_P", "trib_mu_P", "trib_sd_P",
                               "flow_mu_P", "flow_sd_P"), # for new efficient model
                      chains = SRJPEmodel::bt_spas_x_bayes_params$number_chains,
                      iter = SRJPEmodel::bt_spas_x_bayes_params$number_mcmc,
                      seed = 84735)

  return(pcap)
}

#' Generate lt_pCap_U values from the pCap model object
#' @details This function prepares the data list for the abundance STAN model based on what the model name is.
#' @param abundance_inputs A list containing the inputs for the site-specific abundance model including
#' `model_name`, `Nstrata`, `catch_flow`, `use_trib`, `Ntribs`, `Nmr`, `Nwmr`, `Nwomr`, `Uind_wMR`, `Uind_woMR`,
#' `ind_pCap`. This is created by calling `prepare_abundance_inputs()`.
#' @param pCap_model_object A STANfit object with the output of the pCap hierarchical model. This is created
#' by calling `fit_pCap_model()`.
#' @returns a named list with the simulated `lt_pCap_mu` and `lt_pCap_sd` values for
#' `Nstrata`.
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


#' Fit abundance model in BUGS.
#' @details This function runs the BUGS abundance model.
#' @param abundance_inputs the object produced by `prepare_abundance_inputs()`
#' @param lt_pCap_Us the object produced by `generate_lt_pCap_Us()` containing values of
#' length `Nstrata` for `lt_pCap_mu` and `lt_pCap_sd`.
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



#' Fit abundance model in STAN
#' @details This function calls a STAN model for abundance. It is currently in development.
#' @param input A list containing the inputs for the abundance model. This is created by calling
#' `prepare_abundance_inputs()`.
#' @param pCap_fit A STANfit object resulting from running `fit_pCap_model()`.
#' @param lt_pCap_Us A named list produced by running `generate_lt_pCap_Us()`.
#' @returns a STANfit object containing abundance model fits.
#' @md
fit_abundance_model_STAN <- function(input, pCap_fit, lt_pCap_Us) {

  # TODO this is in progress

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

#' Extract parameter estimates from abundance BUGS model object
#' @details This function extracts parameter estimates from the abundance BUGS model and returns them in a tidy data frame format.
#' @param site the site fit.
#' @param run_year the run year fit.
#' @param abundance_inputs a list of inputs for the abundance BUGS model, generated by running `prepare_abundance_inputs()`
#' @param model_object the BUGS object produced by running `fit_abundance_model_BUGS()`.
#' @returns A table with the format:
#' * **model_name**
#' * **site**
#' * **run_year**
#' * **week_fit**
#' * **location_fit**
#' * **parameter**
#' * **statistic**
#' * **value**
#' * **srjpedata_version**
#' @export
#' @md
extract_abundance_estimates <- function(site, run_year,
                                        abundance_inputs,
                                        model_object) {

  # link to actual weeks
  # TODO do we want week formatted as MONTH-DATE ? can do easily
  week_lookup <- tibble("week_fit" = abundance_inputs$weeks_fit) |>
    mutate(week_index = row_number())

  formatted_table <- model_object$summary |>
    as.data.frame() |>
    tibble::rownames_to_column("parameter") |>
    mutate(week_index = ifelse(str_detect(parameter, "b0_pCap|b_flow"), NA,
                               suppressWarnings(readr::parse_number(parameter))),
           site = site,
           run_year = run_year,
           model_name = abundance_inputs$model_name,
           parameter = gsub("[0-9]+|\\[|\\]", "", parameter),
           srjpedata_version = as.character(packageVersion("SRJPEdata"))) |>
    left_join(week_lookup, by = "week_index") |>
    rename(rhat = Rhat) |>
    # now clean up statistics
    pivot_longer(mean:rhat,
                 values_to = "value",
                 names_to = "statistic") |>
    mutate(statistic = str_remove_all(statistic, "\\%"),
           location_fit = NA) |> # this will be necessary for the pCap model fit object
    select(model_name, site, run_year, week_fit, location_fit, parameter, statistic, value, srjpedata_version)

  return(formatted_table)
}


#' Extract parameter estimates from the pCap BUGS model object
#' @details This function extracts parameter estimates from the pCap BUGS model and returns them in a tidy data frame format.
#' @param pCap_inputs a list of inputs for the abundance BUGS model, generated by running `prepare_pCap_inputs()`
#' @param model_object the BUGS object produced by running `fit_pCap_model_BUGS()`.
#' @returns A table with the format:
#' * **model_name**
#' * **site**
#' * **run_year**
#' * **week_fit**
#' * **location_fit**
#' * **parameter**
#' * **statistic**
#' * **value**
#' * **srjpedata_version**
#' @export
#' @md
extract_pCap_estimates <- function(model_object, pCap_inputs) {

  site_lookup <- tibble("location_fit" = pCap_inputs$sites_fit) |>
    mutate(site_index = row_number())

  formatted_table <- rstan::summary(model_object)$summary |>
    as.data.frame() |>
    tibble::rownames_to_column("parameter") |>
    mutate(week_index = NA, # no weekly estimates that are relevant
           site_index = ifelse(str_detect(parameter, "b0_pCap|b_flow"),
                               suppressWarnings(readr::parse_number(substr(parameter, 3, length(parameter)))),
                               NA),
           site = NA,
           run_year = NA,
           model_name = "pCap_all",
           parameter = gsub("[0-9]+|\\[|\\]", "", parameter),
           srjpedata_version = as.character(packageVersion("SRJPEdata"))) |>
    left_join(site_lookup, by = "site_index") |>
    rename(rhat = Rhat) |>
    # now clean up statistics
    pivot_longer(mean:rhat,
                 values_to = "value",
                 names_to = "statistic") |>
    mutate(statistic = str_remove_all(statistic, "\\%"),
           week_fit = NA) |> # this won't be reported for the pCap model
    select(model_name, site, run_year, week_fit, location_fit, parameter, statistic, value, srjpedata_version)

  return(formatted_table)
}

#' Run all JPE sites.
#' @details This function automates running BT-SPAS-X for all JPE sites/run years.
#' @param sites_to_run
#' @param run_pCap
#' @param pCap_model_object_filepath
#' @param bugs_model_file
#' @param bugs_directory
#' @returns A table with the format:
#' * **model_name**
#' * **site**
#' * **run_year**
#' * **week_fit**
#' * **location_fit**
#' * **parameter**
#' * **statistic**
#' * **value**
#' * **srjpedata_version**
#' @export
#' @md
run_bt_spas_x_JPE_sites <- function(sites_to_run,
                                    run_pCap = FALSE,
                                    pCap_model_object_filepath = NULL,
                                    bugs_model_file,
                                    bugs_directory) {

  # run pCap model if necessary
  if(run_pCap) {
    pCap_inputs <- prepare_pCap_inputs(mainstem = FALSE)
    pCap <- fit_pCap_model(pCap_inputs$inputs)
  } else {
    pCap <- readRDS(pCap_model_object_filepath)
  }

  # prep inputs as vectors
  sites_to_run_inputs <- sites_to_run |>
    mutate(bugs_model_file = bugs_model_file,
           bugs_directory = bugs_directory)

  # now run abundance workflow
  SRJPE_fits_table <- purrr::pmap(list(sites_to_run_inputs$site,
                                       sites_to_run_inputs$run_year,
                                       sites_to_run_inputs$bugs_model_file,
                                       sites_to_run_inputs$bugs_directory),
                                  run_abundance_workflow,
                                  .progress = TRUE)

  all_JPE_sites_clean <- SRJPE_fits_table |>
    bind_rows()

  # extract clean table
  return(all_JPE_sites_clean)

}

#' Run abundance BUGS workflow.
#' @details This function runs the abundance BUGS workflow.
#' @param site
#' @param run_year
#' @param pCap
#' @param bugs_model_file
#' @param bugs_directory
#' @returns A table.
#' @export
#' @md
run_abundance_workflow <- function(site,
                                   run_year,
                                   pCap,
                                   bugs_model_file,
                                   bugs_directory) {

  abundance_inputs <- prepare_abundance_inputs(site, run_year, effort_adjust = T)

  lt_pCap_Us <- generate_lt_pCap_Us(abundance_inputs, pCap)

  abundance <- tryCatch({fit_abundance_model_BUGS(abundance_inputs, lt_pCap_Us,
                                        bugs_model_file,
                                        bugs_directory)
  },
    error = function(e) return(tibble("error" = TRUE))
  )

  clean_table <- tryCatch({extract_abundance_estimates(site, run_year,
                                                       abundance_inputs, abundance)
  },
  error = function(e) return(tibble("error" = TRUE))
  )

  return(clean_table)
}



#' BT SPAS X diagnostic plots
#' @details This function produces a plot with data and results of fitting the pCap and abundance models for
#' a given site and run year.
#' @param site_arg The site being fit
#' @param run_year_arg The run year being fit
#' @param abundance_model The STANfit model object produced by running `fit_abundance_model_BUGS()`
#' @returns A plot.
#' @export
#' @md
generate_diagnostic_plot_juv <- function(site_arg, run_year_arg,
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

  # TODO replace with extract() function
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
#' @details This function produces a plot with data for a site and run year.
#' @param site The site being fit
#' @param run_year The run year being fit
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

