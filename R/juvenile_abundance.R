#' Prepare inputs for pCap STAN model
#' @details This function prepares data for input into a pCap STAN model.
#' @param model_type either `one_site` or `all_sites`
#' @param skew either `TRUE` or `FALSE`. Only applies for `one_site` models.
#' @param site_selection if fitting the `one_site` model, the tributary you want to fit.
#' @param exclude_from_hyper a vector of sites you want to drop from being estimated in the hyper distribution
#' parameters `b0_pCap` and `b_flow`. Default is NULL, meaning no sites are dropped. Must be site names
#' in `unique(SRJPEdata::weekly_juvenile_abundance_efficiency_data$site)`
#' @param input_catch_data Optional argument for weekly catch data.
#' Defaults to `SRJPEdata::weekly_juvenile_abundance_catch_data`. If
#' passed in, structure of data frame must match that of the default.
#' @param input_efficiency_data Optional argument for weekly efficiency data.
#' Defaults to `SRJPEdata::SRJPEdata::weekly_juvenile_abundance_efficiency_data`. If
#' passed in, structure of data frame must match that of the default.
#' @returns a list:
#' * **pCap_inputs** a list of data and inits for input into the pCap model.
#' * **sites_fit** a list of site names associated with `ind_trib`.
#' @family Prepare Model Inputs
#' @export
#' @md
prepare_pCap_inputs <- function(model_type = c("all_sites", "one_site"),
                                skew = NULL,
                                site_selection = NULL,
                                exclude_from_hyper = NULL,
                                input_catch_data = NULL,
                                input_efficiency_data = NULL) {

  # check arguments are logical
  if(!is.null(skew) & model_type == "all_sites") {
    cli::cli_warn("Skew model cannot be run for the all_sites model type.")
    return(invisible(NULL))
  }
  if(any(!exclude_from_hyper %in% unique(SRJPEdata::weekly_juvenile_abundance_efficiency_data$site))){
    cli::cli_warn("Non-valid tributary passed to exclude_from_hyper.")
    return(invisible(NULL))
  }
  if(model_type == "one_site" & length(site_selection) != 1){
    cli::cli_warn("One site name must be passed if running one_site model type.")
    return(invisible(NULL))
  }
  if(model_type == "one_site" & is.null(skew)){
    cli::cli_warn("skew must be specified as TRUE/FALSE when running a one_site model type.")
    return(invisible(NULL))
  }

  # default to the SRJPEdata objects
  if(missing(input_catch_data)) {
    input_catch_data <- SRJPEdata::weekly_juvenile_abundance_catch_data
  }
  if(missing(input_efficiency_data)) {
    input_efficiency_data <- SRJPEdata::weekly_juvenile_abundance_efficiency_data
  }

  # available sites for fitting
  available_sites <- c("ubc", "okie dam", "lcc", "ucc", "deer creek", "eye riffle",
                       "gateway riffle", "herringer riffle", "live oak",
                       "steep riffle", "sunset pumps", "mill creek", "hallwood")
  mainstem_sites <- c("red bluff diversion dam", "knights landing", "tisdale")

  # one_site data filtering and parameter names
  if(model_type == "one_site") {
    filtered_efficiency_data <- input_efficiency_data |>
      filter(site == site_selection)
    if(skew) {
      # include alpha
      pCap_parameters <- c("logit_pCap", "b0_pCap", "b_flow","pro_sd_P","yr_sd_P","yr_re", "alpha") #, "b_eff")
      model_name <- "pCap_one_site_skew"
    } else {
      # include alpha
      pCap_parameters <- c("logit_pCap", "b0_pCap", "b_flow","pro_sd_P","yr_sd_P","yr_re") #, "b_eff")
      model_name <- "pCap_one_site"
    }

  } else {
    # filter to tributary sites (can't include mainstem sites in this)
    filtered_efficiency_data <- input_efficiency_data |>
      filter(site %in% available_sites)
    pCap_parameters <- c("logit_pCap","trib_mu_P", "trib_sd_P", "flow_mu_P", "flow_sd_P", "pro_sd_P",
                         "b0_pCap", "b_flow","yr_sd_P","yr_re") # , "b_eff")
    model_name <- "pCap_all_sites"
  }

  # prepare efficiency dataset post-filtering
  mark_recapture_data <- filtered_efficiency_data |>
    # grab standardized_flow
    left_join(input_catch_data |>
                select(year, week, stream, site, run_year, standardized_flow),
              by = c("year", "week", "run_year", "stream", "site")) |>
    # or do we want to filter just no number released?
    dplyr::filter(!is.na(standardized_efficiency_flow),
                  !is.na(number_released),
                  !is.na(number_recaptured)) |>
    # right now there's lifestage in the dataset, so we have to do dplyr::distinct()
    # TODO Check that we need this line
    # dplyr::distinct(site, run_year, week, number_released, number_recaptured, .keep_all = TRUE) |>
    # sort by site
    group_by(site) |>
    left_join(site_order_north_south, by = "site") |>
    arrange(ns_order, year, week) |>
    ungroup()

  if(any(mark_recapture_data$number_recaptured > mark_recapture_data$number_released)) {
    problem_data <- mark_recapture_data |>
      dplyr::filter(number_recaptured > number_released)

    cli::cli_alert_info(paste0(nrow(problem_data), " rows of your data have more Recaptures than Releases for
                          a given week. Filtering out the problematic data for now.
                          Please check your data."))

    mark_recapture_data <- mark_recapture_data |>
      dplyr::filter(number_recaptured <= number_released)
  }

  # prepare effort
  effort <- mark_recapture_data$hours_fished / mark_recapture_data$average_hours_fished_during_efficiency_trials
  effort[effort == 0] <- 0.001
  effort[is.na(effort)] <- 0.001

  # now prepare indexes for use in the model

  # first assign 1:Ntribs to the unique sites in the dataset
  site_lookup <- mark_recapture_data |>
    dplyr::distinct(site) |>
    mutate(ID = row_number())

  # pull the order of those sites
  sites_fit <- site_lookup |>
    dplyr::pull(site)

  # year lookup
  mr_year_lookup <- mark_recapture_data |>
    arrange(ns_order, run_year) |>
    distinct(site, run_year) |>
    mutate(site_run_year_id = row_number())

  site_year_fit <- unique(mr_year_lookup[c("site", "run_year")])
  years_fit <- mr_year_lookup |>
    dplyr::pull(run_year)

  sd_yr_ind <- site_year_fit |>
    mutate(sd_yr_ind = as.integer(factor(site, levels = unique(site)))) |>
    pull(sd_yr_ind)

  # assign the IDs to the sites in the mark-recapture dataset
  mark_recapture_data <- mark_recapture_data |>
    left_join(site_lookup, by = "site") |>
    left_join(mr_year_lookup, by = c("site", "run_year"))

  # get indexing for "mark recap" dataset (pCap model)
  Ntribs <- length(sites_fit) # number of sites (for pCap calculations)
  number_efficiency_experiments <- unique(mark_recapture_data[c("site", "run_year", "week")]) |>
    nrow() # number of efficiency experiments completed, nrow(mark_recapture_data) this depends on whether you have lifestage or not

  # test mr_flow replace
  clean_mr_flow <- mark_recapture_data |>
    group_by(site) |>
    mutate(mean_eff_flow = mean(flow_cfs, na.rm = T),
           sd_eff_flow = sd(flow_cfs, na.rm = T),
           new_mr_flow = (flow_cfs - mean_eff_flow) / sd_eff_flow) |>
    ungroup()

  # prepare data and inits for one_site and all_sites separately
  if(model_type == "all_sites") {
    # drop sites we don't want used in efficiency estimates
    use_in_hyper <- as.integer(!sites_fit %in% exclude_from_hyper)

    data <- list("Nmr" = number_efficiency_experiments,
                 "Ntribs" = Ntribs,
                 "use_in_hyper" = use_in_hyper,
                 "ind_trib" = mark_recapture_data$ID,
                 "Releases" = mark_recapture_data$number_released,
                 "Recaptures" = mark_recapture_data$number_recaptured,
                 "effort" = effort,
                 "mr_flow" = clean_mr_flow$new_mr_flow,
                 "ind_yr" = mark_recapture_data$site_run_year_id,
                 "Nyr_re" = length(unique(mark_recapture_data$site_run_year_id)),
                 "sd_yr_ind" = sd_yr_ind) # liz updated to use dplyr to create

    # initial parameter values
    # TODO Check this
    ini_b0_pCap <- with(mark_recapture_data, {
      recap_sums    <- tapply(number_recaptured, ID, sum)
      released_sums <- tapply(number_released,  ID, sum)
      p             <- qlogis(recap_sums / released_sums)
      p[is.nan(p) | is.infinite(p)] <- -5
      p
    })

    pCap_mu_prior <- qlogis(sum(mark_recapture_data$number_recaptured) /
                              sum(mark_recapture_data$number_released))

    init_list <- list(trib_mu_P = pCap_mu_prior,
                      b0_pCap = ini_b0_pCap,
                      flow_mu_P = 0,
                      b_flow = rep(0, Ntribs),
                      trib_sd_P = 1,
                      flow_sd_P = 1,
                      pro_sd_P = rep(1, Ntribs),
                      yr_sd_P = rep(1, Ntribs))

  } else {
    # one_site
    data <- list("Nmr" = number_efficiency_experiments,
                 "Releases" = mark_recapture_data$number_released,
                 "Recaptures" = mark_recapture_data$number_recaptured,
                 "effort" = effort,
                 "mr_flow" = mark_recapture_data$standardized_efficiency_flow,
                 "ind_yr" = mark_recapture_data$site_run_year_id,
                 "Nyr_re" = length(unique(mark_recapture_data$site_run_year_id)))

    # initial parameter values
    ini_b0_pCap <- qlogis(sum(mark_recapture_data$number_recaptured) /
                            sum(mark_recapture_data$number_released))

    if(is.nan(ini_b0_pCap) | is.infinite(ini_b0_pCap)) {
      # -Inf happens when number recaptured == 0, logit of 0 is -Inf
      ini_b0_pCap <- -5
    }

    init_list <- list(b0_pCap = ini_b0_pCap,
                      b_flow = 0,
                      pro_sd_P = 1,
                      yr_sd_P = 1)
  }


  # check data list for NaNs and Infs
  invisible(lapply(names(data), function(x) {
    if (any(is.nan(data[[x]])) || any(is.infinite(data[[x]]))) {
      cli::cli_warn("NaNs or Infs detected in {.var {x}}. Please check your input data.")
      return(invisible(NULL))
    }
  }))

  invisible(lapply(names(init_list), function(x) {
    if(any(is.nan(init_list[[x]])) | any(is.infinite(init_list[[x]]))) {
      cli::cli_abort(paste0("NaNs detected in ", x, ". Please check your input data."))
    }
  }))

  inits <- list(init_list, init_list, init_list)

  # create list of inputs
  inputs <- list(data = data,
                 inits = inits,
                 parameters = pCap_parameters)

  return(list("inputs" = inputs,
              "sites_fit" = sites_fit,
              "years_fit" = years_fit,
              "sites_dropped" = exclude_from_hyper,
              "site_year_fit" = site_year_fit,
              "model_name" = model_name,
              "skew" = skew,
              "site_selection" = site_selection))

}


#' Prepare inputs for abundance model
#' @details This function prepares data for input into an abundance model.
#' @param site site for which you want to fit the model
#' @param run_year run year for which you want to fit the model
#' @param input_catch_data Optional argument for weekly catch data.
#' Defaults to `SRJPEdata::weekly_juvenile_abundance_catch_data`. If
#' passed in, structure of data frame must match that of the default.
#' @param pCap_model_type One of `all_sites`, `one_site_skew`, or `one_site`. Describes your input
#' for parameter `pCap_model_object`
#' @param min_pCap a parameter that sets the minimum pCap value for setting maximum and initial
#' values for abundance (N). Can be tweaked by min_pCap_mult. Default value is `0.0005`
#' @param min_pCap_mult a parameter that lets you tweak the min pCap value for setting
#' maximum and initial values for abundance (N). Gets multiplied by min_pCap. Default value is `1.0`
#' @param pCap_model_object the STAN model file produced by running `prepare_pCap_inputs()`
#' #' @param input_catch_data Optional argument for weekly catch data.
#' Defaults to `SRJPEdata::SRJPEdata::weekly_juvenile_abundance_catch_data`. If
#' passed in, structure of data frame must match that of the default.
#' @param input_efficiency_data Optional argument for weekly efficiency data.
#' Defaults to `SRJPEdata::SRJPEdata::weekly_juvenile_abundance_efficiency_data`. If
#' passed in, structure of data frame must match that of the default.
#' @returns a list:
#' * **abundance_inputs** a list of data and inits for input into the abundance model.
#' * **model_name** The version of the pCap logic to use for generating lt_pCap_U values. Either
#' `all_mark_recap`, `missing_mark_recap`, `no_mark_recap`, or `no_mark_recap_no_trib`.
#' * **site** The site run.
#' * **run_year** The run year.
#' * **weeks_fit** The weeks fit for the abundance model
#' * **weeks_date** Associated dates for the weeks fit for the abundance model.
#' * **lt_pCap_Us** A list containing output of `generate_lt_pCap_Us`: containing values of
#' length `Nstrata` for `lt_pCap_mu` and `lt_pCap_sd`.
#' @export
#' @md
prepare_abundance_inputs <- function(site, run_year,
                                     pCap_model_type = c("one_site", "one_site_skew", "all_sites"),
                                     min_pCap = 0.0005,
                                     min_pCap_mult = 1.0,
                                     pCap_model_object,
                                     input_catch_data = NULL,
                                     input_efficiency_data = NULL) {
  # checks and default
  if(missing(input_catch_data)) {
    input_catch_data <- SRJPEdata::weekly_juvenile_abundance_catch_data
  }
  if(missing(input_efficiency_data)) {
    input_efficiency_data <- SRJPEdata::weekly_juvenile_abundance_efficiency_data
  }

  catch_data <- input_catch_data |>
    filter(run_year == !!run_year,
           site == !!site,
           week %in% c(seq(45, 53), seq(1, 22))) |>
    group_by(year, week, stream, site, run_year) |>
    # keep NAs in count columns
    summarise(count = if(all(is.na(count))) NA_real_ else sum(count, na.rm = TRUE),
              mean_fork_length = mean(mean_fork_length, na.rm = T),
              hours_fished = mean(hours_fished, na.rm = T),
              flow_cfs = mean(flow_cfs, na.rm = T),
              average_hours_fished_during_efficiency_trials = mean(average_hours_fished_during_efficiency_trials, na.rm = T),
              standardized_flow = mean(standardized_flow, na.rm = T)
              #lgN_prior = mean(lgN_prior, na.rm = T)
              ) |>
    ungroup() |>
    left_join(site_order_north_south, by = "site") |>
    arrange(ns_order) |>
    select(-ns_order) |>
    mutate(count = round(count, 0),
           # change all NaNs to NAs
           # across(mean_fork_length:lgN_prior, ~ifelse(is.nan(.x), NA, .x)),
           # calculate effort
           effort = hours_fished / average_hours_fished_during_efficiency_trials)

  # data checks
  no_data <- nrow(catch_data) == 0
  all_na <- !no_data && all(is.na(catch_data$count))

  if (no_data || all_na) {
    reason <- if (no_data) "no catch data" else "all count values are NA"
    cli::cli_warn(c(
      "Skipping site {.val {site}}, run year {.val {run_year}}:",
      "i" = "{reason}"
    ))
    return(invisible(NULL))
  }

  # Calculate lincoln peterson abundance
  lp_data <- catch_data |>
    left_join(input_efficiency_data |>
                select(-flow_cfs),
              by = c("year", "run_year", "week", "stream", "site")) |>
    mutate(# plot things
           lincoln_peterson_abundance = count * (number_released / number_recaptured),
           lincoln_peterson_efficiency = number_recaptured / number_released,
           sampled = ifelse(is.na(count), FALSE, TRUE),
           efficiency_trial = ifelse(is.na(lincoln_peterson_efficiency), FALSE, TRUE)) |>
    left_join(SRJPEmodel::julian_week_to_date_lookup, by = c("week" = "Jwk")) |>
    left_join(site_order_north_south, by = "site") |>
    arrange(ns_order) |>
    select(-ns_order) |>
    mutate(year = ifelse(week >= 43, run_year - 1, run_year),
           date = factor(date, levels = date),
           week_index = row_number()) |>
    select(week, count, lincoln_peterson_abundance:date, number_released, number_recaptured)

  # min pCap calculations
  min_pCap_new <- input_catch_data |>
    filter(site == !!site) |>
    left_join(input_efficiency_data |>
                select(-flow_cfs),
              by = c("year", "run_year", "week", "stream", "site")) |>
    mutate(lincoln_peterson_efficiency = number_recaptured / number_released) |>
    filter(!is.na(lincoln_peterson_efficiency),
           lincoln_peterson_efficiency > 0)

  if(nrow(min_pCap_new) == 0) {
    min_pCap_new <- min_pCap
  } else {
    min_pCap_new <- min_pCap_new |>
      pull(lincoln_peterson_efficiency) |>
      min(na.rm = T)
  }

  min_pCap_new = min_pCap_new * min_pCap_mult #reduce by a 0-1 factor min_pCap_mult

  # indexing values for BUGS code
  number_weeks_catch <- nrow(catch_data) # for looping through the catch dataset
  indices_with_catch <- which(!is.na(catch_data$count)) # indices of weeks with catch data
  number_weeks_with_catch <- length(indices_with_catch) # how many weeks actually have catch

  # analyze efficiency trials for all relevant sites (do not filter to site)
  # set up filter - if it's a tributary-based model, we cannot use efficiencies from KDL, TIS, RBDD
  if(!site %in% c("knights landing", "tisdale", "red bluff diversion dam")) {
    # drop sites from arguments and also remove mainstem
    remove_sites <- c("knights landing", "tisdale", "red bluff diversion dam")

    # prepare "mark recapture" dataset - all mark-recap trials in the system
    mark_recapture_data <- input_efficiency_data |>
      # grab standardized_flow
      left_join(input_catch_data |>
                  select(year, week, stream, site, run_year, standardized_flow),
                by = c("year", "week", "run_year", "stream", "site")) |>
      # or do we want to filter just no number released?
      dplyr::filter(!site %in% remove_sites &
                      !is.na(standardized_efficiency_flow),
                    !is.na(number_released) &
                      !is.na(number_recaptured)) |>
      # right now there's lifestage in the dataset, so we have to do dplyr::distinct()
      dplyr::distinct(site, run_year, week, number_released, number_recaptured, .keep_all = TRUE) |>
      group_by(site) |>
      left_join(site_order_north_south, by = "site") |>
      arrange(ns_order, year, week) |>
      ungroup()

  } else {
    # prepare "mark recapture" dataset but filter to only the mainstem site
    mark_recapture_data <- input_efficiency_data |>
      # grab standardized_flow
      left_join(input_catch_data |>
                  select(year, week, stream, site, run_year, standardized_flow),
                by = c("year", "week", "run_year", "stream", "site")) |>
      # or do we want to filter just no number released?
      dplyr::filter(site == !!site &
                      !is.na(standardized_efficiency_flow),
                    !is.na(number_released) &
                      !is.na(number_recaptured)) |>
      # right now there's lifestage in the dataset, so we have to do dplyr::distinct()
      dplyr::distinct(site, run_year, week, number_released, number_recaptured, .keep_all = TRUE) |>
      group_by(site) |>
      left_join(site_order_north_south, by = "site") |>
      arrange(ns_order, year, week) |>
      ungroup()
  }


  # bring together efficiency and catch data so that we can get the indices of
  # catch data (hence left join) that correspond to certain efficiency trial
  # information.
  all_data_for_indexing <- left_join(catch_data, mark_recapture_data,
                                     by = c("year", "week", "stream",
                                            "site", "run_year", "flow_cfs",
                                            "standardized_flow")) |>
    select(-ns_order) |>
    left_join(site_order_north_south, by = "site") |>
    arrange(ns_order) |>
    select(-ns_order)

  # get use_trib, sites_fit, and ind_trib indexing
  # first assign 1:Ntribs to the unique sites in the dataset
  site_lookup <- mark_recapture_data |>
    dplyr::distinct(site) |>
    mutate(ID = row_number())

  # pull the order of those sites
  sites_fit <- site_lookup |>
    dplyr::pull(site)

  # assign the IDs to the sites in the mark-recapture dataset
  indices_site_mark_recapture <- mark_recapture_data |>
    left_join(site_lookup, by = "site") |>
    dplyr::pull(ID)

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

  # plotting vectors for josh
  efficiency_plotting_vectors <- all_data_for_indexing |>
    arrange(week) |>
    #arrange(year, week) |>
    filter(week %in% weeks_with_mark_recapture) |>
    select(week, number_released, number_recaptured)

  # set up b-spline basis matrix
  # this corresponds to line 148-153 in josh_original_model_code.R
  if (number_weeks_catch < 4) {
    cli::cli_warn(c(
      "Skipping site {.val {site}}, run year {.val {run_year}}:",
      "i" = "Fewer than 4 weeks of catch data available ({number_weeks_catch} found).",
      "i" = "Spline parameters require at least 4 data points."
    ))
    return(invisible(NULL))
  }

  # spline parameter calculation
  spline_data <- SRJPEmodel::build_spline_data(number_weeks_catch, k_int = 4) # rule of thumb is 1 knot for every 4 data points for a cubic spline (which has 4 parameters)

  # should be same dimensions as catch_flow
  effort <- catch_data$effort  # TODO check should be all Nstrata, and if NA, that's okay because it won't be used but set to 0 for BUGS
  effort[effort == 0] <- 0.001 # catch for 0 values in effort
  effort[is.na(effort)] <- 0.001
  # pass in catch
  weekly_catch_data <- catch_data$count[indices_with_catch]

  # Set prior for log N and ini values based on mean catch across wks trapped fish for all weeks
  # This will provide values for weeks trap wasn't fith
  lgN_max = rep(log(0.001 * (mean(weekly_catch_data) + 1) / min_pCap_new), number_weeks_catch)
  ini_lgN = rep(log(0.001 * (min(weekly_catch_data) + 1) / (min_pCap_new * 2)), number_weeks_catch)

  # Then overide these values for weeks trap fish based on the weekly catch
  for(j in 1:number_weeks_with_catch){
    index_with_catch <- indices_with_catch[j]
    lgN_max[index_with_catch] = log(0.001 * (weekly_catch_data[j] + 1) / min_pCap_new)
    ini_lgN[index_with_catch] = log(0.001 * (weekly_catch_data[j] + 1) / (min_pCap_new * 2))
  }


  # build data list
  data <- list("Nstrata" = number_weeks_catch,
               "Nstrata_wc" = number_weeks_with_catch,
               "u" = weekly_catch_data,
               "K" = spline_data$K,
               "ZP" = spline_data$b_spline_matrix,
               "Uwc_ind" = indices_with_catch,
               "lgN_max" = lgN_max,
               "effort" = effort)

  # data needed for generating lt_pCap_Us
  # also plotting data needs (sorted) for Josh's code
  lt_pCap_U_data <- list("catch_flow" = catch_data$standardized_flow,
                         "use_trib" = indices_sites_pCap,
                         "Nwmr" = number_weeks_with_mark_recapture,
                         "Nwomr" = number_weeks_without_mark_recapture,
                         "Uind_wMR" = indices_with_mark_recapture,
                         "Uind_woMR" = indices_without_mark_recapture,
                         "releases_sort" = efficiency_plotting_vectors$number_released, # this is for josh's plots
                         "recaptures_sort" = efficiency_plotting_vectors$number_recaptured, # for josh's plots
                         "ind_pCap" = indices_pCap)

  # use number of experiments at site to determine which model to call
  number_experiments_at_site <- mark_recapture_data |>
    dplyr::distinct(site, run_year, week) |>
    filter(site == !!site) |>
    nrow()

  if (number_experiments_at_site == 1) {
    cli::cli_abort(c(
      "No model is available for sites with only one experiment.",
      "i" = "Site {.val {site}} has {number_experiments_at_site} experiment recorded.",
      "i" = "Multi-experiment inter-trial variance handling is not yet implemented."
    ))
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
  invisible(lapply(names(data), function(x) {
    if (any(is.nan(data[[x]])) || any(is.infinite(data[[x]]))) {
      cli::cli_abort(c(
        "Invalid values detected in {.var {x}}.",
        "i" = "NaNs or infinite values cannot be used as initial values.",
        "i" = "Please check your input data."
      ))
    }
  }))

  parameters <- c("tau_N", "tau_Ne", "b_sp", "lg_N", "lt_pCap_U",
                  "N", "Ntot", "lg_CumN")

  # inits
  init_list <- list(b_sp = rep(1, spline_data$K),lg_N = ini_lgN)

  invisible(lapply(names(init_list), function(x) {
    if (any(is.nan(init_list[[x]])) || any(is.infinite(init_list[[x]]))) {
      cli::cli_abort(c(
        "Invalid initial values detected in {.var {x}}.",
        "i" = "NaNs or infinite values cannot be used as initial values.",
        "i" = "Please check your input data."
      ))
    }
  }))

  inits <- list(init_list, init_list, init_list)

  inputs_for_abundance <- list(data = data,
                               inits = inits,
                               parameters = parameters)

  weeks_fit <- tibble("Jwk" = catch_data$week) |>
    left_join(SRJPEmodel::julian_week_to_date_lookup, by = "Jwk")

  # run_year_id to index for yr_re
  # and year_sd_id to index for yr_sd_P
  run_year_id_lookup <- mark_recapture_data |>
    arrange(ns_order, run_year) |>
    distinct(site, run_year) |>
    mutate(site_run_year_id = as.integer(factor(paste(site, run_year),
                                                levels = unique(paste(site, run_year)))),
           year_sd_id = as.integer(factor(site, levels = unique(site))))

  run_year_id <- run_year_id_lookup$site_run_year_id[which(run_year_id_lookup$site == site & run_year_id_lookup$run_year == run_year)]
  year_sd_id <- unique(run_year_id_lookup$year_sd_id)[which(unique(run_year_id_lookup$site) == site)]

  abundance_inputs <- list("inputs" = inputs_for_abundance,
                           "lt_pCap_U_data" = lt_pCap_U_data,
                           "model_name" = model_name,
                           "pCap_model_type" = pCap_model_type,
                           "site" = site,
                           "run_year" = run_year,
                           "catch_flow_raw" = catch_data$flow_cfs,
                           "mr_flow_raw" = mark_recapture_data$flow_cfs,
                           "weeks_fit" = weeks_fit$Jwk,
                           "week_date" = weeks_fit$date,
                           "sites_fit" = sites_fit,
                           "run_year_id" = run_year_id,
                           "year_sd_id" = year_sd_id,
                           "min_pCap" = min_pCap_new)

  # generate lt pcap Us based on inputs
  lt_pCap_Us <- generate_lt_pCap_Us(abundance_inputs, pCap_model_object)

  final_abundance_inputs <- modifyList(abundance_inputs,
                                       list("lt_pCap_Us" = lt_pCap_Us,
                                            "lp_data" = lp_data))

  return(final_abundance_inputs)

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
#' @family Fit model
#' @export
#' @md
fit_pCap_model <- function(input) {

  # call the correct model based on pCap_inputs model name (either one_site_skew, one_site, or all_sites)
  stan_model <- eval(parse(text = paste0("SRJPEmodel::bt_spas_x_model_code$", input$model_name)))

  # call model
  options(mc.cores=parallel::detectCores())

  pcap <- rstan::stan(model_name = input$model_name,
                      model_code = stan_model,
                      data = input$inputs$data,
                      init = input$inputs$inits,
                      pars = input$inputs$parameters,
                      chains = SRJPEmodel::bt_spas_x_bayes_params$number_chains,
                      iter = 10000,
                      seed = 84735,
                      control = list(max_treedepth = 15))

  return(pcap)

}

#' Generate lt_pCap_U values from the pCap model object
#' @details This function prepares the data list for the abundance STAN model based on what the model name is.
#' @param abundance_inputs A list containing the inputs for the site-specific abundance model including
#' `model_name`, `Nstrata`, `catch_flow`, `use_trib`, `Ntribs`, `Nmr`, `Nwmr`, `Nwomr`, `Uind_wMR`, `Uind_woMR`,
#' `ind_pCap`. This is created by calling `prepare_abundance_inputs()`.
#' @returns a named list with the simulated `lt_pCap_mu` and `lt_pCap_sd` values for
#' `Nstrata`.
#' @export
#' @md
generate_lt_pCap_Us <- function(abundance_inputs, pCap_model_object){

  pCap_model_type <- ifelse(abundance_inputs$pCap_model_type == "all_sites", "all_sites", "one_site")
  skew <- ifelse(pCap_model_type == "one_site" & abundance_inputs$pCap_model_type == "one_site_skew", TRUE, FALSE)

  # if(any(abundance_inputs$sites_fit %in% c("knights landing", "tisdale", "red bluff diversion dam"))) {
  if(pCap_model_type == "one_site") {
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
    pars_to_extract <- c("logit_pCap", "b0_pCap", "b_flow", "pro_sd_P", "yr_re", "yr_sd_P") # , "b_eff")
    if(skew) {
      pars_to_extract <- c(pars_to_extract, "alpha")
    }

    samples <- rstan::extract(pCap_model_object, pars = pars_to_extract,
                              permuted = TRUE)

    Ntrials <- dim(samples$logit_pCap)[1] # of saved posterior samples from pCap model in stan
    logit_pCap <- samples$logit_pCap # logit_pCap[1:Ntrials,1:Nmr] # The estimated logit pCap posterior for each efficiency trial
    b0_pCap <- samples$b0_pCap # b0_pCap[1:Ntrials,1] #mean logit pCap for the mainstem site
    b_flow <- samples$b_flow # b_flow[1:Ntrials,1] #flow effect for the mainstem site
    # b_eff <- samples$b_eff
    pro_sd_P <- samples$pro_sd_P # pro_sd_P[1:Ntrials]        #process error (sd)
    alpha <- samples$alpha #will be null if Skew==F
    run_year <- abundance_inputs$run_year

    #run_year_id =lookup run_year in mr_year_lookup table for KL or Tis and return run_year_id
    yr_re <- samples$yr_re[, abundance_inputs$run_year_id]#need to know the year abundance model running for and determine which element of yr_re it represents
    yr_sd_P <- samples$yr_sd_P

    # calculations
    lt_pCap_U = matrix(nrow = Ntrials, ncol = Nstrata)
    lt_pCap_mu = matrix(nrow = Nstrata, ncol = Ntrials) #function needs to return this
    lt_pCap_sd = lt_pCap_mu                        #function needs to return this

    if(ModelName == "all_mark_recap" ){#stays as is

      for(i in 1:Nstrata){
        lt_pCap_U[,i] = logit_pCap[,ind_pCap[i]];
      }

    } else if (ModelName == "missing_mark_recap"){

      for(i in 1:Nwmr){
        #Assign estimated pCaps for strata with efficiency data
        lt_pCap_U[,Uind_wMR[i]] = logit_pCap[,ind_pCap[i]];
      }
      for (i in 1:Nwomr) {
        #for weeks without efficiency trials
        if(skew){
          sim_pro_dev=brms::rskew_normal(n=Ntrials, mu = 0, sigma = pro_sd_P, alpha = alpha, xi = NULL, omega = NULL)
        } else {
          sim_pro_dev=rnorm(n=Ntrials, mean=0,sd=pro_sd_P)
        }
        lt_pCap_U[,Uind_woMR[i]] = b0_pCap + b_flow * catch_flow[Uind_woMR[i]] + sim_pro_dev + yr_re + log(abundance_inputs$inputs$data$effort[Uind_woMR[i]])
      }

    } else if (ModelName == "no_mark_recap"){

      for (i in 1:Nwomr) {
          if(skew){
            sim_pro_dev=rskew_normal(n=Ntrials, mu = 0, sigma = pro_sd_P, alpha = alpha, xi = NULL, omega = NULL)
          } else {
            sim_prod_dev=rnorm(n=Ntrials, mean=0,sd=pro_sd_P)
          }
          sim_yr_dev = rnorm(n=Ntrials, mean=0, sd=yr_sd_P)
          lt_pCap_U[,Uind_woMR[i]] = b0_pCap + b_flow * catch_flow[Uind_woMR[i]] + sim_pro_dev + sim_yr_dev + log(abundance_inputs$inputs$data$effort[Uind_woMR[i]])
      }

    } else if (ModelName == "no_mark_recap_no_trib"){

       logit_b0 = rnorm(n=Ntrials, mean=trib_mu_P, sd=trib_sd_P)
       logit_bflow = rnorm(n=Ntrials, mean=flow_mu_P, sd=flow_sd_P)


      for (i in 1:Nwomr) {
        if(skew){
          sim_pro_dev=rskew_normal(n=Ntrials, mu = 0, sigma = pro_sd_P, alpha = alpha, xi = NULL, omega = NULL)
        } else {
          sim_prod_dev=rnorm(n=Ntrials, mean=0,sd=pro_sd_P)
        }
        sim_yr_dev = rnorm(n=Ntrials, mean=0, sd=yr_sd_P)#this won't work since we don't have yr_sd_P for a site with no mr data
        lt_pCap_U[,Uind_woMR[i]] = b0_pCap + b_flow * catch_flow[Uind_woMR[i]] + sim_pro_dev + sim_yr_dev + log(abundance_inputs$inputs$data$effort[Uind_woMR[i]])
      }

    }#end if on ModelName


    #Calculate mean and sd for each lt_pCap_U and return these from function
    for(i in 1:Nstrata){
      lt_pCap_mu[i,]=mean(lt_pCap_U[,i])
      lt_pCap_sd[i,]=sd(lt_pCap_U[,i])
    }

  } else { # if all_sites
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
                                                          "flow_mu_P", "flow_sd_P","yr_re","yr_sd_P"), #, "b_eff"),
                              permuted = TRUE)
    Ntrials <- dim(samples$logit_pCap)[1] # of saved posterior samples from pCap model in stan
    logit_pCap <- samples$logit_pCap # logit_pCap[1:Ntrials,1:Nmr] # The estimated logit pCap posterior for each efficiency trial
    b0_pCap <- samples$b0_pCap # b0_pCap[1:Ntrials,1:Ntribs] #mean logit pCap for each site (at mean discharge)
    b_flow <- samples$b_flow # b_flow[1:Ntrials,1:Ntribs] #flow effect for each site
    # b_eff <- samples$b_eff # b_eff
    pro_sd_P <- samples$pro_sd_P #process error (sd)
    trib_mu_P <- samples$trib_mu_P # trib_mu_P[1:Ntrials] #hyper mean for b0_pCap
    trib_sd_P <- samples$trib_sd_P # trib_sd_P[1:Ntrials] #hyper sd for b0_pCap
    flow_mu_P <- samples$flow_mu_P # flow_mu_P[1:Ntrials] #hyper mean for b_flow
    flow_sd_P <- samples$flow_sd_P # flow_sd_P[1:Ntrials] #hyper sd for b_flow

    site = abundance_inputs$site
    run_year = abundance_inputs$run_year

    #run_year_id =lookup run_year in mr_year_lookup table given site and run year above return run_year_id
    yr_re <- samples$yr_re[, abundance_inputs$run_year_id]#need to know the year abundance model running for and determine which element of yr_re it represents

    #get sd_yr ind from mr_year_lookup based on site and run_year
    yr_sd_P <- samples$yr_sd_P

    # calculations
    lt_pCap_U=matrix(nrow=Ntrials,ncol=Nstrata)
    sim_pro_dev=vector(length=Ntrials)
    lt_pCap_mu=matrix(nrow=Nstrata,ncol=Ntrials) #function needs to return this
    lt_pCap_sd=lt_pCap_mu                        #function needs to return this
    sim_yr_dev=vector(length=Ntrials)

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
        # need to account for effort for weeks without mark recap trials
        sim_pro_dev = rnorm(n=Ntrials, mean=0,sd=pro_sd_P[,use_trib]);
        lt_pCap_U[,Uind_woMR[i]] = b0_pCap[,use_trib] + b_flow[,use_trib] * catch_flow[Uind_woMR[i]] + yr_re + sim_pro_dev + log(abundance_inputs$inputs$data$effort[Uind_woMR[i]])
        # b_eff[,use_trib] * abundance_inputs$inputs$data$effort[Uind_woMR[i]] TODO double check indexing
      }

    } else if (ModelName=="no_mark_recap"){

      for (i in 1:Nwomr) {
        sim_pro_dev = rnorm(n=Ntrials, mean=0, sd=pro_sd_P[,use_trib])
        sim_yr_dev = rnorm(n=Ntrials, mean=0, sd=yr_sd_P[,use_trib])

        lt_pCap_U[,Uind_woMR[i]] = b0_pCap[,use_trib] + b_flow[,use_trib] * catch_flow[Uind_woMR[i]] + sim_yr_dev + sim_pro_dev + log(abundance_inputs$inputs$data$effort[Uind_woMR[i]])
      }

    } else if (ModelName=="no_mark_recap_no_trib"){

      logit_b0 = rnorm(n=Ntrials, mean=trib_mu_P, sd=trib_sd_P)
      logit_bflow = rnorm(n=Ntrials, mean=flow_mu_P, sd=flow_sd_P)


      for (i in 1:Nwomr) {
        sim_pro_dev = rnorm(n=Ntrials, mean=0, sd=rowMeans(pro_sd_P))
        sim_yr_dev = rnorm(n=Ntrials, mean=0, sd=rowMeans(yr_sd_P))
        lt_pCap_U[,Uind_woMR[i]] = logit_b0 + logit_bflow * catch_flow[Uind_woMR[i]] + sim_yr_dev + sim_pro_dev + log(abundance_inputs$inputs$data$effort[Uind_woMR[i]])
      }

    }#end if on ModelName


    #Calculate mean and sd for each lt_pCap_U and return these from function
    for(i in 1:Nstrata){
      lt_pCap_mu[i,]=mean(lt_pCap_U[,i])
      lt_pCap_sd[i,]=sd(lt_pCap_U[,i])
    }
  }

  return(list("lt_pCap_mu" = lt_pCap_mu |> rowMeans(),
              "lt_pCap_sd" = lt_pCap_sd |> rowMeans()))
}


#' Fit abundance model in BUGS.
#' @details This function runs the BUGS abundance model.
#' @param abundance_inputs the object produced by `prepare_abundance_inputs()`
#' @param bugs_directory the filepath pointing to where your WinBUGS14/ directory is
#' @returns a BUGS object.
#' @export
#' @md
fit_abundance_model_BUGS <- function(abundance_inputs,
                                     bugs_directory) {
  parameters <- c("lt_pCap_U", "pCap_U", "Usp", "N", "Ntot", "sd.N", "sd.Ne", "lg_CumN")
  Nmcmc = 2000
  Nburnin = 500
  Nthin = 2
  Nchains = 3

  data <- abundance_inputs$inputs$data
  data$lt_pCap_mu <- abundance_inputs$lt_pCap_Us$lt_pCap_mu
  data$lt_pCap_tau <- 1/abundance_inputs$lt_pCap_Us$lt_pCap_sd^2
  inits_with_lt_pCap_U <- abundance_inputs$inputs$inits[[1]]
  inits_with_lt_pCap_U$lt_pCap_U <- data$lt_pCap_mu
  inits_with_lt_pCap_U$tau.N <- 1
  inits_with_lt_pCap_U$tau.Ne <- 1
  new_inits <- list(inits_with_lt_pCap_U,
                    inits_with_lt_pCap_U,
                    inits_with_lt_pCap_U)

  data$effort <- NULL # not used in abundance.bug any more, can't be in data list if not called

  # Write model code from package object to a temp file
  model_file <- tempfile(fileext = ".bug")
  writeLines(SRJPEmodel::bt_spas_x_model_code$abundance_BUGS, con = model_file)
  on.exit(unlink(model_file))

  abundance <- bugs(data,
                    new_inits,
                    parameters,
                    model_file,
                    debug = F,
                    n.chains = Nchains, n.burnin = Nburnin, n.thin = Nthin, n.iter = Nmcmc,
                    codaPkg = F, DIC = T, clearWD = T,
                    bugs.directory = bugs_directory)

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
extract_abundance_estimates <- function(abundance_inputs,
                                        model_object) {

  # link to actual weeks
  # TODO do we want week formatted as MONTH-DATE ? can do easily
  week_lookup <- tibble("week_fit" = abundance_inputs$weeks_fit) |>
    mutate(week_index = row_number())

  stream_lookup <- SRJPEdata::site_lookup |>
    dplyr::distinct(stream, site)

  formatted_table <- model_object$summary |>
    as.data.frame() |>
    tibble::rownames_to_column("parameter") |>
    mutate(week_index = ifelse(str_detect(parameter, "b0_pCap|b_flow"), NA,
                               suppressWarnings(readr::parse_number(parameter))),
           site = abundance_inputs$site,
           run_year = abundance_inputs$run_year,
           model_name = abundance_inputs$model_name,
           parameter = gsub("[0-9]+|\\[|\\]", "", parameter),
           srjpedata_version = as.character(packageVersion("SRJPEdata"))) |>
    left_join(week_lookup, by = "week_index") |>
    left_join(stream_lookup, by = "site") |>
    # now clean up statistics
    pivot_longer(mean:Rhat,
                 values_to = "value",
                 names_to = "statistic") |>
    mutate(statistic = str_remove_all(statistic, "\\%")) |>
    select(model_name, site, stream, run_year, week_fit, parameter, statistic, value, srjpedata_version)

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

  site_lookup <- tibble("site" = pCap_inputs$sites_fit) |>
    mutate(site_index = row_number())

  stream_lookup <- SRJPEdata::site_lookup |>
    dplyr::distinct(stream, site)

  formatted_table <- rstan::summary(model_object)$summary |>
    as.data.frame() |>
    tibble::rownames_to_column("parameter") |>
    mutate(week_index = NA, # no weekly estimates that are relevant
           site_index = ifelse(str_detect(parameter, "b0_pCap|b_flow"),
                               suppressWarnings(readr::parse_number(substr(parameter, 3, length(parameter)))),
                               NA),
           run_year = NA,
           model_name = model_object@model_name,
           parameter = gsub("[0-9]+|\\[|\\]", "", parameter),
           srjpedata_version = as.character(packageVersion("SRJPEdata"))) |>
    left_join(site_lookup, by = "site_index") |>
    left_join(stream_lookup, by = "site") |>
    mutate(stream = ifelse(is.na(site), NA, stream)) |>
    # now clean up statistics
    pivot_longer(mean:Rhat,
                 values_to = "value",
                 names_to = "statistic") |>
    mutate(statistic = str_remove_all(statistic, "\\%"),
           week_fit = NA) |> # this won't be reported for the pCap model
    select(model_name, site, stream, run_year, week_fit, parameter, statistic, value, srjpedata_version)

  if(pCap_inputs$model_name == "pCap_mainstem_skew_re") {
    formatted_table <- formatted_table |>
      select(-c(site, stream)) |>
      mutate(site = pCap_inputs$site) |>
      left_join(stream_lookup, by = "site")
  }

  return(formatted_table)
}

#' Run all JPE sites.
#' @details This function automates running BT-SPAS-X for all JPE sites/run years.
#' @param sites_to_run
#' @param run_pCap
#' @param mainstem
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
                                    mainstem,
                                    pCap_model_object_filepath,
                                    bugs_model_file,
                                    bugs_directory) {

  # run pCap model if necessary
  if(run_pCap) {
    pCap_inputs <- prepare_pCap_inputs(mainstem = mainstem)
    pCap <- fit_pCap_model(pCap_inputs$inputs)
    saveRDS(pCap, pCap_model_object_filepath)
  }

  # prep inputs as vectors
  sites_to_run_inputs <- sites_to_run

  # now run abundance workflow
  SRJPE_fits_table <- purrr::pmap(list(sites_to_run_inputs$site,
                                       sites_to_run_inputs$run_year,
                                       pCap_model_object_filepath,
                                       bugs_model_file,
                                       bugs_directory),
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
                                   pCap_model_object_filepath,
                                   bugs_model_file,
                                   bugs_directory) {

  cli::cli_bullets(paste0("Running abundance model for ", site, " for run year ", run_year))

  results <- tryCatch({

    abundance_inputs <- prepare_abundance_inputs(site, run_year, effort_adjust = T)

    pCap <- readRDS(pCap_model_object_filepath)

    lt_pCap_Us <- generate_lt_pCap_Us(abundance_inputs, pCap)

    abundance <- fit_abundance_model_BUGS(abundance_inputs,
                                          bugs_model_file,
                                          bugs_directory)
    clean_table <- extract_abundance_estimates(abundance_inputs, abundance)
    return(clean_table)

  },
  error = function(e) return(tibble("site" = site,
                                    "run_year" = run_year,
                                    "error" = TRUE)))

  return(results)
}



#' BT SPAS X diagnostic plots
#' @details This function produces a plot with data and results of fitting the pCap and abundance models for
#' a given site and run year.
#' @param inputs a list of inputs generated by running `prepare_abundance_inputs()`
#' @param results_df a table of parameter estimates generated by running `extract_abundance_estimates()` or
#' by pulling and filtering from the database using `get_most_recent_model_results()`.
#' @returns A plot.
#' @export
#' @md
generate_diagnostic_plot_juv <- function(inputs, results_df) {

  params <- results_df |>
    filter(model_name == "bt_spas_x",
           site == inputs$site,
           year == inputs$run_year) |>
    select(-c(id, model_run_id, model_name)) |>
    pivot_wider(names_from = statistic,
                values_from = value)

  pCap_estimates <- params |>
    filter(parameter == "lt_pCap_U") |>
    select(week = week_fit,
           c("97.5", "50", "mean", "75", "sd", "25", "2.5"))

  N_estimates <- params |>
    filter(parameter == "N") |>
    select(week = week_fit,
           c("97.5", "50", "mean", "75", "sd", "25", "2.5"))

  abundance_plot <- N_estimates |>
    left_join(inputs$lp_data,
              by = "week") |>
    mutate(count_label = ifelse(is.na(count), "", count)) |>
    ggplot(aes(x = date, y = `50`)) +
    geom_bar(stat = "identity", fill = "grey", width = .75) +
    geom_errorbar(aes(x = date, ymin = `2.5`, ymax = `97.5`), width = 0.2) +
    geom_point(aes(x = date, y = lincoln_peterson_abundance),
               shape = 1, color = "blue" ,size = 3) +
    geom_point(aes(x = date, y = Inf, color = sampled),
               size = 3) +
    geom_text(aes(x = date, y = Inf,
                  label = paste(count_label),
                  angle = 90),
              hjust = 1,
              size = 3) +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "#ec5858")) +
    theme_minimal() +
    labs(x = "",
         #x = "Date",
         y = "Abundance",
         title = paste(inputs$site, inputs$run_year)) +
    #theme(axis.text.x=element_blank()) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "")

  # efficiency
  efficiency_plot <- pCap_estimates |>
    left_join(inputs$lp_data,
              by = "week") |>
    mutate(across(c(mean, `50`, `2.5`, `97.5`), plogis),
           number_released_label = ifelse(is.na(number_released), "", number_released),
           number_recaptured_label = ifelse(is.na(number_recaptured), "", number_recaptured)) |>
    ggplot(aes(x = date, y = `50`, fill = efficiency_trial)) +
    geom_bar(stat = "identity", width = .75) +
    geom_errorbar(aes(x = date, ymin = `2.5`, ymax = `97.5`), width = 0.2) +
    geom_point(aes(x = date, y = lincoln_peterson_efficiency),
               shape = 1, color = "blue", size = 3) +
    geom_text(aes(x = date, y = Inf,
                  label = paste(number_released_label, number_recaptured_label),
                  angle = 90),
              hjust = 1,
              size = 3) +
    scale_fill_manual(values = c("TRUE" = "grey", "FALSE" = "#ec5858")) +
    theme_minimal() +
    labs(x = "Date", y = "Weekly Efficiency") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "")

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
    filter(run_year == !!run_year,
           site == !!site,
           week %in% c(seq(45, 53), seq(1, 22))) |>
    group_by(year, week, stream, site, run_year) |>
    # keep NAs in count columns
    summarise(count = if(all(is.na(count))) NA_real_ else sum(count, na.rm = TRUE),
              mean_fork_length = mean(mean_fork_length, na.rm = T),
              hours_fished = mean(hours_fished, na.rm = T),
              catch_flow_cfs = mean(flow_cfs, na.rm = T),
              average_hours_fished_during_efficiency_trials = mean(average_hours_fished_during_efficiency_trials, na.rm = T),
              standardized_flow = mean(standardized_flow, na.rm = T)
              #lgN_prior = mean(lgN_prior, na.rm = T)
              ) |>
    ungroup() |>
    left_join(SRJPEdata::weekly_juvenile_abundance_efficiency_data,
              by = c("year", "run_year", "week", "stream", "site")) |>
    mutate(count = round(count, 0),
           # change all NaNs to NAs
           # across(mean_fork_length:lgN_prior, ~ifelse(is.nan(.x), NA, .x)),
           # plot things
           lincoln_peterson_abundance = count * (number_released / number_recaptured),
           lincoln_peterson_efficiency = number_recaptured / number_released) |>
    left_join(SRJPEmodel::julian_week_to_date_lookup, by = c("week" = "Jwk")) |>
    left_join(site_order_north_south, by = "site") |>
    arrange(ns_order) |>
    select(-ns_order) |>
    mutate(date = factor(date, levels = date),
           week_index = row_number())

  abundance_plot <- data |>
    ggplot(aes(x = date, y = count)) +
    geom_bar(stat = "identity", fill = "grey", width = .75) +
    geom_text(aes(x = date, y = Inf,
                  label = paste(count),
                  angle = 90),
              hjust = 1,
              size = 3) +
    theme_minimal() +
    labs(x = "",
         y = "Abundance",
         title = paste(site, run_year)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  efficiency_plot <- data |>
    mutate(number_released_label = ifelse(is.na(number_released), "", number_released),
           number_recaptured_label = ifelse(is.na(number_recaptured), "", number_recaptured)) |>
    ggplot(aes(x = date, y = lincoln_peterson_efficiency)) +
    geom_point(shape = 1, color = "blue") +
    # geom_bar(stat = "identity", fill = "grey", width = 4) +
    geom_text(aes(x = date, y = Inf,
                  label = paste(number_released_label, number_recaptured_label),
                  angle = 90),
              hjust = 1,
              size = 3) +
    theme_minimal() +
    labs(x = "Date", y = "Weekly Efficiency") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  # arrange together
  gridExtra::grid.arrange(abundance_plot, efficiency_plot)

}

