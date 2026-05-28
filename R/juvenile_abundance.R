#' Prepare inputs for the pCap Stan model
#'
#' @details Prepares the data list, initial values, and indexing structures
#'   needed to fit the capture probability (pCap) Stan model. Supports two
#'   model types: a hierarchical model fit across all tributary sites
#'   (`"all_sites"`) and a single-site model for mainstem sites
#'   (`"one_site"`), optionally with a skew-normal process-error distribution.
#'
#' @param model_type Character. Either `"all_sites"` or `"one_site"`.
#'   The `"all_sites"` model fits a hierarchical model across all available
#'   tributary sites. The `"one_site"` model fits a single mainstem site
#'   specified by `site_selection`.
#' @param skew Logical (`TRUE`/`FALSE`) or `NULL`. Whether to use a
#'   skew-normal distribution for process error. Only applies when
#'   `model_type = "one_site"`; ignored (with a warning) for `"all_sites"`.
#'   Must be specified explicitly for `"one_site"` models.
#' @param site_selection Character. The site name to fit when
#'   `model_type = "one_site"`. Must be exactly one value. Ignored for
#'   `"all_sites"`.
#' @param exclude_from_hyper Character vector of site names to exclude from
#'   the hyper-distribution over `b0_pCap` and `b_flow`. Default `NULL`
#'   (no sites excluded). Values must be present in
#'   `unique(SRJPEdata::weekly_juvenile_abundance_efficiency_data$site)`.
#' @param input_catch_data Optional data frame of weekly catch data.
#'   Defaults to `SRJPEdata::weekly_juvenile_abundance_catch_data`.
#'   If supplied, the structure must match that of the default.
#' @param input_efficiency_data Optional data frame of weekly efficiency
#'   (mark-recapture) data. Defaults to
#'   `SRJPEdata::weekly_juvenile_abundance_efficiency_data`.
#'   If supplied, the structure must match that of the default.
#'
#' @returns A named list with the following elements:
#' * **inputs** A list with three sub-elements passed directly to
#'   `fit_pCap_model()`: `data` (Stan data list), `inits` (initial values
#'   replicated across 3 chains), and `parameters` (character vector of
#'   parameters to monitor).
#' * **sites_fit** Character vector of site names in the order they are
#'   indexed by `ind_trib` in the Stan data.
#' * **years_fit** Numeric vector of run years in the order they are indexed
#'   by `ind_yr`.
#' * **sites_dropped** The value passed to `exclude_from_hyper`; `NULL` if
#'   none were excluded.
#' * **site_year_fit** Data frame with columns `site` and `run_year` giving
#'   the unique site-year combinations indexed by `ind_yr`.
#' * **model_name** Character. The Stan model name selected
#'   (`"pCap_all_sites"`, `"pCap_one_site"`, or `"pCap_one_site_skew"`).
#' * **skew** The value passed to the `skew` argument.
#' * **site_selection** The value passed to the `site_selection` argument.
#'
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


#' Prepare inputs for the abundance model
#'
#' @details Prepares the data list, initial values, spline basis, indexing
#'   structures, and logit-scale capture probability priors (`lt_pCap_mu`,
#'   `lt_pCap_sd`) needed to fit the BT-SPAS-X abundance BUGS model for a
#'   single site and run year. Internally calls `generate_lt_pCap_Us()` using
#'   the supplied pCap model object.
#'
#' @param site Character. The site name to prepare inputs for (e.g.
#'   `"ucc"`, `"knights landing"`).
#' @param run_year Numeric. The run year to prepare inputs for.
#' @param pCap_model_type Character. The type of pCap model whose output is
#'   supplied in `pCap_model_object`. One of `"all_sites"`,
#'   `"one_site_skew"`, or `"one_site"`. Used to determine how to extract
#'   posterior samples and whether a skew-normal process-error distribution
#'   was used.
#' @param min_pCap Numeric. Baseline minimum pCap value used to compute
#'   upper bounds and initial values for abundance (`N`). Overridden by the
#'   minimum observed Lincoln-Peterson efficiency at the site if that is
#'   available and non-zero. Default `0.0005`.
#' @param min_pCap_mult Numeric. Multiplier applied to `min_pCap` (or to
#'   the observed minimum efficiency) to allow sensitivity testing. Values
#'   less than 1 relax the upper bound on `N`. Default `1.0`.
#' @param pCap_model_object A fitted Stan model object produced by
#'   `fit_pCap_model()`. The model type must match `pCap_model_type`.
#' @param input_catch_data Optional data frame of weekly catch data.
#'   Defaults to `SRJPEdata::weekly_juvenile_abundance_catch_data`.
#'   If supplied, the structure must match that of the default.
#' @param input_efficiency_data Optional data frame of weekly efficiency
#'   (mark-recapture) data. Defaults to
#'   `SRJPEdata::weekly_juvenile_abundance_efficiency_data`.
#'   If supplied, the structure must match that of the default.
#'
#' @returns A named list with the following elements:
#' * **inputs** A list with `data`, `inits`, and `parameters` for passing
#'   to `fit_abundance_model_BUGS()`.
#' * **lt_pCap_U_data** Intermediate indexing and flow values used by
#'   `generate_lt_pCap_Us()` to compute the logit-scale pCap priors.
#' * **lt_pCap_Us** Named list with `lt_pCap_mu` and `lt_pCap_sd` vectors
#'   of length `Nstrata`, the logit-scale pCap mean and SD passed into the
#'   BUGS data list.
#' * **lp_data** Data frame of weekly catch and efficiency data used for
#'   plotting, including Lincoln-Peterson abundance and efficiency estimates.
#' * **model_name** Character. The abundance model variant selected based
#'   on mark-recapture data availability: `"all_mark_recap"`,
#'   `"missing_mark_recap"`, `"no_mark_recap"`, or
#'   `"no_mark_recap_no_trib"`.
#' * **pCap_model_type** The value passed to `pCap_model_type`.
#' * **site** The value passed to `site`.
#' * **run_year** The value passed to `run_year`.
#' * **catch_flow_raw** Numeric vector of raw flow (cfs) for each week in
#'   the modelled period.
#' * **mr_flow_raw** Numeric vector of raw flow (cfs) at the time of each
#'   mark-recapture trial used to fit the pCap model.
#' * **weeks_fit** Integer vector of Julian weeks included in the model.
#' * **week_date** Character vector of calendar dates associated with
#'   `weeks_fit`.
#' * **sites_fit** Character vector of sites whose pCap posteriors are
#'   available from the supplied pCap model object.
#' * **run_year_id** Integer. Index into the `yr_re` vector in the pCap
#'   model corresponding to this site–run year combination.
#' * **year_sd_id** Integer. Index into `yr_sd_P` in the pCap model
#'   corresponding to this site.
#' * **min_pCap** Numeric. The effective minimum pCap value after applying
#'   `min_pCap_mult`.
#'
#' @family Prepare Model Inputs
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


#' Build B-spline basis matrix for the BT-SPAS-X abundance model
#'
#' @details Constructs the B-spline basis matrix and associated parameters
#'   used to model the smooth temporal trend in weekly abundance. Knots are
#'   placed at evenly-spaced positions between the second and second-to-last
#'   week, keeping the tails free of knots to reduce boundary effects.
#'   Called internally by `prepare_abundance_inputs()`.
#'
#' @param number_weeks_catch Integer. Total number of weekly strata in the
#'   modelled period (including weeks with no catch).
#' @param k_int Integer. Number of data points per knot. A cubic spline has
#'   4 parameters, so the rule of thumb is `k_int = 4` (one knot per 4
#'   observations). Passed directly from `prepare_abundance_inputs()`.
#'
#' @returns A named list with:
#' * **K** Integer. Total number of columns in the B-spline basis matrix
#'   (number of knots + polynomial degree).
#' * **b_spline_matrix** Numeric matrix of dimensions
#'   `number_weeks_catch` × `K`. Each row is the basis vector for one
#'   weekly stratum.
#' * **knot_positions** Numeric vector of knot positions on the 1-to-`number_weeks_catch`
#'   index scale.
#'
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


#' Fit the pCap Stan model
#'
#' @details Selects the appropriate Stan model code from
#'   `SRJPEmodel::bt_spas_x_model_code` based on the model name stored in
#'   `input`, then calls `rstan::stan()` with fixed MCMC settings (10,000
#'   iterations, seed 84735, max tree depth 15, using all available cores).
#'
#' @param input A named list produced by `prepare_pCap_inputs()`. Must
#'   contain:
#'   * `model_name` — character string identifying which Stan model to run
#'     (`"pCap_all_sites"`, `"pCap_one_site"`, or `"pCap_one_site_skew"`).
#'   * `inputs$data` — the Stan data list.
#'   * `inputs$inits` — initial values (list of 3 chains).
#'   * `inputs$parameters` — character vector of parameters to monitor.
#'
#' @returns A `stanfit` object containing the posterior samples for the
#'   pCap model. Pass this object to `generate_lt_pCap_Us()` via
#'   `prepare_abundance_inputs()`.
#'
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

#' Generate logit-scale pCap priors for the abundance model
#'
#' @details Extracts posterior samples from the fitted pCap Stan model and
#'   uses them to compute, for each weekly stratum, the mean and standard
#'   deviation of the logit-scale capture probability distribution
#'   (`lt_pCap_mu`, `lt_pCap_sd`). These are passed into the BUGS data list
#'   as informative priors for `lt_pCap_U`.
#'
#'   The computation differs by abundance model variant:
#'   * `"all_mark_recap"` — uses the estimated `logit_pCap` posterior
#'     directly for every stratum.
#'   * `"missing_mark_recap"` — uses `logit_pCap` for strata with
#'     efficiency trials and simulates from the flow-regression + process
#'     error for strata without.
#'   * `"no_mark_recap"` — simulates all strata from the flow-regression
#'     with both process error and year random effect.
#'   * `"no_mark_recap_no_trib"` — as above, but draws intercept and slope
#'     from the hyper-distribution (`trib_mu_P`, `trib_sd_P`,
#'     `flow_mu_P`, `flow_sd_P`) because no site-specific pCap was
#'     estimated.
#'
#'   Effort adjustment (`log(effort)`) is added to modelled strata in all
#'   variants except `"all_mark_recap"`.
#'
#' @param abundance_inputs A named list produced by
#'   `prepare_abundance_inputs()`. The relevant sub-elements are
#'   `model_name`, `pCap_model_type`, `inputs$data`, `lt_pCap_U_data`,
#'   `run_year_id`, and `year_sd_id`.
#' @param pCap_model_object A fitted `stanfit` object produced by
#'   `fit_pCap_model()`. Must match the `pCap_model_type` recorded in
#'   `abundance_inputs`.
#'
#' @returns A named list with two numeric vectors, each of length `Nstrata`:
#' * **lt_pCap_mu** Per-stratum posterior mean of `lt_pCap_U` on the logit
#'   scale.
#' * **lt_pCap_sd** Per-stratum posterior standard deviation of `lt_pCap_U`
#'   on the logit scale.
#'
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


#' Fit the BT-SPAS-X abundance model in WinBUGS
#'
#' @details Augments the data list from `prepare_abundance_inputs()` with the
#'   logit-scale pCap priors (`lt_pCap_mu`, `lt_pCap_tau`), writes the BUGS
#'   model code from `SRJPEmodel::bt_spas_x_model_code$abundance_BUGS` to a
#'   temporary file, and calls `R2WinBUGS::bugs()` with fixed MCMC settings
#'   (3 chains, 2000 iterations, 500 burn-in, thinning by 2).
#'
#' @param abundance_inputs A named list produced by
#'   `prepare_abundance_inputs()`. Must contain `inputs$data`, `inputs$inits`,
#'   and `lt_pCap_Us` (the output of `generate_lt_pCap_Us()`).
#' @param bugs_directory Character. File path to the WinBUGS 1.4 executable
#'   directory (e.g. `"C:/Program Files/WinBUGS14/"`).
#'
#' @returns A `bugs` object (from `R2WinBUGS`) containing posterior summaries
#'   and MCMC samples for the parameters `lt_pCap_U`, `pCap_U`, `Usp`, `N`,
#'   `Ntot`, `sd.N`, `sd.Ne`, and `lg_CumN`. Pass this to
#'   `extract_abundance_estimates()`.
#'
#' @family Fit model
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


#' Extract parameter estimates from the abundance BUGS model object
#'
#' @details Tidies the summary table from a fitted `bugs` object into a long
#'   data frame, attaches week dates and stream names, and records the
#'   `SRJPEdata` package version used. The resulting table is in a format
#'   suitable for writing to the database or passing to downstream summaries.
#'
#' @param abundance_inputs A named list produced by
#'   `prepare_abundance_inputs()`, containing at minimum `site`, `run_year`,
#'   `model_name`, and `weeks_fit`.
#' @param model_object A `bugs` object produced by
#'   `fit_abundance_model_BUGS()`.
#'
#' @returns A tidy data frame with one row per parameter–statistic
#'   combination, containing the following columns:
#' * **model_name** Character. The abundance model variant
#'   (e.g. `"missing_mark_recap"`).
#' * **site** Character. The site name.
#' * **stream** Character. The stream associated with the site.
#' * **run_year** Numeric. The run year.
#' * **week_fit** Integer. The Julian week associated with the parameter
#'   (NA for scalars such as `Ntot`).
#' * **parameter** Character. The parameter name with index notation
#'   removed (e.g. `"N"`, `"lt_pCap_U"`, `"Ntot"`).
#' * **statistic** Character. The posterior summary statistic
#'   (e.g. `"mean"`, `"sd"`, `"2.5"`, `"50"`, `"97.5"`, `"Rhat"`).
#' * **value** Numeric. The value of the statistic.
#' * **srjpedata_version** Character. The version of the `SRJPEdata` package
#'   used when producing these estimates.
#'
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


#' Extract parameter estimates from the pCap Stan model object
#'
#' @details Tidies the summary table from a fitted `stanfit` pCap model into
#'   a long data frame, attaches site and stream name lookups, and records the
#'   `SRJPEdata` package version. For `one_site` (mainstem skew) models, the
#'   site name is taken from `pCap_inputs$site` rather than the parameter
#'   index. The resulting table is in the same format as
#'   `extract_abundance_estimates()`.
#'
#' @param model_object A `stanfit` object produced by `fit_pCap_model()`.
#' @param pCap_inputs A named list produced by `prepare_pCap_inputs()`,
#'   containing at minimum `sites_fit` and `model_name`.
#'
#' @returns A tidy data frame with one row per parameter–statistic
#'   combination, containing the following columns:
#' * **model_name** Character. The pCap model name
#'   (e.g. `"pCap_all_sites"`, `"pCap_one_site_skew"`).
#' * **site** Character. The site name associated with the parameter, or
#'   `NA` for hyper-parameters not tied to a specific site.
#' * **stream** Character. The stream associated with the site, or `NA`.
#' * **run_year** Always `NA` for the pCap model (pCap is estimated
#'   across all years).
#' * **week_fit** Always `NA` for the pCap model (no weekly parameters
#'   are extracted).
#' * **parameter** Character. The parameter name with index notation
#'   removed (e.g. `"b0_pCap"`, `"b_flow"`, `"logit_pCap"`).
#' * **statistic** Character. The posterior summary statistic
#'   (e.g. `"mean"`, `"sd"`, `"2.5"`, `"50"`, `"97.5"`, `"Rhat"`).
#' * **value** Numeric. The value of the statistic.
#' * **srjpedata_version** Character. The version of the `SRJPEdata` package
#'   used when producing these estimates.
#'
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

#' Run BT-SPAS-X for all JPE sites and run years
#'
#' @details Iterates over all site–run year combinations in `sites_to_run`
#'   and calls `run_abundance_workflow()` for each, with progress reporting.
#'   Optionally re-fits the pCap model before running abundance. Results from
#'   all sites are row-bound into a single tidy data frame.
#'
#' @param sites_to_run Data frame with (at minimum) columns `site` and
#'   `run_year` specifying which site–run year combinations to model.
#' @param run_pCap Logical. If `TRUE`, the pCap model is re-fit before
#'   running the abundance models and saved to `pCap_model_object_filepath`.
#'   Default `FALSE`.
#' @param mainstem Passed to `prepare_pCap_inputs()` when `run_pCap = TRUE`.
#' @param pCap_model_object_filepath Character. File path (`.rds`) from which
#'   the pCap model object is read (and to which it is written if
#'   `run_pCap = TRUE`).
#' @param bugs_model_file Character. File path to the WinBUGS model code
#'   (`.bug` file) for the abundance model.
#' @param bugs_directory Character. File path to the WinBUGS 1.4 executable
#'   directory.
#'
#' @returns A tidy data frame combining the output of
#'   `extract_abundance_estimates()` across all site–run year combinations,
#'   with columns:
#' * **model_name**, **site**, **stream**, **run_year**, **week_fit**,
#'   **parameter**, **statistic**, **value**, **srjpedata_version**
#'
#'   Rows where the model errored are included with only `site`, `run_year`,
#'   and `error = TRUE` columns (from `run_abundance_workflow()`).
#'
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

#' Run the full abundance estimation workflow for a single site and run year
#'
#' @details Executes the complete BT-SPAS-X pipeline for one site–run year
#'   combination in sequence: prepares abundance inputs, loads the pCap model
#'   object, generates logit-scale pCap priors, fits the BUGS abundance model,
#'   and extracts tidy results. Errors are caught and returned as a one-row
#'   tibble with an `error` flag so that `run_bt_spas_x_JPE_sites()` can
#'   continue to the next site without stopping.
#'
#' @param site Character. The site name.
#' @param run_year Numeric. The run year.
#' @param pCap_model_object_filepath Character. File path to the saved pCap
#'   `stanfit` object (`.rds`), read with `readRDS()`.
#' @param bugs_model_file Character. File path to the WinBUGS model code
#'   (`.bug` file) for the abundance model. Passed to
#'   `fit_abundance_model_BUGS()`.
#' @param bugs_directory Character. File path to the WinBUGS 1.4 executable
#'   directory. Passed to `fit_abundance_model_BUGS()`.
#'
#' @returns On success, the tidy data frame returned by
#'   `extract_abundance_estimates()`. On error, a one-row tibble with columns
#'   `site`, `run_year`, and `error = TRUE`.
#'
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
