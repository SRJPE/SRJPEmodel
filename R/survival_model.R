
#' Prepare inputs for survival model
#' @details This function prepares the data, initial values, and parameter list needed to call the
#' survival STAN model. It will prepare these data based on what model you want to run. This
#' iteration of the model code runs the model using `max_flow` and `fork_length` as covariates. A future
#' version of the model code will allow for the full suite of covariates and biological effects.
#' @returns a named list containing the inputs to the `run_survival_model()` function:
#' * **data_inputs** a named list of data inputs for passing to the STAN model.
#' * **inits** a list of initial values for passing to the STAN model.
#' * **parameters** a vector of parameter names to save from calling the STAN model.
#' @family Prepare Model Inputs
#' @export
#'
#' @md
prepare_survival_inputs <- function(){

  cli::cli_bullets("Preparing inputs to run the model for maximum flow and fork length.")

  # prepare data for sac- and feather/butte-specific models
  sac_data_list <- create_sac_survival_data()
  fb_data_list <- create_butte_feather_survival_data(sac_data_list) # use mean/sd sac fork lengths to standardize butte/feather fork lengths
  # prepare data for travel time models
  tt_data_list <- create_tt_data()
  # workflow for running MaxFlow_FL version of the model
  model_name <- "CovIndCont"
  par_save_list <- c("TT_bCov", "TT_bCovT", "TT_RE", "TT_RET", "TTRE_sd", "TTRE_sdT",
                     "Pro_sd", "Pro_sdT", "T_bSz", "TT_reach", "TT_RelSac", "TT_reachT",
                     "P_b", "muPb", "sdPb", "S_bReach", "S_bTrib", "T_bReach", "T_bTrib", "S_bCov", "S_bCovT", "S_bSz",
                     "RE_sd", "RE_sdT", "S_RE", "S_REt",'log_lik',
                     "pred_surv", "SurvRelSac", "SurvRelSacSz", "pred_survT", "pred_survTSz",
                     "pred_surv_per100", "pred_survT_per100", "pred_pcap",
                     "SurvForecast", "SurvForecastSz", "SurvForecast_nore", "SurvForecastSz_nore",
                     "TribSurvForecast", "TribSurvForecastSz", "TribSurvForecast_nore", "TribSurvForecastSz_nore",
                     "SurvForecastSz_rst", "TribSurvForecastSz_rst", "TTForecastSz_rst", "TTForecastTSz_rst",
                     "TTForecast", "TTForecastT", "TTForecastSz", "TTForecastTSz")

  inits <- list(P_b = matrix(data = 2.2, nrow = sac_data_list$n_years, ncol = sac_data_list$n_reaches),
                muPb = rep(0, sac_data_list$n_reaches),
                sdPb = rep(1.0, sac_data_list$n_reaches),
                S_bReach = c(1.5, 0.5, 0.5, 0.5))

  # prepare data list
  model_data <- list(
    # general variables
    "Nyrs" = sac_data_list$n_years,
    "rch_covind" = c(1, 2, 3, 4),
    # sac-specific survival variables
    "Nind" = sac_data_list$n_ind,
    "Nreaches" = sac_data_list$n_reaches,
    "Ndetlocs" = sac_data_list$n_detection_locations,
    "Rmult" = sac_data_list$Rmult,
    "RmultSac" = sac_data_list$RmultSac,
    "RmultTis" = sac_data_list$RmultTis,
    "Rmultrst" = sac_data_list$Rmultrst,
    "Nrst" = sac_data_list$Nrst,
    "ReachKMrst" = sac_data_list$ReachKMrst,
    "Nrg" = sac_data_list$n_release_groups,
    "CH" = sac_data_list$capture_history,
    "yrind" = sac_data_list$year_index,
    "rgind" = sac_data_list$release_group_index,
    "rgwy_ind" = rep(0, sac_data_list$n_release_groups),
    "firstCap" = sac_data_list$first_capture,
    "lastCap" = sac_data_list$last_capture,
    "CovX" = sac_data_list$cov_x,
    "Sz" = sac_data_list$sz,
    "Xvec" = sac_data_list$x_vec,
    # feather/butte-specific survival variables
    "RmultTrst" = fb_data_list$RmultTrst,
    "NrstT" = fb_data_list$Nrst,
    "TribForRST" = fb_data_list$TribForRST,
    "ReachKMTrst" = fb_data_list$ReachKMTrst,
    "Ntribs" = fb_data_list$n_tribs,
    "NindT" = fb_data_list$n_ind,
    "NreachesT" = fb_data_list$n_reaches,
    "NrgT" = fb_data_list$n_release_groups,
    "trib_ind" = fb_data_list$trib_ind,
    "firstCapT" = fb_data_list$first_capture,
    "lastCapT" = fb_data_list$last_capture,
    "RmultT" = fb_data_list$Rmult,
    "RmultTrib" = fb_data_list$RmultTrib,
    "yrindT" = fb_data_list$year_index,
    "rgindT" = fb_data_list$release_group_index,
    "CHT" = fb_data_list$capture_history,
    "trib_rg" = fb_data_list$trib_index_for_release_groups,
    "rgwy_indT" = rep(0, fb_data_list$n_release_groups),
    "CovXT" = fb_data_list$cov_x,
    "SzT" = fb_data_list$sz,
    "XvecT" = fb_data_list$x_vec,
    # covariate/effect specifying variables
    "UseSizeEffect" = 1,
    "NS_bCov" = 0,
    "Nsz" = 25,
    "Xsz" = fb_data_list$x_sz,
    "NsX" = 25,
    # sac-specific travel time variables
    "Nobs" = tt_data_list$Nobs,
    "ObsTT" = tt_data_list$ObsTT,
    "TTind" = tt_data_list$TTind,
    "ReachKM" = tt_data_list$ReachKM,
    "ReachKM_ind" = tt_data_list$ReachKM_ind,
    # feather-butte travel time variables
    "NobsT" = tt_data_list$NobsT,
    "ObsTTT" = tt_data_list$ObsTTT,
    "TTindT" = tt_data_list$TTindT,
    "ReachKMT" = tt_data_list$ReachKMT,
    "ReachKMT_ind" = tt_data_list$ReachKMT_ind
  )

  inits <- list(inits, inits, inits)

  return(list("inputs" = list("data" = model_data,
                              "inits" = inits,
                              "parameter" = par_save_list),
              "model_name" = model_name))
}

# Create data list for Sacramento model
#' @details This prepares Sacramento-specific inputs within the `prepare_survival_inputs()` function.
#' @returns A named list of some of the data inputs required for the STAN model.
#' @export
#' @md
create_sac_survival_data <- function() {

  # general variables across sac and feather/butte
  all_study_years <- tibble(year = sort(unique(SRJPEdata::survival_model_inputs$year))) |>  # all sac and tribs years combined
    mutate(year_index = row_number())

  # prepare flow covariates
  exceedence_flows <- SRJPEdata::forecast_covariates |>
    filter(name == "3_category_flow_exceedance_year_type",
           stream %in% c("sacramento river", "feather river", "butte creek")) |>
    mutate(value = case_when(text_value == "Wet" ~ 2,
                             text_value == "Average" ~ 1,
                             TRUE ~ 0)) |>
    # Flora's code calls year the same thing as water year (assumption)
    select(stream, year = water_year, month, value) |>
    arrange(stream, year, month)

  max_flows <- SRJPEdata::forecast_covariates |>
    filter(name == "monthly_max_flow",
           stream %in% c("sacramento river", "feather river", "butte creek")) |>
    arrange(stream, year, month) |>
    select(stream, year, month, value)

  # sacramento data
  sac_data <- SRJPEdata::survival_model_inputs |>
    filter(is.na(trib_ind)) |>
    mutate(month = month(as.Date(fish_release_date,format="%m/%d/%Y"))) |>
    left_join(exceedence_flows |>
                filter(stream == "sacramento river") |>
                select(year, exceedence = value),
              by = "year") |>
    left_join(max_flows |>
                filter(stream == "sacramento river") |>
                select(year, month, max_flow = value),
              by = c("month", "year")) |>
    # TODO confirm order of scaling
    mutate(max_flow_standardized = as.numeric(scale(max_flow)))



  data_list <- list(
    # general variables
    n_years = nrow(all_study_years),
    # survival and travel time variables (sacramento-specific)
    n_ind = nrow(sac_data),
    n_reaches = 4, # 1 = release-woodson, 2 = woodson-butte, 3 = butte-sac, 4 = sac-delta
    n_detection_locations = 5,
    n_release_groups = length(unique(sac_data$study_id)),

    Rmult = sac_data |>
      select(dist_rlwoodson_standardized, dist_woodsonbutte_standardized,
             dist_buttesac_standardized, dist_sacdelta_standardized) |>
      as.matrix() |>
      unname(),

    RmultSac = 0.43,
    RmultTis = 1.23,
    ReachKMrst = c(124.9512976, 93.2391355490999, 26.1273401272, 18.3412666561999),
    Rmultrst = c(124.9512976, 93.2391355490999, 26.1273401272, 18.3412666561999) / 100,
    Nrst = 4,

    capture_history = sac_data |>
      separate_wider_position(ch, widths = c("1" = 1, "2" = 1, "3" = 1, "4" = 1, "5" = 1)) |>
      select(`1`, `2`, `3`, `4`, `5`) |>
      mutate_all(as.numeric) |>
      as.matrix() |>
      unname(),

    year_index = sac_data |>
      left_join(all_study_years, by = "year") |>
      dplyr::pull(year_index),

    release_group_index = sac_data |>
      arrange(year) |>
      group_by(study_id) |>
      mutate(index = cur_group_id()) |>
      ungroup() |>
      dplyr::pull(index),

    fish_id_max_flow_index = sac_data |>
      mutate(max_flow_standardized = scale(max_flow)) |>
      group_by(fish_id) |>
      dplyr::reframe(ind = unique(max_flow_standardized),
                     year = unique(year)) |>
      arrange(year) |>
      ungroup() |>
      pull(ind),

    release_group_flow_exceed_index = sac_data |>
      group_by(study_id) |>
      dplyr::reframe(ind = unique(exceedence),
                     year = unique(year)) |>
      arrange(year) |>
      ungroup() |>
      pull(ind),

    release_group_water_year_3_index = sac_data |>
      arrange(year) |>
      group_by(study_id) |>
      summarise(ind = unique(wy_three_categories)) |>
      ungroup() |>
      dplyr::pull(ind),

    release_group_water_year_2_index = sac_data |>
      arrange(year) |>
      group_by(study_id) |>
      summarise(ind = unique(wy_two_categories)) |>
      ungroup() |>
      dplyr::pull(ind),

    first_capture = sac_data$first_capture,
    last_capture = sac_data$last_capture,

    # covariates
    fork_length = sac_data$fish_length,
    weight = sac_data$fish_weight,
    condition = sac_data$fish_k,

    cov_x = data.frame(cbind(sac_data$max_flow_standardized, sac_data$max_flow_standardized,
                             sac_data$max_flow_standardized, sac_data$max_flow_standardized)), # flora's code says max flow sac x 3 and then maxflowdelta, but maxflowdelta is assigned to maxflowsac
    mu_x = mean(sac_data$max_flow, na.rm = T),
    sd_x = sd(sac_data$max_flow, na.rm = T),
    mu_sz = mean(sac_data$fish_length, na.rm = T),
    sd_sz = sd(sac_data$fish_length, na.rm = T),
    sz = (sac_data$fish_length - mean(sac_data$fish_length, na.rm = T)) / sd(sac_data$max_flow, na.rm = T),
    n_sz = 25,
    x_vec = (seq(from = 4000, to = 40000, length.out = 25) - mean(sac_data$fish_length, na.rm = T) / sd(sac_data$fish_length, na.rm = T))
    )

  return(data_list)

}

# Create data list for Butte/Feather model
#' @details This prepares Butte/Feather-specific inputs within the `prepare_survival_inputs()` function.
#' @returns A named list of some of the data inputs required for the STAN model.
#' @export
#' @md
create_butte_feather_survival_data <- function(sac_data_list) {

  # general variables across sac and feather/butte
  all_study_years <- tibble(year = sort(unique(SRJPEdata::survival_model_inputs$year))) |>  # all sac and tribs years combined
    mutate(year_index = row_number())

  # prepare flow covariates
  exceedence_flows <- SRJPEdata::forecast_covariates |>
    filter(name == "3_category_flow_exceedance_year_type",
           stream %in% c("sacramento river", "feather river", "butte creek")) |>
    mutate(value = case_when(text_value == "Wet" ~ 2,
                             text_value == "Average" ~ 1,
                             TRUE ~ 0)) |>
    # Flora's code calls year the same thing as water year (assumption)
    select(stream, year = water_year, month, value) |>
    arrange(stream, year, month)

  max_flows <- SRJPEdata::forecast_covariates |>
    filter(name == "monthly_max_flow",
           stream %in% c("sacramento river", "feather river", "butte creek")) |>
    arrange(stream, year, month) |>
    select(stream, year, month, value)

  feather_butte_data <- SRJPEdata::survival_model_inputs |>
    filter(!is.na(trib_ind)) |>
    mutate(across(fish_weight, ~ replace_na(., mean(., na.rm=TRUE))),
           across(fish_k, ~ replace_na(., mean(., na.rm=TRUE)))) |>
    mutate(month = month(as.Date(fish_release_date,format="%m/%d/%Y"))) |>
    left_join(exceedence_flows |>
                filter(stream != "sacramento river") |>
                pivot_wider(names_from = stream, values_from = value) |>
                select(year, exceedence_but = `butte creek`,
                       exceedence_fea = `feather river`),
              by = "year") |>
    left_join(max_flows |>
                filter(stream != "sacramento river") |>
                pivot_wider(names_from = stream, values_from = value) |>
                select(year, month, max_flow_but = `butte creek`,
                       max_flow_fea = `feather river`),
              by = c("month", "year")) |>
    mutate(exceedence = ifelse(rl %in% c("1B","2B","3B","4B","5B","6B"), exceedence_but, exceedence_fea),
           max_flow = ifelse(rl %in% c("1B","2B","3B","4B","5B","6B"), max_flow_but, max_flow_fea)) |>
    select(-c(exceedence_but, exceedence_fea, max_flow_but, max_flow_fea)) |>
    # TODO confirm order of scaling
    mutate(max_flow_standardized = as.numeric(scale(max_flow)))

  data_list <- list(
    # feather/butte
    # Upper Butte 2019 does not have weight information so use average values across all release group instead.
    n_ind = nrow(feather_butte_data),
    n_tribs = 2,
    n_reaches = 2,
    n_release_groups = length(unique(feather_butte_data$study_id)),
    trib_ind = feather_butte_data$trib_ind,

    Nrst = 3,
    Rmult = feather_butte_data$dist_rlsac_standardized,
    RmultTrib = c(1.17, 0.92),
    ReachKMTrst = c(201.583633469999,89.5817927391999, 131.0731),
    RmultTrst = c(201.583633469999,89.5817927391999, 131.0731) / 100,
    TribForRST = c(1,2,2),

    capture_history = feather_butte_data |>
      separate_wider_position(ch, widths = c("1" = 1, "2" = 1, "3" = 1)) |>
      select(`1`, `2`, `3`) |>
      mutate_all(as.numeric) |>
      as.matrix() |>
      unname(),

    trib_index_for_release_groups = feather_butte_data |>
      arrange(year) |>
      dplyr::distinct(study_id, trib_ind) |>
      dplyr::pull(trib_ind),

    # add month for trib model
    fish_id_max_flow_index = feather_butte_data |>
      mutate(max_flow_standardized = scale(max_flow)) |>
      group_by(fish_id) |>
      dplyr::reframe(ind = unique(max_flow_standardized),
                     year = unique(year),
                     month = unique(month)) |>
      arrange(year, month) |>
      ungroup() |>
      pull(ind),

    year_index = feather_butte_data |>
      left_join(all_study_years, by = "year") |>
      dplyr::pull(year_index),

    release_group_index = feather_butte_data |>
      arrange(year) |>
      group_by(study_id) |>
      mutate(index = cur_group_id()) |>
      ungroup() |>
      dplyr::pull(index),

    release_group_flow_exceed_index = feather_butte_data |>
      arrange(year) |>
      group_by(study_id) |>
      dplyr::reframe(ind = unique(exceedence),
                     year = unique(year)) |>
      arrange(year) |>
      ungroup() |>
      pull(ind),

    release_group_water_year_3_index = feather_butte_data |>
      arrange(year) |>
      group_by(study_id) |>
      summarise(ind = unique(wy_three_categories)) |>
      ungroup() |>
      dplyr::pull(ind),

    release_group_water_year_2_index = feather_butte_data |>
      arrange(year) |>
      group_by(study_id) |>
      summarise(ind = unique(wy_two_categories)) |>
      ungroup() |>
      dplyr::pull(ind),

    first_capture = feather_butte_data$first_capture,
    last_capture = feather_butte_data$last_capture,

    # covariates

    fork_length = feather_butte_data$fish_length,
    weight = feather_butte_data$fish_weight,
    condition = feather_butte_data$fish_k,

    cov_x = data.frame(cbind(feather_butte_data$max_flow_standardized, feather_butte_data$max_flow_standardized)), # flora's code calls the second object deltaflowT but it is same as max flow t
    mu_x = mean(feather_butte_data$max_flow, na.rm = T),
    sd_x = sd(feather_butte_data$max_flow, na.rm = T),
    sz = (feather_butte_data$fish_length - sac_data_list$mu_sz) / sac_data_list$sd_sz,
    n_sX = 25,
    x_sz = (seq(from = 10, to = 150, length.out = 25) - sac_data_list$mu_sz) / sac_data_list$sd_sz,
    x_vec = (seq(from = 150, to = 30500, length.out = 25) - mean(feather_butte_data$fish_length, na.rm = T) / sd(feather_butte_data$fish_length, na.rm = T))
  )

  return(data_list)

}


# Create data list for travel time model
#' @details This prepares travel time inputs within the `prepare_survival_inputs()` function.
#' @returns A named list of some of the data inputs required for the STAN model.
#' @export
#' @md
create_tt_data <- function() {

  # create data
  sac_data <- SRJPEdata::survival_model_inputs |>
    filter(is.na(trib_ind))

  fb_data <- SRJPEdata::survival_model_inputs |>
    filter(!is.na(trib_ind))

  # sacramento
  reach_km_ind <- data.frame(cbind(sac_data$dist_rlwoodson,
                                   sac_data$dist_woodsonbutte,
                                   sac_data$dist_buttesac,
                                   sac_data$dist_sacdelta))

  #identify records with one or more detections after release and get their FishID.
  #Fish not seen after release provide no data for travel time
  TTfR <- SRJPEdata::detection_history_sacramento |>
    select(-c(stream, fish_id))

  Nobs_df <- TTfR |>
    mutate(Obs = 4 - rowSums(is.na(across(everything()))))

  ObsTT <- TTfR |>
    rowwise() |>
    mutate(row_vec = list(c_across(everything()))) |>
    pull(row_vec) |>
    unlist() |>
    na.omit()

  vec <- as.vector(t(TTfR))                    # row-wise flattening
  counter <- seq_len(sum(!is.na(vec)))       # counter for non-NA values
  vec[!is.na(vec)] <- counter                # replace only non-NA

  # Rebuild dataframe with same shape and column names
  TTind <- matrix(vec, nrow = nrow(TTfR), byrow = TRUE) |>
    as.data.frame() |>
    mutate(across(everything(), ~ replace_na(., 0)))

  # feather/butte
  ReachKMT <- data.frame(rbind(cbind(117, 110), cbind(92, 110)))
  ReachKMT_ind <- data.frame(cbind(fb_data$dist_rlsac,
                                   fb_data$dist_sacdelta))
  ReachKMTrst <- c(201.583633469999,89.5817927391999, 131.0731)

  TTfRT <- SRJPEdata::detection_history_feather_butte |>
    select(-c(stream, fish_id))

  NobsT_df <- TTfRT |>
    mutate(Obs = 2 - rowSums(is.na(across(everything()))))

  ObsTTT <- TTfRT |>
    rowwise() |>
    mutate(row_vec = list(c_across(everything()))) |>
    pull(row_vec) |>
    unlist() |>
    na.omit()

  vecT <- as.vector(t(TTfRT))                    # row-wise flattening
  counterT <- seq_len(sum(!is.na(vecT)))       # counter for non-NA values
  vecT[!is.na(vecT)] <- counterT                # replace only non-NA

  # Rebuild dataframe with same shape and column names
  TTindT <- matrix(vecT, nrow = nrow(TTfRT), byrow = TRUE) |>
    as.data.frame() |>
    mutate(across(everything(), ~ replace_na(., 0)))

  data_list <- list(
    # sac-specific travel time variables
    "Nobs" = sum(Nobs_df$Obs),
    "ObsTT" = ObsTT,
    "TTind" = TTind,
    "ReachKM" = c(45, 88, 170, 110),
    "ReachKM_ind" = reach_km_ind,
    # feather-butte travel time variables
    "NobsT" = sum(NobsT_df$Obs),
    "ObsTTT" = ObsTTT,
    "TTindT" = TTindT,
    "ReachKMT" = data.frame(rbind(cbind(117, 110), cbind(92, 110))),
    "ReachKMT_ind" = ReachKMT_ind
    )

    return(data_list)
}

#' Fit Survival Model
#' @details This model calls a survival STAN model that estimates survival in the Sacramento and Feather/Butte systems.
#' @param survival_inputs The named list object generated by running `prepare_survival_inputs()`.
#' @returns a STANfit object with the survival model fit.
#' @family Fit model
#' @export
#' @md
fit_survival_model <- function(survival_inputs) {

  # select model to run (either covariate or no covariate)
  if(survival_inputs$model_name == "CovIndCont") {
    stan_model <- eval(parse(text = "SRJPEmodel::survival_model_code$survival_CovIndCont"))
  } else if(survival_inputs$model_name == "CovIndWY") {
    stan_model <- eval(parse(text = "SRJPEmodel::survival_model_code$survival_CovIndWY"))
  } else {
    stan_model <- eval(parse(text = "SRJPEmodel::survival_model_code$survival_NoCov"))
  }


  options(mc.cores = parallel::detectCores())

  cli::cli_process_start("Fitting STAN survival model")
  fit <- rstan::stan(model_code = stan_model,
                     data = survival_inputs$inputs$data,
                     init = survival_inputs$inputs$inits,
                     chains = 3,
                     iter = 2000,
                     include = T,
                     pars = survival_inputs$inputs$parameter,
                     seed = 1234)
  cli::cli_process_done("STAN survival model fitting complete")

  return(fit)
}

#' Extract parameter estimates from the survival STAN model object
#' @details This function extracts parameter estimates from the survival STAN model,
#' assuming that the model was run with water year type (either 2 or 3) and
#' biological effect of fork length (based on current model structure). It returns them
#' in a tidy data frame format.
#' @param model_object the STAN object produced by running `fit_survival_model()`.
#' @returns A table with the format:
#' * **parameter**
#' * **year**
#' * **reach_name**
#' * **release_group**
#' * **statistic**
#' * **value**
#' * **srjpedata_version**
#' @export
#' @md
extract_survival_estimates <- function(model_object) {
  # get year lookup
  year_lookup <- tibble(year = unique(SRJPEdata::survival_model_inputs$year) |>
                          sort()) |>  # all sac and tribs years combined
    mutate(year_index = row_number())

  # get reach lookup - hard-coded
  reach_lookup_sac <- tibble("sac_reach_index" = 1:4,
                             "reach_name_sac" = c("Woodson", "Butte", "Sacramento", "Delta")) # hard-coded in GetData_flora.R

  reach_lookup_trib <- tibble("trib_reach_index" = 1:2,
                              "reach_name_trib" = c("Sacramento", "Delta"))# taken from powerpoint in sharepoint

  # release group lookup
  release_group_lookup_sac <- SRJPEdata::survival_model_inputs |>
    filter(is.na(trib_ind)) |>
    arrange(year) |>
    group_by(study_id) |>
    mutate(release_group_index_sac = cur_group_id()) |>
    ungroup() |>
    dplyr::distinct(release_group_sac = study_id, release_group_index_sac)

  release_group_lookup_trib <- SRJPEdata::survival_model_inputs |>
    filter(!is.na(trib_ind)) |>
    arrange(year) |>
    group_by(study_id) |>
    mutate(release_group_index_trib = cur_group_id()) |>
    ungroup() |>
    dplyr::distinct(release_group_trib = study_id, release_group_index_trib)

  formatted_table <- rstan::summary(model_object)$summary |>
    as.data.frame() |>
    tibble::rownames_to_column("parameter") |>
    # get indices to match to year, reach, etc. for parameters
    separate_wider_delim(parameter, delim = "[", names = c("parameter", "indices"),
                         too_few = "align_start") |>
    mutate(indices = str_remove_all(indices, "\\]")) |>
    separate_wider_delim(indices, ",", names = c("index_1", "index_2", "index_3"),
                         too_few = "align_start") |>
    mutate(across(index_1:index_3, as.numeric),
           # year-structured parameters
           year_index = ifelse(parameter %in% c("P_b", "pred_pcap"), index_1, NA),
           # TODO aligning for water_year_type for now, but can be any covariate, not just water year type per Flora
           water_year_type = ifelse(parameter %in% c("SurvRelSacSz", "SurvWoodSacSz", "SurvForecastSz",
                                                     "pred_SurvTSz"), index_2, NA),
           environmental_covariate_dimension = ifelse(parameter %in% c("S_bCov", "s_bCovT"), index_1, NA),
           forecast_index = ifelse(parameter == "SurvForecast", index_1, NA), # this is "covariate dimension" per Flora
           # reach-structured parameters
           sac_reach_index = case_when(parameter == "P_b" ~ index_2,
                                       parameter %in% c("muPb", "sdPb", "S_bReach") ~ index_1,
                                       parameter %in% c("pred_surv", "pred_pcap") ~ index_2,
                                       TRUE ~ NA),
           # trib-structured parameters
           trib_reach_index = case_when(parameter == "S_bTrib" ~ index_1,
                                        parameter == "pred_survT" ~ index_2,
                                        parameter %in% c("TribSurvForecastSz", "pred_SurvTSz") ~ index_3,
                                        TRUE ~ NA),
           # release group structured parameters
           # TODO aligning with release group for now, building out for fork length-based model
           # SurvRelSac is release group OR individual depending on covariate used
           # pred_surv is release group OR individual depending on covariate used
           release_group_index_sac = ifelse(parameter %in% c("SurvRelSac", "pred_surv",
                                                             "S_RE"), index_1, NA),
           # TODO aligning with release group for now, building out fork length-based model
           # pred_survT is release group OR individual depending on covariate used
           release_group_index_trib = ifelse(parameter %in% c("S_REt", "pred_survT", "SurvWoodSac"), index_1, NA),
           pred_forecast_index_trib = ifelse(parameter == "TribSurvForecast", index_1, NA),
           size_class_index = ifelse(parameter %in% c("pred_survTSz", "TribSurvForecastSz",
                                                      "SurvRelSacSz", "SurvWoodSacSz",
                                                      "SurvForecastSz", "pred_SurvTSz"), index_1, NA),
           covariate_dimension_index = ifelse(parameter %in% c("pred_survTSz", "TribSurvForecastSz"), index_2 - 1, NA)) |>
    # now join in lookups
    left_join(year_lookup, by = "year_index") |>
    left_join(reach_lookup_sac, by = "sac_reach_index") |>
    left_join(reach_lookup_trib, by = "trib_reach_index") |>
    left_join(release_group_lookup_sac, by = "release_group_index_sac") |>
    left_join(release_group_lookup_trib, by = "release_group_index_trib") |>
    # now combine reach lookups into one column
    mutate(reach_name = ifelse(is.na(reach_name_trib), reach_name_sac, reach_name_trib),
           release_group = ifelse(is.na(release_group_trib), release_group_sac, release_group_trib),
           parameter = ifelse(parameter == "TribSurvForecast",
                              paste0(parameter, "_", pred_forecast_index_trib),
                              parameter),
           parameter == ifelse(parameter == "SurvForecast",
                               paste0(parameter, "_", forecast_index),
                               parameter)) |>
    select(-c(year_index, sac_reach_index, trib_reach_index, release_group_index_sac,
              release_group_index_trib, reach_name_sac, reach_name_trib, release_group_sac,
              release_group_trib, pred_forecast_index_trib, forecast_index,
              index_1, index_2, index_3)) |>
    pivot_longer(mean:Rhat,
                 values_to = "value",
                 names_to = "statistic") |>
    mutate(srjpedata_version = as.character(packageVersion("SRJPEdata")))

  return(formatted_table)
}

#' Plot smolt survival rate
#' @details This function produces a plot showing the parameter estimates for `SurvRelSac` and `SurvWoodSac`, which
#' are survival rates for Release to Woodson Bridge and Woodson Bridge to Sacramento.
#' @param survival_estimates A table that is produced from running `extract_survival_estimates()`.
#' @returns A plot.
#' @export
#' @md
generate_survival_rate_plot <- function(survival_estimates) {
  # TODO this should be the parameter "surv", which right now is not reported...or is it correct?
  survival_estimates |>
    filter(parameter %in% c("SurvRelSac", "SurvWoodSac")) |>
    mutate(Reach = ifelse(parameter == "SurvRelSac",
                          "Release to Woodson Bridge",
                          "Woodson Bridge to Sacramento")) |>
    pivot_wider(names_from = "statistic",
                values_from = "value") |>
    ggplot() +
    geom_errorbar(aes(x = release_group, ymin = `2.5%`, ymax = `97.5%`),
                  color = "blue", width = 0.2) +
    geom_point(aes(x = release_group, y = `50%`), size = 2) +
    facet_wrap(~Reach, scales = "free_x") +
    labs(x = "Release Group",
         y = "Median Survival Rate") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1))

}
