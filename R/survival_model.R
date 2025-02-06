# refactor of Call_Model.R
# library(rstan)
# options(mc.cores = parallel::detectCores())
# rstan_options(auto_write = TRUE)

#' Prepare inputs for survival model
#' @details Runs the survival model for both Sacramento releases and Feather/Butte releases.
#' @param use_covariate either `TRUE` or `FALSE`. If `TRUE`, will call a version of the survival
#' model that fits a covariate effect (water year type). If `FALSE`, it will call a version of the survival model
#' with no covariate effect.
#' @param number_of_water_year_types either `2` or `3`
#' @param effect one of `no_size_effect`, `fork_length_effect`, `weight_effect`, `condition`,
#' @returns TODO
#' @export
#' @md
prepare_survival_inputs <- function(use_covariate,
                                    number_of_water_year_types,
                                    effect){

  # TODO replace with updated SRJPEdata
  # TODO ensure these are sorted, standardized, etc. per PrepData_flora.R
  sac_data <- read.csv(here("data-raw", "survival_model", "Sac_data.csv"),  stringsAsFactors = F)
  fb_data <-  read.csv(here("data-raw", "survival_model", "FeaBut_data.csv"),  stringsAsFactors = F)

  sac_data_list <- get_survival_data_list("sacramento", sac_data, fb_data)
  fb_data_list <- get_survival_data_list("feather_butte", sac_data, fb_data)


  if(use_covariate) {

    parameters <- c("P_b", "muPb", "sdPb", "S_bReach", "S_bTrib", "S_bCov", "S_bCovT", "S_bSz",
                    "RE_sd", "RE_sdT", "S_RE","S_REt", "pred_surv", "SurvRelSac", "SurvWoodSac",
                    "SurvForecast", "SurvRelSacSz", "SurvWoodSacSz", "SurvForecastSz", "pred_pcap",
                    "pred_survT", "pred_survTSz", "TribSurvForecast", "TribSurvForecastSz")


    inits <- list(P_b = matrix(data = 2.2, nrow = Nyrs, ncol = Nreaches),
                  muPb = rep(0, Nreaches),
                  sdPb = rep(1.0, Nreaches),
                  S_bReach = rep(0, 3),
                  S_RE = rep(-3, Nrg),
                  RE_sd = 0.5)

    UseSizeEffect <- 0
    NS_bCov <- 0
    CovX <- rep(0, sac_data_list$n_ind)
    CovXT <- rep(0, fb_data_list$n_ind)
    rgwy_ind <- rep(0, sac_data_list$n_release_groups)
    rgwy_indT <- rep(0, fb_data_list$n_release_groups)
    Sz <- rep(0, sac_data_list$n_ind)
    SzT <- rep(0, fb_data_list$n_ind)
    Xsz <- rep(0, sac_data_list$n_size_classes)

  } else {

    parameters <- c("P_b", "muPb", "sdPb", "S_bReach", "S_bTrib", "S_bCov", "S_bCovT", "S_bSz", "RE_sd",
                    "RE_sdT", "S_RE", "S_REt", "pred_surv", "SurvRelSac", "SurvWoodSac", "SurvForecast",
                    "SurvRelSacSz", "SurvWoodSacSz", "SurvForecastSz", "pred_pcap", "pred_survT",
                    "pred_survTSz", "TribSurvForecast", "TribSurvForecastSz")

    inits <- list(P_b = matrix(data = 2.2, nrow = Nyrs, ncol = Nreaches),
                  muPb = rep(0, Nreaches),
                  sdPb = rep(1.0, Nreaches),
                  S_bReach = rep(0, 3),
                  S_RE = rep(-3, Nrg),
                  RE_sd = 0.5)

    # now step through all possible covariate versions - water year type (either 2 or 3) and
    # effects (no effect, fork length, weight, and condition)

    # set water year type-specific variables
    if(number_of_water_year_types == 2) {

      CovX <- sac_data_list$water_year_2
      CovXT <- fb_data_list$water_year_2
      rgwy_ind <- sac_data_list$release_group_water_year_2_index
      rgwy_indT <- fb_data_list$release_group_water_year_2_index
    } else if (number_of_water_year_types == 3) {

      CovX <- sac_data_list$water_year_3
      CovXT <- fb_data_list$water_year_3
      rgwy_ind <- sac_data_list$release_group_water_year_3_index
      rgwy_indT <- fb_data_list$release_group_water_year_3_index
    }

    # set effect-specific variables
    if (effect == "no_size_effect") {

      Sz <- rep(0, sac_data_list$n_ind)
      SzT <- rep(0, fb_data_list$n_ind)
      Xsz <- rep(0, sac_data_list$n_size_classes)

    } else if (effect == "fork_length_effect") {

      all_fork_lengths <- c(sac_data_list$fork_length, fb_data_list$fork_length)
      mean_fork_length <- mean(all_fork_lengths, na.rm = TRUE)
      sd_fork_length <- sd(all_fork_lengths, na.rm = TRUE)
      Sz <- (sac_data_list$fork_length - mean_fork_length) / sd_fork_length
      SzT <- (fb_data_list$fork_length - mean_fork_length) / sd_fork_length
      Xsz  <- seq(from = min(c(Sz, SzT), na.rm = TRUE),
                  to = max(c(Sz, SzT), na.rm = TRUE),
                  length.out = sac_data_list$n_size_classes)

    } else if (effect == "weight_effect") {

      all_weights <- c(sac_data_list$weight, fb_data_list$weight)
      mean_weight <- mean(all_weights, na.rm = TRUE)
      sd_weight <- sd(all_weights, na.rm = TRUE)
      Sz <- (sac_data_list$weight - mean_weight) / sd_weight
      SzT <- (fb_data_list$weight - mean_weight) / sd_weight
      Xsz  <- seq(from = min(c(Sz, SzT), na.rm = TRUE),
                  to = max(c(Sz, SzT), na.rm = TRUE),
                  length.out = sac_data_list$n_size_classes)

    } else if (effect == "condition") {

      all_conditions <- c(sac_data_list$condition, fb_data_list$condition)
      mean_condition <- mean(all_conditions, na.rm = TRUE)
      sd_condition <- sd(all_conditions, na.rm = TRUE)
      Sz <- (sac_data_list$condition - mean_condition) / sd_condition
      SzT <- (fb_data_list$condition - mean_condition) / sd_condition
      Xsz  <- seq(from = min(c(Sz, SzT), na.rm = TRUE),
                  to = max(c(Sz, SzT), na.rm = TRUE),
                  length.out = sac_data_list$n_size_classes)
    }

  }

  # prepare data list
  model_data <- list(
                     # sac-specific data variables
                     Nind = sac_data_list$n_ind,
                     Nreaches = sac_data_list$n_reaches,
                     Ndetlocs = sac_data_list$n_detection_locations,
                     Rmult = sac_data_list$Rmult,
                     Nyrs = sac_data_list$n_years,
                     Nrg = sac_data_list$n_release_groups,
                     CH = sac_data_list$capture_history,
                     yrind = sac_data_list$year_index,
                     rgind = sac_data_list$release_group_index,

                     # variables set based on number_water_year_types
                     rgwy_ind = rgwy_ind,
                     NS_bCov = ifelse(number_water_year_types == 2, 1, 2),
                     CovX = CovX,

                     # variables set based on effect
                     UseSizeEffect = ifelse(effect == "no_size_effect", 0, 1),

                     firstCap = firstCap,
                     lastCap = lastCap,

                     # feather/butte-specific variables
                     NindT = fb_data_list$n_ind,
                     NrgT = fb_data_list$n_release_groups,
                     trib_ind = fb_data_list$trib_ind,
                     firstCapT = fb_data_list$first_capture,
                     lastCapT = fb_data_list$last_capture,
                     RmultT = fb_data_list$Rmult,
                     RmultTrib = RmultTrib, # TODO stopped here
                     yrindT = yrindT,
                     rgindT = rgindT,
                     CHT = CHT,
                     trib_rg = trib_rg,
                     CovXT = CovXT,
                     rgwy_indT = rgwy_indT,
                     Sz = Sz,
                     SzT = SzT,
                     Nsz = Nsz,
                     Xsz = Xsz

                     # hard-coded variables
                     Ntribs <- 2, # only used for feather-butte part, never changes
                     RmultSac = 0.43, # distance for prediction model corresponds to the average distance from all release locations to Woodson Bridge
                     rch_covind = c(1, 2, 3, 3), # index pointing to covariate effect for each reach (note Butte-Sac and Sac-Delta have same fixed effect index)
                     )

  inits <- list(inits, inits, inits)


}

# Get data list
#' @details TODO
#' @returns TODO
#' @export
#' @md
get_survival_data_list <- function(version, sac_data, feather_butte_data) {

  # some variables are the same for both sacramento version / feather butte version
  all_study_years <- tibble(year = sort(unique(c(sac_data$year, feather_butte_data$year)))) |>  # all sac and tribs years combined
    mutate(year_index = row_number())

  # now we set the data frame to pull from sac-specific or feather/butte-specific
  if(version == "sacramento") {
    data <- sac_data
    n_reaches <- 4 # 1 = release-woodson, 2 = woodson-butte, 3 = butte-sac, 4 = sac-delta
    n_detection_locations = 5

    Rmult <- data |>
      select(dist_rlwoodson.z, dist_woodsonbutte.z, dist_buttesac.z, dist_sacdelta.z) |>
      as.matrix() |>
      unname()

    capture_history <- data |>
      separate_wider_position(ch, widths = c("1" = 1, "2" = 1, "3" = 1, "4" = 1, "5" = 1)) |>
      select(`1`, `2`, `3`, `4`, `5`) |>
      mutate_all(as.numeric) |>
      as.matrix() |>
      unname()

    # specific to feather/butte
    trib_ind = NULL
    trib_rg = NULL

  } else {
    # Upper Butte 2019 does not have weight information so use average values across all release group instead.
    data <- feather_butte_data |>
      mutate(across(fish_weight, ~ replace_na(., mean(., na.rm=TRUE))),
             across(fish_k, ~ replace_na(., mean(., na.rm=TRUE))))
    n_tribs <- 2
    n_reaches <- 2
    n_detection_locations <- 3
    trib_ind <- data$trib_ind

    Rmult <- data$dist_rlsac.z
    RmultTrib <- c(1.17, 0.92) # distances for prediction model correspond to the average distance from all release locations for Butte and Feather respectively

    capture_history <- data |>
      separate_wider_position(ch, widths = c("1" = 1, "2" = 1, "3" = 1)) |>
      select(`1`, `2`, `3`) |>
      mutate_all(as.numeric) |>
      as.matrix() |>
      unname()
  }

  year_index <- data |>
    left_join(all_study_years) |>
    pull(year_index)

  release_group_index <- data |>
    group_by(StudyID) |>
    mutate(index = cur_group_id()) |>
    ungroup() |>
    pull(index)

  # TODO clean this up...I can't handle this one right now
  if(version == "feather_butte") {
    trib_rg <- vector(length = length(unique(feather_butte_data$StudyID))) # the tributary index for each release groups
    for(irg in 1:length(unique(feather_butte_data$StudyID))){
      irecs <- which(release_group_index == irg) #identify all records with the current release group index irg
      trib_rg[irg] <- feather_butte_data$trib_ind[irecs[1]] #Get the tributary index. Only need first records as all individuals with same irg will be from same trib
    }
  }

  release_group_water_year_3_index <- data |>
    arrange(year) |> # TODO confirm
    group_by(StudyID) |>
    summarise(ind = unique(WY3)) |>
    ungroup() |>
    pull(ind)

  release_group_water_year_2_index <- data |>
    arrange(year) |> # TODO confirm
    group_by(StudyID) |>
    summarise(ind = unique(WY2)) |>
    ungroup() |>
    pull(ind)

  fork_length <- data$fish_length
  weight <- data$fish_weight
  condition <- data$fish_k

  # prep the full list
  data_list <- list(n_ind = nrow(data),
                    n_reaches = n_reaches,
                    n_detection_locations = n_detection_locations,
                    n_size_classes = 25, # of size classes to plot size effect over
                    n_years = length(all_study_years), # all sac and tribs years combined
                    release_group = unique(data$StudyID),
                    n_release_groups = length(unique(data$StudyID)),
                    first_capture = data$firstCap,
                    last_capture = data$lastCap,
                    water_year_2 = data$WY2,
                    water_year_3 = data$WY3,
                    capture_history = capture_history,
                    Rmult = Rmult, # reach distances standardized per 100km
                    year_index = year_index,
                    release_group_index = release_group_index,
                    release_group_water_year_3_index = release_group_water_year_3_index,
                    release_group_water_year_2_index = release_group_water_year_2_index,
                    fork_length = fork_length,
                    weight = weight,
                    condition = condition,
                    trib_ind = trib_ind,
                    trib_rg = trib_rg)
  return(data_list)
}


#' Call Survival Model
#' @details TODO
#' @returns TODO
#' @export
#' @md
run_survival_model <- function(use_covariate) {
  # function(survival_model_data, number_detection_locations, number_reaches) {

  # TODO confirm CovWY3_Reach model is the best fitting model (the one we want to use)
  model_name <- "survival_model_STAN" # "CovWY3_Reach" reach-specific 3 water year type covariate effect effect on survival (C, D/BN, W)
  parameters <- c("P_b", "muPb", "sdPb", "S_bReach", "S_bCov", "RE_sd",
                  "S_RE", "pred_surv", "SurvWoodSac", "SurvForecast", "pred_pcap")

  number_years <- length(unique(survival_model_data$year))
  number_release_groups <- length(unique(survival_model_data$study_id))

  initial_parameter_values <- list(P_b = matrix(data = 2.2, nrow = number_years, ncol = number_reaches),
                                   muPb = rep(0, number_reaches),
                                   sdPb = rep(1.0, number_reaches),
                                   S_bCov = matrix(data = 1, nrow = 2, number_reaches - 1),
                                   S_bReach = rep(0, number_reaches - 1),
                                   S_RE = rep(-3, number_release_groups),
                                   RE_sd = 0.5)

  # build out model data
  standard_reach_length <- 100 #Survival calculated for a standardized reach length of 100 km, then converted to reach specific survival in stan models
  reach_lengths <- c(40, 88, 170, 110)
  Rmult <- reach_lengths / standard_reach_length

  prepare_data_variables <- survival_model_data |>
    mutate(firstCap = unname(str_locate(ch, "1")[1, 1]), # first station individual was detected at
           WY3 = case_when(year %in% c(2015, 2021) ~ 0,
                           year %in% c(2013, 2016, 2018, 2020) ~ 1,
                           TRUE ~ 2)) |> # 3 water year type groupings (C, D-BN, W) to used as a fixed effect) |>
    group_by(year) |>
    mutate(yrind = cur_group_id()) |> # year id index for for loop in STAN code
    ungroup() |>
    group_by(study_id) |>
    mutate(rgind = cur_group_id()) |> # study id index for for loop in STAN code
    ungroup()

  # index for S_bCov for each release group - year (0=critical, 1=dry/BN, 2 = wet)
  rgwy_ind <- prepare_data_variables |>
    distinct(study_id, WY3) |>
    mutate(WY3_index = WY3 ) |>
    pull(WY3_index)

  if(length(rgwy_ind) != number_release_groups) {
    cli::cli_alert_danger("rgwy_ind is longer than number_release_groups")
    stop()
  }

  # get detection location of last capture
  lastCap <- sapply(survival_model_data$ch, function(x) {
    max(unlist(str_locate_all(x, "1")))
  })

  # convert to a matrix
  CH <- survival_model_data |>
    separate_wider_position(ch, widths = c("1" = 1, "2" = 1, "3" = 1, "4" = 1, "5" = 1)) |>
    select(`1`, `2`, `3`, `4`, `5`) |>
    mutate_all(as.numeric) |>
    as.matrix() |>
    unname()

  if(dim(CH)[2] != number_detection_locations) {
    cli::cli_abort("capture history matrix must have the dimensions of number individuals x number detection locations")
  }

  model_data <- list(Nind = length(unique(survival_model_data$fish_id)),
                     Nreaches = number_reaches, # for right now all capture histories are length = 4 but can change in future
                     Ndetlocs = number_detection_locations, # for right now using 5 detection locations
                     Rmult = Rmult,
                     Nyrs = number_years,
                     Nrg = number_release_groups,
                     CH = CH,
                     yrind = prepare_data_variables$yrind,
                     rgind = prepare_data_variables$rgind,
                     rch_covind = c(1, 2, 3, 3), # index pointing to covariate effect for each reach (note Butte-Sac and Sac-Delta have same fixed effect index),
                     CovX = prepare_data_variables$WY3, # best fitting model (for now)
                     firstCap = prepare_data_variables$firstCap,
                     lastCap = unname(lastCap),
                     rgwy_ind = rgwy_ind)

  inits <- list(initial_parameter_values, initial_parameter_values, initial_parameter_values)

  cli::cli_process_start("Fitting STAN survival model")
  stan_model <- eval(parse(text = "SRJPEmodel::survival_model_code"))
  fit <- rstan::stan(model_code = stan_model,
                     model_name = model_name,
                     data = model_data,
                     init = inits, chains = 3, iter = 1500, include = T,
                     pars = parameters)
  cli::cli_process_done("STAN survival model fitting complete")

  return("full_object" = fit)
}

