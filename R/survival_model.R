# refactor of Call_Model.R
# library(rstan)
# options(mc.cores = parallel::detectCores())
# rstan_options(auto_write = TRUE)

#' Prepare inputs for survival model
#' @details Runs the survival model for both Sacramento releases and Feather/Butte releases.
#' @returns TODO
#' @export
#' @md
prepare_survival_inputs <- function(){

  data <- SRJPEdata::survival_model_inputs |>
    glimpse()

  # create a sorted sac
  # created a sorted

}

#' Call Survival Model
#' @details TODO
#' @param use_covariate either `TRUE` or `FALSE`. If `TRUE`, will call a version of the survival
#' model that fits a covariate effect (water year type). If `FALSE`, it will call a version of the survival model
#' with no covariate effect.
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

