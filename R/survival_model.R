
#' Prepare inputs for survival model
#' @details This function prepares the data, initial values, and parameter list needed to call the
#' survival STAN model. It will prepare these data based on what model you want to run.
#' @param number_of_water_year_types either `2`, `3`, or leave empty. If empty, the function will call a
#' version of the survival with no covariate effect.
#' @param effect one of `no_biological_effect`, `fork_length_effect`, `weight_effect`, or `condition`.
#' If supplying an effect, you must supply an argument for `number_of_water_year_types`
#' @returns a named list containing the inputs to the `run_survival_model()` function:
#' * **data_inputs** a named list of data inputs for passing to the STAN model.
#' * **inits** a list of initial values for passing to the STAN model.
#' * **parameters** a vector of parameter names to save from calling the STAN model.
#' @export
#' @md
prepare_survival_inputs <- function(number_of_water_year_types = NULL,
                                    effect = NULL){

  # check args
  if(is.null(number_of_water_year_types) & !is.null(effect)) {
    cli::cli_abort("You must specify a number of water year types if you are supplying a biological effect.")
  }

  null_effect = ifelse(is.null(effect), TRUE, FALSE)

  if(is.null(effect)) {
    effect <- "no_biological_effect"
    cli::cli_bullets("No biological effect was supplied, running the no biological effect model")
  }

  if(!effect %in% c("no_biological_effect", "fork_length", "weight", "condition")) {
    cli::cli_abort("Effect argument must be one of no_biological_effect, fork_length, weight, or condition")
  }

  sac_data <- SRJPEdata::survival_model_inputs |>
    filter(is.na(trib_ind))

  fb_data <- SRJPEdata::survival_model_inputs |>
    filter(!is.na(trib_ind))

  sac_data_list <- get_survival_data_list("sacramento", sac_data, fb_data)
  fb_data_list <- get_survival_data_list("feather_butte", sac_data, fb_data)

  if(is.null(number_of_water_year_types) & null_effect) {

    parameters <- c("P_b", "muPb", "sdPb", "S_bReach", "S_bTrib", "RE_sd",
                    "RE_sdT", "S_RE", "S_REt", "pred_surv", "SurvRelSac", "SurvWoodSac",
                    "SurvForecast", "pred_pcap", "pred_survT", "TribSurvForecast")


    inits <- list(P_b = matrix(data = 2.2,
                               nrow = sac_data_list$n_years,
                               ncol = sac_data_list$n_reaches),
                  muPb = rep(0, sac_data_list$n_reaches),
                  sdPb = rep(1.0, sac_data_list$n_reaches),
                  S_bReach = rep(0, 3),
                  S_RE = rep(-3, sac_data_list$n_release_groups),
                  RE_sd = 0.5)

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

    inits <- list(P_b = matrix(data = 2.2,
                               nrow = sac_data_list$n_years,
                               ncol = sac_data_list$n_reaches),
                  muPb = rep(0, sac_data_list$n_reaches),
                  sdPb = rep(1.0, sac_data_list$n_reaches),
                  S_bReach = rep(0, 3),
                  S_RE = rep(-3, sac_data_list$n_release_groups),
                  RE_sd = 0.5)

    # now step through all possible covariate versions - water year type (either 2 or 3) and
    # effects (no effect, fork length, weight, and condition)

    # set water year type-specific variables
    if(number_of_water_year_types == 2) {
      NS_bCov <- 1
      CovX <- sac_data_list$water_year_2
      CovXT <- fb_data_list$water_year_2
      rgwy_ind <- sac_data_list$release_group_water_year_2_index
      rgwy_indT <- fb_data_list$release_group_water_year_2_index
    } else if (number_of_water_year_types == 3) {
      NS_bCov <- 2
      CovX <- sac_data_list$water_year_3
      CovXT <- fb_data_list$water_year_3
      rgwy_ind <- sac_data_list$release_group_water_year_3_index
      rgwy_indT <- fb_data_list$release_group_water_year_3_index
    }

    # set effect-specific variables

    # if no effect supplied, run no_biological_effect

    if (effect == "no_biological_effect") {

      Sz <- rep(0, sac_data_list$n_ind)
      SzT <- rep(0, fb_data_list$n_ind)
      Xsz <- rep(0, sac_data_list$n_size_classes)

    } else if (effect == "fork_length") {

      all_fork_lengths <- c(sac_data_list$fork_length, fb_data_list$fork_length)
      mean_fork_length <- mean(all_fork_lengths, na.rm = TRUE)
      sd_fork_length <- sd(all_fork_lengths, na.rm = TRUE)

      Sz <- (sac_data_list$fork_length - mean_fork_length) / sd_fork_length
      SzT <- (fb_data_list$fork_length - mean_fork_length) / sd_fork_length
      Xsz  <- seq(from = min(c(Sz, SzT), na.rm = TRUE),
                  to = max(c(Sz, SzT), na.rm = TRUE),
                  length.out = sac_data_list$n_size_classes)

    } else if (effect == "weight") {

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
                     firstCap = sac_data_list$first_capture,
                     lastCap = sac_data_list$last_capture,

                     # variables set based on number_water_year_types
                     rgwy_ind = rgwy_ind,
                     rgwy_indT = rgwy_indT,
                     CovX = CovX,
                     CovXT = CovXT,

                     # variables set based on effect
                     UseSizeEffect = ifelse(effect == "no_biological_effect", 0, 1),
                     NS_bCov = NS_bCov,
                     Sz = Sz,
                     SzT = SzT,
                     Xsz = Xsz,

                     # feather/butte-specific variables
                     NindT = fb_data_list$n_ind,
                     NrgT = fb_data_list$n_release_groups,
                     Ntribs = fb_data_list$n_tribs,
                     trib_ind = fb_data_list$trib_ind,
                     firstCapT = fb_data_list$first_capture,
                     lastCapT = fb_data_list$last_capture,
                     RmultT = fb_data_list$Rmult,
                     yrindT = fb_data_list$year_index,
                     rgindT = fb_data_list$release_group_index,
                     CHT = fb_data_list$capture_history,
                     trib_rg = fb_data_list$trib_rg,

                     # hard-coded variables
                     Nsz = 25, # number of size classes to plot size effect over
                     RmultTrib = c(1.17, 0.92), # distances for prediction model correspond to the average distance from all release locations for Butte and Feather respectively
                     Ntribs <- 2, # only used for feather-butte part, never changes
                     RmultSac = 0.43, # distance for prediction model corresponds to the average distance from all release locations to Woodson Bridge
                     rch_covind = c(1, 2, 3, 3) # index pointing to covariate effect for each reach (note Butte-Sac and Sac-Delta have same fixed effect index)
                     )

  inits <- list(inits, inits, inits)

  return(list("data_inputs" = model_data,
              "inits" = inits,
              "parameters" = parameters))


}

# Get data list
#' @details This prepares data within the `prepare_survival_inputs()` function for use in the survival
#' STAN model based on either the Sacramento dataset or Feather/Butte dataset.
#' @returns A named list of some of the data inputs required for the STAN model.
#' @export
#' @md
get_survival_data_list <- function(version, sac_data, feather_butte_data) {

  # some variables are the same for both sacramento version / feather butte version
  all_study_years <- tibble(year = sort(unique(c(sac_data$year, feather_butte_data$year)))) |>  # all sac and tribs years combined
    mutate(year_index = row_number())

  # now we set the data frame to pull from sac-specific or feather/butte-specific
  if(version == "sacramento") {
    data <- sac_data

    # hard-coded value, update if code updates
    n_reaches <- 4 # 1 = release-woodson, 2 = woodson-butte, 3 = butte-sac, 4 = sac-delta
    n_detection_locations = 5

    Rmult <- data |>
      select(dist_rlwoodson_standardized, dist_woodsonbutte_standardized,
             dist_buttesac_standardized, dist_sacdelta_standardized) |>
      as.matrix() |>
      unname()

    capture_history <- data |>
      separate_wider_position(ch, widths = c("1" = 1, "2" = 1, "3" = 1, "4" = 1, "5" = 1)) |>
      select(`1`, `2`, `3`, `4`, `5`) |>
      mutate_all(as.numeric) |>
      as.matrix() |>
      unname()

    # specific to feather/butte
    n_tribs <- NULL
    trib_ind <- NULL
    trib_index_for_release_groups <- NULL

  } else {
    # Upper Butte 2019 does not have weight information so use average values across all release group instead.
    data <- feather_butte_data |>
      mutate(across(fish_weight, ~ replace_na(., mean(., na.rm=TRUE))),
             across(fish_k, ~ replace_na(., mean(., na.rm=TRUE))))

    # hard-coded value, update if code updates
    n_tribs <- 2
    n_reaches <- 2
    n_detection_locations <- 3
    trib_ind <- data$trib_ind

    Rmult <- data$dist_rlsac_standardized

    capture_history <- data |>
      separate_wider_position(ch, widths = c("1" = 1, "2" = 1, "3" = 1)) |>
      select(`1`, `2`, `3`) |>
      mutate_all(as.numeric) |>
      as.matrix() |>
      unname()

    trib_index_for_release_groups <- data |>
      arrange(year) |>
      distinct(study_id, trib_ind) |>
      pull(trib_ind)
  }

  year_index <- data |>
    left_join(all_study_years, by = "year") |>
    pull(year_index)

  release_group_index <- data |>
    arrange(year) |>
    group_by(study_id) |>
    mutate(index = cur_group_id()) |>
    ungroup() |>
    pull(index)

  release_group_water_year_3_index <- data |>
    arrange(year) |>
    group_by(study_id) |>
    summarise(ind = unique(wy_three_categories)) |>
    ungroup() |>
    pull(ind)

  release_group_water_year_2_index <- data |>
    arrange(year) |>
    group_by(study_id) |>
    summarise(ind = unique(wy_two_categories)) |>
    ungroup() |>
    pull(ind)

  fork_length <- data$fish_length
  weight <- data$fish_weight
  condition <- data$fish_k

  # prep the full list
  data_list <- list(n_ind = nrow(data),
                    n_tribs = n_tribs,
                    n_reaches = n_reaches,
                    n_detection_locations = n_detection_locations,
                    n_size_classes = 25, # of size classes to plot size effect over
                    n_years = nrow(all_study_years), # all sac and tribs years combined
                    release_group = unique(data$study_id),
                    n_release_groups = length(unique(data$study_id)),
                    first_capture = data$first_capture,
                    last_capture = data$last_capture,
                    water_year_2 = data$wy_two_categories,
                    water_year_3 = data$wy_three_categories,
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
                    trib_rg = trib_index_for_release_groups)
  return(data_list)
}


#' Fit Survival Model
#' @details This model calls a survival STAN model that estimates survival in the Sacramento and Feather/Butte systems.
#' @param survival_inputs The named list object generated by running `prepare_survival_inputs()`.
#' @returns a STANfit object with the survival model fit.
#' @export
#' @md
fit_survival_model <- function(survival_inputs) {

  # check which model we're running by checking NS_bCov, which is 0 for the no covariate model
  if(survival_inputs$data_inputs$NS_bCov == 0) {
    stan_model <- eval(parse(text = "SRJPEmodel::survival_model_code$survival_NoCov"))
    cli::cli_bullets("Running NoCov survival model")
  } else {
    stan_model <- eval(parse(text = "SRJPEmodel::survival_model_code$survival_CovWY"))
    cli::cli_bullets("Running CovWY survival model")
  }

  options(mc.cores = parallel::detectCores())

  cli::cli_process_start("Fitting STAN survival model")
  fit <- rstan::stan(model_code = stan_model,
                     data = survival_inputs$data_inputs,
                     init = survival_inputs$inits,
                     chains = 3,
                     iter = 2000,
                     include = T,
                     pars = survival_inputs$parameters,
                     seed=1234)
  cli::cli_process_done("STAN survival model fitting complete")

  return("fit" = fit)
}

#' Extract parameter estimates from the survival STAN model object
#' @details This function extracts parameter estimates from the survival STAN model and returns them in a tidy data frame format.
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
    distinct(release_group_sac = study_id, release_group_index_sac)

  release_group_lookup_trib <- SRJPEdata::survival_model_inputs |>
    filter(!is.na(trib_ind)) |>
    arrange(year) |>
    group_by(study_id) |>
    mutate(release_group_index_trib = cur_group_id()) |>
    ungroup() |>
    distinct(release_group_trib = study_id, release_group_index_trib)

  # get release group lookup

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
         water_year_type = ifelse(parameter %in% c("SurvRelSacSz", "SurvWoodSacSz", "SurvForecastSz",
                                                   "pred_SurvTSz"), index_2, NA),
         environmental_covarariate_dimension = ifelse(parameter %in% c("S_bCov", "s_bCovT"), index_1, NA),
         forecast_index = ifelse(parameter == "SurvForecast", index_1, NA),
         # reach-structured parameters
         sac_reach_index = case_when(parameter == "P_b" ~ index_2,
                                     parameter %in% c("muPb", "sdPb", "S_bReach") ~ index_1,
                                     parameter %in% c("pred_surv", "pred_pcap") ~ index_2,
                                     TRUE ~ NA),
         # trib-structured parameters
         trib_reach_index = case_when(parameter == "S_bTrib" ~ index_1,
                                      parameter == "pred_survT" ~ index_2, # TODO confirm this, it's hard coded as c(1, 2), it's the trib model reaches?
                                      parameter %in% c("pred_survTSz", "TribSurvForecastSz", "pred_SurvTSz") ~ index_3,
                                      TRUE ~ NA),
         # release group structured parameters
         release_group_index_sac = ifelse(parameter %in% c("SurvRelSac","SurvWoodSac", "pred_surv",
                                                           "S_RE"), index_1, NA),
         release_group_index_trib = ifelse(parameter %in% c("S_REt", "pred_survT"), index_1, NA),
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
                              parameter)) |>
    select(-c(year_index, sac_reach_index, trib_reach_index, release_group_index_sac,
              release_group_index_trib, reach_name_sac, reach_name_trib, release_group_sac,
              release_group_trib, pred_forecast_index_trib, index_1, index_2, index_3)) |>
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
    facet_wrap(~Reach) +
    labs(x = "Release Group",
         y = "Median Survival Rate") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1))

}
