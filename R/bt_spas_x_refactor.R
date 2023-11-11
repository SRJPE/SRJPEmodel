# refactor of josh_original_model_code.R and MultiRun_SpecialPriors.R

run_bt_spas_x <- function(number_mcmc, number_burnin, number_thin, number_chains,
                          bt_spas_x_input_data, effort_adjust, multi_run_mode,
                          special_priors_data) {
  if(!multi_run_mode) {
    # TODO one trib, site, and year
  } else if (multi_run_mode) {
    # TODO loop through tribs, sites, and years

  }
}

run_single_bt_spas_x <- function(number_mcmc, number_burnin, number_thin, number_chains,
                                 bt_spas_x_input_data, site, run_year,
                                 effort_adjust = c(T, F), trib = c(T, F), special_priors_data,
                                 bugs_directory) {

  # TODO cache bayesian specs as params?

  # TODO get rid of this
  load("~/Desktop/weekly_model_data.rda")
  special_priors_data <- read.csv(here::here("data-raw", "first_draft", "Special_Priors.csv")) |>
    mutate(site = ifelse(Stream_Site == "battle creek_ubc", "ubc", NA)) |>
    select(site, run_year = RunYr, week = Jweek, special_prior = lgN_max)
  bt_spas_x_input_data <- weekly_model_data |> ungroup()
  trib = T
  site_selection <- "ubc"
  run_year_selection <- 2009

  # set up filter - if it's a tributary-based model, we cannot use efficiencies from KDL, TIS, RBDD
  # TODO fill in logic for when it is a mainstem-based model
  if(trib) {
    remove_sites <- c("knights landing", "tisdale", "red bluff diversion dam")
  } else {
    remove_sites <- NA
  }

  # filter datasets to match site, run_year, and weeks
  # catch_flow is average for julian week, mr_flow is average over recapture days (< 1 week)
  input_data <- bt_spas_x_input_data |>
    filter(run_year == run_year_selection,
           site == site_selection)
    # filter(run_year == !!run_year,
    #        site == !!site)

  # analyze efficiency trials for all relevant sites (do not filter to site)
  # TODO remove - moving to SRJPEdata repo as mainstem and trib efficiency versions
  mark_recapture_data <- bt_spas_x_input_data |>
    filter(!site %in% remove_sites,
           !is.na(flow_cfs), # TODO change to efficiency flow when available
           !is.na(number_released),
           !is.na(number_recaptured)) |>
    group_by(stream) |>
    mutate(standardized_flow = (flow_cfs - mean(flow_cfs, na.rm = T)) / sd(flow_cfs, na.rm = T)) # TODO change first instance of flow to efficiency_flow once available

  number_efficiency_experiments <- nrow(mark_recapture_data)
  years_with_efficiency_experiments <- unique(mark_recapture_data$run_year)


  # TODO imagine we are pulling in efficiency data based on argument "trib"

  # prep data for abundance component of pCap model
  # effort adjustment
  # TODO if effort_adjust == T, use standardized_catch_effort, otherwise catch

  # priors for upper limit on log abundance for any week
  # TODO keep this here?
  # first, assign special prior (if relevant), else set to default, then fill in for weeks without catch
  data_with_priors <- input_data |>
    left_join(special_priors_data, by = c("run_year", "week", "site")) |>
    mutate(lgN_prior = ifelse(!is.na(special_prior), special_prior, log((count / 1000) + 1) / 0.025)) |> # maximum possible value for log N across strata
    select(-special_prior)

  # TODO use standardized catch flow

  #Identify the elements in 1:Nstrata (unmarked catch set) without and with pCap and corresponding flow data
  weeks_with_recap <- data_with_priors |>
    filter(!is.na(number_released),
           !is.na(standardized_efficiency_flow))
  weeks_without_recap <- data_with_priors |>
    filter(is.na(number_released),
           is.na(standardized_efficiency_flow))
  # TODO josh has something for BUGS-specific code here
  # if(nrow(weeks_without_recap) == 1 | nrow(weeks_with_recap) == 1)

  # write _data.out file
  # TODO do we need this?


  # set up b-spline basis matrix
  k_int <- 4 # rule of thumb is 1 knot for every 4 data points for a cubic spline (which has 4 parameters)
  number_knots <- round(nrow(input_data) / k_int, 0)
  first_knot_position <- 2
  final_knot_position <- nrow(input_data) - 1 # keep first and/or last knot positions away from tails if there are intervals with no sampling on the tails
  knot_positions <- seq(first_knot_position, final_knot_position, length.out = number_knots) # define position of b-spline knots using even interval if no missing data
  b_spline_matrix <- splines2::bSpline(x = 1:nrow(input_data), knots = knot_positions, deg = 3, intercept = T) # bspline basis matrix. One row for each data point (1:Nstrata), and one column for each term in the cubic polynomial function (4) + number of knots
  ncol_b_spline_matrix <- ncol(b_spline_matrix) # TODO fix naming

  # prep objects for passing to BUGS
  number_experiments_at_site <- mark_recapture_data |>
    group_by(site) |>
    tally() |>
    filter(site == site_selection) |>
    pull(n)

  # set up models and data to pass to bugs
  # if efficiency trials occurred in the site
  if(number_experiments_at_site > 1) {

    if(nrow(weeks_without_recap) == 0) { # all weeks have efficiency trials
      data <- list("Nmr", "Ntribs" ,"ind_trib", "Releases", "Recaptures", "Nstrata", "u", "K", "ZP",
                   "ind_pCap", "Nstrata_wc", "Uwc_ind", "mr_flow", "lgN_max")
      model_name <- pasteo("estN_allMR", ".bug")
    } else { # some or all strata don't have efficiency trials

      if(nrow(weeks_with_recap) > 0) { # some weeks have efficiency trials
        data <- list("Nmr", "Ntribs", "ind_trib", "Releases", "Recaptures", "Nstrata", "u", "K", "ZP",
                     "ind_pCap", "Nwmr", "Nwomr", "Uind_wMR", "Uind_woMR", "use_trib",
                     "Nstrata_wc", "Uwc_ind", "mr_flow", "catch_flow", "lgN_max")
        model_name <- paste0("estN_missMR", ".bug")

      } else if(nrow(weeks_with_recap) == 0) { # no weeks have efficiency trials
        data <- list("Nmr", "Ntribs", "ind_trib", "Releases", "Recaptures", "Nstrata", "u", "K", "ZP", "Nwomr",
                     "Uind_woMR", "use_trib", "Nstrata_wc", "Uwc_ind", "mr_flow", "catch_flow", "lgN_max")
        model_name <- paste0("estN_noMR", ".bug")
      }
    }
  } else if(number_experiments_at_site == 0) { # no efficiency trials were performed at that site
    data <- list("Nmr", "Ntribs", "ind_trib", "Releases", "Recaptures", "Nstrata", "u", "K", "ZP", "Nwomr",
                 "Uind_woMR", "Nstrata_wc", "Uwc_ind", "mr_flow", "catch_flow", "lgN_max")
    model_name <- paste0("estN_noMR_notrib", ".bug")

  }

  parameters <- c("trib_mu.P", "trib_sd.P", "flow_mu.P", "flow_sd.P", "pro_sd.P", "b0_pCap", "b_flow",
                  "pCap_U", "N", "Ntot", "sd.N", "sd.Ne")

  # initialize parameters
  ini_b0_pCap <- mark_recapture_data |>
    filter(site == site_selection,
           run_year == run_year_selection) |>
    # TODO confirm logit is the same as qlogis
    mutate(ini_b0_pCap = stats::qlogis(sum(number_recaptured) / sum(number_released))) |>
    pull(ini_b0_pCap)

  # TODO double check these
  ini_lgN <- data_with_priors |>
    mutate(ini_lgN = log(catch_standardized_by_effort / 1000 + 2),
           ini_lgN = ifelse(is.na(ini_lgN), log(2 / 1000), ini_lgN)) |>
    pull(ini_lgN)

  trib_mu.P <- mark_recapture_data |>
    mutate(trib_mu.P = stats::qlogic(sum(number_recaptured) / sum(number_released))) |>
    pull(trib_mu.P)

  init_list <- list(trib_mu.P = trib_mu.P,
                    b0_pCap = ini_b0_pCap,
                    flow_mu.P = 0,
                    b_flow = rep(0, length(unique(mark_recapture_data$site))),# TODO double check this
                    trib_tau.P = 1,
                    flow_tau.P = 1,
                    pro_tau.P = 1,
                    b_sp = rep(1, ncol_b_spline_matrix),
                    lg_N = ini_lgN)

  inits <- list(inits1 = init_list, inits2 = init_list, inits3 = init_list)

  # run bugs model
  model_results <- bugs(data, inits, parameters, model_name, n.chains = number_chains,
                        n.burnin = number_burnin, n.thin = number_thin, n.iter = number_mcmc,
                        debug = FALSE, codaPkg = FALSE, DIC = TRUE, clearWD = TRUE,
                        bugs.directory = paste0(bugs_directory))

  posterior_output <- model_results$sims.list
  summary_output <- round(model_results$summary, 3)
  dic_output <- c(model_results$pD, model_results$DIC)
  knots_output <- knot_positions

  return(list("posterior_output" = posterior_output,
              "summary_output" = summary_output,
              "dic_output" = dic_output,
              "knots_output" = knots_output))

}
