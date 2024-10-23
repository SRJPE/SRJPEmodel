# compare BT SPAS X run on WinBUGS and on STAN
library(tidyverse)
library(SRJPEdata)
library(SRJPEmodel)

# data prep ---------------------------------------------------------------

# do not run for yearling
weekly_catch <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
  filter(life_stage %in% c("fry", "smolt"))

weekly_efficiency <- SRJPEdata::weekly_juvenile_abundance_efficiency_data |>
  filter(number_released > 0) # TODO temporary fix

trials_to_fit <- weekly_catch |>
  distinct(site, run_year, life_stage)

nrow(trials_to_fit)

# function for running multiple -------------------------------------------

run_multiple_bugs <- function(site, run_year, life_stage) {

  cli::cli_bullets(paste("Running WinBUGS on", site, "for run year", run_year, "and life stage", life_stage))

  results <- tryCatch({run_single_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                            weekly_catch,
                                            weekly_efficiency,
                                            site = site,
                                            run_year = run_year,
                                            lifestage = life_stage,
                                            effort_adjust = F,
                                            bugs_directory = here::here("data-raw", "WinBUGS14"),
                                            debug_mode = FALSE,
                                            no_cut = F)},
                      error = function(e) return("error")
                      )
  return(results)
}

run_multiple_stan <- function(site, run_year, life_stage) {

  cli::cli_bullets(paste("Running STAN on", site, "for run year", run_year, "and life stage", life_stage))

  # TODO do we want to add a function for processing these results?
  results <- tryCatch({run_single_bt_spas_x_stan(SRJPEmodel::bt_spas_x_bayes_params,
                                                 weekly_catch,
                                                 weekly_efficiency,
                                                 site = site,
                                                 run_year = run_year,
                                                 lifestage = life_stage,
                                                 effort_adjust = F)},
                      error = function(e) return("error")
  )
  return(results)
}


# run for BUGS ------------------------------------------------------------

bugs_results <- purrr::pmap(list(trials_to_fit$site,
                                 trials_to_fit$run_year,
                                 trials_to_fit$life_stage),
                            run_multiple_bugs)


# run for STAN ------------------------------------------------------------

stan_results <- purrr::pmap(list(trials_to_fit$site,
                                 trials_to_fit$run_year,
                                 trials_to_fit$life_stage),
                            run_multiple_stan)


# results -----------------------------------------------------------------

# extract_results <- function(element) {
#   if(!is.list(element)) {
#     return(list("error" = TRUE))
#   } else {
#     return(element$final_results)
#   }
# }
#
# clean_results <- lapply(battle_results, extract_results) |>
#   bind_rows() |>
#   filter(error %in% c(FALSE, NA)) |>
#   select(-error) |>
#   glimpse()

