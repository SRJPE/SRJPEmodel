# compare BT SPAS X run on WinBUGS and on STAN
library(tidyverse)
library(SRJPEdata)
library(SRJPEmodel)

# data prep ---------------------------------------------------------------

# do not run for yearling
weekly_catch <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
  filter(life_stage %in% c("fry", "smolt"))

weekly_efficiency <- SRJPEdata::weekly_juvenile_abundance_efficiency_data

trials_to_fit <- weekly_catch |>
  distinct(site, run_year, life_stage)

nrow(trials_to_fit)

trials_to_fit_battle <- trials_to_fit |>
  filter(site %in% c("ubc", "lbc"))

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

bugs_results <- purrr::pmap(list(trials_to_fit_battle$site,
                                 trials_to_fit_battle$run_year,
                                 trials_to_fit_battle$life_stage),
                            run_multiple_bugs,
                            .progress = TRUE)
saveRDS(bugs_results, here::here("data-raw", "juvenile_abundance",
                                 "battle_results_BUGS_oct_2024.rds"))


# run for STAN ------------------------------------------------------------

options(mc.cores=parallel::detectCores())

stan_results_battle <- purrr::pmap(list(trials_to_fit_battle$site[1:2],
                                        trials_to_fit_battle$run_year[1:2],
                                        trials_to_fit_battle$life_stage[1:2]),
                            run_multiple_stan,
                            .progress = TRUE)
saveRDS(stan_results_battle, here::here("data-raw", "juvenile_abundance",
                                         "battle_results_STAN_oct_2024.rds"))


# results -----------------------------------------------------------------

extract_stan_results <- function(element) {

  if(class(element) == "list") {
    results <- element$summary_table |>
      filter(parameter %in% c("pCap_U", "N"))
    if(nrow(results) == 0) {
      results <- data.frame("error" = TRUE)
    }
    return(results)
  }
}

extract_bugs_results <- function(element) {
  par_list <- c("N", "pCap_U")

  if(class(element) == "list") {
    results <- element$final_results |>
      filter(str_detect(parameter, "N\\[") | str_detect(parameter, "pCap_U\\["))
  } else {
    results <- data.frame("error" = TRUE)
  }
}

bugs <- readRDS("data-raw/juvenile_abundance/battle_creek_results_10-23-2024.rds")
stan <- readRDS("data-raw/juvenile_abundance/battle_results_STAN_oct_2024.rds")

battle_stan_results_clean <- lapply(stan,
                                    extract_stan_results) |>
  bind_rows()

battle_bugs_results_clean <- lapply(bugs,
                                    extract_bugs_results) |>
  bind_rows()


# convergence -------------------------------------------------------------

# bugs
battle_bugs_results_clean |>
  filter(is.na(error)) |>
  distinct(site, run_year, life_stage) |>
  tally()
battle_bugs_results_clean |>
  filter(error) |>
  tally()

# converged
15/25

# stan




# plot --------------------------------------------------------------------

# TODO run stan to get clean results and then plot against each other
# and add a 1:1 line
# abundance
battle_bugs_results_clean |>
  filter(is.na(error)) |>
  mutate(parameter = gsub("[0-9]+|\\[|\\]", "", parameter)) |>
  filter(parameter == "N") |>
  pivot_wider(names_from = "statistic", values_from = "value") |>
  ggplot(aes(x = week_fit, y = `50`)) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = 2) +
  geom_point() +
  facet_wrap(~site + run_year + life_stage)

# efficiency
battle_bugs_results_clean |>
  filter(is.na(error)) |>
  mutate(parameter = gsub("[0-9]+|\\[|\\]", "", parameter)) |>
  filter(parameter == "pCap_U") |>
  pivot_wider(names_from = "statistic", values_from = "value") |>
  ggplot(aes(x = week_fit, y = `50`)) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = 2) +
  geom_point() +
  facet_wrap(~site + run_year + life_stage)





