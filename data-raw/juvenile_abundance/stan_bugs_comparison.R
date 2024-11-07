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

nrow(trials_to_fit_battle)

trials_to_fit_butte <- trials_to_fit |>
  filter(site == "okie dam")

nrow(trials_to_fit_butte)

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
bugs_results_butte <- purrr::pmap(list(trials_to_fit_battle$site,
                                 trials_to_fit_battle$run_year,
                                 trials_to_fit_battle$life_stage),
                            run_multiple_bugs,
                            .progress = TRUE)
saveRDS(bugs_results_butte, here::here("data-raw", "juvenile_abundance",
                                 "battle_results_BUGS_oct_2024.rds"))

# run for STAN ------------------------------------------------------------

options(mc.cores=parallel::detectCores())

stan_results_battle <- purrr::pmap(list(trials_to_fit_battle$site,
                                        trials_to_fit_battle$run_year,
                                        trials_to_fit_battle$life_stage),
                            run_multiple_stan,
                            .progress = TRUE)
saveRDS(stan_results_battle, here::here("data-raw", "juvenile_abundance",
                                         "battle_results_STAN_nov_2024.rds"))

stan_results_butte <- purrr::pmap(list(trials_to_fit_battle$site,
                                        trials_to_fit_battle$run_year,
                                        trials_to_fit_battle$life_stage),
                                   run_multiple_stan,
                                   .progress = TRUE)
saveRDS(stan_results_butte, here::here("data-raw", "juvenile_abundance",
                                        "battle_results_STAN_nov_2024.rds"))

# results -----------------------------------------------------------------

extract_stan_results <- function(element) {

  if(class(element) == "list") {
    results <- element$summary_table |>
      filter(parameter %in% c("pCap_U", "N", "Ntot"))
  } else {
    results <- data.frame("error" = TRUE)
  }
  return(results)
}

extract_bugs_results <- function(element) {
  par_list <- c("N", "pCap_U", "Ntot")

  if(class(element) == "list") {
    results <- element$final_results |>
      filter(str_detect(parameter, "N\\[") | str_detect(parameter, "pCap_U\\["))
  } else {
    results <- data.frame("error" = TRUE)
  }
  return(results)
}

# bugs <- readRDS("data-raw/juvenile_abundance/battle_creek_results_10-23-2024.rds")
bugs <- readRDS("data-raw/juvenile_abundance/battle_results_BUGS_nov_2024.rds")
stan <- readRDS("data-raw/juvenile_abundance/battle_results_STAN_nov_2024.rds")

battle_stan_results_clean <- lapply(stan,
                                    extract_stan_results) |>
  bind_rows()

battle_bugs_results_clean <- lapply(bugs,
                                    extract_bugs_results) |>
  bind_rows()


# convergence -------------------------------------------------------------

# bugs
lapply(bugs, is.list) |> unlist() |> sum() # 44/69

# stan
lapply(stan, is.list) |> unlist() |> sum() # 37/69


# streams in both ---------------------------------------------------------

bugs_fit <- battle_bugs_results_clean |>
  filter(is.na(error)) |>
  distinct(run_year, life_stage, site) |>
  mutate(id = paste(run_year, life_stage, site, sep = "-"))

stan_fit <- battle_stan_results_clean |>
  filter(is.na(error)) |>
  distinct(run_year, life_stage, site) |>
  mutate(id = paste(run_year, life_stage, site, sep = "-"))

# fit in bugs but not stan
bugs_fit$id[!bugs_fit$id %in% stan_fit$id]
# fit in stan but not bugs
stan_fit$id[!stan_fit$id %in% bugs_fit$id]
# fit in both
both_fit <- bugs_fit$id[stan_fit$id %in% bugs_fit$id]

# plot --------------------------------------------------------------------

# and add a 1:1 line
# abundance
battle_bugs_results_clean |>
  filter(is.na(error)) |>
  mutate(parameter = gsub("[0-9]+|\\[|\\]", "", parameter),
         model = "BUGS") |>
  bind_rows(battle_stan_results_clean |>
              filter(is.na(error)) |>
            mutate(model = "STAN",
                   statistic = str_remove_all(statistic, "\\%"))) |>
  mutate(id = paste(run_year, life_stage, site, sep = "-")) |>
  filter(parameter == "Ntot",
         statistic == "50",
         id %in% both_fit) |>
  pivot_wider(names_from = "model",
              values_from = "value") |>
  mutate(diff = STAN-BUGS) |>
  ggplot(aes(x = week_fit, y = diff)) +
  geom_hline(aes(yintercept = 0), color = "red", linetype = "dashed") +
  # ggplot(aes(x = BUGS, y = STAN)) +
  #ggplot(aes(x = week_fit, y = value, color = model)) +
  #geom_abline(intercept = 0, slope = 1, color = "red", linetype = 2) +
  geom_point() +
  facet_wrap(~id, scales = "free_y")

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



# try again with a subset -------------------------------------------------


# run for BUGS ------------------------------------------------------------
ubc_trials <- trials_to_fit_battle |>
  filter(site == "lbc") |>
  glimpse()

bugs_results <- purrr::pmap(list(trials_to_fit_battle$site,
                                 trials_to_fit_battle$run_year,
                                 trials_to_fit_battle$life_stage),
                            run_multiple_bugs,
                            .progress = TRUE)
saveRDS(bugs_results, here::here("data-raw", "juvenile_abundance",
                                 "battle_results_BUGS_nov_2024.rds"))


# run for STAN ------------------------------------------------------------

options(mc.cores=parallel::detectCores())

stan_results_battle <- purrr::pmap(list(trials_to_fit_battle$site,
                                        trials_to_fit_battle$run_year,
                                        trials_to_fit_battle$life_stage),
                                   run_multiple_stan,
                                   .progress = TRUE)
saveRDS(stan_results_battle, here::here("data-raw", "juvenile_abundance",
                                        "battle_results_STAN_nov_2024.rds"))

