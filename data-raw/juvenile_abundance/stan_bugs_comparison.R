# compare BT SPAS X run on WinBUGS and on STAN
library(tidyverse)
library(SRJPEdata)
library(SRJPEmodel)
library(scales)

# data prep ---------------------------------------------------------------

# do not run for yearling
weekly_catch <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
  filter(life_stage %in% c("fry", "smolt"))

SRJPEdata::weekly_juvenile_abundance_efficiency_data |>
  filter(number_released == 0)

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
                                            effort_adjust = T,
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
                                                 effort_adjust = T)},
                      error = function(e) return("error")

  )
  return(results)
}


# run for BUGS ------------------------------------------------------------
# battle
bugs_results <- purrr::pmap(list(trials_to_fit_battle$site,
                                 trials_to_fit_battle$run_year,
                                 trials_to_fit_battle$life_stage),
                            run_multiple_bugs,
                            .progress = TRUE)
saveRDS(bugs_results, here::here("data-raw", "juvenile_abundance",
                                 "battle_results_BUGS_oct_2024.rds"))
# butte
bugs_results_butte <- purrr::pmap(list(trials_to_fit_battle$site,
                                 trials_to_fit_battle$run_year,
                                 trials_to_fit_battle$life_stage),
                            run_multiple_bugs,
                            .progress = TRUE)
saveRDS(bugs_results_butte, here::here("data-raw", "juvenile_abundance",
                                 "battle_results_BUGS_oct_2024.rds"))

# run for STAN ------------------------------------------------------------

options(mc.cores=parallel::detectCores())
# battle
stan_results_battle <- purrr::pmap(list(trials_to_fit_battle$site,
                                        trials_to_fit_battle$run_year,
                                        trials_to_fit_battle$life_stage),
                            run_multiple_stan,
                            .progress = TRUE)
saveRDS(stan_results_battle, here::here("data-raw", "juvenile_abundance",
                                         "battle_results_STAN_nov_2024.rds"))
# butte
stan_results_butte <- purrr::pmap(list(trials_to_fit_butte$site,
                                        trials_to_fit_butte$run_year,
                                        trials_to_fit_butte$life_stage),
                                   run_multiple_stan,
                                   .progress = TRUE)
saveRDS(stan_results_butte, here::here("data-raw", "juvenile_abundance",
                                        "butte_results_STAN_nov_2024.rds"))

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

bugs_battle <- readRDS("data-raw/juvenile_abundance/battle_results_BUGS_nov_2024.rds")
bugs_butte <- readRDS("data-raw/juvenile_abundance/butte_results_BUGS_nov_2024.rds")
stan_battle <- readRDS("data-raw/juvenile_abundance/battle_results_STAN_nov_2024.rds")
stan_butte <- readRDS("data-raw/juvenile_abundance/butte_results_STAN_nov_2024.rds")

stan_battle_results_clean <- lapply(stan_battle,
                                    extract_stan_results) |>
  bind_rows()

stan_butte_results_clean <- lapply(stan_butte,
                                   extract_stan_results) |>
  bind_rows()

bugs_battle_results_clean <- lapply(bugs_battle,
                                    extract_bugs_results) |>
  bind_rows()
bugs_butte_results_clean <- lapply(bugs_butte,
                                   extract_bugs_results) |>
  bind_rows()


# convergence -------------------------------------------------------------

# bugs battle
lapply(bugs_battle, is.list) |> unlist() |> sum() # 44/69
# bugs butte
lapply(bugs_butte, is.list) |> unlist() |> sum() # 21/38
# stan battle
lapply(stan_battle, is.list) |> unlist() |> sum() # 37/69
# stan butte
lapply(stan_butte, is.list) |> unlist() |> sum() # 22/38

# which one didn't run for battle/stan
fit <- bugs_battle_results_clean |>
  filter(is.na(error)) |>
  distinct(run_year, life_stage, site) |>
  mutate(id = paste(run_year, life_stage, site, sep = "-"))

all <- trials_to_fit_battle |>
  mutate(id = paste(run_year, life_stage, site, sep = "-"))

# battle sites that didn't fit in the stan model:
all$id[!all$id %in% fit$id]

test <- run_single_bt_spas_x_stan(SRJPEmodel::bt_spas_x_bayes_params,
                                weekly_catch,
                                weekly_efficiency,
                                site = "ubc",
                                run_year = 2022,
                                lifestage = "fry",
                                effort_adjust = F)

# streams in both ---------------------------------------------------------

# battle
bugs_fit <- bugs_battle_results_clean |>
  filter(is.na(error)) |>
  distinct(run_year, life_stage, site) |>
  mutate(id = paste(run_year, life_stage, site, sep = "-"))

stan_fit <- stan_battle_results_clean |>
  filter(is.na(error)) |>
  distinct(run_year, life_stage, site) |>
  mutate(id = paste(run_year, life_stage, site, sep = "-"))

# fit in bugs but not stan
bugs_fit$id[!bugs_fit$id %in% stan_fit$id]
# fit in stan but not bugs
stan_fit$id[!stan_fit$id %in% bugs_fit$id]
# fit in both
both_fit_battle <- intersect(bugs_fit$id, stan_fit$id)

# butte
bugs_fit <- bugs_butte_results_clean |>
  filter(is.na(error)) |>
  distinct(run_year, life_stage, site) |>
  mutate(id = paste(run_year, life_stage, site, sep = "-"))

stan_fit <- stan_butte_results_clean |>
  filter(is.na(error)) |>
  distinct(run_year, life_stage, site) |>
  mutate(id = paste(run_year, life_stage, site, sep = "-"))

# fit in bugs but not stan
bugs_fit$id[!bugs_fit$id %in% stan_fit$id]
# fit in stan but not bugs
stan_fit$id[!stan_fit$id %in% bugs_fit$id]
# fit in both
both_fit_butte <- intersect(bugs_fit$id, stan_fit$id)

# plot battle --------------------------------------------------------------------

# weekly abundance
bugs_battle_results_clean |>
  filter(is.na(error),
         converged) |>
  mutate(parameter = gsub("[0-9]+|\\[|\\]", "", parameter),
         model = "BUGS") |>
  bind_rows(stan_battle_results_clean |>
              filter(is.na(error),
                     converged) |>
            mutate(model = "STAN",
                   statistic = str_remove_all(statistic, "\\%"))) |>
  mutate(id = paste(run_year, life_stage, site, sep = "-")) |>
  filter(parameter == "N",
         statistic == "mean",
         id %in% both_fit_battle) |>
  ggplot(aes(x = week_fit, y = value, color = model)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~id, scales = "free_y")

# Ntot
bugs_battle_results_clean |>
  filter(is.na(error)) |>
  mutate(parameter = gsub("[0-9]+|\\[|\\]", "", parameter)) |>
  filter(parameter == "N",
         statistic == "50") |>
  group_by(site, run_year, life_stage) |>
  summarise(Ntot = sum(value)) |>
  ungroup() |>
  mutate(model = "BUGS") |>
  bind_rows(stan_battle_results_clean |>
              filter(is.na(error),
                     parameter == "Ntot",
                     statistic == "50%") |>
              rename(Ntot = value) |>
              mutate(model = "STAN")) |>
  select(site, run_year, life_stage, Ntot, model) |>
  mutate(id = paste(run_year, life_stage, site, sep = "-")) |>
  filter(id %in% both_fit_battle) |>
  ggplot(aes(x = run_year, y = Ntot, color = model)) +
  geom_point(aes(shape = model), alpha = 0.9) +
  facet_wrap(~site, scales = "free_y") +
  theme_bw() +
  labs(x = "Run Year",
       y = "Annual abundance") +
  scale_y_continuous(labels = label_comma())

# efficiency
bugs_battle_results_clean |>
  filter(is.na(error)) |>
  mutate(parameter = gsub("[0-9]+|\\[|\\]", "", parameter),
         model = "BUGS") |>
  bind_rows(stan_battle_results_clean |>
              filter(is.na(error)) |>
              mutate(model = "STAN",
                     statistic = str_remove_all(statistic, "\\%"))) |>
  mutate(id = paste(run_year, life_stage, site, sep = "-")) |>
  filter(parameter == "pCap_U",
         statistic == "mean",
         id %in% both_fit_battle) |>
  ggplot(aes(x = week_fit, y = value, color = model)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~id, scales = "free_y")

# TODO
# plot efficiency parameters - they should be the same
data <- trials_to_fit_battle |>
  mutate(id = paste(run_year, life_stage, site, sep = "-"))
which(data$id == "2006-fry-ubc")

full_model_object_2006_fry_ubc <- bugs_battle[[30]]
plot_juvenile_abundance("ubc", 2006, "fry", full_model_object_2006_fry_ubc$final_results)
full_model_object_2006_fry_ubc$final_results


# then - is the abundance part of the model different
# then take robust dataset (few missing weeks) and run it

# plot butte --------------------------------------------------------------

# weekly abundance
bugs_butte_results_clean |>
  filter(is.na(error),
         converged) |>
  mutate(parameter = gsub("[0-9]+|\\[|\\]", "", parameter),
         model = "BUGS") |>
  bind_rows(stan_butte_results_clean |>
              filter(is.na(error),
                     converged) |>
              mutate(model = "STAN",
                     statistic = str_remove_all(statistic, "\\%"))) |>
  mutate(id = paste(run_year, life_stage, site, sep = "-")) |>
  filter(parameter == "N",
         statistic == "mean",
         id %in% both_fit_butte) |>
  ggplot(aes(x = week_fit, y = value, color = model)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~id, scales = "free_y")

# Ntot
bugs_butte_results_clean |>
  filter(is.na(error)) |>
  mutate(parameter = gsub("[0-9]+|\\[|\\]", "", parameter)) |>
  filter(parameter == "N",
         statistic == "50") |>
  group_by(site, run_year, life_stage) |>
  summarise(Ntot = sum(value)) |>
  ungroup() |>
  mutate(model = "BUGS") |>
  bind_rows(stan_butte_results_clean |>
              filter(is.na(error),
                     parameter == "Ntot",
                     statistic == "50%") |>
              rename(Ntot = value) |>
              mutate(model = "STAN")) |>
  select(site, run_year, life_stage, Ntot, model) |>
  mutate(id = paste(run_year, life_stage, site, sep = "-")) |>
  filter(id %in% both_fit_butte) |>
  ggplot(aes(x = run_year, y = Ntot, color = model)) +
  geom_point(aes(shape = model), alpha = 0.9) +
  facet_wrap(~site, scales = "free_y") +
  theme_bw() +
  labs(x = "Run Year",
       y = "Annual abundance") +
  scale_y_continuous(labels = label_comma())

# efficiency
bugs_butte_results_clean |>
  filter(is.na(error)) |>
  mutate(parameter = gsub("[0-9]+|\\[|\\]", "", parameter),
         model = "BUGS") |>
  bind_rows(stan_butte_results_clean |>
              filter(is.na(error)) |>
              mutate(model = "STAN",
                     statistic = str_remove_all(statistic, "\\%"))) |>
  mutate(id = paste(run_year, life_stage, site, sep = "-")) |>
  filter(parameter == "pCap_U",
         statistic == "mean",
         id %in% both_fit_butte) |>
  ggplot(aes(x = week_fit, y = value, color = model)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~id, scales = "free_y")

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


# 11-14-2024 --------------------------------------------------------------

# run 10 sites on clear, all lifestages
# do not run for yearling
weekly_catch <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
  filter(life_stage %in% c("fry", "smolt"))

SRJPEdata::weekly_juvenile_abundance_efficiency_data |>
  filter(number_released == 0)

weekly_efficiency <- SRJPEdata::weekly_juvenile_abundance_efficiency_data

# only fit for site/run year combos with data
trials_to_fit <- weekly_catch |>
  mutate(filter_out = ifelse(is.na(life_stage) & count > 0, TRUE, FALSE)) |> # we do not want to keep NA lifestage associated with counts > 0
  filter(!filter_out,
         site == "lcc",
         week %in% c(seq(45, 53), seq(1, 22))) |>
  mutate(count = round(count, 0),
         catch_standardized_by_hours_fished = round(catch_standardized_by_hours_fished, 0)) |>
  group_by(site, run_year) |>
  tally() |>
  arrange(desc(n)) |>
  head(10)

# bugs
# butte
bugs_results <- purrr::pmap(list(trials_to_fit$site,
                                 trials_to_fit$run_year,
                                 NA),
                                 run_multiple_bugs,
                                 .progress = TRUE)
saveRDS(bugs_results, here::here("data-raw", "juvenile_abundance",
                                 "clear_results_BUGS_nov_2024.rds"))

# stan
options(mc.cores=parallel::detectCores())
stan_results <- purrr::pmap(list(trials_to_fit$site,
                                        trials_to_fit$run_year,
                                        NA),
                                   run_multiple_stan,
                                   .progress = TRUE)
saveRDS(stan_results, here::here("data-raw", "juvenile_abundance",
                                 "clear_results_STAN_nov_2024.rds"))

# plot diagnostics
plot_all <- function(element, model_language) {
  if(!is.list(element)) {
    print("model did not fit")
  } else {
    generate_diagnostic_plot(element$final_results, model_language)
  }
}

purrr::pmap(list(element = bugs_results,
                model_language = "BUGS"),
           plot_all)

generate_diagnostic_plot(bugs_results[[3]]$final_results,
                         "BUGS")
