library(tidyverse)
library(SRJPEdata)
library(rstan)
library(R2WinBUGS)
library(scales)

# bt-spas-x ---------------------------------------------------------------

# single run
bt_spas_x_results <- run_single_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                          bt_spas_x_input_data = SRJPEdata::weekly_juvenile_abundance_model_data,
                                          site = "ubc",
                                          run_year = 2004,
                                          lifestage = "fry",
                                          effort_adjust = F,
                                          # when running on remote computer, we have the WinBUGS14 in this folder
                                          # however, messaging should be included here to say that the
                                          # WinBUGS folder needs to be in a folder where you have write permissions
                                          # because this errors out frequently
                                          bugs_directory = here::here("data-raw", "WinBUGS14"),
                                          debug_mode = FALSE,
                                          # TODO remove once we make a decision
                                          no_cut = T)
results <- readRDS("data-raw/juvenile_abundance/ubc_2004_2024-07-23.rds")

# plot single run results
plot_juvenile_abundance(results$final_results)
plot_weekly_capture_probability(results$final_results)

# test with smaller df
multi_run_results <- run_multiple_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                            SRJPEdata::weekly_juvenile_abundance_model_data,
                                            effort_adjust = T,
                                            bugs_directory = here::here("data-raw", "WinBUGS14"),
                                            debug_mode = F)
readr::write_rds(multi_run_results, paste0("data-raw/juvenile_abundance/multi_run_results_", Sys.Date(), ".rds"))

# explore multi run results
multi_run_results <- read_rds("data-raw/juvenile_abundance/multi_run_results.rds")
non_errored_results <- unlist(lapply(multi_run_results, function(x) is_tibble(x)))
multi_run_df <- multi_run_results[non_errored_results] |>
  bind_rows()

# this plot code is stolen from juvenile_abundance_plots.R, will be functionalized
julian_week_to_date_lookup <- read.table(file = "data-raw/juvenile_abundance/btspas_model_code/Jwk_Dates.txt", header = F) |>
  tibble() |>
  filter(V1 != "Jwk") |>
  mutate(V1 = as.numeric(V1)) |>
  select(Jwk = V1, date = V2)

# plot
multi_run_df |>
  filter(site == "eye riffle",
         str_detect(parameter, "N\\["),
         life_stage == "fry") |>
  left_join(julian_week_to_date_lookup, by = c("week_fit" = "Jwk")) |>
  pivot_wider(id_cols = c("date", "site", "run_year", "life_stage", "week_fit"),
              names_from = "statistic",
              values_from = "value") |>
  mutate(fake_date = ifelse(week_fit > 35, paste0(1969, "-", date),
                            paste0(1970, "-", date)),
           # fake_date = ifelse(week_fit > 35, paste0(run_year - 1, "-", date),
           #                  paste0(run_year, "-", date)),
         fake_date = as.Date(fake_date, format = "%Y-%b-%d"),
         week = factor(week_fit,
                       levels = c(35:53, 1:34))) |>
  ggplot(aes(x = fake_date, y = `50`)) +
  geom_bar(stat = "identity", fill = "grey") +
  # geom_errorbar(aes(x = fake_date, ymin = `2.5`, ymax = `97.5`), width = 0.2) +
  facet_wrap(~run_year, scales = "free_y", nrow = 4) +
  theme_minimal() +
  labs(x = "Date", y = "Abundance",
       title = "Upper Battle Creek") +
  scale_x_date(breaks = "1 week", date_labels = "%b-%d") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# passage to spawner ------------------------------------------------------

P2S_results <- run_passage_to_spawner_model(SRJPEdata::observed_adult_input,
                                            SRJPEdata::adult_model_covariates_standard,
                                            "battle creek", "wy_type", FALSE)
P2S_spawners <- get_predicted_spawners_from_P2S(P2S_results$formatted_pars)

P2S_spawners |>
  ggplot(aes(x = year, y = median_predicted_spawners)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2) +
  geom_line() +
  theme_minimal()


# survival model ----------------------------------------------------------

# explore results for survival model
survival_results <- run_survival_model(SRJPEdata::survival_model_inputs,
                                       number_detection_locations = 5,
                                       number_reaches = 4)

