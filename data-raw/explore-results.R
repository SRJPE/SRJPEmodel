library(tidyverse)
library(SRJPEdata)
library(rstan)
library(R2WinBUGS)

# bt-spas-x ---------------------------------------------------------------
# TODO test for different sites - put in stop if running on a site with no data in catch

updated_input <- SRJPEdata::weekly_juvenile_abundance_model_data |>
  group_by(year, week, stream, site, run_year) |>
  summarise(count = sum(count, na.rm = T),
            mean_fork_length = mean(mean_fork_length, na.rm = T),
            number_released = sum(number_released),
            number_recaptured = sum(number_recaptured),
            hours_fished = mean(hours_fished, na.rm = T),
            flow_cfs = mean(flow_cfs, na.rm = T),
            standardized_flow = mean(standardized_flow, na.rm = T),
            catch_standardized_by_hours_fished = count * hours_fished,
            standardized_efficiency_flow = mean(standardized_efficiency_flow, na.rm = T),
            lgN_prior = mean(lgN_prior, na.rm = T)) |>
  ungroup()

bt_spas_x_results <- run_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                   bt_spas_x_input_data = updated_input,
                                   site = "ubc",
                                   run_year = 2009,
                                   life_stage = "fry",
                                   effort_adjust = T,
                                   multi_run_mode = F, # T
                                   mainstem_version = F,
                                   # when running on remote computer, we have the WinBUGS14 in this folder
                                   # however, messaging should be included here to say that the
                                   # WinBUGS folder needs to be in a folder where you have write permissions
                                   # because this errors out frequently
                                   bugs_directory = here::here("data-raw", "WinBUGS14"))


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

