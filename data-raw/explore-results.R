library(tidyverse)
library(SRJPEdata)
library(rstan)
library(R2WinBUGS)

# bt-spas-x ---------------------------------------------------------------
# TODO test for different sites - put in stop if running on a site with no data in catch
bt_spas_x_results <- run_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                   bt_spas_x_input_data = SRJPEdata::weekly_juvenile_abundance_model_data,
                                   site = "ubc",
                                   run_year = 2009,
                                   effort_adjust = T,
                                   multi_run_mode = F, # T
                                   mainstem_version = F,
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
# TODO model runs with 5 detection locations but not 4
survival_results <- run_survival_model(SRJPEdata::survival_model_inputs,
                                       number_detection_locations = 4,
                                       number_reaches = 4)

