library(tidyverse)
library(SRJPEdata)

# bt-spas-x ---------------------------------------------------------------

bt_spas_x_results <- run_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                   bt_spas_x_input_data = SRJPEdata::weekly_juvenile_abundance_model_data,
                                   site = "ubc",
                                   run_year = 2009,
                                   effort_adjust = T,
                                   multi_run_mode = F,
                                   mainstem_version = F,
                                   bugs_directory = "not working")

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
