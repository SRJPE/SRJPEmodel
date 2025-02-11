library(tidyverse)
library(SRJPEdata)
library(rstan)
library(R2WinBUGS)
library(scales)

# bt-spas-x ---------------------------------------------------------------
# example for upper battle creek (ubc) 2018 (trib)

# run pCap model
pCap_inputs <- prepare_pCap_inputs(mainstem = FALSE)
pCap <- fit_pCap_model(pCap_inputs$inputs)

# run abundance model
abundance_inputs <- prepare_abundance_inputs("ubc", 2018, effort_adjust = T)
lt_pCap_Us <- generate_lt_pCap_Us(abundance_inputs, pCap)

abundance <- fit_abundance_model_BUGS(abundance_inputs, lt_pCap_Us,
                                      # point towards where you store the .bug model
                                      "C:/Users/Liz/Documents/SRJPEmodel/model_files/abundance_model.bug",
                                      # point to where you have WinBUGS
                                      "C:/Users/Liz/Documents/SRJPEmodel/data-raw/WinBUGS14")

abundance_table <- extract_abundance_estimates("ubc", 2018, abundance_inputs, abundance)

generate_diagnostic_plot("ubc", 2018, abundance_table)

# example for a mainstem site (tisdale)
tis_inputs <- prepare_pCap_inputs(mainstem = TRUE,
                                  "tisdale")
tis_fit <- fit_pCap_model(tis_inputs$inputs)
tis_abundance_inputs <- prepare_abundance_inputs("tisdale", 2017, effort_adjust = T)
tis_lt_pCap_Us <- generate_lt_pCap_Us(tis_abundance_inputs, tis_fit)

tis_abundance <- fit_abundance_model_BUGS(tis_abundance_inputs, tis_lt_pCap_Us,
                                          # point towards where you store the .bug model
                                          "C:/Users/Liz/Documents/SRJPEmodel/model_files/abundance_model.bug",
                                          # point to where you have WinBUGS
                                          "C:/Users/Liz/Documents/SRJPEmodel/data-raw/WinBUGS14")

tis_abundance_table <- extract_abundance_estimates("tisdale", 2017,
                                                   tis_abundance_inputs, tis_abundance)



# passage to spawner ------------------------------------------------------

P2S_inputs <- prepare_P2S_inputs("battle creek", "wy_type")
P2S_results <- run_passage_to_spawner_model(P2S_inputs)
P2S_spawners <- extract_P2S_estimates(P2S_results)

P2S_spawners |>
  filter(parameter == "predicted_spawners") |>
  ggplot(aes(x = year, y = `50`)) +
  geom_ribbon(aes(ymin = `2.5`, ymax = `97.5`), alpha = 0.2) +
  geom_line() +
  theme_minimal()


# survival model ----------------------------------------------------------

# explore results for survival model
survival_results <- run_survival_model(SRJPEdata::survival_model_inputs,
                                       number_detection_locations = 5,
                                       number_reaches = 4)

