# script for Josh

remotes::install_github("SRJPE/SRJPEdata") # note redd data is in progress
remotes::install_github("SRJPE/SRJPEmodel@wip")

library(SRJPEdata)
library(SRJPEmodel)
library(tidyverse)

# see all data
data(package = "SRJPEdata")

# RST data for BT_SPAS-X
SRJPEdata::weekly_juvenile_abundance_model_data |>
  glimpse()

# documentation
?SRJPEdata::weekly_juvenile_abundance_model_data

# run bt-spas-x
bt_spas_x_results <- run_single_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                          bt_spas_x_input_data = SRJPEdata::weekly_juvenile_abundance_model_data,
                                          site = "eye riffle",
                                          run_year = 1998,
                                          lifestage = "fry",
                                          effort_adjust = T,
                                          mainstem_version = F,
                                          # put your bugs_directory in here as a filepath
                                          bugs_directory = here::here("data-raw", "WinBUGS14"),
                                          debug_mode = FALSE)

# run P2S
P2S_results <- run_passage_to_spawner_model(SRJPEdata::observed_adult_input,
                                            SRJPEdata::adult_model_covariates_standard,
                                            stream_name = "battle creek",
                                            selected_covariate = "wy_type",
                                            extract_predicted_spawners = TRUE)
