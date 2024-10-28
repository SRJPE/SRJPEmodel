# script for Josh
# October 28, 2024

# install SRJPE packages - modeling and data
install.packages("remotes") # this allows you to install packages from github
library(remotes)

remotes::install_github("SRJPE/SRJPEdata") # note redd data is in progress
remotes::install_github("SRJPE/SRJPEmodel")

library(SRJPEdata)
library(SRJPEmodel)
library(tidyverse)

# see all data
data(package = "SRJPEdata")

# RST data for BT_SPAS-X
# catch data
SRJPEdata::weekly_juvenile_abundance_catch_data |>
  glimpse()
# efficiency data
SRJPEdata::weekly_juvenile_abundance_efficiency_data |>
  glimpse()

# documentation
?SRJPEdata::weekly_juvenile_abundance_catch_data

# run bt-spas-x in WinBUGS
bt_spas_x_results <- run_single_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                          SRJPEdata::weekly_juvenile_abundance_catch_data,
                                          SRJPEdata::weekly_juvenile_abundance_efficiency_data,
                                          site = "ubc",
                                          run_year = 2009,
                                          lifestage = "fry",
                                          effort_adjust = T,
                                          # put your bugs_directory in here as a filepath
                                          bugs_directory = here::here("data-raw", "WinBUGS14"),
                                          debug_mode = FALSE,
                                          # whether you want to run the function with the cut() call
                                          no_cut = FALSE)

# run bt-spas-x in stan
bt_spas_x_results_stan <- run_single_bt_spas_x_stan(SRJPEmodel::bt_spas_x_bayes_params,
                                                    SRJPEdata::weekly_juvenile_abundance_catch_data,
                                                    SRJPEdata::weekly_juvenile_abundance_efficiency_data,
                                                    site = "ubc",
                                                    run_year = 2009,
                                                    lifestage = "fry",
                                                    effort_adjust = T)

# run P2S
P2S_results <- run_passage_to_spawner_model(SRJPEdata::observed_adult_input,
                                            SRJPEdata::p2s_model_covariates_standard,
                                            stream_name = "battle creek",
                                            selected_covariate = "wy_type",
                                            extract_predicted_spawners = FALSE)
