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
