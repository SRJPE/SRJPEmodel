# run passage to spawner
library(SRJPEdata)
library(SRJPEmodel)
library(tidyverse)
library(rstan)

# TODO move this to vignette
# pull in data from SRJPEdata ---------------------------------------------
observed_adult_input <- SRJPEdata::observed_adult_input
adult_model_covariates <- SRJPEdata::adult_model_covariates_standard


battle_results <- run_passage_to_spawner_model(observed_adult_input,
                                               adult_model_covariates,
                                               "battle creek",
                                               "wy_type")

clear_results <- run_passage_to_spawner_model(observed_adult_input,
                                              adult_model_covariates,
                                              "clear creek",
                                              "wy_type")

deer_results <- run_passage_to_spawner_model(observed_adult_input,
                                             adult_model_covariates,
                                             "deer creek",
                                             "wy_type")

mill_results <- run_passage_to_spawner_model(observed_adult_input,
                                             adult_model_covariates,
                                             "mill creek",
                                             "wy_type")

# write model summaries ---------------------------------------------------
P2S_model_fits <- bind_rows(battle_results$formatted_pars,
                            clear_results$formatted_pars,
                            mill_results$formatted_pars,
                            deer_results$formatted_pars) |>
  glimpse()

# these are preliminary!
butte_results <- run_passage_to_spawner_model(observed_adult_input,
                                              adult_model_covariates,
                                              "butte creek",
                                              "gdd_std")


# try in comparison mode --------------------------------------------------
comparison_results <- compare_P2S_model_covariates(observed_adult_input,
                                                   adult_model_covariates)


# save as data object ----------------------------------------------------
# this used to be saved in the cloud as "jpe-model-data/adult-model/P2S_model_fits.csv"
usethis::use_data(P2S_model_fits, overwrite = TRUE)
