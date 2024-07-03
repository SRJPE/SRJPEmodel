# run STAN model
library(tidyverse)
library(SRJPEdata)
library(rstan)
library(R2WinBUGS)

# get all the inputs by running on a mac

updated_input <- SRJPEdata::weekly_juvenile_abundance_model_data |>
  group_by(year, week, stream, site, run_year) |>
  summarise(count = sum(count, na.rm = T),
            mean_fork_length = mean(mean_fork_length, na.rm = T),
            number_released = sum(number_released),
            number_recaptured = sum(number_recaptured),
            hours_fished = mean(hours_fished, na.rm = T),
            flow_cfs = mean(flow_cfs, na.rm = T),
            standardized_flow = mean(standardized_flow, na.rm = T),
            catch_standardized_by_hours_fished = sum(catch_standardized_by_hours_fished),
            standardized_efficiency_flow = mean(standardized_efficiency_flow, na.rm = T),
            lgN_prior = mean(lgN_prior, na.rm = T)) |>
  ungroup()

bt_spas_x_results <- run_single_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                           bt_spas_x_input_data = updated_input,
                                           site = "eye riffle",
                                           run_year = 2005,
                                           effort_adjust = T,
                                           mainstem_version = F,
                                           # when running on remote computer, we have the WinBUGS14 in this folder
                                           # however, messaging should be included here to say that the
                                           # WinBUGS folder needs to be in a folder where you have write permissions
                                           # because this errors out frequently
                                           bugs_directory = here::here("data-raw", "WinBUGS14"),
                                           debug_mode = TRUE)

# now run in STAN
all_mark_recap_STAN_code <- read_file("model_files/all_mark_recap_STAN_chatgtp.txt")

STAN_results <- rstan::stan(model_code = all_mark_recap_STAN_code,
                            data = bt_spas_x_results$results$data,
                            init = bt_spas_x_results$results$inits,
                            thin = SRJPEmodel::bt_spas_x_bayes_params$number_thin,
                            warmup = SRJPEmodel::bt_spas_x_bayes_params$number_burnin,
                            chains = SRJPEmodel::bt_spas_x_bayes_params$number_chains,
                            iter = SRJPEmodel::bt_spas_x_bayes_params$number_mcmc)


# test out no cap ---------------------------------------------------------

no_cap_battle_2004 <-


