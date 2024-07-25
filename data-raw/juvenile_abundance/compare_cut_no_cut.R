library(tidyverse)
library(SRJPEdata)
library(rstan)
library(R2WinBUGS)

# compare cut and no cut

cut_results_ubc_2004 <- run_single_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                             bt_spas_x_input_data = SRJPEdata::weekly_juvenile_abundance_model_data,
                                             site = "ubc",
                                             run_year = 2004,
                                             lifestage = "fry",
                                             effort_adjust = F,
                                             mainstem_version = F,
                                             bugs_directory = here::here("data-raw", "WinBUGS14"),
                                             debug_mode = FALSE,
                                             no_cut = FALSE)

no_cut_results_ubc_2004 <- run_single_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                             bt_spas_x_input_data = SRJPEdata::weekly_juvenile_abundance_model_data,
                                             site = "ubc",
                                             run_year = 2004,
                                             lifestage = "fry",
                                             effort_adjust = F,
                                             mainstem_version = F,
                                             bugs_directory = here::here("data-raw", "WinBUGS14"),
                                             debug_mode = FALSE,
                                             no_cut = TRUE)

cut_results_kdl_2020 <- run_single_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                             bt_spas_x_input_data = SRJPEdata::weekly_juvenile_abundance_model_data,
                                             site = "knights landing",
                                             run_year = 2020,
                                             lifestage = "fry",
                                             effort_adjust = F,
                                             mainstem_version = F,
                                             bugs_directory = here::here("data-raw", "WinBUGS14"),
                                             debug_mode = FALSE,
                                             no_cut = FALSE)

no_cut_results_kdl_2020 <- run_single_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                                bt_spas_x_input_data = SRJPEdata::weekly_juvenile_abundance_model_data,
                                                site = "knights landing",
                                                run_year = 2020,
                                                lifestage = "fry",
                                                effort_adjust = F,
                                                mainstem_version = F,
                                                bugs_directory = here::here("data-raw", "WinBUGS14"),
                                                debug_mode = FALSE,
                                                no_cut = TRUE)

cut_results_lf_2022 <- run_single_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                             bt_spas_x_input_data = SRJPEdata::weekly_juvenile_abundance_model_data,
                                             site = "lower feather river",
                                             run_year = 2022,
                                             lifestage = "fry",
                                             effort_adjust = F,
                                             mainstem_version = F,
                                             bugs_directory = here::here("data-raw", "WinBUGS14"),
                                             debug_mode = FALSE,
                                             no_cut = FALSE)

no_cut_results_lf_2022 <- run_single_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                                bt_spas_x_input_data = SRJPEdata::weekly_juvenile_abundance_model_data,
                                                site = "lower feather river",
                                                run_year = 2022,
                                                lifestage = "fry",
                                                effort_adjust = F,
                                                mainstem_version = F,
                                                bugs_directory = here::here("data-raw", "WinBUGS14"),
                                                debug_mode = FALSE,
                                                no_cut = TRUE)

cut_results_deer_1996 <- run_single_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                            bt_spas_x_input_data = SRJPEdata::weekly_juvenile_abundance_model_data,
                                            site = "deer creek",
                                            run_year = 1996,
                                            lifestage = "fry",
                                            effort_adjust = F,
                                            mainstem_version = F,
                                            bugs_directory = here::here("data-raw", "WinBUGS14"),
                                            debug_mode = FALSE,
                                            no_cut = FALSE)

no_cut_results_deer_1996 <- run_single_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                               bt_spas_x_input_data = SRJPEdata::weekly_juvenile_abundance_model_data,
                                               site = "deer creek",
                                               run_year = 1996,
                                               lifestage = "fry",
                                               effort_adjust = F,
                                               mainstem_version = F,
                                               bugs_directory = here::here("data-raw", "WinBUGS14"),
                                               debug_mode = FALSE,
                                               no_cut = TRUE)

cut_comparison_all <- bind_rows(cut_results_ubc_2004$final_results %>%
                                  mutate(cut = T),
                                no_cut_results_ubc_2004$final_results %>%
                                  mutate(cut = F),
                                cut_results_kdl_2020$final_results %>%
                                  mutate(cut = T),
                                no_cut_results_kdl_2020$final_results %>%
                                  mutate(cut = F),
                                cut_results_lf_2022$final_results %>%
                                  mutate(cut = T),
                                no_cut_results_lf_2022$final_results %>%
                                  mutate(cut = F),
                                cut_results_deer_1996$final_results %>%
                                  mutate(cut = T),
                                no_cut_results_deer_1996$final_results %>%
                                  mutate(cut = F))

saveRDS(cut_comparison_all, file = "data-raw/juvenile_abundance/cut_comparison_results_2024-07-25.rds")



compare_ubc_2004 <- cut_results_ubc_2004$final_results %>%
  left_join(no_cut_results_ubc_2004$final_results %>%
              rename(no_cut_value = value)) %>%
  filter(statistic == "mean") %>%
  mutate(cut_diff = value - no_cut_value) %>%
  glimpse()

compare_ubc_2004 %>%
  ggplot(aes(x = cut_diff, fill = parameter)) +
  geom_histogram() +
  facet_wrap(~parameter, scales = "free") +
  theme_minimal() +
  theme(legend.position = "none")
