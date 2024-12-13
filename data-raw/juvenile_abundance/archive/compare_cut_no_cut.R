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
                                             bugs_directory = here::here("data-raw", "WinBUGS14"),
                                             debug_mode = FALSE,
                                             no_cut = FALSE)

no_cut_results_ubc_2004 <- run_single_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                             bt_spas_x_input_data = SRJPEdata::weekly_juvenile_abundance_model_data,
                                             site = "ubc",
                                             run_year = 2004,
                                             lifestage = "fry",
                                             effort_adjust = F,
                                             bugs_directory = here::here("data-raw", "WinBUGS14"),
                                             debug_mode = FALSE,
                                             no_cut = TRUE)

cut_results_kdl_2020 <- run_single_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                             bt_spas_x_input_data = SRJPEdata::weekly_juvenile_abundance_model_data,
                                             site = "knights landing",
                                             run_year = 2020,
                                             lifestage = "fry",
                                             effort_adjust = F,
                                             bugs_directory = here::here("data-raw", "WinBUGS14"),
                                             debug_mode = FALSE,
                                             no_cut = FALSE)

no_cut_results_kdl_2020 <- run_single_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                                bt_spas_x_input_data = SRJPEdata::weekly_juvenile_abundance_model_data,
                                                site = "knights landing",
                                                run_year = 2020,
                                                lifestage = "fry",
                                                effort_adjust = F,
                                                bugs_directory = here::here("data-raw", "WinBUGS14"),
                                                debug_mode = FALSE,
                                                no_cut = TRUE)

cut_results_lf_2022 <- run_single_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                             bt_spas_x_input_data = SRJPEdata::weekly_juvenile_abundance_model_data,
                                             site = "lower feather river",
                                             run_year = 2022,
                                             lifestage = "fry",
                                             effort_adjust = F,
                                             bugs_directory = here::here("data-raw", "WinBUGS14"),
                                             debug_mode = FALSE,
                                             no_cut = FALSE)

no_cut_results_lf_2022 <- run_single_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                                bt_spas_x_input_data = SRJPEdata::weekly_juvenile_abundance_model_data,
                                                site = "lower feather river",
                                                run_year = 2022,
                                                lifestage = "fry",
                                                effort_adjust = F,
                                                bugs_directory = here::here("data-raw", "WinBUGS14"),
                                                debug_mode = FALSE,
                                                no_cut = TRUE)

cut_results_deer_1996 <- run_single_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                            bt_spas_x_input_data = SRJPEdata::weekly_juvenile_abundance_model_data,
                                            site = "deer creek",
                                            run_year = 1996,
                                            lifestage = "fry",
                                            effort_adjust = F,
                                            bugs_directory = here::here("data-raw", "WinBUGS14"),
                                            debug_mode = FALSE,
                                            no_cut = FALSE)

no_cut_results_deer_1996 <- run_single_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                               bt_spas_x_input_data = SRJPEdata::weekly_juvenile_abundance_model_data,
                                               site = "deer creek",
                                               run_year = 1996,
                                               lifestage = "fry",
                                               effort_adjust = F,
                                               bugs_directory = here::here("data-raw", "WinBUGS14"),
                                               debug_mode = FALSE,
                                               no_cut = TRUE)

cut_comparison_all <- bind_rows(cut_results_ubc_2004$final_results |>
                                  mutate(cut = T),
                                no_cut_results_ubc_2004$final_results |>
                                  mutate(cut = F),
                                cut_results_kdl_2020$final_results |>
                                  mutate(cut = T),
                                no_cut_results_kdl_2020$final_results |>
                                  mutate(cut = F),
                                cut_results_lf_2022$final_results |>
                                  mutate(cut = T),
                                no_cut_results_lf_2022$final_results |>
                                  mutate(cut = F),
                                cut_results_deer_1996$final_results |>
                                  mutate(cut = T),
                                no_cut_results_deer_1996$final_results |>
                                  mutate(cut = F))

saveRDS(cut_comparison_all, file = "data-raw/juvenile_abundance/cut_comparison_results_2024-07-25.rds")

palette <- c("#798E87", "#C27D38", "#CCC591", "#29211F")

compare <- readRDS("data-raw/juvenile_abundance/cut_comparison_results_2024-07-25.rds")

compare |>
  filter(str_detect(parameter, "N\\[")) |>
  pivot_wider(id_cols = c(site:srjpedata_version, cut),
              values_from = "value",
              names_from = "statistic") |>
  filter(rhat <= 1.05) |>
  ggplot(aes(x = week_fit, y = `50`, color = cut)) +
  geom_point(position = "jitter") +
  geom_errorbar(aes(x = week_fit, ymin = `2.5`, ymax = `97.5`, color = cut),
                position = "jitter", width = 0.2) +
  facet_wrap(~site, scales = "free") +
  scale_color_manual(values = palette) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "Week", y = "Median weekly abundance",
       title = "Cut vs. no cut model runs for weekly abundance N[]")

compare |>
  filter(str_detect(parameter, "b0_pCap\\[")) |>
  mutate(value = exp(value) / (1 + exp(value))) |>
  mutate(index = substr(parameter, 3, length(parameter)),
         index = readr::parse_number(index)) |>
  pivot_wider(id_cols = c(site:srjpedata_version, cut, index),
              values_from = "value",
              names_from = "statistic") |>
  ggplot(aes(x = index, y = `50`, color = cut)) +
  geom_point(position = "jitter") +
  geom_errorbar(aes(x = index, ymin = `2.5`, ymax = `97.5`, color = cut),
                position = "jitter", width = 0.2) +
  facet_wrap(~site, scales = "free") +
  scale_color_manual(values = palette) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "Trib index", y = "b0_pCap for all tribs",
       title = "Cut vs. no cut model runs for b0_pCap")

compare |>
  filter(str_detect(parameter, "trib_mu")) |>
  mutate(value = exp(value) / (1 + exp(value))) |>
  pivot_wider(id_cols = c(site:srjpedata_version, cut),
              values_from = "value",
              names_from = "statistic") |>
  ggplot(aes(x = 1, y = `50`, color = cut)) +
  geom_point() +
  geom_errorbar(aes(x = 1, ymin = `2.5`, ymax = `97.5`, color = cut), width = 0.2) +
  facet_wrap(~site, scales = "free") +
  scale_color_manual(values = palette) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "", y = "trib_mu",
       title = "Cut vs. no cut model runs for trib_mu")

compare |>
  filter(str_detect(parameter, "pCap_U\\[")) |>
  mutate(index = substr(parameter, 3, length(parameter)),
         index = readr::parse_number(index)) |>
  pivot_wider(id_cols = c(site:srjpedata_version, cut, index),
              values_from = "value",
              names_from = "statistic") |>
  filter(rhat <= 1.05) |>
  ggplot(aes(x = index, y = `50`, color = cut)) +
  geom_point(position = "jitter") +
  geom_errorbar(aes(x = index, ymin = `2.5`, ymax = `97.5`, color = cut),
                width = 0.2, position = "jitter") +
  facet_wrap(~site, scales = "free") +
  scale_color_manual(values = palette) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "Week", y = "pCap_U for all tribs",
       title = "Cut vs. no cut model runs for pCap_U")

compare |>
  filter(str_detect(parameter, "Ntot")) |>
  pivot_wider(id_cols = c(site:srjpedata_version, cut),
              values_from = "value",
              names_from = "statistic") |>
  mutate(run_year = as.factor(run_year)) |>
  ggplot(aes(x = run_year, y = `50`, color = cut)) +
  geom_point(position = "jitter") +
  geom_errorbar(aes(x = run_year, ymin = `2.5`, ymax = `97.5`, color = cut),
                position = "jitter", width = 0.2) +
  facet_wrap(~site, scales = "free") +
  scale_color_manual(values = palette) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "Run year", y = "Median total abundance",
       title = "Cut vs. no cut model runs for total abundance N[]")




# total abundance
compare |>
  filter(site == "lower feather river", run_year == 2022, life_stage == "fry",
         cut) |>
  plot_juvenile_abundance()

compare |>
  filter(site == "deer creek", run_year == 1996, life_stage == "fry",
         cut) |>
  plot_juvenile_abundance()

compare |>
  filter(site == "ubc", run_year == 2004, life_stage == "fry",
         cut) |>
  plot_juvenile_abundance()

compare |>
  filter(site == "knights landing", run_year == 2020, life_stage == "fry",
         cut) |>
  plot_juvenile_abundance()

# capture probability
compare |>
  filter(site == "lower feather river", run_year == 2022, life_stage == "fry",
         cut) |>
  plot_weekly_capture_probability()

compare |>
  filter(site == "deer creek", run_year == 1996, life_stage == "fry",
         cut) |>
  plot_weekly_capture_probability()

compare |>
  filter(site == "ubc", run_year == 2004, life_stage == "fry",
         cut) |>
  plot_weekly_capture_probability()

compare |>
  filter(site == "knights landing", run_year == 2020, life_stage == "fry", cut) |>
  plot_weekly_capture_probability()


# all sites ---------------------------------------------------------------

# run
multi_run_results_cut <- run_multiple_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                            SRJPEdata::weekly_juvenile_abundance_model_data,
                                            effort_adjust = T,
                                            bugs_directory = here::here("data-raw", "WinBUGS14"),
                                            debug_mode = F,
                                            no_cut = F)
readr::write_rds(multi_run_results_cut, paste0("data-raw/juvenile_abundance/multi_run_results_", Sys.Date(), "_cut.rds"))

sites_that_ran <- multi_run_results_cut |>
  mutate(compare_sites = paste0(site, "_", run_year, "_", life_stage))

sites_to_compare <- SRJPEdata::weekly_juvenile_abundance_model_data |>
  mutate(compare_sites = paste0(site, "_", run_year, "_", life_stage)) |>
  filter(compare_sites %in% sites_that_ran$compare_sites) |>
  distinct(site, run_year, life_stage)

multi_run_results_no_cut <- run_multiple_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                                   SRJPEdata::weekly_juvenile_abundance_model_data,
                                                sites_to_compare,
                                                effort_adjust = T,
                                                bugs_directory = here::here("data-raw", "WinBUGS14"),
                                                debug_mode = F,
                                                no_cut = T)
readr::write_rds(multi_run_results_no_cut, paste0("data-raw/juvenile_abundance/multi_run_results_", Sys.Date(), "_no_cut.rds"))

# the above code didn't work (bug in writing results)

# run with fewer sites for smaller analysis

test_sites <- sites_to_compare |>
  arrange(desc(run_year)) |>
  head(20)

analysis_cut <- run_multiple_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                               SRJPEdata::weekly_juvenile_abundance_model_data,
                               test_sites,
                               effort_adjust = T,
                               bugs_directory = here::here("data-raw", "WinBUGS14"),
                               debug_mode = F,
                               no_cut = F)

analysis_no_cut <- run_multiple_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                       SRJPEdata::weekly_juvenile_abundance_model_data,
                                       test_sites,
                                       effort_adjust = T,
                                       bugs_directory = here::here("data-raw", "WinBUGS14"),
                                       debug_mode = F,
                                       no_cut = T)

write_rds(analysis_cut, file = "data-raw/juvenile_abundance/analysis_cut_20_sites.rds")
write_rds(analysis_no_cut, file = "data-raw/juvenile_abundance/analysis_no_cut_20_sites.rds")

analysis_cut <- read_rds("data-raw/juvenile_abundance/analysis_cut_20_sites.rds")
analysis_no_cut <- read_rds("data-raw/juvenile_abundance/analysis_no_cut_20_sites.rds")

# compare
analysis_cut |>
  filter(statistic == "error") |>
  distinct(site, run_year, life_stage)

analysis_no_cut |>
  filter(statistic == "error") |>
  distinct(site, run_year, life_stage)

# tisdale ran with no cut, but not with cut function
data_to_plot <- bind_rows(analysis_cut |>
                            mutate(cut = T),
                          analysis_no_cut |>
                            mutate(cut = F)) |>
  filter(statistic != "error")

combos <- data_to_plot |>
  distinct(site, run_year, life_stage)

for(i in 1:nrow(combos)) {
  cut_abundance_plot <- data_to_plot |>
    filter(site == combos$site[i],
           run_year == combos$run_year[i],
           life_stage == combos$life_stage[i],
           cut) |>
    plot_juvenile_abundance()

   no_cut_abundance_plot <- data_to_plot |>
    filter(site == combos$site[i],
           run_year == combos$run_year[i],
           life_stage == combos$life_stage[i],
           !cut) |>
    plot_juvenile_abundance()

  cut_pcap_plot <- data_to_plot |>
    filter(site == combos$site[i],
           run_year == combos$run_year[i],
           life_stage == combos$life_stage[i],
           cut) |>
    plot_weekly_capture_probability()

  no_cut_pcap_plot <- data_to_plot |>
    filter(site == combos$site[i],
           run_year == combos$run_year[i],
           life_stage == combos$life_stage[i],
           !cut) |>
    plot_weekly_capture_probability()

  gridExtra::grid.arrange(cut_abundance_plot, cut_pcap_plot,
               no_cut_abundance_plot, no_cut_pcap_plot,
               ncol = 2)

}
