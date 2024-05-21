# run deer and mill 2023 through bt spas

library(tidyverse)
library(SRJPEdata)

# prep data ---------------------------------------------------------------
# full

deer_2023 <- SRJPEdata::weekly_juvenile_abundance_model_data_mill_deer_2022_2023 |>
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

deer_mill_2023_input <- SRJPEdata::weekly_juvenile_abundance_model_data |>
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
  ungroup() |>
  bind_rows(deer_2023)

# split by lifestage
yearling <- SRJPEdata::weekly_juvenile_abundance_model_data |>
  bind_rows(SRJPEdata::weekly_juvenile_abundance_model_data_mill_deer_2022_2023) |>
  filter(life_stage == "yearling")

YOY <- SRJPEdata::weekly_juvenile_abundance_model_data |>
  bind_rows(SRJPEdata::weekly_juvenile_abundance_model_data_mill_deer_2022_2023) |>
  filter(life_stage %in% c("fry", "smolt", NA)) |>
  # group by and sum across lifestages so that there are no repeat weeks
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
  ungroup() |>
  glimpse()

# deer --------------------------------------------------------------------

# full



# deer_2023_btspas <- run_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
#                                    bt_spas_x_input_data = deer_mill_2023_input,
#                                    site = "deer creek",
#                                    run_year = 2023,
#                                    effort_adjust = T,
#                                    multi_run_mode = F,
#                                    mainstem_version = F,
#                                    bugs_directory = here::here("data-raw", "WinBUGS14"),
#                                    debug_mode = TRUE)

deer_2023 <- readRDS(here::here("data-raw", "juvenile_abundance",
                                "deer_2023_model_fits.rds"))

test <- extract_bt_spas_x_results(deer_2023$results$model_results)
get_total_juvenile_abundance(test$summary_output)
# how to get the weeks modeled?
get_weekly_juvenile_abundance(test$summary_output) |>
  ggplot(aes(x = week_index, y = median_abundance)) +
  geom_line() +
  geom_ribbon(aes(ymin = lcl_97_5, ymax = ucl_97_5),
              alpha = 0.2) +
  theme_minimal()

# YOY vs Yearling
deer_yearling <- run_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                               bt_spas_x_input_data = yearling,
                               site = "deer creek",
                               run_year = 2023,
                               effort_adjust = F,
                               multi_run_mode = F,
                               mainstem_version = F,
                               bugs_directory = here::here("data-raw", "WinBUGS14"),
                               debug_mode = FALSE)

readr::write_rds(deer_yearling, "data-raw/juvenile_abundance/deer_2023_yearling_model_fits.rds")

deer_YOY <- run_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                          bt_spas_x_input_data = YOY,
                          site = "deer creek",
                          run_year = 2023,
                          effort_adjust = F,
                          multi_run_mode = F,
                          mainstem_version = F,
                          bugs_directory = here::here("data-raw", "WinBUGS14"),
                          debug_mode = FALSE)

readr::write_rds(deer_YOY, "data-raw/juvenile_abundance/deer_2023_YOY_model_fits.rds")




# mill --------------------------------------------------------------------
# full
# mill_2023_btspas <- run_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
#                                   bt_spas_x_input_data = deer_mill_2023_input,
#                                   site = "mill creek",
#                                   run_year = 2023,
#                                   effort_adjust = F,
#                                   multi_run_mode = F,
#                                   mainstem_version = F,
#                                   bugs_directory = here::here("data-raw", "WinBUGS14"),
#                                   debug_mode = TRUE)

mill_2023_btspas <- readRDS(here::here("data-raw", "juvenile_abundance",
                                       "mill_2023_model_fits.rds"))

mill_2023 <- extract_bt_spas_x_results(mill_2023_btspas$results$model_results)
get_total_juvenile_abundance(mill_2023$summary_output)
get_weekly_juvenile_abundance(mill_2023$summary_output) |>
  ggplot(aes(x = week_index, y = median_abundance)) +
  geom_line() +
  geom_ribbon(aes(ymin = lcl_97_5, ymax = ucl_97_5),
              alpha = 0.2) +
  theme_minimal()

# YOY vs yearling

mill_yearling <- run_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                               bt_spas_x_input_data = yearling,
                               site = "mill creek",
                               run_year = 2023,
                               effort_adjust = F,
                               multi_run_mode = F,
                               mainstem_version = F,
                               bugs_directory = here::here("data-raw", "WinBUGS14"),
                               debug_mode = FALSE)

readr::write_rds(mill_yearling, "data-raw/juvenile_abundance/mill_2023_yearling_model_fits.rds")

mill_YOY <- run_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                          bt_spas_x_input_data = YOY,
                          site = "mill creek",
                          run_year = 2023,
                          effort_adjust = F,
                          multi_run_mode = F,
                          mainstem_version = F,
                          bugs_directory = here::here("data-raw", "WinBUGS14"),
                          debug_mode = FALSE)

readr::write_rds(mill_YOY, "data-raw/juvenile_abundance/mill_2023_YOY_model_fits.rds")


# analysis and plotting ---------------------------------------------------
palette <- c("#453F3C", "#C2847A", "#63A088", "#255F85", "#EEE0CB")
big_palette <- colorRampPalette(palette)(20)

# deer
deer_yearling <- readRDS("data-raw/juvenile_abundance/deer_2023_yearling_model_fits.rds")
deer_YOY <- readRDS("data-raw/juvenile_abundance/deer_2023_YOY_model_fits.rds")

get_weekly_juvenile_abundance(deer_yearling$results$model_results$summary) |>
  mutate(life_stage = "yearling") |>
  bind_rows(get_weekly_juvenile_abundance(deer_YOY$results$model_results$summary) |>
              mutate(life_stage = "YOY")) |>
  ggplot(aes(x = week_index, y = median_abundance)) +
  geom_line(aes(color = life_stage), linewidth = 1) +
  geom_ribbon(aes(ymin = lcl_97_5, ymax = ucl_97_5, group = life_stage),
              alpha = 0.2) +
  theme_minimal() +
  facet_wrap(~life_stage, nrow = 2) +
  scale_color_manual(values = palette) +
  theme(legend.position = "bottom") +
  labs(x = "Week",
       y = "Abundance")


# mill
mill_yearling <- readRDS("data-raw/juvenile_abundance/mill_2023_yearling_model_fits.rds")
mill_YOY <- readRDS("data-raw/juvenile_abundance/mill_2023_YOY_model_fits.rds")

get_weekly_juvenile_abundance(mill_yearling$results$model_results$summary) |>
  mutate(life_stage = "yearling") |>
  bind_rows(get_weekly_juvenile_abundance(mill_YOY$results$model_results$summary) |>
              mutate(life_stage = "YOY")) |>
  ggplot(aes(x = week_index, y = median_abundance)) +
  geom_line(aes(color = life_stage), linewidth = 1) +
  geom_ribbon(aes(ymin = lcl_97_5, ymax = ucl_97_5, group = life_stage),
              alpha = 0.2) +
  theme_minimal() +
  facet_wrap(~life_stage, scales = "free_y", nrow = 2) +
  scale_color_manual(values = palette) +
  theme(legend.position = "bottom") +
  labs(x = "Week",
       y = "Abundance")
