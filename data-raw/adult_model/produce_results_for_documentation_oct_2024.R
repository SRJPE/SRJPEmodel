# code to run model results for P2S
# 10-16-2024

library(ggplot2)
library(tidyverse)
library(bayesplot)
library(wesanderson)
library(SRJPEdata)
library(SRJPEmodel)
library(here)

# exclude Deer and Mill Creek pending redd data troubleshooting
streams_to_use <- c("battle creek", "clear creek")

# prep data ------------------------------------------------------------
observed_adult_input <- SRJPEdata::observed_adult_input
adult_model_covariates <- SRJPEdata::p2s_model_covariates_standard

observed_adult_input_wide <- observed_adult_input |>
  select(-reach) |> # empty
  group_by(year, stream, data_type) |>
  summarise(count = sum(count, na.rm = T)) |> # count adipose clipped, run together
  ungroup() |>
  pivot_wider(id_cols = c("year", "stream"),
              names_from = data_type, values_from = count) |>
  glimpse()

# get details on sample sizes
observed_adult_input_wide |>
  filter(!stream %in% c("butte creek", "yuba river", "feather river")) |>
  mutate(spawner_count = ifelse(stream == "deer creek",
                                holding_count, redd_count),
         keep = ifelse(is.na(spawner_count) | is.na(upstream_estimate),
                       FALSE, TRUE)) |>
  filter(keep) |>
  group_by(stream) |>
  tally() |>
  clipr::write_clip()

# plot
observed_adult_input_wide |>
  filter(!stream %in% c("butte creek", "yuba river", "feather river")) |>
  mutate(stream = str_to_title(stream),
         spawner_count = ifelse(stream == "Deer Creek",
                                holding_count, redd_count)) |>
  ggplot(aes(x = upstream_estimate, y = spawner_count)) +
  geom_point() +
  facet_wrap(~stream, scales = "free") +
  theme_minimal() +
  labs(x = "Upstream Estimate",
       y = "Spawner Count",
       title = "Observed Adult Data")

# plot
observed_adult_input_wide |>
  filter(stream %in% streams_to_use) |>
  mutate(stream = str_to_title(stream)) |>
  rename(spawner_count = redd_count) |>
  ggplot(aes(x = year, y = spawner_count, color = "Spawner count (redd)")) +
  geom_line() +
  geom_line(aes(x = year, y = upstream_estimate, color = "Upstream passage")) +
  facet_wrap(~stream, scales = "free_y",
             nrow = 2) +
  xlab("Year") + ylab("Abundance") +
  ggtitle("Spawner survey data and upstream passage by stream") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        legend.position = "bottom") +
  scale_color_manual("Adult data type", values = wes_palette("GrandBudapest1")[2:3])


# do covariate selection --------------------------------------------------
compare_covariates <- compare_P2S_model_covariates()
P2S_comparison_results <- compare_covariates$covariate_comparison_results

# write for table 3
P2S_comparison_results |>
  arrange(stream, par_names) |>
  filter(par_names != "R2_data") |>
  mutate(par_names = case_when(par_names == "mean_redds_per_spawner" ~ "$\\mu_{\\delta}$",
                               par_names == "sigma_redds_per_spawner" ~ "$\\sigma_{\\delta}$",
                               par_names == "b1_survival" ~ "$\\beta_{1}$",
                               par_names == "spawner_abundance_forecast[1]" ~ "Forecasted Spawner Abundance - average",
                               par_names == "spawner_abundance_forecast[2]" ~ "Forecast Spawner Abundance - average + 1 sd",
                               par_names == "R2_fixed" ~ "$R^{2}$ of fixed effects"),
         stream = str_to_title(stream),
         mean = round(mean, 2),
         median = round(median, 2),
         sd = round(sd, 2),
         lcl = round(lcl, 2),
         ucl = round(ucl, 2),
         covar_considered = case_when(covar_considered == "wy_type" ~ "Water year type",
                                      covar_considered == "max_flow_std" ~ "Maximum flow",
                                      covar_considered == "gdd_std" ~ "Growing degree days",
                                      covar_considered == "null_covar" ~ "Null",
                                      covar_considered == "passage_index" ~ "Passage index",
                                      TRUE ~ covar_considered)) |> clipr::write_clip()

# get best covariates
# highest R2 fixed
P2S_comparison_results |>
  filter(par_names == "R2_fixed") |>
  group_by(stream) |>
  slice_max(median) |>
  select(stream, covar_considered)

# highest b1 estimate
P2S_comparison_results |>
  filter(par_names == "b1_survival") |>
  group_by(stream) |>
  slice_max(median) |>
  select(stream, covar_considered, median)

# least variance around b1 estimate
P2S_comparison_results |>
  filter(par_names == "b1_survival") |>
  group_by(stream) |>
  slice_min(sd) |>
  select(stream, covar_considered, sd)

# least variance in prediction
P2S_comparison_results |>
  filter(par_names == "spawner_abundance_forecast[1]") |>
  group_by(stream) |>
  slice_min(sd) |>
  select(stream, covar_considered, sd)

# run model ---------------------------------------------------------------

battle_P2S_results <- run_passage_to_spawner_model("battle creek",
                                                   "wy_type",
                                                   FALSE)

clear_P2S_results <- run_passage_to_spawner_model("clear creek",
                                                  "wy_type",
                                                  FALSE)
#
# deer_P2S_results <- run_passage_to_spawner_model("deer creek",
#                                                  "wy_type",
#                                                  FALSE)
#
# mill_P2S_results <- run_passage_to_spawner_model(adult_model_covariates,
#                                                  "mill creek",
#                                                  "wy_type",
#                                                  FALSE)

# join model summaries
P2S_model_fits <- bind_rows(battle_P2S_results$formatted_pars,
                            clear_P2S_results$formatted_pars)

# generate Table 2
diagnostic_pars <- c("log_mean_redds_per_spawner", "sigma_redds_per_spawner",
                     "b1_survival", "R2_fixed")
P2S_model_fits |>
  filter(par_names %in% diagnostic_pars) |>
  mutate(#across(mean:Rhat, round(digits = 2)),
         stream = str_to_title(stream)) |>
  select(Parameter = par_names, Stream = stream,
         Mean = mean, `Standard Error (mean)` = se_mean,
         `Standard Deviation` = sd, `2.5%`, `50%`, `97.5%`) |>
  clipr::write_clip()


# plots -------------------------------------------------------------------

# observed vs predicted
pred_spawners <- get_predicted_spawners_from_P2S(P2S_model_fits) |>
  glimpse()

obsv_data <- observed_adult_input_wide |>
  mutate(obsv_spawner_count = ifelse(stream == "deer creek", holding_count, redd_count)) |>
  select(year, stream, obsv_spawner_count, obsv_upstream = upstream_estimate) |>
  glimpse()

pred_with_year <- pred_spawners |>
  rename(pred_spawner_count = median_predicted_spawners) |>
  left_join(obsv_data, by = c("year", "stream")) |>
  glimpse()

pred_with_year |>
  mutate(stream = str_to_title(stream)) |>
  ggplot(aes(x = obsv_spawner_count, y = pred_spawner_count)) +
  geom_smooth(color = "#FD6467", method = "lm") +
  geom_point(alpha = 0.8) +
  facet_wrap(~stream, scales = "free", nrow = 2) +
  theme_minimal() +
  xlab("Observed Spawner Count") + ylab("Predicted Spawner Count") +
  ggtitle("Predicted vs. Observed Spawner Count")

spawner_lm <- lm(pred_spawner_count ~ obsv_spawner_count, data = pred_with_year)
summary(spawner_lm)$adj.r.squared

# conversion rate
wy_lookup <- waterYearType::water_year_indices |>
  filter(location == "Sacramento Valley") |>
  mutate(wy_type = ifelse(Yr_type %in% c("Wet", "Above Normal"), "Wet", "Dry")) |>
  select(wy_type, year = WY)

conversion_rates <- P2S_model_fits |>
  filter(str_detect(par_names, "conversion_rate")) |>
  mutate(lcl = `2.5%`, ucl = `97.5%`) |>
  select(median = `50%`, sd, lcl, ucl, stream, year) |>
  rename(pred_conversion_rate = median) |>
  left_join(wy_lookup, by = "year") |>
  mutate(stream = str_to_title(stream)) |>
  glimpse()

conversion_rates |>
  ggplot(aes(x = year, y = pred_conversion_rate)) +
  geom_errorbar(aes(x = year, ymax = ucl, ymin = lcl),
                width = 0.2, alpha = 0.7) +
  geom_line() +
  geom_point(aes(x = year, y = pred_conversion_rate,
                 color = wy_type),
             size = 4) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  facet_wrap(~stream, scales = "free", nrow = 2) +
  theme_minimal() +
  labs(x = "Year",
       y = "Predicted Conversion Rate",
       title = "Conversion rate of passage to spawners") +
  labs(color = "Water Year Type") +
  scale_color_manual(values = wes_palette("GrandBudapest1")[3:2]) +
  theme(plot.title = element_text(size = 15),
        legend.position = "bottom",
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

# model diagnostics (josh's plot)
rps <- P2S_model_fits |>
  filter(str_detect(par_names, "log_mean_redds_per_spawner")) |>
  mutate(rps = exp(`50%`),
         lcl = exp(`2.5%`),
         ucl = exp(`97.5%`),
         log_rps = `50%`)

annual_random_effects <- P2S_model_fits |>
  filter(str_detect(par_names, "log_redds_per_spawner")) |>
  mutate(annual_random_effect = `50%`) |>
  select(stream, year_index, annual_random_effect) |>
  select(-year_index)

R2_estimates <- P2S_model_fits |>
  filter(par_names %in% c("R2_data", "R2_fixed")) |>
  pivot_wider(id_cols = stream,
              names_from = par_names,
              values_from = mean)

b1_table <- P2S_model_fits |>
  filter(par_names %in% c("b1_survival")) |>
  select(stream, parameter = par_names, `2.5%`,
         `50%`, `97.5%`) |>
  left_join(R2_estimates,
            by = "stream")

b1_effect <- obsv_data |>
  drop_na(obsv_upstream, obsv_spawner_count) |>
  left_join(b1_table |>
              select(stream, b1 = `50%`), by = "stream") |>
  left_join(wy_lookup) |>
  drop_na(wy_type) |>
  mutate(covar = ifelse(wy_type == "Wet", 1, 0),
         b1_effect = (b1 * covar)) |>
  glimpse()


obsv_data |>
  filter(stream %in% streams_to_use) |>
  drop_na(obsv_upstream, obsv_spawner_count) |>
  left_join(b1_effect |>
              select(year, stream, b1_effect),
            by = c("year", "stream")) |>
  left_join(rps |>
              select(stream, rps, log_rps, lcl, ucl), by = c("stream")) |>
  mutate(pred_effect = obsv_upstream * exp(log_rps + b1_effect),
         stream = str_to_title(stream)) |>
  ggplot(aes(x = obsv_upstream, y = obsv_spawner_count)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = rps)) +
  geom_ribbon(aes(ymin = (0 + obsv_upstream * lcl),
                  ymax = (0 + obsv_upstream * ucl)),
              fill = "grey70", alpha = 0.3) +
  geom_errorbar(aes(ymin = pmin((rps * obsv_upstream), pred_effect),
                    ymax = pmax((rps * obsv_upstream), pred_effect)),
                color = "#FD6467") +
  facet_wrap(~stream, scale = "free") +
  theme_minimal() +
  xlab("Observed upstream passage") +
  ylab("Observed spawner count (redd or holding)") +
  ggtitle("Model diagnostics")

# forecast plot
forecasts <- comparison_results |>
  filter(stream %in% streams_to_use,
         str_detect(par_names, "abundance_forecast")) |>
  mutate(year_index = readr::parse_number(par_names),
         forecast_level = ifelse(year_index == 1, "Avg", "High"),
         data_type = "forecast",
         covar_considered = case_when(covar_considered == "wy_type" ~ "WY type",
                                      covar_considered == "max_flow_std" ~ "Flow",
                                      covar_considered == "gdd_std" ~ "Temp",
                                      covar_considered == "passage_index" ~ "Total Passage",
                                      covar_considered == "null_covar" ~ "Null"),
         covar_considered_f = factor(covar_considered,
                                     levels = c("Null", "WY type", "Temp", "Flow", "Total Passage"))) |>
  filter(covar_considered != "Null") |>
  select(stream, covar_considered_f, forecast_level, adult_count = median, lcl, ucl, year_index) |>
  glimpse()

forecasts |>
  mutate(stream = str_to_title(stream)) |>
  ggplot(aes(x = forecast_level, y = adult_count)) +
  geom_errorbar(aes(x = forecast_level, ymin = lcl, ymax = ucl),
                width = 0.3, alpha = 0.7) +
  geom_point(aes(x = forecast_level, y = adult_count,
                 color = forecast_level),
             size = 4) +
  facet_wrap(~stream + covar_considered_f,# scales = "free",
             nrow = 2, ncol = 4) +
  xlab("Forecast Level") + ylab("Predicted Spawner Count at Across-year Mean Passage") +
  scale_color_manual("Forecast Level",
                     values = wes_palette("GrandBudapest1")[2:3]) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = "bottom",
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45))


# diagnostics -------------------------------------------------------------

# battle creek
mcmc_pairs(battle_P2S_results$full_object,
           pars = c("b1_survival", "mean_redds_per_spawner"))

mcmc_trace(battle_P2S_results$full_object,
           pars = c("b1_survival", "mean_redds_per_spawner"))
# mcmc_areas(battle_P2S_results$full_object,
#            pars = c("b1_survival", "mean_redds_per_spawner"))
# mcmc_rhat(rhat(battle_P2S_results$full_object))

# clear creek
mcmc_pairs(clear_P2S_results$full_object,
           pars = c("b1_survival", "mean_redds_per_spawner"))

mcmc_trace(clear_P2S_results$full_object,
           pars = c("b1_survival", "mean_redds_per_spawner"))


