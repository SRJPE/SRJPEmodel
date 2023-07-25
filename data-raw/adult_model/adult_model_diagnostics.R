# diagnostics
# libraries and read in predicted model data ------------------------------
library(ggplot2)
library(googleCloudStorageR)
library(tidyverse)
library(wesanderson)


# pull in data from google cloud ------------------------------------------
gcs_auth(json_file = Sys.getenv("GCS_AUTH_FILE"))
# Set global bucket
gcs_global_bucket(bucket = Sys.getenv("GCS_DEFAULT_BUCKET"))

# download (refresh) input data

gcs_get_object(object_name = "jpe-model-data/adult-model/P2S_model_fits.csv",
               bucket = gcs_get_global_bucket(),
               saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
                                       "P2S_model_fits.csv"),
               overwrite = TRUE)
gcs_get_object(object_name = "jpe-model-data/adult-model/adult_data_input_raw.csv",
               bucket = gcs_get_global_bucket(),
               saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
                                       "adult_data_input_raw.csv"),
               overwrite = TRUE)
# download covariate data
gcs_get_object(object_name = "jpe-model-data/adult-model/adult_model_covariates_standard.csv",
               bucket = gcs_get_global_bucket(),
               saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
                                       "adult_model_covariates_standard.csv"),
               overwrite = TRUE)

# read data from csvs -----------------------------------------------------

adult_data_input_raw <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                            "adult_data_input_raw.csv"))

adult_model_covariates <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                              "adult_model_covariates_standard.csv"))

full_data_for_input <- full_join(adult_data_input_raw,
                                 adult_model_covariates,
                                 by = c("year", "stream")) |>
  filter(!is.na(data_type)) |>
  pivot_wider(id_cols = c(year, stream, wy_type, max_flow_std, gdd_std,
                          median_passage_timing_std),
              names_from = data_type,
              values_from = count) |>
  glimpse()

P2S_model_fits <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                      "P2S_model_fits.csv")) |>
  separate(par_names, into = c("par_names", "year_index"), sep = "\\[") |>
  mutate(year_index = str_remove(year_index, "\\]")) |>
  glimpse()


# plot predicted spawner estimates against observed spawner values
pred_spawners <- P2S_model_fits |>
  filter(str_detect(par_names, "predicted_spawners")) |>
  select(par_names, mean, year_index, lcl = X2.5., ucl = X97.5., sd, stream) |>
  glimpse()

obsv_data <- full_data_for_input |>
  mutate(obsv_spawner_count = ifelse(stream == "deer creek", holding_count, redd_count)) |>
  select(year, stream, obsv_spawner_count, obsv_upstream = upstream_estimate) |>
  glimpse()

years_to_join <- P2S_model_fits |>
  filter(str_detect(par_names, "year")) |>
  select(year = mean, stream, year_index) |>
  glimpse()

pred_with_year <- pred_spawners |>
  left_join(years_to_join, by = c("year_index", "stream")) |>
  select(-c(par_names, year_index)) |>
  rename(pred_spawner_count = mean) |>
  left_join(obsv_data, by = c("year", "stream")) |>
  glimpse()

pred_with_year |> ggplot(aes(x = obsv_spawner_count, y = pred_spawner_count)) +
  geom_point() + geom_smooth(method = "lm") + facet_wrap(~stream, scales = "free") +
  theme_minimal() + xlab("Observed Spawner Count") + ylab("Predicted Spawner Count") +
  ggtitle("Predicted vs. Observed Spawner Count")

# look at R2
summary(lm(pred_spawner_count ~ obsv_spawner_count, data = pred_with_year))

# look at estimated value of sigma redds_per_spawner (mean and sd)
P2S_model_fits |>
  filter(str_detect(par_names, "log_redds_per_spawner") |
         str_detect(par_names, "conversion_rate")) |>
  mutate(par_name = ifelse(str_detect(par_names, "log_redds_per_spawner"), "log_rps", "survival")) |>
  left_join(years_to_join, by = c("year_index", "stream")) |>
  ggplot(aes(x = year, y = mean, color = par_name)) +
  geom_point() +
  scale_color_discrete(name = "Parameter name",
                       labels = c("Log redds/spawner (RE)",
                                  "Predicted conversion rate")) +
  facet_wrap(~stream, scales = "free") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  xlab("Year") + ylab("Mean estimated value") + theme_minimal()  +
  theme(legend.position = "bottom") +
  ggtitle("Year-specific random effects vs. predicted conversion rate")


# diagnostic plots --------------------------------------------------------
rps <- P2S_model_fits |>
  filter(str_detect(par_names, "log_mean_redds_per_spawner")) |>
  mutate(rps = exp(mean),
         lcl = exp(X2.5.),
         ucl = exp(X97.5.))

# write b1_survival stats -------------------------------------------------
R2_estimates <- P2S_model_fits |>
  filter(par_names %in% c("R2_data", "R2_fixed")) |>
  pivot_wider(id_cols = stream,
              names_from = par_names,
              values_from = mean)

b1_table <- P2S_model_fits |>
  filter(par_names %in% c("b1_survival")) |>
  select(stream, parameter = par_names, `2.5` = X2.5.,
         `50` = mean, `97.5` = X97.5.) |>
  left_join(R2_estimates,
            by = "stream")


# plot contribution of random effects -------------------------------------
b1_effect <- obsv_data |>
  drop_na(obsv_upstream, obsv_spawner_count) |>
  left_join(b1_table |>
              select(stream, b1 = `50`), by = "stream") |>
  left_join(full_data_for_input |>
              mutate(covar = wy_type) |>
              # mutate(covar = case_when(stream %in% c("battle creek", "clear creek", "mill creek") ~ gdd_std,
              #                          stream == "deer creek" ~ max_flow_std)) |>
              select(year, stream, covar)) |>
  drop_na(covar) |>
  #mutate(b1_effect = (b1 * covar))
  mutate(b1_effect = obsv_upstream * (b1 * covar))

obsv_data |>
  filter(!stream %in% c("feather river", "butte creek", "yuba river")) |>
  drop_na(obsv_upstream) |>
  left_join(b1_effect |>
              select(year, stream, b1_effect),
            by = c("year", "stream")) |>
  left_join(rps |>
              select(stream, rps, lcl, ucl), by = c("stream")) |>
  mutate(b1_effect_pred = (rps * obsv_upstream + b1_effect),
         random_effect_pred = (rps * obsv_upstream)) |>
  ggplot(aes(x = obsv_upstream, y = obsv_spawner_count)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = rps)) +
  geom_ribbon(aes(ymin = (0 + obsv_upstream * lcl),
                  ymax = (0 + obsv_upstream * ucl)),
              fill = "grey70", alpha = 0.3) +
  # geom_errorbar(aes(ymin = pmin(b1_effect_pred, random_effect_pred),
  #                   ymax = pmax(b1_effect_pred, random_effect_pred)),
  #               color = "red") +
  facet_wrap(~stream, scale = "free") +
  theme_minimal() +
  xlab("Observed upstream passage") +
  ylab("Observed spawner count (redd or holding)") +
  ggtitle("Model diagnostics")


# plot conversion rates ---------------------------------------------------
conversion_rates <- P2S_model_fits |>
  filter(str_detect(par_names, "conversion_rate")) |>
  mutate(lcl = `X2.5.`, ucl = `X97.5.`) |>
  select(par_names, mean, sd, lcl, ucl, stream, year_index) |>
  glimpse()

conversion_rates_with_year <- conversion_rates |>
  left_join(years_to_join, by = c("year_index", "stream")) |>
  select(-c(par_names, year_index)) |>
  rename(pred_conversion_rate = mean) |>
  left_join(full_data_for_input |>
              mutate(covar_used = wy_type,
                     covar_type = "Water year type") |>
              # mutate(covar_used = ifelse(stream == "deer creek", max_flow_std, gdd_std),
              #        covar_type = ifelse(stream == "deer creek", "Maximum flow", "Temperature index")) |>
              select(year, stream, covar_used, covar_type),
            by = c("year", "stream")) |>
  mutate(stream = str_to_title(stream),
         stream = ifelse(stream == "Deer Creek", paste0(stream, "- Holding"), paste(stream, "- Redd"))) |>
  glimpse()

conversion_rate_plot <- conversion_rates_with_year |>
  ggplot(aes(x = year, y = pred_conversion_rate)) +
  geom_line() +
  geom_line(aes(x = year, y = covar_used, color = covar_type), linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  facet_wrap(~stream, scales = "free") +
  theme_minimal() + xlab("Year") +
  ylab("Predicted Conversion Rate") +
  ggtitle("Conversion Rate of Passage to Spawners") +
  labs(color = "Covariate Type") +
  scale_color_manual(values = wes_palette("GrandBudapest1")[2:3]) +
  theme(plot.title = element_text(size = 15),
        legend.position = "bottom",
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

ggsave(filename = here::here("data-raw", "adult_model",
                             "adult_model_plots", "conversion_rate_plot.jpg"),
       plot = conversion_rate_plot, width = 12, height = 6)


# standardized conversion rate plot ---------------------------------------

conversion_rate_plot_battle <- conversion_rates_with_year |>
  filter(stream == "Battle Creek - Redd") |>
  mutate(covar_used_scaled = scales::rescale(covar_used, to = c(0, 0.5))) |>
  ggplot(aes(x = year, y = pred_conversion_rate)) +
  geom_line() +
  geom_line(aes(x = year, y = covar_used_scaled, color = covar_type), linetype = "dashed") +
  theme_minimal() + xlab("Year") +
  ylab("Predicted Conversion Rate") +
  ggtitle("Conversion Rate of Passage to Spawners\nBattle Creek") +
  scale_color_manual("Covariate type", values = wes_palette("GrandBudapest1")[2:3]) +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = "bottom",
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

ggsave(filename = here::here("data-raw", "adult_model",
                             "adult_model_plots", "conversion_rate_plot_battle.jpg"),
       plot = conversion_rate_plot_battle, width = 12, height = 7)


# forecasts ---------------------------------------------------------------

yuba_data <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                 "yuba_data.csv"))
butte_data <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                  "butte_data.csv"))
feather_data <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                    "feather_data.csv"))


forecasts <- P2S_model_fits |>
  filter(str_detect(par_names, "abundance_forecast")) |>
  rename(lcl = `X2.5.`, ucl = `X97.5.`) |>
  mutate(par_names = ifelse(year_index == "1",
                            "forecast_dry_year",
                            "forecast_wet_year"),
  # mutate(par_names = ifelse(par_names == "spawner_abundance_forecast[1]",
  #                           "forecast_avg_conditions",
  #                           "forecast_high_conditions"),
         forecast_type = ifelse(par_names == "forecast_dry_year",
                                "Dry", "Wet")) |>
         # forecast_type = ifelse(par_names == "forecast_avg_conditions",
         #                        "Average", "High")) |>
  select(stream, forecast_type, adult_count = mean, lcl, ucl) |> glimpse()
  # pivot_wider(id_cols = c("stream"),
  #             names_from = "par_names",
  #             values_from = c("mean", "lcl", "ucl"))

all_data_sources <- P2S_model_fits |>
  filter(str_detect(par_names, "predicted_spawner")) |>
  left_join(years_to_join, by = c("year_index", "stream")) |>
  select(stream, adult_count = mean,
         lcl = `X2.5.`, ucl = `X97.5.`, year) |>
  mutate(data_type = "P2S spawners estimates") |>
  bind_rows(yuba_data |>
              rename(adult_count = passage_estimate) |>
              mutate(data_type = "passage estimate",
                     stream = "yuba river"),
            butte_data |>
              rename(adult_count = spawner_estimate) |>
              mutate(data_type = "CJS spawner estimate",
                     stream = "butte creek"),
            feather_data |>
              rename(adult_count = spawner_estimate) |>
              mutate(data_type = "CJS spawner estimate",
                     stream = "feather river")) |>
  glimpse()

adult_data_source_plot <- all_data_sources |>
  mutate(stream = str_to_title(stream)) |>
  ggplot(aes(x = year, y = adult_count)) +
  geom_line(aes(color = data_type)) +
  labs(color = "Data Type") +
  geom_ribbon(aes(x = year, ymin = lcl, ymax = ucl), alpha = 0.2) +
  facet_wrap(~ stream, scales = "free", nrow = 2) +
  theme_minimal() +
  xlab("Year") + ylab("Spawner Count") +
  ggtitle("Spawner Counts by Method") +
  scale_color_manual("Covariate type", values = wes_palette("GrandBudapest1")) +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = "bottom",
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(angle = 45))

ggsave(filename = here::here("data-raw", "adult_model",
                             "adult_model_plots", "all_adult_data_sources.jpg"),
       plot = adult_data_source_plot, width = 12, height = 7)


# forecasting plot for one stream -----------------------------------------
forecasts_battle <- forecasts |>
  filter(stream == "battle creek") |>
  mutate(year = 2021,
         data_type = "forecast")

forecast_plot <- all_data_sources |>
  filter(stream == "battle creek") |>
  ggplot(aes(x = year, y = adult_count)) +
  geom_line() +
  geom_point(aes(x = year, y = adult_count,
                 shape = forecast_type, color = forecast_type,
                 size = 0.3),
             data = forecasts_battle) +
  geom_ribbon(aes(x = year, ymin = lcl, ymax = ucl), alpha = 0.2) +
  geom_errorbar(aes(x = year, ymin = lcl, ymax = ucl,
                    color = forecast_type),
                data = forecasts_battle) +
  xlab("Year") + ylab("Predicted Spawner Count") +
  scale_color_manual("Forecast type",
                     values = wes_palette("GrandBudapest1")[2:3]) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = "bottom",
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45)) +
  ggtitle("Forecasted Spawners - Battle Creek")

ggsave(filename = here::here("data-raw", "adult_model",
                             "adult_model_plots", "forecast_plot.jpg"),
       plot = forecast_plot, width = 12, height = 7)

