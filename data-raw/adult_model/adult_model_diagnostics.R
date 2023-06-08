# diagnostics
# libraries and read in predicted model data ------------------------------
library(ggplot2)
library(googleCloudStorageR)
library(tidyverse)


# pull in data from google cloud ------------------------------------------
gcs_auth(json_file = Sys.getenv("GCS_AUTH_FILE"))
# Set global bucket
gcs_global_bucket(bucket = Sys.getenv("GCS_DEFAULT_BUCKET"))

# download input data
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

gcs_get_object(object_name = "jpe-model-data/adult-model/model_fit_summaries.csv",
               bucket = gcs_get_global_bucket(),
               saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
                                       "model_fit_summaries.csv"),
               overwrite = TRUE)
gcs_get_object(object_name = "jpe-model-data/adult-model/model_fit_diagnostic_pars.csv",
               bucket = gcs_get_global_bucket(),
               saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
                                       "model_fit_diagnostic_pars.csv"),
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

report_pars <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                   "model_fit_summaries.csv"))
other_pars <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                       "model_fit_diagnostic_pars.csv"))



# plot predicted spawner estimates against observed spawner values

pred_spawners <- report_pars |>
  filter(par_names != "b1_survival") |>
  select(par_names, mean, lcl = X2.5., ucl = X97.5., sd, stream) |>
  glimpse()

obsv_data <- full_data_for_input |>
  mutate(obsv_spawner_count = ifelse(stream == "deer creek", holding_count, redd_count)) |>
  select(year, stream, obsv_spawner_count, obsv_upstream = upstream_estimate) |>
  glimpse()

years_to_join <- c(full_data_for_input |> filter(stream == "battle creek") |> drop_na(upstream_estimate, redd_count, wy_type) |> pull(year),
                   full_data_for_input |> filter(stream == "clear creek") |> drop_na(upstream_estimate, redd_count, max_flow_std) |> pull(year),
                   full_data_for_input |> filter(stream == "deer creek") |> drop_na(upstream_estimate, holding_count, wy_type) |> pull(year),
                   full_data_for_input |> filter(stream == "mill creek") |> drop_na(upstream_estimate, redd_count, gdd_std) |> pull(year))

pred_with_year <- pred_spawners |>
  arrange(stream) |>
  mutate(year = years_to_join) |>
  select(-par_names) |>
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
other_pars |>
  filter(str_detect(par_names, "log_redds_per_spawner") |
         str_detect(par_names, "conversion_rate")) |>
  mutate(par_name = ifelse(str_detect(par_names, "log_redds_per_spawner"), "log_rps", "survival")) |>
  arrange(par_name) |>
  mutate(year = rep(years_to_join, 2)) |>
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

rps <- other_pars |>
  filter(str_detect(par_names, "log_mean_redds_per_spawner")) |>
  mutate(rps = exp(mean),
         lcl = exp(X2.5.),
         ucl = exp(X97.5.))

obsv_data |>
  filter(!stream %in% c("feather river", "butte creek", "yuba river")) |>
  drop_na(obsv_upstream) |>
  left_join(rps |>
              select(stream, rps, lcl, ucl), by = c("stream")) |>
  ggplot(aes(x = obsv_upstream, y = obsv_spawner_count)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_abline(aes(intercept = 0, slope = rps)) +
  geom_ribbon(aes(ymin = (0 + obsv_upstream * lcl),
                  ymax = (0 + obsv_upstream * ucl)),
              fill = "grey70", alpha = 0.3) +
  facet_wrap(~stream, scale = "free") +
  theme_minimal() +
  xlab("Observed upstream passage") +
  ylab("Observed spawner count (redd or holding)") +
  ggtitle("Model diagnostics")

# write b1_survival stats -------------------------------------------------
b1_table <- report_pars |>
  filter(par_names %in% c("b1_survival")) |>
  select(stream, parameter = par_names, `2.5` = X2.5.,
         `50` = mean, `97.5` = X97.5.) |>
  left_join(other_pars |>
              filter(par_names %in% c("R2_data", "R2_fixed")) |>
              pivot_wider(id_cols = stream,
                          names_from = par_names,
                          values_from = mean),
            by = "stream")

b1_effect <- obsv_data |>
  drop_na(obsv_upstream, obsv_spawner_count) |>
  left_join(b1_table |>
              select(stream, b1 = `50`), by = "stream") |>
  left_join(full_data_for_input |>
              mutate(covar = case_when(stream %in% c("battle creek", "deer creek") ~ wy_type,
                                       stream == "mill creek" ~ gdd_std,
                                       stream == "clear creek" ~ max_flow_std)) |>
              select(year, stream, covar)) |>
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
  geom_errorbar(aes(ymin = pmin(b1_effect_pred, random_effect_pred),
                    ymax = pmax(b1_effect_pred, random_effect_pred)),
                color = "red") +
  facet_wrap(~stream, scale = "free") +
  theme_minimal() +
  xlab("Observed upstream passage") +
  ylab("Observed spawner count (redd or holding)") +
  ggtitle("Model diagnostics")

