# diagnostics
# libraries and read in predicted model data ------------------------------
library(ggplot2)
library(googleCloudStorageR)
library(tidyverse)


# pull in data from google cloud ------------------------------------------
gcs_auth(json_file = Sys.getenv("GCS_AUTH_FILE"))
# Set global bucket
gcs_global_bucket(bucket = Sys.getenv("GCS_DEFAULT_BUCKET"))

gcs_get_object(object_name = "jpe-model-data/adult-model/survival_model_data_raw.csv",
                bucket = gcs_get_global_bucket(),
                saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
                                        "survival_model_data_raw.csv"),
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

survival_model_data_raw <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                               "survival_model_data_raw.csv"))

full_data_for_input <- survival_model_data_raw |>
  mutate(wy_type_std = ifelse(water_year_type == "dry", 0, 1),
         wy_type_std = as.vector(scale(wy_type_std)), # TODO double check this. was producing negative b1_surv estimates for wet year types
         max_flow_std = as.vector(scale(max_flow)),
         gdd_std = as.vector(scale(gdd_total)),
         min_passage_timing_std = as.vector(scale(min_passage_timing))) |>
  select(year, stream, upstream_count, redd_count, holding_count,
         wy_type_std, max_flow_std, gdd_std, min_passage_timing_std) |>
  glimpse()

report_pars <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                   "model_fit_summaries.csv"))
diagnostic_pars <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                       "model_fit_diagnostic_pars.csv"))



# plot predicted spawner estimates against observed spawner values

diagnostics <- report_pars |>
  filter(par_names != "b1_survival") |>
  select(par_names, mean, sd, stream) |>
  glimpse()

obsv_data <- full_data_for_input |>
  mutate(obsv_spawner_count = ifelse(stream == "deer creek", holding_count, redd_count)) |>
  select(year, stream, obsv_spawner_count) |>
  glimpse()

years_to_join <- c(full_data_for_input |> filter(stream == "battle creek") |> drop_na(wy_type_std) |> pull(year),
                   full_data_for_input |> filter(stream == "clear creek") |> drop_na(max_flow_std) |> pull(year),
                   full_data_for_input |> filter(stream == "deer creek") |> drop_na(wy_type_std) |> pull(year),
                   full_data_for_input |> filter(stream == "mill creek") |> drop_na(gdd_std) |> pull(year)) |>
  glimpse()

diagnostics_with_year <- diagnostics |>
  arrange(stream) |>
  mutate(year = years_to_join) |>
  select(-par_names) |>
  rename(pred_spawner_count = mean) |>
  left_join(obsv_data, by = c("year", "stream")) |>
  glimpse()

diagnostics_with_year |> ggplot(aes(x = obsv_spawner_count, y = pred_spawner_count)) +
  geom_point() + geom_smooth(method = "lm") + facet_wrap(~stream, scales = "free") +
  theme_minimal() + xlab("Observed Spawner Count") + ylab("Predicted Spawner Count")

# look at R2
summary(lm(pred_spawner_count ~ obsv_spawner_count, data = diagnostics_with_year))

# look at estimated value of sigma redds_per_spawner (mean and sd)
diagnostic_pars |>
  filter(str_detect(par_names, "log_redds_per_spawner") |
         str_detect(par_names, "survival_rate")) |>
  mutate(par_name = ifelse(str_detect(par_names, "log_redds_per_spawner"), "log_rps", "survival")) |>
  arrange(par_name) |>
  mutate(year = rep(years_to_join, 2)) |>
  ggplot(aes(x = year, y = mean, color = par_name)) +
  geom_point() +
  scale_color_discrete(name = "Parameter name",
                       labels = c("Log redds/spawner (RE)",
                                  "Predicted survival rate")) +
  facet_wrap(~stream, scales = "free") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  xlab("Year") + ylab("Mean estimated value") + theme_minimal()  +
  theme(legend.position = "bottom")

# plot year-specific random effects (log_redds_per_spawner[y])
# a big portion of the variation in predicted survival_rate[y] comes from b1_survival * environmental_covar[y]
# and not random effects


# scratch -----------------------------------------------------------------

get_other_pars <- function(model_fit, stream_name) {
  par_results <- summary(model_fit)$summary
  results_tibble <- as.data.frame(par_results) |>
    rownames_to_column("par_names") |>
    filter(!str_detect(par_names, "predicted_spawners"),
           par_names != "b1_survival") |>
    mutate(stream = stream_name)

  return(results_tibble)
}
