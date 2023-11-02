# put together adult data table for josh

# TODO add in uncertainty here!
# libraries and read in predicted model data ------------------------------
library(ggplot2)
library(tidyverse)
library(tidybayes)
library(googleCloudStorageR)


# pull in data from google cloud ------------------------------------------
gcs_auth(json_file = Sys.getenv("GCS_AUTH_FILE"))
# Set global bucket
gcs_global_bucket(bucket = Sys.getenv("GCS_DEFAULT_BUCKET"))
# get data and save as xlsx
gcs_get_object(object_name = "jpe-model-data/adult-model/survival_model_data.csv",
                                      bucket = gcs_get_global_bucket(),
                                      saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
                                                              "survival_model_data.csv"),
                                      overwrite = TRUE)
gcs_get_object(object_name = "jpe-model-data/adult-model/yuba_data.csv",
                            bucket = gcs_get_global_bucket(),
                            saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
                                                    "yuba_data.csv"),
                            overwrite = TRUE)
gcs_get_object(object_name = "jpe-model-data/adult-model/butte_data.csv",
                             bucket = gcs_get_global_bucket(),
                             saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
                                                     "butte_data.csv"),
                             overwrite = TRUE)
gcs_get_object(object_name = "jpe-model-data/adult-model/feather_data.csv",
                               bucket = gcs_get_global_bucket(),
                               saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
                                                       "feather_data.csv"),
                               overwrite = TRUE)

gcs_get_object(object_name = "jpe-model-data/adult-model/adult_data_input_raw.csv",
               bucket = gcs_get_global_bucket(),
               saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
                                       "adult_data_input_raw.csv"),
               overwrite = TRUE)

gcs_get_object(object_name = "jpe-model-data/adult-model/P2S_model_fits.csv",
               bucket = gcs_get_global_bucket(),
               saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
                             "P2S_model_fits.csv"),
               overwrite = TRUE)

survival_model_data <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                               "survival_model_data.csv"))
adult_data_input_raw <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                            "adult_data_input_raw.csv"))
# we want all data, not just the "best"
# yuba_data <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
#                                                "yuba_data.csv"))
# butte_data <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
#                                                "butte_data.csv"))
# feather_data <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
#                                                "feather_data.csv"))

P2S_model_fits <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                      "P2S_model_fits.csv")) |>
  separate(par_names, into = c("par_names", "year_index"), sep = "\\[") |>
  mutate(year_index = str_remove(year_index, "\\]")) |>
  glimpse()

years_to_join <- P2S_model_fits |>
  filter(str_detect(par_names, "year")) |>
  select(year = mean, stream, year_index) |>
  glimpse()

P2S_model_fits_with_year <- P2S_model_fits |>
  filter(str_detect(par_names, "predicted_spawners")) |>
  select(P2S_median_count = X50., year_index, P2S_lcl = X2.5., P2S_ucl = X97.5., stream) |>
  left_join(years_to_join, by = c("year_index", "stream")) |>
  select(-c(year_index)) |>
  relocate(year, .before = P2S_median_count) |>
  glimpse()

# pull in predicted values from model -------------------------------------

adult_data_all_streams <- bind_rows(P2S_model_fits |>
                                      filter(str_detect(par_names, "predicted_spawners"),
                                             stream == "battle creek") |>
                                      select(-stream) |>
                                      mutate(year = battle_years,
                                             data_type = "modeled_spawners",
                                             stream = "battle creek") |>
                                      select(year, spawner_count = mean, data_type, stream),
                                    P2S_model_fits |>
                                      filter(str_detect(par_names, "predicted_spawners"),
                                             stream == "clear creek") |>
                                      select(-stream) |>
                                      mutate(year = clear_years,
                                             data_type = "modeled_spawners",
                                             stream = "clear creek") |>
                                      select(year, spawner_count = mean, data_type, stream),
                                    P2S_model_fits |>
                                      filter(str_detect(par_names, "predicted_spawners"),
                                             stream == "mill creek") |>
                                      select(-stream) |>
                                      mutate(year = mill_years,
                                             data_type = "modeled_spawners",
                                             stream = "mill creek") |>
                                      select(year, spawner_count = mean, data_type, stream),
                                    P2S_model_fits |>
                                      filter(str_detect(par_names, "predicted_spawners"),
                                             stream == "deer creek") |>
                                      select(-stream) |>
                                      mutate(year = deer_years,
                                             data_type = "modeled_spawners",
                                             stream = "deer creek") |>
                                      select(year, spawner_count = mean, data_type, stream),
                                    yuba_data |>
                                        mutate(stream = "yuba river",
                                               data_type = "upstream_passage_estimate") |>
                                        rename(spawner_count = passage_estimate),
                                    feather_data |>
                                      mutate(stream = "feather river",
                                             data_type = "carcass_cjs_estimate") |>
                                      rename(spawner_count = spawner_estimate),
                                    butte_data |>
                                        mutate(stream = "butte creek",
                                               data_type = "carcass_cjs_estimate") |>
                                        rename(spawner_count = spawner_estimate)) |>
  glimpse()


# upload to google cloud --------------------------------------------------

f <- function(input, output) write_csv(input, file = output)

gcs_upload(adult_data_all_streams,
           object_function = f,
           type = "csv",
           name = "jpe-model-data/adult-model/adult_data_for_spawn_recruit.csv")

write.csv(adult_data_all_streams, here::here("data-raw", "adult_model", "adult_model_data",
                                             "adult_data_for_spawn_recruit.csv"),
          row.names = FALSE)
