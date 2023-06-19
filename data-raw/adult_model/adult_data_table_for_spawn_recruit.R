# put together adult data table for josh


# libraries and read in predicted model data ------------------------------
library(ggplot2)
library(tidyverse)
library(tidybayes)
library(googleCloudStorageR)


# pull in data from google cloud ------------------------------------------
gcs_auth(json_file = Sys.getenv("GCS_AUTH_FILE"))
# Set global bucket
gcs_global_bucket(bucket = Sys.getenv("GCS_DEFAULT_BUCKET"))
# git data and save as xlsx
survival_model_data_raw <- gcs_get_object(object_name = "jpe-model-data/adult-model/survival_model_data_raw.csv",
                                          bucket = gcs_get_global_bucket(),
                                          saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
                                                                  "survival_model_data_raw.csv"),
                                          overwrite = TRUE)
survival_model_data <- gcs_get_object(object_name = "jpe-model-data/adult-model/survival_model_data.csv",
                                      bucket = gcs_get_global_bucket(),
                                      saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
                                                              "survival_model_data.csv"),
                                      overwrite = TRUE)
battle_data <- gcs_get_object(object_name = "jpe-model-data/adult-model/battle_data.csv",
                              bucket = gcs_get_global_bucket(),
                              saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
                                                      "battle_data.csv"),
                              overwrite = TRUE)
clear_data <- gcs_get_object(object_name = "jpe-model-data/adult-model/clear_data.csv",
                             bucket = gcs_get_global_bucket(),
                             saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
                                                     "clear_data.csv"),
                             overwrite = TRUE)
mill_data <- gcs_get_object(object_name = "jpe-model-data/adult-model/mill_data.csv",
                            bucket = gcs_get_global_bucket(),
                            saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
                                                    "mill_data.csv"),
                            overwrite = TRUE)
deer_data <- gcs_get_object(object_name = "jpe-model-data/adult-model/deer_data.csv",
                            bucket = gcs_get_global_bucket(),
                            saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
                                                    "deer_data.csv"),
                            overwrite = TRUE)
yuba_data <- gcs_get_object(object_name = "jpe-model-data/adult-model/yuba_data.csv",
                            bucket = gcs_get_global_bucket(),
                            saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
                                                    "yuba_data.csv"),
                            overwrite = TRUE)
butte_data <- gcs_get_object(object_name = "jpe-model-data/adult-model/butte_data.csv",
                             bucket = gcs_get_global_bucket(),
                             saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
                                                     "butte_data.csv"),
                             overwrite = TRUE)
feather_data <- gcs_get_object(object_name = "jpe-model-data/adult-model/feather_data.csv",
                               bucket = gcs_get_global_bucket(),
                               saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
                                                       "feather_data.csv"),
                               overwrite = TRUE)

# model_fit_summaries <- gcs_get_object(object_name = "jpe-model-data/adult-model/model_fit_summaries.csv",
#                                  bucket = gcs_get_global_bucket(),
#                                  saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
#                                                          "model_fit_summaries.csv"),
#                                  overwrite = TRUE)

gcs_get_object(object_name = "jpe-model-data/adult-model/model_fit_summaries.csv",
               bucket = gcs_get_global_bucket(),
               saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
                             "P2S_model_fits.csv"),
               overwrite = TRUE)

survival_model_data_raw <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                               "survival_model_data_raw.csv"))
survival_model_data <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                               "survival_model_data.csv"))
battle_data <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                               "battle_data.csv"))
clear_data <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                               "clear_data.csv"))
deer_data <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                               "deer_data.csv"))
mill_data <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                               "mill_data.csv"))
yuba_data <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                               "yuba_data.csv"))
butte_data <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                               "butte_data.csv"))
feather_data <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                               "feather_data.csv"))
# model_fit_summaries <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
#                                     "model_fit_summaries.csv"))

P2S_model_fits <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                      "P2S_model_fits.csv")) |>
  glimpse()


# pull in predicted values from model -------------------------------------
battle_years <- survival_model_data |>
  filter(stream == "battle creek") |>
  drop_na(gdd_total) |>
  pull(year)
clear_years <- survival_model_data |>
  filter(stream == "clear creek",
         year >= 2003) |>
  drop_na(gdd_total) |>
  pull(year)
mill_years <- survival_model_data |>
  filter(stream == "mill creek") |>
  drop_na(gdd_total) |>
  pull(year)
deer_years <- survival_model_data |>
  filter(stream == "deer creek") |>
  drop_na(max_flow) |>
  pull(year)

adult_data_all_streams <- bind_rows(P2S_model_fits |>
                                      filter(par_names != "b1_survival",
                                             stream == "battle creek") |>
                                      select(-stream) |>
                                      mutate(year = battle_years,
                                             data_type = "modeled_spawners",
                                             stream = "battle creek") |>
                                      select(year, spawner_count = mean, data_type, stream),
                                    P2S_model_fits |>
                                      filter(par_names != "b1_survival",
                                             stream == "clear creek") |>
                                      select(-stream) |>
                                      mutate(year = clear_years,
                                             data_type = "modeled_spawners",
                                             stream = "clear creek") |>
                                      select(year, spawner_count = mean, data_type, stream),
                                    P2S_model_fits |>
                                      filter(par_names != "b1_survival",
                                             stream == "mill creek") |>
                                      select(-stream) |>
                                      mutate(year = mill_years,
                                             data_type = "modeled_spawners",
                                             stream = "mill creek") |>
                                      select(year, spawner_count = mean, data_type, stream),
                                    P2S_model_fits |>
                                      filter(par_names != "b1_survival",
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
