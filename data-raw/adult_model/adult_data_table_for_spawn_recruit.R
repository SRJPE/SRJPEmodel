# put together adult data table for josh


# libraries and read in predicted model data ------------------------------
library(ggplot2)
library(tidyverse)
library(tidybayes)

load(here::here("data-raw", "adult_model", "adult_data_objects.Rdata"))
load(here::here("data-raw", "adult_model", "best_adult_models.Rdata"))
load(here::here("data-raw", "adult_model", "model_fit_summaries.Rdata"))


# pull in predicted values from model -------------------------------------
battle_years <- adult_data_objects$survival_model_data |>
  filter(stream == "battle creek") |>
  drop_na(water_year_type) |>
  pull(year)
clear_years <- adult_data_objects$survival_model_data |>
  filter(stream == "clear creek") |>
  drop_na(max_flow) |>
  pull(year)
mill_years <- adult_data_objects$survival_model_data |>
  filter(stream == "mill creek") |>
  drop_na(gdd_total) |>
  pull(year)
deer_years <- adult_data_objects$survival_model_data |>
  filter(stream == "deer creek") |>
  drop_na(water_year_type) |>
  pull(year)

adult_data_all_streams <- bind_rows(model_fit_summaries$battle_pred |>
                                      filter(par_names != "b1_survival") |>
                                      mutate(year = battle_years,
                                             data_type = "modeled_spawners",
                                             stream = "battle creek") |>
                                      select(year, spawner_count = mean, data_type, stream),
                                    model_fit_summaries$clear_pred |>
                                      filter(par_names != "b1_survival") |>
                                      mutate(year = clear_years,
                                             data_type = "modeled_spawners",
                                             stream = "clear creek") |>
                                      select(year, spawner_count = mean, data_type, stream),
                                    model_fit_summaries$mill_pred |>
                                      filter(par_names != "b1_survival") |>
                                      mutate(year = mill_years,
                                             data_type = "modeled_spawners",
                                             stream = "mill creek") |>
                                      select(year, spawner_count = mean, data_type, stream),
                                    model_fit_summaries$deer_pred |>
                                      filter(par_names != "b1_survival") |>
                                      mutate(year = deer_years,
                                             data_type = "modeled_spawners",
                                             stream = "deer creek") |>
                                      select(year, spawner_count = mean, data_type, stream),
                                    adult_data_objects$yuba_data |>
                                        mutate(stream = "yuba river",
                                               data_type = "upstream_passage_estimate") |>
                                        rename(spawner_count = passage_estimate),
                                    adult_data_objects$feather_data |>
                                      mutate(stream = "feather river",
                                             data_type = "carcass_cjs_estimate") |>
                                      rename(spawner_count = spawner_estimate),
                                    adult_data_objects$butte_data |>
                                        mutate(stream = "butte creek",
                                               data_type = "carcass_cjs_estimate") |>
                                        rename(spawner_count = spawner_estimate)) |>
  glimpse()

write.csv(adult_data_all_streams, here::here("data-raw", "adult_model", "adult_data_all_streams.csv"),
          row.names = FALSE)
