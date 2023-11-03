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
# get data and save as xlsx
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

# read in data and format correctly ---------------------------------------

# observed data
adult_data_input_raw <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                            "adult_data_input_raw.csv")) |>
  mutate(count = as.numeric(count))

# predicted data from P2S
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
  pivot_longer(c(P2S_median_count, P2S_lcl, P2S_ucl),
               names_to = "data_type", values_to = "count") |> # TODO to make it long, do it here
  glimpse()

# bind together the two groups of data -------------------------------------

# pivot wider so that we have a column for each data type
table_for_spawn_recruit <- bind_rows(P2S_model_fits_with_year,
                                     adult_data_input_raw) |>
  mutate(count = round(count, 0)) |>
  pivot_wider(id_cols = c(year, stream),
              names_from = data_type,
              values_from = count) |>
  mutate(lcl = case_when(!is.na(P2S_lcl) & !is.na(P2S_median_count) ~ P2S_lcl,
                         !is.na(carcass_90_lcl) & !is.na(carcass_estimate) ~ carcass_90_lcl,
                         TRUE ~ NA),
         ucl = case_when(!is.na(P2S_ucl) & !is.na(P2S_median_count) ~ P2S_ucl,
                         !is.na(carcass_90_ucl) & !is.na(carcass_estimate) ~ carcass_90_ucl,
                         TRUE ~ NA)) |>
  select(-c(P2S_lcl, P2S_ucl, carcass_90_lcl, carcass_90_ucl)) |>
  arrange(stream, year) |>
  glimpse()

table_for_spawn_recruit_long <- table_for_spawn_recruit |>
  mutate(cl_type = case_when(!is.na(lcl | ucl) & !is.na(P2S_median_count) ~ "P2S_modeled_count",
                             !is.na(lcl | ucl) & !is.na(carcass_estimate) ~ "carcass_estimate",
                             TRUE ~ NA)) |>
  rename(P2S_modeled_count = P2S_median_count) |>
  pivot_longer(c(P2S_modeled_count, upstream_estimate, redd_count, holding_count,
                 carcass_estimate),
               names_to = "data_type", values_to = "count") |>
  filter(!is.na(count)) |>
  mutate(lcl = ifelse(cl_type != data_type, NA, lcl),
         ucl = ifelse(cl_type != data_type, NA, ucl)) |>
  relocate(c(data_type, count), .before = lcl) |>
    select(-c(cl_type)) |>
  glimpse()


# upload to google cloud --------------------------------------------------

f <- function(input, output) write_csv(input, file = output)

gcs_upload(table_for_spawn_recruit_long,
           object_function = f,
           type = "csv",
           name = "jpe-model-data/adult-model/adult_data_for_spawn_recruit.csv")

write.csv(table_for_spawn_recruit_long, here::here("data-raw", "adult_model", "adult_model_data",
                                             "adult_data_for_spawn_recruit.csv"),
          row.names = FALSE)
write.csv(table_for_spawn_recruit, here::here("data-raw", "adult_model", "adult_model_data",
                                              "table_for_spawn_recruit_wide.csv"),
          row.names = F)
