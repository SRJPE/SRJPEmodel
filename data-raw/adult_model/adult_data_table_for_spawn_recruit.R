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
               names_to = "data_type", values_to = "count") |>
  glimpse()

# bind together the two groups of data -------------------------------------

# pivot wider so that we have a column for each data type
table_for_spawn_recruit <- bind_rows(P2S_model_fits_with_year,
                                     adult_data_input_raw) |>
  mutate(count = round(count, 0)) |>
  pivot_wider(id_cols = c(year, stream),
              names_from = data_type,
              values_from = count) |>
  arrange(stream, year) |>
  glimpse()

# upload to google cloud --------------------------------------------------

f <- function(input, output) write_csv(input, file = output)

gcs_upload(table_for_spawn_recruit,
           object_function = f,
           type = "csv",
           name = "jpe-model-data/adult-model/adult_data_for_spawn_recruit.csv")

write.csv(table_for_spawn_recruit, here::here("data-raw", "adult_model", "adult_model_data",
                                             "adult_data_for_spawn_recruit.csv"),
          row.names = FALSE)
