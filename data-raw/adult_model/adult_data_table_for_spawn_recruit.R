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

# add bt-spas-x estimates -------------------------------------------------

# TODO code for pulling these are in a different branch
bt_spas_estimates <- read.csv("~/Desktop/bt_spas_x_outputs.csv") |>
  mutate(lcl = round(lcl_total_abundance, 0),
         ucl = round(ucl_total_abundance, 0),
         count = round(mean_total_abundance, 0),
         data_type = "bt_spas_estimates") |>
  filter(site != "ucc", # TODO using lcc
         !site %in% c("herringer riffle")) |> # TODO decide on feather to use
  select(-c(filepath, mean_total_abundance, lcl_total_abundance, ucl_total_abundance)) |>
  glimpse()

spawn_recruit_table_full <- bind_rows(table_for_spawn_recruit_long,
                                      bt_spas_estimates) |>
  arrange(stream, year) |>
  glimpse()

full_spawn_recruit_table_wide <- spawn_recruit_table_full |>
  select(-c(ucl, lcl, site)) |>
  pivot_wider(id_cols = c(year, stream),
              names_from = data_type,
              values_from = count) |>
  glimpse()

modelable_years <- full_spawn_recruit_table_wide |>
  mutate(adult_data = ifelse(!is.na(upstream_estimate) |
                              !is.na(carcass_estimate) |
                              !is.na(P2S_modeled_count) |
                              !is.na(redd_count) |
                              !is.na(holding_count),
                            TRUE,
                            FALSE),
         juv_data = ifelse(!is.na(bt_spas_estimates), TRUE, FALSE),
         can_model = ifelse(adult_data & juv_data, TRUE, FALSE)) |>
  filter(can_model == TRUE)

spawn_recruit_table_full |>
  left_join(modelable_years |>
              select(year, stream, can_model),
            by = c("year", "stream")) |>
  mutate(juv = ifelse(data_type == "bt_spas_estimates", TRUE, FALSE)) |>
  filter(can_model == TRUE) |>
  ggplot(aes(x = year, y = count, color = juv)) +
  geom_line() +
  facet_wrap(~stream)

# scale with redds_below_rsts ---------------------------------------------

# TODO are we using this? not for right now
# draft redds below rsts from Skyler analysis in JPE-datasets
rbr <- read.csv(here::here("data-raw", "adult_model", "redds_by_rst_catchment_summary.csv")) |>
  glimpse()
feather_sites <- read.csv(here::here("data-raw", "adult_model", "feather_annual_site_selection.csv")) |>
  filter(year %in% seq(2015, 2017))

sites_to_scale <- rbr |>
  mutate(stream_year = paste0(stream, "_", date_year),
         scale_factor = 1 + percent_below_rst) |>
  filter(!stream_year %in% c("feather river_2012", "feather river_2013", "feather river_2014",
                             "feather river_2018", "feather river_NA")) |>
  select(-stream_year) |>
  mutate(site_to_use = case_when(stream == "battle creek" ~ "ubc",
                                 stream == "clear creek" ~ "lcc",
                                 stream == "feather river" ~ "herringer riffle"),
         supplementary_site_to_use = case_when(stream == "battle creek" ~ "ubc",
                                               stream == "clear creek" ~ "ucc",
                                               stream == "feather river" ~ "gateway riffle")) |>
  select(stream, year = date_year, percent_below_rst, scale_factor)

scaled_table_for_spawn_recruit <- table_for_spawn_recruit_long |>
  left_join(sites_to_scale, by = c("stream", "year")) |>
  mutate(scaled_count = round(count * scale_factor, 0)) |>
  glimpse()



