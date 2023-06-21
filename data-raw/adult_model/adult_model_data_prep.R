# pull in and prep data for adult model
# see adult-modeling.Rmd for details and rationale

# libraries ---------------------------------------------------------------
library(tidyverse)
library(googleCloudStorageR)
library(waterYearType)

# pull adult data & process ----------------------------------------------------------------
gcs_auth(json_file = Sys.getenv("GCS_AUTH_FILE"))
gcs_global_bucket(bucket = Sys.getenv("GCS_DEFAULT_BUCKET"))

# upstream passage - use these data for passage timing calculations
upstream_passage <- read_csv(gcs_get_object(object_name = "standard-format-data/standard_adult_upstream_passage.csv",
                                            bucket = gcs_get_global_bucket())) |>
  filter(!is.na(date)) |>
  mutate(stream = tolower(stream),
         year = year(date)) |>
  filter(run %in% c("spring", NA, "not recorded")) |>
  group_by(year, passage_direction, stream) |>
  summarise(count = sum(count, na.rm = T)) |>
  ungroup() |>
  pivot_wider(names_from = passage_direction, values_from = count) |>
  # calculate upstream passage for streams where passage direction is recorded
  mutate(down = ifelse(is.na(down), 0, down),
         up = case_when(stream %in% c("deer creek", "mill creek") ~ `NA`,
                        !stream %in% c("deer creek", "mill creek") & is.na(up) ~ 0,
                        TRUE ~ up)) |>
  select(-`NA`) |>
  group_by(year, stream) |>
  summarise(count = round(up - down), 0) |>
  select(year, count, stream) |>
  ungroup() |>
  glimpse()

# pull in passage estimates and use these for upstream_count
upstream_passage_estimates <- read_csv(gcs_get_object(object_name = "standard-format-data/standard_adult_passage_estimate.csv",
                                                      bucket = gcs_get_global_bucket())) |>
  mutate(upstream_count = round(passage_estimate, 0)) |>
  glimpse()

# holding
holding <- read_csv(gcs_get_object(object_name = "standard-format-data/standard_holding.csv",
                                   bucket = gcs_get_global_bucket())) |>
  group_by(year, stream) |>
  summarise(count = sum(count, na.rm = T)) |>
  ungroup() |>
  glimpse()

# redd
redd <- read_csv(gcs_get_object(object_name = "standard-format-data/standard_annual_redd.csv",
                                bucket = gcs_get_global_bucket())) |>
  filter(run %in% c("spring", "not recorded")) |>
  # redds in these reaches are likely fall, so set to 0 for battle & clear
  mutate(max_yearly_redd_count = case_when(reach %in% c("R6", "R6A", "R6B", "R7") &
                                             stream %in% c("battle creek", "clear creek") ~ 0,
                                           TRUE ~ max_yearly_redd_count)) |>
  group_by(year, stream) |>
  summarise(count = sum(max_yearly_redd_count, na.rm = T)) |>
  ungroup() |>
  select(year, stream, count) |>
  glimpse()

# raw carcass
carcass <- read_csv(gcs_get_object(object_name = "standard-format-data/standard_carcass.csv",
                                   bucket = gcs_get_global_bucket())) |>
  filter(run %in% c("spring", NA, "unknown")) |>
  group_by(year(date), stream) |>
  summarise(count = sum(count, na.rm = T)) |>
  ungroup() |>
  select(year = `year(date)`, stream, count) |>
  glimpse()

# estimates from CJS model (carcass survey)
carcass_estimates <- read_csv(gcs_get_object(object_name = "standard-format-data/standard_carcass_cjs_estimate.csv",
                                             bucket = gcs_get_global_bucket())) |>
  rename(carcass_spawner_estimate = spawner_abundance_estimate) |>
  glimpse()

# join all together for raw input -----------------------------------------
adult_model_input_raw <- full_join(upstream_passage_estimates |>
                                     select(year, stream,
                                            upstream_estimate = upstream_count),
                                   redd |>
                                     rename(redd_count = count),
                                   by = c("year", "stream")) |>
  full_join(holding |>
              rename(holding_count = count),
            by = c("year", "stream")) |>
  full_join(carcass_estimates |>
              rename(carcass_estimate = carcass_spawner_estimate) |>
              select(-c(lower, upper, confidence_interval)),
            by = c("year", "stream")) |>
  pivot_longer(c(upstream_estimate, redd_count, holding_count, carcass_estimate),
               values_to = "count",
               names_to = "data_type") |>
  filter(!is.na(count)) |>
  arrange(stream, year) |>
  glimpse()



# temperature -------------------------------------------------------------

# threshold
# https://www.noaa.gov/sites/default/files/legacy/document/2020/Oct/07354626766.pdf
threshold <- 20

# migratory temps - sac, months = 3:5
# holding temps - trib specific; 5-7

standard_temp <- read_csv(gcs_get_object(object_name = "standard-format-data/standard_temperature.csv",
                                         bucket = gcs_get_global_bucket()))

# temperature covariates: migratory temperature (march - may in sacramento river)
migratory_temp <- standard_temp |>
  filter(stream == "sacramento river") |>
  filter(month(date) %in% 3:5) |>
  group_by(year(date)) |>
  mutate(above_threshold = ifelse(mean_daily_temp_c > threshold, TRUE, FALSE)) |>
  summarise(prop_days_exceed_threshold = round(sum(above_threshold, na.rm = T)/length(above_threshold), 2)) |>
  ungroup() |>
  mutate(prop_days_below_threshold = 1 - prop_days_exceed_threshold,
         prop_days_below_threshold = ifelse(prop_days_below_threshold == 0, 0.001, prop_days_below_threshold)) |>
  rename(year = `year(date)`) |>
  select(year, prop_days_exceed_threshold_migratory = prop_days_exceed_threshold) |>
  glimpse()

# temperature covariates: migratory temperature (may - july by tributary)
holding_temp <- standard_temp |>
  filter(month(date) %in% 5:7) |>
  group_by(year(date), stream) |>
  mutate(above_threshold = ifelse(mean_daily_temp_c > threshold, TRUE, FALSE)) |>
  summarise(prop_days_exceed_threshold = round(sum(above_threshold, na.rm = T)/length(above_threshold), 2)) |>
  ungroup() |>
  mutate(prop_days_below_threshold = 1 - prop_days_exceed_threshold,
         prop_days_below_threshold = ifelse(prop_days_below_threshold == 0, 0.001, prop_days_below_threshold)) |>
  rename(year = `year(date)`) |>
  select(prop_days_exceed_threshold_holding = prop_days_exceed_threshold,
         stream, year) |>
  glimpse()

# degree days 20 (called gdd here)
gdd_base_sac <- 20 # https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0204274
gdd_base_trib <- 20

gdd_sac <- standard_temp |>
  filter(month(date) %in% 3:5, stream == "sacramento river") |>
  mutate(gdd_sac = mean_daily_temp_c - gdd_base_sac,
         gdd_sac = ifelse(gdd_sac < 0, 0, gdd_sac)) |>
  group_by(year(date)) |>
  summarise(gdd_sac = sum(gdd_sac, na.rm = T)) |>
  rename(year = `year(date)`) |>
  ungroup()

gdd_trib <- standard_temp |>
  filter(month(date) %in% 5:8 & stream != "sacramento river") |>
  mutate(gdd_trib = mean_daily_temp_c - gdd_base_trib,
         gdd_trib = ifelse(gdd_trib < 0, 0, gdd_trib)) |>
  group_by(year(date), stream) |>
  summarise(gdd_trib = sum(gdd_trib, na.rm = T)) |>
  rename(year = `year(date)`) |>
  ungroup()

gdd <- left_join(gdd_trib, gdd_sac,
                 by = c("year")) |>
  mutate(gdd_sac = ifelse(is.na(gdd_sac), 0, gdd_sac),
         gdd_total = round(gdd_sac + gdd_trib, 2)) |>
  glimpse()

# https://www.rdocumentation.org/packages/pollen/versions/0.82.0/topics/gdd
# https://www.researchgate.net/publication/279930331_Fish_growth_and_degree-days_I_Selecting_a_base_temperature_for_a_within-population_study

# flow --------------------------------------------------------------------

standard_flow <- read_csv(gcs_get_object(object_name = "standard-format-data/standard_flow.csv",
                                         bucket = gcs_get_global_bucket())) |>
  filter(month(date) %in% 3:8) |>
  mutate(year = year(date)) |>
  group_by(stream, year) |>
  summarise(mean_flow = mean(flow_cfs, na.rm = T),
            max_flow = max(flow_cfs, na.rm = T)) |>
  glimpse()

# prespawn survival -------------------------------------------------------

carcass_streams <- c("butte creek", "feather river")

prespawn_survival <- left_join(upstream_passage_estimates |>
                                 select(-passage_estimate),
                               redd |>
                                 rename(redd_count = count),
                               by = c("year", "stream")) |>
  left_join(holding |>
              rename(holding_count = count),
            by = c("year", "stream")) |>
  left_join(carcass_estimates |>
              rename(carcass_estimate = carcass_spawner_estimate) |>
              select(-c(upper, lower, confidence_interval)),
            by = c("year", "stream")) |>
  mutate(female_upstream = upstream_count * 0.5,
         prespawn_survival = case_when(stream == "deer creek" ~ holding_count / upstream_count,
                                       stream %in% carcass_streams ~ carcass_estimate / upstream_count,
                                       TRUE ~ redd_count / female_upstream)) |>
  filter(prespawn_survival != Inf) |>
#         stream != "butte creek") |>
  glimpse()


# passage timing ----------------------------------------------------------
upstream_passage_timing <- read_csv(gcs_get_object(object_name = "standard-format-data/standard_adult_upstream_passage.csv",
                                            bucket = gcs_get_global_bucket())) |>
  filter(!is.na(date)) |>
  mutate(stream = tolower(stream),
         year = year(date),
         week = week(date)) |>
  filter(run %in% c("spring","not recorded")) |>
  group_by(year, stream) |>
  summarise(count = sum(count, na.rm = T),
            median_passage_timing = median(week, na.rm = T),
            mean_passage_timing = mean(week, na.rm = T),
            min_passage_timing = min(week, na.rm = T)) |>
  ungroup() |> # TODO look at up-down
  select(-c(count)) |> glimpse()

# water year --------------------------------------------------------------

water_year_data <- waterYearType::water_year_indices |>
  mutate(water_year_type = case_when(Yr_type %in% c("Wet", "Above Normal") ~ "wet",
                               Yr_type %in% c("Dry", "Below Normal", "Critical") ~ "dry",
                               TRUE ~ Yr_type)) |>
  filter(location == "Sacramento Valley") |>
  dplyr::select(WY, water_year_type) |>
  glimpse()


# standardized covariates -------------------------------------------------
adult_model_covariates_standard <- full_join(standard_flow,
                                             gdd,
                                             by = c("year", "stream")) |>
  full_join(upstream_passage_timing,
            by = c("year", "stream")) |>
  full_join(water_year_data,
            by = c("year" = "WY")) |>
  filter(!is.na(stream),
         stream != "sacramento river") |>
  select(-c(mean_flow, mean_passage_timing, min_passage_timing,
            gdd_trib, gdd_sac)) |>
  mutate(wy_type = ifelse(water_year_type == "dry", 0, 1),
         max_flow_std = as.vector(scale(max_flow)),
         gdd_std = as.vector(scale(gdd_total)),
         median_passage_timing_std = as.vector(scale(median_passage_timing))) |>
  select(year, stream, wy_type, max_flow_std, gdd_std, median_passage_timing_std) |>
  arrange(stream, year) |>
  glimpse()

# combine -----------------------------------------------------------------

survival_model_data_raw <- left_join(prespawn_survival |>
                                       drop_na(prespawn_survival),
                                migratory_temp,
                                by = c("year")) |>
  left_join(holding_temp,
            by = c("year", "stream")) |>
  mutate(total_prop_days_exceed_threshold = ifelse(is.na(prop_days_exceed_threshold_migratory), prop_days_exceed_threshold_holding,
                                                   (prop_days_exceed_threshold_migratory + prop_days_exceed_threshold_holding) / 2)) |>
  left_join(standard_flow,
            by = c("year", "stream")) |>
  left_join(upstream_passage_timing,
            by = c("year", "stream")) |>
  left_join(water_year_data,
            by = c("year" = "WY")) |>
  left_join(gdd,
            by = c("year", "stream")) |>
  glimpse()


# yuba data ---------------------------------------------------------------
yuba_data <- upstream_passage_estimates |>
  filter(stream == "yuba river") |>
  select(year, passage_estimate) |>
  glimpse()


# butte data --------------------------------------------------------------
butte_data <- carcass_estimates |>
  filter(stream == "butte creek") |>
  select(year, spawner_estimate = carcass_spawner_estimate) |>
  glimpse()


# feather data ------------------------------------------------------------

feather_data <- carcass_estimates |>
  filter(stream == "feather river") |>
  select(year, spawner_estimate = carcass_spawner_estimate) |>
  glimpse()

# write data objects to bucket -------------------------------------------------------
f <- function(input, output) write_csv(input, file = output)

gcs_upload(adult_model_input_raw,
           object_function = f,
           type = "csv",
           name = "jpe-model-data/adult-model/adult_data_input_raw.csv")
gcs_upload(adult_model_covariates_standard,
           object_function = f,
           type = "csv",
           name = "jpe-model-data/adult-model/adult_model_covariates_standard.csv")
gcs_upload(survival_model_data_raw,
           object_function = f,
           type = "csv",
           name = "jpe-model-data/adult-model/survival_model_data_raw.csv")
gcs_upload(yuba_data,
           object_function = f,
           type = "csv",
           name = "jpe-model-data/adult-model/yuba_data.csv")
gcs_upload(feather_data,
           object_function = f,
           type = "csv",
           name = "jpe-model-data/adult-model/feather_data.csv")
gcs_upload(butte_data,
           object_function = f,
           type = "csv",
           name = "jpe-model-data/adult-model/butte_data.csv")

