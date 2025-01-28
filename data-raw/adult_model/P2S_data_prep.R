# pull in and prep data for adult model
# see adult-modeling.Rmd for details and rationale

# libraries ---------------------------------------------------------------
library(tidyverse)
library(googleCloudStorageR)
library(waterYearType)

# pull adult data & process ----------------------------------------------------------------

# pull in passage estimates and use these for upstream_count
upstream_passage_estimates <- SRJPEdata::upstream_passage_estimates |>
  mutate(upstream_estimate = round(passage_estimate, 0)) |>
  filter(upstream_estimate > 0) |>
  select(year, stream, upstream_estimate) |>
  glimpse()

# holding
holding <- SRJPEdata::holding |>
  group_by(year = year(date), stream) |>
  summarise(holding_count = sum(count, na.rm = T)) |>
  ungroup() |>
  glimpse()

# redd
redd <- SRJPEdata::redd |>
  filter(run %in% c("spring", "not recorded")) |>
  mutate(count = 1) |>
  # redds in these reaches are likely fall, so set to 0 for battle & clear
  mutate(count = case_when(reach %in% c("R6", "R6A", "R6B", "R7") &
                                             stream %in% c("battle creek", "clear creek") ~ 0,
                                           TRUE ~ count)) |>
  group_by(year = year(date), stream) |>
  summarise(redd_count = sum(count, na.rm = T)) |>
  ungroup() |>
  glimpse()

# estimates from CJS model (carcass survey)
carcass_estimates <- SRJPEdata::carcass_estimates |>
  select(year, stream, carcass_estimate) |>
  glimpse()

# join all together for raw input -----------------------------------------
adult_model_input_raw <- full_join(upstream_passage_estimates,
                                   redd,
                                   by = c("year", "stream")) |>
  full_join(holding,
            by = c("year", "stream")) |>
  full_join(carcass_estimates,
            by = c("year", "stream")) |>
  pivot_longer(upstream_estimate:carcass_estimate,
               values_to = "count",
               names_to = "data_type") |>
  filter(!is.na(count)) |>
  arrange(stream, year) |>
  glimpse()


# temperature -------------------------------------------------------------

# threshold
# https://www.noaa.gov/sites/default/files/legacy/document/2020/Oct/07354626766.pdf

# TODO add new thresholds
threshold <- 20

# migratory temps - sac, months = 3:5
# holding temps - trib specific; 5-7

standard_temp <- SRJPEdata::environmental_data |>
  filter(parameter == "temperature",
         statistic == "mean")

# temperature covariates: migratory temperature (march - may in sacramento river)
migratory_temp <- standard_temp |>
  filter(stream == "sacramento river",
         month %in% 3:5) |>
  group_by(year) |>
  mutate(above_threshold = ifelse(value > threshold, TRUE, FALSE)) |>
  summarise(prop_days_exceed_threshold = round(sum(above_threshold, na.rm = T)/length(above_threshold), 2)) |>
  ungroup() |>
  mutate(prop_days_below_threshold = 1 - prop_days_exceed_threshold,
         prop_days_below_threshold = ifelse(prop_days_below_threshold == 0, 0.001, prop_days_below_threshold)) |>
  select(year, prop_days_exceed_threshold_migratory = prop_days_exceed_threshold) |>
  glimpse()

# temperature covariates: migratory temperature (may - july by tributary)
holding_temp <- standard_temp |>
  filter(month %in% 5:7) |>
  group_by(year, stream) |>
  mutate(above_threshold = ifelse(value > threshold, TRUE, FALSE)) |>
  summarise(prop_days_exceed_threshold = round(sum(above_threshold, na.rm = T)/length(above_threshold), 2)) |>
  ungroup() |>
  mutate(prop_days_below_threshold = 1 - prop_days_exceed_threshold,
         prop_days_below_threshold = ifelse(prop_days_below_threshold == 0, 0.001, prop_days_below_threshold)) |>
  select(prop_days_exceed_threshold_holding = prop_days_exceed_threshold,
         stream, year) |>
  glimpse()

# degree days 20 (called gdd here)
gdd_base_sac <- 20 # https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0204274
gdd_base_trib <- 20

gdd_sac <- standard_temp |>
  filter(month %in% 3:5,
         stream == "sacramento river") |>
  mutate(gdd_sac = value - gdd_base_sac,
         gdd_sac = ifelse(gdd_sac < 0, 0, gdd_sac)) |>
  group_by(year) |>
  summarise(gdd_sac = sum(gdd_sac, na.rm = T)) |>
  ungroup()

gdd_trib <- standard_temp |>
  filter(month %in% 5:8 & stream != "sacramento river") |>
  mutate(gdd_trib = value - gdd_base_trib,
         gdd_trib = ifelse(gdd_trib < 0, 0, gdd_trib)) |>
  group_by(year, stream) |>
  summarise(gdd_trib = sum(gdd_trib, na.rm = T)) |>
  ungroup()

gdd <- left_join(gdd_trib, gdd_sac,
                 by = c("year")) |>
  mutate(gdd_sac = ifelse(is.na(gdd_sac), 0, gdd_sac),
         gdd_total = round(gdd_sac + gdd_trib, 2)) |>
  glimpse()

# https://www.rdocumentation.org/packages/pollen/versions/0.82.0/topics/gdd
# https://www.researchgate.net/publication/279930331_Fish_growth_and_degree-days_I_Selecting_a_base_temperature_for_a_within-population_study

# flow --------------------------------------------------------------------

standard_flow <- SRJPEdata::environmental_data |>
  filter(parameter == "flow",
         statistic == "mean") |>
  filter(month %in% 3:8) |>
  group_by(stream, year) |>
  summarise(mean_flow = mean(value, na.rm = T),
            max_flow = max(value, na.rm = T)) |>
  glimpse()

# prespawn survival -------------------------------------------------------

carcass_streams <- c("butte creek", "feather river")

# left join because we need passage estimates
prespawn_survival <- left_join(upstream_passage_estimates,
                               redd,
                               by = c("year", "stream")) |>
  left_join(holding,
            by = c("year", "stream")) |>
  left_join(carcass_estimates,
            by = c("year", "stream")) |>
  mutate(prespawn_survival = case_when(stream == "deer creek" ~ holding_count / upstream_estimate,
                                       stream %in% carcass_streams ~ carcass_estimate / upstream_estimate,
                                       # assume that for every redd counted, it represents two fish (sex ratio of 0.5)
                                       TRUE ~ (redd_count * 2) / upstream_estimate)) |>
  glimpse()


# passage timing ----------------------------------------------------------
upstream_passage_timing <- SRJPEdata::upstream_passage |>
  mutate(year = year(date),
         week = week(date)) |>
  filter(run %in% c("spring", "unknown", NA),
         direction %in% c("up", "not recorded")) |>
  group_by(year, stream) |>
  summarise(count = sum(count, na.rm = T),
            median_passage_timing = median(week, na.rm = T),
            mean_passage_timing = mean(week, na.rm = T),
            min_passage_timing = min(week, na.rm = T)) |>
  ungroup() |>
  select(-c(count)) |>
  glimpse()

# water year --------------------------------------------------------------

water_year_data <- waterYearType::water_year_indices |>
  mutate(water_year_type = case_when(Yr_type %in% c("Wet", "Above Normal") ~ "wet",
                                     Yr_type %in% c("Dry", "Below Normal", "Critical") ~ "dry",
                                     TRUE ~ Yr_type)) |>
  filter(location == "Sacramento Valley") |>
  select(WY, water_year_type) |>
  glimpse()

# total passage as index --------------------------------------------------
upstream_passage_index <- upstream_passage_estimates |>
  mutate(passage_index = upstream_estimate) |>
  select(year, stream, passage_index)

# standardized covariates -------------------------------------------------
scale_covar <- function(x) {
  as.vector(scale(x))
}

adult_model_covariates_standard <- full_join(standard_flow,
                                             gdd,
                                             by = c("year", "stream")) |>
  full_join(upstream_passage_timing,
            by = c("year", "stream")) |>
  full_join(water_year_data,
            by = c("year" = "WY")) |>
  full_join(upstream_passage_index,
            by = c("year", "stream")) |>
  filter(!is.na(stream),
         stream != "sacramento river") |>
  rename(wy_type = water_year_type) |>
  mutate(wy_type = ifelse(wy_type == "dry", 0, 1),
         across(mean_flow:passage_index, scale_covar)) |>
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

