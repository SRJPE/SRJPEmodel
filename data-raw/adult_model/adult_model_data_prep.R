# pull in and prep data for adult model
# see adult-modeling.Rmd for details and rationale

# libraries ---------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(googleCloudStorageR)
library(rstan)
library(rstanarm)
library(bayesplot)
library(GGally) # pairs plot
library(waterYearType)
library(car) # vif
library(glmulti)
library(tidybayes)

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
                                   bucket = gcs_get_global_bucket()))|>
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
  # TODO use carcass counts or carcass estimates for prespawn survival exploration?
  left_join(carcass |>
              rename(carcass_count = count),
            by = c("year", "stream")) |>
  mutate(female_upstream = upstream_count * 0.5,
         prespawn_survival = case_when(stream == "deer creek" ~ holding_count / upstream_count,
                                       stream %in% carcass_streams ~ carcass_count / upstream_count,
                                       TRUE ~ redd_count / female_upstream)) |>
  filter(prespawn_survival != Inf,
         stream != "butte creek") |>
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


# plots to inspect different covars ---------------------------------------
survival_model_data_raw |>
  filter(!stream %in% c("mill creek", "butte creek")) |> # mill & butte only have 3 data points; is skewing plot
  ggplot(aes(x = total_prop_days_exceed_threshold, y = prespawn_survival, fill = stream)) +
  geom_point(aes(color = stream)) + geom_smooth(method = "lm") +
  theme_minimal() + ggtitle("Prespawn survival and temperature by stream") +
  xlab("Proportion of days exceeding threshold temperature") +
  ylab("Prespawn survival")

survival_model_data_raw |>
  filter(!stream %in% c("mill creek", "butte creek")) |> # mill & butte only have 3 data points; is skewing plot
  ggplot(aes(x = gdd_total, y = prespawn_survival, fill = stream)) +
  geom_point(aes(color = stream)) + geom_smooth(method = "lm") +
  theme_minimal() + ggtitle("Prespawn survival and GDD by stream") +
  xlab("Growing degree days over 20 (GDD)") +
  ylab("Prespawn survival")

survival_model_data_raw |>
  filter(!stream %in% c("mill creek", "yuba river")) |>
  ggplot(aes(x = mean_flow, y = prespawn_survival, fill = stream)) +
  geom_point(aes(color = stream)) + geom_smooth(method = "lm")  +
  theme_minimal() + ggtitle("Prespawn survival and mean flow by stream") +
  xlab("Mean flow (cfs)") +
  ylab("Prespawn survival")

survival_model_data_raw |>
  filter(stream != "mill creek") |>
  ggplot(aes(x = max_flow, y = prespawn_survival, fill = stream)) +
  geom_point(aes(color = stream)) + geom_smooth(method = "lm")  +
  theme_minimal() + ggtitle("Prespawn survival and max flow by stream") +
  xlab("Max flow (cfs)") +
  ylab("Prespawn survival")

survival_model_data_raw |>
  filter(stream == "mill creek") |>
  ggplot(aes(x = max_flow, y = prespawn_survival)) +
  geom_point(aes(color = stream))

survival_model_data_raw |>
  #filter(stream != "mill creek") |>
  ggplot(aes(x = min_passage_timing, y = prespawn_survival, fill = stream)) +
  geom_point(aes(color = stream)) + geom_smooth(method = "lm")   +
  theme_minimal() + ggtitle("Prespawn survival and minimum passage time by stream") +
  xlab("Minimum passage time (weeks)") +
  ylab("Prespawn survival")

survival_model_data_raw |>
  ggplot(aes(x = total_prop_days_exceed_threshold, y = prespawn_survival, fill = stream)) +
  geom_point(aes(color = stream)) + geom_smooth(method = "lm") +
  facet_wrap(~water_year_type, scales = "free") +
  theme_minimal() + ggtitle("Prespawn survival and temperature by stream and water year type") +
  xlab("Proportion of days exceeding temperature threshold") +
  ylab("Prespawn survival")


# remove variables with no relationship -----------------------------------

# decide on one passage timing variable - remove minimum passage timing; does not
# represent bulk of the fish moving through the system.
# remove median passage timing; plotting difference between mean and median does not show
# difference greater than +- 2 and previous VIF analyses and AIC model comparisons showed
# median of greater importance

# decide on one gdd variable - total contains the most complete information about the
# lifecycle

# remove covariates that either show no relationship (mean flow) or could be
# better represented by other metrics (all prop_days temperature variables)
# also, remove raw counts
variables_to_remove <- c("mean_flow", "prop_days_exceed_threshold_holding",
                         "prop_days_exceed_threshold_migratory",
                         "total_prop_days_exceed_threshold",
                         "min_passage_timing", "mean_passage_timing",
                         "gdd_sac", "gdd_trib",
                         "upstream_count", "redd_count", "female_upstream",
                         "holding_count", "carcass_count")

survival_model_data <- survival_model_data_raw |>
  select(-all_of(variables_to_remove)) |>
  glimpse()

# check for collinearity for each stream --------------------------------
# function to print Pearsons correlations
print_cors <- function(data, cor_threshold) {

  if("water_year_type" %in% names(data)) {
    new_dat <- data |>
      select(-c(prespawn_survival, water_year_type)) |> # can't calculate cor for a categorical variable
      drop_na() |>
      cor() |>
      as.matrix()
  } else {
    new_dat <- data |>
      select(-prespawn_survival) |> # can't calculate cor for a categorical variable
      drop_na() |>
      cor() |>
      as.matrix()
  }

  new_dat[lower.tri(new_dat)] <- NA # get rid of duplicates (lower tri of matrix)

  new_dat |>
    as.data.frame() |>
    rownames_to_column(var = "variable") |>
    pivot_longer(cols = -variable, names_to = "variable_2",
                 values_to = "correlation") |>
    filter(!is.na(correlation),
           abs(correlation) >= cor_threshold,
           variable != variable_2)
}

# steps for each stream:
# 1. look at ggpairs
# 2. look at Pearson's correlations above threshold 0.65
# 3. use VIF to identify correlated variables (threshold is 5)
# https://online.stat.psu.edu/stat462/node/180/ (VIF)
# 4. use steps 1:3 to select ONE passage timing variable and ONE temperature variable -
# you can't use the glmulti() function with all these variables - it won't converge
# theory - GDD is the standard. stronger than total_prop_days
# 5. use glmuli to look for the best model (by AIC) - including interactions

battle_data_full <- survival_model_data |>
  filter(stream == "battle creek")
battle_data <- battle_data_full |>
  select(-c(year, stream))
ggpairs(battle_data |> drop_na())
print_cors(battle_data, 0.65) # gdd total correlated with each other; mean passage timing correlated with each other
vif(lm(prespawn_survival ~ ., data = battle_data))

# remove variables with highest VIF values
battle_variables_remove <- c()


# now look for interactions using glmulti
best_battle_model <- glmulti(y = "prespawn_survival",
                             xr = battle_data |> select(-c("prespawn_survival",
                                                           all_of(battle_variables_remove))) |>
                               names(),
                             intercept = TRUE,
                             method = "h",
                             maxsize = 1,
                             level = 1,
                             data = battle_data,
                             fitfunction = "lm")
summary(best_battle_model)$bestmodel

# clear
# clear had lots of points where prespawn_survival > 1 (filtered earlier). This gets rid
# of any years with different water year types, so we exclude this as a consideration although
# as data availability increases we may want to include water year type
clear_data_full <- survival_model_data |>
  filter(stream == "clear creek")
clear_data <- clear_data_full |>
  select(-c(stream,  year))
ggpairs(clear_data |> drop_na()) # lots of NAs for median_passage
print_cors(clear_data, 0.65)
vif(lm(prespawn_survival ~ ., data = clear_data |> select(-c(median_passage_timing))))
vif(lm(prespawn_survival ~ ., data = clear_data))


# remove variables with highest VIF values
clear_variables_remove <- c("median_passage_timing") # NAs for many years
clear_variables_remove <- c()
# now look for interactions using glmulti
best_clear_model <- glmulti(y = "prespawn_survival",
                            xr = clear_data |> select(-c("prespawn_survival",
                                               all_of(clear_variables_remove))) |>
                              names(),
                            intercept = TRUE,
                            method = "h",
                            maxsize = 1,
                            level = 1,
                            data = clear_data,
                            fitfunction = "lm")
summary(best_clear_model)$bestmodel

# mill
# mill had lots of points where prespawn_survival > 1 (filtered earlier). This gets rid
# of any years with different water year types, so we exclude this as a consideration although
# as data availability increases we may want to include water year type.
# currently only 3 data points for mill...this is not a lot of statistical power
mill_data_full <- survival_model_data |>
  filter(stream == "mill creek")
mill_data <- mill_data_full |>
  select(-c(year, stream))

ggpairs(mill_data)
print_cors(mill_data, 0.65)
vif(lm(prespawn_survival ~ ., data = mill_data |> select(-c(median_passage_timing))))
vif(lm(prespawn_survival ~ ., data = mill_data |> select(-c(gdd_total)))) # median passage timing increases VIF


# remove variables with highest VIF values
mill_variables_remove <- c("median_passage_timing")

# now look for interactions using glmulti
best_mill_model <- glmulti(y = "prespawn_survival",
                           xr = mill_data |> select(-c("prespawn_survival",
                                                       all_of(mill_variables_remove))) |>
                             names(),
                           intercept = TRUE,
                           method = "h",
                           maxsize = 1,
                           level = 1,
                           data = mill_data,
                           fitfunction = "lm")
summary(best_mill_model)$bestmodel

# deer
deer_data_full <- survival_model_data |>
  filter(stream == "deer creek")
deer_data <- deer_data_full |>
  select(-c(year, stream))
ggpairs(deer_data)
print_cors(deer_data, 0.65) # gdd trib is correlated w/ max flow & median passage
vif(lm(prespawn_survival ~ ., data = deer_data |> select(-c(median_passage_timing))))
vif(lm(prespawn_survival ~ ., data = deer_data |> select(-c(max_flow))))

# remove variables with highest VIF values
deer_variables_remove <- c("median_passage_timing")

# now look for interactions using glmulti
best_deer_model <- glmulti(y = "prespawn_survival",
                           xr = deer_data |> select(-c("prespawn_survival",
                                                       all_of(deer_variables_remove))) |>
                             names(),
                           intercept = TRUE,
                           method = "h",
                           level = 1,
                           maxsize = 1,
                           data = deer_data,
                           fitfunction = "lm")
summary(best_deer_model)$bestmodel


# plot best models and get estimates of coefficients
# this is important because these are used to create the environmental index
# used in the dauphin bayesian models
summary(best_battle_model)$bestmodel
best_battle_lm <- lm(prespawn_survival ~ 1 + water_year_type,
                     data = battle_data)

summary(best_clear_model)$bestmodel
best_clear_lm <- lm(prespawn_survival ~ 1 + max_flow,
                     data = clear_data)

summary(best_mill_model)$bestmodel
best_mill_lm <- lm(prespawn_survival ~ 1 + gdd_total,
                     data = mill_data)

summary(best_deer_model)$bestmodel
best_deer_lm <- lm(prespawn_survival ~ 1 + water_year_type,
                   data = deer_data)


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

gcs_upload(survival_model_data_raw,
           object_function = f,
           type = "csv",
           name = "jpe-model-data/adult-model/survival_model_data_raw.csv")
gcs_upload(survival_model_data,
           object_function = f,
           type = "csv",
           name = "jpe-model-data/adult-model/survival_model_data.csv")
gcs_upload(battle_data_full,
           object_function = f,
           type = "csv",
           name = "jpe-model-data/adult-model/battle_data.csv")
gcs_upload(clear_data_full,
           object_function = f,
           type = "csv",
           name = "jpe-model-data/adult-model/clear_data.csv")
gcs_upload(deer_data_full,
           object_function = f,
           type = "csv",
           name = "jpe-model-data/adult-model/deer_data.csv")
gcs_upload(mill_data_full,
           object_function = f,
           type = "csv",
           name = "jpe-model-data/adult-model/mill_data.csv")
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

# save data objects -------------------------------------------------------

adult_data_objects <- list("survival_model_data_raw" = survival_model_data_raw,
                           "survival_model_data" = survival_model_data,
                           "battle_data" = battle_data_full,
                           "clear_data" = clear_data_full,
                           "mill_data" = mill_data_full,
                           "deer_data" = deer_data_full,
                           "yuba_data" = yuba_data,
                           "butte_data" = butte_data,
                           "feather_data" = feather_data)

save(adult_data_objects, file = here::here("data-raw", "adult_model", "adult_data_objects.Rdata"))

best_adult_models <- list("best_battle_lm" = best_battle_lm,
                          "best_clear_lm" = best_clear_lm,
                          "best_mill_lm" = best_mill_lm,
                          "best_deer_lm" = best_deer_lm)
save(best_adult_models, file = here::here("data-raw", "adult_model", "best_adult_models.Rdata"))

