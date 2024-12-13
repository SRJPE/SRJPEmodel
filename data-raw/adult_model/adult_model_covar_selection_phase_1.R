# This script compares linear regressions of prespawn survival on
# different environmental covariates to identify the covariate with
# the most predictive power for each stream.

# see adult-modeling.Rmd for details

# libraries ---------------------------------------------------------------
library(tidyverse)
library(googleCloudStorageR)
library(GGally) # pairs plot
library(waterYearType)
library(car) # vif
library(glmulti)

# pull in prepped data ----------------------------------------------------------------
# source(here::here("here", "data-raw", "adult_model_data_prep.R"))
gcs_get_object(object_name = "jpe-model-data/adult-model/survival_model_data_raw.csv",
               bucket = gcs_get_global_bucket(),
               saveToDisk = here::here("data-raw", "adult_model", "adult_model_data",
                                       "survival_model_data_raw.csv"),
               overwrite = TRUE)

survival_model_data_raw <- read.csv(here::here("data-raw", "adult_model", "adult_model_data",
                                              "survival_model_data_raw.csv"))



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
  geom_point(aes(color = stream)) +
  theme_minimal()

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
  facet_wrap(~wy_type, scales = "free") +
  theme_minimal() + ggtitle("Prespawn survival and temperature by stream and water year type") +
  xlab("Proportion of days exceeding temperature threshold") +
  ylab("Prespawn survival")


# remove variables with no relationship -----------------------------------

# decide on one passage timing variable, one temperature index variable,
# one flow variable, and one water year type variable

# remove covariates and raw counts
variables_to_remove <- c("mean_flow", "prop_days_exceed_threshold_holding",
                         "prop_days_exceed_threshold_migratory",
                         "total_prop_days_exceed_threshold",
                         "min_passage_timing", "mean_passage_timing",
                         "gdd_sac", "gdd_trib",
                         "upstream_count", "redd_count", "female_upstream",
                         "holding_count", "carcass_estimate")

survival_model_data <- survival_model_data_raw |>
  select(-all_of(variables_to_remove)) |>
  mutate(wy_type = ifelse(water_year_type == "dry", 0, 1),
         max_flow_std = as.vector(scale(max_flow)),
         gdd_std = as.vector(scale(gdd_total)),
         median_passage_timing_std = as.vector(scale(median_passage_timing))) |>
  select(year, stream, prespawn_survival,
         wy_type, max_flow_std, gdd_std, median_passage_timing_std) |>
  glimpse()

# check for collinearity for each stream --------------------------------
# function to print Pearsons correlations
print_cors <- function(data, cor_threshold) {

  if("wy_type" %in% names(data)) {
    new_dat <- data |>
      select(-c(prespawn_survival, wy_type)) |> # can't calculate cor for a categorical variable
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
clear_data_full <- survival_model_data |>
  filter(stream == "clear creek")
clear_data <- clear_data_full |>
  select(-c(stream,  year))
ggpairs(clear_data |> drop_na()) # lots of NAs for median_passage
print_cors(clear_data, 0.65)
vif(lm(prespawn_survival ~ ., data = clear_data |> select(-c(median_passage_timing_std))))
vif(lm(prespawn_survival ~ ., data = clear_data))


# remove variables with highest VIF values
clear_variables_remove <- c("median_passage_timing_std") # NAs for many years

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
mill_data_full <- survival_model_data |>
  filter(stream == "mill creek")
mill_data <- mill_data_full |>
  select(-c(year, stream))

ggpairs(mill_data)
print_cors(mill_data, 0.65)
vif(lm(prespawn_survival ~ ., data = mill_data |> select(-c(median_passage_timing_std))))
vif(lm(prespawn_survival ~ ., data = mill_data |> select(-c(gdd_std)))) # median passage timing increases VIF


# remove variables with highest VIF values
mill_variables_remove <- c("median_passage_timing_std")

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
vif(lm(prespawn_survival ~ ., data = deer_data |> select(-c(median_passage_timing_std))))
vif(lm(prespawn_survival ~ ., data = deer_data |> select(-c(max_flow_std))))

# remove variables with highest VIF values
deer_variables_remove <- c("median_passage_timing_std")

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
best_battle_lm <- lm(prespawn_survival ~ 1 + wy_type,
                     data = battle_data)

summary(best_clear_model)$bestmodel
best_clear_lm <- lm(prespawn_survival ~ 1 + max_flow_std,
                    data = clear_data)

summary(best_mill_model)$bestmodel
best_mill_lm <- lm(prespawn_survival ~ 1 + gdd_std,
                   data = mill_data)

summary(best_deer_model)$bestmodel
best_deer_lm <- lm(prespawn_survival ~ 1 + wy_type,
                   data = deer_data)


# write data objects to bucket -------------------------------------------------------
f <- function(input, output) write_csv(input, file = output)

gcs_upload(survival_model_data,
           object_function = f,
           type = "csv",
           name = "jpe-model-data/adult-model/survival_model_data.csv")

# save data objects -------------------------------------------------------
best_adult_models <- list("best_battle_lm" = best_battle_lm,
                          "best_clear_lm" = best_clear_lm,
                          "best_mill_lm" = best_mill_lm,
                          "best_deer_lm" = best_deer_lm)
save(best_adult_models, file = here::here("data-raw", "adult_model", "best_adult_models.Rdata"))

