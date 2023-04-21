# pull in and prep data for adult model

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

# upstream passage total counts
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

# TODO holding time by stream and year
# TODO look at hatchery/natural fish
# TODO run identification


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

# growing degree days
gdd_base_sac <- 0 # TODO confirm this
gdd_base_trib <- 0

gdd_sac <- standard_temp |>
  filter(month(date) %in% 3:5, stream == "sacramento river") |>
  mutate(gdd_sac = round(mean_daily_temp_c - gdd_base_sac, 3)) |>
  group_by(year(date)) |>
  summarise(gdd_sac = round(sum(gdd_sac, na.rm = T), 2)) |>
  rename(year = `year(date)`) |>
  ungroup()

gdd_trib <- standard_temp |>
  filter(month(date) %in% 5:8 & stream != "sacramento river") |>
  mutate(gdd_trib = round(mean_daily_temp_c - gdd_base_trib), 2) |>
  group_by(year(date), stream) |>
  summarise(gdd_trib = round(sum(gdd_trib, na.rm = T), 2)) |>
  rename(year = `year(date)`) |>
  ungroup()

gdd <- left_join(gdd_trib, gdd_sac,
                 by = c("year")) |>
  mutate(gdd_sac = ifelse(is.na(gdd_sac), 0, gdd_sac)) |>
  group_by(year, stream) |>
  summarise(gdd = sum(gdd_trib, gdd_sac)) |>
  ungroup() |>
  mutate(gdd = ifelse(gdd < 0, 0, gdd)) |>
  glimpse()

# https://www.rdocumentation.org/packages/pollen/versions/0.82.0/topics/gdd
# https://www.researchgate.net/publication/279930331_Fish_growth_and_degree-days_I_Selecting_a_base_temperature_for_a_within-population_study

# flow --------------------------------------------------------------------

standard_flow <- read_csv(gcs_get_object(object_name = "standard-format-data/standard_flow.csv",
                                         bucket = gcs_get_global_bucket())) |>
  mutate(year = year(date)) |>
  group_by(stream, year) |>
  summarise(mean_flow = mean(flow_cfs, na.rm = T)) |>
  glimpse()

# prespawn survival -------------------------------------------------------

streams <- c("battle creek", "clear creek", "deer creek", "mill creek")

prespawn_survival <- inner_join(upstream_passage |>
                                  rename(upstream_count = count),
                                redd |>
                                  rename(redd_count = count),
                                by = c("year", "stream")) |>
  mutate(female_upstream = upstream_count * 0.5,
         prespawn_survival = redd_count / female_upstream,
         prespawn_survival = ifelse(prespawn_survival > 1, 1, prespawn_survival)) |>
  filter(stream %in% streams) |>
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
  ungroup() |># TODO worry about up-down
  select(-c(count)) |> glimpse()



# water year --------------------------------------------------------------

water_year_data <- waterYearType::water_year_indices |>
  mutate(water_year_type = case_when(Yr_type %in% c("Wet", "Above Normal") ~ "wet",
                               Yr_type %in% c("Dry", "Below Normal", "Critical") ~ "dry",
                               TRUE ~ Yr_type)) |>
  filter(location == "Sacramento Valley") |>
  dplyr::select(WY, water_year_type) |>
  glimpse()

later_years <- tibble(WY = 2018:2021,
                      water_year_type = c("dry", "wet", "dry", "dry"))

water_year_data <- rbind(water_year_data, later_years)

# combine -----------------------------------------------------------------

survival_model_data_raw <- left_join(prespawn_survival,
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
  drop_na() |>
  glimpse()

survival_model_data <- survival_model_data_raw |>
  dplyr::select(-c(upstream_count, redd_count, female_upstream)) |>
  glimpse()


survival_model_data |>
  ggplot(aes(x = total_prop_days_exceed_threshold, y = prespawn_survival, fill = stream)) +
  geom_point(aes(color = stream)) + geom_smooth(method = "lm") +
  theme_minimal() + ggtitle("Prespawn survival and temperature by stream") +
  xlab("Proportion of days exceeding threshold temperature") +
  ylab("Prespawn survival")

# TODO anomalous point for 2013 clear creek
survival_model_data |>
  ggplot(aes(x = gdd, y = prespawn_survival, fill = stream)) +
  geom_point(aes(color = stream)) + geom_smooth(method = "lm") +
  theme_minimal() + ggtitle("Prespawn survival and GDD by stream") +
  xlab("Growing degree days (GDD)") +
  ylab("Prespawn survival")

survival_model_data |>
  ggplot(aes(x = mean_flow, y = prespawn_survival, fill = stream)) +
  geom_point(aes(color = stream)) + geom_smooth(method = "lm")  +
  theme_minimal() + ggtitle("Prespawn survival and mean flow by stream") +
  xlab("Mean flow (cfs)") +
  ylab("Prespawn survival")

survival_model_data |>
  filter(stream == "mill creek") |>
  ggplot(aes(x = mean_flow, y = prespawn_survival)) +
  geom_point(aes(color = stream))

survival_model_data |>
  filter(stream != "mill creek") |>
  ggplot(aes(x = min_passage_timing, y = prespawn_survival, fill = stream)) +
  geom_point(aes(color = stream)) + geom_smooth(method = "lm")   +
  theme_minimal() + ggtitle("Prespawn survival and minimum passage time by stream") +
  xlab("Minimum passage time (weeks)") +
  ylab("Prespawn survival")

survival_model_data |>
  ggplot(aes(x = total_prop_days_exceed_threshold, y = prespawn_survival, fill = stream)) +
  geom_point(aes(color = stream)) + geom_smooth(method = "lm") +
  facet_wrap(~water_year_type, scales = "free") +
  theme_minimal() + ggtitle("Prespawn survival and temperature by stream and water year type") +
  xlab("Proportion of days exceeding temperature threshold") +
  ylab("Prespawn survival")

survival_model_data |>
  ggplot(aes(x = year, y = prespawn_survival, fill = stream)) +
  geom_point(aes(color = stream)) + geom_smooth(method = "lm") +
  theme_minimal() + ggtitle("Prespawn survival and year")


# check for collinearity for each stream --------------------------------
# function to print Pearsons correlations
print_cors <- function(data, cor_threshold) {
  new_dat <- data |>
    select(-c(prespawn_survival, water_year_type)) |> # can't calculate cor for a categorical variable
    drop_na() |>
    cor() |>
    as.matrix()

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
# 4. use steps 1:3 to select ONE passage timing variable and ONE temperature variable -
# you can't use the glmulti() function with all these variables - it won't converge
# theory - GDD is the standard. stronger than total_prop_days
# 5. use glmuli to look for the best model (by AIC) - including interactions

battle_data <- survival_model_data |>
  filter(stream == "battle creek") |>
  select(-c(year, stream, prop_days_exceed_threshold_migratory,
            prop_days_exceed_threshold_holding))
ggpairs(battle_data)
print_cors(battle_data, 0.65) # mean passage timing is correlated with other passage timing
vif(lm(prespawn_survival ~ ., data = battle_data)) # remove median and mean passage timing
vif(lm(prespawn_survival ~ ., data = battle_data |> select(-c(mean_passage_timing, median_passage_timing)))) # remove total_prop_days_exceed_threshold

# remove variables with highest VIF values
battle_variables_remove <- c("mean_passage_timing", "median_passage_timing", "total_prop_days_exceed_threshold")

# now look for interactions using glmulti
best_battle_model <- glmulti(y = "prespawn_survival",
                             xr = battle_data |> select(-c("prespawn_survival",
                                                           all_of(battle_variables_remove))) |>
                               names(),
                             intercept = TRUE,
                             method = "h",
                             level = 2,
                             data = battle_data,
                             fitfunction = "lm")

# clear
clear_data <- survival_model_data |>
  filter(stream == "clear creek") |>
  select(-c(year, stream, prop_days_exceed_threshold_migratory,
            prop_days_exceed_threshold_holding))
ggpairs(clear_data)
print_cors(clear_data, 0.65) # mean passage timing is correlated with other passage timing
vif(lm(prespawn_survival ~ ., data = clear_data |> select(-c(mean_passage_timing, median_passage_timing,
                                                             water_year_type)))) # variation inflaction factor - if over 5, concerning
vif(lm(prespawn_survival ~ ., data = clear_data |> select(-c(min_passage_timing, median_passage_timing,
                                                             water_year_type)))) # not as good as above
vif(lm(prespawn_survival ~ ., data = clear_data |> select(-c(min_passage_timing, mean_passage_timing,
                                                             water_year_type))))
vif(lm(prespawn_survival ~ ., data = clear_data |> select(-c(min_passage_timing, mean_passage_timing,
                                                             water_year_type, total_prop_days_exceed_threshold))))
vif(lm(prespawn_survival ~ ., data = clear_data |> select(-c(min_passage_timing, mean_passage_timing,
                                                             water_year_type, gdd)))) # gdd is better

# remove variables with highest VIF values
clear_variables_remove <- c("median_passage_timing", "mean_passage_timing", "water_year_type",
                            "total_prop_days_exceed_threshold")

# now look for interactions using glmulti
best_clear_model <- glmulti(y = "prespawn_survival",
                            xr = clear_data |> select(-c("prespawn_survival",
                                                         all_of(clear_variables_remove))) |>
                              names(),
                            intercept = TRUE,
                            method = "h",
                            level = 2,
                            data = clear_data,
                            fitfunction = "lm")

# mill
mill_data <- survival_model_data |>
  filter(stream == "mill creek") |>
  select(-c(year, stream, prop_days_exceed_threshold_migratory,
            prop_days_exceed_threshold_holding,
            median_passage_timing, mean_passage_timing, min_passage_timing)) # mill does not have passage timing

ggpairs(mill_data)
print_cors(mill_data, 0.65) # mean passage timing is correlated with other passage timing
vif(lm(prespawn_survival ~ ., data = mill_data))
vif(lm(prespawn_survival ~ ., data = mill_data |> select(-c(mean_flow))))
vif(lm(prespawn_survival ~ ., data = mill_data |> select(-c(total_prop_days_exceed_threshold))))

summary(lm(prespawn_survival ~ mean_flow, data = mill_data))$r.squared # mean flow is worse
summary(lm(prespawn_survival ~ total_prop_days_exceed_threshold, data = mill_data))$r.squared
summary(lm(prespawn_survival ~ gdd, data = mill_data))$r.squared

# remove variables with highest VIF values
mill_variables_remove <- c("total_prop_days_exceed_threshold")

# now look for interactions using glmulti
best_mill_model <- glmulti(y = "prespawn_survival",
                           xr = mill_data |> select(-c("prespawn_survival",
                                                       all_of(mill_variables_remove))) |>
                             names(),
                           intercept = TRUE,
                           method = "h",
                           level = 2,
                           data = mill_data,
                           fitfunction = "lm")
# plot best models
summary(best_battle_model)$bestmodel
best_battle_lm <- lm(prespawn_survival ~ 1 + min_passage_timing + gdd:mean_flow,
                     data = battle_data)
avPlots(best_battle_lm)

summary(best_clear_model)$bestmodel
best_clear_lm <- lm(prespawn_survival ~ 1 + mean_flow + min_passage_timing + gdd +
                      min_passage_timing:mean_flow + gdd:min_passage_timing, #+ gdd:mean_flow + ,
                     data = clear_data)
avPlots(best_clear_lm)

summary(best_mill_model)$bestmodel
best_mill_lm <- lm(prespawn_survival ~ 1 + mean_flow + water_year_type:mean_flow,
                     data = mill_data)
avPlots(best_mill_lm)




# build bayesian models for each stream --------------------------------------------------------------

battle_bayes <- stan_lm(prespawn_survival ~ min_passage_timing + mean_flow * gdd,
                        data = battle_data,
                        prior = R2(0.1))

battle_bayes |>
  gather_draws(`(Intercept)`, min_passage_timing, mean_flow, gdd, `mean_flow:gdd`, sigma) |>
  median_qi()

clear_bayes <- stan_lm(prespawn_survival ~ mean_flow + min_passage_timing + gdd +
                         min_passage_timing * mean_flow,
                       data = clear_data,
                       prior = R2(0.1))

clear_bayes |>
  gather_draws(`(Intercept)`, min_passage_timing, mean_flow, gdd, `mean_flow:min_passage_timing`, sigma) |>
  median_qi()

mill_bayes <- stan_lm(prespawn_survival ~ mean_flow + water_year_type * mean_flow,
                      data = mill_data,
                      prior = R2(0.1))
mill_bayes |>
  gather_draws(`(Intercept)`, mean_flow, water_year_typewet, `mean_flow:water_year_typewet`) |>
  median_qi()


# try basic bayesian model w all streams ----------------------------------
all_streams_data <- survival_model_data |>
  drop_na() |>
  mutate(prespawn_survival = round(prespawn_survival, 2)) |>
  select(prespawn_survival, gdd, mean_flow, stream, min_passage_timing) |>
  #filter(stream != "mill creek") |> # to try passage timing
  glimpse()

# TODO try out different covariates here. Can we get this to work?
all_streams_bayes <- stan_glmer(
    prespawn_survival ~ (1 | stream),
    data = all_streams_data, family = gaussian,
    chains = 4, iter = 5000*2, seed = 84735
)

mcmc_areas(all_streams_bayes)
mcmc_dens_overlay(all_streams_bayes)
neff_ratio(all_streams_bayes)# should be >0.1
mcmc_acf(all_streams_bayes)
rhat(all_streams_bayes)

all_streams_bayes |>
  spread_draws(`(Intercept)`, sigma) |>
  median_qi()

all_streams_bayes |>
  spread_draws(b[, stream]) |>
  median_qi()

all_streams_bayes |>
  spread_draws(b[, stream]) |>
  summarise_draws()

hierarchical_bayes_results <- all_streams_bayes |>
  spread_draws(`(Intercept)`, b[, stream]) |>
  mutate(condition_mean = `(Intercept)` + b) |>
  median_qi(condition_mean) |>
  select(stream, condition_mean, lower = .lower, upper = .upper) |>
  glimpse()

# predictions and further documentation
# https://cran.r-project.org/web/packages/tidybayes/vignettes/tidy-rstanarm.html

all_streams_bayes |>
  posterior_predict(draws = 1000) |>
  ppc_stat_grouped(y = all_streams_data$prespawn_survival,
                   group = all_streams_data$stream,
                   stat = "median")



# modified dauphin et al --------------------------------------------------

single_stream_dauphin <- "
  data {
    int N;
    int upstream_count[N];
    int redd_count[N];
    real prespawn_survival[N];
    real environmental_index[N];
  }
  parameters {
    real <lower = 0> mu_k;
    real <lower = 0> sigma_k;
    real mu_a;
    real <lower = 0> sigma_a;
  }
  model {
    // priors
    mu_k ~ gamma(3,1);
    sigma_k ~ gamma(1,2);
    mu_a ~ uniform(0, 1000);
    sigma_a ~ gamma(0.001,0.001);

    real alpha[N];
    real beta;

    beta = mu_k * sigma_k;

    vector[N] lambda;

    // calibration between upstream adults and redd count
    for(i in 1:N) {
      alpha[i] = mu_k * beta * environmental_index[N];
      upstream_count[i] ~ lognormal(mu_a, sigma_a);
      prespawn_survival[i] ~ gamma(alpha[i], beta);
      lambda[i] = upstream_count[i] * prespawn_survival[i];
      redd_count[i] ~ poisson(lambda[i]);
    }

  }"

prep_data_dauphin <- function(data, stream_name) {
  new_dat <- data |>
    drop_na() |>
    filter(stream == stream_name) |>
    mutate(standard_passage_timing = min_passage_timing / sum(min_passage_timing),
           standard_gdd = gdd / sum(gdd),
           standard_flow = mean_flow / sum(mean_flow)) |>
    select(year, upstream_count, redd_count, prespawn_survival,
           standard_passage_timing, standard_gdd, standard_flow, mean_flow)

  return(list(N = length(new_dat$year),
              redd_count = new_dat$redd_count,
              upstream_count = new_dat$upstream_count,
              prespawn_survival = new_dat$prespawn_survival,
              environmental_index = new_dat$mean_flow))
}


battle_dauphin <- stan(model_code = single_stream_dauphin,
                 data = prep_data_dauphin(survival_model_data_raw, "battle creek"),
                 chains = 4, iter = 5000*2, seed = 84735)

mcmc_trace(battle_dauphin, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a"))
mcmc_areas(battle_dauphin, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a"))
mcmc_dens_overlay(battle_dauphin, pars = c("mu_k", "sigma_k",  "mu_a", "sigma_a")) # should be indistinguishable
neff_ratio(battle_dauphin, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should be >0.1
mcmc_acf(battle_dauphin, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should drop to be low
rhat(battle_dauphin, c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should be close to 1


# scratch -----------------------------------------------------------------



# use mean flow (based on m5 r squared > m6 ?)
# recommended models - variation among streams


# custom function for trying different models (replaced by step) ----------

# model 1: temperature
# model 2: flow
# model 3: median passage timing
# model 4: water year type
# model 5: temperature + flow
# model 6: temperature + median passage timing
# model 7: temperature + water year type
# model 8: flow + median passage timing
# model 9: flow + water year type
# model 10: temperature + flow + median passage timing
# model 11: temperature + flow + median passage timing + water year type
# model 12: temperature + median passage timing + water year type
# model 13: temperature + flow + water year type
# TODO: interactions? random effects? water year type?
# TODO model 13

ggpairs(battle_model_data)

compare_models <- function(stream_name) {
  r_squared_tibble <- tibble("model" = seq(1:13),
                             "r_squared" = rep(NA, 13),
                             "stream" = stream_name)

  stream_model_data <- survival_model_data |>
    filter(stream == stream_name)


  m1 <- lm(prespawn_survival ~ total_prop_days_exceed_threshold,
           data = stream_model_data)
  r_squared_tibble$r_squared[1] <- summary(m1)$r.squared

  m2 <- lm(prespawn_survival ~ mean_flow,
           data = stream_model_data)
  r_squared_tibble$r_squared[2] <- summary(m2)$r.squared

  m3 <- lm(prespawn_survival ~ median_passage_timing,
           data = stream_model_data)
  r_squared_tibble$r_squared[3] <- summary(m3)$r.squared

  m4 <- lm(prespawn_survival ~ water_year_type,
           data = stream_model_data)
  r_squared_tibble$r_squared[4] <- summary(m4)$r.squared

  m5 <- lm(prespawn_survival ~ total_prop_days_exceed_threshold + mean_flow,
           data = stream_model_data)
  r_squared_tibble$r_squared[5] <- summary(m5)$r.squared

  m6 <- lm(prespawn_survival ~ total_prop_days_exceed_threshold + median_passage_timing,
           data = stream_model_data)
  r_squared_tibble$r_squared[6] <- summary(m6)$r.squared

  m7 <- lm(prespawn_survival ~ total_prop_days_exceed_threshold + water_year_type,
           data = stream_model_data)
  r_squared_tibble$r_squared[7] <- summary(m7)$r.squared

  m8 <- lm(prespawn_survival ~ mean_flow + median_passage_timing,
           data = stream_model_data)
  r_squared_tibble$r_squared[8] <- summary(m8)$r.squared

  m9 <- lm(prespawn_survival ~ mean_flow + water_year_type,
           data = stream_model_data)
  r_squared_tibble$r_squared[9] <- summary(m9)$r.squared

  m10 <- lm(prespawn_survival ~ total_prop_days_exceed_threshold + mean_flow + median_passage_timing,
            data = stream_model_data)
  r_squared_tibble$r_squared[10] <- summary(m10)$r.squared

  m11 <- lm(prespawn_survival ~ total_prop_days_exceed_threshold + mean_flow + median_passage_timing + water_year_type,
            data = stream_model_data)
  r_squared_tibble$r_squared[11] <- summary(m11)$r.squared

  m12 <- lm(prespawn_survival ~ total_prop_days_exceed_threshold + median_passage_timing + water_year_type,
            data = stream_model_data)
  r_squared_tibble$r_squared[12] <- summary(m12)$r.squared

  m13 <- lm(prespawn_survival ~ total_prop_days_exceed_threshold + mean_flow + water_year_type,
            data = stream_model_data)
  r_squared_tibble$r_squared[13] <- summary(m13)$r.squared

  # for(i in 1:11) {
  #   model_name <- get((paste0("m",i)))
  #   r_squared_tibble$r_squared[i] <- summary(model_name)$r.squared
  # }
  #
  # if(stream == "mill creek"){
  #   stream_model_data$median_passage_timing = NA
  # }

  return(arrange(r_squared_tibble, desc(r_squared)))
}

streams <- c("battle creek",
             "clear creek",
             "mill creek")

all_streams <- purrr::map(streams, compare_models) |>
  reduce(bind_rows)

compare_models("battle creek")




# simple linear regression: prespawn mortality vs temperature

# provide exploratory plots to evaluate the variation across streams

# plots showing relationships between redd and passage, ratio to temp (or redd to temp, passage to temp)
# cumulative normal plots of passage
# additional thinking and exploration of year effect


# more scratch ------------------------------------------------------------

# use step function: https://www.statology.org/multiple-linear-regression-r/


# decide on which param should represent temp & passage ---------------------------------------------


# BATTLE
intercept_only <- lm(prespawn_survival ~ 1, data = battle_data)
all <- lm(prespawn_survival ~ ., data = battle_data |>
            filter(!is.na(mean_flow)))
# forward <- step(intercept_only,
#                 direction = "forward",
#                 scope = formula(all))
#
# backward <- step(all,
#                  direction = "backward",
#                  scope = formula(all))

both_directions <- step(intercept_only,
                        direction = "both",
                        scope = formula(all))
both_directions$anova
both_directions$coefficients

# min passage timing

battle_model <- lm(prespawn_survival ~ min_passage_timing,
                   data = battle_data)
summary(battle_model)
battle_data |>
  ggplot(aes(x = min_passage_timing, y = prespawn_survival)) + geom_point() +
  geom_smooth(method = "lm") + theme_minimal() +
  ggtitle("Linear regression of minimum passage week on prespawn survival - Battle Creek") +
  xlab("Minimum passage timing (julian week)") +
  ylab("Prespawn survival")

# CLEAR
intercept_only <- lm(prespawn_survival ~ 1, data = clear_data)
all <- lm(prespawn_survival ~ ., data = clear_data)

both_directions <- step(intercept_only,
                        direction = "both",
                        scope = formula(all))
both_directions$anova
both_directions$coefficients

summary(lm(prespawn_survival ~ gdd, data = clear_data))$r.squared
summary(lm(prespawn_survival ~ mean_flow, data = clear_data))$r.squared

# gdd
clear_model <- lm(prespawn_survival ~ gdd,
                  data = clear_data)
summary(clear_model)
clear_data |>
  ggplot(aes(x = gdd, y = prespawn_survival)) + geom_point() +
  geom_smooth(method = "lm") + theme_minimal() +
  ggtitle("Linear regression of GDD on prespawn survival - Clear Creek") +
  xlab("Growing degree days") +
  ylab("Prespawn survival")

# MILL
intercept_only <- lm(prespawn_survival ~ 1, data = mill_data |> drop_na())
all <- lm(prespawn_survival ~ ., data = mill_data |> drop_na())

both_directions <- step(intercept_only,
                        direction = "both",
                        scope = formula(all))
both_directions$anova
both_directions$coefficients

summary(lm(prespawn_survival ~ 1, data = mill_data |> drop_na()))$r.squared
summary(lm(prespawn_survival ~ mean_flow + gdd, data = mill_data |> drop_na()))$r.squared

# according to backward, water_year_type and mean_flow
mill_model <- lm(prespawn_survival ~ water_year_type + mean_flow,
                 data = mill_data)
summary(mill_model)
avPlots(mill_model)
mill_data |>
  ggplot(aes(x = mean_flow, y = prespawn_survival)) + geom_point() +
  geom_smooth(method = "lm") + theme_minimal() +
  ggtitle("Linear regression of mean flow on prespawn survival - Clear Creek") +
  xlab("Mean flow (cfs)") +
  ylab("Prespawn survival")


# scratch for temp index -----------------------------------------------------------------

# linear model of prespawn mortality (upstream - redd) to get intercept
upstream_passage_clear <- upstream_passage |>
  filter(stream %in% c("battle creek", "clear creek", "yuba river")) |>
  rename(upstream_count = count)

redd_clear <- redd |>
  filter(stream %in% c("battle creek", "clear creek", "yuba river")) |>
  rename(redd_count = count)

clear_prespawn <- upstream_passage_clear |>
  left_join(redd_clear, by = c("year", "stream")) |>
  mutate(female_upstream = upstream_count * 0.5,
         prespawn_survival = redd_count / female_upstream,
         prespawn_survival = ifelse(prespawn_survival > 1, 1, prespawn_survival)) |>
  glimpse()

temp_prespawn <- left_join(clear_prespawn, migratory_temp |>
                             select(year, prop_days_exceed_threshold_migratory),
                           by = "year") |>
  left_join(holding_temp |>
              select(year, stream, prop_days_exceed_threshold_holding)) |>
  mutate(total_prop_days_exceed_threshold = ifelse(is.na(prop_days_exceed_threshold_holding), prop_days_exceed_threshold_migratory,
                                                   prop_days_exceed_threshold_holding + prop_days_exceed_threshold_migratory)) |>
  glimpse()

temp_prespawn |>
  ggplot(aes(x = total_prop_days_exceed_threshold, y = prespawn_survival, fill = stream)) +
  geom_point(aes(color = stream)) + geom_smooth(method = "lm")

m <- lm(data = temp_prespawn, prespawn_survival ~ total_prop_days_exceed_threshold + stream)
summary(m)

# TODO try median 7-day maximum temps (pull from gauges)

surv_factor <- coef(m)[1] # TODO find estimates

temp_prespawn_scaled <- temp_prespawn |>
  select(year, stream, prespawn_survival, total_prop_days_exceed_threshold) |>
  mutate(scaled_prop_days_exceed = total_prop_days_exceed_threshold * surv_factor) |>
  glimpse()

temp_index <- temp_prespawn_scaled |>
  group_by(year) |>
  summarise(temp_index = mean(scaled_prop_days_exceed, na.rm = T)) |>
  ungroup() |>
  mutate(temp_index = ifelse(is.na(temp_index), 1, temp_index))

