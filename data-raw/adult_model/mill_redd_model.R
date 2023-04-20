# libraries
library(tidyverse)
library(lubridate)
library(googleCloudStorageR)
library(rstan)
library(bayesplot)
library(shinystan)
library(rstanarm)
library(bayesforecast)

# get data ----------------------------------------------------------------

gcs_auth(json_file = Sys.getenv("GCS_AUTH_FILE"))
gcs_global_bucket(bucket = Sys.getenv("GCS_DEFAULT_BUCKET"))

# upstream passage
upstream_passage <- read_csv(gcs_get_object(object_name = "standard-format-data/standard_adult_upstream_passage.csv",
                                            bucket = gcs_get_global_bucket())) |> 
  mutate(stream = tolower(stream)) |> 
  filter(run %in% c("spring", NA, "not recorded")) |> 
  glimpse()

# holding
holding <- read_csv(gcs_get_object(object_name = "standard-format-data/standard_holding.csv",
                                   bucket = gcs_get_global_bucket()))|> 
  group_by(year) |> 
  summarise(count = sum(count, na.rm = T)) |> 
  glimpse()

# redd
redd <- read_csv(gcs_get_object(object_name = "standard-format-data/standard_annual_redd.csv",
                                bucket = gcs_get_global_bucket())) |> 
  group_by(year) |>
  rename(count = max_yearly_redd_count) |> glimpse()

redd_daily <- read_csv(gcs_get_object(object_name = "standard-format-data/standard_daily_redd.csv",
                                      bucket = gcs_get_global_bucket())) |> 
  filter(run %in% c("spring", NA, "not recorded")) |> 
  mutate(count = 1) |> 
  group_by(year, stream) |> 
  summarise(count = sum(count)) |> glimpse()

# temperature

# threshold
# https://www.noaa.gov/sites/default/files/legacy/document/2020/Oct/07354626766.pdf 
threshold <- 21

temp <- read_csv(gcs_get_object(object_name = "standard-format-data/standard_temperature.csv",
                                bucket = gcs_get_global_bucket())) |> 
  filter(stream == "sacramento river") |> 
  filter(month(date) %in% 6:8) |> 
  group_by(year(date)) |> 
  mutate(above_threshold = ifelse(mean_daily_temp_c > threshold, TRUE, FALSE)) |> 
  summarise(prop_days_exceed_threshold = round(sum(above_threshold, na.rm = T)/length(above_threshold), 2)) |> 
  ungroup() |> 
  mutate(prop_days_below_threshold = 1 - prop_days_exceed_threshold) |> 
  rename(year = `year(date)`) |> 
  glimpse()




# join (1999:2020)
data <- holding |> 
  left_join(temp, by = "year") |> 
  filter(!is.na(prop_days_exceed_threshold)) |> 
  glimpse()

# explore
data |> 
  ggplot(aes(x = prop_days_exceed_threshold, y = count)) +
  geom_point()


holding_model <- stan_glm(count ~ prop_days_exceed_threshold, data = data,
                       family = gaussian,
                       #prior_intercept = normal(22, 2),
                       #prior = normal(1, 0.5),
                       #prior_aux = exponential(0.0008),
                       chains = 4, iter = 5000*2, seed = 84735)

launch_shinystan(holding_model)

# poisson: https://mc-stan.org/rstanarm/articles/count.html
# also: https://stats.libretexts.org/Bookshelves/Introductory_Statistics/Book%3A_Introductory_Statistics_(OpenStax)/04%3A_Discrete_Random_Variables/4.07%3A_Poisson_Distribution

holding_model_poisson <- stan_glm(count ~ prop_days_exceed_threshold, data = data,
                                  family = poisson,
                                  #prior_intercept = normal(22, 2),
                                  #prior = normal(1, 0.5),
                                  #prior_aux = exponential(0.0008),
                                  chains = 4, iter = 5000*2, seed = 84735)

launch_shinystan(holding_model_poisson)


# redd bayesian modeling --------------------------------------------------

redd_mill <- redd |> 
  filter(stream == "mill creek") |>
  group_by(year) |> 
  mutate(count = sum(count, na.rm = T)) |> 
  ungroup() |> 
  glimpse()

upstream_mill <- upstream_passage |> 
  filter(stream == "mill creek") |>
  mutate(year = year(date)) |> 
  group_by(year) |> 
  summarise(count = round(sum(count, na.rm = T), 0)) |> 
  filter(!is.na(year)) |> glimpse()

range(redd_mill$year)
range(upstream_mill$year)



# predicted spawning population (count) based on temperature and uncertainty

# TODO go through chapter 9 in bayes book with an eye toward fitting these data
# TODO try poisson distribution in stan_glm()

# prep data
N <- (max(upstream_mill$year) - min(upstream_mill$year) + 1)
redd_count <- redd_mill |> 
  filter(year %in% min(upstream_mill$year):max(upstream_mill$year)) |> 
  distinct(year, count) |> 
  pull(count) |> as.integer()
years <- min(upstream_mill$year):max(upstream_mill$year)
upstream_count <- upstream_mill |> pull(count) |> as.integer()
ratio_k <- upstream_count/redd_count
temp_subset <- temp |> filter(year %in% years) |> pull(prop_days_below_threshold)
temp_subset <- append(temp_subset, 0.15, 2) # add in fake value for 2013 (filtered out)

# write model
mill_upstream_redd_model <- "
  data {
    int N;
    int upstream_count[N];
    int redd_count[N];
    real ratio_k[N];
    real temp[N];
  }
  parameters {
    real mu_k;
    real <lower = 0> sigma_k;
    real mu_a;
    real <lower = 0> sigma_a;
  }
  model {
    vector[N] lambda;
    real alpha;
    real beta;
    beta = mu_k * sigma_k;
    alpha = mu_k * beta;
    for(i in 1:N) {
      upstream_count[i] ~ lognormal(mu_a, sigma_a);
      ratio_k[i] ~ gamma(alpha, beta);
      lambda[i] = upstream_count[i] * ratio_k[i] * temp[N];
      redd_count[i] ~ poisson(lambda[i]);
    }

  }"

mu <- 0.5
sigma <- 0.2
beta <- mu * sigma
alpha <- mu * beta
summary(rgamma(10000, alpha, beta))

init_list <- list()
Nchains <- 4
for(i in 1:Nchains){
  init_list[[i]] <- list(mu_k = rgamma(1, 3, 1),
                         sigma_k = rgamma(1, 1, 2),
                         mu_a = runif(1, 3, 5),
                         sigma_a = rgamma(1, 0.6, 1))
}


fit <- stan(model_code = mill_upstream_redd_model, 
            data = list(N = N,
                        upstream_count = upstream_count,
                        redd_count = redd_count,
                        ratio_k = ratio_k, 
                        temp = temp_subset), 
            #init = init_list,
            chains = 4, iter = 5000*2, seed = 84735)

# diagnostics
mcmc_trace(fit, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a"))
mcmc_areas(fit, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a"))
mcmc_dens_overlay(fit, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should be indistinguishable
neff_ratio(fit, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should be >0.1
mcmc_acf(fit, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should drop to be low
rhat(fit, c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should be close to 1

#launch_shinystan(fit)

# prediction using rstanarm - need to fit model in rstanarm for this to work
# TODO wrapper for passing stan code to rstan?

# prediction with R code and custom function
# get results
pars <- as.data.frame(fit, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a"))
predict_redd <- function(mu_k, sigma_k, mu_a, sigma_a) {
  beta <- mu_k * sigma_k
  alpha <- mu_k * beta
  
  upstream_count <- rlnorm(1, mu_a, sigma_a)
  ratio_k <- rgamma(1, alpha, beta)
  lambda <- upstream_count * ratio_k
  redd_count <- rpois(1, lambda)
  return(list("redd_count" = redd_count,
              "upstream_count" = upstream_count,
              "ratio_k" = ratio_k))
}

pred_redd <- numeric()
pred_adult <- numeric()
pred_ratio <- numeric()
year <- numeric()
for(i in 1:100){
  predict <- predict_redd(sample(pars$mu_k, 1), sample(pars$sigma_k, 1),
                               sample(pars$mu_a, 1), sample(pars$sigma_a, 1))
  year[i] <- i
  pred_redd[i] <- predict$redd_count
  pred_adult[i] <- predict$upstream_count
  pred_ratio[i] <- predict$ratio_k
}

pred_dat <- tibble(year = year,
                   redd = pred_redd,
                   adult = pred_adult,
                   ratio = pred_ratio)
ggplot(pred_dat) +
  geom_line(data = pred_dat, aes(x = year, y = redd)) + 
  geom_line(data = pred_dat, aes(x = year, y = adult),  col = "blue")

ggplot(pred_dat) +
  geom_line(aes(x = year, y = pred_ratio))
# TODO could covariate explain years where ratio is below 1 ? upstream passage has 
# high flow? 

# try predicting with bayesforecast
bayesforecast::forecast(fit, h = 20) # not working
