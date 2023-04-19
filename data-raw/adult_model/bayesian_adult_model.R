# implement modified Dauphin et al model

# libraries ---------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(rstan)
library(bayesplot)

source(here::here("analysis", "adult_model_data_prep.R"))


# TODO input is prop_days_over_threshold * coefficient * ratio_k

# write base model in stan --------------------------------------------------

mill_model_temp_prop <- "
  data {
    int N;
    int upstream_count[N];
    int redd_count[N];
    real ratio_k[N];
    real temp_index[N];
    real prop_surveyed[N];
  }
  parameters {
    real <lower = 0> mu_k;
    real <lower = 0> sigma_k;
    real mu_a;
    real <lower = 0> sigma_a;
  }
  model {
    vector[N] lambda;
    real alpha[N];
    real beta;
    beta = mu_k * sigma_k;
    for(i in 1:N) { 
      alpha[i] = mu_k * beta * temp_index[i];
      upstream_count[i] ~ lognormal(mu_a, sigma_a);
      ratio_k[i] ~ gamma(alpha[i], beta);
      lambda[i] = upstream_count[i] * ratio_k[i] * prop_surveyed[i];
      redd_count[i] ~ poisson(lambda[i]);
    }

  }"

# prep data ---------------------------------------------------------------

all_mill <- redd |> 
  filter(stream == "mill creek") |>
  rename(redd_count = count) |> 
  full_join(upstream_passage |> 
              filter(stream == "mill creek") |> 
              rename(upstream_count = count) |> 
              select(-stream),
            by = "year") |> 
  glimpse()

years <- all_mill |> filter(!is.na(redd_count) & !is.na(upstream_count)) |> 
                              pull(year)
N <- length(years)
redd_count <- all_mill |> 
  filter(year %in% years) |> 
  pull(redd_count)
upstream_count <- all_mill |> 
  filter(year %in% years) |> 
  pull(upstream_count)
ratio_k <- upstream_count / redd_count
temp_index <- temp_index |> 
  filter(year %in% years) |> 
  pull(temp_index)
prop_surveyed <- rep(1.0, N) # assume all of spawning population is surveyed



# fit model ---------------------------------------------------------------
fit <- stan(model_code = mill_model_temp_prop, 
            data = list(N = N,
                        
                        upstream_count = upstream_count,
                        redd_count = redd_count,
                        ratio_k = ratio_k,
                        temp_index = temp_index,
                        prop_surveyed = prop_surveyed), 
            #init = init_list,
            chains = 4, iter = 5000*2, seed = 84735)

mcmc_trace(fit, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a"))
mcmc_areas(fit, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a"))
mcmc_dens_overlay(fit, pars = c("mu_k", "sigma_k",  "mu_a", "sigma_a")) # should be indistinguishable
neff_ratio(fit, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should be >0.1
mcmc_acf(fit, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should drop to be low
rhat(fit, c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should be close to 1


# fix bayesian logic -------------------------------------------------

test <- "
  data {
    int N;
    int upstream_count[N];
    int redd_count[N];
    real ratio_k[N];
    real temp_index[N];
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
      alpha[i] = mu_k * beta * temp_index[i];
      upstream_count[i] ~ lognormal(mu_a, sigma_a);
      ratio_k[i] ~ gamma(alpha[i], beta);
      lambda[i] = upstream_count[i] * ratio_k[i];
      redd_count[i] ~ poisson(lambda[i]);
    }

  }"

test_fit <- stan(model_code = test, 
                  data = list(N = N,
                              upstream_count = upstream_count,
                              redd_count = redd_count,
                              ratio_k = ratio_k,
                              temp_index = temp_index), 
                  #init = init_list,
                  chains = 4, iter = 5000*2, seed = 84735)

mcmc_trace(test_fit, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a"))
mcmc_areas(test_fit, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a"))
mcmc_dens_overlay(test_fit, pars = c("mu_k", "sigma_k",  "mu_a", "sigma_a")) # should be indistinguishable
neff_ratio(test_fit, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should be >0.1
mcmc_acf(test_fit, pars = c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should drop to be low
rhat(test_fit, c("mu_k", "sigma_k", "mu_a", "sigma_a")) # should be close to 1


# try missing data --------------------------------------------------------
test_missing <- "
  data {
    int N;
    int Y;
    int upstream_count[N];
    int redd_count[N];
    int redd_missing[Y];
    real ratio_k[N];
    real temp_index[N];
    real temp_index_old[Y];
  }
  parameters {
    real <lower = 0> mu_k;
    real <lower = 0> sigma_k;
    real mu_a;
    real <lower = 0> sigma_a;
    real upstream_missing[Y];
  }
  model {
    mu_k ~ gamma(3,1);
    sigma_k ~ gamma(1,2);
    mu_a ~ uniform(0, 1000);
    sigma_a ~ gamma(0.001,0.001);
    
    for(i in 1:Y) {
      upstream_missing[i] ~ uniform(0, 1000);
    }
    
    real alpha[N];
    real beta;
    
    beta = mu_k * sigma_k;

    vector[N] lambda;
    vector[Y] lambda_old;
    vector[Y] ratio_k_old;
    vector[Y] alpha_old;
    
    for(i in 1:N) { 
      alpha[i] = mu_k * beta * temp_index[i];
      upstream_count[i] ~ lognormal(mu_a, sigma_a);
      ratio_k[i] ~ gamma(alpha[i], beta);
      lambda[i] = upstream_count[i] * ratio_k[i];
      redd_count[i] ~ poisson(lambda[i]);
    }
    
    for(i in 1:Y){
      alpha_old[i] = mu_k * beta * temp_index_old[i];
      upstream_missing[i] ~ lognormal(mu_a, sigma_a);
      ratio_k_old[i] ~ gamma(alpha_old[i], beta);
      lambda_old[i] = upstream_missing[i] * ratio_k_old[i];
      redd_missing[i] ~ poisson(lambda_old[i]);
      
    }

  }"

years_missing <- all_mill |> 
  filter(is.na(upstream_count)) |> 
  pull(year)
Y <- length(years_missing)
redd_missing <- all_mill |> 
  filter(year %in% years_missing) |> 
  pull(redd_count)

temp_index_old <- temp_prespawn_scaled |> 
  group_by(year) |> 
  summarise(temp_index = mean(scaled_prop_days_exceed, na.rm = T)) |> 
  ungroup() |> 
  mutate(temp_index = ifelse(is.na(temp_index), 1, temp_index)) |> 
  filter(year %in% years_missing) |> 
  pull(temp_index) |> 
  append(1, 0)

test_missing_fit <- stan(model_code = test_missing, 
                         data = list(N = N,
                                     Y = Y,
                                     upstream_count = upstream_count,
                                     redd_count = redd_count,
                                     redd_missing = redd_missing,
                                     ratio_k = ratio_k,
                                     temp_index = temp_index,
                                     temp_index_old = temp_index_old), 
                         #init = init_list,
                         chains = 4, iter = 5000*2, seed = 84735)


# dauphin et al model code ------------------------------------------------

model{
  
  # Prior distributions of the model parameters
  tau.p ~dgamma(0.001,0.001)
  mu.p~dbeta(1,1)I(0.001,)
  L.mu.p <- logit(mu.p)
  
  mu.kappa ~ dgamma(1,0.001)
  tau.kappa ~ dgamma(0.001,0.001)I(0.001,)
  alpha <- mu.kappa * beta
  beta <- mu.kappa * tau.kappa
  
  for (i in 1:3){
    mu.A[i]~dunif(0,500000)
    L.mu.A[i]<-log(mu.A[i])
    tau.A[i]~dgamma(0.001,0.001)
  }
  
  # calibration between adult counts and redd counts
  # T = 8 years, 2001 to 2008; i = rivers (1=Faughan; 2=Finn; 3=Roe)
  for (t in 1:T){
    for (i in 1:3){						
      L.p[t,i] ~ dnorm(L.mu.p, tau.p)
      logit(p[t,i]) <- L.p[t,i] 
      A[t,i] ~ dlnorm(L.mu.A[i],tau.A[i])
      S[t,i] <- A[t,i] - C[t,i]
      kappa[t,i] ~ dgamma(alpha, beta)I(0.001,)
      lambda[t,i] <- (A[t,i] - C[t,i]) * p[t,i] * kappa[t,i]					
      R[t,i] ~dpois(lambda[t,i])
    }					
  }
  # time series
  # Y = 42 years, 1959 to 2000; i = rivers (1=Faughan; 2=Finn; 3=Roe)
  for (y in 1:Y){
    for (i in 1:3){
      L.p.old[y,i] ~ dnorm(L.mu.p, tau.p)
      logit(p.old[y,i]) <- L.p.old[y,i] 
      
      A.old[y,i] ~ dlnorm(L.mu.A[i],tau.A[i])
      
      kappa.old[y,i] ~ dgamma(alpha, beta)I(0.001,)			
      lambda.redds.old[y,i] <- (A.old[y,i] – C[t,i]) * p.old[y,i] * kappa.old[y,i]			
      Redds.old[y,i] ~dpois(lambda.redds.old[y,i])   
    }
    
  }
}
# try hierarchical structure ----------------------------------------------

# prep data
all_tribs_all_data <- upstream_passage |> 
  rename(upstream_count = count) |> 
  full_join(redd |> rename(redd = count),
            by = c("year", "stream")) |> 
  full_join(holding |> rename(holding = count),
            by = c("year", "stream")) |> 
  select(year, stream, upstream = upstream_count, redd, holding) |> 
  glimpse()

# just try redd and upstream for battle & yuba
upstream_hierarchical <- all_tribs_all_data |> 
  filter(stream %in% c("yuba river", "battle creek"),
         year %in% 2011:2019) |> 
  select(year, upstream, stream) |> 
  pivot_wider(names_from = stream, values_from = upstream) |> 
  select(-year) |> 
  as.matrix() |> 
  unname()

redd_hierarchical <- all_tribs_all_data |> 
  filter(stream %in% c("yuba river", "battle creek"),
         year %in% 2011:2019) |>
  select(year, redd, stream) |> 
  pivot_wider(names_from = stream, values_from = redd) |> 
  select(-year) |> 
  as.matrix() |> 
  unname()

ratio_hierarchical <- upstream_hierarchical / redd_hierarchical

S <- dim(redd_hierarchical)[2]
N <- dim(redd_hierarchical)[1]

mill_model_hierarchical <- "
  data {
    int N;
    int S;
    int upstream_count[N,S];
    int redd_count[N,S];
    real ratio_k[N,S];
  }
  parameters {
    real <lower = 0> mu_k;
    real <lower = 0> sigma_k[N];
    real mu_a;
    real <lower = 0> sigma_a;
  }
  model {
    matrix[N,S] lambda;
    real alpha[N];
    real beta[N];
    for(i in 1:N) {
      beta[i] = mu_k * sigma_k[i];
      alpha[i] = mu_k * beta[i];
      for(j in 1:S) {
        upstream_count[i,j] ~ lognormal(mu_a, sigma_a);
        ratio_k[i,j] ~ gamma(alpha[i], beta[i]);
        lambda[i,j] = upstream_count[i,j] * ratio_k[i,j];
        redd_count[i,j] ~ poisson(lambda[i,j]);
      }
    }
  }"

fit_hierarchical <- stan(model_code = mill_model_hierarchical, 
                          data = list(N = N,
                                      S = S,
                                      upstream_count = upstream_hierarchical,
                                      redd_count = redd_hierarchical,
                                      ratio_k = ratio_hierarchical), 
                          chains = 4, iter = 5000*2, seed = 84735)

# diagnostics -------------------------------------------------------------

mcmc_trace(fit_hierarchical, pars = c("mu_k", "mu_a", "sigma_a"))
mcmc_areas(fit_hierarchical, pars = c("mu_k", "mu_a", "sigma_a"))
mcmc_dens_overlay(fit_hierarchical, pars = c("mu_k",  "mu_a", "sigma_a")) # should be indistinguishable
neff_ratio(fit_hierarchical, pars = c("mu_k", "mu_a", "sigma_a")) # should be >0.1
mcmc_acf(fit_hierarchical, pars = c("mu_k",  "mu_a", "sigma_a")) # should drop to be low
rhat(fit_hierarchical, c("mu_k",  "mu_a", "sigma_a")) # should be close to 1


# scratch -----------------------------------------------------------------

# prediction
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



# scratch -----------------------------------------------------------------



mu <- 0.5
sigma <- 0.2
beta <- mu * sigma
alpha <- mu * beta
hist(rgamma(10000, alpha, beta))

init_list <- list()
Nchains <- 4
for(i in 1:Nchains){
  init_list[[i]] <- list(mu_k = rgamma(1, 3, 1),
                         sigma_k = rgamma(1, 1, 2),
                         mu_a = runif(1, 3, 5),
                         sigma_a = rgamma(1, 0.6, 1))
}


# try predicting with bayesforecast
bayesforecast::forecast(fit, h = 20) # not working
