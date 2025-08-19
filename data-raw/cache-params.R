# create parameters to pass to BUGS for BT-SPAS-X
bt_spas_x_bayes_params <- list(number_mcmc = 10000,
                               number_burnin = 5000,
                               number_thin = 5,
                               number_chains = 3)

usethis::use_data(bt_spas_x_bayes_params, overwrite = TRUE)

# JPE forecasting
forecast_seed <- 1234
usethis::use_data(forecast_seed, overwrite = TRUE)

forecast_sites <- c("ubc", "ucc", "mill creek", "deer creek", "okie dam")
usethis::use_data(forecast_sites, overwrite = TRUE)
