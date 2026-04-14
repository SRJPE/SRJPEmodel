# create parameters to pass to BUGS for BT-SPAS-X
bt_spas_x_bayes_params <- list(number_mcmc = 10000,
                               number_burnin = 5000,
                               number_thin = 5,
                               number_chains = 3)

usethis::use_data(bt_spas_x_bayes_params, overwrite = TRUE)

# JPE forecasting
forecast_seed <- 1234
usethis::use_data(forecast_seed, overwrite = TRUE)

forecast_sites <- c("ubc", "lcc", "mill creek", "deer creek", "okie dam")
usethis::use_data(forecast_sites, overwrite = TRUE)

# site order
site_order_north_south <- tibble("ns_order" = 1:15,
                                 "site" = c("lcc", "ucc", "ubc",
                                            "mill creek", "deer creek",
                                            "okie dam",
                                            "steep riffle", "eye riffle",
                                            "gateway riffle", "herringer riffle",
                                            "live oak", "sunset pumps",
                                            "hallwood", "tisdale", "knights landing"))
usethis::use_data(site_order_north_south, overwrite = TRUE)
