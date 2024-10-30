# create parameters to pass to BUGS for BT-SPAS-X
bt_spas_x_bayes_params <- list(number_mcmc = 10000,
                               number_burnin = 5000,
                               number_thin = 5,
                               number_chains = 3)

usethis::use_data(bt_spas_x_bayes_params, overwrite = TRUE)
