# test bt spas x missing mark recap on chatgpt STAN code

ubc_2009_inputs <- readRDS("data-raw/juvenile_abundance/ubc_2009_2024-09-6.rds")
ubc_2009_inputs$data

inits <- ubc_2009_inputs$inits[[1]]
names(inits) <- c("trib_mu_P", "b0_pCap", "flow_mu_P",
                  "b_flow", "trib_tau_P", "flow_tau_P",
                  "pro_tau_P", "b_sp", "lg_N")
# inits$b0_pCap <- c(-4.25399808352841, -4.37002062092314, -4.87083612652123, -4.75584830999987,
#                    -5, -5, -2.78255068642205, -4.72254190733127, -3.81263454572515)

inits <- list(inits, inits, inits)

# subset_data <- inputs$data[c("Ntribs", "Nmr", "Nstrata", "Nstrata_wc",
#                              "ind_trib", "mr_flow", "Recaptures", "Releases",
#                              "Uwc_ind", "Nstrata_wc", "u", "K",
#                              "lgN_max", "ind_pCap", "ZP")]

stan_model <- read_file("model_files/missing_mark_recap.stan")

options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)

test <- rstan::stan(model_name = "ubc_2009_stan_test",
                    model_code = stan_model,
                    data = ubc_2009_inputs$data,
                    init = inits,
                    chains = 3, iter = 1000, seed = 84735)

parlist <- c("trib_mu_P", "trib_sd_P", "flow_mu_P", "pro_sd_P", "flow_sd_P")

print(rstan::summary(test, pars = "b0_pCap")$summary)

test <- run_single_bt_spas_x_stan(SRJPEmodel::bt_spas_x_bayes_params,
                                  weekly_juvenile_abundance_catch_data,
                                  weekly_juvenile_abundance_efficiency_data,
                                  "ubc",
                                  run_year = 2009,
                                  lifestage = "fry",
                                  effort_adjust = F)
