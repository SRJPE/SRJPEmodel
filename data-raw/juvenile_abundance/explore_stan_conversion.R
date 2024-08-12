# test bt spas x all mark recap on chatgpt STAN code

inputs <- readRDS("~/Downloads/inputs.rds")
inputs$data

inits <- inputs$inits[[1]]
names(inits) <- c("trib_mu_P", "b0_pCap", "flow_mu_P",
                  "b_flow", "trib_tau_P", "flow_tau_P",
                  "pro_tau_P", "b_sp", "lg_N")
inits$b0_pCap <- c(-4.25399808352841, -4.37002062092314, -4.87083612652123, -4.75584830999987,
                   -5, -5, -2.78255068642205, -4.72254190733127, -3.81263454572515)

inits <- list(inits, inits, inits)

# subset_data <- inputs$data[c("Ntribs", "Nmr", "Nstrata", "Nstrata_wc",
#                              "ind_trib", "mr_flow", "Recaptures", "Releases",
#                              "Uwc_ind", "Nstrata_wc", "u", "K",
#                              "lgN_max", "ind_pCap", "ZP")]

stan_model <- read_file("model_files/no_cut_bt_spas/missing_mark_recap.stan")
stan_model <- read_file("model_files/no_cut_bt_spas/missing_mark_recap_pCap_only.stan")

options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)

test <- rstan::stan(model_name = "ubc_2004_stan_test",
                    model_code = stan_model,
                    data = inputs$data,
                    init = inits,
                    chains = 3, iter = 1000, seed = 84735)

parlist <- c("trib_mu_P", "trib_sd_P", "flow_mu_P", "pro_sd_P", "flow_sd_P")

print(rstan::summary(test, pars = "b0_pCap")$summary)
