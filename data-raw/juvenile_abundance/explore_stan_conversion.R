# test bt spas x all mark recap on chatgpt STAN code

inputs <- readRDS("~/Downloads/inputs.rds")
inputs$data

# subset_data <- inputs$data[c("Ntribs", "Nmr", "Nstrata", "Nstrata_wc",
#                              "ind_trib", "mr_flow", "Recaptures", "Releases",
#                              "Uwc_ind", "Nstrata_wc", "u", "K",
#                              "lgN_max", "ind_pCap", "ZP")]

stan_model <- read_file("model_files/no_cut_bt_spas/missing_mark_recap.stan")

test <- rstan::stan(model_name = "ubc_2004_stan_test",
                    model_code = stan_model,
                    data = inputs$data,
                    chains = 3, iter = 10, seed = 84735)
