# data
input <- run_single_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                              SRJPEdata::weekly_juvenile_abundance_catch_data,
                              SRJPEdata::weekly_juvenile_abundance_efficiency_data,
                              "lcc", 2018, NA, T, here::here("data-raw", "WinBUGS14"),
                              no_cut =F)

input$data |>
  names()

inits <- input$inits[[1]][-8]
inits <- inits[-8]
inits <- list(inits, inits, inits)

options(mc.cores=parallel::detectCores())
pcap_missing <- rstan::stan(model_code = read_file(here::here("model_files",
                                                            "pCap_missing_mark_recap.stan")),
                          data = input$data,
                          init = inits,
                          chains = SRJPEmodel::bt_spas_x_bayes_params$number_chains,
                          iter = SRJPEmodel::bt_spas_x_bayes_params$number_mcmc,
                          seed = 84735)

print(rstan::summary(pcap_missing,pars=c("trib_mu_P", "trib_sd_P", "flow_mu_P",
                                                  "flow_sd_P", "pro_sd_P",
                                                  "b_flow"))$summary)
rstan::summary(pcap_missing, pars = c("pCap_U"))$summary

bayesplot::mcmc_trace(pcap_missing, pars = "pCap_U[3]")

# pars - need mean and sd
inits <- list("b_sp" = input$inits[[1]]$b_sp,
              "lg_N" = input$inits[[1]]$lg_N)
inits <- list(inits, inits, inits)

pCap_model_results <- rstan::summary(pcap_missing, pars = c("pCap_U"))$summary |>
  data.frame() |>
  summarise(all_mean = mean(mean),
            all_sd = mean(sd))

input$data$lt_pCap_mu <- gtools::logit(pCap_model_results$all_mean)
input$data$lt_pCap_sd <- gtools::logit(pCap_model_results$all_sd)

abundance <- rstan::stan(model_code = read_file(here::here("model_files",
                                                           "abundance_model.stan")),
                         data = input$data,
                         init = inits,
                         chains = SRJPEmodel::bt_spas_x_bayes_params$number_chains,
                         iter = SRJPEmodel::bt_spas_x_bayes_params$number_mcmc,
                         seed = 84735)
