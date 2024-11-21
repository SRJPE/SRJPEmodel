# data
input <- run_single_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                              SRJPEdata::weekly_juvenile_abundance_catch_data,
                              SRJPEdata::weekly_juvenile_abundance_efficiency_data,
                              "lcc", 2018, NA, T, here::here("data-raw", "WinBUGS14"),
                              no_cut =F)
input$data$lgN_max <- input$data$lgN_max + 0.5
input <- run_single_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                              SRJPEdata::weekly_juvenile_abundance_catch_data,
                              SRJPEdata::weekly_juvenile_abundance_efficiency_data,
                              "ubc", 2009, NA, T, here::here("data-raw", "WinBUGS14"),
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
rstan::summary(pcap_missing, pars = c("lt_pCap_U"))$summary

bayesplot::mcmc_trace(pcap_missing, pars = "pCap_U[3]")

# pars - need mean and sd
logit_pCaps <- rstan::summary(pcap_missing, pars = c("lt_pCap_U"))$summary |>
  data.frame()

input$data$lt_pCap_mu <- logit_pCaps$mean
input$data$lt_pCap_sd <- logit_pCaps$sd

inits <- list("b_sp" = input$inits[[1]]$b_sp,
              "lg_N" = input$inits[[1]]$lg_N,
              "lt_pCap_U" = rnorm(input$data$Nstrata_wc,
                                  mean(logit_pCaps$mean),
                                  mean(logit_pCaps$sd)))
inits <- list(inits, inits, inits)

abundance <- rstan::stan(model_code = read_file(here::here("model_files",
                                                           "abundance_model.stan")),
                         data = input$data,
                         init = inits,
                         chains = SRJPEmodel::bt_spas_x_bayes_params$number_chains,
                         iter = SRJPEmodel::bt_spas_x_bayes_params$number_mcmc * 3,
                         seed = 84735)

rstan::summary(abundance,pars=c("lt_pCap_U"))$summary
rstan::summary(abundance,pars=c("N"))$summary
rstan::traceplot(abundance, pars = c("lt_pCap_U"))
rstan::traceplot(pcap_missing, pars = c("pCap_U"))

pcap_efficiency <- rstan::summary(pcap_missing,pars=c("lt_pCap_U"))$summary |>
  data.frame() |>
  select(mean, sd)
abundance_efficiency <- rstan::summary(abundance,pars=c("lt_pCap_U"))$summary |>
  data.frame() |>
  select(mean, sd)

all_efficiency <- bind_rows(pcap_efficiency |>
                              mutate(model = "pcap",
                                     index = row_number(),
                                     mean = plogis(mean)) |>
                              pivot_longer(mean:sd,
                                           names_to = "parameter",
                                           values_to = "value"),
                            abundance_efficiency |>
                              mutate(model = "abundance",
                                     index = row_number(),
                                     mean = plogis(mean)) |>
                              pivot_longer(mean:sd,
                                           names_to = "parameter",
                                           values_to = "value"))

all_efficiency |>
  ggplot(aes(x = index, y = value, color = model)) +
  geom_point(size = 2) +
  geom_hline(aes(yintercept = 0)) +
  facet_wrap(~parameter, scales = "free_y",
             nrow = 2) +
  theme_bw()

SRJPEdata::weekly_juvenile_abundance_catch_data |>
  filter(site == "lcc", run_year == 2018) |>
  group_by(week, site, run_year) |> # lifestages
  summarise(count = sum(count)) |>
  ggplot(aes(x = week, y = count)) +
  geom_point() +
  theme_bw() +
  labs(title = "count data for lcc 2018")
