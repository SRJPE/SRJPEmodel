

# run lcc, feather, ubc ---------------------------------------------------

# function to run the split models

run_split_models <- function(site_arg, run_year_arg,
                             prior_scalar) {

  input <- run_single_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
                                SRJPEdata::weekly_juvenile_abundance_catch_data |>
                                  filter(life_stage %in% c("fry", "smolt"),
                                         hours_fished > 0),
                                SRJPEdata::weekly_juvenile_abundance_efficiency_data,
                                site_arg, run_year_arg, NA, T,
                                here::here("data-raw", "WinBUGS14"),
                                no_cut =F)
  if(site_arg == "lcc") {
    input$data$lgN_max <- input$data$lgN_max + prior_scalar
  }

  inits <- input$inits[[1]][-8]
  inits <- inits[-8]
  inits <- list(inits, inits, inits)

  model_name <- paste0("pCap_", input$model_name, ".stan")
  print(model_name)

  options(mc.cores=parallel::detectCores())
  cli::cli_alert("running pCap model")
  pcap <- rstan::stan(model_code = read_file(here::here("model_files",
                                                        model_name)),
                      data = input$data,
                      init = inits,
                      chains = SRJPEmodel::bt_spas_x_bayes_params$number_chains,
                      iter = SRJPEmodel::bt_spas_x_bayes_params$number_mcmc,
                      seed = 84735)

  # catch any rhat values
  pcap_estimates <- rstan::summary(pcap, pars = c("lt_pCap_U"))$summary |>
    data.frame()

  pcap_estimates |>
    filter(Rhat > 1.05) |>
    row.names()

  # pars - need mean and sd
  logit_pCaps <- rstan::summary(pcap, pars = c("lt_pCap_U"))$summary |>
    data.frame()

  input$data$lt_pCap_mu <- logit_pCaps$mean
  input$data$lt_pCap_sd <- logit_pCaps$sd

  inits <- list("b_sp" = input$inits[[1]]$b_sp,
                "lg_N" = input$inits[[1]]$lg_N,
                "lt_pCap_U" = rnorm(input$data$Nstrata_wc,
                                    mean(logit_pCaps$mean),
                                    mean(logit_pCaps$sd)))
  inits <- list(inits, inits, inits)

  cli::cli_alert("running abundance model")

  abundance <- rstan::stan(model_code = read_file(here::here("model_files",
                                                             "abundance_model.stan")),
                           data = input$data,
                           init = inits,
                           chains = SRJPEmodel::bt_spas_x_bayes_params$number_chains,
                           iter = SRJPEmodel::bt_spas_x_bayes_params$number_mcmc * 3,
                           seed = 84735)

  abundance_estimates <- rstan::summary(abundance, pars = c("N"))$summary |>
    data.frame()

  abundance_estimates |>
    filter(Rhat > 1.05) |>
    row.names()

  return(list("pCap" = pcap,
              "abundance" = abundance,
              "estimates" = list("pcap_estimates" = pcap_estimates,
                                 "abundance_estimates" = abundance_estimates)
              )
  )

}

plot_raw_data <- function(site_arg, run_year_arg) {
  julian_week_to_date_lookup <- read.table(file = "data-raw/juvenile_abundance/archive/btspas_model_code/Jwk_Dates.txt", header = F) |>
    tibble() |>
    filter(V1 != "Jwk") |>
    mutate(V1 = as.numeric(V1)) |>
    select(Jwk = V1, date = V2)

  data <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
    filter(run_year == run_year_arg,
           site == site_arg) |> #,
    #week %in% c(seq(45, 53), seq(1, 22))) |>
    group_by(year, week, stream, site, run_year) |>
    summarise(count = sum(count, na.rm = T),
              mean_fork_length = mean(mean_fork_length, na.rm = T),
              hours_fished = mean(hours_fished, na.rm = T),
              flow_cfs = mean(flow_cfs, na.rm = T),
              average_stream_hours_fished = mean(average_stream_hours_fished, na.rm = T),
              standardized_flow = mean(standardized_flow, na.rm = T),
              catch_standardized_by_hours_fished = sum(catch_standardized_by_hours_fished, na.rm = T),
              lgN_prior = mean(lgN_prior, na.rm = T)) |>
    ungroup() |>
    left_join(SRJPEdata::weekly_juvenile_abundance_efficiency_data,
              by = c("year", "run_year", "week", "stream", "site")) |>
    mutate(count = round(count, 0),
           catch_standardized_by_hours_fished = round(catch_standardized_by_hours_fished, 0),
           # plot things
           lincoln_peterson_abundance = count * (number_released / number_recaptured),
           lincoln_peterson_efficiency = number_recaptured / number_released) |>
    left_join(julian_week_to_date_lookup, by = c("week" = "Jwk")) |>
    mutate(year = ifelse(week >= 43, run_year - 1, run_year),
           fake_date = ymd(paste0(year, "-01-01")),
           final_date = fake_date + weeks(week - 1),
           date = format(final_date, "%b-%d"))

  data |>
    ggplot(aes(x = final_date, y = count)) +
    geom_point() +
    #scale_x_date(date_breaks = "1 week", date_labels = "%b %d") +
    scale_x_date(date_breaks = "1 week", date_labels = "%V") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(x = "Date",
         y = "Raw catch",
         title = paste0("Catch at ", site_arg, " for run year ",
                        run_year_arg))
}

# diagnostic plot function
diagnostic_plots_split <- function(site_arg, run_year_arg,
                                   model_object) {


  julian_week_to_date_lookup <- read.table(file = "data-raw/juvenile_abundance/archive/btspas_model_code/Jwk_Dates.txt", header = F) |>
    tibble() |>
    filter(V1 != "Jwk") |>
    mutate(V1 = as.numeric(V1)) |>
    select(Jwk = V1, date = V2)

  data <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
    filter(run_year == run_year_arg,
           site == site_arg) |> #,
    #week %in% c(seq(45, 53), seq(1, 22))) |>
    group_by(year, week, stream, site, run_year) |>
    summarise(count = sum(count, na.rm = T),
              mean_fork_length = mean(mean_fork_length, na.rm = T),
              hours_fished = mean(hours_fished, na.rm = T),
              flow_cfs = mean(flow_cfs, na.rm = T),
              average_stream_hours_fished = mean(average_stream_hours_fished, na.rm = T),
              standardized_flow = mean(standardized_flow, na.rm = T),
              catch_standardized_by_hours_fished = sum(catch_standardized_by_hours_fished, na.rm = T),
              lgN_prior = mean(lgN_prior, na.rm = T)) |>
    ungroup() |>
    left_join(SRJPEdata::weekly_juvenile_abundance_efficiency_data,
              by = c("year", "run_year", "week", "stream", "site")) |>
    mutate(count = round(count, 0),
           catch_standardized_by_hours_fished = round(catch_standardized_by_hours_fished, 0),
           # plot things
           lincoln_peterson_abundance = count * (number_released / number_recaptured),
           lincoln_peterson_efficiency = number_recaptured / number_released) |>
    left_join(julian_week_to_date_lookup, by = c("week" = "Jwk")) |>
    mutate(year = ifelse(week >= 43, run_year - 1, run_year),
           fake_date = ymd(paste0(year, "-01-01")),
           final_date = fake_date + weeks(week - 1),
           date = format(final_date, "%b-%d"),
           week_index = row_number())

  pCap_estimates <- rstan::summary(model_object$abundance,pars=c("lt_pCap_U"))$summary |>
    data.frame() |>
    tibble::rownames_to_column("parameter") |>
    mutate(week_index = readr::parse_number(parameter)) |>
    select(week_index, mean, sd, `50` = X50., `2.5` = X2.5., `97.5` = X97.5.) |>
    mutate(parameter = "lt_pCap_U")

  N_estimates <- rstan::summary(model_object$abundance,pars=c("N"))$summary |>
    data.frame() |>
    tibble::rownames_to_column("parameter") |>
    mutate(week_index = readr::parse_number(parameter)) |>
    select(week_index, mean, sd, `50` = X50., `2.5` = X2.5., `97.5` = X97.5.) |>
    mutate(parameter = "N")


  data |>
    ggplot(aes(x = final_date, y = catch_standardized_by_hours_fished)) +
    geom_bar(stat = "identity", fill = "grey", width = 5) +
    scale_x_date(date_breaks = "1 week", date_labels = "%b %d") +
    #scale_x_date(date_breaks = "1 week", date_labels = "%V") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(x = "Date",
         y = "Raw catch",
         title = paste0("Catch at ", site_arg, " for run year ",
                        run_year_arg))

  # abundance plot
  abundance_plot <- data |>
    left_join(N_estimates) |>
    #mutate(across(c(mean, `50`, `2.5`, `97.5`), plogis)) |>
    ggplot(aes(x = final_date, y = `50`)) +
    geom_bar(stat = "identity", fill = "grey", width = 5) +
    geom_errorbar(aes(x = final_date, ymin = `2.5`, ymax = `97.5`), width = 0.2) +
    geom_point(aes(x = final_date, y = lincoln_peterson_abundance),
               shape = 1, color = "blue") +
    # geom_point(aes(x = fake_date, y = Inf, color = sampled),
    #            size = 3) +
    geom_text(aes(x = final_date, y = Inf,
                  label = paste(count),
                  angle = 90),
              hjust = 1,
              size = 3) +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "red")) +
    theme_minimal() +
    labs(x = "",
         #x = "Date",
         y = "Abundance",
         title = paste(site_arg, run_year_arg)) +
    theme(axis.text.x=element_blank())
    # scale_x_date(date_breaks = "1 week", date_labels = "%b %d") +
    # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  # efficiency
  efficiency_plot <- data |>
    left_join(pCap_estimates) |>
    mutate(across(c(mean, `50`, `2.5`, `97.5`), plogis),
           number_released_label = ifelse(is.na(number_released), "", number_released),
           number_recaptured_label = ifelse(is.na(number_recaptured), "", number_recaptured)) |>
    ggplot(aes(x = final_date, y = `50`)) +
    geom_bar(stat = "identity", fill = "grey", width = 4) +
    geom_errorbar(aes(x = final_date, ymin = `2.5`, ymax = `97.5`), width = 0.2) +
    geom_point(aes(x = final_date, y = lincoln_peterson_efficiency),
               shape = 1, color = "blue") +
    geom_text(aes(x = final_date, y = Inf,
                  label = paste(number_released_label, number_recaptured_label),
                  angle = 90),
              hjust = 1,
              size = 3) +
    # geom_point(aes(x = fake_date, y = Inf, color = sampled),
    #            size = 3) +
    theme_minimal() +
    labs(x = "Date", y = "Weekly Efficiency") +
    scale_x_date(date_breaks = "1 week", date_labels = "%b %d") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  # arrange together
  gridExtra::grid.arrange(abundance_plot, efficiency_plot)
}

# run on clear creek, no yearling
weekly_catch <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
  filter(life_stage %in% c("fry", "smolt"))

weekly_efficiency <- SRJPEdata::weekly_juvenile_abundance_efficiency_data

# only fit for site/run year combos with data
trials_to_fit <- weekly_catch |>
  mutate(filter_out = ifelse(is.na(life_stage) & count > 0, TRUE, FALSE)) |> # we do not want to keep NA lifestage associated with counts > 0
  filter(!filter_out,
         site == "lcc",
         week %in% c(seq(45, 53), seq(1, 22))) |>
  mutate(count = round(count, 0),
         catch_standardized_by_hours_fished = round(catch_standardized_by_hours_fished, 0)) |>
  group_by(site, run_year) |>
  tally() |>
  arrange(desc(run_year))


# run ---------------------------------------------------------------------

# 2024
plot_raw_data("lcc", 2024)
lcc_2024 <- run_split_models("lcc", 2024, 1) # update scalar
lcc_2024$estimates
diagnostic_plots_split("lcc", 2024, lcc_2024)

# 2023
plot_raw_data("lcc", 2023)
lcc_2023 <- run_split_models("lcc", 2023, 1) # update scalar
lcc_2023$estimates
diagnostic_plots_split("lcc", 2023, lcc_2023)

# 2022
# no efficiency data - not converged
plot_raw_data("lcc", 2022)
lcc_2022 <- run_split_models("lcc", 2022, 12) # update scalar
lcc_2022$estimates
diagnostic_plots_split("lcc", 2022, lcc_2022)

# 2021
plot_raw_data("lcc", 2021)
lcc_2021 <- run_split_models("lcc", 2021, 15) # update scalar
lcc_2021$estimates
diagnostic_plots_split("lcc", 2021, lcc_2021)

# 2020
plot_raw_data("lcc", 2020)
lcc_2020 <- run_split_models("lcc", 2020, 10) # update scalar
lcc_2020$estimates |>
  bind_rows() |>
  filter(Rhat > 1.05)
diagnostic_plots_split("lcc", 2020, lcc_2020)

# 2019
plot_raw_data("lcc", 2019)
lcc_2019 <- run_split_models("lcc", 2019, 9) # update scalar
lcc_2019$estimates |>
  bind_rows() |>
  filter(Rhat > 1.05)
diagnostic_plots_split("lcc", 2019, lcc_2019)

# 2018
plot_raw_data("lcc", 2018)
lcc_2018 <- run_split_models("lcc", 2018, 10) # update scalar
lcc_2018$estimates |>
  bind_rows() |>
  filter(Rhat > 1.05)
diagnostic_plots_split("lcc", 2018, lcc_2018)

# 2017
plot_raw_data("lcc", 2017)
lcc_2017 <- run_split_models("lcc", 2017, 10) # update scalar
lcc_2017$estimates |>
  bind_rows() |>
  filter(Rhat > 1.05)
diagnostic_plots_split("lcc", 2017, lcc_2017)

# 2016
plot_raw_data("lcc", 2016)
lcc_2016 <- run_split_models("lcc", 2016, 12) # update scalar
lcc_2016$estimates |>
  bind_rows() |>
  filter(Rhat > 1.05)
diagnostic_plots_split("lcc", 2016, lcc_2016)

# 2015
plot_raw_data("lcc", 2015)
lcc_2015 <- run_split_models("lcc", 2015, 11) # update scalar
lcc_2015$estimates |>
  bind_rows() |>
  filter(Rhat > 1.05)
diagnostic_plots_split("lcc", 2015, lcc_2015)

# 2014
plot_raw_data("lcc", 2014)
lcc_2014 <- run_split_models("lcc", 2014, 13) # update scalar
lcc_2014$estimates |>
  bind_rows() |>
  filter(Rhat > 1.05)
diagnostic_plots_split("lcc", 2014, lcc_2014)

# 2013
plot_raw_data("lcc", 2013)
lcc_2013 <- run_split_models("lcc", 2013, 12) # update scalar
lcc_2013$estimates |>
  bind_rows() |>
  filter(Rhat > 1.05)
diagnostic_plots_split("lcc", 2013, lcc_2013)

# 2011
plot_raw_data("lcc", 2011)
lcc_2011 <- run_split_models("lcc", 2011, 12) # update scalar
lcc_2011$estimates |>
  bind_rows() |>
  filter(Rhat > 1.05)

saveRDS(lcc_2011, file = here::here("data-raw", "juvenile_abundance",
                                    "split_lcc_results", "lcc_2011.rds"))


cc_results <- purrr::pmap(list(trials_to_fit$site[11:26],
                               trials_to_fit$run_year[11:26],
                               10),
                          run_split_models,
                          .progress = TRUE)

saveRDS(bugs_results, here::here("data-raw", "juvenile_abundance",
                                 "clear_results_BUGS_nov_2024.rds"))



# older analyses ----------------------------------------------------------


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
