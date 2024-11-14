# test out mokelumne data
library(tidyverse)

weekly_juvenile_abundance_catch_data <- read_csv("~/Downloads/weekly_juvenile_abundance_catch_data.csv") # |>
  #mutate(lgN_prior = ifelse(is.na(lgN_prior), log(count/1000 + 1)/0.025, lgN_prior))

priors <- weekly_juvenile_abundance_catch_data |>
  filter(!is.na(lgN_prior)) |>
  group_by(site) |>
  summarise(avg = mean(lgN_prior)) |>
  ungroup() |>
  glimpse()

weekly_juvenile_abundance_catch_data <- weekly_juvenile_abundance_catch_data |>
  mutate(lgN_prior = ifelse(is.na(lgN_prior), mean(priors$avg)/2, lgN_prior))

weekly_juvenile_abundance_efficiency_data <- read_csv("~/Downloads/weekly_juvenile_abundance_efficiency_data.csv")

american <- run_single_bt_spas_x_stan(SRJPEmodel::bt_spas_x_bayes_params,
                                       weekly_juvenile_abundance_catch_data,
                                       weekly_juvenile_abundance_efficiency_data,
                                       "golf",
                                       run_year = 2005,
                                       lifestage = NA,
                                       effort_adjust = F)

american <- run_single_bt_spas_x_stan(SRJPEmodel::bt_spas_x_bayes_params,
                                      weekly_juvenile_abundance_catch_data,
                                      weekly_juvenile_abundance_efficiency_data,
                                      "ubc",
                                      run_year = 2009,
                                      lifestage = "fry",
                                      effort_adjust = F)
