library(tidyverse)
library(SRJPEdata)
library(rstan)
library(R2WinBUGS)
library(scales)
library(DBI) # for stock-recruit
library(RPostgres) # for stock-recruit

# bt-spas-x ---------------------------------------------------------------
# example for upper battle creek (ubc) 2018 (trib)

# run pCap model
pCap_inputs <- prepare_pCap_inputs(mainstem = FALSE, drop_trib_sites = TRUE,
                                   sites_to_drop = "okie dam")
pCap <- fit_pCap_model(pCap_inputs)

# run abundance model
abundance_inputs <- prepare_abundance_inputs("ubc", 2018, effort_adjust = T, pCap)
abundance <- fit_abundance_model_BUGS(abundance_inputs,
                                      # point towards where you store the .bug model
                                      "C:/Users/Liz/Documents/SRJPEmodel/model_files/abundance_model.bug",
                                      # point to where you have WinBUGS
                                      "C:/Users/Liz/Documents/SRJPEmodel/data-raw/WinBUGS14")

abundance_table <- extract_abundance_estimates("ubc", 2018, abundance_inputs, abundance)

plot_juv_data("ubc", 2018)
generate_diagnostic_plot_juv("ubc", 2018, abundance_table)

# example for a mainstem site (tisdale)
tis_inputs <- prepare_pCap_inputs(mainstem = TRUE,
                                  "tisdale")
tis_fit <- fit_pCap_model(tis_inputs$inputs)
tis_abundance_inputs <- prepare_abundance_inputs("tisdale", 2017, effort_adjust = T)
tis_lt_pCap_Us <- generate_lt_pCap_Us(tis_abundance_inputs, tis_fit)

tis_abundance <- fit_abundance_model_BUGS(tis_abundance_inputs, tis_lt_pCap_Us,
                                          # point towards where you store the .bug model
                                          "C:/Users/Liz/Documents/SRJPEmodel/model_files/abundance_model.bug",
                                          # point to where you have WinBUGS
                                          "C:/Users/Liz/Documents/SRJPEmodel/data-raw/WinBUGS14")

tis_abundance_table <- extract_abundance_estimates("tisdale", 2017,
                                                   tis_abundance_inputs, tis_abundance)



# passage to spawner ------------------------------------------------------

P2S_inputs <- prepare_P2S_inputs("battle creek", "wy_type")
P2S_results <- fit_passage_to_spawner_model(P2S_inputs)
P2S_spawners <- extract_P2S_estimates(P2S_results)

P2S_spawners |>
  pivot_wider(names_from = "statistic", values_from = "value") |>
  filter(parameter == "predicted_spawners") |>
  ggplot(aes(x = year, y = `50`)) +
  geom_ribbon(aes(ymin = `2.5`, ymax = `97.5`), alpha = 0.2) +
  geom_line() +
  theme_minimal()


# stock-recruit model -----------------------------------------------------

# get connection
cfg <- config::get()

con <- DBI::dbConnect(RPostgres::Postgres(),
                      dbname = cfg$db_name,
                      host = cfg$db_host,
                      port = cfg$db_port,
                      user = cfg$db_user,
                      password = cfg$db_password)

sr_inputs <- prepare_stock_recruit_inputs(con, stream = "battle creek", adult_data_type = "redd",
                                          covariate = "spawning_above_13_temp_week")
battle_sr <- fit_stock_recruit_model(sr_inputs)
battle_sr_params <- extract_stock_recruit_estimates(sr_inputs, battle_sr) |>
  glimpse()
generate_diagnostic_plot_sr(sr_inputs, battle_sr)
generate_results_plot_sr(sr_inputs, battle_sr)


# survival model ----------------------------------------------------------

# explore results for survival model
survival_inputs <- prepare_survival_inputs(number_of_water_year_types = 3,
                                           effect = "fork_length_effect")
survival_results <- fit_survival_model(survival_inputs)
survival_estimates <- extract_survival_estimates(survival_results)
generate_survival_rate_plot(survival_estimates)


# in season ---------------------------------------------------------------

inseason_inputs <- prepare_inseason_inputs(con, "battle creek", "ubc",
                                             covariate_effect = FALSE,
                                             autocorrelation = FALSE)
inseason_fit <- fit_inseason_model(inseason_inputs)
inseason_estimates <- extract_inseason_estimates(inseason_inputs, inseason_fit)
