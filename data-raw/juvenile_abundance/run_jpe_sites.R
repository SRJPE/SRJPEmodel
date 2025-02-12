library(rstan)
library(tidyverse)
library(SRJPEmodel)
library(SRJPEdata)
library(R2WinBUGS)


# run for one site/run year -----------------------------------------------

# pCap model
pCap_inputs <- prepare_pCap_inputs(mainstem = FALSE)
pCap <- fit_pCap_model(pCap_inputs$inputs)

# save, fit 01-09-2025
saveRDS(pCap, "C:/Users/Liz/Downloads/pCap_model_2025-02-10.rds")

pCap <- readRDS("C:/Users/Liz/Downloads/pCap_model_2025-01-09.rds")

abundance_inputs <- prepare_abundance_inputs("ubc", 2018, effort_adjust = T)
lt_pCap_Us <- generate_lt_pCap_Us(abundance_inputs, pCap)

abundance <- fit_abundance_model_BUGS(abundance_inputs, lt_pCap_Us,
                                      "C:/Users/Liz/Documents/SRJPEmodel/model_files/abundance_model.bug",
                                      "C:/Users/Liz/Documents/SRJPEmodel/data-raw/WinBUGS14")

abundance_table <- extract_abundance_estimates("ubc", 2018, abundance_inputs, abundance)

saveRDS(abundance, "ubc_2018_abundance_model.rds")
saveRDS(abundance_inputs, "ubc_2018_abundance_inputs.rds")
saveRDS(abundance_table, "ubc_2018_abundance_fit_table.rds")


# run for all SRJPE sites -------------------------------------------------

# tribs

sites_to_run <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
  distinct(site, run_year) |>
  filter(!site %in% c("knights landing", "red bluff diversion dam",
                      "tisdale")) |>
  arrange(site, run_year)

all_JPE_sites_clean <- run_bt_spas_x_JPE_sites(sites_to_run = sites_to_run, run_pCap = FALSE,
                                               mainstem = FALSE,
                                               pCap_model_object_filepath = "C:/Users/Liz/Downloads/pCap_model_2025-01-09.rds",
                                               bugs_model_file = "C:/Users/Liz/Documents/SRJPEmodel/model_files/abundance_model.bug",
                                               bugs_directory = "C:/Users/Liz/Documents/SRJPEmodel/data-raw/WinBUGS14")

write_csv(all_JPE_sites_clean, "C:/Users/Liz/Downloads/all_jpe_sites_fit.csv")
saveRDS(all_JPE_sites_clean, "C:/Users/Liz/Downloads/all_JPE_sites_clean.rds")

# mainstem

mainstem_sites <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
  distinct(site, run_year) |>
  filter(site %in% c("knights landing", "red bluff diversion dam",
                      "tisdale")) |>
  arrange(site, run_year)

# tisdale
tisdale <- run_bt_spas_x_JPE_sites(sites_to_run = mainstem_sites |>
                                             filter(site == "tisdale"),
                                           run_pCap = FALSE,
                                           pCap_model_object_filepath = "C:/Users/Liz/Downloads/pCap_tisdale_fit.rds",
                                           bugs_model_file = "C:/Users/Liz/Documents/SRJPEmodel/model_files/abundance_model.bug",
                                           bugs_directory = "C:/Users/Liz/Documents/SRJPEmodel/data-raw/WinBUGS14")

write_csv(tisdale, "C:/Users/Liz/Downloads/tisdale_fits.csv")
saveRDS(tisdale, "C:/Users/Liz/Downloads/tisdale_fits.rds")



# knights landing
knights_landing <- run_bt_spas_x_JPE_sites(sites_to_run = mainstem_sites |>
                                             filter(site == "knights landing"),
                                           run_pCap = FALSE,
                                           pCap_model_object_filepath = "C:/Users/Liz/Downloads/pCap_knights_landing_fit.rds",
                                           bugs_model_file = "C:/Users/Liz/Documents/SRJPEmodel/model_files/abundance_model.bug",
                                           bugs_directory = "C:/Users/Liz/Documents/SRJPEmodel/data-raw/WinBUGS14")
write_csv(knights_landing, "C:/Users/Liz/Downloads/knights_landing_fits.csv")
saveRDS(knights_landing, "C:/Users/Liz/Downloads/knights_landing_fits.rds")

# red bluff diversion dam
rbdd <- run_bt_spas_x_JPE_sites(sites_to_run = mainstem_sites |>
                                             filter(site == "red bluff diversion dam"),
                                           run_pCap = FALSE,
                                           pCap_model_object_filepath = "C:/Users/Liz/Downloads/pCap_rbdd_fit.rds",
                                           bugs_model_file = "C:/Users/Liz/Documents/SRJPEmodel/model_files/abundance_model.bug",
                                           bugs_directory = "C:/Users/Liz/Documents/SRJPEmodel/data-raw/WinBUGS14")
write_csv(rbdd, "C:/Users/Liz/Downloads/rbdd_fits.csv")
saveRDS(rbdd, "C:/Users/Liz/Downloads/rbdd_fits.rds")


# generate data for Noble -------------------------------------------------

# data to run PLAD that applies filter
# sent 1-14-2025

data_for_PLAD <- SRJPEdata::rst_catch |>
  mutate(run_year = ifelse(julian_week >= 45, julian_year + 1, julian_year)) |>
  left_join(SRJPEdata::chosen_site_years_to_model |> # need to make sure to filter out years that have been excluded
              select(run_year = monitoring_year, stream, site) |>
              mutate(include = T)) |>
  filter(include == T) |>
  select(-include) |>
  mutate(exclude = ifelse(site == "mill creek" & run_year == 2025, TRUE, FALSE)) |>
  filter(!exclude) |>
  select(-exclude)

write_csv(data_for_PLAD, "~/Downloads/rst_catch_filtered_for_PLAD.csv")

