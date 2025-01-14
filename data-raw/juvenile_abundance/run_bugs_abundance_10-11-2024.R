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
saveRDS(pCap, "C:/Users/Liz/Downloads/pCap_model_2025-01-09.rds")

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

sites_to_run <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
  distinct(site, run_year) |>
  arrange(site, run_year)

all_JPE_sites_clean <- run_bt_spas_x_JPE_sites(sites_to_run = sites_to_run, run_pCap = FALSE,
                                               pCap_model_object_filepath = "C:/Users/Liz/Downloads/pCap_model_2025-01-09.rds",
                                               bugs_model_file = "C:/Users/Liz/Documents/SRJPEmodel/model_files/abundance_model.bug",
                                               bugs_directory = "C:/Users/Liz/Documents/SRJPEmodel/data-raw/WinBUGS14")

write_csv(all_JPE_sites_clean, "C:/Users/Liz/Downloads/all_jpe_sites_fit.csv")
saveRDS(all_JPE_sites_clean, "C:/Users/Liz/Downloads/all_JPE_sites_clean.rds")

