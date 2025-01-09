library(rstan)
library(tidyverse)
library(SRJPEmodel)
library(SRJPEdata)
library(R2WinBUGS)

# pCap model
pCap_inputs <- prepare_pCap_inputs(mainstem = FALSE)
pCap <- fit_pCap_model(pCap_inputs$inputs)

# save, fit 01-09-2025
saveRDS(pCap, "C:/Users/Liz/Downloads/pCap_model_2025-01-09.rds")

# read in
# pCap <- readRDS("C:/Users/Liz/Downloads/pCap_model.rds")

abundance_inputs <- prepare_abundance_inputs("ubc", 2018, effort_adjust = T)
lt_pCap_Us <- generate_lt_pCap_Us(abundance_inputs, pCap)

# BUGS --------------------------------------------------------------------
abundance <- fit_abundance_model_BUGS(abundance_inputs, lt_pCap_Us,
                                      "C:/Users/Liz/Documents/SRJPEmodel/model_files/abundance_model.bug",
                                      "C:/Users/Liz/Documents/SRJPEmodel/data-raw/WinBUGS14")

abundance_table <- extract_abundance_estimates("ubc", 2018, abundance_inputs, abundance)

saveRDS(abundance, "ubc_2018_abundance_model.rds")
saveRDS(abundance_inputs, "ubc_2018_abundance_inputs.rds")
saveRDS(abundance_table, "ubc_2018_abundance_fit_table.rds")



# run for all SRJPE sites -------------------------------------------------


fit_all_sites_run_years <- function(site, run_year){

  abundance_inputs <- prepare_abundance_inputs(site, run_year, effort_adjust = T)
  lt_pCap_Us <- generate_lt_pCap_Us(abundance_inputs, pCap)
  abundance <- fit_abundance_model_BUGS(abundance_inputs, lt_pCap_Us,
                                        "C:/Users/Liz/Documents/SRJPEmodel/model_files/abundance_model.bug",
                                        "C:/Users/Liz/Documents/SRJPEmodel/data-raw/WinBUGS14")
  clean_table <- extract_abundance_estimates(site, run_year,
                                             abundance_inputs, abundance)
  return(clean_table)
}

# get trials to fit
site_years_to_exclude <- SRJPEdata::years_to_exclude_rst_data |>
  mutate(site_year = paste(site, year, sep = "_")) |>
  pull(site_year)
trials_to_fit <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
  mutate(site_year = paste(site, year, sep = "_")) |>
  filter(!site_year %in% site_years_to_exclude,
         life_stage != "yearling",
         stream != "sacramento river",
         week %in% c(seq(45, 53), seq(1, 22))) |>
  group_by(year, week, stream, site, run_year) |>
  # keep NAs in count columns
  summarise(count = if(all(is.na(count))) NA_real_ else sum(count, na.rm = TRUE),
            mean_fork_length = mean(mean_fork_length, na.rm = T),
            hours_fished = mean(hours_fished, na.rm = T),
            flow_cfs = mean(flow_cfs, na.rm = T),
            average_stream_hours_fished = mean(average_stream_hours_fished, na.rm = T),
            standardized_flow = mean(standardized_flow, na.rm = T),
            catch_standardized_by_hours_fished = if(all(is.na(catch_standardized_by_hours_fished))) NA_real_ else sum(catch_standardized_by_hours_fished, na.rm = TRUE),
            lgN_prior = mean(lgN_prior, na.rm = T)) |>
  ungroup() |>
  distinct(site, run_year)
  # group_by(run_year, site) |>
  # summarise(count = sum(count, na.rm = T)) |>
  # ungroup() |>
  # filter(!is.na(count)) # filter out those sites where we have no catch data at all (? ubc 2015)

write_csv(trials_to_fit, "C:/Users/Liz/Downloads/bt_spas_x_sites_run_years.csv")

SRJPE_fits_table <- purrr::pmap(list(trials_to_fit$site,
                            trials_to_fit$run_year),
                            fit_all_sites_run_years,
                            .progress = TRUE)

SRJPE_clean_table <- purrr::pmap()
