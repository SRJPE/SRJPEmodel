library(tidyverse)
library(DBI)
library(SRJPEdata)
library(R2WinBUGS)
library(AzureStor)
library(RPostgres)

# fit all the battle creek abundance models and push to cloud
battle_run_years <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
  filter(site == "ubc") |>
  distinct(run_year) |>
  pull(run_year)

pcap_object <- readRDS("C:/Users/Liz/Downloads/pCap_model_2025-03-25.rds")

# set up db connection
cfg <- config::get()

con <- DBI::dbConnect(RPostgres::Postgres(),
                      dbname = cfg$db_name,
                      host = cfg$db_host,
                      port = cfg$db_port,
                      user = cfg$db_user,
                      password = cfg$db_password)

# fit battle creek and push to cloud
for(i in battle_run_years) {
  inputs <- prepare_abundance_inputs("ubc", i, effort_adjust = T)
  lt_pcap_us <- generate_lt_pCap_Us(inputs, pcap_object)
  print(paste("Fitting model for", i))
  fit <- fit_abundance_model_BUGS(inputs, lt_pcap_us,
                                  "C:/Users/Liz/Documents/SRJPEmodel/model_files/abundance_model.bug",
                                  "C:/Users/Liz/Documents/SRJPEmodel/data-raw/WinBUGS14")

  # push to cloud
  print(paste("Uploading to cloud for", i))
  store_model_fit(con,
                  storage_account = "jpemodelresults",
                  container_name = "model-results",
                  access_key = Sys.getenv("AZ_CONTAINER_ACCESS_KEY"),
                  model_fit_object = fit,
                  results_name = "juvenile_abundance",
                  site = "ubc",
                  run_year = i,
                  description = paste("ubc", i, "model fit object from auto-run tests"))

}

