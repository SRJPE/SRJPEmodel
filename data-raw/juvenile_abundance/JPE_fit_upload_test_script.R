library(tidyverse)
library(DBI)
library(SRJPEdata)
library(R2WinBUGS)
library(AzureStor)
library(RPostgres)

# fit all the JPE abundance models and push to cloud
JPE_sites <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
  distinct(site,run_year) |>
  arrange(site, run_year) |>
  filter(site != "ubc") # already fit

pcap_object <- readRDS("C:/Users/Liz/Downloads/pCap_model_2025-03-25.rds")

# set up db connection
cfg <- config::get()

con <- DBI::dbConnect(RPostgres::Postgres(),
                      dbname = cfg$db_name,
                      host = cfg$db_host,
                      port = cfg$db_port,
                      user = cfg$db_user,
                      password = cfg$db_password)

# fit sites and push to cloud
# TODO issues with mill creek 2021
JPE_sites_trib <- JPE_sites |>
  filter(!site %in% c("knights landing", "tisdale", "red bluff diversion dam"))

for(i in 1:nrow(JPE_sites_trib)) {

  print(paste("Prepping inputs for", JPE_sites_trib$site[i], JPE_sites_trib$run_year[i]))
  inputs <- prepare_abundance_inputs(JPE_sites_trib$site[i], JPE_sites_trib$run_year[i], effort_adjust = T)
  lt_pcap_us <- generate_lt_pCap_Us(inputs, pcap_object)
  print(paste("Fitting model for", JPE_sites_trib$site[i], JPE_sites_trib$run_year[i]))
  fit <- fit_abundance_model_BUGS(inputs, lt_pcap_us,
                                  "C:/Users/Liz/Documents/SRJPEmodel/model_files/abundance_model.bug",
                                  "C:/Users/Liz/Documents/SRJPEmodel/data-raw/WinBUGS14")

  # push to cloud
  print(paste("Uploading to cloud for", JPE_sites_trib$site[i], JPE_sites_trib$run_year[i]))
  store_model_fit(con,
                  storage_account = "jpemodelresults",
                  container_name = "model-results",
                  access_key = Sys.getenv("AZ_CONTAINER_ACCESS_KEY"),
                  model_fit_object = fit,
                  results_name = "juvenile_abundance",
                  site = JPE_sites_trib$site[i],
                  run_year = JPE_sites_trib$run_year[i],
                  description = paste(JPE_sites_trib$site[i], JPE_sites_trib$run_year[i], "model fit object from auto-run tests"))

}

saveRDS(trib, "C:/Users/Liz/Downloads/jpe_fit_objects_03-31-2025.rds")

# knights landing
kdl_pcap_inputs <- prepare_pCap_inputs(mainstem = TRUE, mainstem_site = "knights landing")
kdl_pCap <- fit_pCap_model(kdl_pcap_inputs$inputs)

store_model_fit(con,
                storage_account = "jpemodelresults",
                container_name = "model-results",
                access_key = Sys.getenv("AZ_CONTAINER_ACCESS_KEY"),
                model_fit_object = kdl_pCap,
                results_name = "pcap_mainstem",
                site = "knights landing",
                description = "knights landing pcap 03-31-2025")

kdl_years <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
  filter(site == "knights landing") |>
  distinct(site, run_year) |>
  pull(run_year)

for(i in kdl_years) {
  print(paste("fitting knights landing", i))
  inputs <- prepare_abundance_inputs("knights landing", i, effort_adjust = T)
  lt_pcap_us <- generate_lt_pCap_Us(inputs, kdl_pCap)
  abundance <- fit_abundance_model_BUGS(inputs, lt_pcap_us,
                                        "C:/Users/Liz/Documents/SRJPEmodel/model_files/abundance_model.bug",
                                        "C:/Users/Liz/Documents/SRJPEmodel/data-raw/WinBUGS14")
  store_model_fit(con,
                  storage_account = "jpemodelresults",
                  container_name = "model-results",
                  access_key = Sys.getenv("AZ_CONTAINER_ACCESS_KEY"),
                  model_fit_object = abundance,
                  results_name = "juvenile_abundance",
                  site = "knights landing",
                  run_year = i,
                  description = paste("knights landing", i))

}

# tisdale
tis_pcap_inputs <- prepare_pCap_inputs(mainstem = TRUE, mainstem_site = "tisdale")
tis_pCap <- fit_pCap_model(tis_pcap_inputs$inputs)

store_model_fit(con,
                storage_account = "jpemodelresults",
                container_name = "model-results",
                access_key = Sys.getenv("AZ_CONTAINER_ACCESS_KEY"),
                model_fit_object = tis_pCap,
                results_name = "pcap_mainstem",
                site = "tisdale",
                description = "tisdale pcap 03-31-2025")

tis_years <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
  filter(site == "tisdale") |>
  distinct(site, run_year) |>
  pull(run_year)

for(i in tis_years) {
  print(paste("fitting tisdale", i))
  inputs <- prepare_abundance_inputs("tisdale", i, effort_adjust = T)
  lt_pcap_us <- generate_lt_pCap_Us(inputs, tis_pCap)
  abundance <- fit_abundance_model_BUGS(inputs, lt_pcap_us,
                                        "C:/Users/Liz/Documents/SRJPEmodel/model_files/abundance_model.bug",
                                        "C:/Users/Liz/Documents/SRJPEmodel/data-raw/WinBUGS14")

  store_model_fit(con,
                  storage_account = "jpemodelresults",
                  container_name = "model-results",
                  access_key = Sys.getenv("AZ_CONTAINER_ACCESS_KEY"),
                  model_fit_object = abundance,
                  results_name = "juvenile_abundance",
                  site = "tisdale",
                  run_year = i,
                  description = paste("tisdale", i))

}


