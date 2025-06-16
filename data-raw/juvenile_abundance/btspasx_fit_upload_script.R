library(tidyverse)
library(R2WinBUGS)

# check you are running on a PC
operating_system <- ifelse(grepl("Mac", Sys.info()['nodename']) | grepl("MBP", Sys.info()['nodename']), "mac", "pc")

if(operating_system != "pc") {
  stop("You must be using a PC to run winBUGS")
}
# set up db connection
cfg <- config::get()

con <- DBI::dbConnect(RPostgres::Postgres(),
                      dbname = cfg$db_name,
                      host = cfg$db_host,
                      port = cfg$db_port,
                      user = cfg$db_user,
                      password = cfg$db_password)


# fit pcap model ----------------------------------------------------------

# trib
pcap_inputs <- prepare_pCap_inputs(mainstem = FALSE)
model_fit_pcap <- fit_pCap_model(pcap_inputs)
store_model_fit(con,
                model_fit_object = model_fit_pcap,
                model_inputs = pcap_inputs,
                results_name = "pcap_all",
                description = paste("pcap model fit", Sys.Date()))



# abundance models --------------------------------------------------------

# trib
JPE_sites_trib <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
  dplyr::distinct(site, run_year) |>
  filter(!site %in% c("knights landing", "tisdale", "red bluff diversion dam"))


for(i in 1:nrow(JPE_sites_trib)) {

  print(paste("Prepping inputs for", JPE_sites_trib$site[i], JPE_sites_trib$run_year[i]))
  inputs <- prepare_abundance_inputs(JPE_sites_trib$site[i], JPE_sites_trib$run_year[i],
                                     effort_adjust = T, model_fit_pcap)

  print(paste("Fitting model for", JPE_sites_trib$site[i], JPE_sites_trib$run_year[i]))
  fit <- fit_abundance_model_BUGS(inputs,
                                  "C:/Users/Liz/Documents/SRJPEmodel/model_files/abundance_model.bug",
                                  "C:/Users/Liz/Documents/SRJPEmodel/data-raw/WinBUGS14")

  # push to cloud
  print(paste("Uploading to cloud for", JPE_sites_trib$site[i], JPE_sites_trib$run_year[i]))
  store_model_fit(con,
                  model_fit_object = fit,
                  model_inputs = inputs,
                  results_name = inputs$model_name,
                  description = paste(JPE_sites_trib$site[i], JPE_sites_trib$run_year[i], "model fit object from auto-run tests", Sys.Date()))

}



# knights landing
kdl_pcap_inputs <- prepare_pCap_inputs(mainstem = TRUE, mainstem_site = "knights landing")
kdl_pCap <- fit_pCap_model(kdl_pcap_inputs)

store_model_fit(con,
                model_fit_object = kdl_pCap,
                model_inputs = kdl_pcap_inputs,
                results_name = "pcap_mainstem",
                description = paste("knights landing pcap", Sys.Date()))

kdl_years <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
  filter(site == "knights landing") |>
  dplyr::distinct(site, run_year) |>
  pull(run_year)

for(i in kdl_years) {
  print(paste("fitting knights landing", i))
  inputs <- prepare_abundance_inputs("knights landing", i,
                                     effort_adjust = T, kdl_pCap)
  abundance <- fit_abundance_model_BUGS(inputs,
                                        "C:/Users/Liz/Documents/SRJPEmodel/model_files/abundance_model.bug",
                                        "C:/Users/Liz/Documents/SRJPEmodel/data-raw/WinBUGS14")
  store_model_fit(con,
                  model_fit_object = abundance,
                  model_inputs = inputs,
                  results_name = inputs$model_name,
                  description = paste("knights landing", i, "fit", Sys.Date()))

}

# tisdale
tis_pcap_inputs <- prepare_pCap_inputs(mainstem = TRUE, mainstem_site = "tisdale")
tis_pCap <- fit_pCap_model(tis_pcap_inputs)

store_model_fit(con,
                model_fit_object = tis_pCap,
                model_inputs = tis_pcap_inputs,
                results_name = "pcap_mainstem",
                description = paste("tisdale pcap", Sys.Date()))

tis_years <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
  filter(site == "tisdale") |>
  dplyr::distinct(site, run_year) |>
  pull(run_year)

for(i in tis_years) {
  print(paste("fitting tisdale", i))
  inputs <- prepare_abundance_inputs("tisdale", i,
                                     effort_adjust = T, tis_pCap)
  abundance <- fit_abundance_model_BUGS(inputs,
                                        "C:/Users/Liz/Documents/SRJPEmodel/model_files/abundance_model.bug",
                                        "C:/Users/Liz/Documents/SRJPEmodel/data-raw/WinBUGS14")
  store_model_fit(con,
                  model_fit_object = abundance,
                  model_inputs = inputs,
                  results_name = inputs$model_name,
                  description = paste("tisdale", i, "fit", Sys.Date()))

}


