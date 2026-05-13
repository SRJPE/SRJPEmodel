
library(tidyverse)
library(SRJPEdata)
library(SRJPEmodel)
library(rstan)
library(R2WinBUGS)

# stan options
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# file path for bugs
bugs_directory  <- "C:/Users/Liz/Documents/SRJPEmodel/data-raw/WinBUGS14"

# connect to db
cfg <- config::get()
con <- DBI::dbConnect(
  RPostgres::Postgres(),
  dbname   = cfg$db_name,
  host     = cfg$db_host,
  port     = cfg$db_port,
  user     = cfg$db_user,
  password = cfg$db_password
)
on.exit(DBI::dbDisconnect(con), add = TRUE)


# fit the 3 pCap models ---------------------------------------------------

for(ii in 1:3){#1:4
  if(ii==1){
    description <- paste0("pCap all sites test ", Sys.Date())
    pCap_inputs <- prepare_pCap_inputs(model_type = "all_sites")

  } else if(ii==2) {
    description <- paste0("pCap one site skew test ", Sys.Date())
    pCap_inputs <- prepare_pCap_inputs(model_type = "one_site",
                                       skew = T,
                                       site_selection = "knights landing")

  } else if(ii==3) {
    description <- paste0("pCap one site skew test ", Sys.Date())
    pCap_inputs <- prepare_pCap_inputs(model_type = "one_site",
                                       skew = T,
                                       site_selection = "tisdale")
  }

  pCap_fit <- fit_pCap_model(pCap_inputs) #automates some of arguments, but code below better for testing

  store_model_fit(
    con = con,
    model_fit_object = pCap_fit,
    model_inputs = pCap_inputs,
    description = description
  )

}


# fit all abundance models ------------------------------------------------

combinations_to_fit <- SRJPEdata::years_to_include_rst_data |>
  filter(!site %in% c("lbc", "red bluff diversion dam")) |>
  arrange(site)

pcap_tisdale <- get_model_fit("pcap_one_site",
                                          con = con,
                                          site_selection = "tisdale")
pcap_knights_landing <- get_model_fit("pcap_one_site",
                                          con = con,
                                          site_selection = "knights landing")
pcap_all_sites <- get_model_fit("pcap_all_sites",
                                          con = con)

for(i in nrow(combinations_to_fit)) {
  site_to_fit <- combinations_to_fit$site[i]
  run_year_to_fit <- combinations_to_fit$run_year[i]
  message(site_to_fit)
  message(run_year_to_fit)

  # load pcap
  if(site_to_fit == "knights landing") {
    abundance_inputs <- prepare_abundance_inputs(
      site = site_to_fit,
      run_year = run_year_to_fit,
      pCap_model_type = "one_site",
      pCap_model_object = pcap_knights_landing,
      min_pCap_mult = 0.5
    )

  } else if(site_to_fit == "tisdale") {
    abundance_inputs <- prepare_abundance_inputs(
      site = site_to_fit,
      run_year = run_year_to_fit,
      pCap_model_type = "one_site",
      pCap_model_object = pcap_tisdale,
      min_pCap_mult = 0.5
    )
  } else {
    abundance_inputs <- prepare_abundance_inputs(
      site = site_to_fit,
      run_year = run_year_to_fit,
      pCap_model_type = "all_sites",
      pCap_model_object = pcap_all_sites,
      min_pCap_mult = 0.5
    )
  }

  # fit abundance
  abundance_fit <- fit_abundance_model_BUGS(
    abundance_inputs,
    bugs_directory
  )

  store_model_fit(
    con              = con,
    model_fit_object = abundance_fit,
    model_inputs     = abundance_inputs,
    description      = paste("TEST —", site_to_fit, run_year_to_fit, "abundance", Sys.Date())
  )
}


