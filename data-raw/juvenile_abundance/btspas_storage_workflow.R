library(tidyverse)
library(SRJPEmodel)
library(SRJPEdata)

# this is what the workflow looks like to produce the objects

# fit the pCap model
pCap_inputs <- prepare_pCap_inputs(mainstem = FALSE)
pCap <- fit_pCap_model(pCap_inputs$inputs) # pCap is a large STANfit object that needs to be stored

# get inputs for every site and run year combination
sites_to_run <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
  distinct(site, run_year) |>
  filter(!site %in% c("knights landing", "red bluff diversion dam",
                      "tisdale")) |>
  arrange(site, run_year)

sites_to_run <- sites_to_run |>
  mutate(filter_out = ifelse(site == "mill creek" & run_year == 2025, TRUE, FALSE)) |>
  filter(!filter_out) |>
  select(-filter_out)

all_abundance_inputs <- purrr::pmap(list(sites_to_run$site,
                                         sites_to_run$run_year,
                                         T), # effort adjust = T
                                    prepare_abundance_inputs,
                                    .progress = TRUE)

# these are the data objects to store
# you can find them here: https://netorg629193-my.sharepoint.com/personal/avizek_flowwest_com/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Favizek%5Fflowwest%5Fcom%2FDocuments%2Fsrjpemodel%5Foutput%5Ffiles&ga=1
object.size(pCap) # example file in sharepoint is pCap_model_2024-01-09.rds
object.size(all_abundance_inputs) # example file in sharepoint is abundance_model_inputs.rds
