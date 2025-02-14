library(readr)
library(DBI)
library(dplyr)
library(tidyverse)
cfg <- config::get(config = "default")

con <- DBI::dbConnect(RPostgres::Postgres(),
                      dbname = cfg$db_name,
                      host = cfg$db_host,
                      port = cfg$db_post,
                      user = cfg$db_user,
                      password = cfg$db_password)


storage_account = "geneticsedidata"
container_name = "model-fits"
model_fits <- readRDS("data/pCap_model_2025-01-09.rds")
store_model_fit(con,
                storage_account = storage_account,
                container_name = container_name,
                # access_key = access_key,
                data = model_fits,
                results_name = "new_model_pcap_test_v6",
                description = "test upload v6",
                mainstem = FALSE)

model_run <- search_model_run(con, keyword=NULL, model_run_id=28, view_all=TRUE)
model_object <- get_model_object(con, keyword = NULL, model_run_id = 28)
model_results <- get_model_results_parameters(con, keyword = NULL, model_run_id=28)
# x <- load_model_fit(con, "missing_mark_recap.bug")
# new_model_fits <- readRDS("data/pCap_model_2025-01-09.rds")
# abundance_fit <- readRDS("data/abundance_model_inputs.rds")
# old_model_fits <- readRDS("data/ubc_2004_2024-07-23.rds")


