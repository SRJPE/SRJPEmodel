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


storage_account = "jpemodelresults"
container_name = "model-results"
model_fits <- readRDS("data/example_american_river_2020.rds")
blob_url <- SRJPEmodel::store_model_fit(con,
                            storage_account = "jpemodelresults",
                            container_name = "model-results",
                            access_key = Sys.getenv("AZ_CONTAINER_ACCESS_KEY"),
                            model_fit_object = model_fits,
                            results_name = "juvenile_abundance",
                            site = "lcc",
                            run_year = 2020,
                            description = "american river 2020 bugs model fit")

model_run <- search_model_run(con, keyword=NULL, model_run_id=NULL, view_all=TRUE)
model_object <- SRJPEmodel::get_model_object(con, keyword = NULL, model_run_id = 10)
model_results <- get_model_results_parameters(con, keyword = NULL, model_run_id=28)
# x <- load_model_fit(con, "missing_mark_recap.bug")
# new_model_fits <- readRDS("data/pCap_model_2025-01-09.rds")
# abundance_fit <- readRDS("data/abundance_model_inputs.rds")
# old_model_fits <- readRDS("data/ubc_2004_2024-07-23.rds")


