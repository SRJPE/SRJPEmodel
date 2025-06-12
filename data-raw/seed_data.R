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
model_fits <- readRDS("data/example_model_fit_06-12-2025.rds")
model_inputs <- readRDS("data/example_model_inputs_06-11-2025.rds")
results_name <- "p2s"
site <- "battle creek"
description <- "test upload battle creek model"
blob_url <- SRJPEmodel::store_model_fit(con,
                            storage_account = "jpemodelresults",
                            container_name = "model-results",
                            access_key = Sys.getenv("AZ_CONTAINER_ACCESS_KEY"),
                            model_fit_object = model_fits,
                            model_inputs = model_inputs,
                            results_name = results_name,
                            site = site,
                            description = description)

model_run <- SRJPEmodel::search_model_run(con, keyword=NULL, model_run_id=NULL, view_all=TRUE)
model_plot_object <- SRJPEmodel::get_model_object(con, model_component = "model_plot", model_run_id = 319)
model_results <- SRJPEmodel::get_model_results_parameters(con, keyword = NULL, model_run_id=28)
most_recent_model_input <- SRJPEmodel::get_most_recent_model_objects(con, model_component = "model_input")
# x <- load_model_fit(con, "missing_mark_recap.bug")
# new_model_fits <- readRDS("data/pCap_model_2025-01-09.rds")
# abundance_fit <- readRDS("data/abundance_model_inputs.rds")
# old_model_fits <- readRDS("data/ubc_2004_2024-07-23.rds")


