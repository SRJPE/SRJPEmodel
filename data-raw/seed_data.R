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


join_lookup <- function(df, db_table, model_lookup_column, db_lookup_column, final_column){
  db_lookup_column_sym <- sym(db_lookup_column)
  model_lookup_column_sym <- sym(model_lookup_column)
  final_column_sym <- sym(final_column)

  all_lookup_val <- unique(df[[model_lookup_column]])
  all_val_id <- tbl(con, db_table) |>
    filter(!!db_lookup_column_sym %in% all_lookup_val) |>
    select(id, !!db_lookup_column_sym) |>
    as_tibble()
  df_with_id <- merge(df, all_val_id, by.x = model_lookup_column, by.y = db_lookup_column, all.x = TRUE)
  colnames(df_with_id)[colnames(df_with_id) == "id"] <- final_column
  final_df <- df_with_id |>
    select(-model_lookup_column)

  return(final_df)
}
# insert_model_name <- function(con, model){
#   try({
#     model_final_results <- model$final_results
#     model_name <- unique(model_final_results$model_name)
#
#     if(length(model_name) != 1){
#       stop("Error: Model name does not contain exactly one unique value.")
#     }
#     query <- glue::glue_sql(
#       "INSERT INTO model_run (
#           model_name,
#         ) VALUES (
#           UNNEST(ARRAY[{model_name*}])
#         );",
#       .con = con
#     )
#     res <- DBI::dbExecute(con, query)
#     return(res)
#   }, silent = TRUE)
# }

insert_model_run <- function(con, model, blob_url, description){
  model_final_results <- data
  try({
    model_name_id <- join_lookup(model_final_results, "model_name", "model_name", "name", "model_name_id") |>
      select(model_name_id) |> unique()
    blob_storage_url <- blob_url
    srjpedata_version <- model_final_results$srjpedata_version |> unique()
    model_run_description <- description

    if(length(model_name_id) != 1){
      stop("Error: There are multiple model name being uploaded. Please only upload one.")
    }
    else if (length(srjpedata_version) != 1){
      stop("Error: There are multiple SRJPE data version being uploaded. Please only upload one.")
    }
    query <- glue::glue_sql(
      "INSERT INTO model_run (
          blob_storage_url,
          model_name_id,
          srjpedata_version,
          description
        ) VALUES (
          UNNEST(ARRAY[{blob_storage_url*}]),
          UNNEST(ARRAY[{model_name_id*}]),
          UNNEST(ARRAY[{srjpedata_version*}]),
          UNNEST(ARRAY[{model_run_description*}])
        );",
      .con = con
    )
    res <- DBI::dbExecute(con, query)
    return(res)
  })
}

insert_model_parameters <- function(con, model, blob_url) {

  model_final_results <- model$final_results
  model_fit_filename <- stringr::str_extract(model$full_object$model.file, "[^/]+$")
  model_final_results$blob_url <- blob_url
  model_final_results$model_fit_filename <- model_fit_filename
  model_final_results <- join_lookup(model_final_results, "model_run", "blob_url", "blob_storage_url", "model_run_id")
  model_final_results <- join_lookup(model_final_results, "model_location", "site", "site", "location_id")
  model_final_results <- join_lookup(model_final_results, "statistic", "statistic", "definition", "statistic_id")
  model_final_results <- join_lookup(model_final_results, "lifestage", "life_stage", "definition", "lifestage_id")
  model_final_results <- join_lookup(model_final_results, "parameter", "parameter", "definition", "parameter_id")

  model_final_results <-  model_final_results |>
    select(model_run_id, location_id, run_year, week_fit, lifestage_id, model_fit_filename, parameter_id, statistic_id, value)

  query <- glue::glue_sql(
    "INSERT INTO model_parameters (
          model_run_id,
          location_id,
          run_year,
          week_fit,
          lifestage_id,
          model_fit_filename,
          parameter_id,
          statistic_id,
          value
        ) VALUES (
          UNNEST(ARRAY[{model_final_results$model_run_id*}]),
          UNNEST(ARRAY[{model_final_results$location_id*}]),
          UNNEST(ARRAY[{model_final_results$run_year*}]),
          UNNEST(ARRAY[{model_final_results$week_fit*}]),
          UNNEST(ARRAY[{model_final_results$lifestage_id*}]),
          UNNEST(ARRAY[{model_final_results$model_fit_filename*}]),
          UNNEST(ARRAY[{model_final_results$parameter_id*}]),
          UNNEST(ARRAY[{model_final_results$statistic_id*}]),
          UNNEST(ARRAY[{model_final_results$value*}])
        );",
    .con = con
  )
  res <- DBI::dbExecute(con, query)

  return(res)
}

#TODO connect to jpedatabase, query for model name (make it unique?), pull model parameters pull lookup tables
load_model_fit <- function(con, model_name){
  model_id <- tbl(con, "model_name") |>
    filter(name == model_name) |>
    select(id) |>
    pull()

  model_run_id <- tbl(con, "model_run") |>
    filter(model_name_id == model_id) |>
    select(id) |>
    pull()

  model_parameters <- tbl(con, "model_parameters") |>
    filter(model_run_id == model_run_id) |>
    as_tibble()

  stat_lookup <- tbl(con, "statistic") |>
    as_tibble() |>
    rename("statistic_id" = "id")

  model_parameters <- model_parameters |> left_join(stat_lookup, by = "statistic_id") |>
    select(-c("statistic_id", "updated_at.x","updated_at.y", "description")) |>
    rename("statistic" = "definition")

  location_lookup <- tbl(con, "model_location") |>
    as_tibble() |>
    rename("location_id" = "id")

  model_parameters <- model_parameters |> left_join(location_lookup, by = "location_id") |>
    select(-c("location_id", "updated_at", "site", "description")) |>
    rename("location" = "stream")

  lifestage_lookup <- tbl(con, "lifestage") |>
    as_tibble() |>
    rename("lifestage_id" = "id")

  model_parameters <- model_parameters |> left_join(lifestage_lookup, by = "lifestage_id") |>
    select(-c("lifestage_id", "updated_at", "description")) |>
    rename("lifestage" = "definition")

  parameter_lookup <- tbl(con, "parameter") |>
    as_tibble() |>
    rename("parameter_id" = "id")

  model_parameters <- model_parameters |> left_join(parameter_lookup, by = "parameter_id") |>
    select(-c("parameter_id", "updated_at", "description")) |>
    rename("parameter" = "definition")

  return(model_parameters)
}

storage_account = "geneticsedidata"
container_name = "model-fits"
# url <- model_pin_board(storage_account, container_name)
# file_path <- "data/pCap_model_2025-01-09.rds"
model_fits <- readRDS("data/pCap_model_2025-01-09.rds")
store_model_fit(con,
                storage_account = storage_account,
                container_name = container_name,
                # access_key = access_key,
                data = model_fits,
                results_name = "new_model_pcap_test_v2",
                description = "test upload")
# insert_model_name(model_fits, con)
insert_model_run(con, model_fits, blob_url, "test upload")
# insert_model_parameters(con, model_fits, "www.azure.com/3")

# x <- load_model_fit(con, "missing_mark_recap.bug")
new_model_fits <- readRDS("data/pCap_model_2025-01-09.rds")
abundance_fit <- readRDS("data/abundance_model_inputs.rds")
old_model_fits <- readRDS("data/ubc_2004_2024-07-23.rds")


