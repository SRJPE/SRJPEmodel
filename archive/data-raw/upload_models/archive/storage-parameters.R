#' @title Get Most Recent Model Results
#' @description This function retrieves the most recent model results for each model name, site, year, stream, etc.
#' @param con A connection object to the database.
#' @return A tibble containing model parameters.
#'
#' @examples
#' \dontrun{
#' con <- dbConnect(RPostgres::Postgres(),
#'                  dbname = "your_db_name",
#'                  host = "your_host",
#'                  port = 5432,
#'                  user = "your_username",
#'                  password = "your_password")
#'
#' Example: get model parameters from Azure JPE Data Storage
#' model_results <- get_model_results_parameters(con)
#'
#' dbDisconnect(con)
#' }
#' @export
get_most_recent_model_results <- function(con) {

  results <- tbl(con, "model_parameters") |>
    # get site (location id), and pull matching stream from db
    left_join(tbl(con, "trap_location") |>
                dplyr::distinct(id, site, stream) |>
                select(location_id = id, site, stream),
              by = "location_id") |>
    left_join(tbl(con, "parameter") |>
                select(parameter_id = id, parameter = definition),
              by = "parameter_id") |>
    left_join(tbl(con, "statistic") |>
                select(statistic_id = id, statistic = definition),
              by = "statistic_id") |>
    left_join(tbl(con, "model_run") |>
                select(model_run_id = id, model_name_id),
              by = "model_run_id") |>
    left_join(tbl(con, "model_name") |>
                select(model_name = name, model_name_id = id),
              by = "model_name_id") |>
    # TODO join any other lookups needed for other model types (i.e. survival)
    collect() |>
    mutate(model_name = case_when(model_name %in% c("no_mark_recap", "no_mark_recap_no_trib",
                                                    "missing_mark_recap", "all_mark_recap") ~ "bt_spas_x",
                                  model_name %in% c("beta_dev_hbmrt", "beta_dv_hbmrt_lag1") ~ "inseason",
                                  model_name %in% c("survival_cov_wy", "survival_no_cov") ~ "survival",
                                  TRUE ~ model_name)) |>
    group_by(model_name, stream, site, year, updated_at) |>
    mutate(recent = cur_group_id()) |>
    ungroup() |>
    group_by(model_name, stream, site, year) |>
    mutate(most_recent_by_year = max(recent)) |>
    ungroup() |>
    # filter to only keep the most recent run for each year
    filter(recent == most_recent_by_year) |>
    select(-c(location_id, parameter_id, statistic_id,
              updated_at, location_fit_id,
              model_name_id, recent, most_recent_by_year)) |>
    dplyr::distinct(year, week_fit, value, stream, site, parameter, statistic, .keep_all = T) # TODO this is an error in the model parameter upload. we are joining on location_id but there are multiple matches for knights landing, so we are getting duplicates of all the parameter esitimates.

  return(results)
}

#' @title Get Model Parameter Results
#' @description This function retrieves the model parameter results from the JPE Database using a keyword(s) from the model run description or a model run ID.
#'
#' @param con A connection object to the database.
#' @param keyword An optional string used to search for the model run description in the database. If provided, it is used to find the corresponding blob URL.
#' @param model_run_id An optional ID for a specific model run. If provided, it is used to find the corresponding blob URL. If neither `keyword` nor `model_run_id` is provided, the latest model run is used.
#' @return A tibble containing model parameters.
#'
#' @examples
#' \dontrun{
#' con <- dbConnect(RPostgres::Postgres(),
#'                  dbname = "your_db_name",
#'                  host = "your_host",
#'                  port = 5432,
#'                  user = "your_username",
#'                  password = "your_password")
#'
#' # Example: get model parameters from Azure JPE Data Storage using a keyword
#' model_results <- get_model_results_parameters(con, keyword = "model description")
#' print(model_results)
#'
#' Example: get model parameters from Azure JPE Data Storage using a model run ID
#' model_results <- get_model_results_parameters(con, model_run_id = 12)
#'
#' dbDisconnect(con)
#' }
#' @export
get_model_results_parameters <- function(con, keyword=NULL, model_run_id=NULL){

  model_run_id <- search_model_run(con, keyword, model_run_id) |>
    dplyr::pull(id)

  model_parameters <- tbl(con, "model_parameters") |>
    filter(model_run_id == !!model_run_id) |>
    as_tibble()

  stat_lookup <- tbl(con, "statistic") |>
    as_tibble() |>
    rename("statistic_id" = "id")

  model_parameters <- model_parameters |> left_join(stat_lookup, by = "statistic_id") |>
    select(-c("statistic_id", "updated_at.x","updated_at.y", "description")) |>
    rename("statistic" = "definition")

  location_lookup <- tbl(con, "trap_location") |>
    as_tibble() |>
    rename("location_id" = "id")

  model_parameters <- model_parameters |> left_join(location_lookup, by = "location_id") |>
    select(-c("location_id", "updated_at", "description")) |>
    rename("location" = "site")


  parameter_lookup <- tbl(con, "parameter") |>
    as_tibble() |>
    rename("parameter_id" = "id")

  model_parameters <- model_parameters |> left_join(parameter_lookup, by = "parameter_id") |>
    select(-c("parameter_id", "updated_at", "description")) |>
    rename("parameter" = "definition")

  return(model_parameters)
}

#' @keywords internal
insert_model_parameters <- function(con, model, blob_url_list, results_name, inputs) {

  # call extract function based on the model name
  if(results_name %in% c("all_mark_recap", "no_mark_recap", "missing_mark_recap", "no_mark_recap_no_trib")) {
    model_final_results <- extract_abundance_estimates(inputs, model) |>
      rename(year = run_year)
  } else if(results_name %in% c("pcap_all", "pcap_mainstem")) {
    model_final_results <- extract_pCap_estimates(model, inputs) |>
      rename(year = run_year)
  } else if(results_name == "p2s") {
    model_final_results <- extract_P2S_estimates(inputs,
                                                 model)
  } else if(results_name == "stock_recruit") {
    model_final_results <- extract_stock_recruit_estimates(inputs,
                                                           model)
  } else if(results_name %in% c("beta_dev_hbmrt", "beta_dv_hbmrt_lag1")) {
    model_final_results <- extract_inseason_estimates(inputs,
                                                      model)
  } else if(results_name %in% c("survival_cov_wy", "survival_no_cov")) {
    model_final_results <- extract_survival_estimates(model)
  }
  model_final_results$blob_url <- blob_url_list$model_fit_url
  model_final_results <- join_lookup(model_final_results, "model_run", "blob_url", "blob_fit_storage_url", "model_run_id")
  model_final_results <- join_lookup(model_final_results, "statistic", "statistic", "definition", "statistic_id")
  model_final_results <- join_lookup(model_final_results, "parameter", "parameter", "definition", "parameter_id")

  # location_id is ALWAYS site, which determines stream.
  # TODO need to do this for survival model
  model_final_results <- model_final_results |>
    left_join(tbl(con, "trap_location") |>
                collect() |>
                dplyr::distinct(site, stream, .keep_all = T) |>
                rename(location_id = id) |>
                select(location_id, site, stream),
              by = c("stream", "site")) |>
    mutate(location_fit_id = NA_integer_)


  model_final_results <-  model_final_results |>
    select(model_run_id, location_id, year, week_fit, parameter_id, statistic_id, value, location_fit_id) |>
    mutate(
      year = as.integer(year),
      week_fit = as.integer(week_fit)
      )


  query <- glue::glue_sql(
    "INSERT INTO model_parameters (
          model_run_id,
          location_id,
          year,
          week_fit,
          parameter_id,
          statistic_id,
          value,
          location_fit_id
        ) VALUES (
          UNNEST(ARRAY[{model_final_results$model_run_id*}]),
          UNNEST(ARRAY[{model_final_results$location_id*}]),
          UNNEST(ARRAY[{model_final_results$year*}]),
          UNNEST(ARRAY[{model_final_results$week_fit*}]),
          UNNEST(ARRAY[{model_final_results$parameter_id*}]),
          UNNEST(ARRAY[{model_final_results$statistic_id*}]),
          UNNEST(ARRAY[{model_final_results$value*}]),
          UNNEST(ARRAY[{model_final_results$location_fit_id*}])
        );",
    .con = con
  )
  res <- DBI::dbExecute(con, query)

  return(res)
}