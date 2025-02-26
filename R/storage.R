#' @title Store Model Fits
#' @description This function stores a model fit object in an Azure Blob Container and the data into JPE database.
#' It uploads the model data to Azure Blob Storage and updates the database with relevant information about the model run and parameters.
#'
#' @param con A connection object to the database.
#' @param storage_account A string specifying the Azure storage account name, typically `jpemodelresults`.
#' @param container_name A string specifying the container name in the Azure storage account, typically `model-results`.
#' @param access_key A string specifying the Azure storage access key with write permissions.
#' @param model_fit_object The model result object that needs to be stored, of class `stanfit` or `bugs`.
#' @param results_name A string specifying a name to identify the model results in Azure Blob Storage. One of `all_mark_recap`,
#'`missing_mark_recap`, `no_mark_recap`,  `no_mark_recap_no_trib`, `pcap_all`, `pcap_mainstem`, or `p2s`.
#' @param site If uploading an abundance model output, or a pCap mainstem object, you must supply the site for which you fit the model.
#' @param run_year If uploading an abundance model output, you must supply the run_year for which you fit the model.
#' @param ... Additional named arguments to be passed as metadata to the blob storage.
#'
#' @return A string representing the URL of the blob in Azure Blob Storage where the model results are stored.
#'
#' @examples
#' \dontrun{
#' con <- DBI::dbConnect(RPostgres::Postgres(),
#'                        dbname = cfg$db_name,
#'                        host = cfg$db_host,
#'                        port = cfg$db_port,
#'                        user = cfg$db_user,
#'                        password = cfg$db_password)
#'
#' blob_url <- store_model_fit(con,
#'                             storage_account = "my_storage_account",
#'                             container_name = "my_container",
#'                             access_key = "my_access_key",
#'                             model_fit_object = model_results,
#'                             results_name = "model_name",
#'                             description = "model_description")
#'
#' print(blob_url)
#'
#' dbDisconnect(con)
#'}
#' @export
store_model_fit <- function(con, storage_account, container_name, access_key, model_fit_object, results_name, site = NULL, run_year = NULL, description, ...){

  # checks
  # check that they supply an approved results name
  if(!results_name %in% c("all_mark_recap", "no_mark_recap", "missing_mark_recap", "no_mark_recap_no_trib",
                          "pcap_all", "pcap_mainstem", "p2s")) {
    cli::cli_abort("You must supply an approved value to the results_name argument.")
  }

  # check that the model objects align with the model type (bugs or stanfit)
  if(results_name %in% c("all_mark_recap", "no_mark_recap", "missing_mark_recap", "no_mark_recap_no_trib") & class(model_fit_object) != "bugs") {
    cli::cli_abort("For models all_mark_recap, no_mark_recap, missing_mark_recap, and no_mark_recap_no_trib, model object must be of class 'bugs'.")
  } else {
    if(class(model_fit_object) != "stanfit") {
      cli::cli_abort("For models pcap_all, pcap_mainstem, and p2s, model object must be of class 'stanfit'.")
    }
  }

  # check that they supply a site and run year if supplying an abundance model fit or pCap mainstem
  if(results_name %in% c("all_mark_recap", "no_mark_recap", "missing_mark_recap", "no_mark_recap_no_trib")) {
    if(is.null(site) | is.null(run_year)){
      cli::cli_abort("You must supply a site and run year if uploading an abundance model result.")
    }
  }
  if(results_name == "pcap_mainstem" & is.null(site)) {
    cli::cli_abort("You must supply a site (either knights landing or tisdale) if uploading a pcap_mainstem model result.")
  }

  model_board <- model_pin_board(storage_account, container_name)

  blob_url <- pin_model_data(
    model_board,
    model_fit_object,
    name = results_name,
    ...
  )

  total_run_rows <- insert_model_run(con, model_fit_object, blob_url, description, results_name, site, run_year)
  message(glue::glue("Inserted new model run into database."))
  total_rows <- insert_model_parameters(con, model_fit_object, blob_url, results_name, site, run_year)
  message(glue::glue("Inserted {total_rows} into database. Uploaded model fit results to {blob_url}."))

  return(blob_url)
}

#' @title Azure Blob Setup
#' @description
#' Creates an `AzureStor` storage container object to be used by the `pins` package.
#'
#' @param storage_account string storage account to create connection to
#' @param access_key the access key that provides for the blob access. By default an env variable
#' with the name `AZ_CONTAINER_ACCESS_KEY` is retrived, optionally you can directly pass in a value.
#'
#' @returns An `AzureStore` storage container object
#' @keywords internal
setup_azure_blob_backend <- function(storage_account, access_key=Sys.getenv("AZ_CONTAINER_ACCESS_KEY")) {
  if (access_key == "") {
    stop("access key is required in order to write and read from the azure board", call. = FALSE)
  }
  storage_endpoint <- glue::glue("https://{storage_account}.blob.core.windows.net")
  store <- AzureStor::storage_endpoint(storage_endpoint, key = access_key)
  return(store)
}

#' @title Model Data Pin Board
#' @description Creates a pin board for model data storage in Azure Blob Storage.
#'
#' @param storage_account The name of the Azure Storage account.
#' @param container The name of the blob container within the storage account.
#' @param ... Additional arguments passed to `setup_azure_blob_backend()`.
#'
#' @return A pins board object connected to the specified Azure Blob Storage container.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' board <- model_pin_board("my_storage_account", "my_container")
#' }
#'
#' @seealso
#' \code{\link[AzureStor]{blob_container}}
#' \code{\link[pins]{board_azure}}
#' @keywords internal
model_pin_board <- function(storage_account, container, ...) {
  storage_client <- setup_azure_blob_backend(storage_account, ...)

  all_containers <- AzureStor::list_blob_containers(storage_client)
  container_names <- names(all_containers)

  if (!(container %in% container_names)) {
    stop(glue::glue("the blob '{container}' does not exist in your azure storage account"), call. = FALSE)
  }

  blob_client <- AzureStor::blob_container(storage_client, container)
  blob_board <- pins::board_azure(container = blob_client, "model-fits/")

  return(blob_board)
}

#' @title Model Data Pin
#' @description
#' Pin model data to provided pin board.
#'
#' @param board a board created with `model_pin_board`
#' @param data data to pin
#' @param name the name assigned to this data on blob, if wanting to add a new version of a previous result, re-use the name
#' @param title title for data
#' @param ... additional named arguments passed is as metadata to the pin
#'
#' @examples
#' \dontrun{
#'
#' # first create a board connected to the azure account you with to use
#' model_board <- model_pin_board("storage-account-name", "model-results")
#'
#' # next run a model
#' bt_spas_x_results <- run_single_bt_spas_x(SRJPEmodel::bt_spas_x_bayes_params,
#'   bt_spas_x_input_data = SRJPEdata::weekly_juvenile_abundance_model_data,
#'   site = "ubc",
#'   run_year = 2004,
#'   lifestage = "fry",
#'   effort_adjust = F,
#'   bugs_directory = here::here("data-raw", "WinBUGS14"),
#'   debug_mode = FALSE, no_cut = T
#' )
#'
#' # now pin the model
#' pin_model_data(
#'   model_board,
#'   bt_spas_x_results,
#'   name = "BT Spas at Butte",
#'   title = "bt spas results",
#'   description = "model results dataframe and model object"
#' )
#'
#' # run again
#' bt_spas_x_results <- run_single_bt_spas_x(...)
#'
#' # pin with same name to create a new version of the data, all previous versions are stored
#' pin_model_data(
#'   model_board,
#'   bt_spas_x_results,
#'   name = "BT Spas at Butte",
#'   title = "bt spas results",
#'   description = "model results dataframe and model object"
#' )
#'
#' # or assign it a new name to create a new object on the board
#' pin_model_data(
#'   model_board,
#'   bt_spas_x_results,
#'   name = "BT Spas at Butte V2",
#'   title = "bt spas results",
#'   description = "model results dataframe and model object"
#' )
#'
#' # search storeage for pins
#' model_board |> pins::pin_search()
#'
#' # A tibble: 3 × 6
#' #   name                 type  title                created               file_size meta
#' #   <chr>                <chr> <chr>                <dttm>              <fs::bytes> <list>
#' # 1 full_model_object    rds   Model with new value 2024-07-31 14:02:13       2.42M <pins_met>
#' # 2 full_model_object_V2 rds   Model with new value 2024-07-31 14:26:42       2.42M <pins_met>
#' # 3 some numbers         rds   some numbers         2024-07-31 14:23:15          61 <pins_met>
#'
#' # retrieve a pin
#' model_data <- pins::pin_read(model_board, "full_model_object")
#' }
#'
#'
#' @export
#' @md
pin_model_data <- function(board, data, name, title = NULL, description = NULL, ...) {
  pin_metadata <- list(...)
  data_name <- pins::pin_write(board,
                               data,
                               name = name,
                               title = title,
                               description = description,
                               metadata = pin_metadata, type = "rds"
  )

  latest_version_df <- board |> pins::pin_versions(data_name)
  latest_version <- latest_version_df$version[1]

  data_url <- glue::glue("{board$container$endpoint$url}/{board$path}/{data_name}/{latest_version}/{data_name}.rds")

  return(data_url)
}

#' @title Search Model Run
#' @description This function searches for model runs in the JPE database using criteria, such as keywords, model run IDs, or an option to view all model runs.
#'
#' @param con A database connection object.
#' @param keyword An optional keyword to search within the model run descriptions.
#' @param model_run_id An optional ID to search for a specific model run. If both `keyword` and `model_run_id` are NULL, the search will pull the latest model_run.
#' @param view_all A boolean indicating if all model runs should be displayed, ignoring other filters (default is FALSE).
#'
#' @return A data frame containing the search results. If `view_all` is set to TRUE, it returns all model runs. If a keyword or model_run_id is provided, it filters the results accordingly.
#'
#' @examples
#' \dontrun{
#' # Example: Search for the latest model run
#' latest_model_run <- search_model_run(con)
#'
#' # Example: Search model runs by keyword
#' keyword_search_results <- search_model_run(con, keyword = "regression")
#'
#' # Example: Search for a specific model run by ID
#' specific_model_run <- search_model_run(con, model_run_id = 5)
#'
#' # Example: View all model runs
#' all_model_runs <- search_model_run(con, view_all = TRUE)
#' }
#' @export
search_model_run <- function(con, keyword=NULL, model_run_id=NULL, view_all=FALSE){

  model_name_table <- tbl(con, "model_name")
  model_run_table <- tbl(con, "model_run")
  model_run_table <- model_run_table |>
    left_join(model_name_table, by = c("model_name_id" = "id"), suffix=c("", ".y")) |>
    select(-ends_with(".y")) |>
    rename(model_name = name) |>
    select(-model_name_id) |>
    collect()

  if (view_all){
    all_model_run <- model_run_table |>
      arrange(desc("updated_at"))

    return(all_model_run)
  } else{
    if (is.null(keyword) && is.null(model_run_id)){
      latest_model_run <- model_run_table |>
        arrange(desc("updated_at")) |>
        head(1)

      return(latest_model_run)

    }else if(view_all){


    }else if(!is.null(keyword)){
      model_run_id <- model_run_table |>
        filter(stringr::str_detect(tolower(description), tolower(keyword))) |>
        collect()
      if (nrow(model_run_id)>1){

        cli::cli_alert_info("Multiple model runs found with the keyword.")
        cli::cli_ul()
        cli::cli_li("{.emph id | description | blob_storage_url | srjpedata_version | updated_at }")
        for (i in seq_len(nrow(model_run_id))){
          cli::cli_li("{model_run_id$id[i]} | {model_run_id$description[i]} | {model_run_id$blob_storage_url[i]} | {model_run_id$srjpedata_version[i]} | {model_run_id$updated_at[i]}")
        }
        cli::cli_end()

        stop("Search term resulted in more than 1 model run. Please see results in the details above and use the corresponding id to pull the specific model you want.")
      }else if (nrow(model_run_id) == 0) {
        stop("No model runs found with the given keyword. Please try another keyword.")
      }

      return(model_run_id)

    }else{
      model_run <- model_run_table |>
        filter(id == model_run_id)

      if (nrow(model_run) == 0){
        stop("There is no model run with this ID. Please enter another model run ID.")
      }

      return(model_run)

    }
  }
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
    filter(model_run_id == model_run_id) |>
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
    select(-c("location_id", "updated_at", "site", "description")) |>
    rename("location" = "stream")


  parameter_lookup <- tbl(con, "parameter") |>
    as_tibble() |>
    rename("parameter_id" = "id")

  model_parameters <- model_parameters |> left_join(parameter_lookup, by = "parameter_id") |>
    select(-c("parameter_id", "updated_at", "description")) |>
    rename("parameter" = "definition")

  return(model_parameters)
}

#' @title Get Model Object
#' @description This function retrieves the model object (.Rds file) from Azure Blob Storage using a keyword from the model run description or a model run ID.
#'
#' @param con A connection object to the database.
#' @param access_key A string specifying the Azure storage access key with read permissions. By default, it retrieves from the environment variable `AZ_CONTAINER_ACCESS_KEY`.
#' @param keyword An optional string used to search for the model run description in the database. If provided, it is used to find the corresponding blob URL.
#' @param model_run_id An optional ID for a specific model run. If provided, it is used to find the corresponding blob URL. If neither `keyword` nor `model_run_id` is provided, the latest model run is used.
#'
#' @return The model object retrieved from the Azure Blob Storage.
#' #' @examples
#' \dontrun{
# Example: Get a model object from Azure Blob Storage using a keyword
#' model_object <- get_model_object(con = db_connection, keyword = "model description")
#'
#' # Example: Get a model object from Azure Blob Storage using a model run ID
#' model_object <- get_model_object(con = db_connection, model_run_id = 12)
#'
#' print(model_object$model.file)
#' }
#' @export
get_model_object <- function(con, keyword=NULL, model_run_id=NULL, access_key=Sys.getenv("AZ_CONTAINER_ACCESS_KEY")){
  if (access_key == "") {
    stop("Access key is required to read from Azure Blob Storage.", call. = FALSE)
  }
  model_run_url <- search_model_run(con, keyword, model_run_id) |>
    dplyr::pull(blob_storage_url)

  storage_account <- sub("https://(.+?)\\.blob\\.core\\.windows\\.net.*", "\\1", model_run_url)
  container_name <- sub("https://.+\\.blob\\.core\\.windows\\.net/(.+?)/.*", "\\1", model_run_url)
  blob_path <- tools::file_path_sans_ext(sub("^.*/", "", model_run_url))

  model_board <- model_pin_board(storage_account, container_name)
  model_object <- pins::pin_read(model_board, blob_path)

  return(model_object)
}

#' @keywords internal
join_lookup <- function(df, db_table, model_lookup_column, db_lookup_column, final_column){
  db_lookup_column_sym <- rlang::sym(db_lookup_column)
  model_lookup_column_sym <- rlang::sym(model_lookup_column)
  final_column_sym <- rlang::sym(final_column)

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

#' @keywords internal
insert_model_run <- function(con, model, blob_url, description, results_name, site = NULL, run_year = NULL){

  # call extract function based on the model name
  if(results_name %in% c("all_mark_recap", "no_mark_recap", "missing_mark_recap", "no_mark_recap_no_trib")) {
    model_final_results <- extract_abundance_estimates(site, run_year,
                                                       prepare_abundance_inputs(site, run_year, effort_adjust = T), model)
  } else if(results_name == "pcap_all") {
    model_final_results <- extract_pCap_estimates(model, prepare_pCap_inputs(mainstem = FALSE))
  } else if(results_name == "pcap_mainstem") {
    model_final_results <- extract_pCap_estimates(model, prepare_pCap_inputs(mainstem = TRUE, mainstem_site = site))
  } else if(results_name == "p2s") {
    model_final_results <- extract_P2S_estimates(model)
  }

  # model_final_results <- model$final_results
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

#' @keywords internal
insert_model_parameters <- function(con, model, blob_url, results_name, site = NULL, run_year = NULL) {

  # call extract function based on the model name
  if(results_name %in% c("all_mark_recap", "no_mark_recap", "missing_mark_recap", "no_mark_recap_no_trib")) {
    model_final_results <- extract_abundance_estimates(site, run_year,
                                                       prepare_abundance_inputs(site, run_year, effort_adjust = T), model)
  } else if(results_name == "pcap_all") {
    model_final_results <- extract_pCap_estimates(model, prepare_pCap_inputs(mainstem = FALSE))
  } else if(results_name == "pcap_mainstem") {
    model_final_results <- extract_pCap_estimates(model, prepare_pCap_inputs(mainstem = TRUE, mainstem_site = site))
  } else if(results_name == "p2s") {
    model_final_results <- extract_P2S_estimates(model)
  }
  # model_fit_filename <- stringr::str_extract(model$full_object$model.file, "[^/]+$")
  # model_final_results$model_fit_filename <- file_name
  model_final_results$blob_url <- blob_url


  model_final_results <- join_lookup(model_final_results, "model_run", "blob_url", "blob_storage_url", "model_run_id")
  model_final_results <- join_lookup(model_final_results, "trap_location", "site", "site", "location_id")
  model_final_results <- join_lookup(model_final_results, "statistic", "statistic", "definition", "statistic_id")
  # model_final_results <- join_lookup(model_final_results, "lifestage", "life_stage", "definition", "lifestage_id")
  model_final_results <- join_lookup(model_final_results, "parameter", "parameter", "definition", "parameter_id")
  model_final_results <- join_lookup(model_final_results, "trap_location", "location_fit", "site", "location_fit_id")

  model_final_results <-  model_final_results |>
    select(model_run_id, location_id, run_year, week_fit, parameter_id, statistic_id, value, location_fit_id) |>
    mutate(
      run_year = as.integer(run_year),
      week_fit = as.integer(week_fit)
      )


  query <- glue::glue_sql(
    "INSERT INTO model_parameters (
          model_run_id,
          location_id,
          run_year,
          week_fit,
          parameter_id,
          statistic_id,
          value,
          location_fit_id
        ) VALUES (
          UNNEST(ARRAY[{model_final_results$model_run_id*}]),
          UNNEST(ARRAY[{model_final_results$location_id*}]),
          UNNEST(ARRAY[{model_final_results$run_year*}]),
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
