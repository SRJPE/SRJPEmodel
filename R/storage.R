#' @title Store Model Fits
#' @description This function stores a model fit object in an Azure Blob Container and the data into JPE database.
#' It uploads the model data to Azure Blob Storage and updates the database with relevant information about the model run and parameters.
#'
#' @param con A connection object to the database.
#' @param storage_account A string specifying the Azure storage account name, typically `jpemodelresults`.
#' @param container_name A string specifying the container name in the Azure storage account, typically `model-results`.
#' @param access_key A string specifying the Azure storage access key with write permissions. Default is stored in your R environment under `AZ_CONTAINER_ACCESS_KEY`.
#' @param model_fit_object The model result object that needs to be stored, of class `stanfit` or `bugs`.
#' @param model_inputs The inputs used to fit the `model_fit_object`.
#' @param results_name A string specifying a name to identify the model results in Azure Blob Storage. One of `bt_spas_x`,
#' `pcap_all`, `pcap_mainstem`, `p2s`, `stock_recruit`, `inseason`, `survival`.
#' @param description A description of the model fit you are uploading.
#' @param ... Additional named arguments to be passed as metadata to the blob storage.
#'
#' @return A string representing the URL of the blob in Azure Blob Storage where the model fit object, diagnostic plots, and inputs are stored.
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
#'                             model_inputs = model_inputs,
#'                             results_name = "model_name",
#'                             description = "model_description")
#'
#' print(blob_url)
#'
#' dbDisconnect(con)
#'}
#' @export
store_model_fit <- function(con, storage_account = "jpemodelresults", container_name = "model-results", access_key = Sys.getenv("AZ_CONTAINER_ACCESS_KEY"),
                            model_fit_object, model_inputs, results_name, description, ...){


  # extracts correct submodel name from "abundance" (assuming user does not know the specifics)
  if(results_name == "bt_spas_x") {
    results_name <- model_inputs$model_name
  }
  if(results_name == "inseason") {
    results_name <- model_inputs$model_name
  }
  if(results_name == "survival") {
    results_name <- model_inputs$model_name
  }
  # check that they supply an approved results name
  if(!results_name %in% c("pcap_all", "pcap_mainstem", "p2s", "stock_recruit",
                          "bt_spas_x", "inseason", "survival",
                          "all_mark_recap", "no_mark_recap", "missing_mark_recap", "no_mark_recap_no_trib",
                          "survival_cov_wy", "survival_no_cov")) {
    cli::cli_abort("You must supply an approved value to the results_name argument.")
  }

  # check that the model objects align with the model type (bugs or stanfit)
  if(results_name %in% c("all_mark_recap", "no_mark_recap", "missing_mark_recap", "no_mark_recap_no_trib") & class(model_fit_object) != "bugs") {
    cli::cli_abort("For models all_mark_recap, no_mark_recap, missing_mark_recap, and no_mark_recap_no_trib, model object must be of class 'bugs'.")
  }
  if(results_name %in% c("pcap_all", "pcap_mainstem", "p2s", "stock_recruit",
                         "beta_dev_hbmrt", "beta_dv_hbmrt_lag1", "survival_cov_wy", "survival_no_cov") & class(model_fit_object) != "stanfit") {
      cli::cli_abort("For models pcap_all, pcap_mainstem, p2s, and stock_recruit, model object must be of class 'stanfit'.")
  }

  # generate diagnostic plot
  model_plot <- generate_diagnostic_plot(model_inputs, model_fit_object)
  # print(model_plot)
  # storage workflow
  model_board <- model_pin_board(storage_account, container_name, results_name)

  blob_url_list <- pin_model_data(
    model_board,
    model_fit_object,
    model_inputs,
    model_plot,
    name = results_name,
    ...
  )

  tryCatch({
    total_run_rows <- insert_model_run(con, model_fit_object, blob_url_list, description, results_name,
                                       model_inputs)

    if (!is.null(total_run_rows) && total_run_rows == 1) {
      message(glue::glue("Inserted new model run into database."))

      total_rows <- insert_model_parameters(con, model_fit_object, blob_url_list, results_name, model_inputs)
      message(glue::glue("Inserted {total_rows} into database. Uploaded model fit results to {blob_url_list$model_fit}."))
    } else {
      message(glue::glue("⚠️ No new model run inserted into database. Skipping parameter insert."))
    }
  }, error = function(e) {
    message(glue::glue("❌ Error during model run insertion: {e$message}"))
    stop(e)
  })

  return(blob_url_list)
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
model_pin_board <- function(storage_account, container, model_name, ...) {
  storage_client <- setup_azure_blob_backend(storage_account, ...)

  all_containers <- AzureStor::list_blob_containers(storage_client)
  container_names <- names(all_containers)

  if (!(container %in% container_names)) {
    stop(glue::glue("the blob '{container}' does not exist in your azure storage account"), call. = FALSE)
  }

  blob_client <- AzureStor::blob_container(storage_client, container)
  blob_board <- pins::board_azure(container = blob_client, paste0("model-fits/", model_name))

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
pin_model_data <- function(board, data, model_input, model_plot, name, title = NULL, description = NULL, ...) {
  pin_metadata <- list(...)
  model_data_name <- pins::pin_write(board,
                               data,
                               name = paste0(name, "_fit"),
                               title = title,
                               description = description,
                               metadata = pin_metadata, type = "rds"
  )
  model_input_name <- pins::pin_write(board,
                                      model_input,
                                      name = paste0(name, "_input"),
                                     # title = title,
                                     # description = description,
                                     metadata = pin_metadata, type = "rds"
  )
  model_plot_name <- pins::pin_write(board,
                                     model_plot,
                                     name = paste0(name, "_plot"),
                                     # title = title,
                                     # description = description,
                                     metadata = pin_metadata, type = "rds"
  )

  fit_latest_version_df <- board |> pins::pin_versions(model_data_name) |>
    arrange(desc(created))
  fit_latest_version <- fit_latest_version_df$version[1]
  print(paste0("model fit version:", fit_latest_version))

  input_latest_version_df <- board |> pins::pin_versions(model_input_name) |>
    arrange(desc(created))
  input_latest_version <- input_latest_version_df$version[1]
  print(paste0("model input version:", input_latest_version))

  plot_latest_version_df <- board |> pins::pin_versions(model_plot_name) |>
    arrange(desc(created))
  plot_latest_version <- plot_latest_version_df$version[1]
  print(paste0("model plot version:", plot_latest_version))

  model_fit_url <- glue::glue("{board$container$endpoint$url}/{board$path}/{model_data_name}/{fit_latest_version}/{model_data_name}.rds")
  model_input_url <- glue::glue("{board$container$endpoint$url}/{board$path}/{model_input_name}/{input_latest_version}/{model_input_name}.rds")
  model_plot_url <- glue::glue("{board$container$endpoint$url}/{board$path}/{model_plot_name}/{plot_latest_version}/{model_plot_name}.rds")

  url_list <- list(
    "model_fit_url"=model_fit_url,
    "model_input_url"=model_input_url,
    "model_plot_url"=model_plot_url
  )
  return(url_list)
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

#' @title Get Model Object
#' @description This function retrieves the model object (.Rds file) from Azure Blob Storage using a keyword from the model run description or a model run ID.
#'
#' @param con A connection object to the database.
#' @param access_key A string specifying the Azure storage access key with read permissions. By default, it retrieves from the environment variable `AZ_CONTAINER_ACCESS_KEY`.
#' @param model_component Select one of fit, input or plot object to pull.
#' @param keyword An optional string used to search for the model run description in the database. If provided, it is used to find the corresponding blob URL.
#' @param model_run_id An optional ID for a specific model run. If provided, it is used to find the corresponding blob URL. If neither `keyword` nor `model_run_id` is provided, the latest model run is used.
#'
#' @return The model object retrieved from the Azure Blob Storage.
#' #' @examples
#' \dontrun{
# Example: Get a model object from Azure Blob Storage using a keyword
#' model_plot_object <- get_model_object(con = db_connection, model_component = "plot", keyword = "model description")
#'
#' # Example: Get a model object from Azure Blob Storage using a model run ID
#' model_input_object <- get_model_object(con = db_connection, model_component = "input", model_run_id = 12)
#'
#' print(model_object$model.file)
#' }
#' @export
get_model_object <- function(con, model_component="model_fit", model_run_id=NULL, keyword=NULL, access_key=Sys.getenv("AZ_CONTAINER_ACCESS_KEY")){
  if (access_key == "") {
    stop("Access key is required to read from Azure Blob Storage.", call. = FALSE)
  }
  if(!model_component %in% c("model_fit", "model_input", "model_plot")) {
    cli::cli_abort("Model component must be model_fit, model_input, or model_plot.")
  }
  if (model_component == "model_fit"){
    model_run_url <- search_model_run(con, keyword, model_run_id) |>
      dplyr::pull(blob_fit_storage_url)
  }else if (model_component == "model_input"){
    model_run_url <- search_model_run(con, keyword, model_run_id) |>
      dplyr::pull(blob_input_storage_url)
  } else {
    model_run_url <- search_model_run(con, keyword, model_run_id) |>
      dplyr::pull(blob_plot_storage_url)
  }

  storage_account <- sub("https://(.+?)\\.blob\\.core\\.windows\\.net.*", "\\1", model_run_url)
  container_name <- "model-results"
  model_id <- sub("^.*/model-fits/([^/]+)/.*$", "\\1", model_run_url)
  model_name <- tools::file_path_sans_ext(sub("^.*/", "", model_run_url))

  model_board <- model_pin_board(storage_account, container_name, model_id)
  version <- sub("^.*/([^/]+)/[^/]+$", "\\1", model_run_url)
  model_object <- pins::pin_read(model_board, model_name, version = version)

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
insert_model_run <- function(con, model, blob_url_list, description, results_name, inputs){

  # call extract function based on the model name
  if(results_name %in% c("all_mark_recap", "no_mark_recap", "missing_mark_recap", "no_mark_recap_no_trib")) {
    model_final_results <- extract_abundance_estimates(inputs, model)
  } else if(results_name == "pcap_all") {
    model_final_results <- extract_pCap_estimates(model, inputs)
  } else if(results_name == "pcap_mainstem") {
    model_final_results <- extract_pCap_estimates(model, inputs)
  } else if(results_name == "p2s") {
    model_final_results <- extract_P2S_estimates(inputs,
                                                 model)
  } else if(results_name == "stock_recruit") {
    model_final_results <- extract_stock_recruit_estimates(inputs,
                                                           model)
  } else if(results_name %in% c("beta_dev_hbmrt", "beta_dv_hbmrt_lag1")) {
    model_final_results <- extract_inseason_estimates(inputs,
                                                      model)
  } else if(results_name %in% c("survival_no_cov", "survival_cov_wy")) {
    model_final_results <- extract_survival_estimates(model)
  }
  # TODO check that model names match user input model names
  model_final_results$model_name <- tolower(model_final_results$model_name)
  # model_final_results <- model$final_results
  try({
    model_name_id <- join_lookup(model_final_results, "model_name", "model_name", "name", "model_name_id") |>
      select(model_name_id) |> unique()
    blob_fit_storage_url <- blob_url_list$model_fit_url
    blob_input_storage_url <- blob_url_list$model_input_url
    blob_plot_storage_url <- blob_url_list$model_plot_url
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
          blob_fit_storage_url,
          blob_input_storage_url,
          blob_plot_storage_url,
          model_name_id,
          srjpedata_version,
          description
        ) VALUES (
          UNNEST(ARRAY[{blob_fit_storage_url*}]),
          UNNEST(ARRAY[{blob_input_storage_url*}]),
          UNNEST(ARRAY[{blob_plot_storage_url*}]),
          UNNEST(ARRAY[{model_name_id*}]),
          UNNEST(ARRAY[{srjpedata_version*}]),
          UNNEST(ARRAY[{model_run_description*}])
        );",
      .con = con
    )
    result <- DBI::dbExecute(con, query)
    return(result)
  })
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

#' @title Get Most Recent Model Objects
#' @description This function retrieves the most recent model objects for each model name, site, year, stream, etc. You can
#' specify model_object (the full fit object, either `BUGS` or `stanfit`, the model plot (a posterior predictive check plot), or
#' the inputs used to fit the associated model object.
#' @param con A connection object to the database.
#' @param model_component A choice of `model_fit`, `model_input` or `model_plot` to pull.
#' @param model_name Can be left empty. If provided, must be one of `pcap_all`, `bt_spas_x`,
#' `pcap_mainstem`, `p2s`, `inseason`, or `stock_recruit`
#' @param stream Can be left empty. If provided, must be one of `battle creek`,
#' `clear creek`, `deer creek`, `mill creek`, `feather river`, `butte creek`,
#' `sacramento river`, or `yuba river`.
#' @return A the most recent model fit objects, inputs, or plots. The format will be a named list,
#' where each element is named by the `model_run_id`, `model_name`, `site`, `stream`, and `year`.
#' @export
get_most_recent_model_objects <- function(con, model_component="model_fit",
                                          access_key=Sys.getenv("AZ_CONTAINER_ACCESS_KEY"),
                                          model_name = NULL,
                                          stream = NULL) {

  recent_models <- SRJPEmodel::get_most_recent_model_results(con)

  if(!missing(model_name)) {

    if(!model_name %in% c("inseason", "pcap_all", "bt_spas_x",
                          "pcap_mainstem", "p2s", "stock_recruit", "survival")) {
      cli::abort("You must provide an approved model name.")
    }
    recent_models <- recent_models |>
      filter(model_name == !!model_name)
  }
  if(!missing(stream)) {
    recent_models <- recent_models |>
      filter(stream == !!stream)
  }

  most_recent_ids <- recent_models |>
    dplyr::distinct(model_run_id) |>
    dplyr::pull(model_run_id)

  model_ids <- recent_models |>
    dplyr::distinct(model_run_id, .keep_all = T) |>
    mutate(model_id = paste(model_run_id, model_name, site, stream, year, sep = "_")) |>
    dplyr::pull(model_id)

  model_object_list <- lapply(most_recent_ids, function(x) {
    cli::cli_bullets(paste("Pulling", which(most_recent_ids == x), "of", length(most_recent_ids), "objects"))
    result <- get_model_object(con, model_component, model_run_id = x)
  })
  names(model_object_list) <- model_ids

  return(model_object_list)
}

#' Diagnostic plots
#' @details This function produces a posterior predictive check plot.
#' @param model_results_name The model type name associated with options in the database: either
#' `all_mark_recap`, `no_mark_recap`, `missing_mark_recap`, `no_mark_recap_no_trib`,
#' `pcap_all`, `pcap_mainstem`, `p2s`, or `stock_recruit`.
#' @param inputs inputs for the fit model, created by running the appropriate `prepare_inputs` function.
#' @param fit the model fit
#' @returns Generates a diagnostic posterior predictive plot for a model fit.
#' @export
#' @md
generate_diagnostic_plot <- function(inputs, fit) {

  dark_JPE <- c("#F5CAC2", "#6E9881", "#9A8723", "#2D4755", "#869AA0")

  model_results_name = inputs$model_name
  print(model_results_name)

  # observed variable is more complicated, need some basic calcs
  if(model_results_name %in% c("pcap_mainstem", "pcap_all")) {
    basic_pcap <- inputs$inputs$data$Recaptures / inputs$inputs$data$Releases
    obsv_variable <- qlogis(basic_pcap)
    pred_variable <- "logit_pCap"
    pred_variable_name <- "Proportion captured (logit)"
    location <- inputs$location
  } else if (model_results_name %in% c("all_mark_recap", "no_mark_recap",
                                       "missing_mark_recap", "no_mark_recap_no_trib")) {
    pCap_mu <- plogis(inputs$lt_pCap_Us$lt_pCap_mu)
    obsv_variable <- inputs$inputs$data$u / pCap_mu[inputs$inputs$data$Nstrata_wc]
    pred_variable <- "N"
    pred_variable_name <- "Weekly juvenile abundance"
    location <- paste(inputs$site, inputs$run_year, sep = "-")
  } else if (model_results_name == "p2s") {
    obsv_variable <- inputs$inputs$data$observed_spawners
    pred_variable <- "predicted_spawners"
    pred_variable_name <- "Spawner abundance"
    location <- inputs$stream
  } else if (model_results_name == "stock_recruit") {
    obsv_variable <- inputs$inputs$data$mu_obslgRS
    pred_variable <- "pred_lgRS"
    pred_variable_name <- "Mean recruits-per-spawner (log scale)"
    location <- inputs$stream
  } else if(model_results_name %in% c("beta_dev_hbmrt", "beta_dv_hbmrt_lag1")) {
    obsv_variable <- inputs$inputs$data$Nx_mu |>
      rowMeans()
    pred_variable <- "pred_pNx"
    pred_variable_name <- "Weekly cumulative juvenile abundance (averaged across years)"
    location <- inputs$stream
  }# TODO add survival

  # Extract posterior samples for predicted values
  n_posterior_samples <- 100

  if(model_results_name %in% c("all_mark_recap", "no_mark_recap",
                               "missing_mark_recap", "no_mark_recap_no_trib")) {
    extract_preds <- fit$sims.list$N[ , inputs$inputs$data$Uwc_ind] # predicted catch for weeks with obsv catch
    y_rep <- extract_preds[sample(nrow(extract_preds), n_posterior_samples), ]
  } else if(model_results_name %in%c("beta_dev_hbmrt", "beta_dv_hbmrt_lag1")) {
    posterior_samples <- rstan::extract(fit)$pred_pNx
    # take weeks estimated and average across years
    y_rep <- posterior_samples[sample(nrow(posterior_samples), n_posterior_samples), , ] |>
      # average across years
      apply(c(1, 2), mean, na.rm = T)
    y_rep[is.nan(y_rep)] <- as.numeric(0) # set NaNs to 0
  } else {
    posterior_samples <- rstan::extract(fit)
    extract_preds <- eval(parse(text = paste0("posterior_samples$", pred_variable)))
    y_rep <- extract_preds[sample(nrow(extract_preds), n_posterior_samples), ]
  }

  # Basic posterior predictive check
  ppc_plot <- bayesplot::ppc_dens_overlay(obsv_variable, y_rep) +
    theme_minimal() +
    labs(
      title = "Posterior Predictive Check",
      subtitle =
      "Is this model representative of the data? If the light blue lines (yrep = simulated data) look like the \ndark blue line (y = observed data) then the model is making predictions similar to the observed data.",
      x = pred_variable_name,
      y = "Density"
    )

  return(ppc_plot)

  # 2. Observed vs Predicted Plot
  # Get the mean prediction for each observation
  # pred_mean <- colMeans(extract_preds)
  #
  # # Create dataframe for prediction vs observed plot
  # obs_pred_df <- data.frame(
  #   Observed = observed_data,
  #   Predicted = pred_mean
  # )
  #
  # # Plot observed vs predicted with 1:1 line
  # obs_pred_plot <- ggplot(obs_pred_df, aes(x = Observed, y = Predicted)) +
  #   geom_point(alpha = 0.7) +
  #   geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  #   theme_minimal() +
  #   labs(
  #     title = "Observed vs. Predicted",
  #     subtitle = "Points closer to the line indicate better fit",
  #     x = "Observed Spawners",
  #     y = "Predicted Spawners"
  #   )
  #
  # gridExtra::grid.arrange(ppc_plot, obs_pred_plot)
  # ggsave(plot = ppc_plot, filename = paste0(local_folder, "/", stream, "_", "p2s_ppc_plot.png"), width = 10, height = 8)
  # ggsave(plot = obs_pred_plot, filename = paste0(local_folder, "/", stream, "_", "p2s_pred_plot.png"), width = 10, height = 8)
}

