store_model_fit <- function(con, storage_account, container_name, access_key, data, results_name, ...){

  model_board <- model_pin_board(storage_account, "model-results")

  blob_url <- pin_model_data(
    model_board,
    data,
    name = results_name,
    ...
  )

  total_rows <- insert_model_parameters(con, data, blob_url)
  message(glue::glue("Inserted {total_rows} into database. Uploaded results to {blob_url}."))

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
#' @export
setup_azure_blob_backend <- function(storage_account, access_key=Sys.getenv("AZ_CONTAINER_ACCESS_KEY")) {
  if (access_key == "") {
    stop("access key is required in order to write and read from the azure boad", call. = FALSE)
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
#' @export
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
                  metadata = pin_metadata
  )

  latest_version_df <- board |> pins::pin_versions(data_name)
  latest_version <- latest_version_df$version[1]

  data_url <- glue::glue("{board$container$endpoint$url}/{board$path}/{data_name}/{latest_version}/{data_name}.rds")

  return(data_url)
}

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

insert_model_parameters <- function(con, model, blob_url) {

  model_final_results <- model$final_results
  model_fit_filename <- stringr::str_extract(model$full_object$model.file, "[^/]+$")
  model_final_results$model_fit_filename <- model_fit_filename
  model_final_results$blob_url <- blob_url
  #TODO: how to extract model_name
  model_final_results$model_name <- "BTSPASX"

  model_final_results <- join_lookup(model_final_results, "model_name", "model_name", "definition", "model_name_id")
  model_final_results <- join_lookup(model_final_results, "model_location", "site", "site", "location_id")
  model_final_results <- join_lookup(model_final_results, "statistic", "statistic", "definition", "statistic_id")
  model_final_results <- join_lookup(model_final_results, "lifestage", "life_stage", "definition", "lifestage_id")
  model_final_results <- join_lookup(model_final_results, "parameter", "parameter", "definition", "parameter_id")

  model_final_results <-  model_final_results |>
    select(model_name_id, location_id, run_year, week_fit, lifestage_id, srjpedata_version, model_fit_filename, parameter_id, statistic_id, value, blob_url)

  query <- glue::glue_sql(
    "INSERT INTO model_parameters (
          model_name_id,
          location_id,
          run_year,
          week_fit,
          lifestage_id,
          srjpedata_version,
          model_fit_filename,
          parameter_id,
          statistic_id,
          value,
          blob_url
        ) VALUES (
          UNNEST(ARRAY[{model_final_results$model_name_id*}]),
          UNNEST(ARRAY[{model_final_results$location_id*}]),
          UNNEST(ARRAY[{model_final_results$run_year*}]),
          UNNEST(ARRAY[{model_final_results$week_fit*}]),
          UNNEST(ARRAY[{model_final_results$lifestage_id*}]),
          UNNEST(ARRAY[{model_final_results$srjpedata_version*}]),
          UNNEST(ARRAY[{model_final_results$model_fit_filename*}]),
          UNNEST(ARRAY[{model_final_results$parameter_id*}]),
          UNNEST(ARRAY[{model_final_results$statistic_id*}]),
          UNNEST(ARRAY[{model_final_results$value*}]),
          UNNEST(ARRAY[{model_final_results$blob_url*}])
        );",
    .con = con
  )
  res <- DBI::dbExecute(con, query)

  return(res)
}


