
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




