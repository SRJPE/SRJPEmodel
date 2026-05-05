# storage.R
# Functions for storing and retrieving model fit objects (.rds) in Azure Blob Storage,
# and writing model run metadata to the JPE database.
#
# ── Public API ────────────────────────────────────────────────────────────────
#   store_model_fit()      Upload a model fit object and record it in the DB.
#   get_model_fit()        Download the most recent (or a specific) model fit.
#   list_model_versions()  List all stored versions for a given model type.
# ── Internal helpers ──────────────────────────────────────────────────────────
#   setup_azure_blob_backend()   Build an AzureStor endpoint object.
#   model_pin_board()            Build a pins board scoped to one model type.
#   insert_model_run()           Write one model run record to the database.
# ─────────────────────────────────────────────────────────────────────────────

# Approved model name values — must match the `name` column of the `model_name`
# table exactly. Extend this vector (and add a row to model_name) when adding
# new model families.
.approved_model_names <- c(
  # juvenile abundance (BUGS)
  "all_mark_recap", "no_mark_recap", "missing_mark_recap", "no_mark_recap_no_trib",
  # pCap (Stan)
  "pcap_all_sites", "pcap_one_site", "pcap_one_site_skew",
  # other model families
  "p2s", "stock_recruit",
  "beta_dev_hbmrt", "beta_dv_hbmrt_lag1",
  "survival_cov_wy", "survival_no_cov"
)

.bugs_models <- c(
  "all_mark_recap", "no_mark_recap", "missing_mark_recap", "no_mark_recap_no_trib"
)

.stan_models <- setdiff(.approved_model_names, .bugs_models)


# ── store_model_fit() ─────────────────────────────────────────────────────────

#' @title Store a Model Fit Object
#' @description
#' Uploads a model fit object as an `.rds` file to Azure Blob Storage and
#' inserts a corresponding record into the `model_run` database table.
#' Each call creates a new *version* of the pin, so all historical fits are
#' preserved and the latest is always the default when retrieving.
#'
#' Metadata is extracted automatically from `model_inputs`, so the same
#' function works for pCap and abundance models without extra arguments.
#' Fields that do not apply to a given model type (e.g. `skew` for abundance
#' models) are stored as `NULL` in the database.
#'
#' The `model_inputs` list is expected to carry the following slots depending
#' on model type:
#'
#' | Slot | DB column | Model type |
#' |---|---|---|
#' | `$model_name` | — | All (required) |
#' | `$site` | `site` | Abundance |
#' | `$run_year` | `run_year` | Abundance |
#' | `$skew` | `skew` | pCap one_site |
#' | `$site_selection` | `site_selection` | pCap one_site |
#' | `$sites_dropped` | `sites_excluded` | pCap all_sites |
#' | `$srjpedata_version` | `srjpedata_version` | All (optional override) |
#'
#' If `$srjpedata_version` is not present on `model_inputs`, the installed
#' version of `SRJPEdata` is used automatically.
#'
#' @param con A database connection object (e.g. from [DBI::dbConnect()]).
#' @param model_fit_object The fitted model object to store. Must be of class
#'   `stanfit` (Stan models) or `bugs` (WinBUGS / OpenBUGS models).
#' @param model_inputs The inputs list returned by the relevant
#'   `prepare_*_inputs()` function.
#' @param description A human-readable description of this model run (e.g.
#'   `"pCap all-sites fit, data through run-year 2024"`).
#' @param storage_account Azure storage account name. Defaults to
#'   `"jpemodelresults"`.
#' @param container_name Azure blob container name. Defaults to
#'   `"model-results"`.
#' @param access_key Azure storage access key with **write** permissions.
#'   Defaults to the `AZ_CONTAINER_ACCESS_KEY` environment variable.
#'
#' @return Invisibly returns the blob URL of the stored `.rds` file.
#'
#' @examples
#' \dontrun{
#' con <- DBI::dbConnect(RPostgres::Postgres(),
#'                       dbname = cfg$db_name, host = cfg$db_host,
#'                       port = cfg$db_port, user = cfg$db_user,
#'                       password = cfg$db_password)
#'
#' # --- pCap all-sites (Stan) ---
#' pCap_inputs <- prepare_pCap_inputs(model_type = "all_sites")
#' pCap_fit    <- fit_pCap_model(pCap_inputs)
#'
#' store_model_fit(
#'   con              = con,
#'   model_fit_object = pCap_fit,
#'   model_inputs     = pCap_inputs,
#'   description      = "pCap all-sites fit, run-year 2024"
#' )
#'
#' # --- pCap one-site with skew (Stan) ---
#' pCap_inputs <- prepare_pCap_inputs(model_type = "one_site",
#'                                    skew = TRUE, site_selection = "ubc")
#' pCap_fit    <- fit_pCap_model(pCap_inputs)
#'
#' store_model_fit(
#'   con              = con,
#'   model_fit_object = pCap_fit,
#'   model_inputs     = pCap_inputs,
#'   description      = "pCap one-site skew fit for ubc, run-year 2024"
#' )
#'
#' # --- Abundance (BUGS) ---
#' abund_inputs <- prepare_abundance_inputs(site = "ubc", run_year = 2024,
#'                                          pCap_model_type = "all_sites",
#'                                          pCap_model_object = pCap_fit)
#' abund_fit    <- fit_abundance_model_BUGS(abund_inputs, bugs_model_file, bugs_dir)
#'
#' store_model_fit(
#'   con              = con,
#'   model_fit_object = abund_fit,
#'   model_inputs     = abund_inputs,
#'   description      = "Abundance fit – ubc 2024"
#' )
#'
#' DBI::dbDisconnect(con)
#' }
#' @export
store_model_fit <- function(con,
                            model_fit_object,
                            model_inputs,
                            description     = NULL,
                            storage_account = "jpemodelresults",
                            container_name  = "model-results",
                            access_key      = Sys.getenv("AZ_CONTAINER_ACCESS_KEY")) {

  results_name <- str_to_lower(model_inputs$model_name)

  # ── Validate results_name ──────────────────────────────────────────────────
  if (!results_name %in% .approved_model_names) {
    cli::cli_abort(c(
      "{.arg model_inputs$model_name} must be one of the approved model names.",
      "i" = "Approved names: {.val {(.approved_model_names)}}",
      "x" = "Got: {.val {results_name}}"
    ))
  }

  # ── Validate model object class ────────────────────────────────────────────
  if (results_name %in% .bugs_models && !inherits(model_fit_object, "bugs")) {
    cli::cli_abort(
      "{.val {results_name}} expects a {.cls bugs} object; got {.cls {class(model_fit_object)}}."
    )
  }
  if (results_name %in% .stan_models && !inherits(model_fit_object, "stanfit")) {
    cli::cli_abort(
      "{.val {results_name}} expects a {.cls stanfit} object; got {.cls {class(model_fit_object)}}."
    )
  }

  # ── Upload to Azure Blob Storage ───────────────────────────────────────────
  board <- model_pin_board(storage_account, container_name, results_name,
                           access_key = access_key)

  pins::pin_write(
    board,
    model_fit_object,
    name        = results_name,
    type        = "rds",
    description = description
  )

  latest_version <- board |>
    pins::pin_versions(results_name) |>
    dplyr::arrange(dplyr::desc(created)) |>
    dplyr::pull(version) |>
    dplyr::first()

  blob_url <- glue::glue(
    "{board$container$endpoint$url}/{board$path}/{results_name}",
    "/{latest_version}/{results_name}.rds"
  )

  cli::cli_alert_success(
    "Stored {.val {results_name}} (version {.val {latest_version}}) to Azure Blob Storage."
  )

  # ── Write metadata record to database ─────────────────────────────────────
  tryCatch({
    insert_model_run(con, blob_url, results_name, description, model_inputs)
    cli::cli_alert_success("Inserted model run record into database.")
  }, error = function(e) {
    cli::cli_alert_danger(
      "Model fit uploaded to blob successfully, but DB insert failed: {e$message}"
    )
    cli::cli_alert_info(
      "The fit is safe in blob storage. Blob URL: {blob_url}"
    )
    stop(e)
  })

  invisible(blob_url)
}


# ── get_model_fit() ───────────────────────────────────────────────────────────

#' @title Retrieve a Model Fit Object
#' @description
#' Downloads a model fit object from Azure Blob Storage. By default the
#' *most recent* version is returned; supply `version` to retrieve a specific
#' historical fit (use [list_model_versions()] to browse available versions).
#'
#' @param results_name The model type name used when the fit was stored. Must
#'   be one of `.approved_model_names` (e.g. `"pcap_all_sites"`,
#'   `"all_mark_recap"`).
#' @param storage_account Azure storage account name. Defaults to
#'   `"jpemodelresults"`.
#' @param container_name Azure blob container name. Defaults to
#'   `"model-results"`.
#' @param access_key Azure storage access key with **read** permissions.
#'   Defaults to the `AZ_CONTAINER_ACCESS_KEY` environment variable.
#' @param version Optional. A specific version string as shown by
#'   [list_model_versions()]. When `NULL` (default) the latest version is
#'   returned.
#'
#' @return The model fit object (class `stanfit` or `bugs`).
#'
#' @examples
#' \dontrun{
#' # Latest pCap all-sites fit
#' pCap_fit  <- get_model_fit("pcap_all_sites")
#'
#' # Latest abundance fit
#' abund_fit <- get_model_fit("all_mark_recap")
#'
#' # A specific historical version
#' old_fit   <- get_model_fit("pcap_all_sites", version = "20240315T102233Z-abc12")
#' }
#' @export
get_model_fit <- function(results_name,
                          storage_account = "jpemodelresults",
                          container_name  = "model-results",
                          access_key      = Sys.getenv("AZ_CONTAINER_ACCESS_KEY"),
                          version         = NULL) {

  if (!results_name %in% .approved_model_names) {
    cli::cli_abort(c(
      "{.arg results_name} must be one of the approved model names.",
      "i" = "Approved names: {.val {(.approved_model_names)}}",
      "x" = "Got: {.val {results_name}}"
    ))
  }

  board <- model_pin_board(storage_account, container_name, results_name,
                           access_key = access_key)

  if (is.null(version)) {
    version <- board |>
      pins::pin_versions(results_name) |>
      dplyr::arrange(dplyr::desc(created)) |>
      dplyr::pull(version) |>
      dplyr::first()

    cli::cli_alert_info(
      "No version specified — downloading latest: {.val {version}}"
    )
  }

  model_object <- pins::pin_read(board, results_name, version = version)

  cli::cli_alert_success(
    "Retrieved {.val {results_name}} version {.val {version}} from Azure Blob Storage."
  )

  model_object
}


# ── list_model_versions() ─────────────────────────────────────────────────────

#' @title List Available Versions of a Model Fit
#' @description
#' Returns a data frame of all stored versions for a given model type, newest
#' first. Use this to inspect history or to supply a `version` argument to
#' [get_model_fit()].
#'
#' @inheritParams get_model_fit
#'
#' @return A tibble with columns `version`, `created`, and `hash`.
#'
#' @examples
#' \dontrun{
#' list_model_versions("pcap_all_sites")
#' list_model_versions("all_mark_recap")
#' }
#' @export
list_model_versions <- function(results_name,
                                storage_account = "jpemodelresults",
                                container_name  = "model-results",
                                access_key      = Sys.getenv("AZ_CONTAINER_ACCESS_KEY")) {

  if (!results_name %in% .approved_model_names) {
    cli::cli_abort(c(
      "{.arg results_name} must be one of the approved model names.",
      "i" = "Approved names: {.val {(.approved_model_names)}}",
      "x" = "Got: {.val {results_name}}"
    ))
  }

  board <- model_pin_board(storage_account, container_name, results_name,
                           access_key = access_key)

  board |>
    pins::pin_versions(results_name) |>
    dplyr::arrange(dplyr::desc(created))
}


# ── Internal helpers ──────────────────────────────────────────────────────────

#' @title Insert a Model Run Record
#' @description
#' Writes one row to the `model_run` table capturing the blob URL and all
#' model metadata extracted from the `model_inputs` list. Fields that do not
#' apply to a given model type are stored as `NULL`.
#'
#' Current `model_run` schema:
#' ```
#'  id                    SERIAL PRIMARY KEY
#'  blob_fit_storage_url  TEXT
#'  model_name_id         INTEGER REFERENCES model_name(id)
#'  updated_at            TIMESTAMP WITHOUT TIME ZONE NOT NULL
#'  srjpedata_version     CHARACTER VARYING
#'  description           TEXT
#'  site                  TEXT        -- abundance models only
#'  run_year              INTEGER     -- abundance models only
#'  skew                  BOOLEAN     -- pcap one_site only
#'  site_selection        TEXT        -- pcap one_site only
#'  sites_excluded        TEXT[]      -- pcap all_sites only
#'  created_at            TIMESTAMPTZ
#' ```
#'
#' @param con A database connection object.
#' @param blob_url Character. Full Azure Blob URL for the stored `.rds`.
#' @param results_name Character. Model type name (must be in `.approved_model_names`).
#' @param description Character or NULL. Human-readable run description.
#' @param model_inputs List. The inputs object from `prepare_*_inputs()`.
#'
#' @return Invisibly returns the number of rows inserted (always 1 on success).
#' @keywords internal
insert_model_run <- function(con, blob_url, results_name, description, model_inputs) {

  # ── Resolve model_name_id from lookup table ────────────────────────────────
  model_name_id <- dplyr::tbl(con, "model_name") |>
    dplyr::filter(name == results_name) |>
    dplyr::pull(id)

  if (length(model_name_id) == 0) {
    cli::cli_abort(c(
      "No matching row found in the {.val model_name} table for {.val {results_name}}.",
      "i" = "Add this model name to the lookup table before storing fits."
    ))
  }

  # ── srjpedata_version ──────────────────────────────────────────────────────
  # Prefer an explicit field on model_inputs; fall back to installed package version.
  srjpedata_version <- if (!is.null(model_inputs$srjpedata_version)) {
    as.character(model_inputs$srjpedata_version)
  } else {
    tryCatch(
      as.character(utils::packageVersion("SRJPEdata")),
      error = function(e) NA_character_
    )
  }

  # ── Model-specific metadata (all nullable) ─────────────────────────────────
  # Abundance models
  site     <- model_inputs$site     %||% NA_character_
  run_year <- model_inputs$run_year %||% NA_integer_

  # pCap one_site models
  skew           <- model_inputs$skew           %||% NA
  site_selection <- model_inputs$site_selection %||% NA_character_

  # pCap all_sites: sites excluded from the hyper-distribution.
  # prepare_pCap_inputs() stores this as $sites_dropped (mirrors exclude_from_hyper arg).
  # Stored as a Postgres TEXT[] — NULL if nothing was excluded.
  sites_excluded <- model_inputs$sites_dropped

  # ── INSERT ─────────────────────────────────────────────────────────────────
  # Use dbAppendTable() so RPostgres handles type serialisation natively —
  # correctly maps NULL, BOOLEAN, INTEGER, and TEXT[] without glue_sql
  # vector-length issues. sites_excluded is wrapped in list() so a character
  # vector is stored as a single array value rather than expanding into rows.
  new_row <- data.frame(
    blob_fit_storage_url = blob_url,
    model_name_id        = model_name_id,
    srjpedata_version    = srjpedata_version,
    description          = description %||% NA_character_,
    site                 = site,
    run_year             = run_year,
    skew                 = skew,
    site_selection       = site_selection,
    updated_at           = Sys.time(),
    stringsAsFactors     = FALSE
  )
  new_row$sites_excluded <- list(sites_excluded)

  rows_inserted <- DBI::dbAppendTable(con, "model_run", new_row)
  invisible(rows_inserted)
}


#' Null-coalescing operator
#' @keywords internal
`%||%` <- function(x, y) if (!is.null(x)) x else y


#' @title Azure Blob Backend
#' @description Creates an `AzureStor` storage endpoint object.
#'
#' @param storage_account Azure storage account name.
#' @param access_key Azure storage access key.
#'
#' @return An `AzureStor` storage endpoint object.
#' @keywords internal
setup_azure_blob_backend <- function(storage_account,
                                     access_key = Sys.getenv("AZ_CONTAINER_ACCESS_KEY")) {
  if (access_key == "") {
    cli::cli_abort(
      "An access key is required to read/write from Azure Blob Storage. \\
       Set {.envvar AZ_CONTAINER_ACCESS_KEY} in your R environment."
    )
  }
  storage_endpoint <- glue::glue("https://{storage_account}.blob.core.windows.net")
  AzureStor::storage_endpoint(storage_endpoint, key = access_key)
}


#' @title Model Pin Board
#' @description
#' Creates a `pins` board scoped to a specific model type's subfolder within
#' the Azure blob container (`model-fits/<model_name>/`).
#'
#' @param storage_account Azure storage account name.
#' @param container Azure blob container name.
#' @param model_name Model type name, used as the blob subfolder.
#' @param access_key Azure storage access key.
#'
#' @return A `pins` board object connected to Azure Blob Storage.
#' @keywords internal
model_pin_board <- function(storage_account, container, model_name,
                            access_key = Sys.getenv("AZ_CONTAINER_ACCESS_KEY")) {

  storage_client <- setup_azure_blob_backend(storage_account, access_key)

  container_names <- names(AzureStor::list_blob_containers(storage_client))
  if (!container %in% container_names) {
    cli::cli_abort(
      "Container {.val {container}} does not exist in storage account {.val {storage_account}}."
    )
  }

  blob_client <- AzureStor::blob_container(storage_client, container)
  pins::board_azure(container = blob_client, paste0("model-fits/", model_name))
}
