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

# ── Model name constants ───────────────────────────────────────────────────────
# .approved_model_names: must match the `name` column of the `model_name` table.
# .abundance_model_types / .pcap_one_site_model_types: the specific variants
#   stored in the `model_type` column of `model_run`. These map to a single
#   consolidated model_name in the DB ("abundance" and "pcap_one_site").
# Extend these vectors (and update model_name / model_type accordingly) when
# adding new model families.

.approved_model_names <- c(
  "abundance",
  "pcap_all_sites", "pcap_one_site",
  "p2s", "stock_recruit",
  "beta_dev_hbmrt", "beta_dv_hbmrt_lag1",
  "survival_cov_wy", "survival_no_cov"
)

# Abundance model type variants (stored in model_run.model_type)
.abundance_model_types <- c(
  "all_mark_recap", "no_mark_recap", "missing_mark_recap", "no_mark_recap_no_trib"
)

# pCap one-site model type variants (stored in model_run.model_type)
.pcap_one_site_model_types <- c("standard", "skew")

# Maps a model_inputs$model_name to the consolidated DB model_name
.model_name_lookup <- c(
  all_mark_recap        = "abundance",
  no_mark_recap         = "abundance",
  missing_mark_recap    = "abundance",
  no_mark_recap_no_trib = "abundance",
  pcap_one_site         = "pcap_one_site",
  pcap_one_site_skew    = "pcap_one_site"
)

# All valid values for model_inputs$model_name (inputs-level names, pre-consolidation)
.approved_input_model_names <- c(
  .abundance_model_types,
  "pcap_all_sites", "pcap_one_site", "pcap_one_site_skew",
  "p2s", "stock_recruit",
  "beta_dev_hbmrt", "beta_dv_hbmrt_lag1",
  "survival_cov_wy", "survival_no_cov"
)

# Model class lookup for validation
.bugs_input_models <- .abundance_model_types
.stan_input_models <- setdiff(.approved_input_model_names, .bugs_input_models)


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

  input_model_name <- str_to_lower(model_inputs$model_name)

  # ── Validate input model name ──────────────────────────────────────────────
  if (!input_model_name %in% .approved_input_model_names) {
    cli::cli_abort(c(
      "{.arg model_inputs$model_name} must be one of the approved input model names.",
      "i" = "Approved names: {.val {(.approved_input_model_names)}}",
      "x" = "Got: {.val {input_model_name}}"
    ))
  }

  # Resolve the consolidated DB model name (e.g. "all_mark_recap" -> "abundance")
  results_name <- unname(.model_name_lookup[input_model_name])
  if (is.na(results_name)) results_name <- input_model_name

  # ── Validate model object class ────────────────────────────────────────────
  if (input_model_name %in% .bugs_input_models && !inherits(model_fit_object, "bugs")) {
    cli::cli_abort(
      "{.val {input_model_name}} expects a {.cls bugs} object; got {.cls {class(model_fit_object)}}."
    )
  }
  if (input_model_name %in% .stan_input_models && !inherits(model_fit_object, "stanfit")) {
    cli::cli_abort(
      "{.val {input_model_name}} expects a {.cls stanfit} object; got {.cls {class(model_fit_object)}}."
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
#' Downloads a model fit object from Azure Blob Storage.
#'
#' By default the most recent version is returned. For model types where
#' multiple fits share the same `results_name` (e.g. `pcap_one_site` fits for
#' different sites, or abundance fits for different site × run_year
#' combinations), supply `con` along with the relevant filter arguments to
#' resolve the correct fit via the database:
#'
#' * **pCap one_site / one_site_skew** — supply `con` + `site_selection`
#' * **Abundance models** — supply `con` + `site` + `run_year`
#' * **pCap all_sites / other single fits** — no filters needed; omit `con`
#'
#' When `con` and filters are provided the database is queried for the most
#' recent matching `model_run` record and the blob URL stored there is used
#' to resolve the exact version. The `version` argument can still be used to
#' override this and pin down a specific historical fit.
#'
#' @param results_name The model type name used when the fit was stored. Must
#'   be one of `.approved_model_names` (e.g. `"pcap_one_site"`,
#'   `"all_mark_recap"`).
#' @param con Optional. A database connection object (e.g. from
#'   [DBI::dbConnect()]). Required when using `site`, `run_year`, or
#'   `site_selection` filters.
#' @param site Optional. Site name to filter on (abundance models).
#' @param run_year Optional. Run year to filter on (abundance models).
#' @param site_selection Optional. Site name to filter on for pCap one_site
#'   models (matches the `site_selection` column in `model_run`).
#' @param storage_account Azure storage account name. Defaults to
#'   `"jpemodelresults"`.
#' @param container_name Azure blob container name. Defaults to
#'   `"model-results"`.
#' @param access_key Azure storage access key with **read** permissions.
#'   Defaults to the `AZ_CONTAINER_ACCESS_KEY` environment variable.
#' @param version Optional. A specific version string as shown by
#'   [list_model_versions()]. Overrides DB-based resolution when supplied.
#'
#' @return The model fit object (class `stanfit` or `bugs`).
#'
#' @examples
#' \dontrun{
#' # pCap all-sites — no filters needed
#' pCap_fit <- get_model_fit("pcap_all_sites")
#'
#' # pCap one-site — resolve by site_selection via DB
#' tisdale_fit <- get_model_fit("pcap_one_site_skew", con = con,
#'                               site_selection = "tisdale")
#' kdl_fit     <- get_model_fit("pcap_one_site_skew", con = con,
#'                               site_selection = "knights landing")
#'
#' # Abundance — resolve by site + run_year via DB
#' abund_fit <- get_model_fit("all_mark_recap", con = con,
#'                             site = "ubc", run_year = 2022)
#'
#' # A specific historical version (bypasses DB lookup)
#' old_fit <- get_model_fit("pcap_all_sites", version = "20240315T102233Z-abc12")
#' }
#' @export
get_model_fit <- function(results_name,
                          con             = NULL,
                          site            = NULL,
                          run_year        = NULL,
                          site_selection  = NULL,
                          storage_account = "jpemodelresults",
                          container_name  = "model-results",
                          access_key      = Sys.getenv("AZ_CONTAINER_ACCESS_KEY"),
                          version         = NULL) {

  if (!results_name %in% .approved_model_names) {
    cli::cli_abort(c(
      "{.arg results_name} must be one of the approved DB model names.",
      "i" = "Approved names: {.val {(.approved_model_names)}}",
      "i" = "For abundance fits use {.val {'abundance'}}; for pCap one-site use {.val {'pcap_one_site'}}.",
      "x" = "Got: {.val {results_name}}"
    ))
  }

  # ── Validate filter usage ──────────────────────────────────────────────────
  using_filters <- !is.null(site) || !is.null(run_year) || !is.null(site_selection)
  if (using_filters && is.null(con)) {
    cli::cli_abort(
      "A database connection {.arg con} is required when using {.arg site}, \\
       {.arg run_year}, or {.arg site_selection} filters."
    )
  }

  board <- model_pin_board(storage_account, container_name, results_name,
                           access_key = access_key)

  # ── Resolve version ────────────────────────────────────────────────────────
  if (is.null(version)) {
    if (using_filters) {
      # Query DB for the most recent matching model_run record
      query <- dplyr::tbl(con, "model_run") |>
        dplyr::inner_join(dplyr::tbl(con, "model_name"),
                          by = c("model_name_id" = "id")) |>
        dplyr::filter(name == results_name)

      if (!is.null(site))           query <- dplyr::filter(query, site           == !!site)
      if (!is.null(run_year))       query <- dplyr::filter(query, run_year       == !!run_year)
      if (!is.null(site_selection)) query <- dplyr::filter(query, site_selection == !!site_selection)

      matched <- query |>
        dplyr::select(blob_fit_storage_url, created_at) |>
        dplyr::collect() |>
        dplyr::arrange(dplyr::desc(created_at)) |>
        dplyr::slice(1)

      if (nrow(matched) == 0) {
        cli::cli_abort(
          "No model run found in the database matching the supplied filters."
        )
      }

      # Parse version token from stored blob URL:
      # .../model-fits/<model_name>/<model_name>/<version>/<model_name>.rds
      version <- basename(dirname(matched$blob_fit_storage_url))
      cli::cli_alert_info("Resolved version from database: {.val {version}}")

    } else {
      # No filters — fall back to most recent pin version
      version <- board |>
        pins::pin_versions(results_name) |>
        dplyr::arrange(dplyr::desc(created)) |>
        dplyr::pull(version) |>
        dplyr::first()

      cli::cli_alert_info(
        "No filters supplied — downloading latest version: {.val {version}}"
      )
    }
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


# ── get_many_model_fits() ─────────────────────────────────────────────────────

#' @title Retrieve Multiple Model Fit Objects
#' @description
#' Downloads the most recent model fit for each site × run_year combination,
#' querying the database to identify which runs to fetch and then pulling each
#' `.rds` from Azure Blob Storage. Results are returned as a named list keyed
#' by `"<site>_<run_year>"`.
#'
#' All filter arguments are optional and combinable. With no filters the
#' function returns the most recent fit for every site × run_year combination
#' present in `model_run` for the given `model_name`.
#'
#' @param con A database connection object (e.g. from [DBI::dbConnect()]).
#' @param model_name The model type to retrieve. Must be one of
#'   `.approved_model_names` (e.g. `"all_mark_recap"`, `"no_mark_recap"`).
#' @param sites Optional character vector of sites to include (e.g.
#'   `c("ubc", "lcc")`). When `NULL` all sites are returned.
#' @param run_years Optional integer vector of run years to include (e.g.
#'   `2020:2024`). When `NULL` all run years are returned.
#' @param storage_account Azure storage account name. Defaults to
#'   `"jpemodelresults"`.
#' @param container_name Azure blob container name. Defaults to
#'   `"model-results"`.
#' @param access_key Azure storage access key with **read** permissions.
#'   Defaults to the `AZ_CONTAINER_ACCESS_KEY` environment variable.
#'
#' @return A named list of model fit objects. Names are formatted as
#'   `"<site>_<run_year>"` (e.g. `"ubc_2020"`, `"lcc_2021"`). Any site ×
#'   run_year that fails to download is returned as `NULL` with a warning
#'   rather than aborting the whole batch.
#'
#' @param model_name The consolidated DB model name to retrieve. Must be one of
#'   `.approved_model_names` (e.g. `"abundance"`, `"pcap_one_site"`). The
#'   specific model variant used for each fit is stored in the `model_type`
#'   column of `model_run` and is not used for filtering here.
#'
#' @examples
#' \dontrun{
#' # All abundance fits for every site × run_year
#' fits <- get_many_model_fits(con, model_name = "abundance")
#'
#' # Filter to specific sites or run years
#' fits <- get_many_model_fits(con, model_name = "abundance",
#'                             sites = c("ubc", "lcc", "mill creek"))
#'
#' fits <- get_many_model_fits(con, model_name = "abundance",
#'                             run_years = 2020:2024)
#'
#' # Access a single result from the list
#' fits[["ubc_2022"]]
#' }
#' @export
get_many_model_fits <- function(con,
                                model_name,
                                sites           = NULL,
                                run_years       = NULL,
                                storage_account = "jpemodelresults",
                                container_name  = "model-results",
                                access_key      = Sys.getenv("AZ_CONTAINER_ACCESS_KEY")) {

  # ── Validate model name ────────────────────────────────────────────────────
  if (!model_name %in% .approved_model_names) {
    cli::cli_abort(c(
      "{.arg model_name} must be one of the approved model names.",
      "i" = "Approved names: {.val {(.approved_model_names)}}",
      "x" = "Got: {.val {model_name}}"
    ))
  }

  # ── Query DB for the most recent run per site × run_year ───────────────────
  runs <- dplyr::tbl(con, "model_run") |>
    dplyr::inner_join(dplyr::tbl(con, "model_name"), by = c("model_name_id" = "id")) |>
    dplyr::filter(name == model_name) |>
    dplyr::filter(!is.na(site), !is.na(run_year))

  if (!is.null(sites)) {
    runs <- dplyr::filter(runs, site %in% !!sites)
  }
  if (!is.null(run_years)) {
    runs <- dplyr::filter(runs, run_year %in% !!run_years)
  }

  # Keep the most recent model_run record per site × run_year, across all
  # model types in the group (e.g. for a given site/year the best available
  # model type is used, as determined by created_at)
  runs <- runs |>
    dplyr::group_by(site, run_year) |>
    dplyr::slice_max(created_at, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(site, run_year, name, blob_fit_storage_url, created_at) |>
    dplyr::arrange(site, run_year) |>
    dplyr::collect()

  if (nrow(runs) == 0) {
    cli::cli_warn("No model runs found in the database matching the supplied filters.")
    return(list())
  }

  cli::cli_alert_info(
    "Downloading {nrow(runs)} fit{?s} for {.val {model_name}} from Azure Blob Storage."
  )

  # ── Download each fit from blob ────────────────────────────────────────────
  board <- model_pin_board(storage_account, container_name, model_name,
                           access_key = access_key)

  fits <- vector("list", nrow(runs))
  names(fits) <- paste0(runs$site, "_", runs$run_year)

  for (i in seq_len(nrow(runs))) {
    key <- names(fits)[i]

    fits[[key]] <- tryCatch({
      version <- basename(dirname(runs$blob_fit_storage_url[i]))
      pins::pin_read(board, model_name, version = version)
    }, error = function(e) {
      cli::cli_warn("Failed to download {.val {key}}: {e$message}")
      NULL
    })

    cli::cli_progress_message("  Downloaded {i}/{nrow(runs)}: {key}")
  }

  n_ok   <- sum(!vapply(fits, is.null, logical(1)))
  n_fail <- nrow(runs) - n_ok

  cli::cli_alert_success("Downloaded {n_ok}/{nrow(runs)} fit{?s} successfully.")
  if (n_fail > 0) {
    cli::cli_alert_warning("{n_fail} fit{?s} failed — returned as NULL in the list.")
  }

  fits
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

  # ── Resolve consolidated DB model_name and model_type ─────────────────────
  # input_model_name is the specific variant (e.g. "all_mark_recap",
  # "pcap_one_site_skew"). results_name is the consolidated DB name
  # (e.g. "abundance", "pcap_one_site"). model_type stores the variant.
  input_model_name <- str_to_lower(model_inputs$model_name)
  results_name     <- unname(.model_name_lookup[input_model_name])
  if (is.na(results_name)) results_name <- input_model_name
  model_type       <- if (input_model_name != results_name) input_model_name else NA_character_

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
    model_type           = model_type,
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
