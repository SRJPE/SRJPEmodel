# canonical.R
# Functions for designating and retrieving the canonical ("production") model
# run for a given model type + scope (site / run_year / site_selection,
# where applicable). Builds on top of storage.R: a canonical designation
# always points at a specific model_run row via model_run_id, so canonical
# fits are downloaded through the same get_model_fit()/version machinery as
# any other stored fit — this file only adds the "which version is
# production" bookkeeping.
#
# canonical_model_run (see db/migrations/0001_create_canonical_model_run.sql)
# is an append-only log, not an update-in-place table: set_canonical_model_run()
# always inserts a new row rather than overwriting an existing one. The
# "current" canonical version for a scope is simply the row with the latest
# set_at for that scope — mirroring the newest-wins pattern store_model_fit()
# already uses for pin versions. This gives a free audit trail with no
# separate is_current flag or transactional flip.
#
# Deliberately not part of the scope: model_type (the mark-recapture variant
# for bt-spas-x, e.g. "all_mark_recap"). A canonical designation answers
# "what's the production run for this site/run_year", not "...for this
# specific variant" — the variant is an implementation detail of how a given
# model_run was fit, not something that should further subdivide which run
# counts as canonical. It's still shown in list_canonical_model_runs()'s
# output (joined from model_run), just not part of the lookup key.
#
# ── Public API ────────────────────────────────────────────────────────────────
#   set_canonical_model_run()      Record a model_run as canonical for a scope.
#   list_canonical_model_runs()    Current canonical designation per scope.
#   get_canonical_history()        Full audit trail for one scope.
#   get_canonical_model_fit()      Download the current canonical fit for a scope.
#   get_many_canonical_model_fits() Download canonical fits for many site × run_year scopes.
# ─────────────────────────────────────────────────────────────────────────────


# ── set_canonical_model_run() ───────────────────────────────────────────────

#' @title Mark a Model Run as Canonical
#' @description
#' Records `model_run_id` as the canonical ("production") fit for a given
#' model type + scope (site / run_year / site_selection, where applicable).
#' Does not move, copy, or modify the underlying fit — it only inserts an
#' audit row pointing at it. Because `canonical_model_run` is append-only,
#' calling this again for the same scope creates a new "current" designation
#' and preserves the previous one in history (see [get_canonical_history()]).
#'
#' @param con A database connection object.
#' @param results_name The consolidated model name (e.g. `"bt-spas-x"`,
#'   `"pcap_one_site"`). Must be one of `.approved_model_names`.
#' @param model_run_id The `model_run.id` of the fit to mark canonical.
#' @param set_by Who is making this designation (e.g. an email address).
#' @param site Optional. Site, for bt-spas-x scopes.
#' @param run_year Optional. Run year, for bt-spas-x scopes.
#' @param site_selection Optional. Site selection, for pcap_one_site scopes.
#' @param note Optional free-text note (e.g. why this version was picked).
#'
#' @return Invisibly returns the number of rows inserted (always 1 on success).
#'
#' @examples
#' \dontrun{
#' set_canonical_model_run(
#'   con = con, results_name = "bt-spas-x", model_run_id = 482,
#'   set_by = "avizek@flowwest.com", site = "ubc", run_year = 2024,
#'   note = "Best convergence diagnostics after QA review"
#' )
#' }
#' @export
set_canonical_model_run <- function(con, results_name, model_run_id, set_by,
                                     site = NULL, run_year = NULL,
                                     site_selection = NULL, note = NULL) {

  if (!results_name %in% .approved_model_names) {
    cli::cli_abort(c(
      "{.arg results_name} must be one of the approved DB model names.",
      "i" = "Approved names: {.val {(.approved_model_names)}}",
      "x" = "Got: {.val {results_name}}"
    ))
  }

  model_name_id <- dplyr::tbl(con, "model_name") |>
    dplyr::filter(name == results_name) |>
    dplyr::pull(id)

  if (length(model_name_id) == 0) {
    cli::cli_abort(
      "No matching row found in the {.val model_name} table for {.val {results_name}}."
    )
  }

  run_model_name_id <- dplyr::tbl(con, "model_run") |>
    dplyr::filter(id == !!model_run_id) |>
    dplyr::pull(model_name_id)

  if (length(run_model_name_id) == 0) {
    cli::cli_abort("No {.val model_run} row found with id {.val {model_run_id}}.")
  }
  if (run_model_name_id != model_name_id) {
    cli::cli_abort(
      "model_run {.val {model_run_id}} belongs to a different model type than {.val {results_name}}."
    )
  }

  new_row <- data.frame(
    model_name_id  = model_name_id,
    site           = site %||% NA_character_,
    run_year       = run_year %||% NA_integer_,
    site_selection = site_selection %||% NA_character_,
    model_run_id   = model_run_id,
    set_by         = set_by,
    note           = note %||% NA_character_,
    set_at         = Sys.time(),
    stringsAsFactors = FALSE
  )

  rows_inserted <- DBI::dbAppendTable(con, "canonical_model_run", new_row)

  cli::cli_alert_success(
    "Marked model_run {.val {model_run_id}} as canonical for {.val {results_name}}."
  )
  invisible(rows_inserted)
}


# ── list_canonical_model_runs() ─────────────────────────────────────────────

#' @title List Current Canonical Model Runs
#' @description
#' Returns the current canonical designation (most recent `set_at`) for each
#' distinct model type + scope, joined to the underlying `model_run` /
#' `model_name` for display.
#'
#' @param con A database connection object.
#' @param results_name Optional. Filter to one consolidated model name.
#'
#' @return A tibble with one row per scope: `model_name`, `model_type`
#'   (the variant of the canonical run, informational only — not part of the
#'   scope), `site`, `run_year`, `site_selection`, `model_run_id`,
#'   `description`, `blob_fit_storage_url`, `set_by`, `set_at`, `note`.
#'
#' @examples
#' \dontrun{
#' list_canonical_model_runs(con)
#' list_canonical_model_runs(con, results_name = "bt-spas-x")
#' }
#' @export
list_canonical_model_runs <- function(con, results_name = NULL) {

  scope_cols <- c("model_name_id", "site", "run_year", "site_selection")

  canonical <- dplyr::tbl(con, "canonical_model_run")

  if (!is.null(results_name)) {
    model_name_id <- dplyr::tbl(con, "model_name") |>
      dplyr::filter(name == results_name) |>
      dplyr::pull(id)
    canonical <- dplyr::filter(canonical, model_name_id == !!model_name_id)
  }

  canonical |>
    dplyr::group_by(dplyr::across(dplyr::all_of(scope_cols))) |>
    dplyr::slice_max(set_at, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::inner_join(
      dplyr::tbl(con, "model_name") |> dplyr::select(id, name),
      by = c("model_name_id" = "id")
    ) |>
    dplyr::rename(model_name = name) |>
    dplyr::inner_join(
      dplyr::tbl(con, "model_run") |>
        dplyr::select(id, model_type, description, blob_fit_storage_url),
      by = c("model_run_id" = "id")
    ) |>
    dplyr::select(model_name, model_type, site, run_year, site_selection,
                   model_run_id, description, blob_fit_storage_url,
                   set_by, set_at, note) |>
    dplyr::arrange(model_name, site, run_year) |>
    dplyr::collect()
}


# ── get_canonical_history() ─────────────────────────────────────────────────

#' @title History of Canonical Designations for One Scope
#' @description
#' Returns every canonical designation ever recorded for one model type +
#' scope, newest first — the audit trail [set_canonical_model_run()] builds
#' by always inserting rather than updating.
#'
#' @inheritParams set_canonical_model_run
#' @return A tibble with columns `model_run_id`, `set_by`, `set_at`, `note`,
#'   ordered newest first.
#'
#' @examples
#' \dontrun{
#' get_canonical_history(con, "bt-spas-x", site = "ubc", run_year = 2024)
#' }
#' @export
get_canonical_history <- function(con, results_name, site = NULL,
                                   run_year = NULL, site_selection = NULL) {

  model_name_id <- dplyr::tbl(con, "model_name") |>
    dplyr::filter(name == results_name) |>
    dplyr::pull(id)

  if (length(model_name_id) == 0) {
    cli::cli_abort(
      "No matching row found in the {.val model_name} table for {.val {results_name}}."
    )
  }

  q <- dplyr::tbl(con, "canonical_model_run") |>
    dplyr::filter(model_name_id == !!model_name_id)

  q <- if (is.null(site)) dplyr::filter(q, is.na(site)) else dplyr::filter(q, site == !!site)
  q <- if (is.null(run_year)) dplyr::filter(q, is.na(run_year)) else dplyr::filter(q, run_year == !!run_year)
  q <- if (is.null(site_selection)) {
    dplyr::filter(q, is.na(site_selection))
  } else {
    dplyr::filter(q, site_selection == !!site_selection)
  }

  q |>
    dplyr::select(model_run_id, set_by, set_at, note) |>
    dplyr::arrange(dplyr::desc(set_at)) |>
    dplyr::collect()
}


# ── get_canonical_model_fit() ───────────────────────────────────────────────

#' @title Retrieve the Canonical Model Fit
#' @description
#' Downloads the current canonical fit for a model type + scope, resolved via
#' [list_canonical_model_runs()] rather than "most recent version" — use this
#' instead of [get_model_fit()] wherever "the production run" is wanted
#' rather than "whatever was uploaded last".
#'
#' @inheritParams set_canonical_model_run
#' @inheritParams get_model_fit
#' @return The model fit object (class `stanfit` or `bugs`).
#'
#' @examples
#' \dontrun{
#' get_canonical_model_fit(con, "bt-spas-x", site = "ubc", run_year = 2024)
#' }
#' @export
get_canonical_model_fit <- function(con, results_name, site = NULL,
                                     run_year = NULL, site_selection = NULL,
                                     storage_account = "jpemodelresults",
                                     container_name  = "model-results",
                                     access_key = Sys.getenv("AZ_CONTAINER_ACCESS_KEY")) {

  canonical <- list_canonical_model_runs(con, results_name = results_name)

  if (!is.null(site))           canonical <- dplyr::filter(canonical, site == !!site)
  if (!is.null(run_year))       canonical <- dplyr::filter(canonical, run_year == !!run_year)
  if (!is.null(site_selection)) canonical <- dplyr::filter(canonical, site_selection == !!site_selection)

  if (nrow(canonical) == 0) {
    cli::cli_abort(
      "No canonical model run found for {.val {results_name}} matching the supplied filters."
    )
  }
  if (nrow(canonical) > 1) {
    cli::cli_abort(c(
      "Multiple canonical runs matched {.val {results_name}}.",
      "i" = "Supply more filters ({.arg site}, {.arg run_year}, {.arg site_selection}) to narrow to one."
    ))
  }

  version <- .version_from_blob_url(canonical$blob_fit_storage_url)

  get_model_fit(
    results_name    = results_name,
    version         = version,
    storage_account = storage_account,
    container_name  = container_name,
    access_key      = access_key
  )
}


#' Parses the pins version token out of a stored blob URL:
#' .../model-fits/<model_name>/<model_name>/<version>/<model_name>.rds
#' @keywords internal
.version_from_blob_url <- function(url) {
  basename(dirname(url))
}


# ── get_many_canonical_model_fits() ─────────────────────────────────────────

#' @title Retrieve Multiple Canonical Model Fit Objects
#' @description
#' Downloads the current canonical fit for every scope matching `model_name`
#' (and any supplied filters), resolved via [list_canonical_model_runs()]
#' rather than "most recent version" — the batch analog of
#' [get_canonical_model_fit()]. Unlike [get_many_model_fits()] (which is
#' specific to the site × run_year scope of abundance models), this works
#' across all scope shapes used in `canonical_model_run`:
#'
#' * **bt-spas-x / plad_btspasx_results** — scoped by site × run_year (PLAD
#'   fits omit run_year and are scoped by site alone).
#' * **pcap_one_site** — scoped by `site_selection` (e.g. `"tisdale"`,
#'   `"knights landing"`).
#' * **pcap_all_sites / other unscoped model types** — a single canonical fit
#'   with no site/run_year/site_selection scope.
#'
#' A single call only spans scopes for one `model_name` at a time — e.g. call
#' once with `"bt-spas-x"` and once with `"pcap_one_site"` to get both
#' families.
#'
#' @param con A database connection object (e.g. from [DBI::dbConnect()]).
#' @param model_name The consolidated DB model name to retrieve. Must be one of
#'   `.approved_model_names` (e.g. `"bt-spas-x"`, `"pcap_one_site"`).
#' @param sites Optional character vector of sites to include (e.g.
#'   `c("ubc", "lcc")`). Applies to site-scoped model types (bt-spas-x,
#'   plad_btspasx_results). When `NULL` all sites are returned.
#' @param run_years Optional integer vector of run years to include (e.g.
#'   `2020:2024`). Applies to bt-spas-x. When `NULL` all run years are
#'   returned.
#' @param site_selections Optional character vector of `site_selection`
#'   values to include (e.g. `c("tisdale", "knights landing")`). Applies to
#'   pcap_one_site. When `NULL` all site selections are returned.
#' @param storage_account Azure storage account name. Defaults to
#'   `"jpemodelresults"`.
#' @param container_name Azure blob container name. Defaults to
#'   `"model-results"`.
#' @param access_key Azure storage access key with **read** permissions.
#'   Defaults to the `AZ_CONTAINER_ACCESS_KEY` environment variable.
#'
#' @return A named list of model fit objects. Names depend on the scope of
#'   `model_name`: `"<site>_<run_year>"` for bt-spas-x (e.g. `"ubc_2020"`),
#'   `"<site>"` for site-only scopes like plad_btspasx_results,
#'   `"<site_selection>"` for pcap_one_site (e.g. `"tisdale"`), or
#'   `model_name` itself for unscoped model types like pcap_all_sites. Any
#'   scope that fails to download is returned as `NULL` with a warning rather
#'   than aborting the whole batch.
#'
#' @examples
#' \dontrun{
#' # All canonical abundance (BT-SPAS-X) fits for every site × run_year
#' fits <- get_many_canonical_model_fits(con, model_name = "bt-spas-x")
#'
#' # Filter to specific sites or run years
#' fits <- get_many_canonical_model_fits(con, model_name = "bt-spas-x",
#'                                       sites = c("ubc", "lcc", "mill creek"))
#'
#' fits <- get_many_canonical_model_fits(con, model_name = "bt-spas-x",
#'                                       run_years = 2020:2024)
#'
#' # Canonical pCap one-site fits for tisdale and knights landing
#' fits <- get_many_canonical_model_fits(con, model_name = "pcap_one_site",
#'                                       site_selections = c("tisdale", "knights landing"))
#'
#' # The single canonical pCap all-sites fit
#' fits <- get_many_canonical_model_fits(con, model_name = "pcap_all_sites")
#' }
#' @export
get_many_canonical_model_fits <- function(con,
                                          model_name,
                                          sites           = NULL,
                                          run_years       = NULL,
                                          site_selections = NULL,
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

  # ── Resolve the current canonical designation for every matching scope ─────
  canonical <- list_canonical_model_runs(con, results_name = model_name)

  if (!is.null(sites))           canonical <- dplyr::filter(canonical, site %in% sites)
  if (!is.null(run_years))       canonical <- dplyr::filter(canonical, run_year %in% run_years)
  if (!is.null(site_selections)) canonical <- dplyr::filter(canonical, site_selection %in% site_selections)

  if (nrow(canonical) == 0) {
    cli::cli_warn("No canonical model runs found matching the supplied filters.")
    return(list())
  }

  cli::cli_alert_info(
    "Downloading {nrow(canonical)} canonical fit{?s} for {.val {model_name}} from Azure Blob Storage."
  )

  # ── Download each fit from blob ────────────────────────────────────────────
  board <- model_pin_board(storage_account, container_name, model_name,
                           access_key = access_key)

  keys <- .canonical_fit_keys(canonical, model_name)

  fits <- vector("list", nrow(canonical))
  names(fits) <- keys

  for (i in seq_len(nrow(canonical))) {
    key <- names(fits)[i]

    fits[[key]] <- tryCatch({
      version <- .version_from_blob_url(canonical$blob_fit_storage_url[i])
      pins::pin_read(board, model_name, version = version)
    }, error = function(e) {
      cli::cli_warn("Failed to download {.val {key}}: {e$message}")
      NULL
    })

    cli::cli_progress_message("  Downloaded {i}/{nrow(canonical)}: {key}")
  }

  n_ok   <- sum(!vapply(fits, is.null, logical(1)))
  n_fail <- nrow(canonical) - n_ok

  cli::cli_alert_success("Downloaded {n_ok}/{nrow(canonical)} canonical fit{?s} successfully.")
  if (n_fail > 0) {
    cli::cli_alert_warning("{n_fail} fit{?s} failed — returned as NULL in the list.")
  }

  fits
}


#' Builds a result-list key per row of a `list_canonical_model_runs()` tibble,
#' based on whichever scope columns are actually populated for that row (they
#' vary by model type): `"<site>_<run_year>"` for bt-spas-x, `"<site>"` for
#' site-only scopes (e.g. plad_btspasx_results), `"<site_selection>"` for
#' pcap_one_site, and `model_name` itself when no scope column is populated
#' (e.g. pcap_all_sites).
#' @keywords internal
.canonical_fit_keys <- function(canonical, model_name) {
  dplyr::case_when(
    !is.na(canonical$site) & !is.na(canonical$run_year) ~
      paste0(canonical$site, "_", canonical$run_year),
    !is.na(canonical$site_selection) ~ canonical$site_selection,
    !is.na(canonical$site)           ~ canonical$site,
    TRUE                             ~ model_name
  )
}
