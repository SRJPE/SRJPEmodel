# test_storage_workflow.R
# Minimal script to test the upload and download functions for a pCap all-sites
# fit and a single abundance model fit.

library(tidyverse)
library(SRJPEdata)
library(SRJPEmodel)
library(rstan)
library(R2WinBUGS)

# ── STAN options ───────────────────────────────────────────────────────────────
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# ── File paths ─────────────────────────────────────────────────────────────────
bugs_model_file <- "model_files/abundance_model.bug"
bugs_directory  <- "" # e.g. "C:/Users/User/Documents/WinBUGS14"

# ── Database connection ────────────────────────────────────────────────────────
cfg <- config::get()
con <- DBI::dbConnect(
  RPostgres::Postgres(),
  dbname   = cfg$db_name,
  host     = cfg$db_host,
  port     = cfg$db_port,
  user     = cfg$db_user,
  password = cfg$db_password
)
on.exit(DBI::dbDisconnect(con), add = TRUE)

# ── Test site and year ─────────────────────────────────────────────────────────
test_site     <- "ubc"
test_run_year <- 2020

# =============================================================================
# UPLOAD
# =============================================================================

# ── 1. Fit and upload pCap all-sites ──────────────────────────────────────────
message("── Fitting pCap all-sites ────────────────────────────────────────────────")

pCap_inputs <- prepare_pCap_inputs(model_type = "all_sites")
pCap_fit    <- fit_pCap_model(pCap_inputs)

store_model_fit(
  con              = con,
  model_fit_object = pCap_fit,
  model_inputs     = pCap_inputs,
  description      = paste("TEST — pCap all-sites", Sys.Date())
)

# ── 2. Fit and upload one abundance model ─────────────────────────────────────
message(glue::glue("── Fitting abundance: {test_site} {test_run_year} ────────────────────────────"))

abundance_inputs <- prepare_abundance_inputs(
  site              = test_site,
  run_year          = test_run_year,
  pCap_model_type   = "all_sites",
  pCap_model_object = pCap_fit,
  min_pCap_mult     = 0.5
)

abundance_fit <- fit_abundance_model_BUGS(
  abundance_inputs,
  bugs_model_file,
  bugs_directory
)

store_model_fit(
  con              = con,
  model_fit_object = abundance_fit,
  model_inputs     = abundance_inputs,
  description      = paste("TEST —", test_site, test_run_year, "abundance", Sys.Date())
)

# =============================================================================
# DOWNLOAD
# =============================================================================

# ── 3. Retrieve the pCap fit and check it came back as a stanfit ──────────────
message("── Retrieving pCap all-sites fit ────────────────────────────────────────")

pCap_retrieved <- get_model_fit("pcap_all_sites")
stopifnot("stanfit" %in% class(pCap_retrieved))
message("  pCap download OK — class: ", paste(class(pCap_retrieved), collapse = ", "))

# ── 4. Retrieve the abundance fit and check it came back as a bugs object ─────
message(glue::glue("── Retrieving abundance fit: {test_site} {test_run_year} ──────────────────────"))

abundance_retrieved <- get_model_fit(abundance_inputs$model_name)
stopifnot("bugs" %in% class(abundance_retrieved))
message("  Abundance download OK — class: ", paste(class(abundance_retrieved), collapse = ", "))

# ── 5. Check the model_run table recorded both entries ────────────────────────
message("── Checking database records ─────────────────────────────────────────────")

recent_runs <- DBI::dbGetQuery(con, "
  SELECT mr.id, mn.name AS model_name, mr.site, mr.run_year,
         mr.description, mr.created_at
  FROM model_run mr
  JOIN model_name mn ON mn.id = mr.model_name_id
  ORDER BY mr.created_at DESC
  LIMIT 5
")
print(recent_runs)

message("── Test complete ─────────────────────────────────────────────────────────")
