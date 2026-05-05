# btspasx_fit_upload_script.R
# Fits and uploads all pCap and abundance (BT-SPAS-X) model objects to Azure
# Blob Storage and records each run in the JPE database.
#
# pCap models:
#   - all_sites       : one fit covering all tributary sites
#   - one_site_skew   : individual fits for knights landing and tisdale
#
# Abundance models:
#   - all tributary sites (uses pcap_all_sites fit)
#   - knights landing and tisdale (each uses their own one_site_skew pCap fit)
#
# Prerequisites:
#   - Must be run on a PC (WinBUGS requirement)
#   - AZ_CONTAINER_ACCESS_KEY set in your .Renviron
#   - config.yml present with database credentials

library(tidyverse)
library(SRJPEdata)
library(SRJPEmodel)
library(rstan)
library(R2WinBUGS)

# ── STAN options ───────────────────────────────────────────────────────────────
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# ── Check OS (WinBUGS requires Windows) ───────────────────────────────────────
operating_system <- ifelse(
  grepl("Mac", Sys.info()["nodename"]) | grepl("MBP", Sys.info()["nodename"]),
  "mac", "pc"
)
if (operating_system != "pc") {
  stop("WinBUGS requires a PC. Use a Windows machine to run this script.")
}

# ── File paths ─────────────────────────────────────────────────────────────────
bugs_model_file <- "model_files/abundance_model.bug"
bugs_directory  <- "" # point to your local WinBUGS14 directory,
# e.g. "C:/Users/User/Documents/WinBUGS14"

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

# ── Site order (tribs first, mainstem last) ────────────────────────────────────
trib_sites     <- c("ucc", "lcc", "ubc", "mill creek", "deer creek", "okie dam",
                    "steep riffle", "gateway riffle", "eye riffle", "herringer riffle",
                    "live oak", "sunset pumps", "hallwood")
mainstem_sites <- c("tisdale", "knights landing")


# ── 1. pCap — all tributary sites ─────────────────────────────────────────────
message("── Fitting pCap all-sites model ──────────────────────────────────────────")

pCap_inputs_all <- prepare_pCap_inputs(model_type = "all_sites")
pCap_all        <- fit_pCap_model(pCap_inputs_all)

store_model_fit(
  con              = con,
  model_fit_object = pCap_all,
  model_inputs     = pCap_inputs_all,
  description      = paste("pCap all-sites fit", Sys.Date())
)


# ── 2. pCap — mainstem sites (one_site_skew per site) ─────────────────────────
pCap_mainstem <- list() # store fits in memory for use in abundance loop below

for (site in mainstem_sites) {
  message(glue::glue("── Fitting pCap one-site-skew for {site} ──────────────────────────────────"))

  pCap_inputs_site    <- prepare_pCap_inputs(model_type    = "one_site",
                                             skew          = TRUE,
                                             site_selection = site)
  pCap_mainstem[[site]] <- fit_pCap_model(pCap_inputs_site)

  store_model_fit(
    con              = con,
    model_fit_object = pCap_mainstem[[site]],
    model_inputs     = pCap_inputs_site,
    description      = paste("pCap one-site-skew fit —", site, Sys.Date())
  )
}


# ── 3. Abundance — all tributary sites ────────────────────────────────────────
message("── Fitting abundance models for tributary sites ───────────────────────────")

JPE_sites_trib <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
  dplyr::distinct(site, run_year) |>
  dplyr::filter(site %in% trib_sites) |>
  # exclude hallwood 2024 (known data issue)
  dplyr::filter(!(site == "hallwood" & run_year == 2024)) |>
  dplyr::arrange(match(site, trib_sites), run_year)

for (i in seq_len(nrow(JPE_sites_trib))) {

  .site     <- JPE_sites_trib$site[i]
  .run_year <- JPE_sites_trib$run_year[i]

  message(glue::glue("  Fitting abundance: {.site} {.run_year}"))

  abundance_inputs <- prepare_abundance_inputs(
    site             = .site,
    run_year         = .run_year,
    pCap_model_type  = "all_sites",
    pCap_model_object = pCap_all,
    min_pCap_mult    = 0.5
  )

  if (is.null(abundance_inputs)) {
    message(glue::glue("  Skipping {.site} {.run_year} — no usable input data."))
    next
  }

  abundance <- fit_abundance_model_BUGS(
    abundance_inputs,
    bugs_model_file,
    bugs_directory
  )

  store_model_fit(
    con              = con,
    model_fit_object = abundance,
    model_inputs     = abundance_inputs,
    description      = paste(.site, .run_year, "abundance fit", Sys.Date())
  )
}


# ── 4. Abundance — mainstem sites ─────────────────────────────────────────────
message("── Fitting abundance models for mainstem sites ────────────────────────────")

for (site in mainstem_sites) {

  site_years <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
    dplyr::filter(site == !!site) |>
    dplyr::distinct(run_year) |>
    dplyr::arrange(run_year) |>
    dplyr::pull(run_year)

  for (run_year in site_years) {

    message(glue::glue("  Fitting abundance: {site} {run_year}"))

    abundance_inputs <- prepare_abundance_inputs(
      site              = site,
      run_year          = run_year,
      pCap_model_type   = "one_site_skew",
      pCap_model_object = pCap_mainstem[[site]],
      min_pCap_mult     = 0.5
    )

    if (is.null(abundance_inputs)) {
      message(glue::glue("  Skipping {site} {run_year} — no usable input data."))
      next
    }

    abundance <- fit_abundance_model_BUGS(
      abundance_inputs,
      bugs_model_file,
      bugs_directory
    )

    store_model_fit(
      con              = con,
      model_fit_object = abundance,
      model_inputs     = abundance_inputs,
      description      = paste(site, run_year, "abundance fit", Sys.Date())
    )
  }
}

message("── All models fitted and uploaded successfully ────────────────────────────")
