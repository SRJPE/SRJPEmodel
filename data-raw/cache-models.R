

# bt-spas-x ----------------------------------------------------------

# stan (pCap and abundance separate)
pCap_all_sites <- readr::read_file(here::here("model_files", "pCap_all_sites.stan"))
pCap_one_site <- readr::read_file(here::here("model_files", "pCap_one_site.stan"))
pCap_one_site_skew <- readr::read_file(here::here("model_files", "pCap_one_site_skew_re.stan"))
abundance_BUGS <- readr::read_file(here::here("model_files", "abundance_model.bug"))

bt_spas_x_model_code <- list(pCap_all_sites = pCap_all_sites,
                             pCap_one_site = pCap_one_site,
                             pCap_one_site_skew = pCap_one_site_skew,
                             abundance_BUGS = abundance_BUGS)

usethis::use_data(bt_spas_x_model_code, overwrite = T)

# read in Jwk lookup # TODO integrate into SRJPEdata package
julian_week_to_date_lookup <- read.table(file = "data-raw/juvenile_abundance/archive/btspas_model_code/Jwk_Dates.txt", header = F) |>
  tibble() |>
  filter(V1 != "Jwk") |>
  mutate(V1 = as.numeric(V1)) |>
  select(Jwk = V1, date = V2)

usethis::use_data(julian_week_to_date_lookup, overwrite = T)

# passage to spawner ------------------------------------------------------

p2s_model_code <- readr::read_file(here::here("model_files", "passage_to_spawner.txt"))
usethis::use_data(p2s_model_code, overwrite = T)

# stock recruit -----------------------------------------------------------

stock_recruit_model_code <- readr::read_file(here::here("model_files", "ricker_stock_recruit.stan"))
usethis::use_data(stock_recruit_model_code, overwrite = T)

# survival ----------------------------------------------------------------
survival_CovWY <- readr::read_file(here::here("model_files", "survival_CovWY.stan"))
survival_NoCov <- readr::read_file(here::here("model_files", "survival_NoCov.stan"))

survival_model_code <- list(survival_CovWY = survival_CovWY,
                            survival_NoCov = survival_NoCov)
usethis::use_data(survival_model_code, overwrite = T)


# in-season ---------------------------------------------------------------

in_season_autocorrelation <- readr::read_file(here::here("model_files", "BetaDevHBMRT_lag1.stan"))
in_season_no_autocorrelation <- readr::read_file(here::here("model_files", "BetaDevHBMRT.stan"))

in_season_model_code <- list(autocorrelation = in_season_autocorrelation,
                             no_autocorrelation = in_season_no_autocorrelation)

usethis::use_data(in_season_model_code, overwrite = T)
