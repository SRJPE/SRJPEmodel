

# bt-spas-x ----------------------------------------------------------

# stan (pCap and abundance separate)
pCap_all <- readr::read_file(here::here("model_files", "pCap_all.stan"))
abundance <- readr::read_file(here::here("model_files", "abundance_model.stan"))
abundance_BUGS <- readr::read_file(here::here("model_files", "abundance_model.bug"))

bt_spas_x_model_code <- list(pCap_all = pCap_all,
                             abundance = abundance,
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


# survival ----------------------------------------------------------------
survival_CovIndWY <- readr::read_file(here::here("model_files", "survival_CovIndWY.stan"))
survival_CovWY <- readr::read_file(here::here("model_files", "survival_CovWY.stan"))
survival_NoCov <- readr::read_file(here::here("model_files", "survival_NoCov.stan"))

survival_model_code <- list(survival_CovIndWY = survival_CovIndWY,
                            survival_CovWY = survival_CovWY,
                            survival_NoCov = survival_NoCov)
usethis::use_data(survival_model_code, overwrite = T)
