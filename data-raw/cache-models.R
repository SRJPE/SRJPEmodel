library(dplyr)

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

usethis::use_data(bt_spas_x_model_code, overwrite = TRUE)

# read in Jwk lookup # TODO integrate into SRJPEdata package
julian_week_to_date_lookup <- read.table(here::here("data-raw", "juvenile_abundance", "Jwk_Dates.txt"), header = FALSE) |>
  tibble() |>
  filter(V1 != "Jwk") |>
  mutate(V1 = as.numeric(V1)) |>
  select(Jwk = V1, date = V2)

usethis::use_data(julian_week_to_date_lookup, overwrite = TRUE)
