

# bt-spas-x ----------------------------------------------------------

# stan (pCap and abundance separate)
# pCap_all_mark_recap <- readr::read_file(here::here("model_files", "pCap_all_mark_recap.stan"))
# pCap_missing_mark_recap <- readr::read_file(here::here("model_files", "pCap_missing_mark_recap.stan"))
# pCap_no_mark_recap <- readr::read_file(here::here("model_files", "pCap_no_mark_recap.stan"))
# pCap_no_mark_recap_no_trib <- readr::read_file(here::here("model_files", "pCap_no_mark_recap_no_trib.stan"))
pCap_all <- readr::read_file(here::here("model_files", "pCap_all.stan"))
abundance <- readr::read_file(here::here("model_files", "abundance_model.stan"))
abundance_BUGS <- readr::read_file(here::here("model_files", "abundance_model.bug"))

bt_spas_x_model_code <- list(#pCap_all_mark_recap = pCap_all_mark_recap,
                             #pCap_missing_mark_recap = pCap_missing_mark_recap,
                             #pCap_no_mark_recap = pCap_no_mark_recap,
                             #pCap_no_mark_recap_no_trib = pCap_no_mark_recap_no_trib,
                             pCap_all = pCap_all,
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

# stan (pCap and abundance combined)
all_mark_recap <- readr::read_file(here::here("model_files", "all_mark_recap.stan"))
missing_mark_recap <- readr::read_file(here::here("model_files", "missing_mark_recap.stan"))
no_mark_recap <- readr::read_file(here::here("model_files", "no_mark_recap.stan"))
no_mark_recap_no_trib <- readr::read_file(here::here("model_files", "no_mark_recap_no_trib.stan"))

stan <- list(all_mark_recap = all_mark_recap,
                              missing_mark_recap = missing_mark_recap,
                              no_mark_recap = no_mark_recap,
                              no_mark_recap_no_trib = no_mark_recap_no_trib)

# winBUGS no cut
all_mark_recap <- readr::read_file(here::here("model_files", "model_files_for_analysis",
                                              "no_cut_all_mark_recap.bug"))
missing_mark_recap <- readr::read_file(here::here("model_files", "model_files_for_analysis",
                                                  "no_cut_missing_mark_recap.bug"))
no_mark_recap <- readr::read_file(here::here("model_files", "model_files_for_analysis",
                                             "no_cut_no_mark_recap.bug"))
no_mark_recap_no_trib <- readr::read_file(here::here("model_files", "model_files_for_analysis",
                                                     "no_cut_no_mark_recap_no_trib.bug"))

no_cut <- list(all_mark_recap = all_mark_recap,
               missing_mark_recap = missing_mark_recap,
               no_mark_recap = no_mark_recap,
               no_mark_recap_no_trib = no_mark_recap_no_trib)

# winBUGS cut
all_mark_recap <- readr::read_file(here::here("model_files", "all_mark_recap.bug"))
missing_mark_recap <- readr::read_file(here::here("model_files", "missing_mark_recap.bug"))
no_mark_recap <- readr::read_file(here::here("model_files", "no_mark_recap.bug"))
no_mark_recap_no_trib <- readr::read_file(here::here("model_files", "no_mark_recap_no_trib.bug"))

cut <- list(all_mark_recap = all_mark_recap,
            missing_mark_recap = missing_mark_recap,
            no_mark_recap = no_mark_recap,
            no_mark_recap_no_trib = no_mark_recap_no_trib)

winbugs <- list(cut = cut,
                no_cut = no_cut)

# bt_spas_x_model_code <- list(stan = stan,
#                              winbugs = winbugs)
#
# usethis::use_data(bt_spas_x_model_code, overwrite = T)

# passage to spawner ------------------------------------------------------

p2s_model_code <- readr::read_file(here::here("model_files", "passage_to_spawner.stan"))
usethis::use_data(p2s_model_code, overwrite = T)
