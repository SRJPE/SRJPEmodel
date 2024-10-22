

# bt-spas-x stan ----------------------------------------------------------

all_mark_recap <- readr::read_file(here::here("model_files", "all_mark_recap.stan"))
missing_mark_recap <- readr::read_file(here::here("model_files", "missing_mark_recap.stan"))
no_mark_recap <- readr::read_file(here::here("model_files", "no_mark_recap.stan"))
no_mark_recap_no_trib <- readr::read_file(here::here("model_files", "no_mark_recap_no_trib.stan"))

bt_spas_x_model_code <- list(all_mark_recap = all_mark_recap,
                             missing_mark_recap = missing_mark_recap,
                             no_mark_recap = no_mark_recap,
                             no_mark_recap_no_trib = no_mark_recap_no_trib)
usethis::use_data(bt_spas_x_model_code, overwrite = T)
