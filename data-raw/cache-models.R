

# bt-spas-x ----------------------------------------------------------

# stan
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

bt_spas_x_model_code <- list(stan = stan,
                             winbugs = winbugs)

usethis::use_data(bt_spas_x_model_code, overwrite = T)

# passage to spawner ------------------------------------------------------

p2s_model_code <- readr::read_file(here::here("model_files", "passage_to_spawner.stan"))
usethis::use_data(p2s_model_code, overwrite = T)
