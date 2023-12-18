# refactor of GetData.R

# TODO determine where this came from - from Flora's file?
detections_with_covariates_sac <- read_csv(file=here::here("data-raw", "survival_model_STAN",
                                                           "SacInp_withcov.csv"))

detections <- detections_with_covariates_sac |>
  # reorder release groups from upstream --> downstream
  mutate(release_group_order = case_when(StudyID == "ColemanFall_2013" ~ 1,
                                         StudyID == "MillCk_Wild_CHK_2013" ~ 2,
                                         StudyID == "MillCk_Wild_CHK_2015" ~ 3,
                                         StudyID == "ColemanFall_2016" ~ 4,
                                         StudyID == "ColemanFall_2017" ~ 5,
                                         StudyID == "RBDD_2017" ~ 6,
                                         StudyID == "MillCk_Wild_CHK_2017" ~ 7,
                                         StudyID == "RBDD_2018" ~ 8,
                                         StudyID == "DeerCk_Wild_CHK_2018" ~ 9,
                                         StudyID == "CNFH_FMR_2019" ~ 10,
                                         StudyID == "CNFH_FMR_2020" ~ 11,
                                         StudyID == "DeerCk_Wild_CHK_2020" ~ 12,
                                         StudyID == "CNFH_FMR_2021" ~ 13)) |>
  glimpse()

n_reaches <- 4
reach_names <- c("release-woodson", "woodson-butte", "butte-sac", "sac-delta")
n_detection_locations <- 5
detection_location_names <- c("release", "woodson", "butte", "sac", "delta")
reach_river_km <- c(40, 88, 170, 110) # TODO reference josh file Summary.xlsx
n_individuals <- nrow(detections)
years <- unique(detections$year) |>
  sort()
n_years <- length(years)
release_groups <- detections |>
  distinct(StudyID) |>
  pull()
n_release_groups <- length(release_groups)

reach_covariate_index <- c(1, 2, 3, 3) # TODO when is this used? note Butte-Sac and Sac-Delta have same fixed effect index


capture_history_matrix <- detections |>
  select(FishID, ch, StudyID, year) |>
  rowwise() |>
  mutate(first_capture = min(unlist(str_locate_all(ch, "1"))),
         last_capture = max(unlist(str_locate_all(ch, "1")))) |>
  ungroup() |>
  mutate(water_year_2 = ifelse(year %in% c(2013, 2015, 2020, 2021), 0, 1),
         water_year_3 = case_when(year %in% c(2015, 2021) ~ 0,
                                  year %in% c(2013, 2016, 2018, 2020) ~ 1,
                                  TRUE ~ 2)) |>
  separate_wider_position(ch, c("1" = 1, "2" = 1, "3" = 1, "4" = 1, "5" = 1)) |>
  rename_at(vars(`1`:`5`), ~ detection_location_names) |> # rename to detection location
  glimpse()

# this is sumout
data_summary <- capture_history_matrix |>
  pivot_longer(release:delta, names_to = "detection_location", values_to = "detections") |>
  mutate(detections = as.numeric(detections)) |>
  group_by(StudyID, detection_location) |>
  summarise(n = sum(detections)) |>
  pivot_wider(id_cols = "StudyID", names_from = detection_location, values_from = n) |>
  left_join(detections |> distinct(StudyID, release_group_order)) |>
  arrange(release_group_order) |>
  select(-c(release_group_order))

# IndStats
individual_summaries <- detections |>
  arrange(release_group_order) |>
  group_by(StudyID) |>
  summarise(mean_fork_length = mean(fish_length, na.rm = T),
            var_fork_length = sd(fish_length, na.rm = T)/mean_fork_length,
            mean_weight = mean(fish_weight, na.rm = T),
            var_weight = sd(fish_weight, na.rm = T)/mean_weight,
            mean_condition = mean(fish_k, na.rm = T),
            var_condition = sd(fish_k, na.rm = T)/mean_condition,
            n = n())


# modeling prep -----------------------------------------------------------
n_covariate_inc <- 25 # number of increments to calculate and plot covariate relationship over
standard_reach_lengths <- 1000 # survival calculated for a standardized reach length of 100 km, then converted to reach specific survival in stan models
r_mult <- reach_river_km / standard_reach_lengths
