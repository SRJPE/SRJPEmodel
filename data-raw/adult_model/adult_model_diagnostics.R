# diagnostics

# plot predicted spawner estimates against observed spawner values

diagnostics <- bind_rows(get_report_pars(battle_pred_spawners, "battle creek"),
                         get_report_pars(clear_pred_spawners, "clear creek"),
                         get_report_pars(deer_pred_spawners, "deer creek"),
                         get_report_pars(mill_pred_spawners, "mill creek")) |>
  filter(par_names != "b1_survival") |>
  select(par_names, mean, sd, stream) |>
  glimpse()

obsv_data <- full_data_for_input |>
  mutate(obsv_spawner_count = ifelse(stream == "deer creek", holding_count, redd_count)) |>
  select(year, stream, obsv_spawner_count) |>
  glimpse()

years_to_join <- c(full_data_for_input |> filter(stream == "battle creek") |> drop_na(wy_type_std) |> pull(year),
                   full_data_for_input |> filter(stream == "clear creek") |> drop_na(max_flow_std) |> pull(year),
                   full_data_for_input |> filter(stream == "deer creek") |> drop_na(wy_type_std) |> pull(year),
                   full_data_for_input |> filter(stream == "mill creek") |> drop_na(gdd_std) |> pull(year)) |>
  glimpse()

diagnostics_with_year <- diagnostics |>
  arrange(stream) |>
  mutate(year = years_to_join) |>
  select(-par_names) |>
  rename(pred_spawner_count = mean) |>
  left_join(obsv_data, by = c("year", "stream")) |>
  glimpse()

diagnostics_with_year |> ggplot(aes(x = obsv_spawner_count, y = pred_spawner_count)) +
  geom_point() + geom_smooth(method = "lm") + facet_wrap(~stream, scales = "free") +
  theme_minimal() + xlab("Observed Spawner Count") + ylab("Predicted Spawner Count")

# look at R2
summary(lm(pred_spawner_count ~ obsv_spawner_count, data = diagnostics_with_year))

# look at estimated value of sigma redds_per_spawner (mean and sd)
get_other_pars <- function(model_fit, stream_name) {
  par_results <- summary(model_fit)$summary
  results_tibble <- as.data.frame(par_results) |>
    rownames_to_column("par_names") |>
    filter(!str_detect(par_names, "predicted_spawners"),
           par_names != "b1_survival") |>
    mutate(stream = stream_name)

  return(results_tibble)
}

other_pars <- bind_rows(get_other_pars(battle_pred_spawners, "battle creek"),
                        get_other_pars(clear_pred_spawners, "clear creek"),
                        get_other_pars(deer_pred_spawners, "deer creek"),
                        get_other_pars(mill_pred_spawners, "mill creek"))

other_pars |>
  filter(str_detect(par_names, "log_redds_per_spawner") |
         str_detect(par_names, "survival_rate")) |>
  mutate(par_name = ifelse(str_detect(par_names, "log_redds_per_spawner"), "log_rps", "survival")) |>
  arrange(par_name) |>
  mutate(year = rep(years_to_join, 2)) |>
  ggplot(aes(x = year, y = mean, color = par_name)) +
  geom_point() +
  scale_color_discrete(name = "Parameter name",
                       labels = c("Log redds/spawner (RE)",
                                  "Predicted survival rate")) +
  facet_wrap(~stream, scales = "free") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  xlab("Year") + ylab("Mean estimated value") + theme_minimal()  +
  theme(legend.position = "bottom")

# plot year-specific random effects (log_redds_per_spawner[y])
# a big portion of the variation in predicted survival_rate[y] comes from b1_survival * environmental_covar[y]
# and not random effects
