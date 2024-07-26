#' Plot weekly estimated juvenile abundance
#' @details This will generate a bar plot for output of a single run of the bt spas x model.
#' @param bt_spas_x_results_clean: the `final_results` element from the object produced by running
#' `run_single_bt_spas_x()`
#' @returns a plot showing weekly estimated abundance from the model (with confidence intervals)
#' compared to lincoln-peterson estimates
#' @export
#' @md
plot_juvenile_abundance <- function(bt_spas_x_results_clean) {

  site_arg <- unique(bt_spas_x_results_clean$site)
  run_year_arg <- unique(bt_spas_x_results_clean$run_year)
  life_stage_arg <- unique(bt_spas_x_results_clean$life_stage)
  lcl <- 0.025
  ucl <- 0.975

  # don't need to use the "sims list" because summary table in BUGS calculates quantiles
  julian_week_to_date_lookup <- tibble("Jwk" = 1:53,
                                       "fake_year" = 1990) |>
    mutate(fake_date = ymd(paste0(fake_year, "-01-01")),
           final_date = fake_date + weeks(Jwk - 1),
           date = format(final_date, "%b-%d")) |>
    select(Jwk, date)

  # calculate lincoln peterson
  input_data <- SRJPEdata::weekly_juvenile_abundance_model_data |>
    filter(site == site_arg,
           run_year == run_year_arg,
           life_stage == life_stage_arg,
           week %in% c(seq(45, 53), seq(1, 22))) |>
    mutate(lincoln_peterson_abundance = count * (number_released / number_recaptured),
           lincoln_peterson_abundance = ifelse(is.infinite(lincoln_peterson_abundance), NA, lincoln_peterson_abundance)) # TODO sd calculation is what

  if(all(unique(input_data$lincoln_peterson_abundance) %in% c(NA, 0))) {
    cli::cli_bullets("The RST data for the site, run year, and life stage selected
                   do not have sufficient recapture data to produce a Lincoln Peterson
                   abundance estimate.")
  }

  summary_output <- bt_spas_x_results_clean |>
    select(-c(model_name, srjpedata_version)) |>
    pivot_wider(id_cols = site:parameter,
                values_from = value,
                names_from = statistic) |>
    mutate(cv = round(100 * (sd / mean), digits = 0)) # TODO confirm we use cv for more than just Ntot plot

  total_abundance <- summary_output |>
    filter(parameter == "Ntot")

  plot_data <- summary_output |>
    filter(str_detect(parameter, "N\\[")) |>
    left_join(julian_week_to_date_lookup, by = c("week_fit" = "Jwk")) |>
    mutate(fake_date = ifelse(week_fit > 35, paste0(run_year - 1, "-", date),
                              paste0(run_year, "-", date)),
           fake_date = as.Date(fake_date, format = "%Y-%b-%d")) |>
    left_join(input_data |>
                select(week, site, run_year, lincoln_peterson_abundance),
              by = c("site", "run_year", "week_fit" = "week"))

  # abundance bar plot
  plot_title <- paste0(str_to_title(site_arg), " ", run_year_arg, " predicted ", life_stage_arg, " annual abundance = ",
                       prettyNum(round(total_abundance$`50`, 0), big.mark = ","), " (",
                       prettyNum(round(total_abundance$`2.5`, 0), big.mark = ","), "-",
                       prettyNum(round(total_abundance$`97.5`, 0), big.mark = ","), ")")
  plot <- plot_data |>
    ggplot(aes(x = fake_date, y = `50`)) +
    geom_bar(stat = "identity", fill = "grey", width = 5) +
    geom_errorbar(aes(x = fake_date, ymin = `2.5`, ymax = `97.5`), width = 0.2) +
    geom_point(aes(x = fake_date, y = lincoln_peterson_abundance),
               shape = 1, color = "blue") +
    theme_minimal() +
    labs(x = "Date", y = "Abundance",
         title = plot_title) +
    scale_y_continuous(labels = label_comma())

  suppressWarnings(print(plot))

}
