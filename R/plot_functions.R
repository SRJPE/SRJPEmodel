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
    scale_y_continuous(labels = label_comma()) +
    scale_x_date(date_breaks = "1 week", date_labels = "%b %d") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  suppressWarnings(print(plot))

}

#' Plot weekly estimated capture probability
#' @details This will generate a bar plot for output of a single run of the bt spas x model.
#' @param bt_spas_x_results_clean: the `final_results` element from the object produced by running
#' `run_single_bt_spas_x()`
#' @returns a plot showing estimated capture probability
#' @export
#' @md
plot_weekly_capture_probability <- function(bt_spas_x_results_clean) {

  site_arg <- unique(bt_spas_x_results_clean$site)
  run_year_arg <- unique(bt_spas_x_results_clean$run_year)
  life_stage_arg <- unique(bt_spas_x_results_clean$life_stage)

  # to get trib-specific pCap, we need to get the indices for each trib
  # that was fit in the model.
  # TODO have bt spas x function return an index (like week fit, trib_fit) for these
  # parameters [1:Ntribs]
  if(!site_arg %in% c("knights landing", "tisdale", "red bluff diversion dam")) {
    remove_sites <- c("knights landing", "tisdale", "red bluff diversion dam")
  } else {
    remove_sites <- c("deer creek", "eye riffle", "live oak",
                      "okie dam", "mill creek", "yuba river", "herringer riffle", "ubc",
                      "lcc", "ucc", "hallwood", "steep riffle", "sunset pumps", "shawn's beach",
                      "gateway riffle", "lower feather river")
  }

  trib_index_data <- SRJPEdata::weekly_juvenile_abundance_model_data |>
    dplyr::filter(!site %in% remove_sites &
                  !is.na(standardized_flow),
                  !is.na(number_released) &
                  !is.na(number_recaptured)) |>
    select(-c(year, mean_fork_length, count, hours_fished, flow_cfs,
              catch_standardized_by_hours_fished, lgN_prior)) |>
    distinct(site)

  trib_index_arg <- which(unique(trib_index_data$site) == site_arg)

  if(length(trib_index_arg) == 0) {
    no_mark_recap_on_trib = T
    trib_index_arg <- -9999
    cli::cli_warn(paste0("No mark recapture was performed at ", site_arg))
  } else {
    no_mark_recap_on_trib = F
  }

  julian_week_to_date_lookup <- tibble("Jwk" = 1:53,
                                       "fake_year" = 1990) |>
    mutate(fake_date = ymd(paste0(fake_year, "-01-01")),
           final_date = fake_date + weeks(Jwk - 1),
           date = format(final_date, "%b-%d")) |>
    select(Jwk, date)

  transform_pCap <- function(estimate) {
    exp(estimate) / (1 + exp(estimate))
  }

  results_wide <- bt_spas_x_results_clean |>
    # transform proportion capture parameters to correct scale
    mutate(value = ifelse(str_detect(parameter, "trib_mu"), transform_pCap(value), value),
           value = ifelse(str_detect(parameter, "b0_pCap\\["), transform_pCap(value), value)) |>
    select(-c(model_name, srjpedata_version)) |>
    pivot_wider(id_cols = site:parameter,
                values_from = value,
                names_from = statistic)

  mean_capture_probability_across_tribs <- results_wide |>
    filter(parameter == "trib_mu.P")

  mean_capture_probability_trib_specific <- results_wide |>
    filter(str_detect(parameter, "b0_pCap\\[")) |>
    mutate(trib_index = substr(parameter, 3, length(parameter)),
           trib_index = readr::parse_number(trib_index)) |>
    filter(trib_index == trib_index_arg)

  plot_data <- results_wide |>
    filter(str_detect(parameter, "pCap_U\\[")) |>
    left_join(julian_week_to_date_lookup, by = c("week_fit" = "Jwk")) |>
    mutate(fake_date = ifelse(week_fit > 35, paste0(run_year - 1, "-", date),
                              paste0(run_year, "-", date)),
           fake_date = as.Date(fake_date, format = "%Y-%b-%d"))

  if(no_mark_recap_on_trib) {
    plot_title <- paste0("Average across-trib capture probability = ",
                         round(mean_capture_probability_across_tribs$`50`, 3), " (",
                         round(mean_capture_probability_across_tribs$`2.5`, 3), "-",
                         round(mean_capture_probability_across_tribs$`97.5`, 3), ")")
  } else {
    plot_title <- paste0(str_to_title(site_arg), " ", run_year_arg, " average ", life_stage_arg, " weekly capture probability = ",
                         round(mean_capture_probability_trib_specific$`50`, 3), " (",
                         round(mean_capture_probability_trib_specific$`2.5`, 3), "-",
                         round(mean_capture_probability_trib_specific$`97.5`, 3), ")",
                         "\n",
                         "Average across-trib capture probability = ",
                         round(mean_capture_probability_across_tribs$`50`, 3), " (",
                         round(mean_capture_probability_across_tribs$`2.5`, 3), "-",
                         round(mean_capture_probability_across_tribs$`97.5`, 3), ")")
  }

  plot <- plot_data |>
    ggplot(aes(x = fake_date, y = `50`)) +
    geom_errorbar(aes(x = fake_date, ymin = `2.5`, ymax = `97.5`), width = 0.2) +
    geom_point(color = "#A42820", size = 2) +
    theme_minimal() +
    labs(x = "Date", y = "Capture Probability",
         title = plot_title) +
    scale_y_continuous(labels = label_comma()) +
    scale_x_date(date_breaks = "1 week", date_labels = "%b %d") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  suppressWarnings(print(plot))

}
