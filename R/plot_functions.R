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
    select(site, run_year, life_stage, statistic, week_fit,
           parameter, value) |>
    #select(-c(model_name, srjpedata_version)) |>
    pivot_wider(id_cols = c(site, run_year, life_stage, week_fit, parameter),
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
                       "\n",
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


#' Plot model diagnostics
#' @param model_fit_summary_object The model summary object produced by a model fit
#' @returns A two-panel plot showing abundance estimates and weekly efficiency estimates compared
#' to lincoln-peterson estimates and raw count and efficiency trial data.
#' @export
#' @md
generate_diagnostic_plot <- function(model_fit_summary_object,
                                     model_language) {

  site_arg <- unique(model_fit_summary_object$site)
  run_year_arg <- unique(model_fit_summary_object$run_year)

  julian_week_to_date_lookup <- read.table(file = "data-raw/juvenile_abundance/archive/btspas_model_code/Jwk_Dates.txt", header = F) |>
    tibble() |>
    filter(V1 != "Jwk") |>
    mutate(V1 = as.numeric(V1)) |>
    select(Jwk = V1, date = V2)

  # input data all lifestages
  input_data <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
    filter(run_year == run_year_arg,
           site == site_arg) |> #,
           #week %in% c(seq(45, 53), seq(1, 22))) |>
    group_by(year, week, stream, site, run_year) |>
    summarise(count = sum(count, na.rm = T),
              mean_fork_length = mean(mean_fork_length, na.rm = T),
              hours_fished = mean(hours_fished, na.rm = T),
              flow_cfs = mean(flow_cfs, na.rm = T),
              average_stream_hours_fished = mean(average_stream_hours_fished, na.rm = T),
              standardized_flow = mean(standardized_flow, na.rm = T),
              catch_standardized_by_hours_fished = sum(catch_standardized_by_hours_fished, na.rm = T),
              lgN_prior = mean(lgN_prior, na.rm = T)) |>
    ungroup() |>
    left_join(SRJPEdata::weekly_juvenile_abundance_efficiency_data,
              by = c("year", "run_year", "week", "stream", "site")) |>
    mutate(count = round(count, 0),
           catch_standardized_by_hours_fished = round(catch_standardized_by_hours_fished, 0),
           # plot things
           lincoln_peterson_abundance = count * (number_released / number_recaptured),
           lincoln_peterson_efficiency = number_recaptured / number_released)

  # create week lookup for year to tell you if it was sampled or not
  sampled <- tibble("week" = c(45:53, 1:22),
                    "sampled" = TRUE,
                    "run_year" = run_year_arg) |>
    mutate(year = ifelse(week >= 43, run_year - 1, run_year),
           fake_date = ymd(paste0(year, "-01-01")),
           final_date = fake_date + weeks(week - 1),
           date = format(final_date, "%b-%d")) |>
    select(date, run_year, week, sampled)

  input_data <- full_join(input_data, sampled,
                          by = c("week", "run_year")) |>
    mutate(sampled = ifelse(is.na(count), FALSE, sampled),
           site = ifelse(is.na(site), site_arg, site),
           run_year = ifelse(is.na(run_year), run_year_arg, run_year)) |>
    filter(week %in% c(45:53, 1:22))

  if(model_language == "STAN") {
    summary_output <- model_fit_summary_object |>
      select(-c(model_name, srjpedata_version, converged, life_stage)) |>
      mutate(statistic = str_remove_all(statistic, "\\%")) |>
      filter(!parameter %in% c("logit_pCap", "pro_dev_P")) |>
      pivot_wider(id_cols = c(site, run_year, week_fit, site_fit_hierarchical, parameter),
                  values_from = value,
                  names_from = statistic) |>
      mutate(cv = round(100 * (sd / mean), digits = 0))

    total_abundance <- summary_output |>
      filter(parameter == "Ntot")

  } else if(model_language == "BUGS") {
    summary_output <- model_fit_summary_object |>
      select(-c(model_name, srjpedata_version, converged, life_stage)) |>
      mutate(parameter = gsub("[0-9]+|\\[|\\]", "", parameter)) |>
      filter(!parameter %in% c("logit_pCap", "pro_dev_P")) |>
      pivot_wider(id_cols = c(site, run_year, week_fit, site_fit_hierarchical, parameter),
                  values_from = value,
                  names_from = statistic) |>
      mutate(cv = round(100 * (sd / mean), digits = 0))

    total_abundance <- model_fit_summary_object |>
      mutate(parameter = gsub("[0-9]+|\\[|\\]", "", parameter)) |>
      filter(parameter == "N",
             statistic == "50") |>
      group_by(site, run_year) |> # not lifestage
      summarise(Ntot = sum(value)) |>
      ungroup()
  }

  # plot abundance only
  # set up cv, lcl, and ucl and make sure to scale to 1000s

  abundance_plot_data <- summary_output |>
    filter(parameter == "N") |>
    left_join(julian_week_to_date_lookup, by = c("week_fit" = "Jwk")) |>
    full_join(input_data |>
              select(week, site, run_year, count, lincoln_peterson_abundance, sampled),
              by = c("site", "run_year", "week_fit" = "week")) |>
    # fill in date for weeks not sampled
    mutate(year = ifelse(week_fit >= 43, run_year - 1, run_year),
           fake_date = ymd(paste0(year, "-01-01")),
           final_date = fake_date + weeks(week_fit - 1),
           date = format(final_date, "%b-%d")) |>
    select(-c(year, fake_date, final_date)) |>
    mutate(fake_date = ifelse(week_fit > 35, paste0(run_year - 1, "-", date),
                              paste0(run_year, "-", date)),
           fake_date = as.Date(fake_date, format = "%Y-%b-%d"))

  # abundance bar plot
  if(model_language == "STAN") {
    abundance_plot_title <- paste0(str_to_title(site_arg), " ", run_year_arg, " predicted annual abundance = ",
                                   prettyNum(round(total_abundance$`50`, 0), big.mark = ","), " (",
                                   prettyNum(round(total_abundance$`2.5`, 0), big.mark = ","), "-",
                                   prettyNum(round(total_abundance$`97.5`, 0), big.mark = ","), ")")
  } else {
    # TODO get uncertainty around total estimate
    abundance_plot_title <- paste0(str_to_title(site_arg), " ", run_year_arg, " predicted annual abundance = ",
                                   prettyNum(round(total_abundance$Ntot, 0), big.mark = ","))
  }

  abundance_plot <- abundance_plot_data |>
    ggplot(aes(x = fake_date, y = `50`)) +
    geom_bar(stat = "identity", fill = "grey") +
    geom_errorbar(aes(x = fake_date, ymin = `2.5`, ymax = `97.5`), width = 0.2) +
    geom_point(aes(x = fake_date, y = lincoln_peterson_abundance),
               shape = 1, color = "blue") +
    geom_point(aes(x = fake_date, y = Inf, color = sampled),
               size = 3) +
    geom_text(aes(x = fake_date, y = Inf,
                  label = paste(count),
                  angle = 90),
              hjust = 1,
              size = 3) +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "red")) +
    theme_minimal() +
    labs(x = "Date", y = "Abundance",
         title = abundance_plot_title) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "")

  efficiency_plot_data <- summary_output |>
    filter(parameter == "pCap_U") |>
    left_join(julian_week_to_date_lookup, by = c("week_fit" = "Jwk")) |>
    full_join(input_data |>
                select(week, site, run_year, number_recaptured, number_released,
                       lincoln_peterson_efficiency, sampled),
              by = c("site", "run_year", "week_fit" = "week")) |>
    # fill in date for weeks not sampled
    mutate(year = ifelse(week_fit >= 43, run_year - 1, run_year),
           fake_date = ymd(paste0(year, "-01-01")),
           final_date = fake_date + weeks(week_fit - 1),
           date = format(final_date, "%b-%d")) |>
    select(-c(year, fake_date, final_date)) |>
    mutate(fake_date = ifelse(week_fit > 35, paste0(run_year - 1, "-", date),
                              paste0(run_year, "-", date)),
           fake_date = as.Date(fake_date, format = "%Y-%b-%d"),
           number_released_label = ifelse(is.na(number_released), "", number_released),
           number_recaptured_label = ifelse(is.na(number_recaptured), "", number_recaptured))

  efficiency_plot <- efficiency_plot_data |>
    ggplot(aes(x = fake_date, y = `50`)) +
    geom_bar(stat = "identity", fill = "grey") +
    geom_errorbar(aes(x = fake_date, ymin = `2.5`, ymax = `97.5`), width = 0.2) +
    geom_point(aes(x = fake_date, y = lincoln_peterson_efficiency),
               shape = 1, color = "blue") +
    geom_text(aes(x = fake_date, y = Inf,
                  label = paste(number_released_label, number_recaptured_label),
                  angle = 90),
              hjust = 1,
              size = 3) +
    # geom_point(aes(x = fake_date, y = Inf, color = sampled),
    #            size = 3) +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "red")) +
    theme_minimal() +
    labs(x = "Date", y = "Weekly Efficiency") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "")

  # arrange together
  gridExtra::grid.arrange(abundance_plot, efficiency_plot)

}
