library(gridExtra)

# plot code for diagnostic plots
options(mc.cores=parallel::detectCores())
ubc_2002 <- run_single_bt_spas_x_stan(SRJPEmodel::bt_spas_x_bayes_params,
                                      SRJPEdata::weekly_juvenile_abundance_catch_data,
                                      SRJPEdata::weekly_juvenile_abundance_efficiency_data,
                                      site = "ubc",
                                      run_year = 2002,
                                      lifestage = NA, # group
                                      effort_adjust = T)

# inputs
model_fit_summary_object <- ubc_2002$summary_table
site <- "ubc"
run_year <- 2002

julian_week_to_date_lookup <- read.table(file = "data-raw/juvenile_abundance/archive/btspas_model_code/Jwk_Dates.txt", header = F) |>
  tibble() |>
  filter(V1 != "Jwk") |>
  mutate(V1 = as.numeric(V1)) |>
  select(Jwk = V1, date = V2)

# input data all lifestages
input_data <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
  filter(run_year == 2002,
         site == "ubc",
         week %in% c(seq(45, 53), seq(1, 22))) |>
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
         lincoln_peterson_abundance = count * (number_released / number_recaptured))

# create week lookup for year to tell you if it was sampled or not
# TODO how to catch weeks with no sampling?
sampled <- c(45:53, 1:22) %in% input_data$week

input_data$sampled <- sampled

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

# plot abundance only
# set up cv, lcl, and ucl and make sure to scale to 1000s

abundance_plot_data <- summary_output |>
  filter(parameter == "N") |>
  left_join(julian_week_to_date_lookup, by = c("week_fit" = "Jwk")) |>
  mutate(fake_date = ifelse(week_fit > 35, paste0(run_year - 1, "-", date),
                            paste0(run_year, "-", date)),
         fake_date = as.Date(fake_date, format = "%Y-%b-%d")) |>
  left_join(input_data |>
              select(week, site, run_year, count, lincoln_peterson_abundance),
            by = c("site", "run_year", "week_fit" = "week"))

# abundance bar plot
# TODO add sample (T/F)
# TODO add catch numbers
abundance_plot_title <- paste0(str_to_title(site), " ", run_year, " predicted annual abundance = ",
                     prettyNum(round(total_abundance$`50`, 0), big.mark = ","), " (",
                     prettyNum(round(total_abundance$`2.5`, 0), big.mark = ","), "-",
                     prettyNum(round(total_abundance$`97.5`, 0), big.mark = ","), ")")
abundance_plot <- abundance_plot_data |>
  ggplot(aes(x = fake_date, y = `50`)) +
  geom_bar(stat = "identity", fill = "grey") +
  geom_errorbar(aes(x = fake_date, ymin = `2.5`, ymax = `97.5`), width = 0.2) +
  geom_point(aes(x = fake_date, y = lincoln_peterson_abundance),
             shape = 1, color = "blue") +
  geom_text(aes(x = fake_date, y = Inf,
                label = paste(count),
                angle = 90),
            hjust = 1,
            size = 3) +
  geom_point(aes(x = fake_date, y = -Inf, color = sampled)) +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "red")) +
  theme_minimal() +
  labs(x = "Date", y = "Abundance",
       title = abundance_plot_title) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "")

efficiency_plot_data <- summary_output |>
  filter(parameter == "pCap_U") |>
  left_join(julian_week_to_date_lookup, by = c("week_fit" = "Jwk")) |>
  mutate(fake_date = ifelse(week_fit > 35, paste0(run_year - 1, "-", date),
                            paste0(run_year, "-", date)),
         fake_date = as.Date(fake_date, format = "%Y-%b-%d")) |>
  left_join(input_data |>
              select(week, site, run_year, number_recaptured, number_released,
                     lincoln_peterson_abundance),
            by = c("site", "run_year", "week_fit" = "week"))

efficiency_plot <- efficiency_plot_data |>
  ggplot(aes(x = fake_date, y = `50`)) +
  geom_bar(stat = "identity", fill = "grey") +
  geom_errorbar(aes(x = fake_date, ymin = `2.5`, ymax = `97.5`), width = 0.2) +
  geom_point(aes(x = fake_date, y = lincoln_peterson_abundance),
             shape = 1, color = "blue") +
  geom_text(aes(x = fake_date, y = Inf,
                label = paste(number_released, number_recaptured),
                angle = 90),
            hjust = 1,
            size = 3) +
  geom_point(aes(x = fake_date, y = -Inf, color = sampled)) +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "red")) +
  theme_minimal() +
  labs(x = "Date", y = "Weekly Efficiency") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "")

gridExtra::grid.arrange(abundance_plot, efficiency_plot)
