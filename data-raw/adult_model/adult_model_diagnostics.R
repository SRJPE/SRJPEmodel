# diagnostics
# libraries and read in predicted model data ------------------------------
library(ggplot2)
library(googleCloudStorageR)
library(tidyverse)
library(wesanderson)
library(latex2exp)
library(SRJPEdata)


# join observed adult input and covariates --------------------------------
full_data_for_input <- SRJPEdata::observed_adult_input |>
  select(-reach) |> # empty
  group_by(year, stream, data_type, ) |>
  summarise(count = sum(count, na.rm = T)) |> # count adipose clipped, run together
  ungroup() |>
  full_join(SRJPEdata::p2s_model_covariates_standard,
            by = c("year", "stream")) |>
  filter(!is.na(data_type)) |>
  pivot_wider(id_cols = c(year, stream, wy_type, max_flow_std, gdd_std,
                          median_passage_timing_std, passage_index),
              names_from = data_type,
              values_from = count) |>
  glimpse()


# plot predicted spawner estimates against observed spawner values
pred_spawners <- get_predicted_spawners_from_P2S(SRJPEmodel::P2S_model_fits) |>
  glimpse()

obsv_data <- full_data_for_input |>
  mutate(obsv_spawner_count = ifelse(stream == "deer creek", holding_count, redd_count)) |>
  select(year, stream, obsv_spawner_count, obsv_upstream = upstream_estimate) |>
  glimpse()

pred_with_year <- pred_spawners |>
  rename(pred_spawner_count = median_predicted_spawners) |>
  left_join(obsv_data, by = c("year", "stream")) |>
  glimpse()

obsv_vs_pred_plot <- pred_with_year |>
  mutate(stream = str_to_title(stream)) |>
  ggplot(aes(x = obsv_spawner_count, y = pred_spawner_count)) +
  geom_smooth(color = "#FD6467", method = "lm") +
  geom_point(alpha = 0.8) +
  facet_wrap(~stream, scales = "free") +
  theme_minimal() +
  xlab("Observed Spawner Count") + ylab("Predicted Spawner Count") +
  ggtitle("Predicted vs. Observed Spawner Count")

ggsave(here::here("data-raw", "adult_model", "adult_model_plots",
                  "obsv_vs_pred_plot.jpg"),
       obsv_vs_pred_plot,
       width = 8, height = 6)

# look at R2
summary(lm(pred_spawner_count ~ obsv_spawner_count, data = pred_with_year |>
             filter(stream == "battle creek")))$r.squared

# look at estimated value of sigma redds_per_spawner (mean and sd)
fixed_random_effects <- SRJPEmodel::P2S_model_fits |>
  filter(str_detect(par_names, "log_redds_per_spawner") |
         str_detect(par_names, "conversion_rate")) |>
  mutate(par_name = ifelse(str_detect(par_names, "log_redds_per_spawner"), "log_rps", "survival")) |>
  rename(median = `50%`) |>
  mutate(stream = str_to_title(stream)) |>
  ggplot(aes(x = year, y = median, color = par_name)) +
  geom_point() +
  geom_line(alpha = 0.5) +
  facet_wrap(~stream, scales = "free") +
  scale_color_manual(values = wes_palette("GrandBudapest1")[3:2],
                     name = "Parameter",
                     labels = c(TeX("$log(\\mu_{\\delta})$"),
                                TeX("$R_{y}$"))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlab("Year") + ylab("Median estimated value") + theme_minimal()  +
  theme(legend.position = "bottom") +
  ggtitle("Year-specific random effects vs. predicted conversion rate")

ggsave(here::here("data-raw", "adult_model", "adult_model_plots",
                  "random_effects_vs_conversion_rate.jpg"),
       fixed_random_effects,
       width = 8, height = 6)


# diagnostic plots --------------------------------------------------------
rps <- SRJPEmodel::P2S_model_fits |>
  filter(str_detect(par_names, "log_mean_redds_per_spawner")) |>
  mutate(rps = exp(`50%`),
         lcl = exp(`2.5%`),
         ucl = exp(`97.5%`),
         log_rps = `50%`)

annual_random_effects <- SRJPEmodel::P2S_model_fits |>
  filter(str_detect(par_names, "log_redds_per_spawner")) |>
  mutate(annual_random_effect = `50%`) |>
  select(stream, year_index, annual_random_effect) |>
  select(-year_index)

# write b1_survival stats -------------------------------------------------
R2_estimates <- SRJPEmodel::P2S_model_fits |>
  filter(par_names %in% c("R2_data", "R2_fixed")) |>
  pivot_wider(id_cols = stream,
              names_from = par_names,
              values_from = mean)

b1_table <- SRJPEmodel::P2S_model_fits |>
  filter(par_names %in% c("b1_survival")) |>
  select(stream, parameter = par_names, `2.5%`,
         `50%`, `97.5%`) |>
  left_join(R2_estimates,
            by = "stream")


# plot contribution of random effects -------------------------------------
b1_effect <- obsv_data |>
  drop_na(obsv_upstream, obsv_spawner_count) |>
  left_join(b1_table |>
              select(stream, b1 = `50%`), by = "stream") |>
  left_join(full_data_for_input |>
              mutate(covar = wy_type) |>
              select(year, stream, covar)) |>
  drop_na(covar) |>
  mutate(b1_effect = (b1 * covar))
  #mutate(b1_effect = obsv_upstream * (b1 * covar))

create_josh_plot <- function(obsv_data, b1_effect, rps, annual_random_effects,
                             stream_choice_arg = c("ALL", "battle creek", "clear creek",
                                               "deer creek", "mill creek")) {

  if(stream_choice_arg == "ALL") {
    stream_choice <- c("battle creek", "clear creek",
                       "deer creek", "mill creek")
  } else {
    stream_choice <- stream_choice_arg
  }

  if(stream_choice_arg == "ALL") {
    josh_plots <- obsv_data |>
      filter(stream %in% stream_choice) |>
      #filter(!stream %in% c("feather river", "butte creek", "yuba river")) |>
      drop_na(obsv_upstream, obsv_spawner_count) |>
      left_join(b1_effect |>
                  select(year, stream, b1_effect),
                by = c("year", "stream")) |>
      left_join(rps |>
                  select(stream, rps, log_rps, lcl, ucl), by = c("stream")) |>
      mutate(pred_effect = obsv_upstream * exp(log_rps + b1_effect)) |>
      ggplot(aes(x = obsv_upstream, y = obsv_spawner_count)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      geom_abline(aes(intercept = 0, slope = rps)) +
      geom_ribbon(aes(ymin = (0 + obsv_upstream * lcl),
                      ymax = (0 + obsv_upstream * ucl)),
                  fill = "grey70", alpha = 0.3) +
      geom_errorbar(aes(ymin = pmin((rps * obsv_upstream), pred_effect),
                        ymax = pmax((rps * obsv_upstream), pred_effect)),
                    color = "#FD6467") +
      facet_wrap(~stream, scale = "free") +
      theme_minimal() +
      xlab("Observed upstream passage") +
      ylab("Observed spawner count (redd or holding)") +
      ggtitle("Model diagnostics")
  } else {
    josh_plots <- obsv_data |>
      filter(stream %in% stream_choice) |>
      drop_na(obsv_upstream, obsv_spawner_count) |>
      left_join(b1_effect |>
                  select(year, stream, b1_effect),
                by = c("year", "stream")) |>
      left_join(rps |>
                  select(stream, rps, log_rps, lcl, ucl), by = c("stream")) |>
      mutate(pred_effect = obsv_upstream * exp(log_rps + b1_effect)) |>
      ggplot(aes(x = obsv_upstream, y = obsv_spawner_count)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      geom_abline(aes(intercept = 0, slope = rps)) +
      geom_ribbon(aes(ymin = (0 + obsv_upstream * lcl),
                      ymax = (0 + obsv_upstream * ucl)),
                  fill = "grey70", alpha = 0.3) +
      geom_errorbar(aes(ymin = pmin((rps * obsv_upstream), pred_effect),
                          ymax = pmax((rps * obsv_upstream), pred_effect)),
                      color = "#FD6467") +
      theme_minimal() +
      xlab("Observed upstream passage") +
      ylab("Observed spawner count (redd or holding)") +
      ggtitle("Model diagnostics")
  }

  ggsave(here::here("data-raw", "adult_model", "adult_model_plots",
                    paste0("josh_plots_", stream_choice_arg, ".png")),
         josh_plots,
         width = 8, height = 6)
}

create_josh_plot(obsv_data, b1_effect, rps, annual_random_effects,
                 "ALL")




# plot conversion rates ---------------------------------------------------
conversion_rates <- SRJPEmodel::P2S_model_fits |>
  filter(str_detect(par_names, "conversion_rate")) |>
  mutate(lcl = `2.5%`, ucl = `97.5%`) |>
  select(par_names, median = `50%`, sd, lcl, ucl, stream, year) |>
  glimpse()

conversion_rates_with_year <- conversion_rates |>
  select(-c(par_names)) |>
  rename(pred_conversion_rate = median) |>
  left_join(full_data_for_input |>
              mutate(covar_used = wy_type,
                     covar_type = "Water year type") |>
              # mutate(covar_used = ifelse(stream == "deer creek", max_flow_std, gdd_std),
              #        covar_type = ifelse(stream == "deer creek", "Maximum flow", "Temperature index")) |>
              select(year, stream, covar_used, covar_type),
            by = c("year", "stream")) |>
  mutate(stream_capital = str_to_title(stream),
         stream_label = ifelse(stream_capital == "Deer Creek", paste0(stream_capital, "- Holding"), paste(stream_capital, "- Redd")),
         water_year_type_label = ifelse(covar_used == 1, "Wet", "Dry")) |>
  glimpse()

plot_conversion_rate <- function(conversion_rates_with_year,
                                 stream_choice_arg = c("ALL", "battle creek", "clear creek",
                                                       "deer creek", "mill creek")) {

  if(stream_choice_arg == "ALL") {
    stream_choice <- c("battle creek", "clear creek",
                       "deer creek", "mill creek")
  } else {
    stream_choice <- stream_choice_arg
  }

  if(stream_choice_arg == "ALL") {
    conversion_rate_plot <- conversion_rates_with_year |>
      filter(stream %in% stream_choice) |>
      ggplot(aes(x = year, y = pred_conversion_rate)) +
      geom_errorbar(aes(x = year, ymax = ucl, ymin = lcl),
                    width = 0.2, alpha = 0.7) +
      geom_line() +
      geom_point(aes(x = year, y = pred_conversion_rate,
                     color = water_year_type_label),
                 size = 4) +
      #geom_line(aes(x = year, y = covar_used, color = covar_type), linetype = "dashed") +
      geom_hline(yintercept = 1, linetype = "dotted") +
      facet_wrap(~stream_label, scales = "free") +
      theme_minimal() + xlab("Year") +
      ylab("Predicted Conversion Rate") +
      ggtitle(paste0("Conversion Rate of Passage to Spawners - ", stream_choice_arg)) +
      labs(color = "Covariate Type") +
      scale_color_manual(values = wes_palette("GrandBudapest1")[3:2]) +
      theme(plot.title = element_text(size = 15),
            legend.position = "bottom",
            strip.text = element_text(size = 12),
            axis.text = element_text(size = 12),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10))

  } else {
    conversion_rate_plot <- conversion_rates_with_year |>
      filter(stream %in% stream_choice) |>
      ggplot(aes(x = year, y = pred_conversion_rate)) +
      geom_errorbar(aes(x = year, ymax = ucl, ymin = lcl),
                    width = 0.2, alpha = 0.7) +
      geom_line() +
      geom_point(aes(x = year, y = pred_conversion_rate,
                     color = water_year_type_label),
                 size = 4) +
      #geom_line(aes(x = year, y = covar_used, color = covar_type), linetype = "dashed") +
      geom_hline(yintercept = 1, linetype = "dotted") +
      theme_minimal() + xlab("Year") +
      ylab("Predicted Conversion Rate") +
      ggtitle(paste0("Conversion Rate of Passage to Spawners - ", stream_choice_arg)) +
      labs(color = "Covariate Type") +
      scale_color_manual(values = wes_palette("GrandBudapest1")[3:2]) +
      theme(plot.title = element_text(size = 15),
            legend.position = "bottom",
            strip.text = element_text(size = 12),
            axis.text = element_text(size = 12),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10))
  }

  ggsave(filename = here::here("data-raw", "adult_model",
                             "adult_model_plots",
                             paste0("conversion_rate_plot_", stream_choice_arg, ".jpg")),
       plot = conversion_rate_plot, width = 12, height = 6)
}

plot_conversion_rate(conversion_rates_with_year,
                     stream_choice_arg = "mill creek")



# standardized conversion rate plot ---------------------------------------

conversion_rate_plot_battle <- conversion_rates_with_year |>
  filter(stream == "Battle Creek - Redd") |>
  mutate(covar_used_scaled = scales::rescale(covar_used, to = c(0, 0.5))) |>
  ggplot(aes(x = year, y = pred_conversion_rate)) +
  geom_line() +
  geom_line(aes(x = year, y = covar_used_scaled, color = covar_type), linetype = "dashed") +
  theme_minimal() + xlab("Year") +
  ylab("Predicted Conversion Rate") +
  ggtitle("Conversion Rate of Passage to Spawners\nBattle Creek") +
  scale_color_manual("Covariate type", values = wes_palette("GrandBudapest1")[2:3]) +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = "bottom",
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

ggsave(filename = here::here("data-raw", "adult_model",
                             "adult_model_plots", "conversion_rate_plot_battle.jpg"),
       plot = conversion_rate_plot_battle, width = 12, height = 7)


# forecasts ---------------------------------------------------------------

yuba_data <- SRJPEdata::observed_adult_input |>
  filter(stream == "yuba river")

butte_data <- SRJPEdata::observed_adult_input |>
  filter(stream == "butte creek")

feather_data <- SRJPEdata::observed_adult_input |>
  filter(stream == "feather river")

all_streams <- SRJPEmodel::P2S_comparison_results |>
  mutate(year_index = suppressWarnings(readr::parse_number(par_names))) |>
  glimpse()

forecasts <- all_streams |>
  filter(str_detect(par_names, "abundance_forecast")) |>
  mutate(par_names = ifelse(year_index == "1", "forecast_dry_year", "forecast_wet_year"),
         forecast_type = ifelse(par_names == "forecast_dry_year", "Dry", "Wet")) |>
  select(stream, forecast_type, adult_count = median, lcl, ucl, covar_considered, year_index) |>
  glimpse()


all_data_sources <- SRJPEmodel::P2S_model_fits |>
  get_predicted_spawners_from_P2S() |>
  select(stream, adult_count = median_predicted_spawners, lcl, ucl, year) |>
  mutate(data_type = "P2S spawners estimates") |>
  bind_rows(yuba_data |>
              filter(data_type == "upstream_estimate") |>
              rename(adult_count = count) |>
              mutate(stream = "yuba river"),
            butte_data |>
              filter(data_type == "carcass_estimate") |>
              rename(adult_count = count) |>
              mutate(stream = "butte creek"),
            feather_data |>
                filter(data_type == "carcass_estimate") |>
                rename(adult_count = count) |>
                mutate(stream = "feather river")) |>
  glimpse()

adult_data_source_plot <- all_data_sources |>
  mutate(stream = str_to_title(stream),
         fake_date = as.Date(paste0(year, "-01-01"))) |>
  ggplot(aes(x = fake_date, y = adult_count)) +
  geom_line(aes(color = data_type)) +
  labs(color = "Data Type") +
  geom_ribbon(aes(x = fake_date, ymin = lcl, ymax = ucl), alpha = 0.2) +
  facet_wrap(~ stream, scales = "free", nrow = 2) +
  theme_minimal() +
  xlab("Year") + ylab("Spawner Count") +
  ggtitle("Spawner Counts by Method") +
  scale_x_date(date_labels = "%Y") +
  scale_color_manual("Covariate type", values = wes_palette("GrandBudapest1")) +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = "bottom",
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(angle = 45))

ggsave(filename = here::here("data-raw", "adult_model",
                             "adult_model_plots", "all_adult_data_sources.jpg"),
       plot = adult_data_source_plot, width = 12, height = 7)


# unmodeled data ----------------------------------------------------------

plot_raw_spawners <- function(all_data_sources, stream_name_arg) {

  type_for_plot <- unique(all_data_sources |>
                            filter(stream == stream_name_arg) |>
                            pull(data_type))
  plot <- all_data_sources |>
    mutate(year = as.integer(year)) |>
    filter(stream == stream_name_arg) |>
    left_join(waterYearType::water_year_indices |>
                filter(location == "Sacramento Valley") |>
                mutate(wy_type = ifelse(Yr_type %in% c("Wet", "Above Normal"), "Wet", "Dry")) |>
                select(wy_type, year = WY),
              by = "year") |>
    ggplot(aes(x = year, y = adult_count)) +
    geom_line() +
    geom_point(aes(x = year, y = adult_count, color = wy_type),
               size = 4) +
    geom_errorbar(aes(x = year, ymin = lcl_90, ymax = ucl_90),
                alpha = 0.7, width = 0.3) +
    theme_minimal() +
    scale_color_manual(name = "water year type",
                       values = wes_palette("GrandBudapest1")[3:2]) +
    ggtitle(paste0("Raw Spawner Counts - ",
                   str_to_title(stream_name_arg),
                   " (", type_for_plot, ")")) +
    xlab("Year") + ylab("Spawner Count")

  ggsave(filename = here::here("data-raw", "adult_model",
                               "adult_model_plots",
                               paste0(stream_name_arg, "_raw_spawners_plot.jpg")),
         plot = plot, width = 12, height = 7)


}
plot_raw_spawners(all_data_sources, "feather river")



# alternative forecasting plot for one stream -----------------------------------------

alternative_forecast_plot <- function(forecasts, stream_name_arg) {

  if(stream_name_arg == "ALL") {
    stream_name_choice <- c("battle creek", "clear creek", "deer creek", "mill creek")
  } else {
    stream_name_choice <- stream_name_arg
  }

  forecasts_stream <- forecasts |>
    filter(stream %in% stream_name_choice,
           !adult_count %in% c(0, Inf)) |>
    rename(forecast_level = forecast_type) |>
    mutate(data_type = "forecast",
           covar_considered = case_when(covar_considered == "wy_type" ~ "WY type",
                                        covar_considered == "max_flow_std" ~ "Flow",
                                        covar_considered == "gdd_std" ~ "Temp",
                                        covar_considered == "passage_index" ~ "Total Passage",
                                        covar_considered == "null_covar" ~ "Null"),
           forecast_level = ifelse(forecast_level == "Dry", "Low", "High"),
           covar_considered_f = factor(covar_considered,
                                       levels = c("Null", "WY type", "Temp", "Flow", "Total Passage"))) |>
    filter(!(covar_considered == "Null" & forecast_level == "High"))

  if(stream_name_arg == "ALL") {
    forecast_plot <- forecasts_stream |>
      mutate(stream = str_to_title(stream)) |>
      ggplot(aes(x = forecast_level, y = adult_count)) +
      geom_errorbar(aes(x = forecast_level, ymin = lcl, ymax = ucl),
                    width = 0.3, alpha = 0.7) +
      geom_point(aes(x = forecast_level, y = adult_count,
                     color = forecast_level),
                 size = 4) +
      facet_wrap(~stream + covar_considered_f, scales = "free") +
      xlab("Forecast Level") + ylab("Predicted Spawner Count at Across-year Mean Passage") +
      scale_color_manual("Forecast Level",
                         values = wes_palette("GrandBudapest1")[2:3]) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 15),
            legend.position = "bottom",
            strip.text = element_text(size = 12),
            axis.text = element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12),
            axis.text.x = element_text(angle = 45)) +
      ggtitle(paste0("Forecasted Spawners - ", str_to_title(stream_name_arg)))
  } else {
    forecast_plot <- forecasts_stream |>
      mutate(stream = str_to_title(stream)) |>
      ggplot(aes(x = forecast_level, y = adult_count)) +
      geom_errorbar(aes(x = forecast_level, ymin = lcl, ymax = ucl),
                    width = 0.3, alpha = 0.7) +
      geom_point(aes(x = forecast_level, y = adult_count,
                     color = forecast_level),
                 size = 4) +
      facet_wrap(~covar_considered_f, nrow = 1) +
      xlab("Forecast Level") + ylab("Predicted Spawner Count at Across-year Mean Passage") +
      scale_color_manual("Forecast Level",
                         values = wes_palette("GrandBudapest1")[2:3]) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 15),
            legend.position = "bottom",
            strip.text = element_text(size = 12),
            axis.text = element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12),
            axis.text.x = element_text(angle = 45)) +
      ggtitle(paste0("Forecasted Spawners - ", str_to_title(stream_name_arg)))
  }

  return(forecast_plot)

  # ggsave(filename = here::here("data-raw", "adult_model",
  #                              "adult_model_plots",
  #                              paste0(stream_name_arg, "_alternative_forecast_plot.jpg")),
  #        plot = forecast_plot, width = 12, height = 8)
}

alternative_forecast_plot(forecasts, "battle creek")


# report tables ------------------------------------------------------------

compare <- SRJPEmodel::P2S_comparison_results |>
  mutate(covar_considered = case_when(covar_considered == "wy_type" ~ "Water year type",
                                      covar_considered == "max_flow_std" ~ "Flow",
                                      covar_considered == "gdd_std" ~ "Temp",
                                      covar_considered == "passage_index" ~ "Total Passage",
                                      covar_considered == "null_covar" ~ "Null model"),
         mean = round(mean, 3),
         sd = round(sd, 3)) |>
  mutate(report_vars = paste0(mean, " (", sd, ")"),
         stream = str_to_title(stream)) |>
  filter(par_names %in% c("R2_fixed", "mean_redds_per_spawner",
                          "sigma_redds_per_spawner", "b1_survival")) |>
  pivot_wider(id_cols = c(stream, covar_considered),
              names_from = par_names,
              values_from = report_vars) |>
  select(Stream = stream, `Covariate` = covar_considered,
         `Survival (b1)` = b1_survival,
         `Spawner per passage (mean)` = mean_redds_per_spawner,
         `Spawner per passage (sd)` = sigma_redds_per_spawner,
         `R2 (fixed)` = R2_fixed) |>
  arrange(Stream, Covariate)
  # select(Stream = stream, Parameter = par_names,
  #        `Covariate Considered` = covar_considered,
  #        `Est. Mean` = mean, `Est. SD` = sd,
  #        Converged = convergence_metric) |>
  # filter(!str_detect(Parameter, "spawner_abundance_forecast")) |>
  # arrange(Stream, Parameter)

clipr::write_clip(compare)

#trib, covariate, mean, sd, R2, LOOic, sigma redds per spawner, R2_data, R2_fixed

# forecasting plot for one stream - not used -----------------------------------------

forecast_plot <- function(forecasts, all_data_sources, stream_name_arg) {


  forecasts_stream <- forecasts |>
    filter(stream == stream_name_arg,
           !adult_count %in% c(0, Inf)) |>
    rename(forecast_level = forecast_type) |>
    mutate(year = 2021,
           data_type = "forecast",
           year = case_when(covar_considered == "wy_type" ~ 2021,
                            covar_considered == "max_flow_std" ~ 2022,
                            covar_considered == "gdd_std" ~ 2023,
                            covar_considered == "null_covar" ~ 2024),
           covar_considered = case_when(covar_considered == "wy_type" ~ "WY type",
                                        covar_considered == "max_flow_std" ~ "Flow",
                                        covar_considered == "gdd_std" ~ "Temp",
                                        covar_considered == "null_covar" ~ "Null"))

  forecast_plot <- all_data_sources |>
    filter(stream == stream_name_arg) |>
    ggplot(aes(x = year, y = adult_count)) +
    geom_line() +
    geom_point(aes(x = year, y = adult_count,
                   color = covar_considered,
                   shape = forecast_level),
               size = 4,
               alpha = 0.8,
               data = forecasts_stream) +
    geom_ribbon(aes(x = year, ymin = lcl, ymax = ucl), alpha = 0.2) +
    geom_errorbar(aes(x = year, ymin = lcl, ymax = ucl,
                      color = covar_considered,
                      linetype = forecast_level),
                  data = forecasts_stream) +
    xlab("Year") + ylab("Predicted Spawner Count") +
    scale_color_manual("Forecast type",
                       values = wes_palette("GrandBudapest1")[2:6]) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 15),
          legend.position = "bottom",
          strip.text = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          axis.text.x = element_text(angle = 45)) +
    ggtitle(paste0("Forecasted Spawners - ", str_to_title(stream_name_arg)))


  ggsave(filename = here::here("data-raw", "adult_model",
                               "adult_model_plots",
                               paste0(stream_name_arg, "_forecast_plot.jpg")),
         plot = forecast_plot, width = 12, height = 7)
}

forecast_plot(forecasts, all_data_sources, "battle creek")
