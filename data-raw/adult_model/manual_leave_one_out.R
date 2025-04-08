library(tidyverse)
library(SRJPEdata)

# helper function for manually predicting out of sample
predict_out_of_sample <- function(fit, obsv_spawners, obsv_passage, obsv_env, year) {
  log_redds_per_spawner <- rstan::summary(fit, pars = c("log_mean_redds_per_spawner"))$summary |>
    as.data.frame() |>
    select(low = `2.5%`, median = `50%`, high = `97.5%`, mean)
  b1_survival <- rstan::summary(fit, pars = c("b1_survival"))$summary |>
    as.data.frame() |>
    select(low = `2.5%`, median = `50%`, high = `97.5%`, mean)

  conversion_rate_median <- exp(log_redds_per_spawner$median + b1_survival$median * obsv_env)
  conversion_rate_mean <- exp(log_redds_per_spawner$mean + b1_survival$mean * obsv_env)
  conversion_rate_low <- exp(log_redds_per_spawner$low + b1_survival$low * obsv_env)
  conversion_rate_high <- exp(log_redds_per_spawner$high + b1_survival$high * obsv_env)

  pred_spawners_median <- obsv_passage * conversion_rate_median
  pred_spawners_mean <- obsv_passage * conversion_rate_mean
  pred_spawners_low <- obsv_passage * conversion_rate_low
  pred_spawners_high <- obsv_passage * conversion_rate_high

  diff_median <- obsv_spawners - pred_spawners_median
  diff_mean <- obsv_spawners - pred_spawners_median

  results <- tibble(pred_spawners_median = pred_spawners_median,
                    pred_spawners_mean = pred_spawners_mean,
                    pred_spawners_low = pred_spawners_low,
                    pred_spawners_high = pred_spawners_high,
                    diff_median = diff_median,
                    diff_mean = diff_median,
                    year = year)

  return(results)
}


# function for workflow
manual_LOO <- function(stream, covariate) {
  # prepare full inputs
  inputs <- prepare_P2S_inputs(stream, covariate, FALSE)
  results_tibble <- list()

  for(i in 1:inputs$N) {

    # drop a year
    new_inputs <- inputs
    new_inputs$N <- inputs$N - 1
    new_inputs$input_years <- inputs$input_years[-i]
    new_inputs$observed_passage <- inputs$observed_passage[-i]
    new_inputs$observed_spawners <- inputs$observed_spawners[-i]
    new_inputs$environmental_covar <- inputs$environmental_covar[-i]

    fit <- fit_passage_to_spawner_model(new_inputs)

    results_tibble[[i]] <- predict_out_of_sample(fit, inputs$observed_spawners[i],
                                                 inputs$observed_passage[i],
                                                 inputs$environmental_covar[i],
                                                 inputs$input_years[i])




  }

  final_results <- results_tibble |>
    bind_rows() |>
    mutate(stream = stream,
           covariate = covariate,
           observed_spawners = inputs$observed_spawners)

  return(final_results)
}

test <- manual_LOO("battle creek", "wy_type")
saveRDS(test, "~/Downloads/battle_wy_type_manual_loo_results_04-05-2025.rds")

test |>
  mutate(observed_spawners = inputs$observed_spawners) |>
  ggplot() +
  geom_point(aes(x = year, y = pred_spawners_median), color = "red") +
  geom_errorbar(aes(x = year, ymin = pred_spawners_low,
                    ymax = pred_spawners_high),
                color = "red", width = 0.3) +
  geom_point(aes(x = year, y = observed_spawners))

# now run for all
to_run <- expand.grid(stream = c("battle creek", "clear creek"),
                      covar = c("wy_type", "max_flow_std", "gdd_std", "passage_index",
                                "median_passage_timing_std", "null_covar"))

all_results <- purrr::pmap(list(to_run$stream,
                                to_run$covar),
                           manual_LOO,
                           .progress = TRUE) |>
  bind_rows()

saveRDS(all_results, "~/Downloads/all_p2s_manual_loo_results_04-08-2025.rds")

all_results |>
  mutate(covariate = case_when(covariate == "wy_type" ~ "Water year type",
                               covariate == "max_flow_std" ~ "Maximum flow",
                               covariate == "gdd_std" ~ "Growing degree days",
                               covariate == "passage_index" ~ "Passage index",
                               covariate == "median_passage_timing_std" ~ "Passage timing",
                               covariate == "null_covar" ~ "Null",
                               TRUE ~ covariate),
         stream = str_to_title(stream),
         year = factor(year)) |> # for x axis ticks
  arrange(stream, covariate, year) |>
  ggplot() +
  geom_point(aes(x = year, y = pred_spawners_median,
                 color = "Predicted spawners")) +
  geom_errorbar(aes(x = year, ymin = pred_spawners_low,
                    ymax = pred_spawners_high,
                    color = "Predicted spawners"), width = 0.3) +
  geom_point(aes(x = year, y = observed_spawners,
                 color = "Observed spawners"),
             shape = 1) +
  facet_wrap(~covariate + stream, nrow = 2 ,scales = "free") +
  theme_minimal() +
  labs(x = "Spawner count",
       y = "Year") +
  scale_color_manual(name = "",
                     breaks = c("Predicted spawners", "Observed spawners"),
                     values = c("Predicted spawners" = "black",
                                "Observed spawners" = "red")) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1))

all_results |>
  filter(covariate == "wy_type") |>
  mutate(covariate == "Water year type",
         stream = str_to_title(stream),
         year = factor(year)) |>
  ggplot() +
  geom_point(aes(x = year, y = pred_spawners_median,
                 color = "Predicted spawners")) +
  geom_errorbar(aes(x = year, ymin = pred_spawners_low,
                    ymax = pred_spawners_high,
                    color = "Predicted spawners"), width = 0.1) +
  geom_point(aes(x = year, y = observed_spawners,
                 color = "Observed spawners"),
             shape = 1) +
  facet_wrap(~stream, nrow = 2 ,scales = "free") +
  theme_minimal() +
  labs(x = "Spawner count",
       y = "Year") +
  scale_color_manual(name = "",
                     breaks = c("Predicted spawners", "Observed spawners"),
                     values = c("Predicted spawners" = "black",
                                "Observed spawners" = "red")) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1))






