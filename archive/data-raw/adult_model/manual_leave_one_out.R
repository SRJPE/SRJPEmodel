library(tidyverse)
library(SRJPEdata)

# helper function for manually predicting out of sample
predict_out_of_sample <- function(fit, obsv_spawners, obsv_passage, obsv_env, year) {

  # check for fit
  error <- ifelse(sum(rstan::summary(fit)$summary[,"Rhat"] > 1.05, na.rm = T) > 0,
                  TRUE, FALSE)

  posteriors <- as.data.frame(fit, pars = c("log_mean_redds_per_spawner", "b1_survival"))

  conversion_rate <- exp(posteriors$log_mean_redds_per_spawner + posteriors$b1_survival * obsv_env)
  conversion_rate_low <- quantile(conversion_rate, 0.025)
  conversion_rate_high <- quantile(conversion_rate, 0.975)

  pred_spawners <- obsv_passage * conversion_rate
  pred_spawners_low <- quantile(pred_spawners, 0.025)
  pred_spawners_median <- quantile(pred_spawners, 0.5)
  pred_spawners_high <- quantile(pred_spawners, 0.975)

  mean_absolute_difference <- mean(abs(pred_spawners - obsv_spawners))
  mean_absolute_difference_rel <- mean(abs(pred_spawners - obsv_spawners) / obsv_spawners)

  results <- tibble(mean_absolute_difference = mean_absolute_difference,
                    mean_absolute_difference_rel = mean_absolute_difference_rel,
                    pred_spawners_median = pred_spawners_median,
                    pred_spawners_low = pred_spawners_low,
                    pred_spawners_high = pred_spawners_high,
                    year = year,
                    rhat_error = error)

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

# now run for all
to_run <- expand.grid(covar = c("wy_type", "max_flow_std", "gdd_std", "passage_index",
                                "median_passage_timing_std", "null_covar"),
                      stream = c("battle creek", "clear creek"))

all_results <- purrr::pmap(list(to_run$stream,
                                to_run$covar),
                           manual_LOO,
                           .progress = TRUE) |>
  bind_rows()

saveRDS(all_results, here::here("data-raw", "adult_model", "all_p2s_manual_loo_results_04-08-2025.rds"))


# process results ---------------------------------------------------------

all_results <- readRDS(here::here("data-raw", "adult_model", "all_p2s_manual_loo_results_04-08-2025.rds"))

all_results |>
  group_by(stream, covariate) |>
  summarise(MAE = median(mean_absolute_difference)) |>
  ungroup() |>
  mutate(MAE = abs(round(MAE, 4))) |>
  arrange(stream, MAE) |>
  View()

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
  facet_wrap(~covariate + stream, nrow = 2 ,scales = "free_y") +
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

all_results <- all_results |>
  mutate(mean_absolute_difference = abs)

