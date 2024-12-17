# LOO analysis for P2S
library(loo)
library(SRJPEdata)
library(SRJPEmodel)
library(tidyverse)


# prep data
data <- SRJPEdata::observed_adult_input |>
  filter(stream %in% c("battle creek", "clear creek")) |>
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
  select(year, stream, wy_type, redd_count, upstream_estimate) |>
  drop_na(redd_count, upstream_estimate, wy_type)

# get sample sizes
data |>
  group_by(stream) |>
  tally()

# run battle
battle_data <- data |>
  filter(stream == "battle creek")

# calculate ss tot inits
ss_tot <- 0
N <- length(battle_data$year)
for(i in 1:N) {
  ss_tot <- ss_tot + (battle_data$redd_count[i] - mean(battle_data$redd_count))^2
}

battle_data_list <- list("N" = length(unique(battle_data$year)),
                        "input_years" = unique(battle_data$year),
                        "observed_passage" = battle_data$upstream_estimate,
                        "observed_spawners" = battle_data$redd_count,
                        "percent_female" = 0.5,
                        "environmental_covar" = battle_data$wy_type,
                        "ss_total" = ss_tot,
                        "average_upstream_passage" = mean(battle_data$upstream_estimate, na.rm = TRUE))

p2s_model <- eval(parse(text = "SRJPEmodel::p2s_model_code"))

battle_fit <- rstan::stan(model_name = "battle_p2s",
                          model_code = p2s_model,
                          data = battle_data_list,
                          chains = 3, iter = 20000*2, seed = 84735)

log_lik_battle <- extract_log_lik(battle_fit, merge_chains = FALSE)
r_eff_battle <- relative_eff(exp(log_lik_battle), cores = 2)

loo_battle <- loo(log_lik_battle, r_eff = r_eff_battle, cores = 2)
print(loo_battle) # Pareto k diagnostics are bad, but not very bad

pointwise <- tibble("year" = battle_data$year,
                    "elpd" = loo_battle$pointwise[,1],
                    "elpd_mcse" = loo_battle$pointwise[,2])

# Clear
clear_data <- data |>
  filter(stream == "clear creek",
         upstream_estimate > 0)

# calculate ss tot inits
ss_tot <- 0
N <- length(clear_data$year)
for(i in 1:N) {
  ss_tot <- ss_tot + (clear_data$redd_count[i] - mean(clear_data$redd_count))^2
}

clear_data_list <- list("N" = length(unique(clear_data$year)),
                         "input_years" = unique(clear_data$year),
                         "observed_passage" = clear_data$upstream_estimate,
                         "observed_spawners" = clear_data$redd_count,
                         "percent_female" = 0.5,
                         "environmental_covar" = clear_data$wy_type,
                         "ss_total" = ss_tot,
                         "average_upstream_passage" = mean(clear_data$upstream_estimate, na.rm = TRUE))


clear_fit <- rstan::stan(model_name = "clear_p2s",
                          model_code = p2s_model,
                          data = clear_data_list,
                          chains = 3, iter = 20000*2, seed = 84735)

log_lik_clear <- extract_log_lik(clear_fit, merge_chains = FALSE)
r_eff_clear <- relative_eff(exp(log_lik_clear), cores = 2)

loo_clear <- loo(log_lik_clear, r_eff = r_eff_clear, cores = 2)
print(loo_clear) # Pareto k diagnostics are bad, but not very bad

clear_pointwise <- tibble("year" = clear_data$year,
                          "elpd" = loo_clear$pointwise[,1],
                          "elpd_mcse" = loo_clear$pointwise[,2])


# run for all environmental covars ----------------------------------------
# this doesn't actually make statistical sense - we can only compare
# if fit to the same dataset

run_loo <- function(stream_arg) {

  data <- SRJPEdata::observed_adult_input |>
    filter(stream == stream_arg) |>
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
    drop_na(wy_type:redd_count) |>
    mutate(null = 0) |>
    relocate(null, .before = upstream_estimate)

  # calculate ss tot inits
  ss_tot <- 0
  N <- length(data$year)
  for(i in 1:N) {
    ss_tot <- ss_tot + (data$redd_count[i] - mean(data$redd_count))^2
  }

  results <- list()

  for(i in names(data)[3:8]) {
    print(i)
    covar_clean <- data[i] |>
      as.vector() |> unlist() |> unname()

    print(covar_clean)

    data_list <- list("N" = length(unique(data$year)),
                      "input_years" = unique(data$year),
                      "observed_passage" = data$upstream_estimate,
                      "observed_spawners" = data$redd_count,
                      "percent_female" = 0.5,
                      "environmental_covar" = covar_clean,
                      "ss_total" = ss_tot,
                      "average_upstream_passage" = mean(data$upstream_estimate, na.rm = TRUE))

    p2s_model <- eval(parse(text = "SRJPEmodel::p2s_model_code"))

    fit <- rstan::stan(model_name = "p2s",
                       model_code = p2s_model,
                       data = data_list,
                       chains = 3, iter = 20000*2, seed = 84735)

    log_lik <- extract_log_lik(fit, merge_chains = FALSE)
    r_eff <- relative_eff(exp(log_lik), cores = 2)

    loo <- loo(log_lik, r_eff = r_eff, cores = 2)

    results[[i]] <- loo$estimates |>
      data.frame() |>
      tibble::rownames_to_column("diagnostic") |>
      mutate(covar = i)
  }

  return(results)
}

# run for battle
battle <- run_loo("battle creek")
# run for clear
clear <- run_loo("clear creek")

# get LOOIC to compare
bind_rows(battle) |>
  filter(diagnostic == "looic") |>
  arrange(Estimate) |>
  clipr::write_clip()
bind_rows(clear) |>
  filter(diagnostic == "looic") |>
  arrange(Estimate) |>
  clipr::write_clip()

# get p_loo to compare
bind_rows(battle) |>
  filter(diagnostic == "p_loo") |>
  clipr::write_clip()
bind_rows(clear) |>
  filter(diagnostic == "p_loo") |>
  clipr::write_clip()
