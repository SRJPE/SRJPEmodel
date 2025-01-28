# code to compare environmental covariates with the STAN model
# 1-26-2024
# old code ----------------------------------------------------------------

# stan model --------------------------------------------------------------

predicted_spawners <- "
  data {
    int N;
    int input_years[N];
    int observed_passage[N];
    int observed_spawners[N]; // this is redd count but also holding count
    real percent_female;
    real environmental_covar[N];
    real ss_total;
    real average_upstream_passage;
  }

  parameters {
    real log_mean_redds_per_spawner; //no lower constraint of 0 since in log space
    real <lower = 0> sigma_redds_per_spawner;
    real b1_survival;
    real log_redds_per_spawner[N];
  }

  transformed parameters{

    real conversion_rate[N];
    real redds_per_spawner[N];
    real predicted_spawners[N];
    real mean_redds_per_spawner = exp(log_mean_redds_per_spawner) * percent_female;

    for(i in 1:N) {
      conversion_rate[i] = exp(log_redds_per_spawner[i] + b1_survival * environmental_covar[i]);

      // predicted redds is product of observed passage * conversion rate
      predicted_spawners[i] = observed_passage[i] * conversion_rate[i];
      redds_per_spawner[i] = exp(log_redds_per_spawner[i]); //for ouptut only
    }
  }

  model {
    log_redds_per_spawner ~ normal(log_mean_redds_per_spawner, sigma_redds_per_spawner);
    observed_spawners ~ poisson(predicted_spawners);
  }

  generated quantities {

    // report years for ease later
    int years[N];

    // calculate R2 fixed effects
    real y_pred[N];
    real RE_draws[N]; // vector for storing random draws for RE
    real var_fit = 0;
    real var_res = 0.0;
    real ss_res = 0.0;

    for(i in 1:N) {
      years[i] = input_years[i];
      RE_draws[i] = normal_rng(log_mean_redds_per_spawner, sigma_redds_per_spawner);
      y_pred[i] = observed_passage[i] * exp(RE_draws[i]);
    }

    real y_pred_mu = mean(y_pred[]);
    real R2_data;
    real R2_fixed;

    // calculate R2_data and R2_fixed
    for(i in 1:N) {
      var_fit = var_fit + (y_pred[i] - y_pred_mu)^2;
      // TODO var_res is actually based off y_pred with full env effects?
      // TODO and then subtract y_pred from that?
      var_res = var_res + (observed_spawners[i] - y_pred[i])^2;
      ss_res = ss_res + (observed_spawners[i] - predicted_spawners[i])^2;
    }

    R2_data = 1.0 - (ss_res / ss_total);
    R2_fixed = 1.0 - (var_fit / (var_fit + var_res));

    // forecast error in spawner abundance at average upstream passage abundance
    real spawner_abundance_forecast[2];
    real rps_forecast = normal_rng(log_mean_redds_per_spawner, sigma_redds_per_spawner);// random draw from hyper

    spawner_abundance_forecast[1] = average_upstream_passage * exp(rps_forecast); //for dummy variable = 0 condition, or at average covariate conditions
    spawner_abundance_forecast[2] = average_upstream_passage * exp(rps_forecast + b1_survival); //for dummy variable=1 condition, or at 1 SD from mean covariate condition
  }
"

# test covariates with bayesian model -------------------------------------
compare_covars <- function(data, stream_name, truncate_data) {
  final_results <- tibble("par_names" = "DELETE",
                          "stream" = "DELETE",
                          "mean" = 0.0,
                          "median" = 0.0,
                          "sd" = 0.0,
                          "covar_considered" = "DELETE")

  covars_to_consider = c("wy_type", "max_flow_std", "gdd_std",
                         "null_covar", "passage_index") # no more "median_passage_timing_std" bc of sample size

  # set percent female variable to 1 for carcass and holding surveys, 0.5 for redd
  percent_female <- case_when(stream_name %in% c("battle creek", "clear creek", "mill creek") ~ 0.5,
                              stream_name %in% c("yuba river", "feather river", "butte creek", "deer creek") ~ 1)

  # loop through covariates to test
  for(i in 1:length(covars_to_consider)) {

    selected_covariate <- covars_to_consider[i]

    if(truncate_data == TRUE) {
      # drop NAs for all covariates being considered
      # makes sure you are using the same dataset for each cycle
      stream_data <- data |>
        mutate(observed_spawners = ifelse(stream == "deer creek", holding_count, redd_count)) |>
        filter(upstream_estimate > 0) |>
        filter(stream == stream_name) |>
        drop_na(wy_type, max_flow_std, gdd_std, observed_spawners, upstream_estimate) |>
        mutate(null_covar = 0)
    } else if(truncate_data == FALSE) {
      stream_data <- data |>
        mutate(observed_spawners = ifelse(stream == "deer creek", holding_count, redd_count)) |>
        filter(upstream_estimate > 0) |>
        filter(stream == stream_name) |>
        mutate(null_covar = 0) |>
        drop_na(observed_spawners, upstream_estimate, all_of(selected_covariate))
    }

    print(selected_covariate)
    print(unique(stream_data$stream))

    covar <- stream_data |>
      pull(all_of(selected_covariate))

    stream_data_list <- list("N" = length(unique(stream_data$year)),
                             "input_years" = unique(stream_data$year),
                             "observed_passage" = stream_data$upstream_estimate,
                             "observed_spawners" = stream_data$observed_spawners,
                             "percent_female" = percent_female,
                             "environmental_covar" = covar,
                             "ss_total" = calculate_ss_tot(stream_data),
                             "average_upstream_passage" = mean(stream_data$upstream_estimate,
                                                               na.rm = TRUE))

    print(stream_data_list)

    stream_model_fit <- stan(model_code = predicted_spawners,
                             data = stream_data_list,
                             chains = 3, iter = 20000*2, seed = 84735)

    # check rhat
    rhat_diagnostic <- rhat(stream_model_fit) |> unname()
    rhat_diagnostic <- rhat_diagnostic[rhat_diagnostic != "NaN"]
    rhats_over_threshold <- sum(rhat_diagnostic > 1.05) > 0

    if(rhats_over_threshold){
      convergence_rhat <- FALSE
    } else {
      convergence_rhat <- TRUE
    }

    stream_results <- get_all_pars(stream_model_fit, stream_name) |>
      filter(par_names %in% c("R2_data", "R2_fixed", "mean_redds_per_spawner",
                              "b1_survival", "sigma_redds_per_spawner") |
               str_detect(par_names, "spawner_abundance_forecast")) |> # TODO add pspawn after it's added to STAN code
      select(par_names, stream, mean, median = `50%`, sd, lcl = `2.5%`, ucl = `97.5%`) |>
      #select(par_names, stream, mean, sd) |>
      mutate(covar_considered = selected_covariate,
             convergence_metric = convergence_rhat)

    final_results <- bind_rows(final_results, stream_results)
  }

  final_results <- final_results |>
    filter(par_names != "DELETE")

  return(final_results)
}

all_streams <- tibble(par_names = "DELETE",
                      stream = "DELETE",
                      mean = 0.0,
                      median = 0.0,
                      sd = 0.0,
                      lcl = 0.0,
                      ucl = 0.0,
                      covar_considered = "DELETE",
                      convergence_metric = NA)


for(i in c("battle creek", "clear creek", "deer creek", "mill creek")) {
  stream_results <- compare_covars(full_data_for_input, i, 84735, TRUE)
  all_streams <- bind_rows(all_streams, stream_results)
}

write.csv(all_streams |> filter(par_names != "DELETE"),
          here::here("data-raw", "adult_model", "adult_model_data",
                     "covar_compare_with_null.csv"),
          row.names = FALSE)

# identify best covariate for each trib -----------------------------------
# The best covariate will have the
# highest R2_fixed, the
# largest absolute value of b1_surv, the
# most precise b1_surv, and the
# most precise forecast error.

get_rankings <- function(results, stream_name) {
  results <- results |>
    filter(stream == stream_name) |>
    # b1_survival and spawner_forecast[2] not valid for null model
    mutate(mean = ifelse(covar_considered == "null_covar" &
                           par_names %in% c("b1_survival", "spawner_abundance_forecast[2]"),
                         NA, mean),
           sd = ifelse(covar_considered == "null_covar" &
                         par_names %in% c("b1_survival", "spawner_abundance_forecast[2]"),
                       NA, sd)) |>
    mutate(mean = case_when(covar_considered == "null_covar" &
                              par_names %in% c("b1_survival", "spawner_abundance_forecast[2]") ~ NA,
                            covar_considered != "null_covar" &
                              par_names == "spawner_abundance_forecast[1]" ~ NA,
                            TRUE ~ mean),
           sd = case_when(covar_considered == "null_covar" &
                            par_names %in% c("b1_survival", "spawner_abundance_forecast[2]") ~ NA,
                          covar_considered != "null_covar" &
                            par_names == "spawner_abundance_forecast[1]" ~ NA,
                          TRUE ~ sd)) |>
    drop_na(mean, sd) |>
    mutate(par_names = case_when(par_names == "spawner_abundance_forecast[1]" ~ "spawner_abundance_forecast",
                                 par_names == "spawner_abundance_forecast[2]" ~ "spawner_abundance_forecast",
                                 TRUE ~ par_names))

  # assign rankings
  max_rankings <- results |>
    filter(par_names %in% c("b1_survival")) |>
    #filter(par_names %in% c("R2_fixed", "b1_survival")) |>
    mutate(mean = ifelse(par_names == "R2_fixed", mean, abs(mean))) |>
    group_by(par_names) |>
    mutate(rank = order(order(mean, decreasing = TRUE))) |>
    arrange(rank)

  min_rankings <- results |>
    filter(par_names %in% c("b1_survival", "spawner_abundance_forecast")) |>
    mutate(par_names = ifelse(par_names == "b1_survival", "b1_survival_sd", par_names)) |>
    filter(mean != Inf) |>
    group_by(par_names) |>
    mutate(rank = order(order(sd, decreasing = FALSE))) |>
    arrange(rank)

  results_with_rankings <- bind_rows(max_rankings, min_rankings)

  results_with_rankings_wide <- results_with_rankings |>
    pivot_wider(id_cols = c(covar_considered, stream),
                names_from = par_names,
                values_from = rank) |>
    rowwise() |>
    mutate(mean_rank = mean(c(b1_survival, b1_survival_sd, spawner_abundance_forecast), na.rm = T)) |>
    arrange(mean_rank)

  best_covar <- results_with_rankings_wide[1, ] |>
    pull(covar_considered)

  second_best_covar <- results_with_rankings_wide[2, ] |>
    pull(covar_considered)

  return(list("best_model" = best_covar,
              "compare_rankings" = results_with_rankings_wide,
              "best_model_values" = results_with_rankings |>
                filter(covar_considered == best_covar) |>
                select(-c(stream)),
              "second_best_model_values" = results_with_rankings |>
                filter(covar_considered == second_best_covar) |>
                select(-c(stream))))
}

battle_rankings <- all_streams |> get_rankings("battle creek") # gdd_std
clear_rankings <- all_streams |> get_rankings("clear creek") # gdd_std model
deer_rankings <- all_streams |> get_rankings("deer creek") # null model - but doesn't converge. second best is max_flow_std
mill_rankings <- all_streams |> get_rankings("mill creek") # gdd_std
