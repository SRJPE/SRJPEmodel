# Function to create diagnostic plots for Bayesian models from stanfit objects
create_diagnostic_plots <- function(stanfit, observed_data, n_posterior_samples = 100) {
  library(rstan)
  library(ggplot2)
  library(bayesplot)
  library(dplyr)
  library(tidyr)

  pp_plots <- function(model_type, inputs, fit, n_posterior_samples = 100) {
    model_type <- "p2s"

    # Extract posterior samples for predicted values
    # Assuming your generated quantities block has R_pred for predictions
    posterior_samples <- rstan::extract(fit)

    if(model_type == "p2s") {
      pred_variable <- "predicted_spawners"
      observed_data <- inputs$observed_spawners
      key_params <- c("b_survival", "conversion_rate")
    } else if(model_type == "stock_recruit") {
      pred_variable <- "pred_R"
      observed_data <- inputs$data$R
      key_params <- c("alpha", "beta", "gamma")
    } else if(model_type == "survival") {
      pred_variable <- "pred_surv"
      observed_data <- NULL # TODO
      key_params <- NULL # TODO
    } else if(model_type == "bt_spas_x") {
      pred_variable <- "Ntot"
      observed_data <- NULL # TODO fix to be the lincoln-peterson expanded ?
      key_params <- NULL # TODO
    }

    # extract predicted variable
    extract_preds <- eval(parse(text = paste0("posterior_samples$", pred_variable)))
    y_rep <- extract_preds[sample(nrow(extract_preds), n_posterior_samples), ]

    # Basic posterior predictive check
    ppc_plot <- ppc_dens_overlay(observed_data, y_rep) +
      theme_minimal() +
      labs(
        title = "Posterior Predictive Check",
        subtitle = "Distribution of observed data vs. posterior predictions",
        x = pred_variable,
        y = "Density"
      )

    # 2. Observed vs Predicted Plot
    # -----------------------------
    # Get the mean prediction for each observation
    pred_mean <- colMeans(extract_preds)

    # Create dataframe for prediction vs observed plot
    obs_pred_df <- data.frame(
      Observed = observed_data,
      Predicted = pred_mean
    )

    # Plot observed vs predicted with 1:1 line
    obs_pred_plot <- ggplot(obs_pred_df, aes(x = Observed, y = Predicted)) +
      geom_point(alpha = 0.7) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      theme_minimal() +
      labs(
        title = "Observed vs. Predicted",
        subtitle = "Points closer to the line indicate better fit",
        x = "Observed",
        y = "Predicted"
      )

    # 3. Residual Plot
    # ----------------
    residuals_df <- obs_pred_df %>%
      mutate(
        Residual = Observed - Predicted,
        Order = 1:n()
      )

    residual_plot <- ggplot(residuals_df, aes(x = Predicted, y = Residual)) +
      geom_point(alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      theme_minimal() +
      labs(
        title = "Residual Plot",
        subtitle = "Should show random scatter around zero",
        x = "Predicted",
        y = "Residual (Observed - Predicted)"
      )

    # 4. Parameter Posterior Distributions
    # -----------------------------------
    # Extract posterior for key parameters

    samples <- posterior_samples |>
      select(all_of(key_params))

    key_param_samples_1 <- eval(parse(text = "posterior_samples$", key_params[1]))
    sigma_samples <- posterior_samples$sigma

    # Create dataframe for parameter posteriors
    params_df <- data.frame(
      alpha = alpha_samples,
      sigma = sigma_samples
    ) %>%
      pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Value")

    param_plot <- ggplot(params_df, aes(x = Value)) +
      geom_density(fill = "skyblue", alpha = 0.5) +
      facet_wrap(~ Parameter, scales = "free") +
      theme_minimal() +
      labs(
        title = "Posterior Distributions",
        subtitle = "Key model parameters",
        x = "Parameter Value",
        y = "Density"
      )

    # 5. Trace plots for convergence
    # -----------------------------
    trace_plot <- mcmc_trace(
      as.array(stanfit),
      pars = c("alpha", "sigma"),
      facet_args = list(ncol = 1)
    ) +
      labs(
        title = "MCMC Trace Plots",
        subtitle = "Should show good mixing and convergence"
      )

    # Combine plots into a list
    plots <- list(
      posterior_predictive = ppc_plot,
      observed_vs_predicted = obs_pred_plot,
      residuals = residual_plot,
      parameter_posteriors = param_plot,
      trace_plots = trace_plot
    )

    return(plots)

  }


  # 1. Posterior Predictive Check Plot
  # ----------------------------------
  # Get a subset of posterior predictions for plotting
  y_rep <- posterior_samples$R_pred[sample(nrow(posterior_samples$R_pred), n_posterior_samples), ]

  # Basic posterior predictive check
  ppc_plot <- ppc_dens_overlay(observed_data, y_rep) +
    theme_minimal() +
    labs(
      title = "Posterior Predictive Check",
      subtitle = "Distribution of observed data vs. posterior predictions",
      x = "Recruitment",
      y = "Density"
    )

  # 2. Observed vs Predicted Plot
  # -----------------------------
  # Get the mean prediction for each observation
  pred_mean <- colMeans(posterior_samples$R_pred)

  # Create dataframe for prediction vs observed plot
  obs_pred_df <- data.frame(
    Observed = observed_data,
    Predicted = pred_mean
  )

  # Plot observed vs predicted with 1:1 line
  obs_pred_plot <- ggplot(obs_pred_df, aes(x = Observed, y = Predicted)) +
    geom_point(alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(
      title = "Observed vs. Predicted",
      subtitle = "Points closer to the line indicate better fit",
      x = "Observed Recruitment",
      y = "Predicted Recruitment"
    )

  # 3. Residual Plot
  # ----------------
  residuals_df <- obs_pred_df %>%
    mutate(
      Residual = Observed - Predicted,
      Order = 1:n()
    )

  residual_plot <- ggplot(residuals_df, aes(x = Predicted, y = Residual)) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(
      title = "Residual Plot",
      subtitle = "Should show random scatter around zero",
      x = "Predicted Recruitment",
      y = "Residual (Observed - Predicted)"
    )

  # 4. Parameter Posterior Distributions
  # -----------------------------------
  # Extract posterior for key parameters
  alpha_samples <- posterior_samples$alpha
  sigma_samples <- posterior_samples$sigma

  # Create dataframe for parameter posteriors
  params_df <- data.frame(
    alpha = alpha_samples,
    sigma = sigma_samples
  ) %>%
    pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Value")

  param_plot <- ggplot(params_df, aes(x = Value)) +
    geom_density(fill = "skyblue", alpha = 0.5) +
    facet_wrap(~ Parameter, scales = "free") +
    theme_minimal() +
    labs(
      title = "Posterior Distributions",
      subtitle = "Key model parameters",
      x = "Parameter Value",
      y = "Density"
    )

  # 5. Trace plots for convergence
  # -----------------------------
  trace_plot <- mcmc_trace(
    as.array(stanfit),
    pars = c("alpha", "sigma"),
    facet_args = list(ncol = 1)
  ) +
    labs(
      title = "MCMC Trace Plots",
      subtitle = "Should show good mixing and convergence"
    )

  # Combine plots into a list
  plots <- list(
    posterior_predictive = ppc_plot,
    observed_vs_predicted = obs_pred_plot,
    residuals = residual_plot,
    parameter_posteriors = param_plot,
    trace_plots = trace_plot
  )

  return(plots)
}

# Example usage
# ------------
# Assuming you have:
# - stanfit = your fitted Stan model object
# - R_observed = vector of observed recruitment values

# Run the function
# diagnostic_plots <- create_diagnostic_plots(stanfit, R_observed)

# Display individual plots
# diagnostic_plots$posterior_predictive
# diagnostic_plots$observed_vs_predicted
# diagnostic_plots$residuals
# diagnostic_plots$parameter_posteriors
# diagnostic_plots$trace_plots

# Or display all plots together (requires gridExtra)
# grid.arrange(
#   diagnostic_plots$posterior_predictive,
#   diagnostic_plots$observed_vs_predicted,
#   diagnostic_plots$residuals,
#   diagnostic_plots$parameter_posteriors,
#   ncol = 2
# )
