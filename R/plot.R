#' Plot pCap diagnostics for all sites (tributary) model
#'
#' @description Produces three ggplot objects from a fitted `pCap_all_sites`
#'   Stan model:
#'   1. Observed vs. predicted capture probability (scatter, one point per trial).
#'   2. Mean pCap per site with 95 % CI and hyper-distribution overlay.
#'   3. Flow–pCap relationship per site (one facet each).
#'   4. Observed vs. predicted pCap histograms per site (one facet each).
#'
#' @param pcap A fitted Stan model object for the `all_sites` pCap model.
#' @param lb Numeric. Lower quantile bound. Default `0.025`.
#' @param ub Numeric. Upper quantile bound. Default `0.975`.
#'
#' @returns A named list with four ggplot objects:
#'   \describe{
#'     \item{`obs_vs_pred_plot`}{Observed vs. predicted pCap scatter.}
#'     \item{`site_mean_plot`}{Per-site mean pCap with CI and hyper-distribution.}
#'     \item{`flow_plot`}{Flow–pCap curves faceted by site.}
#'     \item{`hist_plot`}{Observed vs. predicted histograms faceted by site.}
#'   }
#'
#' @export
plot_pCap_all_sites <- function(pcap, lb = 0.025, ub = 0.975) {

  inv_logit <- function(x) exp(x) / (1 + exp(x))

  # ---- Prepare inputs -------------------------------------------------------
  pCap_inputs <- prepare_pCap_inputs(model_type = "all_sites")

  Nmr        <- pCap_inputs$inputs$data$Nmr
  Ntribs     <- pCap_inputs$inputs$data$Ntribs
  ind_trib   <- pCap_inputs$inputs$data$ind_trib
  mr_flow    <- pCap_inputs$inputs$data$mr_flow
  sites      <- pCap_inputs$sites_fit           # character vector of site names
  Releases   <- pCap_inputs$inputs$data$Releases
  Recaptures <- pCap_inputs$inputs$data$Recaptures
  pCap_obs   <- Recaptures / Releases

  # ---- 1. Observed vs predicted scatter -------------------------------------
  dp_lp <- as.data.frame(pcap, pars = "logit_pCap")
  pCap_pred <- vapply(seq_len(Nmr), function(i) {
    mean(inv_logit(dp_lp[[paste0("logit_pCap[", i, "]")]]))
  }, numeric(1))

  scatter_df <- data.frame(obs = pCap_obs, pred = pCap_pred)
  xy_lim     <- range(c(pCap_obs, pCap_pred))

  obs_vs_pred_plot <- ggplot2::ggplot(scatter_df, ggplot2::aes(x = obs, y = pred)) +
    ggplot2::geom_point(shape = 19, size = 1.5, alpha = 0.7) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::coord_equal(xlim = xy_lim, ylim = xy_lim) +
    ggplot2::labs(
      x = "Observed Capture Probability",
      y = "Predicted Capture Probability",
      title = "All Sites: Observed vs. Predicted pCap"
    ) +
    ggplot2::theme_bw()

  # ---- 2. Per-site mean pCap with CI ----------------------------------------
  dp_b0 <- as.data.frame(pcap, pars = c("b0_pCap", "pro_sd_P"))
  Nsims  <- nrow(dp_b0)

  site_stats <- lapply(seq_len(Ntribs), function(i) {
    b0   <- dp_b0[[paste0("b0_pCap[", i, "]")]]
    pC   <- inv_logit(b0)
    q    <- quantile(pC, probs = c(lb, 0.5, ub))
    irecs <- which(ind_trib == i)
    obs_mean <- sum(Recaptures[irecs]) / sum(Releases[irecs])
    n_samps  <- length(irecs)
    data.frame(
      site      = paste0(sites[i], " (", n_samps, ")"),
      site_ord  = i,
      mean_pred = mean(pC),
      lo        = q[[1]],
      hi        = q[[3]],
      obs_mean  = obs_mean
    )
  })
  site_df <- do.call(rbind, site_stats)
  site_df$site <- factor(site_df$site, levels = site_df$site[order(site_df$site_ord)])

  # Hyper-distribution summary
  dp_mu <- as.data.frame(pcap, pars = "trib_mu_P")
  dp_sd <- as.data.frame(pcap, pars = "trib_sd_P")
  hmu   <- inv_logit(mean(dp_mu[, 1]))
  hci   <- inv_logit(quantile(dp_mu[, 1], probs = c(lb, ub)))

  # Hyper-distribution density overlay (logit-normal -> probability scale)
  xmax_dens <- 0.3
  pvec       <- seq(0, xmax_dens, length.out = 200)
  lpvec      <- log(pvec / (1 - pvec))
  dens_vals  <- dnorm(lpvec, mean = mean(dp_mu[, 1]), sd = mean(dp_sd[, 1]))
  hyper_df   <- data.frame(pCap = pvec, density = dens_vals)

  site_mean_plot <- ggplot2::ggplot(site_df) +
    ggplot2::annotate(
      "rect",
      xmin = hci[[1]], xmax = hci[[2]],
      ymin = -Inf,     ymax = Inf,
      fill = "darkgray", alpha = 0.25
    ) +
    ggplot2::geom_vline(xintercept = hmu, linetype = "dashed") +
    ggplot2::geom_errorbarh(
      ggplot2::aes(y = site, xmin = lo, xmax = hi),
      height = 0.3
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = mean_pred, y = site),
      shape = 19, size = 2
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = obs_mean, y = site),
      shape = 21, colour = "red", size = 2.5
    ) +
    ggplot2::scale_x_continuous(limits = c(0, xmax_dens)) +
    ggplot2::labs(
      x     = "Mean Trap Efficiency",
      y     = NULL,
      title = "All Sites: Per-Site Capture Probability"
    ) +
    ggplot2::theme_bw()

  # ---- 3. Flow–pCap curves per site -----------------------------------------
  Qs    <- seq(-2, 6, by = 0.1)
  Nq    <- length(Qs)

  flow_rows <- lapply(seq_len(Ntribs), function(i) {
    b0nm  <- paste0("b0_pCap[", i, "]")
    bfnm  <- paste0("b_flow[", i, "]")
    dp1   <- as.data.frame(pcap, pars = b0nm)
    dp2   <- as.data.frame(pcap, pars = bfnm)
    b0    <- dp1[[b0nm]]
    bf    <- dp2[[bfnm]]

    pred_mat <- sapply(Qs, function(q) {
      p <- inv_logit(b0 + bf * q)
      c(mean = mean(p), lo = quantile(p, lb), hi = quantile(p, ub))
    })

    irecs <- which(ind_trib == i)
    obs_rows <- data.frame(
      site  = sites[i],
      flow  = mr_flow[irecs],
      pCap  = Recaptures[irecs] / Releases[irecs],
      type  = "observed"
    )
    curve_rows <- data.frame(
      site      = sites[i],
      flow      = Qs,
      mean_pred = pred_mat["mean", ],
      lo        = pred_mat["lo.2.5%", ],
      hi        = pred_mat["hi.97.5%", ]
    )
    list(obs = obs_rows, curve = curve_rows)
  })

  obs_flow_df   <- do.call(rbind, lapply(flow_rows, `[[`, "obs"))
  curve_flow_df <- do.call(rbind, lapply(flow_rows, `[[`, "curve"))

  flow_plot <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = curve_flow_df,
      ggplot2::aes(x = flow, ymin = lo, ymax = hi, group = site),
      alpha = 0.25
    ) +
    ggplot2::geom_line(
      data = curve_flow_df,
      ggplot2::aes(x = flow, y = mean_pred, colour = site),
      linewidth = 1
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    ggplot2::geom_point(
      data = obs_flow_df,
      ggplot2::aes(x = flow, y = pCap, colour = site),
      shape = 19, size = 1.5, alpha = 0.8
    ) +
    ggplot2::facet_wrap(~ site, scales = "free_y") +
    ggplot2::labs(
      x     = "Standardized Discharge",
      y     = "Capture Probability",
      title = "All Sites: Flow–pCap Relationships"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")

  # ---- 4. Observed vs predicted histograms per site -------------------------
  tp <- 0.25
  dp_hist <- as.data.frame(pcap, pars = c("b0_pCap", "yr_sd_P", "pro_sd_P"))
  Nsims_h <- nrow(dp_hist)

  hist_rows <- lapply(seq_len(Ntribs), function(i) {
    b0      <- dp_hist[[paste0("b0_pCap[",  i, "]")]]
    yr_sd   <- dp_hist[[paste0("yr_sd_P[",  i, "]")]]
    pro_sd  <- dp_hist[[paste0("pro_sd_P[", i, "]")]]
    pdev    <- rnorm(Nsims_h, 0, pro_sd)
    pC      <- inv_logit(b0 + pdev)

    irecs   <- which(ind_trib == i)
    pCap_o  <- pCap_obs[irecs]
    xmax    <- max(pCap_o)

    obs_cut <- pCap_o[pCap_o <= xmax]
    gen_cut <- pC[pC <= xmax]

    data.frame(
      site   = sites[i],
      value  = c(obs_cut, gen_cut),
      source = c(rep("Observed",  length(obs_cut)),
                 rep("Predicted", length(gen_cut)))
    )
  })
  hist_df <- do.call(rbind, hist_rows)

  mean_rows <- lapply(seq_len(Ntribs), function(i) {
    irecs <- which(ind_trib == i)
    b0    <- dp_hist[[paste0("b0_pCap[",  i, "]")]]
    pro_sd <- dp_hist[[paste0("pro_sd_P[", i, "]")]]
    pdev  <- rnorm(Nsims_h, 0, pro_sd)
    pC    <- inv_logit(b0 + pdev)
    data.frame(
      site   = sites[i],
      source = c("Observed", "Predicted"),
      xint   = c(mean(pCap_obs[irecs]), mean(pC))
    )
  })
  mean_df <- do.call(rbind, mean_rows)

  hist_plot <- ggplot2::ggplot(hist_df, ggplot2::aes(x = value, fill = source)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      position = "identity", alpha = 0.4, bins = 40
    ) +
    ggplot2::geom_vline(
      data     = mean_df,
      ggplot2::aes(xintercept = xint, colour = source,
                   linetype = source),
      linewidth = 0.8
    ) +
    ggplot2::scale_fill_manual(values   = c("Observed" = "red",  "Predicted" = "blue")) +
    ggplot2::scale_colour_manual(values = c("Observed" = "red",  "Predicted" = "blue")) +
    ggplot2::scale_linetype_manual(values = c("Observed" = "solid", "Predicted" = "dashed")) +
    ggplot2::facet_wrap(~ site, scales = "free") +
    ggplot2::labs(
      x      = "Capture Probability",
      y      = "Density",
      title  = "All Sites: Observed vs. Predicted pCap Distributions",
      fill   = NULL,
      colour = NULL,
      linetype = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top")

  list(
    obs_vs_pred_plot = obs_vs_pred_plot,
    site_mean_plot   = site_mean_plot,
    flow_plot        = flow_plot,
    hist_plot        = hist_plot
  )
}

#' Plot pCap flow relationship and predicted vs observed for mainstem sites
#'
#' @description Produces two ggplot panels per site:
#'   1. Flow–capture probability relationship with observed and predicted points
#'      and uncertainty ribbons (with and without process error).
#'   2. Overlapping histograms of observed vs. predicted capture probability.
#'
#' @param site_name Character. One of the mainstem sites supported by the
#'   one-site skew model (e.g. `"knights landing"`, `"tisdale"`).
#' @param pcap A fitted Stan model object (loaded from an `.Rdata` file via
#'   `load()`) corresponding to `site_name`.
#' @param lb Numeric. Lower quantile bound. Default `0.025`.
#' @param ub Numeric. Upper quantile bound. Default `0.975`.
#'
#' @returns A named list with two ggplot objects:
#'   \describe{
#'     \item{`flow_plot`}{Flow–pCap relationship plot.}
#'     \item{`hist_plot`}{Observed vs. predicted pCap histogram.}
#'   }
#'
#' @export
plot_pCap_main <- function(site_name, pcap, lb = 0.025, ub = 0.975) {

  inv_logit <- function(x) exp(x) / (1 + exp(x))

  # ---- Prepare inputs -------------------------------------------------------
  pCap_inputs <- prepare_pCap_inputs(
    model_type     = "one_site",
    skew           = TRUE,
    site_selection = site_name
  )

  Releases    <- pCap_inputs$inputs$data$Releases
  Recaptures  <- pCap_inputs$inputs$data$Recaptures
  obs_pCap    <- Recaptures / Releases
  flow        <- pCap_inputs$inputs$data$mr_flow
  Nobs        <- length(obs_pCap)

  dp <- as.data.frame(pcap, pars = c("logit_pCap", "b0_pCap", "b_flow",
                                     "pro_sd_P", "alpha"))

  # ---- Predicted pCap per trial ---------------------------------------------
  pstats <- matrix(nrow = Nobs, ncol = 3)
  for (i in seq_len(Nobs)) {
    col_nm  <- paste0("logit_pCap[", i, "]")
    samples <- inv_logit(dp[[col_nm]])
    pstats[i, ] <- quantile(samples, probs = c(lb, 0.5, ub))
    pstats[i, 2] <- mean(samples)
  }

  # ---- Flow–pCap relationship -----------------------------------------------
  b0_pCap  <- dp[["b0_pCap"]]
  b_flow   <- dp[["b_flow"]]
  pro_sd_P <- dp[["pro_sd_P"]]
  alpha    <- dp[["alpha"]]

  dp_yr <- as.data.frame(pcap, pars = "yr_sd_P")
  yr_sd_P <- dp_yr[, 1]

  Nsims <- length(pro_sd_P)
  pdev  <- brms::rskew_normal(
    n = Nsims, mu = 0, sigma = pro_sd_P, alpha = alpha, xi = NULL, omega = NULL
  )

  Nq   <- 50
  qseq <- seq(min(flow) * 1.1, max(flow) * 1.1, length.out = Nq)

  flow_df <- data.frame(
    flow      = qseq,
    mean_pred = NA_real_,
    lo        = NA_real_,
    hi        = NA_real_,
    lo_pe     = NA_real_,
    hi_pe     = NA_real_
  )

  for (i in seq_len(Nq)) {
    pred    <- inv_logit(b0_pCap + b_flow * qseq[i])
    pred_pe <- inv_logit(b0_pCap + b_flow * qseq[i] + pdev)
    flow_df$mean_pred[i] <- mean(pred)
    flow_df$lo[i]        <- quantile(pred,    lb)
    flow_df$hi[i]        <- quantile(pred,    ub)
    flow_df$lo_pe[i]     <- quantile(pred_pe, lb)
    flow_df$hi_pe[i]     <- quantile(pred_pe, ub)
  }

  obs_df <- data.frame(
    flow       = flow,
    obs_pCap   = obs_pCap,
    pred_mean  = pstats[, 2],
    pred_lo    = pstats[, 1],
    pred_hi    = pstats[, 3]
  )

  flow_plot <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = flow_df,
      ggplot2::aes(x = flow, ymin = lo_pe, ymax = hi_pe),
      fill = "grey", alpha = 0.3
    ) +
    ggplot2::geom_ribbon(
      data = flow_df,
      ggplot2::aes(x = flow, ymin = lo, ymax = hi),
      fill = "grey", alpha = 0.6
    ) +
    ggplot2::geom_line(
      data = flow_df,
      ggplot2::aes(x = flow, y = mean_pred),
      linewidth = 1
    ) +
    ggplot2::geom_point(
      data = obs_df,
      ggplot2::aes(x = flow, y = obs_pCap),
      shape = 21, colour = "blue", size = 2.5
    ) +
    ggplot2::geom_point(
      data = obs_df,
      ggplot2::aes(x = flow, y = pred_mean),
      shape = 19, size = 1.5
    ) +
    ggplot2::geom_linerange(
      data = obs_df,
      ggplot2::aes(x = flow, ymin = pred_lo, ymax = pred_hi),
      linewidth = 0.5
    ) +
    ggplot2::labs(
      title = site_name,
      x     = "Standardized Discharge",
      y     = "Capture Probability"
    ) +
    ggplot2::theme_bw()

  # ---- Histogram: observed vs predicted -------------------------------------
  dp2 <- as.data.frame(pcap)

  b0_h   <- dp2[["b0_pCap"]]
  pro_h  <- dp2[["pro_sd_P"]]
  alpha_h <- dp2[["alpha"]]
  Nsims2  <- nrow(dp2)

  yr_sd_h <- as.data.frame(pcap, pars = "yr_sd_P")[, 1]
  yrdev   <- rnorm(Nsims2, mean = 0, sd = yr_sd_h)
  pdev_h  <- brms::rskew_normal(
    n = Nsims2, mu = 0, sigma = pro_h, alpha = alpha_h, xi = NULL, omega = NULL
  )
  genP_wyr <- inv_logit(b0_h + yrdev + pdev_h)

  xmax    <- 0.05
  obs_cut <- obs_pCap[obs_pCap <= xmax]
  gen_cut <- genP_wyr[genP_wyr <= xmax]

  hist_df <- data.frame(
    value  = c(obs_cut, gen_cut),
    source = c(rep("Observed", length(obs_cut)), rep("Predicted", length(gen_cut)))
  )

  mean_df <- data.frame(
    source = c("Observed", "Predicted"),
    xint   = c(mean(obs_pCap), mean(genP_wyr))
  )

  hist_plot <- ggplot2::ggplot(hist_df, ggplot2::aes(x = value, fill = source)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      position = "identity", alpha = 0.4, bins = 40
    ) +
    ggplot2::geom_vline(
      data     = mean_df,
      ggplot2::aes(xintercept = xint, colour = source),
      linetype = c("solid", "dashed"),
      linewidth = 1
    ) +
    ggplot2::scale_fill_manual(values   = c("Observed" = "red",  "Predicted" = "blue")) +
    ggplot2::scale_colour_manual(values = c("Observed" = "red",  "Predicted" = "blue")) +
    ggplot2::labs(
      title = site_name,
      x     = "Capture Probability",
      y     = "Density",
      fill  = NULL,
      colour = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top")

  list(flow_plot = flow_plot, hist_plot = hist_plot)
}

#' Plot weekly abundance and capture probability for a site–year combination
#'
#' @description Produces two ggplot panels for a given site and run year:
#'   1. Weekly estimated abundance (bars) with 95 % CI, Lincoln-Peterson point
#'      estimates, discharge overlay, and catch labels.
#'   2. Weekly estimated capture probability (bars) with 95 % CI, observed
#'      efficiency trials, and discharge overlay.
#'
#' @param site Character. Site name (e.g. `"ucc"`, `"knights landing"`).
#' @param run_year Numeric. Run year to plot.
#' @param pcap A fitted Stan model object for the pCap model appropriate for
#'   `site`.  The model type is inferred automatically:
#'   `"one_site"` (skew) for mainstem sites (`"tisdale"`, `"knights landing"`),
#'   `"all_sites"` for all other sites.
#' @param pCap_inputs Optional. Pre-computed output of `prepare_pCap_inputs()`
#'   for the `all_sites` model.  If `NULL` (default) and `site` is not a
#'   mainstem site, the function calls `prepare_pCap_inputs()` internally.
#' @param abundance A fitted BUGS abundance model object (from
#'   `fit_abundance_model_BUGS()` or loaded from `.Rdata`).
#' @param min_pCap_mult Numeric. Multiplier passed to `prepare_abundance_inputs()`.
#'   Default `0.5`.
#' @param Nscale Numeric. Abundance scaling factor (result displayed in units of
#'   `1 / Nscale` fish). Default `0.001` (thousands of fish).
#' @param lb Numeric. Lower quantile bound. Default `0.025`.
#' @param ub Numeric. Upper quantile bound. Default `0.975`.
#'
#' @returns A named list with two ggplot objects:
#'   \describe{
#'     \item{`abundance_plot`}{Weekly abundance.}
#'     \item{`pcap_plot`}{Weekly capture probability.}
#'   }
#'
#' @export
plot_pCap_site_year <- function(site,
                                run_year,
                                pcap,
                                pCap_inputs  = NULL,
                                abundance,
                                min_pCap_mult = 0.5,
                                Nscale        = 0.001,
                                lb            = 0.025,
                                ub            = 0.975) {

  mainstem_sites <- c("tisdale", "knights landing")

  # Determine model type and pCap_inputs
  if (site %in% mainstem_sites) {
    model_type  <- "one_site"
    pCap_inputs <- prepare_pCap_inputs(
      model_type     = "one_site",
      skew           = TRUE,
      site_selection = site
    )
  } else {
    model_type <- "all_sites"
    if (is.null(pCap_inputs)) {
      pCap_inputs <- prepare_pCap_inputs(model_type = "all_sites")
    }
  }

  abundance_inputs <- prepare_abundance_inputs(
    site              = site,
    run_year          = run_year,
    pCap_model_type   = model_type,
    min_pCap_mult     = min_pCap_mult,
    pCap_model_object = pcap
  )

  Nstrata  <- abundance_inputs$inputs$data$Nstrata
  Jdt      <- abundance_inputs$week_date        # date labels per stratum
  Uwc_ind  <- abundance_inputs$inputs$data$Uwc_ind
  u        <- abundance_inputs$inputs$data$u
  effort   <- abundance_inputs$inputs$data$effort
  Flow     <- abundance_inputs$catch_flow / 1000 # kcfs

  # Which strata have efficiency data?
  irecs_N <- which(
    !is.na(abundance_inputs$lp_data$number_released) &
      !is.na(abundance_inputs$lp_data$number_recaptured)
  )

  Releases   <- NA
  Recaptures <- NA
  pN         <- numeric(0)

  if (length(irecs_N) > 0) {
    Releases   <- abundance_inputs$lp_data$number_released[irecs_N]
    Recaptures <- abundance_inputs$lp_data$number_recaptured[irecs_N]
    irecs_u    <- which(!is.na(match(Uwc_ind, irecs_N)))
    pN         <- (u[irecs_u] / ((Recaptures + 0.01) / Releases)) * Nscale
    pN[Recaptures == 0] <- NA
  }

  # ---- Weekly N estimates ---------------------------------------------------
  N_sims <- as.data.frame(abundance$sims.list$N)
  N_mat  <- matrix(nrow = Nstrata, ncol = 3)
  for (i in seq_len(Nstrata)) {
    v           <- N_sims[, i] * Nscale
    N_mat[i, ]  <- quantile(v, probs = c(lb, 0.5, ub))
    N_mat[i, 2] <- mean(v)
  }

  # Total N label for title
  Ntot_sims <- abundance$sims.list$Ntot
  Ntot       <- quantile(Ntot_sims, probs = c(lb, 0.5, ub)) * Nscale
  Ntot_cv    <- round(100 * sd(Ntot_sims) / mean(Ntot_sims), 0)
  NtotLab    <- paste0(
    round(Ntot[2], 0), " (", round(Ntot[1], 0), " – ", round(Ntot[3], 0),
    ") cv=", Ntot_cv, "%"
  )

  bar_colour <- ifelse(seq_len(Nstrata) %in% Uwc_ind, "grey70", "#e05555")

  abund_df <- data.frame(
    stratum    = seq_len(Nstrata),
    date_label = factor(Jdt, levels = Jdt),
    mean_N     = N_mat[, 2],
    lo_N       = N_mat[, 1],
    hi_N       = N_mat[, 3],
    bar_colour = bar_colour,
    flow       = Flow,
    sampled    = seq_len(Nstrata) %in% Uwc_ind
  )

  # Catch counts per sampled stratum (same order as Uwc_ind)
  abund_df$catch_label <- NA_character_
  abund_df$effort_label <- NA_character_
  k <- 0
  for (i in seq_len(Nstrata)) {
    if (abund_df$sampled[i]) {
      k <- k + 1
      abund_df$catch_label[i]  <- as.character(u[k])
      abund_df$effort_label[i] <- as.character(round(effort[i], 2))
    }
  }

  # Lincoln-Peterson abundance points
  lp_df <- if (length(irecs_N) > 0) {
    data.frame(stratum = irecs_N, pN = pN)
  } else {
    data.frame(stratum = integer(0), pN = numeric(0))
  }

  # Normalise flow for secondary axis
  flow_range <- range(Flow, na.rm = TRUE)
  N_range    <- range(c(N_mat, pN), na.rm = TRUE)
  scale_flow <- function(f) {
    N_range[1] + (f - flow_range[1]) /
      diff(flow_range) * diff(N_range)
  }

  abundance_plot <- ggplot2::ggplot(abund_df,
                                    ggplot2::aes(x = date_label)) +
    ggplot2::geom_bar(
      ggplot2::aes(y = mean_N, fill = sampled),
      stat = "identity", width = 0.75
    ) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "grey70", "FALSE" = "#e05555"),
      guide  = "none"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lo_N, ymax = hi_N),
      width = 0.25
    ) +
    ggplot2::geom_point(
      data = lp_df,
      ggplot2::aes(x = abund_df$date_label[stratum], y = pN),
      shape = 21, colour = "blue", size = 2.5, inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(y = max(N_mat) * 1.05, label = catch_label),
      angle = 90, hjust = 0, vjust = 0.5, size = 2.5, na.rm = TRUE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(y = max(N_mat) * 0.85, label = effort_label),
      angle = 90, hjust = 0, vjust = 0.5, size = 2.2, colour = "grey40",
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = scale_flow(flow), group = 1),
      colour = "darkgreen", linewidth = 1
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = scale_flow(flow)),
      colour = "darkgreen", size = 1.2
    ) +
    ggplot2::scale_y_continuous(
      name      = "Abundance ('000s)",
      sec.axis  = ggplot2::sec_axis(
        transform = ~ flow_range[1] + (. - N_range[1]) / diff(N_range) * diff(flow_range),
        name      = "Discharge (kcfs)"
      )
    ) +
    ggplot2::labs(
      x     = "",
      title = paste0(site, " ", run_year, "  Ntot = ", NtotLab)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x       = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      axis.title.y.right = ggplot2::element_text(colour = "darkgreen"),
      axis.text.y.right  = ggplot2::element_text(colour = "darkgreen")
    )

  # ---- Weekly pCap estimates ------------------------------------------------
  dp_pcap  <- abundance$sims.list$pCap_U
  pCap_mat <- matrix(nrow = Nstrata, ncol = 3)
  for (i in seq_len(Nstrata)) {
    irecs_min <- which(dp_pcap[, i] >= abundance_inputs$min_pCap)
    if (length(irecs_min) == 0) irecs_min <- seq_len(nrow(dp_pcap))
    v              <- dp_pcap[irecs_min, i]
    pCap_mat[i, ]  <- quantile(v, probs = c(lb, 0.5, ub))
    pCap_mat[i, 2] <- mean(v)
  }

  pcap_df <- data.frame(
    stratum    = seq_len(Nstrata),
    date_label = factor(Jdt, levels = Jdt),
    mean_pCap  = pCap_mat[, 2],
    lo_pCap    = pCap_mat[, 1],
    hi_pCap    = pCap_mat[, 3],
    has_trial  = seq_len(Nstrata) %in% irecs_N,
    flow       = Flow
  )

  # Observed efficiencies at trial strata
  obs_eff_df <- if (length(irecs_N) > 0 && !all(is.na(Releases))) {
    data.frame(
      stratum     = irecs_N,
      date_label  = factor(Jdt[irecs_N], levels = Jdt),
      obs_eff     = Recaptures / Releases,
      Recaptures  = Recaptures,
      Releases    = Releases
    )
  } else {
    NULL
  }

  pCap_range <- range(pCap_mat, na.rm = TRUE)
  scale_flow_p <- function(f) {
    pCap_range[1] + (f - flow_range[1]) /
      diff(flow_range) * diff(pCap_range)
  }

  pcap_plot <- ggplot2::ggplot(pcap_df, ggplot2::aes(x = date_label)) +
    ggplot2::geom_bar(
      ggplot2::aes(y = mean_pCap, fill = has_trial),
      stat = "identity", width = 0.75
    ) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "grey70", "FALSE" = "#e05555"),
      guide  = "none"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lo_pCap, ymax = hi_pCap),
      width = 0.25
    ) +
    {
      if (!is.null(obs_eff_df)) {
        list(
          ggplot2::geom_point(
            data = obs_eff_df,
            ggplot2::aes(x = date_label, y = obs_eff),
            shape = 21, colour = "blue", size = 2.5, inherit.aes = FALSE
          ),
          ggplot2::geom_text(
            data = obs_eff_df,
            ggplot2::aes(x = date_label,
                         y = max(pCap_mat) * 0.97,
                         label = Recaptures),
            angle = 90, hjust = 1, size = 2.2, inherit.aes = FALSE
          ),
          ggplot2::geom_text(
            data = obs_eff_df,
            ggplot2::aes(x = date_label,
                         y = max(pCap_mat) * 0.84,
                         label = Releases),
            angle = 90, hjust = 1, size = 2.2, inherit.aes = FALSE
          )
        )
      }
    } +
    ggplot2::geom_line(
      ggplot2::aes(y = scale_flow_p(flow), group = 1),
      colour = "darkgreen", linewidth = 1
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = scale_flow_p(flow)),
      colour = "darkgreen", size = 1.2
    ) +
    ggplot2::scale_y_continuous(
      name     = "Capture Probability",
      sec.axis = ggplot2::sec_axis(
        transform = ~ flow_range[1] + (. - pCap_range[1]) / diff(pCap_range) * diff(flow_range),
        name      = "Discharge (kcfs)"
      )
    ) +
    ggplot2::labs(
      x     = "First Date of Week",
      title = paste0(site, " ", run_year, ": Weekly Capture Probability")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x        = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      axis.title.y.right = ggplot2::element_text(colour = "darkgreen"),
      axis.text.y.right  = ggplot2::element_text(colour = "darkgreen")
    )

  list(abundance_plot = abundance_plot, pcap_plot = pcap_plot)
}

#' Plot year random effects and site-year capture probabilities
#'
#' @description Produces two ggplot objects from a fitted pCap model:
#'   1. All year random effects on a single panel, colour-coded by site.
#'   2. Per-site panels showing estimated pCap by year, plus a "no efficiency
#'      trial" (no_et) reference point in red.
#'
#' @param pcap A fitted Stan model object for the pCap model.
#' @param model_type Character.  Either `"all_sites"` or `"one_site"`.
#' @param site_selection Character. Required when `model_type = "one_site"`.
#'   Site name passed to `prepare_pCap_inputs()`.
#' @param use_pdev Logical. Include process error when computing predicted
#'   pCap? Default `FALSE`.
#' @param use_effort Logical. Include effort adjustment when computing
#'   predicted pCap? Default `FALSE`.
#' @param lb Numeric. Lower quantile bound. Default `0.025`.
#' @param ub Numeric. Upper quantile bound. Default `0.975`.
#'
#' @returns A named list with two ggplot objects:
#'   \describe{
#'     \item{`yr_re_plot`}{All year random effects.}
#'     \item{`site_year_plot`}{Per-site pCap by year (faceted, or single panel
#'       for `one_site`).}
#'   }
#'
#' @export
plot_year_re_with_effort <- function(pcap,
                                     model_type     = c("all_sites", "one_site"),
                                     site_selection = NULL,
                                     use_pdev       = FALSE,
                                     use_effort     = FALSE,
                                     lb             = 0.025,
                                     ub             = 0.975) {

  model_type <- match.arg(model_type)
  inv_logit  <- function(x) exp(x) / (1 + exp(x))

  # ---- Prepare inputs & posterior samples -----------------------------------
  if (model_type == "all_sites") {
    pCap_inputs <- prepare_pCap_inputs(model_type = "all_sites")
    xlabels     <- paste(
      strtrim(pCap_inputs$site_year_fit$site, 3),
      pCap_inputs$site_year_fit$run_year,
      sep = "-"
    )
    dp <- as.data.frame(pcap, pars = c("b0_pCap", "yr_re", "yr_sd_P", "pro_sd_P"))

  } else {
    if (is.null(site_selection)) {
      stop("`site_selection` must be provided when `model_type = 'one_site'`.")
    }
    pCap_inputs <- prepare_pCap_inputs(
      model_type     = "one_site",
      skew           = TRUE,
      site_selection = site_selection
    )
    xlabels <- pCap_inputs$years_fit
    dp      <- as.data.frame(pcap, pars = c("b0_pCap", "yr_re", "yr_sd_P",
                                            "pro_sd_P", "alpha"))
  }

  Nsims    <- nrow(dp)
  Ntribs   <- pCap_inputs$inputs$data$Ntribs
  Nyr_re   <- pCap_inputs$inputs$data$Nyr_re
  obs_pCap <- pCap_inputs$inputs$data$Recaptures /
    pCap_inputs$inputs$data$Releases
  ind_trib <- pCap_inputs$inputs$data$ind_trib

  # ---- 1. Year random effects plot ------------------------------------------
  yr_re_mat <- matrix(nrow = Nyr_re, ncol = 3)
  for (i in seq_len(Nyr_re)) {
    col_nm        <- paste0("yr_re[", i, "]")
    yr_re_mat[i, ] <- quantile(dp[[col_nm]], probs = c(lb, 0.5, ub))
  }

  # Colour-code by site (all_sites only)
  if (model_type == "all_sites") {
    site_prefix <- strtrim(xlabels, 3)
    uniq_prefix <- unique(site_prefix)
    palette     <- rep(c("black", "red", "goldenrod2", "purple", "hotpink",
                         "forestgreen"), length.out = length(uniq_prefix))
    pt_col <- palette[match(site_prefix, uniq_prefix)]
    site_label <- site_prefix
  } else {
    pt_col     <- rep("black", Nyr_re)
    site_label <- rep(site_selection, Nyr_re)
  }

  yr_re_df <- data.frame(
    idx       = seq_len(Nyr_re),
    label     = xlabels,
    mean_re   = yr_re_mat[, 2],
    lo_re     = yr_re_mat[, 1],
    hi_re     = yr_re_mat[, 3],
    pt_colour = pt_col,
    site      = site_label
  )

  yr_re_plot <- ggplot2::ggplot(yr_re_df,
                                ggplot2::aes(x = idx, y = mean_re)) +
    ggplot2::geom_vline(
      xintercept = seq_len(Nyr_re),
      colour = "grey80", linetype = "dotted"
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lo_re, ymax = hi_re),
      width = 0.3, colour = "grey40"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = site),
      shape = 19, size = 2
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq_len(Nyr_re),
      labels = xlabels
    ) +
    ggplot2::labs(
      x      = "",
      y      = "Random Year Effect",
      title  = "Year Random Effects",
      colour = "Site"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 7),
      legend.position = if (model_type == "all_sites") "right" else "none"
    )

  # ---- 2. Per-site pCap by year ---------------------------------------------
  # Compute pCap stats for each site-year random effect
  pCap_stats <- matrix(nrow = Nyr_re, ncol = 3)
  for (i in seq_len(Nyr_re)) {
    yr_re_col <- dp[[paste0("yr_re[", i, "]")]]

    if (model_type == "all_sites") {
      j        <- which(pCap_inputs$sites_fit ==
                          pCap_inputs$site_year_fit$site[i])
      b0       <- dp[[paste0("b0_pCap[",  j, "]")]]
      pro_sd   <- dp[[paste0("pro_sd_P[", j, "]")]]
      yr_sd    <- dp[[paste0("yr_sd_P[",  j, "]")]]
      pdev     <- rnorm(Nsims, 0, pro_sd) * use_pdev
    } else {
      b0     <- dp[["b0_pCap"]]
      pro_sd <- dp[["pro_sd_P"]]
      yr_sd  <- dp[["yr_sd_P"]]
      alpha  <- dp[["alpha"]]
      pdev   <- brms::rskew_normal(Nsims, mu = 0, sigma = pro_sd,
                                   alpha = alpha, xi = NULL, omega = NULL) *
        use_pdev
    }

    irecs  <- which(pCap_inputs$inputs$data$ind_yr == i)
    effort <- pCap_inputs$inputs$data$effort[irecs]
    effort_adj <- if (use_effort) log(mean(effort)) else 0

    p1           <- inv_logit(b0 + yr_re_col + pdev + effort_adj)
    q            <- quantile(p1, probs = c(lb, 0.5, ub))
    pCap_stats[i, ] <- c(q[[1]], mean(p1), q[[3]])
  }

  # Build per-site "no efficiency trial" (no_et) reference
  site_year_rows <- lapply(seq_len(Ntribs), function(i) {

    if (model_type == "all_sites") {
      b0     <- dp[[paste0("b0_pCap[",  i, "]")]]
      pro_sd <- dp[[paste0("pro_sd_P[", i, "]")]]
      yr_sd  <- dp[[paste0("yr_sd_P[",  i, "]")]]
      yrdev  <- rnorm(Nsims, 0, yr_sd)
      pdev   <- rnorm(Nsims, 0, pro_sd) * use_pdev
      site_nm <- pCap_inputs$sites_fit[i]
    } else {
      b0     <- dp[["b0_pCap"]]
      pro_sd <- dp[["pro_sd_P"]]
      yr_sd  <- dp[["yr_sd_P"]]
      alpha  <- dp[["alpha"]]
      yrdev  <- rnorm(Nsims, 0, yr_sd)
      pdev   <- brms::rskew_normal(Nsims, 0, pro_sd, alpha, NULL, NULL) *
        use_pdev
      site_nm <- site_selection
    }

    irecs2 <- which(pCap_inputs$inputs$data$ind_trib == i)
    effort <- pCap_inputs$inputs$data$effort[irecs2]
    effort_adj <- if (use_effort) log(mean(effort)) else 0

    forP <- inv_logit(b0 + yrdev + pdev + effort_adj)
    forP_q <- quantile(forP, probs = c(lb, 0.5, ub))

    # Which site-year indices belong to this site?
    irecs_sy <- which(pCap_inputs$inputs$data$sd_yr_ind == i)
    years_sy  <- pCap_inputs$site_year_fit$run_year[irecs_sy]
    obs_mean  <- if (length(irecs2) > 0) mean(obs_pCap[irecs2]) else NA

    # Per-year pCap rows
    per_yr <- if (length(irecs_sy) > 0) {
      data.frame(
        site     = site_nm,
        year_lbl = as.character(years_sy),
        x_ord    = seq_along(irecs_sy),
        mean_pCap = pCap_stats[irecs_sy, 2],
        lo_pCap  = pCap_stats[irecs_sy, 1],
        hi_pCap  = pCap_stats[irecs_sy, 3],
        is_no_et = FALSE,
        ref_line = NA_real_
      )
    } else {
      NULL
    }

    # no_et reference row
    n_yr <- length(irecs_sy)
    no_et <- data.frame(
      site      = site_nm,
      year_lbl  = "no_et",
      x_ord     = n_yr + 1L,
      mean_pCap = forP_q[[2]],
      lo_pCap   = forP_q[[1]],
      hi_pCap   = forP_q[[3]],
      is_no_et  = TRUE,
      ref_line  = forP_q[[2]]
    )

    rbind(per_yr, no_et)
  })

  sy_df <- do.call(rbind, site_year_rows)

  # Make x_ord a factor with year_lbl as labels, nested within site so each
  # facet gets its own independent axis — avoids the breaks/labels length mismatch.
  sy_df$x_fac <- factor(
    paste(sy_df$site, sy_df$x_ord, sep = "__"),
    levels = paste(sy_df$site, sy_df$x_ord, sep = "__")
  )

  # Horizontal reference lines (forP mean) per site
  ref_df <- sy_df[sy_df$is_no_et, c("site", "ref_line")]

  y_label <- if (use_effort) {
    "Effort-Adjusted Capture Probability"
  } else {
    "Capture Probability"
  }

  site_year_plot <- ggplot2::ggplot(sy_df,
                                    ggplot2::aes(x = x_fac, y = mean_pCap)) +
    ggplot2::geom_hline(
      data     = ref_df,
      ggplot2::aes(yintercept = ref_line),
      linetype = "dotted", colour = "red"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lo_pCap, ymax = hi_pCap,
                   colour = is_no_et),
      width = 0.25
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = is_no_et),
      shape = 19, size = 2
    ) +
    ggplot2::scale_colour_manual(
      values = c("FALSE" = "black", "TRUE" = "red"),
      guide  = "none"
    ) +
    ggplot2::scale_x_discrete(
      labels = setNames(sy_df$year_lbl, sy_df$x_fac)
    ) +
    ggplot2::facet_wrap(~ site, scales = "free_x") +
    ggplot2::labs(
      x     = "Run Year",
      y     = y_label,
      title = "Site-Year Capture Probability"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 7)
    )

  list(yr_re_plot = yr_re_plot, site_year_plot = site_year_plot)
}


#' BT SPAS X diagnostic plots
#' @details This function produces a plot with data and results of fitting the pCap and abundance models for
#' a given site and run year.
#' @param inputs a list of inputs generated by running `prepare_abundance_inputs()`
#' @param results_df a table of parameter estimates generated by running `extract_abundance_estimates()` or
#' by pulling and filtering from the database using `get_most_recent_model_results()`.
#' @returns A plot.
#' @export
#' @md
generate_diagnostic_plot_juv <- function(inputs, results_df) {

  params <- results_df |>
    filter(model_name == "bt_spas_x",
           site == inputs$site,
           year == inputs$run_year) |>
    select(-c(id, model_run_id, model_name)) |>
    pivot_wider(names_from = statistic,
                values_from = value)

  pCap_estimates <- params |>
    filter(parameter == "lt_pCap_U") |>
    select(week = week_fit,
           c("97.5", "50", "mean", "75", "sd", "25", "2.5"))

  N_estimates <- params |>
    filter(parameter == "N") |>
    select(week = week_fit,
           c("97.5", "50", "mean", "75", "sd", "25", "2.5"))

  abundance_plot <- N_estimates |>
    left_join(inputs$lp_data,
              by = "week") |>
    mutate(count_label = ifelse(is.na(count), "", count)) |>
    ggplot(aes(x = date, y = `50`)) +
    geom_bar(stat = "identity", fill = "grey", width = .75) +
    geom_errorbar(aes(x = date, ymin = `2.5`, ymax = `97.5`), width = 0.2) +
    geom_point(aes(x = date, y = lincoln_peterson_abundance),
               shape = 1, color = "blue" ,size = 3) +
    geom_point(aes(x = date, y = Inf, color = sampled),
               size = 3) +
    geom_text(aes(x = date, y = Inf,
                  label = paste(count_label),
                  angle = 90),
              hjust = 1,
              size = 3) +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "#ec5858")) +
    theme_minimal() +
    labs(x = "",
         #x = "Date",
         y = "Abundance",
         title = paste(inputs$site, inputs$run_year)) +
    #theme(axis.text.x=element_blank()) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "")

  # efficiency
  efficiency_plot <- pCap_estimates |>
    left_join(inputs$lp_data,
              by = "week") |>
    mutate(across(c(mean, `50`, `2.5`, `97.5`), plogis),
           number_released_label = ifelse(is.na(number_released), "", number_released),
           number_recaptured_label = ifelse(is.na(number_recaptured), "", number_recaptured)) |>
    ggplot(aes(x = date, y = `50`, fill = efficiency_trial)) +
    geom_bar(stat = "identity", width = .75) +
    geom_errorbar(aes(x = date, ymin = `2.5`, ymax = `97.5`), width = 0.2) +
    geom_point(aes(x = date, y = lincoln_peterson_efficiency),
               shape = 1, color = "blue", size = 3) +
    geom_text(aes(x = date, y = Inf,
                  label = paste(number_released_label, number_recaptured_label),
                  angle = 90),
              hjust = 1,
              size = 3) +
    scale_fill_manual(values = c("TRUE" = "grey", "FALSE" = "#ec5858")) +
    theme_minimal() +
    labs(x = "Date", y = "Weekly Efficiency") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "")

  # arrange together
  gridExtra::grid.arrange(abundance_plot, efficiency_plot)
}


#' BT SPAS X raw data plots
#' @details This function produces a plot with data for a site and run year.
#' @param site The site being fit
#' @param run_year The run year being fit
#' @returns A plot
#' @export
#' @md
plot_juv_data <- function(site, run_year) {
  data <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
    filter(run_year == !!run_year,
           site == !!site,
           week %in% c(seq(45, 53), seq(1, 22))) |>
    group_by(year, week, stream, site, run_year) |>
    # keep NAs in count columns
    summarise(count = if(all(is.na(count))) NA_real_ else sum(count, na.rm = TRUE),
              mean_fork_length = mean(mean_fork_length, na.rm = T),
              hours_fished = mean(hours_fished, na.rm = T),
              catch_flow_cfs = mean(flow_cfs, na.rm = T),
              average_hours_fished_during_efficiency_trials = mean(average_hours_fished_during_efficiency_trials, na.rm = T),
              standardized_flow = mean(standardized_flow, na.rm = T)
              #lgN_prior = mean(lgN_prior, na.rm = T)
    ) |>
    ungroup() |>
    left_join(SRJPEdata::weekly_juvenile_abundance_efficiency_data,
              by = c("year", "run_year", "week", "stream", "site")) |>
    mutate(count = round(count, 0),
           # change all NaNs to NAs
           # across(mean_fork_length:lgN_prior, ~ifelse(is.nan(.x), NA, .x)),
           # plot things
           lincoln_peterson_abundance = count * (number_released / number_recaptured),
           lincoln_peterson_efficiency = number_recaptured / number_released) |>
    left_join(SRJPEmodel::julian_week_to_date_lookup, by = c("week" = "Jwk")) |>
    left_join(site_order_north_south, by = "site") |>
    arrange(ns_order) |>
    select(-ns_order) |>
    mutate(date = factor(date, levels = date),
           week_index = row_number())

  abundance_plot <- data |>
    ggplot(aes(x = date, y = count)) +
    geom_bar(stat = "identity", fill = "grey", width = .75) +
    geom_text(aes(x = date, y = Inf,
                  label = paste(count),
                  angle = 90),
              hjust = 1,
              size = 3) +
    theme_minimal() +
    labs(x = "",
         y = "Abundance",
         title = paste(site, run_year)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  efficiency_plot <- data |>
    mutate(number_released_label = ifelse(is.na(number_released), "", number_released),
           number_recaptured_label = ifelse(is.na(number_recaptured), "", number_recaptured)) |>
    ggplot(aes(x = date, y = lincoln_peterson_efficiency)) +
    geom_point(shape = 1, color = "blue") +
    # geom_bar(stat = "identity", fill = "grey", width = 4) +
    geom_text(aes(x = date, y = Inf,
                  label = paste(number_released_label, number_recaptured_label),
                  angle = 90),
              hjust = 1,
              size = 3) +
    theme_minimal() +
    labs(x = "Date", y = "Weekly Efficiency") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  # arrange together
  gridExtra::grid.arrange(abundance_plot, efficiency_plot)

}


