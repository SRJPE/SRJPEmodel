#' @title Passage to Spawner Parameter Estimates
#' @name P2S_model_fits
#' @description The passage to spawner model, documented in \code{vignette("passage_to_spawner_submodel", package = "SRJPEmodel")},
#' estimates the relationship between upstream passage counts and spawner counts (either `redd` or `holding`) for
#' several streams. When fit, the model produces an object `P2S_model_fits` that contains the parameter names,
#' mean parameter estimate, confidence intervals, and diagnostics. Currently, this model is only fit for
#' `battle creek`, `clear creek`, `deer creek`, and `mill creek`.
#' @format A tibble with 343 rows and 12 columns
#' \itemize{
#'   \item \code{par_names}: Parameter name
#'   \item \code{mean}: Mean of the posterior distribution for a parameter
#'   \item \code{se_mean}: Monte Carlo standard error for summary of all chains merged (see [details](https://mc-stan.org/rstan/reference/stanfit-method-summary.html))
#'   \item \code{sd}: Standard deviation of the posterior distribution for a parameter
#'   \item \code{`2.5%`}: 2.5% quantile of posterior distribution for a parameter
#'   \item \code{`25%`}: 25% quantile of posterior distribution for a parameter
#'   \item \code{`50%`}: 50% quantile of posterior distribution for a parameter
#'   \item \code{`75%`}: 75% quantile of posterior distribution for a parameter
#'   \item \code{`97.5%`}: 97.5% quantile of posterior distribution for a parameter
#'   \item \code{n_eff}: Effective sample size for a parameter
#'   \item \code{Rhat}: Split Rhats for a parameter
#'   \item \code{stream}: Stream name for the parameter estimate
#'   }
'P2S_model_fits'


#' @title Passage to Spawner Covariate Comparison results
#' @name P2S_comparison_results
#' @description Output from running `SRJPEmodel::compare_P2S_model_covariates()`
#' @format A tibble with 140 rows and 9 columns
#' \itemize{
#'   \item \code{par_names}: Parameter name
#'   \item \code{stream}: Stream for which the model was run
#'   \item \code{mean}: Mean of the posterior distribution for a parameter
#'   \item \code{median}: Median the posterior distribution for a parameter
#'   \item \code{sd`}: Standard deviation of the posterior distribution for a parameter
#'   \item \code{lcl}: 2.5% lower confidence limit
#'   \item \code{ucl}: 97.5 upper confidence limit
#'   \item \code{covar_considered}: The covariate associated with the results
#'   \item \code{convergence_metric}: Whether or not the convergence metric was reached
#'   }
'P2S_model_fits'
