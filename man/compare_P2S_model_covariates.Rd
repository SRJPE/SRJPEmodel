% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/passage_to_spawner_model.R
\name{compare_P2S_model_covariates}
\alias{compare_P2S_model_covariates}
\title{Compare Environmental Covariates in Passage to Spawner (P2S) Model}
\usage{
compare_P2S_model_covariates(stream)
}
\value{
A table containing the following variables:
\itemize{
\item \strong{parameter} Parameter name
\item \strong{mean} Mean of the posterior distribution for a parameter
\item \strong{se_mean} Monte Carlo standard error for summary of all chains merged (see \href{https://mc-stan.org/rstan/reference/stanfit-method-summary.html}{details})
\item \strong{sd} Standard deviation of the posterior distribution for a parameter
\item \strong{\code{2.5}} 2.5\% quantile of posterior distribution for a parameter.
\item \strong{\code{25}} 25\% quantile of posterior distribution for a parameter.
\item \strong{\code{50}} 50\% quantile of posterior distribution for a parameter.
\item \strong{\verb{7\%}} 75\% quantile of posterior distribution for a parameter.
\item \strong{\code{97.5}} 97.5\% quantile of posterior distribution for a parameter.
\item \strong{n_eff} Effective sample size for a parameter
\item \strong{Rhat} Split Rhats for a parameter
\item \strong{year} Year observed data came from
}
}
\description{
This function fits the Passage to Spawner (\code{fit_passage_to_spawner_model()}) for a given stream to all
environmental covariates. For a stream, the dataset is truncated to only include years for which all environmental
covariates are available. Parameters that are helpful for selecting a covariate are: \code{R2_data}, \code{R2_fixed}, \code{mean_redds_per_spawner},
\code{b1_survival}, \code{sigma_redds_per_spawner}, and \code{spawner_abundance_forecast}.
}
\seealso{
Other passage_to_spawner: 
\code{\link{extract_P2S_estimates}()},
\code{\link{fit_passage_to_spawner_model}()}
}
\concept{passage_to_spawner}
