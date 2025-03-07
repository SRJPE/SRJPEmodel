% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survival_model.R
\name{fit_survival_model}
\alias{fit_survival_model}
\title{Fit Survival Model}
\usage{
fit_survival_model(survival_inputs)
}
\arguments{
\item{survival_inputs}{The named list object generated by running \code{prepare_survival_inputs()}.}
}
\value{
a STANfit object with the survival model fit.
}
\description{
Fit Survival Model
}
\details{
This model calls a survival STAN model that estimates survival in the Sacramento and Feather/Butte systems.
}
