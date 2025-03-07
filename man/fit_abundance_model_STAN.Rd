% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/juvenile_abundance.R
\name{fit_abundance_model_STAN}
\alias{fit_abundance_model_STAN}
\title{Fit abundance model in STAN}
\usage{
fit_abundance_model_STAN(input, pCap_fit, lt_pCap_Us)
}
\arguments{
\item{input}{A list containing the inputs for the abundance model. This is created by calling
\code{prepare_abundance_inputs()}.}

\item{pCap_fit}{A STANfit object resulting from running \code{fit_pCap_model()}.}

\item{lt_pCap_Us}{A named list produced by running \code{generate_lt_pCap_Us()}.}
}
\value{
a STANfit object containing abundance model fits.
}
\description{
Fit abundance model in STAN
}
\details{
This function calls a STAN model for abundance. It is currently in development.
}
