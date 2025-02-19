% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/juvenile_abundance.R
\name{prepare_abundance_inputs}
\alias{prepare_abundance_inputs}
\title{Prepare inputs for abundance model}
\usage{
prepare_abundance_inputs(
  site,
  run_year,
  effort_adjust = c(T, F),
  input_catch_data = NULL,
  input_efficiency_data = NULL
)
}
\arguments{
\item{site}{site for which you want to fit the model}

\item{run_year}{run year for which you want to fit the model}

\item{effort_adjust}{whether or not you want to use catch adjusted by effort.}

\item{input_catch_data}{Optional argument for weekly catch data.
Defaults to \code{SRJPEdata::weekly_juvenile_abundance_catch_data}. If
passed in, structure of data frame must match that of the default.}

\item{input_efficiency_data}{Optional argument for weekly efficiency data.
Defaults to \verb{SRJPEdata::SRJPEdata::weekly_juvenile_abundance_efficiency_data}. If
passed in, structure of data frame must match that of the default.}
}
\value{
a list:
\itemize{
\item \strong{abundance_inputs} a list of data and inits for input into the abundance model.
\item \strong{model_name} The version of the pCap logic to use for generating lt_pCap_U values. Either
\code{all_mark_recap}, \code{missing_mark_recap}, \code{no_mark_recap}, or \code{no_mark_recap_no_trib}.
\item \strong{weeks_fit} The weeks fit for the abundance model
\item \strong{weeks_date} Associated dates for the weeks fit for the abundance model.
}
}
\description{
Prepare inputs for abundance model
}
\details{
This function prepares data for input into an abundance model.
}
