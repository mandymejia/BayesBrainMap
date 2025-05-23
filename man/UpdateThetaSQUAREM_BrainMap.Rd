% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EM_BrainMap.R
\name{UpdateThetaSQUAREM_BrainMap}
\alias{UpdateThetaSQUAREM_BrainMap}
\title{Update theta SQUAREM}
\usage{
UpdateThetaSQUAREM_BrainMap(
  theta_vec,
  prior_mean,
  prior_var,
  meshes,
  BOLD,
  C_diag,
  H,
  Hinv,
  s0_vec,
  D,
  Dinv_s0,
  update_nu0sq,
  verbose
)
}
\arguments{
\item{theta_vec}{Vector of initial parameter values}

\item{prior_mean}{Passed to UpdateTheta_BrainMap function}

\item{prior_var}{Passed to UpdateTheta_BrainMap function}

\item{meshes}{Passed to UpdateTheta_BrainMap function}

\item{BOLD}{Passed to UpdateTheta_BrainMap function}

\item{C_diag}{Passed to UpdateTheta_BrainMap function}

\item{s0_vec}{Passed to UpdateTheta_BrainMap function}

\item{D}{Passed to UpdateTheta_BrainMap function}

\item{Dinv_s0}{Passed to UpdateTheta_BrainMap function}

\item{update_nu0sq}{For non-spatial model: updating \code{nu0sq} is recommended
if dimension reduction was not performed, and is not recommended if it was.}

\item{verbose}{Passed to UpdateTheta_BrainMap function}
}
\value{
Vector of updated parameter values
}
\description{
Helper function for SQUAREM for estimating parameters
}
\keyword{internal}
