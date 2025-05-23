% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_prior.R
\name{estimate_prior_from_DR}
\alias{estimate_prior_from_DR}
\title{Estimate prior from DR}
\usage{
estimate_prior_from_DR(DR, LV = NULL)
}
\arguments{
\item{DR}{the test/retest dual regression estimates, as an array with
dimensions \eqn{M \times N \times (L \times V)}, where \eqn{M} is the number
of visits (2), \eqn{N} is the number of subjects, \eqn{L} is the number of
brain networks, and \eqn{V} is the number of data locations.

(\eqn{L} and \eqn{V} are collapsed because they are treated equivalently
in the context of calculating the variance decomposition and priors).}

\item{LV}{A length-two integer vector giving the dimensions \eqn{L} and
\eqn{V} to reshape the result. Default: \code{NULL} (do not reshape the
result).}
}
\value{
List of two elements: the priors and the variance decomposition.

There are two version of the variance prior: \code{varUB} gives the
unbiased variance estimate, and \code{varNN} gives the upwardly-biased
non-negative variance estimate. Values in \code{varUB} will need to be
clamped above zero before using in \code{BrainMap}.
}
\description{
Estimate variance decomposition and priors from DR estimates.
}
\keyword{internal}
