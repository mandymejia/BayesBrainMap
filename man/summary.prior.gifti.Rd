% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_prior.methods.R
\name{summary.prior.gifti}
\alias{summary.prior.gifti}
\alias{print.summary.prior.gifti}
\alias{print.prior.gifti}
\title{Summarize a \code{"prior.gifti"} object}
\usage{
\method{summary}{prior.gifti}(object, ...)

\method{print}{summary.prior.gifti}(x, ...)

\method{print}{prior.gifti}(x, ...)
}
\arguments{
\item{object}{Object of class \code{"prior.gifti"}.}

\item{...}{further arguments passed to or from other methods.}

\item{x}{The prior from \code{estimate_prior.gifti}}
}
\value{
A list summarizing the prior: data dimensions, options used for
prior estimation, etc.

Nothing, invisibly.

Nothing, invisibly.
}
\description{
Summary method for class \code{"prior.gifti"}
}
