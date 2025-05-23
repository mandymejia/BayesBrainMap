% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_prior.methods.R
\name{plot.prior.gifti}
\alias{plot.prior.gifti}
\title{Plot prior}
\usage{
\method{plot}{prior.gifti}(
  x,
  stat = c("both", "mean", "sd", "var"),
  var_method = c("non-negative", "unbiased"),
  ...
)
}
\arguments{
\item{x}{The prior from \code{estimate_prior.gifti}}

\item{stat}{\code{"mean"}, \code{"sd"}, or \code{"both"} (default). By
default the square root of the variance prior is shown; another option is
\code{stat="var"} to instead display the variance prior directly.}

\item{var_method}{\code{"non-negative"} (default) or \code{"unbiased"}}

\item{...}{Additional arguments to \code{view_xifti}}
}
\value{
The plot
}
\description{
Plot prior
}
