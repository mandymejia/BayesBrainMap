% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BrainMap.methods.R
\name{plot.bMap.nifti}
\alias{plot.bMap.nifti}
\title{Plot prior}
\usage{
\method{plot}{bMap.nifti}(
  x,
  stat = c("mean", "se"),
  plane = c("axial", "sagittal", "coronal"),
  n_slices = 9,
  slices = NULL,
  ...
)
}
\arguments{
\item{x}{The result of \code{BrainMap} with NIFTI data}

\item{stat}{\code{"mean"} (default), \code{"se"}}

\item{plane, n_slices, slices}{Anatomical plane and which slice indices to show.
Default: 9 axial slices.}

\item{...}{Additional arguments}
}
\value{
The plot
}
\description{
Plot prior
}
