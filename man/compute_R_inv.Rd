% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EM_BrainMap.R
\name{compute_R_inv}
\alias{compute_R_inv}
\title{Compute SPDE and prior precision matrices for S}
\usage{
compute_R_inv(meshes, kappa, C1 = 1/(4 * pi), rm_extra = FALSE)
}
\arguments{
\item{meshes}{\code{NULL} for spatial independence model, otherwise a list of
objects of class "BrainMap_mesh" containing the triangular mesh (see
\code{\link{make_mesh}}) for each brain structure.}

\item{kappa}{Current estimates of SPDE parameter kappa for each latent field}

\item{C1}{Constant, equal to \eqn{1/(4*pi)} for a 2-dimensional field with alpha=2}

\item{rm_extra}{If \code{TRUE}. remove extra (non-data) vertices from the mesh for greater computational efficiency}
}
\value{
A list containing R inverse and SPDE matrices
}
\description{
Compute SPDE matrices (F, G and GFinvG) and prior precision matrix for S
}
\keyword{internal}
