#' TR
#'
#' @param TR The temporal resolution of the data, i.e. the time between volumes,
#'  in seconds. \code{TR} is required for detrending with \code{hpf}.
#'
#' @name TR_param
#' @keywords internal
NULL

#' hpf
#'
#' @param hpf The frequency at which to apply a highpass filter to the data
#'  during pre-processing, in Hertz. Default: \code{NULL} (disabled). If the
#'  data has not already been highpass filtered, a recommended filter value is
#'  \code{.01} Hz.

#'  The highpass filter serves to detrend the data, since low-frequency
#'  variance is associated with noise. Highpass filtering is accomplished by
#'  nuisance regression of discrete cosine transform (DCT) bases.
#'
#'  Note the \code{TR} argument is required for highpass filtering. If
#'  \code{TR} is not provided, an error will be raised if \code{hpf} is set.
#'
#' @name hpf_param
#' @keywords internal
NULL

#' lpf
#'
#' @param lpf The frequency at which to apply a lowpass filter to the data
#'  during pre-processing, in Hertz. Default: \code{NULL} (disabled). Though
#'  this argument is offered, we do not recommend using a a LPF. 

#'  Note the \code{TR} argument is required for lowpass filtering. If
#'  \code{TR} is not provided, an error will be raised if \code{lpf} is set.
#'
#' @name lpf_param
#' @keywords internal
NULL

#' hpf_lpf
#' 
#' @param hpf,lpf The frequencies at which to apply temporal filtering to the
#'  data during pre-processing, in Hertz. Set either to \code{NULL} to disable.
#'  Default: \code{0.01} Hz highpass filter, and \code{NULL} for the lowpass
#'  filter (disabled). Highpass filtering is accomplished by nuisance regression
#'  of discrete cosine transform (DCT) bases. Lowpass filtering is accomplished
#'  by a non-linear filter: \code{\link[fMRItools]{fsl_bptf}}.
#' 
#'  Note the \code{TR} argument is required for temporal filtering. If
#'  \code{TR} is not provided, an error will be raised if \code{hpf} or
#'  \code{lpf} is set.
#' 
#' @name hpf_lpf_param
#' @keywords internal
NULL

#' varTol
#'
#' @param varTol Tolerance for variance of each data location. For each scan,
#'  locations which do not meet this threshold are masked out of the analysis.
#'  Default: \code{1e-6}. Variance is calculated on the original data, before
#'  any normalization. Set to \code{0} to avoid removing locations due to
#'  low variance.
#'
#' @name varTol_Param
#' @keywords internal
NULL

#' scale_by
#'
#' @param scale_by Scale the BOLD at each voxel based on either its 
#'  \code{"mean"} (default), or its \code{"sd"}. Mean scaling cannot be used if
#'  the \code{BOLD} have already been de-meaned. 
#' 
#' @name scale_by_Param
#' @keywords internal
NULL

#' scale_sm_FWHM
#'
#' @param scale_sm_FWHM Full width at half maximum (FWHM) for smoothing the
#'  estimates of scale across brain locations (see \code{scale_by}), to reduce
#'  the variance of the estimates. Set to \code{0} to disable smoothing, or
#'  \code{Inf} for "global" smoothing (estimate and use one measure of scale 
#'  across the entire brain). Otherwise, for "local" smoothing, this should be a
#'  positive number. Note that local smoothing is only available for surface 
#'  data input (CIFTI or GIFTI \code{BOLD}). Default: local smoothing with a 
#'  FWHM of \code{4}.
#' 
#' @name scale_sm_FWHM_Param
#' @keywords internal
NULL

#' GSR
#'
#' @param GSR Center BOLD across columns (each image)? This
#'  is equivalent to performing global signal regression. Default:
#'  \code{FALSE}.
#'
#' @name GSR_Param
#' @keywords internal
NULL
