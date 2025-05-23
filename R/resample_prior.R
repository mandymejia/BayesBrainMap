#' Resample CIFTI prior
#'
#' Resample a CIFTI prior to a new spatial resolution.
#'
#' @param x The \code{"prior.cifti"} object.
#' @param resamp_res The new resampling resolution.
#' @param verbose Give occasional updates? Default: \code{FALSE}.
#'
#' @return The resampled \code{"prior.cifti"} object.
#' @keywords internal
#'
resample_prior <- function(x, resamp_res, verbose=FALSE){
  stopifnot(inherits(x, "prior.cifti"))

  if (!requireNamespace("ciftiTools", quietly = TRUE)) {
    stop("Package \"ciftiTools\" needed to work with CIFTI data. Please install it.", call. = FALSE)
  }

  tm_struct_mask <- !(names(x$prior) %in% c("FC", "FC_Chol")) # mean, varUB, varNN

  # Resample the data.
  if (verbose) { cat("Resampling priors.\n") }
  x$prior[tm_struct_mask] <- lapply(x$prior[tm_struct_mask], function(y){
    as.matrix(ciftiTools::resample_xifti(ciftiTools::newdata_xifti(x$dat_struct, y), resamp_res=resamp_res, verbose=FALSE))
  })

  if (verbose) { cat("Resampling variance decomposition.\n") }
  x$var_decomp <- lapply(x$var_decomp, function(y){
    as.matrix(ciftiTools::resample_xifti(ciftiTools::newdata_xifti(x$dat_struct, y), resamp_res=resamp_res, verbose=FALSE))
  })
  if (!is.null(x$sigma_sq0)) {
    x$sigma_sq0 <- c(as.matrix(ciftiTools::resample_xifti(ciftiTools::newdata_xifti(x$dat_struct, x$sigma_sq), resamp_res=resamp_res, verbose=FALSE)))
  }

  if (verbose) { cat("Formatting new prior object.\n") }
  # Replace `NaN` values with NA values.
  x$prior[tm_struct_mask] <- lapply(x$prior[tm_struct_mask], function(y){y[] <- ifelse(is.nan(y), NA, y)})
  x$var_decomp <- lapply(x$var_decomp, function(y){y[] <- ifelse(is.nan(y), NA, y)})
  if (!is.null(x$sigma_sq0)) {
    x$sigma_sq0 <- ifelse(is.nan(x$sigma_sq0), NA, x$sigma_sq0)
  }

  # Get new `dat_struct` and mask.
  x$dat_struct <- ciftiTools::resample_xifti(x$dat_struct, resamp_res=resamp_res)
  x$mask <- !is.na(x$prior$mean[,1])

  x
}
