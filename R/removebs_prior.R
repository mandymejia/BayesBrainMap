#' Remove brain structure from CIFTI prior
#'
#' Remove a brain structure from a CIFTI prior
#'
#' @param x The \code{"prior.cifti"} object.
#' @param remove \code{"cortex_left"}, \code{"cortex_right"}, and/or \code{"subcortical"}.
#'
#' @keywords internal
removebs_prior <- function(x, remove=NULL){
  stopifnot(inherits(x, "prior.cifti"))
  remove <- match.arg(remove, c("cortex_left", "cortex_right", "subcortical"), several.ok=TRUE)

  # Remove brain structure(s) from data.
  x$prior[1:3] <- lapply(x$prior[1:3], function(y){
    as.matrix(ciftiTools::remove_xifti(ciftiTools::newdata_xifti(x$dat_struct, y), remove=remove))
  })

  x$var_decomp <- lapply(x$var_decomp, function(y){
    as.matrix(ciftiTools::remove_xifti(ciftiTools::newdata_xifti(x$dat_struct, y), remove=remove))
  })
  if (!is.null(x$sigma_sq0)) {
    x$sigma_sq0 <- as.matrix(ciftiTools::remove_xifti(
      ciftiTools::newdata_xifti(x$dat_struct, x$sigma_sq0),
      remove=remove
    ))
  }

  # Get new `dat_struct` and mask.
  x$dat_struct <- ciftiTools::remove_xifti(x$dat_struct, remove=remove)
  x$mask <- !is.na(x$prior$mean[,1])

  x
}
