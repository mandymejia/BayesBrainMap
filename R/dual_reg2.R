#' Dual Regression wrapper
#'
#' Wrapper to \code{dual_reg} used by `estimate_prior`. The format of `BOLD`
#'  (and `BOLD2`) must be provided, and `template` must be vectorized if applicable.
#'
#' @param BOLD,BOLD2 Subject-level fMRI data in one of the following formats:
#'  a CIFTI file path, a \code{"xifti"} object, a NIFTI file path, a \code{"nifti"} object, or
#'  \eqn{V \times T} numeric matrices, where \eqn{V} is the number of data locations and
#'  \eqn{T} is the number of timepoints.
#'
#'  If \code{BOLD2} is provided it must be in the same format as \code{BOLD};
#'  \code{BOLD} will be the test data and \code{BOLD2} will be the retest data.
#'
#'  If \code{BOLD2} is not provided, \code{BOLD} will be split in half;
#'  the first half will be the test data and the second half will be the retest data.
#' @param template The group ICA map or parcellation as a (vectorized) numeric
#'  matrix (\eqn{V \times Q}). If it's an ICA map, its columns will be centered.
#' @param template_parc_table If the template is a parcellation, provide the
#'  parcellation table here. Default: \code{NULL}.
#' @param keepA Keep the resulting \strong{A} matrices, or only return the
#'  \strong{S} matrices (default)?
#' @param drop_first (Optional) Number of volumes to drop from the start of each
#'  BOLD session. Default: \code{0}.
#' @param nuisance (Optional) Nuisance matrix to regress from the BOLD data.
#'  If \code{BOLD2} is provided, should be a length-2 list with the first entry
#'  corresponding to \code{BOLD} and the second to \code{BOLD2}. If \code{NULL},
#'  do not remove any nuisance signals.
#'
#'  Nuisance regression is performed in a simultaneous regression with any spike
#'  regressors from \code{scrub} and DCT bases from \code{hpf}.
#'
#'  Note that the nuisance matrices should be provided with timepoints matching
#'  the original \code{BOLD} and \code{BOLD2} irregardless of \code{drop_first}.
#'  Nuisance matrices will be truncated automatically if \code{drop_first>0}.
#' @param scrub (Optional) Numeric vector of integers giving the indices
#'  of volumes to scrub from the BOLD data. (List the volumes to remove, not the
#'  ones to keep.) If \code{BOLD2} is provided, should be a length-two list with
#'  the first entry corresponding to \code{BOLD} and the second to \code{BOLD2}.
#'
#'  Scrubbing is performed within a nuisance regression by adding a spike
#'  regressor to the nuisance design matrix for each volume to scrub.
#'
#'  Note that indices are counted beginning with the first index in the
#'  \code{BOLD} session irregardless of \code{drop_first}. The indices will be
#'  adjusted automatically if \code{drop_first>0}.
#' @inheritParams TR_param
#' @inheritParams hpf_param
#' @inheritParams lpf_param
#' @inheritParams GSR_Param
#' @inheritParams scale_by_Param
#' @inheritParams scale_sm_FWHM_Param
#' @param scale_sm_surfL,scale_sm_surfR Required only for "local" smoothing
#'  (see \code{scale_sm_FWHM}). To smooth the scale estimates, provide the surface
#'  geometries along which to smooth, as GIFTI geometry files or 
#'  \code{ciftiTools} \code{"surf"} objects. The resolutions should match with
#'  those of the \code{BOLD} data.
#'
#'  To create a \code{"surf"} object from data, see
#'  \code{\link[ciftiTools]{make_surf}}.
#' 
#'  If not provided, the fs_LR "midthickness" surfaces will be used.
#' @param brainstructures Only applies if the entries of \code{BOLD} are CIFTI file paths.
#'  Character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{c("all")}.
#' @param resamp_res Only applies if the entries of \code{BOLD} are CIFTI file paths.
#'  Resample the data upon reading it in? Default: \code{NULL} (no resampling).
#' @param mask Required if and only if the entries of \code{BOLD} are NIFTI file paths or
#'  \code{"nifti"} objects. This is a brain map formatted as a binary array of the same
#'  size as the fMRI data, with \code{TRUE} corresponding to in-mask voxels.
#' @param format Expected format of \code{BOLD} and \code{BOLD2}. Should be one
#'  of the following: a \code{"CIFTI"} file path, a \code{"xifti"} object, a
#'  \code{"NIFTI"} file path, a \code{"nifti"} object, or a \code{"data"} matrix.
#' @param Q2,Q2_max Obtain dual regression estimates after denoising? Denoising is
#'  based on modeling and removing nuisance ICs. It may result in a cleaner
#'  estimate for smaller datasets, but it may be unnecessary (and time-consuming)
#'  for larger datasets.
#'
#'  Set \code{Q2} to control denoising: use a positive integer to specify the
#'  number of nuisance ICs, \code{NULL} to have the number of nuisance ICs
#'  estimated by PESEL, or zero (default) to skip denoising.
#'
#'  If \code{is.null(Q2)}, use \code{Q2_max} to specify the maximum number of
#'  nuisance ICs that should be estimated by PESEL. \code{Q2_max} must be less
#'  than \eqn{T * .75 - Q} where \eqn{T} is the minimum number of timepoints in
#'  each fMRI scan and \eqn{Q} is the number of networks in \code{template}. If \code{NULL}
#'  (default), \code{Q2_max} will be set to \eqn{T * .50 - Q}, rounded.
#' @param FC_updateA_path Where to save the BOLD if \code{FC_updateA}.
#' @inheritParams varTol_Param
#' @param maskTol Tolerance for number of locations masked out due to low
#'  variance or missing values. If more than this many locations are masked out,
#'  this subject is skipped without calculating dual regression. \code{maskTol}
#'  can be specified either as a proportion of the number of locations (between
#'  zero and one), or as a number of locations (integers greater than one).
#'  Default: \code{.1}, i.e. up to 10\% of locations can be masked out.
#'
#'  If \code{BOLD2} is provided, masks are calculated for each scan and then
#'  the intersection of the masks is used.
#' @param verbose Display progress updates? Default: \code{TRUE}.
#'
#' @return The dual regression \strong{S} matrices, or both the \strong{S}
#'  and \strong{A} matrices if \code{keepA}, or \code{NULL} if dual
#'  regression was skipped due to too many masked data locations.
#'
#' @importFrom fMRItools dual_reg norm_BOLD mask_BOLD
#'
#' @keywords internal
dual_reg2 <- function(
  BOLD, BOLD2=NULL,
  format=c("CIFTI", "xifti", "GIFTI", "gifti", "NIFTI", "nifti", "RDS", "data"),
  template, template_parc_table=NULL,
  mask=NULL,
  keepA=FALSE,
  drop_first=0, nuisance=NULL, scrub=NULL,
  TR=NULL, hpf=NULL, lpf=NULL,
  GSR=FALSE,
  scale_by=c("mean", "sd", "none"),
  scale_sm_FWHM=4,
  scale_sm_surfL=NULL,
  scale_sm_surfR=NULL,
  Q2=0, Q2_max=NULL,
  # NA_limit=.1,
  brainstructures="all", resamp_res=NULL,
  FC_updateA_path=NULL,
  varTol=1e-6, maskTol=.1,
  verbose=TRUE){

  if (verbose) { extime <- Sys.time() }

  keepA <- as.logical(keepA); stopifnot(length(keepA)==1)

  scale_by <- match.arg(scale_by, c("mean", "sd", "none"))
  if (scale_by == "none") { scale_sm_FWHM <- 0 } # avoid unnecessary work
  stopifnot(fMRItools::is_1(scale_sm_FWHM, "numeric"))
  scale_sm <- switch(
    as.character(scale_sm_FWHM), 
    "0"="none", "Inf"="global", "local"
  )
  if (scale_sm=="local") { stopifnot(scale_sm_FWHM > 0) }

  if (identical(TR, "from_xifti_metadata")) {
    if ((!is.null(hpf) && hpf != 0) || (!is.null(lpf) && is.finite(lpf))) {
      stop("`hpf` or `lpf` was requested, but `TR` was neither provided nor able to be inferred from the data. Please provide `TR`.")
    }
    TR <- NULL  # no temporal filtering requested, so TR isn't needed
  }

  if (!is.null(BOLD2)) {
    if (!is.null(nuisance)) { stopifnot(is.list(nuisance) && length(nuisance)==2) }
    if (!is.null(scrub)) { stopifnot(is.list(scrub) && length(scrub)==2) }
  }

  do_denoise <- !((!is.null(Q2) && Q2==0) || (!is.null(Q2_max) && Q2_max==0))

  # No other arg checks: check them before calling this function.

  # For `"xifti"` data for handling the medial wall and smoothing.
  xii1 <- NULL

  # Prepare output.
  out <- list(test = NULL, retest = NULL)

  # Load helper variables.
  retest <- !is.null(BOLD2)
  format <- match.arg(format, c("CIFTI", "xifti", "GIFTI", "gifti", "NIFTI", "nifti", "RDS", "data"))
  FORMAT <- get_FORMAT(format)
  check_req_ifti_pkg(FORMAT)

  template_parc <- !is.null(template_parc_table)
  nQ <- if (template_parc) { nrow(template_parc_table) } else { ncol(template) } # unused

  if (FORMAT=="NIFTI") { stopifnot(!is.null(mask)) } 
  if (is.null(mask)) {
    nI <- nV <- nrow(template)
  } else if (FORMAT=="NIFTI") {
    nI <- dim(drop(mask))
    nV <- sum(mask)
  } else {
    nI <- length(mask); nV <- sum(mask)
  }

  # Get `BOLD` (and `BOLD2`) as a data matrix or array.  -----------------------
  if (verbose) { cat("\tReading in data... ") }
  if (FORMAT == "CIFTI") {
    if (is.character(BOLD)) { BOLD <- ciftiTools::read_cifti(BOLD, brainstructures=brainstructures, resamp_res=resamp_res) }
    if (ciftiTools::is.xifti(BOLD)) {
      if (scale_sm == "local") {
        xii1 <- ciftiTools::convert_xifti(ciftiTools::select_xifti(BOLD, 1), "dscalar") * 0 # Extract surface and mwall from input xifti for scale smoothing
      }
      BOLD <- as.matrix(BOLD)
    }
    stopifnot(is.matrix(BOLD))
    if (retest) {
      if (is.character(BOLD2)) { BOLD2 <- ciftiTools::read_cifti(BOLD2, brainstructures=brainstructures, resamp_res=resamp_res) }
      if (ciftiTools::is.xifti(BOLD2)) { BOLD2 <- as.matrix(BOLD2) }
      stopifnot(is.matrix(BOLD2))
    }
  } else if (FORMAT == "GIFTI") {
    if (is.character(BOLD)) { BOLD <- gifti::readgii(BOLD) }
    stopifnot(gifti::is.gifti(BOLD))
    ghemi <- BOLD$file_meta["AnatomicalStructurePrimary"]
    if (!(ghemi %in% c("CortexLeft", "CortexRight"))) {
      stop("AnatomicalStructurePrimary metadata missing or invalid for template.")
    }
    ghemi <- switch(ghemi, CortexLeft="left", CortexRight="right")
    if (scale_sm == "local") {
      if (ghemi == "left") {
        xii1 <- ciftiTools::select_xifti(ciftiTools::as.xifti(cortexL=do.call(cbind, BOLD$data)), 1) * 0
      } else if (ghemi == "right") {
        xii1 <- ciftiTools::select_xifti(ciftiTools::as.xifti(cortexR=do.call(cbind, BOLD$data)), 1) * 0
      } else { stop() }
      xii1$meta$cifti$intent <- 3006
    }
    BOLD <- do.call(cbind, BOLD$data)

    stopifnot(is.matrix(BOLD))
    if (retest) {
      if (is.character(BOLD2)) { BOLD2 <- gifti::readgii(BOLD2) }
      if (inherits(BOLD2, "gifti")) { BOLD2 <- do.call(cbind, BOLD2$data) }
      stopifnot(is.matrix(BOLD2))
    }
    nI <- nV <- nrow(template)
  } else if (FORMAT == "NIFTI") {
    if (is.character(BOLD)) { BOLD <- RNifti::readNifti(BOLD) }
    stopifnot(length(dim(BOLD)) > 1)
    if (retest) {
      if (is.character(BOLD2)) { BOLD2 <- RNifti::readNifti(BOLD2) }
      stopifnot(length(dim(BOLD2)) > 1)
    }
  } else if (FORMAT == "MATRIX") {
    if (is.character(BOLD)) { BOLD <- readRDS(BOLD) }
    stopifnot(is.matrix(BOLD))
    if (retest) {
      if (is.character(BOLD2)) { BOLD2 <- readRDS(BOLD2) }
      stopifnot(is.matrix(BOLD2))
    }
    nI <- nV <- nrow(template)
  } else { stop() }

  dBOLD <- dim(BOLD)
  ldB <- length(dim(BOLD))

  # If `retest`, ensure that spatial dimensions of `BOLD2` match with `BOLD`.
  if (retest) {
    stopifnot(length(dim(BOLD)) == length(dim(BOLD2)))
    stopifnot(all(dBOLD[seq(ldB-1)] == dim(BOLD2)[seq(ldB-1)]))
  }

  # Check BOLD (and BOLD2) dimensions correspond with `template` and `mask`.
  if(!(ldB-1 == length(nI))) { stop("`template` and BOLD spatial dimensions do not match.") }
  if(!all(dBOLD[seq(ldB-1)] == nI)) { stop("`template` and BOLD spatial dimensions do not match.") }

  # Vectorize `BOLD` (and `BOLD2`). --------------------------------------------
  if (FORMAT=="NIFTI") {
    BOLD <- matrix(BOLD, nrow=prod(nI))[as.logical(mask),,drop=FALSE]
    stopifnot(nrow(BOLD) == nV)
    if (retest) {
      BOLD2 <- matrix(BOLD2, nrow=prod(nI))[as.logical(mask),,drop=FALSE]
      stopifnot(nrow(BOLD2) == nV)
    }
  } else if (!is.null(mask)) {
    # Mask out the locations.
    BOLD <- BOLD[mask,,drop=FALSE]
    if (!is.null(xii1)) {
      xiitmp <- as.matrix(xii1)
      xiitmp[!mask,] <- NA
      xii1 <- ciftiTools::move_to_mwall(ciftiTools::newdata_xifti(xii1, xiitmp))
    }
    nV <- nrow(BOLD)
    if (retest) {
      BOLD2 <- BOLD2[mask,,drop=FALSE]
      stopifnot(nrow(BOLD2)==nV)
    }
  }

  # `drop_first` ---------------------------------------------------------------
  # Do here, before NA values check.
  stopifnot(fMRItools::is_posNum(drop_first, zero_ok=TRUE))
  if (drop_first > 0) {
    stopifnot(drop_first < ncol(BOLD) - 2)
    # Drop columns from BOLD; drop rows from nuisance; adjust `scrub`.
    # (Already done for scrubbing.)
    if (!retest) {
      BOLD <- BOLD[,-seq(drop_first),drop=FALSE]
      if (!is.null(scrub)) { scrub <- scrub[scrub > drop_first] - drop_first }
      if (!is.null(nuisance)) { nuisance <- nuisance[-seq(drop_first),,drop=FALSE] }
    } else {
      stopifnot(drop_first < ncol(BOLD2) - 2)
      BOLD <- BOLD[,-seq(drop_first),drop=FALSE]
      BOLD2 <- BOLD2[,-seq(drop_first),drop=FALSE]
      if (!is.null(scrub[[1]])) { scrub[1] <- list(scrub[[1]][scrub[[1]] > drop_first] - drop_first) }
      if (!is.null(scrub[[2]])) { scrub[2] <- list(scrub[[2]][scrub[[2]] > drop_first] - drop_first) }
      if (!is.null(nuisance[[1]])) { nuisance[1] <- list(nuisance[[1]][-seq(drop_first),,drop=FALSE]) }
      if (!is.null(nuisance[[2]])) { nuisance[2] <- list(nuisance[[2]][-seq(drop_first),,drop=FALSE]) }
    }

    # Do not do this again.
    drop_first <- 0
  }

  # Check for missing values. --------------------------------------------------
  mask2 <- fMRItools::mask_BOLD(BOLD, varTol=varTol)
  if (retest) { mask2 <- mask2 & fMRItools::mask_BOLD(BOLD2, varTol=varTol) }
  use_mask2 <- !all(mask2)
  if (use_mask2) {
    # Coerce `maskTol` to number of locations.
    stopifnot(is.numeric(maskTol) && length(maskTol)==1 && maskTol >= 0)
    if (maskTol < 1) { maskTol <- maskTol * nV }
    # Skip this scan if `maskTol` is surpassed.
    if (sum(!mask2) > maskTol) { return(NULL) }
    # Mask out the locations.
    BOLD <- BOLD[mask2,,drop=FALSE]
    template <- template[mask2,,drop=FALSE]
    if (retest) { BOLD2 <- BOLD2[mask2,,drop=FALSE] }
    if (!is.null(xii1)) {
      xiitmp <- as.matrix(xii1)
      xiitmp[!mask2,] <- NA
      xii1 <- ciftiTools::move_to_mwall(ciftiTools::newdata_xifti(xii1, xiitmp))
    }
    nV <- nrow(BOLD)

    # [TO DO]: replace with fMRIscrub::unmask_mat(..., mask_dim=2)
    # For later
    unmask <- function(S, mask) {
      S2 <- matrix(NA, nrow=nrow(S), ncol=length(mask))
      S2[,mask] <- S
      S2
    }
    unmask_vec <- function(vec, mask) {
      vec2 <- rep(NA, length(mask))
      vec2[mask] <- vec
      vec2
    }
  }

  # Prep for dual regression ---------------------------------------------------
  if (is.null(xii1) && scale_sm=="local") { 
    message("No surface data: skipping smoothing of scale estimates.")
    scale_sm_FWHM <- 0
    scale_sm <- "none" 
  }

  # Add surfaces to `xii1`
  if (!is.null(xii1) && scale_sm=="local") {
    xii1 <- ciftiTools::add_surf(xii1, surfL=scale_sm_surfL, surfR=scale_sm_surfR)
  }

  ### Define helper functions ---

  # Do the big regression. Do not center and do not scale.
  big_nreg_BOLD <- function(B) { norm_BOLD (
    BOLD=B,
    nuisance=nuisance, scrub=scrub,
    TR=TR, hpf=hpf, lpf=lpf,
    center_rows=FALSE, center_cols=FALSE,
    scale_by="none", scale_sm_FWHM=0
  ) }

  # Center and scale. Do not do the big regression again.
  center_scale_BOLD <- function(B) { norm_BOLD(
    BOLD=B,
    TR=TR, hpf=NULL, lpf=NULL,
    scale_by=scale_by, scale_sm_FWHM=scale_sm_FWHM, scale_sm_xifti=xii1,
    center_rows=TRUE, center_cols=GSR
  ) }

  # Handle continuous vs. discrete prior
  DR_FUN <- if (template_parc) {
    function(template, ...) { fMRItools::dual_reg_parc(parc=template, ...) }
  } else {
    function(template, parc_vals, ...) { fMRItools::dual_reg(GICA=template, ...) }
  }

  # Dual regression without any norm_BOLD stuff
  DR_noNorm <- function(B) { DR_FUN(
    B, template=template, parc_vals=template_parc_table$Key,
    # Disable norm stuff that's enabled by default
    hpf=0, scale_by="none"#, GSR=FALSE
  ) }

  # Get the first dual regression results. -------------------------------------
  if (verbose) { cat("\n\tDual regression... ") }

  if (!retest) {
    # 1. Normalizing. ---
    # Do the big nuisance regression (first half of `norm_BOLD`)
    #   (everything but centering and scaling).
    BOLD <- big_nreg_BOLD(BOLD)
    # Get `nT` (was updated by scrubbing and `drop_first`)
    nT <- ncol(BOLD)
    # Split BOLD in half.
    part1 <- seq(round(nT/2))
    part2 <- setdiff(seq(nT), part1)
    # Center and scale. (No nuisance regression, temporal filtering, etc.)
    BOLDh1 <- center_scale_BOLD(BOLD[, part1, drop=FALSE]) #first half of data
    BOLDh2 <- center_scale_BOLD(BOLD[, part2, drop=FALSE]) #second half of data
    
    # 2. Two DR's. ---
    out$test <- DR_noNorm(BOLDh1)
    out$retest <- DR_noNorm(BOLDh2)

  } else {
    # 1. Normalizing. ---
    BOLD <- norm_BOLD(
      BOLD, 
      nuisance=nuisance[[1]], scrub=scrub[[1]],
      TR=TR, hpf=hpf, lpf=lpf,
      center_rows=TRUE, center_cols=GSR,
      scale_by=scale_by, scale_sm_FWHM=scale_sm_FWHM, 
      scale_sm_xifti=xii1
    )
    BOLD2 <- norm_BOLD(
      BOLD2, 
      nuisance=nuisance[[2]], scrub=scrub[[2]],
      TR=TR, hpf=hpf, lpf=lpf,
      center_rows=TRUE, center_cols=GSR,
      scale_by=scale_by, scale_sm_FWHM=scale_sm_FWHM, 
      scale_sm_xifti=xii1
    )

    # 2. Two DR's. ---
    out$test <- DR_noNorm(BOLD)
    out$retest <- DR_noNorm(BOLD2)
  }

  BOLDss <- list(
    test = if (!retest) { BOLDh1 } else { BOLD },
    retest = if (!retest) { BOLDh2 } else { BOLD2 }
  )

  if (!retest && !do_denoise) rm(BOLD)

  # Get `sigma_sq` -------------------------------------------------------------
  # part inside colSums() is TxV
  calc_sigma_sq <- function(DR, B) colSums((DR$A %*% DR$S - t(B))^2) / ncol(B)
  for (sess in c("test","retest")) {
    out[[sess]]$sigma_sq <- calc_sigma_sq(out[[sess]], BOLDss[[sess]])
  }
  rm(BOLDss)

  # Return these DR results if denoising is not needed. ------------------------
  if (!do_denoise) {

    if (!is.null(FC_updateA_path)) {
      BOLDkeep <- list(
        test = if (!retest) { BOLDh1 } else { BOLD },
        retest = if (!retest) { BOLDh2 } else { BOLD2 }
      )
      if (use_mask2) { BOLDkeep$mask2 <- mask2 }
      saveRDS(BOLDkeep, file.path(FC_updateA_path, "BOLDkeep.rds"))
    }

    if (retest) { rm(BOLD, BOLD2) } else { rm(BOLDh1, BOLDh2) }

    for (sess in c("test", "retest")) {
      if (use_mask2) { out[[sess]]$sigma_sq <- unmask_vec(out[[sess]]$sigma_sq, mask2) }
      if (!keepA) { out[[sess]]$A <- out[[sess]]$A2 <- NULL }
      if (use_mask2) { out[[sess]]$S <- unmask(out[[sess]]$S, mask2) }
    }

    if (verbose) { cat(" Done!\n") }
    if (verbose) { print(Sys.time() - extime) }
    return(out)
  }

  if (!retest) { rm(BOLDh1, BOLDh2) }

  # Estimate and deal with nuisance ICs. ---------------------------------------
  if (verbose) { cat(" Denoising... ") }
  # If !retest, we prefer to estimate nuisance ICs across the full scan
  # and then halve it after.
  if (!retest) {
    # Note: up to this line, `BOLD` has gone thru the big regression
    #   but has not been centered or scaled.
    BOLD <- center_scale_BOLD(BOLD)
    # Get initial estimate of networks for full scan.
    BOLD_DR <- DR_noNorm(BOLD)
    BOLD <- rm_nuisIC(BOLD, DR=BOLD_DR[c("A", "S")], Q2=Q2, Q2_max=Q2_max, verbose=verbose)
    rm(BOLD_DR)
    BOLD2 <- BOLD[, part2, drop=FALSE]
    BOLD <- BOLD[, part1, drop=FALSE]
  } else {
    BOLD <- rm_nuisIC(BOLD, DR=out$test[c("A", "S")], Q2=Q2, Q2_max=Q2_max, verbose=verbose)
    BOLD2 <- rm_nuisIC(BOLD2, DR=out$retest[c("A", "S")], Q2=Q2, Q2_max=Q2_max, verbose=verbose)
  }
  
  # Center `BOLD` and `BOLD2` (again).
  BOLD <- BOLD - rowMeans(BOLD)
  BOLD2 <- BOLD2 - rowMeans(BOLD2)

  if (!is.null(FC_updateA_path)) {
    BOLDkeep <- list(
      test = BOLD,
      retest = BOLD2
    )
    saveRDS(BOLDkeep, file.path(FC_updateA_path, "BOLDkeep.rds"))
  }

  # Do DR again. ---------------------------------------------------------------
  if (verbose) { cat("\n\tDual regression again... ") }
  
  out$test_preclean <- out$test
  out$test <- DR_noNorm(BOLD)
  out$retest_preclean <- out$retest
  out$retest <- DR_noNorm(BOLD2)

  out$test$sigma_sq <- calc_sigma_sq(out$test, BOLD)
  out$retest$sigma_sq <- calc_sigma_sq(out$retest, BOLD2)
  rm(BOLD, BOLD2); gc()

  for (sess in c("test", "retest", "test_preclean", "retest_preclean")) {
    if (use_mask2) { out[[sess]]$sigma_sq <- unmask_vec(out[[sess]]$sigma_sq, mask2) }
    if (!keepA) { out[[sess]]$A <- out[[sess]]$A2 <- NULL }
    if (use_mask2) { out[[sess]]$S <- unmask(out[[sess]]$S, mask2) }
  }

  if (verbose) { cat(" Done!\n") }
  if (verbose) { print(Sys.time() - extime) }
  out
}
