#' Engagements of (spatial) Bayesian brain mapping
#'
#' Identify areas of engagement in each network from the result of (spatial) 
#'  Bayesian brain mapping.
#'
#' @param bMap Fitted (spatial) Bayesian brain map from \code{\link{BrainMap}}.
#' @param u,z Set a threshold value for engagement? A threshold value can be
#'  specified directly with \code{u}, or a z-score-like threshold in terms of
#'  standard deviations (the SD of values in the mean prior) can be specified
#'  with \code{z}. Only one type of threshold can be used. Default: \code{NULL}
#'  (do not use a threshold). Either argument can also be a vector to test
#'  multiple thresholds at once, as long as \code{type} is not \code{"!="}
#'  (to ensure the engagement regions are successive subsets).
#' @param alpha Significance level for hypothesis testing. Default: \code{0.01}.
#' @param type Type of region: \code{">"} (default), \code{"abs >"}, \code{"<"},
#'  or \code{"!="}. \code{"abs >"} tests for magnitude by taking the absolute
#'  value and then testing if they are greater than... .
#' @param method_p If the input is a \code{"bMap.[format]"} model object, the 
#'  type of multiple comparisons correction to use for p-values, or \code{NULL}
#'  for no correction. See \code{help(p.adjust)}. Default: \code{"BH"} 
#'  (Benjamini & Hochberg, i.e. the false discovery rate). Note that multiple
#'  comparisons will account for data locations, but not networks.
#' @param verbose If \code{TRUE}, display progress of algorithm. Default:
#'  \code{FALSE}.
#' @param which.nets Indices of networks for which to identify engagements. If
#'  \code{NULL} (default), use all networks.
#' @param deviation If \code{TRUE} identify significant deviations from the
#'  prior mean, rather than significant areas of engagement. Default:
#'  \code{FALSE}.
#'
#' @return A list containing engagement maps for each network, the joint and
#'  marginal PPMs for each network, and the parameters used for computing 
#'  engagement. If the input represented CIFTI- or NIFTI-format data, then the
#'  engagements maps will be formatted accordingly.
#'
#'  Use \code{summary} to obtain information about the engagements results.
#'  For CIFTI-format engagements, use \code{plot} to visualize the engagement
#'  maps.
#'
#' @export
#'
#' @importFrom fMRItools is_1 is_integer is_posNum
#' @importFrom stats pnorm p.adjust
#' @importFrom utils packageVersion
#' @importFrom grDevices colorRamp rgb
#'
#' @examples
#' \dontrun{
#'  engagements(bMap_result, alpha=.05, deviation=TRUE)
#' }
engagements <- function(
  bMap, u=NULL, z=NULL, alpha=0.01,
  type=c(">", "abs >", "<", "!="),
  method_p='BH',
  verbose=FALSE, which.nets=NULL, deviation=FALSE){

  # Setup ----------------------------------------------------------------------
  is_bMap <- inherits(bMap, "bMap.matrix") || inherits(bMap, "bMap.cifti") || inherits(bMap, "bMap.nifti")
  is_sbMap <- inherits(bMap, "sbMap.matrix") || inherits(bMap, "sbMap.cifti") || inherits(bMap, "sbMap.nifti")
  if (!(xor(is_bMap, is_sbMap))) { stop("bMap must be of class sbMap or bMap") }

  if (is_sbMap) {
    if (!requireNamespace("excursions", quietly = TRUE)) {
      stop("Package \"excursions\" needed for spatial Bayesian brain map engagements. Please install.", call. = FALSE)
    }
  }

  FORMAT <- class(bMap)[grepl("bMap", class(bMap))]
  if (length(FORMAT) != 1) { stop("Not a bMap.") }
  FORMAT <- switch(FORMAT,
    bMap.cifti = "CIFTI",
    bMap.gifti = "GIFTI",
    bMap.nifti = "NIFTI",
    bMap.matrix = "DATA",
    sbMap.cifti = "CIFTI",
    sbMap.gifti = "GIFTI",
    sbMap.nifti = "NIFTI",
    sbMap.matrix = "DATA"
  )

  if (FORMAT == "CIFTI") {
    if (!requireNamespace("ciftiTools", quietly = TRUE)) {
      stop("Package \"ciftiTools\" needed to read NIFTI data. Please install it.", call. = FALSE)
    }
  }

  # Simple argument checks
  if ((!is.null(u)) && (!is.null(z))) { stop("Set only one of `u` or `z`.") }
  stopifnot(is.null(u) || is.numeric(u))
  stopifnot(is.null(z) || is.numeric(z))
  stopifnot(fMRItools::is_posNum(alpha))
  stopifnot(alpha <= 1)
  type <- match.arg(type, c(">", "abs >", "<", "!="))
  if(alpha <= 0 | alpha >= 1) stop('alpha must be between 0 and 1')
  stopifnot(fMRItools::is_1(verbose, "logical"))
  stopifnot(fMRItools::is_1(deviation, "logical"))

  nL <- ncol(bMap$A)
  if(is.null(which.nets)) which.nets <- seq(nL)
  stopifnot(is.numeric(which.nets))
  stopifnot(which.nets == round(which.nets))
  if(min((which.nets) %in% (1:nL))==0) stop('Invalid entries in which.nets')
  nI <- length(which.nets)

  if (deviation) {
    if (!is.null(z)) { stop("`z` not compatible with `deviation==TRUE`.") }
    if (u != 0) { warning("`u != 0` not advised for `deviation==TRUE`. Proceeding anyway.") }
  }

  # Get needed metadata from `bMap`.
  Q <- bMap$omega
  mask <- bMap[["mask"]] # avoid grabbing mask_nii

  # Vectorize data.
  if (FORMAT == "CIFTI") {
    xii1 <- ciftiTools::newdata_xifti(ciftiTools::select_xifti(bMap$subjNet_mean,1), 0)
    bMap$subjNet_mean <- as.matrix(bMap$subjNet_mean)
    bMap$subjNet_se <- as.matrix(bMap$subjNet_se)
  } else if (FORMAT == "NIFTI") {
    mask_nii <- bMap$mask_nii
    bMap$subjNet_mean <- matrix(bMap$subjNet_mean[rep(mask_nii, nL)], ncol=nL)
    bMap$subjNet_se <- matrix(bMap$subjNet_se[rep(mask_nii, nL)], ncol=nL)
  }

  bMap <- bMap[c("prior_mean", "prior_var", "subjNet_mean", "subjNet_se")]
  names(bMap) <- c("t_mean", "t_var", "s_mean", "s_se")

  # Apply data mask.
  use_mask <- (!is.null(mask)) && (!all(mask))
  if (use_mask) {
    bMap$s_mean <- bMap$s_mean[mask,]
    bMap$s_se <- bMap$s_se[mask,]
  }

  nV <- nrow(bMap$s_mean)

  # Convert `z` to `u`.
  # Make `u` a nU x nL matrix (cutoffs by networks).
  u_og <- u
  if (!is.null(z)) {
    stopifnot(!is.matrix(z))
    z <- sort(z, decreasing=type=="<")
    u_mat <- outer(z, sqrt(colVars(bMap$t_mean[,which.nets,drop=FALSE])))
  } else if (!is.null(u)) {
    stopifnot(!is.matrix(u))
    u <- sort(u, decreasing=type=="<")
    u_mat <- outer(u, rep(1, nL))
  } else {
    u_mat <- matrix(0, nrow=1, ncol=nL)
  }
  nU <- nrow(u_mat)

  if (type == "!=" && nU > 1) {
    stop("Multiple u/z not compatible with '!=' test.")
  }

  if (type == "abs >") {
    bMap$s_mean <- abs(bMap$s_mean)
    type <- ">"
  }

  eng_name <- format_engagement_name(
    u=u, z=z, type=type, deviation=deviation, collapse=FALSE
  )

  # Loop over `u` to compute engagements. --------------------------------------
  out <- vector("list", nU)
  names(out) <- eng_name
  for (uu in seq(nU)) {
    if (verbose) { cat(eng_name[uu], ".\n") }
    uu_mat <- matrix(u_mat[uu,], nrow=nV, ncol=nL, byrow=TRUE)

    # Spatial Bayesian brain map engagements -----------------------------------
    if (is_sbMap) {
      if(verbose) cat('Determining areas of engagements based on joint posterior distribution of latent fields\n')

      #identify areas of engagement in each network
      engaged <- jointPPM <- marginalPPM <- vars <- matrix(NA, nrow=nV, ncol=nL)

      for(q in which.nets){
        if(verbose) cat(paste0('.. network ',q,' (',which(which.nets==q),' of ',length(which.nets),') \n'))
        inds_q <- (1:nV) + (q-1)*nV
        if(deviation){
          Dinv_mu_s <-  (as.vector(bMap$s_mean) - as.vector(bMap$t_mean) - uu_mat)/as.vector(sqrt(bMap$t_var))
        } else {
          Dinv_mu_s <- (as.vector(bMap$s_mean) - uu_mat)/as.vector(sqrt(bMap$t_var))
        }

        if(q==which.nets[1]) {
          #we scale mu by D^(-1) to use Omega for precision (actual precision of s|y is D^(-1) * Omega * D^(-1) )
          #we subtract u first since rescaling by D^(-1) would affect u too
          #save rho from first time running excursions, pass into excursions for other networks
          tmp <- system.time(res_q <- excursions::excursions(alpha = alpha, mu = Dinv_mu_s, Q = Q, type = type, u = u_mat[1], ind = inds_q)) #I think u does not matter, should check because may differ across fields
          if(verbose) print(tmp)
          rho <- res_q$rho
        } else {
          tmp <- system.time(res_q <- excursions::excursions(alpha = alpha, mu = Dinv_mu_s, Q = Q, type = type, u = u_mat[1], ind = inds_q, rho=rho))
          if(verbose) print(tmp)
        }
        engaged[,q] <- res_q$E[inds_q]
        jointPPM[,q] <- res_q$F[inds_q]
        marginalPPM[,q] <- res_q$rho[inds_q]
        vars[,q] <- res_q$vars[inds_q]
      }

      if (length(unique(u))==1) { u <- u[1] }
      if (length(unique(z))==1) { z <- z[1] }
      out[[uu]] <- list(
        engaged=engaged, jointPPM=jointPPM, marginalPPM=marginalPPM, vars=vars
      )
    }

    # Bayesian brain mapping engagements -------------------------------------------------
    if (is_bMap) {
      if(verbose) cat('\tDetermining areas of engagements based on hypothesis testing at each location\n')

      nL <- ncol(bMap$s_mean)
      if(deviation){
        t_stat <- (as.matrix(bMap$s_mean) - bMap$t_mean - uu_mat) / as.matrix(bMap$s_se)
      } else {
        t_stat <- (as.matrix(bMap$s_mean) - uu_mat) / as.matrix(bMap$s_se)
      }

      if(type=='>') pvals <- 1-pnorm(t_stat)
      if(type=='<') pvals <- pnorm(t_stat)
      if(type=='!=') pvals <- 2*(1-pnorm(abs(t_stat)))

      if(verbose) cat(paste0('\tCorrecting for multiple comparisons with method ', method_p, "\n"))

      pvals_adj <- engaged <- matrix(NA, nrow=nV, ncol=nL)
      if(is.null(method_p)) method_p <- 'none'
      for(q in which.nets){
        pvals_adj[,q] <- p.adjust(pvals[,q], method=method_p)
        engaged[,q] <- (pvals_adj[,q] < alpha)
      }

      out[[uu]] <- list(
        engaged = engaged,
        pvals = pvals, pvals_adj = pvals_adj,
        se = bMap$s_se, tstats = t_stat
      )
    }
  }

  # Format result. -------------------------------------------------------------
  engaged <- rowSums(abind(lapply(out, "[[", "engaged"), along=3), dims=2)
  engaged[] <- as.numeric(engaged)
  dimnames(engaged) <- NULL

  # Unmask data.
  if (use_mask) { engaged <- fMRItools::unmask_mat(engaged, mask) }

  # Get colors.
  if (FORMAT == "CIFTI") {
    engaged_colors <- rev(ciftiTools::make_color_pal("plasma")$color)[round(seq(nU)/nU*256)]
    # ac_rgb <- grDevices::colorRamp(c("white", "red"), space="Lab")(seq(0, nU)/(nU))
    # engaged_colors <- vector("character", nU+1)
    # for (uu in seq(nU+1)) {
    #   engaged_colors[uu] <- rgb(ac_rgb[uu,1], ac_rgb[uu,2], ac_rgb[uu,3], 255, maxColorValue=255)
    # }
  }

  # Un-vectorize data.
  if (FORMAT == "CIFTI") {
    engaged <- ciftiTools::newdata_xifti(xii1, engaged)
    engaged <- ciftiTools::move_from_mwall(engaged, -1)
    engaged <- ciftiTools::convert_xifti(
      engaged, "dlabel",
      levels_old=c(-1, 0, seq(nU)),
      levels=c(-1, 0, seq(nU)),
      labels=c("Medial Wall", "Not engaged", paste("Engaged:", eng_name)),
      colors=c("#888888", "white", engaged_colors),
      add_white=FALSE
    )
  }

  result <- c(
    list(engaged=engaged),
    out,
    list(params=list(
      alpha=alpha, method_p=method_p, type=type, u=u, z=z, deviation=deviation
    ))
  )

  if (FORMAT == "CIFTI") {
    class(result) <- "bMap_eng.cifti"
  } else if (FORMAT == "NIFTI") {
    engaged <- fMRItools::unvec_vol(engaged, mask_nii, fill=NA)
    class(result) <- "bMap_eng.nifti"
  } else {
    class(result) <- "bMap_eng.matrix"
  }

  result
}

