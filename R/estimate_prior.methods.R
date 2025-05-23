#' Summarize a \code{"prior.cifti"} object
#'
#' Summary method for class \code{"prior.cifti"}
#'
#' @param object Object of class \code{"prior.cifti"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A list summarizing the prior: data dimensions, options used for
#'  prior estimation, etc.
#' @method summary prior.cifti
summary.prior.cifti <- function(object, ...) {
  tmean <- struct_prior(object$prior$mean, "CIFTI", object$mask_input, 
    object$params, object$dat_struct, object$template_parc_table)
  tparams <- lapply(
    object$params,
    function(q) {
      if (is.null(q)) { q <- "NULL"};
      paste0(as.character(q), collapse=" ")
    }
  )
  x <- c(
    summary(tmean),
    list(has_DR="DR" %in% names(object)),
    tparams
  )

  class(x) <- "summary.prior.cifti"
  return(x)
}

#' Summarize a \code{"prior.gifti"} object
#'
#' Summary method for class \code{"prior.gifti"}
#'
#' @param object Object of class \code{"prior.gifti"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A list summarizing the prior: data dimensions, options used for
#'  prior estimation, etc.
#' @method summary prior.gifti
summary.prior.gifti <- function(object, ...) {
  tparams <- lapply(
    object$params,
    function(q) {
      if (is.null(q)) { q <- "NULL"};
      paste0(as.character(q), collapse=" ")
    }
  )

  x <- c(
    list(
      nV=nrow(object$prior$mean),
      nL=ncol(object$prior$mean),
      hemisphere=object$dat_struct$hemisphere,
      hasDR="DR" %in% names(object)
    ),
    tparams
  )

  class(x) <- "summary.prior.gifti"
  return(x)
}

#' Summarize a \code{"prior.nifti"} object
#'
#' Summary method for class \code{"prior.nifti"}
#'
#' @param object Object of class \code{"prior.nifti"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A list summarizing the prior: data dimensions, options used for
#'  prior estimation, etc.
#' @method summary prior.nifti
summary.prior.nifti <- function(object, ...) {
  tparams <- lapply(
    object$params,
    function(q) {
      if (is.null(q)) { q <- "NULL"};
      paste0(as.character(q), collapse=" ")
    }
  )

  x <- c(
    list(
      mask_dims=dim(object$dat_struct),
      nV=nrow(object$prior$mean),
      nL=ncol(object$prior$mean),
      hasDR="DR" %in% names(object)
    ),
    tparams
  )

  class(x) <- "summary.prior.nifti"
  return(x)
}

#' Summarize a \code{"prior.matrix"} object
#'
#' Summary method for class \code{"prior.matrix"}
#'
#' @param object Object of class \code{"prior.matrix"}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A list summarizing the prior: data dimensions, options used for
#'  prior estimation, etc.
#' @method summary prior.matrix
summary.prior.matrix <- function(object, ...) {
  tparams <- lapply(
    object$params,
    function(q) {
      if (is.null(q)) { q <- "NULL"};
      paste0(as.character(q), collapse=" ")
    }
  )

  x <- c(
    list(
      nV=nrow(object$prior$mean),
      nL=ncol(object$prior$mean),
      hasDR="DR" %in% names(object)
    ),
    tparams
  )

  class(x) <- "summary.prior.matrix"
  return(x)
}

#' @rdname summary.prior.cifti
#' @export
#'
#' @param x The prior from \code{estimate_prior.cifti}
#' @param ... further arguments passed to or from other methods.
#' @return Nothing, invisibly.
#' @method print summary.prior.cifti
print.summary.prior.cifti <- function(x, ...) {
  # Get TR
  the_TR <- if (x$TR == "NULL") {
    "not provided"
  } else {
    paste("TR=", x$TR, "s.")
  }
  # Get highpass filter
  the_hpf <- if (x$hpf == "0") {
    "not used"
  } else {
    paste(as.character(x$hpf), "Hz")
  }

  cat("====PRIOR INFO=======================\n")
  cat("# Subjects:      ", x$num_subjects, "\n")
  cat("Temporal Res.:   ", the_TR, "\n")
  cat("Highpass filter: ", the_hpf, "\n")
  cat("Spatial scaling: ", x$scale, "\n")
  cat("Q2 and Q2_max:   ", paste0(x$Q2, ", ", x$Q2_max), "\n")
  cat("Pseudo retest:   ", x$pseudo_retest, "\n")
  cat("FC prior:        ", x$FC, "\n")
  cat("\n")

  class(x) <- "summary.xifti"
  print(x)
  invisible(NULL)
}

#' @rdname summary.prior.gifti
#' @export
#'
#' @param x The prior from \code{estimate_prior.gifti}
#' @param ... further arguments passed to or from other methods.
#' @return Nothing, invisibly.
#' @method print summary.prior.gifti
print.summary.prior.gifti <- function(x, ...) {
  # Get TR
  the_TR <- if (x$TR == "NULL") {
    "not provided"
  } else {
    paste("TR=", x$TR, "s.")
  }
  # Get highpass filter
  the_hpf <- if (x$hpf == "0") {
    "not used"
  } else {
    paste(as.character(x$hpf), "Hz")
  }

  cat("====PRIOR INFO=======================\n")
  cat("# Subjects:      ", x$num_subjects, "\n")
  cat("Temporal Res.:   ", the_TR, "\n")
  cat("Highpass filter: ", the_hpf, "\n")
  cat("Spatial scaling: ", x$scale, "\n")
  cat("Q2 and Q2_max:   ", paste0(x$Q2, ", ", x$Q2_max), "\n")
  cat("Pseudo retest:   ", x$pseudo_retest, "\n")
  cat("-------------------------------------\n")
  cat("# Locations:     ", x$nV, "\n")
  cat("# Networks:      ", x$nL, "\n")
  cat("Hemisphere:      ", x$hemisphere, "\n")
  cat("\n")

  invisible(NULL)
}

#' @rdname summary.prior.nifti
#' @export
#'
#' @param x The prior from \code{estimate_prior.nifti}
#' @param ... further arguments passed to or from other methods.
#' @return Nothing, invisibly.
#' @method print summary.prior.nifti
print.summary.prior.nifti <- function(x, ...) {
  # Get TR
  the_TR <- if (x$TR == "NULL") {
    "not provided"
  } else {
    paste("TR=", x$TR, "s.")
  }
  # Get highpass filter
  the_hpf <- if (x$hpf == "0") {
    "not used"
  } else {
    paste(as.character(x$hpf), "Hz")
  }

  cat("====PRIOR INFO=======================\n")
  cat("# Subjects:      ", x$num_subjects, "\n")
  cat("Temporal Res.:   ", the_TR, "\n")
  cat("Highpass filter: ", the_hpf, "\n")
  cat("Spatial scaling: ", x$scale, "\n")
  cat("Q2 and Q2_max:   ", paste0(x$Q2, ", ", x$Q2_max), "\n")
  cat("Pseudo retest:   ", x$pseudo_retest, "\n")
  cat("-------------------------------------\n")
  cat("Mask dims:       ", paste0(x$mask_dims, collapse=" x "), "\n")
  cat("Vectorized dims:\n")
  cat("# Locations:     ", x$nV, "\n")
  cat("# Networks:      ", x$nL, "\n")
  cat("\n")

  invisible(NULL)
}

#' @rdname summary.prior.matrix
#' @export
#'
#' @param x The prior from \code{estimate_prior.cifti}
#' @param ... further arguments passed to or from other methods.
#' @return Nothing, invisibly.
#' @method print summary.prior.matrix
print.summary.prior.matrix <- function(x, ...) {
  # Get TR
  the_TR <- if (x$TR == "NULL") {
    "not provided"
  } else {
    paste("TR=", x$TR, "s.")
  }
  # Get highpass filter
  the_hpf <- if (x$hpf == "0") {
    "not used"
  } else {
    paste(as.character(x$hpf), "Hz")
  }

  cat("====PRIOR INFO=======================\n")
  cat("# Subjects:      ", x$num_subjects, "\n")
  cat("Temporal Res.:   ", the_TR, "\n")
  cat("Highpass filter: ", the_hpf, "\n")
  cat("Spatial scaling: ", x$scale, "\n")
  cat("Q2 and Q2_max:   ", paste0(x$Q2, ", ", x$Q2_max), "\n")
  cat("Pseudo retest:   ", x$pseudo_retest, "\n")
  cat("-------------------------------------\n")
  cat("Dimensions:      \n")
  cat("# Locations:     ", x$nV, "\n")
  cat("# Networks:      ", x$nL, "\n")
  cat("\n")

  invisible(NULL)
}

#' @rdname summary.prior.cifti
#' @export
#'
#' @return Nothing, invisibly.
#' @method print prior.cifti
print.prior.cifti <- function(x, ...) {
  print.summary.prior.cifti(summary(x))
}

#' @rdname summary.prior.gifti
#' @export
#'
#' @return Nothing, invisibly.
#' @method print prior.gifti
print.prior.gifti <- function(x, ...) {
  print.summary.prior.gifti(summary(x))
}

#' @rdname summary.prior.nifti
#' @export
#'
#' @return Nothing, invisibly.
#' @method print prior.nifti
print.prior.nifti <- function(x, ...) {
  print.summary.prior.nifti(summary(x))
}

#' @rdname summary.prior.matrix
#' @export
#'
#' @return Nothing, invisibly.
#' @method print prior.matrix
print.prior.matrix <- function(x, ...) {
  print.summary.prior.matrix(summary(x))
}

#' Plot prior
#'
#' @param x The prior from \code{estimate_prior.cifti}
#' @param stat \code{"mean"}, \code{"sd"}, or \code{"both"} (default). By
#'  default the square root of the variance prior is shown; another option is
#'  \code{stat="var"} to instead display the variance prior directly.
#' @param var_method \code{"non-negative"} (default) or \code{"unbiased"}
#' @param ... Additional arguments to \code{view_xifti}
#' @return The plot
#' @export
#' @method plot prior.cifti
plot.prior.cifti <- function(x, stat=c("both", "mean", "sd", "var"),
  var_method=c("non-negative", "unbiased"), ...) {
  stopifnot(inherits(x, "prior.cifti"))

  if (!requireNamespace("ciftiTools", quietly = TRUE)) {
    stop("Package \"ciftiTools\" needed to read NIFTI data. Please install it.", call. = FALSE)
  }

  var_method <- match.arg(var_method, c("non-negative", "unbiased"))

  # Check `...`
  args <- list(...)
  has_title <- "title" %in% names(args)
  has_idx <- "idx" %in% names(args)
  has_fname <- "fname" %in% names(args)

  # Check `stat`
  stat <- tolower(stat)
  if (has_idx && length(args$idx)>1 && !("fname" %in% names(args))) {
    if (identical(stat, c("both", "mean", "sd", "var"))) {
      stat <- "mean"
    } else {
      stat <- match.arg(stat, c("both", "mean", "sd", "var"))
    }
    if (stat == "both") {
      if (!("fname" %in% names(args))) {
        warning(
          "For multiple `idx`, use one call to plot() ",
          "for the mean prior, ",
          "and a separate call for the variance prior. ",
          "Showing the mean prior now."
        )
        stat <- "mean"
      }
    }
  }
  stat <- match.arg(stat, c("both", "mean", "sd", "var"))

  # Print message saying what's happening.
  msg1 <- ifelse(has_idx,
    "Plotting the",
    "Plotting the first component's"
  )
  msg2 <- switch(stat,
    both="mean and sqrt(variance) prior.",
    mean="mean prior.",
    sd="sqrt(variance) prior.",
    var="variance prior."
  )
  cat(msg1, msg2, "\n")

  # Plot
  out <- list(mean=NULL, var=NULL)
  if (stat == "both") { stat <- c("mean", "sd") }
  for (ss in stat) {
    ssname <- if (ss == "mean") {
      ss
    } else if (var_method=="non-negative") {
      "varNN"
    } else {
      "varUB"
    }
    if (ss=="var" && var_method=="unbiased") { x$prior[[ssname]][] <- pmax(0, x$prior[[ssname]]) }
    if (ss=="sd") {
      x$prior[[ssname]] <- sqrt(x$prior[[ssname]])
    }
    tss <- struct_prior(x$prior[[ssname]], "CIFTI", x$mask_input, x$params, x$dat_struct, x$template_parc_table)
    if (ss=="sd") {
      ssname <- paste0("sqrt ", ssname)
    }

    args_ss <- args
    # Handle title and idx
    ### No title: use the component names if available, and the indices if not.
    if (!has_title && !has_idx) {
      args_ss$title <- if (!is.null(tss$meta$cifti$names)) {
        tss$meta$cifti$names[1]
      } else {
        "First component"
      }
    } else if (!has_title) {
      args_ss$title <- if (!is.null(tss$meta$cifti$names)) {
        tss$meta$cifti$names[args$idx]
      } else {
        paste("Component", args$idx)
      }
    }
    ### Append the statistic name.
    args_ss$title <- paste0(args_ss$title, " (", ssname, ")")
    # Handle fname
    if (has_fname) {
      fext <- if (grepl("html$", args_ss$fname[1])) {
        "html"
      } else if (grepl("pdf$", args_ss$fname[1])) {
        "pdf"
      } else {
        "png"
      }
      args_ss$fname <- gsub(paste0(".", fext), "", args_ss$fname, fixed=TRUE)
      args_ss$fname <- paste0(args_ss$fname, "_", ss, ".", fext)
    }
    out[[ss]] <- do.call(
      ciftiTools::view_xifti, c(list(tss), args_ss)
    )
  }

  invisible(out)
}

#' Plot prior
#'
#' @param x The prior from \code{estimate_prior.gifti}
#' @param stat \code{"mean"}, \code{"sd"}, or \code{"both"} (default). By
#'  default the square root of the variance prior is shown; another option is
#'  \code{stat="var"} to instead display the variance prior directly.
#' @param var_method \code{"non-negative"} (default) or \code{"unbiased"}
#' @param ... Additional arguments to \code{view_xifti}
#' @return The plot
#' @export
#' @method plot prior.gifti
plot.prior.gifti <- function(x, stat=c("both", "mean", "sd", "var"),
  var_method=c("non-negative", "unbiased"), ...) {
  stopifnot(inherits(x, "prior.gifti"))

  if (!requireNamespace("ciftiTools", quietly = TRUE)) {
    stop("Package \"ciftiTools\" needed to read NIFTI data. Please install it.", call. = FALSE)
  }

  if (x$dat_struct$hemisphere == "left")  {
    y <- ciftiTools::as_cifti(cortexL=x$prior$mean[,1,drop=FALSE] * 0)
  } else {
    y <- ciftiTools::as_cifti(cortexR=x$prior$mean[,1,drop=FALSE] * 0)
  }
  y <- ciftiTools::move_from_mwall(y)
  x$dat_struct <- y; class(x) <- "prior.cifti"
  plot.prior.cifti(x, stat, var_method, ...)
}

#' Plot prior
#'
#' Based on \code{oro.nifti::image}.
#'
#' Consider using \code{struct_prior} to obtain the 3D volumes to plot with a different
#'  viewer function (e.g. from \code{oro.nifti}) if desired.
#'
#' @param x The prior from \code{estimate_prior.nifti}
#' @param stat \code{"mean"} (default), \code{"sd"}, or \code{"var"}.
#'  (\code{"sd"} will show the square root of the variance prior.)
#' @param var_method \code{"non-negative"} (default) or \code{"unbiased"}
#' @param plane,n_slices,slices Anatomical plane and which slice indices to
#'  show.
#'  Default: 9 axial slices.
#' @param ... Additional arguments to \code{oro.nifti::image}
#' @return The plot
#' @export
#' @method plot prior.nifti
plot.prior.nifti <- function(x, stat=c("mean", "sd", "var"),
  plane=c("axial", "sagittal", "coronal"), n_slices=9, slices=NULL,
  var_method=c("non-negative", "unbiased"), ...) {
  stopifnot(inherits(x, "prior.nifti"))

  if (!requireNamespace("oro.nifti", quietly = TRUE)) {
    stop("Package \"oro.nifti\" needed to read NIFTI data. Please install it.", call. = FALSE)
  }

  var_method <- match.arg(var_method, c("non-negative", "unbiased"))

  # Check `...`
  args <- list(...)
  has_title <- "title" %in% names(args)
  has_idx <- "idx" %in% names(args)
  has_fname <- "fname" %in% names(args)

  # Check `idx`
  if (has_idx) {
    stopifnot(length(args$idx)==1)
    stopifnot(is.numeric(args$idx) && args$idx==round(args$idx))
    stopifnot(args$idx %in% seq(ncol(x$prior$mean)))
  } else {
    args$idx <- 1
  }
  idx <- args$idx; args$idx <- NULL

  # Check `stat`
  stat <- tolower(stat)
  if (has_idx && length(args$idx)>1 && !("fname" %in% names(args)) && identical(stat, c("mean", "sd", "var"))) {
    stat <- "mean"
  }
  stat <- match.arg(stat, c("mean", "sd", "var"))

  # Print message saying what's happening.
  msg1 <- ifelse(has_idx,
    "Plotting the",
    "Plotting the first component's"
  )
  msg2 <- switch(stat,
    mean="mean prior.",
    sd="sqrt(variance) prior.",
    var="variance prior."
  )
  cat(msg1, msg2, "\n")

  # Plot
  out <- list(mean=NULL, var=NULL)

  plane <- match.arg(plane, c("axial", "sagittal", "coronal"))
  args$plane <- plane
  plane_dim <- switch(plane, axial=3, coronal=2, sagittal=1)
  if (is.null(slices)) {
    if (is.null(n_slices)) { warning("Using 9 slices."); n_slices <- 9 }
    n_slices <- as.numeric(n_slices)
    if (length(n_slices) > 1) { warning("Using the first entry of `slice`."); n_slices <- n_slices[1] }
    # Pick slices that are spaced out, and with many voxels.
    mask_count <- apply(x$dat_struct, plane_dim, sum)
    ns_all <- length(mask_count)
    slices <- seq(ns_all)
    # Remove slices with zero voxels.
    slices <- slices[mask_count != 0]
    mask_count <- mask_count[mask_count != 0]
    ns_all <- length(mask_count)
    if (n_slices > length(slices)) {
      warning(
        "`n_slices` is larger than the number of non-empty slices (",
        length(slices), "). Showing all non-empty slices."
      )
      n_slices <- length(slices)
    }
    # Remove slices with few voxels.
    if (n_slices < (ns_all / 2)) {
      slices <- slices[mask_count > quantile(mask_count, .33)]
    }
    slices <- slices[round(seq(1, length(slices), length.out=n_slices))]
  } else {
    slices <- as.numeric(slices)
    stopifnot(all(slices %in% seq(dim(x$dat_struct)[plane_dim])))
  }

  ssname <- if (stat == "mean") {
    stat
  } else if (var_method=="non-negative") {
    "varNN"
  } else {
    "varUB"
  }
  if (stat=="var" && var_method=="unbiased") { x$prior[[ssname]][] <- pmax(0, x$prior[[ssname]]) }
  tss <- struct_prior(x$prior[[ssname]], "NIFTI", x$mask_input, x$params, x$dat_struct, x$template_parc_table)
  tss <- tss[,,,idx]

  if (plane=="axial") {
    tss <- tss[,,slices,drop=FALSE]
  } else if (plane=="coronal") {
    tss <- tss[,slices,,drop=FALSE]
  } else if (plane=="sagittal") {
    tss <- tss[slices,,,drop=FALSE]
  } else { stop() }

  if (stat=="sd") {
    tss <- sqrt(tss)
    ssname <- paste0("sqrt ", ssname)
  }

  args_ss <- args
  args_ss$plane <- plane
  # Handle title and idx
  if (!has_title && !has_idx) {
    c1name <- "First component"
  }
  if (has_title) { stop("Not supported yet.") }
  if (has_fname) { stop("Not supported yet. Call `pdf` or `png` beforehand, and then `dev.off`.") }
  do.call(
    oro.nifti::image,
    c(list(oro.nifti::as.nifti(tss)), args_ss)
  )
}

#' Plot prior
#'
#' @param x The prior from \code{estimate_prior.matrix}
#' @param ... Additional arguments
#' @return The plot
#' @export
#' @method plot prior.matrix
plot.prior.matrix <- function(x, ...) {
  stop("Not supported yet.")
}
