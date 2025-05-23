#' Apply data structure to priors
#'
#' @param prior The prior
#' @param FORMAT "CIFTI", "GIFTI", "NIFTI", or "DATA"
#' @param dat_struct The data structure
#' @param mask_input The input mask
#' @param params The params
#'
#' @keywords internal
struct_prior <- function(prior, FORMAT, mask_input, params, dat_struct, template_parc_table){
  # Un-apply the input mask.
  if (!is.null(mask_input)) {
    if (FORMAT=="NIFTI") {
      prior <- fMRItools::unvec_vol(prior, drop(mask_input))
    } else {
      prior <- fMRItools::unmask_mat(prior, mask_input)
    }
  }

  # Add metadata.
  if (FORMAT == "CIFTI") {
    if (!requireNamespace("ciftiTools", quietly = TRUE)) {
      stop("Package \"ciftiTools\" needed to work with CIFTI data. Please install it.", call. = FALSE)
    }
    prior <- ciftiTools::newdata_xifti(dat_struct, prior)
    prior$meta$cifti$names <- if (!is.null(template_parc_table)) {
      rownames(template_parc_table)
    } else {
      paste("network", params$inds)
    }

  } else if (FORMAT == "GIFTI") {
    if (!requireNamespace("ciftiTools", quietly = TRUE)) {
      stop("Package \"ciftiTools\" needed to work with GIFTI data. Please install it.", call. = FALSE)
    }

    prior <- ciftiTools:::as.metric_gifti(
      prior, hemisphere=dat_struct$hemisphere
    )

  } else if (FORMAT == "NIFTI") {
    prior <- RNifti::asNifti(prior, reference=mask_input)
  }

  prior
}

#' Export prior
#'
#' Export the priors (mean, variance, and FC) as separate files for
#'  visualization or processing outside of \code{BayesBrainMap}.
#'
#' @param x The result of \code{estimate_prior}
#' @param out_fname Use \code{NULL} (default) to just return the prior
#'  objects directly. Otherwise, use a character vector of length 3 or 4 of file
#'  path(s) to save the output to:
#'  the mean prior, the variance prior, the variance decomposition, and
#'  the FC prior if present, in that order. If one file name is provided,
#'  it will be appended with
#'  \code{"_mean.[file_ext]"} for the prior mean map,
#'  \code{"_var.[file_ext]"} for the prior variance map,
#'  \code{"_varDecomp.rds"} for the variance decomposition, and
#'  \code{"_FC.rds"} where \code{[file_ext]}
#'  will be \code{"dscalar.nii"} for CIFTI input, \code{"nii"} for NIFTI input,
#'  and \code{"rds"} for data input.
#' @param var_method \code{"non-negative"} (default) or \code{"unbiased"}
#'
#' @return If \code{is.null(out_fname)}, the priors in data matrix,
#'  \code{"xifti"}, or \code{"nifti"} format, to match the format of the
#'  original BOLD data. Otherwise, the paths to the new files specified by
#'  \code{out_fname}. If prior includes functional connectivity components,
#'  the FC prior and its mean and variance will be included.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'  tm <- estimate_prior(cii1_fnames, cii2_fnames, gICA_fname)
#'  export_prior(tm, out_fname="my_prior", var_method="unbiased")
#' }
#'
export_prior <- function(x, out_fname=NULL, var_method=c("non-negative", "unbiased")){

  # Check prior format.
  FORMAT <- class(x)[grepl("prior", class(x))]
  if (length(FORMAT) != 1) { stop("Not a prior.") }
  FORMAT <- switch(FORMAT,
    prior.cifti = "CIFTI",
    prior.gifti = "GIFTI",
    prior.nifti = "NIFTI",
    prior.matrix = "MATRIX"
  )
  FORMAT_extn <- switch(FORMAT, CIFTI=".dscalar.nii", GIFTI=".func.gii", NIFTI=".nii", MATRIX=".rds")

  var_method <- match.arg(var_method, c("non-negative", "unbiased"))
  var_name <- switch(var_method, `non-negative`="varNN", unbiased="varUB")

  x$prior$varUB[] <- pmax(0, x$prior$varUB)

  FC <- "FC" %in% names(x$prior)

  # Fix for old version
  if (FORMAT == "CIFTI" && !is.null(x$dat_struct$meta$subcort$labels)) {
    if (!requireNamespace("ciftiTools", quietly = TRUE)) {
      stop("Package \"ciftiTools\" needed to export prior. Please install it.", call. = FALSE)
    }
    sub_levs <- levels(x$dat_struct$meta$subcort$labels)
    if (length(sub_levs) != length(ciftiTools::substructure_table()$ciftiTools_Name)) {
      x$dat_struct$meta$subcort$labels <- factor(
        x$dat_struct$meta$subcort$labels,
        levels = ciftiTools::substructure_table()$ciftiTools_Name
      )
      stopifnot(ciftiTools:::is.subcort_labs(x$dat_struct$meta$subcort$labels))
    }
  }

  # `out_fname` ----------------------------------------------------------------
  if (!is.null(out_fname)) {
    out_fname <- as.character(out_fname)
    if (!all(dir.exists(dirname(out_fname)))) { stop('Directory part of `out_fname` does not exist.') }
    if (length(out_fname) == 1) {
      if (!endsWith(out_fname, FORMAT_extn)) { out_fname <- paste0(out_fname, FORMAT_extn) }
      out_fname <- c(
        gsub(paste0(FORMAT_extn, "$"), paste0("_mean", FORMAT_extn), out_fname),
        gsub(paste0(FORMAT_extn, "$"), paste0("_var", FORMAT_extn), out_fname),
        gsub(paste0(FORMAT_extn, "$"), paste0("_varDecomp.rds"), out_fname),
        gsub(paste0(FORMAT_extn, "$"), paste0("_FC.rds"), out_fname)
      )
      if (!FC) { out_fname <- out_fname[seq(3)] }
    } else if (length(out_fname) == 3 + as.numeric(FC)) {
      if (!all(endsWith(out_fname[seq(2)], FORMAT_extn))) {
        out_fname[seq(2)] <- paste0(out_fname[seq(2)], FORMAT_extn)
      }
      if (!endsWith(out_fname[3], ".rds")) {
        out_fname[3] <- paste0(out_fname[3], ".rds")
      }
      if (FC && !endsWith(out_fname[4], ".rds")) {
        out_fname[4] <- paste0(out_fname[4], ".rds")
      }
    } else {
      stop(
        "`out_fname` should be a length 1 or 3/4 character vector giving the ",
        "names for:\n\tThe mean prior,\n\tthe variance prior,",
        "\n\tthe variance decomposition, and \n\tthe FC prior.\n"
      )
    }
  }

  tm_struct_mask <- !(names(x$prior) %in% c("FC", "FC_Chol"))
  x$prior[tm_struct_mask] <- lapply(
    x$prior[tm_struct_mask], struct_prior,
    FORMAT, x$mask_input, x$params, x$dat_struct, x$template_parc_table
  )

  # Select the chosen variance decomposition.
  x$prior <- list(
    mean = x$prior$mean,
    var = x$prior[[var_name]],
    FC = x$prior$FC
    #, FC_Chol=FC_Chol # [TO DO]: want to export anything in `FC_Chol`?
  )

  #compute mean and variance of FC
  if(FC){
    Q <- nrow(x$prior$FC$psi)
    FC_mean <- x$prior$FC$psi/(x$prior$FC$nu - Q - 1)
    FC_var <- FC_mean*0
    for(q1 in 1:Q){
      for(q2 in 1:Q){
        FC_var[q1,q2] <- IW_var(x$prior$FC$nu, Q, FC_mean[q1,q2], FC_mean[q1,q1], FC_mean[q2,q2])
      }
    }
    x$prior$FC$mean <- FC_mean
    x$prior$FC$var <- FC_var
  }

  # Add params to `"xifti"` metadata; resample it.
  if (FORMAT == "CIFTI") {
    x$params <- lapply(
      x$params,
      function(q) {
        if (is.null(q)) { q <- "NULL"};
        paste0(as.character(q), collapse=" ")
      }
    )
    x$prior$mean$meta$cifti$misc <- c(list(prior="mean"), x$params)
    x$prior$var$meta$cifti$misc <- c(list(prior="var"), x$params)
  }

  # Save
  if (!is.null(out_fname)) {
    if (FORMAT == "CIFTI") {
      ciftiTools::write_cifti(x$prior$mean, out_fname[1])
      ciftiTools::write_cifti(x$prior$var, out_fname[2])
    } else if (FORMAT == "GIFTI") {
      if (!requireNamespace("gifti", quietly = TRUE)) {
        stop("Package \"gifti\" needed to write NIFTI data. Please install it.", call. = FALSE)
      }
      gifti::writegii(x$prior$mean, out_fname[1])
      gifti::writegii(x$prior$var, out_fname[2])
    } else if (FORMAT == "NIFTI") {
      if (!requireNamespace("RNifti", quietly = TRUE)) {
        stop("Package \"RNifti\" needed to write NIFTI data. Please install it.", call. = FALSE)
      }
      RNifti::writeNifti(x$prior$mean, out_fname[1])
      RNifti::writeNifti(x$prior$var, out_fname[2])
    } else {
      saveRDS(x$prior$mean, out_fname[1])
      saveRDS(x$prior$var, out_fname[2])
    }
    saveRDS(x$var_decomp, out_fname[3])
    if (FC) { saveRDS(x$prior$FC, out_fname[4]) }
  }

  if (is.null(out_fname)) {
    return(x$prior)
  } else {
    return(invisible(out_fname))
  }
}



