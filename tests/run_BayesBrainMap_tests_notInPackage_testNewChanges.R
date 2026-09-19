# Build --> Install and Restart

# Setup ------------------------------------------------------------------------
# ciftiTools
library(ciftiTools)
print(packageVersion("ciftiTools"))
ciftiTools.setOption("wb_path", "~/Applications")
library(fMRItools)

# BayesBrainMap
library(BayesBrainMap)
# roxygen2::roxygenize("../../BayesBrainMap")
print(packageVersion("BayesBrainMap"))

library(RNifti)
library(gifti)
library(rgl)

# file paths
data_dir <- "data_notInPackage"
subjects <- c(100307, 100408, 100610)
cii_fnames <- c(
  paste0(data_dir, "/", subjects, "_rfMRI_REST1_LR_Atlas.dtseries.nii"),
  paste0(data_dir, "/", subjects, "_rfMRI_REST2_LR_Atlas.dtseries.nii")
)
giiL_fnames <- gsub("dtseries.nii", "sep.L.func.gii", cii_fnames, fixed=TRUE)
giiL_ROI_fnames <- gsub("dtseries.nii", "sep.ROI_L.func.gii", cii_fnames, fixed=TRUE)
nii_fnames <- gsub("_Atlas.dtseries.nii", ".nii.gz", cii_fnames, fixed=TRUE)
rds_fnames <-gsub("dtseries.nii", "rds", cii_fnames, fixed=TRUE)

template_fname <- c(
  cii = file.path(data_dir, "melodic_IC_100.4k.dscalar.nii"),
  gii = file.path(data_dir, "melodic_IC_100.4k.sep.L.func.gii"),
  nii = file.path(data_dir, "melodic_IC_sum.nii.gz"),
  rds = file.path(data_dir, "melodic_IC_100.4k.rds")
)

xii1 <- select_xifti(read_cifti(template_fname["cii"]), 1) * 0

# basic -----
pr_cii <- estimate_prior(
  cii_fnames[seq(3)], template = template_fname["cii"], TR=.72, FC=FALSE,
  brainstructures=c("left", "right")
)
pr_cii

# `hpf` -----
# `hpf` arg, HPF in `nuisance`, and manually doing the HPF give the same result.
pr_ciiA <- estimate_prior(
  cii_fnames[seq(3)], template = template_fname["cii"], TR=.72, FC=FALSE,
  brainstructures=c("left", "right"), hpf=.01
)

bold <- lapply(cii_fnames[seq(3)], read_cifti, brainstructures=c("left", "right"))
nreg <- lapply(bold, function(q){
  fMRItools::dct_bases(
    ncol(q),
    round(fMRItools::dct_convert(ncol(q), TR=.72, f=.01))
  )
})
pr_ciiB <- estimate_prior(
  bold, nuisance = nreg, template = template_fname["cii"], TR=.72, FC=FALSE,
  brainstructures=c("left", "right")
)

for (bb in seq(length(bold))) {
  bold[[bb]] <- newdata_xifti(bold[[bb]], nuisance_regression(bold[[bb]], cbind(1, nreg[[bb]])))
}
pr_ciiC <- estimate_prior(
  bold, template = template_fname["cii"], TR=.72, FC=FALSE,
  brainstructures=c("left", "right")
)

testthat::expect_equal(pr_ciiA$prior, pr_ciiB$prior)
testthat::expect_equal(pr_ciiB$prior, pr_ciiC$prior)


# `drop_first` -----
# `drop_first` arg, scrubbing, and manually removing give the same result.
## with hpf -----
pr_ciiA <- estimate_prior(
  cii_fnames[seq(3)], template = template_fname["cii"], TR=.72, FC=TRUE,
  brainstructures="left", hpf=.01, drop_first=234, Q2=2,
  FC_nPivots=10, FC_nSamp=1000
)

bMap_cii <- fit_BBM(
  cii_fnames[4], pr_ciiA, brainstructures="left", maxiter=5, TR="prior", resamp_res=2000
)

bold <- lapply(cii_fnames[seq(3)], read_cifti, brainstructures="left")
bold <- lapply(bold, select_xifti, seq(235, 1200))
pr_ciiB <- estimate_prior(
  bold, template = template_fname["cii"], TR=.72, FC=TRUE,
  brainstructures="left", hpf=.01, Q2=2,
  FC_nPivots=10, FC_nSamp=1000
)

## w/o hpf -----
pr_cii0 <- estimate_prior(
  cii_fnames[seq(2)], template = template_fname["cii"], TR=.72, FC=FALSE,
  brainstructures="right", drop_first=0
)

pr_ciiA <- estimate_prior(
  cii_fnames[seq(2)], template = template_fname["cii"], TR=.72, FC=FALSE,
  brainstructures="right", drop_first=2
)

pr_ciiB <- estimate_prior(
  cii_fnames[seq(2)], template = template_fname["cii"], TR=.72, FC=FALSE,
  brainstructures="right",
  scrub=list(seq(2), seq(2))
)

pr_ciiC <- estimate_prior(
  cii_fnames[seq(2)], template = template_fname["cii"], TR=.72, FC=FALSE,
  brainstructures="right", drop_first=1,
  scrub=list(seq(2), seq(2))
)

## w/o hpf v2 -----
pr_ciiA <- estimate_prior(
  cii_fnames[seq(3)], template = template_fname["cii"], TR=.72, FC=FALSE,
  brainstructures="right", drop_first=1, scale_sm_FWHM=Inf, GSR=TRUE,
  scrub=list(c(7, 17, 777), NULL, 5)
)

pr_ciiB <- estimate_prior(
  cii_fnames[seq(3)], template = template_fname["cii"], TR=.72, FC=FALSE,
  brainstructures="right", scale_sm_FWHM=Inf, GSR=TRUE,
  scrub=list(c(1, 7, 17, 777), 1, c(1,5))
)

bold <- lapply(cii_fnames[seq(3)], read_cifti, brainstructures="right")
my_idx <- seq(2, 1200)[!(seq(2, 1200) %in% c(7, 17, 777))]
bold[[1]] <- select_xifti(bold[[1]], my_idx)
bold[[2]] <- select_xifti(bold[[2]], seq(2, 1200))
bold[[3]] <- select_xifti(bold[[3]], c(seq(2,4), seq(6, 1200)))
pr_ciiC <- estimate_prior(
  bold, template = template_fname["cii"], TR=.72, FC=FALSE,
  brainstructures="right", scale_sm_FWHM=Inf, GSR=TRUE
)

# w/ BOLD2 -----
pr_ciiA <- estimate_prior(
  cii_fnames[seq(3)], BOLD2=cii_fnames[c(4,5,6)],
  template = template_fname["cii"], TR=.72, FC=TRUE,
  brainstructures="right", hpf=.01, Q2=NULL,
  FC_nPivots=10, FC_nSamp=1000
)

pr_ciiB <- estimate_prior(
  cii_fnames[seq(4,6)], BOLD2=cii_fnames[seq(3)],
  template = template_fname["cii"], TR=.72, FC=TRUE,
  brainstructures="right", hpf=.01, Q2=NULL,
  FC_nPivots=10, FC_nSamp=1000
)

stopifnot(
  max(as.matrix(pr_ciiA$prior$mean) - as.matrix(pr_ciiB$prior$mean)) < 1e-8
)

bMap_cii <- fit_BBM(
  cii_fnames[4], pr_ciiA, maxiter=5, TR="prior", resamp_res=2000
)
