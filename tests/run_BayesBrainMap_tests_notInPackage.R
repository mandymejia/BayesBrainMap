# Build --> Install and Restart

# Setup ------------------------------------------------------------------------
# ciftiTools
library(ciftiTools)
print(packageVersion("ciftiTools"))
ciftiTools.setOption("wb_path", "~/Applications")

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

# Quick little check of the three main functions, w/ CIFTI ---------------------
pr_cii <- estimate_prior(
  cii_fnames[seq(3)], template = template_fname["cii"], TR=.72, FC=FALSE,
  brainstructures=c("left", "right")
)
# # w/ mask
# pr2_cii <- estimate_prior(
#   cii_fnames[seq(3)], template = template_fname["cii"], TR=.72, FC=FALSE, varTol=5000,
#   brainstructures=c("left", "right"), mask=c(rep(FALSE, 144), rep(TRUE, 7100), rep(FALSE, 100))
# )
#
# pr_cii <- estimate_prior(
#   cii_fnames[seq(3)], template = template_fname["cii"], TR=.72, FC=FALSE, varTol=4000,
#   brainstructures=c("left", "right")
# )
# # w/ mask
# pr2_cii <- estimate_prior(
#   cii_fnames[seq(3)], template = template_fname["cii"], TR=.72, FC=FALSE, varTol=4000,
#   brainstructures=c("left", "right"), mask=c(rep(FALSE, 144), rep(TRUE, 7100), rep(FALSE, 100))
# )

bMap_cii <- BrainMap(
  cii_fnames[4], pr_cii, brainstructures="left", maxiter=5, TR="prior", resamp_res=2000
)

bMapA_cii <- BrainMap(
  cii_fnames[4], pr_cii, brainstructures="left", maxiter=5, TR="prior", resamp_res=2000,
  scrub=c(seq(7), 10, 51, 78), hpf=0
)

bMapB_cii <- BrainMap(
  cii_fnames[4], pr_cii, brainstructures="left", maxiter=5, TR="prior", resamp_res=2000,
  drop_first=6, scrub=c(3, 7, 10, 51, 78), hpf=0
)

z <- read_cifti(cii_fnames[4], resamp_res=2000)
bMap_cii2 <- BrainMap(
  z, pr_cii, brainstructures="left", maxiter=5, TR="prior", resamp_res=2000
)

engMap_cii <- engagements(bMap_cii)
engMap_cii <- engagements(bMap_cii, z=c(0, .1, 3, 11))

engMap_cii
plot(engMap_cii)
close3d()

pr_cii2 <- estimate_prior(
  cii_fnames[seq(3)], template = template_fname["cii"], TR=.72, FC=TRUE,
  brainstructures=c("left", "right")
)
bMap_cii2 <- BrainMap(
  cii_fnames[4], pr_cii2, brainstructures="left", maxiter=5, TR="prior", resamp_res=2000
) # ***summary TRUE

# `estimate_prior`: check for same result w/ different file types -----------
### Test 1: basic ----
pr_cii <- estimate_prior(
  cii_fnames[seq(5)], brainstructures="left", template = template_fname["cii"],
  FC=TRUE, TR=.72, scale_sm_FWHM=0
)
# pr_gii <- estimate_prior(
#   giiL_fnames[seq(5)], template = template_fname["gii"],
#   FC=TRUE, TR=.72, scale_sm_FWHM=0
# )
# pr_rds <- estimate_prior(
#   rds_fnames[seq(5)], template = template_fname["rds"],
#   FC=TRUE, TR=.72, scale_sm_FWHM=0
# )
# testthat::expect_equal(
#   lapply(pr_cii$prior[seq(3)], fMRItools::unmask_mat, pr_cii$dat_struct$meta$cortex$medial_wall_mask$left),
#   pr_gii$prior[seq(3)]
# )
# testthat::expect_equal(pr_cii$prior, pr_rds$prior)

pr_cii#; pr_gii; pr_rds
plot(pr_cii)#; plot(pr_gii)

### Test 2: with various parameters changed ----
pr_cii <- estimate_prior(
  cii_fnames[seq(3)], cii_fnames[seq(4,6)], template = template_fname["cii"],
  inds=c(2,7,11,90), scale="global", scale_sm_FWHM=5,
  maskTol=.9, brainstructures="left", wb_path="~/Applications",
  usePar=TRUE, FC=TRUE, varTol=10000
)
# pr_gii <- estimate_prior(
#   giiL_fnames[seq(3)], giiL_fnames[seq(4,6)], template = template_fname["gii"],
#   inds=c(2,7,11,90), scale="global", scale_sm_FWHM=5,
#   maskTol=.9, wb_path="~/Applications",
#   usePar=TRUE, FC=TRUE, varTol=10000
# )
# testthat::expect_equal(
#   lapply(pr_cii$prior[seq(3)], fMRItools::unmask_mat, pr_cii$dat_struct$meta$cortex$medial_wall_mask$left),
#   pr_gii$prior[seq(3)]
# )
# pr_gii <- estimate_prior(
#   giiL_fnames[seq(3)], giiL_fnames[seq(4,6)], template = template_fname["gii"],
#   inds=seq(5), scale="none", Q2=5,
#   maskTol=.9, wb_path="~/Applications",
#   usePar=TRUE
# )
# pr_rds <- estimate_prior(
#   lapply(rds_fnames[seq(4,6)], readRDS), lapply(rds_fnames[seq(3)], readRDS), template = template_fname["rds"],
#   inds=seq(5), scale="none",  Q2=5,
#   maskTol=.9, usePar=TRUE
# )
# testthat::expect_equal(
#   pr_gii$prior[seq(3)],
#   lapply(pr_rds$prior[seq(3)], fMRItools::unmask_mat, pr_cii$dat_struct$meta$cortex$medial_wall_mask$left),
# )

close3d(); close3d(); close3d(); close3d()

# `export_prior` and `BrainMap`: check for same result w/ different file types -----------------
pr_cii <- estimate_prior(
  cii_fnames[seq(3)], brainstructures="left", template = template_fname["cii"], inds=seq(3),
  scale="global"
)
# pr_gii <- estimate_prior(
#   giiL_fnames[seq(3)], template = template_fname["gii"], inds=seq(3),
#   scale="global"
# )
# pr_rds <- estimate_prior(
#   rds_fnames[seq(3)], template = template_fname["rds"], inds=seq(3),
#   scale="global"
# )

### `export_prior` ----
out_fname=export_prior(pr_cii, tempfile())
pr_cii2 <- list(read_cifti(out_fname[1]), read_cifti(out_fname[2]), readRDS(out_fname[3]))
#out_fname=export_prior(pr_gii, tempfile())
# pr_gii2 <- list(readgii(out_fname[1]), readgii(out_fname[2]), readRDS(out_fname[3]))
# out_fname=export_prior(pr_rds, tempfile())
# pr_rds2 <- lapply(out_fname, readRDS)

### `BrainMap` ----
bMap_cii <- BrainMap(cii_fnames[4], brainstructures="left", pr_cii, maxiter=20, Q2=0, TR=.72)
#bMap_gii <- BrainMap(giiL_fnames[4], pr_gii, Q2=0, maxiter=20, TR=.72)
#bMap_rds <- BrainMap(rds_fnames[4], pr_rds, Q2=0, maxiter=20, TR=.72)
bMap_cii#; bMap_gii; bMap_rds
#testthat::expect_equal(bMap_gii$A, bMap_rds$A)
#engMap_rds <- engagements(bMap_rds)
engMap_cii <- engagements(bMap_cii)
plot(engagements(bMap_cii))#; plot(engagements(bMap_gii))
close3d(); close3d()

gamma <- 2
gamma_scaled <- gamma*sqrt(matrixStats::colVars(bMap_cii$prior_mean))
eng2 <- engagements(bMap_cii, u=gamma_scaled)
eng3 <- engagements(bMap_cii, z=gamma)
# breaks
eng2 <- cbind(eng2[[2]]$engaged[,1], eng2[[3]]$engaged[,2], eng2[[4]]$engaged[,3])
eng3 <- eng3[[2]]$engaged
testthat::expect_equal(eng2, eng3)

### Without FC ----
pr_cii <- estimate_prior(
  cii_fnames[seq(3)], brainstructures="left", template = template_fname["cii"], inds=seq(3),
  scale="none", TR=.72, FC=FALSE
)
# pr_rds <- estimate_prior(
#   rds_fnames[seq(3)], template = template_fname["rds"], inds=seq(3),
#   scale="none", TR=.72, FC=FALSE
# )
bMap_cii <- BrainMap(cii_fnames[4], brainstructures="left", pr_cii, maxiter=20, Q2=5, TR=.72)
#bMap_rds <- BrainMap(rds_fnames[4], pr_rds, Q2=5, maxiter=20, TR=.72)
#testthat::expect_equal(bMap_cii$theta_MLE, bMap_rds$theta_MLE)

# CIFTI ------------------------------------------------------------------------

### With subcortex ----
tm <- estimate_prior(
  cii_fnames[seq(4)], template=template_fname["cii"], scale=FALSE#, FC=TRUE
)
tm
plot(tm)
close3d(); close3d()

cii <- read_cifti(cii_fnames[5], brainstructures="all")
cii$data$cortex_left[33,] <- mean(cii$data$cortex_left[33,])
bMap <- BrainMap(cii, tm, scale=FALSE, miniter=2, maxiter=3, Q2=0)
plot(bMap)
close3d()

engMap <- engagements(bMap)
engMap_fname <- paste0(tempfile(), ".dlabel.nii")
engMap$engaged <- move_to_mwall(engMap$engaged)
engMap$engaged$meta$cortex$medial_wall_mask$right <- rep(TRUE, 4002)
write_cifti(engMap$engaged, engMap_fname)
engMap2 <- read_cifti(engMap_fname, brainstructures="all")
engMap2$meta$cortex$medial_wall_mask$right <- rep(TRUE, 4002)
plot(engMap); plot(engMap2)
close3d(); close3d()


rm(cii); rm(z); rm(engMap); rm(engMap2);
rm(tm); rm(pr_cii); rm(pr_cii2); rm(bMap); rm(bMap_cii); rm(bMap_cii2)
gc()

tm <- estimate_prior(
  cii_fnames[seq(3)], cii_fnames[seq(4, 6)],
  template=template_fname["cii"], scale="local",
  brainstructures="right", varTol=1, verbose=FALSE
)
tm
plot(tm, "var")
close3d()

tm2 <- estimate_prior(
  cii_fnames[seq(3)], cii_fnames[seq(4, 6)],
  template=template_fname["cii"], scale="local", scale_sm_FWHM=20,
  brainstructures="right", varTol=1, verbose=FALSE
)

cii <- lapply(cii_fnames[seq(4)], read_xifti, brainstructures="left")
cii[[1]]$data$cortex_left[3,17] <- NA
cii[[1]]$data$cortex_left[11,5] <- NA
cii[[2]]$data$cortex_left[11,seq(10)] <- NA
cii[[3]]$data$cortex_left[11,] <- NA
cii[[1]]$data$cortex_left[78,5] <- NA
cii[[2]]$data$cortex_left[78,seq(10)] <- NA
cii[[3]]$data$cortex_left[78,] <- NA
cii[[4]]$data$cortex_left[,] <- NA
tm <- estimate_prior(
  cii, template=read_cifti(template_fname["cii"], brainstructures="left"),
  FC=FALSE,
  scale="global", inds=c(1,4,7,11), maskTol = .5, missingTol=.5
)
tm
rm(cii)
plot(tm, idx=3)
close3d(); close3d()
cii <- read_cifti(cii_fnames[5], brainstructures="left")
#squarem1 error ... ?
# BrainMap(cii, tm, brainstructures="left", scale="global", maxiter=7, Q2=0, spatial_model = TRUE)
cii$data$cortex_left[33,] <- mean(cii$data$cortex_left[33,])
bMap <- testthat::expect_error( # Not supported yet: flat or NA voxels in data, after applying prior mask, with spatial model.
  BrainMap(cii, tm, brainstructures="left", scale="global", maxiter=7, Q2=0, spatial_model = TRUE, TR=.72)
)

cii <- lapply(cii_fnames[seq(4)], read_xifti, brainstructures="right")
cii0 <- lapply(cii, as.matrix)
cii0f <- paste0(c(tempfile(), tempfile(), tempfile(), tempfile()), ".rds")
for (ii in seq(4)) { saveRDS(cii0[[ii]], cii0f[ii]) }
tm <- estimate_prior(
  cii0f,
  template=as.matrix(read_cifti(template_fname["cii"], brainstructures="right")),
  scale=FALSE, inds=c(1,4,7,11)
)


# CIFTI pseudo retest vs data true retest: should get same results.
tm <- estimate_prior(
  cii,
  template=as.matrix(read_cifti(template_fname["cii"], brainstructures="right")),
  scale="global", inds=c(1,4,7,11),
)
tm2 <- estimate_prior(
  lapply(cii, function(x){as.matrix(x)[,seq(600)]}),
  lapply(cii, function(x){as.matrix(x)[,seq(601,1200)]}),
  template=as.matrix(read_cifti(template_fname["cii"], brainstructures="right")),
  scale="global", inds=c(1,4,7,11),
)
stopifnot(
  max(abs(do.call(c, tm$var_decomp) - do.call(c, tm2$var_decomp)), na.rm=TRUE) < 1e-8
)
rm(tm2)

# NIFTI ------------------------------------------------------------------------
rm(xii1)

# Load NIFTI group IC
ngIC_fname <- file.path(data_dir, "melodic_IC_sum.nii.gz")
ngIC <- readNifti(ngIC_fname)
nmask <- apply(ngIC!=0, seq(3), all)
mask_fname <- file.path(data_dir, "mask.nii.gz")
RNifti::writeNifti(nmask, mask_fname)

# mask erosion?
tm <- estimate_prior(
  nii_fnames[seq(4)], template=ngIC_fname, #scale=FALSE,
  mask=mask_fname, varTol = 500, maskTol=.3, missingTol=.9
)
tm
bMap <- BrainMap(
  nii_fnames[2], tm, scale=FALSE,
  miniter=1, maxiter=1, mask=mask_fname, Q2=0, TR=.72
)
bMap
eng <- engagements(bMap, u=.2)
eng
