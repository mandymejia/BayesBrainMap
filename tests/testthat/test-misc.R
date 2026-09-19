test_that("Data normalization completes", {
  nv <- 1000
  nt <- 170
  x <- fMRItools::norm_BOLD(matrix(rnorm(nv*nt), nrow=nv), scale_by="sd", scale_sm_FWHM=0)
  testthat::expect_true(all(dim(x) == c(nv, nt)))
})
