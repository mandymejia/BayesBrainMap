# ---- Simulation --------------------------------------------------------------
set.seed(42)

nV <- 500
nT <- 150
nQs <- 5     # signal networks (in template)
nQn <- 2     # nuisance ICs (not in template)
nN  <- 6     # subjects used to build the prior
drop_first <- 3
spikes <- c(10, 40, 90)

sim_GICA <- function(nV, nQs, nQn){
  nQ <- nQs + nQn
  M <- matrix(rnorm(nV*nQ, sd=.3), nrow=nV)
  blk <- floor(nV/nQ)
  for (q in seq(nQ)) {
    idx <- ((q-1)*blk + 1):(q*blk)
    M[idx, q] <- M[idx, q] + 2
  }
  list(S = M[, seq(nQs), drop=FALSE], N = M[, nQs + seq(nQn), drop=FALSE])
}

sim_tc <- function(nQ, nT){
  t(sapply(seq(nQ), function(q){
    as.numeric(stats::filter(rnorm(nT + 10), rep(1/3, 3), sides=1))[-seq(10)]
  }))
}

sim_BOLD <- function(nV, nT, gica, nuis=NULL, spikes=NULL,
                     sd_noise=1, amp_nuisIC=2, amp_trend=3){
  nQs <- ncol(gica$S); nQn <- ncol(gica$N)
  X <- gica$S %*% sim_tc(nQs, nT) +
    amp_nuisIC * gica$N %*% sim_tc(nQn, nT) +
    matrix(rnorm(nV*nT, sd=sd_noise), nrow=nV)
  if (!is.null(nuis)) {
    nuis <- nuis[seq(nT), , drop=FALSE]
    X <- X + matrix(rnorm(nV*ncol(nuis), sd=2), nrow=nV) %*% t(nuis)
  }
  tt <- (seq(nT) - .5) / nT
  trends <- rbind(cos(pi*tt), cos(2*pi*tt))
  X <- X + matrix(rnorm(nV*2, sd=amp_trend), nrow=nV) %*% trends
  if (!is.null(spikes)) {
    X[, spikes] <- X[, spikes] + matrix(rnorm(nV*length(spikes), sd=15), nrow=nV)
  }
  X + 100
}

gica <- sim_GICA(nV, nQs, nQn)
tmp  <- gica$S

# Subject-specific spatial maps -> nonzero between-subject variance in the prior.
sim_group <- function(nN){
  B1 <- B2 <- vector("list", nN)
  for (i in seq(nN)) {
    g_i <- list(S = gica$S + matrix(rnorm(nV*nQs, sd=.5), nrow=nV), N = gica$N)
    B1[[i]] <- sim_BOLD(nV, nT, g_i, amp_trend=0, amp_nuisIC=1)
    B2[[i]] <- sim_BOLD(nV, nT, g_i, amp_trend=0, amp_nuisIC=1)
  }
  list(B1=B1, B2=B2)
}
grp <- sim_group(nN)

# The prior used throughout: matrix format, no FC (kept fast), TR/hpf fixed so
# that fit_BBM's "prior" sentinel arguments resolve without needing xifti
# metadata.
prior0 <- estimate_prior(
  grp$B1, grp$B2, template=tmp,
  scale_by="mean", scale_sm_FWHM=0,
  TR=1, hpf=0,
  FC=FALSE, verbose=FALSE
)

# A small-FC prior, used only by the FC-specific tests (kept cheap).
prior_FC <- estimate_prior(
  grp$B1, grp$B2, template=tmp,
  scale_by="mean", scale_sm_FWHM=0,
  TR=1, hpf=0,
  FC=TRUE, FC_nPivots=4, FC_nSamp=40, verbose=FALSE
)

# A prior with a nontrivial `mask` (one subject has a flat location).
bd_B1 <- grp$B1; bd_B1[[2]][100, ] <- 100  # flat voxel for subject 2 only
prior_masked <- estimate_prior(
  bd_B1, grp$B2, template=tmp,
  scale_by="mean", scale_sm_FWHM=0,
  TR=1, hpf=0,
  FC=FALSE, missingTol=.1, verbose=FALSE
)

# New (held-out) subject data to run `fit_BBM` on.
g_new <- list(S = gica$S + matrix(rnorm(nV*nQs, sd=.5), nrow=nV), N = gica$N)
nuis_new <- scale(matrix(cumsum(rnorm(nT*2)), nrow=nT))

Bnew  <- sim_BOLD(nV, nT, g_new, nuis=nuis_new, amp_trend=0, amp_nuisIC=1)
Bnew2 <- sim_BOLD(nV, nT, g_new, amp_trend=0, amp_nuisIC=1)  # a 2nd session, no nuisance

# Quiet wrapper. Defaults are fast (no dimension reduction, no denoising, no
# FC) so most tests run quickly; override anything via `...`.
run_fit <- function(BOLD, prior=prior0, ...){
  args <- list(
    BOLD=BOLD, prior=prior,
    scale_by="prior", scale_sm_FWHM="prior",
    drop_first="prior",
    hpf="prior", TR="prior", GSR="prior",
    Q2=0, Q2_max=NULL,
    reduce_dim=FALSE, method_FC="none",
    usePar=FALSE, verbose=FALSE
  )
  dots <- list(...)
  args[names(dots)] <- dots
  if (!all(args[c("scale_by", "scale_sm_FWHM", "drop_first", "hpf", "TR", "GSR")] == "prior")) {
    suppressWarnings(do.call(BayesBrainMap::fit_BBM, args))
  } else {
    do.call(BayesBrainMap::fit_BBM, args)
  }
}

# ---- Basic structure -----------------------------------------------------

test_that("single-session matrix input: structure and dimensions", {
  res <- run_fit(Bnew)
  expect_true(inherits(res, "bMap.matrix"))
  expect_equal(dim(res$subjNet_mean), c(nV, nQs))
  expect_equal(dim(res$subjNet_se),   c(nV, nQs))
  expect_true(all(is.finite(res$subjNet_mean)))
  expect_true(all(res$subjNet_se >= 0))
  expect_equal(dim(res$result_DR$A), c(nT, nQs))
  expect_equal(dim(res$result_DR$S), c(nQs, nV))
  expect_named(res$comptime, c("DR", "bMap", "FC-bMap", "sbMap"))
  expect_true(res$comptime["DR"] >= 0)
  expect_true(res$comptime["bMap"] >= res$comptime["DR"])
})

test_that("doMLE=TRUE returns an MLE map; doMLE=FALSE omits it", {
  res  <- run_fit(Bnew, doMLE=TRUE)
  res2 <- run_fit(Bnew, doMLE=FALSE)
  expect_equal(dim(res$MLE), c(nV, nQs))
  expect_null(res2$MLE)
})

test_that("estimated subjNet_mean correlates with the true (simulated) subject maps", {
  res <- run_fit(Bnew, scale_by="none")
  r <- diag(cor(res$subjNet_mean, g_new$S))
  expect_true(all(r > .3), info=paste("r =", paste(round(r, 2), collapse=", ")))
})

test_that("results are reproducible given identical inputs", {
  a <- run_fit(Bnew, drop_first=drop_first, nuisance=nuis_new, scrub=spikes)
  b <- run_fit(Bnew, drop_first=drop_first, nuisance=nuis_new, scrub=spikes)
  expect_equal(a$subjNet_mean, b$subjNet_mean)
  expect_equal(a$subjNet_se,   b$subjNet_se)
})

# ---- var_method / scale_by / GSR ------------------------------------------

test_that("var_method changes the prior variance used, and thus the estimates", {
  res_nn <- run_fit(Bnew, var_method="non-negative")
  res_ub <- run_fit(Bnew, var_method="unbiased")
  expect_false(isTRUE(all.equal(res_nn$subjNet_mean, res_ub$subjNet_mean)))
  expect_equal(res_nn$params$var_method, "non-negative")
  expect_equal(res_ub$params$var_method, "unbiased")
})

test_that("scale_by='prior' resolves to the value used for the prior", {
  res <- run_fit(Bnew, scale_by="prior")
  expect_equal(res$params$scale_by, prior0$params$scale_by)
})

test_that("scale_by options run and change results", {
  res_mean <- run_fit(Bnew, scale_by="mean")
  res_sd   <- run_fit(Bnew, scale_by="sd")
  res_none <- run_fit(Bnew, scale_by="none")
  for (r in list(res_mean, res_sd, res_none)) {
    expect_equal(dim(r$subjNet_mean), c(nV, nQs))
  }
  expect_false(isTRUE(all.equal(res_mean$subjNet_mean, res_none$subjNet_mean)))
})

test_that("GSR=TRUE runs", {
  res <- run_fit(Bnew, GSR=TRUE)
  expect_equal(dim(res$subjNet_mean), c(nV, nQs))
})

test_that("mismatched parameters relative to the prior trigger a warning", {
  expect_warning(run_fit(Bnew, varTol=1), "varTol")
})

# ---- reduce_dim -------------------------------------------------------------

test_that("reduce_dim TRUE and FALSE both run and give correlated results", {
  res_full <- run_fit(Bnew, reduce_dim=FALSE)
  res_red  <- run_fit(Bnew, reduce_dim=TRUE)
  expect_equal(dim(res_red$subjNet_mean), c(nV, nQs))
  r <- diag(cor(res_full$subjNet_mean, res_red$subjNet_mean))
  expect_true(all(r > .5))
})

# ---- Q2 denoising ------------------------------------------------------------

test_that("Q2 denoising runs and its estimate is recorded", {
  res <- run_fit(Bnew, Q2=nQn)
  expect_equal(res$params$Q2_est, nQn)
  expect_equal(dim(res$subjNet_mean), c(nV, nQs))
})

test_that("Q2=NULL (PESEL) runs", {
  res <- run_fit(Bnew, Q2=NULL, Q2_max=nQn + 2)
  expect_equal(dim(res$subjNet_mean), c(nV, nQs))
})

test_that("Q2 denoising with multiple BOLD sessions re-normalizes every session", {
  res <- run_fit(list(Bnew, Bnew2), Q2=nQn)
  expect_equal(dim(res$subjNet_mean), c(nV, nQs))
})

# ---- drop_first / nuisance / scrub ------------------------------------------

test_that("drop_first drops the requested number of volumes", {
  res <- run_fit(Bnew, drop_first=drop_first)
  expect_equal(nrow(res$result_DR$A), nT - drop_first)
  expect_equal(res$params$drop_first, drop_first)
})

test_that("drop_first with no nuisance/scrub supplied does not error", {
  expect_no_error(run_fit(Bnew, drop_first=drop_first))
})

test_that("drop_first combined with scrub: correct number of volumes retained", {
  sc <- spikes[spikes > drop_first]
  res <- run_fit(Bnew, drop_first=drop_first, scrub=sc)
  expect_equal(nrow(res$result_DR$A), nT - drop_first - length(sc))
})

test_that("a single (non-list) nuisance/scrub applies correctly to a single-session run", {
  res_plain <- run_fit(Bnew, nuisance=nuis_new, scrub=spikes)
  res_list  <- run_fit(Bnew, nuisance=list(nuis_new), scrub=list(spikes))
  expect_equal(res_plain$subjNet_mean, res_list$subjNet_mean, tolerance=1e-6)
})

test_that("multiple BOLD sessions run with no nuisance/scrub (defaults)", {
  res <- run_fit(list(Bnew, Bnew2))
  expect_equal(dim(res$subjNet_mean), c(nV, nQs))
  expect_equal(nrow(res$result_DR$A), nT * 2)
})

test_that("multiple sessions: per-session nuisance regressors are applied to their own session", {
  # If fit_BBM normalized only the last session (stale loop index), swapping
  # which nuisance matrix goes with which session would have NO effect on the
  # result, since only nuisance[[2]] would ever be used either way.
  nuis1 <- scale(matrix(cumsum(rnorm(nT*2)), nrow=nT))
  nuis2 <- scale(matrix(cumsum(rnorm(nT*2)), nrow=nT))
  res_both <- run_fit(list(Bnew, Bnew2), nuisance=list(nuis1, nuis2))
  res_swap <- run_fit(list(Bnew, Bnew2), nuisance=list(nuis2, nuis1))
  expect_false(isTRUE(all.equal(res_both$subjNet_mean, res_swap$subjNet_mean)))
})

# ---- hpf / TR -----------------------------------------------------------------

test_that("hpf requested without a usable TR errors", {
  expect_error(run_fit(Bnew, hpf=.02, TR=NULL))
})

test_that("TR='prior' resolves to the TR used for the prior", {
  res <- run_fit(Bnew, TR="prior")
  expect_equal(res$params$TR, prior0$params$TR)
})

# ---- covariates ---------------------------------------------------------------

test_that("covariates provided but none used for the prior: errors", {
  expect_error(run_fit(Bnew, covariates=c(age=30)))
})

# ---- masking ------------------------------------------------------------------

test_that("flat locations in new BOLD are masked out and don't affect other locations", {
  B <- Bnew
  B[1:3, ] <- 100   # flat -> masked by varTol
  bad <- 1:3
  res <- run_fit(B)
  expect_true(all(is.na(res$subjNet_mean[bad, ])))
  expect_false(anyNA(res$subjNet_mean[-bad, ]))
})

test_that("prior$mask locations are propagated through and unmasked in the output", {
  expect_false(is.null(prior_masked$mask))
  res <- run_fit(Bnew, prior=prior_masked)
  expect_equal(dim(res$subjNet_mean), c(nV, nQs))
  expect_true(all(is.na(res$subjNet_mean[!prior_masked$mask, ])))
  expect_false(anyNA(res$subjNet_mean[prior_masked$mask, ]))
})

# ---- FC -----------------------------------------------------------------------

test_that("method_FC='VB1' runs using an FC-enabled prior", {
  res <- expect_warning(
    run_fit(Bnew, prior=prior_FC, method_FC="VB1", miniter=2, maxiter=3),
    "Failed to converge"
  )
  expect_true(res$params$FC)
  expect_false(is.null(res$result_bMap))
  expect_equal(dim(res$subjNet_mean), c(nV, nQs))
})

test_that("method_FC requested but prior lacks FC info: warns and falls back to 'none'", {
  expect_warning(res <- run_fit(Bnew, prior=prior0, method_FC="VB1"), "FC information")
  expect_false(res$params$FC)
})

# ---- spatial_model (matrix format: meshes required explicitly) ---------------

test_that("`spatial_model=TRUE` for matrix data requires explicit meshes", {
  expect_error(run_fit(Bnew, spatial_model=TRUE), "meshes|INLA")
})

# ---- Bad inputs -----------------------------------------------------------

test_that("bad inputs error", {
  expect_error(run_fit(Bnew[-1, ]))              # wrong nV vs. prior
  expect_error(run_fit(Bnew, var_method="nonsense"))
  expect_error(run_fit(Bnew, method_FC="nonsense"))
})

test_that("BOLD/prior format mismatch errors", {
  fake_prior <- prior0
  class(fake_prior) <- "prior.cifti"
  expect_error(run_fit(Bnew, prior=fake_prior), "format")
})
