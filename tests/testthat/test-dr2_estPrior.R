# written w/ Claude

# ---- Simulation --------------------------------------------------------------
set.seed(0)

nV <- 900
nT <- 200
nQs <- 7   # signal networks (in template)
nQn <- 3   # nuisance ICs (not in template)
drop_first <- 3
spikes <- c(10, 55, 120)  # all > drop_first

# Returns list(S = nV x nQs signal maps [the template], N = nV x nQn nuisance maps)
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

# Smooth-ish random timecourses (Q x T)
sim_tc <- function(nQ, nT){
  t(sapply(seq(nQ), function(q){
    as.numeric(stats::filter(rnorm(nT + 10), rep(1/3, 3), sides=1))[-seq(10)]
  }))
}

# nT x k nuisance matrix
nuis <- scale(matrix(cumsum(rnorm(nT*2)), nrow=nT))

sim_BOLD <- function(nV, nT, gica, nuis=NULL, spikes=NULL,
                     sd_noise=1, amp_nuisIC=2, amp_trend=3){
  nQs <- ncol(gica$S); nQn <- ncol(gica$N)
  # Signal + nuisance ICs + white noise
  X <- gica$S %*% sim_tc(nQs, nT) +
    amp_nuisIC * gica$N %*% sim_tc(nQn, nT) +
    matrix(rnorm(nV*nT, sd=sd_noise), nrow=nV)
  # Nuisance regressors (voxel-wise loadings)
  if (!is.null(nuis)) {
    nuis <- nuis[seq(nT), , drop=FALSE]
    X <- X + matrix(rnorm(nV*ncol(nuis), sd=2), nrow=nV) %*% t(nuis)
  }
  # Slow trends
  tt <- (seq(nT) - .5) / nT
  trends <- rbind(cos(pi*tt), cos(2*pi*tt))
  X <- X + matrix(rnorm(nV*2, sd=amp_trend), nrow=nV) %*% trends
  # Spike volumes
  if (!is.null(spikes)) {
    X[, spikes] <- X[, spikes] + matrix(rnorm(nV*length(spikes), sd=15), nrow=nV)
  }
  # Baseline mean (so mean-scaling is valid)
  X + 100
}

gica <- sim_GICA(nV, nQs, nQn)
tmp <- gica$S

B1s <- list(
  sim_BOLD(nV, nT, gica, nuis=nuis),
  sim_BOLD(nV, nT, gica, nuis=nuis, spikes=spikes),
  sim_BOLD(nV, nT, gica, spikes=spikes),
  sim_BOLD(nV, nT, gica)
)
B2s <- list(
  sim_BOLD(nV, nT, gica, nuis=nuis),
  sim_BOLD(nV, nT, gica, nuis=nuis, spikes=spikes),
  sim_BOLD(nV, nT, gica, spikes=spikes),
  sim_BOLD(nV, nT, gica)
)

nB <- length(B1s)

# Wrapper: quiet, sensible defaults (no scale smoothing since no surfaces).
run <- function(B1, B2=NULL, ...){
  args <- list(
    BOLD=B1, BOLD2=B2, format="data", template=tmp,
    TR=1, scale_by="sd", scale_sm_FWHM=0, verbose=FALSE
  )
  dots <- list(...)
  args[names(dots)] <- dots
  do.call(BayesBrainMap:::dual_reg2, args)
}

# ---- Tests -------------------------------------------------------------------

test_that("retest, no denoising: structure and dimensions", {
  dr <- run(B1s[[1]], B2s[[1]], drop_first=drop_first, keepA=TRUE,
            nuisance=list(nuis, nuis), hpf=.01)
  expect_named(dr, c("test", "retest"))
  for (s in c("test", "retest")) {
    expect_equal(dim(dr[[s]]$S), c(nQs, nV))
    expect_equal(dim(dr[[s]]$A), c(nT - drop_first, nQs))
    expect_length(dr[[s]]$sigma_sq, nV)
    expect_true(all(is.finite(dr[[s]]$sigma_sq)))
    expect_true(all(dr[[s]]$sigma_sq > 0))
  }
})

test_that("keepA=FALSE drops A and A2", {
  dr <- run(B1s[[4]], B2s[[4]], keepA=FALSE)
  for (s in c("test", "retest")) {
    expect_null(dr[[s]]$A)
    expect_null(dr[[s]]$A2)
    expect_false(is.null(dr[[s]]$S))
  }
})

test_that("scrub + drop_first: correct number of volumes retained", {
  # Spike at index 2 is inside the dropped window -> not counted.
  sc <- c(2, spikes)
  dr <- run(B1s[[2]], B2s[[2]], drop_first=drop_first, keepA=TRUE,
            nuisance=list(nuis, nuis), scrub=list(sc, sc))
  expect_equal(nrow(dr$test$A), nT - drop_first - length(spikes))
  expect_equal(nrow(dr$retest$A), nT - drop_first - length(spikes))
})

test_that("retest: NULL list entries for nuisance/scrub are allowed", {
  dr <- run(B1s[[1]], B2s[[3]], drop_first=drop_first, keepA=TRUE,
            nuisance=list(nuis, NULL), scrub=list(NULL, spikes))
  expect_equal(nrow(dr$test$A), nT - drop_first)
  expect_equal(nrow(dr$retest$A), nT - drop_first - length(spikes))
})

test_that("retest scans may have different numbers of timepoints", {
  nT2 <- 150
  B2 <- sim_BOLD(nV, nT2, gica, nuis=nuis)
  dr <- run(B1s[[1]], B2, drop_first=drop_first, keepA=TRUE,
            nuisance=list(nuis, nuis[seq(nT2), ]))
  expect_equal(nrow(dr$test$A), nT - drop_first)
  expect_equal(nrow(dr$retest$A), nT2 - drop_first)
})

test_that("pseudo-retest: BOLD is split in half", {
  dr <- run(B1s[[2]], NULL, drop_first=drop_first, keepA=TRUE,
            nuisance=nuis, scrub=spikes, hpf=.01)
  n <- nT - drop_first - length(spikes)
  expect_equal(nrow(dr$test$A), round(n/2))
  expect_equal(nrow(dr$retest$A), n - round(n/2))
  expect_equal(dim(dr$test$S), c(nQs, nV))
})

test_that("matches manual norm_BOLD + dual_reg (retest, no denoising)", {
  args <- list(drop_first=drop_first, nuisance=nuis, scrub=spikes,
               TR=1, hpf=.01, scale_by="sd", scale_sm_FWHM=0)
  manual <- function(B){
    Bn <- do.call(fMRItools::norm_BOLD, c(list(BOLD=B), args))
    fMRItools::dual_reg(BOLD=Bn, GICA=tmp, scale_by="none", hpf=0)
  }
  dr <- run(B1s[[2]], B2s[[2]], drop_first=drop_first, hpf=.01,
            nuisance=list(nuis, nuis), scrub=list(spikes, spikes))
  expect_equal(dr$test$S, manual(B1s[[2]])$S, tolerance=1e-6)
  expect_equal(dr$retest$S, manual(B2s[[2]])$S, tolerance=1e-6)
})

test_that("scale_by options run and change results", {
  d_sd   <- run(B1s[[4]], B2s[[4]], scale_by="sd")
  d_mean <- run(B1s[[4]], B2s[[4]], scale_by="mean")
  d_none <- run(B1s[[4]], B2s[[4]], scale_by="none")
  d_glob <- run(B1s[[4]], B2s[[4]], scale_by="sd", scale_sm_FWHM=Inf)
  for (d in list(d_sd, d_mean, d_none, d_glob)) {
    expect_equal(dim(d$test$S), c(nQs, nV))
  }
  expect_false(isTRUE(all.equal(d_sd$test$S, d_none$test$S)))
  expect_false(isTRUE(all.equal(d_sd$test$S, d_glob$test$S)))
})

test_that("local smoothing without surface data messages and continues", {
  expect_message(
    run(B1s[[4]], B2s[[4]], scale_by="sd", scale_sm_FWHM=4),
    "No surface"
  )
})

test_that("hpf without TR errors; TR sentinel handled", {
  expect_error(run(B1s[[4]], B2s[[4]], TR=NULL, hpf=.02))
  expect_error(run(B1s[[4]], B2s[[4]], TR="from_xifti_metadata", hpf=.01))
  expect_no_error(run(B1s[[4]], B2s[[4]], TR="from_xifti_metadata", hpf=NULL))
})

test_that("low-variance and NA locations are masked and unmasked", {
  B1 <- B1s[[4]]; B2 <- B2s[[4]]
  B1[1:3, ] <- 100            # constant voxels in scan 1
  B2[10, 5] <- NA             # missing value in scan 2
  bad <- c(1:3, 10)
  dr <- run(B1, B2, keepA=TRUE)
  for (s in c("test", "retest")) {
    expect_equal(dim(dr[[s]]$S), c(nQs, nV))
    expect_true(all(is.na(dr[[s]]$S[, bad])))
    expect_false(anyNA(dr[[s]]$S[, -bad]))
    expect_length(dr[[s]]$sigma_sq, nV)
    expect_true(all(is.na(dr[[s]]$sigma_sq[bad])))
    expect_false(anyNA(dr[[s]]$sigma_sq[-bad]))
  }
  # Same as running on the clean subset directly.
  keep <- setdiff(seq(nV), bad)
  dr2 <- run(B1[keep, ], B2[keep, ], template=tmp[keep, ])
  expect_equal(dr$test$S[, keep], dr2$test$S, tolerance=1e-6)
  expect_equal(dr$retest$sigma_sq[keep], dr2$retest$sigma_sq, tolerance=1e-6)
})

test_that("too many masked locations returns NULL", {
  B1 <- B1s[[4]]
  B1[1:(nV*.2), ] <- 100   # 20% masked > default maskTol (10%)
  expect_null(run(B1, B2s[[4]]))
  expect_false(is.null(run(B1, B2s[[4]], maskTol=.3)))
})

test_that("`mask` argument: template is pre-masked, output matches subsetting", {
  mask <- rep(TRUE, nV); mask[1:50] <- FALSE
  dr <- run(B1s[[4]], B2s[[4]], mask=mask, template=tmp[mask, ])
  expect_equal(dim(dr$test$S), c(nQs, sum(mask)))
  dr2 <- run(B1s[[4]][mask, ], B2s[[4]][mask, ], template=tmp[mask, ])
  expect_equal(dr$test$S, dr2$test$S, tolerance=1e-6)
  # Mismatched template/mask
  expect_error(run(B1s[[4]], B2s[[4]], mask=mask))
})

test_that("bad inputs error", {
  expect_error(run(B1s[[4]], B2s[[4]], nuisance=nuis))          # not a length-2 list
  expect_error(run(B1s[[4]], B2s[[4]], scrub=spikes))           # not a length-2 list
  expect_error(run(B1s[[4]], B2s[[4]], drop_first=nT))          # drops everything
  expect_error(run(B1s[[4]], B2s[[4]], format="nonsense"))
  expect_error(run(B1s[[4]], B2s[[4]], template=tmp[-1, ]))     # wrong nV
  expect_error(run(B1s[[4]], B2s[[4]][-1, ]))                   # BOLD2 spatial mismatch
})

test_that("denoising (Q2): structure, dimensions, preclean fields", {
  dr <- run(B1s[[1]], B2s[[1]], drop_first=drop_first, keepA=TRUE,
            nuisance=list(nuis, nuis), Q2=nQn)
  expect_true(all(c("test", "retest", "test_preclean", "retest_preclean") %in% names(dr)))
  for (s in c("test", "retest", "test_preclean", "retest_preclean")) {
    expect_equal(dim(dr[[s]]$S), c(nQs, nV))
    expect_length(dr[[s]]$sigma_sq, nV)
    expect_true(all(is.finite(dr[[s]]$sigma_sq)))
  }
  # Removing nuisance ICs should reduce residual variance overall.
  expect_lt(mean(dr$test$sigma_sq), mean(dr$test_preclean$sigma_sq))
})

test_that("denoising, pseudo-retest and Q2=NULL (PESEL) run", {
  dr1 <- run(B1s[[2]], NULL, drop_first=drop_first, nuisance=nuis,
             scrub=spikes, Q2=nQn, keepA=TRUE)
  n <- nT - drop_first - length(spikes)
  expect_equal(nrow(dr1$test$A), round(n/2))
  expect_equal(nrow(dr1$retest$A), n - round(n/2))
  dr2 <- run(B1s[[4]], B2s[[4]], Q2=NULL)
  expect_equal(dim(dr2$test$S), c(nQs, nV))
})

test_that("denoising + masking: outputs unmasked to full length", {
  B1 <- B1s[[4]]; B2 <- B2s[[4]]
  B1[1:3, ] <- 100
  for (s in c("test", "retest", "test_preclean", "retest_preclean")) {
    dr <- run(B1, B2, Q2=nQn)
    expect_equal(dim(dr[[s]]$S), c(nQs, nV))
    expect_true(all(is.na(dr[[s]]$sigma_sq[1:3])))
  }
})

test_that("Q2=0 and Q2_max=0 both skip denoising", {
  expect_named(run(B1s[[4]], B2s[[4]], Q2=0), c("test", "retest"))
  expect_named(run(B1s[[4]], B2s[[4]], Q2=NULL, Q2_max=0), c("test", "retest"))
})

test_that("estimated S recover true signal maps (clean data)", {
  set.seed(1)
  # Cleaner data: 2x signal maps, lower noise, milder nuisance ICs, no slow trends.
  # (Correlation with the template is invariant to the 2x scaling.)
  gica_c <- list(S = gica$S * 2, N = gica$N)
  Bc1 <- sim_BOLD(nV, nT, gica_c, sd_noise=.5, amp_nuisIC=1, amp_trend=0)
  Bc2 <- sim_BOLD(nV, nT, gica_c, sd_noise=.5, amp_nuisIC=1, amp_trend=0)

  # No scaling: SD scaling flattens the block contrast between signal/non-signal voxels.
  dr <- run(Bc1, Bc2, scale_by="none", Q2=nQn)

  for (s in c("test", "retest")) {
    r <- diag(cor(t(dr[[s]]$S), tmp))  # row q of S <-> template column q
    expect_true(all(r > .8), info=paste(s, "r =", paste(round(r, 2), collapse=", ")))
  }
})

test_that("GSR=TRUE runs", {
  dr <- run(B1s[[4]], B2s[[4]], GSR=TRUE)
  expect_equal(dim(dr$test$S), c(nQs, nV))
})

test_that("FC_updateA_path: BOLDkeep saved with mask2 (both paths)", {
  B1 <- B1s[[4]]; B1[1:3, ] <- 100   # ensures mask2 is non-trivial
  for (Q2 in c(0, nQn)) {
    d <- file.path(tempdir(), paste0("fcua_", Q2)); dir.create(d, showWarnings=FALSE)
    run(B1, B2s[[4]], Q2=Q2, FC_updateA_path=d)
    keep <- readRDS(file.path(d, "BOLDkeep.rds"))
    expect_true(all(c("test", "retest", "mask2") %in% names(keep)))
    expect_length(keep$mask2, nV)
    expect_equal(nrow(keep$test), sum(keep$mask2))
    expect_equal(nrow(keep$retest), sum(keep$mask2))
  }
})

test_that("no-mask case also saves mask2 (all TRUE)", {
  d <- file.path(tempdir(), "fcua_nomask"); dir.create(d, showWarnings=FALSE)
  run(B1s[[4]], B2s[[4]], FC_updateA_path=d)
  keep <- readRDS(file.path(d, "BOLDkeep.rds"))
  expect_true("mask2" %in% names(keep))
  expect_true(all(keep$mask2))
})

test_that("results are reproducible (Q2=0)", {
  a <- run(B1s[[1]], B2s[[1]], nuisance=list(nuis, nuis))
  b <- run(B1s[[1]], B2s[[1]], nuisance=list(nuis, nuis))
  expect_equal(a, b)
})

# ---- Masking stress tests ----------------------------------------------------
# Steps in dual_reg2 being exercised:
#   (1) input `mask`  -> BOLD[mask,], nV <- sum(mask), template pre-masked
#   (2) drop_first    -> happens BEFORE the variance/NA check
#   (3) mask2 = mask_BOLD(BOLD) & mask_BOLD(BOLD2)  (raw data, pre-regression)
#   (4) maskTol       -> compared against nV AFTER the input mask
#   (5) unmask        -> outputs have length(sum(mask)), NA at mask2 positions

test_that("input mask + mask2: NA positions are in masked space, intersection across scans", {
  mask <- rep(TRUE, nV); mask[1:50] <- FALSE
  B1 <- B1s[[4]]; B2 <- B2s[[4]]
  B1[60, ] <- 100          # constant in scan 1 only
  B2[70, 5] <- NA          # NA in scan 2 only
  # Original voxels 60 and 70 -> positions in masked space:
  bad <- match(c(60, 70), which(mask))   # 10, 20

  dr <- run(B1, B2, mask=mask, template=tmp[mask, ], keepA=TRUE)
  for (s in c("test", "retest")) {
    expect_equal(dim(dr[[s]]$S), c(nQs, sum(mask)))      # NOT nV, NOT sum(mask)-2
    expect_length(dr[[s]]$sigma_sq, sum(mask))
    expect_true(all(is.na(dr[[s]]$S[, bad])))            # masked in BOTH scans
    expect_false(anyNA(dr[[s]]$S[, -bad]))
    expect_true(all(is.na(dr[[s]]$sigma_sq[bad])))
    expect_false(anyNA(dr[[s]]$sigma_sq[-bad]))
  }
  # Identical to running on the clean subset directly.
  keep <- setdiff(which(mask), c(60, 70))
  dr2 <- run(B1[keep, ], B2[keep, ], template=tmp[keep, ])
  expect_equal(dr$test$S[, -bad], dr2$test$S, tolerance=1e-6)
  expect_equal(dr$retest$sigma_sq[-bad], dr2$retest$sigma_sq, tolerance=1e-6)
})

test_that("maskTol proportion is relative to post-input-mask nV (boundary)", {
  mask <- seq_len(nV) <= 300            # nV becomes 300 -> .1 allows 30 masked
  ok  <- B1s[[4]]; ok[1:30, ]  <- 100   # 30 masked: 30 > 30 is FALSE -> runs
  bad <- B1s[[4]]; bad[1:31, ] <- 100   # 31 masked -> NULL
  # (31 is only 3.4% of 900, so a wrong base would let this through.)
  expect_false(is.null(run(ok,  B2s[[4]], mask=mask, template=tmp[mask, ])))
  expect_null(        run(bad, B2s[[4]], mask=mask, template=tmp[mask, ]))
})

test_that("maskTol >= 1 is a count of locations", {
  B <- B1s[[4]]; B[1:3, ] <- 100
  expect_null(run(B, B2s[[4]], maskTol=2))
  expect_false(is.null(run(B, B2s[[4]], maskTol=3)))   # 3 > 3 is FALSE
})

test_that("drop_first happens before the variance/NA check", {
  B1 <- B1s[[4]]
  B1[7, ] <- 100; B1[7, seq(drop_first)] <- rnorm(drop_first, 100, 10)  # varies only in dropped vols
  B1[8, 1] <- NA                                                        # NA only in a dropped vol
  dr <- run(B1, B2s[[4]], drop_first=drop_first)
  expect_true(all(is.na(dr$test$S[, 7])))     # constant after dropping -> masked
  expect_false(anyNA(dr$test$S[, 8]))         # NA was dropped -> kept
  expect_true(is.na(dr$test$sigma_sq[7]))
  expect_false(is.na(dr$test$sigma_sq[8]))
})

# [TO DO]: revisit these
# # ---- Likely to FAIL: mask2 is computed on raw data, before regression/scrub/split.
# # A voxel can pass the raw variance check but have ~zero variance afterward,
# # which makes SD scaling stop with "zero or negative scaling measures".
# # Desired behavior: no error (voxel masked/NA, or handled). Failing = real gap.
#
# test_that("voxel constant except at a scrubbed volume does not crash SD scaling", {
#   B1 <- B1s[[4]]; B1[11, ] <- 100; B1[11, 55] <- 300
#   expect_no_error(dr <- run(B1, B2s[[4]], scrub=list(55, 55), scale_by="sd"))
#   expect_equal(dim(dr$test$S), c(nQs, nV))
# })
#
# test_that("voxel fully explained by nuisance regressors does not crash SD scaling", {
#   B1 <- B1s[[4]]; B1[12, ] <- 100 + 5 * nuis[, 1]
#   expect_no_error(dr <- run(B1, B2s[[4]], nuisance=list(nuis, nuis), scale_by="sd"))
#   expect_equal(dim(dr$test$S), c(nQs, nV))
# })
#
# test_that("pseudo-retest: voxel constant in only one half does not crash SD scaling", {
#   B1 <- B1s[[4]]; B1[9, seq(nT/2 + 1, nT)] <- 100   # varies in half 1, constant in half 2
#   expect_no_error(dr <- run(B1, NULL, scale_by="sd"))
#   expect_equal(dim(dr$retest$S), c(nQs, nV))
# })

# ---- Group simulation ---------------------------------------------------------
set.seed(77)
nN <- 6

# Subject-specific spatial maps (shared by that subject's two scans) so that the
# between-subject variance is > 0.
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

# Per-subject nuisance/scrub (different for each subject and scan)
nu1 <- lapply(seq(nN), function(i) matrix(rnorm(nT*2), nrow=nT))
nu2 <- lapply(seq(nN), function(i) matrix(rnorm(nT*2), nrow=nT))
sc1 <- lapply(seq(nN), function(i) c(10 + i, 55))
sc2 <- lapply(seq(nN), function(i) c(20 + i, 90))

# Input mask: drop first 50 locations. Positions below are in MASKED space.
mask <- rep(TRUE, nV); mask[1:50] <- FALSE
pos <- function(v) match(v, which(mask))   # orig voxel -> masked-space index

# Quiet wrapper
ep <- function(B1, B2=NULL, ...){
  args <- list(
    BOLD=B1, BOLD2=B2, template=tmp, scale_by="none", scale_sm_FWHM=0,
    FC=FALSE, verbose=FALSE
  )
  dots <- list(...)
  args[names(dots)] <- dots
  do.call(estimate_prior, args)
}

# Expected (rescaled) prior from raw S (2 x N x L x V), at locations `v`, using subjects `subj`.
exp_prior <- function(S, sigma_sq0, subj, v){
  rs <- sqrt(sigma_sq0) / mean(sqrt(sigma_sq0), na.rm=TRUE)
  S <- S[, subj, , v, drop=FALSE]
  m <- (S[1,,,,drop=FALSE] + S[2,,,,drop=FALSE]) / 2; dim(m) <- dim(S)[-1]
  d <- (S[1,,,,drop=FALSE] - S[2,,,,drop=FALSE]);     dim(d) <- dim(S)[-1]
  mu  <- t(apply(S, c(3, 4), mean))
  vnn <- t(apply(m, c(2, 3), var))
  vub <- vnn - t(apply(d, c(2, 3), var)) / 4
  list(mean = mu / rs[v], varNN = vnn / rs[v]^2, varUB = vub / rs[v]^2)
}

# ---- estimate_prior_from_DR (unit) -------------------------------------------

test_that("estimate_prior_from_DR matches hand-computed mean/varNN/varUB (M=2)", {
  set.seed(1)
  N <- 10; L <- 3; V <- 4
  DR <- array(rnorm(2*N*L*V), c(2, N, L*V))
  DR[2, , ] <- DR[1, , ] + matrix(rnorm(N*L*V, sd=.3), nrow=N)
  res <- BayesBrainMap:::estimate_prior_from_DR(DR, LV=c(L, V))
  m <- (DR[1, , ] + DR[2, , ]) / 2
  d <- DR[1, , ] - DR[2, , ]
  vnn <- apply(m, 2, var)
  expect_equal(dim(res$prior$mean), c(L, V))
  expect_equal(c(res$prior$mean),  colMeans(m), tolerance=1e-8)
  expect_equal(c(res$prior$varNN), vnn, tolerance=1e-8)
  # If this fails, my two-way-ANOVA formula may differ from var_decomp's; compare
  # against estimate_prior_from_DR_two.
  expect_equal(c(res$prior$varUB), vnn - apply(d, 2, var)/4, tolerance=1e-8)
  expect_error(BayesBrainMap:::estimate_prior_from_DR(DR, LV=c(5, 5)))
})

# ---- Basic structure ---------------------------------------------------------

test_that("retest, matrix input: structure and dimensions", {
  res <- ep(grp$B1, grp$B2)
  expect_s3_class(res, "prior.matrix")
  expect_true(all(c("prior", "var_decomp", "sigma_sq0", "mask_input", "mask",
                    "params", "dat_struct", "template_parc_table") %in% names(res)))
  expect_named(res$prior, c("mean", "varUB", "varNN"))
  for (x in res$prior) expect_equal(dim(x), c(nV, nQs))
  for (x in res$var_decomp) expect_equal(dim(x), c(nV, nQs))
  expect_length(res$sigma_sq0, nV)
  expect_null(res$mask)         # nothing masked
  expect_null(res$mask_input)
  expect_equal(res$params$num_subjects, nN)
  expect_false(res$params$pseudo_retest)
  expect_true(all(is.finite(res$prior$mean)))
  expect_true(all(res$prior$varNN >= 0))
})

test_that("pseudo-retest (BOLD2 = NULL) runs", {
  res <- ep(grp$B1, NULL)
  expect_true(res$params$pseudo_retest)
  expect_equal(dim(res$prior$mean), c(nV, nQs))
})

test_that("prior aggregation matches hand computation (keep_S)", {
  res <- ep(grp$B1, grp$B2, keep_S=TRUE)
  expect_equal(dim(res$S), c(2, nN, nQs, nV))
  e <- exp_prior(res$S, res$sigma_sq0, subj=seq(nN), v=seq(nV))
  expect_equal(res$prior$mean,  e$mean,  tolerance=1e-8)
  expect_equal(res$prior$varNN, e$varNN, tolerance=1e-8)
  expect_equal(res$prior$varUB, e$varUB, tolerance=1e-8)
})

test_that("keep_S as a file path saves rather than returns", {
  f <- tempfile(fileext=".rds")
  res <- ep(grp$B1, grp$B2, keep_S=f)
  expect_null(res$S)
  expect_true(file.exists(f))
  expect_equal(dim(readRDS(f)), c(2, nN, nQs, nV))
  expect_error(ep(grp$B1, grp$B2, keep_S=list(1)))
})

# ---- Argument handling -------------------------------------------------------

test_that("scale smoothing is switched off (with warning) for matrix data", {
  expect_warning(res <- ep(grp$B1, grp$B2, scale_sm_FWHM=4), "Scale smoothing")
  expect_equal(res$params$scale_sm_FWHM, 0)
})

test_that("hpf without TR: error, or message for hpf=.01", {
  expect_error(ep(grp$B1, grp$B2, hpf=.02), "TR")
  expect_message(res <- ep(grp$B1, grp$B2, hpf=.01), "hpf=0")
  expect_equal(res$params$hpf, 0)
})

test_that("`inds` subsets networks AFTER dual regression", {
  full <- ep(grp$B1, grp$B2)
  sub  <- ep(grp$B1, grp$B2, inds=c(2, 5))
  expect_equal(sub$params$inds, c(2, 5))
  expect_equal(sub$prior$mean,  full$prior$mean[, c(2, 5)])
  expect_equal(sub$prior$varNN, full$prior$varNN[, c(2, 5)])
  expect_equal(sub$sigma_sq0, full$sigma_sq0)
  expect_error(ep(grp$B1, grp$B2, inds=99), "inds")
})

test_that("Q2 denoising works through estimate_prior and lowers residual variance", {
  r0 <- ep(grp$B1, grp$B2, Q2=0)
  r3 <- ep(grp$B1, grp$B2, Q2=nQn)
  expect_lt(mean(r3$sigma_sq0), mean(r0$sigma_sq0))
})

test_that("per-subject nuisance/scrub are mapped to the right subject (retest)", {
  res <- ep(grp$B1, grp$B2, drop_first=drop_first, TR=1, hpf=.01,
            nuisance=list(nu1, nu2), scrub=list(sc1, sc2), keep_S=TRUE)
  for (i in c(1, 3, nN)) {   # includes first and last subject
    d <- run(grp$B1[[i]], grp$B2[[i]], drop_first=drop_first, hpf=.01, scale_by="none",
             nuisance=list(nu1[[i]], nu2[[i]]), scrub=list(sc1[[i]], sc2[[i]]))
    expect_equal(res$S[1, i, , ], d$test$S,   tolerance=1e-8, ignore_attr=TRUE)
    expect_equal(res$S[2, i, , ], d$retest$S, tolerance=1e-8, ignore_attr=TRUE)
  }
})

test_that("per-subject nuisance/scrub (pseudo-retest)", {
  res <- ep(grp$B1, NULL, drop_first=drop_first, TR=1, hpf=.01,
            nuisance=nu1, scrub=sc1, keep_S=TRUE)
  for (i in c(1, nN)) {
    d <- run(grp$B1[[i]], NULL, drop_first=drop_first, hpf=.01, scale_by="none",
             nuisance=nu1[[i]], scrub=sc1[[i]])
    expect_equal(res$S[1, i, , ], d$test$S, tolerance=1e-8, ignore_attr=TRUE)
  }
  expect_error(ep(grp$B1, NULL, nuisance=nu1[-1]))                 # wrong length
  expect_error(ep(grp$B1, grp$B2, nuisance=nu1))                   # retest needs list of 2
})

test_that("missing retest file: pair is excluded, nuisance/scrub stay aligned", {
  d <- tempfile(); dir.create(d)
  f1 <- file.path(d, sprintf("b1_%d.rds", seq(nN)))
  f2 <- file.path(d, sprintf("b2_%d.rds", seq(nN)))
  for (i in seq(nN)) {
    saveRDS(grp$B1[[i]], f1[i])
    if (i != 2) saveRDS(grp$B2[[i]], f2[i])   # subject 2 has no retest file
  }
  args <- list(drop_first=drop_first, TR=1, hpf=.01, keep_S=TRUE)
  expect_warning(
    res <- do.call(ep, c(list(f1, f2, nuisance=list(nu1, nu2), scrub=list(sc1, sc2)), args)),
    "excluded"
  )
  ref <- do.call(ep, c(list(grp$B1[-2], grp$B2[-2],
                            nuisance=list(nu1[-2], nu2[-2]), scrub=list(sc1[-2], sc2[-2])), args))
  expect_equal(res$params$num_subjects, nN - 1)
  expect_equal(res$S, ref$S, tolerance=1e-8)
  expect_equal(res$prior$mean, ref$prior$mean, tolerance=1e-8)
})

# ---- Masking -----------------------------------------------------------------

test_that("input `mask` == pre-subsetting data and template", {
  res <- ep(grp$B1, grp$B2, mask=mask)
  ref <- ep(lapply(grp$B1, function(x) x[mask, ]),
            lapply(grp$B2, function(x) x[mask, ]), template=tmp[mask, ])
  expect_equal(nrow(res$prior$mean), sum(mask))
  expect_equal(res$mask_input, mask)
  expect_equal(res$prior$mean,  ref$prior$mean,  tolerance=1e-8)
  expect_equal(res$prior$varUB, ref$prior$varUB, tolerance=1e-8)
  expect_equal(res$sigma_sq0,   ref$sigma_sq0,   tolerance=1e-8)
})

test_that("numeric 0/1 mask is coerced; wrong-length mask errors", {
  expect_message(res <- ep(grp$B1, grp$B2, mask=as.numeric(mask)), "Coercing")
  expect_equal(res$mask_input, mask)
  expect_error(ep(grp$B1, grp$B2, mask=mask[-1]))
})

# Bad locations, in masked space:  pos(60)=10, pos(70)=20, pos(80)=30
#   subj 2: constant voxel 60 in TEST scan only     -> subject-level NA (both scans)
#   subj 3: NA at voxel 70 in RETEST scan only      -> subject-level NA (both scans)
#   subj 1 and 2: bad at voxel 80                   -> 2 subjects missing
#   subj 4: constant in the INPUT-MASKED-OUT region -> must have no effect
bad_data <- function(){
  B1 <- grp$B1; B2 <- grp$B2
  B1[[2]][60, ] <- 100
  B2[[3]][70, 5] <- NA
  B1[[1]][80, ] <- 100; B2[[2]][80, ] <- 100
  B1[[4]][1:50, ] <- 100
  list(B1=B1, B2=B2)
}

test_that("missingTol: location masked when >= tol subjects lack data (masked-space indexing)", {
  bd <- bad_data()

  # Default missingTol=.1 -> 0.6 subjects: any single bad subject masks the location.
  res0 <- ep(bd$B1, bd$B2, mask=mask, keep_S=TRUE)
  expect_equal(which(!res0$mask), sort(pos(c(60, 70, 80))))
  expect_true(all(is.na(res0$prior$mean[!res0$mask, ])))
  expect_false(anyNA(res0$prior$mean[res0$mask, ]))
  expect_equal(which(is.na(res0$sigma_sq0)), sort(pos(c(60, 70, 80))))

  # missingTol=2: locations missing in 1 subject are kept; 2 subjects -> masked.
  res2 <- ep(bd$B1, bd$B2, mask=mask, keep_S=TRUE, missingTol=2)
  expect_equal(which(!res2$mask), pos(80))
  S <- res2$S
  expect_equal(dim(S), c(2, nN, nQs, sum(mask)))
  expect_true(all(is.na(S[, 2, , pos(60)])))          # subject-level NA is for BOTH visits
  expect_true(all(is.na(S[, 3, , pos(70)])))          # (retest-only NA masks the test scan too)
  expect_false(anyNA(S[, c(1, 3:nN), , pos(60)]))
  expect_false(anyNA(S[, 4, , -pos(80)]))             # input-masked-out region had no effect
  expect_true(all(is.na(res2$prior$mean[pos(80), ])))
  expect_false(is.na(res2$sigma_sq0[pos(60)]))        # averaged over the subjects that have data

  # missingTol=3: nothing masked at the group level.
  res3 <- ep(bd$B1, bd$B2, mask=mask, missingTol=3)
  expect_null(res3$mask)
  expect_false(anyNA(res3$prior$mean))
})

test_that("prior mean at a location with a missing subject uses only available subjects", {
  bd <- bad_data()
  res <- ep(bd$B1, bd$B2, mask=mask, keep_S=TRUE, missingTol=2)
  e <- exp_prior(res$S, res$sigma_sq0, subj=c(1, 3:nN), v=pos(60))
  expect_equal(res$prior$mean[pos(60), ], c(e$mean), tolerance=1e-8)
})

# LIKELY TO FAIL if estimate_prior_from_DR uses nN = dim(DR)[2] (includes the NA subject)
# for the divisor, instead of the number of subjects with data at that location.
test_that("prior VARIANCE at a location with a missing subject uses n_available - 1", {
  bd <- bad_data()
  res <- ep(bd$B1, bd$B2, mask=mask, keep_S=TRUE, missingTol=2)
  e <- exp_prior(res$S, res$sigma_sq0, subj=c(1, 3:nN), v=pos(60))
  expect_equal(res$prior$varNN[pos(60), ], c(e$varNN), tolerance=1e-8)
  expect_equal(res$prior$varUB[pos(60), ], c(e$varUB), tolerance=1e-8)
})

test_that("a subject skipped by maskTol: default missingTol errors, lenient missingTol works", {
  B1 <- grp$B1
  B1[[3]][1:200, ] <- 100      # 200/900 > maskTol=.1 -> dual_reg2 returns NULL for subject 3
  expect_error(ep(B1, grp$B2), "No locations")
  res <- ep(B1, grp$B2, missingTol=.2, keep_S=TRUE)     # .2*6 = 1.2 subjects allowed
  expect_true(all(is.na(res$S[, 3, , ])))
  expect_false(anyNA(res$S[, -3, , ]))
  expect_false(anyNA(res$prior$mean))
})

# Dropping the skipped subject up front should give the SAME prior.
# Mean should pass; variance LIKELY FAILS if the n-1 divisor counts the skipped subject.
test_that("skipped subject == subject removed from the input", {
  B1 <- grp$B1
  B1[[3]][1:200, ] <- 100
  a <- ep(B1, grp$B2, missingTol=.2)
  b <- ep(B1[-3], grp$B2[-3], missingTol=.2)
  expect_equal(a$sigma_sq0,   b$sigma_sq0,   tolerance=1e-8)
  expect_equal(a$prior$mean,  b$prior$mean,  tolerance=1e-8)
  expect_equal(a$prior$varNN, b$prior$varNN, tolerance=1e-8)
  expect_equal(a$prior$varUB, b$prior$varUB, tolerance=1e-8)
})

test_that("error on the first subject stops; skipped-for-masking first subject does not", {
  bad <- grp$B1; bad[[1]] <- bad[[1]][-1, ]      # wrong number of locations
  expect_error(ep(bad, grp$B2), "first subject")
  B1 <- grp$B1; B1[[1]][1:200, ] <- 100
  expect_no_error(ep(B1, grp$B2, missingTol=.2))
})

# ---- FC ----------------------------------------------------------------------

test_that("FC prior: shapes and validity", {
  res <- ep(grp$B1, grp$B2, FC=TRUE, FC_nPivots=4, FC_nSamp=40, keep_FC=TRUE)
  expect_equal(dim(res$FC), c(2, nN, nQs, nQs))
  for (m in 1:2) for (i in seq(nN)) {
    expect_equal(diag(res$FC[m, i, , ]), rep(1, nQs), tolerance=1e-6)
    expect_equal(res$FC[m, i, , ], t(res$FC[m, i, , ]), tolerance=1e-8)
  }
  fc <- res$prior$FC
  expect_named(fc, c("empirical", "IW", "Chol"))
  expect_equal(dim(fc$empirical$mean), c(nQs, nQs))
  expect_equal(dim(fc$Chol$mean), c(nQs, nQs))
  expect_equal(diag(fc$Chol$mean), rep(1, nQs), tolerance=1e-6)
  expect_length(fc$Chol$pivots, 4)
  expect_true(all(is.na(diag(fc$IW$mean))))
  expect_error(ep(grp$B1, grp$B2, FC=TRUE, FC_nPivots=3, FC_nSamp=40), "multiple")
})

test_that("FC_nPivots=0 skips the Cholesky prior", {
  res <- ep(grp$B1, grp$B2, FC=TRUE, FC_nPivots=0)
  expect_null(res$prior$FC$Chol)
  expect_false(is.null(res$prior$FC$IW))
})

# Trips: FC_updateA + `inds` subset + denoising + a subject with a bad location
# that is masked at the group level (prior$mean has NA rows) + input mask.
test_that("FC_updateA with inds subset, Q2, input mask, and masked locations", {
  B1 <- grp$B1; B1[[2]][60, ] <- 100
  inds <- c(1, 3, 5)
  res <- ep(B1, grp$B2, mask=mask, inds=inds, Q2=nQn,
            FC=TRUE, FC_updateA=TRUE, FC_nPivots=0, keep_FC=TRUE)
  expect_equal(which(!res$mask), pos(60))
  expect_equal(dim(res$prior$mean), c(sum(mask), length(inds)))
  expect_equal(dim(res$FC), c(2, nN, length(inds), length(inds)))
  expect_false(anyNA(res$FC))
  for (m in 1:2) for (i in seq(nN)) {
    expect_equal(diag(res$FC[m, i, , ]), rep(1, length(inds)), tolerance=1e-6)
  }
  # IW mean should have the same sign as the empirical mean
  # (fails if the IW denominator nu - nQ - 1 goes negative when nQ > length(inds)).
  emp <- res$prior$FC$empirical$mean
  iw  <- res$prior$FC$IW$mean
  expect_true(all(sign(iw) == sign(emp), na.rm=TRUE))
})

make_dr <- function(M = 2, N = 12, V = 6, seed = 1) {
  set.seed(seed)
  subj  <- matrix(rnorm(N * V, sd = 2), N, V)
  visit <- matrix(rnorm(M * V, sd = 0.5), M, V)
  x <- array(NA_real_, c(M, N, V))
  for (m in seq_len(M))
    x[m, , ] <- 10 + subj + rep(visit[m, ], each = N) + matrix(rnorm(N * V), N, V)
  x
}

test_that("prior: per-location df, NA where nS < 2 (catches the vd$nN bug)", {
  x <- make_dr(N = 8, V = 6)
  x[, , 1] <- NA; x[, -2, 2] <- NA; x[2, c(1, 4), 3] <- NA
  res <- BayesBrainMap:::estimate_prior_from_DR(x)
  expect_length(res$prior$varUB, 6)
  expect_true(all(is.na(res$prior$varUB[1:2])))
  expect_false(anyNA(res$prior$varUB[3:6]))
  for (v in 3:6) {   # == running the location alone on its complete cases
    keep <- colSums(is.na(matrix(x[, , v], nrow = 2))) == 0
    one <- BayesBrainMap:::estimate_prior_from_DR(x[, keep, v, drop = FALSE])
    for (nm in c("mean", "varUB", "varNN"))
      expect_equal(unname(res$prior[[nm]][v]), unname(one$prior[[nm]]), info = nm)
  }
})

test_that("prior agrees with mean_squares for M = 2, 3", {
  for (M in c(2, 3)) {
    x <- make_dr(M = M, N = 10, V = 4); x[c(3, 50, 77)] <- NA
    ms <- fMRItools::mean_squares(var_decomp(x)); res <- BayesBrainMap:::estimate_prior_from_DR(x)
    expect_equal(unname(res$prior$varNN), unname(ms$MSB / M))
    expect_equal(unname(res$prior$varUB), unname((ms$MSB - ms$MSR) / M))
  }
})

test_that("LV reshaping and cleanup", {
  x <- make_dr(V = 6)
  res <- BayesBrainMap:::estimate_prior_from_DR(x, LV = c(2, 3))
  expect_equal(dim(res$prior$varUB), c(2, 3))
  expect_equal(dim(res$var_decomp$SSR), c(2, 3))
  expect_equal(dim(res$var_decomp$nS), c(2, 3))
  expect_null(res$var_decomp$nM)          # catches the dropped cleanup line
  expect_error(BayesBrainMap:::estimate_prior_from_DR(x, LV = c(2, 2)))
})
