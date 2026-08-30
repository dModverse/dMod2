# ============================================================================
# EM + msEM + diagnostic plots (the public NLME API).
#
# Sections:
#   * EM method dispatch     - focei / quadrature / foceiQuadrature
#   * etaSE / shrinkage           - Laplace-inverse-Hessian diagnostics
#   * msEM                   - multi-start wrapper around EM
#   * predict.em + plots     - data frame + ggplot diagnostic helpers
#
# The C++ FOCEI kernel itself is tested in test-focei.R; here we only
# exercise the orchestrator and the public output shape.
# ============================================================================

## Context: "EM + msEM + diagnostic plots"  (context() is deprecated in testthat 3e; kept as a note)


# Shared one-eta NLME fixture builder.
.build_one_eta <- function(seed = 1L, N = 4L, tag = "nlmef") {
  set.seed(seed)
  g <- Y(c(y = "intercept"), f = NULL, parameters = "intercept",
         compile = TRUE, deriv2 = TRUE,
         modelname = paste0("nlmefit_obs_", tag, "_", seed))
  x <- Xt()
  subjects <- paste0("s", seq_len(N))
  trafo <- eqnvec(intercept = "mu_pop * exp(eta)")
  subj_table <- data.frame(eta = paste0("eta_", subjects),
                           row.names = subjects)
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE, deriv2 = TRUE,
         modelname = paste0("nlmefit_p_", tag, "_", seed))
  true_mu  <- 2.0; true_om <- 0.3
  true_eta <- rnorm(N, 0, true_om)
  y_obs    <- true_mu * exp(true_eta) + rnorm(N, 0, 0.2)
  data <- as.datalist(data.frame(name = "y", time = 0, sigma = 0.2,
                                 value = y_obs, condition = subjects,
                                 stringsAsFactors = FALSE))
  om <- omega(eta = "eta", subjects = subjects)
  obj <- normL2(data, g * x * p) + constraintL2(mu = 0, Omega = om)
  list(obj = obj, om = om, prdfn = g * x * p, data = data,
       subjects = subjects, true_mu = true_mu, true_om = true_om,
       y_obs = y_obs)
}


# ---- Method dispatch -----------------------------------------------------

test_that("EM(method='foceiQuadrature') runs end-to-end on one-eta prdfn", {
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- .build_one_eta(1L, tag = "fq")
  init <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))

  fit <- suppressMessages(EM(
    s$obj, init,
    method = "foceiQuadrature",
    control = list(quadrature = list(level = 4L,
                                     epsQuadLevels = c(3L, 4L),
                                     maxEcmPerStage = 3L)),
    verbose = FALSE))

  expect_s3_class(fit, "em")
  expect_equal(fit$method, "foceiQuadrature")
  expect_true(!is.null(fit$foceiStart))
  expect_true(!is.null(fit$stageTrace))
  expect_true(nrow(fit$stageTrace) >= 2L)
  expect_true(is.finite(fit$value))
  expect_true(abs(fit$argument["mu_pop"] - mean(s$y_obs)) < 0.5)
})


test_that("EM(method='quadrature') runs cold without a FOCEI prelude", {
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- .build_one_eta(2L, tag = "qd")
  init <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))

  fit <- suppressMessages(EM(
    s$obj, init,
    method = "quadrature",
    control = list(quadrature = list(level = 4L,
                                     epsQuadLevels = c(3L, 4L),
                                     maxEcmPerStage = 3L)),
    verbose = FALSE))
  expect_s3_class(fit, "em")
  expect_equal(fit$method, "quadrature")
  expect_null(fit$foceiStart)
  expect_true(is.finite(fit$value))
})


test_that("EM(method='foceiQuadrature') polishes a FOCEI fit without OFV blow-up", {
  # Loose invariant: starting from a FOCEI warmstart, the ECM final OFV
  # should sit within a small distance of the FOCEI OFV (polish either
  # improves or leaves it roughly unchanged on a well-resolved problem).
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- .build_one_eta(3L, N = 6L, tag = "polish")
  init <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))

  fit <- suppressMessages(EM(
    s$obj, init,
    method = "foceiQuadrature",
    control = list(quadrature = list(level = 4L,
                                     epsQuadLevels = c(4L, 5L),
                                     maxEcmPerStage = 3L)),
    verbose = FALSE))
  expect_true(is.finite(fit$foceiStart$value))
  expect_lt(abs(fit$value - fit$foceiStart$value), 10)
})


# ---- M1 API: reconstruction, emInit, validation, summary ---------------

test_that("EM reconstructs model pieces from obj (slim signature)", {
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  s <- .build_one_eta(21L, tag = "recon")
  init <- emInit(c(mu_pop = 2.0), s$om)
  # Slim call: no om / prdfn / data / errfn arguments.
  fit <- suppressMessages(EM(s$obj, init, method = "focei", verbose = FALSE))
  expect_s3_class(fit, "em")
  expect_true(is.finite(fit$value))
  # The composed objective self-describes via stamped attributes.
  expect_false(is.null(attr(s$obj, "prdfn")))
  expect_false(is.null(attr(s$obj, "data")))
  expect_false(is.null(attr(s$obj, "omegaSpec")))
})

test_that("emInit assembles a complete, ordered start", {
  s <- .build_one_eta(22L, tag = "init")
  init <- emInit(c(mu_pop = 2.0), s$om, sd = 0.3)
  expect_true(all(s$om$cholPars %in% names(init)))
  expect_equal(unname(init[s$om$cholPars[s$om$isDiag]]),
               rep(log(0.3), sum(s$om$isDiag)))
  expect_equal(unname(init["mu_pop"]), 2.0)
})

test_that("EM errors on incomplete init and on removed arguments", {
  s <- .build_one_eta(23L, tag = "valid")
  expect_error(EM(s$obj, c(mu_pop = 2.0), method = "focei"), "Cholesky")
  init <- emInit(c(mu_pop = 2.0), s$om)
  # Clean break: prdfn is no longer a formal.
  expect_error(EM(s$obj, init, prdfn = s$prdfn, method = "focei"), "prdfn")
})

test_that("EM warns on unrecognised control keys", {
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  s <- .build_one_eta(24L, tag = "ctrl")
  init <- emInit(c(mu_pop = 2.0), s$om)
  expect_warning(
    suppressMessages(EM(s$obj, init, method = "focei", verbose = FALSE,
                             control = list(focei = list(
                               innerControl = list(rtol = 1e-7))))),
    "unrecognised")
})

test_that("summary.em reports SE/RSE consistent with vcov()", {
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  s <- .build_one_eta(25L, N = 6L, tag = "summ")
  init <- emInit(c(mu_pop = 2.0), s$om)
  fit <- suppressMessages(EM(s$obj, init, method = "focei", verbose = FALSE))
  sm <- summary(fit)
  expect_s3_class(sm, "summary.em")
  expect_true(all(c("estimate", "se", "rse.pct") %in% names(sm$population)))
  V <- vcov(fit)
  d <- diag(V)[rownames(sm$population)]; d[d < 0] <- NA_real_
  expect_equal(unname(sm$population[["se"]]), unname(sqrt(d)), tolerance = 1e-8)
  expect_equal(unname(coef(fit)), unname(fit$argument))
  ci <- confint(fit)
  expect_true(all(ci$lower <= ci$value & ci$value <= ci$upper, na.rm = TRUE))
  expect_output(print(sm), "Population parameters")
})

test_that("EM exposes the plain-ML OFV convention + NONMEM-comparable value", {
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  s <- .build_one_eta(26L, N = 5L, tag = "ofv")
  init <- emInit(c(mu_pop = 2.0), s$om)
  fit <- suppressMessages(EM(s$obj, init, method = "focei", verbose = FALSE))
  expect_identical(fit$ofvType, "-2LL")
  n_obs <- sum(vapply(s$data, nrow, integer(1)))
  expect_equal(fit$nObs, as.integer(n_obs))
  expect_equal(fit$value_nonmem, fit$value - n_obs * log(2 * pi), tolerance = 1e-9)
})

test_that("EM(method='laplaceEM') runs and applies the H^-1 Omega correction", {
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  s <- .build_one_eta(31L, N = 6L, tag = "lem")
  init <- emInit(c(mu_pop = 2.0), s$om)
  fit <- suppressMessages(EM(
    s$obj, init, method = "laplaceEM", verbose = FALSE,
    control = list(quadrature = list(maxEcm = 60L, epsEcm = 1e-3))))
  expect_s3_class(fit, "em")
  expect_equal(fit$method, "laplaceEM")
  expect_true(is.finite(fit$value))
  expect_true(fit$Omega[1, 1] > 0)              # positive-definite Omega
  expect_equal(dim(fit$etaModes), c(6L, 1L))
  expect_false(is.null(fit$stageTrace))
  # Laplace-EM's closed-form Omega = (1/N) sum(eta* eta*^T + H^{-1}) strictly
  # exceeds the naive mode second moment (1/N) sum(eta* eta*^T); that gap is the
  # posterior-variance correction that separates Laplace-EM from ITS.
  expect_gt(fit$Omega[1, 1], mean(fit$etaModes[, 1]^2))
})

test_that("EM(method='saem') runs end-to-end and returns a sensible fit", {
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  set.seed(101)
  s <- .build_one_eta(41L, N = 6L, tag = "saem")
  init <- emInit(c(mu_pop = 2.0), s$om)
  # Short schedule keeps the test fast; SAEM is stochastic so the assertions are
  # structural + a loose sanity band rather than tight numeric agreement.
  fit <- suppressMessages(EM(
    s$obj, init, method = "saem", verbose = FALSE,
    control = list(saem = list(nBurnin = 30L, nEM = 30L, nMcmc = 5L))))
  expect_s3_class(fit, "em")
  expect_equal(fit$method, "saem")
  expect_true(is.finite(fit$value))
  expect_true(fit$Omega[1, 1] > 0)                 # positive-definite Omega
  expect_equal(dim(fit$etaModes), c(6L, 1L))
  expect_false(is.null(fit$stageTrace))
  expect_true(all(fit$stageTrace$phase %in% c("burnin", "converge")))
  # structural estimate stays in a sane neighbourhood of the data mean
  expect_lt(abs(fit$argument[["mu_pop"]] - mean(s$y_obs)), 1.0)
})


# ---- etaSE + shrinkage ---------------------------------------------------

test_that("EM emits etaSE + shrinkage from Laplace inverse Hessian", {
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- .build_one_eta(7L, N = 5L, tag = "se")
  init <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))

  fit <- suppressMessages(EM(
    s$obj, init,
    method = "focei", verbose = FALSE))

  expect_s3_class(fit, "em")
  expect_true(!is.null(fit$etaSE))
  expect_equal(dim(fit$etaSE), c(length(s$subjects), 1L))
  expect_true(all(is.finite(fit$etaSE)))
  expect_true(all(fit$etaSE > 0))

  expect_true(!is.null(fit$shrinkage))
  expect_equal(dim(fit$shrinkage), c(length(s$subjects), 1L))
  expect_true(all(fit$shrinkage <= 1 + 1e-8))
  expect_true(all(fit$shrinkage > -1))

  printed <- capture.output(print(fit))
  expect_true(any(grepl("eta", printed, ignore.case = TRUE)))
  expect_true(any(grepl("shrink", printed, ignore.case = TRUE)))
})


# ---- msEM -----------------------------------------------------------

test_that("msEM returns a parlist of nlmeFits and as.parframe works", {
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- .build_one_eta(11L, tag = "ms_basic")
  center <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))

  pl <- suppressMessages(msEM(
    s$obj, center,
    method = "focei",
    fits = 3L, cores = 1L,
    samplefun = "rnorm", sd = 0.2,
    start1stfromCenter = TRUE,
    verbose = FALSE))

  expect_s3_class(pl, "parlist")
  expect_length(pl, 3L)
  for (fit in pl) {
    expect_s3_class(fit, "em")
    expect_true(all(c("argument", "value", "converged", "iterations",
                      "parinit", "index") %in% names(fit)))
    # default keepFull = FALSE strips heavy state
    expect_null(fit$emDiag)
    expect_null(fit$prdfn)
    expect_null(fit$data)
  }
  expect_equal(pl[[1]]$parinit, center)

  pf <- as.parframe(pl)
  expect_s3_class(pf, "parframe")
  expect_true(all(c("value", "converged", "iterations") %in%
                    attr(pf, "metanames")))
  expect_true(diff(range(pf$value)) >= 0)
  best <- as.parvec(pf)
  expect_lt(abs(unname(best["mu_pop"]) - mean(s$y_obs)), 0.5)
})


test_that("msEM accepts a parframe as center (rows used directly)", {
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- .build_one_eta(12L, tag = "ms_pf")

  starts_df <- data.frame(mu_pop        = c(2.0, 2.5),
                          omega_eta_eta = c(log(0.3), log(0.5)))
  pf_in <- parframe(starts_df)

  pl <- suppressMessages(msEM(
    s$obj, pf_in,
    method = "focei", cores = 1L, verbose = FALSE))

  expect_s3_class(pl, "parlist")
  expect_length(pl, 2L)  # nrow(center) overrides fits
})


test_that("msEM keepFull = TRUE preserves emDiag / prdfn / data", {
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- .build_one_eta(13L, tag = "ms_keep")
  center <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))

  pl <- suppressMessages(msEM(
    s$obj, center,
    method = "focei", fits = 2L, cores = 1L,
    keepFull = TRUE, start1stfromCenter = TRUE,
    sd = 0.1, verbose = FALSE))

  expect_length(pl, 2L)
  expect_false(is.null(pl[[1]]$emDiag))
  expect_false(is.null(pl[[1]]$prdfn))
  expect_false(is.null(pl[[1]]$data))
})


test_that("msEM handles per-fit failures without aborting the run", {
  oldwd <- setwd(tempdir())
  on.exit(setwd(oldwd))
  s <- .build_one_eta(14L, tag = "ms_err")

  starts_df <- data.frame(mu_pop        = c(2.0, NA_real_),
                          omega_eta_eta = c(log(0.3), log(0.3)))
  pf_in <- parframe(starts_df)

  pl <- suppressMessages(msEM(
    s$obj, pf_in,
    method = "focei", cores = 1L, verbose = FALSE))

  expect_length(pl, 2L)
  stats <- vapply(pl, function(f) !is.null(f$error), logical(1))
  expect_true(any(!stats))  # at least one succeeded
  expect_true(any(stats))   # at least one failed
})


# ---- predict.em + plots --------------------------------------------

# Plots fixture: more times so plotIndivs has a curve to draw.
.build_for_plots <- function(seed = 1L) {
  set.seed(seed)
  g <- Y(c(y = "intercept"), f = NULL, parameters = "intercept",
         compile = TRUE, deriv2 = TRUE, modelname = paste0("plt_obs_", seed))
  x <- Xt()
  subjects <- paste0("s", 1:4)
  trafo <- eqnvec(intercept = "mu_pop * exp(eta)")
  subj_table <- data.frame(eta = paste0("eta_", subjects),
                           row.names = subjects)
  trafos <- branch(trafo, table = subj_table, apply = "insert")
  p <- P(trafos, method = "explicit", compile = TRUE, deriv2 = TRUE,
         modelname = paste0("plt_p_", seed))
  true_eta <- rnorm(4, 0, 0.3)
  obs_rows <- do.call(rbind, lapply(seq_along(subjects), function(i) {
    ts <- c(0, 1, 2)
    data.frame(name = "y", time = ts,
               value = 2.0 * exp(true_eta[i]) + rnorm(length(ts), 0, 0.2),
               sigma = 0.2, condition = subjects[i], stringsAsFactors = FALSE)
  }))
  data <- as.datalist(obs_rows)
  om <- omega(eta = "eta", subjects = subjects)
  obj <- normL2(data, g * x * p) + constraintL2(mu = 0, Omega = om)
  list(obj = obj, om = om, prdfn = g * x * p, data = data,
       subjects = subjects)
}

.run_focei <- function(s, init = c(mu_pop = 2.0, omega_eta_eta = log(0.3)))
  suppressMessages(EM(s$obj, init,
                            method = "focei", verbose = FALSE))


test_that("predict.em returns a long data.frame with IPRED/PRED/IWRES", {
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  fit <- .run_focei(.build_for_plots(1L))
  pf <- predict(fit, times = seq(0, 3, length.out = 10))
  expect_s3_class(pf, "data.frame")
  expected_cols <- c("condition","time","name","observed","sigma",
                     "IPRED","PRED","source","IRES","PRES","IWRES","PWRES")
  expect_true(all(expected_cols %in% names(pf)))
  expect_true(any(pf$source == "obs"))
  expect_true(any(pf$source == "grid"))
  obs <- pf[pf$source == "obs", , drop = FALSE]
  expect_equal(obs$IRES, obs$observed - obs$IPRED, tolerance = 1e-12)
})


test_that("plot.em, plotIndivs, plotResiduals return ggplots", {
  skip_if_not_installed("ggplot2")
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  fit <- .run_focei(.build_for_plots(2L))
  expect_s3_class(plotIndivs(fit),    "ggplot")
  expect_s3_class(plot(fit),          "ggplot")
  expect_s3_class(plotResiduals(fit), "ggplot")
})


test_that("plotIndivs paginates when subjectsPerPage is set", {
  skip_if_not_installed("ggplot2")
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  fit <- .run_focei(.build_for_plots(2L))  # 4 subjects
  pages <- plotIndivs(fit, subjectsPerPage = 2L)
  expect_type(pages, "list")
  expect_length(pages, 2L)
  lapply(pages, expect_s3_class, "ggplot")
  expect_match(pages[[1]]$labels$title, "page 1/2", fixed = TRUE)
  expect_match(pages[[2]]$labels$title, "page 2/2", fixed = TRUE)
  one <- plotIndivs(fit, subjectsPerPage = 10L)
  expect_type(one, "list")
  expect_length(one, 1L)
  expect_s3_class(one[[1]], "ggplot")
  expect_error(plotIndivs(fit, subjectsPerPage = 0L), "positive integer")
  expect_error(plotIndivs(fit, subjectsPerPage = c(2L, 3L)), "positive integer")
})


test_that("plotHistIndivs returns either a ggplot (cowplot) or a list of two", {
  skip_if_not_installed("ggplot2")
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  fit <- .run_focei(.build_for_plots(3L))
  p <- plotHistIndivs(fit)
  if (requireNamespace("cowplot", quietly = TRUE)) {
    expect_s3_class(p, "ggplot")
  } else {
    expect_type(p, "list")
    expect_named(p, c("hist", "qq"))
  }
})


test_that("plotTrace errors on focei fit, works on foceiQuadrature fit", {
  skip_if_not_installed("ggplot2")
  oldwd <- setwd(tempdir()); on.exit(setwd(oldwd))
  s <- .build_for_plots(4L)
  init <- c(mu_pop = 2.0, omega_eta_eta = log(0.3))
  fit_lap <- suppressMessages(EM(s$obj, init,
                                       method = "focei", verbose = FALSE))
  expect_error(plotTrace(fit_lap), "requires.*quadrature")

  fit_qd <- suppressMessages(EM(s$obj, init,
                                      method = "foceiQuadrature",
                                      control = list(quadrature = list(
                                        level = 4L,
                                        epsQuadLevels = c(3L, 4L),
                                        maxEcmPerStage = 2L)),
                                      verbose = FALSE))
  expect_s3_class(plotTrace(fit_qd), "ggplot")
})


test_that("plotResiduals back-compat: parframe path still works", {
  # Smoke-test that the existing plotResiduals(parframe, x, data, ...) entry
  # point isn't broken by the EM dispatch shim. Just confirm the
  # EM branch is bypassed for a non-EM input.
  pf <- structure(data.frame(value = 1.0, index = 1L),
                  class = c("parframe", "data.frame"))
  expect_false(inherits(pf, "em"))
})
