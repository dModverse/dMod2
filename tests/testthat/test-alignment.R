# alignmentModel() and align(): free expressions, normalization, derivatives, fits.

.alignmentSim <- function(nGel = 6, seed = 1) {
  set.seed(seed)
  times <- c(0, 10, 30, 60, 120)
  d <- expand.grid(time = times, gel = sprintf("G%02d", seq_len(nGel)), stringsAsFactors = FALSE)
  truth <- c(0.02, 1, 4, 2, 0.5)[match(d$time, times)]
  s <- rlnorm(nGel, 7, 0.5)[match(d$gel, unique(d$gel))]
  b <- rlnorm(nGel, 4, 0.3)[match(d$gel, unique(d$gel))]
  d$value <- (s * truth + b) * rlnorm(nrow(d), 0, 0.1)
  d$time <- as.character(d$time)
  d
}

# every symbol on the log scale
.logs <- function(...) paste0(c(...), " <- exp(", c(...), ")")

.am <- local({
  cache <- list()
  function(key, ...) {
    if (is.null(cache[[key]])) cache[[key]] <<- alignmentModel(...)
    cache[[key]]
  }
})

.rel <- function(d, ...) alignmentModel(d, model = "yi * sj + bj", errmodel = "sigmaR * prediction",
                                        biological = yi ~ time, scaling = sj ~ gel, error = bj + sigmaR ~ 1,
                                        trafo = .logs("yi", "sj", "bj", "sigmaR"), ...)

test_that("the values are normalized per block of connected experiments and the outputs have exact derivatives", {
  d <- .alignmentSim(4)
  m <- .am("bg", d, model = "yi * sj + bj", errmodel = "sqrt(sigmaA^2 + (sigmaR * prediction)^2)",
           biological = yi ~ time, scaling = sj + bj ~ gel, error = sigmaA + sigmaR ~ 1,
           trafo = .logs("yi", "sj", "bj", "sigmaA", "sigmaR"))
  expect_s3_class(m, "alignmentmodel")
  expect_equal(colnames(m$design$sj$X), c("sj", "sj_gelG02", "sj_gelG03", "sj_gelG04"))
  # one block, the mean of the predictors of yi over its five values
  nz <- m$normalization
  expect_equal(nz$symbol, "yi")
  expect_equal(unique(nz$block), 1L)
  expect_equal(unname(nz$A[1, ]), c(1, rep(0.2, 4)))
  # experiments without a shared value form blocks of their own
  d2 <- rbind(d, transform(d, gel = paste0(gel, "b"), time = paste0(time, "b")))
  m2 <- alignmentModel(d2, model = "yi * sj", errmodel = "sigmaR * prediction", biological = yi ~ time,
                       scaling = sj ~ gel, error = sigmaR ~ 1, trafo = .logs("yi", "sj", "sigmaR"))
  expect_equal(sort(unique(m2$normalization$block)), 1:2)
  expect_equal(m2$normalization$block, rep(1:2, each = nrow(d)))
  expect_equal(dMod2:::.connectedBlocks(c("a", "b", "b", "c"), c("x", "x", "y", "z")), c(1L, 1L, 1L, 2L))
  ap <- dMod2:::.alignmentPars(m)
  set.seed(2)
  p <- setNames(runif(length(ap$fixed), -0.3, 0.3), ap$fixed)
  rows <- 1:7
  ev <- dMod2:::.alignmentEval(m, rows, p, deriv2 = TRUE)
  pn <- dimnames(ev$J)[[3]]
  f <- function(q, r, o) { pp <- p; pp[pn] <- q; dMod2:::.alignmentEval(m, rows, pp, deriv = FALSE)$y[r, o] }
  for (o in 1:2) for (r in c(1, 6)) {
    expect_equal(unname(ev$J[r, o, ]), numDeriv::grad(f, unname(p[pn]), r = r, o = o), tolerance = 1e-7)
    g <- function(q) { pp <- p; pp[pn] <- q; dMod2:::.alignmentEval(m, rows, pp)$J[r, o, ] }
    expect_equal(unname(ev$H[r, o, , ]), numDeriv::jacobian(g, unname(p[pn])), tolerance = 1e-7)
  }
})

test_that("align() holds the geometric mean of the values per block at 1", {
  d <- .alignmentSim(4)
  d2 <- rbind(d, transform(d, gel = paste0(gel, "b"), time = paste0(time, "b"), value = 3 * value))
  m <- alignmentModel(d2, model = "yi * sj", errmodel = "sigmaR * prediction", biological = yi ~ time,
                      scaling = sj ~ gel, error = sigmaR ~ 1, trafo = .logs("yi", "sj", "sigmaR"))
  fit <- align(m, fits = 1)
  a <- fit$aligned
  second <- grepl("b$", a$time)
  expect_equal(mean(a$estimate[!second]), 0, tolerance = 1e-4)
  expect_equal(mean(a$estimate[second]), 0, tolerance = 1e-4)
  expect_equal(fit$npar, length(fit$fit$argument) - 2L)
  expect_true(all(is.finite(a$se)))
})

test_that("start values and trafos of any scale and name", {
  expect_equal(dMod2:::.trafoInverse(list(trafo = c(yi = "yi")), "yi", 1000), 1000, tolerance = 1e-10)
  expect_equal(dMod2:::.trafoInverse(list(trafo = c(yi = "exp(yi)")), "yi", 1e60), log(1e60), tolerance = 1e-10)
  expect_true(is.na(dMod2:::.trafoInverse(list(trafo = c(yi = "exp(yi)")), "yi", -1)))
  # biological symbols with different formulas; a symbol with "_"; one lambda for the symbols of a formula
  d <- .alignmentSim(4)
  m <- alignmentModel(d, model = "y_i * sj + zi", errmodel = "sigmaR * prediction",
                      biological = list(y_i ~ time, zi ~ 1), scaling = sj ~ gel, error = sigmaR ~ 1,
                      trafo = .logs("y_i", "sj", "zi", "sigmaR"))
  fit <- suppressMessages(align(m, fits = 1))
  expect_equal(unique(fit$aligned$symbol), c("y_i", "zi"))
  expect_true(is.na(fit$aligned$time[fit$aligned$symbol == "zi"]))
  ml <- alignmentModel(d, model = "yi * sj + bj", errmodel = "sigmaR * prediction", biological = yi ~ time,
                       scaling = sj ~ gel, error = bj + sigmaR ~ 1 + (1 | gr(gel, dist = laplace())),
                       trafo = .logs("yi", "sj", "bj", "sigmaR"))
  expect_length(dMod2:::.alignmentRandomEffects(ml)$blocks, 1L)
})

test_that("a numeric column on the right gives levels, I() a covariate, and the formulas use lm() coding", {
  d <- .alignmentSim(3)
  d$t <- as.numeric(d$time)
  m <- alignmentModel(d, model = "yi * sj", errmodel = "sigmaR * prediction", biological = yi ~ t,
                      scaling = sj ~ gel, error = sigmaR ~ 1, trafo = .logs("yi", "sj", "sigmaR"))
  expect_equal(ncol(m$design$yi$X), 5L)
  expect_true(is.numeric(m$data$t))
  m <- alignmentModel(d, model = "yi * sj", errmodel = "sigmaR * prediction", biological = yi ~ I(t),
                      scaling = sj ~ gel, error = sigmaR ~ 1, trafo = .logs("yi", "sj", "sigmaR"))
  expect_equal(ncol(m$design$yi$X), 2L)
  m <- alignmentModel(d, model = "yi * sj", errmodel = "sigmaR * prediction", biological = yi ~ time,
                      scaling = sj ~ gel, error = sigmaR ~ 1, trafo = .logs("yi", "sj", "sigmaR"))
  expect_equal(ncol(m$design$yi$X), 5L)
  expect_error(alignmentModel(d, model = "yi * sj", errmodel = "s", biological = yi ~ time,
                              scaling = sj ~ gel, error = s + u ~ 1), "not in the expressions")
  expect_error(alignmentModel(d, model = "yi * sj", errmodel = "s", biological = yi ~ time,
                              scaling = sj ~ gel, error = s ~ 1, trafo = "zz <- exp(zz)"), "not a symbol")
})

test_that("align() fits the values and a finite error for each value", {
  d <- .alignmentSim(6)
  m <- .am("rel", d, model = "yi * sj + bj", errmodel = "sigmaR * prediction", biological = yi ~ time,
           scaling = sj ~ gel, error = bj + sigmaR ~ 1, trafo = .logs("yi", "sj", "bj", "sigmaR"))
  fit <- align(m, fits = 3)
  expect_s3_class(fit, "alignment")
  expect_equal(nrow(fit$aligned), 5L)
  expect_true(all(is.finite(fit$aligned$se)))
  rel <- fit$aligned$value / fit$aligned$value[fit$aligned$time == "10"]
  expect_equal(rel[fit$aligned$time == "30"], 4, tolerance = 0.2)
  expect_equal(rel[fit$aligned$time == "60"], 2, tolerance = 0.2)
  expect_equal(fit$aligned$sigma, fit$aligned$value * fit$aligned$se)
  expect_equal(fit$scaled$residual, (d$value - fit$scaled$prediction) / fit$scaled$sigma)
  # the scaled reading is the value of the biological symbol at which the prediction meets the reading
  lin <- dMod2:::.alignmentPredictors(fit$model, seq_len(nrow(d)), fit$parameters)
  sj <- exp(lin[, "sj"]); bj <- exp(lin[, "bj"])
  ex <- unname((d$value - bj) / sj)
  expect_equal(fit$scaled$scaled, ifelse(ex > 0, ex, NaN), tolerance = 1e-8)
  ok <- ex > 0
  expect_equal(fit$scaled$scaledSigma[ok], unname(fit$scaled$sigma / sj)[ok], tolerance = 1e-8)
  dl <- as.datalist(fit)
  expect_equal(dl[[1]]$value, fit$aligned$value)
  expect_equal(unique(dl[[1]]$name), "yi")
  expect_s3_class(plotAligned(fit), "ggplot")
  expect_equal(nrow(attr(plotValues(fit), "data")), nrow(fit$multistart))
  expect_s3_class(plotScaled(fit), "ggplot")
  expect_s3_class(plotResiduals(fit), "ggplot")
  # the fit is the trust optimum: zero gradient of the objective
  expect_lt(max(abs(fit$objfn(fit$fit$argument)$gradient)), 1e-2)
})

test_that("an error model without symbols takes error = NULL, and init sets starts", {
  d <- .alignmentSim(4)
  expect_error(alignmentModel(d, model = "yi * sj", errmodel = "0.1 * prediction", biological = yi ~ time,
                              scaling = sj ~ gel), "error = NULL")
  expect_error(alignmentModel(d, model = "yi * sj", errmodel = "0.1 * prediction", biological = yi ~ time,
                              scaling = NULL, error = NULL), "need formulas")
  m <- alignmentModel(d, model = "yi * sj", errmodel = "0.1 * prediction", biological = yi ~ time,
                      scaling = sj ~ gel, error = NULL, trafo = .logs("yi", "sj"))
  expect_false("error" %in% m$symbols$role)
  fit <- suppressMessages(align(m, fits = 1, init = c(yi = 1)))
  expect_true(all(is.finite(fit$aligned$se)))
  expect_equal(unname(fit$scaled$sigma), unname(0.1 * fit$scaled$prediction))
  expect_error(suppressMessages(align(m, fits = 1, init = c(yi = NA))), "give `init`")
  # a standard deviation from a data column: the reverse sweep evaluates the prediction alone
  d$sdev <- 0.1 * d$value
  m0 <- alignmentModel(d, model = "yi * sj", errmodel = "sdev", biological = yi ~ time,
                       scaling = sj ~ gel, error = NULL, trafo = .logs("yi", "sj"))
  fit0 <- suppressMessages(align(m0, fits = 1, init = c(yi = 1)))
  p <- fit0$fit$argument + 0.05
  expect_equal(fit0$objfn(p, sweep = "reverse")$gradient[names(p)], fit0$objfn(p)$gradient[names(p)])
  shared <- environment(dMod2:::.fnLeaves(attr(fit0$objfn, "prdfn"))[[1]]$kernel)$shared
  expect_false(shared$sigmaVjp)
  expect_length(shared$vjpEval, 0L)
  expect_null(shared$vjpEvalMany)
})

test_that("without random terms every experiment is a condition and trust() runs sparse", {
  d <- .alignmentSim(5)
  m <- .rel(d)
  expect_setequal(unique(dMod2:::.alignmentConditions(m)), unique(d$gel))
  # the scale of every experiment starts at the median of the scales that meet its readings
  st <- dMod2:::.alignmentStart(m, dMod2:::.alignmentPars(m)$fixed)
  lin <- dMod2:::.alignmentPredictors(m, seq_len(nrow(d)), st)
  expect_equal(tapply(lin[, "sj"], d$gel, function(x) x[1]),
               tapply(log((d$value - 1) / exp(st[["yi"]])), d$gel, median), tolerance = 1e-8)
  fit <- withr::with_options(list(dMod2.sparseAuto = list(parameters = 0, fill = 1)),
                              suppressMessages(align(m, fits = 1)))
  expect_s4_class(fit$fit$hessian, "sparseMatrix")
  dense <- suppressMessages(align(m, fits = 1))
  expect_true(is.matrix(dense$fit$hessian))
  expect_equal(dense$value, fit$value, tolerance = 1e-6)
  expect_equal(dense$aligned$se, fit$aligned$se, tolerance = 1e-4)
  # all conditions in one batch as the conditions one by one, forward and backward
  o <- dMod2:::.alignmentObjective(m)
  p <- fit$fit$argument + 0.05
  all <- o(p, deriv2 = TRUE)
  one <- Reduce(`+`, lapply(unique(dMod2:::.alignmentConditions(m)), function(cn) o(p, deriv2 = TRUE, conditions = cn)))
  nm <- names(p)
  expect_equal(all$value, one$value, tolerance = 1e-12)
  expect_equal(all$gradient[nm], one$gradient[nm], tolerance = 1e-12)
  expect_equal(all$hessian[nm, nm], one$hessian[nm, nm], tolerance = 1e-12)
  expect_equal(o(p, sweep = "reverse")$gradient[nm], all$gradient[nm], tolerance = 1e-10)
})

test_that("a prior enters the objective, and the values of a symbol without trafo have the mean 0", {
  d <- .alignmentSim(4)
  d$value <- log2(d$value)
  add <- function(...) alignmentModel(d, model = "yi + sj", errmodel = "sigma", biological = yi ~ time,
                                      error = sigma ~ 1, trafo = "sigma <- exp(sigma)", ...)
  m <- add(scaling = sj ~ gel, prior = constraintL2(mu = c(sigma = log(0.5)), sigma = 1e-4))
  fit <- align(m, fits = 1)
  expect_equal(unname(fit$parameters["sigma"]), log(0.5), tolerance = 1e-4)
  expect_equal(mean(fit$aligned$value), 0, tolerance = 1e-4)
  expect_s3_class(add(scaling = sj ~ 1 + (1 | gel), prior = constraintL2(mu = c(omega_sj_sj = 0))),
                  "alignmentmodel")
  expect_error(add(scaling = sj ~ gel, prior = constraintL2(mu = c(omega_sj_sj = 0))), "no parameter")
  expect_error(add(scaling = sj ~ gel, prior = c(sigma = 0)), "constraintL2")
})

test_that("a random background per gel fits through EM, censored readings by M3", {
  d <- .alignmentSim(8)
  d$lloq <- NA_real_
  low <- d$time == "0" & d$gel %in% c("G01", "G02")
  d$lloq[low] <- d$value[low] * 1.2
  m <- alignmentModel(d, model = "yi * sj + bj", errmodel = "sigmaR * prediction", biological = yi ~ time,
                      scaling = sj ~ gel, error = bj + sigmaR ~ 1 + (1 | gel),
                      trafo = .logs("yi", "sj", "bj", "sigmaR"))
  expect_equal(unique(dMod2:::.alignmentConditions(m)), unique(d$gel))
  fit <- suppressMessages(align(m, fits = 2))
  expect_true(is.finite(fit$value))
  expect_named(fit$multistart, c("value", "converged", "iterations"))
  expect_equal(nrow(fit$multistart), 2L)
  expect_true(all(c("eta_bj_G01", "eta_sigmaR_G01") %in% names(fit$parameters)))
  expect_equal(sort(unique(fit$omega$eta1)), c("bj", "sigmaR"))
  expect_equal(fit$scaled$bloq, low)
  expect_equal(fit$aligned$lloq[fit$aligned$time == "0"], min(fit$scaled$scaled[low]))
  expect_true(all(is.na(fit$aligned$lloq[fit$aligned$time != "0"])))
  # FOCEI on the sparse Hessian: the same marginal, gradient and Hessian, and standard errors
  # beside the optimum, where the gradient is not of the size of the convergence of the modes
  ctl <- list(trustControl = list(iterlim = 0), innerControl = list(ftol = 1e-12, mtol = 1e-12, iterlim = 200))
  point <- fit$fit$argument + 0.05
  at <- function(sparse) em(fit$objfn, fit$randomEffects, point, control = ctl, verbose = FALSE,
                            sparse = sparse)
  dense <- at(FALSE)
  sp <- at(TRUE)
  nm <- names(dense$argument)
  # forward over reverse through the vjp of the alignment kernels, against forward over forward
  expect_identical(dense$derivMode, "forward-reverse")
  ff <- em(fit$objfn, fit$randomEffects, point, control = c(ctl, list(derivMode = "forward-forward")),
           verbose = FALSE)
  expect_equal(dense$value, ff$value, tolerance = 1e-10)
  expect_equal(dense$gradient, ff$gradient, tolerance = 1e-6)
  expect_equal(dense$hessian, ff$hessian, tolerance = 1e-6)
  # the evaluation the two outputs of a sweep share is gone after it
  shared <- environment(dMod2:::.fnLeaves(attr(fit$objfn, "prdfn"))[[1]]$kernel)$shared
  # requests of the groups with etas of their own: one batch as the conditions one by one
  conds <- names(shared$rowsOf)[1:3]
  pl <- setNames(lapply(seq_along(conds), function(j) {
    lp <- shared$local[[conds[j]]]$pars
    p <- point[intersect(names(point), lp)]
    e <- setdiff(lp, names(point))
    c(p, setNames(0.1 * j + seq_along(e) / 10, e))
  }), conds)
  tm <- lapply(conds, function(cn) as.numeric(shared$rowsOf[[cn]]))
  batch <- dMod2:::.alignmentOutputsMany(shared, tm, pl, rep(list(NULL), 3), TRUE, FALSE, 1L, as.list(conds))
  for (j in seq_along(conds))
    expect_equal(batch[[j]], dMod2:::.alignmentOutputs(shared, shared$rowsOf[[conds[j]]], conds[j], tm[[j]],
                                                       pl[[j]], NULL, TRUE, FALSE), tolerance = 1e-12)
  expect_true(shared$sigmaVjp)
  expect_length(shared$vjpEval, 0L)
  expect_null(shared$vjpEvalMany)
  expect_equal(sp$value, dense$value, tolerance = 1e-8)
  expect_equal(sp$gradient[nm], dense$gradient[nm], tolerance = 1e-6)
  expect_equal(as.matrix(sp$hessian)[nm, nm], dense$hessian[nm, nm], tolerance = 1e-6)
  fs <- withr::with_options(list(dMod2.sparseAuto = list(parameters = 0, fill = 1)),
                            suppressMessages(align(m, fits = 1)))
  expect_s4_class(fs$fit$hessian, "sparseMatrix")
  # the standard errors of the sparse solve are those of the dense one on the same Hessian
  free <- intersect(names(fs$fit$argument), rownames(fs$fit$hessian))
  X <- m$design$yi$X[!duplicated(m$data$time), , drop = FALSE]
  expect_equal(fs$aligned$se, sqrt(dMod2:::.alignmentVariance(as.matrix(fs$fit$hessian), free)(X)))
  expect_equal(sqrt(dMod2:::.alignmentVariance(sp$hessian, nm)(X)),
               sqrt(dMod2:::.alignmentVariance(dense$hessian, nm)(X)), tolerance = 1e-5)
})

test_that("a Student-t alignment uses normT", {
  d <- .alignmentSim(5)
  d$value[3] <- d$value[3] * 5
  m <- .am("t", d, model = "yi * sj", errmodel = "sigmaR * prediction", biological = yi ~ time,
           scaling = sj ~ gel, error = sigmaR ~ 1, trafo = .logs("yi", "sj", "sigmaR"),
           distribution = "t", df = 3)
  fit <- suppressMessages(align(m, fits = 2))
  expect_true(is.finite(fit$value))
  expect_error(alignmentModel(d, model = "yi * sj", errmodel = "s * prediction", biological = yi ~ time,
                              scaling = sj ~ gel, error = s ~ 1, distribution = "t"), "df")
})

test_that("a covariate switches the background off for selected readings", {
  d <- .alignmentSim(4)
  d$wb <- rep(c(1, 0), length.out = nrow(d))
  m <- alignmentModel(d, model = "yi * sj + wb * bj", errmodel = "sigmaR * prediction",
                      biological = yi ~ time, scaling = sj ~ gel, error = bj + sigmaR ~ 1,
                      trafo = .logs("yi", "sj", "bj", "sigmaR"))
  expect_equal(m$covariates, "wb")
  fit <- suppressMessages(align(m, fits = 1))
  off <- d$wb == 0
  lin <- dMod2:::.alignmentPredictors(fit$model, seq_len(nrow(d)), fit$parameters)
  expect_equal(fit$scaled$prediction[off], unname(exp(lin[, "yi"] + lin[, "sj"]))[off], tolerance = 1e-10)
})

test_that("censored readings without random terms have the exact Hessian of M3", {
  d <- .alignmentSim(5)
  d$lloq <- NA_real_
  low <- d$time == "0" & d$gel %in% c("G01", "G02", "G03")
  d$lloq[low] <- d$value[low] * 1.3
  fit <- suppressWarnings(suppressMessages(align(.rel(d), fits = 1, iterlim = 5)))
  p <- fit$fit$argument + 0.05
  v <- function(q) fit$objfn(setNames(q, names(p)), deriv = FALSE)$value
  gr <- function(q) unname(fit$objfn(setNames(q, names(p)))$gradient[names(p)])
  out <- fit$objfn(p, deriv2 = TRUE)
  expect_equal(unname(out$gradient[names(p)]), numDeriv::grad(v, unname(p)), tolerance = 1e-6)
  expect_equal(unname(out$hessian[names(p), names(p)]), numDeriv::jacobian(gr, unname(p)), tolerance = 1e-6)
})


# ---- EM machinery through alignment models -----------------------------------

.foceiAt <- function(fit, p, iterlim = 0) {
  ctrl <- list(trustControl = list(iterlim = iterlim),
               innerControl = list(ftol = 1e-12, mtol = 1e-12, iterlim = 200))
  em(fit$objfn, fit$randomEffects, p, method = "focei", verbose = FALSE, control = ctrl)
}

.rbg <- function(d, ...) alignmentModel(d, model = "yi * sj + bj", errmodel = "sigmaR * prediction",
                                        biological = yi ~ time, scaling = list(sj ~ gel, bj ~ 1 + (1 | gel)),
                                        error = sigmaR ~ 1, trafo = .logs("yi", "sj", "bj", "sigmaR"), ...)

test_that("the outer model of FOCEI learns the curvature that leads off a saddle of Gauss-Newton", {
  m <- .rbg(.alignmentSim(8))
  obj <- dMod2:::.alignmentObjective(m) + dMod2:::.alignmentNormalization(m)
  st <- dMod2:::.alignmentStart(m, dMod2:::.alignmentPars(m)$fixed)
  fit <- em(obj, dMod2:::.alignmentRandomEffects(m), st, control = list(trustControl = list(mtol = 1e-4)),
            sparse = FALSE, verbose = FALSE)
  expect_lt(fit$value, 472)
})

test_that("the FOCEI gradient and Hessian are those of the FOCEI objective, censored data too", {
  d <- .alignmentSim(8)
  d$lloq <- NA_real_
  low <- d$time == "0" & d$gel %in% c("G01", "G02", "G03", "G04")
  d$lloq[low] <- d$value[low] * 1.2
  fit <- suppressMessages(align(.rbg(d), fits = 2))
  obj <- fit$objfn
  sel <- c("bj", "omega_bj_bj", "sigmaR")
  th <- fit$fit$argument
  p1 <- th + 0.05
  ofv <- function(q, at, nm = sel) { p <- at; p[nm] <- q; .foceiAt(fit, p)$value }
  num <- numDeriv::grad(ofv, unname(p1[sel]), at = p1, method.args = list(eps = 1e-4, d = 1e-3, r = 2))
  expect_equal(unname(.foceiAt(fit, p1)$gradient[sel]), num, tolerance = 1e-4)
  # the Gauss-Newton FOCEI Hessian leaves out third-order terms; in the Omega direction they are
  # large when sigma depends on eta, so only the structural diagonal is compared
  st <- c("bj", "sigmaR")
  H <- numDeriv::hessian(ofv, unname(th[st]), at = th, nm = st, method.args = list(eps = 1e-3, d = 1e-3))
  expect_equal(unname(diag(fit$fit$hessian[st, st])), diag(H), tolerance = 0.2)
  l2 <- normL2(attr(obj, "data"), attr(obj, "prdfn"), attr(obj, "errfn"), optBLOQ = "M1")
  expect_error(em(l2, fit$randomEffects, th, method = "focei", verbose = FALSE), "M3")
})

test_that("the FOCEI gradient with a random effect in the error model is the slope of the objective", {
  set.seed(5)
  d <- .alignmentSim(8)
  sc <- exp(rnorm(8, 0, 0.7))[match(d$gel, unique(d$gel))]
  d$value <- d$value * exp(rnorm(nrow(d), 0, 0.1 * sc))
  m <- alignmentModel(d, model = "yi * sj + bj", errmodel = "sigmaR * prediction", biological = yi ~ time,
                      scaling = sj ~ gel, error = bj + sigmaR ~ 1 + (1 | gel),
                      trafo = .logs("yi", "sj", "bj", "sigmaR"))
  fit <- suppressMessages(align(m, fits = 1))
  sel <- c("bj", "omega_bj_bj", "omega_sigmaR_sigmaR", "sigmaR")
  th <- fit$fit$argument + 0.05
  ofv <- function(q) { p <- th; p[sel] <- q; .foceiAt(fit, p)$value }
  num <- numDeriv::grad(ofv, unname(th[sel]), method.args = list(eps = 1e-4, d = 1e-3, r = 2))
  expect_equal(unname(.foceiAt(fit, th)$gradient[sel]), num, tolerance = 5e-3)
})

test_that("quadrature has the exact Hessian of its objective, Omega included", {
  d <- .alignmentSim(6)
  fo <- suppressMessages(align(.rbg(d), fits = 2))
  qu <- suppressMessages(align(.rbg(d), fits = 2, method = "quadrature"))
  expect_equal(qu$value, fo$value, tolerance = 0.05)
  # away from the optimum, where the gradient is not of the size of the differencing noise
  th <- qu$fit$argument + 0.05
  q <- dMod2:::.quadratureObjective(dMod2:::.emDesign(qu$objfn, qu$randomEffects), th, NULL, list())
  q$rebuild(th, 3L)
  st <- names(th)
  v <- function(x) q$value(setNames(x, st), deriv = FALSE)$value
  gr <- function(x) unname(q$value(setNames(x, st), deriv2 = TRUE)$gradient[st])
  o <- q$value(th, deriv2 = TRUE)
  expect_equal(unname(o$gradient[st]), numDeriv::grad(v, unname(th)), tolerance = 1e-6)
  expect_equal(unname(o$hessian[st, st]), numDeriv::jacobian(gr, unname(th)), tolerance = 1e-5)
  expect_true(all(is.finite(qu$fit$etaSE)) && all(qu$fit$etaSE > 0))
})

test_that("absolute readings set the units of their block, the other blocks stay normalized", {
  set.seed(4)
  times <- c("0", "10", "30", "60")
  truth <- c(0.5, 2, 4, 1)
  wb <- expand.grid(time = times, gel = sprintf("G%02d", 1:4), stringsAsFactors = FALSE)
  wb$value <- truth[match(wb$time, times)] * rlnorm(4, 3, 0.5)[match(wb$gel, unique(wb$gel))] *
    rlnorm(nrow(wb), 0, 0.05)
  ms <- data.frame(time = times, gel = NA, value = truth * rlnorm(4, 0, 0.05))
  other <- expand.grid(time = paste0(times, "b"), gel = c("H01", "H02"), stringsAsFactors = FALSE)
  other$value <- rep(truth, 2) * rlnorm(2, 3, 0.5)[match(other$gel, c("H01", "H02"))] * rlnorm(8, 0, 0.05)
  d <- rbind(cbind(wb, absolute = FALSE), cbind(ms, absolute = TRUE), cbind(other, absolute = FALSE))
  am <- function(...) alignmentModel(d, model = "yi * sj", errmodel = "sigmaR * prediction", biological = yi ~ time,
                                     scaling = sj ~ gel, error = sigmaR ~ 1, trafo = .logs("yi", "sj", "sigmaR"), ...)
  m <- am()
  expect_equal(nrow(m$normalization$A), 1L)
  expect_equal(sum(abs(m$coef$sj[d$absolute, ])), 0)
  expect_output(print(m), "4 absolute readings, model yi: 1 of 2 blocks on their scale")
  fit <- suppressMessages(align(m, fits = 1))
  a <- fit$aligned
  first <- a[a$time %in% times, ]
  expect_equal(first$value[match(times, first$time)], truth, tolerance = 0.1)
  expect_true(all(is.finite(first$se)))
  expect_equal(mean(a$estimate[!a$time %in% times]), 0, tolerance = 1e-4)
  expect_equal(fit$scaled$scaled[d$absolute], d$value[d$absolute], tolerance = 1e-8)
  p <- fit$fit$argument + 0.05
  expect_equal(fit$objfn(p, sweep = "reverse")$gradient[names(p)], fit$objfn(p)$gradient[names(p)])
  # the units of the values by the expression
  fit3 <- suppressMessages(align(am(absolute = "1000 * yi"), fits = 1))
  first3 <- fit3$aligned[fit3$aligned$time %in% times, ]
  expect_equal(first3$value[match(times, first3$time)], truth / 1000, tolerance = 0.1)
  expect_equal(fit3$value, fit$value, tolerance = 1e-4)
  expect_error(am(absolute = "yi * sj"), "sj scaling")
  expect_error(alignmentModel(d[!d$absolute, ], model = "yi * sj", errmodel = "sigmaR * prediction",
                              biological = yi ~ time, scaling = sj ~ gel, error = sigmaR ~ 1, absolute = "yi"),
               "column `absolute`")
  d$absolute[1] <- NA
  expect_error(am(), "TRUE or FALSE")
})

test_that("absolute readings of a model without trafo are their values", {
  set.seed(8)
  times <- c("0", "10", "30", "60")
  truth <- c(0.5, 2, 4, 1)
  wb <- expand.grid(time = times, gel = sprintf("G%02d", 1:4), stringsAsFactors = FALSE)
  wb$value <- truth[match(wb$time, times)] * runif(4, 5, 20)[match(wb$gel, unique(wb$gel))] *
    rlnorm(nrow(wb), 0, 0.05)
  ms <- data.frame(time = times, gel = NA, value = truth * rlnorm(4, 0, 0.05))
  d <- rbind(cbind(wb, absolute = FALSE), cbind(ms, absolute = TRUE))
  m <- alignmentModel(d, model = "yi * sj", errmodel = "sigmaR * prediction", biological = yi ~ time,
                      scaling = sj ~ gel, error = sigmaR ~ 1, trafo = "sigmaR <- exp(sigmaR)")
  fit <- suppressMessages(align(m, fits = 1))
  expect_equal(fit$aligned$value[match(times, fit$aligned$time)], truth, tolerance = 0.1)
})

test_that("a background model with absolute readings recovers the values in their units", {
  set.seed(9)
  times <- c(0, 10, 30, 60, 120)
  conc <- c(2, 40, 85, 60, 25)
  gels <- sprintf("G%02d", 1:8)
  wb <- expand.grid(time = times, gel = gels, stringsAsFactors = FALSE)
  bg <- setNames(rnorm(8, 1, 0.3), gels)
  shift <- setNames(rnorm(8, 3, 1), gels)
  wb$value <- log2(conc[match(wb$time, times)] + 2^bg[wb$gel]) + shift[wb$gel] + rnorm(nrow(wb), 0, 0.1)
  ms <- data.frame(time = rep(c(0, 30, 120), each = 2), gel = NA)
  ms$value <- log2(conc[match(ms$time, times)]) + rnorm(nrow(ms), 0, 0.1)
  d <- rbind(cbind(wb, method = "WB", absolute = FALSE), cbind(ms, method = "MS", absolute = TRUE))
  m <- alignmentModel(d, model = "log2(2^yi + 2^(b + bj)) + sj", errmodel = "sigma",
                      biological = yi ~ time, scaling = list(sj ~ gel, b ~ 1, bj ~ 0 + (1 | gel)),
                      error = sigma ~ method, trafo = "sigma <- exp(sigma)",
                      prior = constraintL2(mu = c(omega_bj_bj = log(0.3), b = 0), sigma = c(omega_bj_bj = 0.3, b = 2)))
  expect_true(all(is.na(m$design$bj$random[[1]]$groups[d$absolute])))
  # a condition per gel, named by it, and one of the absolute readings
  expect_setequal(unique(dMod2:::.alignmentConditions(m)), c(gels, "absolute"))
  # derivatives of absolute and relative rows
  ap <- dMod2:::.alignmentPars(m)
  set.seed(2)
  p <- setNames(runif(length(ap$fixed) + length(ap$eta), -0.3, 0.3), c(ap$fixed, ap$eta))
  rows <- c(which(d$absolute)[1:2], 1L, 7L)
  ev <- dMod2:::.alignmentEval(m, rows, p, deriv2 = TRUE)
  pn <- dimnames(ev$J)[[3]]
  f <- function(q, r, o) { pp <- p; pp[pn] <- q; dMod2:::.alignmentEval(m, rows, pp, deriv = FALSE)$y[r, o] }
  for (o in 1:2) for (r in seq_along(rows)) {
    expect_equal(unname(ev$J[r, o, ]), numDeriv::grad(f, unname(p[pn]), r = r, o = o), tolerance = 1e-7)
    g <- function(q) { pp <- p; pp[pn] <- q; dMod2:::.alignmentEval(m, rows, pp)$J[r, o, ] }
    expect_equal(unname(ev$H[r, o, , ]), numDeriv::jacobian(g, unname(p[pn])), tolerance = 1e-7)
  }
  fit <- suppressMessages(align(m, fits = 1))
  a <- fit$aligned[order(as.numeric(fit$aligned$time)), ]
  expect_lt(max(abs(a$estimate - log2(conc))), 0.3)
  expect_true(all(is.finite(fit$scaled$scaled[d$absolute])))
})

test_that("random terms of two crossed groupings recover both standard deviations", {
  set.seed(11)
  ab <- c("A1", "A2"); times <- sprintf("t%d", 1:6); gels <- sprintf("G%02d", 1:10)
  truth <- matrix(rnorm(12), 2, dimnames = list(ab, times))
  d <- expand.grid(time = times, antibody = ab, gel = gels, stringsAsFactors = FALSE)
  d$lane <- d$time
  scale <- setNames(rnorm(20), paste(rep(gels, each = 2), ab))
  lane <- setNames(rnorm(60, 0, 0.3), paste(rep(gels, each = 6), times))
  d$value <- truth[cbind(d$antibody, d$time)] + scale[paste(d$gel, d$antibody)] + lane[paste(d$gel, d$lane)] +
    rnorm(nrow(d), 0, 0.05)
  m <- alignmentModel(d, model = "yi + sj + lj", errmodel = "sigma", biological = yi ~ antibody * time,
                      scaling = list(sj ~ antibody + (1 | gel:antibody), lj ~ 0 + (1 | gel:lane)),
                      error = sigma ~ 1, trafo = "sigma <- exp(sigma)")
  # a condition per gel and antibody, the lanes within it; a block per antibody
  expect_length(unique(dMod2:::.alignmentConditions(m)), 20L)
  expect_equal(nrow(m$normalization$A), 2L)
  fit <- suppressMessages(align(m, fits = 1))
  expect_true(fit$fit$converged)
  # the standard deviations of the drawn effects, the scales within each antibody
  sds <- setNames(fit$omega$estimate, fit$omega$eta1)
  expect_equal(sds[["sj"]], sqrt(mean(tapply(scale, rep(ab, 10), var))), tolerance = 0.15)
  expect_equal(sds[["lj"]], sd(lane), tolerance = 0.15)
  expect_setequal(unique(fit$omega$grouping), c("gel:antibody", "gel:lane"))
  expect_setequal(unique(fit$random$grouping), c("gel:antibody", "gel:lane"))
  expect_equal(nrow(fit$random), 80L)
  expect_true(inherits(plotRandomEffects(fit$fit), c("gg", "list")))
  expect_error(alignmentModel(d, model = "yi + sj + lj", errmodel = "sigma", biological = yi ~ antibody * time,
                              scaling = list(sj ~ antibody + (1 | gr(gel, dist = laplace())),
                                             lj ~ 0 + (1 | gel:lane)),
                              error = sigma ~ 1), "several groupings are normal")
})

test_that("align() runs the fits of a multistart in parallel as one by one", {
  d <- .alignmentSim(5)
  m <- .rel(d)
  expect_error(align(m, cores = c(fits = 2, conditions = 2)), "number of fits")
  run <- function(model, cores) { set.seed(3); suppressMessages(align(model, fits = 3, cores = cores)) }
  one <- run(m, 1)
  two <- run(m, 2)
  expect_equal(two$multistart, one$multistart, tolerance = 1e-10)
  expect_equal(two$parameters, one$parameters, tolerance = 1e-10)
  r <- .rbg(.alignmentSim(6))
  one <- run(r, 1)
  two <- run(r, 2)
  expect_equal(two$multistart, one$multistart, tolerance = 1e-10)
  expect_equal(two$parameters, one$parameters, tolerance = 1e-10)
})
