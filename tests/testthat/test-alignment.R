# alignmentModel() and align(): free expressions, normalization, derivatives, fits.

withr::local_options(dMod2.outdir = tempdir())

# readings of nGel gels at five times: `spread` the log standard deviation of the backgrounds
# between the gels, `deviation` the log deviations of the readings, drawn where NULL; attribute
# `truth` holds the value, scale and background of every reading
.alignmentSim <- function(nGel = 6, seed = 1, deviation = NULL, spread = 0.3) {
  set.seed(seed)
  times <- c(0, 10, 30, 60, 120)
  d <- expand.grid(time = times, gel = sprintf("G%02d", seq_len(nGel)), stringsAsFactors = FALSE)
  truth <- c(0.02, 1, 4, 2, 0.5)[match(d$time, times)]
  s <- rlnorm(nGel, 7, 0.5)[match(d$gel, unique(d$gel))]
  b <- rlnorm(nGel, 4, spread)[match(d$gel, unique(d$gel))]
  d$value <- (s * truth + b) * if (is.null(deviation)) rlnorm(nrow(d), 0, 0.1) else exp(deviation)
  d$time <- as.character(d$time)
  attr(d, "truth") <- data.frame(yi = truth, sj = s, bj = b)
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
  # one background in all gels, as the model has it: a background that differs between the gels
  # can put the smallest value at 0, where the errors are infinite
  d <- .alignmentSim(6, spread = 0)
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
  # the values and scales are at the optimum of the objective for the error parameters of the fit
  mean <- dMod2:::.alignmentMean(m)
  expect_lt(max(abs(fit$objfn(fit$fit$argument)$gradient[mean])), 1e-2)
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
  fit <- expect_no_warning(suppressMessages(align(m, fits = 1, init = c(yi = 1))))
  # without error parameters and random terms the restricted likelihood is the likelihood
  expect_identical(fit$remlRounds, 0L)
  expect_identical(attr(fit$parameters, "estimation"), "ML")
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

test_that("without random terms the standard errors come from the exact Hessian", {
  # the Hessian of the fit, without the second derivatives of the outputs, can be indefinite at a
  # minimum with large residuals; the standard errors take the exact one
  d <- data.frame(time = c("0", "1", "0", "1"), gel = c("A", "A", "B", "B"), value = c(100, 100, 100, 400))
  m <- alignmentModel(d, model = "yi * sj", errmodel = "0.1 * prediction", biological = yi ~ time,
                      scaling = sj ~ gel, error = NULL, trafo = .logs("yi", "sj"))
  fit <- suppressMessages(align(m, fits = 1))
  nm <- names(fit$fit$argument)
  expect_false(dMod2:::.positiveDefinite(fit$fit$hessian[nm, nm]))
  H <- numDeriv::hessian(function(q) fit$objfn(setNames(q, nm))$value, fit$fit$argument)
  V <- solve(H / 2)
  dimnames(V) <- list(nm, nm)
  X <- as.matrix(m$design$yi$X[!duplicated(d$time), , drop = FALSE])
  expect_equal(fit$aligned$se, unname(sqrt(diag(X %*% V[colnames(X), colnames(X)] %*% t(X)))), tolerance = 1e-5)
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
  # a random effect the fit sets to 0 at the boundary leaves the variances, its etas held at 0;
  # whether the error scale, drawn without a spread, ends there depends on the draw
  dropped <- sub("^omega_(.+)_\\1$", "\\1", fit$fit$boundary, perl = TRUE)
  expect_setequal(c(fit$omega$eta1, dropped), c("bj", "sigmaR"))
  expect_false(any(dropped %in% fit$omega$eta1))
  for (s in dropped)
    expect_equal(unname(fit$parameters[grep(paste0("^eta_", s, "_"), names(fit$parameters))]), rep(0, 8))
  expect_equal(fit$scaled$bloq, low)
  expect_equal(fit$aligned$lloq[fit$aligned$time == "0"], min(fit$scaled$scaled[low]))
  expect_true(all(is.na(fit$aligned$lloq[fit$aligned$time != "0"])))
  # FOCEI on the sparse Hessian: the same marginal, gradient and Hessian, and standard errors
  # beside the optimum, where the gradient is not of the size of the convergence of the modes
  ctl <- list(trustControl = list(iterlim = 0), innerControl = list(ftol = 1e-12, mtol = 1e-12, iterlim = 200))
  point <- fit$fit$argument + 0.05
  at <- function(sparse) em(fit$objfn, fit$randomEffects, point, fixed = fit$fit$fixed, control = ctl,
                            verbose = FALSE, sparse = sparse)
  dense <- at(FALSE)
  sp <- at(TRUE)
  nm <- names(dense$argument)
  # forward over reverse through the vjp of the alignment kernels, against forward over forward
  expect_identical(dense$derivMode, "forward-reverse")
  ff <- em(fit$objfn, fit$randomEffects, point, fixed = fit$fit$fixed,
           control = c(ctl, list(derivMode = "forward-forward")), verbose = FALSE)
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
                            suppressMessages(align(m, fits = 1, reml = FALSE)))
  expect_s4_class(fs$fit$hessian, "sparseMatrix")
  # the standard errors of the sparse solve are those of the dense one on the same Hessian
  free <- intersect(names(fs$fit$argument), rownames(fs$fit$hessian))
  X <- m$design$yi$X[!duplicated(m$data$time), , drop = FALSE]
  H <- if (identical(fs$seSource, "Hessian")) fs$fit$hessian
       else dMod2:::.alignmentGaussNewton(fs$model, fs$objfn, fs$randomEffects, fs$fit)
  expect_s4_class(H, "sparseMatrix")
  expect_equal(fs$aligned$se, sqrt(dMod2:::.alignmentVariance(as.matrix(H), free)(X)))
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

test_that("technical replicates are averaged per sample, and the fit is that of the means taken by hand", {
  d <- .alignmentSim(4)
  set.seed(12)
  lanes <- d[rep(seq_len(nrow(d)), each = 3), ]
  lanes$lane <- rep(1:3, nrow(d))
  lanes$value <- lanes$value * rlnorm(nrow(lanes), 0, 0.15)
  lanes <- lanes[sample(nrow(lanes)), ]
  # a relative error model: the geometric mean
  m <- .rel(lanes, technical = "lane")
  hand <- aggregate(value ~ time + gel, lanes, function(v) exp(mean(log(v))))
  expect_equal(m$technical$readings, nrow(lanes))
  expect_false("lane" %in% names(m$data))
  expect_equal(nrow(m$data), nrow(d))
  key <- function(x) paste(x$gel, x$time)
  hand <- hand[match(key(m$data), key(hand)), ]
  expect_equal(m$data$value, hand$value)
  fit <- suppressMessages(align(m, fits = 1))
  ref <- suppressMessages(align(.rel(hand), fits = 1))
  expect_equal(fit$value, ref$value, tolerance = 1e-10)
  expect_equal(fit$parameters, ref$parameters, tolerance = 1e-10)
  expect_equal(fit$aligned, ref$aligned, tolerance = 1e-10)
  # an absolute error model: the arithmetic mean; a sample with censored readings is censored at
  # the mean of its readings with the limits in place of the censored ones
  lanes$value <- log2(lanes$value)
  lanes$lloq <- NA_real_
  low <- lanes$gel == "G01" & lanes$time %in% c("0", "10")
  # every reading of the two samples is censored but the first lane at time 10
  lanes$lloq[low] <- max(lanes$value[low & lanes$lane != 1]) + 0.01
  lanes$lloq[low & lanes$time == "10" & lanes$lane == 1] <- -Inf
  add <- alignmentModel(lanes, model = "yi + sj", errmodel = "sigma", biological = yi ~ time, scaling = sj ~ gel,
                        error = sigma ~ 1, trafo = "sigma <- exp(sigma)", technical = "lane")
  hand <- aggregate(value ~ time + gel, lanes, mean)
  at <- match(key(add$data), key(hand))
  cens <- add$data$gel == "G01" & add$data$time %in% c("0", "10")
  expect_equal(add$data$value[!cens], hand$value[at][!cens])
  expect_true(all(is.na(add$data$value[cens])))
  expect_equal(add$censored, cens)
  b <- lanes[low, ]
  bound <- tapply(ifelse(b$value <= b$lloq, b$lloq, b$value), b$time, mean)
  expect_equal(add$data$lloq[cens], as.vector(bound[add$data$time[cens]]))
  expect_error(.rel(lanes, technical = "gel"), "technical column gel")
  expect_error(.rel(lanes, technical = "plate"), "names columns")
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
  d <- .alignmentSim(8)
  m <- .rbg(d)
  obj <- dMod2:::.alignmentObjective(m) + dMod2:::.alignmentNormalization(m)
  re <- dMod2:::.alignmentRandomEffects(m)
  st <- dMod2:::.alignmentStart(m, dMod2:::.alignmentPars(m)$fixed)
  fit <- em(obj, re, st, control = list(trustControl = list(mtol = 1e-4)), sparse = FALSE, verbose = FALSE)
  # the parameters the readings were drawn from, the values at a geometric mean of 1; -2 log L
  # there lies above the optimum and below the saddle
  tr <- attr(d, "truth")
  shift <- mean(log(unique(tr$yi)))
  target <- list(yi = log(tr$yi) - shift, sj = log(tr$sj) + shift, bj = rep(4, nrow(d)),
                 sigmaR = rep(log(0.1), nrow(d)))
  truth <- c(unlist(lapply(names(target), function(k) {
    X <- as.matrix(m$design[[k]]$X)
    setNames(qr.solve(X, target[[k]]), colnames(X))
  })), omega_bj_bj = log(0.3))
  at <- em(obj, re, truth[names(fit$argument)], control = list(trustControl = list(iterlim = 0)), sparse = FALSE,
           verbose = FALSE)
  expect_lt(fit$value, at$value)
})

test_that("the FOCEI gradient and Hessian are those of the FOCEI objective, censored data too", {
  d <- .alignmentSim(8)
  d$lloq <- NA_real_
  low <- d$time == "0" & d$gel %in% c("G01", "G02", "G03", "G04")
  d$lloq[low] <- d$value[low] * 1.2
  fit <- suppressMessages(align(.rbg(d), fits = 2, reml = FALSE))
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
  # error scales that spread by a factor of 20 between the gels: the fit keeps their random effect
  sc <- rep(exp(seq(-1.5, 1.5, length.out = 8)), each = 5)
  set.seed(5)
  dev <- rnorm(40, 0, 0.1 * sc)
  d <- .alignmentSim(8, deviation = dev)
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
  fo <- suppressMessages(align(.rbg(d), fits = 2, reml = FALSE))
  qu <- suppressMessages(align(.rbg(d), fits = 2, method = "quadrature"))
  expect_identical(attr(qu$parameters, "estimation"), "ML")
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
  # an optimised BLAS rounds by the memory layout of its arguments, so these fits agree to their
  # stopping tolerance only
  expect_equal(two$multistart$value, one$multistart$value, tolerance = 1e-8)
  expect_equal(two$parameters, one$parameters, tolerance = 1e-6)
})

test_that("loadDLL() finds the compiled expressions of an alignment from any working directory", {
  m <- .rel(.alignmentSim(5))
  obj <- dMod2:::.alignmentObjective(m)
  files <- withr::with_dir(withr::local_tempdir(), withr::with_options(list(dMod2.outdir = NULL),
                                                                       dMod2:::.dllFiles(list(obj))))
  expect_true(paste0(attr(m$fn, "modelname"), .Platform$dynlib.ext) %in% basename(files))
})

test_that("the restricted likelihood of a model linear in the mean parameters gives sigma^2 = RSS / (n - p)", {
  set.seed(7)
  times <- c("0", "10", "30", "60"); gels <- c("G1", "G2", "G3")
  d <- expand.grid(time = times, gel = gels, stringsAsFactors = FALSE)
  d$value <- c(0, 1, 2, 1)[match(d$time, times)] + rnorm(3)[match(d$gel, gels)] + rnorm(12, 0, 0.3)
  m <- alignmentModel(d, model = "yi + sj", errmodel = "sigma", biological = yi ~ time, scaling = sj ~ gel,
                      error = sigma ~ 1, trafo = "sigma <- exp(sigma)")
  ml <- align(m, fits = 1, reml = FALSE)
  re <- align(m, fits = 1)
  expect_identical(attr(re$parameters, "estimation"), "REML")
  expect_identical(attr(ml$parameters, "estimation"), "ML")
  rss <- function(fit) sum((fit$scaled$value - fit$scaled$prediction)^2)
  # p counts the values and the scales, less one for the normalization
  expect_equal(exp(2 * re$parameters[["sigma"]]), rss(re) / (12 - 6), tolerance = 1e-4)
  expect_equal(exp(2 * ml$parameters[["sigma"]]), rss(ml) / 12, tolerance = 1e-4)
  expect_equal(re$aligned$estimate, ml$aligned$estimate, tolerance = 1e-6)
  expect_lte(re$remlRounds, 3L)
  expect_identical(ml$remlRounds, 0L)
})

test_that("the restricted likelihood through FOCEI: sigma near the drawn noise, its value, gradient and sparse path", {
  set.seed(1)
  times <- c(0, 10, 30, 60, 120, 240); conc <- c(1, 40, 85, 60, 25, 10)
  gels <- sprintf("G%02d", 1:10)
  shift <- setNames(rnorm(10, 3, 1.5), gels)
  background <- setNames(1.5 + rnorm(10, 0, 0.3), gels)
  wb <- do.call(rbind, lapply(gels, function(g) data.frame(method = "WB", gel = g, time = sort(sample(times, 4)))))
  noiseWB <- rnorm(nrow(wb), 0, 0.4)
  wb$value <- log2(conc[match(wb$time, times)] + 2^background[wb$gel]) + shift[wb$gel] + noiseWB
  ms <- data.frame(method = "MS", gel = NA, time = rep(c(0, 30, 60, 240), each = 2))
  ms$value <- log2(conc[match(ms$time, times)]) + rnorm(nrow(ms), 0, 0.35)
  d <- rbind(cbind(wb, absolute = FALSE), cbind(ms, absolute = TRUE))
  m <- alignmentModel(d, model = "log2(2^yi + 2^(b + bj)) + sj", errmodel = "sigma",
                      biological = yi ~ time, scaling = list(sj ~ gel, b ~ 1, bj ~ 0 + (1 | gel)),
                      error = sigma ~ method, trafo = "sigma <- exp(sigma)",
                      prior = constraintL2(mu = c(omega_bj_bj = log(0.3), b = 0), sigma = c(omega_bj_bj = 0.3, b = 2)))
  ml <- suppressMessages(align(m, fits = 1, reml = FALSE))
  re <- suppressMessages(align(m, fits = 1))
  sigmas <- function(fit) c(exp(fit$parameters[["sigma"]]), exp(fit$parameters[["sigma"]] + fit$parameters[["sigma_methodWB"]]))
  # restricted likelihood removes the downward bias of the maximum likelihood estimate
  expect_true(all(sigmas(re) > sigmas(ml)))
  # the relative method near the standard deviation of its noise as drawn; the eight readings of
  # the absolute method leave its sigma far from that in some draws
  expect_lt(abs(sigmas(re)[2] / sqrt(mean(noiseWB^2)) - 1), 0.35)
  expect_gte(re$remlRounds, 2L)
  # the term of a round: log|S / 2| - p log(2 pi) at the parameters and modes of the last round, the
  # error parameters and the hyperparameter replaced; FOCEI with it has the exact gradient, dense
  # and sparse
  design <- dMod2:::.emDesign(re$objfn, re$randomEffects)
  joint <- design$objfn + dMod2:::.normalPrior(design)
  mean <- intersect(dMod2:::.alignmentMean(m), names(re$fit$argument))
  theta <- setdiff(names(re$fit$argument), mean)
  expect_setequal(theta, c("sigma", "sigma_methodWB", "omega_bj_bj"))
  anchor <- c(re$fit$argument, re$fit$eta)
  term <- dMod2:::.remlAnchor(design, anchor, mean, theta)
  point <- re$fit$argument + 0.05
  full <- anchor; full[theta] <- point[theta]
  # sigma without the values and scales: the information is the Gauss-Newton Hessian of the joint
  H <- as.matrix(joint(pars = full, deriv = TRUE)$hessian)
  V <- c(mean, re$randomEffects$etaNames$gel[, 1])
  ee <- function(e) as.numeric(determinant(H[e, e, drop = FALSE])$modulus)
  expect_equal(term(point)$value, ee(V) - sum(vapply(V[-seq_along(mean)], ee, 0)) - length(mean) * log(4 * pi),
               tolerance = 1e-10)
  expect_equal(unname(term(point)$gradient[setdiff(names(point), theta)]), rep(0, length(mean)))
  ctl <- list(trustControl = list(iterlim = 0), innerControl = list(ftol = 1e-12, mtol = 1e-12, iterlim = 200))
  at <- function(p, sparse = FALSE)
    suppressMessages(em(re$objfn + term, re$randomEffects, p, control = ctl, verbose = FALSE, sparse = sparse))
  r0 <- at(point)
  sel <- c(theta, "b", "yi_time30", "sj_gelG03")
  ofv <- function(q) { p <- point; p[sel] <- q; at(p)$value }
  num <- numDeriv::grad(ofv, unname(point[sel]), method.args = list(eps = 1e-4, d = 1e-3, r = 2))
  expect_equal(unname(r0$gradient[sel]), num, tolerance = 1e-4)
  sp <- at(point, sparse = TRUE)
  expect_equal(sp$value, r0$value, tolerance = 1e-8)
  expect_equal(sp$gradient[names(point)], r0$gradient[names(point)], tolerance = 1e-6)
})

test_that("the gradient of the restricted term is exact with censored readings and sigma of the prediction", {
  d <- .alignmentSim(6)
  d$lloq <- NA_real_
  low <- d$time == "0" & d$gel %in% c("G01", "G02")
  d$lloq[low] <- d$value[low] * 1.2
  m <- .rbg(d)
  obj <- dMod2:::.alignmentObjective(m) + dMod2:::.alignmentNormalization(m)
  design <- dMod2:::.emDesign(obj, dMod2:::.alignmentRandomEffects(m))
  joint <- design$objfn + dMod2:::.normalPrior(design)
  ap <- dMod2:::.alignmentPars(m)
  set.seed(3)
  # distinct values: with all values equal the background and the scales are confounded
  full <- setNames(c(dMod2:::.alignmentStart(m, ap$fixed) + rnorm(length(ap$fixed), 0, 0.5),
                     rnorm(length(ap$eta), 0, 0.2)), c(ap$fixed, ap$eta))
  full <- c(full, omega_bj_bj = log(0.3))
  mean <- dMod2:::.alignmentMean(m)
  theta <- c("sigmaR", "omega_bj_bj")
  term <- function(f, deriv = TRUE) dMod2:::.remlTerm(design, f, NULL, mean, theta, deriv = deriv)
  g <- term(full)$gradient
  h <- 1e-4
  num <- vapply(theta, function(k) {
    e <- setNames(numeric(length(full)), names(full)); e[k] <- h
    (term(full + e, FALSE)$value - term(full - e, FALSE)$value) / (2 * h)
  }, 0)
  expect_equal(g[theta], num, tolerance = 1e-7)
})

test_that("alignment models build in the output directory of the session, by forks of the same model at once", {
  skip_on_os("windows")
  dir <- file.path(tempdir(), "alignOutdir")
  dir.create(dir, showWarnings = FALSE)
  withr::local_options(dMod2.outdir = dir)
  d <- .alignmentSim(3)
  build <- function() alignmentModel(d, model = "1.5 * yi * sj", errmodel = "sigmaR * prediction",
                                     biological = yi ~ time, scaling = sj ~ gel, error = sigmaR ~ 1,
                                     trafo = .logs("yi", "sj", "sigmaR"))
  value <- function(m) {
    ap <- dMod2:::.alignmentPars(m)$fixed
    sum(dMod2:::.alignmentEval(m, 1:3, setNames(rep(0.1, length(ap)), ap), deriv = FALSE)$y)
  }
  res <- parallel::mclapply(1:2, function(i) tryCatch(value(build()), error = conditionMessage), mc.cores = 2)
  expect_true(all(vapply(res, is.numeric, TRUE)))
  expect_equal(res[[1]], res[[2]])
  expect_true(any(grepl("align_.*\\.so$", list.files(dir, recursive = TRUE))))
})

test_that("a fit reports the modes at its parameters, a positive definite Hessian and finite standard errors", {
  times <- c(0, 10, 30, 60, 120, 240, 480, 1440, 2880)
  signal <- c(-2.6, 1, 2.3, 2.5, 2, 1.2, 0.8, 0.2, -0.2)
  gels <- sprintf("G%02d", 1:50)
  course <- setNames(rep(c("short", "long"), 25), gels)
  drop <- rep(rep_len(c(0, 0, 1, 0, 1, 2, 0, 1, 0, 1), 25), each = 2)
  lanes <- lapply(setNames(seq_along(gels), gels), function(k) {
    tk <- if (course[k] == "short") times[1:6] else times[c(1, 4, 6:9)]
    tk[seq_along(tk) > drop[k]]
  })
  set.seed(27)
  shift <- setNames(ifelse(course == "short", 2, -2) + rnorm(50, 0, 1.5), gels)
  background <- setNames(-1.8 + rnorm(50, 0, 0.3), gels)
  d <- do.call(rbind, lapply(gels, function(g) data.frame(gel = g, time = lanes[[g]])))
  d$value <- log2(2^signal[match(d$time, times)] + 2^background[d$gel]) + shift[d$gel] + rnorm(nrow(d), 0, 0.4)
  m <- alignmentModel(d, model = "log2(2^yi + 2^(b + bj)) + sj", errmodel = "sigma",
                      biological = yi ~ time, scaling = list(sj ~ gel, b ~ 1, bj ~ 0 + (1 | gel)),
                      error = sigma ~ 1, trafo = "sigma <- exp(sigma)",
                      prior = constraintL2(mu = c(omega_bj_bj = log(0.3), b = -2), sigma = c(omega_bj_bj = 0.3, b = 1)))
  fit <- suppressMessages(align(m, fits = 1, reml = FALSE))
  # the modes the fit reports are those at its parameters, and its Hessian there is positive definite
  ctl <- list(trustControl = list(iterlim = 0), innerControl = list(ftol = 1e-12, mtol = 1e-12, iterlim = 200))
  at <- suppressMessages(em(fit$objfn, fit$randomEffects, fit$fit$argument, control = ctl, verbose = FALSE))
  expect_equal(fit$fit$etaModes, at$etaModes, tolerance = 1e-4)
  H <- as.matrix(fit$fit$hessian)
  expect_gt(min(eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values), 0)
  expect_identical(fit$seSource, "Hessian")
  expect_true(all(is.finite(fit$aligned$se)))
  # the Gauss-Newton Schur complement, which takes over where the Hessian is not positive
  # definite, gives finite standard errors of the same order
  gn <- dMod2:::.alignmentGaussNewton(fit$model, fit$objfn, fit$randomEffects, fit$fit)
  on <- names(fit$fit$argument)
  X <- m$design$yi$X[!duplicated(m$data$time), , drop = FALSE]
  se <- sqrt(dMod2:::.alignmentVariance(as.matrix(gn), on)(X))
  expect_true(all(is.finite(se)))
  expect_lt(max(abs(log(se / fit$aligned$se))), log(2))
})

test_that("the Gauss-Newton fallback is the expected information, positive definite where the joint Hessian is not", {
  times <- c("0", "1", "2")
  d <- expand.grid(time = times, gel = c("A", "B", "C"), stringsAsFactors = FALSE)
  d$value <- 100 * exp(c(-0.4, 0, 0.4))[match(d$gel, c("A", "B", "C"))] * c(1, 2, 4)[match(d$time, times)]
  d$value[d$time == "2" & d$gel == "A"] <- 4 * d$value[d$time == "2" & d$gel == "A"]
  m <- alignmentModel(d, model = "yi * sj", errmodel = "0.1 * prediction", biological = yi ~ time,
                      scaling = sj ~ 1 + (1 | gel), error = NULL, trafo = .logs("yi", "sj"))
  fit <- suppressMessages(align(m, fits = 1, reml = FALSE))
  f <- fit$fit
  eta <- dMod2:::.etaVector(f)
  on <- names(f$argument)
  e <- names(eta)
  free <- c(on, e)
  # the Schur complement of the first-order Hessian of the joint at the modes is indefinite
  design <- dMod2:::.emDesign(fit$objfn, fit$randomEffects)
  joint <- design$objfn + dMod2:::.normalPrior(design)
  H <- as.matrix(joint(pars = c(f$argument, eta), fixed = f$fixed, hessian = "sparse")$hessian)[free, free]
  first <- H[on, on] - H[on, e] %*% solve(H[e, e], H[e, on])
  expect_lt(min(eigen(first, symmetric = TRUE, only.values = TRUE)$values), 0)
  gn <- as.matrix(dMod2:::.alignmentGaussNewton(m, fit$objfn, fit$randomEffects, f))
  expect_gt(min(eigen(gn, symmetric = TRUE, only.values = TRUE)$values), 0)
  # the information of the readings, 2 / sigma^2 in the prediction and 4 / sigma^2 in sigma
  n <- nrow(d)
  out <- function(q) as.vector(dMod2:::.alignmentEval(m, seq_len(n), c(setNames(q, free), f$fixed), deriv = FALSE)$y)
  q <- unname(c(f$argument, eta)[free])
  J <- numDeriv::jacobian(out, q)
  s <- out(q)[n + seq_len(n)]
  expect_equal(unname(as.matrix(dMod2:::.alignmentFisher(m, c(setNames(q, free), f$fixed), free))),
               2 * crossprod(J[seq_len(n), ] / s) + 4 * crossprod(J[n + seq_len(n), ] / s), tolerance = 1e-6)
})
