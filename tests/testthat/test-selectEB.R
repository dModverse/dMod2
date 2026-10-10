# selectEB(): score, g update, search, thresholds, refits and diagnostics.

.eb_env <- function(obj, spec, rule = "pip", prior = list(a = 1, g = 3)) {
  ctl <- list(trust = list(rinit = 0.1, rmax = 10, iterlim = 200), tol = 1e-3, kappa = 0.01,
              window = 2)
  dMod2:::.ebEnv(obj, spec, prior, ctl, rule, 0.05)
}

test_that("the score is the exact marginal likelihood of a linear model under the g-prior", {
  set.seed(1)
  n <- 30; X <- matrix(rnorm(n * 4), n, 4, dimnames = list(NULL, paste0("b", 1:4)))
  sigma <- 0.5
  y <- drop(X %*% c(1, 0, -0.5, 0)) + rnorm(n, 0, sigma)
  obj  <- .lsq_obj(.linear_model(X), y, sigma)
  spec <- dMod2:::.ebSpec(NULL, colnames(X), NULL)
  env  <- .eb_env(obj, spec, prior = list(a = 1, b = 4, g = 0))
  env$D0 <- sum((y / sigma)^2)
  start <- stats::setNames(rep(0.1, 4), colnames(X))
  for (z in list(c(TRUE, FALSE, TRUE, FALSE), c(TRUE, TRUE, TRUE, TRUE), c(FALSE, TRUE, FALSE, FALSE))) {
    s <- dMod2:::.ebEval(env, z, start)
    Xz <- X[, z, drop = FALSE]
    H  <- Xz %*% solve(crossprod(Xz), t(Xz))
    for (g in c(0.5, 4, 50)) {
      C <- sigma^2 * (diag(n) + g * H)
      m2 <- as.numeric(determinant(C)$modulus) + sum(y * solve(C, y)) + n * log(2 * pi)
      G  <- dMod2:::.ebScore(env, s, g) - dMod2:::.ebStructPrior(z, spec, env$prior)
      expect_equal(G - m2, -n * log(2 * pi * sigma^2), tolerance = 1e-6)
    }
  }
})

test_that("the update of g minimises the score at fixed structure", {
  s <- list(D = 40, p = 3, z = c(TRUE, TRUE, TRUE, FALSE))
  spec <- dMod2:::.ebSpec(NULL, paste0("b", 1:4), NULL)
  env <- .eb_env(function(...) NULL, spec)
  env$D0 <- 400
  num <- stats::optimize(function(g) dMod2:::.ebScore(env, s, g), c(0, 1e4), tol = 1e-10)$minimum
  expect_equal(dMod2:::.ebGhat(env, s), num, tolerance = 1e-5)
  env$D0 <- 41
  expect_identical(dMod2:::.ebGhat(env, s), 0)
})

test_that("the structure prior is beta-binomial per family and favours sparse structures", {
  spec <- dMod2:::.ebSpec(list(s_a = "la", s_b = "lb"), c("r1", "r2", "r3"), NULL)
  z <- c(TRUE, FALSE, TRUE, TRUE, FALSE)
  pr <- dMod2:::.ebStructPrior(z, spec, list(a = 1))
  expect_equal(pr, -2 * (lbeta(2, 3) - lbeta(1, 2)) - 2 * (lbeta(3, 4) - lbeta(1, 3)))
  zs <- lapply(0:3, function(n) c(FALSE, FALSE, seq_len(3) <= n))
  pz <- vapply(zs, dMod2:::.ebStructPrior, 0, spec = spec, prior = list(a = 1))
  expect_true(all(diff(pz) >= 0))
  expect_gt(pz[2], pz[1])
})

test_that("the search reaches the minimum over all structures", {
  f <- .gated_lin_problem(2, sigma = 0.3, k = c(0.6, 0.4, 1))
  fit <- selectEB(f$obj, f$center, zero = f$zero, fits = 6,
                  control = list(seed = 2, swaps = FALSE, fullFits = 3))
  spec <- fit$spec
  env <- .eb_env(f$obj, spec)
  env$D0 <- fit$empty$D
  J <- length(spec$cand)
  G <- vapply(0:(2^J - 1), function(i) {
    z <- as.logical(bitwAnd(i, 2^(0:(J - 1))))
    s <- dMod2:::.ebEval(env, z, f$center)
    dMod2:::.ebScore(env, s, dMod2:::.ebGhat(env, s))
  }, 0)
  expect_equal(fit$best$G, min(G), tolerance = 1e-4)
  expect_true(all(diff(fit$best$trace) < 0))
})

test_that("a gated linear model is selected with thresholds consistent with the score", {
  f <- .gated_lin_problem(1)
  fit <- selectEB(f$obj, f$center, zero = f$zero, fits = 4, control = list(seed = 1))
  expect_s3_class(fit, "selectEB")
  expect_identical(fit$best$key, f$truth)
  expect_identical(fit$selected, f$truth)
  tm <- fit$terms
  expect_identical(tm$on, tm$gain > tm$threshold)
  expect_identical(tm$on, tm$pip > 0.5)
  g <- fit$hyper$g
  off <- !tm$on
  # an absent candidate's threshold at the estimated g, from the score
  sp <- function(n) -2 * (lbeta(n + 1, 8 - n + 8) - lbeta(1, 8))
  expect_equal(tm$threshold[off], rep((1 + g) / g * (log1p(g) + sp(4) - sp(3)), sum(off)),
               tolerance = 1e-8)
  expect_equal(fit$hyper$penalty, (1 + g) / g * log1p(g))
  expect_true(all(fit$structures$G >= fit$best$G - 1e-8))
  expect_output(print(fit), "selected after refits")
  for (type in c("terms", "waterfall", "trace")) expect_s3_class(plot(fit, type = type), "ggplot")
})

test_that("rule alpha thresholds every present parameter at the chi-square quantile", {
  f <- .gated_lin_problem(1)
  fit <- selectEB(f$obj, f$center, zero = f$zero, rule = "alpha", alpha = 0.01, fits = 3,
                  control = list(seed = 1, swaps = FALSE))
  expect_equal(fit$terms$threshold, rep(qchisq(0.99, 1), 8), tolerance = 1e-8)
  expect_identical(fit$terms$on, fit$terms$gain > qchisq(0.99, 1))
  expect_true(all(is.na(fit$terms$pip)))
  expect_identical(fit$best$key, f$truth)
})

test_that("candidates that replace each other are reported", {
  t <- seq(0.2, 6, by = 0.2)
  m <- .paths_model(t)
  truth <- c(la = log10(0.5), lb = log10(0.5), lc = -1, s_a = 1, s_b = 0, s_c = 0)
  set.seed(3)
  y <- m(truth)$y + rnorm(length(t), 0, 0.01)
  obj <- .lsq_obj(m, y, 0.01)
  zero <- list(s_a = "la", s_b = "lb", s_c = "lc")
  fit <- selectEB(obj, c(la = -0.5, lb = -0.5, lc = -1), zero = zero, fits = 4,
                  control = list(seed = 3))
  expect_equal(sum(fit$best$z[1:2]), 1)
  expect_false(fit$best$z[3])
  expect_setequal(fit$alternatives$key[1:2], c("-s_b -s_c", "-s_a -s_c"))
  expect_lt(abs(diff(fit$alternatives$G[1:2])), 1e-3)
  # both paths present: the rates are collinear
  spec <- fit$spec
  env <- .eb_env(obj, spec)
  s <- dMod2:::.ebEval(env, c(TRUE, TRUE, FALSE), c(la = log10(0.3), lb = log10(0.2), lc = -1))
  grp <- dMod2:::.ebGroups(env, list(key = s$key, z = s$z))
  expect_length(grp, 1L)
  expect_setequal(grp[[1]]$candidates, c("s_a", "s_b"))
  expect_lt(grp[[1]]$eigenvalue, 1e-6)
})

test_that("differences of cell lines are selected as reference parameters", {
  lines <- c("R", "B", "C", "D")
  times <- rep(c(0.5, 1, 2, 4, 8), 4)
  line  <- rep(lines, each = 5)
  m <- .lines_model(line, times, lines)
  truth <- c(la = 0, lk = -0.5, ra_B = 0, ra_C = 0.3, ra_D = 0, rk_B = 0, rk_C = 0, rk_D = -0.4)
  set.seed(4)
  y <- m(truth)$y + rnorm(length(times), 0, 0.02)
  obj <- .lsq_obj(m, y, 0.02)
  ref <- c(paste0("ra_", lines[-1]), paste0("rk_", lines[-1]))
  fit <- selectEB(obj, c(la = 0, lk = 0, stats::setNames(rep(0, 6), ref)), reference = ref,
                  fits = 4, control = list(seed = 4))
  expect_identical(fit$terms$candidate[fit$terms$selected], c("ra_C", "rk_D"))
  expect_identical(unname(fit$fit[c("ra_B", "ra_D", "rk_B", "rk_C")]), rep(0, 4))
})

test_that("the tests use the data term of the objective", {
  f <- .gated_lin_problem(1)
  withPrior <- function(pars, fixed = NULL, ...) {
    out <- f$obj(pars, fixed = fixed, ...)
    attr(out, "data") <- out$value
    out$value <- out$value + 7
    out
  }
  a <- selectEB(withPrior, f$center, zero = f$zero, fits = 2, control = list(seed = 1, swaps = FALSE))
  b <- selectEB(withPrior, f$center, zero = f$zero, fits = 2, lrt = "total",
                control = list(seed = 1, swaps = FALSE))
  expect_equal(a$refits$data, a$refits$value - 7)
  expect_equal(b$refits$data, b$refits$value)
  expect_equal(a$refits$stat, b$refits$stat)
})

test_that("selectEB checks its arguments", {
  f <- .gated_lin_problem(1)
  expect_error(selectEB(f$obj, f$center), "no candidates")
  expect_error(selectEB(f$obj, f$center[-2], zero = f$zero), "center lacks")
  expect_error(selectEB(f$obj, f$center, zero = f$zero, class = c(foo = "x")), "named by candidates")
  expect_error(selectEB(f$obj, f$center, zero = f$zero, rule = "bic"))
})
