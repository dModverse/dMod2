# scanL1(ssl = ): spike-and-slab lasso path, E-step, penalty and selection.

test_that("the E-step is the posterior share of the slab at its fixed point", {
  ssl <- list(lambda1 = 1, a = 1, b = 4)
  d   <- c(0, 0.01, 0.3, 1)
  es  <- dMod2:::.sslEstep(d, 200, ssl)
  psi <- function(d, l) l / 4 * exp(-l * d / 2)
  th  <- es$theta
  expect_equal(es$p, th * psi(d, 1) / (th * psi(d, 1) + (1 - th) * psi(d, 200)))
  expect_equal(th, (sum(es$p) + ssl$a - 1) / (length(d) + ssl$a + ssl$b - 2))
})

test_that("the derivative of the penalty in |d| is the EM weight", {
  ssl <- list(lambda1 = 2, a = 1, b = 3)
  d   <- c(0.05, 0.2, 0.7)
  th  <- 0.3
  es  <- dMod2:::.sslEstep(d, 50, ssl, theta = th, maxit = 0L)
  h   <- 1e-6
  num <- vapply(seq_along(d), function(j) {
    e <- replace(numeric(3), j, h)
    (dMod2:::.sslPenalty(d + e, th, 50, ssl) - dMod2:::.sslPenalty(d - e, th, 50, ssl)) / (2 * h)
  }, 0)
  expect_equal(num, es$p * 2 + (1 - es$p) * 50, tolerance = 1e-6)
})

test_that("pairwise terms and fuse weights follow the same order", {
  th <- c(r_B = 0.1, r_C = 0.4, r_D = 0.4, s_k = 0.5)
  g  <- list(list(pars = c("r_B", "r_C", "r_D"), anchor = 0))
  d  <- dMod2:::.sslTerms(th, "s_k", g)
  expect_identical(names(d), c("s_k", "r_B:r_C", "r_B:r_D", "r_C:r_D",
                               "r_B:anchor", "r_C:anchor", "r_D:anchor"))
  expect_equal(unname(d), c(0.5, 0.3, 0.3, 0, 0.1, 0.4, 0.4))
})

test_that("the spike-and-slab terms of a block are the gaps of sorted neighbours", {
  th <- c(r_B = 0.1, r_C = 0.4, r_D = -0.2, s_k = 0.5)
  g  <- list(list(pars = c("r_B", "r_C", "r_D"), anchor = 0))
  ed <- dMod2:::.sslEdges(th, "s_k", g)
  expect_identical(names(ed$d), c("s_k", "r_D:anchor", "anchor:r_B", "r_B:r_C"))
  expect_equal(unname(ed$d), c(0.5, 0.2, 0.1, 0.3))
  expect_equal(dMod2:::.sslGaps(th, "s_k", g, ed$edges), unname(ed$d))
})

test_that("a merge puts spike terms onto their kink, transitively", {
  th <- c(s_k = 0.01, s_m = 0.8, r_B = 0.004, r_C = 0.5, r_D = 0.52, r_E = 0.53)
  g  <- list(list(pars = c("r_B", "r_C", "r_D", "r_E"), anchor = 0))
  d  <- dMod2:::.sslEdges(th, c("s_k", "s_m"), g)$d
  p  <- stats::setNames(rep(1, length(d)), names(d))
  p[c("s_k", "anchor:r_B", "r_C:r_D", "r_D:r_E")] <- 0
  m  <- dMod2:::.sslMerge(th, p, c("s_k", "s_m"), g)
  expect_identical(unname(m[c("s_k", "s_m", "r_B")]), c(0, 0.8, 0))
  expect_equal(unname(m[c("r_C", "r_D", "r_E")]), rep(mean(c(0.5, 0.52, 0.53)), 3))
  # nothing to merge once every spike term sits on its kink
  p2 <- stats::setNames(c(0, 1, 0, 1, 0, 0), names(dMod2:::.sslEdges(m, c("s_k", "s_m"), g)$d))
  expect_null(dMod2:::.sslMerge(m, p2, c("s_k", "s_m"), g))
})

test_that("the spike-and-slab path fuses and anchors from one start per lambda", {
  set.seed(2)
  type  <- rep(c("R", "B", "C", "D"), each = 8)
  model <- .type_model(type)
  truth <- c(mu = 1, r_B = 0, r_C = 0.3, r_D = 0.3)
  obj   <- .lsq_obj(model, model(truth)$y + rnorm(length(type), 0, 0.05), 0.05)
  start <- c(mu = 0, r_B = 0, r_C = 0, r_D = 0)
  block <- list(r = list(pars = c("r_B", "r_C", "r_D"), anchor = 0))
  grid  <- 10^seq(0, 4, length.out = 13)

  fit <- scanL1(obj, start, groups = block, lambda = grid, fits = 3, pathFits = 1,
                ssl = list(lambda1 = 1), select = "plateau")
  expect_identical(fit$selected, "{r_B=0 | r_C,r_D}")
  # the last third of the grid holds one structure
  expect_true(all(tail(fit$path$key, 4) == fit$selected))
  # each lambda runs the EM from the previous mode and from the sparse point,
  # the downward pass from the mode above
  expect_true(all(fit$path$starts %in% 1:2) && any(fit$path$starts == 2L))
  # on the plateau only the weak slab shrinks: the estimate is close to the refit
  last <- fit$coefficients[nrow(fit$coefficients), c("r_B", "r_C", "r_D")]
  expect_equal(unname(last), unname(fit$fit[c("r_B", "r_C", "r_D")]), tolerance = 0.01)
  expect_identical(length(fit$full$values), 3L)
  expect_s3_class(plot(fit, type = "waterfall"), "ggplot")
  expect_identical(colnames(fit$inclusion)[-1],
                   c("r_B:r_C", "r_B:r_D", "r_C:r_D", "r_B:anchor", "r_C:anchor",
                     "r_D:anchor"))
  expect_true(all(fit$inclusion[nrow(fit$inclusion), c("r_B:r_C", "r_C:anchor")] > 0.99))
  expect_true(all(fit$inclusion[nrow(fit$inclusion), c("r_C:r_D", "r_B:anchor")] < 0.01))
  expect_s3_class(plot(fit, type = "clusters"), "ggplot")
  expect_s3_class(plot(fit, type = "inclusion"), "ggplot")
  expect_error(scanL1(obj, start, groups = block, q = 0.8, ssl = list()), "two different")
})

test_that("the spike-and-slab path removes the gate of an unneeded term", {
  set.seed(1)
  X     <- matrix(runif(60), 20, 3)
  model <- .gated_model(X)
  truth <- c(log10_k1 = log10(0.5), log10_k2 = -1, log10_k3 = log10(0.2),
             s_k1 = 1, s_k2 = 0, s_k3 = 1)
  obj   <- .lsq_obj(model, model(truth)$y + rnorm(20, 0, 0.01), 0.01)
  gates <- list(s_k1 = "log10_k1", s_k2 = "log10_k2", s_k3 = "log10_k3")
  fit <- scanL1(obj, c(log10_k1 = 0, log10_k2 = -0.5, log10_k3 = 0), zero = gates,
                lambda = 10^seq(0, 4, length.out = 9), fits = 3, pathFits = 1,
                ssl = list(lambda1 = 1), select = "plateau")
  expect_identical(fit$selected, "-s_k2")
  expect_identical(fit$fit[["s_k2"]], 0)
})

test_that("the toy decay of Hauber et al. clusters the two mutated cell types", {
  set.seed(3)
  times <- rep(seq(0, 100, by = 10), 3)
  type  <- rep(c("c1", "c2", "c3"), each = 11)
  model <- .decay_model(type, times)
  truth <- c(lk = -1.5, r_c2 = 0.2, r_c3 = 0.3)
  sigma <- 10^-1.3
  obj   <- .lsq_obj(model, model(truth)$y + rnorm(length(times), 0, sigma), sigma)
  block <- list(r = list(pars = c("r_c2", "r_c3"), anchor = 0))
  fit <- scanL1(obj, c(lk = -1, r_c2 = 0, r_c3 = 0), groups = block,
                lambda = 10^seq(0, 4, length.out = 13), fits = 5, pathFits = 1,
                ssl = list(lambda1 = 1), select = "plateau")
  expect_identical(fit$selected, "{r_c2,r_c3}")
})

test_that("spike terms within the snap distance are put onto their kink", {
  th <- c(r_B = 1e-9, r_C = 0.5, r_D = 0.5 + 1e-8)
  g  <- list(list(pars = c("r_B", "r_C", "r_D"), anchor = 0))
  p  <- c(0, 1, 0)
  m  <- dMod2:::.sslMerge(th, p, character(0), g, gap = 1e-6)
  expect_identical(m[["r_B"]], 0)
  expect_identical(m[["r_C"]], m[["r_D"]])
  expect_null(dMod2:::.sslMerge(c(r_B = 0.01, r_C = 0.5, r_D = 0.6), p, character(0), g, 1e-6))
})

test_that("a merge joins through spike edges that already sit on their kink", {
  th <- c(r_B = -2e-22, r_C = -1e-22, r_D = -1e-22, r_E = 0)
  g  <- list(list(pars = c("r_B", "r_C", "r_D", "r_E"), anchor = 0))
  m  <- dMod2:::.sslMerge(th, rep(0, 4), character(0), g, gap = 1e-6)
  expect_identical(unname(m), rep(0, 4))
})
