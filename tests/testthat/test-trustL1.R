# trustL1() on quadratic objectives, whose L1-penalized optimum is the
# soft-thresholding solution on both sides of the kink.


# Quadratic objfn centered at `target` with unit "sigma". Returns an
# objlist so trustL1's `+.objlist` accumulation works correctly.
.quadratic_objfn <- function(target) {
  d <- length(target)
  function(p, ...) {
    gr <- 2 * (p - target); names(gr) <- names(p)
    hs <- 2 * diag(d); dimnames(hs) <- list(names(p), names(p))
    objlist(value = sum((p - target)^2), gradient = gr, hessian = hs)
  }
}


## ---- Soft-threshold above the kink: |x0 - mu| > lambda/2 ---------------

test_that("trustL1 lands on soft(x0, lambda/2) when above the L1 kink", {
  x0 <- 1.5
  lambda <- 0.4
  obj <- .quadratic_objfn(c(x = x0))

  fit <- trustL1(obj, parinit = c(x = 0), mu = c(x = 0),
                 lambda = lambda, rinit = 0.5, rmax = 5, iterlim = 100)
  # x0 lies outside the dead zone, so the soft-threshold is active.
  expected <- x0 - sign(x0 - 0) * lambda / 2
  expect_equal(unname(fit$argument[["x"]]), expected, tolerance = 1e-4)
})


## ---- Multivariate: soft-threshold acts elementwise --------------------

test_that("trustL1 with diagonal quadratic applies soft-threshold elementwise", {
  x0 <- c(a = 1.0, b = 0.05, c = -0.8)
  lambda <- 0.3
  obj <- .quadratic_objfn(x0)
  mu_L1 <- c(a = 0, b = 0, c = 0)

  fit <- trustL1(obj, parinit = c(a = 0, b = 0, c = 0),
                 mu = mu_L1, lambda = lambda,
                 rinit = 0.5, rmax = 5, iterlim = 200)

  # Each coordinate: soft(x0_i, lambda/2) w.r.t. mu_i = 0.
  soft <- function(z, t) if (abs(z) <= t) 0 else z - sign(z) * t
  expected <- vapply(x0, soft, t = lambda / 2, FUN.VALUE = 0.0)
  expect_equal(unname(fit$argument[names(x0)]),
               unname(expected), tolerance = 5e-3)
})


## ---- Lambda = 0 reduces to plain trust ---------------------------------

test_that("trustL1 with lambda = 0 converges to the unpenalized minimum", {
  x0 <- c(a = 0.7, b = -1.2)
  obj <- .quadratic_objfn(x0)
  fit <- trustL1(obj, parinit = c(a = 0, b = 0),
                 mu = c(a = 0, b = 0), lambda = 0,
                 rinit = 0.5, rmax = 5, iterlim = 100)
  expect_equal(unname(fit$argument[names(x0)]), unname(x0),
               tolerance = 1e-5)
})


## ---- Sub-kink coords are pinned to mu exactly -------------------------

test_that("trustL1 pins sub-kink coords to mu exactly via the active set", {
  # x0 places 'b' and 'd' inside the soft-threshold dead zone, 'a' and 'c' outside.
  x0     <- c(a = 1.5, b = 0.10, c = -0.8, d = -0.05)
  lambda <- 0.5
  obj    <- .quadratic_objfn(x0)
  mu_L1  <- c(a = 0, b = 0, c = 0, d = 0)

  fit <- trustL1(obj, parinit = c(a = 0, b = 0, c = 0, d = 0),
                 mu = mu_L1, lambda = lambda,
                 rinit = 0.5, rmax = 5, iterlim = 200)

  # b and d must land at mu = 0 exactly (active-set pinning, not just
  # numerically close).
  expect_identical(unname(fit$argument[["b"]]), 0)
  expect_identical(unname(fit$argument[["d"]]), 0)

  # a and c follow the analytic soft-threshold solution.
  expect_equal(unname(fit$argument[["a"]]),
               1.5 - lambda / 2, tolerance = 1e-6)
  expect_equal(unname(fit$argument[["c"]]),
               -0.8 + lambda / 2, tolerance = 1e-6)

  # Reported value equals smooth-part value at theta plus L1 contribution.
  smooth_at_theta <- sum((fit$argument - x0)^2)
  l1_at_theta    <- lambda * sum(abs(fit$argument - mu_L1))
  expect_equal(fit$value, smooth_at_theta + l1_at_theta, tolerance = 1e-8)
})


## ---- Sub-kink coords stay pinned when parinit starts off the kink -----

test_that("trustL1 pulls sub-kink coords to mu when started away from it", {
  # Same dead-zone setup, but initialize far from mu so the optimizer has
  # to *cross* the kink to pin the parameter.
  x0     <- c(a = 1.5, b = 0.10)
  lambda <- 0.5
  obj    <- .quadratic_objfn(x0)
  mu_L1  <- c(a = 0, b = 0)

  fit <- trustL1(obj, parinit = c(a = 2.0, b = 1.0),
                 mu = mu_L1, lambda = lambda,
                 rinit = 0.5, rmax = 5, iterlim = 200)

  expect_identical(unname(fit$argument[["b"]]), 0)
  expect_equal(unname(fit$argument[["a"]]),
               1.5 - lambda / 2, tolerance = 1e-6)
})


## ---- One-sided penalty acts below mu only -------------------------------

test_that("trustL1 one-sided penalty acts below mu only", {
  # lambda * max(0, mu - p): below mu the soft threshold, at mu where the
  # objective pulls down by at most lambda, above mu nothing
  x0     <- c(a = -1.0, b = 0.5, c = -0.1, d = 0.1)
  lambda <- 0.4
  obj    <- .quadratic_objfn(x0)
  mu_L1  <- c(a = 0, b = 0, c = 0, d = 0)

  fit <- trustL1(obj, parinit = c(a = 0.5, b = 0.5, c = 0.5, d = 0),
                 mu = mu_L1, one.sided = TRUE, lambda = lambda,
                 rinit = 0.5, rmax = 5, iterlim = 200)

  expect_equal(unname(fit$argument[["a"]]), -1 + lambda / 2, tolerance = 1e-6)
  expect_equal(unname(fit$argument[["b"]]), 0.5, tolerance = 1e-6)
  expect_identical(unname(fit$argument[["c"]]), 0)
  # started on the kink with the objective pulling up, d leaves it
  expect_equal(unname(fit$argument[["d"]]), 0.1, tolerance = 1e-6)
})


## ---- Box bounds alongside the L1 kink ---------------------------------
# The kink active set is independent of the box and pins to mu bit-exactly.

test_that("trustL1 combines a box bound with an L1 kink", {
  # a is capped by the upper bound; b lies inside the dead zone and stays at mu.
  obj <- function(p, ...) {
    d <- as.numeric(p - c(3, 0.05))
    list(value = sum(d^2), gradient = 2 * d, hessian = 2 * diag(2))
  }
  fit <- trustL1(obj, c(a = 0.5, b = 1), mu = c(a = 0, b = 0), lambda = 1,
                 rinit = 1, rmax = 10, iterlim = 100,
                 parupper = c(a = 1, b = Inf))

  expect_true(fit$converged)
  expect_identical(unname(fit$argument[["b"]]), 0)
  expect_equal(unname(fit$argument[["a"]]), 1, tolerance = 1e-6)
  expect_lt(fit$argument[["a"]], 1)
  expect_true(fit$atBound[["a"]])
})


test_that("trustL1 rejects a kink outside the box", {
  obj <- function(p, ...) list(value = sum(p^2), gradient = 2 * p, hessian = 2 * diag(2))
  expect_error(
    trustL1(obj, c(a = 0.5, b = 0.5), mu = c(a = 2), lambda = 1,
            rinit = 1, rmax = 10, parupper = c(a = 1, b = Inf)),
    "mu must lie strictly inside")
})


test_that("trustL1 boundary = 'clip' reproduces the historical pinning", {
  obj <- function(p, ...) {
    d <- as.numeric(p - c(2, 0.05))
    list(value = sum(d^2), gradient = 2 * d, hessian = 2 * diag(2))
  }
  refl <- trustL1(obj, c(a = 1, b = 1), mu = c(a = 0, b = 0), lambda = 1,
                  rinit = 1, rmax = 10, boundary = "reflective")
  clip <- trustL1(obj, c(a = 1, b = 1), mu = c(a = 0, b = 0), lambda = 1,
                  rinit = 1, rmax = 10, boundary = "clip")
  expect_identical(unname(refl$argument[["b"]]), 0)
  expect_identical(unname(clip$argument[["b"]]), 0)
  expect_equal(unname(refl$argument), unname(clip$argument), tolerance = 1e-8)
})


## ---- gates and fusion blocks --------------------------------------------

test_that("a gate stays on the kink unless the data pull it upward", {
  for (case in list(c(y = -1, s = 0), c(y = 0.3, s = 0), c(y = 2, s = 1.5))) {
    fit <- trustL1(.quadratic_objfn(c(s = case[["y"]])), parinit = c(s = 1),
                   mu = c(s = 0), lambda = 1, gate = "s")
    expect_identical(fit$argument[["s"]] == 0, case[["s"]] == 0)
    expect_equal(fit$argument[["s"]], case[["s"]], tolerance = 1e-6)
  }
  start_below <- trustL1(.quadratic_objfn(c(s = -1)), parinit = c(s = -0.5),
                         mu = c(s = 0), lambda = 1, gate = "s")
  expect_identical(start_below$argument[["s"]], 0)
})

test_that("a fusion block shrinks toward fused values and fuses exactly", {
  y   <- c(a = 0, b = 1, c = 5)
  obj <- .quadratic_objfn(y)
  fit <- function(lambda, start = c(a = 0.5, b = -1, c = 2))
    trustL1(obj, parinit = start, fuse = list(names(y)), lambda = lambda,
            rinit = 0.5, rmax = 5, iterlim = 200)$argument

  # No fusion: each value moves by lambda/2 per neighbor above minus below.
  expect_equal(unname(fit(0.8)), c(0.8, 1, 4.2), tolerance = 1e-6)

  # a and b fuse and stay fused; c keeps a separate value.
  p <- fit(1.2)
  expect_identical(p[["a"]], p[["b"]])
  expect_equal(unname(p), c(1.1, 1.1, 3.8), tolerance = 1e-6)

  # Strong penalty: everything at the mean, also from a fused start.
  expect_equal(unname(fit(10)), rep(2, 3), tolerance = 1e-6)
  p <- fit(1.2, start = c(a = 2, b = 2, c = 2))
  expect_equal(unname(p), c(1.1, 1.1, 3.8), tolerance = 1e-6)
})

test_that("an anchored block soft-thresholds toward the anchor", {
  y <- c(a = 0.3, b = -2)
  w <- matrix(1, 3, 3); w[1:2, 1:2] <- 0
  fit <- trustL1(.quadratic_objfn(y), parinit = c(a = 0, b = 0),
                 fuse = list(list(pars = names(y), anchor = 0, weights = w)),
                 lambda = 1, rinit = 0.5, rmax = 5)
  expect_identical(fit$argument[["a"]], 0)
  expect_equal(fit$argument[["b"]], -1.5, tolerance = 1e-6)
})

test_that("an anchored fused block releases only the members the data pull away", {
  y <- c(a = 0.1, b = 0.2, c = 3)
  w <- matrix(1, 4, 4); w[4, 1:2] <- w[1:2, 4] <- 3
  fit <- trustL1(.quadratic_objfn(y), parinit = c(a = 0, b = 0, c = 0),
                 fuse = list(list(pars = names(y), anchor = 0, weights = w)),
                 lambda = 0.3, rinit = 0.5, rmax = 5, iterlim = 200)
  expect_identical(fit$argument[["a"]], 0)
  expect_identical(fit$argument[["b"]], 0)
  # the released member balances the anchor and the pinned members
  expect_equal(fit$argument[["c"]], 3 - 0.45, tolerance = 1e-6)
})

test_that("gates and fuse blocks reject the clip boundary and overlap with mu", {
  obj <- .quadratic_objfn(c(a = 1, b = 2))
  expect_error(trustL1(obj, c(a = 0, b = 0), fuse = list(c("a", "b")),
                       boundary = "clip"), "reflective")
  expect_error(trustL1(obj, c(a = 0, b = 0), mu = c(a = 0),
                       fuse = list(c("a", "b"))), "both in mu")
})


test_that("trustL1 maximizing subtracts the penalty", {
  x0 <- c(a = 3, b = -0.2, c = 0.1)
  neg <- function(pars, ...) {
    o <- .quadratic_objfn(x0)(pars, ...)
    o$value <- -o$value; o$gradient <- -o$gradient; o$hessian <- -o$hessian
    o
  }
  fit <- trustL1(neg, parinit = c(a = 0, b = 0, c = 0), mu = 0 * x0, lambda = 1,
                 minimize = FALSE, rinit = 0.5, rmax = 5, iterlim = 200)
  expect_equal(unname(fit$argument), c(2.5, 0, 0), tolerance = 1e-4)
  ref <- trustL1(.quadratic_objfn(x0), parinit = c(a = 0, b = 0, c = 0), mu = 0 * x0,
                 lambda = 1, rinit = 0.5, rmax = 5, iterlim = 200)
  expect_equal(fit$value, -ref$value, tolerance = 1e-6)
})


test_that("trustL1 takes the tolerances in tolControl and refuses them flat", {
  obj <- .quadratic_objfn(c(a = 1, b = 0.5))
  init <- c(a = 0, b = 0)
  expect_error(trustL1(obj, init, gtol = 1e-8), "tolControl\\$gtol")
  expect_error(trustL1(obj, init, theta.max = 0.9), "stepControl\\$thetaMax")
  expect_error(trustL1(obj, init, tolControl = list(ctol = 1)), "unknown entries ctol")
  fit <- trustL1(obj, init, lambda = 0.1, tolControl = list(gtol = 1e-10),
                 stepControl = list(thetaMax = 0.9))
  expect_true(fit$converged)
})


test_that("trustL1 stops for convergence only where the full Newton step predicts at most mtol", {
  # a coordinate held at the kink and a flat free one, minimum 0 at (0, 9)
  obj <- function(p, ...) {
    h <- c(x = 1, y = 1e-3)
    d <- p - c(x = 0, y = 9)
    H <- diag(unname(h))
    dimnames(H) <- list(names(p), names(p))
    objlist(value = 0.5 * sum(h * d^2), gradient = h * d, hessian = H)
  }
  tol <- list(ftol = 1e-4, mtol = 1e-4, gtol = 1e-2)
  grad <- trustL1(obj, c(x = 0, y = 0), mu = c(x = 0), lambda = 0.1, rinit = 1, rmax = 10,
                  tolControl = tol)
  expect_identical(grad$stopReason, "gradient")
  expect_identical(unname(grad$argument[["x"]]), 0)
  expect_lte(grad$value, 1e-4)
  val <- trustL1(obj, c(x = 0, y = 0), mu = c(x = 0), lambda = 0.1, rinit = 0.05, rmax = 0.05,
                 tolControl = modifyList(tol, list(ftol = 1e-3, gtol = 1e-6)))
  expect_identical(val$stopReason, "fvalue")
  expect_lte(val$value, 1e-4)
})

test_that("trustL1 ends on a stop held back by the Newton decrease when the radius collapses", {
  # the gradient predicts a decrease along the flat direction that the noisy value does not follow
  h <- c(1, 1e-4)
  obj <- function(p, ...) {
    objlist(0.5 * sum(h * p^2) + 5e-5 * sin(1e6 * p[["y"]]), h * p + c(0, 0.01),
            diag(h) |> `dimnames<-`(list(names(p), names(p))))
  }
  fit <- trustL1(obj, c(x = 1, y = 1), mu = c(x = 0, y = 0), lambda = 0, rinit = 1, rmax = 10,
                 iterlim = 200, tolControl = list(ftol = 1e-4, mtol = 1e-4, gtol = 1e-6, rmin = 1e-3))
  expect_true(fit$converged)
  expect_identical(fit$stopReason, "fvalue")
})

test_that("trustL1 ends on stagnation unconverged where the Newton decrease exceeds mtol", {
  # a value that follows no step: the gradient predicts a decrease the run never finds
  flat <- function(p, ...) objlist(value = 1, gradient = c(x = 1, y = 1),
                                   hessian = diag(2) |> `dimnames<-`(list(names(p), names(p))))
  fit <- trustL1(flat, c(x = 0, y = 0), mu = c(x = 5), lambda = 0.1, iterlim = 500)
  expect_identical(fit$stopReason, "stagnation")
  expect_false(fit$converged)
})
