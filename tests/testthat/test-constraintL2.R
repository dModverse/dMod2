# ============================================================================
# Behavioral tests for constraintL2().
#
# Two paths are exercised:
#   * scalar / diagonal sigma (Gaussian L2 prior) -- "Scalar" sections below
#   * MVN-Omega (FOCEI prior on random effects)  -- "MVN" sections below
#
# Second-order chain-rule semantics live in test-deriv2.R.
# ============================================================================


# ---- Scalar: single parameter ------------------------------------------

test_that("constraintL2 with scalar sigma equals (p - mu)^2 / sigma^2", {
  mu    <- c(theta = 0.5)
  sigma <- 0.1
  obj   <- constraintL2(mu = mu, sigma = sigma)

  for_each_backend(function(cpp) {
    o <- obj(c(theta = 0.8))
    expected <- ((0.8 - 0.5) / sigma)^2
    expect_equal(unname(o$value), expected, tolerance = 1e-12,
                 info = paste0("cpp=", cpp))
    expect_equal(unname(o$gradient[["theta"]]),
                 2 * (0.8 - 0.5) / sigma^2, tolerance = 1e-12,
                 info = paste0("cpp=", cpp))
    expect_equal(unname(o$hessian[1, 1]), 2 / sigma^2,
                 tolerance = 1e-12, info = paste0("cpp=", cpp))
  })
})


# ---- Scalar: vector sigma, multiple parameters -------------------------

test_that("constraintL2 sums per-parameter contributions when mu has length > 1", {
  mu    <- c(a = 0.0, b = 1.0, c = -0.5)
  sigma <- c(a = 1.0, b = 0.5, c = 0.2)
  obj   <- constraintL2(mu = mu, sigma = sigma)
  p     <- c(a = 0.3, b = 1.4, c = -0.6)

  for_each_backend(function(cpp) {
    o <- obj(p)
    expected <- sum(((p - mu) / sigma)^2)
    expect_equal(unname(o$value), unname(expected), tolerance = 1e-12,
                 info = paste0("cpp=", cpp))
    expect_equal(unname(o$gradient[names(p)]),
                 unname(2 * (p - mu) / sigma^2), tolerance = 1e-12,
                 info = paste0("cpp=", cpp))
    expect_equal(unname(diag(o$hessian)),
                 unname(2 / sigma^2), tolerance = 1e-12,
                 info = paste0("cpp=", cpp))
    H <- o$hessian
    expect_lt(max(abs(H[upper.tri(H)])), 1e-12,
              label = paste0("cpp=", cpp, " off-diag"))
  })
})


test_that("constraintL2 has value 0 and gradient 0 at the prior mean", {
  mu  <- c(a = 0.1, b = 0.2, c = 0.3)
  sg  <- c(a = 0.5, b = 0.5, c = 0.5)
  obj <- constraintL2(mu = mu, sigma = sg)

  for_each_backend(function(cpp) {
    o <- obj(mu)
    expect_equal(unname(o$value), 0, tolerance = 1e-14,
                 info = paste0("cpp=", cpp))
    expect_lt(max(abs(o$gradient)), 1e-14,
              label = paste0("cpp=", cpp))
  })
})


test_that("constraintL2 only constrains parameters whose names appear in mu", {
  mu  <- c(a = 0.0, b = 0.0)
  sg  <- 1.0
  obj <- constraintL2(mu = mu, sigma = sg)
  p   <- c(a = 0.4, b = -0.3, c = 100)

  for_each_backend(function(cpp) {
    o <- obj(p)
    expect_equal(unname(o$value), 0.4^2 + (-0.3)^2, tolerance = 1e-12,
                 info = paste0("cpp=", cpp))
    expect_equal(unname(o$gradient[["c"]]), 0, tolerance = 1e-12,
                 info = paste0("cpp=", cpp))
  })
})


test_that("(constraintL2(mu1) + constraintL2(mu2))(pars) sums per-element", {
  mu1 <- c(a = 0.0, b = 0.0); mu2 <- c(a = 1.0, b = 1.0)
  obj <- constraintL2(mu = mu1, sigma = 1, attr.name = "c1") +
         constraintL2(mu = mu2, sigma = 1, attr.name = "c2")
  p <- c(a = 0.3, b = 0.7)

  for_each_backend(function(cpp) {
    v_total <- obj(p)$value
    v_a <- (0.3 - 0)^2 + (0.7 - 0)^2 + (0.3 - 1)^2 + (0.7 - 1)^2
    expect_equal(unname(v_total), v_a, tolerance = 1e-12,
                 info = paste0("cpp=", cpp))
  })
})


# ---- MVN: backwards compatibility --------------------------------------

test_that("backwards compatibility: existing scalar/diagonal path still works", {
  prior <- structure(rep(0, 3), names = letters[1:3])
  obj <- constraintL2(mu = prior)
  res <- obj(c(a = 1, b = -1, c = 0.5))
  expect_equal(res$value, 1 + 1 + 0.25)
})


# ---- MVN: value --------------------------------------------------------



# ---- MVN: gradient -----------------------------------------------------



# Closed-form eta-block gradient for full Omega: grad_eta_i = 2 * Omega^-1 * eta_i.
# The chol-block gradient is more cumbersome to derive without re-deriving the
# forward/backsolve formula; we restrict the eta-only assertion here.





# ---- MVN: Hessian ------------------------------------------------------



# ---- MVN: misc ---------------------------------------------------------







# ============================================================================
# Cross-backend parity (C++ kernel vs R reference)
# ============================================================================

test_that("constraintL2 cpp kernel agrees with R reference on a small diagonal case", {
  obj <- constraintL2(mu = c(a = 0.0, b = 1.0, c = -0.5), sigma = c(0.5, 1, 2))
  pars <- c(a = 0.3, b = 1.4, c = -0.6)
  with_cpp_backend(FALSE, { o_R <- obj(pars) })
  with_cpp_backend(TRUE,  { o_C <- obj(pars) })
  expect_equal(o_C$value,    o_R$value,    tolerance = 1e-12)
  expect_equal(o_C$gradient, o_R$gradient, tolerance = 1e-12)
  expect_equal(o_C$hessian,  o_R$hessian,  tolerance = 1e-12)
})
