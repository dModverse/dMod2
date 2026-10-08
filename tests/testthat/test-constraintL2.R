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
  obj <- constraintL2(mu = mu1, sigma = 1, attrName = "c1") +
         constraintL2(mu = mu2, sigma = 1, attrName = "c2")
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
# Composition with a parameter transformation
# ============================================================================

# Log trafos for one and for two conditions, compiled together on first use.
.cl2_p <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)
    tr <- function(k) eqnvec(A = "exp(A_log)", k = paste0("exp(", k, ")"))
    p1 <- P(tr("k_log"), condition = "C1", modelname = "test_cl2_p1",
            compile = FALSE)
    p2 <- P(list(C1 = tr("k_log"), C2 = tr("k2_log")), modelname = "test_cl2_p2",
            compile = FALSE)
    compile(p1, p2, output = "test_cl2_p", cores = test_cores())
    cache <<- list(p1 = p1, p2 = p2)
    cache
  }
})

test_that("constraintL2 * P takes fixed outer parameters", {
  skip_on_cran()
  withr::local_options(dMod.batch.check = TRUE)
  fx <- .cl2_p()
  prior <- constraintL2(mu = c(A = 1, k = 0.5), sigma = 0.5)
  all <- c(A_log = 0.2, k_log = -0.4, k2_log = 0.1)

  for (p in list(fx$p1, fx$p2)) {
    obj <- prior * p
    pars <- all[getParameters(p)]
    free <- setdiff(names(pars), "A_log")
    full <- obj(pars)
    part <- obj(pars[free], fixed = pars["A_log"])
    expect_equal(part$value, full$value)
    expect_equal(part$gradient, full$gradient[free])
    expect_equal(part$hessian, full$hessian[free, free, drop = FALSE])
  }

  fit <- trust(prior * fx$p1, all["k_log"], rinit = 1, rmax = 10,
               fixed = all["A_log"])
  expect_true(fit$converged)
  expect_equal(fit$argument[["k_log"]], log(0.5), tolerance = 1e-6)
})


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

test_that("constraintL2 refuses numbers among sigma parameter names", {
  expect_error(constraintL2(c(a = 0, b = 0), sigma = c(a = "s_a", b = 2)),
               "all numbers or all parameter names")
})
