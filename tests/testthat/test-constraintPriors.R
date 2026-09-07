# ============================================================================
# Behavioral tests for the soft constraints beyond L2: constraintL1,
# constraintCauchy, constraintGamma, constraintExponential, constraintChisq
# and constraintRayleigh.
#
# Each is checked against R's own density on the -2 log scale, and its
# derivatives against finite differences. Composition with a parfn lives at
# the bottom.
# ============================================================================

.priorCases <- function() {
  list(
    list(name = "constraintL1",
         fn   = constraintL1(c(k = 3), sigma = 5),
         ld   = function(x) -log(2 * 5) - abs(x - 3) / 5),
    list(name = "constraintCauchy",
         fn   = constraintCauchy(c(k = 3), sigma = 5),
         ld   = function(x) stats::dcauchy(x, 3, 5, log = TRUE)),
    list(name = "constraintGamma",
         fn   = constraintGamma(c(k = 3), scale = 5),
         ld   = function(x) stats::dgamma(x, shape = 3, scale = 5, log = TRUE)),
    list(name = "constraintExponential",
         fn   = constraintExponential(c(k = 3)),
         ld   = function(x) stats::dexp(x, rate = 1 / 3, log = TRUE)),
    list(name = "constraintChisq",
         fn   = constraintChisq(c(k = 4)),
         ld   = function(x) stats::dchisq(x, 4, log = TRUE)),
    list(name = "constraintRayleigh",
         fn   = constraintRayleigh(c(k = 3)),
         ld   = function(x) log(x) - 2 * log(3) - x^2 / (2 * 3^2)))
}


test_that("each constraint is -2 log of its density", {
  for (case in .priorCases())
    expect_equal(unname(case$fn(pars = c(k = 5))$value), -2 * case$ld(5),
                 tolerance = 1e-12, info = case$name)
})


test_that("gradient and Hessian match finite differences", {
  h <- 1e-5
  for (case in .priorCases()) {
    f <- case$fn
    o <- f(pars = c(k = 5))
    up <- f(pars = c(k = 5 + h))$value
    dn <- f(pars = c(k = 5 - h))$value
    expect_equal(unname(o$gradient[["k"]]), (up - dn) / (2 * h),
                 tolerance = 1e-4, info = case$name)
    expect_equal(unname(o$hessian[1, 1]), (up - 2 * o$value + dn) / h^2,
                 tolerance = 1e-3, info = case$name)
  }
})


test_that("the positive families are Inf outside their support", {
  for (nm in c("constraintGamma", "constraintChisq", "constraintRayleigh"))
    expect_equal(unname(.priorCases()[[
      which(vapply(.priorCases(), `[[`, "", "name") == nm)]]$fn(
        pars = c(k = -1))$value), Inf, info = nm)
  expect_equal(unname(constraintExponential(c(k = 3))(pars = c(k = -1))$value), Inf)
})


test_that("constraints sum over their parameters and skip unnamed ones", {
  obj <- constraintCauchy(c(a = 0, b = 1), sigma = c(a = 2, b = 3))
  one <- constraintCauchy(c(a = 0), sigma = 2)(pars = c(a = 0.5))$value
  two <- constraintCauchy(c(b = 1), sigma = 3)(pars = c(b = 4))$value
  o   <- obj(pars = c(a = 0.5, b = 4, c = 9))
  expect_equal(unname(o$value), unname(one + two), tolerance = 1e-12)
  expect_equal(unname(o$gradient[["c"]]), 0)
})


test_that("a fixed parameter contributes to the value but not the gradient", {
  obj <- constraintCauchy(c(k = 3, m = 0), sigma = 5)
  o   <- obj(pars = c(k = 5), fixed = c(m = 1))
  expect_equal(names(o$gradient), "k")
  expect_equal(unname(o$value),
               unname(constraintCauchy(c(k = 3), sigma = 5)(pars = c(k = 5))$value +
                      constraintCauchy(c(m = 0), sigma = 5)(pars = c(m = 1))$value),
               tolerance = 1e-12)
})


test_that("an estimated scale is rejected with a pointer to constraintL2", {
  expect_error(constraintCauchy(c(k = 3), sigma = "s"), "constraintL2")
})


test_that("a duplicated parameter name is rejected", {
  expect_error(constraintCauchy(c(k = 1, k = 2), sigma = 1), "Duplicated")
})


test_that("composing with a parfn carries the chain rule", {
  withr::local_dir(tempdir())
  p <- P(eqnvec(k = "exp(logk)"), condition = NULL, compile = TRUE,
         modelname = "constraintChain", deriv2 = TRUE)
  f    <- constraintCauchy(c(k = 3), sigma = 5)
  comp <- f * p
  lk   <- log(5)
  h    <- 1e-6

  o <- comp(pars = c(logk = lk), deriv = TRUE, deriv2 = TRUE)
  up <- comp(pars = c(logk = lk + h), deriv = FALSE)$value
  dn <- comp(pars = c(logk = lk - h), deriv = FALSE)$value

  expect_equal(unname(o$value), unname(f(pars = c(k = exp(lk)))$value),
               tolerance = 1e-12)
  expect_equal(unname(o$gradient[["logk"]]), (up - dn) / (2 * h),
               tolerance = 1e-5)
  # deriv2 is opt-in: only then does the second-order trafo term enter.
  expect_equal(unname(o$hessian[1, 1]),
               (up - 2 * comp(pars = c(logk = lk), deriv = FALSE)$value + dn) / h^2,
               tolerance = 1e-3)
})
