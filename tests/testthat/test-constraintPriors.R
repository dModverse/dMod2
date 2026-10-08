# The soft constraints beyond L2, each against its own density on the -2 log
# scale with closed-form derivatives, and composed with a parameter trafo.

.priorCases <- function() {
  list(
    list(name = "constraintL1",
         fn   = constraintL1(c(k = 3), sigma = 5),
         ld   = function(x) -log(2 * 5) - abs(x - 3) / 5,
         d1   = function(x) 2 * sign(x - 3) / 5,
         d2   = function(x) 0),
    list(name = "constraintCauchy",
         fn   = constraintCauchy(c(k = 3), sigma = 5),
         ld   = function(x) stats::dcauchy(x, 3, 5, log = TRUE),
         d1   = function(x) 4 * (x - 3) / (25 + (x - 3)^2),
         d2   = function(x) 4 * (25 - (x - 3)^2) / (25 + (x - 3)^2)^2),
    list(name = "constraintGamma",
         fn   = constraintGamma(c(k = 3), scale = 5),
         ld   = function(x) stats::dgamma(x, shape = 3, scale = 5, log = TRUE),
         d1   = function(x) -2 * (2 / x - 1 / 5),
         d2   = function(x) 4 / x^2),
    list(name = "constraintExponential",
         fn   = constraintExponential(c(k = 3)),
         ld   = function(x) stats::dexp(x, rate = 1 / 3, log = TRUE),
         d1   = function(x) 2 / 3,
         d2   = function(x) 0),
    list(name = "constraintChisq",
         fn   = constraintChisq(c(k = 4)),
         ld   = function(x) stats::dchisq(x, 4, log = TRUE),
         d1   = function(x) -2 * (1 / x - 1 / 2),
         d2   = function(x) 2 / x^2),
    list(name = "constraintRayleigh",
         fn   = constraintRayleigh(c(k = 3)),
         ld   = function(x) log(x) - 2 * log(3) - x^2 / (2 * 3^2),
         d1   = function(x) -2 * (1 / x - x / 9),
         d2   = function(x) 2 * (1 / x^2 + 1 / 9)))
}


test_that("each constraint is -2 log of its density", {
  for (case in .priorCases())
    expect_equal(unname(case$fn(pars = c(k = 5))$value), -2 * case$ld(5),
                 tolerance = 1e-12, info = case$name)
})


test_that("gradient and Hessian are the derivatives of -2 log density", {
  for (case in .priorCases()) {
    o <- case$fn(pars = c(k = 5))
    expect_equal(unname(o$gradient[["k"]]), case$d1(5), tolerance = 1e-12,
                 info = case$name)
    expect_equal(unname(o$hessian[1, 1]), case$d2(5), tolerance = 1e-12,
                 info = case$name)
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


# An R-level log trafo k = exp(logk), with its own Jacobian and Hessian.
.priorLogTrafo <- function() {
  p2p <- function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE) {
    k <- exp(pars[["logk"]])
    as.parvec(c(k = k),
              deriv  = if (deriv) matrix(k, 1, 1, dimnames = list("k", "logk")) else FALSE,
              deriv2 = if (deriv2) array(k, c(1, 1, 1),
                                         dimnames = list("k", "logk", "logk")) else FALSE)
  }
  parfn(p2p, "logk", NULL)
}

test_that("composing with a parfn applies the chain rule", {
  case <- .priorCases()[[2]]
  comp <- case$fn * .priorLogTrafo()
  k <- 5
  o <- comp(pars = c(logk = log(k)), deriv = TRUE, deriv2 = TRUE)
  gn <- comp(pars = c(logk = log(k)), deriv = TRUE)

  expect_equal(unname(o$value), unname(case$fn(pars = c(k = k))$value),
               tolerance = 1e-12)
  expect_equal(unname(o$gradient[["logk"]]), case$d1(k) * k, tolerance = 1e-12)
  # deriv2 is opt-in: only then does the second-order trafo term enter.
  expect_equal(unname(gn$hessian[1, 1]), case$d2(k) * k^2, tolerance = 1e-12)
  expect_equal(unname(o$hessian[1, 1]), case$d2(k) * k^2 + case$d1(k) * k,
               tolerance = 1e-12)
})
