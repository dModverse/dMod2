# Argument checks and return types across the API.

test_that("an argument a function does not take is an error", {
  expect_error(constraintGamma(c(a = 1), bogus = 1), "unused argument")
  expect_error(normL2(NULL, NULL, attr.name = "p"), "unused argument")
  expect_error(dMod2:::.pimplPTC(list(rtol = 1e-9)), "unknown controlsPTC entry rtol")
  expect_error(trust(constraintL2(c(a = 0)), c(a = 1), ftol = 1e-3),
               "tolControl\\$ftol")
})

test_that("the snake-case names of the symmetry controls are derived", {
  expect_identical(dMod2:::.symControlSnake(c("perPrimeCap", "minSupportCandCap",
                                              "degreeCap")),
                   c("perprime_cap", "minsupport_cand_cap", "degree_cap"))
})

test_that("rref names its result", {
  A <- rbind(c(1, 2, 3), c(2, 4, 7))
  r <- rref(A)
  expect_named(r, c("rref", "pivots"))
  expect_identical(r[[1]], r$rref)
})

test_that("repar takes the transformation first and refuses the other order", {
  tr <- c(a = "a", b = "b")
  expect_identical(repar(tr, "x ~ exp(x)", x = "a"), c(a = "exp(a)", b = "b"))
  expect_error(repar("x ~ exp(x)", tr, x = "a"), "repar\\(trafo, expr\\)")
  expect_equal(unclass(repar(NULL, "x ~ x", x = "a")), c(a = "a"))
})

test_that("a transformation started from NULL keeps the right-hand side", {
  expect_equal(unclass(define(NULL, "x ~ exp(x)", x = c("a", "b"))),
               c(a = "exp(a)", b = "exp(b)"))
  expect_equal(unclass(insert(NULL, "x ~ exp(x)", x = "a")), c(a = "exp(a)"))
  expect_s3_class(define(NULL, "x ~ x", x = "a"), "eqnvec")
})

test_that("trust returns a trustfit; vcov and profile are stats methods", {
  obj <- constraintL2(mu = c(a = 0.5, b = -0.3), sigma = 1)
  fit <- trust(obj, c(a = 0, b = 0), rinit = 1, rmax = 10)
  expect_s3_class(fit, "trustfit")
  expect_true(is.list(fit))
  expect_identical(vcov(fit), stats::vcov(fit))
  expect_false(exists("profile", envir = asNamespace("dMod2"), inherits = FALSE))

  expect_s3_class(profile(obj, c(a = 0.5, b = -0.3), "a", limits = c(-1, 1),
                          cores = 1), "parframe")
  plain <- function(p, ...) obj(p, ...)
  expect_s3_class(profile(plain, c(a = 0.5, b = -0.3), "a", limits = c(-1, 1),
                          cores = 1), "parframe")

  fits <- mstrust(obj, c(a = 0, b = 0), fits = 2, cores = 1, sd = 1,
                  samplefun = stats::rnorm)
  expect_s3_class(fits[[1]], "trustfit")
})

test_that("msParframe leaves the global RNG alone", {
  set.seed(7)
  before <- .Random.seed
  msParframe(c(a = 0, b = 1), n = 3, seed = 1)
  expect_identical(.Random.seed, before)
})

test_that("plotFluxes takes its legend title", {
  times <- 0:5
  grid <- data.frame(name = "A", time = times, row.names = paste0("A", times))
  x <- Xd(grid)
  pars <- structure(exp(-times / 2), names = getParameters(x))
  p <- plotFluxes(pars, x, seq(0, 5, 0.5), c(prod = "0.2"), legendTitle = "F")
  expect_s3_class(p, "ggplot")
  expect_error(symmetryDetection(eqnvec(x = "-k*x")), "give the observables")
})
