# Behavioral tests for Y() (observation function).
#
# Verifies:
#   * value: observable g(states) evaluates correctly
#   * composition: (Y * Xs)(...) equals Y applied to Xs output
#   * derivMode: "reverse" and "forward" builds agree on the value
#   * attach.input: pass-through of inputs alongside outputs
#   * gradient: analytic chain rule on y = A^2 (no numDeriv)
#
# Second-order chain rule is covered by test-deriv2-Y.R.

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("cppDE")
  testthat::skip_on_cran()
}


# Every observation function the file needs, compiled into one shared object on
# first use. Prediction and trafo of the decay chain come from the shared fixture.
.y_fx <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    bench <- fx_decay_compiled()
    oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

    g_sq <- Y(c(y = "A^2"), f = bench$xfn, condition = NULL,
              attach.input = FALSE, modelname = "test_Y_sq", compile = FALSE)
    g_rev <- Y(c(y = "A^2"), f = bench$xfn, condition = NULL,
               attach.input = FALSE, derivMode = "reverse",
               modelname = "test_Y_dm_rev", compile = FALSE)
    g_dual <- Y(c(y = "A^2"), f = bench$xfn, condition = NULL,
                attach.input = FALSE, derivMode = "forward",
                modelname = "test_Y_dm_dual", compile = FALSE)
    g_attach <- Y(c(y = "A"), f = bench$xfn, condition = NULL,
                  attach.input = TRUE, modelname = "test_Y_attach", compile = FALSE)

    x_np <- Xs(odemodel(as.eqnvec(c(A = "-k*A")), modelname = "noparam_y_ode",
                        compile = FALSE, backend = "cppDE"))
    g_np <- Y(c(y1 = "1.0"), f = NULL, states = c("A"),
              parameters = character(0), derivMode = "forward",
              compile = FALSE, modelname = "noparam_y_obs")

    compile(g_sq, g_rev, g_dual, g_attach, x_np, g_np,
            output = "test_Y_all", cores = 4L)
    cache <<- list(bench = bench, g_sq = g_sq, g_rev = g_rev, g_dual = g_dual,
                   g_attach = g_attach, x_np = x_np, g_np = g_np)
    cache
  }
})


## ---- Value: linear observable ------------------------------------------

test_that("Y(y = A) on linear decay matches A(t) directly", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  prd_obs <- bench$gfn * bench$xfn * bench$pfn_id

  times <- c(0, 1, 2, 5)
  pars  <- c(A = 1.7, k = 0.42)
  out <- prd_obs(times = times, pars = pars, deriv = FALSE)

  closed <- pars[["A"]] * exp(-pars[["k"]] * times)
  expect_equal(out$C1[, "y"], closed, tolerance = 1e-5)
})


## ---- Value: nonlinear observable ---------------------------------------

test_that("Y(y = A^2) evaluates the closed-form (A(t))^2", {
  skip_if_no_compile()
  fx <- .y_fx()
  bench <- fx$bench
  prd_sq <- fx$g_sq * bench$xfn * bench$pfn_id

  times <- c(0, 1, 2, 5)
  pars  <- c(A = 1.4, k = 0.3)
  out <- prd_sq(times = times, pars = pars, deriv = FALSE)

  closed <- (pars[["A"]] * exp(-pars[["k"]] * times))^2
  expect_equal(out$C1[, "y"], closed, tolerance = 1e-5)
})


## ---- derivMode parity --------------------------------------------------

test_that("Y derivMode 'reverse' and 'forward' agree on a nonlinear observable", {
  skip_if_no_compile()
  fx <- .y_fx()
  bench <- fx$bench

  prd_sym  <- fx$g_rev  * bench$xfn * bench$pfn_id
  prd_dual <- fx$g_dual * bench$xfn * bench$pfn_id

  times <- c(0, 1, 2, 5)
  pars  <- c(A = 1.4, k = 0.3)
  o_sym  <- prd_sym (times = times, pars = pars, deriv = TRUE)
  o_dual <- prd_dual(times = times, pars = pars, deriv = TRUE)

  expect_equal(o_sym$C1[, "y"], o_dual$C1[, "y"], tolerance = 1e-12)
  # The reverse build carries no forward entries, the forward one does.
  expect_null(attr(o_sym$C1, "deriv"))
  expect_false(is.null(attr(o_dual$C1, "deriv")))
})


## ---- attach.input ------------------------------------------------------

test_that("Y with attach.input = TRUE returns inputs and outputs", {
  skip_if_no_compile()
  fx <- .y_fx()
  bench <- fx$bench
  prd_in <- fx$g_attach * bench$xfn * bench$pfn_id

  times <- c(0, 1, 2)
  pars  <- c(A = 1.0, k = 0.5)
  out <- prd_in(times = times, pars = pars, deriv = FALSE)

  cn <- colnames(out$C1)
  expect_true("y" %in% cn)
  expect_true("A" %in% cn)
  # When attach.input is TRUE, the observable column should equal the input
  # column for the identity observable.
  expect_equal(out$C1[, "y"], out$C1[, "A"], tolerance = 1e-12)
})


## ---- Gradient: analytic chain rule on y = A^2 --------------------------

test_that("Y gradient on y = A^2 follows the analytic chain rule dy/dtheta = 2 A * dA/dtheta", {
  skip_if_no_compile()
  fx <- .y_fx()
  bench <- fx$bench
  prd_sq <- fx$g_sq * bench$xfn * bench$pfn_id

  times <- c(0, 1, 2, 5)
  pars  <- c(A = 1.0, k = 0.5)
  out <- prd_sq(times = times, pars = pars, deriv = TRUE)
  d <- attr(out$C1, "deriv")  # [time, var, par]

  # Closed form: A(t) = A0 * exp(-k * t), y(t) = A(t)^2.
  #   dy/dA0 = 2 * A * (dA/dA0) = 2 * A0 * exp(-k*t) * exp(-k*t) = 2 * A0 * exp(-2 k t)
  #   dy/dk  = 2 * A * (dA/dk)  = 2 * A0 * exp(-k*t) * (-t * A0 * exp(-k*t))
  #                              = -2 * t * A0^2 * exp(-2 k t)
  A0 <- pars[["A"]]; k <- pars[["k"]]
  ref_dA <- 2 * A0 * exp(-2 * k * times)
  ref_dk <- -2 * times * A0^2 * exp(-2 * k * times)
  expect_equal(d[, "y", "A"], ref_dA, tolerance = 1e-5)
  expect_equal(d[, "y", "k"], ref_dk, tolerance = 1e-5)
})


# ============================================================================
# Edge case: Y with pure-numeric observable (no outer parameters)
# ============================================================================

test_that("Y with pure-numeric observable composes with an Xs prediction", {
  fx <- .y_fx()
  out <- (fx$g_np * fx$x_np)(seq(0, 5, length.out = 3), c(A = 1.0, k = 0.1))
  pred <- out[[1]]
  expect_true(all(pred[, "y1"] == 1.0))
  expect_equal(pred[, "time"], c(0, 2.5, 5))
})
