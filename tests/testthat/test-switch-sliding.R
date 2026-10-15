# A state that slides along a switching surface, through the whole chain.

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("cppDE")
  testthat::skip_on_cran()
}

fx_register(extra = function(bench) {
  xfn <- Xs(odemodel(eqnvec(x = "kin - v*(x > 0)"), modelname = "slide_ode",
                     compile = FALSE))
  list(slide_x = xfn,
       slide_g = Y(c(xo = "x"), f = xfn, attachInput = FALSE,
                   modelname = "slide_obs", compile = FALSE),
       slide_p = P(eqnvec(x = "x", kin = "kin", v = "v"), condition = "C1",
                   modelname = "slide_p", compile = FALSE))
})

test_that("normL2 of a sliding model has the closed-form gradient and Hessian", {
  skip_if_no_compile()
  x <- fx_extra()
  prd <- x$slide_g * x$slide_x * x$slide_p
  tt <- c(0.5, 1, 2, 3, 4)
  sigma <- 0.05
  value <- c(0.72, 0.38, 0.01, -0.02, 0)
  data <- as.datalist(data.frame(name = "xo", time = tt, value = value,
                                 sigma = sigma, condition = "C1"),
                      splitBy = "condition")
  p <- c(x = 1, kin = 0.2, v = 0.8)
  on <- tt < 1 / 0.6
  res <- (ifelse(on, 1 - 0.6 * tt, 0) - value) / sigma
  J <- cbind(x = on, kin = tt * on, v = -tt * on) / sigma
  out <- normL2(data, prd)(p)
  expect_equal(unname(out$gradient[names(p)]), unname(2 * colSums(res * J)),
               tolerance = 1e-6)
  expect_equal(unname(out$hessian[names(p), names(p)]), unname(2 * crossprod(J)),
               tolerance = 1e-6)
})
