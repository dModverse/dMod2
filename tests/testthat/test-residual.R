# res() and datapointL2() on hand-built predictions.
# End-to-end coverage is in test-normL2.R.


# Hand-built decay prediction with analytical sensitivities.
.make_prdframe <- function(times = c(0, 1, 2)) {
  prdf <- matrix(c(times, exp(-0.5 * times)),
                 nrow = length(times), ncol = 2,
                 dimnames = list(NULL, c("time", "A")))
  d_arr <- array(0, c(length(times), 1, 2),
                 dimnames = list(NULL, "A", c("A_par", "k_par")))
  d_arr[, "A", "A_par"] <- exp(-0.5 * times)
  d_arr[, "A", "k_par"] <- -times * 1.0 * exp(-0.5 * times)
  attr(prdf, "deriv") <- d_arr
  prdf
}


# ---- res ----------------------------------------------------------------

test_that("res computes weighted.residual = (pred - value) / sigma per row", {
  prdf <- .make_prdframe()
  dat <- data.frame(time = c(1, 2), name = "A",
                    value = c(0.55, 0.30), sigma = c(0.1, 0.2),
                    lloq = c(-Inf, -Inf),
                    stringsAsFactors = FALSE)
  r <- res(dat, prdf)
  pred <- c(exp(-0.5), exp(-1.0))
  expect_equal(r$prediction, pred, tolerance = 1e-12)
  expect_equal(r$weighted.residual, (pred - dat$value) / dat$sigma,
               tolerance = 1e-12)
  expect_equal(r$weighted.0, pred / dat$sigma, tolerance = 1e-12)
  expect_true(all(r$bloq == FALSE))
})


test_that("res sets value = pmax(value_raw, lloq) and bloq mask correctly", {
  prdf <- .make_prdframe()
  dat <- data.frame(time = c(1, 2), name = "A",
                    value = c(0.55, 0.05),  # 0.05 < lloq -> censored
                    sigma = c(0.1, 0.1),
                    lloq  = c(-Inf, 0.10),
                    stringsAsFactors = FALSE)
  r <- res(dat, prdf)
  expect_equal(r$value, c(0.55, 0.10), tolerance = 1e-12)
  expect_equal(r$bloq, c(FALSE, TRUE))
  expect_equal(r$weighted.residual[2], (exp(-1) - 0.10) / 0.1,
               tolerance = 1e-12)
})


test_that("res fills sigma from errmodel matrix where data$sigma is NA", {
  prdf <- .make_prdframe()
  err <- matrix(c(0, 1, 2,
                  0.2, 0.2, 0.2),
                nrow = 3, ncol = 2,
                dimnames = list(NULL, c("time", "A")))

  dat <- data.frame(time = c(1, 2), name = "A",
                    value = c(0.55, 0.30),
                    sigma = c(NA, 0.1),
                    lloq  = c(-Inf, -Inf),
                    stringsAsFactors = FALSE)
  r <- res(dat, prdf, err = err)
  expect_equal(r$sigma, c(0.2, 0.1), tolerance = 1e-12)
  pred <- c(exp(-0.5), exp(-1.0))
  expect_equal(r$weighted.residual,
               (pred - dat$value) / r$sigma, tolerance = 1e-12)
})


test_that("res 'deriv' attribute has shape [n_rows, n_params] with expected names", {
  prdf <- .make_prdframe()
  dat <- data.frame(time = c(1, 2), name = "A",
                    value = c(0.55, 0.30), sigma = c(0.1, 0.1),
                    lloq = c(-Inf, -Inf),
                    stringsAsFactors = FALSE)
  r <- res(dat, prdf)
  d <- attr(r, "deriv")
  expect_equal(dim(d), c(2, 2))
  expect_equal(colnames(d), c("A_par", "k_par"))
  expect_equal(d[, "A_par"], c(exp(-0.5), exp(-1.0)), tolerance = 1e-12)
  expect_equal(d[, "k_par"], c(-1 * exp(-0.5), -2 * exp(-1.0)),
               tolerance = 1e-12)
})


test_that("res errors if an observable in data is missing from the prediction", {
  prdf <- .make_prdframe()
  dat <- data.frame(time = 1, name = "B",  # B does not exist in prdf
                    value = 0.5, sigma = 0.1, lloq = -Inf,
                    stringsAsFactors = FALSE)
  expect_error(res(dat, prdf), regexp = "Observable not found")
})


# ---- datapointL2 --------------------------------------------------------

# datapointL2 penalises (prediction[name, t] - target_par) / sigma in L2.
# It reads prediction from the shared `env` set by an upstream normL2.

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("cppDE")
  testthat::skip_on_cran()
}


test_that("datapointL2 value equals ((pred - target) / sigma)^2 at a known point", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.1)

  obj_main <- normL2(data, bench$prd_id)
  obj_val  <- datapointL2(name = "y", time = 2.0, parameter = "newpoint",
                          sigma = 0.05, condition = "C1")
  obj <- obj_main + obj_val

  pars <- c(bench$outerpars_id, newpoint = 0.5)
  o    <- obj(pars)
  o_main <- obj_main(bench$outerpars_id)

  pred_at_t <- unname(bench$prd_id(times = c(0, 2), pars = bench$outerpars_id,
                                   deriv = FALSE)$C1[2, "y"])
  expected_val <- ((pred_at_t - pars[["newpoint"]]) / 0.05)^2
  expect_equal(unname(o$value - o_main$value), expected_val, tolerance = 1e-3)
})


test_that("normL2 + datapointL2 value equals the sum of the parts", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.1)

  obj_main <- normL2(data, bench$prd_id)
  obj_val  <- datapointL2(name = "y", time = 5.0, parameter = "vtarget",
                          sigma = 0.1, condition = "C1")
  obj <- obj_main + obj_val
  pars <- c(bench$outerpars_id, vtarget = 0.08)

  # +.objfn evaluates main first, then datapointL2 reuses env$prediction.
  v_main <- obj_main(bench$outerpars_id)
  env_main <- attr(v_main, "env")
  v_pt   <- obj_val(pars = pars, env = env_main)

  v_total <- obj(pars)
  expect_equal(v_total$value, v_main$value + v_pt$value, tolerance = 1e-10)
})


test_that("datapointL2 contribution to the combined gradient follows the closed form", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.1)

  obj_main <- normL2(data, bench$prd_id)
  obj_val  <- datapointL2(name = "y", time = 3.0, parameter = "target",
                          sigma = 0.05, condition = "C1")
  obj <- obj_main + obj_val
  pars <- c(A = 1.0, k = 0.5, target = 0.2)

  t_pt <- 3.0; sigma_pt <- 0.05
  pred_pt <- pars[["A"]] * exp(-pars[["k"]] * t_pt)
  r <- pred_pt - pars[["target"]]
  contrib <- c(A      = 2 * r * exp(-pars[["k"]] * t_pt) / sigma_pt^2,
               k      = 2 * r * (-t_pt * pars[["A"]] *
                                   exp(-pars[["k"]] * t_pt)) / sigma_pt^2,
               target = -2 * r / sigma_pt^2)

  g_combined <- obj(pars)$gradient
  g_main     <- obj_main(pars[c("A", "k")])$gradient

  expect_equal(unname(g_combined[["A"]]),
               unname(g_main[["A"]] + contrib[["A"]]), tolerance = 1e-3)
  expect_equal(unname(g_combined[["k"]]),
               unname(g_main[["k"]] + contrib[["k"]]), tolerance = 1e-3)
  expect_equal(unname(g_combined[["target"]]), unname(contrib[["target"]]),
               tolerance = 1e-3)
})
