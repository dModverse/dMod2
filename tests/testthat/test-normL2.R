# ============================================================================
# Behavioral tests for normL2() against mathematical ground truth.
#
# Every value / gradient claim is checked against a closed-form analytic
# formula (Gaussian log-likelihood, decay sensitivities composed by chain
# rule, errmodel propagation). No finite-difference reference is used.
#
# Hessian and second-order chain-rule semantics live in test-deriv2.R.
#
# All blocks run under both objfn backends (R reference and C++ kernel)
# via for_each_backend(); the info= tag distinguishes failures.
# ============================================================================

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("cppDE")
  testthat::skip_on_cran()
}


# Error-model chains (sigma = sigma_y, sigma = srel * y) and a second-condition
# trafo on the shared decay fixture, in one shared object built on first use.
# The trafos pass the error parameters through, so normL2's errmodel call sees them.
.nl2_fx <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    bench <- fx_decay_compiled()
    oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

    e_const <- Y(c(y = "sigma_y"), f = bench$gfn, attach.input = FALSE,
                 condition = "C1", modelname = "nl2_err_const", compile = FALSE)
    p_sig <- P(eqnvec(A = "A", k = "k", sigma_y = "sigma_y"), condition = "C1",
               modelname = "nl2_p_sig", compile = FALSE)
    e_prop <- Y(c(y = "srel * y"), f = bench$gfn, attach.input = FALSE,
                condition = "C1", modelname = "nl2_err_prop", compile = FALSE)
    p_prop <- P(eqnvec(A = "A", k = "k", srel = "srel"), condition = "C1",
                modelname = "nl2_p_prop", compile = FALSE)
    p_C2 <- P(eqnvec(A = "A", k = "k"), condition = "C2",
              modelname = "nl2_p_id", compile = FALSE)
    compile(e_const, p_sig, e_prop, p_prop, p_C2, output = "nl2_all", cores = 4L)

    cache <<- list(
      const = list(prd = bench$gfn * bench$xfn * p_sig,  e = e_const),
      prop  = list(prd = bench$gfn * bench$xfn * p_prop, e = e_prop),
      pfn_C2 = p_C2)
    cache
  }
})


# ---- Basics: value / gradient -------------------------------------------

test_that("normL2 value equals sum(wr^2) + sum(log(2*pi*sigma^2)) at a known point", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.1)

  prd_at <- bench$prd_id(times = data$C1$time, pars = bench$outerpars_id)
  closed <- truth_nll_aloq(prd_at$C1[, "y"], data$C1$value, data$C1$sigma)

  for_each_backend(function(cpp) {
    obj <- normL2(data, bench$prd_id)
    o <- obj(bench$outerpars_id)
    expect_equal(o$value, closed, tolerance = 1e-10,
                 info = paste0("cpp=", cpp))
  })
})


test_that("normL2 value scales as 1/sigma^2 when sigma is rescaled", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()

  c1 <- 0.1; c2 <- 0.2
  data1 <- fx_decay_data(sigma = c1)
  data2 <- data1
  data2$C1$value <- data1$C1$value
  data2$C1$sigma <- c2

  for_each_backend(function(cpp) {
    o1 <- normL2(data1, bench$prd_id)(bench$outerpars_id)
    o2 <- normL2(data2, bench$prd_id)(bench$outerpars_id)
    chi1 <- o1$value - sum(log(2 * pi * data1$C1$sigma^2))
    chi2 <- o2$value - sum(log(2 * pi * data2$C1$sigma^2))
    expect_equal(chi2 * (c2 / c1)^2, chi1, tolerance = 1e-10,
                 info = paste0("cpp=", cpp))
  })
})


test_that("normL2 gradient equals 2 * Jt * (pred - y) / sigma^2 (analytic decay sens)", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.1)
  pars <- bench$outerpars_id + c(A = 0.15, k = -0.1)

  t <- data$C1$time
  obs <- data$C1$value
  sigma_vec <- data$C1$sigma
  pred <- pars[["A"]] * exp(-pars[["k"]] * t)
  wr <- (pred - obs) / sigma_vec
  J_A <- exp(-pars[["k"]] * t)
  J_k <- -t * pars[["A"]] * exp(-pars[["k"]] * t)
  g_ref <- c(A = 2 * sum(wr * J_A / sigma_vec),
             k = 2 * sum(wr * J_k / sigma_vec))

  for_each_backend(function(cpp) {
    obj <- normL2(data, bench$prd_id)
    g_ana <- obj(pars)$gradient
    expect_equal(unname(g_ana[names(g_ref)]), unname(g_ref),
                 tolerance = 1e-4, info = paste0("cpp=", cpp))
  })
})


# ---- Sigma source equivalence -------------------------------------------

test_that("sigma from data column == sigma from errmodel, constant case", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  sigma_const <- 0.1

  data_col <- fx_decay_data(sigma = sigma_const)
  data_em  <- data_col
  data_em$C1$sigma <- NA_real_

  ec <- .nl2_fx()$const
  pars_em <- c(bench$outerpars_id, sigma_y = sigma_const)

  for_each_backend(function(cpp) {
    o_col <- normL2(data_col, ec$prd)(pars_em)
    o_em  <- normL2(data_em,  ec$prd, errmodel = ec$e)(pars_em)
    expect_equal(o_em$value, o_col$value, tolerance = 1e-10,
                 info = paste0("cpp=", cpp))
  })
})


# ---- Multi-condition aggregation ----------------------------------------

test_that("normL2 sums per-condition contributions across two conditions", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()

  prd_multi <- bench$gfn * bench$xfn * (bench$pfn_id + .nl2_fx()$pfn_C2)
  data_multi <- fx_decay_data_multi(
    parslist = list(C1 = c(A = 1.0, k = 0.5), C2 = c(A = 1.0, k = 1.0)),
    sigma = 0.1)

  pars <- c(A = 1.0, k = 0.7)

  for_each_backend(function(cpp) {
    o_joint <- normL2(data_multi, prd_multi)(pars)
    o_C1 <- normL2(data_multi["C1"], prd_multi)(pars)
    o_C2 <- normL2(data_multi["C2"], prd_multi)(pars)
    expect_equal(o_joint$value, o_C1$value + o_C2$value, tolerance = 1e-10,
                 info = paste0("cpp=", cpp))
  })
})


# ---- BLOQ: closed-form value (M3) ---------------------------------------

test_that("normL2 with BLOQ rows adds -2 * sum(log Phi(-wr_bloq)) (M3) over ALOQ value", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()

  data <- fx_decay_data_bloq(sigma = 0.05, lloq = 0.1,
                             times = seq(0, 10, by = 1))
  pars <- bench$outerpars_id

  prd_at <- bench$prd_id(times = data$C1$time, pars = pars,
                         deriv = FALSE)$C1
  pred <- prd_at[, "y"]
  sigma_vec <- data$C1$sigma
  lloq_vec  <- data$C1$lloq
  val_post  <- pmax(data$C1$value, lloq_vec)
  is_bloq   <- val_post <= lloq_vec

  closed_aloq <- truth_nll_aloq(pred[!is_bloq], data$C1$value[!is_bloq],
                                sigma_vec[!is_bloq])
  closed_bloq <- truth_nll_bloq_m3(pred[is_bloq], lloq_vec[is_bloq],
                                   sigma_vec[is_bloq])
  expected <- closed_aloq + closed_bloq

  for_each_backend(function(cpp) {
    obj <- normL2(data, bench$prd_id)
    o   <- obj(pars)
    expect_equal(o$value, expected, tolerance = 1e-4,
                 info = paste0("cpp=", cpp))
  })
})


test_that("adding a BLOQ row to data strictly increases the normL2 value", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()

  raw_pars <- bench$outerpars_id
  data_no <- fx_decay_data(pars = raw_pars, sigma = 0.05,
                           times = seq(0, 10, by = 1), seed = 7L)
  data_yes <- data_no
  data_yes$C1$lloq <- 0.1

  for_each_backend(function(cpp) {
    o_no  <- normL2(data_no,  bench$prd_id)(raw_pars)
    o_yes <- normL2(data_yes, bench$prd_id)(raw_pars)
    expect_gt(o_yes$value, o_no$value,
              label = paste0("cpp=", cpp, " bloq monotonicity"))
  })
})


test_that("normL2 gradient on a BLOQ dataset equals analytic ALOQ + M3 BLOQ contributions", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data_bloq(sigma = 0.05, lloq = 0.1,
                              times = seq(0, 10, by = 1))
  pars <- bench$outerpars_id + c(A = 0.2, k = -0.1)

  t <- data$C1$time; obs <- data$C1$value; sigma_vec <- data$C1$sigma
  lloq_vec <- data$C1$lloq
  pred <- pars[["A"]] * exp(-pars[["k"]] * t)
  J_A  <- exp(-pars[["k"]] * t)
  J_k  <- -t * pars[["A"]] * exp(-pars[["k"]] * t)
  val_post <- pmax(obs, lloq_vec)
  is_bloq  <- val_post <= lloq_vec
  wr_aloq <- (pred[!is_bloq] - obs[!is_bloq]) / sigma_vec[!is_bloq]
  s_aloq  <- sigma_vec[!is_bloq]
  grad_aloq_A <- 2 * sum(wr_aloq * J_A[!is_bloq] / s_aloq)
  grad_aloq_k <- 2 * sum(wr_aloq * J_k[!is_bloq] / s_aloq)
  if (any(is_bloq)) {
    wr_b <- (pred[is_bloq] - lloq_vec[is_bloq]) / sigma_vec[is_bloq]
    G    <- exp(stats::dnorm(-wr_b, log = TRUE) -
                  stats::pnorm(-wr_b, log.p = TRUE))
    s_b  <- sigma_vec[is_bloq]
    grad_bloq_A <- 2 * sum(G * J_A[is_bloq] / s_b)
    grad_bloq_k <- 2 * sum(G * J_k[is_bloq] / s_b)
  } else {
    grad_bloq_A <- 0; grad_bloq_k <- 0
  }
  g_ref <- c(A = grad_aloq_A + grad_bloq_A, k = grad_aloq_k + grad_bloq_k)

  for_each_backend(function(cpp) {
    obj <- normL2(data, bench$prd_id)
    g_ana <- obj(pars)$gradient
    expect_equal(unname(g_ana[names(g_ref)]), unname(g_ref),
                 tolerance = 1e-3, info = paste0("cpp=", cpp))
  })
})


test_that("R and C++ backends give the same value on a BLOQ dataset", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data_bloq(sigma = 0.05, lloq = 0.1,
                              times = seq(0, 10, by = 1))
  pars <- bench$outerpars_id

  with_cpp_backend(FALSE, {
    v_R <- normL2(data, bench$prd_id)(pars)$value
  })
  with_cpp_backend(TRUE, {
    v_C <- normL2(data, bench$prd_id)(pars)$value
  })
  expect_equal(v_C, v_R, tolerance = 1e-9)
})


test_that("normL2(opt.BLOQ = ...) selects the BLOQ method on both backends", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data_bloq(sigma = 0.05, lloq = 0.1,
                              times = seq(0, 10, by = 1))
  pars  <- bench$outerpars_id

  prd_at    <- bench$prd_id(times = data$C1$time, pars = pars,
                            deriv = FALSE)$C1
  pred      <- prd_at[, "y"]
  sigma_vec <- data$C1$sigma
  lloq_vec  <- data$C1$lloq
  val_post  <- pmax(data$C1$value, lloq_vec)
  is_bloq   <- val_post <= lloq_vec

  aloq_value <- truth_nll_aloq(pred[!is_bloq], data$C1$value[!is_bloq],
                               sigma_vec[!is_bloq])
  bloq_m3    <- truth_nll_bloq_m3(pred[is_bloq], lloq_vec[is_bloq],
                                  sigma_vec[is_bloq])
  bloq_m4    <- truth_nll_bloq_m4(pred[is_bloq], lloq_vec[is_bloq],
                                  sigma_vec[is_bloq])
  w0_aloq <- pred[!is_bloq] / sigma_vec[!is_bloq]
  m4beal_aloq_correction <- 2 * sum(stats::pnorm(w0_aloq, log.p = TRUE))

  expected <- c(
    M1     = aloq_value,
    M3     = aloq_value + bloq_m3,
    M4NM   = aloq_value + bloq_m4,
    M4BEAL = aloq_value + m4beal_aloq_correction + bloq_m4
  )

  for_each_backend(function(cpp) {
    for (mode in names(expected)) {
      obj <- normL2(data, bench$prd_id, opt.BLOQ = mode)
      o   <- obj(pars)
      expect_equal(o$value, expected[[mode]], tolerance = 1e-3,
                   info = paste0("cpp=", cpp, " mode=", mode))
    }
  })

  for (mode in c("M1", "M3", "M4NM", "M4BEAL")) {
    with_cpp_backend(FALSE, {
      g_R <- normL2(data, bench$prd_id, opt.BLOQ = mode)(pars)$gradient
    })
    with_cpp_backend(TRUE, {
      g_C <- normL2(data, bench$prd_id, opt.BLOQ = mode)(pars)$gradient
    })
    expect_equal(g_C, g_R, tolerance = 1e-9,
                 info = paste0("gradient parity, mode=", mode))
  }
})


test_that("normL2 rejects unknown opt.BLOQ values", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data_bloq(sigma = 0.05, lloq = 0.1,
                              times = seq(0, 10, by = 1))
  # match.arg() rejects the unknown value; assert on the listed choices rather
  # than the "should be one of" prefix, which match.arg translates per locale.
  expect_error(normL2(data, bench$prd_id, opt.BLOQ = "M2"),
               "M4BEAL")
})


# ---- errmodel: proportional sigma ---------------------------------------

test_that("normL2 with sigma = srel*y matches the proportional-error log-likelihood", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  ec <- .nl2_fx()$prop

  pars <- c(A = 1.0, k = 0.5, srel = 0.1)
  data <- fx_decay_data(pars = pars[c("A", "k")], sigma = 0.05)
  data$C1$sigma <- NA_real_

  prd_at <- ec$prd(times = data$C1$time, pars = pars, deriv = FALSE)$C1
  pred <- prd_at[, "y"]
  sigma_pred <- pars[["srel"]] * pred
  expected <- truth_nll_aloq(pred, data$C1$value, sigma_pred)

  for_each_backend(function(cpp) {
    obj <- normL2(data, ec$prd, errmodel = ec$e)
    o   <- obj(pars)
    expect_equal(o$value, expected, tolerance = 1e-4,
                 info = paste0("cpp=", cpp))
  })
})


test_that("normL2 gradient with proportional errmodel follows the analytic closed form", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  ec <- .nl2_fx()$prop

  data <- fx_decay_data(pars = c(A = 1.0, k = 0.5), sigma = 0.05)
  data$C1$sigma <- NA_real_
  pars <- c(A = 1.2, k = 0.45, srel = 0.08)

  t <- data$C1$time; obs <- data$C1$value
  A <- pars[["A"]]; k <- pars[["k"]]; srel <- pars[["srel"]]
  pred <- A * exp(-k * t)
  sigma_vec <- srel * pred
  wr <- (pred - obs) / sigma_vec
  dpred <- cbind(A = exp(-k * t), k = -t * A * exp(-k * t), srel = 0)
  dsigma <- srel * dpred
  dsigma[, "srel"] <- pred
  dlog_sigma <- dsigma / sigma_vec
  dwr <- (dpred * sigma_vec - (pred - obs) * dsigma) / sigma_vec^2
  g_ref <- colSums(2 * wr * dwr) + colSums(2 * dlog_sigma)

  for_each_backend(function(cpp) {
    obj <- normL2(data, ec$prd, errmodel = ec$e)
    g_ana <- obj(pars)$gradient
    expect_equal(unname(g_ana[names(g_ref)]), unname(g_ref),
                 tolerance = 1e-3, info = paste0("cpp=", cpp))
  })
})


test_that("rows with explicit sigma keep it; NA rows fall through to errmodel", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  ec <- .nl2_fx()$prop

  pars <- c(A = 1.0, k = 0.5, srel = 0.1)
  data <- fx_decay_data(pars = pars[c("A", "k")], sigma = 0.05)
  n <- nrow(data$C1)
  data$C1$sigma <- ifelse(seq_len(n) <= n %/% 2, NA_real_, 0.05)

  for_each_backend(function(cpp) {
    obj <- normL2(data, ec$prd, errmodel = ec$e)
    o   <- obj(pars)
    expect_true(is.finite(o$value),
                label = paste0("cpp=", cpp, " mixed sigma finite"))
  })

  data_na <- data; data_na$C1 <- data_na$C1[is.na(data$C1$sigma), ]
  data_ex <- data; data_ex$C1 <- data_ex$C1[!is.na(data$C1$sigma), ]
  # Use a common time grid across all three so the adaptive integrator
  # produces bit-identical predictions at the shared data times.
  all_times <- sort(unique(data$C1$time))
  with_cpp_backend(FALSE, {
    v_full <- normL2(data,    ec$prd, errmodel = ec$e, times = all_times)(pars)$value
    v_na   <- normL2(data_na, ec$prd, errmodel = ec$e, times = all_times)(pars)$value
    v_ex   <- normL2(data_ex, ec$prd, errmodel = ec$e, times = all_times)(pars)$value
  })
  expect_equal(v_full, v_na + v_ex, tolerance = 1e-9)
})


test_that("getParameters(normL2(..., errmodel = ec$e)) includes errmodel pars", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  ec <- .nl2_fx()$prop
  data <- fx_decay_data()
  data$C1$sigma <- NA_real_
  obj <- normL2(data, ec$prd, errmodel = ec$e)
  expect_true("srel" %in% attr(obj, "parameters"))
})


# ---- BLOQ + error model: end-to-end wiring ------------------------------
# The errmodel-derived sigma (and its parameter derivatives) must reach the
# BLOQ partition of the kernel, not only the ALOQ rows. Validate the M3 value
# against the closed form built with sigma = srel * pred on both partitions,
# and the full gradient against finite differences (this exercises dsigma
# propagation into the BLOQ rows).
test_that("normL2 BLOQ M3 + proportional errmodel: value and gradient match", {
  skip_if_no_compile()
  ec <- .nl2_fx()$prop

  pars <- c(A = 1.0, k = 0.5, srel = 0.1)
  data <- fx_decay_data_bloq(pars = pars[c("A", "k")], sigma = 0.05,
                             lloq = 0.1, times = seq(0, 10, by = 1))
  data$C1$sigma <- NA_real_   # both ALOQ and BLOQ rows draw sigma from errmodel

  d   <- data$C1
  obj <- normL2(data, ec$prd, errmodel = ec$e, opt.BLOQ = "M3")
  o <- obj(pars)

  # Use the exact prediction normL2 integrated (its env), so the closed-form
  # reference is not perturbed by a separate integrator-grid call.
  env_pred <- attr(o, "env")$prediction[["C1"]]
  pred <- env_pred[match(d$time, env_pred[, "time"]), "y"]
  sigma_pred <- pars[["srel"]] * pred
  bloq <- d$value <= d$lloq
  expect_true(any(bloq) && any(!bloq))   # fixture has both partitions

  expected <- truth_nll_aloq(pred[!bloq], d$value[!bloq], sigma_pred[!bloq]) +
              truth_nll_bloq_m3(pred[bloq], d$lloq[bloq], sigma_pred[bloq])
  expect_equal(o$value, expected, tolerance = 1e-6)

  skip_if_not_installed("numDeriv")
  g_fd <- numDeriv::grad(
    function(p) obj(setNames(p, names(pars)), deriv = FALSE)$value, pars)
  expect_equal(unname(o$gradient[names(pars)]), unname(g_fd), tolerance = 1e-3)
})


test_that("normL2 BLOQ M4 + proportional errmodel: value and gradient match", {
  skip_if_no_compile()
  ec <- .nl2_fx()$prop

  pars <- c(A = 1.0, k = 0.5, srel = 0.1)
  data <- fx_decay_data_bloq(pars = pars[c("A", "k")], sigma = 0.05,
                             lloq = 0.1, times = seq(0, 10, by = 1))
  data$C1$sigma <- NA_real_

  d   <- data$C1
  obj <- normL2(data, ec$prd, errmodel = ec$e, opt.BLOQ = "M4NM")
  o <- obj(pars)

  env_pred <- attr(o, "env")$prediction[["C1"]]
  pred <- env_pred[match(d$time, env_pred[, "time"]), "y"]
  sigma_pred <- pars[["srel"]] * pred
  bloq <- d$value <= d$lloq

  expected <- truth_nll_aloq(pred[!bloq], d$value[!bloq], sigma_pred[!bloq]) +
              truth_nll_bloq_m4(pred[bloq], d$lloq[bloq], sigma_pred[bloq])
  expect_equal(o$value, expected, tolerance = 1e-6)

  skip_if_not_installed("numDeriv")
  g_fd <- numDeriv::grad(
    function(p) obj(setNames(p, names(pars)), deriv = FALSE)$value, pars)
  expect_equal(unname(o$gradient[names(pars)]), unname(g_fd), tolerance = 1e-3)
})


# ============================================================================
# Cross-backend parity (C++ kernel vs R reference)
# ============================================================================

test_that("normL2 cpp kernel agrees with R reference on the linear-decay fixture", {
  testthat::skip_if_not_installed("cppDE")
  testthat::skip_on_cran()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.05)
  pars  <- bench$outerpars_id

  with_cpp_backend(FALSE, {
    o_R <- normL2(data, bench$prd_id)(pars)
  })
  with_cpp_backend(TRUE, {
    o_C <- normL2(data, bench$prd_id)(pars)
  })
  expect_equal(o_C$value,    o_R$value,    tolerance = 1e-9)
  expect_equal(o_C$gradient, o_R$gradient, tolerance = 1e-8)
  expect_equal(o_C$hessian,  o_R$hessian,  tolerance = 1e-8)
})


# ============================================================================
# chi2 attribute
# ============================================================================

test_that("normL2 reports the sum of squares as a chi2 attribute", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.1)

  o <- normL2(data, bench$prd_id)(bench$outerpars_id)
  n <- sum(vapply(data, nrow, 0L))
  # value = chi2 + sum(log(2 pi sigma^2)), with one sigma over the fixture
  expect_equal(unname(attr(o, "chi2")),
               o$value - n * log(2 * pi * 0.1^2), tolerance = 1e-9)
  expect_equal(names(attr(o, "chi2")), "data")
})


test_that("terms sharing an attr.name pool their chi2, others split", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.1)

  one  <- normL2(data, bench$prd_id)
  chi1 <- unname(attr(one(bench$outerpars_id), "chi2"))

  same <- (one + normL2(data, bench$prd_id))(bench$outerpars_id)
  expect_equal(unname(attr(same, "chi2")), 2 * chi1, tolerance = 1e-9)
  expect_null(attr(same, "chi2_data"))

  split <- (one + normL2(data, bench$prd_id, attr.name = "validation"))(bench$outerpars_id)
  expect_null(attr(split, "chi2"))
  expect_equal(unname(attr(split, "chi2_data")), chi1, tolerance = 1e-9)
  expect_equal(unname(attr(split, "chi2_validation")), chi1, tolerance = 1e-9)

  # a third term folds back into the contribution it belongs to
  three <- (one + normL2(data, bench$prd_id, attr.name = "validation") +
              normL2(data, bench$prd_id))(bench$outerpars_id)
  expect_equal(unname(attr(three, "chi2_data")), 2 * chi1, tolerance = 1e-9)
  expect_equal(unname(attr(three, "chi2_validation")), chi1, tolerance = 1e-9)
})


# ============================================================================
# Printing
# ============================================================================

test_that("print.objlist skips the blocks a deriv = FALSE call does not carry", {
  o <- structure(list(value = -480.3, gradient = NULL, hessian = NULL),
                 class = "objlist")
  attr(o, "data") <- -480.3
  attr(o, "chi2") <- c(data = 541)
  out <- capture.output(print(o))
  expect_true(any(grepl("^value", out)))
  expect_false(any(grepl("^(gradient|hessian)\\[", out)))
  # the attribute block carries names and numbers, not storage modes
  expect_true(any(grepl("^ chi2 +541", out)))
  expect_false(any(grepl("num |chr |List of", out)))

  g <- setNames(c(1, 2, 3), letters[1:3])
  h <- matrix(0, 3, 3, dimnames = list(letters[1:3], letters[1:3]))
  o2 <- structure(list(value = 1.5, gradient = g, hessian = h), class = "objlist")
  out2 <- capture.output(print(o2))
  expect_true(any(grepl("gradient\\[1:3\\]", out2)))
  expect_true(any(grepl("hessian\\[1:3,1:3\\]", out2)))
})


test_that("normL2 gradient and Hessian follow the order of the parameter vector", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.1)

  p1 <- bench$outerpars_id
  p2 <- p1[rev(names(p1))]

  for_each_backend(function(cpp) {
    obj <- normL2(data, bench$prd_id)
    o1 <- obj(p1)
    o2 <- obj(p2)

    expect_identical(names(o1$gradient), names(p1), info = paste0("cpp=", cpp))
    expect_identical(names(o2$gradient), names(p2), info = paste0("cpp=", cpp))
    expect_equal(o1$value, o2$value, info = paste0("cpp=", cpp))
    expect_equal(o1$gradient[names(p2)], o2$gradient, info = paste0("cpp=", cpp))
    expect_equal(o1$hessian[names(p2), names(p2)], o2$hessian,
                 info = paste0("cpp=", cpp))
  })
})


## ---- hessian = FALSE: gradient without the Hessian ---------------------

test_that("hessian = FALSE returns value and gradient but no Hessian", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(pars = c(A = 1.0, k = 0.5), sigma = 0.05)
  obj   <- normL2(data, bench$prd_id)
  p     <- bench$outerpars_id

  full   <- obj(p)
  nohess <- obj(p, hessian = FALSE)
  expect_true(is.matrix(full$hessian))
  expect_null(nohess$hessian)
  expect_equal(nohess$value, full$value)
  expect_equal(nohess$gradient, full$gradient)

  # The NULL Hessian survives objective composition (normL2 + constraintL2).
  comp <- obj + constraintL2(setNames(rep(0, length(p)), names(p)), sigma = 4)
  expect_true(is.matrix(comp(p)$hessian))
  expect_null(comp(p, hessian = FALSE)$hessian)
  expect_equal(comp(p, hessian = FALSE)$gradient, comp(p)$gradient)
  expect_equal(comp(p, hessian = FALSE)$value,    comp(p)$value)
})


test_that("a prediction cut short is an error, not a read past its rows", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data <- fx_decay_data()
  obj  <- normL2(data, bench$prd_id)
  pars <- bench$outerpars_id
  full <- bench$prd_id(sort(unique(data[[1]]$time)), pars)
  v <- obj(pars, .prediction = full)$value

  # what the solver returns when it stops early: the first rows only
  cut <- full
  pr  <- full[[1]]
  keep <- seq_len(2L)
  at <- attributes(pr)
  short <- unclass(pr)[keep, , drop = FALSE]
  for (a in setdiff(names(at), c("dim", "dimnames"))) attr(short, a) <- at[[a]]
  attr(short, "deriv") <- attr(pr, "deriv")[keep, , , drop = FALSE]
  cut[[1]] <- short
  expect_error(obj(pars, .prediction = cut), "stopped early")
  # and the full prediction afterwards still reads its own rows
  expect_equal(obj(pars, .prediction = full)$value, v)
})
