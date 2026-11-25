# normL2() against closed-form likelihoods, gradients and error models, BLOQ
# treatments, chi2 bookkeeping and argument checks. Hessians live in test-deriv2.R.

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("cppDE")
  testthat::skip_on_cran()
}


# Error-model chains (absolute and relative) and a second-condition
# trafo on the shared decay fixture, in one shared object built on first use.
# The trafos pass the error parameters through, so normL2's errmodel call sees them.
.nl2_fx <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    bench <- fx_decay_compiled()
    oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

    e_const <- Y(c(y = "sigma_y"), f = bench$gfn, attachInput = FALSE,
                 condition = "C1", modelname = "nl2_err_const", compile = FALSE)
    p_sig <- P(eqnvec(A = "A", k = "k", sigma_y = "sigma_y"), condition = "C1",
               modelname = "nl2_p_sig", compile = FALSE)
    e_prop <- Y(c(y = "srel * y"), f = bench$gfn, attachInput = FALSE,
                condition = "C1", modelname = "nl2_err_prop", compile = FALSE)
    p_prop <- P(eqnvec(A = "A", k = "k", srel = "srel"), condition = "C1",
                modelname = "nl2_p_prop", compile = FALSE)
    p_C2 <- P(eqnvec(A = "A", k = "k"), condition = "C2",
              modelname = "nl2_p_id", compile = FALSE)
    # Two sources from 0, so the observed ratio is 0/0 at t = 0 only.
    re <- eqnlist() |> addReaction("", "A", "k1") |> addReaction("", "B", "k2")
    x_nan <- Xs(odemodel(re, modelname = "nl2nan_ode", compile = FALSE), condition = "C1")
    g_nan <- Y(c(frac = "A/(A+B)"), re, modelname = "nl2nan_obs", compile = FALSE)
    p_nan <- P(eqnvec(A = "0", B = "0", k1 = "k1", k2 = "k2"), condition = "C1",
               modelname = "nl2nan_p", compile = FALSE)
    compile(e_const, p_sig, e_prop, p_prop, p_C2, x_nan, g_nan, p_nan,
            output = "nl2_all", cores = test_cores())

    cache <<- list(
      const = list(prd = bench$gfn * bench$xfn * p_sig,  e = e_const),
      prop  = list(prd = bench$gfn * bench$xfn * p_prop, e = e_prop),
      pfn_C2 = p_C2,
      nan = g_nan * x_nan * p_nan)
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

  obj <- normL2(data, bench$prd_id)
  o <- obj(bench$outerpars_id)
  expect_equal(o$value, closed, tolerance = 1e-10)
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

  obj <- normL2(data, bench$prd_id)
  g_ana <- obj(pars)$gradient
  expect_equal(unname(g_ana[names(g_ref)]), unname(g_ref),
               tolerance = 1e-4)
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

  o_col <- normL2(data_col, ec$prd)(pars_em)
  o_em  <- normL2(data_em,  ec$prd, errmodel = ec$e)(pars_em)
  expect_equal(o_em$value, o_col$value, tolerance = 1e-10)
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

  o_joint <- normL2(data_multi, prd_multi)(pars)
  o_C1 <- normL2(data_multi["C1"], prd_multi)(pars)
  o_C2 <- normL2(data_multi["C2"], prd_multi)(pars)
  expect_equal(o_joint$value, o_C1$value + o_C2$value, tolerance = 1e-10)
})


# ---- BLOQ: closed-form value (M3) ---------------------------------------

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

  obj <- normL2(data, bench$prd_id)
  g_ana <- obj(pars)$gradient
  expect_equal(unname(g_ana[names(g_ref)]), unname(g_ref),
               tolerance = 1e-3)
})


test_that("normL2(optBLOQ = ...) selects the BLOQ method", {
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

  for (mode in names(expected)) {
    obj <- normL2(data, bench$prd_id, optBLOQ = mode)
    o   <- obj(pars)
    expect_equal(o$value, expected[[mode]], tolerance = 1e-3, info = mode)
  }
  # M3 is the default
  expect_equal(normL2(data, bench$prd_id)(pars)$value, expected[["M3"]], tolerance = 1e-3)
})


test_that("normL2 rejects unknown optBLOQ values", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data_bloq(sigma = 0.05, lloq = 0.1,
                              times = seq(0, 10, by = 1))
  # match.arg() rejects the unknown value; assert on the listed choices rather
  # than the "should be one of" prefix, which match.arg translates per locale.
  expect_error(normL2(data, bench$prd_id, optBLOQ = "M2"),
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

  obj <- normL2(data, ec$prd, errmodel = ec$e)
  o   <- obj(pars)
  expect_equal(o$value, expected, tolerance = 1e-4)
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

  obj <- normL2(data, ec$prd, errmodel = ec$e)
  g_ana <- obj(pars)$gradient
  expect_equal(unname(g_ana[names(g_ref)]), unname(g_ref),
               tolerance = 1e-3)
})


test_that("hessian = \"expected\" gives the expected information, positive semidefinite with an error model", {
  skip_if_no_compile()
  ec <- .nl2_fx()$prop
  data <- fx_decay_data(pars = c(A = 1.0, k = 0.5), sigma = 0.05)
  data$C1$sigma <- NA_real_
  pars <- c(A = 1.2, k = 0.45, srel = 0.08)
  t <- data$C1$time
  A <- pars[["A"]]; k <- pars[["k"]]; srel <- pars[["srel"]]
  pred <- A * exp(-k * t)
  sigma <- srel * pred
  dpred <- cbind(A = exp(-k * t), k = -t * A * exp(-k * t), srel = 0)
  dsigma <- srel * dpred
  dsigma[, "srel"] <- pred
  h_ref <- crossprod(dpred, 2 / sigma^2 * dpred) + crossprod(dsigma, 4 / sigma^2 * dsigma)
  obj <- normL2(data, ec$prd, errmodel = ec$e)
  h <- obj(pars, hessian = "expected")$hessian[colnames(h_ref), colnames(h_ref)]
  expect_equal(unname(h), unname(h_ref), tolerance = 1e-4)
  expect_gt(min(eigen(h, symmetric = TRUE, only.values = TRUE)$values), 0)
  # the Gauss-Newton Hessian keeps the exact curvature of the loss, value and gradient are shared
  expect_equal(obj(pars, hessian = "expected")$gradient, obj(pars)$gradient)
  expect_false(isTRUE(all.equal(obj(pars)$hessian, obj(pars, hessian = "expected")$hessian)))
})


test_that("rows with explicit sigma keep it; NA rows fall through to errmodel", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  ec <- .nl2_fx()$prop

  pars <- c(A = 1.0, k = 0.5, srel = 0.1)
  data <- fx_decay_data(pars = pars[c("A", "k")], sigma = 0.05)
  n <- nrow(data$C1)
  data$C1$sigma <- ifelse(seq_len(n) <= n %/% 2, NA_real_, 0.05)

  data_na <- data; data_na$C1 <- data_na$C1[is.na(data$C1$sigma), ]
  data_ex <- data; data_ex$C1 <- data_ex$C1[!is.na(data$C1$sigma), ]
  # Use a common time grid across all three so the adaptive integrator
  # produces bit-identical predictions at the shared data times.
  all_times <- sort(unique(data$C1$time))
  v_full <- normL2(data,    ec$prd, errmodel = ec$e, times = all_times)(pars)$value
  v_na   <- normL2(data_na, ec$prd, errmodel = ec$e, times = all_times)(pars)$value
  v_ex   <- normL2(data_ex, ec$prd, errmodel = ec$e, times = all_times)(pars)$value
  expect_equal(v_full, v_na + v_ex, tolerance = 1e-9)
})


test_that("a sigma column of NA alone is numeric and falls through to errmodel", {
  skip_if_no_compile()
  ec <- .nl2_fx()$prop
  pars <- c(A = 1.0, k = 0.5, srel = 0.1)
  d <- fx_decay_data(pars = pars[c("A", "k")])$C1
  d$sigma <- NA
  lgl <- as.datalist(list(C1 = d))
  d$sigma <- NA_real_
  dbl <- as.datalist(list(C1 = d))

  expect_type(lgl$C1$sigma, "double")
  expect_equal(normL2(lgl, ec$prd, errmodel = ec$e)(pars),
               normL2(dbl, ec$prd, errmodel = ec$e)(pars))
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
# The errmodel sigma and the sigma derivatives reach the BLOQ rows as well as the
# ALOQ rows: M3 value against the closed form, gradient against differences.
test_that("normL2 BLOQ M3 + proportional errmodel: value and gradient match", {
  skip_if_no_compile()
  ec <- .nl2_fx()$prop

  pars <- c(A = 1.0, k = 0.5, srel = 0.1)
  data <- fx_decay_data_bloq(pars = pars[c("A", "k")], sigma = 0.05,
                             lloq = 0.1, times = seq(0, 10, by = 1))
  data$C1$sigma <- NA_real_   # both ALOQ and BLOQ rows draw sigma from errmodel

  d   <- data$C1
  obj <- normL2(data, ec$prd, errmodel = ec$e, optBLOQ = "M3")
  o <- obj(pars)

  # Use the exact prediction normL2 integrated (the normL2 env), so the closed-form
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


test_that("the native chain gives the objective of the R path, error model and censored rows included", {
  skip_if_no_compile()
  ec <- .nl2_fx()$prop
  pars <- c(A = 1.0, k = 0.5, srel = 0.1)
  data <- fx_decay_data_bloq(pars = pars[c("A", "k")], sigma = 0.05, lloq = 0.1, times = seq(0, 10, by = 1))
  data$C1$sigma[seq(1L, nrow(data$C1), by = 2L)] <- NA_real_
  obj <- normL2(data, ec$prd, errmodel = ec$e, optBLOQ = "M3")
  run <- function(chain, ...) withr::with_options(list(dMod2.evaluation = if (chain) "C++" else "R"), obj(...))
  for (args in list(list(pars), list(pars, deriv = FALSE), list(pars[-1], fixed = pars[1]))) {
    a <- do.call(run, c(list(TRUE), args))
    b <- do.call(run, c(list(FALSE), args))
    expect_equal(a$value, b$value, tolerance = 1e-12)
    expect_equal(a$gradient, b$gradient, tolerance = 1e-12)
    expect_equal(a$hessian, b$hessian, tolerance = 1e-12)
    expect_equal(attr(a, "chi2"), attr(b, "chi2"), tolerance = 1e-12)
  }
  entries <- environment(obj)$.chain$entries
  expect_false(isTRUE(get(ls(entries)[1], envir = entries)[[1]]$none))
  # the predictions of the environment, computed when asked for
  expect_equal(attr(run(TRUE, pars), "env")$prediction$C1, attr(run(FALSE, pars), "env")$prediction$C1)

  mc <- fx_decay_multicond_compiled()
  dm <- fx_decay_data_multi(setNames(rep(list(c(A = 1, k = 0.5)), 4), mc$conditions))
  objm <- normL2(dm, mc$prd)
  a <- withr::with_options(list(dMod2.evaluation = "C++", dMod2.cores = 2L), objm(mc$outerpars))
  b <- withr::with_options(list(dMod2.evaluation = "R"), objm(mc$outerpars))
  expect_error(withr::with_options(list(dMod2.evaluation = "native"), objm(mc$outerpars)), "C\\+\\+")
  expect_equal(a$value, b$value, tolerance = 1e-12)
  expect_equal(a$gradient, b$gradient, tolerance = 1e-12)
  expect_equal(a$hessian, b$hessian, tolerance = 1e-12)
})

test_that("hessian = \"sparse\" takes the native chain and gives the sparse Hessian of the R path", {
  skip_if_no_compile()
  ec <- .nl2_fx()$prop
  pars <- c(A = 1.0, k = 0.5, srel = 0.1)
  data <- fx_decay_data_bloq(pars = pars[c("A", "k")], sigma = 0.05, lloq = 0.1, times = seq(0, 10, by = 1))
  data$C1$sigma[seq(1L, nrow(data$C1), by = 2L)] <- NA_real_
  mc <- fx_decay_multicond_compiled()
  dm <- fx_decay_data_multi(setNames(rep(list(c(A = 1, k = 0.5)), 4), mc$conditions))
  calls <- 0L
  chainEval <- .chainEval
  local_mocked_bindings(.chainEval = function(...) {
    calls <<- calls + 1L
    chainEval(...)
  })
  cases <- list(list(normL2(data, ec$prd, errmodel = ec$e, optBLOQ = "M3"), pars),
                list(normL2(dm, mc$prd), mc$outerpars),
                list(normL2(dm, mc$prd), mc$outerpars[-1], fixed = mc$outerpars[1]),
                list(normL2(dm, mc$prd), c(mc$outerpars, unused = 1)))
  for (args in cases) {
    obj <- args[[1]]
    run <- function(chain, ...) withr::with_options(list(dMod2.evaluation = if (chain) "C++" else "R"),
                                                    do.call(obj, c(args[-1], list(...))))
    before <- calls
    a <- run(TRUE, hessian = "sparse")
    expect_equal(calls, before + 1L)
    b <- run(FALSE, hessian = "sparse")
    d <- run(TRUE)
    expect_s4_class(a$hessian, "dgCMatrix")
    expect_null(a$lowrank)
    expect_equal(a$value, b$value, tolerance = 1e-12)
    expect_equal(a$gradient, b$gradient, tolerance = 1e-12)
    expect_equal(as.matrix(a$hessian), as.matrix(b$hessian), tolerance = 1e-12)
    expect_equal(Matrix::nnzero(a$hessian), Matrix::nnzero(b$hessian))
    expect_identical(as.matrix(a$hessian), d$hessian)
    expect_identical(a$gradient, d$gradient)
  }
})

test_that("normL2 BLOQ M4 + proportional errmodel: value and gradient match", {
  skip_if_no_compile()
  ec <- .nl2_fx()$prop

  pars <- c(A = 1.0, k = 0.5, srel = 0.1)
  data <- fx_decay_data_bloq(pars = pars[c("A", "k")], sigma = 0.05,
                             lloq = 0.1, times = seq(0, 10, by = 1))
  data$C1$sigma <- NA_real_

  d   <- data$C1
  obj <- normL2(data, ec$prd, errmodel = ec$e, optBLOQ = "M4NM")
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


# ---- chi2 attribute -----------------------------------------------------

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


test_that("terms sharing an attrName pool their chi2, others split", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.1)

  one  <- normL2(data, bench$prd_id)
  chi1 <- unname(attr(one(bench$outerpars_id), "chi2"))

  same <- (one + normL2(data, bench$prd_id))(bench$outerpars_id)
  expect_equal(unname(attr(same, "chi2")), 2 * chi1, tolerance = 1e-9)
  expect_null(attr(same, "chi2_data"))

  split <- (one + normL2(data, bench$prd_id, attrName = "validation"))(bench$outerpars_id)
  expect_null(attr(split, "chi2"))
  expect_equal(unname(attr(split, "chi2_data")), chi1, tolerance = 1e-9)
  expect_equal(unname(attr(split, "chi2_validation")), chi1, tolerance = 1e-9)

  # a third term folds back into the contribution it belongs to
  three <- (one + normL2(data, bench$prd_id, attrName = "validation") +
              normL2(data, bench$prd_id))(bench$outerpars_id)
  expect_equal(unname(attr(three, "chi2_data")), 2 * chi1, tolerance = 1e-9)
  expect_equal(unname(attr(three, "chi2_validation")), chi1, tolerance = 1e-9)
})


test_that("normL2 gradient and Hessian follow the order of the parameter vector", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.1)

  p1 <- bench$outerpars_id
  p2 <- p1[rev(names(p1))]

  obj <- normL2(data, bench$prd_id)
  o1 <- obj(p1)
  o2 <- obj(p2)

  expect_identical(names(o1$gradient), names(p1))
  expect_identical(names(o2$gradient), names(p2))
  expect_equal(o1$value, o2$value)
  expect_equal(o1$gradient[names(p2)], o2$gradient)
  expect_equal(o1$hessian[names(p2), names(p2)], o2$hessian)
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


test_that("a prediction cut short is an error, not a read past the prediction rows", {
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
  # and the full prediction afterwards still reads the full rows
  expect_equal(obj(pars, .prediction = full)$value, v)
})


test_that("an undefined observable counts only at a data point", {
  skip_if_no_compile()
  prd  <- .nl2_fx()$nan
  pars <- c(k1 = 1, k2 = 3)

  later <- as.datalist(data.frame(name = "frac", time = c(1, 2), value = 0.25,
                                  sigma = 0.1, condition = "C1"))
  v <- normL2(later, prd)(pars)
  expect_true(is.finite(v$value))
  expect_true(all(is.finite(v$gradient)))
  expect_equal(v$value, 2 * log(2 * pi * 0.1^2), tolerance = 1e-8)

  at0 <- as.datalist(data.frame(name = "frac", time = 0, value = 0.25,
                                sigma = 0.1, condition = "C1"))
  expect_error(normL2(at0, prd)(pars), "NaN at data point\\(s\\) of condition 'C1': frac \\(t = 0\\)")
})

test_that("evalConditionResidual keeps fixed error parameters without sensitivities", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data <- fx_decay_data(sigma = 0.1)
  data$C1$sigma <- NA_real_
  ec <- .nl2_fx()$const
  pars <- bench$outerpars_id
  fixed <- c(sigma_y = 0.1)
  pr <- ec$prd(times = data$C1$time, pars = pars, fixed = fixed, deriv = FALSE)[["C1"]]
  r <- evalConditionResidual(data$C1, pr, errfn = ec$e, fixed = fixed, cn = "C1",
                             deriv = FALSE)
  ref <- normL2(data, ec$prd, errmodel = ec$e)(pars, fixed = fixed, deriv = FALSE)
  expect_equal(r$value, ref$value, tolerance = 1e-10)
})


test_that("an error parameter outside the prediction's derivatives has a gradient", {
  x <- Xd(data.frame(name = "y", time = 0:3, row.names = paste0("a", 0:3)), condition = "C1")
  X2Y <- function(out, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE) {
    p <- c(unclass(pars), unclass(fixed))
    sig <- rep(exp(p[["ls"]]), nrow(out))
    jac <- if (deriv) array(sig, c(nrow(out), 1, 1), dimnames = list(NULL, "y", "ls"))
    prdframe(cbind(time = out[, "time"], y = sig), deriv = jac, parameters = p)
  }
  e <- obsfn(X2Y, "ls", "C1")
  data <- as.datalist(data.frame(name = "y", time = 0:3, value = c(1, 2, 2.5, 4), sigma = NA,
                                 condition = "C1"))
  obj <- normL2(data, x, errmodel = e)
  pars <- c(a0 = 1.2, a1 = 1.7, a2 = 2.9, a3 = 3.6, ls = -0.3)
  num <- numDeriv::grad(function(q) obj(setNames(q, names(pars)), deriv = FALSE)$value, unname(pars))
  expect_equal(unname(obj(pars)$gradient[names(pars)]), num, tolerance = 1e-6)
  pr <- x(times = 0:3, pars = pars, conditions = "C1")[["C1"]]
  r <- evalConditionResidual(data$C1, pr, errfn = e, cn = "C1")
  expect_equal(unname(r$gradient[names(pars)]), num, tolerance = 1e-6)
})


# A prediction and an error model with exact second derivatives written out, so the Hessians of
# the objectives need no compiled model.
.nl2_handmade <- function() {
  P2X <- function(times, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE) {
    p <- c(unclass(pars), unclass(fixed)); an <- paste0("a", 0:3)
    jac <- array(0, c(length(times), 1, 4), dimnames = list(NULL, "y", an))
    for (i in seq_along(times)) jac[i, 1, times[i] + 1] <- 1
    jac2 <- if (deriv2) array(0, c(length(times), 1, 4, 4), dimnames = list(NULL, "y", an, an))
    prdframe(cbind(time = times, y = p[an][times + 1]), deriv = jac, deriv2 = jac2, parameters = p)
  }
  x <- prdfn(P2X, paste0("a", 0:3), "C1")
  X2Y <- function(out, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE) {
    p <- c(unclass(pars), unclass(fixed))
    sig <- exp(p[["ls"]]) * (1 + 0.1 * out[, "y"]^2)
    jac <- jac2 <- NULL
    if (deriv || deriv2) {
      jac <- array(0, c(nrow(out), 1, 5), dimnames = list(NULL, "y", c(paste0("a", 0:3), "ls")))
      for (i in seq_len(nrow(out))) {
        jac[i, 1, i] <- exp(p[["ls"]]) * 0.2 * out[i, "y"]
        jac[i, 1, 5] <- sig[i]
      }
    }
    if (deriv2) {
      jac2 <- array(0, c(nrow(out), 1, 5, 5), dimnames = list(NULL, "y", dimnames(jac)[[3]], dimnames(jac)[[3]]))
      for (i in seq_len(nrow(out))) {
        jac2[i, 1, i, i] <- exp(p[["ls"]]) * 0.2
        jac2[i, 1, i, 5] <- jac2[i, 1, 5, i] <- exp(p[["ls"]]) * 0.2 * out[i, "y"]
        jac2[i, 1, 5, 5] <- sig[i]
      }
    }
    prdframe(cbind(time = out[, "time"], y = sig), deriv = jac, deriv2 = jac2, parameters = p)
  }
  list(x = x, e = obsfn(X2Y, "ls", "C1"))
}


test_that("normT has the Student-t likelihood with exact gradient and Hessian, censored rows too", {
  hm <- .nl2_handmade()
  x <- hm$x
  e <- hm$e
  data <- as.datalist(data.frame(name = "y", time = 0:3, value = c(1, 2, 2.5, 0.2), sigma = NA,
                                 lloq = c(-Inf, -Inf, -Inf, 0.5), condition = "C1"))
  obj <- normT(data, x, errmodel = e, df = 4)
  pars <- c(a0 = 1.2, a1 = 1.7, a2 = 2.9, a3 = 0.9, ls = -0.3)
  v <- function(q) obj(setNames(q, names(pars)), deriv = FALSE)$value
  gr <- function(q) unname(obj(setNames(q, names(pars)))$gradient[names(pars)])
  # the likelihood itself, rows above the limit and the censored one
  sg <- exp(-0.3) * (1 + 0.1 * pars[1:4]^2)
  ref <- -2 * sum(dt((pars[1:3] - c(1, 2, 2.5)) / sg[1:3], 4, log = TRUE) - log(sg[1:3])) -
    2 * pt((0.5 - pars[[4]]) / sg[[4]], 4, log.p = TRUE)
  expect_equal(v(unname(pars)), unname(ref), tolerance = 1e-10)
  out <- obj(pars, deriv2 = TRUE)
  expect_equal(unname(out$gradient[names(pars)]), numDeriv::grad(v, unname(pars)), tolerance = 1e-7)
  expect_equal(unname(out$hessian[names(pars), names(pars)]), numDeriv::jacobian(gr, unname(pars)),
               tolerance = 1e-6)
  expect_error(normT(data, x, df = 0), "positive")
})


test_that("normL2 has the exact Hessian with censored rows under every BLOQ method", {
  hm <- .nl2_handmade()
  data <- as.datalist(data.frame(name = "y", time = 0:3, value = c(1, 2, 2.5, 0.2), sigma = NA,
                                 lloq = c(-Inf, -Inf, 0.4, 0.5), condition = "C1"))
  pars <- c(a0 = 1.2, a1 = 1.7, a2 = 2.9, a3 = 0.9, ls = -0.3)
  for (mode in c("M3", "M4NM", "M4BEAL")) {
    obj <- normL2(data, hm$x, errmodel = hm$e, optBLOQ = mode)
    v <- function(q) obj(setNames(q, names(pars)), deriv = FALSE)$value
    gr <- function(q) unname(obj(setNames(q, names(pars)))$gradient[names(pars)])
    out <- obj(pars, deriv2 = TRUE)
    expect_equal(unname(out$gradient[names(pars)]), numDeriv::grad(v, unname(pars)),
                 tolerance = 1e-7, info = mode)
    expect_equal(unname(out$hessian[names(pars), names(pars)]), numDeriv::jacobian(gr, unname(pars)),
                 tolerance = 1e-6, info = mode)
  }
})


test_that("the row loss derivatives of the Student-t likelihood match the closed-form expressions", {
  f <- function(x, y, cens, nu) if (cens) -2 * pt((y - x[1]) / x[2], nu, log.p = TRUE) else
    (nu + 1) * log1p(((x[1] - y) / x[2])^2 / nu) + 2 * log(x[2])
  for (cens in c(FALSE, TRUE)) {
    D <- dMod2:::.rowLossDerivs(1.3, 0.7, 0.9, cens, df = 4)
    d3 <- numDeriv::jacobian(function(x) {
      E <- dMod2:::.rowLossDerivs(x[1], x[2], 0.9, cens, df = 4); c(E[1, 3, 1], E[1, 2, 2], E[1, 1, 3]) }, c(1.3, 0.7))
    h <- numDeriv::hessian(function(x) f(x, 0.9, cens, 4), c(1.3, 0.7))
    expect_equal(c(D[1, 3, 1], D[1, 2, 2], D[1, 1, 3]), c(h[1, 1], h[1, 2], h[2, 2]), tolerance = 1e-7)
    expect_equal(c(D[1, 4, 1], D[1, 3, 2], D[1, 2, 3], D[1, 1, 4]), c(d3[1, 1], d3[1, 2], d3[2, 2], d3[3, 2]),
                 tolerance = 1e-6)
  }
  # normal rows above the limit in closed form, as the series gives them
  p <- c(1.3, -0.2, 4); s <- c(0.7, 1.5, 0.2); y <- c(0.9, 0.3, 4.1)
  expect_equal(dMod2:::.rowLossDerivs(p, s, y, rep(FALSE, 3)),
               dMod2:::.rowLossTaylor(p, s, y, rep(FALSE, 3)), tolerance = 1e-12)
})
