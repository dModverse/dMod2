# normL2() against closed-form likelihoods, gradients and error models, BLOQ
# treatments, chi2 bookkeeping and argument checks. Hessians live in test-deriv2.R.

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
  obj <- normL2(data, ec$prd, errmodel = ec$e, optBLOQ = "M3")
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


# ============================================================================
# Printing
# ============================================================================

test_that("print.objlist skips the blocks a deriv = FALSE call lacks", {
  o <- structure(list(value = -480.3, gradient = NULL, hessian = NULL),
                 class = "objlist")
  attr(o, "data") <- -480.3
  attr(o, "chi2") <- c(data = 541)
  out <- capture.output(print(o))
  expect_true(any(grepl("^value", out)))
  expect_false(any(grepl("^(gradient|hessian)\\[", out)))
  # the attribute block shows names and numbers, not storage modes
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
