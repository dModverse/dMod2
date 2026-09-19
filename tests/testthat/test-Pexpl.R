# Behavioral tests for the explicit parameter transformation Pexpl() / P().
#
# Verifies:
#   * identity trafo round-trips (value + Jacobian)
#   * log trafo gives Jacobian = diag(exp(theta)) = diag(p)
#   * a mixed nonlinear trafo's Jacobian matches the algebraic derivative
#   * derivMode "reverse" and "forward" agree on the value
#   * getParameters() consistency through composition (Y * Xs * P)
#
# Second-order chain rule is covered by test-deriv2-Pexpl.R.

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("cppDE")
  testthat::skip_on_cran()
}


# Every model the file compiles itself, in one shared object built on first use.
.pexpl_fx <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

    p_mix <- P(eqnvec(A = "a^2", k = "a * b"), condition = "C1",
               modelname = "test_P_mix", compile = FALSE)
    p_rev <- P(eqnvec(A = "exp(a)", k = "exp(b)"), condition = "C1",
               method = "explicit", derivMode = "reverse",
               modelname = "test_P_dm_rev", compile = FALSE)
    p_fwd <- P(eqnvec(A = "exp(a)", k = "exp(b)"), condition = "C1",
               method = "explicit", derivMode = "forward",
               modelname = "test_P_dm_fwd", compile = FALSE)

    const <- c(A = "1.0", B = "2.5")
    p_const_val <- Pexpl(const, deriv = FALSE, compile = FALSE,
                         modelname = "noparam_pexpl_val")
    p_const_fwd <- Pexpl(const, derivMode = "forward", compile = FALSE,
                         modelname = "noparam_pexpl_fwd")

    x_full <- Xs(odemodel(as.eqnvec(c(A = "-k*A")), modelname = "noparam_full_ode",
                          compile = FALSE, backend = "cppDE"))
    p_full <- Pexpl(c(A = "1.0", k = "0.5"), derivMode = "forward",
                    compile = FALSE, modelname = "noparam_full_p")
    # attach.input keeps the state alongside the observable, which is what the
    # full-chain test compares.
    g_full <- Y(c(y1 = "A"), f = NULL, states = c("A"),
                parameters = character(0), attach.input = TRUE,
                derivMode = "forward", compile = FALSE,
                modelname = "noparam_full_g")

    compile(p_mix, p_rev, p_fwd, p_const_val, p_const_fwd, x_full, p_full, g_full,
            output = "test_P_all", cores = 4L)
    cache <<- list(p_mix = p_mix, p_rev = p_rev, p_fwd = p_fwd,
                   p_const_val = p_const_val, p_const_fwd = p_const_fwd,
                   x_full = x_full, p_full = p_full, g_full = g_full)
    cache
  }
})


## ---- Identity transformation -------------------------------------------

test_that("Pexpl identity trafo round-trips and has identity Jacobian", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  outer <- c(A = 1.5, k = 0.3)
  inner <- bench$pfn_id(outer, deriv = TRUE)
  # parlist[[C1]] is a parvec carrying value + "deriv" attr.
  pv <- inner$C1
  expect_equal(as.numeric(pv), as.numeric(outer))
  J  <- attr(pv, "deriv")
  expect_equal(unname(J), diag(2))
  expect_setequal(rownames(J), c("A", "k"))
  expect_setequal(colnames(J), c("A", "k"))
})


## ---- Log transformation -----------------------------------------------

test_that("Pexpl log trafo maps theta -> exp(theta) with Jacobian diag(exp(theta))", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  outer <- c(A_log = log(1.5), k_log = log(0.3))
  inner <- bench$pfn_log(outer, deriv = TRUE)
  pv <- inner$C1
  expect_equal(as.numeric(pv), c(1.5, 0.3))

  # J[i, j] = d (inner_i) / d (outer_j). Diagonal entries are exp(theta) = p.
  J <- attr(pv, "deriv")
  expect_equal(unname(diag(J)), c(1.5, 0.3))
  expect_equal(J[upper.tri(J)], rep(0, sum(upper.tri(J))))
  expect_equal(J[lower.tri(J)], rep(0, sum(lower.tri(J))))
})


## ---- Mixed nonlinear trafo: analytical Jacobian -----------------------

test_that("Pexpl Jacobian on a mixed nonlinear trafo equals the algebraic derivative", {
  skip_if_no_compile()

  # Mixed trafo: A = a^2, k = a * b
  # Analytical Jacobian: J[1,] = (2a, 0), J[2,] = (b, a).
  pfn <- .pexpl_fx()$p_mix

  outer <- c(a = 1.3, b = 0.7)
  inner <- pfn(outer, deriv = TRUE)$C1
  J <- attr(inner, "deriv")

  J_ref <- rbind(
    A = c(a = 2 * outer[["a"]], b = 0),
    k = c(a = outer[["b"]],     b = outer[["a"]]))
  expect_equal(unname(J), unname(J_ref), tolerance = 1e-8)
})


## ---- derivMode parity --------------------------------------------------

test_that("Pexpl derivMode 'reverse' and 'forward' agree on the value", {
  skip_if_no_compile()
  fx <- .pexpl_fx()
  pfn_rev <- fx$p_rev
  pfn_fwd <- fx$p_fwd

  outer <- c(a = 0.3, b = -0.5)
  i_rev <- pfn_rev(outer, deriv = TRUE)$C1
  i_fwd <- pfn_fwd(outer, deriv = TRUE)$C1
  expect_equal(as.numeric(i_rev), as.numeric(i_fwd), tolerance = 1e-12)
  J <- attr(i_fwd, "deriv")
  expect_equal(unname(J[c("A", "k"), c("a", "b")]), diag(exp(c(0.3, -0.5))),
               tolerance = 1e-12)
})


## ---- Parameter set propagation through composition --------------------

test_that("getParameters(Y * Xs * P) equals getParameters(P) (outer-pars view)", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  expect_setequal(getParameters(bench$prd_id),  getParameters(bench$pfn_id))
  expect_setequal(getParameters(bench$prd_log), getParameters(bench$pfn_log))
})


# ============================================================================
# Edge case: Pexpl with pure-numeric trafo (no outer parameters)
# ============================================================================

test_that("Pexpl with pure-numeric trafo evaluates (values and forward)", {
  fx <- .pexpl_fx()
  out_val <- fx$p_const_val(c(dummy = 1.0))
  expect_equal(unclass(out_val[[1]])[c("A", "B")], c(A = 1.0, B = 2.5))

  out_fwd <- fx$p_const_fwd(c(dummy = 1.0))
  expect_equal(unclass(out_fwd[[1]])[c("A", "B")], c(A = 1.0, B = 2.5))
})


test_that("Full g*x*p chain with constant-only Pexpl evaluates", {
  fx <- .pexpl_fx()
  out <- (fx$g_full * fx$x_full * fx$p_full)(seq(0, 5, length.out = 3), c(dummy = 1.0))
  pred <- out[[1]]
  expect_equal(unname(pred[, "y1"]), unname(pred[, "A"]))
  expect_equal(unname(pred[1, "A"]), 1.0, tolerance = 1e-8)
  expect_equal(unname(pred[3, "A"]), exp(-0.5 * 5), tolerance = 1e-4)
})

test_that("an uncompiled Pexpl asks for compile()", {
  withr::local_dir(tempdir())
  p <- Pexpl(c(A = "exp(a)"), compile = FALSE, modelname = "uncompiled_pexpl")
  expect_error(p(c(a = 1)), "is not compiled; call compile\\(\\)")
  expect_error(p(c(a = 1), deriv = FALSE), "is not compiled; call compile\\(\\)")
})
