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
    # attachInput keeps the state alongside the observable, which is what the
    # full-chain test compares.
    g_full <- Y(c(y1 = "A"), f = NULL, states = c("A"),
                parameters = character(0), attachInput = TRUE,
                derivMode = "forward", compile = FALSE,
                modelname = "noparam_full_g")

    # The mixed trafo again, passing every input it does not map through.
    p_thru <- P(eqnvec(A = "a^2", k = "a * b"), condition = "C1",
                attachInput = TRUE, deriv2 = TRUE,
                modelname = "test_P_thru", compile = FALSE)

    compile(p_mix, p_rev, p_fwd, p_const_val, p_const_fwd, x_full, p_full, g_full,
            p_thru, output = "test_P_all", cores = test_cores())
    cache <<- list(p_mix = p_mix, p_rev = p_rev, p_fwd = p_fwd,
                   p_const_val = p_const_val, p_const_fwd = p_const_fwd,
                   x_full = x_full, p_full = p_full, g_full = g_full,
                   p_thru = p_thru)
    cache
  }
})


## ---- Identity transformation -------------------------------------------

test_that("Pexpl identity trafo round-trips and has identity Jacobian", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  outer <- c(A = 1.5, k = 0.3)
  inner <- bench$pfn_id(outer, deriv = TRUE)
  # parlist[[C1]] is a parvec with value + "deriv" attr.
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


## ---- attachInput: inputs passed through -------------------------------

# A = a^2, k = a * b, and every input that is not A or k handed on untouched.
# An input without a derivative row counts as fixed downstream, so a missing
# row zeroes the forward gradient along it.

test_that("Pexpl(attachInput = TRUE) gives every input it passes through its own row", {
  skip_if_no_compile()
  p <- .pexpl_fx()$p_thru
  outer <- c(a = 1.3, b = 0.7, s = 2, u = -1)
  out <- p(outer, deriv2 = TRUE)$C1

  thru <- c("a", "b", "s", "u")
  expect_identical(names(out), c("A", "k", thru))
  expect_equal(unclass(out)[thru], outer[thru], ignore_attr = TRUE)
  expect_null(attr(out, "fixed"))

  # s and u are read by nothing, and are directions of their own all the same.
  J <- attr(out, "deriv")
  expect_identical(dimnames(J), list(names(out), thru))
  expect_equal(unname(J["A", ]), c(2 * 1.3, 0, 0, 0))
  expect_equal(unname(J["k", ]), c(0.7, 1.3, 0, 0))
  expect_equal(unname(J[thru, ]), diag(4))

  # The identity has no curvature; the trafo's own rows keep theirs.
  H <- attr(out, "deriv2")
  expect_identical(dimnames(H), list(names(out), thru, thru))
  expect_equal(H[thru, , ], array(0, c(4, 4, 4)), ignore_attr = TRUE)
  expect_equal(H["A", "a", "a"], 2)
  expect_equal(H["k", "a", "b"], 1)
  expect_equal(H["k", "b", "a"], 1)

  # The value path hands on the same inputs, in the same order.
  expect_identical(names(p(outer, deriv = FALSE)$C1), names(out))
})

test_that("an input passed through that the caller fixed stays fixed", {
  skip_if_no_compile()
  p <- .pexpl_fx()$p_thru

  # One the trafo does not read is still handed on, as a constant.
  out <- p(c(a = 1.3, b = 0.7, u = -1), fixed = c(s = 2), deriv2 = TRUE)$C1
  expect_true("s" %in% names(out))
  expect_equal(unclass(out)[["s"]], 2)
  expect_identical(attr(out, "fixed"), "s")
  J <- attr(out, "deriv")
  expect_false("s" %in% rownames(J))
  expect_false("s" %in% colnames(J))
  expect_equal(unname(J["u", ]), c(0, 0, 1))
  expect_false("s" %in% dimnames(attr(out, "deriv2"))[[1L]])

  # One the trafo reads: A = a^2 no longer moves, k = a * b only along b.
  out <- p(c(b = 0.7, s = 2), fixed = c(a = 1.3))$C1
  expect_setequal(names(out), c("A", "k", "a", "b", "s"))
  expect_setequal(attr(out, "fixed"), c("A", "a"))
  J <- attr(out, "deriv")
  expect_setequal(colnames(J), c("b", "s"))
  expect_equal(J["k", "b"], 1.3)
  expect_equal(unname(J[c("b", "s"), c("b", "s")]), diag(2))
})

test_that("behind another trafo an input passed through keeps what it brought", {
  skip_if_no_compile()
  p <- .pexpl_fx()$p_thru
  J0 <- rbind(a = c(x = 1, y = 0), b = c(x = 0.5, y = 2), s = c(x = 3, y = -1))
  H0 <- array(0, c(3, 2, 2), dimnames = list(rownames(J0), c("x", "y"), c("x", "y")))
  H0["s", "x", "y"] <- H0["s", "y", "x"] <- 0.25
  H0["b", "x", "x"] <- 0.5
  pin <- as.parvec(c(a = 1.3, b = 0.7, s = 2), deriv = J0, deriv2 = H0)
  out <- p(pin, deriv2 = TRUE)$C1

  J <- attr(out, "deriv")
  H <- attr(out, "deriv2")
  expect_identical(colnames(J), c("x", "y"))
  expect_equal(J[c("a", "b", "s"), ], J0)
  expect_equal(H[c("a", "b", "s"), , ], H0)
  # The chain rule on the trafo's own rows, for reference.
  expect_equal(J["k", ], 0.7 * J0["a", ] + 1.3 * J0["b", ])
  expect_equal(H["k", , ], outer(J0["a", ], J0["b", ]) + outer(J0["b", ], J0["a", ]) +
                 1.3 * H0["b", , ])
})

test_that("the batched entry passes the inputs through as the kernel does", {
  skip_if_no_compile()
  kernel <- attr(.pexpl_fx()$p_thru, "mappings")$C1
  bf <- attr(kernel, "batchfn")
  pars  <- list(as.parvec(c(a = 1.3, b = 0.7, u = -1)),
                as.parvec(c(a = 0.4, b = -2, u = 3)))
  fixed <- list(as.parvec(c(s = 2)), as.parvec(c(s = 5)))

  for (d2 in c(FALSE, TRUE)) {
    res <- bf(pars, fixed, deriv = TRUE, deriv2 = d2,
              conditions = list("C1", "C1"), cores = 1L)
    for (i in 1:2) {
      ref <- kernel(pars[[i]], fixed[[i]], deriv = TRUE, deriv2 = d2)
      expect_equal(res[[i]], ref, tolerance = 1e-14, info = paste(d2, i))
      expect_true("u" %in% rownames(attr(res[[i]], "deriv")))
      expect_identical(attr(res[[i]], "fixed"), "s")
    }
  }
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
