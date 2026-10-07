# scanL1(): L1 selection over a lambda grid, gates and anchored fusion blocks.

# Weighted least squares over a model y = f(theta) with Jacobian, as an
# objective of the shape trust() and trustL1() call: value, gradient and
# Gauss-Newton Hessian over the free parameters.
.lsq_obj <- function(model, data, sigma) {
  function(pars, fixed = NULL, ...) {
    all <- c(pars, fixed)
    m   <- model(all)
    r   <- (m$y - data) / sigma
    J   <- m$J[, names(pars), drop = FALSE] / sigma
    objlist(value = sum(r^2), gradient = setNames(drop(2 * crossprod(J, r)), names(pars)),
            hessian = 2 * crossprod(J))
  }
}

# y_i = sum_j 10^log10_kj * s_kj * x_ij: the gated linear model.
.gated_model <- function(X) function(p) {
  k  <- paste0("k", seq_len(ncol(X)))
  a  <- 10^p[paste0("log10_", k)]; s <- p[paste0("s_", k)]
  J  <- cbind(sweep(X, 2, log(10) * a * s, `*`), sweep(X, 2, a, `*`))
  colnames(J) <- c(paste0("log10_", k), paste0("s_", k))
  list(y = drop(X %*% (a * s)), J = J)
}

# Four cell types around one level mu, fold changes r_<type> to type R.
.type_model <- function(type) function(p) {
  r <- c(R = 0, p[c("r_B", "r_C", "r_D")])
  names(r) <- c("R", "B", "C", "D")
  J <- cbind(mu = 1, sapply(c("B", "C", "D"), function(t) as.numeric(type == t)))
  colnames(J) <- c("mu", "r_B", "r_C", "r_D")
  list(y = p[["mu"]] + r[type], J = J)
}


test_that("a tied objective pulls gradient and Hessian back through the tie", {
  f <- function(pars, fixed = NULL, deriv = TRUE, ...) {
    p <- c(pars, fixed)[c("a", "b", "c")]
    w <- c(1, 2, 3)
    objlist(value = sum(w * p^2) + p[["a"]] * p[["c"]],
            gradient = c(a = 2 * p[["a"]] + p[["c"]], b = 4 * p[["b"]],
                         c = 6 * p[["c"]] + p[["a"]])[names(pars)],
            hessian = matrix(c(2, 0, 1, 0, 4, 0, 1, 0, 6), 3,
                             dimnames = list(c("a", "b", "c"), c("a", "b", "c")))[
                               names(pars), names(pars)])
  }
  ft <- dMod2:::.l1Tie(f, list(a = c("a", "c")))
  out <- ft(c(a = 0.5, b = 1))
  expect_equal(out$value, 0.5^2 * (1 + 3) + 2 + 0.25)
  expect_equal(out$gradient[["a"]], 2 * 0.5 + 0.5 + 6 * 0.5 + 0.5)
  expect_equal(out$hessian["a", "a"], 2 + 1 + 1 + 6)
  expect_equal(out$hessian["a", "b"], 0)
})

test_that("a structure key names removed parameters and groups by exact equality", {
  st <- dMod2:::.l1Structure(c(s_k = 0, s_m = 0.4, r_A = 0, r_B = 0.2, r_C = 0.2),
                             gates = c("s_k", "s_m"), reference = NULL,
                             groups = list(list(pars = c("r_A", "r_B", "r_C"),
                                                anchor = 0)))
  expect_identical(st$removed, "s_k")
  expect_identical(st$key, "-s_k {r_A=0 | r_B,r_C}")
  st0 <- dMod2:::.l1Structure(c(r_A = 0.1, r_B = 0.2), NULL, NULL,
                              list(list(pars = c("r_A", "r_B"))))
  expect_identical(st0$key, "{r_A | r_B}")
})

test_that("scanL1 removes the gate of a term the data do not need", {
  set.seed(1)
  X     <- matrix(runif(60), 20, 3)
  model <- .gated_model(X)
  truth <- c(log10_k1 = log10(0.5), log10_k2 = -1, log10_k3 = log10(0.2),
             s_k1 = 1, s_k2 = 0, s_k3 = 1)
  obj <- .lsq_obj(model, model(truth)$y + rnorm(20, 0, 0.01), 0.01)
  gates <- list(s_k1 = "log10_k1", s_k2 = "log10_k2", s_k3 = "log10_k3")

  fit <- scanL1(obj, c(log10_k1 = 0, log10_k2 = -0.5, log10_k3 = 0), zero = gates,
                lambda = 10^seq(-1, 3, length.out = 5), fits = 3,
                control = list(ndata = 20))
  expect_identical(fit$selected, "-s_k2")
  expect_identical(fit$fit[["s_k2"]], 0)
  expect_identical(unname(fit$fit[c("s_k1", "s_k3")]), c(1, 1))
  expect_gt(fit$refits$p[fit$refits$key == "-s_k2"], 0.05)
  expect_s3_class(plot(fit), "ggplot")
  expect_s3_class(plot(fit, type = "test"), "ggplot")
})

test_that("scanL1 fuses equal fold changes and anchors the shared ones", {
  set.seed(2)
  type  <- rep(c("R", "B", "C", "D"), each = 8)
  model <- .type_model(type)
  truth <- c(mu = 1, r_B = 0, r_C = 0.3, r_D = 0.3)
  obj   <- .lsq_obj(model, model(truth)$y + rnorm(length(type), 0, 0.05), 0.05)
  start <- c(mu = 0, r_B = 0, r_C = 0, r_D = 0)
  block <- list(list(pars = c("r_B", "r_C", "r_D"), anchor = 0))
  grid  <- 10^seq(-1, 3, length.out = 9)

  fit <- scanL1(obj, start, groups = block, lambda = grid, fits = 3,
                control = list(ndata = length(type)))
  expect_identical(fit$selected, "{r_B=0 | r_C,r_D}")
  expect_identical(fit$fit[["r_B"]], 0)
  expect_identical(fit$fit[["r_C"]], fit$fit[["r_D"]])
  expect_lt(abs(fit$fit[["r_C"]] - 0.3), 0.05)
  expect_true(all(c(4, 2) %in% fit$refits$nfree))

  bic <- scanL1(obj, start, groups = block, lambda = grid, fits = 3,
                select = "bic", control = list(ndata = length(type)))
  expect_identical(bic$selected, "{r_B=0 | r_C,r_D}")

  q08 <- scanL1(obj, start, groups = block, lambda = grid, fits = 3, q = 0.8)
  expect_identical(q08$selected, "{r_B=0 | r_C,r_D}")

  expect_error(scanL1(obj, start), "at least one")
  expect_error(scanL1(obj, start, reference = "r_X"), "not among")
  expect_error(scanL1(obj, start, groups = block, select = "bic"), "ndata")
})

test_that("scanL1 runs on a compiled dMod objective", {
  fx  <- fx_decay_compiled()
  obj <- normL2(fx_decay_data(), fx$prd_log)
  fit <- scanL1(obj, c(A_log = 0.3, k_log = 0), reference = "A_log",
                lambda = 10^seq(0, 2, length.out = 3), fits = 2)
  expect_identical(fit$selected, "-A_log")
  expect_identical(fit$fit[["A_log"]], 0)
  expect_lt(abs(fit$fit[["k_log"]] - log(0.5)), 0.1)
})
