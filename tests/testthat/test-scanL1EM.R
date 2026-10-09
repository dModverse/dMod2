# scanL1(lambda = "em"): posterior of a term, strength update and selection.

test_that("the posterior of a term matches quadrature on a fine grid", {
  for (cs in list(list(b = 0.3, g = 0, c = 400, l = 2, q = 1, one = FALSE),
                  list(b = 0.05, g = 1, c = 50, l = 3, q = 0.8, one = FALSE),
                  list(b = 0.8, g = -0.5, c = 20, l = 1.5, q = 0.8, one = TRUE))) {
    r <- dMod2:::.l1EmTerm(cs$b, cs$g, cs$c, cs$l, cs$q, cs$one)
    x <- seq(if (cs$one) 0 else -6, 6, length.out = 4e5)
    lf <- -(cs$g * (x - cs$b) + 0.5 * cs$c * (x - cs$b)^2) - cs$l * abs(x)^cs$q
    h <- exp(lf)
    expect_equal(r$logZ, log(sum(h) * diff(x)[1]), tolerance = 1e-4)
    expect_equal(r$Eq, sum(abs(x)^cs$q * h) / sum(h), tolerance = 1e-4)
  }
})

test_that("a term without data curvature has the prior as posterior", {
  r <- dMod2:::.l1EmTerm(0.01, 1e-30, 1e-40, 2, 0.8, TRUE)
  expect_equal(r$Eq, 1 / (0.8 * 2), tolerance = 1e-6)
  expect_equal(r$logZ, 0.5 * dMod2:::.l1EmLogNorm(2, 0.8, TRUE), tolerance = 1e-6)
})

test_that("the strength removes the gate of a term the data do not need", {
  set.seed(1)
  X     <- matrix(runif(60), 20, 3)
  model <- .gated_model(X)
  truth <- c(log10_k1 = log10(0.5), log10_k2 = -1, log10_k3 = log10(0.2),
             s_k1 = 1, s_k2 = 0, s_k3 = 1)
  obj <- .lsq_obj(model, model(truth)$y + rnorm(20, 0, 0.01), 0.01)
  gates <- list(s_k1 = "log10_k1", s_k2 = "log10_k2", s_k3 = "log10_k3")
  for (q in c(1, 0.8)) {
    fit <- scanL1(obj, c(log10_k1 = 0, log10_k2 = -0.5, log10_k3 = 0), zero = gates,
                  lambda = "em", q = q, fits = 3)
    expect_identical(fit$selected, "-s_k2")
    expect_identical(fit$fit[["s_k2"]], 0)
    tm <- fit$em$terms
    lam <- fit$em$lambda[["gate"]]
    expect_equal(lam, (3 / q) / sum(tm$Eq), tolerance = 0.01)
    expect_true(all(diff(fit$em$trace[, "value"]) < 1e-6 | abs(diff(fit$em$trace[, "gate"])) < 0.05))
  }
  expect_output(print(fit), "lambda by EM")
})

test_that("fold changes are selected relative to their full estimate", {
  type  <- rep(c("R", "B", "C", "D"), each = 8)
  model <- .type_model(type)
  start <- c(mu = 0, r_B = 0, r_C = 0, r_D = 0)
  ref   <- c("r_B", "r_C", "r_D")
  set.seed(2)
  obj <- .lsq_obj(model, model(c(mu = 1, r_B = 0, r_C = 0.3, r_D = -0.5))$y +
                    rnorm(length(type), 0, 0.05), 0.05)
  fit <- scanL1(obj, start, reference = ref, lambda = "em", q = 0.8, fits = 3)
  expect_identical(fit$selected, "-r_B")
  expect_identical(fit$fit[["r_B"]], 0)
  set.seed(2)
  obj0 <- .lsq_obj(model, model(c(mu = 1, r_B = 0, r_C = 0, r_D = 0))$y +
                     rnorm(length(type), 0, 0.05), 0.05)
  fit0 <- scanL1(obj0, start, reference = ref, lambda = "em", q = 0.8, fits = 3)
  expect_gt(fit0$em$lambda[["reference"]], fit$em$lambda[["reference"]])
  # every step is a likelihood ratio test of one term at level alpha
  st <- fit0$steps
  expect_identical(st$accepted, ifelse(st$step == "remove", st$p >= 0.05, st$p < 0.05))
  expect_true(all(st$df == 1))
  expect_identical(sort(c(st$term)), sort(ref))
})

test_that("lambda = \"em\" takes no groups", {
  obj <- function(...) NULL
  expect_error(scanL1(obj, c(a = 0, b = 0), groups = list(c("a", "b")), lambda = "em"),
               "no groups")
})
