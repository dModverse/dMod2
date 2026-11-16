# Sparse Hessians: normL2(hessian = "sparse"), sums of objectives, trust(sparse = TRUE).

.sparseFixture <- function() {
  grid <- function(cn) data.frame(name = "y", time = 0:3, row.names = paste0(cn, "_a", 0:3))
  x <- Xd(grid("C1"), condition = "C1") + Xd(grid("C2"), condition = "C2")
  set.seed(3)
  data <- as.datalist(data.frame(name = "y", time = rep(0:3, 2), value = rnorm(8, 2), sigma = 0.3,
                                 condition = rep(c("C1", "C2"), each = 4)))
  pars <- setNames(rnorm(8, 2), c(paste0("C1_a", 0:3), paste0("C2_a", 0:3)))
  list(obj = normL2(data, x), pars = pars)
}

test_that("the sparse Hessian of normL2 equals the dense one, also in a sum", {
  f <- .sparseFixture()
  obj <- f$obj + constraintL2(c(C1_a0 = 1, C2_a3 = 2), sigma = 1)
  dense <- obj(f$pars)
  sp <- obj(f$pars, hessian = "sparse")
  expect_s4_class(sp$hessian, "sparseMatrix")
  expect_equal(as.matrix(sp$hessian)[names(f$pars), names(f$pars)],
               dense$hessian[names(f$pars), names(f$pars)])
  expect_equal(sp$gradient, dense$gradient)
})

test_that("trust(sparse = TRUE) reaches the dense optimum", {
  f <- .sparseFixture()
  obj <- f$obj + constraintL2(c(C1_a0 = 1, C2_a3 = 2), sigma = 1)
  d <- trust(obj, f$pars, rinit = 1, rmax = 10)
  s <- trust(obj, f$pars, rinit = 1, rmax = 10, sparse = TRUE)
  expect_true(s$converged)
  expect_equal(s$value, d$value, tolerance = 1e-8)
  expect_equal(s$argument[names(d$argument)], d$argument, tolerance = 1e-5)
})

test_that("sparse Hessians of the constraints equal the dense ones, also through a chain", {
  outer <- c(u = 0.3, v = -0.2, w = 1.1, z = 0.5)
  inner <- c(a = 0.3 * -0.2, b = exp(1.1), c = 0.5, s = 0.1)
  dP <- matrix(0, 4, 5, dimnames = list(names(inner), c(names(outer), "s")))
  dP["a", c("u", "v")] <- c(-0.2, 0.3); dP["b", "w"] <- exp(1.1); dP["c", "z"] <- 1; dP["s", "s"] <- 1
  dP2 <- array(0, c(4, 5, 5), dimnames = list(names(inner), colnames(dP), colnames(dP)))
  dP2["a", "u", "v"] <- dP2["a", "v", "u"] <- 1; dP2["b", "w", "w"] <- exp(1.1)
  chained <- as.parvec(inner, deriv = dP, deriv2 = dP2)
  plain <- as.parvec(inner)
  terms <- list(constraintL2(c(a = 0.1, b = 2), sigma = c(a = "s", b = "s")),
                constraintL2(c(a = 0.1, c = 1), sigma = 0.5),
                constraintCauchy(c(b = 2, c = 0), sigma = 1.5))
  for (f in terms) for (p in list(chained, plain)) {
    dense <- f(p, deriv2 = TRUE)
    sp <- f(p, deriv2 = TRUE, hessian = "sparse")
    expect_s4_class(sp$hessian, "sparseMatrix")
    nm <- names(dense$gradient)
    expect_equal(sp$gradient[nm], dense$gradient)
    expect_equal(as.matrix(sp$hessian)[nm, nm], dense$hessian[nm, nm])
  }
})

test_that("vcov of a sparse Hessian equals the dense one, also for a block and when singular", {
  f <- .sparseFixture()
  obj <- f$obj + constraintL2(c(C1_a0 = 1, C2_a3 = 2), sigma = 1)
  dense <- obj(f$pars)
  sp <- obj(f$pars, hessian = "sparse")
  nm <- names(f$pars)
  expect_equal(vcov(sp)[nm, nm], vcov(dense)[nm, nm])
  blk <- c("C2_a1", "C1_a2")
  expect_equal(vcov(sp, which = blk, fixed = "C1_a0"), vcov(dense, fixed = "C1_a0")[blk, blk])
  fit <- trust(obj, f$pars, rinit = 1, rmax = 10, sparse = TRUE)
  expect_equal(vcov(fit, which = blk), vcov(obj(fit$argument), which = blk), tolerance = 1e-6)
  sing <- objlist(0, setNames(numeric(2), c("a", "b")),
                  Matrix::sparseMatrix(i = c(1, 1, 2, 2), j = c(1, 2, 1, 2), x = 2, dims = c(2, 2),
                                       dimnames = list(c("a", "b"), c("a", "b"))))
  expect_equal(unname(vcov(sing)), matrix(0.25, 2, 2))
  # semidefinite matrices that sparse Cholesky factorizes with a pivot at rounding level
  for (seed in 1:3) {
    set.seed(seed)
    A <- matrix(rnorm(12), 4, 3)
    H <- A %*% t(A)
    nm4 <- paste0("p", 1:4)
    semi <- objlist(0, setNames(numeric(4), nm4), Matrix::Matrix(H, sparse = TRUE, dimnames = list(nm4, nm4)))
    expect_equal(unname(vcov(semi)), MASS::ginv(H / 2), tolerance = 1e-8)
  }
})

test_that("mstrust and profile run on sparse Hessians", {
  f <- .sparseFixture()
  obj <- f$obj + constraintL2(c(C1_a0 = 1, C2_a3 = 2), sigma = 1)
  d <- trust(obj, f$pars, rinit = 1, rmax = 10)
  set.seed(1)
  ms <- mstrust(obj, f$pars, fits = 2, cores = 1, sparse = TRUE, sd = 0.5)
  best <- as.parframe(ms)
  expect_equal(best$value[1], d$value, tolerance = 1e-6)

  for (W in c("hessian", "identity")) {
    dense <- profile(obj, d$argument, "C1_a2", limits = c(-1, 1), algoControl = list(W = W))
    sp <- profile(obj, d$argument, "C1_a2", limits = c(-1, 1), algoControl = list(W = W),
                  hessian = "sparse")
    expect_equal(sp$value, dense$value, tolerance = 1e-6)
    expect_equal(sp$constraint, dense$constraint, tolerance = 1e-6)
  }
})

test_that("a singular sparse path equation gives the profile and the warning of the dense one", {
  set.seed(18)
  J <- matrix(rnorm(24), 6, 4, dimnames = list(NULL, c("a", "b", "c", "d")))
  J[, "b"] <- J[, "a"] * exp(rnorm(1))
  r <- rnorm(6)
  obj <- function(pars, fixed = NULL, hessian = NULL, ...) {
    res <- as.vector(J %*% c(pars, fixed)[colnames(J)] - r)
    free <- names(pars)
    g <- setNames(2 * as.vector(crossprod(J, res)), colnames(J))
    H <- (2 * crossprod(J))[free, free, drop = FALSE]
    objlist(sum(res^2), g[free], if (identical(hessian, "sparse")) Matrix::Matrix(H, sparse = TRUE) else H)
  }
  p0 <- trust(obj, c(a = 0.5, b = 0, c = 0, d = 0), rinit = 1, rmax = 10)$argument
  run <- function(...) {
    msg <- character()
    pr <- withCallingHandlers(profile(obj, p0, "d", limits = c(-1, 1), cores = 1, ...),
                              warning = function(w) {
                                msg <<- c(msg, conditionMessage(w))
                                invokeRestart("muffleWarning")
                              })
    list(pr = pr, msg = msg)
  }
  dense <- run(sparse = FALSE)
  sp <- run(hessian = "sparse")
  expect_match(dense$msg, "singular at")
  expect_equal(sp$msg, dense$msg)
  expect_equal(sp$pr$constraint, dense$pr$constraint, tolerance = 1e-8)
  expect_equal(sp$pr$value, dense$pr$value, tolerance = 1e-8)

  # above 5000 parameters the singular path equation is an error that names the pair
  nm <- c(paste0("x", 1:5000), "a", "b")
  big <- function(pars, fixed = NULL, hessian = NULL, ...) {
    p <- c(pars, fixed)[nm]
    u <- p[["a"]] + p[["b"]] - 1
    g <- c(2 * (p[1:5000] - 1), a = 2.6 * u, b = 2.6 * u)
    H <- Matrix::bdiag(Matrix::Diagonal(5000, 2), Matrix::Matrix(2.6, 2, 2))
    dimnames(H) <- list(nm, nm)
    free <- names(pars)
    objlist(sum((p[1:5000] - 1)^2) + 1.3 * u^2, g[free], H[free, free])
  }
  expect_error(profile(big, setNames(c(rep(1, 5000), 0.4, 0.6), nm), "x1", limits = c(-1, 1),
                       cores = 1, hessian = "sparse"),
               "non-identifiable parameters: [ab]\\.")
})

test_that("trust(sparse = \"auto\") takes the sparse path above the thresholds, from the first evaluation", {
  f <- .sparseFixture()
  obj <- f$obj + constraintL2(c(C1_a0 = 1, C2_a3 = 2), sigma = 1)
  calls <- 0
  counted <- function(...) { calls <<- calls + 1; obj(...) }
  dense <- trust(obj, f$pars, rinit = 1, rmax = 10)
  expect_true(is.matrix(dense$hessian))
  forced <- trust(counted, f$pars, rinit = 1, rmax = 10, sparse = TRUE)
  n <- calls; calls <- 0
  auto <- withr::with_options(list(dMod2.sparseAuto = list(parameters = 0, fill = 1)),
                              trust(counted, f$pars, rinit = 1, rmax = 10))
  expect_s4_class(auto$hessian, "sparseMatrix")
  expect_equal(calls, n)
  expect_equal(auto$value, dense$value, tolerance = 1e-8)
  # a fill above the threshold or a quasi-Newton Hessian keep the dense path; bounds take the
  # same steps on both
  withr::with_options(list(dMod2.sparseAuto = list(parameters = 0, fill = 0.01)), {
    expect_true(is.matrix(trust(obj, f$pars, rinit = 1, rmax = 10)$hessian))
  })
  withr::with_options(list(dMod2.sparseAuto = list(parameters = 0, fill = 1)), {
    bounded <- trust(obj, f$pars, parupper = c(C1_a0 = 10))
    expect_s4_class(bounded$hessian, "sparseMatrix")
    expect_equal(bounded$argument, trust(obj, f$pars, parupper = c(C1_a0 = 10), sparse = FALSE)$argument,
                 tolerance = 1e-10)
    expect_true(is.matrix(trust(obj, f$pars, hessianMethod = "bfgs")$hessian))
    pr <- profile(obj, dense$argument, "C1_a2", limits = c(-1, 1))
    expect_equal(pr$value, suppressWarnings(profile(obj, dense$argument, "C1_a2", limits = c(-1, 1),
                                                    sparse = FALSE))$value, tolerance = 1e-6)
  })
  plain <- function(x, ...) list(value = sum(x^2), gradient = 2 * x, hessian = diag(2, length(x)))
  expect_error(trust(plain, f$pars, sparse = TRUE), "dense Hessian")
  expect_error(trust(obj, f$pars, sparse = "yes"), "TRUE, FALSE")
})

test_that("normT, censored rows and datapointL2 have the sparse Hessian of the dense one", {
  f <- .sparseFixture()
  grid <- function(cn) data.frame(name = "y", time = 0:3, row.names = paste0(cn, "_a", 0:3))
  x <- Xd(grid("C1"), condition = "C1") + Xd(grid("C2"), condition = "C2")
  set.seed(3)
  data <- data.frame(name = "y", time = rep(0:3, 2), value = rnorm(8, 2), sigma = 0.3,
                     condition = rep(c("C1", "C2"), each = 4))
  censored <- as.datalist(transform(data, lloq = 1.8))
  data <- as.datalist(data)
  pars <- c(f$pars, target = 1.5)
  terms <- function(d) list(normT(d, x, df = 4), normL2(d, x),
                            normL2(d, x) + datapointL2("y", 2, "target", sigma = 0.2, condition = "C1"))
  for (obj in terms(censored)) {
    dense <- obj(pars, deriv = TRUE)
    sp <- obj(pars, deriv = TRUE, hessian = "sparse")
    expect_s4_class(sp$hessian, "sparseMatrix")
    nm <- names(dense$gradient)
    expect_equal(sp$value, dense$value)
    expect_equal(sp$gradient[nm], dense$gradient)
    expect_equal(as.matrix(sp$hessian)[nm, nm], dense$hessian[nm, nm])
  }
  # the optima without censored rows, whose parameters of their own run off to -Inf
  for (obj in terms(data)) {
    start <- if ("target" %in% attr(obj, "parameters")) pars else f$pars
    a <- trust(obj, start, rinit = 1, rmax = 10, sparse = FALSE)
    b <- trust(obj, start, rinit = 1, rmax = 10, sparse = TRUE)
    expect_equal(b$value, a$value, tolerance = 1e-8)
    expect_equal(b$argument[names(a$argument)], a$argument, tolerance = 1e-5)
  }
})

test_that("a sum of objectives holds the terms, not the frame it is built in", {
  f <- .sparseFixture()
  prior <- constraintL2(c(C1_a0 = 1), sigma = 1)
  build <- function() {
    ballast <- numeric(1e6)
    f$obj + prior
  }
  obj <- build()
  size <- length(serialize(obj, NULL))
  expect_lt(size, length(serialize(f$obj, NULL)) + length(serialize(prior, NULL)) + 1e5)
})

# One step of both curvature models of trust() on H plus the low-rank part `lr`.
.modelSteps <- function(H, g, r, lr = NULL) {
  nm <- paste0("p", seq_along(g))
  H <- Matrix::Matrix(as.matrix(H), sparse = TRUE, dimnames = list(nm, nm))
  if (!is.null(lr)) rownames(lr$U) <- nm
  k <- dMod2:::.kernelHessian(list(hessian = H, lowrank = lr), nm)$kernel
  full <- as.matrix(H) + if (is.null(lr)) 0 else lr$U %*% (lr$d * t(lr$U))
  list(sparse = dMod2:::trust_model_step(k, full, g, r, TRUE), dense = dMod2:::trust_model_step(k, full, g, r, FALSE),
       B = full)
}

test_that("both models solve the subproblem of an indefinite sum with a low-rank part, Newton without it", {
  H <- Matrix::Diagonal(x = c(1, 2, 3, 4))
  lr <- list(U = cbind(c(1, 1, 0, 0), c(0, 0.5, 1, 0)), d = c(-2, 0.5))
  g <- c(1, -1, 0.5, 2)
  st <- .modelSteps(H, g, 0.7, lr)
  expect_lt(min(eigen(st$B, symmetric = TRUE)$values), 0)
  expect_equal(sqrt(sum(st$sparse$p^2)), 0.7, tolerance = 1e-10)
  expect_equal(st$sparse$p, st$dense$p, tolerance = 1e-9)
  expect_equal(st$sparse$pred, st$dense$pred, tolerance = 1e-9)
  plain <- .modelSteps(H, g, 10)
  expect_true(plain$sparse$newton && plain$dense$newton)
  expect_equal(plain$sparse$p, -g / c(1, 2, 3, 4), tolerance = 1e-12)
})

test_that("both models step to the boundary along a flat direction, zero at a zero gradient, and take tiny gradients", {
  H <- Matrix::Diagonal(x = c(1, 1e-200))
  st <- .modelSteps(H, c(1, 1), 1)
  expect_equal(sqrt(sum(st$sparse$p^2)), 1, tolerance = 1e-10)
  expect_equal(st$sparse$p, st$dense$p, tolerance = 1e-9)
  expect_equal(.modelSteps(H, c(0, 0), 1)$sparse$p, c(0, 0))
  for (g in list(c(1e-160, 1e-160), c(1e-170, 1e-170))) {
    tiny <- .modelSteps(Matrix::Diagonal(x = c(1, 2)), g, 1e-170)
    expect_equal(dMod2:::.norm2(tiny$sparse$p) / 1e-170, 1, tolerance = 1e-10)
  }
})

test_that("the sparse model takes a low-rank part whose weights differ by many orders", {
  setTimeLimit(elapsed = 10, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf))
  st <- .modelSteps(Matrix::Diagonal(x = c(1, 2, 3)), c(1, 1, 1), 10, list(U = cbind(c(1, 0, 0), c(0, 1, 0)), d = c(1, -1e-13)))
  expect_true(st$sparse$newton)
  expect_equal(st$sparse$p, -1 / c(2, 2 - 1e-13, 3), tolerance = 1e-12)
})

test_that("the sparse and the dense model take the same steps, indefinite and with a low-rank part", {
  # an objective indefinite away from its minimum, with a sparse Hessian and a dense part of rank one
  n <- 30
  f <- function(x, deriv = TRUE, hessian = TRUE, deriv2 = FALSE, ...) {
    a <- x[-n]; b <- x[-1]; r <- b - a^2
    v <- sum(r^2) + sum((1 - a)^2) + 0.1 * sum(x^4) - sum(x^2) + 50 * (sum(x) - 3)^2
    g <- numeric(n)
    g[-n] <- g[-n] - 4 * a * r - 2 * (1 - a); g[-1] <- g[-1] + 2 * r
    g <- setNames(g + 0.4 * x^3 - 2 * x + 100 * (sum(x) - 3), names(x))
    d <- numeric(n); d[-n] <- d[-n] + 12 * a^2 - 4 * b + 2; d[-1] <- d[-1] + 2; d <- d + 1.2 * x^2 - 2
    S <- Matrix::sparseMatrix(i = c(seq_len(n), seq_len(n - 1), seq_len(n - 1) + 1),
                              j = c(seq_len(n), seq_len(n - 1) + 1, seq_len(n - 1)),
                              x = c(d, -4 * a, -4 * a), dims = c(n, n), dimnames = list(names(x), names(x)))
    if (!identical(hessian, "sparse")) return(objlist(v, g, as.matrix(S) + 100))
    out <- objlist(v, g, S)
    out$lowrank <- list(U = matrix(1, n, 1, dimnames = list(names(x), NULL)), d = 100)
    out
  }
  x0 <- setNames(seq(-1, 1, length.out = n), paste0("x", seq_len(n)))
  for (bounds in list(NULL, list(lower = setNames(rep(-0.5, n), names(x0)), upper = setNames(rep(0.8, n), names(x0))))) {
    run <- function(sparse) suppressWarnings(trust(f, x0, rinit = 1, rmax = 10, sparse = sparse, blather = TRUE,
                                                   parlower = bounds$lower, parupper = bounds$upper))
    d <- run(FALSE); s <- run(TRUE)
    expect_identical(s$steptype, d$steptype)
    expect_equal(s$valpath, d$valpath, tolerance = 1e-10)
    expect_equal(s$argument, d$argument, tolerance = 1e-10)
  }
})
