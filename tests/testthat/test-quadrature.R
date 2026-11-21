# Quadrature stack, bottom-up from the sparse grids to the quadrature objective and the certificate.


# ---- sparse grids --------------------------------------------------------

# Analytic moments of the weight exp(-z^2): int z^k exp(-z^2) dz = 0 if k odd,
# (2m)!/(4^m m!) * sqrt(pi) if k = 2m.
gh_moment_1d <- function(k) {
  if (k %% 2L == 1L) return(0)
  m <- k / 2L
  factorial(2 * m) / (4^m * factorial(m)) * sqrt(pi)
}

test_that("sparseGridGK is exact to total degree 2 (level - K) + 1 with mass pi^(K/2)", {
  for (K in 1:3) for (level in K:(K + 3L)) {
    g <- dMod2:::sparseGridGK(as.integer(K), as.integer(level))
    expect_equal(sum(g$weights), pi^(K / 2), tolerance = 1e-12)
    deg <- 2L * (level - K) + 1L
    a <- as.matrix(expand.grid(rep(list(0:deg), K)))
    a <- a[rowSums(a) <= deg, , drop = FALSE]
    for (r in seq_len(nrow(a))) {
      est <- sum(g$weights * apply(g$nodes^matrix(a[r, ], nrow(g$nodes), K, byrow = TRUE), 1, prod))
      expect_equal(est, prod(vapply(a[r, ], gh_moment_1d, 0)), tolerance = 1e-10,
                   info = sprintf("K = %d, level = %d, a = %s", K, level, paste(a[r, ], collapse = ",")))
    }
  }
  # level K + 1 is the 2K + 1 point rule
  expect_equal(nrow(dMod2:::sparseGridGK(4L, 5L)$nodes), 9L)
  expect_error(dMod2:::sparseGridGK(3L, 2L), "below dimension")
})

test_that("the nodes integrate a normal of the posterior precision exactly", {
  set.seed(42)
  K <- 3L
  mode <- c(0.4, -0.2, 0.7)
  A <- matrix(rnorm(K * K), K, K)
  P <- crossprod(A) + diag(K)
  Sig <- solve(P)
  for (level in (K + 1L):5L) {
    qn <- dMod2:::.quadratureNodes(mode, P, level)
    W <- qn$signs * exp(qn$logW)
    d <- sweep(qn$eta, 2, mode)
    dens <- exp(-K / 2 * log(2 * pi) + as.numeric(determinant(P)$modulus) / 2 - rowSums((d %*% P) * d) / 2)
    expect_equal(sum(W * dens), 1, tolerance = 1e-8)
    expect_equal(colSums(W * dens * qn$eta), mode, tolerance = 1e-8)
    expect_equal(crossprod(d, W * dens * d), Sig, tolerance = 1e-7)
  }
  # pruning drops the nodes of smallest share
  qn <- dMod2:::.quadratureNodes(mode, P, K + 3L, pruneTol = 8)
  d <- sweep(qn$eta, 2, mode)
  dens <- exp(-K / 2 * log(2 * pi) + as.numeric(determinant(P)$modulus) / 2 - rowSums((d %*% P) * d) / 2)
  expect_lt(nrow(qn$eta), nrow(dMod2:::.quadratureNodes(mode, P, K + 3L)$eta))
  expect_equal(sum(qn$signs * exp(qn$logW) * dens), 1, tolerance = 1e-2)
})


# A one-eta model y = exp(mu + eta) per group, observed once per condition.
.quad_models <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    withr::local_preserve_seed()
    dir <- file.path(tempdir(), "quadrature_models")
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    withr::local_dir(dir)
    g <- Y(c(y = "intercept"), f = NULL, parameters = "intercept", compile = FALSE,
           modelname = "quad_obs", deriv2 = TRUE)
    groups <- c("s1", "s2", "s3")
    tab <- data.frame(ID = rep(groups, each = 2), w = rep(c(-1, 0, 1.5), each = 2),
                      row.names = paste0(rep(groups, each = 2), c("a", "b")))
    tr <- branch(define(NULL, "intercept <- exp(mu)"), table = tab)
    tr <- insert(tr, mu ~ 1 + (1 | ID))
    p <- P(tr, method = "explicit", compile = FALSE, modelname = "quad_p", deriv2 = TRUE)
    # a group covariate on mu and the parameter s of an error model sigma = exp(s)
    tw <- branch(define(NULL, "x <- x", x = c("intercept", "s")), table = tab)
    tw <- insert(tw, "intercept <- exp(mu)")
    tw <- insert(tw, mu ~ w + (1 | ID))
    pw <- P(tw, method = "explicit", compile = FALSE, modelname = "quad_pw", deriv2 = TRUE)
    err <- Y(c(y = "exp(s)"), g, attachInput = FALSE, compile = FALSE, modelname = "quad_err",
             deriv2 = TRUE)
    compile(g, p, pw, err, output = "quadrature_models", cores = 2)
    values <- c(2.0, 2.3, 1.6, 1.5, 2.6, 2.4)
    data <- as.datalist(data.frame(name = "y", time = 0, sigma = 0.3, value = values,
                                   condition = rownames(tab)))
    dataE <- as.datalist(data.frame(name = "y", time = 0, sigma = NA, value = values,
                                    condition = rownames(tab)))
    cache <<- list(g = g, p = p, pw = pw, data = data, obj = normL2(data, g * Xt() * p),
                   objE = normL2(dataE, g * Xt() * pw, errmodel = err))
    cache
  }
})

test_that("evalConditionResidual sums to normL2", {
  m <- .quad_models()
  pars <- c(mu = 0.7, eta_mu_s1 = 0.1, eta_mu_s2 = -0.1, eta_mu_s3 = 0.2)
  prediction <- (m$g * Xt() * m$p)(times = 0, pars = pars, deriv = TRUE, conditions = names(m$data))
  agg <- Reduce(`+`, lapply(names(m$data), function(cn)
    evalConditionResidual(dataI = m$data[[cn]], predictionI = prediction[[cn]],
                          cn = cn, deriv = TRUE)))
  res <- m$obj(pars)
  expect_equal(agg$value, res$value, tolerance = 1e-10)
  expect_equal(agg$gradient, res$gradient, tolerance = 1e-10)
  expect_equal(agg$hessian, res$hessian, tolerance = 1e-10)
})

test_that("the quadrature of a group over the group conditions recovers the integral", {
  m <- .quad_models()
  d <- dMod2:::.emDesign(m$obj, randomEffects(m$p))
  psi <- c(mu = 0.7, omega_mu_mu = log(0.3))
  truth <- -2 * sum(vapply(d$groups, function(g) {
    y <- vapply(d$conditions[[g]], function(cn) m$data[[cn]]$value, 0)
    f <- function(e) vapply(e, function(e0) exp(sum(dnorm(y, exp(0.7 + e0), 0.3, log = TRUE)) +
                                                dnorm(e0, 0, 0.3, log = TRUE)), 0)
    log(integrate(f, -Inf, Inf, rel.tol = 1e-12)$value)
  }, 0))
  # one eta: the levels of the 3, 9, 19 and 35 point rules
  q <- dMod2:::.quadratureObjective(d, psi, NULL, list())
  for (k in 1:4) {
    q$rebuild(psi, c(2, 4, 9, 16)[k])
    v <- q$value(psi["mu"], deriv = FALSE)$value
    expect_equal(v, truth, tolerance = c(1e-2, 1e-4, 1e-7, 1e-10)[k])
  }
  # by tolerance: each group at the lowest level whose next one changes -2 log L by tol / G
  q <- dMod2:::.quadratureObjective(d, psi, NULL, list())
  q$rebuild(psi, 1, tol = 1e-6)
  expect_true(all(q$state$levelError <= 1e-6 / 3))
  expect_true(all(q$state$levels == 9L))
  expect_lt(abs(q$value(psi["mu"], deriv = FALSE)$value - truth), 1e-6)
  q$rebuild(psi, 1, tol = 1)
  expect_true(all(q$state$levels == 9L))
})

test_that("on frozen nodes the quadrature has the exact gradient and Hessian, Omega included", {
  m <- .quad_models()
  d <- dMod2:::.emDesign(m$obj, randomEffects(m$p))
  psi <- c(mu = 0.7, omega_mu_mu = log(0.3))
  q <- dMod2:::.quadratureObjective(d, psi, NULL, list())
  q$rebuild(psi, 5)
  at <- q$value(psi, deriv2 = TRUE)
  f <- function(v) q$value(setNames(v, names(psi)), deriv = FALSE)$value
  gr <- function(v) q$value(setNames(v, names(psi)), deriv2 = TRUE)$gradient
  expect_equal(unname(at$gradient), numDeriv::grad(f, unname(psi)), tolerance = 1e-7)
  expect_equal(unname(at$hessian), numDeriv::jacobian(gr, unname(psi)), tolerance = 1e-6)
})

test_that("on nodes frozen in phi the referenced fixed effects, the error model and Omega have exact derivatives", {
  m <- .quad_models()
  d <- dMod2:::.emDesign(m$objE, randomEffects(m$pw))
  expect_equal(dMod2:::.referenceShift(d)$pars, c("mu", "mu_w"))
  psi <- c(mu = 0.7, mu_w = 0.1, s = log(0.3), omega_mu_mu = log(0.3))
  q <- dMod2:::.quadratureObjective(d, psi, NULL, list())
  q$rebuild(psi, 4)
  at <- psi + c(0.05, -0.03, 0.02, 0.04)
  ev <- q$value(at, deriv2 = TRUE)
  f <- function(v) q$value(setNames(v, names(psi)), deriv = FALSE)$value
  gr <- function(v) q$value(setNames(v, names(psi)), deriv2 = TRUE)$gradient
  expect_equal(unname(ev$gradient), numDeriv::grad(f, unname(at)), tolerance = 1e-7)
  expect_equal(unname(ev$hessian), numDeriv::jacobian(gr, unname(at)), tolerance = 1e-6)
  # at the rebuild point the frozen value is the quadrature there
  expect_equal(q$value(psi, deriv = FALSE)$value,
               -2 * sum(vapply(seq_along(d$groups), function(i) {
                 n <- q$state$nodes[[i]]
                 v <- q$integrand(d$groups[i], list(n$eta), psi)[[1]] + n$logW
                 log(sum(n$signs * exp(v)))
               }, 0)), tolerance = 1e-10)
})

test_that("a fit on capped levels ends where the levels resolve, at the best iterate", {
  m <- .quad_models()
  fit <- suppressMessages(em(m$objE, randomEffects(m$pw), c(mu = 0.7, mu_w = 0.1, s = log(0.3)),
                             method = "quadrature", verbose = FALSE,
                             control = list(quadrature = list(maxLevel = 2, tol = 1e-12, tolDecrement = 1e-12,
                                                              certificate = FALSE))))
  tr <- fit$stageTrace
  expect_true(fit$converged)
  expect_lt(nrow(tr), 10L)
  # every group is held at level 2 with an error above tol / G
  expect_true(all(tr$capped == 3L))
  lv <- attr(tr, "levels")
  expect_equal(lv$level, rep(2L, 3))
  expect_true(all(lv$capped & lv$levelError > 1e-12))
  # the decrement stays above tolDecrement: the stop is a rise on unchanged levels
  expect_identical(fit$stopReason, "resolution")
  expect_gt(tr$OFV[nrow(tr)], tr$OFV[nrow(tr) - 1L])
  # the best iterate, not the last
  expect_equal(fit$value, min(tr$OFV), tolerance = 1e-10)
  expect_lt(fit$value, tr$OFV[nrow(tr)])
})

test_that("a fit returns the best iterate when the rebuilds run out", {
  m <- .quad_models()
  fit <- suppressMessages(em(m$objE, randomEffects(m$pw), c(mu = 0.7, mu_w = 0.1, s = log(0.3)),
                             method = "quadrature", verbose = FALSE,
                             control = list(quadrature = list(maxLevel = 2, tol = 1e-12, tolDecrement = 0,
                                                              maxEcmPerStage = 2, certificate = FALSE))))
  tr <- fit$stageTrace
  expect_false(fit$converged)
  expect_identical(fit$stopReason, "iterlim")
  expect_equal(nrow(tr), 2L)
  expect_equal(fit$value, min(tr$OFV), tolerance = 1e-10)
})

test_that("a fit on uncapped levels is stationary at tolDecrement", {
  m <- .quad_models()
  fit <- suppressMessages(em(m$objE, randomEffects(m$pw), c(mu = 0.7, mu_w = 0.1, s = log(0.3)),
                             method = "quadrature", verbose = FALSE,
                             control = list(quadrature = list(tol = 1e-8, tolDecrement = 1e-10,
                                                              certificate = FALSE))))
  tr <- fit$stageTrace
  expect_true(fit$converged)
  expect_identical(fit$stopReason, "decrement")
  expect_true(all(tr$capped == 0L))
  expect_lt(tr$decrement[nrow(tr)], 1e-10)
  expect_lt(max(abs(fit$gradient)), 1e-5)
})

test_that("levels over the node budget are rejected and the error comes from the level below", {
  # nodes per group at K = 20 for levels K + 1 to K + 3: the default budget ends at K + 2
  expect_equal(vapply(21:23, function(l) dMod2:::.levelNodes(20L, l), 0L), c(41L, 801L, 10041L))
  m <- .quad_models()
  d <- dMod2:::.emDesign(m$obj, randomEffects(m$p))
  psi <- c(mu = 0.7, omega_mu_mu = log(0.3))
  # one eta: levels 2, 4 and 9 have 3, 9 and 19 nodes
  logL <- function(level) {
    q <- dMod2:::.quadratureObjective(d, psi, NULL, list())
    q$rebuild(psi, level)
    vapply(seq_along(d$groups), function(i) {
      n <- q$state$nodes[[i]]
      log(sum(n$signs * exp(q$integrand(d$groups[i], list(n$eta), psi)[[1]] + n$logW)))
    }, 0)
  }
  q <- dMod2:::.quadratureObjective(d, psi, NULL, list(maxNodes = 9))
  q$rebuild(psi, 1, tol = 1e-6)
  st <- q$state
  expect_equal(st$levels, rep(4L, 3))
  expect_true(all(st$capped))
  expect_true(all(vapply(st$nodes, function(n) nrow(n$eta), 0L) <= 9L))
  expect_equal(st$levelError, 2 * abs(logL(4) - logL(2)), tolerance = 1e-8)
  expect_false(q$allowed(9L))
  # a requested level over the budget falls to the highest one within it
  q <- dMod2:::.quadratureObjective(d, psi, NULL, list(maxNodes = 9))
  q$rebuild(psi, 9)
  expect_equal(vapply(q$state$nodes, function(n) nrow(n$eta), 0L), rep(9L, 3))
  expect_true(all(q$state$capped))
  # within the budget the levels are those without it
  q <- dMod2:::.quadratureObjective(d, psi, NULL, list(maxNodes = 35))
  q$rebuild(psi, 1, tol = 1e-6)
  expect_equal(q$state$levels, rep(9L, 3))
  expect_false(any(q$state$capped))
})

test_that("the certificate at the estimate agrees with the quadrature and is reproducible", {
  m <- .quad_models()
  d <- dMod2:::.emDesign(m$objE, randomEffects(m$pw))
  psi <- c(mu = 0.7, mu_w = 0.1, s = log(0.3), omega_mu_mu = log(0.3))
  q <- dMod2:::.quadratureObjective(d, psi, NULL, list())
  q$rebuild(psi, 4)
  cert <- dMod2:::.quadratureCertificate(d, q, psi, NULL, list(points = 256, seed = 3))
  expect_equal(cert$group, d$groups)
  expect_false(any(cert$flag))
  expect_true(all(abs(cert$quadrature - cert$importance) < 3 * cert$se + 1e-3))
  expect_identical(cert, dMod2:::.quadratureCertificate(d, q, psi, NULL, list(points = 256, seed = 3)))
})

test_that("the Newton decrement of a singular Hessian counts the flat directions linearly", {
  # positive definite: g' H^-1 g / 2
  ev <- list(gradient = c(1, 2), hessian = diag(c(2, 4)))
  expect_equal(dMod2:::.newtonDecrement(ev), (1 / 2 + 4 / 4) / 2)
  # a flat direction without gradient, as at a variance that runs to 0, adds nothing
  ev <- list(gradient = c(1, 0), hessian = diag(c(2, 0)))
  expect_equal(dMod2:::.newtonDecrement(ev), 1 / 4)
  # a gradient along a flat direction gains rmax times the gradient size
  ev <- list(gradient = c(1, 1e-3), hessian = diag(c(2, 0)))
  expect_equal(dMod2:::.newtonDecrement(ev, rmax = 10), 1 / 4 + 1e-2)
})
