# Laplace and spike-and-slab random effects: the integrals against numerical integration, the
# exact gradient of the approximate marginal, its exactness on a linear model, and fits.

test_that("the normal-Laplace integral and its moments match a numerical integral", {
  for (case in list(c(a = 3, m = 0.4, lambda = 2), c(a = 0.5, m = -1.2, lambda = 0.7))) {
    a <- case[["a"]]; m <- case[["m"]]; lam <- case[["lambda"]]
    f <- function(e, k = 0) e^k * lam / 2 * exp(-lam * abs(e) - a / 2 * (e - m)^2)
    I <- integrate(f, -Inf, Inf, rel.tol = 1e-12)$value
    nl <- dMod2:::.normalLaplace(a, m, lam)
    expect_equal(nl$logI, log(I), tolerance = 1e-10)
    expect_equal(nl$mean, integrate(f, -Inf, Inf, k = 1, rel.tol = 1e-12)$value / I, tolerance = 1e-8)
    expect_equal(nl$dm, numDeriv::grad(function(x) dMod2:::.normalLaplace(a, x, lam)$logI, m), tolerance = 1e-7)
    expect_equal(nl$da, numDeriv::grad(function(x) dMod2:::.normalLaplace(x, m, lam)$logI, a), tolerance = 1e-7)
    expect_equal(nl$dlambda, numDeriv::grad(function(x) dMod2:::.normalLaplace(a, m, x)$logI, lam), tolerance = 1e-7)
  }
  # a strength far above the curvature: the integral tends to exp(-a m^2 / 2)
  nl <- dMod2:::.normalLaplace(2, 0.3, 1e5)
  expect_equal(nl$logI, -0.09, tolerance = 1e-4)
  expect_equal(nl$dm, numDeriv::grad(function(x) dMod2:::.normalLaplace(2, x, 1e5)$logI, 0.3), tolerance = 1e-7)
  expect_equal(nl$da, numDeriv::grad(function(x) dMod2:::.normalLaplace(x, 0.3, 1e5)$logI, 2), tolerance = 1e-7)
})

# A decay model y = A exp(-k dose) per condition, five groups with two conditions each, and the
# linear model y = a + b dose on the same design with a known sigma.
.lap <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    withr::local_preserve_seed()
    dir <- file.path(tempdir(), "emlaplace_models")
    dir.create(dir, showWarnings = FALSE)
    withr::local_dir(dir)
    # second derivatives forward over forward and forward over reverse
    .dm <- c("forward", "reverse", "forward-reverse")
    groups <- paste0("L", 1:5)
    tab <- data.frame(ID = rep(groups, each = 2), dose = rep(c(0.5, 2), 5),
                      row.names = paste0(rep(groups, each = 2), c("a", "b")))
    base <- function() {
      tr <- branch(define(NULL, "x <- x", x = c("k", "A", "s")), table = tab)
      tr <- define(tr, "dose <- d", d = dose)
      insert(tr, "x <- exp(x)", x = c("k", "A", "s"))
    }
    g <- Y(c(y = "A * exp(-k * dose)"), f = NULL, parameters = c("A", "k", "dose"), compile = FALSE,
           modelname = "lap_g", deriv2 = TRUE, derivMode = .dm)
    err <- Y(c(y = "s"), g, attachInput = FALSE, compile = FALSE, modelname = "lap_err", deriv2 = TRUE, derivMode = .dm)
    pc <- P(insert(base(), k + A ~ 1 + (1 | gr(ID, dist = laplace()))), method = "explicit",
            compile = FALSE, modelname = "lap_centre", deriv2 = TRUE, derivMode = .dm)
    pm <- P(insert(base(), k ~ 1 + (1 | gr(ID, dist = laplace())), A ~ 1 + (1 | ID)), method = "explicit",
            compile = FALSE, modelname = "lap_mixed", deriv2 = TRUE, derivMode = .dm)
    pn <- P(insert(base(), A ~ 1 + (1 | ID)), method = "explicit", compile = FALSE, modelname = "lap_normal",
            deriv2 = TRUE, derivMode = .dm)
    pl <- P(insert(base(), k ~ 1 + (1 | gr(ID, dist = laplace()))), method = "explicit", compile = FALSE,
            modelname = "lap_konly", deriv2 = TRUE, derivMode = .dm)
    ps <- P(insert(base(), k + A ~ 1 + (1 | gr(ID, dist = spikeSlab()))), method = "explicit",
            compile = FALSE, modelname = "lap_slab", deriv2 = TRUE, derivMode = .dm)
    psk <- P(insert(base(), k ~ 1 + (1 | gr(ID, dist = dirichletProcess()))), method = "explicit",
             compile = FALSE, modelname = "lap_dpcentre", deriv2 = TRUE, derivMode = .dm)
    psm <- P(insert(base(), k ~ 1 + (1 | gr(ID, dist = spikeSlab())), A ~ 1 + (1 | ID)), method = "explicit",
             compile = FALSE, modelname = "lap_slabmixed", deriv2 = TRUE, derivMode = .dm)
    pskm <- P(insert(base(), k ~ 1 + (1 | gr(ID, dist = dirichletProcess())), A ~ 1 + (1 | ID)),
              method = "explicit", compile = FALSE, modelname = "lap_slabclmixed", deriv2 = TRUE, derivMode = .dm)
    psk2 <- P(insert(base(), k + A ~ 1 + (1 | gr(ID, dist = dirichletProcess()))), method = "explicit",
              compile = FALSE, modelname = "lap_slabcltwo", deriv2 = TRUE, derivMode = .dm)
    linear <- function(dist, a = a ~ 1 + (1 | ID)) {
      tr <- branch(define(NULL, "x <- x", x = c("a", "b")), table = tab)
      tr <- define(tr, "dose <- d", d = dose)
      insert(tr, a, b ~ 1 + (1 | gr(ID, dist = dist)))
    }
    gl <- Y(c(y = "a + b * dose"), f = NULL, parameters = c("a", "b", "dose"), compile = FALSE,
            modelname = "lap_glin", deriv2 = TRUE, derivMode = .dm)
    pll <- P(linear(laplace()), method = "explicit", compile = FALSE, modelname = "lap_linlap", deriv2 = TRUE, derivMode = .dm)
    pls <- P(linear(spikeSlab()), method = "explicit", compile = FALSE, modelname = "lap_linslab", deriv2 = TRUE, derivMode = .dm)
    pld <- P(linear(dirichletProcess()), method = "explicit", compile = FALSE, modelname = "lap_lindp", deriv2 = TRUE, derivMode = .dm)
    pl2 <- P(linear(laplace(), a ~ 1 + (1 | gr(ID, dist = laplace()))), method = "explicit", compile = FALSE,
             modelname = "lap_lintwo", deriv2 = TRUE, derivMode = .dm)
    compile(g, err, pc, pm, pn, pl, ps, psk, psm, pskm, psk2, gl, pll, pls, pld, pl2, output = "emlaplace_models",
            cores = 2)
    set.seed(4)
    dk <- c(L1 = 0, L2 = 0, L3 = 0.8, L4 = 0.8, L5 = -0.6)
    rows <- lapply(rownames(tab), function(cn) {
      k <- 0.6 * exp(dk[[tab[cn, "ID"]]])
      data.frame(name = "y", time = c(0, 1, 2, 3), value = 10 * exp(-k * tab[cn, "dose"]) + rnorm(4, 0, 0.1),
                 sigma = NA_real_, condition = cn)
    })
    data <- as.datalist(do.call(rbind, rows))
    ea <- setNames(rnorm(5, 0, 0.5), groups)
    eb <- c(L1 = 0, L2 = 0, L3 = 0.4, L4 = -0.3, L5 = 0)
    rows <- lapply(rownames(tab), function(cn) {
      id <- tab[cn, "ID"]
      data.frame(name = "y", time = c(0, 1, 2, 3),
                 value = 1 + ea[[id]] + (0.5 + eb[[id]]) * tab[cn, "dose"] + rnorm(4, 0, 0.3),
                 sigma = 0.3, condition = cn)
    })
    cache <<- list(g = g, err = err, pc = pc, pm = pm, pn = pn, pl = pl, ps = ps, psk = psk, psm = psm,
                   pskm = pskm, psk2 = psk2, pl2 = pl2, data = data, init = c(k = log(0.6), A = log(10), s = log(0.1)),
                   gl = gl, pll = pll, pls = pls, pld = pld, tab = tab, linear = as.datalist(do.call(rbind, rows)))
    cache
  }
})

# The marginal of `p` at `theta` with every Laplace strength `lambda`, slab probability or
# concentration `share` and slab standard deviation `tau`; the decay data unless `obj` is given.
.marginalAt <- function(p, theta, lambda = NA, share = NA, tau = NA, grad = TRUE, obj = NULL,
                        coupling = "expectationConsistent", scores = FALSE, derivMode = NULL) {
  m <- .lap()
  obj <- obj %||% normL2(m$data, m$g * Xt() * p, errmodel = m$err)
  d <- dMod2:::.emDesign(obj, randomEffects(p))
  lay <- dMod2:::.sparseLayout(d)
  hyp <- lapply(list(lambda = lambda, prob = share, concentration = share, tau = tau),
                function(v) setNames(rep(v, nrow(lay)), lay$stem))
  prior <- if (d$normal$K) dMod2:::.normalPrior(d)
  warm <- matrix(0, length(d$groups), ncol(d$etaNames))
  dMod2:::.laplaceMarginal(d, dMod2:::.conditionObjectives(d), prior, theta, hyp, names(theta),
                           dMod2:::.sparseHypers(lay), NULL, warm, grad = grad,
                           control = list(iterlim = 300, coupling = coupling, derivMode = derivMode), scores = scores)
}

# The slope of a marginal: its inner iterations leave noise of about 1e-10 in the value, so the
# difference steps stay above 1e-3.
.slope <- function(f, x) numDeriv::grad(f, x, method.args = list(eps = 1e-2, d = 1e-2, r = 4))

test_that("the approximate marginal has the exact gradient", {
  m <- .lap()
  theta <- m$init + c(0.05, -0.02, 0.1)
  lam <- 3
  r <- .marginalAt(m$pc, theta, lam)
  f <- function(v) .marginalAt(m$pc, setNames(v, names(theta)), lam, grad = FALSE)$value
  expect_equal(unname(r$gradient[names(theta)]), .slope(f, unname(theta)), tolerance = 1e-5)
  fl <- function(l) .marginalAt(m$pc, theta, l, grad = FALSE)$value
  expect_equal(unname(r$gradient[length(theta) + 1L]), .slope(fl, lam), tolerance = 1e-5)
})

test_that("the marginal alone needs no second-order sensitivities under forward over reverse", {
  m <- .lap()
  theta <- m$init + c(0.05, -0.02, 0.1)
  full <- .marginalAt(m$pc, theta, 3, derivMode = "forward-reverse")
  real <- dMod2:::.predictMany
  local_mocked_bindings(.predictMany = function(...) {
    if (isTRUE(list(...)$deriv2)) stop("second-order sensitivities requested")
    real(...)
  }, .package = "dMod2")
  alone <- .marginalAt(m$pc, theta, 3, grad = FALSE, derivMode = "forward-reverse")
  expect_equal(alone$value, full$value, tolerance = 1e-10)
})

test_that("forward over reverse gives the gradient and the model Hessian rows of forward over forward", {
  m <- .lap()
  th <- m$init + c(0.05, -0.02, 0.1)
  cases <- list(list(p = m$pc, theta = th, lambda = 3),
                list(p = m$pm, theta = c(th, omega_A_A = log(0.2)), lambda = 3),
                list(p = m$ps, theta = th, share = 0.4, tau = 0.5),
                list(p = m$psk, theta = th, share = 2, tau = 0.5))
  for (cs in cases) {
    at <- function(mode) .marginalAt(cs$p, cs$theta, cs$lambda %||% NA, cs$share %||% NA, cs$tau %||% NA,
                                     derivMode = mode)
    ff <- at(NULL); fr <- at("forward-reverse")
    expect_equal(fr$value, ff$value, tolerance = 1e-12)
    expect_equal(fr$gradient, ff$gradient, tolerance = 1e-8)
    expect_equal(fr$modeShift, ff$modeShift, tolerance = 1e-8)
  }
})

# -2 log of the exact marginal of the linear model: per group the normal intercept in closed form,
# the sparse slope by integrate() (Laplace) or in closed form (spike and slab).
.linearExact <- function(theta, lambda = NA, prob = NA, tau = NA) {
  m <- .lap()
  ids <- unique(m$tab$ID)
  om <- exp(theta[["omega_a_a"]])
  sum(vapply(ids, function(id) {
    cn <- rownames(m$tab)[m$tab$ID == id]
    y <- unlist(lapply(m$linear[cn], `[[`, "value"))
    dose <- rep(m$tab[cn, "dose"], each = 4)
    C <- diag(0.3^2, length(y)) + om^2
    logN <- function(e, V) {
      r <- y - theta[["a"]] - (theta[["b"]] + e) * dose
      -0.5 * (sum(r * solve(V, r)) + as.numeric(determinant(2 * pi * V)$modulus))
    }
    if (!is.na(prob))
      return(-2 * log((1 - prob) * exp(logN(0, C)) + prob * exp(logN(0, C + tau^2 * tcrossprod(dose)))))
    f <- function(e) vapply(e, function(x) exp(logN(x, C)) * lambda / 2 * exp(-lambda * abs(x)), 0)
    -2 * log(integrate(f, -Inf, 0, rel.tol = 1e-12)$value + integrate(f, 0, Inf, rel.tol = 1e-12)$value)
  }, 0))
}

test_that("with one sparse coordinate per group the marginal of a linear model is exact", {
  m <- .lap()
  theta <- c(a = 1.1, b = 0.4, omega_a_a = log(0.4))
  objOf <- function(p) normL2(m$linear, m$gl * Xt() * p)
  lap <- .marginalAt(m$pll, theta, 3, obj = objOf(m$pll))
  expect_equal(lap$value, .linearExact(theta, lambda = 3), tolerance = 1e-8)
  slab <- .marginalAt(m$pls, theta, share = 0.4, tau = 0.5, obj = objOf(m$pls))
  expect_equal(slab$value, .linearExact(theta, prob = 0.4, tau = 0.5), tolerance = 1e-8)
  f <- function(v) .marginalAt(m$pll, setNames(v, names(theta)), 3, grad = FALSE, obj = objOf(m$pll))$value
  expect_equal(unname(lap$gradient[names(theta)]), .slope(f, unname(theta)), tolerance = 1e-6)
})

# The exact posterior of the linear model with a normal intercept and a Dirichlet-process slope:
# per partition of the groups and the anchor its log prior plus log N(y; mu, Sigma), y normal
# given the partition (sigma^2 I, the intercept per group, tau^2 per cluster), and per group the
# probability of the population value.
.linearPartitions <- function(theta, conc, tau) {
  m <- .lap()
  ids <- unique(m$tab$ID); G <- length(ids)
  cn <- unlist(lapply(ids, function(id) rownames(m$tab)[m$tab$ID == id]))
  y <- unlist(lapply(m$linear[cn], `[[`, "value"))
  owner <- rep(match(m$tab[cn, "ID"], ids), each = 4)
  dose <- rep(m$tab[cn, "dose"], each = 4)
  base <- diag(0.3^2, length(y)) + exp(2 * theta[["omega_a_a"]]) * outer(owner, owner, `==`)
  mu <- theta[["a"]] + theta[["b"]] * dose
  parts <- list(1L)
  for (k in seq_len(G)) parts <- unlist(lapply(parts, function(p) lapply(seq_len(max(p) + 1L), function(j) c(p, j))),
                                       recursive = FALSE)
  terms <- vapply(parts, function(lab) {
    blocks <- split(seq_len(G + 1L), lab)
    V <- base
    for (b in blocks) if (!((G + 1L) %in% b)) { z <- dose * (owner %in% b); V <- V + tau^2 * tcrossprod(z) }
    prior <- (length(blocks) - 1) * log(conc) + sum(lgamma(lengths(blocks))) - sum(log(conc + seq_len(G)))
    prior - 0.5 * (sum((y - mu) * solve(V, y - mu)) + as.numeric(determinant(2 * pi * V)$modulus))
  }, 0)
  w <- exp(terms - max(terms)); w <- w / sum(w)
  zero <- vapply(parts, function(lab) lab[seq_len(G)] == lab[G + 1L], logical(G))
  list(logI = max(terms) + log(sum(exp(terms - max(terms)))), p0 = setNames(drop(zero %*% w), ids))
}

test_that("a Dirichlet-process slope of a linear model: the marginal sums all partitions exactly", {
  m <- .lap()
  theta <- c(a = 1.1, b = 0.4, omega_a_a = log(0.4))
  r <- .marginalAt(m$pld, theta, share = 0.8, tau = 0.5, obj = normL2(m$linear, m$gl * Xt() * m$pld))
  expect_equal(r$value, -2 * .linearPartitions(theta, 0.8, 0.5)$logI, tolerance = 1e-8)
})

test_that("SAEM with Dirichlet-process random effects: the sampler draws the exact posterior", {
  m <- .lap()
  set.seed(5)
  theta <- c(a = 1.1, b = 0.4, omega_a_a = log(0.4))
  fixed <- c(theta, concentration_ID = 0.8, tau_ID = 0.5)
  obj <- normL2(m$linear, m$gl * Xt() * m$pld)
  fit <- em(obj, randomEffects(m$pld), fixed = fixed, method = "saem",
            control = list(saem = list(nEM = 1500L, particles = 2000L)), verbose = FALSE)
  exact <- .linearPartitions(theta, 0.8, 0.5)
  p0 <- setNames(fit$random$p0[fit$random$eta == "b"], fit$random$group[fit$random$eta == "b"])
  expect_lt(max(abs(p0[names(exact$p0)] - exact$p0)), 0.08)
  # the importance-sampling value of the marginal within a few standard errors
  expect_lt(abs(fit$value + 2 * exact$logI), 4 * fit$valueSE + 1e-6)
})

test_that("two coupled sparse coordinates per group: mean field bounds -2 log L, expectation consistent is closer", {
  m <- .lap()
  theta <- c(a = 1.1, b = 0.4)
  lam <- 2
  obj <- normL2(m$linear, m$gl * Xt() * m$pl2)
  r <- .marginalAt(m$pl2, theta, lam, obj = obj, coupling = "meanField")
  ec <- .marginalAt(m$pl2, theta, lam, obj = obj)
  expect_identical(ec$coupling, "expectationConsistent")
  exact <- sum(vapply(unique(m$tab$ID), function(id) {
    cn <- rownames(m$tab)[m$tab$ID == id]
    y <- unlist(lapply(m$linear[cn], `[[`, "value"))
    dose <- rep(m$tab[cn, "dose"], each = 4)
    dens <- function(x, z) exp(sum(dnorm(y, theta[["a"]] + x + (theta[["b"]] + z) * dose, 0.3, log = TRUE))) *
      lam^2 / 4 * exp(-lam * (abs(x) + abs(z)))
    half <- function(f) integrate(f, -Inf, 0, rel.tol = 1e-10)$value + integrate(f, 0, Inf, rel.tol = 1e-10)$value
    inner <- function(x) vapply(x, function(xi) half(function(z) vapply(z, function(zi) dens(xi, zi), 0)), 0)
    -2 * log(half(inner))
  }, 0))
  expect_gte(r$value, exact - 1e-6)
  expect_lt(abs(ec$value - exact), abs(r$value - exact))
})

test_that("expectation consistent for two coupled Dirichlet-process stems against the sum over pairs of partitions", {
  G <- 3; conc <- 1.2; tau <- 0.6
  parts <- list(1L)
  for (k in seq_len(G)) parts <- unlist(lapply(parts, function(p) lapply(seq_len(max(p) + 1L), function(j) c(p, j))),
                                       recursive = FALSE)
  info <- lapply(parts, function(lab) {
    blocks <- split(seq_len(G + 1L), lab)
    free <- Filter(function(b) !((G + 1L) %in% b), blocks)
    Z <- matrix(0, G, length(free))
    for (j in seq_along(free)) Z[free[[j]], j] <- 1
    list(prior = (length(blocks) - 1) * log(conc) + sum(lgamma(lengths(blocks))) - sum(log(conc + seq_len(G))), Z = Z)
  })
  # log of the integral over the cluster values of a pair of partitions, the data in stem-major order
  pairLog <- function(Z, PF, mv) {
    L <- crossprod(Z, PF %*% Z) + diag(1 / tau^2, ncol(Z)); b <- crossprod(Z, PF %*% mv)
    -ncol(Z) * log(tau) - 0.5 * as.numeric(determinant(L)$modulus) - 0.5 * sum(mv * (PF %*% mv)) +
      0.5 * sum(b * solve(L, b))
  }
  exact <- function(P, mm) {
    PF <- matrix(0, 2 * G, 2 * G)
    for (g in seq_len(G)) PF[c(g, G + g), c(g, G + g)] <- P[[g]]
    terms <- unlist(lapply(info, function(a) vapply(info, function(b) {
      Z <- rbind(cbind(a$Z, matrix(0, G, ncol(b$Z))), cbind(matrix(0, G, ncol(a$Z)), b$Z))
      if (!ncol(Z)) return(a$prior + b$prior - 0.5 * sum(c(mm) * (PF %*% c(mm))))
      a$prior + b$prior + pairLog(Z, PF, c(mm))
    }, 0)))
    max(terms) + log(sum(exp(terms - max(terms))))
  }
  lay <- data.frame(stem = c("s1", "s2"), block = 1L, dist = "dirichletProcess", lambdaName = NA_character_,
                    probName = NA_character_, concentrationName = "concentration_ID", tauName = "tau_ID",
                    lambda = NA_real_, prob = NA_real_, concentration = NA_real_, tau = NA_real_,
                    stringsAsFactors = FALSE)
  hyp <- list(lambda = c(s1 = NA, s2 = NA), prob = c(s1 = NA, s2 = NA), concentration = c(s1 = conc, s2 = conc),
              tau = c(s1 = tau, s2 = tau))
  draw <- function(rho) {
    S <- diag(c(0.3, 0.36)) %*% matrix(c(1, rho, rho, 1), 2) %*% diag(c(0.3, 0.36))
    P <- lapply(seq_len(G), function(g) structure(solve(S), dimnames = list(c("s1", "s2"), c("s1", "s2"))))
    mm <- matrix(c(0, 0.8, 0.8, 0, 0.6, 0), G, 2, dimnames = list(paste0("g", seq_len(G)), c("s1", "s2"))) +
      t(replicate(G, drop(t(chol(S)) %*% rnorm(2))))
    list(P = P, m = mm)
  }
  set.seed(4)
  # without coupling both are exact
  d0 <- draw(0)
  expect_equal(dMod2:::.sparseIntegral(d0$P, d0$m, d0$m, lay, hyp, character(0))$logI, exact(d0$P, d0$m),
               tolerance = 1e-10)
  err <- t(replicate(10, {
    d <- draw(0.6)
    ex <- exact(d$P, d$m)
    ec <- dMod2:::.sparseIntegral(d$P, d$m, d$m, lay, hyp, character(0))
    mf <- dMod2:::.sparseIntegral(d$P, d$m, d$m, lay, hyp, character(0), coupling = "meanField")
    c(coupling = ec$coupling == "expectationConsistent", ec = abs(ec$logI - ex), mf = abs(mf$logI - ex))
  }))
  expect_true(all(err[, "coupling"] == 1))
  expect_lt(mean(err[, "ec"]), mean(err[, "mf"]) / 3)
})

# the gradient of a fit on the unbounded scale of its hyperparameters, log or logit
.freeGradient <- function(fit) {
  g <- fit$gradient[names(fit$argument)]
  v <- fit$argument
  pos <- grepl("^(lambda|tau|concentration)_", names(g))
  prob <- grepl("^prob_", names(g))
  g[pos] <- g[pos] * v[pos]
  g[prob] <- g[prob] * v[prob] * (1 - v[prob])
  g
}

test_that("a centre fit keeps the shared parameter shared and estimates its strength", {
  m <- .lap()
  obj <- normL2(m$data, m$g * Xt() * m$pc, errmodel = m$err)
  fit <- em(obj, randomEffects(m$pc), m$init, verbose = FALSE)
  expect_true(fit$converged)
  expect_lt(max(abs(.freeGradient(fit))), 0.02)
  expect_true(is.finite(fit$strengths$estimate) && fit$strengths$estimate > 0)
  expect_lt(max(abs(fit$etaModes[, "eta_A"])), 0.05)
  expect_gt(abs(fit$etaModes["L3", "eta_k"]), 0.3)
})

test_that("a prior beside the data enters a fit of sparse random effects", {
  m <- .lap()
  obj <- normL2(m$data, m$g * Xt() * m$pc, errmodel = m$err)
  plain <- em(obj, randomEffects(m$pc), m$init, verbose = FALSE)
  target <- plain$argument[["A"]] + 0.3
  pr <- constraintL2(mu = c(A = target), sigma = 1e-3)
  fit <- em(obj + pr, randomEffects(m$pc), m$init, verbose = FALSE)
  expect_equal(fit$argument[["A"]], target, tolerance = 1e-3)
  expect_error(em(obj + constraintL2(mu = c(eta_k_L1 = 0), sigma = 1), randomEffects(m$pc), m$init),
               "fixed effects only")
})

test_that("normal and Laplace random effects in one fit: exact gradient and the two limits", {
  m <- .lap()
  theta <- c(m$init + c(0.05, -0.02, 0.1), omega_A_A = log(0.2))
  r <- .marginalAt(m$pm, theta, 3)
  f <- function(v) .marginalAt(m$pm, setNames(v, names(theta)), 3, grad = FALSE)$value
  expect_equal(unname(r$gradient[names(theta)]), .slope(f, unname(theta)), tolerance = 1e-5)
  objOf <- function(p) normL2(m$data, m$g * Xt() * p, errmodel = m$err)
  # a huge strength pins the Laplace etas at 0: the fit of the normal ones alone
  big <- em(objOf(m$pm), randomEffects(m$pm), m$init, fixed = c(lambda_ID = 1e5), verbose = FALSE)
  normal <- em(objOf(m$pn), randomEffects(m$pn), m$init, verbose = FALSE)
  expect_equal(big$value, normal$value, tolerance = 1e-3)
  expect_equal(unname(big$argument[names(normal$argument)]), unname(normal$argument), tolerance = 1e-3)
  # a tiny Omega pins the normal etas at 0: the fit of the Laplace ones alone
  tiny <- em(objOf(m$pm), randomEffects(m$pm), m$init, fixed = c(omega_A_A = log(1e-5)), verbose = FALSE)
  lap <- em(objOf(m$pl), randomEffects(m$pl), m$init, verbose = FALSE)
  expect_equal(tiny$value, lap$value, tolerance = 1e-3)
  # sparse fits have no Hessian of the marginal, their standard errors come from the empirical
  # information of the group scores; CWRES without normal random effects is IWRES
  expect_null(lap$hessian)
  expect_equal(summary(lap)$seMethod, "scores")
  expect_equal(vcov(lap), solve(lap$information), ignore_attr = TRUE)
  expect_true(all(is.finite(summary(lap)$population$se)))
  expect_equal(residuals(lap, "CWRES"), residuals(lap, "IWRES"))
  expect_true(big$strengths$fixed)
  d <- do.call(rbind, lapply(names(m$data), function(cn) transform(m$data[[cn]], condition = cn)))
  d$lloq <- stats::quantile(d$value, 0.2)
  m1 <- normL2(as.datalist(d), m$g * Xt() * m$pl, errmodel = m$err, optBLOQ = "M1")
  expect_error(em(m1, randomEffects(m$pl), m$init, verbose = FALSE), "M3")
})

test_that("the group scores of a sparse marginal sum to its gradient, each the gradient of its group alone", {
  m <- .lap()
  objOf <- function(p, conds = names(m$data)) normL2(m$data[conds], m$g * Xt() * p, errmodel = m$err)
  theta <- m$init + c(0.05, -0.02, 0.1)
  normal <- c(theta, omega_A_A = log(0.2))
  for (cs in list(list(m$pc, theta), list(m$pm, normal), list(m$ps, theta), list(m$psm, normal))) {
    at <- function(obj = NULL, scores = FALSE)
      .marginalAt(cs[[1]], cs[[2]], 3, share = 0.4, tau = 0.5, obj = obj, scores = scores)
    r <- at(scores = TRUE)
    free <- replace(0 * r$gradient, names(r$freeGradient), r$freeGradient)
    expect_equal(colSums(r$scores)[names(r$gradient)] + free, r$gradient, tolerance = 1e-10)
    for (g in c("L1", "L3")) {
      one <- at(objOf(cs[[1]], rownames(m$tab)[m$tab$ID == g]))
      expect_equal(r$scores[g, ], one$gradient[colnames(r$scores)], tolerance = 1e-6)
    }
  }
})

test_that("where expectation consistent fails at the start the fit is by mean field and says so", {
  m <- .lap()
  obj <- normL2(m$data, m$g * Xt() * m$pc, errmodel = m$err)
  mf <- em(obj, randomEffects(m$pc), m$init, verbose = FALSE, control = list(laplace = list(coupling = "meanField")))
  expect_equal(mf$coupling, "meanField")
  expect_error(em(obj, randomEffects(m$pc), m$init, control = list(laplace = list(coupling = "ec"))), "coupling")
  expect_error(em(obj, randomEffects(m$pc), m$init, control = list(laplace = list(ecmSteps = 5))), "unknown")
  local_mocked_bindings(.sparseExpectation = function(...) NULL, .package = "dMod2")
  fb <- em(obj, randomEffects(m$pc), m$init, verbose = FALSE)
  expect_equal(fb$coupling, "meanFieldFallback")
  expect_true(fb$converged)
  expect_equal(fb$value, mf$value, tolerance = 1e-5)
})

test_that("a trial point at which expectation consistent fails is a failed evaluation, the coupling stays", {
  m <- .lap()
  obj <- normL2(m$data, m$g * Xt() * m$pc, errmodel = m$err)
  ec <- em(obj, randomEffects(m$pc), m$init, verbose = FALSE)
  expect_equal(ec$coupling, "expectationConsistent")
  real <- dMod2:::.sparseExpectation
  calls <- 0L
  local_mocked_bindings(.sparseExpectation = function(...) {
    calls <<- calls + 1L
    if (calls %in% 3:4) NULL else real(...)
  }, .package = "dMod2")
  fit <- em(obj, randomEffects(m$pc), m$init, verbose = FALSE)
  expect_gt(calls, 4L)
  expect_equal(fit$coupling, "expectationConsistent")
  expect_true(fit$converged)
  expect_equal(fit$value, ec$value, tolerance = 1e-5)
})


# Every set partition of the groups and the anchor (element G + 1): its Ewens prior with
# concentration `conc`, the block integrals by integrate(), the groups at the anchor and the labels.
.partitionsBrute <- function(a, m, conc, tau) {
  G <- length(a); n <- G + 1L
  parts <- list(1L)
  for (k in seq_len(n - 1L)) parts <- unlist(lapply(parts, function(p)
    lapply(seq_len(max(p) + 1L), function(j) c(p, j))), recursive = FALSE)
  lapply(parts, function(lab) {
    blocks <- split(seq_len(n), lab)
    logw <- (length(blocks) - 1) * log(conc) + sum(lgamma(lengths(blocks))) - sum(log(conc + seq_len(G)))
    for (b in blocks) {
      g <- b[b <= G]
      if (n %in% b) {
        logw <- logw - sum(a[g] * m[g]^2) / 2
        next
      }
      f <- function(v) vapply(v, function(x) exp(-sum(a[g] * (x - m[g])^2) / 2), 0) * dnorm(v, 0, tau)
      logw <- logw + log(integrate(f, -Inf, Inf, rel.tol = 1e-13)$value)
    }
    anchorBlock <- lab[n]
    free <- setdiff(unique(lab[seq_len(G)]), anchorBlock)
    labels <- ifelse(lab[seq_len(G)] == anchorBlock, 0L, match(lab[seq_len(G)], free))
    list(w = exp(logw), zero = lab[seq_len(G)] == anchorBlock, labels = labels)
  })
}

test_that("the Dirichlet-process integral is the Ewens sum over all set partitions, with exact derivatives", {
  a <- c(3, 1.5, 4, 2); m <- c(0.9, -0.2, 0.75, 0.05); conc <- 0.7; tau <- 0.8
  ps <- dMod2:::.slabPartition(a, m, conc, tau)
  all <- .partitionsBrute(a, m, conc, tau)
  w <- vapply(all, `[[`, 0, "w")
  expect_equal(ps$logI, log(sum(w)), tolerance = 1e-9)
  expect_equal(ps$spike, colSums(w * t(vapply(all, `[[`, logical(4), "zero"))) / sum(w), tolerance = 1e-8)
  # without data every block integral is 1 and the prior sums to 1
  expect_equal(dMod2:::.slabPartition(0 * a, m, conc, tau)$logI, 0, tolerance = 1e-12)
  # the prior does not depend on the order of the groups
  o <- c(3, 1, 4, 2)
  expect_equal(dMod2:::.slabPartition(a[o], m[o], conc, tau)$logI, ps$logI, tolerance = 1e-12)
  logI <- function(a, m, conc, tau) dMod2:::.slabPartition(a, m, conc, tau)$logI
  expect_equal(ps$dm, numDeriv::grad(function(x) logI(a, x, conc, tau), m), tolerance = 1e-7)
  expect_equal(ps$da, numDeriv::grad(function(x) logI(x, m, conc, tau), a), tolerance = 1e-7)
  expect_equal(ps$dconcentration, numDeriv::grad(function(x) logI(a, m, x, tau), conc), tolerance = 1e-7)
  expect_equal(ps$dtau, numDeriv::grad(function(x) logI(a, m, conc, x), tau), tolerance = 1e-7)
  # the most probable partition
  expect_equal(ps$labels, all[[which.max(w)]]$labels)
  # the M-step of the concentration: the expected number of slab blocks
  th <- dMod2:::.ewensConcentration(1.3, c(4, 6))
  expect_equal(sum(th / (th + c(1:4, 1:6))), 1.3, tolerance = 1e-10)
})

test_that("the spike-and-slab integral of a group sums the subsets on the full curvature", {
  P <- matrix(c(4, 1.5, 1.5, 2), 2); m <- c(0.6, -0.3); prob <- c(0.3, 0.6); tau <- c(0.8, 0.5)
  ss <- dMod2:::.slabSubsets(P, m, prob, tau)
  quad <- function(e) exp(-sum((e - m) * (P %*% (e - m))) / 2)
  line <- function(k, v) { e <- c(0, 0); e[k] <- v; quad(e) }
  one <- function(k) integrate(function(v) vapply(v, function(x) line(k, x), 0) * dnorm(v, 0, tau[k]),
                               -Inf, Inf, rel.tol = 1e-12)$value
  both <- integrate(function(x) vapply(x, function(xi) integrate(function(y) vapply(y, function(yi)
    quad(c(xi, yi)), 0) * dnorm(y, 0, tau[2]), -Inf, Inf, rel.tol = 1e-11)$value, 0) * dnorm(x, 0, tau[1]),
    -Inf, Inf, rel.tol = 1e-10)$value
  I <- (1 - prob[1]) * (1 - prob[2]) * quad(c(0, 0)) + prob[1] * (1 - prob[2]) * one(1) +
    (1 - prob[1]) * prob[2] * one(2) + prob[1] * prob[2] * both
  expect_equal(ss$logI, log(I), tolerance = 1e-8)
  logI <- function(P, m, prob, tau) dMod2:::.slabSubsets(P, m, prob, tau)$logI
  expect_equal(ss$dm, numDeriv::grad(function(x) logI(P, x, prob, tau), m), tolerance = 1e-7)
  # P is symmetric: the off-diagonal entries move together
  sym <- function(x) matrix(c(x[1], x[2], x[2], x[3]), 2)
  dP <- numDeriv::grad(function(x) logI(sym(x), m, prob, tau), c(P[1, 1], P[1, 2], P[2, 2]))
  expect_equal(c(ss$dP[1, 1], 2 * ss$dP[1, 2], ss$dP[2, 2]), dP, tolerance = 1e-7)
  expect_equal(ss$dprob, numDeriv::grad(function(x) logI(P, m, x, tau), prob), tolerance = 1e-7)
  expect_equal(ss$dtau, numDeriv::grad(function(x) logI(P, m, prob, x), tau), tolerance = 1e-7)
  # one coordinate is the partition of one group with the anchor, c / (1 + c) its slab share
  one <- dMod2:::.slabSubsets(matrix(3), 0.4, 0.3 / 1.3, 0.7)
  expect_equal(one$logI, dMod2:::.slabPartition(3, 0.4, 0.3, 0.7)$logI)
})

test_that("more spike-and-slab coordinates than integrated jointly: one factor each, expectation consistent", {
  K <- 13L; G <- 2L; st <- paste0("s", seq_len(K))
  # the subset sum of .slabSubsets without its bound on K
  exact <- function(P, m, prob, tau) {
    b <- drop(P %*% m); c0 <- sum(m * b)
    terms <- apply(as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), K))), 1, function(F) {
      out <- sum(log(prob[F])) + sum(log1p(-prob[!F])) - c0 / 2
      if (!any(F)) return(out)
      R <- chol(P[F, F, drop = FALSE] + diag(1 / tau[F]^2, sum(F)))
      out + sum(backsolve(R, b[F], transpose = TRUE)^2) / 2 - sum(log(diag(R))) - sum(log(tau[F]))
    })
    max(terms) + log(sum(exp(terms - max(terms))))
  }
  lay <- data.frame(stem = st, block = 1L, dist = "spikeSlab", lambdaName = NA_character_, probName = "prob_ID",
                    concentrationName = NA_character_, tauName = "tau_ID", lambda = NA_real_, prob = NA_real_,
                    concentration = NA_real_, tau = NA_real_, stringsAsFactors = FALSE)
  hyp <- list(lambda = setNames(rep(NA, K), st), prob = setNames(rep(0.3, K), st),
              concentration = setNames(rep(NA, K), st), tau = setNames(rep(0.8, K), st))
  draw <- function(rho) {
    P <- lapply(seq_len(G), function(g) {
      A <- diag(runif(K, 2, 6))
      A[cbind(1:(K - 1), 2:K)] <- A[cbind(2:K, 1:(K - 1))] <- rho * sqrt(diag(A)[-K] * diag(A)[-1])
      structure(A, dimnames = list(st, st))
    })
    m <- matrix(rnorm(G * K, 0, 0.6) * rbinom(G * K, 1, 0.4), G, K, dimnames = list(paste0("g", seq_len(G)), st))
    list(P = P, m = m)
  }
  logI <- function(d) sum(vapply(seq_len(G), function(g) exact(d$P[[g]], d$m[g, ], hyp$prob, hyp$tau), 0))
  set.seed(1)
  d0 <- draw(0)
  expect_equal(lengths(dMod2:::.sparseFactors(d0$P, st, lay)$factors), rep(1L, K))
  # without coupling the factors are exact
  expect_equal(dMod2:::.sparseIntegral(d0$P, d0$m, d0$m, lay, hyp, character(0))$logI, logI(d0), tolerance = 1e-10)
  d <- draw(0.3)
  ec <- dMod2:::.sparseIntegral(d$P, d$m, d$m, lay, hyp, character(0))
  mf <- dMod2:::.sparseIntegral(d$P, d$m, d$m, lay, hyp, character(0), coupling = "meanField")
  expect_identical(ec$coupling, "expectationConsistent")
  expect_lt(abs(ec$logI - logI(d)), abs(mf$logI - logI(d)) / 3)
})

test_that("the spike-and-slab and cluster marginals have the exact gradient, normal etas beside", {
  m <- .lap()
  # the inner modes are solved to 1e-10: central differences with a step above that noise and
  # Richardson extrapolation, error of fourth order in the step
  cdv <- function(f, x, h = 2e-4) vapply(seq_along(x), function(j) {
    d <- function(s) { e <- replace(numeric(length(x)), j, s); (f(x + e) - f(x - e)) / (2 * s) }
    (4 * d(h / 2) - d(h)) / 3
  }, 0)
  for (p in list(m$ps, m$psk, m$psm, m$pskm, m$psk2)) {
    theta <- m$init + c(0.05, -0.02, 0.1)
    if (any(randomEffects(p)$terms$dist == "normal")) theta <- c(theta, omega_A_A = log(0.2))
    r <- .marginalAt(p, theta, share = 0.4, tau = 0.5)
    f <- function(v) .marginalAt(p, setNames(v, names(theta)), share = 0.4, tau = 0.5, grad = FALSE)$value
    ft <- function(x) .marginalAt(p, theta, share = x, tau = 0.5, grad = FALSE)$value
    fs <- function(x) .marginalAt(p, theta, share = 0.4, tau = x, grad = FALSE)$value
    expect_equal(unname(r$gradient[names(theta)]), cdv(f, unname(theta)), tolerance = 1e-5)
    # the hyperparameters with a larger step: the value is noisy at 1e-8 below steps of 1e-4
    expect_equal(unname(r$gradient[grep("^(prob|concentration)_", names(r$gradient))]), cdv(ft, 0.4, 2e-3),
                 tolerance = 1e-5)
    expect_equal(unname(r$gradient[grep("^tau_", names(r$gradient))]), cdv(fs, 0.5, 2e-3), tolerance = 1e-5)
  }
})

test_that("the curvature of the hyperparameters is that of the complete-data prior, its secant the M-step", {
  lay <- data.frame(stem = c("a", "b", "c"), lambdaName = c("lambda_ID", NA, NA),
                    probName = c(NA, "prob_ID", NA), concentrationName = c(NA, NA, "concentration_ID"),
                    tauName = c(NA, "tau_ID", "tau_ID_2"), stringsAsFactors = FALSE)
  kind <- dMod2:::.hyperKind(lay)
  G <- 5
  r <- list(abs = matrix(c(0.3, 0.1, 0.7, 0.2, 0.4), G, 3, dimnames = list(NULL, lay$stem)),
            slab = matrix(c(0, 0, 0, 2.6, G, 1.3, 1.8, G, 0.9), 3, byrow = TRUE,
                          dimnames = list(lay$stem, c("free", "groups", "slab2"))))
  hyper <- c(lambda_ID = 2.1, prob_ID = 0.35, concentration_ID = 1.4, tau_ID = 0.6, tau_ID_2 = 0.8)
  # the expected -2 log prior of every hyperparameter as a function of its unbounded value
  prior <- list(
    lambda_ID = function(u) sum(-2 * log(exp(u) / 2) + 2 * exp(u) * r$abs[, "a"]),
    prob_ID = function(u) { p <- plogis(u); K <- r$slab["b", "free"]; -2 * (K * log(p) + (G - K) * log(1 - p)) },
    concentration_ID = function(u) -2 * (r$slab["c", "free"] * u - sum(log(exp(u) + seq_len(G)))),
    tau_ID = function(u) 2 * r$slab["b", "free"] * u + r$slab["b", "slab2"] * exp(-2 * u),
    tau_ID_2 = function(u) 2 * r$slab["c", "free"] * u + r$slab["c", "slab2"] * exp(-2 * u))
  u <- dMod2:::.hyperToFree(kind, hyper)
  expect_equal(dMod2:::.hyperFromFree(kind, u), hyper)
  second <- vapply(names(hyper), function(h) numDeriv::hessian(prior[[h]], u[[h]])[1, 1], 0)
  expect_equal(dMod2:::.hyperInformation(lay, kind, r, hyper), second, tolerance = 1e-7)
  # the Newton step on the secant curvature lands at the minimum of the expected -2 log prior,
  # and close to the minimum the secant is the curvature
  B <- dMod2:::.hyperSecant(lay, kind, r, hyper)
  expect_true(all(B > 0))
  slope <- vapply(names(hyper), function(h) numDeriv::grad(prior[[h]], u[[h]]), 0)
  minimum <- vapply(names(hyper), function(h) optimize(prior[[h]], u[[h]] + c(-20, 20), tol = 1e-10)$minimum, 0)
  expect_equal(u - slope / B, minimum, tolerance = 1e-6)
  near <- dMod2:::.hyperFromFree(kind, minimum + 1e-5)
  expect_equal(dMod2:::.hyperSecant(lay, kind, r, near), dMod2:::.hyperInformation(lay, kind, r, near),
               tolerance = 1e-4)
})

test_that("multistart perturbs the hyperparameters on their unbounded scale", {
  m <- .lap()
  set.seed(2)
  for (p in list(m$ps, m$psk, m$pc)) {
    d <- dMod2:::.emDesign(normL2(m$data, m$g * Xt() * p, errmodel = m$err), randomEffects(p))
    center <- dMod2:::.emStart(d, m$init)
    draws <- t(replicate(200, dMod2:::.perturbStart(d, center, 3)))
    hyper <- grep("^(lambda|concentration|tau)_", colnames(draws))
    expect_true(all(draws[, hyper] > 0))
    prob <- grep("^prob_", colnames(draws))
    expect_true(all(draws[, prob] > 0 & draws[, prob] < 1))
    expect_gt(sd(draws[, "k"]), 1)
  }
})

test_that("SAEM with Dirichlet-process random effects selects the groups that differ", {
  m <- .lap()
  set.seed(3)
  obj <- normL2(m$data, m$g * Xt() * m$psk, errmodel = m$err)
  fit <- em(obj, randomEffects(m$psk), m$init, method = "saem", control = list(saem = list(nEM = 300L)),
            verbose = FALSE)
  cl <- setNames(fit$random$cluster[fit$random$eta == "k"], fit$random$group[fit$random$eta == "k"])
  expect_identical(cl[["L1"]], cl[["L2"]])
  expect_identical(cl[["L3"]], cl[["L4"]])
  expect_false(cl[["L3"]] == cl[["L1"]])
  expect_false(cl[["L5"]] %in% cl[c("L1", "L3")])
  expect_true(is.finite(fit$value) && fit$valueSE > 0)
  expect_gt(fit$argument[["concentration_ID"]], 0)
  expect_error(vcov(fit), "Dirichlet-process")
})

test_that("spike-and-slab fits select the groups that differ and estimate the threshold", {
  m <- .lap()
  objOf <- function(p) normL2(m$data, m$g * Xt() * p, errmodel = m$err)
  fit <- em(objOf(m$psk), randomEffects(m$psk), m$init, verbose = FALSE)
  expect_true(fit$converged)
  # started at the estimate, the fit finds no lower value
  again <- em(objOf(m$psk), randomEffects(m$psk), fit$argument, verbose = FALSE)
  expect_gt(again$value, fit$value - 1e-4)
  th <- fit$slab$estimate[fit$slab$name == "concentration_ID"]
  expect_true(th > 0)
  e <- fit$etaModes[, "eta_k"]
  cl <- setNames(fit$random$cluster, fit$random$group)
  expect_identical(e[["L3"]], e[["L4"]])
  expect_identical(cl[["L1"]], cl[["L2"]])
  expect_false(cl[["L3"]] == cl[["L1"]])
  expect_false(cl[["L5"]] %in% cl[c("L1", "L3")])
  # with the structural parameters fixed only the concentration and tau are estimated
  only <- em(objOf(m$psk), randomEffects(m$psk), fixed = fit$argument[c("k", "A", "s")],
             control = list(laplace = list(maxOuter = 2L)), verbose = FALSE)
  expect_setequal(names(only$argument), c("concentration_ID", "tau_ID"))
  expect_true(is.finite(only$value))
  expect_lte(only$iterations, 2L)
  expect_error(vcov(fit), "Dirichlet-process")
  expect_equal(summary(fit)$seMethod, "none")
  centre <- em(objOf(m$ps), randomEffects(m$ps), m$init, verbose = FALSE)
  expect_true(centre$converged)
  again <- em(objOf(m$ps), randomEffects(m$ps), centre$argument, verbose = FALSE)
  expect_gt(again$value, centre$value - 1e-5)
  # five groups give scores in four directions, fewer than the five parameters
  expect_error(vcov(centre), "do not span")
  expect_true(all(is.na(summary(centre)$population$se)))
  A <- centre$random[centre$random$eta == "A", ]
  expect_true(all(A$mode == 0))
  expect_true(all(A$p0 > 0.5))
})
