# Spike-and-slab and Dirichlet-process random effects: the integrals against numerical integration,
# the exact gradient of the approximate marginal, exactness on a linear model, and fits.

# The integral over the sparse coordinates in center form, int exp(-(eta - m)' P (eta - m) / 2) prod p:
# that of the natural form with b = P m, less m' P m / 2 per group.
.centered <- function(P, m, ...) {
  b <- t(vapply(seq_along(P), function(i) drop(P[[i]] %*% m[i, ]), numeric(ncol(m))))
  dimnames(b) <- dimnames(m)
  r <- dMod2:::.sparseIntegral(P, b, ...)
  if (is.null(r)) return(r)
  r$logI <- r$logI - sum(vapply(seq_along(P), function(i) sum(m[i, ] * (P[[i]] %*% m[i, ])) / 2, 0))
  r
}

# A decay model y = A exp(-k dose) per condition, five groups with two conditions each, and the
# linear model y = a + b dose on the same design with a known sigma.
.sparse <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    withr::local_preserve_seed()
    dir <- file.path(tempdir(), "emsparse_models")
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
           modelname = "sp_g", deriv2 = TRUE, derivMode = .dm)
    err <- Y(c(y = "s"), g, attachInput = FALSE, compile = FALSE, modelname = "sp_err", deriv2 = TRUE, derivMode = .dm)
    pn <- P(insert(base(), A ~ 1 + (1 | ID)), method = "explicit", compile = FALSE, modelname = "sp_normal",
            deriv2 = TRUE, derivMode = .dm)
    pks <- P(insert(base(), k ~ 1 + (1 | gr(ID, dist = spikeSlab()))), method = "explicit", compile = FALSE,
             modelname = "sp_konly", deriv2 = TRUE, derivMode = .dm)
    ps <- P(insert(base(), k + A ~ 1 + (1 | gr(ID, dist = spikeSlab()))), method = "explicit",
            compile = FALSE, modelname = "sp_slab", deriv2 = TRUE, derivMode = .dm)
    psk <- P(insert(base(), k ~ 1 + (1 | gr(ID, dist = dirichletProcess()))), method = "explicit",
             compile = FALSE, modelname = "sp_dpcenter", deriv2 = TRUE, derivMode = .dm)
    psm <- P(insert(base(), k ~ 1 + (1 | gr(ID, dist = spikeSlab())), A ~ 1 + (1 | ID)), method = "explicit",
             compile = FALSE, modelname = "sp_slabmixed", deriv2 = TRUE, derivMode = .dm)
    pskm <- P(insert(base(), k ~ 1 + (1 | gr(ID, dist = dirichletProcess())), A ~ 1 + (1 | ID)),
              method = "explicit", compile = FALSE, modelname = "sp_slabclmixed", deriv2 = TRUE, derivMode = .dm)
    psk2 <- P(insert(base(), k + A ~ 1 + (1 | gr(ID, dist = dirichletProcess()))), method = "explicit",
              compile = FALSE, modelname = "sp_slabcltwo", deriv2 = TRUE, derivMode = .dm)
    linear <- function(dist, a = a ~ 1 + (1 | ID)) {
      tr <- branch(define(NULL, "x <- x", x = c("a", "b")), table = tab)
      tr <- define(tr, "dose <- d", d = dose)
      insert(tr, a, b ~ 1 + (1 | gr(ID, dist = dist)))
    }
    gl <- Y(c(y = "a + b * dose"), f = NULL, parameters = c("a", "b", "dose"), compile = FALSE,
            modelname = "sp_glin", deriv2 = TRUE, derivMode = .dm)
    pls <- P(linear(spikeSlab()), method = "explicit", compile = FALSE, modelname = "sp_linslab", deriv2 = TRUE, derivMode = .dm)
    pld <- P(linear(dirichletProcess()), method = "explicit", compile = FALSE, modelname = "sp_lindp", deriv2 = TRUE, derivMode = .dm)
    # a second observable z = b + c that only the first two groups measure
    gz <- Y(c(y = "a + b * dose", z = "b + c"), f = NULL, parameters = c("a", "b", "c", "dose"), compile = FALSE,
            modelname = "sp_gz", deriv2 = TRUE, derivMode = .dm)
    reachP <- function(dist, name) {
      tr <- branch(define(NULL, "x <- x", x = c("a", "b", "c")), table = tab)
      tr <- define(tr, "dose <- d", d = dose)
      P(insert(tr, a ~ 1 + (1 | ID), b + c ~ 1 + (1 | gr(ID, dist = dist))), method = "explicit", compile = FALSE,
        modelname = name, deriv2 = TRUE, derivMode = .dm)
    }
    pz <- list(spikeSlab = reachP(spikeSlab(), "sp_zslab"), dirichletProcess = reachP(dirichletProcess(), "sp_zdp"))
    # a state X that a stimulus u drives from 0, u = 0 in the groups L1 and L2, read with a drift w
    # that only a state at 0 multiplies
    tabu <- cbind(tab, u = ifelse(tab$ID %in% c("L1", "L2"), 0, 1))
    xu <- Xs(odemodel(eqnvec(X = "s * u - dg * X", W = "w * X - dg * W"), modelname = "sp_xu", compile = FALSE,
                      deriv2 = TRUE, derivMode = .dm))
    gu <- Y(c(y = "X + W + c0"), xu, compile = FALSE, modelname = "sp_gu", deriv2 = TRUE, derivMode = .dm)
    tru <- branch(define(NULL, "x <- x", x = c("s", "w", "dg", "c0", "u", "X", "W")), table = tabu)
    tru <- define(tru, "u <- v", v = u)
    tru <- insert(insert(tru, "X <- 0"), "W <- 0")
    pu <- P(insert(tru, c0 ~ 1 + (1 | ID), s + w ~ 1 + (1 | gr(ID, dist = spikeSlab()))), method = "explicit",
            compile = FALSE, modelname = "sp_pu", deriv2 = TRUE, derivMode = .dm)
    compile(g, err, pn, pks, ps, psk, psm, pskm, psk2, gl, pls, pld, gz, pz$spikeSlab, pz$dirichletProcess, xu, gu,
            pu, output = "emsparse_models", cores = 2)
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
    linear <- do.call(rbind, rows)
    ec <- c(L1 = 0.5, L2 = -0.4)
    zrows <- lapply(rownames(tab)[tab$ID %in% names(ec)], function(cn) {
      id <- tab[cn, "ID"]
      data.frame(name = "z", time = c(0, 1, 2, 3), value = 0.5 + eb[[id]] + 0.2 + ec[[id]] + rnorm(4, 0, 0.3),
                 sigma = 0.3, condition = cn)
    })
    cache <<- list(g = g, err = err, pn = pn, pks = pks, ps = ps, psk = psk, psm = psm,
                   pskm = pskm, psk2 = psk2, data = data, init = c(k = log(0.6), A = log(10), s = log(0.1)),
                   gl = gl, pls = pls, pld = pld, tab = tab, linear = as.datalist(linear),
                   gz = gz, pz = pz, zdata = as.datalist(rbind(linear, do.call(rbind, zrows))),
                   xu = xu, gu = gu, pu = pu, udata = as.datalist(transform(linear, value = value - 1)))
    cache
  }
})

# The marginal of `p` at `theta` with every slab probability or concentration `share` and slab
# standard deviation `tau`; the decay data unless `obj` is given.
.marginalAt <- function(p, theta, share = NA, tau = NA, grad = TRUE, obj = NULL,
                        scores = FALSE, derivMode = NULL, reach = TRUE, tolControl = NULL,
                        warm = NULL, iterlim = 300L) {
  m <- .sparse()
  obj <- obj %||% normL2(m$data, m$g * Xt() * p, errmodel = m$err)
  d <- dMod2:::.emDesign(obj, randomEffects(p))
  if (!reach) d$reach <- NULL
  lay <- dMod2:::.sparseLayout(d)
  hyp <- lapply(list(prob = share, concentration = share, tau = tau),
                function(v) setNames(rep(v, nrow(lay)), lay$stem))
  prior <- if (d$normal$K) dMod2:::.normalPrior(d)
  warm <- warm %||% matrix(0, length(d$groups), ncol(d$etaNames))
  dMod2:::.sparseMarginal(d, dMod2:::.conditionObjectives(d), prior, theta, hyp, names(theta),
                          dMod2:::.sparseHypers(lay), NULL, warm, grad = grad,
                          control = list(iterlim = iterlim, derivMode = derivMode,
                                         trustControl = list(tolControl = tolControl)), scores = scores)
}

# The slope of a marginal: the inner iterations leave noise of about 1e-10 in the value, so the
# difference steps stay above 1e-3.
.slope <- function(f, x) numDeriv::grad(f, x, method.args = list(eps = 1e-2, d = 1e-2, r = 4))

test_that("the marginal alone needs no second-order sensitivities under forward over reverse", {
  m <- .sparse()
  theta <- m$init + c(0.05, -0.02, 0.1)
  full <- .marginalAt(m$ps, theta, 0.4, 0.5, derivMode = "forward-reverse")
  real <- dMod2:::.predictMany
  local_mocked_bindings(.predictMany = function(...) {
    if (isTRUE(list(...)$deriv2)) stop("second-order sensitivities requested")
    real(...)
  }, .package = "dMod2")
  alone <- .marginalAt(m$ps, theta, 0.4, 0.5, grad = FALSE, derivMode = "forward-reverse")
  # the value alone lacks the first-order correction for the modes, which needs the exact inner Hessian
  expect_equal(alone$value, full$value, tolerance = 1e-8)
})

test_that("forward over reverse gives the gradient and the model Hessian rows of forward over forward", {
  m <- .sparse()
  th <- m$init + c(0.05, -0.02, 0.1)
  cases <- list(list(p = m$ps, theta = th, share = 0.4),
                list(p = m$psm, theta = c(th, omega_A_A = log(0.2)), share = 0.4),
                list(p = m$psk, theta = th, share = 2))
  for (cs in cases) {
    at <- function(mode) .marginalAt(cs$p, cs$theta, cs$share, 0.5, derivMode = mode)
    ff <- at(NULL); fr <- at("forward-reverse")
    expect_equal(fr$value, ff$value, tolerance = 1e-12)
    expect_equal(fr$gradient, ff$gradient, tolerance = 1e-8)
    expect_equal(fr$modeShift, ff$modeShift, tolerance = 1e-8)
  }
})

# -2 log of the exact marginal of the linear model: per group the normal intercept and the
# spike-and-slab slope in closed form.
.linearExact <- function(theta, prob, tau) {
  m <- .sparse()
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
    -2 * log((1 - prob) * exp(logN(0, C)) + prob * exp(logN(0, C + tau^2 * tcrossprod(dose))))
  }, 0))
}

test_that("with one sparse coordinate per group the marginal of a linear model is exact", {
  m <- .sparse()
  theta <- c(a = 1.1, b = 0.4, omega_a_a = log(0.4))
  objOf <- function(p) normL2(m$linear, m$gl * Xt() * p)
  slab <- .marginalAt(m$pls, theta, share = 0.4, tau = 0.5, obj = objOf(m$pls))
  expect_equal(slab$value, .linearExact(theta, prob = 0.4, tau = 0.5), tolerance = 1e-8)
  f <- function(v) .marginalAt(m$pls, setNames(v, names(theta)), 0.4, 0.5, grad = FALSE, obj = objOf(m$pls))$value
  expect_equal(unname(slab$gradient[names(theta)]), .slope(f, unname(theta)), tolerance = 1e-6)
})

# The exact posterior of the linear model with a normal intercept and a Dirichlet-process slope:
# per partition of the groups and the anchor the log prior of the partition plus log N(y; mu, Sigma), y normal
# given the partition (sigma^2 I, the intercept per group, tau^2 per cluster), and per group the
# probability of the population value.
.linearPartitions <- function(theta, conc, tau) {
  m <- .sparse()
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
  m <- .sparse()
  theta <- c(a = 1.1, b = 0.4, omega_a_a = log(0.4))
  r <- .marginalAt(m$pld, theta, share = 0.8, tau = 0.5, obj = normL2(m$linear, m$gl * Xt() * m$pld))
  expect_equal(r$value, -2 * .linearPartitions(theta, 0.8, 0.5)$logI, tolerance = 1e-8)
})

test_that("SAEM with Dirichlet-process random effects: the sampler draws the exact posterior", {
  m <- .sparse()
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
  lay <- data.frame(stem = c("s1", "s2"), block = 1L, dist = "dirichletProcess", probName = NA_character_,
                    concentrationName = "concentration_ID", tauName = "tau_ID", prob = NA_real_,
                    concentration = NA_real_, tau = NA_real_, stringsAsFactors = FALSE)
  hyp <- list(prob = c(s1 = NA, s2 = NA), concentration = c(s1 = conc, s2 = conc), tau = c(s1 = tau, s2 = tau))
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
  expect_equal(.centered(d0$P, d0$m, lay, hyp, character(0))$logI, exact(d0$P, d0$m), tolerance = 1e-10)
  err <- t(replicate(10, {
    d <- draw(0.6)
    ex <- exact(d$P, d$m)
    ec <- .centered(d$P, d$m, lay, hyp, character(0))
    c(coupling = ec$coupling == "expectationConsistent", ec = abs(ec$logI - ex))
  }))
  expect_true(all(err[, "coupling"] == 1))
  expect_lt(mean(err[, "ec"]), 0.05)
})

# the gradient of a fit on the unbounded scale of the hyperparameters, log or logit
.freeGradient <- function(fit) {
  g <- fit$gradient[names(fit$argument)]
  v <- fit$argument
  pos <- grepl("^(tau|concentration)_", names(g))
  prob <- grepl("^prob_", names(g))
  g[pos] <- g[pos] * v[pos]
  g[prob] <- g[prob] * v[prob] * (1 - v[prob])
  g
}

test_that("a prior beside the data enters a fit of sparse random effects", {
  m <- .sparse()
  obj <- normL2(m$data, m$g * Xt() * m$ps, errmodel = m$err)
  plain <- em(obj, randomEffects(m$ps), m$init, verbose = FALSE)
  target <- plain$argument[["A"]] + 0.3
  pr <- constraintL2(mu = c(A = target), sigma = 1e-3)
  fit <- em(obj + pr, randomEffects(m$ps), m$init, verbose = FALSE)
  expect_equal(fit$argument[["A"]], target, tolerance = 1e-3)
  expect_error(em(obj + constraintL2(mu = c(eta_k_L1 = 0), sigma = 1), randomEffects(m$ps), m$init),
               "fixed effects only")
})

test_that("normal and spike-and-slab random effects in one fit: the limit of a tiny Omega", {
  m <- .sparse()
  objOf <- function(p) normL2(m$data, m$g * Xt() * p, errmodel = m$err)
  # a tiny Omega pins the normal etas at 0: the fit of the slab ones alone
  tiny <- em(objOf(m$psm), randomEffects(m$psm), m$init, fixed = c(omega_A_A = log(1e-5), tau_ID = 0.5),
             verbose = FALSE)
  slab <- em(objOf(m$pks), randomEffects(m$pks), m$init, fixed = c(tau_ID = 0.5), verbose = FALSE)
  expect_equal(tiny$value, slab$value, tolerance = 1e-3)
  # sparse fits have no Hessian of the marginal, their standard errors come from the empirical
  # information of the group scores; CWRES without normal random effects is IWRES
  expect_null(slab$hessian)
  expect_equal(summary(slab)$seMethod, "scores")
  expect_equal(vcov(slab), solve(slab$information), ignore_attr = TRUE)
  expect_true(all(is.finite(summary(slab)$population$se)))
  expect_equal(residuals(slab, "CWRES"), residuals(slab, "IWRES"))
  expect_true(slab$slab$fixed[slab$slab$name == "tau_ID"])
  d <- do.call(rbind, lapply(names(m$data), function(cn) transform(m$data[[cn]], condition = cn)))
  d$lloq <- stats::quantile(d$value, 0.2)
  m1 <- normL2(as.datalist(d), m$g * Xt() * m$pks, errmodel = m$err, optBLOQ = "M1")
  expect_error(em(m1, randomEffects(m$pks), m$init, verbose = FALSE), "M3")
})

test_that("the group scores of a sparse marginal sum to the gradient, each the gradient of one group alone", {
  m <- .sparse()
  objOf <- function(p, conds = names(m$data)) normL2(m$data[conds], m$g * Xt() * p, errmodel = m$err)
  theta <- m$init + c(0.05, -0.02, 0.1)
  normal <- c(theta, omega_A_A = log(0.2))
  for (cs in list(list(m$ps, theta), list(m$psm, normal))) {
    at <- function(obj = NULL, scores = FALSE)
      .marginalAt(cs[[1]], cs[[2]], share = 0.4, tau = 0.5, obj = obj, scores = scores)
    r <- at(scores = TRUE)
    free <- replace(0 * r$gradient, names(r$freeGradient), r$freeGradient)
    expect_equal(colSums(r$scores)[names(r$gradient)] + free, r$gradient, tolerance = 1e-10)
    for (g in c("L1", "L3")) {
      one <- at(objOf(cs[[1]], rownames(m$tab)[m$tab$ID == g]))
      expect_equal(r$scores[g, ], one$gradient[colnames(r$scores)], tolerance = 1e-6)
    }
  }
})

test_that("em refuses entries of control$sparseFit, $quadrature, $saem and $innerControl it does not read", {
  m <- .sparse()
  obj <- normL2(m$data, m$g * Xt() * m$ps, errmodel = m$err)
  expect_error(em(obj, randomEffects(m$ps), m$init, control = list(sparseFit = list(coupling = "meanField"))),
               "unknown entries of `control\\$sparseFit`: coupling")
  expect_error(em(obj, randomEffects(m$ps), m$init, control = list(quadrature = list(level = 3))),
               "unknown entries of `control\\$quadrature`: level")
  expect_error(em(obj, randomEffects(m$ps), m$init, control = list(saem = list(nIter = 3))),
               "unknown entries of `control\\$saem`: nIter")
  expect_error(em(obj, randomEffects(m$ps), m$init, control = list(innerControl = list(ftol = 1e-8))),
               "unknown entries of `control\\$innerControl`: ftol")
  expect_error(em(obj, randomEffects(m$ps), m$init, control = list(sparseFit = list(trustControl = list(nope = 1)))),
               "unknown entries nope")
})

test_that("the joint trust region of sparse random effects runs the phases of trust()", {
  m <- .sparse()
  obj <- normL2(m$data, m$g * Xt() * m$ps, errmodel = m$err)
  ph <- list(list(hessianMethod = "gn", iterlim = 2), list(hessianMethod = "bfgs"))
  fit <- em(obj, randomEffects(m$ps), m$init, verbose = FALSE,
            control = list(sparseFit = list(trustControl = list(phases = ph))))
  expect_true(is.finite(fit$value))
  expect_lt(max(abs(.freeGradient(fit))), 0.05)
  expect_warning(em(obj, randomEffects(m$ps), m$init, verbose = FALSE,
                    control = list(sparseFit = list(maxOuter = 1L, trustControl = list(
                      phases = list(list(hessianMethod = "gn"), list(hessianMethod = "nope")))))),
                 "trust region stopped")
})

test_that("an eta that reaches no data of the group rests at the prior, with the marginal and gradient unchanged", {
  m <- .sparse()
  theta <- c(a = 1.1, b = 0.4, c = 0.2, omega_a_a = log(0.4))
  for (dist in names(m$pz)) {
    p <- m$pz[[dist]]
    obj <- normL2(m$zdata, m$gz * Xt() * p)
    d <- dMod2:::.emDesign(obj, randomEffects(p))
    expect_identical(unname(d$reach[, "eta_c"]), c(TRUE, TRUE, FALSE, FALSE, FALSE))
    expect_true(all(d$reach[, c("eta_a", "eta_b")]))
    # structurally unreached: the data do not move with it
    full <- c(theta, setNames(rep(0.1, length(d$etaNames)), d$etaNames))
    gr <- obj(full[setdiff(names(full), "omega_a_a")])$gradient
    expect_equal(unname(gr[d$etaNames[3:5, "eta_c"]]), c(0, 0, 0))
    expect_true(all(abs(gr[d$etaNames[1:2, "eta_c"]]) > 1e-3))
    for (mode in list(NULL, "forward-reverse")) {
      at <- function(reach) .marginalAt(p, theta, share = if (dist == "spikeSlab") 0.4 else 2, tau = 0.5,
                                        obj = obj, derivMode = mode, reach = reach)
      masked <- at(TRUE); all <- at(FALSE)
      expect_equal(masked$value, all$value, tolerance = 1e-8)
      expect_equal(masked$gradient, all$gradient, tolerance = 1e-6)
      expect_equal(unname(masked$modes[3:5, "eta_c"]), c(0, 0, 0))
    }
  }
})

test_that("a state that starts at 0 and that nothing drives passes no random effect on", {
  m <- .sparse()
  obj <- normL2(m$udata, m$gu * m$xu * m$pu)
  d <- dMod2:::.emDesign(obj, randomEffects(m$pu))
  stimulated <- !(d$groups %in% c("L1", "L2"))
  expect_identical(unname(d$reach[, "eta_s"]), stimulated)
  expect_identical(unname(d$reach[, "eta_w"]), stimulated)
  # the data do not move with an unreached eta, at a point away from the modes
  set.seed(2)
  full <- c(s = 0.7, w = 0.5, dg = 0.3, c0 = 0.1, setNames(rnorm(length(d$etaNames), 0, 0.3), d$etaNames))
  gr <- obj(full)$gradient
  expect_equal(unname(gr[d$etaNames[!stimulated, c("eta_s", "eta_w")]]), rep(0, 4))
  expect_true(all(abs(gr[d$etaNames[stimulated, c("eta_s", "eta_w")]]) > 1e-6))
})

test_that("a root of a steady state depends on the inputs the equations and conservation rows reach", {
  env <- list2env(list(
    dependent = c("A", "B", "C", "D", "E"),
    all_exprs = c(A = "k1 - k2 * A", B = "k2 * A - k3 * B", C = "k4 * TNF - k5 * C",
                  D = "-kf * D + kb * E", E = "kf * D - kb * E", T = "(D + E) - T"),
    parms_all = c("k1", "k2", "k3", "k4", "k5", "TNF", "kf", "kb", "T")))
  dep <- dMod2:::.pimplDepends(env, "TNF")
  expect_setequal(dep$k1, c("A", "B"))
  expect_setequal(dep$k3, "B")
  expect_setequal(dep$k5, "C")
  expect_false("k4" %in% names(dep))
  expect_setequal(dep$T, c("D", "E"))
  expect_setequal(dep$kf, c("D", "E"))
  # a subsystem without full structural rank: every root moves
  env$all_exprs <- c(A = "k1 - B", B = "k2 - B", C = "k9 - C")
  env$dependent <- c("A", "B", "C")
  env$parms_all <- c("k1", "k2", "k9")
  dep <- dMod2:::.pimplDepends(env, character(0))
  expect_setequal(dep$k9, c("A", "B", "C"))
  expect_setequal(dep$k2, c("A", "B"))
})

test_that("the trace gradients of the sweeps direction by direction are those along all tangents", {
  m <- .sparse()
  p <- m$pz$spikeSlab
  d <- dMod2:::.emDesign(normL2(m$zdata, m$gz * Xt() * p), randomEffects(p))
  full <- c(a = 1.1, b = 0.4, c = 0.2, setNames(seq(-0.2, 0.3, length.out = length(d$etaNames)), d$etaNames))
  # units of three and of two etas, and per unit a full matrix and one with a single entry
  units <- lapply(d$groups, function(g) list(eta = unname(d$etaNames[g, d$reach[g, ]]), conditions = d$conditions[[g]]))
  set.seed(1)
  Ms <- lapply(units, function(u) {
    k <- length(u$eta)
    A <- matrix(rnorm(k * k), k)
    list(crossprod(A), diag(c(1, rep(0, k - 1)), k))
  })
  wide <- dMod2:::.sweepUnits(d, full, NULL, units, Ms, narrow = FALSE)
  narrow <- dMod2:::.sweepUnits(d, full, NULL, units, Ms, narrow = TRUE)
  for (u in seq_along(units)) {
    nm <- wide[[u]]$reach
    expect_equal(narrow[[u]]$rows[, nm], wide[[u]]$rows[, nm], tolerance = 1e-10)
    for (j in 1:2) expect_equal(narrow[[u]]$traces[[j]][nm], wide[[u]]$traces[[j]][nm], tolerance = 1e-10)
  }
})

test_that("expectation consistent converges beside a coordinate without data and a pair the data tie", {
  # a pair with correlation 0.997, a coordinate the data do not see, large linear terms
  P <- diag(c(4, 4, 2, 1, 0, 3)); P[1, 2] <- P[2, 1] <- 0.997 * 4; P[3, 4] <- P[4, 3] <- 0.5
  b <- matrix(c(30, -28, 2, -1, 0, 1.5), 1)
  # no budget for the exact sum, so that expectation consistent couples the coordinates
  ec <- dMod2:::sparse_integral(list(P), b, rep(0L, 6), rep(0.3, 6), rep(0.8, 6), rep(NA_real_, 6),
                                matrix(TRUE, 1, 6), NULL, NULL, 0)
  expect_true(ec$ok)
  expect_identical(ec$coupling, "expectationConsistent")
  # the exact subset sum over all six coordinates
  exact <- log(sum(apply(as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), 6))), 1, function(F) {
    w <- prod(ifelse(F, 0.3, 0.7))
    if (!any(F)) return(w)
    A <- P[F, F, drop = FALSE] + diag(1 / 0.8^2, sum(F))
    w * exp(sum(b[F] * solve(A, b[F])) / 2) / sqrt(det(A)) / 0.8^sum(F)
  })))
  expect_lt(abs(ec$logI - exact), 1e-3)
  # the coordinate without data keeps the prior
  expect_equal(ec$mean[1, 5], 0, tolerance = 1e-8)
  expect_equal(ec$spike[1, 5], 0.7, tolerance = 1e-8)
})

test_that("a trial point at which expectation consistent fails is a failed evaluation, the coupling stays", {
  m <- .sparse()
  obj <- normL2(m$data, m$g * Xt() * m$psk2, errmodel = m$err)
  ec <- em(obj, randomEffects(m$psk2), m$init, verbose = FALSE)
  expect_equal(ec$coupling, "expectationConsistent")
  real <- dMod2:::.sparseIntegral
  first <- bad <- NULL
  after <- 0L
  # the integral fails at the first trial point whose slab standard deviation is not that of the start
  local_mocked_bindings(.sparseIntegral = function(P, b, lay, hyp, ...) {
    at <- hyp$tau
    if (is.null(first)) first <<- at
    if (is.null(bad) && !isTRUE(all.equal(at, first))) bad <<- at
    if (isTRUE(all.equal(at, bad))) return(NULL)
    if (!is.null(bad)) after <<- after + 1L
    real(P, b, lay, hyp, ...)
  }, .package = "dMod2")
  fit <- em(obj, randomEffects(m$psk2), m$init, verbose = FALSE)
  expect_false(is.null(bad))
  expect_gt(after, 0L)
  expect_equal(fit$coupling, "expectationConsistent")
  expect_true(fit$converged)
  # both fits end at the optimum to their tolerance, a few 1e-4 in -2 log L
  expect_lt(abs(fit$value - ec$value), 1e-3)
})


# Every set partition of the groups and the anchor (element G + 1): the Ewens prior with
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
  # in natural parameters t = a, h = a m; log Z of the center form is less sum(a m^2) / 2
  stem <- function(h, t, conc, tau) dMod2:::sparse_factor(1L, h, t, NA, tau, conc)
  ps <- stem(a * m, a, conc, tau)
  all <- .partitionsBrute(a, m, conc, tau)
  w <- vapply(all, `[[`, 0, "w")
  expect_equal(ps$logZ - sum(a * m^2) / 2, log(sum(w)), tolerance = 1e-9)
  expect_equal(ps$spike, colSums(w * t(vapply(all, `[[`, logical(4), "zero"))) / sum(w), tolerance = 1e-8)
  # without data every block integral is 1 and the prior sums to 1
  expect_equal(stem(0 * a, 0 * a, conc, tau)$logZ, 0, tolerance = 1e-12)
  # the prior does not depend on the order of the groups
  o <- c(3, 1, 4, 2)
  expect_equal(stem((a * m)[o], a[o], conc, tau)$logZ, ps$logZ, tolerance = 1e-12)
  logZ <- function(h, t, conc, tau) stem(h, t, conc, tau)$logZ
  expect_equal(ps$m1, numDeriv::grad(function(x) logZ(x, a, conc, tau), a * m), tolerance = 1e-7)
  expect_equal(-ps$m2 / 2, numDeriv::grad(function(x) logZ(a * m, x, conc, tau), a), tolerance = 1e-7)
  expect_equal(ps$dconc, numDeriv::grad(function(x) logZ(a * m, a, x, tau), conc), tolerance = 1e-7)
  expect_equal(ps$dtau, numDeriv::grad(function(x) logZ(a * m, a, conc, x), tau), tolerance = 1e-7)
  # the most probable partition
  expect_equal(ps$labels, all[[which.max(w)]]$labels)
  # the M-step of the concentration: the expected number of slab blocks
  th <- dMod2:::.ewensConcentration(1.3, c(4, 6))
  expect_equal(sum(th / (th + c(1:4, 1:6))), 1.3, tolerance = 1e-10)
})

test_that("the spike-and-slab integral of a group sums the subsets on the full curvature", {
  P <- matrix(c(4, 1.5, 1.5, 2), 2); m <- c(0.6, -0.3); prob <- c(0.3, 0.6); tau <- c(0.8, 0.5)
  st <- c("s1", "s2")
  lay <- data.frame(stem = st, block = 1:2, dist = "spikeSlab", probName = c("p1", "p2"),
                    concentrationName = NA_character_, tauName = c("t1", "t2"), prob = NA_real_,
                    concentration = NA_real_, tau = NA_real_, stringsAsFactors = FALSE)
  joint <- function(P, b, prob, tau) {
    hyp <- list(prob = setNames(prob, st), concentration = setNames(c(NA, NA), st), tau = setNames(tau, st))
    dMod2:::.sparseIntegral(list(structure(P, dimnames = list(st, st))), matrix(b, 1, dimnames = list("g", st)),
                            lay, hyp, c("p1", "p2", "t1", "t2"))
  }
  b <- drop(P %*% m)
  ss <- joint(P, b, prob, tau)
  expect_null(ss$coupling)
  quad <- function(e) exp(-sum((e - m) * (P %*% (e - m))) / 2)
  line <- function(k, v) { e <- c(0, 0); e[k] <- v; quad(e) }
  one <- function(k) integrate(function(v) vapply(v, function(x) line(k, x), 0) * dnorm(v, 0, tau[k]),
                               -Inf, Inf, rel.tol = 1e-12)$value
  both <- integrate(function(x) vapply(x, function(xi) integrate(function(y) vapply(y, function(yi)
    quad(c(xi, yi)), 0) * dnorm(y, 0, tau[2]), -Inf, Inf, rel.tol = 1e-11)$value, 0) * dnorm(x, 0, tau[1]),
    -Inf, Inf, rel.tol = 1e-10)$value
  I <- (1 - prob[1]) * (1 - prob[2]) * quad(c(0, 0)) + prob[1] * (1 - prob[2]) * one(1) +
    (1 - prob[1]) * prob[2] * one(2) + prob[1] * prob[2] * both
  expect_equal(ss$logI - sum(m * b) / 2, log(I), tolerance = 1e-8)
  # d log I / d b = E eta, d log I / d P = -E[eta eta'] / 2, the symmetric entries moving together
  logI <- function(P, b, prob, tau) joint(P, b, prob, tau)$logI
  expect_equal(unname(drop(ss$mean)), numDeriv::grad(function(x) logI(P, x, prob, tau), b), tolerance = 1e-7)
  sym <- function(x) matrix(c(x[1], x[2], x[2], x[3]), 2)
  dP <- numDeriv::grad(function(x) logI(sym(x), b, prob, tau), c(P[1, 1], P[1, 2], P[2, 2]))
  M2 <- -(ss$cov[[1]] + tcrossprod(drop(ss$mean))) / 2
  expect_equal(c(M2[1, 1], 2 * M2[1, 2], M2[2, 2]), dP, tolerance = 1e-7)
  expect_equal(unname(ss$dhyper[c("p1", "p2")]), numDeriv::grad(function(x) logI(P, b, x, tau), prob), tolerance = 1e-7)
  expect_equal(unname(ss$dhyper[c("t1", "t2")]), numDeriv::grad(function(x) logI(P, b, prob, x), tau), tolerance = 1e-7)
  # one coordinate is the partition of one group with the anchor, c / (1 + c) the slab share
  expect_equal(dMod2:::sparse_factor(0L, 1.2, 3, 0.3 / 1.3, 0.7, NA)$logZ,
               dMod2:::sparse_factor(1L, 1.2, 3, NA, 0.7, 0.3)$logZ)
})

test_that("13 coupled spike-and-slab coordinates: exact as one factor, beyond the budget expectation consistent", {
  K <- 13L; G <- 2L; st <- paste0("s", seq_len(K))
  # the subset sum without the bound on K of the joint factors
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
  lay <- data.frame(stem = st, block = 1L, dist = "spikeSlab", probName = "prob_ID",
                    concentrationName = NA_character_, tauName = "tau_ID", prob = NA_real_,
                    concentration = NA_real_, tau = NA_real_, stringsAsFactors = FALSE)
  hyp <- list(prob = setNames(rep(0.3, K), st), concentration = setNames(rep(NA, K), st),
              tau = setNames(rep(0.8, K), st))
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
  # without coupling the factors are exact, also without budget
  expect_equal(.centered(d0$P, d0$m, lay, hyp, character(0), budget = 0)$logI, logI(d0), tolerance = 1e-10)
  # coupled: the joint sum is exact, without budget expectation consistent is near it
  d <- draw(0.3)
  one <- .centered(d$P, d$m, lay, hyp, character(0))
  expect_null(one$coupling)
  expect_equal(one$logI, logI(d), tolerance = 1e-10)
  ec <- .centered(d$P, d$m, lay, hyp, character(0), budget = 0)
  expect_identical(ec$coupling, "expectationConsistent")
  expect_lt(abs(ec$logI - logI(d)), 0.02)
})

test_that("the spike-and-slab and cluster marginals have the exact gradient, normal etas beside", {
  m <- .sparse()
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

test_that("the marginal is that at the exact modes to first order, whatever the outer tolerances", {
  m <- .sparse()
  theta <- m$init + c(0.05, -0.02, 0.1)
  for (p in list(m$ps, m$psm)) {
    th <- if (any(randomEffects(p)$terms$dist == "normal")) c(theta, omega_A_A = log(0.2)) else theta
    exact <- .marginalAt(p, th, share = 0.4, tau = 0.5)
    loose <- .marginalAt(p, th, share = 0.4, tau = 0.5, tolControl = list(ftol = 1e-4, mtol = 1e-4, gtol = 1e-2))
    expect_equal(loose$value, exact$value, tolerance = 1e-9)
    expect_equal(loose$gradient, exact$gradient, tolerance = 1e-6)
    # modes moved off and held there: the value is off at second order in the shift only
    off <- vapply(c(2e-3, 1e-3), function(h)
      abs(suppressWarnings(.marginalAt(p, th, share = 0.4, tau = 0.5, warm = exact$modes + h, iterlim = 0L))$value -
            exact$value), 0)
    expect_gt(off[1] / off[2], 3.5)
  }
})

test_that("the curvature of the hyperparameters is that of the complete-data prior, the secant the M-step", {
  lay <- data.frame(stem = c("b", "c"), probName = c("prob_ID", NA), concentrationName = c(NA, "concentration_ID"),
                    tauName = c("tau_ID", "tau_ID_2"), stringsAsFactors = FALSE)
  kind <- dMod2:::.hyperKind(lay)
  G <- 5
  r <- list(slab = matrix(c(2.6, G, 1.3, 1.8, G, 0.9), 2, byrow = TRUE,
                          dimnames = list(lay$stem, c("free", "groups", "slab2"))))
  hyper <- c(prob_ID = 0.35, concentration_ID = 1.4, tau_ID = 0.6, tau_ID_2 = 0.8)
  # the expected -2 log prior of every hyperparameter as a function of the unbounded value
  prior <- list(
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
  m <- .sparse()
  set.seed(2)
  for (p in list(m$ps, m$psk)) {
    d <- dMod2:::.emDesign(normL2(m$data, m$g * Xt() * p, errmodel = m$err), randomEffects(p))
    center <- dMod2:::.emStart(d, m$init)
    draws <- t(replicate(200, dMod2:::.perturbStart(d, center, 3)))
    hyper <- grep("^(concentration|tau)_", colnames(draws))
    expect_true(all(draws[, hyper] > 0))
    prob <- grep("^prob_", colnames(draws))
    expect_true(all(draws[, prob] > 0 & draws[, prob] < 1))
    expect_gt(sd(draws[, "k"]), 1)
  }
  # a named sd perturbs the parameters it names, by name
  sd <- dMod2:::.perturbSd(c(s = 0.5, A = 2), center)
  expect_equal(sd[c("k", "A", "s")], c(k = 0, A = 2, s = 0.5))
  expect_error(dMod2:::.perturbSd(c(0.5, 2), center), "names")
  expect_error(dMod2:::.perturbSd(c(B = 1), center), "B")
})

test_that("SAEM with Dirichlet-process random effects selects the groups that differ", {
  m <- .sparse()
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
  m <- .sparse()
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
             control = list(sparseFit = list(maxOuter = 2L)), verbose = FALSE)
  expect_setequal(names(only$argument), c("concentration_ID", "tau_ID"))
  expect_true(is.finite(only$value))
  expect_lte(only$iterations, 2L)
  expect_error(vcov(fit), "Dirichlet-process")
  expect_equal(summary(fit)$seMethod, "none")
  # a center fit keeps the shared parameter shared and estimates the probability of a deviation
  center <- em(objOf(m$ps), randomEffects(m$ps), m$init, verbose = FALSE)
  expect_true(center$converged)
  expect_lt(max(abs(.freeGradient(center))), 0.02)
  share <- center$slab$estimate[center$slab$name == "prob_ID"]
  expect_true(share > 0 && share < 1)
  expect_gt(abs(center$etaModes["L3", "eta_k"]), 0.3)
  # a fit as init starts from the inner modes of that fit and evaluates to the same value
  again0 <- em(objOf(m$ps), randomEffects(m$ps), center, control = list(sparseFit = list(maxOuter = 0)),
               verbose = FALSE)
  expect_equal(again0$value, center$value, tolerance = 1e-8)
  again <- em(objOf(m$ps), randomEffects(m$ps), center$argument, verbose = FALSE)
  expect_gt(again$value, center$value - 1e-5)
  # five groups give scores in four directions, fewer than the five parameters
  expect_error(vcov(center), "do not span")
  expect_true(all(is.na(summary(center)$population$se)))
  A <- center$random[center$random$eta == "A", ]
  expect_true(all(A$mode == 0))
  expect_true(all(A$p0 > 0.5))
})

test_that("spike-and-slab and Dirichlet-process fits on inexact modes end at the nested optimum", {
  m <- .sparse()
  objOf <- function(p) normL2(m$data, m$g * Xt() * p, errmodel = m$err)
  for (p in list(m$ps, m$psk, m$psm)) {
    init <- if (any(randomEffects(p)$terms$dist == "normal")) c(m$init, omega_A_A = log(0.2)) else m$init
    nested <- em(objOf(p), randomEffects(p), init, verbose = FALSE)
    sqp <- em(objOf(p), randomEffects(p), init, verbose = FALSE, control = list(inner = "sqp"))
    expect_true(sqp$converged)
    expect_equal(sqp$value, nested$value, tolerance = 1e-6)
    expect_lt(max(abs(sqp$argument - nested$argument[names(sqp$argument)])), 1e-3)
    expect_equal(sqp$etaModes, nested$etaModes, tolerance = 1e-3)
    expect_true(is.integer(sqp$innerIterations) && sqp$innerIterations > 0L)
    expect_true(sqp$sqpIterations >= 1L && sqp$sqpIterations <= sqp$iterations)
    expect_true(nested$innerIterations > 0L)
  }
})

test_that("NPDE and the random-effect plots of a fit with spike-and-slab random effects alone", {
  m <- .sparse()
  fit <- em(normL2(m$data, m$g * Xt() * m$ps, errmodel = m$err), randomEffects(m$ps), m$init, verbose = FALSE)
  npde <- residuals(fit, type = "NPDE", nsim = 100, seed = 1)
  expect_length(npde, sum(vapply(m$data, nrow, 0L)))
  expect_true(all(is.finite(npde)))
  set.seed(1)
  expect_s3_class(plotRandomEffects(fit), "ggplot")
  d <- plotRandomEffects(fit, type = "inclusion")$data
  expect_equal(d$inclusion, 1 - fit$random$p0)
  expect_true(all(c("prob_ID", "tau_ID") %in% levels(plotTrace(fit)$data$panel)))
})

test_that("importance sampling on the model agrees with the marginal where the approximation is exact", {
  m <- .sparse()
  withr::local_seed(1)
  fit <- em(normL2(m$linear, m$gl * Xt() * m$pls), randomEffects(m$pls), c(a = 1, b = 0.5, omega_a_a = log(0.5)),
            verbose = FALSE)
  l <- logLik(fit, "importance", draws = 2000L)
  expect_equal(attr(l, "approximation"), as.numeric(logLik(fit)), tolerance = 1e-6)
  expect_lt(abs(as.numeric(l) - as.numeric(logLik(fit))), 4 * attr(l, "se"))
  expect_true(all(attr(l, "ess") > 200))
  fit$design <- NULL
  expect_error(logLik(fit, "importance"), "no design")
})
