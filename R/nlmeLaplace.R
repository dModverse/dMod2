# Laplace random effects with density lambda / 2 exp(-lambda |eta|) per coordinate and
# spike-and-slab ones, normal ones beside them; marginal, exact gradient and M-steps as derived
# in the guide chapter Group-specific effects.


## ---- normal-Laplace integral -------------------------------------------------

# sqrt(pi) x erfcx(x) - 1 by its asymptotic series, accurate to machine precision for x > 10
.erfcxSeries <- function(x) {
  t <- 1 / (2 * x^2); term <- 1; s <- 0
  for (k in 1:14) { term <- -term * (2 * k - 1) * t; s <- s + term }
  s
}

# log of erfcx(x) = exp(x^2) erfc(x) and its derivative, stable for every real x
.logErfcx <- function(x) {
  out <- log(2) + x^2 + stats::pnorm(-x * sqrt(2), log.p = TRUE)
  big <- x > 10
  out[big] <- log1p(.erfcxSeries(x[big])) - log(x[big] * sqrt(pi))
  out
}

.dLogErfcx <- function(x, logErfcx) {
  out <- 2 * x - 2 / sqrt(pi) * exp(-logErfcx)
  big <- x > 10
  s <- .erfcxSeries(x[big])
  out[big] <- 2 * x[big] * s / (1 + s)
  out
}

.logSumExp2 <- function(a, b) pmax(a, b) + log1p(exp(-abs(a - b)))

## log I(a, m, lambda), I = int lambda / 2 exp(-lambda |eta|) exp(-a / 2 (eta - m)^2) d eta, its
## partial derivatives and the posterior moments E eta, E|eta| and E (eta - m)^2, vectorised.
.normalLaplace <- function(a, m, lambda) {
  n <- max(length(a), length(m), length(lambda))
  a <- rep_len(a, n); m <- rep_len(m, n); lambda <- rep_len(lambda, n)
  cc <- sqrt(a / 2)
  u <- (lambda / a - m) * cc
  v <- (lambda / a + m) * cc
  lEu <- .logErfcx(u); lEv <- .logErfcx(v)
  lse <- .logSumExp2(lEu, lEv)
  logI <- log(lambda) - log(2) + 0.5 * (log(pi) - log(2 * a)) - a * m^2 / 2 + lse
  wPos <- exp(lEu - lse); wNeg <- exp(lEv - lse)
  dEu <- .dLogErfcx(u, lEu); dEv <- .dLogErfcx(v, lEv)
  dcda <- cc / (2 * a)
  dudm <- -cc; dvdm <- cc
  dudl <- cc / a; dvdl <- cc / a
  duda <- -lambda / a^2 * cc + (lambda / a - m) * dcda
  dvda <- -lambda / a^2 * cc + (lambda / a + m) * dcda
  dm <- -a * m + wPos * dEu * dudm + wNeg * dEv * dvdm
  dl <- 1 / lambda + wPos * dEu * dudl + wNeg * dEv * dvdl
  da <- -1 / (2 * a) - m^2 / 2 + wPos * dEu * duda + wNeg * dEv * dvda
  list(logI = logI, dm = dm, da = da, dlambda = dl,
       mean = m + dm / a, abs = 1 / lambda - dl, cen2 = -2 * da)
}


## ---- design pieces ---------------------------------------------------------------

## The Laplace, spike-and-slab and cluster coordinates: per stem its block, distribution, the
## names of its hyperparameters and their fixed values, NA when estimated.
.sparseLayout <- function(design) {
  blocks <- Filter(function(b) b$dist != "normal", design$blocks)
  name <- function(b, kind) if (kind %in% names(b$hyper)) b$hyper[[kind]] else NA_character_
  do.call(rbind, lapply(seq_along(blocks), function(i) {
    b <- blocks[[i]]
    data.frame(stem = b$stems, block = i, dist = b$dist, lambdaName = name(b, "lambda"),
               probName = name(b, "prob"), concentrationName = name(b, "concentration"),
               tauName = name(b, "tau"), lambda = b$lambda, prob = b$prob, concentration = b$concentration,
               tau = b$tau, stringsAsFactors = FALSE)
  }))
}

## The hyperparameters of every stem, lambda, theta and tau, from `pars` or the fixed values.
.stemHypers <- function(lay, pars) {
  value <- function(fixed, name) setNames(ifelse(is.na(fixed), unname(pars[name]), fixed), lay$stem)
  list(lambda = value(lay$lambda, lay$lambdaName), prob = value(lay$prob, lay$probName),
       concentration = value(lay$concentration, lay$concentrationName), tau = value(lay$tau, lay$tauName))
}

## The estimated hyperparameters of the layout.
.sparseHypers <- function(lay) {
  h <- c(lay$lambdaName[is.na(lay$lambda)], lay$probName[is.na(lay$prob)],
         lay$concentrationName[is.na(lay$concentration)], lay$tauName[is.na(lay$tau)])
  unique(h[!is.na(h)])
}

## The slab precision 2 / tau^2 of the -2 log prior of every spike-and-slab eta, named by eta.
.slabPrecision <- function(E, lay, hyp) {
  st <- lay$stem[lay$dist %in% c("spikeSlab", "dirichletProcess")]
  setNames(rep(2 / hyp$tau[st]^2, each = nrow(E)), as.vector(E[, st, drop = FALSE]))
}

## One data objective per condition, so that the groups evaluate in one batch.
.conditionObjectives <- function(design) {
  grid <- attr(design$objfn, "timesD", exact = TRUE)
  # the objective over all conditions, which .objEvalMany() calls once for many of them
  whole <- .normL2Objective(design$data, design$prdfn, design$errfn, grid, "data", design$bloq, NULL, design$df)
  lapply(setNames(nm = names(design$data)), function(cn) {
    o <- .normL2Objective(design$data[cn], design$prdfn, design$errfn,
                          if (is.list(grid)) grid[cn] else grid, "data", design$bloq, NULL, design$df)
    attr(o, "whole") <- whole
    o
  })
}

## The terms of the objective of `design` besides its data, such as priors and constraints; NULL
## without.
.extraObjective <- function(design) {
  terms <- Filter(function(t) is.null(attr(t, "l2spec", exact = TRUE)), .objTerms(design$objfn))
  if (length(terms)) Reduce(`+`, terms)
}

## The data -2 log L of every group at its own parameter vector, summed over its conditions.
.groupData <- function(design, objs, parsByGroup, deriv = TRUE, deriv2 = FALSE) {
  conds <- unlist(design$conditions, use.names = FALSE)
  owner <- rep(seq_along(design$groups), lengths(design$conditions))
  outs <- .objEvalMany(objs[conds], parsByGroup[owner], deriv = deriv, deriv2 = deriv2)
  lapply(seq_along(design$groups), function(i) Reduce(`+`, outs[owner == i]))
}


## ---- inner mode --------------------------------------------------------------------

## The joint mode of all etas at `theta` by trustL1: data -2 log L, normal prior, 2 lambda |eta|
## per Laplace coordinate and the slab alone, eta^2 / tau^2, per spike-and-slab eta.
.laplaceModes <- function(design, objs, prior, theta, hyp, start, control = list()) {
  lay <- .sparseLayout(design)
  E <- design$etaNames
  names <- as.vector(E)
  slab <- .slabPrecision(E, lay, hyp)
  objective <- function(eta, ...) {
    full <- c(theta, setNames(as.numeric(eta), names))
    outs <- .groupData(design, objs, rep(list(full), length(design$groups)))
    gr <- setNames(numeric(length(names)), names)
    he <- matrix(0, length(names), length(names), dimnames = list(names, names))
    v <- 0
    for (i in seq_along(outs)) {
      en <- E[i, ]
      v <- v + outs[[i]]$value
      gr[en] <- outs[[i]]$gradient[en]
      he[en, en] <- outs[[i]]$hessian[en, en]
    }
    if (!is.null(prior)) {
      pr <- prior(pars = full, deriv = TRUE)
      nn <- intersect(names, names(pr$gradient))
      v <- v + pr$value
      gr[nn] <- gr[nn] + pr$gradient[nn]
      he[nn, nn] <- he[nn, nn] + pr$hessian[nn, nn]
    }
    if (length(slab)) {
      es <- full[names(slab)]
      v <- v + sum(slab * es^2) / 2
      gr[names(slab)] <- gr[names(slab)] + slab * es
      he[cbind(names(slab), names(slab))] <- he[cbind(names(slab), names(slab))] + slab
    }
    list(value = v, gradient = gr, hessian = he)
  }
  laplace <- lay$stem[lay$dist == "laplace"]
  mu <- setNames(numeric(length(laplace) * nrow(E)), as.vector(E[, laplace]))
  args <- list(objfun = objective, parinit = setNames(as.vector(start), names), mu = mu,
               rinit = 1, rmax = 10, iterlim = control$iterlim %||% 100L,
               tolControl = list(ftol = 1e-10, mtol = 1e-10))
  if (length(mu)) args$lambda <- setNames(rep(2 * hyp$lambda[laplace], each = nrow(E)), names(mu))
  fit <- do.call(trustL1, args)
  matrix(fit$argument[names], nrow(E), dimnames = dimnames(E))
}


## ---- curvature -------------------------------------------------------------------------

## The Gauss-Newton Hessian of the data -2 log L of a group in its etas, sum over its rows of
## J' L J with J = d(p, s) / d eta and L the Hessian of the row loss, and its derivatives in the
## parameters `A` and in the etas, [K, K, direction].
.gnCurvature <- function(design, g, preds, eta, A, fixed) {
  frames <- lapply(design$conditions[[g]], function(cn) .volumeFrames(design, cn, preds[[cn]], fixed))
  rw <- .stackRows(Map(function(cn, fr) .volumeRows(design, cn, fr, eta, A, FALSE),
                       design$conditions[[g]], frames))
  lr <- .rowLossDerivs(rw$p, rw$s, rw$y, rw$cens, design$df)
  l <- function(i, j) lr[, i + 1L, j + 1L]
  Lpp <- l(2, 0); Lps <- l(1, 1); Lss <- l(0, 2)
  a <- rw$a; b <- rw$b
  K <- length(eta); tn <- nrow(a)
  sym <- function(X, Y, w) crossprod(X, w * Y) + crossprod(Y, w * X)
  H <- crossprod(a, Lpp * a) + sym(a, b, Lps) + crossprod(b, Lss * b)
  # dH along a direction with first derivatives (dp, ds) of p, s and (Da, Db) of a, b
  along <- function(dp, ds, Da, Db)
    crossprod(a, (l(3, 0) * dp + l(2, 1) * ds) * a) + sym(Da, a, Lpp) +
      sym(a, b, l(2, 1) * dp + l(1, 2) * ds) + sym(Da, b, Lps) + sym(a, Db, Lps) +
      crossprod(b, (l(1, 2) * dp + l(0, 3) * ds) * b) + sym(Db, b, Lss)
  dA <- array(0, c(K, K, length(A)), dimnames = list(eta, eta, A))
  dE <- array(0, c(K, K, K), dimnames = list(eta, eta, eta))
  for (q in seq_along(A))
    dA[, , q] <- along(rw$pA[, q], rw$sA[, q], matrix(rw$Wp[, , q], tn), matrix(rw$Ws[, , q], tn))
  for (j in seq_len(K))
    dE[, , j] <- along(a[, j], b[, j], matrix(rw$Vp[, , j], tn), matrix(rw$Vs[, , j], tn))
  dimnames(H) <- list(eta, eta)
  list(H = H, dA = dA, dE = dE)
}


## ---- sparse block ----------------------------------------------------------------------------

## The integral of the sparse coordinates (precision P[[i]], centre m[i, ] per group) by mean field
## over factors, each integrated exactly given the means of the others; the bound is stationary in
## the factors, so its derivatives in m, P and the hyperparameters hold them fixed.
.sparseMeanField <- function(P, m, start, lay, hyp, hypers, tol = 1e-11, maxit = 1000L) {
  G <- nrow(m); K <- ncol(m); stems <- colnames(m)
  sf <- .sparseFactors(P, stems, lay)
  r <- sf$rows; factors <- sf$factors; within <- sf$within; coupled <- sf$coupled
  mu <- start
  fits <- vector("list", length(factors))
  for (it in seq_len(maxit)) {
    before <- mu
    for (j in seq_along(factors)) {
      f <- factors[[j]]
      o <- setdiff(seq_len(K), f)
      shifted <- matrix(vapply(seq_len(G), function(i) {
        if (!length(o)) return(m[i, f])
        m[i, f] - drop(solve(P[[i]][f, f, drop = FALSE], P[[i]][f, o, drop = FALSE] %*% (mu[i, o] - m[i, o])))
      }, numeric(length(f))), G, length(f), byrow = TRUE)
      fits[[j]] <- .factorIntegral(f, P, shifted, lay[r[f], , drop = FALSE], hyp)
      mu[, f] <- fits[[j]]$mean
    }
    if (!coupled || max(abs(mu - before)) < tol) break
  }
  d <- mu - m
  cov <- lapply(seq_len(G), function(i) matrix(0, K, K, dimnames = list(stems, stems)))
  for (j in seq_along(factors)) for (i in seq_len(G)) cov[[i]][factors[[j]], factors[[j]]] <- fits[[j]]$cov[[i]]
  out <- .factorOutputs(fits, factors, m, lay, hypers)
  logI <- out$logI
  # the coupling between factors: sum_{b < c} d_b' P_bc d_c + 1/2 sum_b r_b' P_bb^-1 r_b
  dm <- matrix(0, G, K, dimnames = dimnames(m))
  dP <- vector("list", G)
  for (i in seq_len(G)) {
    Pc <- P[[i]]; Pc[within] <- 0
    Pw <- P[[i]]; Pw[!within] <- 0
    rr <- drop(Pc %*% d[i, ])
    logI <- logI + sum(d[i, ] * rr) / 2 + sum(rr * solve(Pw, rr)) / 2
    dm[i, ] <- drop(P[[i]] %*% d[i, ])
    dP[[i]] <- -(tcrossprod(d[i, ]) + cov[[i]]) / 2
  }
  c(list(logI = logI, dm = dm, dP = dP, mean = mu, cov = cov, sweeps = it, coupling = "meanField"),
    out[c("dhyper", "hyperGroups", "abs", "values", "clusters", "spike", "slab")])
}

## The sum of log I of the factors and their reported values, clusters, E|eta|, M-step statistics
## and derivatives in the hyperparameters, collected over the stems; `hyperGroups` holds the
## derivatives per group, NA where a factor couples the groups (Dirichlet process).
.factorOutputs <- function(fits, factors, m, lay, hypers) {
  G <- nrow(m); K <- ncol(m)
  values <- clusters <- spike <- matrix(NA_real_, G, K, dimnames = dimnames(m))
  absE <- matrix(0, G, K, dimnames = dimnames(m))
  slab <- matrix(0, nrow(lay), 3, dimnames = list(lay$stem, c("free", "groups", "slab2")))
  dhyper <- setNames(numeric(length(hypers)), hypers)
  hyperGroups <- matrix(0, G, length(hypers), dimnames = list(rownames(m), hypers))
  logI <- 0
  for (j in seq_along(factors)) {
    f <- factors[[j]]; ft <- fits[[j]]
    logI <- logI + ft$logI
    values[, f] <- ft$values; clusters[, f] <- ft$clusters; spike[, f] <- ft$spike
    if (!is.null(ft$abs)) absE[, f] <- ft$abs
    if (!is.null(ft$stats)) slab[rownames(ft$stats), ] <- ft$stats
    for (k in seq_along(ft$hyper))
      if (names(ft$hyper)[k] %in% hypers) dhyper[[names(ft$hyper)[k]]] <- dhyper[[names(ft$hyper)[k]]] + ft$hyper[[k]]
    if (is.null(ft$hyperGroups)) { hyperGroups[] <- NA_real_; next }
    for (k in which(colnames(ft$hyperGroups) %in% hypers)) {
      h <- colnames(ft$hyperGroups)[k]
      hyperGroups[, h] <- hyperGroups[, h] + ft$hyperGroups[, k]
    }
  }
  list(logI = logI, dhyper = dhyper, hyperGroups = hyperGroups, abs = absE, values = values,
       clusters = clusters, spike = spike, slab = slab)
}

## The factors of the sparse coordinates: every Laplace and Dirichlet-process stem alone, the
## spike-and-slab stems of a group jointly up to .slabJoint of them, beyond that each alone;
## `coupled` when the data precision links two factors.
.sparseFactors <- function(P, stems, lay) {
  r <- match(stems, lay$stem)
  joint <- which(lay$dist[r] == "spikeSlab")
  if (length(joint) > .slabJoint) joint <- integer(0)
  factors <- c(as.list(setdiff(seq_along(stems), joint)), if (length(joint)) list(joint))
  owner <- integer(length(stems))
  for (j in seq_along(factors)) owner[factors[[j]]] <- j
  within <- outer(owner, owner, `==`)
  list(factors = factors, rows = r, within = within,
       coupled = any(vapply(P, function(Pi) any(Pi[!within] != 0), TRUE)))
}

## Expectation consistent (Opper and Winther 2005): every factor exact under a diagonal Gaussian
## cavity, the coupling exact under diagonal sites with the same means and variances; stationary in
## the sites, so its derivatives hold them fixed. NULL when the damped iteration does not converge;
## one that stalls within 100 tol counts as converged.
.sparseExpectation <- function(P, m, lay, hyp, hypers, tol = 1e-8, maxit = 500L, start = NULL) {
  G <- nrow(m); K <- ncol(m); stems <- colnames(m)
  sf <- .sparseFactors(P, stems, lay)
  Pm <- matrix(vapply(seq_len(G), function(i) drop(P[[i]] %*% m[i, ]), numeric(K)), G, K, byrow = TRUE)
  Pd <- lapply(P, function(Pi) { Pi <- as.matrix(Pi); storage.mode(Pi) <- "double"; Pi })
  storage.mode(Pm) <- "double"
  cavity <- function(Lam, gam) ec_cavity(Pd, Lam, gam, Pm, stems, dimnames(m))
  integrate <- function(f, cv) .factorIntegral(f, lapply(seq_len(G), function(i) diag(cv$t[i, ], K)),
                                               cv$n[, f, drop = FALSE], lay[sf$rows[f], , drop = FALSE], hyp)
  Lam <- gam <- matrix(0, G, K, dimnames = dimnames(m))
  if (!is.null(start) && identical(dim(start$Lam), dim(Lam)) && all(.sparseCavityOk(P, start$Lam, start$gam))) {
    Lam[] <- start$Lam; gam[] <- start$gam
  }
  damp <- 0.5; last <- Inf; converged <- FALSE
  for (it in seq_len(maxit)) {
    delta <- 0
    for (f in sf$factors) {
      cv <- cavity(Lam, gam)
      ft <- integrate(f, cv)
      v <- matrix(vapply(ft$cov, function(C) diag(as.matrix(C)), numeric(length(f))), G, length(f), byrow = TRUE)
      newL <- 1 / v - cv$t[, f, drop = FALSE]
      newG <- ft$mean / v - cv$t[, f, drop = FALSE] * cv$n[, f, drop = FALSE]
      if (!all(is.finite(newL)) || !all(is.finite(newG))) return(NULL)
      step <- damp
      repeat {
        L2 <- Lam; G2 <- gam
        L2[, f] <- Lam[, f] + step * (newL - Lam[, f])
        G2[, f] <- gam[, f] + step * (newG - gam[, f])
        if (all(cavity(L2, G2)$t > 0)) break
        step <- step / 2
        if (step < 1e-6) return(NULL)
      }
      delta <- max(delta, abs(newL - Lam[, f]) / (1 + abs(Lam[, f])), abs(newG - gam[, f]) / (1 + abs(gam[, f])))
      Lam <- L2; gam <- G2
    }
    damp <- if (delta > last) max(damp / 2, 0.01) else min(damp * 1.2, 0.7)
    last <- delta
    if (delta < tol) { converged <- TRUE; break }
  }
  if (!converged && !(last < 100 * tol)) return(NULL)
  # log I = log Z of the factors under the cavity + log Z of the coupling under the sites
  #         - log Z of the cavity under the sites
  cv <- cavity(Lam, gam)
  fits <- lapply(sf$factors, function(f) integrate(f, cv))
  out <- .factorOutputs(fits, sf$factors, m, lay, hypers)
  logR <- 0
  dm <- matrix(0, G, K, dimnames = dimnames(m)); dP <- vector("list", G)
  for (i in seq_len(G)) {
    A <- P[[i]] + diag(Lam[i, ], K); h <- Pm[i, ] + gam[i, ]
    logR <- logR + K / 2 * log(2 * pi) - 0.5 * as.numeric(determinant(A)$modulus) - 0.5 * sum(m[i, ] * Pm[i, ]) +
      0.5 * sum(h * solve(A, h))
    d <- cv$mu[i, ] - m[i, ]
    dm[i, ] <- drop(P[[i]] %*% d)
    dP[[i]] <- -(tcrossprod(d) + cv$cov[[i]]) / 2
  }
  logS <- sum(0.5 * log(2 * pi) - 0.5 * log(cv$t + Lam) - cv$t * cv$n^2 / 2 +
                (cv$t * cv$n + gam)^2 / (2 * (cv$t + Lam)))
  c(list(logI = out$logI + logR - logS, dm = dm, dP = dP, mean = cv$mu, cov = cv$cov, sweeps = it,
         coupling = "expectationConsistent", sites = list(Lam = Lam, gam = gam)),
    out[c("dhyper", "hyperGroups", "abs", "values", "clusters", "spike", "slab")])
}

## Whether sites leave every cavity precision positive.
.sparseCavityOk <- function(P, Lam, gam) {
  unlist(lapply(seq_along(P), function(i) 1 / diag(solve(P[[i]] + diag(Lam[i, ], ncol(Lam)))) - Lam[i, ] > 0))
}

## The sparse integral: exact where the factors are not coupled, otherwise by expectation
## consistent or mean field as `coupling` says; NULL where expectation consistent does not
## converge. Expectation consistent starts from the sites `sites` where given, otherwise from zero.
.sparseIntegral <- function(P, m, start, lay, hyp, hypers, coupling = "expectationConsistent", sites = NULL) {
  if (identical(coupling, "expectationConsistent") && .sparseFactors(P, colnames(m), lay)$coupled) {
    ec <- if (!is.null(sites)) .sparseExpectation(P, m, lay, hyp, hypers, start = sites)
    if (is.null(ec)) ec <- .sparseExpectation(P, m, lay, hyp, hypers)
    return(ec)
  }
  .sparseMeanField(P, m, start, lay, hyp, hypers)
}

## The integral of Dirichlet-process stems by sequential importance sampling
## (.dirichletImportance), with the standard error of log I; the moments and derivatives of the
## mean-field result are not formed, they stay 0 or NA.
.importanceField <- function(P, m, lay, hyp, hypers, particles) {
  G <- nrow(m); K <- ncol(m); stems <- colnames(m)
  if (!all(lay$dist[match(stems, lay$stem)] == "dirichletProcess"))
    stop("em: importance sampling integrates Dirichlet-process random effects only.", call. = FALSE)
  is <- .dirichletImportance(P, m, hyp$concentration[stems], hyp$tau[stems], particles)
  none <- matrix(NA_real_, G, K, dimnames = dimnames(m))
  list(logI = is$logI, se = is$se, dm = 0 * m, dP = lapply(P, function(x) 0 * x),
       dhyper = setNames(numeric(length(hypers)), hypers), mean = none,
       cov = lapply(seq_len(G), function(i) matrix(NA_real_, K, K, dimnames = list(stems, stems))),
       abs = 0 * m, values = none, clusters = none, spike = none,
       slab = matrix(0, nrow(lay), 3, dimnames = list(lay$stem, c("free", "groups", "slab2"))))
}

## The exact integral of one mean-field factor at the centres `mt` (groups x its stems): log Z,
## per group the posterior mean and covariance, E|eta| of Laplace stems, the Fisher derivatives in
## the hyperparameters, also per group where the factor is a product over the groups, the reported
## values and clusters and the M-step statistics.
.factorIntegral <- function(f, P, mt, rows, hyp) {
  G <- nrow(mt); st <- rows$stem
  a <- vapply(P, function(Pi) Pi[f[1], f[1]], 0)
  one <- function(v) matrix(v, G, 1)
  if (rows$dist[1] == "laplace") {
    nl <- .normalLaplace(a, mt[, 1], hyp$lambda[[st]])
    return(list(logI = sum(nl$logI), mean = one(nl$mean), cov = as.list(pmax(nl$cen2 - (nl$mean - mt[, 1])^2, 0)),
                abs = nl$abs, hyper = setNames(sum(nl$dlambda), rows$lambdaName),
                hyperGroups = matrix(nl$dlambda, G, 1, dimnames = list(NULL, rows$lambdaName)),
                values = one(NA_real_), clusters = one(NA_real_), spike = one(NA_real_)))
  }
  if (rows$dist[1] == "dirichletProcess") {
    ss <- .slabPartition(a, mt[, 1], hyp$concentration[[st]], hyp$tau[[st]])
    return(list(logI = ss$logI, mean = one(ss$mean), cov = as.list(pmax(ss$cen2 - (ss$mean - mt[, 1])^2, 0)),
                hyper = setNames(c(ss$dconcentration, ss$dtau), c(rows$concentrationName, rows$tauName)),
                values = one(ss$values), clusters = one(ss$labels), spike = one(ss$spike),
                stats = matrix(ss$stats, 1, dimnames = list(st, NULL))))
  }
  per <- lapply(seq_len(G), function(i)
    .slabSubsets(P[[i]][f, f, drop = FALSE], mt[i, ], hyp$prob[st], hyp$tau[st]))
  byGroup <- function(field) matrix(vapply(per, `[[`, numeric(length(f)), field), G, length(f), byrow = TRUE)
  incl <- byGroup("incl"); sel <- byGroup("selected")
  # the groups away from 0 of a symbol, numbered
  clusters <- apply(sel, 2, function(x) ifelse(x > 0, cumsum(x > 0), 0))
  stats <- cbind(colSums(incl), G, colSums(byGroup("slab2")), deparse.level = 0)
  rownames(stats) <- st
  list(logI = sum(vapply(per, `[[`, 0, "logI")), mean = byGroup("mean"), cov = lapply(per, `[[`, "cov"),
       hyper = c(setNames(colSums(byGroup("dprob")), rows$probName), setNames(colSums(byGroup("dtau")), rows$tauName)),
       hyperGroups = `colnames<-`(cbind(byGroup("dprob"), byGroup("dtau")), c(rows$probName, rows$tauName)),
       values = byGroup("values"), clusters = matrix(clusters, G, length(f)), spike = byGroup("spike"),
       stats = stats)
}


## ---- marginal ------------------------------------------------------------------------------

## The approximate marginal -2 log L at `theta`: the normal etas on the quadratic model of the joint
## mode, the sparse coordinates on the Schur complement by mean field; its exact gradient in
## `active` and `hypers`, a model Hessian, the moments and the selected clusters. With
## `scores`, the gradient split into the rows of the groups and that of the conditions without
## random effects, with their Hessian.
.laplaceMarginal <- function(design, objs, prior, theta, hyp, active, hypers, fixed, warm,
                             grad = TRUE, control = list(), scores = FALSE) {
  lay <- .sparseLayout(design)
  nl <- design$normal
  E <- design$etaNames
  G <- nrow(E); K <- ncol(E)
  isL <- colnames(E) %in% lay$stem
  iN <- which(!isL); iL <- which(isL)
  modes <- .laplaceModes(design, objs, prior, theta, hyp, warm, control)
  full <- c(theta, setNames(as.vector(modes), as.vector(E)))
  # forward over reverse: the data Hessians Gauss-Newton, the exact rows and the derivatives of the
  # curvature from sweeps along the etas, also for the value alone
  fr <- identical(control$derivMode, "forward-reverse")
  at <- .groupData(design, objs, rep(list(full), G), deriv = TRUE, deriv2 = grad && !fr)
  pr <- if (!is.null(prior)) prior(pars = full, fixed = fixed, deriv = TRUE)
  free <- if (length(design$freeConditions))
    .objEvalMany(objs[design$freeConditions], rep(list(full), length(design$freeConditions)),
                 deriv = TRUE, deriv2 = grad && !fr)
  # the terms of the objective besides the data, priors and constraints, enter as the free conditions
  extra <- .extraObjective(design)
  if (!is.null(extra)) free <- c(free, list(extra(pars = full, fixed = fixed, deriv = TRUE, deriv2 = grad && !fr)))
  A <- active
  units <- lapply(seq_len(G), function(i) list(eta = unname(E[i, ]), conditions = design$conditions[[design$groups[i]]]))
  if (fr) {
    fw <- .sweepForward(design, full, fixed, units)
    curv <- lapply(.unitCurvature(design, fw, units), function(H) list(H = H))
  } else {
    conds <- unlist(design$conditions, use.names = FALSE)
    times <- lapply(conds, function(cn) .emGrid(design, cn))
    preds <- setNames(.predictMany(design$prdfn, times = times, parsList = rep(list(full), length(conds)),
                                   conditions = conds, fixed = fixed, deriv = TRUE, deriv2 = TRUE), conds)
    curv <- lapply(seq_len(G), function(i) .gnCurvature(design, design$groups[i], preds, E[i, ], A, fixed))
  }
  gr <- matrix(t(vapply(seq_len(G), function(i) at[[i]]$gradient[E[i, ]], numeric(K))), G,
               dimnames = dimnames(E))

  # the normal etas: H_NN of every group with the prior precision 2 Omega^-1; the sparse block
  # of `Hs` is its Schur complement, the curvature left once the normal etas are integrated out
  cd <- if (length(iN)) .choleskyDerivatives(nl, .choleskyFactor(nl, c(full, fixed)))
  HN <- lapply(curv, function(cv) if (length(iN)) cv$H[iN, iN, drop = FALSE] + 2 * cd$Oinv)
  volume <- sum(vapply(HN, function(H) if (is.null(H)) 0 else as.numeric(determinant(H / 2)$modulus), 0))
  # the sparse block of every group: its curvature S, the Schur complement of the normal block,
  # a normal factor with precision S / 2 and centre mode - S^-1 g
  S <- lapply(seq_len(G), function(i) {
    H <- curv[[i]]$H
    s <- H[iL, iL, drop = FALSE]
    if (length(iN)) s <- s - H[iL, iN, drop = FALSE] %*% solve(HN[[i]], H[iN, iL, drop = FALSE])
    s + diag(1e-12, length(iL))
  })
  KL <- length(iL)
  bySparse <- function(f) matrix(vapply(seq_len(G), f, numeric(KL)), G, KL, byrow = TRUE,
                                 dimnames = list(rownames(E), colnames(E)[iL]))
  gL <- gr[, iL, drop = FALSE]
  u <- bySparse(function(i) solve(S[[i]], gL[i, ]))
  sp <- if (identical(control$integral, "importance"))
    .importanceField(lapply(S, `/`, 2), modes[, iL, drop = FALSE] - u, lay, hyp, hypers, control$particles %||% 500L)
  else .sparseIntegral(lapply(S, `/`, 2), modes[, iL, drop = FALSE] - u, modes[, iL, drop = FALSE],
                       lay, hyp, hypers, control$coupling %||% "expectationConsistent", control$sites)
  if (is.null(sp)) return(list(infeasible = TRUE))
  value <- sum(vapply(at, `[[`, 0, "value")) + (pr$value %||% 0) + volume - sum(gL * u) / 2 - 2 * sp$logI +
    sum(vapply(free, `[[`, 0, "value"))
  selected <- modes
  selected[, iL] <- ifelse(is.na(sp$values), modes[, iL], sp$values)
  clusters <- spike <- matrix(NA_real_, G, K, dimnames = dimnames(E))
  clusters[, iL] <- sp$clusters; spike[, iL] <- sp$spike
  mean <- absE <- matrix(0, G, K, dimnames = dimnames(E))
  mean[, iL] <- sp$mean; absE[, iL] <- sp$abs
  # the normal etas given the sparse ones are normal with a mean linear in them
  cov <- lapply(seq_len(G), function(i) {
    V <- matrix(0, K, K, dimnames = list(colnames(E), colnames(E)))
    V[iL, iL] <- sp$cov[[i]]
    if (length(iN)) {
      B <- solve(HN[[i]], curv[[i]]$H[iN, iL, drop = FALSE])
      VS <- V[iL, iL, drop = FALSE]
      VS[is.na(VS)] <- 0
      V[iN, iN] <- 2 * solve(HN[[i]]) + B %*% VS %*% t(B)
      V[iN, iL] <- -B %*% VS
      V[iL, iN] <- t(V[iN, iL, drop = FALSE])
    }
    V
  })
  var <- matrix(vapply(cov, diag, numeric(K)), G, K, byrow = TRUE, dimnames = dimnames(E))
  out <- list(value = value, modes = modes, mean = mean, var = var, cov = cov, abs = absE, slab = sp$slab,
              selected = selected, clusters = clusters, spike = spike, valueSE = 2 * (sp$se %||% 0),
              coupling = if (.sparseFactors(S, colnames(E)[iL], lay)$coupled) sp$coupling, sites = sp$sites)
  if (!grad) return(out)

  # derivatives of the value in g, S and the mode on the sparse block, the mean-field factors
  # held fixed, then through the Schur complement on the cross and normal blocks of H
  w <- bySparse(function(i) solve(S[[i]], sp$dm[i, ]))
  Q <- length(A)
  own <- matrix(0, G, Q, dimnames = list(rownames(E), A))
  names <- as.vector(t(E))
  dEta <- setNames(numeric(G * K), names)
  Hex <- matrix(0, G * K, G * K, dimnames = list(names, names))
  Hcross <- matrix(0, G * K, Q, dimnames = list(names, A))
  schur <- matrix(0, Q, Q, dimnames = list(A, A))
  # per group the matrices the derivatives of the curvature are contracted with: MHi on all etas,
  # H_NN^-1 on the normal ones
  contraction <- lapply(seq_len(G), function(i) {
    MHi <- matrix(0, K, K, dimnames = list(E[i, ], E[i, ]))
    MHi[iL, iL] <- tcrossprod(u[i, ]) / 2 - tcrossprod(u[i, ], w[i, ]) - tcrossprod(w[i, ], u[i, ]) - sp$dP[[i]]
    Hi <- NULL
    if (length(iN)) {
      Hi <- solve(HN[[i]])
      B <- curv[[i]]$H[iL, iN, drop = FALSE] %*% Hi
      Ms <- MHi[iL, iL, drop = FALSE]
      MHi[iL, iN] <- -Ms %*% B
      MHi[iN, iL] <- t(MHi[iL, iN, drop = FALSE])
      MHi[iN, iN] <- crossprod(B, Ms %*% B)
    }
    list(MHi = MHi, Hi = Hi)
  })
  # forward over reverse: the exact eta rows of the data Hessian and the trace gradients
  sw <- if (fr) .sweepUnits(design, full, fixed, units, fw = fw, Ms = lapply(contraction, function(ct) {
    M2 <- 0 * ct$MHi
    if (!is.null(ct$Hi)) M2[iN, iN] <- ct$Hi
    list(ct$MHi, M2)
  }))
  pick <- function(v, nm) { out <- setNames(numeric(length(nm)), nm); h <- intersect(nm, names(v)); out[h] <- v[h]; out }
  for (i in seq_len(G)) {
    en <- E[i, ]
    o <- at[[i]]
    hh <- o$hessian
    an <- intersect(A, names(o$gradient))
    cv <- curv[[i]]
    Mgi <- Meta <- numeric(K)
    Mgi[iL] <- 2 * w[i, ] - u[i, ]
    Meta[iL] <- -2 * sp$dm[i, ]
    MHi <- contraction[[i]]$MHi
    Hi <- contraction[[i]]$Hi
    # the exact rows of the group's data Hessian at its etas, and tr(M dH) along A and the etas
    if (fr) {
      rowsAt <- function(cols) {
        out <- matrix(0, K, length(cols), dimnames = list(en, cols))
        h <- intersect(cols, colnames(sw[[i]]$rows))
        out[, h] <- sw[[i]]$rows[, h]
        out
      }
      exactE <- rowsAt(en); exactA <- rowsAt(an)
      dMA <- pick(sw[[i]]$traces[[1]], A); dME <- pick(sw[[i]]$traces[[1]], en)
      dHA <- pick(sw[[i]]$traces[[2]], A); dHE <- pick(sw[[i]]$traces[[2]], en)
    } else {
      exactE <- hh[en, en, drop = FALSE]; exactA <- hh[en, an, drop = FALSE]
      dMA <- vapply(seq_len(Q), function(q) sum(MHi * cv$dA[, , q]), 0)
      dME <- vapply(seq_len(K), function(j) sum(MHi * cv$dE[, , j]), 0)
      if (length(iN)) {
        dHA <- vapply(seq_len(Q), function(q) sum(Hi * cv$dA[iN, iN, q]), 0)
        dHE <- vapply(seq_len(K), function(j) sum(Hi * cv$dE[iN, iN, j]), 0)
      }
    }
    own[i, an] <- own[i, an] + o$gradient[an] + drop(Mgi %*% exactA) + dMA[match(an, A)]
    dEta[en] <- o$gradient[en] + drop(Mgi %*% exactE) + dME + Meta
    if (length(iN)) {
      # d log|H_NN| along the parameters, the Cholesky entries and the etas; the prior precision
      # in H_NN also enters the Schur complement
      own[i, ] <- own[i, ] + dHA
      dEta[en] <- dEta[en] + dHE
      ch <- intersect(A, nl$cholPars)
      MN <- Hi + MHi[iN, iN, drop = FALSE]
      own[i, ch] <- own[i, ch] + vapply(ch, function(cp) 2 * sum(MN * cd$dOinv[[match(cp, nl$cholPars)]]), 0)
    }
    Hex[en, en] <- exactE
    Hcross[en, an] <- exactA
    schur[an, an] <- schur[an, an] + hh[an, an]
  }
  explicit <- colSums(own)
  if (!is.null(pr)) {
    pa <- intersect(A, names(pr$gradient)); pe <- intersect(names, names(pr$gradient))
    explicit[pa] <- explicit[pa] + pr$gradient[pa]
    dEta[pe] <- dEta[pe] + pr$gradient[pe]
    Hex[pe, pe] <- Hex[pe, pe] + pr$hessian[pe, pe]
    Hcross[pe, pa] <- Hcross[pe, pa] + pr$hessian[pe, pa]
    schur[pa, pa] <- schur[pa, pa] + pr$hessian[pa, pa]
  }
  freeGradient <- setNames(numeric(Q), A)
  freeHessian <- matrix(0, Q, Q, dimnames = list(A, A))
  for (o in free) {
    an <- intersect(A, names(o$gradient))
    freeGradient[an] <- freeGradient[an] + o$gradient[an]
    freeHessian[an, an] <- freeHessian[an, an] + o$hessian[an, an]
  }
  explicit <- explicit + freeGradient
  schur <- schur + freeHessian

  # the model Hessian: the Schur complement over all etas, the sparse etas with the effective prior
  # curvature of the integral, the sites of expectation consistent or per factor 2 C^-1 - S of its
  # posterior covariance C; with normal etas alone it is the Gauss-Newton Schur complement
  posterior <- Hex
  factors <- .sparseFactors(lapply(S, `/`, 2), colnames(E)[iL], lay)$factors
  for (i in seq_len(G)) {
    en <- E[i, iL]
    if (!is.null(sp$sites)) {
      posterior[en, en] <- posterior[en, en] + diag(2 * sp$sites$Lam[i, ], KL)
      next
    }
    V <- sp$cov[[i]]
    if (anyNA(V)) next
    for (f in factors) {
      e <- eigen(V[f, f, drop = FALSE], symmetric = TRUE)
      inv <- e$vectors %*% (t(e$vectors) / pmax(e$values, 1e-12 * max(1, e$values)))
      posterior[en[f], en[f]] <- posterior[en[f], en[f]] + 2 * inv - S[[i]][f, f, drop = FALSE]
    }
  }
  model <- if (Q) tryCatch(schur - crossprod(Hcross, solve(posterior, Hcross)), error = function(e) NULL)

  slabPrec <- .slabPrecision(E, lay, hyp)
  Hex[cbind(names(slabPrec), names(slabPrec))] <- Hex[cbind(names(slabPrec), names(slabPrec))] + slabPrec

  # the shift of the mode: normal etas, slab etas and free clusters move with the parameters; the
  # curvature of the mode is block diagonal by groups, so the shift of a group is its own
  C <- .freeClusters(modes, E, lay, hyp)
  gradient <- explicit
  hyperGrad <- -2 * sp$dhyper
  byGroup <- if (scores) {
    hg <- sp$hyperGroups %||% matrix(NA_real_, G, length(hypers), dimnames = list(NULL, hypers))
    cbind(own, -2 * hg[, hypers, drop = FALSE])
  }
  if (scores && !is.null(pr)) byGroup[, A] <- byGroup[, A] + .priorGroups(nl, E, full, fixed, cd, A)
  out$modeShift <- matrix(0, length(names), Q, dimnames = list(names, A))
  if (ncol(C$map)) {
    map <- C$map[names, , drop = FALSE]
    CH <- crossprod(map, Hex %*% map)
    md <- drop(crossprod(map, dEta))
    owner <- outer(seq_len(G), rep(seq_len(G), each = K)[max.col(abs(t(map)), "first")], `==`) * 1
    if (Q) {
      shift <- solve(CH, crossprod(map, Hcross))
      out$modeShift <- -map %*% shift
      gradient <- gradient - drop(md %*% shift)
      if (scores) byGroup[, A] <- byGroup[, A] - owner %*% (md * shift)
      schur <- schur - crossprod(Hcross, map %*% shift)
    }
    for (h in hypers) {
      r <- crossprod(map, ifelse(!is.na(C$hyper[names]) & C$hyper[names] == h, C$sign[names], 0))
      v <- md * drop(solve(CH, r))
      hyperGrad[h] <- hyperGrad[h] - sum(v)
      if (scores) byGroup[, h] <- byGroup[, h] - drop(owner %*% v)
    }
  }
  out$gradient <- c(gradient, hyperGrad)
  out$hessian <- model %||% schur
  if (scores) out[c("scores", "freeGradient", "freeHessian")] <- list(byGroup, freeGradient, freeHessian)
  out
}

## The gradient of the normal prior of every group, eta' Omega^-1 eta + log|Omega|, in the
## Cholesky parameters among `A`, a row per group.
.priorGroups <- function(nl, E, full, fixed, cd, A) {
  ch <- intersect(A, nl$cholPars)
  out <- matrix(0, nrow(E), length(A), dimnames = list(rownames(E), A))
  if (!length(ch)) return(out)
  all <- c(full, fixed)
  for (i in seq_len(nrow(E))) {
    eta <- all[E[i, nl$stems]]
    out[i, ch] <- vapply(match(ch, nl$cholPars), function(m) sum(eta * (cd$dOinv[[m]] %*% eta)), 0) +
      2 * unname(nl$isDiag)[match(ch, nl$cholPars)]
  }
  out
}

## The free variables of the mode as columns over all etas (normal and slab etas, Laplace ones
## away from 0), and per eta the derivative of its prior gradient in its
## hyperparameter with that hyperparameter's name.
.freeClusters <- function(modes, E, lay, hyp, tol = 1e-8) {
  names <- as.vector(t(E))
  sgn <- setNames(numeric(length(names)), names)
  hyper <- setNames(rep(NA_character_, length(names)), names)
  cols <- list()
  unit <- function(at) { v <- setNames(numeric(length(names)), names); v[at] <- 1; v }
  for (k in seq_len(ncol(E))) {
    x <- modes[, k]; en <- E[, k]
    r <- match(colnames(E)[k], lay$stem)
    if (is.na(r)) { cols <- c(cols, lapply(en, unit)); next }
    if (lay$dist[r] %in% c("spikeSlab", "dirichletProcess")) {
      cols <- c(cols, lapply(en, unit))
      hyper[en] <- lay$tauName[r]
      sgn[en] <- -4 * x / hyp$tau[[lay$stem[r]]]^3
      next
    }
    hyper[en] <- lay$lambdaName[r]
    cols <- c(cols, lapply(en[abs(x) > tol], unit))
    sgn[en] <- 2 * sign(x)
  }
  map <- if (length(cols)) do.call(cbind, cols) else matrix(0, length(names), 0)
  rownames(map) <- names
  list(map = map, sign = sgn, hyper = hyper)
}


## ---- fit -------------------------------------------------------------------------------------

## Sparse random effects: one trust region over all parameters, the hyperparameters unbounded, on
## the Schur complement under the posterior covariance of the sparse etas and per hyperparameter
## the curvature whose Newton step is its M-step, with a BFGS correction for the missing information.
.fitLaplace <- function(design, start, fixed, control = list(), verbose = FALSE) {
  .requireM3(design, NULL, "FOCEI")
  lay <- .sparseLayout(design)
  objs <- .conditionObjectives(design)
  prior <- if (design$normal$K) .normalPrior(design)
  hypers <- intersect(.sparseHypers(lay), names(start))
  off <- intersect(attr(.extraObjective(design), "parameters"), c(.sparseHypers(lay), as.vector(design$etaNames)))
  if (length(off))
    stop("em: with Laplace, spike-and-slab or Dirichlet-process random effects the terms beside the data ",
         "take the fixed effects only; ", paste(off, collapse = ", "), " are hyperparameters or random effects.",
         call. = FALSE)
  kind <- .hyperKind(lay)[hypers]
  structural <- setdiff(names(start), hypers)
  E <- design$etaNames
  warm <- matrix(0, nrow(E), ncol(E), dimnames = dimnames(E))
  warmAt <- shift <- NULL
  hyper <- start[hypers]
  par <- start[structural]
  lc <- control$laplace %||% list()
  lc$derivMode <- .secondOrderMode(design, design$objfn, par, fixed, control[setdiff(names(control), "secondOrderCorrection")],
                                   units = lapply(design$groups, function(g)
                                     list(eta = unname(E[g, ]), conditions = design$conditions[[g]])))
  maxOuter <- lc$maxOuter %||% 200L
  jointControl <- .trustControl(list(rinit = 1, rmax = 10, tolControl = list(ftol = 1e-10, mtol = 1e-5)),
                                lc$trustControl, label = "trustControl")
  # the inner start at p: the last modes moved by their first-order shift along the parameters
  predicted <- function(p) {
    if (is.null(shift)) return(warm)
    q <- intersect(colnames(shift), names(p))
    out <- warm
    at <- match(rownames(shift), as.vector(E))
    out[at] <- out[at] + drop(shift[, q, drop = FALSE] %*% (p[q] - warmAt[q]))
    out
  }
  # one coupling of the sparse random effects for the whole fit, so that the objective is one
  # function: expectation consistent where it converges at the start, mean field otherwise; a trial
  # point at which expectation consistent does not converge is a failed evaluation
  requested <- match.arg(lc$coupling %||% "expectationConsistent", c("expectationConsistent", "meanField"))
  coupling <- requested
  sites <- NULL
  # every evaluation starts from the modes and sites of the last accepted point, so that the
  # objective near it does not depend on the trials before
  marginal <- function(p, grad = TRUE, commit = TRUE, scores = FALSE) {
    r <- .laplaceMarginal(design, objs, prior, c(p, fixed), .stemHypers(lay, c(hyper, fixed)), names(p),
                          hypers, fixed, predicted(p), grad = grad, scores = scores,
                          control = modifyList(lc, list(coupling = coupling, sites = sites)))
    if (isTRUE(r$infeasible)) return(r)
    if (commit) accept(p, r)
    r
  }
  accept <- function(p, r) {
    warm <<- r$modes
    warmAt <<- p
    if (!is.null(r$modeShift)) shift <<- r$modeShift
    if (!is.null(r$sites)) sites <<- r$sites
  }
  fellBack <- identical(coupling, "expectationConsistent") && isTRUE(marginal(par, grad = FALSE)$infeasible)
  if (fellBack) coupling <- "meanField"
  rows <- list()
  record <- function(phase, value, dpar)
    rows[[length(rows) + 1L]] <<- data.frame(iter = length(rows) + 1L, phase = phase, OFV = value,
                                             deltaPar = dpar, t(hyper))
  # the iterations are the recorded steps; the start of every trust run is recorded without one
  steps <- function() sum(vapply(rows, function(r) !is.na(r$deltaPar), TRUE))
  # the objective in the structural parameters and the unbounded hyperparameters; its model Hessian
  # is the model Hessian of the marginal and the secant curvature of every hyperparameter to its
  # M-step, plus a BFGS correction from the points that lowered the value
  best <- Inf; bestAt <- bestGrad <- correction <- NULL
  jointObjective <- function(x, ...) {
    p <- x[names(par)]
    hyper <<- .hyperFromFree(kind, x[hypers])
    r <- marginal(p, commit = FALSE)
    if (isTRUE(r$infeasible)) return(.failedEvaluation(x))
    n <- names(x)
    base <- matrix(0, length(n), length(n), dimnames = list(n, n))
    base[names(p), names(p)] <- r$hessian[names(p), names(p)]
    base[cbind(hypers, hypers)] <- .hyperSecant(lay, kind, r, hyper)
    g <- c(r$gradient[names(p)], r$gradient[hypers] * .hyperJacobian(kind, hyper))
    if (is.null(correction)) correction <<- 0 * base
    if (r$value < best) {
      accept(p, r)
      if (!is.null(bestAt)) {
        record("trust", r$value, max(abs(x - bestAt)))
        correction <<- .secantCorrection(base, correction, x - bestAt, g - bestGrad)
      } else record("trust", r$value, NA_real_)
      if (verbose) message(sprintf("em(laplace) iteration %d: -2 log L %.6f", length(rows), r$value))
      best <<- r$value; bestAt <<- x; bestGrad <<- g
    }
    objlist(value = r$value, gradient = g, hessian = base + correction)
  }
  # one run of the trust region from the current point: the trust result, NULL with a warning on
  # an error
  jointRun <- function(rinit = jointControl$rinit) {
    x0 <- c(par, .hyperToFree(kind, hyper))
    best <<- Inf; bestAt <<- bestGrad <<- correction <<- NULL
    jt <- tryCatch(.innerTrust(modifyList(jointControl, list(objfun = jointObjective, parinit = x0, rinit = rinit,
                                                             iterlim = max(maxOuter - steps(), 1L)))),
                   error = function(e) {
                     warning("em: the trust region stopped: ", conditionMessage(e), call. = FALSE)
                     NULL
                   })
    hyper <<- .hyperFromFree(kind, x0[hypers])
    if (is.list(jt)) {
      par <<- jt$argument[names(par)]
      hyper <<- .hyperFromFree(kind, jt$argument[hypers])
    }
    jt
  }

  # a run that stops short starts again with a fresh correction while the runs lower -2 log L
  estimated <- length(par) + length(hyper) > 0L
  conv <- !estimated
  last <- Inf
  mtol <- jointControl$tolControl$mtol %||% 1e-5
  while (estimated && steps() < maxOuter) {
    jt <- jointRun()
    if (is.null(jt)) { conv <- FALSE; break }
    conv <- isTRUE(jt$converged)
    if (conv || !(jt$value < last - mtol)) break
    last <- jt$value
  }
  final <- marginal(par, scores = TRUE)
  if (isTRUE(final$infeasible) || !is.finite(final$value %||% NA_real_))
    stop("em: the integral over the sparse random effects is not finite at the last point of the ",
         "fit; a start nearer the optimum, or another coupling in control$laplace, may avoid it.",
         call. = FALSE)
  coupled <- if (!is.null(final$coupling)) if (fellBack) "meanFieldFallback" else final$coupling
  # the model Hessian approximates the curvature of the integral by Gauss-Newton, so the fit has none; with
  # independent groups the empirical information is G / (G - 1) times the sum of the outer products
  # of their centred scores of log L, the conditions without random effects add half their Hessian
  information <- if (!is.null(final$scores) && !anyNA(final$scores)) {
    G <- nrow(final$scores)
    centred <- sweep(final$scores, 2, colMeans(final$scores)) / 2
    I <- if (G > 1L) G / (G - 1) * crossprod(centred) else 0 * crossprod(centred)
    I[names(par), names(par)] <- I[names(par), names(par)] + final$freeHessian / 2
    I
  }
  list(argument = c(par, hyper), value = final$value, gradient = final$gradient,
       hessian = NULL, information = information, etaModes = final$selected,
       etaCov = setNames(final$cov, design$groups),
       clusters = final$clusters, spike = final$spike, converged = conv, iterations = steps(),
       stageTrace = do.call(rbind, rows), coupling = coupled)
}

## A failed evaluation for trust(): it shrinks the radius, three in a row end the run.
.failedEvaluation <- function(x)
  objlist(value = Inf, gradient = setNames(0 * as.numeric(x), names(x)),
          hessian = matrix(0, length(x), length(x), dimnames = list(names(x), names(x))))

## The BFGS update of B = base + correction for the step s with gradient change y, returned as the
## new correction; a pair with s'y <= 1e-8 |s| |y| is skipped, a correction that leaves B
## indefinite is dropped.
.secantCorrection <- function(base, correction, s, y) {
  B <- base + correction
  Bs <- drop(B %*% s)
  sy <- sum(s * y); sBs <- sum(s * Bs)
  if (sy <= 1e-8 * sqrt(sum(s^2) * sum(y^2)) || sBs <= 0) return(correction)
  updated <- B - tcrossprod(Bs) / sBs + tcrossprod(y) / sy
  if (min(eigen(updated, symmetric = TRUE, only.values = TRUE)$values) <= 0) return(0 * correction)
  updated - base
}

## The kind of every hyperparameter of the layout: lambda, prob, concentration or tau.
.hyperKind <- function(lay) {
  kinds <- c("lambda", "prob", "concentration", "tau")
  out <- unlist(lapply(kinds, function(k) setNames(rep(k, nrow(lay)), lay[[paste0(k, "Name")]])))
  out[!is.na(names(out)) & !duplicated(names(out))]
}

## The hyperparameters on their unbounded scale, logit for a probability and log otherwise, back,
## and the derivative of a hyperparameter in its unbounded value.
.hyperToFree <- function(kind, h) {
  pr <- kind[names(h)] == "prob"
  u <- h; u[!pr] <- log(h[!pr]); u[pr] <- stats::qlogis(h[pr])
  u
}
.hyperFromFree <- function(kind, u) {
  pr <- kind[names(u)] == "prob"
  h <- u; h[!pr] <- exp(u[!pr]); h[pr] <- stats::plogis(u[pr])
  h
}
.hyperJacobian <- function(kind, h) {
  pr <- kind[names(h)] == "prob"
  d <- h; d[pr] <- h[pr] * (1 - h[pr])
  d
}

## The complete-data information of every hyperparameter on its unbounded scale: the second
## derivative of the expected -2 log prior of the random effects at the posterior moments `r`,
## at least 1e-8 so that it seeds a positive definite model.
.hyperInformation <- function(lay, kind, r, hyper) {
  info <- vapply(names(hyper), function(h) {
    v <- hyper[[h]]
    k <- kind[[h]]
    st <- lay$stem[lay[[paste0(k, "Name")]] %in% h]
    S <- r$slab[st, , drop = FALSE]
    switch(k,
      lambda = 2 * v * sum(r$abs[, st]),
      tau = 4 * sum(S[, "slab2"]) / v^2,
      prob = 2 * sum(S[, "groups"]) * v * (1 - v),
      concentration = { i <- unlist(lapply(S[, "groups"], seq_len)); 2 * sum(v * i / (v + i)^2) })
  }, 0)
  pmax(info, 1e-8)
}

## The curvature of every hyperparameter on its unbounded scale whose Newton step is its M-step:
## the slope of the complete-data gradient between the value and the M-step, Q being convex on that
## scale; the complete-data information where the two coincide.
.hyperSecant <- function(lay, kind, r, hyper) {
  info <- .hyperInformation(lay, kind, r, hyper)
  d <- .hyperToFree(kind, hyper) - .hyperToFree(kind, .updateSparse(lay, r, hyper))
  slope <- vapply(names(hyper), function(h) {
    v <- hyper[[h]]
    k <- kind[[h]]
    st <- lay$stem[lay[[paste0(k, "Name")]] %in% h]
    S <- r$slab[st, , drop = FALSE]
    switch(k,
      lambda = 2 * (v * sum(r$abs[, st]) - length(r$abs[, st])),
      tau = 2 * (sum(S[, "free"]) - sum(S[, "slab2"]) / v^2),
      prob = -2 * (sum(S[, "free"]) - sum(S[, "groups"]) * v),
      concentration = { i <- unlist(lapply(S[, "groups"], seq_len)); -2 * (sum(S[, "free"]) - sum(v / (v + i))) })
  }, 0)
  pmax(ifelse(abs(d) > 1e-8 & slope * d > 0, slope / d, info), 1e-8)
}

## M-steps from the posterior moments `r`: lambda = n / sum E|eta|, tau^2 = E[sum of squared slab
## values] / E[slab blocks], prob = E[slab etas] / groups, and the concentration c with
## E[slab blocks] = sum_i c / (c + i) over the groups of every stem.
.updateSparse <- function(lay, r, current) {
  for (h in names(current)) {
    st <- lay$stem[lay$lambdaName %in% h]
    if (length(st)) {
      current[[h]] <- length(r$abs[, st]) / sum(r$abs[, st])
      next
    }
    rows <- r$slab[lay$stem[lay$probName %in% h | lay$concentrationName %in% h | lay$tauName %in% h], ,
                   drop = FALSE]
    S <- colSums(rows)
    if (h %in% lay$tauName) {
      if (S[["free"]] > 1e-8) current[[h]] <- sqrt(S[["slab2"]] / S[["free"]])
    } else if (h %in% lay$probName) {
      current[[h]] <- min(max(S[["free"]] / S[["groups"]], 1e-8), 1 - 1e-8)
    } else {
      current[[h]] <- .ewensConcentration(S[["free"]], rows[, "groups"])
    }
  }
  current
}

## The Ewens concentration c > 0 whose expected number of new blocks over stems with `groups`
## groups each equals `free`: sum_s sum_{i <= G_s} c / (c + i) = free, by Newton on log c (the
## left side is increasing in c); `free` at most half a block below its limit, so c stays finite.
.ewensConcentration <- function(free, groups) {
  i <- unlist(lapply(groups, seq_len))
  free <- min(max(free, 1e-8), length(i) - 0.5)
  lt <- log(max(free, 0.1))
  for (it in 1:100) {
    th <- exp(lt)
    f <- sum(th / (th + i)) - free
    df <- sum(th * i / (th + i)^2)
    step <- f / df
    lt <- lt - max(min(step, 2), -2)
    if (abs(step) < 1e-12) break
  }
  exp(lt)
}
