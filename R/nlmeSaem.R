# SAEM: Metropolis draws of the etas, Robbins-Monro averages of the sufficient statistics, the
# mu-referenced fixed effects by GLS and Omega in closed form, the other parameters on the averaged
# quadratic model; the estimate is the mean of the late iterates, the marginal at the estimate by quadrature.


## -2 log likelihood of data rows with prediction p and sigma s, normal or Student-t, censored
## rows by M3, with the constants of the FOCEI kernel.
.rowLoss <- function(p, s, y, censored, df = Inf) {
  w <- (p - y) / s
  a <- !censored
  out <- if (is.finite(df))
    sum((df + 1) * log1p(w[a]^2 / df) + 2 * log(s[a])) -
      2 * sum(a) * (lgamma((df + 1) / 2) - lgamma(df / 2) - log(df * pi) / 2)
  else sum(w[a]^2 + log(2 * pi) + 2 * log(s[a]))
  if (!any(censored)) return(out)
  out - 2 * sum(if (is.finite(df)) stats::pt(-w[censored], df, log.p = TRUE)
                else stats::pnorm(-w[censored], log.p = TRUE))
}

## -2 log p(y_g, eta_g) of every group at the rows of `eta`, the value of the FOCEI kernel per
## group: the data term of every condition in one evaluation without derivatives, in C++ where the
## model allows it, else from the predictions and the row losses in R.
.groupDensity <- function(design, meta, fixed, cores = 1L) {
  errfn <- .errorModel(design)
  node <- .fnNode(errfn)
  blocks <- unlist(lapply(meta$fast_meta, `[[`, "blocks"), recursive = FALSE)
  group <- rep(seq_along(meta$fast_meta), vapply(meta$fast_meta, function(m) length(m$blocks), 1L))
  conditions <- vapply(blocks, function(b) as.character(b$condition), "")
  times <- lapply(blocks, `[[`, "times")
  etaNames <- lapply(meta$fast_meta, `[[`, "eta_names")
  dataTerm <- function(parsList) {
    d <- .chainLocal(design$chain, conditions, parsList, rep(list(fixed), length(parsList)), FALSE, FALSE, cores,
                     design$bloq, design$df)
    if (is.null(d)) return(NULL)
    v <- vapply(d, function(x) as.numeric(x$value), 0)
    vapply(seq_along(etaNames), function(g) sum(v[group == g]), 0)
  }
  function(full, eta, Oinv, logDet) {
    parsList <- lapply(seq_along(blocks), function(j) {
      p <- full; p[etaNames[[group[j]]]] <- eta[group[j], ]; p
    })
    value <- dataTerm(parsList)
    if (!is.null(value)) return(value + rowSums((eta %*% Oinv) * eta) + logDet)
    preds <- .predictMany(design$prdfn, times = times, parsList = parsList, conditions = conditions,
                          fixed = fixed, deriv = FALSE, cores = cores)
    inner <- lapply(preds, function(q) attr(q, "parameters"))
    errs <- as.list(if (!is.null(node))
      .evalNode(node, .bundle(conds = conditions, out = as.list(preds), pars = inner,
                              fixed = vector("list", length(preds))), FALSE, FALSE, NULL, cores)
      else lapply(seq_along(preds), function(j)
        errfn(out = preds[[j]], pars = inner[[j]], conditions = conditions[j])[[1L]]))
    value <- numeric(length(etaNames))
    for (j in seq_along(blocks)) {
      b <- blocks[[j]]
      err <- errs[[j]] %||% .dataSigmaErrfn(design$data)(out = preds[[j]], conditions = conditions[j])[[1L]]
      p <- preds[[j]][cbind(b$t_idx_in_pred, b$o_idx_in_pred)]
      s <- err[cbind(b$t_idx_in_err, b$o_idx_in_err)]
      recorded <- !is.na(b$sigma_data)
      s[recorded] <- b$sigma_data[recorded]
      value[group[j]] <- value[group[j]] + .rowLoss(p, s, b$y_data, seq_along(p) > b$n_aloq, design$df)
    }
    value + rowSums((eta %*% Oinv) * eta) + logDet
  }
}


## Metropolis steps on all groups at once for the -2 log density `density(eta)` with current
## values `value`: a random walk of size `step` per group or, with `modes` and lower Cholesky
## factors `chol` (NULL where a group has none), the independence proposal N(mode, chol chol').
.metropolisGroups <- function(density, eta, value, n, step = NULL, modes = NULL, chol = NULL) {
  G <- nrow(eta); K <- ncol(eta)
  independent <- !is.null(modes)
  movable <- if (independent) !vapply(chol, is.null, TRUE) else rep(TRUE, G)
  logProposal <- function(x) if (!independent) numeric(G) else vapply(seq_len(G), function(g)
    if (movable[g]) -sum(forwardsolve(chol[[g]], x[g, ] - modes[g, ])^2) / 2 else 0, 0)
  current <- logProposal(eta)
  accepted <- numeric(G)
  for (s in seq_len(n)) {
    z <- matrix(stats::rnorm(G * K), G, K)
    u <- stats::runif(G)
    prop <- if (!independent) eta + step * z else matrix(vapply(seq_len(G), function(g)
      if (movable[g]) modes[g, ] + drop(chol[[g]] %*% z[g, ]) else eta[g, ], numeric(K)),
      G, K, byrow = TRUE)
    pv <- tryCatch(density(prop), error = function(e) rep(Inf, G))
    lp <- logProposal(prop)
    ok <- movable & is.finite(pv) & log(u) < (value - pv) / 2 + current - lp
    eta[ok, ] <- prop[ok, ]
    value[ok] <- pv[ok]
    current[ok] <- lp[ok]
    accepted <- accepted + ok
  }
  list(eta = eta, value = value, accept = accepted / max(n, 1L))
}

## Lower Cholesky factors of scale (H_g / 2)^-1, the conditional covariances of the Laplace
## approximation; NULL for a group whose H_g is not positive definite.
.proposalFactors <- function(H, scale)
  lapply(H, function(h) tryCatch(t(chol(chol2inv(chol(h / 2)) * scale)), error = function(e) NULL))


## t statistics of the slopes of the columns of `x` over the rows, 0 for a constant column.
.slopeStatistics <- function(x) {
  i <- seq_len(nrow(x)) - (nrow(x) + 1) / 2
  b <- colSums(i * x) / sum(i^2)
  res <- sweep(x, 2, colMeans(x)) - outer(i, b)
  se <- sqrt(colSums(res^2) / (nrow(x) - 2) / sum(i^2))
  ifelse(se > 0, b / se, 0)
}

## Monte Carlo covariance of the column means of `x` by ten batch means.
.batchMeansCov <- function(x) {
  size <- nrow(x) %/% 10L
  rows <- utils::tail(seq_len(nrow(x)), 10L * size)
  means <- rowsum(x[rows, , drop = FALSE], rep(seq_len(10L), each = size)) / size
  stats::cov(means) / 10
}

## The expected increase of -2 log L from the Monte Carlo error of SAEM: the covariance `V` of
## the mean of the single-draw M-steps through the EM map, (I - A)^-1 = Iobs^-1 Ic, measured in
## the observed information Iobs.
.monteCarloGap <- function(Iobs, Ic, V) sum(diag(Ic %*% V %*% Ic %*% solve(Iobs)))

## Half the Hessian of -2 log L of the quadrature marginal at level `level`, the inverse
## covariance of the estimate, with absolute eigenvalues where it is not positive definite.
.marginalInformation <- function(q, psi, fixed, level) {
  q$rebuild(psi, level)
  out <- tryCatch(q$value(psi, fixed = fixed, deriv2 = TRUE), error = function(e) q$value(psi, fixed = fixed))
  ev <- eigen((out$hessian + t(out$hessian)) / 4, symmetric = TRUE)
  lambda <- pmax(abs(ev$values), 1e-8 * max(abs(ev$values)))
  structure(ev$vectors %*% (t(ev$vectors) * lambda), dimnames = list(names(psi), names(psi)))
}

## The Robbins-Monro average of the quadratic models 0.5 t'A t - c't of the complete-data objective
## at the drawn random effects, `out` the average value, gradient and Hessian at `theta`.
.averageQuadratic <- function(model, out, theta, gamma) {
  A <- unclass(out$hessian)
  c <- drop(A %*% theta) - out$gradient
  if (is.null(model)) return(list(A = A, c = c))
  list(A = model$A + gamma * (A - model$A), c = model$c + gamma * (c - model$c))
}

## The Newton step at `theta` on the averaged model, at most `radius` long; the fixed point of the step zeroes
## the averaged gradient, the score of the marginal by Fisher's identity.
.quadraticStep <- function(model, theta, radius = 1) {
  e <- eigen(model$A, symmetric = TRUE)
  d <- pmax(abs(e$values), 1e-8 * max(abs(e$values), 1e-300))
  step <- -drop(e$vectors %*% (crossprod(e$vectors, drop(model$A %*% theta) - model$c) / d))
  size <- sqrt(sum(step^2))
  if (size > radius) step <- step * radius / size
  setNames(theta + step, names(theta))
}


## ---- Dirichlet-process random effects ------------------------------------------------------

## One Gibbs sweep (Neal 2000, algorithm 8, `m` auxiliary clusters) over the clusters `s$z` of the
## groups of one Dirichlet-process stem, 0 the population value, values `s$v`; `dataAt(v)` is
## -2 log p(y_g | eta) of every group with the eta of the stem at v[g].
.dpGibbs <- function(s, dataAt, concentration, tau, m = 2L) {
  G <- length(s$z)
  vals <- c(0, s$v)
  L <- matrix(vapply(vals, function(v) dataAt(rep(v, G)), numeric(G)), G)
  aux <- matrix(stats::rnorm(G * m, 0, tau), G, m)
  La <- matrix(vapply(seq_len(m), function(j) dataAt(aux[, j]), numeric(G)), G)
  sizes <- tabulate(s$z + 1L, length(vals))
  sizes[1] <- sizes[1] + 1L
  for (g in sample.int(G)) {
    k <- s$z[g] + 1L
    sizes[k] <- sizes[k] - 1L
    a <- aux[g, ]; la <- La[g, ]
    # a group alone in a cluster offers the value of that cluster as one auxiliary
    if (k > 1L && sizes[k] == 0L) { a[1] <- vals[k]; la[1] <- L[g, k] }
    open <- which(sizes > 0L)
    logw <- c(log(sizes[open]) - L[g, open] / 2, log(concentration / m) - la / 2)
    j <- sample.int(length(logw), 1L, prob = exp(logw - max(logw)))
    if (j <= length(open)) {
      s$z[g] <- open[j] - 1L
      sizes[open[j]] <- sizes[open[j]] + 1L
      next
    }
    v <- a[j - length(open)]
    slot <- which(sizes == 0L & seq_along(sizes) > 1L)[1]
    if (is.na(slot)) { vals <- c(vals, NA); sizes <- c(sizes, 0L); L <- cbind(L, NA); slot <- length(vals) }
    if (!identical(vals[slot], v)) { vals[slot] <- v; L[, slot] <- dataAt(rep(v, G)) }
    s$z[g] <- slot - 1L
    sizes[slot] <- 1L
  }
  keep <- which(sizes[-1] > 0L) + 1L
  s$z <- match(s$z + 1L, c(1L, keep)) - 1L
  s$v <- vals[keep]
  s
}

## Random-walk Metropolis steps on the values of the clusters of one stem, all clusters at once;
## the acceptance rate per step.
.dpValues <- function(s, dataAt, tau, step, n = 3L) {
  K <- length(s$v)
  if (!K) return(list(s = s, accept = NA_real_))
  byCluster <- function(x) vapply(seq_len(K), function(j) sum(x[s$z == j]), 0)
  current <- dataAt(c(0, s$v)[s$z + 1L])
  accepted <- 0
  for (it in seq_len(n)) {
    prop <- s$v + step * stats::rnorm(K)
    pg <- dataAt(c(0, prop)[s$z + 1L])
    delta <- byCluster(pg) - byCluster(current)
    ok <- log(stats::runif(K)) < -delta / 2 - (prop^2 - s$v^2) / (2 * tau^2)
    s$v[ok] <- prop[ok]
    moved <- s$z %in% which(ok)
    current[moved] <- pg[moved]
    accepted <- accepted + mean(ok)
  }
  list(s = s, accept = accepted / n)
}

## The state of a stem in canonical form: cluster 0 for the population value, the others
## numbered in the order of their first group, their values in that order.
.canonicalClusters <- function(s) {
  free <- unique(s$z[s$z > 0L])
  list(z = ifelse(s$z > 0L, match(s$z, free), 0L), v = s$v[free])
}

.fitSaem <- function(design, start, fixed, control = list(), verbose = FALSE) {
  lay <- design$normal
  sc <- control$saem %||% list()
  maxBurnin <- sc$nBurnin %||% 200L
  maxAveraging <- sc$nEM %||% 300L
  nMcmc <- sc$nMcmc %||% 3L
  nLaplace <- sc$nLaplace %||% 3L
  laplaceEvery <- sc$laplaceEvery %||% 5L
  laplaceScale <- sc$laplaceScale %||% 1.5
  window <- sc$window %||% 25L
  mcTol <- sc$mcTol %||% 0.05
  cores <- getOption("dMod2.cores", 1L)
  E <- design$etaNames[, lay$stems, drop = FALSE]
  G <- nrow(E); K <- ncol(E)
  meta <- .foceiMeta(design, start, fixed)
  .requireM3(design, meta, "SAEM")
  quad <- modifyList(control, control$quadrature %||% list())
  cm1 <- .trustControl(list(rinit = 1, rmax = 10, iterlim = 1, tolControl = list(ftol = 1e-6, mtol = 1e-6)),
                       sc$cm1Control, label = "cm1Control")

  # phi = eta + D_g beta held while beta moves, for the random intercepts with X constant in the
  # groups: GLS for those of .muReference(), the trust step for those whose symbol also enters
  # conditions without a group
  grouped <- unlist(design$conditions, use.names = FALSE)
  groupedDesign <- design
  groupedDesign$data <- design$data[grouped]
  ref <- .muReference(groupedDesign)
  gls <- names(.muReference(design))
  refPars <- unlist(lapply(unname(ref), `[[`, "pars"))
  free <- .intersect(unlist(lapply(unname(ref[gls]), `[[`, "pars")), names(start))
  structural <- .setdiff(names(start), lay$cholPars)
  trustPars <- .setdiff(structural, free)
  partial <- .intersect(trustPars, refPars)
  D <- lapply(design$groups, function(g) {
    m <- matrix(0, K, length(refPars), dimnames = list(lay$stems, refPars))
    for (e in ref) m[e$stem, e$pars] <- e$X[g, ]
    m
  })
  shift <- function(beta) matrix(vapply(D, function(m) drop(m %*% beta[refPars]), numeric(K)), G, K,
                                 byrow = TRUE, dimnames = dimnames(E))
  # Omega and the GLS beta that maximize the expected prior given the averages s1 of phi_g and
  # s2 of phi_g phi_g' / G; the diagonal of S at least `minVar` while annealing
  mStep <- function(s1, s2, beta, minVar = NULL) {
    for (it in 1:3) {
      sh <- shift(beta)
      S <- s2 - (crossprod(sh, s1) + crossprod(s1, sh) - crossprod(sh)) / G
      if (!is.null(minVar)) diag(S) <- pmax(diag(S), minVar)
      chol <- .updateNormal(design, list(S), fixed)
      if (!length(free)) break
      Omega <- tcrossprod(.choleskyFactor(lay, c(chol, fixed)))
      beta[free] <- .glsFixedEffects(ref, s1, Omega, beta, free)
    }
    list(beta = beta, chol = chol)
  }

  full <- setNames(numeric(length(meta$pars_full_names)), meta$pars_full_names)
  full[names(start)] <- start
  im <- .innerModes(design, c(start, fixed), matrix(0, G, K, dimnames = dimnames(E)), fixed, meta)
  eta <- im$modes
  density <- .groupDensity(design, meta, fixed, cores)
  objPars <- attr(design$objfn, "parameters")
  dataTerm <- function(p, active, conditions) {
    active <- .intersect(active, objPars)
    if (!length(active) || !length(conditions))
      return(list(value = 0, gradient = .named(NULL, trustPars), hessian = .namedHessian(NULL, trustPars)))
    held <- .setdiff(.intersect(objPars, names(p)), active)
    out <- design$objfn(pars = p[active], fixed = c(fixed, p[held]), deriv = TRUE, conditions = conditions)
    list(value = out$value, gradient = .named(out$gradient, trustPars),
         hessian = .namedHessian(out$hessian, trustPars))
  }
  # the complete-data objective at the drawn phi: the data of the groups depend on the partially
  # referenced beta through phi only, those of conditions without a group directly, and the prior
  cm1Objective <- function(theta, ...) {
    p <- full; p[trustPars] <- theta
    p[as.vector(E)] <- as.vector(eta)
    if (!length(partial)) {
      out <- dataTerm(p, trustPars, names(design$data))
      return(objlist(value = out$value, gradient = out$gradient, hessian = out$hessian))
    }
    atPhi <- p; atPhi[partial] <- beta[partial]
    inGroups <- dataTerm(atPhi, .setdiff(trustPars, partial), grouped)
    outside <- dataTerm(p, trustPars, design$freeConditions)
    b <- beta; b[partial] <- theta[partial]
    r <- phi - shift(b)
    W <- chol2inv(t(.choleskyFactor(lay, c(full, fixed))))
    Dp <- lapply(D, function(m) m[, partial, drop = FALSE])
    gp <- -2 * Reduce(`+`, Map(function(d, g) crossprod(d, W %*% r[g, ]), Dp, seq_len(G)))
    Hp <- 2 * Reduce(`+`, lapply(Dp, function(d) crossprod(d, W %*% d)))
    grad <- inGroups$gradient + outside$gradient
    grad[partial] <- grad[partial] + drop(gp)
    hess <- inGroups$hessian + outside$hessian
    hess[partial, partial] <- hess[partial, partial] + Hp
    objlist(value = inGroups$value + outside$value + sum((r %*% W) * r), gradient = grad, hessian = hess)
  }

  # half the Hessian of the complete-data -2 log p at the averaged statistics, block diagonal
  # in the GLS beta, the Cholesky parameters and the trust parameters
  completeInformation <- function() {
    Ic <- matrix(0, length(start), length(start), dimnames = list(names(start), names(start)))
    L <- .choleskyFactor(lay, c(full, fixed))
    W <- chol2inv(t(L))
    if (length(free))
      Ic[free, free] <- Reduce(`+`, lapply(D, function(m) crossprod(m[, free, drop = FALSE], W %*% m[, free, drop = FALSE])))
    sh <- shift(beta)
    S <- s2 - (crossprod(sh, s1) + crossprod(s1, sh) - crossprod(sh)) / G
    chol <- .intersect(lay$cholPars, names(start))
    Ic[chol, chol] <- G / 2 * .choleskyDerivatives(lay, L, S)$hcc[chol, chol]
    if (length(trustPars)) Ic[trustPars, trustPars] <- cm1Objective(full[trustPars])$hessian / 2
    Ic
  }

  beta <- c(start, fixed)[refPars]
  step <- rep(0.4, G)
  modePhi <- im$modes + shift(beta)
  factors <- .proposalFactors(im$H, laplaceScale)
  s1 <- s2 <- phi <- NULL
  psi <- start
  history <- matrix(NA_real_, 0, length(start), dimnames = list(NULL, names(start)))
  targets <- history
  rows <- list()
  burnin <- TRUE; burnEnd <- NA_integer_; stationary <- FALSE
  info <- NULL; mcGap <- NA_real_; q <- NULL; model <- NULL
  for (k in seq_len(maxBurnin + maxAveraging)) {
    gamma <- if (burnin) 1 else (k - burnEnd)^-0.7
    L <- .choleskyFactor(lay, c(full, fixed))
    Oinv <- chol2inv(t(L)); logDet <- 2 * sum(log(diag(L)))
    target <- function(e) density(full, e, Oinv, logDet)
    # a state where the model fails takes any proposal it can evaluate
    value <- tryCatch(target(eta), error = function(e) rep(Inf, G))
    rw <- .metropolisGroups(target, eta, value, nMcmc, step = step)
    eta <- rw$eta
    if (nLaplace > 0L) {
      # where the modes cannot be refreshed the previous proposal stays, MH is exact either way
      im <- if (k > 1L && (k - 1L) %% laplaceEvery == 0L)
        tryCatch(.innerModes(design, c(full[names(start)], fixed), eta, fixed, meta), error = function(e) NULL)
      if (!is.null(im)) {
        modePhi <- im$modes + shift(beta)
        factors <- .proposalFactors(im$H, laplaceScale)
      }
      lap <- .metropolisGroups(target, eta, rw$value, nLaplace, modes = modePhi - shift(beta), chol = factors)
      eta <- lap$eta
    }
    if (burnin) step <- pmin(pmax(step * exp(0.3 * (rw$accept - 0.4)), 0.02), 5)

    phi <- eta + shift(beta)
    s1 <- if (is.null(s1)) phi else s1 + gamma * (phi - s1)
    s2 <- if (is.null(s2)) crossprod(phi) / G else s2 + gamma * (crossprod(phi) / G - s2)
    upd <- mStep(s1, s2, beta, if (burnin) 0.95 * diag(tcrossprod(L)))
    one <- if (!burnin) mStep(phi, crossprod(phi) / G, beta)
    beta <- upd$beta
    full[free] <- beta[free]
    held <- .intersect(names(upd$chol), names(start))
    full[held] <- upd$chol[held]
    eta <- phi - shift(beta)
    # the other parameters: a trust step while annealing, then the step on the averaged quadratic
    # model of the complete-data objective
    # where the complete-data objective fails at the draws, the trust parameters keep their values
    if (length(trustPars) && burnin) {
      cm <- tryCatch(.innerTrust(modifyList(cm1, list(objfun = cm1Objective, parinit = full[trustPars]))),
                     error = function(e) NULL)
      if (!is.null(cm)) full[trustPars] <- cm$argument
    }
    if (length(trustPars) && !burnin) {
      out <- tryCatch(cm1Objective(full[trustPars]), error = function(e) NULL)
      one$trust <- full[trustPars]
      if (!is.null(out)) {
        model <- .averageQuadratic(if (gamma < 1) model, out, full[trustPars], gamma)
        one$trust <- .quadraticStep(list(A = model$A, c = drop(model$A %*% full[trustPars]) - out$gradient),
                                    full[trustPars], Inf)
        full[trustPars] <- .quadraticStep(model, full[trustPars])
      }
    }
    if (length(trustPars)) {
      beta[partial] <- full[partial]
      eta <- phi - shift(beta)
    }
    prev <- psi
    psi <- full[names(start)]
    history <- rbind(history, psi)

    if (burnin) {
      stationary <- k >= window &&
        all(abs(.slopeStatistics(utils::tail(history, window))) < 2)
      if (stationary || k >= maxBurnin) {
        burnin <- FALSE; burnEnd <- k
        q <- .quadratureObjective(design, psi, fixed, quad)
        q$state$modes <- eta
        info <- .marginalInformation(q, psi, fixed, K + 1L)
        if (verbose) message(sprintf("em(saem) burn-in ends at iteration %d", k))
      }
    } else {
      # the parameters from this iteration's draws alone; the estimate is the mean of the
      # iterates after the first `window` of the averaging
      draw <- psi
      draw[free] <- one$beta[free]
      held <- .intersect(names(one$chol), names(start))
      draw[held] <- one$chol[held]
      if (length(trustPars)) draw[trustPars] <- one$trust
      targets <- rbind(targets, draw)
      n <- k - burnEnd
      Ic <- if (n >= 2L * window && n %% 10L == 0L) tryCatch(completeInformation(), error = function(e) NULL)
      if (!is.null(Ic))
        mcGap <- .monteCarloGap(info, Ic, .batchMeansCov(targets[-seq_len(window), , drop = FALSE]))
    }
    rows[[k]] <- data.frame(iter = k, phase = if (is.na(burnEnd) || k <= burnEnd) "burnin" else "averaging",
                            gamma = gamma, deltaPsi = max(abs(psi - prev)),
                            deltaStruct = if (length(structural)) max(abs(psi[structural] - prev[structural])) else 0,
                            accept = mean(rw$accept), mcGap = mcGap)
    if (verbose && k %% 50L == 0L)
      message(sprintf("em(saem) iteration %d, %s", k, rows[[k]]$phase))
    if (!is.na(mcGap) && mcGap <= mcTol) break
  }
  trace <- do.call(rbind, rows)
  if (!is.na(burnEnd) && nrow(history) > burnEnd + window)
    psi <- colMeans(history[(burnEnd + window + 1L):nrow(history), , drop = FALSE])
  q <- q %||% .quadratureObjective(design, psi, fixed, quad)
  q$state$modes <- eta
  q$rebuild(psi, K + 2L)
  converged <- stationary && !is.na(mcGap) && mcGap <= mcTol
  .quadratureResult(design, q, psi, fixed, nrow(trace), converged, trace)
}


## SAEM with Dirichlet-process random effects, normal ones beside them: Metropolis and Gibbs draws,
## Robbins-Monro statistics with closed-form M-steps, the structural parameters on the averaged
## quadratic model; -2 log L by importance sampling, the clusters the most frequent late partition.
.fitSaemDirichlet <- function(design, start, fixed, control = list(), verbose = FALSE) {
  lay <- design$normal
  sl <- .sparseLayout(design)
  if (any(sl$dist != "dirichletProcess"))
    stop("em: method = \"saem\" takes normal and Dirichlet-process random effects.", call. = FALSE)
  sc <- control$saem %||% list()
  maxBurnin <- sc$nBurnin %||% 200L
  maxAveraging <- sc$nEM %||% 300L
  nMcmc <- sc$nMcmc %||% 3L
  window <- sc$window %||% 25L
  cores <- getOption("dMod2.cores", 1L)
  E <- design$etaNames[, lay$stems, drop = FALSE]
  Dp <- design$etaNames[, sl$stem, drop = FALSE]
  G <- nrow(Dp); K <- ncol(E)
  hypers <- .intersect(.sparseHypers(sl), names(start))
  outer <- start[.setdiff(names(start), hypers)]
  meta <- .foceiMeta(design, c(outer, setNames(numeric(length(Dp)), as.vector(Dp))), fixed)
  .requireM3(design, meta, "SAEM")
  full <- setNames(numeric(length(meta$pars_full_names)), meta$pars_full_names)
  full[names(outer)] <- outer
  hyper <- start[hypers]
  # phi = eta + D_g beta held while the mu-referenced fixed effects move, those by GLS where their
  # symbol enters no condition outside the groups; the other structural parameters by trust
  ref <- .muReference(design)
  free <- .intersect(unlist(lapply(unname(ref), `[[`, "pars")), names(outer))
  refPars <- unlist(lapply(unname(ref), `[[`, "pars"))
  D <- lapply(design$groups, function(g) {
    m <- matrix(0, K, length(refPars), dimnames = list(lay$stems, refPars))
    for (e in ref) m[e$stem, e$pars] <- e$X[g, ]
    m
  })
  shift <- function(beta) matrix(vapply(D, function(m) drop(m %*% beta[refPars]), numeric(K)), G, K,
                                 byrow = TRUE, dimnames = dimnames(E))
  structural <- .setdiff(names(outer), c(lay$cholPars, free))
  density <- .groupDensity(design, meta, fixed, cores)
  objPars <- attr(design$objfn, "parameters")
  cm1 <- .trustControl(list(rinit = 1, rmax = 10, iterlim = 1, tolControl = list(ftol = 1e-6, mtol = 1e-6)),
                       sc$cm1Control, label = "cm1Control")
  # the complete-data -2 log L of the structural parameters at the drawn etas
  cm1Objective <- function(theta, ...) {
    p <- full; p[structural] <- theta
    p[as.vector(E)] <- as.vector(eta)
    active <- .intersect(structural, objPars)
    held <- .setdiff(.intersect(objPars, names(p)), active)
    out <- design$objfn(pars = p[active], fixed = c(fixed, p[held]), deriv = TRUE)
    objlist(value = out$value, gradient = .named(out$gradient, structural),
            hessian = .namedHessian(out$hessian, structural))
  }

  # every group starts in a separate cluster at the group mode under the slab alone, so the sampler
  # merges clusters rather than finding deviations by draws from the slab
  objs <- .conditionObjectives(design)
  prior <- if (K) .normalPrior(design)
  start0 <- matrix(0, G, ncol(design$etaNames), dimnames = dimnames(design$etaNames))
  mode0 <- .sparseModes(design, objs, prior, c(outer, fixed), .stemHypers(sl, c(hyper, fixed)), start0,
                        list(iterlim = 300L))
  state <- lapply(setNames(nm = sl$stem), function(st) list(z = seq_len(G), v = unname(mode0[, st])))
  for (st in sl$stem) full[Dp[, st]] <- mode0[, st]
  eta <- mode0[, lay$stems, drop = FALSE]
  step <- rep(0.4, G)
  valueStep <- setNames(rep(0.3, nrow(sl)), sl$stem)
  beta <- c(outer, fixed)[refPars]
  s1 <- s2 <- NULL
  stats <- NULL
  blocks <- Filter(function(b) b$dist == "dirichletProcess", design$blocks)
  psiOf <- function() c(full[names(outer)], hyper)[names(start)]
  history <- matrix(NA_real_, 0, length(start), dimnames = list(NULL, names(start)))
  draws <- list()
  rows <- list()
  burnin <- TRUE; burnEnd <- NA_integer_; stationary <- FALSE; model <- NULL
  for (k in seq_len(maxBurnin + maxAveraging)) {
    gamma <- if (burnin) 1 else (k - burnEnd)^-0.7
    L <- if (K) .choleskyFactor(lay, c(full, fixed))
    Oinv <- if (K) chol2inv(t(L)) else matrix(0, 0, 0)
    logDet <- if (K) 2 * sum(log(diag(L))) else 0
    hyp <- .stemHypers(sl, c(hyper, fixed))
    accept <- NA_real_
    if (K) {
      rw <- .metropolisGroups(function(e) density(full, e, Oinv, logDet), eta,
                              tryCatch(density(full, eta, Oinv, logDet), error = function(e) rep(Inf, nrow(eta))),
                              nMcmc, step = step)
      eta <- rw$eta
      accept <- mean(rw$accept)
      if (burnin) step <- pmin(pmax(step * exp(0.3 * (rw$accept - 0.4)), 0.02), 5)
    }
    for (st in sl$stem) {
      at <- Dp[, st]
      dataAt <- function(v) { f <- full; f[at] <- v; density(f, eta, Oinv, logDet) }
      state[[st]] <- .dpGibbs(state[[st]], dataAt, hyp$concentration[[st]], hyp$tau[[st]])
      mv <- .dpValues(state[[st]], dataAt, hyp$tau[[st]], valueStep[[st]])
      state[[st]] <- mv$s
      if (burnin && !is.na(mv$accept))
        valueStep[[st]] <- min(max(valueStep[[st]] * exp(0.3 * (mv$accept - 0.4)), 0.01), 5)
      full[at] <- c(0, state[[st]]$v)[state[[st]]$z + 1L]
    }

    # A-step and M-steps: Omega and the GLS fixed effects on the averages of phi, then per block
    # the slab sd and the concentration
    if (K) {
      phi <- eta + shift(beta)
      s1 <- if (is.null(s1)) phi else s1 + gamma * (phi - s1)
      s2 <- if (is.null(s2)) crossprod(phi) / G else s2 + gamma * (crossprod(phi) / G - s2)
      for (it in 1:3) {
        sh <- shift(beta)
        S <- s2 - (crossprod(sh, s1) + crossprod(s1, sh) - crossprod(sh)) / G
        if (burnin) diag(S) <- pmax(diag(S), 0.95 * diag(tcrossprod(L)))
        chol <- .updateNormal(design, list(S), fixed)
        held <- .intersect(names(chol), names(outer))
        full[held] <- chol[held]
        if (!length(free)) break
        beta[free] <- .glsFixedEffects(ref, s1, tcrossprod(.choleskyFactor(lay, c(full, fixed))), beta, free)
      }
      full[free] <- beta[free]
      eta <- phi - shift(beta)
    }
    cur <- vapply(blocks, function(b) c(free = sum(vapply(b$stems, function(st) length(state[[st]]$v), 0)),
                                        slab2 = sum(vapply(b$stems, function(st) sum(state[[st]]$v^2), 0))),
                  numeric(2))
    stats <- if (is.null(stats)) cur else stats + gamma * (cur - stats)
    # while annealing the slab and the concentration change by a factor of at most 0.95 to 1.25
    # per iteration, so one sweep cannot drive them to a limit
    bound <- function(new, old) if (burnin) min(max(new, 0.95 * old), 1.25 * old) else new
    for (j in seq_along(blocks)) {
      b <- blocks[[j]]
      tn <- b$hyper[["tau"]]; cn <- b$hyper[["concentration"]]
      if (tn %in% hypers && stats["free", j] > 1e-8)
        hyper[[tn]] <- bound(sqrt(stats["slab2", j] / stats["free", j]), hyper[[tn]])
      if (cn %in% hypers)
        hyper[[cn]] <- bound(.ewensConcentration(stats["free", j], rep(G, length(b$stems))), hyper[[cn]])
    }
    # where the complete-data objective fails at the draws, the structural parameters keep their values
    if (length(structural) && burnin) {
      cm <- tryCatch(.innerTrust(modifyList(cm1, list(objfun = cm1Objective, parinit = full[structural]))),
                     error = function(e) NULL)
      if (!is.null(cm)) full[structural] <- cm$argument
    }
    out <- if (length(structural) && !burnin) tryCatch(cm1Objective(full[structural]), error = function(e) NULL)
    if (!is.null(out)) {
      model <- .averageQuadratic(if (gamma < 1) model, out, full[structural], gamma)
      full[structural] <- .quadraticStep(model, full[structural])
    }
    prev <- if (nrow(history)) history[nrow(history), ] else psiOf()
    psi <- psiOf()
    history <- rbind(history, psi)
    if (!burnin) draws[[length(draws) + 1L]] <- lapply(state, .canonicalClusters)

    if (burnin) {
      stationary <- k >= window && all(abs(.slopeStatistics(utils::tail(history, window))) < 2)
      if (stationary || k >= maxBurnin) {
        burnin <- FALSE; burnEnd <- k
        if (verbose) message(sprintf("em(saem) burn-in ends at iteration %d", k))
      }
    }
    rows[[k]] <- data.frame(iter = k, phase = if (is.na(burnEnd) || k <= burnEnd) "burnin" else "averaging",
                            gamma = gamma, deltaPsi = if (length(psi)) max(abs(psi - prev)) else 0,
                            deltaStruct = if (length(structural)) max(abs(psi[structural] - prev[structural])) else 0,
                            accept = accept, clusters = sum(vapply(state, function(s) length(s$v), 0)))
    if (verbose && k %% 50L == 0L) message(sprintf("em(saem) iteration %d, %s", k, rows[[k]]$phase))
    if (!burnin && k - burnEnd >= maxAveraging) break
  }
  trace <- do.call(rbind, rows)
  late <- if (!is.na(burnEnd) && nrow(history) > burnEnd + window) (burnEnd + window + 1L):nrow(history)
          else seq_len(nrow(history))
  psi <- colMeans(history[late, , drop = FALSE])
  full[names(outer)] <- psi[names(outer)]
  hyper[] <- psi[hypers]

  # the most frequent partition of every stem in the late draws, the cluster values averaged over
  # the draws with that partition, and the share of draws at the population value per group
  keep <- utils::tail(draws, length(late))
  modes <- matrix(0, G, ncol(design$etaNames), dimnames = dimnames(design$etaNames))
  modes[, lay$stems] <- eta
  clusters <- spike <- matrix(NA_real_, G, ncol(design$etaNames), dimnames = dimnames(design$etaNames))
  for (st in sl$stem) {
    keys <- vapply(keep, function(d) paste(d[[st]]$z, collapse = " "), "")
    top <- names(which.max(table(keys)))
    with <- keep[keys == top]
    z <- with[[1]][[st]]$z
    values <- colMeans(do.call(rbind, lapply(with, function(d) d[[st]]$v)) %||% matrix(0, 1, 0))
    modes[, st] <- c(0, values)[z + 1L]
    clusters[, st] <- z
    spike[, st] <- colMeans(do.call(rbind, lapply(keep, function(d) d[[st]]$z == 0L)))
  }
  r <- .sparseMarginal(design, objs, prior, c(full[names(outer)], fixed), .stemHypers(sl, c(hyper, fixed)),
                       as.character(.intersect(names(outer), names(start))), character(0), fixed, modes, grad = TRUE,
                       control = list(integral = "importance", particles = sc$particles %||% 500L, iterlim = 300L))
  modes[, lay$stems] <- r$modes[, lay$stems]
  list(argument = psi, value = r$value, valueSE = r$valueSE, gradient = NULL, hessian = NULL,
       etaModes = modes, etaCov = setNames(r$cov, design$groups), clusters = clusters, spike = spike,
       converged = stationary, iterations = nrow(trace), stageTrace = trace)
}
