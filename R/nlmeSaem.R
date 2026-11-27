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

## -2 log p(y_g, eta_g) of every group at the rows of `eta`, as the FOCEI kernel per group, all
## conditions in one evaluation without derivatives; with `full` a list of parameter vectors, one
## column per vector.
.groupDensity <- function(design, meta, fixed, cores = 1L) {
  errfn <- .errorModel(design)
  node <- .fnNode(errfn)
  blocks <- unlist(lapply(meta$fast_meta, `[[`, "blocks"), recursive = FALSE)
  group <- rep(seq_along(meta$fast_meta), vapply(meta$fast_meta, function(m) length(m$blocks), 1L))
  conditions <- vapply(blocks, function(b) as.character(b$condition), "")
  times <- lapply(blocks, `[[`, "times")
  etaNames <- lapply(meta$fast_meta, `[[`, "eta_names")
  G <- length(etaNames)
  dataTerm <- function(parsList, n) {
    d <- .chainLocal(design$chain, rep(conditions, n), parsList, rep(list(fixed), length(parsList)), FALSE, FALSE,
                     cores, design$bloq, design$df)
    if (is.null(d)) return(NULL)
    v <- vapply(d, function(x) as.numeric(x$value), 0)
    matrix(rowsum(v, rep(group, n) + G * rep(seq_len(n) - 1L, each = length(group)), reorder = TRUE), G, n)
  }
  one <- function(full, eta) {
    parsList <- lapply(seq_along(blocks), function(j) {
      p <- full; p[etaNames[[group[j]]]] <- eta[group[j], ]; p
    })
    preds <- .predictMany(design$prdfn, times = times, parsList = parsList, conditions = conditions,
                          fixed = fixed, deriv = FALSE, cores = cores)
    inner <- lapply(preds, function(q) attr(q, "parameters"))
    errs <- as.list(if (!is.null(node))
      .evalNode(node, .bundle(conds = conditions, out = as.list(preds), pars = inner,
                              fixed = vector("list", length(preds))), FALSE, FALSE, NULL, cores)
      else lapply(seq_along(preds), function(j)
        errfn(out = preds[[j]], pars = inner[[j]], conditions = conditions[j])[[1L]]))
    value <- numeric(G)
    for (j in seq_along(blocks)) {
      b <- blocks[[j]]
      err <- errs[[j]] %||% .dataSigmaErrfn(design$data)(out = preds[[j]], conditions = conditions[j])[[1L]]
      p <- preds[[j]][cbind(b$t_idx_in_pred, b$o_idx_in_pred)]
      s <- err[cbind(b$t_idx_in_err, b$o_idx_in_err)]
      recorded <- !is.na(b$sigma_data)
      s[recorded] <- b$sigma_data[recorded]
      value[group[j]] <- value[group[j]] + .rowLoss(p, s, b$y_data, seq_along(p) > b$n_aloq, design$df)
    }
    value
  }
  function(full, eta, Oinv, logDet) {
    many <- is.list(full)
    sets <- if (many) full else list(full)
    parsList <- unlist(lapply(sets, function(f) lapply(seq_along(blocks), function(j) {
      p <- f; p[etaNames[[group[j]]]] <- eta[group[j], ]; p
    })), recursive = FALSE)
    value <- dataTerm(parsList, length(sets)) %||% vapply(sets, one, numeric(G), eta = eta)
    value <- matrix(value, G) + rowSums((eta %*% Oinv) * eta) + logDet
    if (many) value else value[, 1]
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
  saem <- .saemSettings(control)
  cores <- getOption("dMod2.cores", 1L)
  E <- design$etaNames[, lay$stems, drop = FALSE]
  G <- nrow(E); K <- ncol(E)
  meta <- .foceiMeta(design, start, fixed)
  .requireM3(design, meta, "SAEM")
  quad <- modifyList(control, control$quadrature %||% list())

  # phi = eta + D_g beta held while beta moves, for the random intercepts with X constant in the
  # groups: GLS for those of .muReference(), the trust step for those whose symbol also enters
  # conditions without a group
  grouped <- unlist(design$conditions, use.names = FALSE)
  groupedDesign <- design
  groupedDesign$data <- design$data[grouped]
  ref <- .muReference(groupedDesign)
  gls <- names(.muReference(design))
  refs <- .referenceMatrices(design, ref, E)
  refPars <- refs$refPars
  D <- refs$D
  shift <- refs$shift
  free <- .intersect(unlist(lapply(unname(ref[gls]), `[[`, "pars")), names(start))
  structural <- .setdiff(names(start), lay$cholPars)
  trustPars <- .setdiff(structural, free)
  partial <- .intersect(trustPars, refPars)
  # Omega and the GLS beta that maximize the expected prior given the averages s1 of phi_g and
  # s2 of phi_g phi_g' / G; the diagonal of S at least `minVar` while annealing
  mStep <- function(s1, s2, beta, minVar = NULL) {
    for (it in 1:3) {
      sh <- shift(beta)
      S <- s2 - (crossprod(sh, s1) + crossprod(s1, sh) - crossprod(sh)) / G
      if (!is.null(minVar)) diag(S) <- pmax(diag(S), minVar)
      cholPars <- .updateNormal(design, list(S), fixed)
      if (!length(free)) break
      Omega <- tcrossprod(.choleskyFactor(lay, c(cholPars, fixed)))
      beta[free] <- .glsFixedEffects(ref, s1, Omega, beta, free)
    }
    list(beta = beta, chol = cholPars)
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
    cp <- .intersect(lay$cholPars, names(start))
    Ic[cp, cp] <- G / 2 * .choleskyDerivatives(lay, L, S)$hcc[cp, cp]
    if (length(trustPars)) Ic[trustPars, trustPars] <- cm1Objective(full[trustPars])$hessian / 2
    Ic
  }

  beta <- c(start, fixed)[refPars]
  step <- rep(0.4, G)
  modePhi <- im$modes + shift(beta)
  factors <- .proposalFactors(im$H, saem$laplaceScale)
  s1 <- s2 <- phi <- NULL
  psi <- start
  history <- matrix(NA_real_, 0, length(start), dimnames = list(NULL, names(start)))
  targets <- history
  rows <- list()
  burnin <- TRUE; burnEnd <- NA_integer_; stationary <- FALSE
  info <- NULL; mcGap <- NA_real_; q <- NULL; model <- NULL
  for (k in seq_len(saem$nBurnin + saem$nEM)) {
    gamma <- if (burnin) 1 else (k - burnEnd)^-0.7
    L <- .choleskyFactor(lay, c(full, fixed))
    Oinv <- chol2inv(t(L)); logDet <- 2 * sum(log(diag(L)))
    target <- function(e) density(full, e, Oinv, logDet)
    # a state where the model fails takes any proposal it can evaluate
    value <- tryCatch(target(eta), error = function(e) rep(Inf, G))
    rw <- .metropolisGroups(target, eta, value, saem$nMcmc, step = step)
    eta <- rw$eta
    if (saem$nLaplace > 0L) {
      # where the modes cannot be refreshed the previous proposal stays, MH is exact either way
      im <- if (k > 1L && (k - 1L) %% saem$laplaceEvery == 0L)
        tryCatch(.innerModes(design, c(full[names(start)], fixed), eta, fixed, meta), error = function(e) NULL)
      if (!is.null(im)) {
        modePhi <- im$modes + shift(beta)
        factors <- .proposalFactors(im$H, saem$laplaceScale)
      }
      lap <- .metropolisGroups(target, eta, rw$value, saem$nLaplace, modes = modePhi - shift(beta), chol = factors)
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
      cm <- tryCatch(.innerTrust(modifyList(saem$cm1, list(objfun = cm1Objective, parinit = full[trustPars]))),
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
      stationary <- k >= saem$window && all(abs(.slopeStatistics(utils::tail(history, saem$window))) < 2)
      if (stationary || k >= saem$nBurnin) {
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
      Ic <- if (n >= 2L * saem$window && n %% 10L == 0L) tryCatch(completeInformation(), error = function(e) NULL)
      if (!is.null(Ic))
        mcGap <- .monteCarloGap(info, Ic, .batchMeansCov(targets[-seq_len(saem$window), , drop = FALSE]))
    }
    rows[[k]] <- data.frame(iter = k, phase = if (is.na(burnEnd) || k <= burnEnd) "burnin" else "averaging",
                            gamma = gamma, deltaPsi = max(abs(psi - prev)),
                            deltaStruct = if (length(structural)) max(abs(psi[structural] - prev[structural])) else 0,
                            accept = mean(rw$accept), mcGap = mcGap)
    if (verbose && k %% 50L == 0L)
      message(sprintf("em(saem) iteration %d, %s", k, rows[[k]]$phase))
    if (!is.na(mcGap) && mcGap <= saem$mcTol) break
  }
  trace <- do.call(rbind, rows)
  if (!is.na(burnEnd) && nrow(history) > burnEnd + saem$window)
    psi <- colMeans(history[(burnEnd + saem$window + 1L):nrow(history), , drop = FALSE])
  q <- q %||% .quadratureObjective(design, psi, fixed, quad)
  q$state$modes <- eta
  q$rebuild(psi, K + 2L)
  converged <- stationary && !is.na(mcGap) && mcGap <= saem$mcTol
  .quadratureResult(design, q, psi, fixed, nrow(trace), converged, trace)
}


## ---- spike-and-slab random effects ---------------------------------------------------------

## The quadratic model of the data -2 log L of every group in the spike-and-slab coordinates `at`
## (G x J names) at the values `ref`: exact gradient and the curvature of the objective at `pars`
## with `at` set to `ref`.
.slabModel <- function(design, objs, pars, at, ref) {
  G <- nrow(at)
  out <- .groupData(design, objs, lapply(seq_len(G), function(g) { p <- pars[[g]]; p[at[g, ]] <- ref[g, ]; p }))
  grad <- t(vapply(seq_len(G), function(g) .named(out[[g]]$gradient, at[g, ]), numeric(ncol(at))))
  if (ncol(at) == 1L) grad <- t(grad)
  if (any(!is.finite(grad))) return(NULL)
  H <- lapply(seq_len(G), function(g) .namedHessian(out[[g]]$hessian, at[g, ]))
  if (any(!is.finite(unlist(H)))) return(NULL)
  # the curvature with negative eigenvalues set to 0, so every subset has a Gaussian
  H <- lapply(H, function(h) {
    e <- eigen((h + t(h)) / 2, symmetric = TRUE)
    structure(e$vectors %*% (t(e$vectors) * pmax(e$values, 0)), dimnames = dimnames(h))
  })
  list(ref = ref, grad = grad, H = H)
}

## The Gaussian of the slab values of group `g` in the subset `s`, the other coordinates at 0, on
## the quadratic model plus the slab prior of mean `mu`: the mean (0 off `s`) and the upper factor R
## of the precision P, the covariance 2 `scale` P^-1.
.slabGaussian <- function(model, g, s, tau, scale, mu = 0 * tau) {
  mean <- numeric(length(s))
  if (!any(s)) return(list(mean = mean, R = NULL, s = s, scale = scale))
  H <- model$H[[g]]
  R <- chol(H[s, s, drop = FALSE] + diag(2 / tau[s]^2, sum(s)))
  rhs <- drop(H %*% model$ref[g, ])[s] - model$grad[g, s] + 2 * mu[s] / tau[s]^2
  mean[s] <- backsolve(R, forwardsolve(t(R), rhs))
  list(mean = mean, R = R, s = s, scale = scale)
}

## Probabilists' Gauss-Hermite nodes and weights of order 3 or 5: sum w f(t) = E f(Z), Z ~ N(0, 1),
## exact for polynomials of degree 2 n - 1.
.hermiteNodes <- function(n) {
  if (!n %in% c(3, 5)) stop("em: `control$saem$nodes` is 3 or 5.", call. = FALSE)
  gh <- .gaussHermite(n)
  list(t = sqrt(2) * gh$x, w = gh$w / sqrt(pi))
}

## The collapsed move on indicator `j[g]` of every group g (NA: none) with spike `spike` and slab N(`mu`, tau^2):
## the slab value integrated out by Gauss-Hermite, independence Metropolis-Hastings on the exact posterior
## (methods-nlme, SAEM); for a spike at -Inf `override` holds the proposal mean and sd per group.
.slabCollapsed <- function(V, Z, value, dataAt, model, prob, tau, j, scale, gh, spike = 0 * tau, mu = 0 * tau,
                           override = NULL, beta = NULL) {
  G <- nrow(V); J <- ncol(V)
  go <- !is.na(j)
  base <- V; dir <- matrix(0, G, J); m <- sd <- rep(NA_real_, G)
  # a group whose local model has no Gaussian along the line does not move
  for (g in which(go)) tryCatch({
    i <- j[g]
    o <- which(Z[g, ] & seq_len(J) != i)
    H <- model$H[[g]]
    S <- c(i, o)
    P <- H[S, S, drop = FALSE] + diag(2 / tau[S]^2, length(S))
    B <- if (length(o) && is.finite(spike[i])) -solve(P[-1, -1, drop = FALSE], P[-1, 1]) else numeric(0)
    xi <- if (Z[g, i]) V[g, i] else spike[i]
    if (length(B)) base[g, o] <- V[g, o] - B * xi
    base[g, i] <- 0
    dir[g, i] <- 1
    if (length(B)) dir[g, o] <- B
    if (!is.null(override) && !is.na(override$m[g])) { m[g] <- override$m[g]; sd[g] <- override$sd[g]; next }
    # coordinates at a spike of -Inf leave the quadratic model at the reference
    dv <- base[g, ] - model$ref[g, ]; dv[!is.finite(dv)] <- 0
    dp <- base[g, ] - mu; dp[!is.finite(dp)] <- 0
    grad <- model$grad[g, ] + drop(H %*% dv) + 2 * dp / tau^2
    s <- P[1, 1] + if (length(B)) sum(P[1, -1] * B) else 0
    if (!(s > 0)) stop("no curvature")
    m[g] <- -(grad[i] + if (length(B)) sum(B * grad[o]) else 0) / s
    sd[g] <- sqrt(2 / s)
  }, error = function(e) go[g] <<- FALSE)
  lineAt <- function(x) base + ifelse(dir == 0, 0, x * dir)
  out <- list(V = V, Z = Z, value = value, birth = rep(NA, G), ok = rep(FALSE, G),
              pinc = rep(NA_real_, G), slab2 = rep(NA_real_, G))
  if (!any(go)) return(out)
  ig <- cbind(seq_len(G), ifelse(go, j, 1L))
  # log p(y, eta) of a group on the line at x from the data term d, with the prior of every
  # coordinate of the group
  lp1 <- matrix(log(prob), G, J, byrow = TRUE); lp0 <- matrix(log1p(-prob), G, J, byrow = TRUE)
  # the coordinates under a beta prior on prob enter by the count of their block instead
  if (!is.null(beta)) lp1[, !is.na(beta$block)] <- lp0[, !is.na(beta$block)] <- 0
  sdSlab <- matrix(tau, G, J, byrow = TRUE)
  muSlab <- matrix(mu, G, J, byrow = TRUE)
  xOff <- ifelse(go, spike[ig[, 2]], 0)
  logF <- function(x, d, slab) {
    z <- Z; z[ig] <- slab
    v <- lineAt(ifelse(go, x, 0))
    out <- rowSums(ifelse(z, lp1 + stats::dnorm(v, muSlab, sdSlab, log = TRUE), lp0)) - d / 2
    if (!is.null(beta)) out <- out + .betaCountPrior(z, beta)
    ifelse(go, out, NA_real_)
  }
  # the data term on the line at the columns of `x`, all in one evaluation; where that fails each
  # column alone, a failed column Inf
  at <- function(x) {
    x <- matrix(x, G)
    x[!go, ] <- 0
    one <- function(k) tryCatch(dataAt(lineAt(x[, k])), error = function(e) rep(Inf, G))
    d <- tryCatch(matrix(dataAt(lapply(seq_len(ncol(x)), function(k) lineAt(x[, k]))), G),
                  error = function(e) vapply(seq_len(ncol(x)), one, numeric(G)))
    matrix(ifelse(is.finite(d), d, Inf), G)
  }
  inSlab <- go & Z[ig] == 1
  nodes <- outer(m, rep(1, length(gh$t))) + outer(sd, gh$t)
  d <- at(cbind(xOff, nodes))
  d0 <- ifelse(inSlab, d[, 1], value)
  lf0 <- logF(xOff, d0, FALSE)
  lw <- vapply(seq_along(gh$t), function(q) {
    x <- nodes[, q]
    log(gh$w[q]) + logF(x, d[, q + 1L], TRUE) - stats::dnorm(x, m, sd, log = TRUE)
  }, numeric(G))
  sp <- sd * sqrt(scale)
  lw <- matrix(lw, G)
  top <- apply(lw, 1, max)
  logI <- top + log(rowSums(exp(lw - top)))
  logI[!is.finite(top)] <- -Inf
  m2 <- rowSums(exp(lw - top) * (nodes - mu[ig[, 2]])^2) / rowSums(exp(lw - top))
  pinc <- ifelse(is.finite(lf0), 1 / (1 + exp(lf0 - logI)), 1)
  pinc[!go] <- NA
  q1 <- 0.01 + 0.98 * pinc
  logW <- function(slab, x, lf) ifelse(slab, lf - log(q1) - stats::dnorm(x, m, sp, log = TRUE), lf - log1p(-q1))
  cur <- logW(inSlab, V[ig], ifelse(inSlab, logF(V[ig], value, TRUE), lf0))
  slabNew <- stats::runif(G) < q1
  xNew <- ifelse(slabNew, m + sp * stats::rnorm(G), xOff)
  dNew <- d0
  if (any(go & slabNew)) dNew[go & slabNew] <- at(xNew)[go & slabNew, 1]
  new <- logW(slabNew, xNew, ifelse(slabNew, logF(xNew, dNew, TRUE), lf0))
  ok <- go & is.finite(new) & (!is.finite(cur) | log(stats::runif(G)) < new - cur)
  newV <- lineAt(ifelse(go, xNew, 0))
  for (g in which(ok)) {
    out$V[g, ] <- newV[g, ]
    out$Z[g, j[g]] <- slabNew[g]
    out$value[g] <- dNew[g]
  }
  # a proposed change of the indicator: TRUE a birth, FALSE a death
  out$birth <- ifelse(go & slabNew != inSlab, slabNew, NA)
  out$ok <- ok
  out$pinc <- pinc
  out$slab2 <- ifelse(go, pinc * ifelse(is.finite(m2), m2, 0), NA)
  out
}

## Per group the log prior of the indicators `z` (groups x stems) under the beta priors of `beta`: per block
## log B(a + n, b + K - n) / B(a, b) with n of the K coordinates of the block that reach data in the slab.
.betaCountPrior <- function(z, beta) {
  out <- numeric(nrow(z))
  for (bl in unique(stats::na.omit(beta$block))) {
    cols <- which(beta$block == bl)
    n <- rowSums(z[, cols, drop = FALSE] & beta$reach[, cols, drop = FALSE])
    K <- rowSums(beta$reach[, cols, drop = FALSE])
    out <- out + lbeta(beta$a[[bl]] + n, beta$b[[bl]] + K - n) - lbeta(beta$a[[bl]], beta$b[[bl]])
  }
  out
}

## A random-walk Metropolis step on all slab values of every group where the indicator is set, of
## size `step` per group times `sd` per coordinate; the acceptance per group.
.slabWalk <- function(V, Z, value, dataAt, tau, step, sd, mu = 0 * tau) {
  G <- nrow(V)
  prop <- ifelse(Z, V + step * sd * matrix(stats::rnorm(length(V)), G), V)
  pv <- tryCatch(dataAt(prop), error = function(e) rep(Inf, G))
  t2 <- matrix(rep(tau^2, each = G), G)
  mm <- matrix(rep(mu, each = G), G)
  ok <- rowSums(Z) > 0 & is.finite(pv) &
    log(stats::runif(G)) < (value - pv) / 2 - rowSums(ifelse(Z, ((prop - mm)^2 - (V - mm)^2) / (2 * t2), 0))
  V[ok, ] <- prop[ok, ]
  value[ok] <- pv[ok]
  list(V = V, value = value, accept = ok)
}


## SAEM with spike-and-slab and Dirichlet-process random effects, normal ones beside them; the
## estimate is the mean of the late iterates, the selected effects and clusters come from the late
## draws.
.fitSaemSparse <- function(design, start, fixed, control = list(), verbose = FALSE) {
  lay <- design$normal
  sl <- .sparseLayout(design)
  saem <- .saemSettings(control)
  gh <- .hermiteNodes(saem$nodes)
  cores <- getOption("dMod2.cores", 1L)
  E <- design$etaNames[, lay$stems, drop = FALSE]
  Dp <- design$etaNames[, sl$stem, drop = FALSE]
  G <- nrow(Dp); K <- ncol(E)
  dpStems <- sl$stem[sl$dist == "dirichletProcess"]
  slabStems <- sl$stem[sl$dist == "spikeSlab"]
  # the spike value and the slab mean of every spike-and-slab stem; other than 0 the final -2 log L is the
  # sum over subsets of em(method = "focei")
  blockOf <- function(st) Filter(function(b) st %in% b$stems, design$blocks)[[1]]
  spikeS <- vapply(slabStems, function(st) blockOf(st)$spike %||% 0, 0)
  meanS <- vapply(slabStems, function(st) blockOf(st)$mean %||% 0, 0)
  # beta priors on prob: per slab stem the block, and the shapes by block
  betaBlock <- vapply(slabStems, function(st) {
    b <- blockOf(st)
    if (is.null(b$probPrior)) NA_integer_ else match(TRUE, vapply(design$blocks, identical, TRUE, b))
  }, 0L)
  betaS <- if (any(!is.na(betaBlock))) list(block = betaBlock,
    a = lapply(design$blocks, function(b) b$probPrior$a), b = lapply(design$blocks, function(b) b$probPrior$b))
  hyperprior <- any(vapply(slabStems, function(st) !is.null(blockOf(st)$probPrior) || !is.null(blockOf(st)$tauPrior),
                           TRUE))
  general <- any(spikeS != 0 | meanS != 0) || hyperprior
  if (general && length(dpStems))
    stop("em: spikeSlab() with `mean` or `spike` other than 0 is fitted by method = \"saem\" beside no ",
         "Dirichlet-process random effects.", call. = FALSE)
  # -2 log L by importance sampling over the partitions where the exact sum beside spike-and-slab
  # random effects is too large or `integral` asks for it, else by the exact sum or EC
  importance <- length(dpStems) > 0L && (!length(slabStems) || G > 14L || identical(saem$integral, "importance"))
  hypers <- .intersect(.sparseHypers(sl), names(start))
  off <- .intersect(attr(.extraObjective(design), "parameters"), c(.sparseHypers(sl), as.vector(design$etaNames)))
  if (length(off))
    stop("em: with spike-and-slab or Dirichlet-process random effects the terms beside the data ",
         "take the fixed effects only; ", paste(off, collapse = ", "), " are hyperparameters or random effects.",
         call. = FALSE)
  outer <- start[.setdiff(names(start), hypers)]
  meta <- .foceiMeta(design, c(outer, setNames(numeric(length(Dp)), as.vector(Dp))), fixed)
  .requireM3(design, meta, "SAEM")
  chain <- list(full = setNames(numeric(length(meta$pars_full_names)), meta$pars_full_names), hyper = start[hypers])
  chain$full[names(outer)] <- outer
  # phi = eta + D_g beta held while the mu-referenced fixed effects move, those by GLS where their
  # symbol enters no condition outside the groups; the other structural parameters by trust
  ref <- .muReference(design)
  free <- .intersect(unlist(lapply(unname(ref), `[[`, "pars")), names(outer))
  refs <- .referenceMatrices(design, ref, E)
  structural <- .setdiff(names(outer), c(lay$cholPars, free))
  density <- .groupDensity(design, meta, fixed, cores)
  objPars <- attr(design$objfn, "parameters")
  # the complete-data -2 log L of the structural parameters at the parameters `full` and the drawn etas `eta`
  completeData <- function(full, eta) function(theta, ...) {
    p <- full; p[structural] <- theta
    p[as.vector(E)] <- as.vector(eta)
    active <- .intersect(structural, objPars)
    held <- .setdiff(.intersect(objPars, names(p)), active)
    out <- design$objfn(pars = p[active], fixed = c(fixed, p[held]), deriv = TRUE)
    objlist(value = out$value, gradient = .named(out$gradient, structural),
            hessian = .namedHessian(out$hessian, structural))
  }

  # every group starts in the slab and in a separate cluster at the group mode under the slab
  # alone, so the sampler removes deviations rather than finding them by draws from the slab; a
  # fit as init gives the start instead: the selected values, 0 in the spike, clusters by value
  objs <- .conditionObjectives(design)
  prior <- if (K) .normalPrior(design)
  En <- design$etaNames
  fitsDesign <- function(m) is.matrix(m) && all(rownames(En) %in% rownames(m)) && all(colnames(En) %in% colnames(m))
  selected0 <- if (fitsDesign(control$selected)) control$selected[rownames(En), colnames(En), drop = FALSE]
  start0 <- if (fitsDesign(control$modes)) control$modes[rownames(En), colnames(En), drop = FALSE]
            else matrix(0, G, ncol(En), dimnames = dimnames(En))
  mode0 <- .sparseModes(design, c(outer, fixed), .stemHypers(sl, c(chain$hyper, fixed)), start0, list(iterlim = 300L))
  if (!is.null(selected0)) mode0[, sl$stem] <- selected0[, sl$stem]
  chain$state <- lapply(setNames(nm = dpStems), function(st) {
    if (is.null(selected0)) return(list(z = seq_len(G), v = unname(mode0[, st])))
    v <- unique(mode0[mode0[, st] != 0, st])
    list(z = match(mode0[, st], v, nomatch = 0L), v = unname(v))
  })
  for (st in sl$stem) chain$full[Dp[, st]] <- mode0[, st]
  Ds <- Dp[, slabStems, drop = FALSE]
  reach <- .reachOf(design)[, slabStems, drop = FALSE]
  # the moves of all groups at once keep the posterior where a beta prior couples no two groups
  for (bl in unique(stats::na.omit(betaBlock)))
    if (sum(rowSums(reach[, which(betaBlock == bl), drop = FALSE]) > 0) > 1L)
      stop("em: method = \"saem\" takes betaPrior() for spike-and-slab random effects that reach the data of one ",
           "group.", call. = FALSE)
  spikeM <- matrix(spikeS, G, length(slabStems), byrow = TRUE)
  # a slab start without a finite slab-alone mode is the slab mean
  chain$full[Ds] <- ifelse(is.finite(chain$full[Ds]), chain$full[Ds], rep(meanS, each = G))
  chain$full[Ds[!reach]] <- spikeM[!reach]
  chain$incl <- if (is.null(selected0)) reach else reach & matrix(chain$full[Ds] != spikeM, G)
  chain$slabRef <- matrix(chain$full[Ds], G, dimnames = dimnames(Ds))
  # an estimated slab starts at least as wide as the largest of those modes: in a narrow slab a coordinate
  # without a deviation costs little, and tau and prob drift together to a slab of width 0 including all
  for (b in Filter(function(b) b$dist == "spikeSlab", design$blocks)) {
    tn <- unname(b$hyper["tau"])
    if (tn %in% hypers && is.null(selected0))
      chain$hyper[[tn]] <- max(chain$hyper[[tn]], abs(chain$slabRef[, match(b$stems, slabStems)] -
                                                       rep(meanS[match(b$stems, slabStems)], each = G)))
  }
  chain$slabStep <- rep(0.5, G)
  chain$eta <- mode0[, lay$stems, drop = FALSE]
  chain$step <- rep(0.4, G)
  chain$valueStep <- setNames(rep(0.3, length(dpStems)), dpStems)
  chain$beta <- c(outer, fixed)[refs$refPars]
  blocks <- Filter(function(b) b$dist != "normal", design$blocks)
  psiOf <- function() c(chain$full[names(outer)], chain$hyper)[names(start)]
  history <- matrix(NA_real_, 0, length(start), dimnames = list(NULL, names(start)))
  draws <- list()
  rows <- list()
  burnin <- TRUE; burnEnd <- NA_integer_; stationary <- FALSE
  for (k in seq_len(saem$nBurnin + saem$nEM)) {
    gamma <- if (burnin) 1 else (k - burnEnd)^-0.7
    L <- if (K) .choleskyFactor(lay, c(chain$full, fixed))
    Oinv <- if (K) chol2inv(t(L)) else matrix(0, 0, 0)
    logDet <- if (K) 2 * sum(log(diag(L))) else 0
    hyp <- .stemHypers(sl, c(chain$hyper, fixed))
    accept <- NA_real_
    moves <- c(birth = NA_real_, death = NA_real_, walk = NA_real_)
    if (K) {
      walk <- .saemNormalWalk(density, chain$full, chain$eta, Oinv, logDet, saem$nMcmc, chain$step, burnin)
      chain$eta <- walk$eta
      chain$step <- walk$step
      accept <- walk$accept
    }
    chain <- .saemClusterSweep(chain, Dp, dpStems, density, Oinv, logDet, hyp, burnin)
    if (length(slabStems)) {
      slabs <- .saemSlabMoves(chain, design, objs, E, Ds, reach, fixed, hyp, slabStems, density, Oinv, logDet,
                              saem, gh, burnin, spikeS, meanS, betaS)
      chain <- slabs$chain
      moves <- slabs$moves
    }
    if (K) chain <- .saemNormalStep(chain, design, ref, free, refs, outer, fixed, L, gamma, burnin)
    chain <- .saemHyperStep(chain, blocks, slabStems, reach, hypers, G, gamma, burnin)
    chain <- .saemStructuralStep(chain, completeData(chain$full, chain$eta), structural, saem$cm1, gamma, burnin)
    prev <- if (nrow(history)) history[nrow(history), ] else psiOf()
    psi <- psiOf()
    history <- rbind(history, psi)
    if (!burnin) draws[[length(draws) + 1L]] <- list(dp = lapply(chain$state, .canonicalClusters), incl = chain$incl,
                                                     pinc = if (length(slabStems)) chain$rbIncl,
                                                     v = matrix(chain$full[Ds], G, dimnames = dimnames(Ds)))
    if (burnin) {
      stationary <- k >= saem$window && all(abs(.slopeStatistics(utils::tail(history, saem$window))) < 2)
      if (stationary || k >= saem$nBurnin) {
        burnin <- FALSE; burnEnd <- k
        if (verbose) message(sprintf("em(saem) burn-in ends at iteration %d", k))
      }
    }
    rows[[k]] <- data.frame(iter = k, phase = if (is.na(burnEnd) || k <= burnEnd) "burnin" else "averaging",
                            gamma = gamma, deltaPsi = if (length(psi)) max(abs(psi - prev)) else 0,
                            deltaStruct = if (length(structural)) max(abs(psi[structural] - prev[structural])) else 0,
                            accept = accept, clusters = sum(vapply(chain$state, function(s) length(s$v), 0)),
                            included = sum(chain$incl), acceptBirth = moves[["birth"]], acceptDeath = moves[["death"]],
                            acceptWalk = moves[["walk"]])
    if (verbose && k %% 50L == 0L) message(sprintf("em(saem) iteration %d, %s", k, rows[[k]]$phase))
    if (!burnin && k - burnEnd >= saem$nEM) break
  }
  trace <- do.call(rbind, rows)
  late <- if (!is.na(burnEnd) && nrow(history) > burnEnd + saem$window) (burnEnd + saem$window + 1L):nrow(history)
          else seq_len(nrow(history))
  psi <- colMeans(history[late, , drop = FALSE])
  chain$full[names(outer)] <- psi[names(outer)]
  chain$hyper[] <- psi[hypers]

  keep <- utils::tail(draws, length(late))
  if (!length(keep)) keep <- list(list(dp = lapply(chain$state, .canonicalClusters), incl = chain$incl, pinc = chain$incl,
                                       v = matrix(chain$full[Ds], G, dimnames = dimnames(Ds))))
  picked <- .saemSelection(design, keep, chain$eta, dpStems, slabStems, reach, hyp, spikeS)
  if (general) {
    # -2 log L as em(method = "focei") defines it for a given spike or slab mean: the sum over subsets
    lay2 <- .spikeLayout(design)
    hypS <- .stemHypers(sl, c(chain$hyper, fixed))
    est <- .spikeEStep(design, objs, lay2, chain$full[.setdiff(names(outer), lay$cholPars)], fixed,
                       list(prob = hypS$prob[lay2$stem], tau = hypS$tau[lay2$stem]), list(), .spikeBases(design),
                       control$sparseFit$termBudget %||% .spikeBudget, control$sparseFit$ec, NULL, list(),
                       .spikeOmega(design, c(chain$full, fixed)))
    extra <- .extraObjective(design)
    # -2 log p(tau) of the scale priors at the estimate, and the posterior of prob under a beta prior
    penalty <- sum(vapply(unique(lay2$block), function(b)
      .spikeTauPenalty(hypS$tau[[lay2$stem[lay2$block == b][1]]], design$blocks[[b]]$tauPrior), 0))
    probPosterior <- Filter(Negate(is.null), setNames(lapply(unique(lay2$block), function(b)
      .spikeProbPosterior(lay2, est, b)), lay2$probName[match(unique(lay2$block), lay2$block)]))
    value <- sum(vapply(est, `[[`, 0, "m2logI")) + penalty +
      if (!is.null(extra)) extra(chain$full[.setdiff(names(outer), lay$cholPars)], fixed = fixed, deriv = FALSE)$value else 0
    cov <- lapply(design$groups, function(g) {
      e <- est[[g]]
      best <- e$parts[[which.max(vapply(e$parts, `[[`, 0, "weight"))]]
      C <- matrix(0, ncol(En), ncol(En), dimnames = list(colnames(En), colnames(En)))
      k <- match(best$on, En[g, ])
      C[k, k] <- best$cov
      C
    })
    return(list(argument = psi, value = value, gradient = NULL, hessian = NULL, etaModes = picked$modes,
                etaCov = setNames(cov, design$groups), clusters = picked$clusters, spike = picked$spike,
                converged = stationary, iterations = nrow(trace), stageTrace = trace, probPosterior = probPosterior,
                hyperprior = if (any(!vapply(design$blocks, function(b) is.null(b$tauPrior), TRUE))) penalty))
  }
  # -2 log L as em(method = "focei") defines it, with the second derivatives the models provide
  reachAll <- .reachOf(design)
  units <- lapply(design$groups, function(g)
    list(eta = unname(design$etaNames[g, reachAll[g, ]]), conditions = design$conditions[[g]]))
  derivMode <- .secondOrderMode(design, design$objfn, chain$full[names(outer)], fixed,
                                control[.setdiff(names(control), "secondOrderCorrection")], units = units)
  integral <- if (importance) list(integral = "importance", particles = saem$particles)
  r <- .sparseMarginal(design, objs, prior, c(chain$full[names(outer)], fixed), .stemHypers(sl, c(chain$hyper, fixed)),
                       as.character(.intersect(names(outer), names(start))), character(0), fixed, picked$modes,
                       grad = FALSE, control = c(integral, list(iterlim = 300L, derivMode = derivMode)))
  modes <- picked$modes
  modes[, lay$stems] <- r$modes[, lay$stems]
  list(argument = psi, value = r$value, valueSE = if (importance) r$valueSE, gradient = NULL, hessian = NULL,
       etaModes = modes, innerModes = r$modes, etaCov = setNames(r$cov, design$groups), clusters = picked$clusters,
       spike = picked$spike, converged = stationary, iterations = nrow(trace), stageTrace = trace)
}

## The settings of SAEM, `control$saem` with the defaults, and in `cm1` the trust settings of the Newton steps
## of the parameters outside the closed-form M-steps.
.saemSettings <- function(control) {
  sc <- control$saem %||% list()
  list(nBurnin = sc$nBurnin %||% 200L, nEM = sc$nEM %||% 300L, nMcmc = sc$nMcmc %||% 3L,
       nLaplace = sc$nLaplace %||% 3L, laplaceEvery = sc$laplaceEvery %||% 5L,
       laplaceScale = sc$laplaceScale %||% 1.5, window = sc$window %||% 25L, mcTol = sc$mcTol %||% 0.05,
       nFlip = sc$nFlip %||% 1L, nodes = sc$nodes %||% 3L, particles = sc$particles %||% 500L,
       integral = sc$integral,
       cm1 = .trustControl(list(rinit = 1, rmax = 10, iterlim = 1, tolControl = list(ftol = 1e-6, mtol = 1e-6)),
                           sc$cm1Control, label = "cm1Control"))
}

## The mu-referenced fixed effects `ref` on the etas `E`: their names, per group the matrix D_g, and
## shift(beta), the G x K matrix of the shifts D_g beta.
.referenceMatrices <- function(design, ref, E) {
  refPars <- unlist(lapply(unname(ref), `[[`, "pars"))
  K <- ncol(E)
  D <- lapply(design$groups, function(g) {
    m <- matrix(0, K, length(refPars), dimnames = list(colnames(E), refPars))
    for (e in ref) m[e$stem, e$pars] <- e$X[g, ]
    m
  })
  shift <- function(beta) matrix(vapply(D, function(m) drop(m %*% beta[refPars]), numeric(K)), nrow(E), K,
                                 byrow = TRUE, dimnames = dimnames(E))
  list(refPars = refPars, D = D, shift = shift)
}

## Random-walk Metropolis steps on the normal etas of all groups, the step sizes adapted while annealing.
.saemNormalWalk <- function(density, full, eta, Oinv, logDet, n, step, burnin) {
  current <- tryCatch(density(full, eta, Oinv, logDet), error = function(e) rep(Inf, nrow(eta)))
  rw <- .metropolisGroups(function(e) density(full, e, Oinv, logDet), eta, current, n, step = step)
  if (burnin) step <- pmin(pmax(step * exp(0.3 * (rw$accept - 0.4)), 0.02), 5)
  list(eta = rw$eta, accept = mean(rw$accept), step = step)
}

## Per Dirichlet-process stem a Gibbs sweep over the clusters and random-walk steps on the cluster values, the
## step sizes adapted while annealing; `chain$full` takes the values of the drawn clusters.
.saemClusterSweep <- function(chain, Dp, stems, density, Oinv, logDet, hyp, burnin) {
  full <- chain$full
  for (st in stems) {
    at <- Dp[, st]
    dataAt <- function(v) { f <- full; f[at] <- v; density(f, chain$eta, Oinv, logDet) }
    chain$state[[st]] <- .dpGibbs(chain$state[[st]], dataAt, hyp$concentration[[st]], hyp$tau[[st]])
    mv <- .dpValues(chain$state[[st]], dataAt, hyp$tau[[st]], chain$valueStep[[st]])
    chain$state[[st]] <- mv$s
    if (burnin && !is.na(mv$accept))
      chain$valueStep[[st]] <- min(max(chain$valueStep[[st]] * exp(0.3 * (mv$accept - 0.4)), 0.01), 5)
    full[at] <- c(0, chain$state[[st]]$v)[chain$state[[st]]$z + 1L]
  }
  chain$full <- full
  chain
}

## The indicators and slab values against proposals from the local models, which stay where they cannot be
## refreshed: per sweep the collapsed birth or death move on every indicator once, then random-walk steps on the
## slab values; the Rao-Blackwellized inclusion and slab second moment per coordinate and the acceptance rates.
.saemSlabMoves <- function(chain, design, objs, E, Ds, reach, fixed, hyp, slabStems, density, Oinv, logDet,
                           saem, gh, burnin, spikeS = 0 * hyp$tau[slabStems], meanS = 0 * hyp$tau[slabStems],
                           betaS = NULL) {
  G <- nrow(Ds)
  tauS <- hyp$tau[slabStems]; probS <- hyp$prob[slabStems]
  groupPars <- function() lapply(seq_len(G), function(g) {
    p <- chain$full; p[names(fixed)] <- fixed; p[E[g, ]] <- chain$eta[g, ]; p
  })
  local <- tryCatch(.slabModel(design, objs, groupPars(), Ds, chain$slabRef), error = function(e) NULL)
  if (!is.null(local)) {
    chain$slabProp <- local
    # the reference of the next model: the mean of the local Gaussian in the slab alone, within five slab
    # standard deviations
    for (g in seq_len(G)) chain$slabRef[g, ] <- tryCatch(
      pmin(pmax(.slabGaussian(local, g, reach[g, ], tauS, 1, meanS)$mean, meanS - 5 * tauS), meanS + 5 * tauS),
      error = function(e) chain$slabRef[g, ])
  }
  chain$slabProp <- chain$slabProp %||%
    list(ref = 0 * chain$slabRef, grad = 0 * chain$slabRef, H = rep(list(diag(0, ncol(Ds))), G))
  model <- chain$slabProp
  V <- matrix(chain$full[Ds], G)
  blockAt <- function(x) {
    set <- function(v) { f <- chain$full; f[Ds] <- v; f }
    density(if (is.list(x)) lapply(x, set) else set(x), chain$eta, Oinv, logDet)
  }
  value <- tryCatch(blockAt(V), error = function(e) rep(Inf, G))
  # the conditional standard deviations of the local model, the scales of the random walk
  walkScale <- t(vapply(model$H, function(h) sqrt(2 / (pmax(diag(h), 0) + 2 / tauS^2)), numeric(ncol(Ds))))
  if (ncol(Ds) == 1L) walkScale <- t(walkScale)
  incl <- chain$incl
  # for a spike at -Inf the proposal of coordinate j[g] is the Laplace approximation of the coordinate given
  # the others at V, by Newton steps on the model itself from the slab reference
  newton <- function(V, j) {
    minf <- !is.na(j) & !is.finite(spikeS[pmax(j, 1L)])
    if (!any(minf)) return(NULL)
    x <- chain$slabRef[cbind(seq_len(G), pmax(j, 1L))]
    curv <- rep(NA_real_, G)
    base <- groupPars()
    for (it in seq_len(3L)) {
      pars <- lapply(seq_len(G), function(g) {
        p <- base[[g]]; p[Ds[g, ]] <- V[g, ]
        if (minf[g]) p[Ds[g, j[g]]] <- x[g]
        p
      })
      out <- tryCatch(.groupData(design, objs, pars), error = function(e) NULL)
      if (is.null(out)) return(NULL)
      for (g in which(minf)) {
        nm <- Ds[g, j[g]]; k <- j[g]
        gr <- out[[g]]$gradient[[nm]] + 2 * (x[g] - meanS[k]) / tauS[k]^2
        h <- as.matrix(out[[g]]$hessian)[nm, nm] + 2 / tauS[k]^2
        if (!is.finite(gr) || !(h > 0)) { minf[g] <- FALSE; next }
        x[g] <- x[g] - max(min(gr / h, 2 * tauS[k]), -2 * tauS[k])
        curv[g] <- h
      }
    }
    list(m = ifelse(minf, x, NA_real_), sd = ifelse(minf, sqrt(2 / curv), NA_real_))
  }
  tried <- c(birth = 0, death = 0); taken <- tried
  rbIncl <- rbSlab2 <- rbCount <- matrix(0, G, ncol(Ds))
  for (s in seq_len(saem$nFlip)) {
    # every indicator of a group once, in random order, one round over all groups per indicator
    visit <- lapply(seq_len(G), function(g) { i <- which(reach[g, ]); i[sample.int(length(i))] })
    for (r in seq_len(max(lengths(visit), 0L))) {
      j <- vapply(visit, function(i) if (r <= length(i)) i[r] else NA_integer_, 0L)
      fl <- .slabCollapsed(V, incl, value, blockAt, model, probS, tauS, j, saem$laplaceScale, gh, spikeS, meanS,
                           newton(V, j), if (!is.null(betaS)) c(betaS, list(reach = reach)))
      V <- fl$V; incl <- fl$Z; value <- fl$value
      tried <- tried + c(sum(fl$birth, na.rm = TRUE), sum(!fl$birth, na.rm = TRUE))
      taken <- taken + c(sum(fl$ok & fl$birth, na.rm = TRUE), sum(fl$ok & !fl$birth, na.rm = TRUE))
      done <- which(!is.na(fl$pinc))
      at <- cbind(done, j[done])
      rbIncl[at] <- rbIncl[at] + fl$pinc[done]
      rbSlab2[at] <- rbSlab2[at] + fl$slab2[done]
      rbCount[at] <- rbCount[at] + 1
    }
  }
  # a coordinate whose collapsed step failed in every sweep enters with the draw
  chain$rbIncl <- ifelse(rbCount > 0, rbIncl / pmax(rbCount, 1), incl)
  chain$rbSlab2 <- ifelse(rbCount > 0, rbSlab2 / pmax(rbCount, 1), ifelse(incl, (V - rep(meanS, each = G))^2, 0))
  walked <- numeric(G)
  for (s in seq_len(saem$nMcmc)) {
    sw <- .slabWalk(V, incl, value, blockAt, tauS, chain$slabStep, walkScale, meanS)
    V <- sw$V; value <- sw$value
    walked <- walked + sw$accept / saem$nMcmc
  }
  moves <- c(taken / tried, walk = sum(walked) / max(sum(rowSums(incl) > 0), 1))
  chain$full[Ds] <- V
  chain$incl <- incl
  if (burnin)
    chain$slabStep <- ifelse(rowSums(incl) > 0, pmin(pmax(chain$slabStep * exp(0.3 * (walked - 0.25)), 0.01), 5),
                             chain$slabStep)
  list(chain = chain, moves = moves)
}

## The A-step on phi = eta + D beta and the M-steps of Omega and the GLS fixed effects on the averages of phi;
## while annealing the diagonal of S stays at least 0.95 of the current variances.
.saemNormalStep <- function(chain, design, ref, free, refs, outer, fixed, L, gamma, burnin) {
  G <- nrow(chain$eta)
  phi <- chain$eta + refs$shift(chain$beta)
  chain$s1 <- if (is.null(chain$s1)) phi else chain$s1 + gamma * (phi - chain$s1)
  chain$s2 <- if (is.null(chain$s2)) crossprod(phi) / G else chain$s2 + gamma * (crossprod(phi) / G - chain$s2)
  for (it in 1:3) {
    sh <- refs$shift(chain$beta)
    S <- chain$s2 - (crossprod(sh, chain$s1) + crossprod(chain$s1, sh) - crossprod(sh)) / G
    if (burnin) diag(S) <- pmax(diag(S), 0.95 * diag(tcrossprod(L)))
    cholPars <- .updateNormal(design, list(S), fixed)
    held <- .intersect(names(cholPars), names(outer))
    chain$full[held] <- cholPars[held]
    if (!length(free)) break
    Omega <- tcrossprod(.choleskyFactor(design$normal, c(chain$full, fixed)))
    chain$beta[free] <- .glsFixedEffects(ref, chain$s1, Omega, chain$beta, free)
  }
  chain$full[free] <- chain$beta[free]
  chain$eta <- phi - refs$shift(chain$beta)
  chain
}

## The averaged statistics of every sparse block (free coordinates, groups, slab second moments) and the M-steps
## of the slab sd, the concentration and the probability of a deviation; while annealing each changes by a
## factor of 0.95 to 1.25 at most per iteration, so that one sweep cannot drive it to a limit.
.saemHyperStep <- function(chain, blocks, slabStems, reach, hypers, G, gamma, burnin) {
  current <- vapply(blocks, function(b) {
    if (b$dist == "dirichletProcess")
      return(c(free = sum(vapply(b$stems, function(st) length(chain$state[[st]]$v), 0)), groups = 0,
               slab2 = sum(vapply(b$stems, function(st) sum(chain$state[[st]]$v^2), 0))))
    j <- match(b$stems, slabStems)
    c(free = sum(chain$rbIncl[, j][reach[, j]]), groups = sum(reach[, j]), slab2 = sum(chain$rbSlab2[, j][reach[, j]]))
  }, numeric(3))
  chain$blockStats <- if (is.null(chain$blockStats)) current else chain$blockStats + gamma * (current - chain$blockStats)
  averaged <- chain$blockStats
  bound <- function(new, old) if (burnin) min(max(new, 0.95 * old), 1.25 * old) else new
  odds <- function(p) p / (1 - p)
  for (j in seq_along(blocks)) {
    b <- blocks[[j]]
    tn <- unname(b$hyper["tau"]); cn <- unname(b$hyper["concentration"]); pn <- unname(b$hyper["prob"])
    if (tn %in% hypers && averaged["free", j] > 1e-8)
      chain$hyper[[tn]] <- bound(.spikeTauMode(averaged["slab2", j], averaged["free", j], b$tauPrior,
                                               chain$hyper[[tn]]), chain$hyper[[tn]])
    if (cn %in% hypers)
      chain$hyper[[cn]] <- bound(.ewensConcentration(averaged["free", j], rep(G, length(b$stems))), chain$hyper[[cn]])
    if (pn %in% hypers && averaged["groups", j] > 0) {
      o <- bound(odds(min(max(averaged["free", j] / averaged["groups", j], 1e-8), 1 - 1e-8)), odds(chain$hyper[[pn]]))
      chain$hyper[[pn]] <- o / (1 + o)
    }
  }
  chain
}

## The structural parameters: a trust step on the complete-data objective at the draws while annealing, then the
## step on the averaged quadratic model; where the objective fails at the draws they keep their values.
.saemStructuralStep <- function(chain, objective, structural, cm1, gamma, burnin) {
  if (!length(structural)) return(chain)
  if (burnin) {
    cm <- tryCatch(.innerTrust(modifyList(cm1, list(objfun = objective, parinit = chain$full[structural]))),
                   error = function(e) NULL)
    if (!is.null(cm)) chain$full[structural] <- cm$argument
    return(chain)
  }
  out <- tryCatch(objective(chain$full[structural]), error = function(e) NULL)
  if (is.null(out)) return(chain)
  chain$model <- .averageQuadratic(if (gamma < 1) chain$model, out, chain$full[structural], gamma)
  chain$full[structural] <- .quadraticStep(chain$model, chain$full[structural])
  chain
}

## The selection from the late draws `keep`: per Dirichlet-process stem the most frequent partition with the
## cluster values averaged over the draws of that partition; per spike-and-slab coordinate the mean conditional inclusion
## probability, selected above 1/2 with the mean slab value.
.saemSelection <- function(design, keep, eta, dpStems, slabStems, reach, hyp, spikeS = rep(0, length(slabStems))) {
  G <- nrow(design$etaNames)
  modes <- matrix(0, G, ncol(design$etaNames), dimnames = dimnames(design$etaNames))
  modes[, design$normal$stems] <- eta
  clusters <- spike <- matrix(NA_real_, G, ncol(design$etaNames), dimnames = dimnames(design$etaNames))
  for (st in dpStems) {
    keys <- vapply(keep, function(d) paste(d$dp[[st]]$z, collapse = " "), "")
    top <- names(which.max(table(keys)))
    same <- keep[keys == top]
    z <- same[[1]]$dp[[st]]$z
    values <- colMeans(do.call(rbind, lapply(same, function(d) d$dp[[st]]$v)) %||% matrix(0, 1, 0))
    modes[, st] <- c(0, values)[z + 1L]
    clusters[, st] <- z
    spike[, st] <- colMeans(do.call(rbind, lapply(keep, function(d) d$dp[[st]]$z == 0L)))
  }
  if (length(slabStems)) {
    Z <- Reduce(`+`, lapply(keep, `[[`, "incl"))
    V <- Reduce(`+`, lapply(keep, function(d) ifelse(d$incl, d$v, 0)))
    inclusion <- Reduce(`+`, lapply(keep, `[[`, "pinc")) / length(keep)
    selected <- inclusion > 0.5
    modes[, slabStems] <- ifelse(selected & Z > 0, V / pmax(Z, 1), matrix(spikeS, G, length(slabStems), byrow = TRUE))
    spike[, slabStems] <- ifelse(reach, 1 - inclusion, 1 - rep(hyp$prob[slabStems], each = G))
    clusters[, slabStems] <- apply(selected, 2, function(s) ifelse(s, cumsum(s), 0))
  }
  list(modes = modes, clusters = clusters, spike = spike)
}
