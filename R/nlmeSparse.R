# Spike-and-slab and Dirichlet-process random effects, normal ones beside them; marginal, exact
# gradient and M-steps as derived in the methods article Sparsity priors.


## ---- design pieces ---------------------------------------------------------------

## The spike-and-slab and cluster coordinates: per stem the block, distribution, the names of the
## hyperparameters and their fixed values, NA when estimated.
.sparseLayout <- function(design) {
  blocks <- Filter(function(b) b$dist != "normal", design$blocks)
  name <- function(b, kind) if (kind %in% names(b$hyper)) b$hyper[[kind]] else NA_character_
  do.call(rbind, lapply(seq_along(blocks), function(i) {
    b <- blocks[[i]]
    data.frame(stem = b$stems, block = i, dist = b$dist, probName = name(b, "prob"),
               concentrationName = name(b, "concentration"), tauName = name(b, "tau"), prob = b$prob,
               concentration = b$concentration, tau = b$tau, stringsAsFactors = FALSE)
  }))
}

## The hyperparameters of every stem, prob, concentration and tau, from `pars` or the fixed values.
.stemHypers <- function(lay, pars) {
  value <- function(fixed, name) setNames(ifelse(is.na(fixed), unname(pars[name]), fixed), lay$stem)
  list(prob = value(lay$prob, lay$probName), concentration = value(lay$concentration, lay$concentrationName),
       tau = value(lay$tau, lay$tauName))
}

## The estimated hyperparameters of the layout.
.sparseHypers <- function(lay) {
  h <- c(lay$probName[is.na(lay$prob)], lay$concentrationName[is.na(lay$concentration)],
         lay$tauName[is.na(lay$tau)])
  unique(h[!is.na(h)])
}

## The etas of every group that reach the group data (.etaReach()), a G x K logical like design$etaNames.
.reachOf <- function(design) {
  E <- design$etaNames
  design$reach %||% matrix(TRUE, nrow(E), ncol(E), dimnames = dimnames(E))
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

## The terms of the objective of `design` besides the data, such as priors and constraints; NULL
## without.
.extraObjective <- function(design) {
  terms <- Filter(function(t) is.null(attr(t, "l2spec", exact = TRUE)), .objTerms(design$objfn))
  if (length(terms)) Reduce(`+`, terms)
}

## The data -2 log L of every group at the group parameter vector, summed over the group conditions.
.groupData <- function(design, objs, parsByGroup, deriv = TRUE, deriv2 = FALSE) {
  conds <- unlist(design$conditions, use.names = FALSE)
  owner <- rep(seq_along(design$groups), lengths(design$conditions))
  outs <- .objEvalMany(objs[conds], parsByGroup[owner], deriv = deriv, deriv2 = deriv2)
  lapply(seq_along(design$groups), function(i) Reduce(`+`, outs[owner == i]))
}


## ---- inner mode --------------------------------------------------------------------

## The modes of the etas at `theta` by the group kernel, the slab alone per sparse eta, etas off the group data
## at 0, each group from the lowest of `start` and `control$optima`; attributes "iterations", the inner
## iterations summed over the groups, and with `control$phi` "phi", the inner values.
.sparseModes <- function(design, theta, hyp, start, control = list()) {
  lay <- .sparseLayout(design)
  E <- design$etaNames
  K <- ncol(E)
  reach <- .reachOf(design)
  isL <- colnames(E) %in% lay$stem
  meta <- control$meta %||% .foceiMeta(design, theta, NULL, colnames(E), reach)
  Pm <- matrix(0, K, K)
  ld <- 0
  if (any(!isL)) {
    L <- .choleskyFactor(design$normal, theta)
    Pm[!isL, !isL] <- chol2inv(t(L))
    ld <- 2 * sum(log(diag(L)))
  }
  diag(Pm)[isL] <- 1 / hyp$tau[colnames(E)[isL]]^2
  predictor <- .guardPredictor(.batchPredictor(design))
  full <- setNames(numeric(length(meta$pars_full_names)), meta$pars_full_names)
  full[names(theta)] <- theta
  iterations <- 0L
  run <- function(warm, iterlim) {
    warm <- matrix(as.vector(warm), nrow(E))
    warm[!reach] <- 0
    k <- nlme_group_kernel(design$prdfn, .errorModel(design), full, warm, meta, Pm, ld, NULL,
                           list(iterlim = iterlim), predict_cb = predictor)
    iterations <<- iterations + sum(k$iterations)
    list(modes = matrix(as.vector(k$eta_modes), nrow(E), dimnames = dimnames(E)), phi = as.numeric(k$value))
  }
  # the starts: per group the row of `start` or of an earlier optimum with the lowest value
  start <- matrix(as.vector(start), nrow(E), dimnames = dimnames(E))
  depth <- max(0L, vapply(control$optima, function(o) if (is.null(o)) 0L else nrow(o), 0L))
  if (depth) {
    low <- run(start, 0L)$phi
    for (r in seq_len(depth)) {
      cand <- start
      for (i in seq_len(nrow(E))) if (!is.null(control$optima[[i]]) && nrow(control$optima[[i]]) >= r)
        cand[i, ] <- control$optima[[i]][r, ]
      phi <- run(cand, 0L)$phi
      lower <- phi < low
      start[lower, ] <- cand[lower, ]
      low[lower] <- phi[lower]
    }
  }
  best <- run(start, control$iterlim %||% 100L)
  modes <- best$modes
  modes[!reach] <- 0
  if (isTRUE(control$phi)) attr(modes, "phi") <- setNames(best$phi, rownames(E))
  attr(modes, "iterations") <- iterations
  modes
}

## The lowest inner minimum per group at `theta` from the modes `warm`, from 0 and from `spread`
## starts around `warm`, each eta moved by `scale` times a deterministic pseudo-normal number: the
## inner problem of a group may have several minima, and every search keeps to the one it starts in.
.exploreModes <- function(design, theta, hyp, warm, control, scale, spread = 2L) {
  ctl <- utils::modifyList(control, list(phi = TRUE, optima = NULL))
  best <- .sparseModes(design, theta, hyp, warm, ctl)
  low <- attr(best, "phi")
  iterations <- attr(best, "iterations")
  n <- length(warm)
  for (k in 0:spread) {
    z <- if (k == 0L) 0 else stats::qnorm(pmin(pmax(((sin(seq_len(n) * k * 12.9898 + 78.233) * 43758.5453) %% 1), 0.01), 0.99))
    start <- if (k == 0L) 0 * warm else warm + scale * matrix(z, nrow(warm))
    m <- .sparseModes(design, theta, hyp, start, ctl)
    iterations <- iterations + attr(m, "iterations")
    lower <- !is.na(attr(m, "phi")) & (is.na(low) | attr(m, "phi") < low - 1e-6)
    best[lower, ] <- m[lower, ]
    low[lower] <- attr(m, "phi")[lower]
  }
  attr(best, "phi") <- NULL
  attr(best, "iterations") <- iterations
  best
}


## ---- curvature -------------------------------------------------------------------------

## The Gauss-Newton Hessian of the data -2 log L of a group in the group etas, sum over the group rows of
## J' L J with J = d(p, s) / d eta and L the Hessian of the row loss, and the derivatives of that Hessian in the
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

## The factors of the sparse coordinates: every Dirichlet-process stem alone, the spike-and-slab
## stems of a group jointly; `coupled` when the data precision links two factors.
.sparseFactors <- function(P, stems, lay) {
  r <- match(stems, lay$stem)
  slabs <- which(lay$dist[r] == "spikeSlab")
  factors <- c(as.list(.setdiff(seq_along(stems), slabs)), if (length(slabs)) list(slabs))
  owner <- integer(length(stems))
  for (j in seq_along(factors)) owner[factors[[j]]] <- j
  within <- outer(owner, owner, `==`)
  list(factors = factors, rows = r, within = within,
       coupled = any(vapply(P, function(Pi) any(Pi[!within] != 0), TRUE)))
}

## I = int exp(-eta' P_i eta / 2 + b_i' eta) prod p(eta) per group (src/sparse_kernel.cpp): exact within `budget`
## nodes where P couples no factors, else by EC from `sites` (`ec` marks those groups); log I, moments, hyper
## derivatives, values, clusters, spike and M-step statistics, NULL where it does not converge.
.sparseIntegral <- function(P, b, lay, hyp, hypers, sites = NULL, reach = NULL, budget = .slabBudget) {
  G <- nrow(b); K <- ncol(b); stems <- colnames(b)
  r <- match(stems, lay$stem)
  kind <- match(lay$dist[r], c("spikeSlab", "dirichletProcess")) - 1L
  reach <- reach %||% matrix(TRUE, G, K)
  num <- function(v) as.numeric(unname(v[stems]))
  sp <- sparse_integral(lapply(P, function(Pi) { Pi <- as.matrix(Pi); storage.mode(Pi) <- "double"; Pi }),
                        unname(b), kind, num(hyp$prob), num(hyp$tau), num(hyp$concentration),
                        unname(reach), sites$Lam, sites$gam, as.numeric(budget))
  if (!isTRUE(sp$ok)) return(NULL)
  dn <- dimnames(b)
  named <- function(x) `dimnames<-`(x, dn)
  # the derivatives in the hyperparameters by name, summed over the stems that share them
  dhyper <- setNames(numeric(length(hypers)), hypers)
  hyperGroups <- matrix(0, G, length(hypers), dimnames = list(rownames(b), hypers))
  add <- function(name, total, groups) {
    for (j in which(!is.na(name) & name %in% hypers)) {
      dhyper[[name[j]]] <<- dhyper[[name[j]]] + total[j]
      hyperGroups[, name[j]] <<- hyperGroups[, name[j]] + if (is.null(groups)) NA_real_ else groups[, j]
    }
  }
  slab <- kind == 0L; dp <- kind == 1L
  add(ifelse(slab, lay$probName[r], NA), sp$dprob, sp$gprob)
  add(ifelse(slab, lay$tauName[r], NA), sp$dtau, sp$gtau)
  add(ifelse(dp, lay$concentrationName[r], NA), sp$dconc, NULL)
  add(ifelse(dp, lay$tauName[r], NA), sp$dtau, NULL)
  # a Dirichlet-process stem couples the groups: no derivatives per group
  if (any(dp)) hyperGroups[] <- NA_real_
  stats <- cbind(free = sp$stats$free, groups = sp$stats$groups, slab2 = sp$stats$slab2)
  slabStats <- matrix(0, nrow(lay), 3, dimnames = list(lay$stem, c("free", "groups", "slab2")))
  slabStats[stems[slab | dp], ] <- stats[slab | dp, , drop = FALSE]
  list(logI = sp$logI, mean = named(sp$mean), cov = lapply(sp$cov, `dimnames<-`, list(stems, stems)),
       values = named(sp$values), clusters = named(sp$clusters),
       spike = named(pmin(pmax(sp$spike, 0), 1)),
       dhyper = dhyper, hyperGroups = hyperGroups, slab = slabStats,
       coupling = if (sp$coupling != "none") sp$coupling, ec = setNames(sp$ec, rownames(b)),
       sites = if (any(sp$ec)) list(Lam = named(sp$Lam), gam = named(sp$gam)),
       sweeps = sp$outer)
}

## The integral of Dirichlet-process stems, spike-and-slab ones beside them, by sequential
## importance sampling (.dirichletImportance), with the standard error of log I; the moments and
## derivatives of the mean-field result are not formed, they stay 0 or NA.
.importanceField <- function(P, m, lay, hyp, hypers, particles) {
  G <- nrow(m); K <- ncol(m); stems <- colnames(m)
  dist <- lay$dist[match(stems, lay$stem)]
  if (!any(dist == "dirichletProcess"))
    stop("em: importance sampling integrates Dirichlet-process random effects.", call. = FALSE)
  is <- .dirichletImportance(P, m, hyp$concentration[stems], hyp$tau[stems], particles,
                             slab = dist == "spikeSlab", prob = hyp$prob[stems])
  none <- matrix(NA_real_, G, K, dimnames = dimnames(m))
  list(logI = is$logI, se = is$se, dm = 0 * m, dP = lapply(P, function(x) 0 * x),
       dhyper = setNames(numeric(length(hypers)), hypers), mean = none,
       cov = lapply(seq_len(G), function(i) matrix(NA_real_, K, K, dimnames = list(stems, stems))),
       values = none, clusters = none, spike = none,
       slab = matrix(0, nrow(lay), 3, dimnames = list(lay$stem, c("free", "groups", "slab2"))))
}

## ---- marginal ------------------------------------------------------------------------------

## The approximate marginal -2 log L at `theta`: the sparse etas integrated on the Schur complement S at the
## joint mode, P = S / 2 and b = (S mode - g) / 2; the exact gradient in `active` and `hypers`, a model Hessian,
## the moments and clusters, with `scores` the gradient per group. Etas off the group data keep the prior.
.sparseMarginal <- function(design, objs, prior, theta, hyp, active, hypers, fixed, warm,
                            grad = TRUE, control = list(), scores = FALSE) {
  modes <- .sparseModes(design, theta, hyp, warm, control)
  innerIterations <- attr(modes, "iterations")
  attr(modes, "iterations") <- NULL
  at <- .sparsePoint(design, theta, modes, identical(control$derivMode, "forward-reverse"))
  terms <- .sparseTerms(design, objs, prior, at, fixed, grad)
  curv <- .sparseCurvature(design, at, active, fixed)
  blocks <- .sparseBlocks(at, curv, terms, fixed)
  importance <- identical(control$integral, "importance")
  integral <- if (importance) .sparseImportanceIntegral(at, blocks, hyp, hypers, control$particles %||% 500L)
              else .sparseIntegral(lapply(blocks$schur, `/`, 2), blocks$b, at$lay, hyp, hypers, control$sites,
                                   at$reach[, at$sparseCols, drop = FALSE])
  if (is.null(integral)) return(list(infeasible = TRUE, innerIterations = innerIterations))
  out <- .sparseMoments(at, terms, curv, blocks, integral, innerIterations)
  if (!grad || importance) return(out)

  d <- .sparseGroupDerivatives(design, at, terms, curv, blocks, integral, active, fixed)
  model <- .sparseModelHessian(at, blocks, integral, d)
  slabPrecision <- .slabPrecision(at$E, at$lay, hyp)
  diagonal <- cbind(names(slabPrecision), names(slabPrecision))
  d$etaHessian[diagonal] <- d$etaHessian[diagonal] + slabPrecision
  out$value <- out$value + .modeCorrection(at, terms, d, slabPrecision)
  outer <- .sparseOuterDerivatives(at, terms, blocks, integral, d, hyp, hypers, active, fixed, scores)
  out$modeShift <- outer$modeShift
  out$gradient <- outer$gradient
  out$hessian <- model %||% outer$schur
  if (scores)
    out[c("scores", "freeGradient", "freeHessian")] <-
      list(outer$byGroup, d$freeGradient, .freeInformation(design, objs, terms, at, active, fixed, d$freeHessian))
  out
}

## The point of .sparseMarginal(): the layout, the normal and sparse columns of the etas, the modes and the
## full parameter vector.
.sparsePoint <- function(design, theta, modes, forwardReverse) {
  lay <- .sparseLayout(design)
  E <- design$etaNames
  isSparse <- colnames(E) %in% lay$stem
  list(lay = lay, normal = design$normal, E = E, G = nrow(E), K = ncol(E), normalCols = which(!isSparse),
       sparseCols = which(isSparse), reach = .reachOf(design), modes = modes,
       full = c(theta, setNames(as.vector(modes), as.vector(E))), forwardReverse = forwardReverse)
}

## The data terms of the groups, the normal prior, the terms without random effects (free conditions,
## priors, constraints) and the data gradient in the etas, 0 off the group data. Under forward over reverse
## the data Hessians are Gauss-Newton.
.sparseTerms <- function(design, objs, prior, at, fixed, grad) {
  deriv2 <- grad && !at$forwardReverse
  groups <- .groupData(design, objs, rep(list(at$full), at$G), deriv = TRUE, deriv2 = deriv2)
  normalPrior <- if (!is.null(prior)) prior(pars = at$full, fixed = fixed, deriv = TRUE)
  free <- if (length(design$freeConditions))
    .objEvalMany(objs[design$freeConditions], rep(list(at$full), length(design$freeConditions)),
                 deriv = TRUE, deriv2 = deriv2)
  extra <- .extraObjective(design)
  if (!is.null(extra)) free <- c(free, list(extra(pars = at$full, fixed = fixed, deriv = TRUE, deriv2 = deriv2)))
  etaGradient <- matrix(t(vapply(seq_len(at$G), function(i) groups[[i]]$gradient[at$E[i, ]], numeric(at$K))), at$G,
                        dimnames = dimnames(at$E))
  etaGradient[!at$reach] <- 0
  list(groups = groups, prior = normalPrior, free = free, extra = extra, etaGradient = etaGradient)
}

## The data curvature H of every group on all etas of the group, 0 in those off the group data: from forward
## sweeps under forward over reverse, else Gauss-Newton with the derivatives dA in `active` and dE in the etas.
.sparseCurvature <- function(design, at, active, fixed) {
  units <- lapply(seq_len(at$G), function(i)
    list(eta = unname(at$E[i, at$reach[i, ]]), conditions = design$conditions[[design$groups[i]]]))
  if (at$forwardReverse) {
    forward <- .sweepForward(design, at$full, fixed, units)
    perUnit <- .unitCurvature(design, forward, units)
    groups <- lapply(seq_len(at$G), function(i) {
      H <- matrix(0, at$K, at$K)
      H[at$reach[i, ], at$reach[i, ]] <- perUnit[[i]]
      list(H = H)
    })
    return(list(groups = groups, units = units, forward = forward))
  }
  conds <- unlist(design$conditions, use.names = FALSE)
  times <- lapply(conds, function(cn) .emGrid(design, cn))
  preds <- setNames(.predictMany(design$prdfn, times = times, parsList = rep(list(at$full), length(conds)),
                                 conditions = conds, fixed = fixed, deriv = TRUE, deriv2 = TRUE), conds)
  groups <- lapply(seq_len(at$G), function(i) .gnCurvature(design, design$groups[i], preds, at$E[i, ], active, fixed))
  list(groups = groups, units = units)
}

## A G x (sparse etas) matrix from `f(i)`, the row of group i.
.bySparse <- function(at, f) {
  n <- length(at$sparseCols)
  matrix(vapply(seq_len(at$G), f, numeric(n)), at$G, n, byrow = TRUE,
         dimnames = list(rownames(at$E), colnames(at$E)[at$sparseCols]))
}

## Per group the normal block H_NN plus the prior precision 2 Omega^-1 with the volume sum log|H_NN / 2|, and the
## sparse block: S, the Schur complement of the normal block, and b, the data term on it being
## exp(-eta' P eta / 2 + b' eta) up to `offset`, mode' S mode / 2 - g' mode.
.sparseBlocks <- function(at, curv, terms, fixed) {
  normal <- at$normalCols; sparse <- at$sparseCols
  cholesky <- if (length(normal)) .choleskyDerivatives(at$normal, .choleskyFactor(at$normal, c(at$full, fixed)))
  normalHessian <- lapply(curv$groups, function(cv) if (length(normal)) cv$H[normal, normal, drop = FALSE] + 2 * cholesky$Oinv)
  volume <- sum(vapply(normalHessian, function(H) if (is.null(H)) 0 else as.numeric(determinant(H / 2)$modulus), 0))
  schur <- lapply(seq_len(at$G), function(i) {
    H <- curv$groups[[i]]$H
    s <- H[sparse, sparse, drop = FALSE]
    if (length(normal)) s <- s - H[sparse, normal, drop = FALSE] %*% solve(normalHessian[[i]], H[normal, sparse, drop = FALSE])
    (s + t(s)) / 2
  })
  sparseGradient <- terms$etaGradient[, sparse, drop = FALSE]
  sparseModes <- at$modes[, sparse, drop = FALSE]
  b <- .bySparse(at, function(i) drop(schur[[i]] %*% sparseModes[i, ] - sparseGradient[i, ]) / 2)
  offset <- sum(vapply(seq_len(at$G), function(i)
    sum(sparseModes[i, ] * (schur[[i]] %*% sparseModes[i, ])) / 2 - sum(sparseGradient[i, ] * sparseModes[i, ]), 0))
  list(cholesky = cholesky, normalHessian = normalHessian, volume = volume, schur = schur, b = b, offset = offset,
       sparseGradient = sparseGradient, sparseModes = sparseModes)
}

## The sparse integral by importance sampling in center form, m = mode - S^-1 g: log I of the natural form is
## that of the center form plus m' P m / 2.
.sparseImportanceIntegral <- function(at, blocks, hyp, hypers, particles) {
  sparse <- at$sparseCols
  center <- blocks$sparseModes - .bySparse(at, function(i) {
    on <- at$reach[i, sparse]
    out <- numeric(length(sparse))
    out[on] <- solve(blocks$schur[[i]][on, on, drop = FALSE], blocks$sparseGradient[i, on])
    out
  })
  sampled <- .importanceField(lapply(blocks$schur, `/`, 2), center, at$lay, hyp, hypers, particles)
  sampled$logI <- sampled$logI + sum(vapply(seq_len(at$G), function(i)
    sum(center[i, ] * (blocks$schur[[i]] %*% center[i, ])) / 4, 0))
  sampled
}

## The value, the selected values, clusters and spike probabilities, and the posterior moments of the etas;
## the normal etas given the sparse ones are normal with a mean linear in them.
.sparseMoments <- function(at, terms, curv, blocks, integral, innerIterations) {
  G <- at$G; K <- at$K; E <- at$E
  normal <- at$normalCols; sparse <- at$sparseCols
  value <- sum(vapply(terms$groups, `[[`, 0, "value")) + (terms$prior$value %||% 0) + blocks$volume -
    2 * integral$logI + blocks$offset + sum(vapply(terms$free, `[[`, 0, "value"))
  selected <- at$modes
  selected[, sparse] <- ifelse(is.na(integral$values), at$modes[, sparse], integral$values)
  clusters <- spike <- matrix(NA_real_, G, K, dimnames = dimnames(E))
  clusters[, sparse] <- integral$clusters; spike[, sparse] <- integral$spike
  etaMean <- matrix(0, G, K, dimnames = dimnames(E))
  etaMean[, sparse] <- integral$mean
  etaCov <- lapply(seq_len(G), function(i) {
    V <- matrix(0, K, K, dimnames = list(colnames(E), colnames(E)))
    V[sparse, sparse] <- integral$cov[[i]]
    if (length(normal)) {
      B <- solve(blocks$normalHessian[[i]], curv$groups[[i]]$H[normal, sparse, drop = FALSE])
      VS <- V[sparse, sparse, drop = FALSE]
      VS[is.na(VS)] <- 0
      V[normal, normal] <- 2 * solve(blocks$normalHessian[[i]]) + B %*% VS %*% t(B)
      V[normal, sparse] <- -B %*% VS
      V[sparse, normal] <- t(V[normal, sparse, drop = FALSE])
    }
    V
  })
  etaVar <- matrix(vapply(etaCov, diag, numeric(K)), G, K, byrow = TRUE, dimnames = dimnames(E))
  list(value = value, modes = at$modes, innerIterations = innerIterations, mean = etaMean, var = etaVar,
       cov = etaCov, slab = integral$slab, selected = selected, clusters = clusters, spike = spike,
       valueSE = 2 * (integral$se %||% 0), coupling = integral$coupling, sites = integral$sites, reach = at$reach,
       local = list(P = lapply(blocks$schur, `/`, 2), b = blocks$b))
}

## Per group the matrices the derivatives of the curvature are contracted with: `weight` on all etas, from
## delta = E eta - mode and Cov eta of the integral, and `normalInverse`, H_NN^-1, on the normal ones.
.curvatureContractions <- function(at, curv, blocks, integral, delta) {
  normal <- at$normalCols; sparse <- at$sparseCols
  lapply(seq_len(at$G), function(i) {
    weight <- matrix(0, at$K, at$K, dimnames = list(at$E[i, ], at$E[i, ]))
    weight[sparse, sparse] <- (tcrossprod(delta[i, ]) + integral$cov[[i]]) / 2
    normalInverse <- NULL
    if (length(normal)) {
      normalInverse <- solve(blocks$normalHessian[[i]])
      B <- curv$groups[[i]]$H[sparse, normal, drop = FALSE] %*% normalInverse
      Ms <- weight[sparse, sparse, drop = FALSE]
      weight[sparse, normal] <- -Ms %*% B
      weight[normal, sparse] <- t(weight[sparse, normal, drop = FALSE])
      weight[normal, normal] <- crossprod(B, Ms %*% B)
    }
    list(weight = weight, normalInverse = normalInverse)
  })
}

## For group i the exact rows of the data Hessian at the group etas `en` and in `an`, and the traces
## tr(weight dH) and tr(H_NN^-1 dH_NN) along `active` and the etas: from the sweeps `sweep` under forward over
## reverse, else from the curvature derivatives `cv`.
.groupTraces <- function(at, term, cv, contraction, sweep, en, an, active) {
  if (at$forwardReverse) {
    rowsAt <- function(cols) {
      out <- matrix(0, at$K, length(cols), dimnames = list(en, cols))
      h <- .intersect(cols, colnames(sweep$rows))
      out[rownames(sweep$rows), h] <- sweep$rows[, h]
      out
    }
    pick <- function(v, nm) { out <- setNames(numeric(length(nm)), nm); h <- .intersect(nm, names(v)); out[h] <- v[h]; out }
    return(list(exactE = rowsAt(en), exactA = rowsAt(an),
                dMA = pick(sweep$traces[[1]], active), dME = pick(sweep$traces[[1]], en),
                dHA = pick(sweep$traces[[2]], active), dHE = pick(sweep$traces[[2]], en)))
  }
  normal <- at$normalCols
  weight <- contraction$weight
  out <- list(exactE = term$hessian[en, en, drop = FALSE], exactA = term$hessian[en, an, drop = FALSE],
              dMA = vapply(seq_along(active), function(q) sum(weight * cv$dA[, , q]), 0),
              dME = vapply(seq_len(at$K), function(j) sum(weight * cv$dE[, , j]), 0))
  if (length(normal)) {
    out$dHA <- vapply(seq_along(active), function(q) sum(contraction$normalInverse * cv$dA[normal, normal, q]), 0)
    out$dHE <- vapply(seq_len(at$K), function(j) sum(contraction$normalInverse * cv$dE[normal, normal, j]), 0)
  }
  out
}

## The derivatives of the value at fixed modes: the explicit gradient in `active` per group, the gradient in
## the etas, the exact eta and cross blocks of the Hessian, the model ones (Gauss-Newton under forward over
## reverse, so that the Schur complement is positive semidefinite), and the block of `active`.
.sparseGroupDerivatives <- function(design, at, terms, curv, blocks, integral, active, fixed) {
  G <- at$G; K <- at$K; E <- at$E
  normal <- at$normalCols; sparse <- at$sparseCols
  nl <- at$normal
  delta <- integral$mean - blocks$sparseModes
  nActive <- length(active)
  etaOrder <- as.vector(t(E))
  groupGradient <- matrix(0, G, nActive, dimnames = list(rownames(E), active))
  etaGradient <- setNames(numeric(G * K), etaOrder)
  etaHessian <- matrix(0, G * K, G * K, dimnames = list(etaOrder, etaOrder))
  crossHessian <- matrix(0, G * K, nActive, dimnames = list(etaOrder, active))
  etaModel <- etaHessian; crossModel <- crossHessian
  activeHessian <- matrix(0, nActive, nActive, dimnames = list(active, active))
  contraction <- .curvatureContractions(at, curv, blocks, integral, delta)
  sweeps <- if (at$forwardReverse)
    .sweepUnits(design, at$full, fixed, curv$units, fw = curv$forward, Ms = lapply(seq_len(G), function(i) {
      ct <- contraction[[i]]
      M2 <- 0 * ct$weight
      if (!is.null(ct$normalInverse)) M2[normal, normal] <- ct$normalInverse
      list(ct$weight[at$reach[i, ], at$reach[i, ], drop = FALSE], M2[at$reach[i, ], at$reach[i, ], drop = FALSE])
    }))
  for (i in seq_len(G)) {
    en <- E[i, ]
    term <- terms$groups[[i]]
    hh <- term$hessian
    an <- .intersect(active, names(term$gradient))
    modeWeight <- shiftWeight <- numeric(K)
    modeWeight[sparse] <- delta[i, ]
    shiftWeight[sparse] <- -drop(blocks$schur[[i]] %*% delta[i, ]) - blocks$sparseGradient[i, ]
    tr <- .groupTraces(at, term, curv$groups[[i]], contraction[[i]], sweeps[[i]], en, an, active)
    groupGradient[i, an] <- groupGradient[i, an] + term$gradient[an] + drop(modeWeight %*% tr$exactA) +
      tr$dMA[match(an, active)]
    etaGradient[en] <- term$gradient[en] + drop(modeWeight %*% tr$exactE) + tr$dME + shiftWeight
    if (length(normal)) {
      # d log|H_NN| along the parameters, the Cholesky entries and the etas; the prior precision in H_NN
      # also enters the Schur complement
      groupGradient[i, ] <- groupGradient[i, ] + tr$dHA
      etaGradient[en] <- etaGradient[en] + tr$dHE
      ch <- .intersect(active, nl$cholPars)
      MN <- contraction[[i]]$normalInverse + contraction[[i]]$weight[normal, normal, drop = FALSE]
      groupGradient[i, ch] <- groupGradient[i, ch] +
        vapply(ch, function(cp) 2 * sum(MN * blocks$cholesky$dOinv[[match(cp, nl$cholPars)]]), 0)
    }
    etaHessian[en, en] <- tr$exactE
    crossHessian[en, an] <- tr$exactA
    etaModel[en, en] <- if (at$forwardReverse) hh[en, en] else tr$exactE
    crossModel[en, an] <- if (at$forwardReverse) hh[en, an] else tr$exactA
    activeHessian[an, an] <- activeHessian[an, an] + hh[an, an]
  }
  explicit <- colSums(groupGradient)
  pr <- terms$prior
  if (!is.null(pr)) {
    pa <- .intersect(active, names(pr$gradient)); pe <- .intersect(etaOrder, names(pr$gradient))
    explicit[pa] <- explicit[pa] + pr$gradient[pa]
    etaGradient[pe] <- etaGradient[pe] + pr$gradient[pe]
    etaHessian[pe, pe] <- etaHessian[pe, pe] + pr$hessian[pe, pe]
    crossHessian[pe, pa] <- crossHessian[pe, pa] + pr$hessian[pe, pa]
    etaModel[pe, pe] <- etaModel[pe, pe] + pr$hessian[pe, pe]
    crossModel[pe, pa] <- crossModel[pe, pa] + pr$hessian[pe, pa]
    activeHessian[pa, pa] <- activeHessian[pa, pa] + pr$hessian[pa, pa]
  }
  freeGradient <- setNames(numeric(nActive), active)
  freeHessian <- matrix(0, nActive, nActive, dimnames = list(active, active))
  for (o in terms$free) {
    an <- .intersect(active, names(o$gradient))
    freeGradient[an] <- freeGradient[an] + o$gradient[an]
    freeHessian[an, an] <- freeHessian[an, an] + o$hessian[an, an]
  }
  list(groupGradient = groupGradient, explicit = explicit + freeGradient, etaGradient = etaGradient,
       etaHessian = etaHessian, crossHessian = crossHessian, etaModel = etaModel, crossModel = crossModel,
       activeHessian = activeHessian + freeHessian, freeGradient = freeGradient, freeHessian = freeHessian,
       etaOrder = etaOrder,
       onData = etaOrder[as.vector(t(at$reach))])
}

## The model Hessian: the Schur complement over all etas, the sparse etas with the effective prior curvature
## of the integral, the sites of EC or per factor 2 C^-1 - S of the posterior covariance C; NULL without active
## parameters or where the solve fails.
.sparseModelHessian <- function(at, blocks, integral, d) {
  posterior <- d$etaModel
  factors <- .sparseFactors(lapply(blocks$schur, `/`, 2), colnames(at$E)[at$sparseCols], at$lay)$factors
  for (i in seq_len(at$G)) {
    en <- at$E[i, at$sparseCols]
    if (isTRUE(integral$ec[i])) {
      posterior[en, en] <- posterior[en, en] + diag(2 * pmax(integral$sites$Lam[i, ], 0), length(at$sparseCols))
      next
    }
    V <- integral$cov[[i]]
    if (anyNA(V)) next
    for (f in factors) {
      e <- eigen(V[f, f, drop = FALSE], symmetric = TRUE)
      inv <- e$vectors %*% (t(e$vectors) / pmax(e$values, 1e-12 * max(1, e$values)))
      posterior[en[f], en[f]] <- posterior[en[f], en[f]] + 2 * inv - blocks$schur[[i]][f, f, drop = FALSE]
    }
  }
  on <- d$onData
  if (ncol(d$crossModel))
    tryCatch(d$activeHessian - crossprod(d$crossModel[on, , drop = FALSE],
                                         solve(posterior[on, on, drop = FALSE], d$crossModel[on, , drop = FALSE])),
             error = function(e) NULL)
}

## The change of the value to the exact modes, to first order: the value depends on the modes by the eta
## gradient, and the inner search leaves a gradient at them that a Newton step on the eta Hessian removes.
.modeCorrection <- function(at, terms, d, slabPrecision) {
  inner <- setNames(numeric(length(d$etaOrder)), d$etaOrder)
  for (i in seq_len(at$G)) inner[at$E[i, ]] <- terms$etaGradient[i, ]
  pr <- terms$prior
  if (!is.null(pr)) {
    pe <- .intersect(d$etaOrder, names(pr$gradient))
    inner[pe] <- inner[pe] + pr$gradient[pe]
  }
  inner[names(slabPrecision)] <- inner[names(slabPrecision)] + slabPrecision * at$full[names(slabPrecision)]
  on <- d$onData
  tryCatch(-sum(d$etaGradient[on] * solve(d$etaHessian[on, on, drop = FALSE], inner[on])), error = function(e) 0)
}

## The total gradient with the shift of the modes along the parameters, block diagonal by groups, the Schur
## complement of the active block, and with `scores` the gradient per group.
.sparseOuterDerivatives <- function(at, terms, blocks, integral, d, hyp, hypers, active, fixed, scores) {
  G <- at$G; K <- at$K
  nActive <- length(active)
  tauDeriv <- .tauDerivative(at$modes, at$E, at$lay, hyp)
  gradient <- d$explicit
  hyperGradient <- -2 * integral$dhyper
  byGroup <- if (scores) {
    hg <- integral$hyperGroups %||% matrix(NA_real_, G, length(hypers), dimnames = list(NULL, hypers))
    cbind(d$groupGradient, -2 * hg[, hypers, drop = FALSE])
  }
  if (scores && !is.null(terms$prior))
    byGroup[, active] <- byGroup[, active] + .priorGroups(at$normal, at$E, at$full, fixed, blocks$cholesky, active)
  modeShift <- matrix(0, length(d$etaOrder), nActive, dimnames = list(d$etaOrder, active))
  owner <- outer(seq_len(G), rep(seq_len(G), each = K), `==`) * 1
  schur <- d$activeHessian
  if (nActive) {
    shift <- solve(d$etaHessian, d$crossHessian)
    modeShift <- -shift
    gradient <- gradient - drop(d$etaGradient %*% shift)
    if (scores) byGroup[, active] <- byGroup[, active] - owner %*% (d$etaGradient * shift)
    schur <- schur - crossprod(d$crossHessian, shift)
  }
  for (h in hypers) {
    v <- d$etaGradient * drop(solve(d$etaHessian, ifelse(!is.na(tauDeriv$hyper) & tauDeriv$hyper == h, tauDeriv$d, 0)))
    hyperGradient[h] <- hyperGradient[h] - sum(v)
    if (scores) byGroup[, h] <- byGroup[, h] - drop(owner %*% v)
  }
  list(modeShift = modeShift, gradient = c(gradient, hyperGradient), schur = schur, byGroup = byGroup)
}

## The expected information in `active` of the terms without random effects, by which they enter the
## information; `freeHessian` gives the layout.
.freeInformation <- function(design, objs, terms, at, active, fixed, freeHessian) {
  info <- 0 * freeHessian
  extra <- terms$extra
  free <- c(lapply(design$freeConditions, function(cn)
    objs[[cn]](at$full, fixed = fixed, deriv = TRUE, hessian = "expected")),
    if (!is.null(extra)) list(extra(pars = at$full, fixed = fixed, deriv = TRUE, hessian = "expected")))
  for (o in free) {
    an <- .intersect(active, rownames(o$hessian))
    info[an, an] <- info[an, an] + o$hessian[an, an]
  }
  info
}

## The gradient of the normal prior of every group, eta' Omega^-1 eta + log|Omega|, in the
## Cholesky parameters among `A`, a row per group.
.priorGroups <- function(nl, E, full, fixed, cd, A) {
  ch <- .intersect(A, nl$cholPars)
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

## Per eta, in the order of t(E), the slab standard deviation the prior gradient of a sparse eta
## depends on, by name, and the derivative of that gradient in it at the modes.
.tauDerivative <- function(modes, E, lay, hyp) {
  names <- as.vector(t(E))
  out <- list(hyper = setNames(rep(NA_character_, length(names)), names), d = setNames(numeric(length(names)), names))
  for (k in which(colnames(E) %in% lay$stem)) {
    r <- match(colnames(E)[k], lay$stem)
    out$hyper[E[, k]] <- lay$tauName[r]
    out$d[E[, k]] <- -4 * modes[, k] / hyp$tau[[lay$stem[r]]]^3
  }
  out
}


## ---- fit -------------------------------------------------------------------------------------

## Sparse random effects: one trust region over all parameters, the hyperparameters unbounded, on
## the Schur complement under the posterior covariance of the sparse etas and per hyperparameter
## the curvature whose Newton step is the M-step, with a BFGS correction for the missing information.
.fitSparse <- function(design, start, fixed, control = list(), verbose = FALSE) {
  .requireM3(design, NULL, "FOCEI")
  lay <- .sparseLayout(design)
  objs <- .conditionObjectives(design)
  prior <- if (design$normal$K) .normalPrior(design)
  hypers <- .intersect(.sparseHypers(lay), names(start))
  off <- .intersect(attr(.extraObjective(design), "parameters"), c(.sparseHypers(lay), as.vector(design$etaNames)))
  if (length(off))
    stop("em: with spike-and-slab or Dirichlet-process random effects the terms beside the data ",
         "take the fixed effects only; ", paste(off, collapse = ", "), " are hyperparameters or random effects.",
         call. = FALSE)
  kind <- .hyperKind(lay)[hypers]
  structural <- .setdiff(names(start), hypers)
  E <- design$etaNames
  reach <- .reachOf(design)
  blind <- Filter(function(h) {
    st <- lay$stem[lay[[paste0(kind[[h]], "Name")]] %in% h]
    !any(reach[, st])
  }, hypers)
  if (length(blind))
    stop("em: the random effects of ", paste(blind, collapse = ", "), " reach the data of no group; ",
         "fix ", if (length(blind) > 1L) "them" else "it", " or drop the effect.", call. = FALSE)
  # the inner modes of a fit given as init where they fit the design, else 0
  warm <- matrix(0, nrow(E), ncol(E), dimnames = dimnames(E))
  given <- control$modes
  if (is.matrix(given) && all(dimnames(E)[[1]] %in% rownames(given)) && all(colnames(E) %in% colnames(given)))
    warm[] <- given[rownames(E), colnames(E)]
  warmAt <- shift <- NULL
  hyper <- start[hypers]
  par <- start[structural]
  lc <- control$sparseFit %||% list()
  units <- lapply(design$groups, function(g) list(eta = unname(E[g, reach[g, ]]), conditions = design$conditions[[g]]))
  lc$derivMode <- .secondOrderMode(design, design$objfn, par, fixed,
                                   control[.setdiff(names(control), "secondOrderCorrection")], units = units)
  lc$meta <- .foceiMeta(design, c(par, fixed), NULL, colnames(E), reach)
  maxOuter <- lc$maxOuter %||% 200L
  jointControl <- .trustControl(list(rinit = 1, rmax = 10, tolControl = list(ftol = 1e-10, mtol = 1e-5)),
                                lc$trustControl, label = "trustControl")
  # control$inner = "sqp": the modes by sqpIterlim inner iterations from the predicted modes, without the earlier
  # optima as starts, until the trust region stops at the tolerances of .sqpTolerance(); then on nested
  nested <- lc
  inexact <- identical(control[["inner"]], "sqp")
  if (inexact) lc$iterlim <- as.integer(control$innerControl$sqpIterlim %||% 2L)
  runControl <- if (inexact) .sqpTolerance(jointControl) else jointControl
  innerIterations <- 0L
  # the inner start at p: the last modes moved by their first-order shift along the parameters
  predicted <- function(p) {
    if (is.null(shift)) return(warm)
    q <- .intersect(colnames(shift), names(p))
    out <- warm
    at <- match(rownames(shift), as.vector(E))
    out[at] <- out[at] + drop(shift[, q, drop = FALSE] %*% (p[q] - warmAt[q]))
    out
  }
  # a trial point at which the integral does not converge is a failed evaluation
  sites <- optima <- NULL
  # every evaluation starts from the modes and sites of the last accepted point, so that the
  # objective near it does not depend on the trials before
  marginal <- function(p, grad = TRUE, commit = TRUE, scores = FALSE) {
    # the solves and sweeps of an evaluation leave freed memory behind that the process would keep
    on.exit(release_heap(), add = TRUE)
    r <- .sparseMarginal(design, objs, prior, c(p, fixed), .stemHypers(lay, c(hyper, fixed)), names(p),
                         hypers, fixed, predicted(p), grad = grad, scores = scores,
                         control = modifyList(lc, list(sites = sites, optima = if (!inexact) optima)))
    innerIterations <<- innerIterations + r$innerIterations
    if (isTRUE(r$infeasible)) return(r)
    if (commit) accept(p, r)
    r
  }
  # the inner minima from several starts, at the start of the fit and every 25 accepted steps
  explored <- NA_integer_
  scale <- matrix(1, nrow(E), ncol(E), dimnames = dimnames(E))
  explore <- function(p) {
    hy <- .stemHypers(lay, c(hyper, fixed))
    st <- lay$stem[lay$dist %in% c("spikeSlab", "dirichletProcess")]
    scale[, st] <- rep(hy$tau[st], each = nrow(E))
    found <- .exploreModes(design, c(p, fixed), hy, warm, nested, scale * .reachOf(design))
    innerIterations <<- innerIterations + attr(found, "iterations")
    optima <<- .modeOptima(optima, found)
    explored <<- steps()
  }
  accept <- function(p, r) {
    warm <<- r$modes
    warmAt <<- p
    if (!inexact) optima <<- .modeOptima(optima, r$modes)
    if (!is.na(explored) && steps() - explored >= 25L) explore(p)
    if (!is.null(r$modeShift)) shift <<- r$modeShift
    if (!is.null(r$sites)) sites <<- r$sites
  }
  rows <- list()
  record <- function(phase, value, dpar) {
    row <- data.frame(iter = length(rows) + 1L, phase = phase, OFV = value, deltaPar = dpar, t(hyper))
    rows[[length(rows) + 1L]] <<- row
  }
  # the iterations are the recorded steps; the start of every trust run is recorded without one
  steps <- function() sum(vapply(rows, function(r) !is.na(r$deltaPar), TRUE))
  # the objective in the structural parameters and the unbounded hyperparameters; the model Hessian
  # is the model Hessian of the marginal and the secant curvature of every hyperparameter to the
  # M-step, plus a BFGS correction from the points that lowered the value
  best <- Inf; bestAt <- bestGrad <- correction <- NULL
  # why evaluations of the trust region failed, and the warnings of the inner searches, by the
  # first line of the message
  failures <- integer(0)
  fail <- function(why) failures[why] <<- if (is.na(failures[why])) 1L else failures[why] + 1L
  jointObjective <- function(x, ...) {
    p <- x[names(par)]
    hyper <<- .hyperFromFree(kind, x[hypers])
    r <- withCallingHandlers(tryCatch(marginal(p, commit = FALSE), error = function(e) {
      fail(sub("\n.*", "", conditionMessage(e)))
      stop(e)
    }), warning = function(w) fail(paste("warning:", sub("\n.*", "", conditionMessage(w)))))
    if (isTRUE(r$infeasible)) {
      fail("the expectation consistent (EC) approximation did not converge on the sparse random effects")
      return(.failedEvaluation(x))
    }
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
      if (verbose) message(sprintf("em(sparse) iteration %d: -2 log L %.6f", length(rows), r$value))
      best <<- r$value; bestAt <<- x; bestGrad <<- g
    }
    objlist(value = r$value, gradient = g, hessian = base + correction)
  }
  # one run of the trust region from the current point: the trust result, NULL with a warning on
  # an error
  jointRun <- function(rinit = runControl$rinit) {
    x0 <- c(par, .hyperToFree(kind, hyper))
    best <<- Inf; bestAt <<- bestGrad <<- correction <<- NULL
    settings <- modifyList(runControl, list(objfun = jointObjective, parinit = x0, rinit = rinit,
                                            iterlim = max(maxOuter - steps(), 1L)))
    jt <- tryCatch(.innerTrust(settings), error = function(e) {
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

  if (length(par) + length(hyper) > 0L) explore(par)
  # a run that stops short starts again with a fresh correction while the runs lower -2 log L
  estimated <- length(par) + length(hyper) > 0L
  runs <- list(jt = list(), converged = !estimated, stopReason = if (estimated) "iterlim" else NA_character_)
  if (estimated) runs <- .sparseRuns(jointRun, steps, maxOuter, runControl$tolControl$mtol %||% 1e-5, runs)
  sqpIterations <- NULL
  if (inexact) {
    # the nested search from the estimates and modes of the inexact one
    inexact <- FALSE
    lc <- nested
    runControl <- jointControl
    if (estimated && !is.null(runs$jt)) {
      sqpIterations <- steps()
      runs <- .sparseRuns(jointRun, steps, maxOuter, jointControl$tolControl$mtol %||% 1e-5, runs)
    }
  }
  conv <- runs$converged
  stopReason <- runs$stopReason
  explored <- NA_integer_
  final <- marginal(par, scores = TRUE)
  if (isTRUE(final$infeasible) || !is.finite(final$value %||% NA_real_))
    stop("em: the integral over the sparse random effects has no finite value at the last point of the ",
         "fit, or the expectation consistent (EC) approximation did not converge there; a start nearer ",
         "the optimum may avoid it.", call. = FALSE)
  coupled <- final$coupling
  list(argument = c(par, hyper), value = final$value, gradient = final$gradient,
       hessian = NULL, information = .empiricalInformation(final, names(par)), etaModes = final$selected,
       etaCov = setNames(final$cov, design$groups),
       clusters = final$clusters, spike = final$spike, converged = conv, iterations = steps(),
       stopReason = stopReason, failures = failures, stageTrace = do.call(rbind, rows), coupling = coupled,
       innerModes = final$modes, innerIterations = innerIterations, sqpIterations = sqpIterations)
}

## Runs of the trust region by `run()` while each lowers -2 log L by more than `mtol`, up to `maxOuter` steps:
## the last run with convergence and stop reason, `previous` where none starts, `jt` NULL after an error.
.sparseRuns <- function(run, steps, maxOuter, mtol, previous) {
  out <- previous
  last <- Inf
  while (steps() < maxOuter) {
    jt <- run()
    if (is.null(jt)) return(list(jt = NULL, converged = FALSE, stopReason = "error"))
    out <- list(jt = jt, converged = isTRUE(jt$converged), stopReason = jt$stopReason %||% NA_character_)
    if (out$converged || !(jt$value < last - mtol)) break
    last <- jt$value
  }
  out
}

## The empirical information of a sparse fit from the scores of `final`: the model Hessian approximates the
## curvature of the integral by Gauss-Newton, so with independent groups it is G / (G - 1) times the sum of the
## outer products of the centered scores, plus half the Hessian of the conditions without random effects.
.empiricalInformation <- function(final, par) {
  if (is.null(final$scores) || anyNA(final$scores)) return(NULL)
  G <- nrow(final$scores)
  centered <- sweep(final$scores, 2, colMeans(final$scores)) / 2
  I <- if (G > 1L) G / (G - 1) * crossprod(centered) else 0 * crossprod(centered)
  I[par, par] <- I[par, par] + final$freeHessian / 2
  I
}


## The inner optima per group, `optima`, with the modes `modes` of an accepted point added where they
## differ from those held by more than 1e-3, the newest first, at most `keep` per group.
.modeOptima <- function(optima, modes, keep = 4L) {
  lapply(seq_len(nrow(modes)), function(i) {
    m <- modes[i, , drop = FALSE]
    known <- optima[[i]]
    if (is.null(known)) return(m)
    if (min(apply(abs(sweep(known, 2, m[1, ])), 1, max)) <= 1e-3) return(known)
    rbind(m, known)[seq_len(min(keep, nrow(known) + 1L)), , drop = FALSE]
  })
}

## A failed evaluation for trust(): it shrinks the radius, three in a row end the run.
.failedEvaluation <- function(x)
  objlist(value = Inf, gradient = setNames(0 * as.numeric(x), names(x)),
          hessian = matrix(0, length(x), length(x), dimnames = list(names(x), names(x))))

## The BFGS update of B = base + correction for the step s with gradient change y, returned as the
## new correction; a pair with 0 < s'y <= 1e-8 |s| |y| is skipped, a correction that a pair with
## s'y <= 0 contradicts or that leaves B indefinite is dropped.
.secantCorrection <- function(base, correction, s, y) {
  B <- base + correction
  Bs <- drop(B %*% s)
  sy <- sum(s * y); sBs <- sum(s * Bs)
  if (sy <= 0) return(0 * correction)
  if (sy <= 1e-8 * sqrt(sum(s^2) * sum(y^2)) || sBs <= 0) return(correction)
  updated <- B - tcrossprod(Bs) / sBs + tcrossprod(y) / sy
  if (min(eigen(updated, symmetric = TRUE, only.values = TRUE)$values) <= 0) return(0 * correction)
  updated - base
}

## The kind of every hyperparameter of the layout: prob, concentration or tau.
.hyperKind <- function(lay) {
  kinds <- c("prob", "concentration", "tau")
  out <- unlist(lapply(kinds, function(k) setNames(rep(k, nrow(lay)), lay[[paste0(k, "Name")]])))
  out[!is.na(names(out)) & !duplicated(names(out))]
}

## The hyperparameters on their unbounded scale, logit for a probability and log otherwise, back,
## and the derivative of a hyperparameter in the unbounded value.
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

## The complete-data information of every hyperparameter on the unbounded scale: the second
## derivative of the expected -2 log prior of the random effects at the posterior moments `r`,
## at least 1e-8 so that it seeds a positive definite model.
.hyperInformation <- function(lay, kind, r, hyper) {
  info <- vapply(names(hyper), function(h) {
    v <- hyper[[h]]
    k <- kind[[h]]
    st <- lay$stem[lay[[paste0(k, "Name")]] %in% h]
    S <- r$slab[st, , drop = FALSE]
    switch(k,
      tau = 4 * sum(S[, "slab2"]) / v^2,
      prob = 2 * sum(S[, "groups"]) * v * (1 - v),
      concentration = { i <- unlist(lapply(S[, "groups"], seq_len)); 2 * sum(v * i / (v + i)^2) })
  }, 0)
  pmax(info, 1e-8)
}

## The curvature of every hyperparameter on the unbounded scale whose Newton step is the M-step:
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
      tau = 2 * (sum(S[, "free"]) - sum(S[, "slab2"]) / v^2),
      prob = -2 * (sum(S[, "free"]) - sum(S[, "groups"]) * v),
      concentration = { i <- unlist(lapply(S[, "groups"], seq_len)); -2 * (sum(S[, "free"]) - sum(v / (v + i))) })
  }, 0)
  pmax(ifelse(abs(d) > 1e-8 & slope * d > 0, slope / d, info), 1e-8)
}

## M-steps from the posterior moments `r` over the coordinates that reach their data: tau^2 = E[sum of squared
## slab values] / E[slab blocks], prob = E[slab etas] / groups, and the concentration c with
## E[slab blocks] = sum_i c / (c + i) over the groups of every stem.
.updateSparse <- function(lay, r, current) {
  for (h in names(current)) {
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
## left side is increasing in c); `free` at most half a block below the limit, so c stays finite.
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


## ---- importance sampling on the model ------------------------------------------------

## -2 log L of a spike-and-slab fit by importance sampling on the model, `draws` per group from the
## mixture of the approximation with spreads widened by `widen`: value, standard error, effective
## sample size per group and the value of the approximation at the same point.
.sparseImportance <- function(object, draws = 1000L, widen = 1.5, chunk = 250L) {
  design <- object$design
  if (is.null(design))
    stop("logLik: the fit has no design; method = \"importance\" evaluates the model.", call. = FALSE)
  lay <- .sparseLayout(design)
  if (!is.null(design$groupings) || is.null(lay) || any(lay$dist != "spikeSlab"))
    stop("logLik: method = \"importance\" takes spike-and-slab random effects of one grouping, normal ones beside.",
         call. = FALSE)
  fixed <- object$fixed
  hypers <- .intersect(.sparseHypers(lay), names(object$argument))
  par <- object$argument[.setdiff(names(object$argument), hypers)]
  theta <- c(par, fixed)
  hyp <- .stemHypers(lay, c(object$argument, fixed))
  objs <- .conditionObjectives(design)
  prior <- if (design$normal$K) .normalPrior(design)
  E <- design$etaNames
  # the inner modes of the approximation at the estimates, of a fit without them the modes of the changes
  warm <- object$innerModes %||% object$etaModes %||% matrix(0, nrow(E), ncol(E), dimnames = dimnames(E))
  reach <- .reachOf(design)
  dm <- .secondOrderMode(design, design$objfn, par, fixed, list(), units = lapply(design$groups, function(g)
    list(eta = unname(E[g, reach[g, ]]), conditions = design$conditions[[g]])))
  r <- .sparseMarginal(design, objs, prior, theta, hyp, names(par), hypers, fixed, warm, grad = FALSE,
                       control = list(derivMode = dm))
  if (isTRUE(r$infeasible)) stop("logLik: the approximation has no value at the estimates.", call. = FALSE)
  G <- nrow(E)
  isL <- colnames(E) %in% lay$stem
  full <- c(theta, setNames(as.vector(r$modes), as.vector(E)))
  if (design$normal$K) {
    Lc <- .choleskyFactor(design$normal, theta)
    Omega <- Lc %*% t(Lc)
  }
  dmvn <- function(X, m, S) {
    R <- chol(S)
    Z <- backsolve(R, t(X) - m, transpose = TRUE)
    -colSums(Z^2) / 2 - sum(log(diag(R))) - nrow(R) * log(2 * pi) / 2
  }
  groups <- lapply(seq_len(G), function(i) {
    on <- r$reach[i, ]
    jL <- which(isL & on); jN <- which(!isL & on)
    X <- matrix(0, draws, ncol(E))
    lq <- lp <- numeric(draws)
    # the subset by the spike probabilities, the slab values of it jointly from the local model
    # exp(-eta' P eta / 2 + b' eta) of the data times the slab prior
    if (length(jL)) {
      kL <- match(jL, which(isL))
      st <- colnames(E)[jL]
      tau <- hyp$tau[st]; prob <- hyp$prob[st]
      p0 <- pmin(pmax(r$spike[i, jL], 0.01), 0.99)
      P <- r$local$P[[i]][kL, kL, drop = FALSE]
      bb <- r$local$b[i, kL]
      Z <- matrix(stats::runif(draws * length(jL)) > rep(p0, each = draws), draws)
      for (k in seq_len(draws)) {
        z <- Z[k, ]
        lq[k] <- sum(ifelse(z, log(1 - p0), log(p0)))
        lp[k] <- sum(ifelse(z, log(prob), log(1 - prob)))
        if (!any(z)) next
        Lam <- P[z, z, drop = FALSE] + diag(1 / tau[z]^2, sum(z))
        R <- chol(Lam / widen^2)
        mu <- backsolve(chol(Lam), backsolve(chol(Lam), bb[z], transpose = TRUE))
        x <- mu + backsolve(R, stats::rnorm(sum(z)))
        X[k, jL[z]] <- x
        u <- R %*% (x - mu)
        lq[k] <- lq[k] + sum(log(diag(R))) - sum(z) * log(2 * pi) / 2 - sum(u^2) / 2
        lp[k] <- lp[k] + sum(stats::dnorm(x, 0, tau[z], log = TRUE))
      }
    }
    if (length(jN)) {
      S <- widen^2 * r$cov[[i]][jN, jN, drop = FALSE]
      m <- r$modes[i, jN]
      Xn <- matrix(stats::rnorm(draws * length(jN)), draws) %*% chol(S) + rep(m, each = draws)
      X[, jN] <- Xn
      lq <- lq + dmvn(Xn, m, S)
      on <- match(colnames(E)[jN], design$normal$stems)
      lp <- lp + dmvn(Xn, rep(0, length(jN)), Omega[on, on, drop = FALSE])
    }
    list(X = X, lw = lp - lq)
  })
  # the data -2 log L of every group at the draws, in chunks of draws over all groups
  conds <- lapply(seq_len(G), function(i) design$conditions[[design$groups[i]]])
  data <- matrix(0, G, draws)
  for (start in seq(1L, draws, by = chunk)) {
    ks <- start:min(draws, start + chunk - 1L)
    jobs <- do.call(rbind, lapply(seq_len(G), function(i)
      expand.grid(group = i, draw = ks, cond = conds[[i]], stringsAsFactors = FALSE)))
    pars <- lapply(seq_len(nrow(jobs)), function(q) {
      p <- full
      p[E[jobs$group[q], ]] <- groups[[jobs$group[q]]]$X[jobs$draw[q], ]
      p
    })
    v <- vapply(.objEvalMany(objs[jobs$cond], pars, deriv = FALSE), `[[`, 0, "value")
    agg <- tapply(v, list(jobs$group, jobs$draw), sum)
    data[as.integer(rownames(agg)), as.integer(colnames(agg))] <- agg
  }
  terms <- vapply(seq_len(G), function(i) {
    lw <- groups[[i]]$lw - data[i, ] / 2
    mx <- max(lw)
    w <- exp(lw - mx)
    c(value = -2 * (mx + log(mean(w))), se = 2 * stats::sd(w) / (sqrt(draws) * mean(w)),
      ess = sum(w)^2 / sum(w^2))
  }, numeric(3))
  free <- if (length(design$freeConditions))
    sum(vapply(.objEvalMany(objs[design$freeConditions], rep(list(full), length(design$freeConditions)),
                            deriv = FALSE), `[[`, 0, "value"))
  extra <- .extraObjective(design)
  beside <- (free %||% 0) + if (!is.null(extra)) extra(pars = full, fixed = fixed, deriv = FALSE)$value else 0
  list(value = sum(terms["value", ]) + beside, se = sqrt(sum(terms["se", ]^2)),
       ess = setNames(terms["ess", ], rownames(E)), approximation = r$value)
}
