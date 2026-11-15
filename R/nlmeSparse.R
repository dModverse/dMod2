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

## The mode of the etas at `theta`, group by group, since the etas of a group enter only the group
## conditions: data -2 log L, normal prior and the slab alone, eta^2 / tau^2, per spike-and-slab
## eta, by trust. The model is differentiated in the etas of the group that reach the group data
## only, the others rest at 0; the Hessian is Gauss-Newton with a BFGS correction from the points
## that lowered the value, which adds the curvature of the residuals. A group starts from the lowest of `start` and the modes in
## `control$optima`, the optima found at earlier points, since the inner problem may have several.
## The tolerances of the outer trust region, which bound the error the solver leaves, end the search.
## With `control$phi` the value of every group at the group mode is attribute "phi".
.sparseModes <- function(design, objs, prior, theta, hyp, start, control = list()) {
  lay <- .sparseLayout(design)
  E <- design$etaNames
  names <- as.vector(E)
  slab <- .slabPrecision(E, lay, hyp)
  modes <- matrix(as.vector(start), nrow(E), dimnames = dimnames(E))
  reach <- .reachOf(design)
  modes[!reach] <- 0
  given <- as.list(control$trustControl$tolControl)
  tol <- utils::modifyList(list(ftol = 1e-10, mtol = 1e-10),
                           given[intersect(c("ftol", "mtol", "gtol"), names(given))])
  phis <- setNames(rep(NA_real_, nrow(E)), rownames(E))
  for (i in seq_len(nrow(E))) {
    en <- E[i, reach[i, ]]
    if (!length(en)) next
    conds <- design$conditions[[i]]
    sl <- slab[intersect(names(slab), en)]
    # the value alone
    phi <- function(eta) {
      full <- c(theta, setNames(as.vector(modes), names))
      full[en] <- as.numeric(eta)
      v <- sum(vapply(.objEvalMany(objs[conds], rep(list(full[en]), length(conds)), deriv = FALSE,
                                   fixed = full[setdiff(names(full), en)]), `[[`, 0, "value"))
      if (!is.null(prior)) v <- v + prior(pars = full, deriv = FALSE)$value
      if (length(sl)) v <- v + sum(sl * full[names(sl)]^2) / 2
      v
    }
    known <- control$optima[[i]]
    if (length(known)) {
      starts <- rbind(modes[i, reach[i, ]], known[, reach[i, ], drop = FALSE])
      modes[i, reach[i, ]] <- starts[which.min(vapply(seq_len(nrow(starts)), function(r) phi(starts[r, ]), 0)), ]
    }
    best <- Inf; bestAt <- bestGrad <- correction <- NULL
    objective <- function(eta, ...) {
      full <- c(theta, setNames(as.vector(modes), names))
      full[en] <- as.numeric(eta)
      out <- Reduce(`+`, .objEvalMany(objs[conds], rep(list(full[en]), length(conds)),
                                      fixed = full[setdiff(names(full), en)]))
      v <- out$value
      gr <- out$gradient[en]
      he <- out$hessian[en, en, drop = FALSE]
      if (!is.null(prior)) {
        pr <- prior(pars = full, deriv = TRUE)
        nn <- intersect(en, names(pr$gradient))
        v <- v + pr$value
        gr[nn] <- gr[nn] + pr$gradient[nn]
        he[nn, nn] <- he[nn, nn] + pr$hessian[nn, nn]
      }
      if (length(sl)) {
        es <- full[names(sl)]
        v <- v + sum(sl * es^2) / 2
        gr[names(sl)] <- gr[names(sl)] + sl * es
        he[cbind(names(sl), names(sl))] <- he[cbind(names(sl), names(sl))] + sl
      }
      if (is.null(correction)) correction <<- 0 * he
      # a correction from an earlier point that leaves the model at this one indefinite is dropped
      if (any(correction != 0) && min(eigen(he + correction, symmetric = TRUE, only.values = TRUE)$values) <= 0)
        correction <<- 0 * he
      if (v < best) {
        if (!is.null(bestAt)) correction <<- .secantCorrection(he, correction, eta - bestAt, gr - bestGrad)
        best <<- v; bestAt <<- eta; bestGrad <<- gr
      }
      list(value = v, gradient = gr, hessian = he + correction)
    }
    args <- list(objfun = objective, parinit = setNames(modes[i, reach[i, ]], en), rinit = 1, rmax = 10,
                 iterlim = control$iterlim %||% 100L, tolControl = tol)
    fit <- do.call(trust, c(args, list(sparse = FALSE)))
    modes[i, reach[i, ]] <- fit$argument[en]
    if (isTRUE(control$phi)) phis[i] <- phi(fit$argument[en])
  }
  if (isTRUE(control$phi)) attr(modes, "phi") <- phis
  modes
}

## The lowest inner minimum per group at `theta` from the modes `warm`, from 0 and from `spread`
## starts around `warm`, each eta moved by `scale` times a deterministic pseudo-normal number: the
## inner problem of a group may have several minima, and every search keeps to the one it starts in.
.exploreModes <- function(design, objs, prior, theta, hyp, warm, control, scale, spread = 2L) {
  ctl <- utils::modifyList(control, list(phi = TRUE, optima = NULL))
  best <- .sparseModes(design, objs, prior, theta, hyp, warm, ctl)
  low <- attr(best, "phi")
  n <- length(warm)
  for (k in 0:spread) {
    z <- if (k == 0L) 0 else stats::qnorm(pmin(pmax(((sin(seq_len(n) * k * 12.9898 + 78.233) * 43758.5453) %% 1), 0.01), 0.99))
    start <- if (k == 0L) 0 * warm else warm + scale * matrix(z, nrow(warm))
    m <- .sparseModes(design, objs, prior, theta, hyp, start, ctl)
    lower <- !is.na(attr(m, "phi")) & (is.na(low) | attr(m, "phi") < low - 1e-6)
    best[lower, ] <- m[lower, ]
    low[lower] <- attr(m, "phi")[lower]
  }
  attr(best, "phi") <- NULL
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
  factors <- c(as.list(setdiff(seq_along(stems), slabs)), if (length(slabs)) list(slabs))
  owner <- integer(length(stems))
  for (j in seq_along(factors)) owner[factors[[j]]] <- j
  within <- outer(owner, owner, `==`)
  list(factors = factors, rows = r, within = within,
       coupled = any(vapply(P, function(Pi) any(Pi[!within] != 0), TRUE)))
}

## The integral over the sparse coordinates, I = int exp(-eta' P_i eta / 2 + b_i' eta) prod p(eta)
## per group (src/sparse_kernel.cpp), group by group where no Dirichlet-process stem joins groups:
## exact where P couples no two factors and the joint spike-and-slab sum of a group stays within
## `budget` nodes, else by the expectation consistent (EC) approximation from the sites `sites`
## where given (`ec` marks those groups). Only the coordinates of `reach` whose row of P and entry
## of b are not all 0 enter, the others keep their prior. NULL where the integral does not
## converge. log I, the moments E eta and Cov eta that the derivatives of log I in b and P are, the
## derivatives in the hyperparameters `hypers`, also per group, the reported values, clusters and
## spike probabilities and the M-step statistics.
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
       values = none, clusters = none, spike = none,
       slab = matrix(0, nrow(lay), 3, dimnames = list(lay$stem, c("free", "groups", "slab2"))))
}

## ---- marginal ------------------------------------------------------------------------------

## The approximate marginal -2 log L at `theta`: the normal etas on the quadratic model of the joint
## mode, the sparse coordinates on the Schur complement in natural parameters, P = S / 2 and
## b = (S mode - g) / 2, integrated by .sparseIntegral(); the exact gradient in `active` and
## `hypers`, a model Hessian, the moments and the selected clusters. An eta that does not reach
## the data of the group rests at 0 with the prior. With
## `scores`, the gradient split into the rows of the groups and that of the conditions without
## random effects, with their Hessian.
.sparseMarginal <- function(design, objs, prior, theta, hyp, active, hypers, fixed, warm,
                            grad = TRUE, control = list(), scores = FALSE) {
  lay <- .sparseLayout(design)
  nl <- design$normal
  E <- design$etaNames
  G <- nrow(E); K <- ncol(E)
  isL <- colnames(E) %in% lay$stem
  iN <- which(!isL); iL <- which(isL)
  reach <- .reachOf(design)
  modes <- .sparseModes(design, objs, prior, theta, hyp, warm, control)
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
  units <- lapply(seq_len(G), function(i)
    list(eta = unname(E[i, reach[i, ]]), conditions = design$conditions[[design$groups[i]]]))
  if (fr) {
    fw <- .sweepForward(design, full, fixed, units)
    # the curvature of a group on all etas of the group, 0 in those that do not reach the group data
    cu <- .unitCurvature(design, fw, units)
    curv <- lapply(seq_len(G), function(i) {
      H <- matrix(0, K, K)
      H[reach[i, ], reach[i, ]] <- cu[[i]]
      list(H = H)
    })
  } else {
    conds <- unlist(design$conditions, use.names = FALSE)
    times <- lapply(conds, function(cn) .emGrid(design, cn))
    preds <- setNames(.predictMany(design$prdfn, times = times, parsList = rep(list(full), length(conds)),
                                   conditions = conds, fixed = fixed, deriv = TRUE, deriv2 = TRUE), conds)
    curv <- lapply(seq_len(G), function(i) .gnCurvature(design, design$groups[i], preds, E[i, ], A, fixed))
  }
  gr <- matrix(t(vapply(seq_len(G), function(i) at[[i]]$gradient[E[i, ]], numeric(K))), G,
               dimnames = dimnames(E))
  gr[!reach] <- 0

  # the normal etas: H_NN of every group with the prior precision 2 Omega^-1; the sparse block
  # of `Hs` is the Schur complement of `Hs`, the curvature left once the normal etas are integrated out
  cd <- if (length(iN)) .choleskyDerivatives(nl, .choleskyFactor(nl, c(full, fixed)))
  HN <- lapply(curv, function(cv) if (length(iN)) cv$H[iN, iN, drop = FALSE] + 2 * cd$Oinv)
  volume <- sum(vapply(HN, function(H) if (is.null(H)) 0 else as.numeric(determinant(H / 2)$modulus), 0))
  # the sparse block of every group: the curvature S, the Schur complement of the normal block;
  # the data term on it is exp(-eta' P eta / 2 + b' eta) up to the constant mode' S mode / 2 - g' mode
  S <- lapply(seq_len(G), function(i) {
    H <- curv[[i]]$H
    s <- H[iL, iL, drop = FALSE]
    if (length(iN)) s <- s - H[iL, iN, drop = FALSE] %*% solve(HN[[i]], H[iN, iL, drop = FALSE])
    (s + t(s)) / 2
  })
  KL <- length(iL)
  bySparse <- function(f) matrix(vapply(seq_len(G), f, numeric(KL)), G, KL, byrow = TRUE,
                                 dimnames = list(rownames(E), colnames(E)[iL]))
  gL <- gr[, iL, drop = FALSE]
  mL <- modes[, iL, drop = FALSE]
  bL <- bySparse(function(i) drop(S[[i]] %*% mL[i, ] - gL[i, ]) / 2)
  offset <- sum(vapply(seq_len(G), function(i) sum(mL[i, ] * (S[[i]] %*% mL[i, ])) / 2 - sum(gL[i, ] * mL[i, ]), 0))
  importance <- identical(control$integral, "importance")
  sp <- if (importance) {
    # in center form, m = mode - S^-1 g, log I of the natural form is that of the center form plus m' P m / 2
    mc <- mL - bySparse(function(i) {
      a <- reach[i, iL]
      out <- numeric(KL)
      out[a] <- solve(S[[i]][a, a, drop = FALSE], gL[i, a])
      out
    })
    is <- .importanceField(lapply(S, `/`, 2), mc, lay, hyp, hypers, control$particles %||% 500L)
    is$logI <- is$logI + sum(vapply(seq_len(G), function(i) sum(mc[i, ] * (S[[i]] %*% mc[i, ])) / 4, 0))
    is
  } else .sparseIntegral(lapply(S, `/`, 2), bL, lay, hyp, hypers, control$sites, reach[, iL, drop = FALSE])
  if (is.null(sp)) return(list(infeasible = TRUE))
  value <- sum(vapply(at, `[[`, 0, "value")) + (pr$value %||% 0) + volume - 2 * sp$logI + offset +
    sum(vapply(free, `[[`, 0, "value"))
  selected <- modes
  selected[, iL] <- ifelse(is.na(sp$values), modes[, iL], sp$values)
  clusters <- spike <- matrix(NA_real_, G, K, dimnames = dimnames(E))
  clusters[, iL] <- sp$clusters; spike[, iL] <- sp$spike
  mean <- matrix(0, G, K, dimnames = dimnames(E))
  mean[, iL] <- sp$mean
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
  out <- list(value = value, modes = modes, mean = mean, var = var, cov = cov, slab = sp$slab,
              selected = selected, clusters = clusters, spike = spike, valueSE = 2 * (sp$se %||% 0),
              coupling = sp$coupling, sites = sp$sites, reach = reach)
  if (!grad || importance) return(out)

  # derivatives of the value in g, S and the mode on the sparse block, from the moments of the
  # integral at the fixed point, delta = E eta - mode and Cov eta, then through the Schur
  # complement on the cross and normal blocks of H
  delta <- sp$mean - mL
  Q <- length(A)
  own <- matrix(0, G, Q, dimnames = list(rownames(E), A))
  names <- as.vector(t(E))
  dEta <- setNames(numeric(G * K), names)
  Hex <- matrix(0, G * K, G * K, dimnames = list(names, names))
  Hcross <- matrix(0, G * K, Q, dimnames = list(names, A))
  # the blocks of the model Hessian: under forward over reverse those of the Gauss-Newton Hessian, as
  # the parameter block is, so that the Schur complement is positive semidefinite; the exact rows
  # enter the gradient and the mode shift
  HexM <- Hex; HcrossM <- Hcross
  schur <- matrix(0, Q, Q, dimnames = list(A, A))
  # per group the matrices the derivatives of the curvature are contracted with: MHi on all etas,
  # H_NN^-1 on the normal ones
  contraction <- lapply(seq_len(G), function(i) {
    MHi <- matrix(0, K, K, dimnames = list(E[i, ], E[i, ]))
    MHi[iL, iL] <- (tcrossprod(delta[i, ]) + sp$cov[[i]]) / 2
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
  sw <- if (fr) .sweepUnits(design, full, fixed, units, fw = fw, Ms = lapply(seq_len(G), function(i) {
    ct <- contraction[[i]]
    M2 <- 0 * ct$MHi
    if (!is.null(ct$Hi)) M2[iN, iN] <- ct$Hi
    list(ct$MHi[reach[i, ], reach[i, ], drop = FALSE], M2[reach[i, ], reach[i, ], drop = FALSE])
  }))
  pick <- function(v, nm) { out <- setNames(numeric(length(nm)), nm); h <- intersect(nm, names(v)); out[h] <- v[h]; out }
  for (i in seq_len(G)) {
    en <- E[i, ]
    o <- at[[i]]
    hh <- o$hessian
    an <- intersect(A, names(o$gradient))
    cv <- curv[[i]]
    Mgi <- Meta <- numeric(K)
    Mgi[iL] <- delta[i, ]
    Meta[iL] <- -drop(S[[i]] %*% delta[i, ]) - gL[i, ]
    MHi <- contraction[[i]]$MHi
    Hi <- contraction[[i]]$Hi
    # the exact rows of the group's data Hessian at the group etas, and tr(M dH) along A and the etas
    if (fr) {
      rowsAt <- function(cols) {
        out <- matrix(0, K, length(cols), dimnames = list(en, cols))
        h <- intersect(cols, colnames(sw[[i]]$rows))
        out[rownames(sw[[i]]$rows), h] <- sw[[i]]$rows[, h]
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
    HexM[en, en] <- if (fr) hh[en, en] else exactE
    HcrossM[en, an] <- if (fr) hh[en, an] else exactA
    schur[an, an] <- schur[an, an] + hh[an, an]
  }
  explicit <- colSums(own)
  if (!is.null(pr)) {
    pa <- intersect(A, names(pr$gradient)); pe <- intersect(names, names(pr$gradient))
    explicit[pa] <- explicit[pa] + pr$gradient[pa]
    dEta[pe] <- dEta[pe] + pr$gradient[pe]
    Hex[pe, pe] <- Hex[pe, pe] + pr$hessian[pe, pe]
    Hcross[pe, pa] <- Hcross[pe, pa] + pr$hessian[pe, pa]
    HexM[pe, pe] <- HexM[pe, pe] + pr$hessian[pe, pe]
    HcrossM[pe, pa] <- HcrossM[pe, pa] + pr$hessian[pe, pa]
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
  # curvature of the integral, the sites of EC or per factor 2 C^-1 - S of the
  # posterior covariance C; with normal etas alone it is the Gauss-Newton Schur complement
  posterior <- HexM
  factors <- .sparseFactors(lapply(S, `/`, 2), colnames(E)[iL], lay)$factors
  for (i in seq_len(G)) {
    en <- E[i, iL]
    if (isTRUE(sp$ec[i])) {
      posterior[en, en] <- posterior[en, en] + diag(2 * pmax(sp$sites$Lam[i, ], 0), KL)
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
  on <- names[as.vector(t(reach))]
  model <- if (Q) tryCatch(schur - crossprod(HcrossM[on, , drop = FALSE],
                                             solve(posterior[on, on, drop = FALSE], HcrossM[on, , drop = FALSE])),
                           error = function(e) NULL)

  slabPrec <- .slabPrecision(E, lay, hyp)
  Hex[cbind(names(slabPrec), names(slabPrec))] <- Hex[cbind(names(slabPrec), names(slabPrec))] + slabPrec

  # the shift of the mode: every eta moves with the parameters; the curvature of the mode is block
  # diagonal by groups, so the shift of a group is separate
  td <- .tauDerivative(modes, E, lay, hyp)
  gradient <- explicit
  hyperGrad <- -2 * sp$dhyper
  byGroup <- if (scores) {
    hg <- sp$hyperGroups %||% matrix(NA_real_, G, length(hypers), dimnames = list(NULL, hypers))
    cbind(own, -2 * hg[, hypers, drop = FALSE])
  }
  if (scores && !is.null(pr)) byGroup[, A] <- byGroup[, A] + .priorGroups(nl, E, full, fixed, cd, A)
  out$modeShift <- matrix(0, length(names), Q, dimnames = list(names, A))
  owner <- outer(seq_len(G), rep(seq_len(G), each = K), `==`) * 1
  if (Q) {
    shift <- solve(Hex, Hcross)
    out$modeShift <- -shift
    gradient <- gradient - drop(dEta %*% shift)
    if (scores) byGroup[, A] <- byGroup[, A] - owner %*% (dEta * shift)
    schur <- schur - crossprod(Hcross, shift)
  }
  for (h in hypers) {
    v <- dEta * drop(solve(Hex, ifelse(!is.na(td$hyper) & td$hyper == h, td$d, 0)))
    hyperGrad[h] <- hyperGrad[h] - sum(v)
    if (scores) byGroup[, h] <- byGroup[, h] - drop(owner %*% v)
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
  hypers <- intersect(.sparseHypers(lay), names(start))
  off <- intersect(attr(.extraObjective(design), "parameters"), c(.sparseHypers(lay), as.vector(design$etaNames)))
  if (length(off))
    stop("em: with spike-and-slab or Dirichlet-process random effects the terms beside the data ",
         "take the fixed effects only; ", paste(off, collapse = ", "), " are hyperparameters or random effects.",
         call. = FALSE)
  kind <- .hyperKind(lay)[hypers]
  structural <- setdiff(names(start), hypers)
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
  lc$derivMode <- .secondOrderMode(design, design$objfn, par, fixed, control[setdiff(names(control), "secondOrderCorrection")],
                                   units = lapply(design$groups, function(g)
                                     list(eta = unname(E[g, reach[g, ]]), conditions = design$conditions[[g]])))
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
  # a trial point at which the integral does not converge is a failed evaluation
  sites <- optima <- NULL
  # every evaluation starts from the modes and sites of the last accepted point, so that the
  # objective near it does not depend on the trials before
  marginal <- function(p, grad = TRUE, commit = TRUE, scores = FALSE) {
    # the solves and sweeps of an evaluation leave freed memory behind that the process would keep
    on.exit(release_heap(), add = TRUE)
    r <- .sparseMarginal(design, objs, prior, c(p, fixed), .stemHypers(lay, c(hyper, fixed)), names(p),
                         hypers, fixed, predicted(p), grad = grad, scores = scores,
                         control = modifyList(lc, list(sites = sites, optima = optima)))
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
    found <- .exploreModes(design, objs, prior, c(p, fixed), hy, warm, lc, scale * .reachOf(design))
    optima <<- .modeOptima(optima, found)
    explored <<- steps()
  }
  accept <- function(p, r) {
    warm <<- r$modes
    warmAt <<- p
    optima <<- .modeOptima(optima, r$modes)
    if (!is.na(explored) && steps() - explored >= 25L) explore(p)
    if (!is.null(r$modeShift)) shift <<- r$modeShift
    if (!is.null(r$sites)) sites <<- r$sites
  }
  rows <- list()
  record <- function(phase, value, dpar)
    rows[[length(rows) + 1L]] <<- data.frame(iter = length(rows) + 1L, phase = phase, OFV = value,
                                             deltaPar = dpar, t(hyper))
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

  if (length(par) + length(hyper) > 0L) explore(par)
  # a run that stops short starts again with a fresh correction while the runs lower -2 log L
  estimated <- length(par) + length(hyper) > 0L
  conv <- !estimated
  stopReason <- if (estimated) "iterlim" else NA_character_
  last <- Inf
  mtol <- jointControl$tolControl$mtol %||% 1e-5
  while (estimated && steps() < maxOuter) {
    jt <- jointRun()
    if (is.null(jt)) { conv <- FALSE; stopReason <- "error"; break }
    conv <- isTRUE(jt$converged)
    stopReason <- jt$stopReason %||% NA_character_
    if (conv || !(jt$value < last - mtol)) break
    last <- jt$value
  }
  explored <- NA_integer_
  final <- marginal(par, scores = TRUE)
  if (isTRUE(final$infeasible) || !is.finite(final$value %||% NA_real_))
    stop("em: the integral over the sparse random effects has no finite value at the last point of the ",
         "fit, or the expectation consistent (EC) approximation did not converge there; a start nearer ",
         "the optimum may avoid it.",
         call. = FALSE)
  coupled <- final$coupling
  # the model Hessian approximates the curvature of the integral by Gauss-Newton, so the fit has none; with
  # independent groups the empirical information is G / (G - 1) times the sum of the outer products
  # of their centered scores of log L, the conditions without random effects add half their Hessian
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
       stopReason = stopReason, failures = failures, stageTrace = do.call(rbind, rows), coupling = coupled,
       innerModes = final$modes)
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

## M-steps from the posterior moments `r` over the coordinates that reach their data:
## tau^2 = E[sum of squared slab values] / E[slab blocks],
## prob = E[slab etas] / groups, and the concentration c with E[slab blocks] = sum_i c / (c + i)
## over the groups of every stem.
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
