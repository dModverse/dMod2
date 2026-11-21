# Adaptive sparse-grid quadrature of the marginal likelihood of every group around the
# conditional mode of the normal random effects of the group, on nested Genz-Keister grids, and the ECM fits
# on it: one trust over every parameter on frozen nodes, the nodes rebuilt between the steps.


# Smolyak grids depend on the dimension and level only.
.sparseGridCache <- new.env(parent = emptyenv())

## The Smolyak grid of nested Genz-Keister rules; `key` identifies a node across levels.
.sparseGrid <- function(K, level) {
  key <- paste0(K, "_", level)
  if (is.null(.sparseGridCache[[key]])) {
    raw <- sparseGridGK(as.integer(K), as.integer(level))
    .sparseGridCache[[key]] <- list(z = raw$nodes, logAbsW = log(abs(raw$weights)),
                                    signs = sign(raw$weights), z2 = rowSums(raw$nodes^2),
                                    key = apply(round(raw$nodes * 1e12), 1, paste, collapse = " "))
  }
  .sparseGridCache[[key]]
}

## The next level above `level` whose grid differs, NA above the largest one.
.nextLevel <- function(K, level) {
  a <- .sparseGrid(K, level)
  for (l in seq(level + 1L, length.out = 25L)) {
    if (l - K + 1L > 26L) return(NA_integer_)
    b <- .sparseGrid(K, l)
    if (!identical(a$key, b$key) || !identical(a$logAbsW, b$logAbsW)) return(l)
  }
  NA_integer_
}

## The highest level below `level`, at least K, whose grid differs, NA if there is none.
.prevLevel <- function(K, level) {
  a <- .sparseGrid(K, level)
  for (l in rev(seq_len(level - K) + K - 1L)) {
    b <- .sparseGrid(K, l)
    if (!identical(a$key, b$key) || !identical(a$logAbsW, b$logAbsW)) return(l)
  }
  NA_integer_
}

## The node count of a level after pruning.
.levelNodes <- function(K, level, pruneTol = Inf) {
  w <- .sparseGrid(K, level)$logAbsW
  sum(w >= max(w) - pruneTol)
}

## Nodes eta = mode + sqrt(2) R^-1 z and log weights w exp(z'z) 2^(K/2) / |R| around the mode
## with posterior precision P = R'R, exact for a Gaussian posterior; nodes whose |w| lies
## `pruneTol` in log below the largest are dropped.
.quadratureNodes <- function(mode, precision, level, pruneTol = Inf) {
  K <- length(mode)
  grid <- .sparseGrid(K, level)
  R <- chol(precision)
  eta <- t(sqrt(2) * backsolve(R, t(grid$z))) + matrix(mode, nrow(grid$z), K, byrow = TRUE)
  colnames(eta) <- names(mode)
  logW <- (K / 2) * log(2) - sum(log(diag(R))) + grid$logAbsW + grid$z2
  keep <- grid$logAbsW >= max(grid$logAbsW) - pruneTol
  list(eta = eta[keep, , drop = FALSE], logW = logW[keep], signs = grid$signs[keep], key = grid$key[keep])
}


## Conditional modes and inner Hessians H_g (of the -2 log density, so the posterior precision
## is H_g / 2) of every group at `psi`, by the inner search of the group kernel.
.innerModes <- function(design, psi, etaInit, fixed, meta, control = list()) {
  lay <- design$normal
  E <- design$etaNames[, lay$stems, drop = FALSE]
  full <- setNames(numeric(length(meta$pars_full_names)), meta$pars_full_names)
  full[meta$outer_names] <- psi[meta$outer_names]
  full[as.vector(E)] <- as.vector(etaInit)
  L <- .choleskyFactor(lay, psi)
  ic <- modifyList(list(rinit = 1, rmax = 10, iterlim = 30, steptol = 1e-8, mtol = 1e-12,
                        eigen_floor_relative = 1e-10), control)
  fr <- nlme_group_kernel(model_cb = design$prdfn, err_cb = .errorModel(design), pars_full = full,
                          eta_warmstart = etaInit, subject_meta = meta,
                          Omega_inv_mat = chol2inv(t(L)), Omega_log_det = 2 * sum(log(diag(L))),
                          fixed = fixed, control = ic, predict_cb = .guardPredictor(.batchPredictor(design)))
  modes <- fr$eta_modes
  dimnames(modes) <- dimnames(etaInit)
  list(modes = modes, H = setNames(fr$H_GN, design$groups))
}


## Per group the matrix D_g that maps a change of the mu-referenced fixed effects to the shift
## of the group etas, eta = phi - D_g beta with phi held, one column per referenced parameter.
.referenceShift <- function(design) {
  lay <- design$normal
  # the conditions without a group are evaluated at psi, so referencing reads the grouped ones
  grouped <- design
  grouped$data <- design$data[unlist(design$conditions, use.names = FALSE)]
  ref <- .muReference(grouped)
  pars <- unlist(lapply(unname(ref), `[[`, "pars"))
  D <- lapply(setNames(nm = design$groups), function(g) {
    m <- matrix(0, lay$K, length(pars), dimnames = list(lay$stems, pars))
    for (e in ref) m[e$stem, e$pars] <- e$X[g, ]
    m
  })
  list(pars = pars, D = D)
}


## -2 log p(y_g | eta) at the nodes of several groups in one batch: `etas` holds per group the
## node etas, `base` every other parameter. Per group the values and, with `grad`, gradient and
## Hessian (exact with `deriv2`) in `active`, parameters of `base` the prediction depends on.
.nodeData <- function(design, groups, etas, base, active, grad, deriv2, cores) {
  E <- design$etaNames[, design$normal$stems, drop = FALSE]
  grad <- grad && length(active) > 0L
  reqGroup <- reqNode <- integer(0); reqCond <- character(0)
  for (i in seq_along(groups)) {
    conds <- design$conditions[[groups[i]]]
    B <- nrow(etas[[i]])
    reqGroup <- c(reqGroup, rep(i, B * length(conds)))
    reqNode <- c(reqNode, rep(seq_len(B), length(conds)))
    reqCond <- c(reqCond, rep(conds, each = B))
  }
  full <- lapply(seq_along(reqGroup), function(r) {
    p <- base
    p[E[groups[reqGroup[r]], ]] <- etas[[reqGroup[r]]][reqNode[r], ]
    p
  })
  rest <- .setdiff(names(base), active)
  kr <- .nodeChain(design, reqCond, full, active, rest, grad, deriv2, cores) %||%
    .nodeR(design, reqCond, full, active, rest, grad, deriv2, cores)
  Q <- length(active)
  lapply(seq_along(groups), function(i) {
    B <- nrow(etas[[i]])
    value <- numeric(B)
    gr <- if (grad) matrix(0, B, Q, dimnames = list(NULL, active))
    he <- if (grad) array(0, c(B, Q, Q), dimnames = list(NULL, active, active))
    for (r in which(reqGroup == i)) {
      loc <- kr$local[[r]]
      b <- reqNode[r]
      value[b] <- value[b] + loc$value
      if (!grad) next
      nm <- .intersect(active, names(loc$gradient))
      gr[b, nm] <- gr[b, nm] + loc$gradient[nm]
      he[b, nm, nm] <- he[b, nm, nm] + loc$hessian[nm, nm]
    }
    list(value = value, gradient = gr, hessian = he)
  })
}

## The data term of every request of .nodeData() in C++: list(local), NULL where the chain does not
## cover the model.
.nodeChain <- function(design, reqCond, full, active, rest, grad, deriv2, cores) {
  loc <- .chainLocal(design$chain, reqCond, if (grad) lapply(full, `[`, active) else full,
                     if (grad) lapply(full, `[`, rest), grad, grad && deriv2, cores, design$bloq, design$df)
  if (!is.null(loc)) list(local = loc)
}

## The data term of every request of .nodeData() through the R evaluation: list(local).
.nodeR <- function(design, reqCond, full, active, rest, grad, deriv2, cores) {
  preds <- .predictMany(design$prdfn, times = lapply(reqCond, function(cn) .emGrid(design, cn)),
                        parsList = if (grad) lapply(full, `[`, active) else full,
                        conditions = reqCond, fixed = if (grad) lapply(full, `[`, rest),
                        deriv = grad, deriv2 = grad && deriv2, cores = cores)
  preds <- unname(as.list(preds))
  errs <- NULL
  if (!is.null(design$errfn)) {
    node <- .fnNode(design$errfn)
    inner <- lapply(preds, attr, "parameters")
    errs <- if (!is.null(node))
      .evalNode(node, .bundle(conds = reqCond, out = preds, pars = inner, fixed = vector("list", length(preds))),
                grad, grad && deriv2, NULL, cores)
    else lapply(seq_along(preds), function(r)
      design$errfn(out = preds[[r]], pars = inner[[r]], deriv = grad, deriv2 = grad && deriv2,
                   conditions = reqCond[r])[[1L]])
    errs <- unname(as.list(errs))
  }
  # the row map of a request depends on the request condition only
  first <- match(unique(reqCond), reqCond)
  meta <- lapply(first, function(r)
    .build_normL2_meta(setNames(list(design$data[[reqCond[r]]]), "c"), list(c = preds[[r]]),
                       if (!is.null(errs)) list(errs[[r]]), "c")[[1]])
  names(meta) <- unique(reqCond)
  pn <- if (grad) unique(unlist(lapply(c(preds, errs), function(f) dimnames(attr(f, "deriv"))[[3]])))
  normL2_kernel(preds, errs, unname(meta[reqCond]), as.character(pn), isTRUE(grad && deriv2), 1L,
                design$bloq, build_hessian = grad, local_blocks = TRUE, t_df = design$df)
}


## log N(eta_b | 0, Omega) of the node etas with, given `grad`, the gradient and Hessian in the
## referenced fixed effects `refPars` (eta_b = phi_b - D beta) and the Cholesky parameters `ch`.
.nodePrior <- function(lay, L, eta, D, refPars, ch, grad) {
  K <- lay$K
  Z <- forwardsolve(L, t(eta))
  value <- -K / 2 * log(2 * pi) - sum(log(diag(L))) - colSums(Z^2) / 2
  if (!grad) return(list(value = value))
  B <- nrow(eta)
  act <- c(refPars, ch)
  gr <- matrix(0, B, length(act), dimnames = list(NULL, act))
  he <- array(0, c(B, length(act), length(act)), dimnames = list(NULL, act, act))
  Oinv <- chol2inv(t(L))
  Dr <- D[, refPars, drop = FALSE]
  if (length(refPars)) {
    gr[, refPars] <- (eta %*% Oinv) %*% Dr
    he[, refPars, refPars] <- rep(-crossprod(Dr, Oinv %*% Dr), each = B)
  }
  if (length(ch)) {
    cd <- .choleskyDerivatives(lay, L)
    ia <- match(ch, lay$cholPars)
    for (a in seq_along(ch)) {
      dO <- cd$dOinv[[ia[a]]]
      gr[, ch[a]] <- -(rowSums((eta %*% dO) * eta) + if (lay$isDiag[[ch[a]]]) 2 else 0) / 2
      if (length(refPars)) he[, refPars, ch[a]] <- he[, ch[a], refPars] <- (eta %*% dO) %*% Dr
      for (k in seq_len(a))
        he[, ch[a], ch[k]] <- he[, ch[k], ch[a]] <-
          -rowSums((eta %*% cd$d2Oinv[ia[a], ia[k], , ]) * eta) / 2
    }
  }
  list(value = value, gradient = gr, hessian = he)
}


## Signed log-sum-exp of the log integrand over the nodes: log L, the weights, the posterior
## moments of the etas and, with per-node gradients and Hessians of the log integrand, those of
## -2 log L (the weighted node Hessians plus the weighted covariance of the node gradients).
.nodeSum <- function(logInt, signs, eta, stems, gr = NULL, he = NULL) {
  M <- max(logInt)
  s <- sum(signs * exp(logInt - M))
  if (!is.finite(s) || s <= 0) return(NULL)
  logL <- M + log(s)
  w <- signs * exp(logInt - logL)
  out <- list(logL = logL, mean = setNames(drop(w %*% eta), stems),
              moment = `dimnames<-`(crossprod(eta, w * eta), list(stems, stems)),
              maxWeight = max(abs(w)), nEff = 1 / sum(w^2))
  if (is.null(gr)) return(out)
  g1 <- drop(w %*% gr)
  Q <- ncol(gr)
  H <- matrix(colSums(w * matrix(he, length(w))), Q, Q) + crossprod(gr, w * gr) - tcrossprod(g1)
  out$gradient <- setNames(-2 * g1, colnames(gr))
  out$hessian <- `dimnames<-`(-2 * H, list(colnames(gr), colnames(gr)))
  out
}


## The quadrature objective: `rebuild(psi, level, tol)` places the nodes around the modes at
## `psi`, each group at the group level (see `rebuild`); `value(pars)` is -2 log L on the frozen nodes,
## frozen in phi = X beta + eta for the mu-referenced fixed effects, in eta for the other etas.
.quadratureObjective <- function(design, start, fixed, control = list()) {
  lay <- design$normal
  E <- design$etaNames[, lay$stems, drop = FALSE]
  meta <- .foceiMeta(design, start, fixed)
  .requireM3(design, meta, "Quadrature")
  cores <- getOption("dMod2.cores", 1L)
  maxLevel <- control$maxLevel %||% (lay$K + 8L)
  maxNodes <- control$maxNodes %||% 1000L
  pruneTol <- control$pruneTol %||% Inf
  shift <- .referenceShift(design)
  state <- new.env()
  state$modes <- matrix(0, nrow(E), lay$K, dimnames = dimnames(E))
  state$stamp <- 0L

  # log integrand of `groups` at their node etas at the rebuild point, forward only
  integrand <- function(groups, etas, psi) {
    L <- .choleskyFactor(lay, psi)
    base <- c(psi, setNames(as.vector(state$modes), as.vector(E)))
    data <- .nodeData(design, groups, etas, base, character(0), FALSE, FALSE, cores)
    Map(function(d, eta) -d$value / 2 + .nodePrior(lay, L, eta, NULL, character(0), character(0), FALSE)$value,
        data, etas)
  }

  # a level within the node budget; it is evaluated, and the tolerance raises a group to it up to
  # maxLevel
  allowed <- function(l) !is.na(l) && .levelNodes(lay$K, l, pruneTol) <= maxNodes
  # the highest level within the budget up to `level`, at least K
  accepted <- function(level) {
    l <- as.integer(level)
    while (l > lay$K && !allowed(l)) l <- l - 1L
    l
  }

  # each group at the lowest level from `level` up, never below the previous one, whose error, the
  # change of the group -2 log L to the next level within the budget, is at most tol / G; a group
  # above tol / G at the highest allowed level is capped
  rebuild <- function(psi, level, tol = Inf, etaInit = state$modes) {
    psi <- c(psi, fixed[.setdiff(names(fixed), names(psi))])
    im <- .innerModes(design, psi, etaInit, fixed, meta, control$innerControl %||% list())
    state$modes <- im$modes
    state$H <- im$H
    state$psi0 <- psi
    state$stamp <- state$stamp + 1L
    state$memo <- NULL
    G <- length(design$groups)
    levels <- pmax(rep(accepted(level), G), state$levels %||% 0L)
    capped <- levels < as.integer(level)
    error <- rep(NA_real_, G)
    nodesAt <- function(i, l) .quadratureNodes(im$modes[i, ], im$H[[i]] / 2, l, pruneTol)
    if (is.finite(tol)) {
      # lock step over the groups, the node values kept by key
      seen <- replicate(G, numeric(0), simplify = FALSE)
      open <- seq_len(G)
      while (length(open)) {
        nxt <- vapply(open, function(i) .nextLevel(lay$K, levels[i]), 0L)
        est <- vapply(nxt, allowed, TRUE)
        up <- est & !is.na(nxt) & nxt <= maxLevel
        ref <- ifelse(est, nxt, vapply(open, function(i) .prevLevel(lay$K, levels[i]), 0L))
        grids <- lapply(seq_along(open), function(j)
          list(at = nodesAt(open[j], levels[open[j]]),
               ref = if (!is.na(ref[j])) nodesAt(open[j], ref[j])))
        need <- lapply(seq_along(open), function(j) {
          gr <- grids[[j]]
          eta <- rbind(gr$at$eta, gr$ref$eta)
          key <- c(gr$at$key, gr$ref$key)
          new <- !duplicated(key) & !(key %in% names(seen[[open[j]]]))
          list(eta = eta[new, , drop = FALSE], key = key[new])
        })
        some <- which(vapply(need, function(n) length(n$key) > 0L, TRUE))
        vals <- if (length(some)) integrand(design$groups[open[some]], lapply(need[some], `[[`, "eta"), psi)
        for (k in seq_along(some)) {
          j <- some[k]
          seen[[open[j]]] <- c(seen[[open[j]]], setNames(vals[[k]], need[[j]]$key))
        }
        done <- logical(length(open))
        for (j in seq_along(open)) {
          i <- open[j]
          at <- grids[[j]]$at; rf <- grids[[j]]$ref
          la <- .nodeSum(at$logW + seen[[i]][at$key], at$signs, at$eta, lay$stems)
          lr <- if (!is.null(rf)) .nodeSum(rf$logW + seen[[i]][rf$key], rf$signs, rf$eta, lay$stems)
          e <- if (is.null(rf)) NA_real_ else if (is.null(la) || is.null(lr)) Inf else 2 * abs(lr$logL - la$logL)
          error[i] <- e
          if (up[j] && e > tol / G) { levels[i] <- nxt[j]; next }
          done[j] <- TRUE
          capped[i] <- !up[j] && !isTRUE(e <= tol / G)
        }
        open <- open[!done]
      }
    }
    state$levels <- levels
    state$levelError <- error
    state$capped <- capped
    state$nodes <- lapply(seq_len(G), function(i) nodesAt(i, levels[i]))
    invisible(state)
  }

  # the nodes of a rebuild, to return to without a new mode search
  fields <- c("modes", "H", "psi0", "levels", "levelError", "capped", "nodes")
  snapshot <- function() mget(fields, envir = state)
  restore <- function(snap) {
    list2env(snap, envir = state)
    state$stamp <- state$stamp + 1L
    state$memo <- NULL
    invisible(state)
  }

  # the data term at the frozen nodes, recomputed only when a parameter it depends on changes or
  # more derivatives are asked for
  dataTerm <- function(base, dataActive, grad, deriv2) {
    key <- list(state$stamp, base[dataActive])
    m <- state$memo
    if (!is.null(m) && identical(m$key, key) && m$grad >= grad && m$deriv2 >= deriv2) return(m$value)
    out <- .nodeData(design, design$groups, lapply(state$nodes, `[[`, "eta"), base, dataActive,
                     grad, deriv2, cores)
    state$memo <- list(key = key, grad = grad, deriv2 = deriv2, value = out)
    out
  }

  value <- function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE, ...) {
    psi <- state$psi0
    psi[names(pars)] <- pars
    psi[names(fixed)] <- fixed
    active <- names(pars)
    refPars <- .intersect(active, shift$pars)
    ch <- .intersect(active, lay$cholPars)
    dataActive <- .setdiff(active, c(refPars, lay$cholPars))
    L <- .choleskyFactor(lay, psi)
    # the prediction sees phi at the rebuild point and the current other parameters
    base <- c(state$psi0, setNames(as.vector(state$modes), as.vector(E)))
    base[dataActive] <- psi[dataActive]
    data <- dataTerm(base, dataActive, deriv, deriv2)
    per <- lapply(seq_along(design$groups), function(i) {
      g <- design$groups[i]
      nodes <- state$nodes[[i]]
      dBeta <- psi[shift$pars] - state$psi0[shift$pars]
      eta <- nodes$eta - matrix(drop(shift$D[[g]] %*% dBeta), nrow(nodes$eta), lay$K, byrow = TRUE)
      pr <- .nodePrior(lay, L, eta, shift$D[[g]], refPars, ch, deriv)
      logInt <- nodes$logW + pr$value - data[[i]]$value / 2
      gr <- he <- NULL
      if (deriv) {
        Q <- length(active)
        gr <- matrix(0, length(logInt), Q, dimnames = list(NULL, active))
        he <- array(0, c(length(logInt), Q, Q), dimnames = list(NULL, active, active))
        pa <- c(refPars, ch)
        gr[, pa] <- pr$gradient
        he[, pa, pa] <- pr$hessian
        if (length(dataActive)) {
          gr[, dataActive] <- -data[[i]]$gradient / 2
          he[, dataActive, dataActive] <- -data[[i]]$hessian / 2
        }
      }
      out <- .nodeSum(logInt, nodes$signs, eta, lay$stems, gr, he)
      if (is.null(out))
        stop("em: the quadrature of group ", g, " is not resolved; raise the level.", call. = FALSE)
      out
    })
    v <- -2 * sum(vapply(per, `[[`, 0, "logL"))
    gr <- if (deriv) Reduce(`+`, lapply(per, `[[`, "gradient"))
    hs <- if (deriv) Reduce(`+`, lapply(per, `[[`, "hessian"))
    # the data of conditions without a group enter as they are
    if (length(design$freeConditions)) {
      free <- design$objfn(pars = psi[.intersect(names(psi), attr(design$objfn, "parameters"))],
                           fixed = fixed, deriv = deriv, deriv2 = deriv2,
                           conditions = design$freeConditions)
      v <- v + free$value
      if (deriv) {
        nm <- .intersect(active, names(free$gradient))
        gr[nm] <- gr[nm] + free$gradient[nm]
        hs[nm, nm] <- hs[nm, nm] + free$hessian[nm, nm]
      }
    }
    out <- objlist(value = v, gradient = gr %||% setNames(numeric(length(pars)), active),
                   hessian = hs %||% matrix(0, length(pars), length(pars)))
    attr(out, "moments") <- lapply(per, `[[`, "moment")
    attr(out, "means") <- lapply(per, `[[`, "mean")
    attr(out, "diagnostics") <- data.frame(group = design$groups,
                                           maxWeight = vapply(per, `[[`, 0, "maxWeight"),
                                           nEff = vapply(per, `[[`, 0, "nEff"))
    out
  }
  class(value) <- c("objfn", "fn")
  attr(value, "parameters") <- names(start)
  list(value = value, rebuild = rebuild, integrand = integrand, state = state, allowed = allowed,
       snapshot = snapshot, restore = restore)
}


## The predicted decrease of -2 log L of a Newton step, g' H^-1 g / 2. Where the Hessian is not
## positive definite, as where the variance of a random effect runs to 0, the directions without
## positive curvature add the linear decrease of a step of `rmax` along their gradient.
.newtonDecrement <- function(ev, rmax = 10) {
  R <- tryCatch(chol(ev$hessian), error = function(e) NULL)
  if (!is.null(R)) return(sum(backsolve(R, ev$gradient, transpose = TRUE)^2) / 2)
  e <- eigen((ev$hessian + t(ev$hessian)) / 2, symmetric = TRUE)
  g <- drop(crossprod(e$vectors, ev$gradient))
  curved <- e$values > sqrt(.Machine$double.eps) * max(abs(e$values))
  sum(g[curved]^2 / e$values[curved]) / 2 + rmax * sqrt(sum(g[!curved]^2))
}

## Whether the prediction and error model have second derivatives; the evaluation stays in the
## memo of the data term.
.hasDeriv2 <- function(q, psi, fixed)
  tryCatch({ q$value(psi, fixed = fixed, deriv2 = TRUE); TRUE }, error = function(e) FALSE)


## The final evaluation over every parameter on the nodes at `psi`, with the exact Hessian where
## the chain has second derivatives, and the conditional covariance of the etas, M - m m' per group.
## The stage trace holds the levels per group as attribute "levels".
.quadratureResult <- function(design, q, psi, fixed, iterations, converged, trace,
                              exact = .hasDeriv2(q, psi, fixed), stopReason = NULL) {
  out <- q$value(psi, fixed = fixed, deriv2 = exact)
  cov <- Map(function(M, m) M - tcrossprod(m), attr(out, "moments"), attr(out, "means"))
  if (!is.null(trace)) attr(trace, "levels") <- .levelTable(design, q$state)
  list(argument = psi, value = out$value, gradient = out$gradient, hessian = out$hessian,
       etaModes = q$state$modes, etaCov = setNames(cov, design$groups),
       converged = converged, iterations = iterations, stopReason = stopReason, stageTrace = trace)
}

## Per group the level, node count, estimated level error and whether the level is capped.
.levelTable <- function(design, st)
  data.frame(group = design$groups, level = st$levels, nodes = vapply(st$nodes, function(n) nrow(n$eta), 0L),
             levelError = st$levelError, capped = st$capped %||% FALSE)


## ECM by quadrature from `start`: per stage the nodes rebuilt at the current parameters and one
## trust over every parameter on them, until the stop "decrement" or "resolution" of em(); the fit
## returns the best iterate.
.fitQuadrature <- function(design, start, fixed, control = list(), verbose = FALSE,
                           etaInit = NULL) {
  lay <- design$normal
  qc <- control$quadrature %||% list()
  levels <- qc$levels %||% (lay$K + 1L)
  tol <- qc$tol %||% 0.1
  nIt <- qc$maxEcmPerStage %||% 10L
  q <- .quadratureObjective(design, start, fixed, modifyList(control, qc))
  if (!is.null(etaInit)) q$state$modes <- etaInit
  psi <- start
  q$rebuild(psi, levels[1], tol)
  exact <- .hasDeriv2(q, psi, fixed)
  step <- .trustControl(list(rinit = 1, rmax = 10, hessianMethod = if (exact) "exact" else "gn",
                             tolControl = list(ftol = 1e-8, mtol = 1e-8)),
                        qc$cm1Control, label = "cm1Control")
  rows <- list(); best <- NULL; reason <- "iterlim"
  for (level in levels) {
    prev <- NULL; reason <- "iterlim"
    for (it in seq_len(nIt)) {
      if (it > 1L || level != levels[1]) q$rebuild(psi, level, tol)
      ev <- q$value(psi, fixed = fixed, deriv2 = exact)
      dec <- .newtonDecrement(ev)
      st <- q$state
      at <- list(psi = psi, value = ev$value, levels = st$levels, snapshot = q$snapshot())
      # values on higher levels supersede those on lower ones
      if (is.null(best) || !identical(best$levels, at$levels) || at$value < best$value) best <- at
      last <- at
      rows[[length(rows) + 1L]] <- data.frame(level = level, ecmIter = it, OFV = ev$value, decrement = dec,
                                              resolution = sum(st$levelError[is.finite(st$levelError)]),
                                              maxLevel = max(st$levels), capped = sum(st$capped),
                                              nodes = sum(vapply(st$nodes, function(n) nrow(n$eta), 0L)),
                                              maxSoftmax = max(attr(ev, "diagnostics")$maxWeight),
                                              nEffMin = min(attr(ev, "diagnostics")$nEff))
      if (verbose) message(sprintf("em(quadrature) level %d ecm %d: -2 log L %.6f, decrement %.2e, capped groups %d",
                                   level, it, ev$value, dec, sum(st$capped)))
      if (dec < (qc$tolDecrement %||% 1e-4)) { reason <- "decrement"; break }
      # the frozen nodes predicted a decrease the rebuilt ones do not show: the step is below
      # what these levels resolve
      if (!is.null(prev) && identical(prev$levels, at$levels) && at$value > prev$value) {
        reason <- "resolution"; break
      }
      if (it == nIt) break
      prev <- at
      cm <- .innerTrust(modifyList(step, list(objfun = q$value, parinit = psi, fixed = fixed,
                                              iterlim = qc$maxCm1Iter %||% 30L)))
      psi <- cm$argument
    }
  }
  if (!identical(best, last)) q$restore(best$snapshot)
  fit <- .quadratureResult(design, q, best$psi, fixed, length(rows), reason != "iterlim",
                           do.call(rbind, rows), exact, reason)
  if (!isFALSE(qc$certificate))
    fit$certificate <- .quadratureCertificate(design, q, best$psi, fixed, qc$certificate %||% list(), tol)
  fit
}


## Laplace-EM: CM-1 on the frozen nodes of level K + 1 (one node would drop the posterior
## variance), Omega from the moments eta* eta*' + (H / 2)^-1, until the extrapolated distance
## to the fixed point is below `epsEcm`; the reported value at level K + 2.
.fitLaplaceEM <- function(design, start, fixed, control = list(), verbose = FALSE) {
  lay <- design$normal
  qc <- control$quadrature %||% list()
  q <- .quadratureObjective(design, start, fixed, modifyList(control, qc))
  structural <- .setdiff(names(start), lay$cholPars)
  psi <- start
  q$rebuild(psi, lay$K + 1L)
  exact <- length(structural) > 0L && .hasDeriv2(q, psi[structural], fixed)
  cm1 <- .trustControl(list(rinit = 1, rmax = 10, hessianMethod = if (exact) "exact" else "gn",
                            tolControl = list(ftol = 1e-8, mtol = 1e-8)),
                       qc$cm1Control, label = "cm1Control")
  rows <- list(); conv <- FALSE; last <- Inf
  for (it in seq_len(qc$maxEcm %||% 200L)) {
    if (it > 1L) q$rebuild(psi, lay$K + 1L)
    before <- psi
    # with every structural parameter fixed only Omega moves
    cm <- if (length(structural))
      .innerTrust(modifyList(cm1, list(objfun = q$value, parinit = psi[structural], fixed = fixed,
                                       iterlim = qc$maxCm1Iter %||% 30L)))
      else list(argument = psi[structural], value = q$value(psi[structural], fixed = c(fixed, psi), deriv = FALSE)$value)
    psi[structural] <- cm$argument
    q$rebuild(psi, lay$K + 1L)
    moments <- lapply(design$groups, function(g)
      tcrossprod(q$state$modes[g, ]) + 2 * chol2inv(chol(q$state$H[[g]])))
    upd <- .updateNormal(design, moments, fixed)
    psi[.intersect(names(psi), names(upd))] <- upd[.intersect(names(psi), names(upd))]
    dpsi <- max(abs(psi - before))
    # EM contracts linearly: the remaining distance is about dpsi rate / (1 - rate)
    rate <- min(dpsi / last, 0.999)
    distance <- if (is.finite(last)) dpsi * rate / (1 - rate) else Inf
    last <- dpsi
    rows[[length(rows) + 1L]] <- data.frame(level = lay$K + 1L, ecmIter = it, OFV = cm$value,
                                            deltaPsi = dpsi, distance = distance)
    if (verbose) message(sprintf("em(laplaceEM) ecm %d: -2 log L %.6f, |delta psi| %.2e", it, cm$value, dpsi))
    if (max(dpsi, distance) < (qc$epsEcm %||% 1e-5)) { conv <- TRUE; break }
  }
  q$rebuild(psi, lay$K + 2L)
  .quadratureResult(design, q, psi, fixed, length(rows), conv, do.call(rbind, rows), exact)
}


## Importance sampling at the estimate: per group a t5 proposal at (mode, (H / 2)^-1), `points`
## draws on 8 randomly shifted rank-1 lattices. A group is flagged when the group quadrature is off by
## more than 3 SE plus max(tol / G, level error), or the Pareto shape exceeds 0.7.
.quadratureCertificate <- function(design, q, psi, fixed, control = list(), tol = Inf) {
  lay <- design$normal
  K <- lay$K
  st <- q$state
  shifts <- 8L
  n <- .nextPrime(ceiling((control$points %||% 256L) / shifts))
  lattice <- .latticePoints(n, K + 1L)
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, globalenv()))
  set.seed(control$seed %||% 1L)
  u <- do.call(rbind, lapply(seq_len(shifts), function(s) (lattice + matrix(stats::runif(K + 1L), n, K + 1L, byrow = TRUE)) %% 1))
  z <- stats::qnorm(u[, 1:K, drop = FALSE])
  scale <- sqrt(stats::qchisq(u[, K + 1L], 5) / 5)
  ests <- st$nodes
  # the quadrature at psi: the rebuilt nodes are those of psi after the final evaluation
  quad <- vapply(seq_along(design$groups), function(i) {
    ev <- q$integrand(design$groups[i], list(ests[[i]]$eta), c(psi, fixed))[[1]]
    s <- .nodeSum(ests[[i]]$logW + ev, ests[[i]]$signs, ests[[i]]$eta, lay$stems)
    if (is.null(s)) NA_real_ else s$logL
  }, 0)
  draws <- lapply(seq_along(design$groups), function(i) {
    # the proposal covariance 2 H^-1 = C C' with C^-1 = chol(H) / sqrt(2), without inverting H,
    # which a collapsed random effect makes ill-conditioned
    Cinv <- chol(st$H[[i]]) / sqrt(2)
    eta <- matrix(st$modes[i, ], nrow(z), K, byrow = TRUE) + t(backsolve(Cinv, t(z))) / scale
    colnames(eta) <- lay$stems
    d <- Cinv %*% (t(eta) - st$modes[i, ])
    logProposal <- lgamma((5 + K) / 2) - lgamma(5 / 2) - K / 2 * log(5 * pi) + sum(log(diag(Cinv))) -
      (5 + K) / 2 * log1p(colSums(d^2) / 5)
    list(eta = eta, logProposal = logProposal)
  })
  logF <- q$integrand(design$groups, lapply(draws, `[[`, "eta"), c(psi, fixed))
  # the error of the levels against the next one where the levels were not chosen by tolerance
  levelError <- st$levelError
  for (i in which(is.na(levelError))) {
    up <- .nextLevel(K, st$levels[i])
    if (!q$allowed(up)) up <- .prevLevel(K, st$levels[i])
    if (is.na(up)) next
    nodes <- .quadratureNodes(st$modes[i, ], st$H[[i]] / 2, up)
    s <- .nodeSum(nodes$logW + q$integrand(design$groups[i], list(nodes$eta), c(psi, fixed))[[1]],
                  nodes$signs, nodes$eta, lay$stems)
    levelError[i] <- if (is.null(s) || is.na(quad[i])) Inf else 2 * abs(s$logL - quad[i])
  }
  rows <- lapply(seq_along(design$groups), function(i) {
    lw <- logF[[i]] - draws[[i]]$logProposal
    mx <- max(lw)
    perShift <- vapply(split(exp(lw - mx), rep(seq_len(shifts), each = n)), mean, 0)
    est <- mean(perShift)
    se <- stats::sd(perShift) / sqrt(shifts) / est
    data.frame(group = design$groups[i], level = st$levels[i], quadrature = -2 * quad[i],
               importance = -2 * (mx + log(est)), se = 2 * se, khat = .paretoShape(lw),
               levelError = levelError[i], capped = isTRUE(st$capped[i]))
  })
  out <- do.call(rbind, rows)
  allowed <- pmax(if (is.finite(tol)) tol / nrow(out) else 0, out$levelError, na.rm = TRUE)
  out$flag <- abs(out$quadrature - out$importance) > 3 * out$se + allowed | out$khat > 0.7
  out
}

## The smallest prime >= n.
.nextPrime <- function(n) {
  n <- max(2L, as.integer(n))
  repeat {
    if (n < 4L || all(n %% 2:floor(sqrt(n)) != 0L)) return(n)
    n <- n + 1L
  }
}

## A Korobov rank-1 lattice of n points in [0, 1)^s whose generator minimizes the P2 criterion.
.latticePoints <- function(n, s) {
  k <- 0:(n - 1L)
  b2 <- function(x) x^2 - x + 1 / 6
  merit <- function(a) {
    z <- (a^(0:(s - 1L))) %% n
    x <- outer(k, z) %% n / n
    mean(apply(1 + 2 * pi^2 * b2(x), 1, prod))
  }
  cand <- 2:max(2L, n - 1L)
  a <- cand[which.min(vapply(cand, merit, 0))]
  outer(k, (a^(0:(s - 1L))) %% n) %% n / n
}

## The Pareto shape estimate of the upper tail of the log weights `lw` (Zhang and Stephens,
## with the weakly informative prior of Vehtari et al.), over the largest min(S / 5, 3 sqrt(S)).
.paretoShape <- function(lw) {
  S <- length(lw)
  M <- ceiling(min(0.2 * S, 3 * sqrt(S)))
  lw <- sort(lw - max(lw))
  cut <- lw[S - M]
  x <- exp(lw[(S - M + 1):S]) - exp(cut)
  if (all(x <= 0)) return(-Inf)
  N <- length(x)
  grid <- 30 + floor(sqrt(N))
  xstar <- x[floor(N / 4 + 0.5)]
  theta <- 1 / x[N] + (1 - sqrt(grid / (seq_len(grid) - 0.5))) / 3 / xstar
  profile <- vapply(theta, function(t) { k <- mean(log1p(-t * x)); log(-t / k) - k - 1 }, 0)
  lt <- N * profile
  w <- exp(lt - max(lt)); w <- w / sum(w)
  k <- mean(log1p(-sum(theta * w) * x))
  (k * N + 5) / (N + 10)
}
