# Spike-and-slab random effects with a given slab mean and spike, -Inf included, beside normal ones: per
# group the sum over the subsets of slab etas, each a Laplace integral or a sparse-grid quadrature (methods
# article Sparsity priors, section Spike at a given value).


## The most Laplace integrals per evaluation, summed over the groups.
.spikeBudget <- 4096L

## Per stem of the spike-and-slab blocks: block, slab mean, spike and the hyperparameter names.
.spikeLayout <- function(design) {
  blocks <- design$blocks
  do.call(rbind, lapply(which(vapply(blocks, function(b) b$dist == "spikeSlab", TRUE)), function(i) {
    b <- blocks[[i]]
    data.frame(stem = b$stems, block = i, mean = b$mean, spike = b$spike, prob = b$prob, tau = b$tau,
               probName = b$hyper[["prob"]], tauName = b$hyper[["tau"]], stringsAsFactors = FALSE)
  }))
}

## The prior of a group: slab means `mu` and standard deviations `tau` of the slab etas, and with normal
## random effects their etas `normal`, Omega^-1 `Oinv` and log|Omega| `logdetO`.
.spikePriorOf <- function(mu, tau, normal = character(0), Oinv = NULL, logdetO = 0)
  list(mu = mu, tau = tau, normal = normal, Oinv = Oinv, logdetO = logdetO)

## -2 log of the prior `pr` without the normalizer at the etas `p` (named, the slab etas of `p` in the
## slab), with gradient and Hessian in `p`.
.spikePriorTerm <- function(pr, p) {
  n <- length(p)
  gr <- setNames(numeric(n), names(p))
  H <- matrix(0, n, n, dimnames = list(names(p), names(p)))
  s <- .intersect(names(p), names(pr$tau))
  d <- p[s] - pr$mu[s]
  v <- sum((d / pr$tau[s])^2)
  gr[s] <- 2 * d / pr$tau[s]^2
  if (length(s)) H[s, s] <- diag(2 / pr$tau[s]^2, length(s))
  if (length(pr$normal)) {
    q <- pr$normal
    Oq <- drop(pr$Oinv %*% p[q])
    v <- v + sum(p[q] * Oq)
    gr[q] <- 2 * Oq
    H[q, q] <- 2 * pr$Oinv
  }
  list(value = v, gradient = gr, hessian = H)
}

## The normalizer of the prior over the etas `on`, -2 log of (2 pi)^(-d/2) |V|^(-1/2) without the 2 pi.
.spikePriorNorm <- function(pr, on)
  sum(log(pr$tau[.intersect(on, names(pr$tau))]^2)) + if (length(pr$normal)) pr$logdetO else 0

## The Laplace integral of group `g` over the etas `on` (the normal etas and the slab etas of the subset,
## the other slab etas at the spike): -2 log Z with the normalizer of the prior, the mode and the
## covariance 2 H^-1; NULL where the inner search fails.
.spikeSubset <- function(design, objs, g, theta, fixed, etas, on, pr, start) {
  conds <- design$conditions[[g]]
  off <- .setdiff(etas, on)
  base <- c(theta, fixed)
  if (!length(on)) {
    out <- Reduce(`+`, .objEvalMany(objs[conds], rep(list(c(base, start[off])), length(conds)), deriv = FALSE))
    if (!is.finite(out$value)) return(NULL)
    return(list(m2logZ = out$value, mode = numeric(0), cov = matrix(0, 0, 0)))
  }
  fx <- c(base, start[off])
  objfun <- function(p, fixed = NULL, deriv = TRUE, ...) {
    outs <- .objEvalMany(objs[conds], rep(list(p), length(conds)), deriv = deriv, fixed = fixed)
    out <- Reduce(`+`, outs)
    prior <- .spikePriorTerm(pr, p)
    gr <- prior$gradient
    if (!is.null(out$gradient)) gr <- gr + out$gradient[names(p)]
    H <- prior$hessian
    if (!is.null(out$hessian)) H <- H + as.matrix(out$hessian)[names(p), names(p)]
    objlist(value = out$value + prior$value, gradient = gr, hessian = H)
  }
  fit <- tryCatch(.innerTrust(list(objfun = objfun, parinit = start[on], fixed = fx, rinit = 1, rmax = 10,
                                   iterlim = 100L)),
                  error = function(e) NULL)
  if (is.null(fit) || !is.finite(fit$value)) return(NULL)
  at <- objfun(fit$argument, fixed = fx)
  H <- (at$hessian + t(at$hessian)) / 2
  ch <- tryCatch(chol(H / 2), error = function(e) NULL)
  if (is.null(ch)) return(NULL)
  list(m2logZ = at$value + .spikePriorNorm(pr, on) + 2 * sum(log(diag(ch))),
       mode = fit$argument, cov = chol2inv(ch), prec = H / 2)
}

## The data -2 log L of group `g` at the parameter vectors `full`, one per node, in one call of the
## native chain (threads by option dMod2.cores only): values and, with `active`, gradients (nodes x
## active) and Gauss-Newton Hessians (nodes x active x active).
.spikeNodes <- function(design, g, full, active = character(0)) {
  conds <- design$conditions[[g]]
  B <- length(full)
  grad <- length(active) > 0L
  reqCond <- rep(conds, each = B)
  reqFull <- rep(full, times = length(conds))
  rest <- .setdiff(names(full[[1L]]), active)
  cores <- getOption("dMod2.cores", 1L)
  kr <- .nodeChain(design, reqCond, reqFull, active, rest, grad, FALSE, cores) %||%
    .nodeR(design, reqCond, reqFull, active, rest, grad, FALSE, cores)
  node <- rep(seq_len(B), times = length(conds))
  Q <- length(active)
  value <- numeric(B)
  gr <- if (grad) matrix(0, B, Q, dimnames = list(NULL, active))
  he <- if (grad) array(0, c(B, Q, Q), dimnames = list(NULL, active, active))
  for (r in seq_along(node)) {
    loc <- kr$local[[r]]
    b <- node[r]
    value[b] <- value[b] + loc$value
    if (!grad) next
    nm <- .intersect(active, names(loc$gradient))
    gr[b, nm] <- gr[b, nm] + loc$gradient[nm]
    he[b, nm, nm] <- he[b, nm, nm] + as.matrix(loc$hessian)[nm, nm]
  }
  list(value = value, gradient = gr, hessian = he)
}

## The subsets `rs` (from .spikeSubset()) of group `g` on adaptive sparse grids (.adaptiveGrids()), a subset of
## relative weight `w` to the tolerance tol / w, from the index sets `starts`: -2 log Z, nodes, signed posterior
## weights, index set, level and whether `tol` was met.
.spikeQuadrature <- function(design, g, theta, fixed, etas, rs, pr, alpha, w, starts, qc = list()) {
  tol <- qc$tol
  K <- vapply(rs, function(r) length(r$on), 0L)
  fx <- lapply(rs, function(r) c(theta, fixed, alpha[.setdiff(etas, r$on)]))
  names0 <- unique(c(names(theta), names(fixed), etas))
  Rs <- lapply(rs, function(r) chol(r$prec))
  etaOf <- function(i, z) {
    eta <- t(sqrt(2) * backsolve(Rs[[i]], t(z))) + matrix(rs[[i]]$mode, nrow(z), K[i], byrow = TRUE)
    `colnames<-`(eta, rs[[i]]$on)
  }
  logf <- function(which, z) {
    etaList <- Map(etaOf, which, z)
    full <- unlist(Map(function(i, eta) lapply(seq_len(nrow(eta)), function(b) c(fx[[i]], eta[b, ])[names0]),
                       which, etaList), recursive = FALSE)
    data <- .spikeNodes(design, g, full)$value
    owner <- rep(seq_along(which), vapply(etaList, nrow, 0L))
    Map(function(k, eta) {
      prior <- vapply(seq_len(nrow(eta)), function(b) .spikePriorTerm(pr, eta[b, ])$value, 0)
      -(data[owner == k] + prior) / 2
    }, seq_along(which), etaList)
  }
  wn <- w / sum(w)
  grids <- .adaptiveGrids(K, K + 1L, starts, logf, tol / pmax(wn, 1e-6), qc$maxNodes %||% 1000L,
                          qc$maxLevel %||% (K + 8L))
  Map(function(r, gr) {
    nd <- .gridNodes(gr, r$mode, r$prec)
    lt <- nd$logW + gr$logValues[match(nd$key, gr$key)]
    top <- max(lt)
    sw <- nd$signs * exp(lt - top)
    if (!all(is.finite(sw)) || sum(sw) <= 0) return(r)
    r$m2logZ <- -2 * (top + log(sum(sw))) + .spikePriorNorm(pr, r$on) + length(r$on) * log(2 * pi)
    r$nodes <- nd$eta
    r$logW <- nd$logW
    r$signs <- nd$signs
    r$omega <- sw / sum(sw)
    r$indices <- gr$indices
    r$level <- gr$level
    r$quadConverged <- !gr$capped
    r
  }, rs, grids)
}

## All subsets of `n` items as logical rows.
.subsets <- function(n) {
  if (!n) return(matrix(FALSE, 1L, 0L))
  as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), n))) > 0
}

## The E-step of every group: per subset of the slab etas that reach the group data the weight, mode and
## covariance; -2 log I per group; the probability of the spike per eta. Groups beyond `budget` subsets,
## the largest first, take EC (.spikeECGroup()) and the subset of inclusion above 1/2 as their one part.
.spikeEStep <- function(design, objs, lay, theta, fixed, hyp, warm, bases, budget, control = list(),
                        quadrature = NULL, gridWarm = list(), omega = NULL) {
  E <- design$etaNames
  reach <- .reachOf(design)
  stems <- lay$stem
  nstems <- design$normal$stems
  lbStem <- setNames(ifelse(stems %in% names(bases), log(bases[stems]), NA_real_), stems)
  sizes <- setNames(rowSums(reach[design$groups, stems, drop = FALSE]), design$groups)
  exact <- design$groups[order(sizes)][cumsum(2^sort(sizes)) <= budget]
  lapply(setNames(nm = design$groups), function(g) {
    slab <- unname(E[g, stems])
    nrm <- unname(E[g, nstems])
    etas <- c(slab, nrm)
    mu <- setNames(lay$mean, slab)
    # the spike of the slab etas, the normal ones are always integrated
    alpha <- setNames(c(lay$spike, rep(0, length(nrm))), etas)
    tauG <- setNames(hyp$tau[stems], slab)
    probG <- setNames(hyp$prob[stems], slab)
    pr <- .spikePriorOf(mu, tauG, nrm, omega$Oinv, omega$logdetO)
    cand <- unname(E[g, stems[reach[g, stems]]])
    keyOf <- function(on) paste(on, collapse = ",")
    startOf <- function(on) {
      v <- alpha
      v[nrm] <- 0
      v[on] <- warm[[g]][[keyOf(on)]] %||% c(mu, setNames(rep(0, length(nrm)), nrm))[on]
      v
    }
    if (!g %in% exact) {
      # the data term first expanded at the mode with every candidate in the slab
      start <- startOf(c(nrm, cand))
      full <- .spikeSubset(design, objs, g, theta, fixed, etas, c(nrm, cand), pr, start)
      if (!is.null(full)) start[c(nrm, cand)] <- full$mode
      ec <- .spikeECGroup(design, objs, g, theta, fixed, etas, cand, nrm, pr, alpha, probG,
                          setNames(lbStem, slab), start, control)
      if (is.null(ec))
        stop("em: the expectation consistent (EC) approximation of group ", g, " did not converge.", call. = FALSE)
      on <- c(nrm, ec$on)
      r <- .spikeSubset(design, objs, g, theta, fixed, etas, on, pr, startOf(on))
      if (is.null(r)) stop("em: the inner search of group ", g, " failed.", call. = FALSE)
      r$on <- on
      r$weight <- 1
      inclusion <- probG
      inclusion[cand] <- ec$inclusion[cand]
      return(list(parts = list(r), m2logI = ec$m2logI, inclusion = inclusion, slab = slab, normal = nrm,
                  etas = etas, cand = cand, mu = mu, alpha = alpha,
                  ec = list(inclusion = ec$inclusion, dev2 = ec$dev2, normalMoment = ec$normalMoment)))
    }
    S <- .subsets(length(cand))
    parts <- lapply(seq_len(nrow(S)), function(i) {
      on <- c(nrm, cand[S[i, ]])
      r <- .spikeSubset(design, objs, g, theta, fixed, etas, on, pr, startOf(on))
      if (is.null(r)) return(NULL)
      r$on <- on
      r$logw <- sum(log(probG[cand[S[i, ]]])) + sum(log1p(-probG[cand[!S[i, ]]])) - r$m2logZ / 2
      r
    })
    parts <- Filter(Negate(is.null), parts)
    if (!length(parts)) stop("em: no subset of the spike-and-slab random effects of group ", g,
                             " has a finite integral.", call. = FALSE)
    lw <- vapply(parts, `[[`, 0, "logw")
    sel <- which(lw >= max(lw) - log(1e8) & vapply(parts, function(p) length(p$on) > 0L, TRUE))
    if (!is.null(quadrature) && length(sel)) {
      # the subsets of non-negligible Laplace weight by quadrature, the others keep the Laplace value
      keys <- vapply(parts[sel], function(p) keyOf(p$on), "")
      qs <- .spikeQuadrature(design, g, theta, fixed, etas, parts[sel], pr, alpha,
                             exp(lw[sel] - max(lw)), lapply(keys, function(k) gridWarm[[g]][[k]]),
                             modifyList(quadrature, list(tol = (quadrature$tol %||% 0.1) / length(design$groups))))
      for (k in seq_along(sel)) {
        qs[[k]]$logw <- qs[[k]]$logw + (parts[[sel[k]]]$m2logZ - qs[[k]]$m2logZ) / 2
        parts[[sel[k]]] <- qs[[k]]
      }
      lw <- vapply(parts, `[[`, 0, "logw")
    }
    logI <- max(lw) + log(sum(exp(lw - max(lw))))
    for (i in seq_along(parts)) parts[[i]]$weight <- exp(parts[[i]]$logw - logI)
    inclusion <- setNames(vapply(slab, function(e) {
      if (!e %in% cand) return(probG[[e]])
      sum(vapply(parts, function(p) if (e %in% p$on) p$weight else 0, 0))
    }, 0), slab)
    list(parts = parts, m2logI = -2 * logI, inclusion = inclusion, slab = slab, normal = nrm, etas = etas,
         cand = cand, mu = mu, alpha = alpha)
  })
}

## E[eta eta'] of the normal etas of every group under the Laplace posteriors of the subsets with their
## weights; EC groups by the EC Gaussian.
.spikeNormalMoments <- function(est) {
  lapply(est, function(e) {
    q <- e$normal
    if (!is.null(e$ec)) return(e$ec$normalMoment)
    Reduce(`+`, lapply(e$parts, function(p) {
      j <- match(q, p$on)
      p$weight * (tcrossprod(p$mode[q]) + p$cov[j, j, drop = FALSE])
    }))
  })
}

## Omega^-1 and log|Omega| from the Cholesky parameters in `pars`.
.spikeOmega <- function(design, pars) {
  if (!design$normal$K) return(NULL)
  L <- .choleskyFactor(design$normal, pars)
  list(Oinv = chol2inv(t(L)), logdetO = 2 * sum(log(diag(L))))
}

## The M-step of the hyperparameters: prob as the mean inclusion, tau^2 as the mean of
## (eta - mean)^2 under the slab, per block; fixed values stay.
.spikeHyperStep <- function(lay, est, hyp, hypers) {
  stems <- lay$stem
  for (b in unique(lay$block)) {
    st <- lay$stem[lay$block == b]
    k <- match(st, stems)
    incl <- unlist(lapply(est, function(e) e$inclusion[k]))
    if (lay$probName[lay$block == b][1] %in% hypers) hyp$prob[st] <- min(max(mean(incl), 1e-6), 1 - 1e-6)
    if (!lay$tauName[lay$block == b][1] %in% hypers) next
    num <- den <- 0
    for (e in est) {
      own <- e$slab[k]
      if (!is.null(e$ec)) {
        on <- .intersect(names(e$ec$inclusion), own)
        num <- num + sum(e$ec$inclusion[on] * e$ec$dev2[on])
        den <- den + sum(e$ec$inclusion[on])
      }
      for (p in if (is.null(e$ec)) e$parts) {
        on <- .intersect(p$on, own)
        if (!length(on)) next
        j <- match(on, p$on)
        num <- num + p$weight * sum((p$mode[on] - e$mu[on])^2 + diag(p$cov)[j])
        den <- den + p$weight * length(on)
      }
      # an eta off the group data keeps the slab of the prior
      blind <- .setdiff(own, e$cand)
      num <- num + sum(hyp$prob[st[match(blind, own)]] * hyp$tau[st[match(blind, own)]]^2)
      den <- den + sum(hyp$prob[st[match(blind, own)]])
    }
    if (den > 0) hyp$tau[st] <- sqrt(num / den)
  }
  hyp
}

## The M-step of the parameters: the weighted subset objectives (data plus prior) minimized jointly in `theta`
## and one copy of the subset etas per subset, every copy then the inner mode at the new `theta`; the log
## determinant of the Laplace integrals enters the weights only.
.spikeM <- function(design, objs, est, extra, lay, hyp, omega) {
  items <- list()
  for (g in names(est)) {
    e <- est[[g]]
    pr <- .spikePriorOf(e$mu, setNames(hyp$tau[lay$stem], e$slab), e$normal, omega$Oinv, omega$logdetO)
    for (p in e$parts) {
      if (p$weight < 1e-8) next
      k <- length(items) + 1L
      copies <- if (length(p$on)) paste0(".s", k, ".", p$on) else character(0)
      items[[k]] <- list(g = g, on = p$on, copies = copies, w = p$weight, mode = p$mode,
                         off = e$alpha[.setdiff(e$etas, p$on)], pr = pr)
    }
  }
  init <- unlist(lapply(items, function(it) setNames(it$mode, it$copies)))
  objfun <- function(pars, fixed = NULL, deriv = TRUE, ...) {
    th <- .setdiff(names(pars), names(init))
    n <- length(pars)
    val <- 0
    gr <- setNames(numeric(n), names(pars))
    H <- matrix(0, n, n, dimnames = list(names(pars), names(pars)))
    for (it in items) {
      conds <- design$conditions[[it$g]]
      local <- c(pars[th], setNames(pars[it$copies], it$on))
      outs <- .objEvalMany(objs[conds], rep(list(local), length(conds)), deriv = deriv, fixed = c(fixed, it$off))
      out <- Reduce(`+`, outs)
      prior <- .spikePriorTerm(it$pr, local[it$on])
      val <- val + it$w * (out$value + prior$value)
      if (!deriv) next
      map <- setNames(c(th, it$copies), c(th, it$on))
      gl <- out$gradient[names(local)]
      gl[it$on] <- gl[it$on] + prior$gradient[it$on]
      gr[map[names(local)]] <- gr[map[names(local)]] + it$w * gl
      Hl <- as.matrix(out$hessian)[names(local), names(local)]
      if (length(it$on)) Hl[it$on, it$on] <- Hl[it$on, it$on] + prior$hessian[it$on, it$on]
      idx <- map[names(local)]
      H[idx, idx] <- H[idx, idx] + it$w * Hl
    }
    out <- objlist(value = val, gradient = gr, hessian = H)
    if (!is.null(extra) && length(th)) {
      ex <- extra(pars[th], fixed = fixed, deriv = deriv)
      out$value <- out$value + ex$value
      if (deriv) {
        out$gradient[th] <- out$gradient[th] + ex$gradient[th]
        out$hessian[th, th] <- out$hessian[th, th] + as.matrix(ex$hessian)[th, th]
      }
    }
    out
  }
  list(objfun = objfun, init = init)
}

## -2 log L on the frozen nodes of every subset (a subset without a grid by the Laplace node at the mode) in
## the fixed effects, Cholesky parameters, log tau and logit prob, gradient by Fisher's and Hessian by Louis'
## identity; mu-referenced normal etas frozen in phi. `back` maps the argument to the parameters.
.spikeFrozen <- function(design, est, lay, theta, omegaPars, hyp, hypers, fixed, extra) {
  nl <- design$normal
  groups <- lapply(names(est), function(g) {
    e <- est[[g]]
    rows <- list(); logW <- numeric(0); signs <- numeric(0); subset <- list()
    for (p in e$parts) {
      if (p$weight < 1e-12) next
      if (!is.null(p$nodes)) {
        for (i in seq_len(nrow(p$nodes))) {
          v <- e$alpha; v[p$on] <- p$nodes[i, p$on]
          rows[[length(rows) + 1L]] <- v
          subset[[length(subset) + 1L]] <- p$on
        }
        logW <- c(logW, p$logW); signs <- c(signs, p$signs)
        next
      }
      v <- e$alpha; v[p$on] <- p$mode
      rows[[length(rows) + 1L]] <- v
      subset[[length(subset) + 1L]] <- p$on
      d <- length(p$on)
      logW <- c(logW, if (d) d / 2 * log(2 * pi) - 0.5 * as.numeric(determinant(p$prec)$modulus) else 0)
      signs <- c(signs, 1)
    }
    list(g = g, eta = do.call(rbind, rows), logW = logW, signs = signs, subset = subset, e = e)
  })
  probFree <- .intersect(lay$probName, hypers)
  tauFree <- .intersect(lay$tauName, hypers)
  cholFree <- names(omegaPars)
  tn <- names(theta)
  shift <- if (nl$K) .referenceShift(design) else list(pars = character(0))
  refAll <- .intersect(shift$pars, tn)
  beta0 <- theta[refAll]
  init <- c(theta, omegaPars, setNames(log(hyp$tau[lay$stem[match(tauFree, lay$tauName)]]), tauFree),
            setNames(stats::qlogis(hyp$prob[lay$stem[match(probFree, lay$probName)]]), probFree))
  back <- function(arg) {
    h <- hyp
    for (n in tauFree) h$tau[lay$stem[lay$tauName == n]] <- exp(arg[[n]])
    for (n in probFree) h$prob[lay$stem[lay$probName == n]] <- stats::plogis(arg[[n]])
    list(theta = arg[tn], omegaPars = arg[cholFree], hyp = h)
  }
  objfun <- function(pars, deriv = TRUE, ...) {
    all <- c(pars, fixed)
    act <- names(pars)
    Q <- length(act)
    total <- 0
    gr <- setNames(numeric(Q), act)
    H <- matrix(0, Q, Q, dimnames = list(act, act))
    h <- back(c(init[.setdiff(names(init), act)], pars))$hyp
    L <- if (nl$K) .choleskyFactor(nl, all)
    thA <- .intersect(tn, act)
    refA <- .intersect(refAll, act)
    dataA <- .setdiff(thA, refAll)
    # the data see phi at the rebuild point: the referenced fixed effects at beta0, the etas frozen
    base <- all[.intersect(names(all), c(tn, names(fixed)))]
    base[refAll] <- beta0
    for (G in groups) {
      e <- G$e
      B <- nrow(G$eta)
      full <- lapply(seq_len(B), function(b) c(base, G$eta[b, ]))
      data <- .spikeNodes(design, G$g, full, if (deriv) dataA else character(0))
      logInt <- G$logW - data$value / 2
      g1 <- if (deriv) matrix(0, B, Q, dimnames = list(NULL, act))
      h1 <- if (deriv) array(0, c(B, Q, Q), dimnames = list(NULL, act, act))
      if (deriv && length(dataA)) {
        g1[, dataA] <- -data$gradient[, dataA, drop = FALSE] / 2
        h1[, dataA, dataA] <- -data$hessian[, dataA, dataA, drop = FALSE] / 2
      }
      if (nl$K) {
        ch <- .intersect(cholFree, act)
        etaN <- G$eta[, e$normal, drop = FALSE]
        D <- if (length(refAll)) shift$D[[G$g]][, refAll, drop = FALSE]
        if (length(refAll)) etaN <- etaN - matrix(drop(D %*% (all[refAll] - beta0)), B, nl$K, byrow = TRUE)
        pa <- c(refA, ch)
        pr <- .nodePrior(nl, L, etaN, D, refA, ch, deriv && length(pa) > 0L)
        logInt <- logInt + pr$value
        if (deriv && length(pa)) {
          g1[, pa] <- g1[, pa] + pr$gradient[, pa, drop = FALSE]
          h1[, pa, pa] <- h1[, pa, pa] + pr$hessian[, pa, pa, drop = FALSE]
        }
      }
      inCand <- e$slab %in% e$cand
      prob <- h$prob[lay$stem]
      for (b in seq_len(B)) {
        on <- .intersect(G$subset[[b]], e$slab)
        st <- lay$stem[match(on, e$slab)]
        tau <- h$tau[st]
        dev2 <- ((G$eta[b, on] - e$mu[on]) / tau)^2
        z <- e$slab %in% G$subset[[b]]
        logInt[b] <- logInt[b] + sum(-0.5 * log(2 * pi) - log(tau) - dev2 / 2) +
          sum(ifelse(z, log(prob), log1p(-prob))[inCand])
        if (!deriv) next
        for (n in .intersect(tauFree, act)) {
          k <- lay$tauName[match(st, lay$stem)] == n
          g1[b, n] <- g1[b, n] + sum(-1 + dev2[k])
          h1[b, n, n] <- h1[b, n, n] - 2 * sum(dev2[k])
        }
        for (n in .intersect(probFree, act)) {
          k <- lay$probName == n & inCand
          r <- unname(prob[k])
          g1[b, n] <- g1[b, n] + sum(z[k] - r)
          h1[b, n, n] <- h1[b, n, n] - sum(r * (1 - r))
        }
      }
      ns <- .nodeSum(logInt, G$signs, matrix(0, B, 1, dimnames = list(NULL, "z")), "z", g1, h1)
      if (is.null(ns)) return(objlist(value = Inf, gradient = gr, hessian = H))
      total <- total - 2 * ns$logL
      if (deriv) { gr <- gr + ns$gradient[act]; H <- H + ns$hessian[act, act] }
    }
    out <- objlist(value = total, gradient = gr, hessian = H)
    if (!is.null(extra) && length(thA)) {
      ex <- extra(pars[thA], fixed = fixed, deriv = deriv)
      out$value <- out$value + ex$value
      if (deriv) {
        n <- names(ex$gradient)
        out$gradient[n] <- out$gradient[n] + ex$gradient
        out$hessian[n, n] <- out$hessian[n, n] + as.matrix(ex$hessian)
      }
    }
    out
  }
  list(objfun = objfun, init = init, back = back)
}

## Fit of a design whose random effects are spike-and-slab ones, beside normal ones, by the sum over
## subsets: EM on the Laplace integrals, then with `quadrature` (control$quadrature) trust steps on
## frozen quadrature nodes.
.fitSpikeSubsets <- function(design, start, fixed, control = list(), verbose = FALSE, quadrature = NULL) {
  .requireM3(design, NULL, "FOCEI")
  if (any(vapply(design$blocks, function(b) b$dist == "dirichletProcess", TRUE)))
    stop("em: spikeSlab() with `mean` or `spike` other than 0, or with method = \"quadrature\", is fitted ",
         "beside normal random effects only.", call. = FALSE)
  lay <- .spikeLayout(design)
  objs <- .conditionObjectives(design)
  extra <- .extraObjective(design)
  hypers <- .intersect(c(lay$probName, lay$tauName), names(start))
  cholPars <- design$normal$cholPars
  given <- c(start, fixed)
  value <- function(col, name) setNames(ifelse(is.na(lay[[col]]), unname(given[lay[[name]]]), lay[[col]]), lay$stem)
  st <- new.env()
  st$theta <- start[.setdiff(names(start), c(hypers, cholPars))]
  st$omegaPars <- start[.intersect(cholPars, names(start))]
  st$hyp <- list(prob = value("prob", "probName"), tau = value("tau", "tauName"))
  st$warm <- list(); st$gridWarm <- list(); st$rows <- list()
  budget <- control$sparseFit$budget %||% .spikeBudget
  bases <- .spikeBases(design)
  omegaOf <- function() .spikeOmega(design, c(st$omegaPars, fixed))
  keyOf <- function(e) vapply(e$parts, function(p) paste(p$on, collapse = ","), "")
  eStep <- function(qc) {
    est <- .spikeEStep(design, objs, lay, st$theta, fixed, st$hyp, st$warm, bases, budget, control$sparseFit$ec, qc,
                       st$gridWarm, omegaOf())
    st$warm <- lapply(est, function(e) setNames(lapply(e$parts, `[[`, "mode"), keyOf(e)))
    st$gridWarm <- lapply(est, function(e) setNames(lapply(e$parts, `[[`, "indices"), keyOf(e)))
    est
  }
  record <- function(stage, k, value, delta) {
    st$rows[[length(st$rows) + 1L]] <- data.frame(stage = stage, ecmIter = k, OFV = value, delta = delta)
    if (verbose) message(sprintf("em(spike subsets) %s %d: -2 log L %.6f, change %.2e", stage, k, value, delta))
  }

  # Laplace stage: EM until no estimate changes by more than epsEcm
  conv <- FALSE
  for (k in seq_len(control$sparseFit$maxEcm %||% 200L)) {
    est <- eStep(NULL)
    before <- c(st$theta, unlist(st$hyp), st$omegaPars)
    st$hyp <- .spikeHyperStep(lay, est, st$hyp, hypers)
    if (length(st$omegaPars))
      st$omegaPars[] <- .updateNormal(design, .spikeNormalMoments(est), fixed)[names(st$omegaPars)]
    if (length(st$theta)) {
      ms <- .spikeM(design, objs, est, extra, lay, st$hyp, omegaOf())
      m <- .innerTrust(list(objfun = ms$objfun, parinit = c(st$theta, ms$init), fixed = fixed, rinit = 1, rmax = 10,
                            iterlim = 100L))
      st$theta <- m$argument[names(st$theta)]
    }
    delta <- max(abs(c(st$theta, unlist(st$hyp), st$omegaPars) - before))
    record("laplace", k, sum(vapply(est, `[[`, 0, "m2logI")), delta)
    conv <- delta < (control$sparseFit$epsEcm %||% 1e-6)
    if (conv) break
  }
  stopReason <- if (conv) "converged" else "iterlim"

  # quadrature stage: a trust step over every parameter on the frozen nodes, until a Newton step on the
  # rebuilt nodes gains less than tolDecrement
  if (!is.null(quadrature)) {
    conv <- FALSE; stopReason <- "iterlim"
    for (k in seq_len(quadrature$maxEcmPerStage %||% 10L)) {
      est <- eStep(quadrature)
      fr <- .spikeFrozen(design, est, lay, st$theta, st$omegaPars, st$hyp, hypers, fixed, extra)
      ev <- fr$objfun(fr$init)
      dec <- .newtonDecrement(ev)
      record("quadrature", k, ev$value, dec)
      if (dec < (quadrature$tolDecrement %||% 1e-4)) { conv <- TRUE; stopReason <- "decrement"; break }
      tc <- modifyList(list(rinit = 1, rmax = 10), quadrature$cm1Control %||% list())
      m <- .innerTrust(modifyList(tc, list(objfun = fr$objfun, parinit = fr$init,
                                           iterlim = quadrature$maxCm1Iter %||% 30L)))
      b <- fr$back(m$argument)
      st$theta <- b$theta; st$omegaPars[] <- b$omegaPars; st$hyp <- b$hyp
    }
  }

  est <- eStep(quadrature)
  theta <- st$theta; hyp <- st$hyp
  total <- sum(vapply(est, `[[`, 0, "m2logI"))
  if (!is.null(extra)) total <- total + extra(theta, fixed = fixed, deriv = FALSE)$value
  E <- design$etaNames
  modes <- clusters <- spike <- matrix(NA_real_, nrow(E), ncol(E), dimnames = dimnames(E))
  cov <- list()
  for (gi in seq_along(design$groups)) {
    e <- est[[gi]]
    best <- e$parts[[which.max(vapply(e$parts, `[[`, 0, "weight"))]]
    v <- e$alpha
    v[best$on] <- best$mode
    modes[gi, match(e$etas, E[gi, ])] <- v
    sc <- match(e$slab, E[gi, ])
    clusters[gi, sc] <- as.numeric(e$slab %in% best$on)
    spike[gi, sc] <- pmin(pmax(1 - e$inclusion, 0), 1)
    C <- matrix(0, ncol(E), ncol(E), dimnames = list(colnames(E), colnames(E)))
    k <- match(best$on, E[gi, ])
    C[k, k] <- best$cov
    cov[[gi]] <- C
  }
  hyperValues <- c(setNames(hyp$prob[lay$stem], lay$probName), setNames(hyp$tau[lay$stem], lay$tauName))
  capped <- !unlist(lapply(est, function(e) lapply(e$parts, function(p) p$quadConverged %||% TRUE)))
  if (!is.null(quadrature) && any(capped))
    warning("em: the quadrature of a subset reached maxLevel before tol.", call. = FALSE)
  list(argument = c(theta, st$omegaPars, hyperValues[hypers]), value = total, gradient = NULL, hessian = NULL,
       information = NULL, etaModes = modes, etaCov = setNames(cov, design$groups), clusters = clusters,
       spike = spike, converged = conv, iterations = length(st$rows), stopReason = stopReason,
       stageTrace = do.call(rbind, st$rows))
}
