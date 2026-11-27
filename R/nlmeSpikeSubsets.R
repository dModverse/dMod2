# Spike-and-slab random effects with a given slab mean and spike, -Inf included, beside normal ones: per
# group the sum over the subsets of slab etas, each a Laplace integral or a sparse-grid quadrature (methods
# article Sparsity priors, section Spike at a given value).


## The most Laplace integrals per evaluation, summed over the groups.
.spikeBudget <- 4096L

## The approximations of the sum over subsets that a phase takes.
.spikeApprox <- c("meanfield", "ec", "exact")

## TRUE for one number of at least `from`, and for `.isCount()` a whole one.
.isNumberFrom <- function(x, from) is.numeric(x) && length(x) == 1L && !is.na(x) && x >= from
.isCount <- function(x) .isNumberFrom(x, 1) && x == round(x)

## The burn-in of control$sparseFit: NULL without, else `start` and `iterlim`.
.spikeBurninOf <- function(sf) {
  b <- sf$burnin %||% FALSE
  if (isFALSE(b)) return(NULL)
  if (isTRUE(b)) b <- list()
  start <- if (is.list(b)) b$start %||% "random"
  iterlim <- if (is.list(b)) b$iterlim %||% 200L
  ok <- is.list(b) && all(names(b) %in% c("start", "iterlim")) && is.character(start) && !anyNA(start) &&
    .isCount(iterlim)
  if (!ok)
    stop("em: `control$sparseFit$burnin` is FALSE, TRUE or a list of `start` (\"random\", \"spike\", \"slab\" or ",
         "names of spike-and-slab random effects) and `iterlim` (a whole number of at least 1).", call. = FALSE)
  list(start = start, iterlim = as.integer(iterlim))
}

## The phases of control$sparseFit in the order they run, each with `approx`, `maxECM` and `epsECM`; epsECM
## defaults to 1e-2, in the last phase to 1e-6.
.spikePhases <- function(sf) {
  ph <- sf$phases %||% .spikeApprox
  if (is.character(ph)) ph <- setNames(rep(list(list()), length(ph)), ph)
  entryOk <- function(p) is.list(p) && all(names(p) %in% c("maxECM", "epsECM")) &&
    (is.null(p$maxECM) || .isCount(p$maxECM)) && (is.null(p$epsECM) || .isNumberFrom(p$epsECM, 0))
  ok <- is.list(ph) && length(ph) > 0L && !is.null(names(ph)) && all(names(ph) %in% .spikeApprox) &&
    !anyDuplicated(names(ph)) &&
    all(vapply(ph, entryOk, TRUE))
  if (!ok)
    stop("em: `control$sparseFit$phases` names each of \"meanfield\", \"ec\" and \"exact\" at most once, as a ",
         "character vector or a list of lists with `maxECM` (a whole number of at least 1) and `epsECM`.",
         call. = FALSE)
  last <- length(ph)
  Map(function(approx, p, i) list(approx = approx, maxECM = as.integer(p$maxECM %||% 200L),
                                  epsECM = p$epsECM %||% if (i == last) 1e-6 else 1e-2),
      names(ph), ph, seq_along(ph))
}

## The spike-and-slab etas in the slab during the burn-in from `start`: each with the probability `prob` of the
## stem, none, every one, or the given ones; the others stay at the spike.
.spikeBurninSet <- function(start, design, lay, prob) {
  E <- design$etaNames[, lay$stem, drop = FALSE]
  etas <- as.vector(E)
  if (identical(start, "random")) return(etas[stats::runif(length(etas)) < prob[rep(lay$stem, each = nrow(E))]])
  if (identical(start, "spike")) return(character(0))
  if (identical(start, "slab")) return(etas)
  bad <- .setdiff(start, etas)
  if (length(bad))
    stop("em: `control$sparseFit$burnin$start` names no spike-and-slab random effect: ", paste(bad, collapse = ", "),
         ".", call. = FALSE)
  start
}

## The error `e` of phase `i` with approximation `approx` as a condition with the phase in field `phase`.
.spikePhaseError <- function(e, i, approx, several) {
  msg <- conditionMessage(e)
  if (several) msg <- sprintf("em: phase %d (%s): %s", i, approx, sub("^em: ", "", msg))
  structure(class = c("dMod2PhaseError", "error", "condition"), list(message = msg, call = NULL, phase = i))
}

## Per stem of the spike-and-slab blocks: block, slab mean, spike, the hyperparameter names and the shapes of a
## beta prior on `prob`, NA without.
.spikeLayout <- function(design) {
  blocks <- design$blocks
  do.call(rbind, lapply(which(vapply(blocks, function(b) b$dist == "spikeSlab", TRUE)), function(i) {
    b <- blocks[[i]]
    data.frame(stem = b$stems, block = i, mean = b$mean, spike = b$spike, prob = b$prob, tau = b$tau,
               probName = b$hyper[["prob"]], tauName = b$hyper[["tau"]], betaA = b$probPrior$a %||% NA_real_,
               betaB = b$probPrior$b %||% NA_real_, stringsAsFactors = FALSE)
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
## covariance 2 H^-1; NULL where the inner search fails. With `laplace = FALSE` the value of data plus prior
## at the mode and the mode alone.
.spikeSubset <- function(design, objs, g, theta, fixed, etas, on, pr, start, laplace = TRUE) {
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
  # a failed inner search stops the fit: every subset stays in the sum over subsets
  fit <- tryCatch(.innerTrust(list(objfun = objfun, parinit = start[on], fixed = fx, rinit = 1, rmax = 10,
                                   iterlim = 100L)),
                  error = function(e) stop("em: the inner search of group ", g, " with the slab etas ",
                                           paste(on, collapse = ", "), " and the others at the spike failed: ",
                                           conditionMessage(e), call. = FALSE))
  if (!is.finite(fit$value)) return(NULL)
  if (!laplace) return(list(m2logZ = fit$value, mode = fit$argument))
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

## The E-step of every group by the approximation `approx`: per subset of the slab etas that reach the group
## data the weight, mode and covariance; -2 log I per group; the probability of the slab per eta; `ref`, the
## candidates of inclusion above 1/2. "exact" sums exactly over the subsets of the groups within `termBudget`
## subsets, the smallest first, and takes EC for the others; "burnin" takes the one subset `ref` at the joint
## mode, without integral. EC and mean field take the subset of inclusion above 1/2 as their one part, and
## expand around `ref` of the step before.
.spikeEStep <- function(design, objs, lay, theta, fixed, hyp, warm, bases, termBudget, control = list(),
                        quadrature = NULL, gridWarm = list(), omega = NULL, approx = "exact", ref = list()) {
  E <- design$etaNames
  reach <- .reachOf(design)
  stems <- lay$stem
  nstems <- design$normal$stems
  lbStem <- setNames(ifelse(stems %in% names(bases), log(bases[stems]), NA_real_), stems)
  sizes <- setNames(rowSums(reach[design$groups, stems, drop = FALSE]), design$groups)
  exact <- if (approx == "exact") design$groups[order(sizes)][cumsum(2^sort(sizes)) <= termBudget] else character(0)
  betaBlocks <- unique(lay$block[!is.na(lay$betaA)])
  out <- lapply(setNames(nm = design$groups), function(g) {
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
    refG <- .intersect(ref[[g]] %||% cand, cand)
    keyOf <- function(on) paste(on, collapse = ",")
    startOf <- function(on) {
      v <- alpha
      v[nrm] <- 0
      v[on] <- warm[[g]][[keyOf(on)]] %||% c(mu, setNames(rep(0, length(nrm)), nrm))[on]
      v
    }
    group <- list(slab = slab, normal = nrm, etas = etas, cand = cand, mu = mu, alpha = alpha)
    if (approx == "burnin") {
      on <- c(nrm, refG)
      r <- .spikeSubset(design, objs, g, theta, fixed, etas, on, pr, startOf(on), laplace = FALSE)
      if (is.null(r)) stop("em: the joint mode of group ", g, " is not finite.", call. = FALSE)
      r$on <- on
      r$weight <- 1
      return(c(list(parts = list(r), m2logI = r$m2logZ, inclusion = probG, ref = refG), group))
    }
    if (approx == "meanfield")
      return(c(.spikeMeanField(design, objs, g, theta, fixed, etas, cand, nrm, pr, probG, refG, startOf), group))
    if (!g %in% exact) {
      # the data term first expanded at the mode with the candidates of `ref` in the slab
      start <- startOf(c(nrm, refG))
      full <- .spikeSubset(design, objs, g, theta, fixed, etas, c(nrm, refG), pr, start)
      if (!is.null(full)) start[c(nrm, refG)] <- full$mode
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
      return(c(list(parts = list(r), m2logI = ec$m2logI, inclusion = inclusion, ref = ec$on,
                    factor = list(base = -ec$m2logI / 2 - sum(ec$logZsite), w0 = ec$massSpike, w1 = ec$massSlab,
                                  logZ = ec$logZsite),
                    perCandidate = list(inclusion = ec$inclusion, dev2 = ec$dev2, normalMoment = ec$normalMoment)),
               group))
    }
    S <- .subsets(length(cand))
    # candidates under a beta prior on `prob` enter by their count per block, the others by the product
    betaBlock <- setNames(ifelse(is.na(lay$betaA), NA_integer_, lay$block), slab)[cand]
    prodC <- is.na(betaBlock)
    parts <- lapply(seq_len(nrow(S)), function(i) {
      on <- c(nrm, cand[S[i, ]])
      r <- .spikeSubset(design, objs, g, theta, fixed, etas, on, pr, startOf(on))
      if (is.null(r)) return(NULL)
      r$on <- on
      r$logw <- sum(log(probG[cand[S[i, ] & prodC]])) + sum(log1p(-probG[cand[!S[i, ] & prodC]])) - r$m2logZ / 2
      r$count <- vapply(betaBlocks, function(b) sum(S[i, ] & betaBlock %in% b), 0L)
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
    # with a beta prior the groups of a block are summed jointly, by .spikeBetaSum()
    if (length(betaBlocks)) return(c(list(parts = parts, probG = probG), group))
    .spikeWeights(parts, -2 * .logSumExp(lw), probG, group)
  })
  if (!length(betaBlocks)) return(out)
  factorized <- names(Filter(function(e) !is.null(e$factor), out))
  if (length(factorized)) out <- .spikeBetaFactorized(out, lay, factorized)
  if (approx != "exact" || !length(exact)) return(out)
  .spikeBetaSum(out, lay, exact)
}

## The weights of the subsets `parts` of a group from their log weight and -2 log I, `m2logI`; the inclusion per
## slab eta, the prior `probG` for an eta off the group data.
.spikeWeights <- function(parts, m2logI, probG, group) {
  for (i in seq_along(parts)) parts[[i]]$weight <- exp(parts[[i]]$logw + m2logI / 2)
  inclusion <- setNames(vapply(group$slab, function(e) {
    if (!e %in% group$cand) return(probG[[e]])
    sum(vapply(parts, function(p) if (e %in% p$on) p$weight else 0, 0))
  }, 0), group$slab)
  c(list(parts = parts, m2logI = m2logI, inclusion = inclusion, ref = group$cand[inclusion[group$cand] > 0.5]),
    group)
}

## The log of the sum of exp(`x`).
.logSumExp <- function(x) {
  top <- max(x)
  if (!is.finite(top)) return(top)
  top + log(sum(exp(x - top)))
}

## Log weights by count vector, named "n1,n2,..." over the beta blocks: the table of the log weights of the sum
## of the counts.
.convolveCounts <- function(t1, t2) {
  out <- list()
  for (k1 in names(t1)) for (k2 in names(t2)) {
    k <- paste(.countsOf(k1) + .countsOf(k2), collapse = ",")
    out[[k]] <- c(out[[k]], t1[[k1]] + t2[[k2]])
  }
  vapply(out, .logSumExp, 0)
}

.countsOf <- function(key) as.integer(strsplit(key, ",", fixed = TRUE)[[1]])

## The exact sum over the subsets of the groups `exact` in `est` (from .spikeEStep()) with `prob` under a beta
## prior integrated out: per beta block the weight B(a + n, b + K - n) / B(a, b) of the count n of slab etas over
## all groups, K the etas of the block that reach data. The groups of a block are summed jointly by convolving
## their sums by count; -2 log I of the joint sum is shared equally by the groups.
.spikeBetaSum <- function(est, lay, exact) {
  blocks <- unique(lay$block[!is.na(lay$betaA)])
  a <- lay$betaA[match(blocks, lay$block)]
  b <- lay$betaB[match(blocks, lay$block)]
  # the slab etas of a group follow the order of the stems in `lay`
  K <- vapply(blocks, function(bl) sum(vapply(est[exact], function(e)
    sum(lay$block[match(e$cand, e$slab)] == bl), 0L)), 0L)
  logW <- function(key) {
    n <- .countsOf(key)
    sum(lbeta(a + n, b + K - n) - lbeta(a, b))
  }
  tables <- lapply(est[exact], function(e) {
    keys <- vapply(e$parts, function(p) paste(p$count, collapse = ","), "")
    vapply(split(vapply(e$parts, `[[`, 0, "logw"), keys), .logSumExp, 0)
  })
  zero <- setNames(0, paste(integer(length(blocks)), collapse = ","))
  total <- Reduce(.convolveCounts, tables, zero)
  logI <- .logSumExp(total + vapply(names(total), logW, 0))
  for (g in exact) {
    others <- Reduce(.convolveCounts, tables[.setdiff(exact, g)], zero)
    e <- est[[g]]
    for (i in seq_along(e$parts)) {
      own <- e$parts[[i]]$count
      joint <- vapply(names(others), function(k) logW(paste(own + .countsOf(k), collapse = ",")), 0)
      e$parts[[i]]$logw <- e$parts[[i]]$logw + .logSumExp(others + joint)
    }
    est[[g]] <- .spikeWeights(e$parts, -2 * logI, e$probG, e[c("slab", "normal", "etas", "cand", "mu", "alpha")])
    est[[g]]$m2logI <- -2 * logI / length(exact)
  }
  est
}

## Mean field over the inclusion of the candidates `cand` of group `g`: log Z of a subset linear in the
## indicators around the reference subset `ref`, from the Laplace integrals of the reference and of every
## subset that differs from it in one candidate, so that the sum over subsets factorizes per candidate. The
## reference moves to the candidates of inclusion above 1/2 until it stays; a subset whose inner search fails
## has weight 0 and is counted in `failed`.
.spikeMeanField <- function(design, objs, g, theta, fixed, etas, cand, nrm, pr, prob, ref, startOf) {
  seen <- new.env()
  # the key of the empty subset is "s", not "", which an environment refuses
  keyOf <- function(on) paste0("s", paste(sort(on), collapse = ","))
  laplaceOf <- function(on) {
    on <- cand[cand %in% on]
    key <- keyOf(on)
    if (!is.null(seen[[key]])) return(seen[[key]]$r)
    r <- tryCatch(.spikeSubset(design, objs, g, theta, fixed, etas, c(nrm, on), pr, startOf(c(nrm, on))),
                  error = function(e) NULL)
    if (!is.null(r)) r$on <- c(nrm, on)
    assign(key, list(r = r), envir = seen)
    r
  }
  visited <- character(0)
  repeat {
    base <- laplaceOf(ref)
    if (is.null(base)) stop("em: the inner search of group ", g, " at the reference subset failed.", call. = FALSE)
    flips <- lapply(cand, function(j) laplaceOf(if (j %in% ref) .setdiff(ref, j) else c(ref, j)))
    # log odds of the slab per candidate: the subsets with and without it, the others as in the reference
    flip <- vapply(flips, function(r) if (is.null(r)) Inf else r$m2logZ, 0)
    inRef <- cand %in% ref
    m2On <- ifelse(inRef, base$m2logZ, flip)
    m2Off <- ifelse(inRef, flip, base$m2logZ)
    logOdds <- log(prob[cand]) - log1p(-prob[cand]) + (m2Off - m2On) / 2
    inclusion <- stats::plogis(logOdds)
    visited <- c(visited, keyOf(ref))
    nxt <- cand[inclusion > 0.5]
    if (setequal(nxt, ref) || keyOf(nxt) %in% visited || length(visited) > length(cand)) break
    ref <- nxt
  }
  failed <- sum(vapply(as.list(seen), function(x) is.null(x$r), TRUE))
  # log I: prior and integral of the reference, times 1 + exp(flip log odds) per candidate
  flipLogit <- ifelse(inRef, -logOdds, logOdds)
  log1pExp <- ifelse(flipLogit > 0, flipLogit + log1p(exp(-flipLogit)), log1p(exp(flipLogit)))
  logPriorRef <- sum(log(prob[cand[inRef]])) + sum(log1p(-prob[cand[!inRef]]))
  logI <- logPriorRef - base$m2logZ / 2 + sum(log1pExp)
  dev2 <- setNames(vapply(seq_along(cand), function(k) {
    r <- if (inRef[k]) base else flips[[k]]
    if (is.null(r)) return(0)
    j <- match(cand[k], r$on)
    (r$mode[[cand[k]]] - pr$mu[[cand[k]]])^2 + r$cov[j, j]
  }, 0), cand)
  q <- match(nrm, base$on)
  normalMoment <- if (length(nrm)) tcrossprod(base$mode[nrm]) + base$cov[q, q, drop = FALSE]
  base$weight <- 1
  inc <- prob
  inc[cand] <- inclusion
  # log Z of a subset differs from that of the reference by `delta` per candidate on the other side
  delta <- ifelse(is.finite(flip), (base$m2logZ - flip) / 2, -Inf)
  w1 <- setNames(ifelse(inRef, 0, delta), cand)
  w0 <- setNames(ifelse(inRef, delta, 0), cand)
  logZ <- .logSumExp2(log(prob[cand]) + w1, log1p(-prob[cand]) + w0)
  factor <- list(base = -base$m2logZ / 2, w0 = w0, w1 = w1, logZ = logZ)
  list(parts = list(base), m2logI = -2 * logI, inclusion = inc, ref = ref, failed = failed, factor = factor,
       perCandidate = list(inclusion = setNames(inclusion, cand), dev2 = dev2, normalMoment = normalMoment))
}

## Mean field or EC over the groups `groups` of `est` with `prob` under a beta prior integrated out. Both factorize
## over the candidates: per group `factor` holds `base`, log I without the candidates, plus per candidate the log
## masses `w0` of the spike and `w1` of the slab and `logZ` of their mixture under the prior, so that log I is
## base + sum(logZ). The candidates under a beta prior take the weight of .spikeBetaSum() by convolving one
## count table per candidate in place of their mixtures; -2 log I is shared equally by the groups.
.spikeBetaFactorized <- function(est, lay, groups) {
  blocks <- unique(lay$block[!is.na(lay$betaA)])
  a <- lay$betaA[match(blocks, lay$block)]
  b <- lay$betaB[match(blocks, lay$block)]
  units <- list(); base <- 0
  for (g in groups) {
    e <- est[[g]]
    f <- e$factor
    block <- lay$block[match(e$cand, e$slab)]
    beta <- !is.na(lay$betaA[match(e$cand, e$slab)])
    base <- base + f$base + sum(f$logZ[e$cand][!beta])
    for (k in which(beta)) {
      unit <- as.integer(blocks == block[k])
      eta <- e$cand[k]
      units[[length(units) + 1L]] <- list(g = g, eta = eta, unit = unit, w1 = f$w1[[eta]],
        table = setNames(c(f$w0[[eta]], f$w1[[eta]]), c(paste(0L * unit, collapse = ","), paste(unit, collapse = ","))))
    }
  }
  K <- vapply(blocks, function(bl) sum(vapply(units, function(u) u$unit[blocks == bl], 0L)), 0L)
  logW <- function(keys) vapply(keys, function(key) {
    n <- .countsOf(key)
    sum(lbeta(a + n, b + K - n) - lbeta(a, b))
  }, 0)
  zero <- setNames(0, paste(integer(length(blocks)), collapse = ","))
  tables <- lapply(units, `[[`, "table")
  total <- Reduce(.convolveCounts, tables, zero)
  logBeta <- .logSumExp(total + logW(names(total)))
  for (i in seq_along(units)) {
    u <- units[[i]]
    rest <- Reduce(.convolveCounts, tables[-i], zero)
    on <- vapply(names(rest), function(k) paste(.countsOf(k) + u$unit, collapse = ","), "")
    incl <- exp(u$w1 + .logSumExp(rest + logW(on)) - logBeta)
    est[[u$g]]$inclusion[[u$eta]] <- incl
    est[[u$g]]$perCandidate$inclusion[[u$eta]] <- incl
  }
  for (g in groups) {
    est[[g]]$m2logI <- -2 * (base + logBeta) / length(groups)
    est[[g]]$ref <- est[[g]]$cand[est[[g]]$inclusion[est[[g]]$cand] > 0.5]
  }
  est
}

## Elementwise log(exp(x) + exp(y)).
.logSumExp2 <- function(x, y) {
  top <- pmax(x, y)
  ifelse(is.finite(top), top + log(exp(x - top) + exp(y - top)), top)
}

## E[eta eta'] of the normal etas of every group under the Laplace posteriors of the subsets with their
## weights; EC and mean-field groups from their per-candidate moments.
.spikeNormalMoments <- function(est) {
  lapply(est, function(e) {
    q <- e$normal
    if (!is.null(e$perCandidate)) return(e$perCandidate$normalMoment)
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

## Under a beta prior on `prob` of block `b`, the shapes of Beta(a + n, b + K - n), n the expected number of slab
## etas among the K etas of the block that reach data; NULL without a beta prior.
.spikeProbPosterior <- function(lay, est, b) {
  rows <- lay$block == b
  if (is.na(lay$betaA[rows][1])) return(NULL)
  n <- K <- 0
  for (e in est) {
    on <- .intersect(e$slab[rows], e$cand)
    n <- n + sum(e$inclusion[on])
    K <- K + length(on)
  }
  c(a = lay$betaA[rows][1] + n, b = lay$betaB[rows][1] + K - n)
}

## The M-step of the hyperparameters: prob as the mean inclusion, under a beta prior the predictive
## probability of the slab (a + n_-k) / (a + b + K - 1) that EC and the burn-in take, tau^2 as the mean of
## (eta - mean)^2 under the slab, under a scale prior the posterior mode of .spikeTauMode(), per block;
## `tauPriors` by block; fixed values stay.
.spikeHyperStep <- function(lay, est, hyp, hypers, tauPriors = list()) {
  stems <- lay$stem
  for (b in unique(lay$block)) {
    st <- lay$stem[lay$block == b]
    k <- match(st, stems)
    incl <- unlist(lapply(est, function(e) e$inclusion[k]))
    if (lay$probName[lay$block == b][1] %in% hypers) hyp$prob[st] <- min(max(mean(incl), 1e-6), 1 - 1e-6)
    # under a beta prior, the predictive probability of the slab given the expected count of the others
    post <- .spikeProbPosterior(lay, est, b)
    if (!is.null(post))
      hyp$prob[st] <- (post[["a"]] - mean(incl, na.rm = TRUE)) / (sum(post) - 1)
    if (!lay$tauName[lay$block == b][1] %in% hypers) next
    num <- den <- 0
    for (e in est) {
      own <- e$slab[k]
      pc <- e$perCandidate
      if (!is.null(pc)) {
        on <- .intersect(names(pc$inclusion), own)
        num <- num + sum(pc$inclusion[on] * pc$dev2[on])
        den <- den + sum(pc$inclusion[on])
      }
      for (p in if (is.null(pc)) e$parts) {
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
    if (den > 0) hyp$tau[st] <- .spikeTauMode(num, den, tauPriors[[b]], hyp$tau[st[1]])
  }
  hyp
}

## The tau that minimizes den log tau^2 + num / tau^2 - 2 log p(tau): without a prior sqrt(num / den); under an
## inverse-gamma prior on tau^2 in closed form; under a half-Cauchy or half-normal prior on tau by Newton steps
## in log tau, from `tau0`, where the objective is convex.
.spikeTauMode <- function(num, den, prior = NULL, tau0 = 1) {
  if (is.null(prior)) return(sqrt(num / den))
  if (prior$dist == "invGamma") return(sqrt((num + 2 * prior$rate) / (den + 2 * prior$shape + 2)))
  if (num <= 0) return(tau0)
  u <- log(unname(tau0))
  for (i in seq_len(100L)) {
    pen <- .spikeTauPenalty(exp(u), prior, derivs = TRUE)
    d1 <- 2 * den - 2 * num * exp(-2 * u) + pen[["d1"]]
    d2 <- 4 * num * exp(-2 * u) + pen[["d2"]]
    step <- max(min(-d1 / d2, 1), -1)
    u <- u + step
    if (abs(step) < 1e-12) break
  }
  exp(u)
}

## -2 log p(tau) of a scale prior, the density in tau for half-Cauchy and half-normal and in tau^2 for the
## inverse gamma; with `derivs` also the first and second derivative in log tau.
.spikeTauPenalty <- function(tau, prior, derivs = FALSE) {
  if (is.null(prior)) return(if (derivs) c(value = 0, d1 = 0, d2 = 0) else 0)
  v <- unname(tau)^2
  out <- switch(prior$dist,
    halfCauchy = {
      r <- v / prior$scale^2
      c(value = -2 * log(2 / (pi * prior$scale)) + 2 * log1p(r), d1 = 4 * r / (1 + r), d2 = 8 * r / (1 + r)^2)
    },
    halfNormal = {
      r <- v / prior$scale^2
      c(value = -2 * log(2 / (prior$scale * sqrt(2 * pi))) + r, d1 = 2 * r, d2 = 4 * r)
    },
    invGamma = c(value = -2 * (prior$shape * log(prior$rate) - lgamma(prior$shape)) + 2 * (prior$shape + 1) * log(v) +
                   2 * prior$rate / v,
                 d1 = 4 * (prior$shape + 1) - 4 * prior$rate / v, d2 = 8 * prior$rate / v))
  if (derivs) out else out[["value"]]
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
      # a part of a node of the integral over tau has the slab width of the node
      prPart <- pr
      if (!is.null(p$tau)) prPart$tau[match(names(p$tau), lay$stem)] <- p$tau
      items[[k]] <- list(g = g, on = p$on, copies = copies, w = p$weight, mode = p$mode,
                         off = e$alpha[.setdiff(e$etas, p$on)], pr = prPart)
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
.spikeFrozen <- function(design, est, lay, theta, omegaPars, hyp, hypers, fixed, extra, tauPriors = list()) {
  # a beta prior on prob weights the count of slab etas of a block, here within one group
  betaBlocks <- unique(lay$block[!is.na(lay$betaA)])
  for (bl in betaBlocks) {
    holders <- Filter(function(e) any(lay$block[match(e$cand, e$slab)] == bl), est)
    if (length(holders) > 1L)
      stop("em: method = \"quadrature\" takes betaPrior() for spike-and-slab random effects that reach the data ",
           "of one group.", call. = FALSE)
  }
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
      betaC <- inCand & !is.na(lay$betaA)
      prodC <- inCand & is.na(lay$betaA)
      K <- vapply(betaBlocks, function(bl) sum(betaC & lay$block == bl), 0L)
      a <- lay$betaA[match(betaBlocks, lay$block)]
      bb <- lay$betaB[match(betaBlocks, lay$block)]
      for (b in seq_len(B)) {
        on <- .intersect(G$subset[[b]], e$slab)
        st <- lay$stem[match(on, e$slab)]
        tau <- h$tau[st]
        dev2 <- ((G$eta[b, on] - e$mu[on]) / tau)^2
        z <- e$slab %in% G$subset[[b]]
        count <- vapply(betaBlocks, function(bl) sum(z & betaC & lay$block == bl), 0L)
        logInt[b] <- logInt[b] + sum(-0.5 * log(2 * pi) - log(tau) - dev2 / 2) +
          sum(ifelse(z, log(prob), log1p(-prob))[prodC]) + sum(lbeta(a + count, bb + K - count) - lbeta(a, bb))
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
    # -2 log p(tau) of the scale priors on the log tau coordinates
    for (n in tauFree) {
      prior <- tauPriors[[lay$block[match(n, lay$tauName)]]]
      if (is.null(prior)) next
      pen <- .spikeTauPenalty(h$tau[[lay$stem[match(n, lay$tauName)]]], prior, derivs = TRUE)
      total <- total + pen[["value"]]
      if (deriv && n %in% act) { gr[n] <- gr[n] + pen[["d1"]]; H[n, n] <- H[n, n] + pen[["d2"]] }
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
## subsets: EM on the Laplace integrals in the phases of control$sparseFit, each from the estimates and inner
## modes of the one before, then with `quadrature` (control$quadrature) trust steps on frozen quadrature nodes.
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
  beta <- !is.na(lay$betaA)
  st$hyp$prob[beta] <- lay$betaA[beta] / (lay$betaA[beta] + lay$betaB[beta])
  # the scale priors of the blocks whose tau is estimated, and -2 log p(tau) of them at the current tau
  tauPriors <- lapply(design$blocks, function(b) if (isTRUE(b$hyper["tau"] %in% hypers)) b$tauPrior)
  # a block whose tau the phase "exact" integrates over log tau, on a Gauss-Hermite grid of mean and sd in log tau
  integrated <- Filter(function(b) isTRUE(tauPriors[[b]]$integrate), unique(lay$block))
  if (length(integrated) > 1L)
    stop("em: `integrate = TRUE` applies to the tau of one spike-and-slab block.", call. = FALSE)
  if (length(integrated) && !is.null(quadrature))
    stop("em: method = \"quadrature\" takes the posterior mode of tau; drop `integrate = TRUE`.", call. = FALSE)
  intStems <- lay$stem[lay$block %in% integrated]
  st$tauGrid <- if (length(integrated)) c(mean = log(st$hyp$tau[[intStems[1]]]), sd = 0.5)
  st$integratedNow <- FALSE
  # -2 log p(tau) of the scale priors, but of the tau the current E-step integrated
  penaltyOf <- function() sum(vapply(unique(lay$block), function(b)
    if (st$integratedNow && b %in% integrated) 0
    else .spikeTauPenalty(st$hyp$tau[lay$stem[lay$block == b][1]], tauPriors[[b]]), 0))
  st$warm <- list(); st$gridWarm <- list(); st$rows <- list(); st$ref <- list()
  sf <- control$sparseFit
  termBudget <- sf$termBudget %||% .spikeBudget
  phases <- .spikePhases(sf)
  bases <- .spikeBases(design)
  burn <- .spikeBurninOf(sf)
  burnSet <- if (!is.null(burn)) .spikeBurninSet(burn$start, design, lay, st$hyp$prob)
  if (!is.null(burnSet))
    st$ref <- lapply(setNames(nm = design$groups),
                     function(g) .intersect(unname(design$etaNames[g, lay$stem]), burnSet))
  omegaOf <- function() .spikeOmega(design, c(st$omegaPars, fixed))
  keyOf <- function(e) vapply(e$parts, function(p) paste(p$on, collapse = ","), "")
  extraOf <- function() (if (is.null(extra)) 0 else extra(st$theta, fixed = fixed, deriv = FALSE)$value) + penaltyOf()
  eStep <- function(qc, approx) {
    if (approx == "exact" && length(integrated)) return(eStepTau(qc))
    eStepAt(qc, approx, st$hyp)
  }
  eStepAt <- function(qc, approx, hyp) {
    st$integratedNow <- FALSE
    est <- .spikeEStep(design, objs, lay, st$theta, fixed, hyp, st$warm, bases, termBudget, sf$ec, qc,
                       st$gridWarm, omegaOf(), approx, st$ref)
    # the inner modes of every subset met so far start the inner searches of later steps and phases
    for (g in names(est)) for (p in est[[g]]$parts)
      if (length(p$on)) st$warm[[g]][[paste(p$on, collapse = ",")]] <- p$mode
    st$gridWarm <- lapply(est, function(e) setNames(lapply(e$parts, `[[`, "indices"), keyOf(e)))
    st$ref <- lapply(est, `[[`, "ref")
    est
  }
  # the exact sum integrated over log tau by 7 Gauss-Hermite nodes on st$tauGrid: the parts of every node weighted
  # by the posterior of the node, which moves the grid to the posterior mean and sd of log tau
  eStepTau <- function(qc) {
    gh <- .gaussHermite(7L)
    prior <- tauPriors[[integrated]]
    u <- st$tauGrid[["mean"]] + sqrt(2) * st$tauGrid[["sd"]] * gh$x
    nodes <- lapply(u, function(ui) {
      h <- st$hyp
      h$tau[intStems] <- exp(ui)
      eStepAt(qc, "exact", h)
    })
    m2 <- vapply(nodes, function(est) sum(vapply(est, `[[`, 0, "m2logI")), 0)
    # the prior density in log tau: of tau, or of tau^2 for the inverse gamma
    logJacobian <- if (prior$dist == "invGamma") log(2) + 2 * u else u
    lw <- -(m2 + vapply(exp(u), .spikeTauPenalty, 0, prior = prior)) / 2 + logJacobian + log(gh$w) + gh$x^2
    logI <- .logSumExp(lw) + log(sqrt(2) * st$tauGrid[["sd"]])
    w <- exp(lw - .logSumExp(lw))
    mean <- sum(w * u)
    st$tauGrid <- c(mean = mean, sd = max(sqrt(sum(w * (u - mean)^2)), 1e-3))
    st$hyp$tau[intStems] <- exp(mean)
    merged <- lapply(setNames(nm = names(nodes[[1]])), function(g) {
      e <- nodes[[1]][[g]]
      e$parts <- unlist(Map(function(est, wi, ui) lapply(est[[g]]$parts, function(p) {
        p$weight <- p$weight * wi
        p$tau <- setNames(rep(exp(ui), length(intStems)), intStems)
        p
      }), nodes, w, u), recursive = FALSE)
      e$inclusion <- Reduce(`+`, Map(function(est, wi) wi * est[[g]]$inclusion, nodes, w))
      e$m2logI <- -2 * logI / length(nodes[[1]])
      e$ref <- e$cand[e$inclusion[e$cand] > 0.5]
      e
    })
    st$ref <- lapply(merged, `[[`, "ref")
    st$integratedNow <- TRUE
    merged
  }
  record <- function(stage, phase, k, value, delta) {
    st$rows[[length(st$rows) + 1L]] <- data.frame(stage = stage, phase = phase, ecmIter = k, OFV = value,
                                                  delta = delta)
    if (verbose) message(sprintf("em(spike subsets) %s %d: -2 log L %.6f, change %.2e", stage, k, value, delta))
  }

  # the phase "exact" beyond termBudget ends the fit at EC, without EC it stops
  approxes <- vapply(phases, `[[`, "", "approx")
  nSubsets <- sum(2^rowSums(.reachOf(design)[design$groups, lay$stem, drop = FALSE]))
  skipExact <- nSubsets > termBudget
  if (skipExact && "exact" %in% approxes && !"ec" %in% approxes)
    stop("em: the exact sum runs over ", format(nSubsets), " subsets, more than `termBudget` (", format(termBudget),
         "); add \"ec\" to `control$sparseFit$phases` or raise `termBudget`.", call. = FALSE)
  phaseRows <- list()
  several <- length(phases) + as.integer(!is.null(burn)) > 1L
  # burn-in: one trust region over the fixed effects and the slab etas of the start subset at the joint mode,
  # the others at the spike, prob, tau and Omega held
  if (!is.null(burn)) {
    t0 <- proc.time()[["elapsed"]]
    m <- tryCatch({
      est <- eStep(NULL, "burnin")
      m <- list(value = sum(vapply(est, `[[`, 0, "m2logI")) + extraOf(), iterations = 0L, converged = TRUE)
      if (length(st$theta)) {
        ms <- .spikeM(design, objs, est, extra, lay, st$hyp, omegaOf())
        m <- .innerTrust(list(objfun = ms$objfun, parinit = c(st$theta, ms$init), fixed = fixed, rinit = 1,
                              rmax = 10, iterlim = burn$iterlim))
        st$theta <- m$argument[names(st$theta)]
        # .spikeM() names the copy of the etas of item k `.s<k>.<eta>`; the burn-in has one item per group
        for (k in seq_along(est)) {
          on <- est[[k]]$parts[[1]]$on
          if (length(on))
            st$warm[[names(est)[k]]][[paste(on, collapse = ",")]] <- setNames(m$argument[paste0(".s", k, ".", on)], on)
        }
      }
      m
    }, error = function(e) stop(.spikePhaseError(e, 1L, "burnin", several)))
    record("burnin", 1L, m$iterations, m$value, NA_real_)
    phaseRows[[1L]] <- data.frame(phase = 1L, approx = "burnin", iterations = as.integer(m$iterations),
                                  value = m$value, converged = isTRUE(m$converged),
                                  seconds = proc.time()[["elapsed"]] - t0, failed = 0L)
  }

  # Laplace stage: per phase ECM until no estimate changes by more than epsECM
  conv <- FALSE
  for (j in seq_along(phases)) {
    ph <- phases[[j]]
    i <- j + as.integer(!is.null(burn))
    if (ph$approx == "exact" && skipExact) {
      phaseRows[[i]] <- data.frame(phase = i, approx = "exact", iterations = 0L, value = NA_real_, converged = NA,
                                   seconds = 0, failed = 0L)
      next
    }
    t0 <- proc.time()[["elapsed"]]
    conv <- FALSE
    k <- 0L
    last <- NA_real_
    failed <- 0L
    tryCatch(for (k in seq_len(ph$maxECM)) {
      est <- eStep(NULL, ph$approx)
      failed <- failed + sum(vapply(est, function(e) e$failed %||% 0L, 0L))
      before <- c(st$theta, unlist(st$hyp), st$omegaPars)
      held <- if (st$integratedNow) lay$tauName[match(integrated, lay$block)]
      st$hyp <- .spikeHyperStep(lay, est, st$hyp, .setdiff(hypers, held), tauPriors)
      if (length(st$omegaPars))
        st$omegaPars[] <- .updateNormal(design, .spikeNormalMoments(est), fixed)[names(st$omegaPars)]
      if (length(st$theta)) {
        ms <- .spikeM(design, objs, est, extra, lay, st$hyp, omegaOf())
        m <- .innerTrust(list(objfun = ms$objfun, parinit = c(st$theta, ms$init), fixed = fixed, rinit = 1,
                              rmax = 10, iterlim = 100L))
        st$theta <- m$argument[names(st$theta)]
      }
      delta <- max(abs(c(st$theta, unlist(st$hyp), st$omegaPars) - before))
      last <- sum(vapply(est, `[[`, 0, "m2logI"))
      record("laplace", i, k, last, delta)
      conv <- delta < ph$epsECM
      if (conv) break
    }, error = function(e) stop(.spikePhaseError(e, i, ph$approx, several)))
    phaseRows[[i]] <- data.frame(phase = i, approx = ph$approx, iterations = k, value = last + extraOf(),
                                 converged = conv, seconds = proc.time()[["elapsed"]] - t0, failed = failed)
  }
  stopReason <- if (conv) "converged" else "iterlim"
  final <- utils::tail(.setdiff(approxes, if (skipExact) "exact"), 1L)

  # quadrature stage: a trust step over every parameter on the frozen nodes, until a Newton step on the
  # rebuilt nodes gains less than tolDecrement
  if (!is.null(quadrature)) {
    conv <- FALSE; stopReason <- "iterlim"
    for (k in seq_len(quadrature$maxECMPerStage %||% 10L)) {
      est <- eStep(quadrature, final)
      fr <- .spikeFrozen(design, est, lay, st$theta, st$omegaPars, st$hyp, hypers, fixed, extra, tauPriors)
      ev <- fr$objfun(fr$init)
      dec <- .newtonDecrement(ev)
      record("quadrature", length(phaseRows), k, ev$value, dec)
      if (dec < (quadrature$tolDecrement %||% 1e-4)) { conv <- TRUE; stopReason <- "decrement"; break }
      tc <- modifyList(list(rinit = 1, rmax = 10), quadrature$cm1Control %||% list())
      m <- .innerTrust(modifyList(tc, list(objfun = fr$objfun, parinit = fr$init,
                                           iterlim = quadrature$maxCm1Iter %||% 30L)))
      b <- fr$back(m$argument)
      st$theta <- b$theta; st$omegaPars[] <- b$omegaPars; st$hyp <- b$hyp
    }
  }

  est <- eStep(quadrature, final)
  theta <- st$theta; hyp <- st$hyp
  total <- sum(vapply(est, `[[`, 0, "m2logI")) + penaltyOf()
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
       stageTrace = do.call(rbind, st$rows), phases = do.call(rbind, phaseRows), burnin = burnSet,
       hyperprior = if (any(!vapply(tauPriors, is.null, TRUE))) penaltyOf(),
       tauPosterior = if (st$integratedNow)
         setNames(list(c(logMean = st$tauGrid[["mean"]], logSd = st$tauGrid[["sd"]])),
                  lay$tauName[match(integrated, lay$block)]),
       probPosterior = Filter(Negate(is.null), setNames(lapply(unique(lay$block), function(b)
         .spikeProbPosterior(lay, est, b)), lay$probName[match(unique(lay$block), lay$block)])))
}
