#' Structure selection by empirical Bayes
#'
#' @description
#' Selects which candidate terms of a model are present: gated parameters
#' (mechanism terms, from [gateL1]) and reference parameters (differences to a
#' reference, e.g. log fold changes of a cell line). Every candidate is either
#' exactly absent or drawn from a slab whose width is estimated with the
#' structure, so the penalty per term is estimated rather than set by hand.
#' One multistart of this single problem replaces a penalty grid with a
#' multistart per strength.
#'
#' @details
#' A structure \eqn{z \in \{0,1\}^J} switches the candidates on or off; an absent
#' gate or reference parameter is fixed at 0. The slab is a g-prior on the
#' effects of the present candidates: normal around 0 with \eqn{g} times the
#' covariance of their estimate, so \eqn{g} is the squared effect size in units
#' of its standard error and does not depend on parametrisation or units. With
#' the likelihood ratio in place of the Wald statistic, the structure scores
#' \deqn{G(z, g) = \frac{g\,D_z + D_0}{1 + g} + (p_z + a_g) \log(1 + g)
#'   - 2 \log p(z),}
#' where \eqn{D_z} is the minimal \eqn{-2 \log L} of structure \eqn{z}, \eqn{D_0}
#' that of the structure without candidates, \eqn{p_z} the number of parameters
#' of the present candidates, \eqn{(1 + g)^{-a_g/2}} the hyper-g prior and
#' \eqn{p(z)} the beta-binomial structure prior per family, i.e. the share of
#' present candidates integrated out under a \eqn{\mathrm{Beta}(a, b)} prior,
#' which corrects for the number of candidates. Minimising over \eqn{g} gives
#' \eqn{1 + \hat g = (D_0 - D_z) / (p_z + a_g)}.
#'
#' Each start draws a structure and parameters, then alternates a search over
#' single switches of candidates (each a fit warm-started from the current one)
#' with the update of \eqn{g}, until neither lowers \eqn{G}. The values of
#' \eqn{G} are comparable across starts.
#'
#' A candidate with \eqn{p_j} parameters is present iff its gain
#' \eqn{\Delta_j = D_{\mathrm{off}} - D_{\mathrm{on}}} exceeds
#' \deqn{\lambda_j = \frac{1 + g}{g} \left(p_j \log(1 + g) - 2 \log
#'   \frac{p(z_{\mathrm{on}})}{p(z_{\mathrm{off}})}\right),}
#' equivalently iff its conditional inclusion probability exceeds 0.5. With
#' `rule = "alpha"`, every present parameter costs `qchisq(1 - alpha, 1)`
#' instead.
#'
#' The final structure comes from refits along the candidates ordered by their
#' gain above the threshold, present ones first: the smallest structure that a
#' likelihood ratio test against the full model does not reject at level
#' `alpha`.
#'
#' Candidates whose parameters are close to collinear at the best structure
#' are reported in `nonidentifiable`, and structures whose score is within
#' `control$window` of the best one in `alternatives`.
#'
#' @param obj Objective function on the \eqn{-2 \log L} scale, called as
#'   `obj(pars, fixed = , deriv = )`.
#' @param center Named start vector of all parameters except the gates.
#' @param zero Named list, gate name to the parameters it multiplies, as
#'   `attr(gateL1(...), "gates")`.
#' @param reference Names of reference parameters; absent means 0.
#' @param class Optional named character vector, candidate to family; each
#'   family has its own structure prior. Default: one family of gates, one of
#'   reference parameters.
#' @param rule `"pip"`: structure of lowest score \eqn{G}. `"alpha"`: penalty
#'   `qchisq(1 - alpha, 1)` per present parameter.
#' @param alpha Level of the likelihood ratio tests and of `rule = "alpha"`.
#' @param lrt `"data"`: the tests use the attribute `"data"` of the objective
#'   value where present. `"total"`: the objective value.
#' @param prior List: `a`, `b` (Beta prior on the share of present candidates,
#'   default \eqn{a = 1} and \eqn{b} the number of candidates of the family,
#'   which favours sparse structures) and `g` (exponent \eqn{a_g} of the
#'   hyper-g prior, default 3).
#' @param fits Number of starts per batch.
#' @param cores Number of forks for the starts.
#' @param sd Standard deviation of the start perturbation of the parameters.
#' @param control List: `trust` (arguments of [trust], default `rinit = 0.1`,
#'   `rmax = 10`, `iterlim = 200`), `hits` and `tolHits` (further batches until
#'   `hits` starts end within `tolHits` of the best score, default 1 and 0.1),
#'   `maxFits` (default `10 * fits`), `refitFits` (random starts per refit,
#'   default 3), `fullFits` (starts of the full and the empty structure,
#'   default `fits`), `kappa` (eigenvalue bound of the correlation of the
#'   present parameters, default 0.01), `window` (default 2), `swaps`
#'   (evaluate the exchange of every present with every absent candidate at
#'   the best structure, default `TRUE`), `tol` (default 1e-3), `seed`.
#' @return Object of class `selectEB`, a list with
#'   \describe{
#'     \item{`terms`}{per candidate: family, `on` at the best structure, `gain`,
#'       `threshold`, `pip` (`NA` with `rule = "alpha"`), `selected` after the
#'       refits.}
#'     \item{`hyper`}{\eqn{\hat g}, the penalty \eqn{(1 + g)/g \log(1 + g)} per
#'       parameter, and per family the posterior mean share of present
#'       candidates.}
#'     \item{`structures`}{one row per start: final structure and score.}
#'     \item{`alternatives`}{evaluated structures within `window` of the best.}
#'     \item{`refits`}{refitted structures with value, test statistic, degrees
#'       of freedom and p-value.}
#'     \item{`selected`, `fit`}{key and parameters of the selected structure.}
#'     \item{`nonidentifiable`}{groups of near-collinear present candidates.}
#'     \item{`full`, `empty`, `best`}{fits of the full and the empty structure,
#'       and the best start.}
#'   }
#' @example inst/examples/selectEB.R
#' @seealso [scanL1], [gateL1]
#' @export
selectEB <- function(obj, center, zero = NULL, reference = NULL, class = NULL,
                     rule = c("pip", "alpha"), alpha = 0.05, lrt = c("data", "total"),
                     prior = list(), fits = 20, cores = 1, sd = 0.5, control = list()) {
  rule <- match.arg(rule)
  lrt  <- match.arg(lrt)
  if (.Platform$OS.type == "windows") cores <- 1L
  prior <- utils::modifyList(list(a = 1, g = 3), prior)
  ctl <- utils::modifyList(list(trust = list(rinit = 0.1, rmax = 10, iterlim = 200),
                                hits = 1, tolHits = 0.1, maxFits = 10 * fits,
                                refitFits = 3, fullFits = fits, kappa = 0.01,
                                window = 2, swaps = TRUE, tol = 1e-3, seed = NULL),
                           control)
  spec <- .ebSpec(zero, reference, class)
  center <- center[setdiff(names(center), names(zero))]
  miss <- setdiff(unlist(spec$pars), names(center))
  if (length(miss))
    stop("selectEB: center lacks ", paste(miss, collapse = ", "), call. = FALSE)
  if (!is.null(ctl$seed)) set.seed(ctl$seed)

  env <- .ebEnv(obj, spec, prior, ctl, rule, alpha)
  J <- length(spec$cand)
  full  <- .ebAnchor(env, rep(TRUE, J), center, sd, cores)
  empty <- .ebAnchor(env, rep(FALSE, J), center, sd, cores)
  env$D0 <- empty$D
  runs <- .ebMultistart(env, center, full$argument, fits, sd, cores)
  best <- runs$best
  if (ctl$swaps) .ebSwaps(env, best)
  terms <- .ebTerms(env, best)
  refit <- .ebRefits(env, best, terms, full, center, sd, cores, alpha, lrt)
  terms$selected <- terms$candidate %in% refit$on

  out <- list(terms = terms, hyper = .ebHyperTable(env, best), structures = runs$table,
              alternatives = .ebAlternatives(env, best), refits = refit$table,
              selected = refit$key, fit = refit$argument,
              nonidentifiable = .ebGroups(env, best), full = full, empty = empty,
              best = best[c("key", "z", "g", "G", "D", "argument", "trace", "hits",
                            "starts")],
              nfit = env$nfit, rule = rule, alpha = alpha, lrt = lrt, prior = prior,
              window = ctl$window, spec = spec)
  class(out) <- "selectEB"
  out
}


# Candidates: gates, then reference parameters, with their parameters and
# families.
.ebSpec <- function(zero, reference, class) {
  gates <- names(zero)
  cand  <- c(gates, reference)
  if (!length(cand))
    stop("selectEB: no candidates, give `zero` or `reference`.", call. = FALSE)
  if (anyDuplicated(cand))
    stop("selectEB: candidate names must be unique.", call. = FALSE)
  type <- stats::setNames(rep(c("gate", "reference"), c(length(gates), length(reference))),
                          cand)
  pars <- c(lapply(zero, as.character), stats::setNames(as.list(reference), reference))
  fam  <- stats::setNames(ifelse(type == "gate", "zero", "reference"), cand)
  if (!is.null(class)) {
    if (is.null(names(class)) || !all(names(class) %in% cand))
      stop("selectEB: `class` must be named by candidates.", call. = FALSE)
    fam[names(class)] <- as.character(class)
  }
  list(cand = cand, type = type, pars = pars, family = fam)
}

# Shared state of one selection: the best fit per structure over all starts,
# and the structures fitted in the current start.
.ebEnv <- function(obj, spec, prior, ctl, rule, alpha) {
  env <- new.env(parent = emptyenv())
  env$obj <- obj; env$spec <- spec; env$prior <- prior; env$ctl <- ctl
  env$rule <- rule; env$pen <- stats::qchisq(1 - alpha, 1)
  env$best <- list(); env$seen <- character(0); env$nfit <- 0L
  env
}

.ebKey <- function(z, spec) {
  off <- spec$cand[!z]
  if (!length(off)) "full" else paste0("-", paste(off, collapse = " -"))
}

.ebUnkey <- function(key, spec) {
  if (key == "full") return(rep(TRUE, length(spec$cand)))
  !spec$cand %in% sub("^-", "", strsplit(key, " ", fixed = TRUE)[[1]])
}

# Fixed values of structure `z`: gates at 0 or 1, the parameters of absent
# gates at `par`, absent reference parameters at 0.
.ebFixed <- function(z, spec, par) {
  g    <- spec$type == "gate"
  offG <- spec$cand[g & !z]
  offR <- spec$cand[!g & !z]
  c(stats::setNames(as.numeric(z[g]), spec$cand[g]),
    par[unique(unlist(spec$pars[offG]))],
    stats::setNames(rep(0, length(offR)), offR))
}

# Fit with the parameters in `fx` fixed; a structure without free parameters
# is evaluated.
.ebTrust <- function(env, start, fx) {
  if (!length(start)) {
    v <- env$obj(start, fixed = fx, deriv = FALSE)
    return(list(value = v$value, argument = start, hessian = matrix(0, 0, 0),
                converged = TRUE))
  }
  do.call(trust, c(list(env$obj, start), .l1TrustArgs(env$ctl$trust), list(fixed = fx)))
}

# Full parameter vector of a fit: `start` with the fitted and the fixed values.
.ebArg <- function(start, fit, fx) {
  start[names(fit$argument)] <- fit$argument
  both <- intersect(names(fx), names(start))
  start[both] <- fx[both]
  start
}

.ebFit <- function(env, z, start) {
  fx   <- .ebFixed(z, env$spec, start)
  free <- setdiff(names(start), names(fx))
  env$nfit <- env$nfit + 1L
  fit <- try(suppressWarnings(.ebTrust(env, start[free], fx)), silent = TRUE)
  if (inherits(fit, "try-error") || !is.finite(fit$value)) return(NULL)
  list(D = fit$value, argument = .ebArg(start, fit, fx), hessian = fit$hessian,
       converged = isTRUE(fit$converged))
}

# Structure `z`, fitted once per start from `start` and scored with the best
# fit over all starts.
.ebEval <- function(env, z, start) {
  key <- .ebKey(z, env$spec)
  if (!key %in% env$seen) {
    env$seen <- c(env$seen, key)
    f <- .ebFit(env, z, start)
    old <- env$best[[key]]
    if (!is.null(f) && (is.null(old) || f$D < old$D)) env$best[[key]] <- f
  }
  b <- env$best[[key]]
  list(z = z, key = key, D = if (is.null(b)) Inf else b$D,
       argument = if (is.null(b)) start else b$argument,
       p = sum(lengths(env$spec$pars[z])))
}

# -2 log beta-binomial probability of the structure, per family.
.ebStructPrior <- function(z, spec, prior) {
  sum(vapply(split(z, spec$family[spec$cand]), function(zf) {
    b <- prior$b %||% length(zf)
    -2 * (lbeta(sum(zf) + prior$a, sum(!zf) + b) - lbeta(prior$a, b))
  }, 0))
}

.ebScore <- function(env, s, g) {
  if (env$rule == "alpha") return(s$D + env$pen * s$p)
  (g * s$D + env$D0) / (1 + g) + (s$p + env$prior$g) * log1p(g) +
    .ebStructPrior(s$z, env$spec, env$prior)
}

.ebGhat <- function(env, s) {
  if (env$rule == "alpha") return(NA_real_)
  max(0, (env$D0 - s$D) / (s$p + env$prior$g) - 1)
}

# One start: switches of single candidates in random order at fixed g until
# none lowers G, then the update of g; both until neither improves.
.ebRun <- function(env, z, start) {
  tol <- env$ctl$tol
  env$seen <- character(0)
  cur <- .ebEval(env, z, start)
  g <- .ebGhat(env, cur)
  cur$G <- .ebScore(env, cur, g)
  trace <- cur$G
  repeat {
    repeat {
      moved <- FALSE
      for (j in sample.int(length(z))) {
        z1 <- cur$z; z1[j] <- !z1[j]
        nb <- .ebEval(env, z1, cur$argument)
        nb$G <- .ebScore(env, nb, g)
        if (nb$G < cur$G - tol) { cur <- nb; moved <- TRUE; trace <- c(trace, cur$G) }
      }
      if (!moved) break
    }
    if (env$rule == "alpha") break
    g1 <- .ebGhat(env, cur)
    G1 <- .ebScore(env, cur, g1)
    if (!(G1 < cur$G - tol)) break
    g <- g1; cur$G <- G1; trace <- c(trace, G1)
  }
  cur$g <- g
  cur$trace <- trace
  cur
}

# Multistart fit of the full or the empty structure, entered into the
# structure table.
.ebAnchor <- function(env, z, center, sd, cores) {
  ctl <- env$ctl
  fx  <- .ebFixed(z, env$spec, center)
  free <- setdiff(names(center), names(fx))
  r <- .l1Multistart(function(st) .ebTrust(env, st, fx), center[free], ctl$fullFits,
                     sd, cores, env$obj,
                     wf = list(hits = ctl$hits, tol = ctl$tolHits, max = ctl$maxFits))
  if (is.null(r))
    stop("selectEB: no start of the ", if (all(z)) "full" else "empty",
         " structure converged.", call. = FALSE)
  env$nfit <- env$nfit + r$starts
  arg <- .ebArg(center, r, fx)
  key <- .ebKey(z, env$spec)
  env$best[[key]] <- list(D = r$value, argument = arg, hessian = r$hessian,
                          converged = isTRUE(r$converged))
  list(D = r$value, argument = arg, values = r$values, starts = r$starts, hits = r$hits)
}

# Starts: all candidates on from the full optimum, all off, then random
# structures and parameters; batches until `hits` starts end within `tolHits`.
.ebMultistart <- function(env, center, fullArg, fits, sd, cores) {
  spec <- env$spec
  ctl  <- env$ctl
  J    <- length(spec$cand)
  draw <- function(i) {
    if (i == 1L) return(list(z = rep(TRUE, J), par = fullArg))
    if (i == 2L) return(list(z = rep(FALSE, J), par = center))
    list(z = stats::runif(J) < 0.5, par = center + stats::rnorm(length(center), 0, sd))
  }
  one <- function(s) {
    try(resetWarmStarts(env$obj, verbose = FALSE), silent = TRUE)
    n0 <- env$nfit
    r <- try(.ebRun(env, s$z, s$par), silent = TRUE)
    if (inherits(r, "try-error")) NULL else list(run = r, best = env$best, nfit = env$nfit - n0)
  }
  runs <- list(); n <- 0L
  repeat {
    batch <- lapply(n + seq_len(max(fits, 1L)), draw)
    out <- if (cores > 1L) parallel::mclapply(batch, one, mc.cores = cores,
                                              mc.preschedule = FALSE)
           else lapply(batch, one)
    n <- n + length(batch)
    for (o in Filter(is.list, out)) {
      runs <- c(runs, list(o$run))
      if (cores > 1L) {
        .ebMerge(env, o$best)
        env$nfit <- env$nfit + o$nfit
      }
    }
    if (!length(runs)) break
    G <- vapply(runs, `[[`, 0, "G")
    hits <- sum(G <= min(G) + ctl$tolHits)
    if (hits >= ctl$hits || n >= ctl$maxFits) break
  }
  if (!length(runs)) stop("selectEB: no start finished.", call. = FALSE)
  G <- vapply(runs, `[[`, 0, "G")
  tab <- data.frame(start = seq_along(runs), key = vapply(runs, `[[`, "", "key"),
                    G = G, g = vapply(runs, `[[`, 0, "g"),
                    size = vapply(runs, function(r) sum(r$z), 0L))
  tab <- tab[order(tab$G), ]
  rownames(tab) <- NULL
  best <- runs[[which.min(G)]]
  best$hits <- hits
  best$starts <- n
  list(best = best, table = tab)
}

.ebMerge <- function(env, best) {
  for (k in names(best)) {
    old <- env$best[[k]]
    if (is.null(old) || best[[k]]$D < old$D) env$best[[k]] <- best[[k]]
  }
}

# Exchange of every present with every absent candidate at the best structure.
.ebSwaps <- function(env, best) {
  env$seen <- character(0)
  for (i in which(best$z)) for (j in which(!best$z)) {
    z1 <- best$z; z1[i] <- FALSE; z1[j] <- TRUE
    .ebEval(env, z1, best$argument)
  }
}

.ebState <- function(env, z, g) {
  key <- .ebKey(z, env$spec)
  b <- env$best[[key]]
  s <- list(z = z, key = key, D = if (is.null(b)) Inf else b$D,
            p = sum(lengths(env$spec$pars[z])))
  s$G <- .ebScore(env, s, g)
  s
}

# Per candidate at the best structure: gain in -2 log L of switching it on,
# the threshold the gain has to exceed, and the conditional inclusion
# probability.
.ebTerms <- function(env, best) {
  spec <- env$spec
  w <- if (env$rule == "alpha") 1 else best$g / (1 + best$g)
  rows <- lapply(seq_along(spec$cand), function(j) {
    z1 <- best$z; z1[j] <- !z1[j]
    nb  <- .ebState(env, z1, best$g)
    inn <- if (best$z[j]) best else nb
    out <- if (best$z[j]) nb else best
    gain <- out$D - inn$D
    data.frame(candidate = spec$cand[j], family = spec$family[[j]], on = best$z[j],
               gain = gain, threshold = gain - (out$G - inn$G) / w,
               pip = if (env$rule == "pip") stats::plogis((out$G - inn$G) / 2) else NA_real_)
  })
  do.call(rbind, rows)
}

.ebHyperTable <- function(env, best) {
  spec <- env$spec
  if (env$rule == "alpha")
    return(data.frame(g = NA_real_, penalty = env$pen))
  g <- best$g
  fam <- do.call(rbind, lapply(unique(spec$family), function(f) {
    zf <- best$z[spec$family == f]
    b <- env$prior$b %||% length(zf)
    data.frame(family = f, candidates = length(zf), on = sum(zf),
               share = (sum(zf) + env$prior$a) / (length(zf) + env$prior$a + b))
  }))
  cbind(g = g, penalty = if (g > 0) (1 + g) / g * log1p(g) else Inf, fam)
}

# Evaluated structures whose score at the estimated g is within `window` of
# the best.
.ebAlternatives <- function(env, best) {
  rows <- lapply(names(env$best), function(k) {
    s <- .ebState(env, .ebUnkey(k, env$spec), best$g)
    data.frame(key = k, size = sum(s$z), D = s$D, G = s$G)
  })
  tab <- do.call(rbind, rows)
  tab <- tab[tab$G <= best$G + env$ctl$window, ]
  tab <- tab[order(tab$G), ]
  rownames(tab) <- NULL
  tab
}

# Groups of present candidates whose parameters are close to collinear at
# the best structure: eigenvectors of the correlation form of the profile
# Hessian with eigenvalue below `kappa`, and parameters without curvature.
.ebGroups <- function(env, best) {
  H <- env$best[[best$key]]$hessian
  I <- intersect(unlist(env$spec$pars[best$z]), rownames(H))
  if (is.null(H) || !length(I)) return(list())
  owner <- stats::setNames(rep(names(env$spec$pars), lengths(env$spec$pars)),
                           unlist(env$spec$pars))
  P <- .ebSchur(H / 2, I, setdiff(rownames(H), I))
  d <- diag(P)
  flat <- I[d <= 1e-10 * max(d, 1e-300)]
  out <- lapply(flat, function(p)
    list(candidates = unname(owner[p]), parameters = p, loadings = stats::setNames(1, p),
         eigenvalue = 0))
  J <- setdiff(I, flat)
  if (length(J) < 2L) return(out)
  s <- 1 / sqrt(d[J])
  e <- eigen(P[J, J, drop = FALSE] * outer(s, s), symmetric = TRUE)
  for (k in which(e$values < env$ctl$kappa)) {
    ld <- e$vectors[, k]
    sel <- abs(ld) > 0.3
    out <- c(out, list(list(candidates = unique(unname(owner[J[sel]])), parameters = J[sel],
                            loadings = stats::setNames(ld[sel], J[sel]),
                            eigenvalue = e$values[k])))
  }
  out
}

# Profile Hessian of `I`: the Schur complement of `A` over `N`, with the
# directions of `N` without curvature left out.
.ebSchur <- function(A, I, N) {
  if (!length(N)) return(A[I, I, drop = FALSE])
  e <- eigen(A[N, N, drop = FALSE], symmetric = TRUE)
  keep <- e$values > 1e-10 * max(abs(e$values), 1e-300)
  V <- e$vectors[, keep, drop = FALSE]
  Ninv <- V %*% (t(V) / e$values[keep])
  P <- A[I, I, drop = FALSE] - A[I, N, drop = FALSE] %*% Ninv %*% A[N, I, drop = FALSE]
  (P + t(P)) / 2
}

# Refits along the candidates ordered by gain above threshold, from the size
# of the best structure down while the test against the full model does not
# reject, up while it does.
.ebRefits <- function(env, best, terms, full, center, sd, cores, alpha, lrt) {
  spec <- env$spec
  ctl  <- env$ctl
  cand <- terms$candidate[order(!terms$on, terms$threshold - terms$gain)]
  J <- length(cand)
  Dfull <- .ebDataValue(env, full$argument, rep(TRUE, J), lrt)
  nFull <- sum(lengths(spec$pars))
  cache <- list()
  refit <- function(k) {
    kk <- as.character(k)
    if (!is.null(cache[[kk]])) return(cache[[kk]])
    z <- spec$cand %in% cand[seq_len(k)]
    key <- .ebKey(z, spec)
    fx  <- .ebFixed(z, spec, full$argument)
    free <- setdiff(names(center), names(fx))
    warm <- env$best[[key]]$argument %||% full$argument
    r <- .l1Multistart(function(st) .ebTrust(env, st, fx), full$argument[free],
                       ctl$refitFits, sd, cores, env$obj, extra = list(warm[free]))
    env$nfit <- env$nfit + r$starts
    arg <- .ebArg(full$argument, r, fx)
    D  <- .ebDataValue(env, arg, z, lrt)
    df <- nFull - sum(lengths(spec$pars[z]))
    stat <- max(0, D - Dfull)
    res <- data.frame(size = k, key = key, value = r$value, data = D, stat = stat, df = df,
                      p = if (df > 0) stats::pchisq(stat, df, lower.tail = FALSE) else 1,
                      starts = r$starts, hits = r$hits)
    res$argument <- list(arg)
    cache[[kk]] <<- res
    res
  }
  k <- sum(best$z)
  if (refit(k)$p >= alpha) {
    while (k > 0L && refit(k - 1L)$p >= alpha) k <- k - 1L
  } else {
    while (k < J && refit(k)$p < alpha) k <- k + 1L
  }
  tab <- do.call(rbind, cache)
  sel <- refit(k)
  tab <- tab[order(tab$size), setdiff(names(tab), "argument")]
  rownames(tab) <- NULL
  list(table = tab, key = sel$key, argument = sel$argument[[1]], on = cand[seq_len(k)])
}

.ebDataValue <- function(env, arg, z, lrt) {
  fx <- .ebFixed(z, env$spec, arg)
  free <- setdiff(names(arg), names(fx))
  v <- env$obj(arg[free], fixed = fx, deriv = FALSE)
  d <- attr(v, "data", exact = TRUE)
  if (lrt == "data" && is.numeric(d)) sum(d) else v$value
}


#' @export
print.selectEB <- function(x, ...) {
  cat(sprintf("selectEB: rule %s, %d candidates, %d starts (%d hits), %d fits\n",
              x$rule, nrow(x$terms), x$best$starts, x$best$hits, x$nfit))
  print(x$hyper, row.names = FALSE, digits = 3)
  cat("\n")
  print(x$terms, row.names = FALSE, digits = 3)
  cat(sprintf("\nbest structure: %s\nselected after refits (alpha = %g): %s\n",
              x$best$key, x$alpha, x$selected))
  if (nrow(x$alternatives) > 1L)
    cat(sprintf("%d further structures within %g of the best score\n",
                nrow(x$alternatives) - 1L, x$window))
  for (g in x$nonidentifiable)
    cat(sprintf("near-collinear (eigenvalue %.2g): %s\n", g$eigenvalue,
                paste(g$candidates, collapse = ", ")))
  invisible(x)
}

utils::globalVariables(c("candidate", "gain", "threshold", "decision", "start", "G", "step"))

#' @param x Object of class `selectEB`.
#' @param type `"terms"`: gain and threshold per candidate. `"waterfall"`: final
#'   score of every start above the best. `"trace"`: score along the best start.
#' @param ... Not used.
#' @return A ggplot.
#' @rdname selectEB
#' @export
plot.selectEB <- function(x, type = c("terms", "waterfall", "trace"), ...) {
  type <- match.arg(type)
  if (type == "waterfall") {
    s <- x$structures
    s$start <- seq_len(nrow(s))
    return(ggplot2::ggplot(s, ggplot2::aes(start, G - min(G), colour = key)) +
             ggplot2::geom_point() +
             ggplot2::scale_y_continuous(trans = "log1p", breaks = c(0, 10^(0:6))) +
             ggplot2::labs(x = "start, sorted", y = "score above the best", colour = NULL) +
             theme_dMod() + ggplot2::theme(legend.position = "bottom", legend.direction = "vertical"))
  }
  if (type == "trace") {
    tr <- x$best$trace
    return(ggplot2::ggplot(data.frame(step = seq_along(tr) - 1L, G = tr), ggplot2::aes(step, G)) +
             ggplot2::geom_step() + ggplot2::geom_point(size = 0.8) +
             ggplot2::labs(x = "accepted move", y = "score G") + theme_dMod())
  }
  tm <- x$terms
  tm$candidate <- factor(tm$candidate, levels = rev(tm$candidate))
  tm$decision <- ifelse(tm$selected, "selected", "not selected")
  ggplot2::ggplot(tm, ggplot2::aes(y = candidate)) +
    ggplot2::geom_point(ggplot2::aes(x = pmax(gain, 0), colour = decision)) +
    ggplot2::geom_point(ggplot2::aes(x = pmin(pmax(threshold, 0), 1e6)), shape = 124, size = 4) +
    ggplot2::scale_x_continuous(trans = "log1p", breaks = c(0, 10^(0:6))) +
    ggplot2::labs(x = "gain in -2 log L, threshold as bar", y = NULL, colour = NULL) +
    theme_dMod()
}
