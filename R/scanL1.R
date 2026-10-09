utils::globalVariables(c("lambda", "parameter", "member", "term"))

#' Gate parameters for an L1 selection toward zero
#'
#' @description Multiplies each named entry of a parameter transformation by a
#' gate parameter `<prefix><par>`, so that `scanL1()` can pull the parameter to
#' exactly zero. In a log parametrisation such as `k = exp10(log10_k)` zero lies
#' at minus infinity and an L1 term on `log10_k` pulls `k` toward one, not zero.
#' The gate `s_k` in `k = exp10(log10_k) * s_k` is linear and dimensionless: at
#' `s_k = 0` the parameter is gone, and one strength applies to all gates alike.
#' The full model fixes every gate at one.
#'
#' @param trafo An [eqnvec], or a list of them (one per condition, as from
#'   [branch]).
#' @param pars Character, names of `trafo` to gate.
#' @param prefix Prefix of the gate names.
#' @return `trafo` with the gates inserted and attribute `"gates"`: a named list,
#'   gate name to the symbols of the original expression. `scanL1()` takes it as
#'   `zero` and fixes these symbols while the gate is free.
#' @seealso [scanL1]
#' @example inst/examples/scanL1.R
#' @export
gateL1 <- function(trafo, pars, prefix = "s_") {
  one <- is.eqnvec(trafo)
  tl  <- if (one) list(trafo) else trafo
  miss <- setdiff(pars, names(tl[[1]]))
  if (length(miss))
    stop("gateL1: not in the transformation: ", paste(miss, collapse = ", "),
         call. = FALSE)
  gates <- stats::setNames(lapply(pars, function(p) getSymbols(tl[[1]][[p]])),
                           paste0(prefix, pars))
  tl <- lapply(tl, function(tr) {
    for (p in pars) tr[[p]] <- paste0("(", tr[[p]], ") * ", prefix, p)
    tr
  })
  out <- if (one) tl[[1]] else tl
  attr(out, "gates") <- gates
  out
}


#' L1 selection over a grid of penalty strengths
#'
#' @description Selects a sparse model structure by L1 regularisation in the
#' manner of Hauber, Rosenblatt and Timmer (2023). For every `lambda` of the
#' grid the penalised objective is minimised with [trustL1] from several
#' starts. The structure the best fit lands on (parameters at zero, fused
#' groups) is then refitted without penalty, and the refits are compared with
#' the full model.
#'
#' Three penalties can be combined:
#' \describe{
#'   \item{`zero`}{gates from [gateL1], penalty `lambda * s`; a gate at zero
#'     removes its parameter.}
#'   \item{`reference`}{parameters pulled toward zero, penalty
#'     `lambda * |r|`, typically log fold changes to a reference cell type.}
#'   \item{`groups`}{blocks with the pairwise penalty
#'     `lambda * sum_{i<j} |r_i - r_j|`; members with equal values form one
#'     group.}
#' }
#'
#' @param obj Objective without penalty, e.g. `normL2(data, g * x * p)`, on the
#'   `-2 log L` scale.
#' @param center Named numeric, start of the full model. Gate names are
#'   ignored.
#' @param zero Named list, gate name to the parameters it multiplies, as in
#'   `attr(gateL1(...), "gates")`. In the selection step these parameters are
#'   fixed at the full fit and the gate is free; in the refit a surviving gate
#'   is fixed at one with its parameters free, a removed gate at zero with its
#'   parameters fixed.
#' @param reference Character, parameters pulled toward zero.
#' @param groups List of blocks. A block is a character vector, or a list with
#'   `pars` and `anchor`: the anchor is a fixed value every member is also
#'   pulled toward, e.g. 0 for the reference cell type.
#' @param lambda Grid of penalty strengths.
#' @param fixed Named numeric, parameters fixed throughout.
#' @param q Exponent of the penalty. `q < 1` (e.g. 0.8) is approximated by
#'   reweighted L1 fits.
#' @param ssl `NULL`, or a list with `lambda1` (slab strength, default 1), `a`
#'   and `b` (Beta prior of the slab share, default 2 and 2, which keeps the
#'   share away from zero). Switches to the spike-and-slab lasso: every
#'   penalised term is drawn from the slab with strength `lambda1` or from the
#'   spike with strength `lambda`, and `lambda` becomes the spike strength of
#'   the path. See Details.
#' @param select `"lrt"`: the largest `lambda` whose refit the likelihood ratio
#'   test against the full model does not reject at level `alpha`. `"bic"`: the
#'   refit with the smallest BIC. `"plateau"`: the structure at the largest
#'   `lambda`, meant for `ssl`, where the structure stops changing once the
#'   spike is strong enough.
#' @param alpha Level of the likelihood ratio test.
#' @param fits Number of starts of the full model.
#' @param pathFits Number of starts per `lambda` and per refit. For a `lambda`
#'   they include the optimum of the next smaller `lambda` and the sparse point
#'   (every gate and reference at zero, every block at its anchor or mean); a
#'   downward pass then restarts each `lambda` from the optimum of the next
#'   larger one and keeps the better fit. A refit starts from the full optimum
#'   and from the penalised optimum of its structure. The remaining starts are
#'   drawn around the full optimum.
#' @param cores Number of forked processes for the starts.
#' @param sd Standard deviation of the start perturbation.
#' @param control List: `trust` (arguments of [trustL1] and [trust], default
#'   `rinit = 0.1, rmax = 10, iterlim = 200`), `nq` (reweighting rounds for
#'   `q < 1`, default 3), `eps` (floor of `|difference|` in the reweighting,
#'   default 0.01), `ndata` (data points for BIC, default from `obj`), `nem`
#'   (iterations of the spike-and-slab EM, default 50), `tolp` (change of the
#'   inclusion probabilities that ends the EM, default 1e-4), `nmerge` (merge
#'   moves per `lambda`, default 5, see Details), `snap` (spike terms closer
#'   to their kink than this are put onto it without a refit, default 1e-6),
#'   `hits` (runs that must reach
#'   the best value within `tolHits`, default 1 and 0.1, for the full model, a
#'   `lambda` and a refit; further batches of starts are added until then, up
#'   to `maxFits`, default ten times the batch).
#' @details With `ssl`, term `j` (a gate, a reference parameter, or in a block
#'   the gap between neighbours of the sorted values, anchor included, as in
#'   Ke, Fan and Wu 2015) has the prior
#'   \deqn{\pi(d_j \mid \theta) = \theta\,\psi_1(d_j) + (1-\theta)\,\psi_0(d_j),
#'   \quad \psi_i(d) = \tfrac{\lambda_i}{4} e^{-\lambda_i |d| / 2},}
#'   so that `-2 log` of each component is `lambda_i * |d|` up to a constant,
#'   the scale of `obj`. The posterior mode at one spike strength `lambda` is
#'   found by EM: the E-step gives the inclusion probability
#'   \eqn{p_j = \theta\psi_1(d_j) / \pi(d_j \mid \theta)} and
#'   \eqn{\theta = (\sum_j p_j + a - 1)/(J + a + b - 2)}, the M-step is an
#'   L1 fit by [trustL1] with weight \eqn{p_j \lambda_1 + (1 - p_j) \lambda}
#'   on term `j`, the neighbours taken at the current values. Terms with \eqn{p_j < 1/2} that the weighted fit leaves
#'   off their kink are then put onto it and the EM is rerun; the merge is kept
#'   if it lowers the `-2 log` posterior. The path is one chain over the
#'   increasing grid, each
#'   `lambda` started from the optimum of the previous one (Rockova and George
#'   2018) and from the sparse point, then one chain back down; each `lambda`
#'   keeps the mode of smaller `-2 log` posterior. `pathFits` then only
#'   applies to the refits.
#' @return Object of class `scanL1` with
#'   \describe{
#'     \item{`path`}{one row per `lambda`: penalised value, structure key,
#'       number of removed parameters, p-value of its refit, convergence of the
#'       penalised fit and, with `ssl`, the number of EM iterations.}
#'     \item{`coefficients`}{penalised estimates per `lambda`.}
#'     \item{`arguments`}{all free parameters of the penalised fit per `lambda`.}
#'     \item{`level`}{one row per `lambda` and structure on the lowest level of
#'       its waterfall, the runs within `tolHits` of the best value: structure
#'       key, best value and number of runs.}
#'     \item{`levelFits`}{per `lambda`, these runs: free parameters and value.}
#'     \item{`refits`}{one row per distinct structure: `-2 log L`, free
#'       parameters, LRT statistic, degrees of freedom, p-value, BIC.}
#'     \item{`full`}{the full fit: value, parameters and the sorted values of
#'       all its starts.}
#'     \item{`selected`, `lambdaSelected`}{the chosen structure and the
#'       `lambda` it was chosen at.}
#'     \item{`structure`}{per key: removed parameters and groups.}
#'     \item{`fit`}{named parameters of the chosen refit, fixed ones included.}
#'     \item{`inclusion`, `theta`}{with `ssl`: inclusion probability of every
#'       gate, reference parameter and pairwise difference at the final slab
#'       share, and that share, per `lambda`.}
#'   }
#' @references Hauber AL, Rosenblatt M, Timmer J (2023). Uncovering specific
#'   mechanisms across cell types in dynamical models. PLoS Comput Biol 19(9):
#'   e1010867.
#'
#'   Ke ZT, Fan J, Wu Y (2015). Homogeneity pursuit. J Am Stat Assoc
#'   110(509):175-194.
#'
#'   Rockova V, George EI (2018). The spike-and-slab LASSO. J Am Stat Assoc
#'   113(521): 431-444.
#' @seealso [gateL1], [trustL1]
#' @example inst/examples/scanL1.R
#' @export
scanL1 <- function(obj, center, zero = NULL, reference = NULL, groups = NULL,
                   lambda = 10^seq(-3, 3, length.out = 25), fixed = NULL,
                   q = 1, ssl = NULL, select = c("lrt", "bic", "plateau"),
                   alpha = 0.05, fits = 10, pathFits = fits, cores = 1, sd = 0.5,
                   control = list()) {
  select <- match.arg(select)
  if (!length(zero) && !length(reference) && !length(groups))
    stop("scanL1: give at least one of zero, reference, groups.", call. = FALSE)
  if (q <= 0 || q > 1) stop("scanL1: q must lie in (0, 1].", call. = FALSE)
  if (!is.null(ssl) && q < 1)
    stop("scanL1: ssl and q < 1 are two different penalties, give one.", call. = FALSE)
  if (identical(lambda, "em") && (length(groups) || !is.null(ssl)))
    stop("scanL1: lambda = \"em\" takes gates and reference parameters, no groups or ssl.",
         call. = FALSE)
  ctl <- utils::modifyList(list(trust = list(rinit = 0.1, rmax = 10, iterlim = 200L),
                                nq = 3L, eps = 0.01, ndata = NULL, nem = 50L,
                                tolp = 1e-4, nmerge = 5L, snap = 1e-6, hits = 1L,
                                tolHits = 0.1, maxFits = NULL,
                                em = list(init = 1, a = 1, b = 0, tol = 1e-3, adaptive = TRUE)),
                           control)
  wf <- if (ctl$hits > 1L) list(hits = ctl$hits, tol = ctl$tolHits, max = ctl$maxFits)
  if (.Platform$OS.type == "windows") cores <- 1L

  gates  <- names(zero)
  groups <- lapply(groups, function(b) if (is.character(b)) list(pars = b) else b)
  center <- center[setdiff(names(center), gates)]
  pen    <- c(reference, unlist(lapply(groups, `[[`, "pars")))
  miss   <- setdiff(c(pen, unlist(zero)), names(center))
  if (length(miss))
    stop("scanL1: not among the free parameters: ", paste(miss, collapse = ", "),
         call. = FALSE)

  ## 1. Full model ------------------------------------------------------------
  fixFull <- c(fixed, stats::setNames(rep(1, length(gates)), gates))
  full <- .l1Multistart(function(st)
    do.call(trust, c(list(obj, st), .l1TrustArgs(ctl$trust), list(fixed = fixFull))),
                        center, fits, sd, cores, obj, wf = wf)
  if (is.null(full)) stop("scanL1: every fit of the full model failed.", call. = FALSE)
  fullPars <- full$argument
  nFull    <- length(fullPars)

  ## 2. Penalised path --------------------------------------------------------
  gatedPars <- unique(unlist(zero))
  fixSel    <- c(fixed, fullPars[gatedPars])
  start0    <- c(fullPars[setdiff(names(fullPars), gatedPars)],
                 stats::setNames(rep(1, length(gates)), gates))
  sparse    <- start0
  sparse[c(gates, reference)] <- 0
  for (b in groups) sparse[b$pars] <- b$anchor %||% mean(start0[b$pars])
  if (identical(lambda, "em"))
    return(.l1EmScan(obj, start0, sparse, gates, reference, fixSel, q, ctl, pathFits,
                     sd, cores, wf, full, fullPars, zero, fixed, alpha))
  lambda    <- sort(lambda)
  if (!is.null(ssl))
    ssl <- utils::modifyList(list(lambda1 = 1, a = 2, b = 2), ssl)
  pen1 <- function(start, l, n, extra, prior = NULL)
    .l1Penalised(obj, start, lambda[l], gates, reference, groups, fixSel, q,
                 ctl, n, sd, cores, extra, ssl, if (n > 1L) wf, prior)
  # Upward from the full optimum, then downward from each larger lambda's
  # optimum; a lambda keeps the better of the two. With ssl both are chains of
  # EM runs, upward from the previous mode and from the sparse point.
  fitsL <- vector("list", length(lambda))
  warm  <- start0
  for (l in seq_along(lambda)) {
    fitsL[[l]] <- if (is.null(ssl)) pen1(start0, l, pathFits, list(warm, sparse))
                  else pen1(warm, l, 2L, list(sparse))
    if (!is.null(fitsL[[l]])) warm <- fitsL[[l]]$argument
  }
  for (l in rev(seq_along(lambda))[-1]) {
    up <- fitsL[[l + 1L]]
    if (is.null(up)) next
    down <- pen1(up$argument, l, 1L, NULL)
    cur  <- fitsL[[l]]
    if (is.null(down) || (!is.null(cur) && down$value >= cur$value)) next
    # A better optimum from above: more starts until the waterfall reaches it.
    new <- if (!is.null(ssl) || is.null(wf) || is.null(cur)) down
      else pen1(start0, l, pathFits, NULL,
                list(fits = list(down), values = cur$values, starts = cur$starts + 1L))
    if (is.null(new)) next
    new$level  <- .l1Level(c(new$level, cur$level), ctl$tolHits)
    fitsL[[l]] <- new
  }
  path <- lapply(seq_along(lambda), function(l) {
    f <- fitsL[[l]]
    if (is.null(f)) return(NULL)
    c(f, list(lambda = lambda[l],
              structure = .l1Structure(f$argument, gates, reference, groups)))
  })
  path <- Filter(Negate(is.null), path)
  if (!length(path)) stop("scanL1: every penalised fit failed.", call. = FALSE)

  ## 3. Refit each distinct structure without penalty ------------------------
  keys    <- vapply(path, function(z) z$structure$key, "")
  structs <- lapply(path, `[[`, "structure")[!duplicated(keys)]
  names(structs) <- unique(keys)
  ndata <- ctl$ndata %||% .l1Ndata(obj)
  pathArg <- lapply(path, `[[`, "argument")[!duplicated(keys)]
  rf <- Map(.l1Refit, structs, pathArg,
            MoreArgs = list(obj = obj, fullPars = fullPars, zero = zero, fixed = fixed,
                            ctl = ctl, fits = pathFits, sd = sd, cores = cores, wf = wf))
  rf <- Filter(Negate(is.null), rf)
  if (!length(rf)) stop("scanL1: every refit failed.", call. = FALSE)
  refitTab <- data.frame(key = names(rf),
                         value = vapply(rf, `[[`, 0, "value"),
                         nfree = vapply(rf, `[[`, 0L, "nfree"),
                         starts = vapply(rf, function(z) z$starts %||% NA_integer_, 0L),
                         hits = vapply(rf, function(z) z$hits %||% NA_integer_, 0L),
                         stringsAsFactors = FALSE)
  refitTab$stat <- pmax(0, refitTab$value - full$value)
  refitTab$df   <- nFull - refitTab$nfree
  refitTab$p    <- ifelse(refitTab$df > 0,
                          stats::pchisq(refitTab$stat, refitTab$df, lower.tail = FALSE), 1)
  refitTab$bic  <- if (is.null(ndata)) NA_real_ else refitTab$value + log(ndata) * refitTab$nfree
  rownames(refitTab) <- NULL

  ## 4. Selection -------------------------------------------------------------
  pathTab <- data.frame(
    lambda  = vapply(path, `[[`, 0, "lambda"),
    value   = vapply(path, `[[`, 0, "value"),
    key     = keys,
    removed = vapply(path, function(z) length(z$structure$removed), 0L),
    stringsAsFactors = FALSE)
  pathTab$p <- refitTab$p[match(pathTab$key, refitTab$key)]
  pathTab$converged <- vapply(path, function(z) isTRUE(z$converged), TRUE)
  pathTab$starts <- vapply(path, function(z) z$starts %||% NA_integer_, 0L)
  pathTab$hits   <- vapply(path, function(z) z$hits %||% NA_integer_, 0L)
  if (!is.null(ssl)) pathTab$em <- vapply(path, function(z) z$emIterations %||% NA_integer_, 0L)
  sel <- .l1Select(pathTab, refitTab, select, alpha)

  coefs <- do.call(rbind, lapply(path, function(z) z$argument[c(gates, pen)]))
  rownames(coefs) <- NULL

  if (!is.null(ssl)) {
    incl <- do.call(rbind, lapply(path, `[[`, "inclusion"))
    out_ssl <- list(inclusion = cbind(lambda = pathTab$lambda, incl),
                    theta = vapply(path, `[[`, 0, "theta"), ssl = ssl)
  } else out_ssl <- NULL

  levelTab <- do.call(rbind, lapply(path, function(z) {
    k <- vapply(z$level, function(f) .l1Structure(f$argument, gates, reference, groups)$key, "")
    v <- vapply(z$level, `[[`, 0, "value")
    data.frame(lambda = z$lambda, key = unique(k),
               value = vapply(unique(k), function(u) min(v[k == u]), 0, USE.NAMES = FALSE),
               fits = vapply(unique(k), function(u) sum(k == u), 0L, USE.NAMES = FALSE),
               stringsAsFactors = FALSE)
  }))

  out <- list(path = pathTab, coefficients = cbind(lambda = pathTab$lambda, coefs),
              arguments = lapply(path, `[[`, "argument"),
              refits = refitTab,
              level = levelTab, levelFits = lapply(path, `[[`, "level"),
              full = list(value = full$value, argument = fullPars, values = full$values,
                          starts = full$starts, hits = full$hits),
              selected = sel$key, lambdaSelected = sel$lambda,
              structure = structs, select = select, alpha = alpha, q = q,
              fit = rf[[sel$key]]$argument, groups = groups, gates = gates,
              reference = reference)
  out <- c(out, out_ssl)
  class(out) <- "scanL1"
  out
}


# `control$trust` as `trust()` takes it: the flat tolerances of `trustL1()`
# go into `tolControl`.
.l1TrustArgs <- function(a) {
  tol <- intersect(names(a), c("ftol", "mtol", "gtol", "xtol", "rmin"))
  if (!length(tol)) return(a)
  a$tolControl <- utils::modifyList(a$tolControl %||% list(), a[tol])
  a[setdiff(names(a), tol)]
}

# Best of `fits` runs of `run(start)` from `center` and around it, with the runs
# within `levelTol` of it as `level`. With `wf`, batches follow until `wf$hits` runs
# lie within `wf$tol` or `wf$max` starts are spent; `prior` adds earlier runs.
.l1Multistart <- function(run, center, fits, sd, cores, obj, extra = NULL,
                          positive = character(0), wf = NULL, prior = NULL,
                          levelTol = wf$tol %||% 0.1) {
  draw <- function(n) lapply(seq_len(max(n, 0L)), function(i) {
    st <- center + stats::rnorm(length(center), 0, sd)
    st[positive] <- abs(st[positive])
    st
  })
  one <- function(st) {
    try(resetWarmStarts(obj, verbose = FALSE), silent = TRUE)
    # A start the solver cannot follow is dropped; its warnings say nothing more.
    f <- try(suppressWarnings(run(st)), silent = TRUE)
    if (inherits(f, "try-error") || !is.finite(f$value)) NULL else f
  }
  maxN   <- if (is.null(wf)) 0L else wf$max %||% (10L * max(fits, 1L))
  starts <- c(list(center), extra)
  starts <- c(starts, draw(max(fits, 1L) - length(starts)))
  res   <- prior$fits %||% list()
  pv    <- prior$values %||% numeric(0)
  total <- prior$starts %||% 0L
  repeat {
    out <- if (cores > 1L) parallel::mclapply(starts, one, mc.cores = cores,
                                              mc.preschedule = FALSE)
           else lapply(starts, one)
    total <- total + length(starts)
    res  <- c(res, Filter(function(f) is.list(f) && !is.null(f$value), out))
    hits <- if (length(res)) {
      v <- c(vapply(res, `[[`, 0, "value"), pv)
      sum(v <= min(v) + (wf$tol %||% 0))
    } else 0L
    if (is.null(wf) || hits >= wf$hits || total >= maxN) break
    starts <- draw(min(max(fits, cores, 1L), maxN - total))
  }
  if (!length(res)) return(NULL)
  vals <- vapply(res, `[[`, 0, "value")
  best <- res[[which.min(vals)]]
  best$values <- sort(c(vals, pv))
  best$starts <- total
  best$hits   <- hits
  best$level  <- .l1Level(lapply(res, `[`, c("argument", "value")), levelTol)
  best
}

# The runs within `tol` of the best value: the lowest level of a waterfall.
.l1Level <- function(fits, tol) {
  v <- vapply(fits, `[[`, 0, "value")
  fits[v <= min(v) + tol]
}

# One penalised fit at `lambda` from `start0`, the `extra` starts and random
# ones; q < 1 reweights the L1 fit of every start and scores by the Lq
# objective, `ssl` runs the EM of the spike-and-slab lasso from every start.
.l1Penalised <- function(obj, start0, lambda, gates, reference, groups,
                         fixed, q, ctl, fits, sd, cores, extra = NULL, ssl = NULL,
                         wf = NULL, prior = NULL) {
  singles <- c(gates, reference)
  mu <- stats::setNames(rep(0, length(singles)), singles)
  wS <- stats::setNames(rep(1, length(singles)), singles)
  wB <- lapply(groups, function(b) matrix(1, length(b$pars) + 1L, length(b$pars) + 1L))
  run <- function(st, wS, wB) {
    fuse <- if (length(groups)) Map(function(b, w)
      list(pars = b$pars, anchor = b$anchor, weights = w, lambda = lambda),
      groups, wB)
    args <- c(list(obj, st), ctl$trust, list(fixed = fixed, fuse = fuse))
    if (length(mu)) args <- c(args, list(mu = mu, lambda = lambda * wS[names(mu)],
                                         gate = if (length(gates)) gates))
    do.call(trustL1, args)
  }
  if (!is.null(ssl)) {
    em <- function(st) .sslEM(run, st, singles, groups, lambda, ssl, ctl)
    return(.l1Multistart(em, start0, fits, sd, cores, obj, extra = extra,
                         positive = gates, levelTol = ctl$tolHits))
  }
  if (q == 1)
    return(.l1Multistart(function(st) run(st, wS, wB), start0, fits, sd, cores, obj,
                         extra = extra, positive = gates, wf = wf, prior = prior,
                         levelTol = ctl$tolHits))
  lq <- function(st) {
    best <- run(st, wS, wB)
    for (k in seq_len(ctl$nq)) {
      th <- best$argument
      wS <- q * pmax(abs(th[singles]), ctl$eps)^(q - 1)
      wB <- lapply(groups, function(b) {
        v <- c(th[b$pars], if (is.null(b$anchor)) 0 else b$anchor)
        q * pmax(abs(outer(v, v, `-`)), ctl$eps)^(q - 1)
      })
      nxt <- try(run(th, wS, wB), silent = TRUE)
      if (inherits(nxt, "try-error")) break
      best <- nxt
    }
    d <- .sslTerms(best$argument, singles, groups)
    best$value <- obj(best$argument, fixed = fixed, deriv = FALSE)$value +
      lambda * sum(d^q)
    best
  }
  .l1Multistart(lq, start0, fits, sd, cores, obj, extra = extra, positive = gates,
                wf = wf, prior = prior, levelTol = ctl$tolHits)
}

# Absolute penalised terms of `th`: the singles, then the upper triangle of
# every block's pairwise differences (anchor last), named "a:b".
.sslTerms <- function(th, singles, groups) {
  d <- abs(th[singles])
  for (b in groups) {
    v <- c(th[b$pars], if (!is.null(b$anchor)) stats::setNames(b$anchor, "anchor"))
    ut <- which(upper.tri(diag(length(v))), arr.ind = TRUE)
    d <- c(d, stats::setNames(abs(v[ut[, 1]] - v[ut[, 2]]),
                              paste0(names(v)[ut[, 1]], ":", names(v)[ut[, 2]])))
  }
  d
}

# Log densities of the slab and spike components on the -2 log L scale:
# psi_i(d) = lambda_i / 4 exp(-lambda_i |d| / 2).
.sslLogPsi <- function(d, l) log(l / 4) - l * d / 2

# Inclusion probabilities of terms `d` and the slab share `theta`, iterated to
# their joint fixed point at fixed `d`.
.sslEstep <- function(d, lambda0, ssl, theta = 0.5, maxit = 100L, tol = 1e-10) {
  l1 <- .sslLogPsi(d, ssl$lambda1)
  l0 <- .sslLogPsi(d, lambda0)
  J <- length(d)
  for (i in seq_len(maxit)) {
    p <- stats::plogis(log(theta) - log1p(-theta) + l1 - l0)
    thNew <- (sum(p) + ssl$a - 1) / (J + ssl$a + ssl$b - 2)
    thNew <- min(max(thNew, 1e-12), 1 - 1e-12)
    done <- abs(thNew - theta) < tol
    theta <- thNew
    if (done) break
  }
  p <- stats::plogis(log(theta) - log1p(-theta) + l1 - l0)
  list(p = p, theta = theta)
}

# -2 log prior of the terms and the slab share, up to a constant.
.sslPenalty <- function(d, theta, lambda0, ssl) {
  l1 <- log(theta) + .sslLogPsi(d, ssl$lambda1)
  l0 <- log1p(-theta) + .sslLogPsi(d, lambda0)
  mx <- pmax(l1, l0)
  -2 * sum(mx + log(exp(l1 - mx) + exp(l0 - mx))) -
    2 * ((ssl$a - 1) * log(theta) + (ssl$b - 1) * log1p(-theta))
}

# Terms of the spike-and-slab lasso: the singles, then per block the gaps
# between neighbours of the sorted values, anchor included. `d` holds the
# gaps, `edges` per block the index pairs into members and anchor.
.sslEdges <- function(th, singles, groups) {
  d <- abs(th[singles])
  edges <- lapply(groups, function(b) {
    v <- c(th[b$pars], if (!is.null(b$anchor)) stats::setNames(b$anchor, "anchor"))
    o <- order(v)
    cbind(o[-length(o)], o[-1L])
  })
  for (i in seq_along(groups)) {
    b <- groups[[i]]
    v <- c(th[b$pars], if (!is.null(b$anchor)) stats::setNames(b$anchor, "anchor"))
    e <- edges[[i]]
    d <- c(d, stats::setNames(abs(v[e[, 2]] - v[e[, 1]]),
                              paste0(names(v)[e[, 1]], ":", names(v)[e[, 2]])))
  }
  list(d = d, edges = edges)
}

# Gaps of `th` on given edges.
.sslGaps <- function(th, singles, groups, edges) {
  d <- abs(th[singles])
  for (i in seq_along(groups)) {
    b <- groups[[i]]
    v <- c(th[b$pars], if (!is.null(b$anchor)) b$anchor)
    d <- c(d, abs(v[edges[[i]][, 2]] - v[edges[[i]][, 1]]))
  }
  unname(d)
}

# EM of the spike-and-slab lasso from one start. Spike terms the weighted fit
# leaves off their kink are merged onto it; the merge is kept when it lowers
# the -2 log posterior.
.sslEM <- function(run, st, singles, groups, lambda0, ssl, ctl) {
  fit <- .sslLoop(run, st, singles, groups, lambda0, ssl, ctl)
  snap <- function(f) {
    th <- .sslMerge(f$argument, f$edgeP, singles, groups, ctl$snap)
    if (!is.null(th)) f$argument <- th
    f
  }
  fit <- snap(fit)
  for (m in seq_len(ctl$nmerge)) {
    st2 <- .sslMerge(fit$argument, fit$edgeP, singles, groups)
    if (is.null(st2)) break
    alt <- try(.sslLoop(run, st2, singles, groups, lambda0, ssl, ctl), silent = TRUE)
    if (inherits(alt, "try-error") || alt$value >= fit$value) break
    alt$merges <- (fit$merges %||% 0L) + 1L
    fit <- snap(alt)
  }
  fit
}

# EM from one start: E-step on the current terms, M-step a trustL1 fit with
# weights p lambda1 + (1 - p) lambda0 on the edges. Returns the last fit with
# the -2 log posterior as `value`, `edgeP`, `inclusion` and `theta`.
.sslLoop <- function(run, st, singles, groups, lambda0, ssl, ctl) {
  th <- st
  ed <- .sslEdges(th, singles, groups)
  es <- .sslEstep(ed$d, lambda0, ssl)
  fit <- NULL
  for (k in seq_len(ctl$nem)) {
    w   <- es$p * ssl$lambda1 / lambda0 + (1 - es$p)
    wS  <- w[singles]
    off <- length(singles)
    wB  <- Map(function(b, e) {
      m <- length(b$pars) + 1L
      W <- matrix(0, m, m)
      W[e] <- w[off + seq_len(nrow(e))]
      W[e[, 2:1, drop = FALSE]] <- W[e]
      off <<- off + nrow(e)
      W
    }, groups, ed$edges)
    fit <- run(th, wS, wB)
    th  <- fit$argument
    pen <- lambda0 * sum(w * .sslGaps(th, singles, groups, ed$edges))
    edNew <- .sslEdges(th, singles, groups)
    esNew <- .sslEstep(edNew$d, lambda0, ssl, es$theta)
    done <- identical(edNew$edges, ed$edges) && max(abs(esNew$p - es$p)) < ctl$tolp
    ed <- edNew
    es <- esNew
    if (done) break
  }
  fit$value <- fit$value - pen + .sslPenalty(ed$d, es$theta, lambda0, ssl)
  fit$emIterations <- k
  fit$edgeP <- es$p
  all <- .sslTerms(th, singles, groups)
  fit$inclusion <- stats::plogis(log(es$theta) - log1p(-es$theta) +
                                   .sslLogPsi(all, ssl$lambda1) - .sslLogPsi(all, lambda0))
  fit$theta <- es$theta
  fit
}

# Start with every spike term (p < 1/2) that is off its kink by less than
# `gap` put onto it: singles to zero, block members joined by spike edges to
# their anchor or to their mean. NULL if there is nothing to merge.
.sslMerge <- function(th, p, singles, groups, gap = Inf) {
  ed <- .sslEdges(th, singles, groups)
  spike <- p < 0.5 & ed$d < gap
  if (!any(spike & ed$d > 0)) return(NULL)
  s <- singles[spike[seq_along(singles)]]
  th[s] <- 0
  off <- length(singles)
  for (i in seq_along(groups)) {
    b  <- groups[[i]]
    v  <- c(b$pars, if (!is.null(b$anchor)) "anchor")
    e  <- ed$edges[[i]]
    sp <- spike[off + seq_len(nrow(e))]
    off <- off + nrow(e)
    lab <- seq_along(v)
    repeat {
      old <- lab
      for (j in which(sp)) lab[e[j, ]] <- min(lab[e[j, ]])
      for (j in seq_along(lab)) lab[j] <- lab[lab[j]]
      if (identical(old, lab)) break
    }
    for (cl in unique(lab)) {
      mem <- v[lab == cl]
      if (length(mem) < 2L) next
      pars <- setdiff(mem, "anchor")
      th[pars] <- if ("anchor" %in% mem) b$anchor else mean(th[pars])
    }
  }
  th
}

# Removed parameters and fused groups of a penalised optimum, as exact
# equalities, plus a key that identifies the structure.
.l1Structure <- function(th, gates, reference, groups) {
  removed <- c(gates[th[gates] == 0], reference[th[reference] == 0])
  cls <- lapply(groups, function(b) {
    v <- th[b$pars]
    split(b$pars, match(v, unique(v)))
  })
  anch <- lapply(groups, function(b) {
    if (is.null(b$anchor)) character(0) else b$pars[th[b$pars] == b$anchor]
  })
  grpKey <- vapply(seq_along(cls), function(i) {
    g <- vapply(cls[[i]], function(z) {
      lab <- paste(sort(z), collapse = ",")
      if (all(z %in% anch[[i]])) paste0(lab, "=", groups[[i]]$anchor) else lab
    }, "")
    paste0("{", paste(sort(g), collapse = " | "), "}")
  }, "")
  rm <- if (length(removed)) paste0("-", sort(removed))
  list(removed = removed, groups = unname(cls), anchored = anch,
       anchor = lapply(groups, `[[`, "anchor"),
       key = paste(c(rm, grpKey), collapse = " "))
}

# Unpenalised refit of one structure from the full optimum and from `start`,
# the penalised optimum. Returns value, number of free parameters and the full
# named parameter vector.
.l1Refit <- function(st, start, obj, fullPars, zero, fixed, ctl, fits, sd, cores,
                     wf = NULL) {
  gates   <- names(zero)
  offGate <- intersect(st$removed, gates)
  fixR <- c(fixed,
            stats::setNames(rep(1, length(setdiff(gates, offGate))), setdiff(gates, offGate)),
            stats::setNames(rep(0, length(offGate)), offGate),
            fullPars[unique(unlist(zero[offGate]))])
  offRef <- setdiff(st$removed, gates)
  fixR <- c(fixR, stats::setNames(rep(0, length(offRef)), offRef))
  ties <- list()
  for (i in seq_along(st$groups)) {
    for (cl in st$groups[[i]]) {
      if (all(cl %in% st$anchored[[i]])) {
        fixR <- c(fixR, stats::setNames(rep(st$anchor[[i]], length(cl)), cl))
        next
      }
      if (length(cl) > 1L) ties[[cl[1]]] <- cl
    }
  }
  tied  <- unlist(lapply(ties, `[`, -1L))
  free  <- setdiff(names(fullPars), c(names(fixR), tied))
  objR  <- .l1Tie(obj, ties)
  warm  <- fullPars
  both  <- intersect(names(start), names(warm))
  warm[both] <- start[both]
  r <- .l1Multistart(function(s)
    do.call(trust, c(list(objR, s), .l1TrustArgs(ctl$trust), list(fixed = fixR))),
                     fullPars[free], fits, sd, cores, obj, extra = list(warm[free]),
                     wf = wf)
  if (is.null(r)) return(NULL)
  arg <- r$argument
  for (rep in names(ties)) arg[ties[[rep]]] <- arg[[rep]]
  list(value = r$value, nfree = length(free), argument = c(arg, fixR),
       starts = r$starts, hits = r$hits)
}

# Objective over one representative per tied group: gradient and Hessian are
# pulled back through the indicator map.
.l1Tie <- function(obj, ties) {
  if (!length(ties)) return(obj)
  function(pars, fixed = NULL, deriv = TRUE, ...) {
    full <- pars
    for (rep in names(ties)) full[ties[[rep]]] <- pars[[rep]]
    out <- obj(full, fixed = fixed, deriv = deriv, ...)
    if (!deriv || is.null(out$gradient)) return(out)
    g  <- out$gradient
    M  <- matrix(0, length(g), length(pars), dimnames = list(names(g), names(pars)))
    for (j in names(pars)) M[intersect(c(j, ties[[j]]), names(g)), j] <- 1
    out$gradient <- drop(crossprod(M, g))
    if (!is.null(out$hessian)) out$hessian <- crossprod(M, out$hessian %*% M)
    out
  }
}

.l1Ndata <- function(obj) {
  d <- attr(obj, "data", exact = TRUE)
  if (is.null(d)) return(NULL)
  sum(vapply(d, NROW, 0L))
}

.l1Select <- function(pathTab, refitTab, select, alpha) {
  if (select == "bic") {
    if (all(is.na(refitTab$bic)))
      stop("scanL1: BIC needs the number of data points, set control$ndata.",
           call. = FALSE)
    key <- refitTab$key[which.min(refitTab$bic)]
    return(list(key = key, lambda = max(pathTab$lambda[pathTab$key == key])))
  }
  if (select == "plateau") {
    l <- which.max(pathTab$lambda)
    return(list(key = pathTab$key[l], lambda = pathTab$lambda[l]))
  }
  ok <- which(!is.na(pathTab$p) & pathTab$p >= alpha)
  l  <- if (length(ok)) max(ok) else which.min(pathTab$lambda)
  list(key = pathTab$key[l], lambda = pathTab$lambda[l])
}


#' @export
print.scanL1 <- function(x, ...) {
  cat(sprintf("scanL1: %d lambdas, %d distinct structures, selection by %s",
              nrow(x$path), nrow(x$refits), toupper(x$select)))
  if (x$select == "lrt") cat(sprintf(" (alpha = %g)", x$alpha))
  if (x$q < 1) cat(sprintf(", q = %g", x$q))
  if (!is.null(x$ssl)) cat(sprintf(", spike-and-slab with lambda1 = %g", x$ssl$lambda1))
  if (!is.null(x$em)) {
    cat(sprintf(", lambda by EM: %s\n\nTerms:\n",
                paste(sprintf("%s %.4g", names(x$em$lambda), x$em$lambda), collapse = ", ")))
    print(x$em$terms, row.names = FALSE, digits = 4)
  }
  rt  <- x$refits
  ids <- paste0("S", seq_len(nrow(rt)))
  cat("\n\nRefits:\n")
  print(data.frame(id = ids, rt[setdiff(names(rt), "key")]), row.names = FALSE, digits = 4)
  cat("\nStructures:\n")
  cat(sprintf("  %s  %s\n", ids, ifelse(nzchar(rt$key), rt$key, "(full model)")), sep = "")
  sel <- if (is.null(x$em)) sprintf(" at lambda = %.4g", x$lambdaSelected) else ""
  cat(sprintf("\nSelected%s: %s\n", sel, ids[match(x$selected, rt$key)]))
  invisible(x)
}

#' Plot an L1 scan
#'
#' @param x A `scanL1` result.
#' @param type `"path"`: penalised estimates over `lambda`. `"test"`: p-value of
#'   each refit over `lambda`, with the level of the test. `"waterfall"`: sorted
#'   values of the starts of the full model. `"clusters"`: per
#'   block, the group every member (and the anchor) belongs to at each
#'   `lambda`. `"inclusion"`: inclusion probabilities of the spike-and-slab
#'   lasso over `lambda`.
#' @param ... Ignored.
#' @return A ggplot.
#' @export
plot.scanL1 <- function(x, type = c("path", "test", "waterfall", "clusters", "inclusion"),
                        ...) {
  type <- match.arg(type)
  vline <- ggplot2::geom_vline(xintercept = x$lambdaSelected, linetype = 2)
  if (type == "waterfall") {
    v <- x$full$values
    return(ggplot2::ggplot(data.frame(index = seq_along(v), value = v - v[1]),
                           ggplot2::aes(index, value)) +
             ggplot2::geom_point() +
             ggplot2::scale_y_continuous(trans = "log1p", breaks = c(0, 10^(0:6))) +
             ggplot2::labs(x = "start, sorted",
                           y = "-2 log L of the full model above the best") +
             theme_dMod())
  }
  if (type == "path") {
    cf <- x$coefficients
    df <- data.frame(lambda = rep(cf[, "lambda"], ncol(cf) - 1L),
                     parameter = rep(colnames(cf)[-1], each = nrow(cf)),
                     value = as.vector(cf[, -1, drop = FALSE]))
    return(ggplot2::ggplot(df, ggplot2::aes(lambda, value, colour = parameter)) +
             ggplot2::geom_line() + ggplot2::geom_point(size = 0.8) +
             ggplot2::scale_x_log10() + vline +
             ggplot2::labs(x = "lambda", y = "penalised estimate") +
             theme_dMod())
  }
  if (type == "clusters") {
    if (!length(x$groups)) stop("plot.scanL1: no blocks to cluster.", call. = FALSE)
    df <- .l1ClusterTable(x)
    br <- seq(floor(log10(min(df$lambda))), ceiling(log10(max(df$lambda))))
    br <- br[seq(1, length(br), by = ceiling(length(br) / 3))]
    return(ggplot2::ggplot(df, ggplot2::aes(log10(lambda), member, fill = cluster)) +
             ggplot2::geom_tile(colour = "white") +
             ggplot2::facet_wrap(~block, scales = "free_y") +
             ggplot2::scale_x_continuous(breaks = br, labels = paste0("1e", br)) +
             ggplot2::geom_vline(xintercept = log10(x$lambdaSelected), linetype = 2) +
             ggplot2::labs(x = "lambda", y = NULL, fill = "group") +
             theme_dMod())
  }
  if (type == "inclusion") {
    if (is.null(x$inclusion)) stop("plot.scanL1: not a spike-and-slab scan.", call. = FALSE)
    ic <- x$inclusion
    df <- data.frame(lambda = rep(ic[, "lambda"], ncol(ic) - 1L),
                     term = rep(colnames(ic)[-1], each = nrow(ic)),
                     p = as.vector(ic[, -1, drop = FALSE]))
    return(ggplot2::ggplot(df, ggplot2::aes(lambda, p, group = term)) +
             ggplot2::geom_line(alpha = 0.6) +
             ggplot2::scale_x_log10() + vline +
             ggplot2::labs(x = "lambda", y = "inclusion probability") +
             theme_dMod())
  }
  ggplot2::ggplot(x$path, ggplot2::aes(lambda, p)) +
    ggplot2::geom_step() + ggplot2::geom_point(size = 0.8) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::geom_hline(yintercept = x$alpha, linetype = 2) + vline +
    ggplot2::labs(x = "lambda", y = "p-value of the refit against the full model") +
    theme_dMod()
}

# Long table of group labels per block, member and lambda. Groups are numbered
# by first appearance with the anchor first, so the anchor group is always 1.
.l1ClusterTable <- function(x) {
  cf <- x$coefficients
  bn <- names(x$groups) %||% paste0("block", seq_along(x$groups))
  bn[!nzchar(bn)] <- paste0("block", which(!nzchar(bn)))
  do.call(rbind, lapply(seq_along(x$groups), function(i) {
    b <- x$groups[[i]]
    do.call(rbind, lapply(seq_len(nrow(cf)), function(l) {
      v <- c(if (!is.null(b$anchor)) c(anchor = b$anchor), cf[l, b$pars])
      data.frame(lambda = unname(cf[l, "lambda"]), block = bn[i], member = names(v),
                 cluster = factor(match(v, unique(v))), stringsAsFactors = FALSE)
    }))
  }))
}
