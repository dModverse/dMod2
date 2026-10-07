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
#' @param select `"lrt"`: the largest `lambda` whose refit the likelihood ratio
#'   test against the full model does not reject at level `alpha`. `"bic"`: the
#'   refit with the smallest BIC.
#' @param alpha Level of the likelihood ratio test.
#' @param fits Number of starts per `lambda` and per refit. For a `lambda`
#'   they include the optimum of the next smaller `lambda` and the sparse point
#'   (every gate and reference at zero, every block at its anchor or mean); a
#'   downward pass then restarts each `lambda` from the optimum of the next
#'   larger one and keeps the better fit. A refit starts from the full optimum.
#'   The remaining starts are drawn around the full optimum.
#' @param cores Number of forked processes for the starts.
#' @param sd Standard deviation of the start perturbation.
#' @param control List: `trust` (arguments of [trustL1] and [trust], default
#'   `rinit = 0.1, rmax = 10, iterlim = 200`), `nq` (reweighting rounds for
#'   `q < 1`, default 3), `eps` (floor of `|difference|` in the reweighting,
#'   default 0.01), `ndata` (data points for BIC, default from `obj`).
#' @return Object of class `scanL1` with
#'   \describe{
#'     \item{`path`}{one row per `lambda`: penalised value, structure key,
#'       number of removed parameters and of free parameters in the refit.}
#'     \item{`coefficients`}{penalised estimates per `lambda`.}
#'     \item{`refits`}{one row per distinct structure: `-2 log L`, free
#'       parameters, LRT statistic, degrees of freedom, p-value, BIC.}
#'     \item{`full`}{the full fit: value and parameters.}
#'     \item{`selected`, `lambdaSelected`}{the chosen structure and the
#'       `lambda` it was chosen at.}
#'     \item{`structure`}{per key: removed parameters and groups.}
#'     \item{`fit`}{named parameters of the chosen refit, fixed ones included.}
#'   }
#' @references Hauber AL, Rosenblatt M, Timmer J (2023). Uncovering specific
#'   mechanisms across cell types in dynamical models. PLoS Comput Biol 19(9):
#'   e1010867.
#' @seealso [gateL1], [trustL1]
#' @example inst/examples/scanL1.R
#' @export
scanL1 <- function(obj, center, zero = NULL, reference = NULL, groups = NULL,
                   lambda = 10^seq(-3, 3, length.out = 25), fixed = NULL,
                   q = 1, select = c("lrt", "bic"), alpha = 0.05,
                   fits = 10, cores = 1, sd = 0.5, control = list()) {
  select <- match.arg(select)
  if (!length(zero) && !length(reference) && !length(groups))
    stop("scanL1: give at least one of zero, reference, groups.", call. = FALSE)
  if (q <= 0 || q > 1) stop("scanL1: q must lie in (0, 1].", call. = FALSE)
  ctl <- utils::modifyList(list(trust = list(rinit = 0.1, rmax = 10, iterlim = 200L),
                                nq = 3L, eps = 0.01, ndata = NULL), control)
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
  full <- .l1Multistart(function(st) do.call(trust, c(list(obj, st), ctl$trust,
                                                      list(fixed = fixFull))),
                        center, fits, sd, cores, obj)
  if (is.null(full)) stop("scanL1: every fit of the full model failed.", call. = FALSE)
  fullPars <- full$argument
  nFull    <- length(fullPars)

  ## 2. Penalised path --------------------------------------------------------
  gatedPars <- unique(unlist(zero))
  fixSel    <- c(fixed, fullPars[gatedPars])
  start0    <- c(fullPars[setdiff(names(fullPars), gatedPars)],
                 stats::setNames(rep(1, length(gates)), gates))
  lambda    <- sort(lambda)
  sparse    <- start0
  sparse[c(gates, reference)] <- 0
  for (b in groups) sparse[b$pars] <- b$anchor %||% mean(start0[b$pars])
  pen1 <- function(start, l, n, extra)
    .l1Penalised(obj, start, lambda[l], gates, reference, groups, fixSel, q,
                 ctl, n, sd, cores, extra)
  # Upward from the full optimum, then downward from each larger lambda's
  # optimum; a lambda keeps the better of the two.
  fitsL <- vector("list", length(lambda))
  warm  <- start0
  for (l in seq_along(lambda)) {
    fitsL[[l]] <- pen1(start0, l, fits, list(warm, sparse))
    if (!is.null(fitsL[[l]])) warm <- fitsL[[l]]$argument
  }
  for (l in rev(seq_along(lambda))[-1]) {
    up <- fitsL[[l + 1L]]
    if (is.null(up)) next
    down <- pen1(up$argument, l, 1L, NULL)
    if (!is.null(down) && (is.null(fitsL[[l]]) || down$value < fitsL[[l]]$value))
      fitsL[[l]] <- down
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
  rf <- lapply(structs, .l1Refit, obj = obj, fullPars = fullPars, zero = zero,
               fixed = fixed, ctl = ctl, fits = fits, sd = sd, cores = cores)
  rf <- Filter(Negate(is.null), rf)
  if (!length(rf)) stop("scanL1: every refit failed.", call. = FALSE)
  refitTab <- data.frame(key = names(rf),
                         value = vapply(rf, `[[`, 0, "value"),
                         nfree = vapply(rf, `[[`, 0L, "nfree"),
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
  sel <- .l1Select(pathTab, refitTab, select, alpha)

  coefs <- do.call(rbind, lapply(path, function(z) z$argument[c(gates, pen)]))
  rownames(coefs) <- NULL

  out <- list(path = pathTab, coefficients = cbind(lambda = pathTab$lambda, coefs),
              refits = refitTab, full = list(value = full$value, argument = fullPars),
              selected = sel$key, lambdaSelected = sel$lambda,
              structure = structs, select = select, alpha = alpha, q = q,
              fit = rf[[sel$key]]$argument)
  class(out) <- "scanL1"
  out
}


# Best of `fits` runs of `run(start)`: start 1 is `center`, the others are
# drawn around it. `abs()` keeps gates non-negative.
.l1Multistart <- function(run, center, fits, sd, cores, obj, extra = NULL,
                          positive = character(0)) {
  starts <- c(list(center), extra)
  while (length(starts) < max(fits, 1L)) {
    st <- center + stats::rnorm(length(center), 0, sd)
    st[positive] <- abs(st[positive])
    starts <- c(starts, list(st))
  }
  one <- function(st) {
    try(resetWarmStarts(obj, verbose = FALSE), silent = TRUE)
    # A start the solver cannot follow is dropped; its warnings say nothing more.
    f <- try(suppressWarnings(run(st)), silent = TRUE)
    if (inherits(f, "try-error") || !is.finite(f$value)) NULL else f
  }
  res <- if (cores > 1L) parallel::mclapply(starts, one, mc.cores = cores,
                                            mc.preschedule = FALSE)
         else lapply(starts, one)
  res <- Filter(function(f) is.list(f) && !is.null(f$value), res)
  if (!length(res)) return(NULL)
  res[[which.min(vapply(res, `[[`, 0, "value"))]]
}

# One penalised fit at `lambda` from `start0`, the `extra` starts and random
# ones around `start0`; q < 1 by reweighting the best L1 fit.
.l1Penalised <- function(obj, start0, lambda, gates, reference, groups,
                         fixed, q, ctl, fits, sd, cores, extra = NULL) {
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
  best <- .l1Multistart(function(st) run(st, wS, wB), start0, fits, sd, cores, obj,
                        extra = extra, positive = gates)
  if (is.null(best) || q == 1) return(best)
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
  best
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

# Unpenalised refit of one structure. Returns value, number of free parameters
# and the full named parameter vector.
.l1Refit <- function(st, obj, fullPars, zero, fixed, ctl, fits, sd, cores) {
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
  r <- .l1Multistart(function(s) do.call(trust, c(list(objR, s), ctl$trust,
                                                  list(fixed = fixR))),
                     fullPars[free], fits, sd, cores, obj)
  if (is.null(r)) return(NULL)
  arg <- r$argument
  for (rep in names(ties)) arg[ties[[rep]]] <- arg[[rep]]
  list(value = r$value, nfree = length(free), argument = c(arg, fixR))
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
  cat("\n\nRefits:\n")
  print(x$refits, row.names = FALSE, digits = 4)
  cat(sprintf("\nSelected at lambda = %.4g:\n  %s\n", x$lambdaSelected,
              if (nzchar(x$selected)) x$selected else "(full model)"))
  invisible(x)
}

#' Plot an L1 scan
#'
#' @param x A `scanL1` result.
#' @param type `"path"`: penalised estimates over `lambda`. `"test"`: p-value of
#'   each refit over `lambda`, with the level of the test.
#' @param ... Ignored.
#' @return A ggplot.
#' @export
plot.scanL1 <- function(x, type = c("path", "test"), ...) {
  type <- match.arg(type)
  if (type == "path") {
    cf <- x$coefficients
    df <- data.frame(lambda = rep(cf[, "lambda"], ncol(cf) - 1L),
                     parameter = rep(colnames(cf)[-1], each = nrow(cf)),
                     value = as.vector(cf[, -1, drop = FALSE]))
    return(ggplot2::ggplot(df, ggplot2::aes(lambda, value, colour = parameter)) +
             ggplot2::geom_line() + ggplot2::geom_point(size = 0.8) +
             ggplot2::scale_x_log10() +
             ggplot2::geom_vline(xintercept = x$lambdaSelected, linetype = 2) +
             ggplot2::labs(x = "lambda", y = "penalised estimate") +
             theme_dMod())
  }
  ggplot2::ggplot(x$path, ggplot2::aes(lambda, p)) +
    ggplot2::geom_step() + ggplot2::geom_point(size = 0.8) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::geom_hline(yintercept = x$alpha, linetype = 2) +
    ggplot2::geom_vline(xintercept = x$lambdaSelected, linetype = 2) +
    ggplot2::labs(x = "lambda", y = "p-value of the refit against the full model") +
    theme_dMod()
}
