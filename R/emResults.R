# The result of a mixed-effects fit and its methods: Omega and the random effects with standard
# errors, log-likelihood, residuals, predictions and confidence intervals.


## The `em` object of a fit of `design`; `fit` holds argument, value, gradient, hessian, the
## modes of the random effects per group and their conditional covariances, or with several
## groupings the modes as one vector and the covariances per component.
.emResult <- function(design, fit, method, fixed) {
  pars <- c(fit$argument, fixed)
  nObs <- sum(vapply(design$data, nrow, 0L))
  # censored rows enter by -2 log Phi, without the 2 pi of a density
  nDensity <- sum(vapply(names(design$data), function(cn) {
    d <- .conditionData(design, cn)
    sum(d$value > d$lloq)
  }, 0))
  out <- list(argument = fit$argument, value = fit$value, gradient = fit$gradient,
              hessian = fit$hessian, information = fit$information, converged = fit$converged,
              iterations = fit$iterations, trustIterations = fit$trustIterations,
              method = method, fixed = fixed,
              nObs = nObs, valueNonmem = if (!is.finite(design$df)) fit$value - nDensity * log(2 * pi) else NA_real_,
              prdfn = design$prdfn, data = design$data, errfn = design$errfn, design = design,
              trace = fit$trace, stageTrace = fit$stageTrace, certificate = fit$certificate,
              valueSE = fit$valueSE, coupling = fit$coupling, derivMode = fit$derivMode,
              hessianExact = fit$hessianExact)
  if (is.null(design$groupings)) {
    gr <- .groupingResult(design, pars, .emCurvature(fit), fit$etaModes, fit$etaCov, fit$clusters, fit$spike)
    out <- c(out, gr[c("Omega", "etaModes", "etaCov", "etaSE", "shrinkage", "omega", "random")])
    out$eta <- .etaVector(list(design = design, etaModes = fit$etaModes))
  } else {
    per <- lapply(design$groupings, function(sd) {
      modes <- matrix(fit$etas[sd$etaNames], nrow(sd$etaNames), dimnames = dimnames(sd$etaNames))
      cov <- lapply(sd$groups, function(g) {
        en <- sd$etaNames[g, ]
        V <- fit$etaCovUnits[[which(vapply(fit$etaCovUnits, function(V) all(en %in% rownames(V)), TRUE))[1]]]
        V[en, en, drop = FALSE]
      })
      .groupingResult(sd, pars, fit$hessian, modes, setNames(cov, sd$groups))
    })
    for (f in c("Omega", "etaModes", "etaCov", "etaSE", "shrinkage")) out[[f]] <- lapply(per, `[[`, f)
    for (f in c("omega", "random")) out[[f]] <- do.call(rbind, Map(function(t, g) cbind(grouping = g, t),
                                                                   lapply(per, `[[`, f), names(per)))
    out$eta <- fit$etas
  }
  out$strengths <- .strengthTable(design, pars, names(fixed))
  out$slab <- .slabTable(design, pars, names(fixed))
  class(out) <- c("em", "list")
  out
}

## Omega, modes, conditional standard errors, shrinkage and the tables of one grouping.
.groupingResult <- function(sd, pars, hessian, etaModes, etaCov, clusters = NULL, spike = NULL) {
  lay <- sd$normal
  Omega <- if (lay$K) tcrossprod(.choleskyFactor(lay, pars))
  info <- .etaInfo(etaCov, Omega, etaModes)
  list(Omega = Omega, etaModes = etaModes, etaCov = etaCov, etaSE = info$se, shrinkage = info$shrinkage,
       omega = .omegaTable(sd, pars, hessian), random = .randomTable(sd, etaModes, info$se, clusters, spike))
}

## The modes of every random effect as one named vector, eta names as in the model.
.etaVector <- function(object) {
  if (!is.null(object[["eta", exact = TRUE]])) return(object[["eta", exact = TRUE]])
  E <- object$design$etaNames[, colnames(object$etaModes), drop = FALSE]
  setNames(as.vector(object$etaModes[rownames(E), , drop = FALSE]), as.vector(E))
}

## Conditional standard errors of the random effects from their conditional covariances per
## group, and the shrinkage per random effect, 1 - SD(modes) / SD(Omega).
.etaInfo <- function(etaCov, Omega, etaModes) {
  if (is.null(etaModes)) return(list(se = NULL, shrinkage = NULL))
  shrinkage <- NULL
  if (!is.null(Omega)) {
    sd <- sqrt(pmax(diag(Omega), 0))[colnames(etaModes)]
    shrinkage <- setNames(1 - apply(etaModes, 2, stats::sd) / sd, colnames(etaModes))
    shrinkage[!is.finite(sd) | sd <= 0] <- NA_real_
  }
  if (is.null(etaCov)) return(list(se = NULL, shrinkage = shrinkage))
  v <- vapply(etaCov, function(C) pmax(diag(C), 0), numeric(ncol(etaModes)))
  se <- matrix(sqrt(v), nrow(etaModes), byrow = TRUE, dimnames = dimnames(etaModes))
  list(se = se, shrinkage = shrinkage)
}


## Omega as standard deviations and correlations with standard errors by the delta method on
## the covariance of the Cholesky parameters.
.omegaTable <- function(design, pars, hessian) {
  lay <- design$normal
  if (!lay$K) return(NULL)
  L <- .choleskyFactor(lay, pars)
  O <- tcrossprod(L)
  K <- lay$K
  sd <- sqrt(diag(O))
  loc <- unname(lay$loc)
  dg <- unname(lay$isDiag)
  dO <- lapply(seq_along(lay$cholPars), function(m) {
    d <- matrix(0, K, K)
    d[loc[m, 1L], loc[m, 2L]] <- if (dg[m]) L[loc[m, , drop = FALSE]] else 1
    d %*% t(L) + L %*% t(d)
  })
  pairs <- which(lower.tri(O, diag = TRUE), arr.ind = TRUE)
  pairs <- pairs[order(pairs[, 1] != pairs[, 2], pairs[, 2], pairs[, 1]), , drop = FALSE]
  k <- pairs[, 1]; l <- pairs[, 2]
  isSd <- k == l
  est <- ifelse(isSd, sd[k], O[cbind(k, l)] / (sd[k] * sd[l]))
  # d sd_k = dO_kk / (2 sd_k); d cor_kl = dO_kl / (sd_k sd_l) - cor_kl (dO_kk / (2 O_kk) + dO_ll / (2 O_ll))
  J <- vapply(dO, function(d) ifelse(isSd, d[cbind(k, k)] / (2 * sd[k]),
                                      d[cbind(k, l)] / (sd[k] * sd[l]) -
                                        est * (d[cbind(k, k)] / (2 * O[cbind(k, k)]) +
                                               d[cbind(l, l)] / (2 * O[cbind(l, l)]))),
              numeric(length(k)))
  J <- matrix(J, length(k))
  se <- rep(NA_real_, length(k))
  # fixed Cholesky entries have no variance
  cp <- intersect(lay$cholPars, rownames(hessian))
  if (!is.null(hessian) && length(cp)) {
    V <- tryCatch(if (.isSparse(hessian)) .vcovSparse(hessian, NULL, cp) else solve(hessian / 2),
                  error = function(e) NULL)
    Jf <- J[, match(cp, lay$cholPars), drop = FALSE]
    if (!is.null(V)) se <- sqrt(pmax(rowSums((Jf %*% V[cp, cp, drop = FALSE]) * Jf), 0))
  }
  keep <- isSd | O[cbind(k, l)] != 0
  data.frame(type = ifelse(isSd, "sd", "cor"), eta1 = lay$short[k], eta2 = lay$short[l],
             estimate = est, se = se, stringsAsFactors = FALSE)[keep, , drop = FALSE]
}

## The strengths of the Laplace blocks, estimated or fixed in the design or by name in `fixed`.
.strengthTable <- function(design, pars, fixed = NULL) {
  blocks <- Filter(function(b) b$dist == "laplace", design$blocks)
  if (!length(blocks)) return(NULL)
  data.frame(name = vapply(blocks, `[[`, "", "hyper"),
             symbols = vapply(blocks, function(b) paste(b$short, collapse = ", "), ""),
             estimate = vapply(blocks, function(b) if (is.na(b$lambda)) pars[[b$hyper]] else b$lambda, 0),
             fixed = !vapply(blocks, function(b) is.na(b$lambda), TRUE) |
               vapply(blocks, function(b) b$hyper, "") %in% fixed, stringsAsFactors = FALSE)
}

## The hyperparameters of the spike-and-slab and cluster blocks: the probability of a deviation
## or the concentration, and the slab standard deviation.
.slabTable <- function(design, pars, fixed = NULL) {
  blocks <- Filter(function(b) b$dist %in% c("spikeSlab", "dirichletProcess"), design$blocks)
  if (!length(blocks)) return(NULL)
  value <- function(b, kind) if (is.na(b[[kind]])) pars[[b$hyper[[kind]]]] else b[[kind]]
  do.call(rbind, lapply(blocks, function(b) data.frame(
    name = unname(b$hyper), distribution = b$dist, symbols = paste(b$short, collapse = ", "),
    estimate = vapply(names(b$hyper), function(k) value(b, k), 0),
    fixed = !is.na(unlist(b[names(b$hyper)])) | unname(b$hyper) %in% fixed, row.names = NULL,
    stringsAsFactors = FALSE)))
}

## The random effects per group: mode and conditional standard error; for spike-and-slab and
## Dirichlet-process ones the selected cluster (0 the fixed effect) and the probability of 0.
.randomTable <- function(design, etaModes, etaSE, clusters = NULL, spike = NULL) {
  if (is.null(etaModes)) return(NULL)
  groups <- rownames(etaModes)
  E <- design$etaNames[groups, colnames(etaModes), drop = FALSE]
  out <- data.frame(group = rep(groups, ncol(etaModes)),
                    eta = rep(sub("^eta_", "", colnames(etaModes)), each = length(groups)),
                    parameter = as.vector(E), mode = as.vector(etaModes),
                    se = if (!is.null(etaSE)) as.vector(etaSE) else NA_real_,
                    stringsAsFactors = FALSE)
  if (is.null(clusters) || all(is.na(clusters))) return(out)
  out$cluster <- as.vector(clusters[groups, colnames(etaModes), drop = FALSE])
  out$p0 <- as.vector(spike[groups, colnames(etaModes), drop = FALSE])
  out
}


#' Print a mixed-effects fit
#'
#' @param x An `em` object, see [em()].
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.em <- function(x, ...) {
  cat(sprintf("EM (%s): -2 log L %.6f, converged %s, %s iterations, %d groups\n",
              x$method, x$value, x$converged, format(x$iterations %||% NA_integer_),
              sum(vapply(.subDesigns(x$design), function(sd) length(sd$groups), 0L))))
  if (!is.null(x$valueSE))
    cat(sprintf("-2 log L by importance sampling, standard error %.3g\n", x$valueSE))
  if (!is.null(x$coupling))
    cat("Coupled sparse random effects: ", switch(x$coupling, expectationConsistent = "expectation consistent",
        meanField = "mean field (an upper bound on -2 log L)",
        meanFieldFallback = "mean field, expectation consistent did not converge at the start"), "\n", sep = "")
  print(x$argument)
  if (length(x$boundary))
    cat("Set to 0 at the boundary: ", paste(x$boundary, collapse = ", "), "\n", sep = "")
  if (!is.null(x$omega)) {
    cat("Random effects:\n")
    print(x$omega, row.names = FALSE)
  }
  if (!is.null(x$strengths)) {
    cat("Laplace strengths:\n")
    print(x$strengths, row.names = FALSE)
  }
  if (!is.null(x$slab)) {
    kinds <- c(spikeSlab = "Spike-and-slab", dirichletProcess = "Dirichlet-process")
    cat(paste(kinds[unique(x$slab$distribution)], collapse = " and "), "hyperparameters:\n")
    print(x$slab, row.names = FALSE)
  }
  if (!is.null(x$certificate)) {
    ce <- x$certificate
    cat(sprintf("Quadrature: levels %d to %d, estimated level error %.3g in -2 log L; importance sampling flags %d of %d groups\n",
                min(ce$level), max(ce$level), sum(ce$levelError, na.rm = TRUE), sum(ce$flag), nrow(ce)))
    if (any(ce$flag)) print(ce[ce$flag, ], row.names = FALSE)
  }
  invisible(x)
}


#' Estimates of a mixed-effects fit
#'
#' @param object An `em` object, see [em()].
#' @param ... Ignored.
#' @return The named estimates of the fixed parameters and hyperparameters.
#' @export
coef.em <- function(object, ...) object$argument


#' Covariance of the estimates of a mixed-effects fit
#'
#' @param object An `em` object, see [em()].
#' @param ... `parupper`, `parlower` and `which`, passed to [vcov()][vcov.trustfit].
#' @return The inverse of half the Hessian of the marginal -2 log L. A fit with
#'   Laplace or spike-and-slab random effects has no Hessian; there it is the
#'   inverse of the empirical information `information` of the group scores,
#'   see [em()], which needs more groups than parameters. A fit
#'   with Dirichlet-process random effects, whose groups are not independent,
#'   has neither: an error, as where the scores do not span every direction.
#' @export
vcov.em <- function(object, ...) {
  H <- .emCurvature(object)
  if (is.null(H) && !is.null(object$information))
    stop("vcov: the scores of the ", length(object$design$groups), " groups do not span every ",
         "direction of the ", length(object$argument), " parameters; confint(method = \"profile\") gives ",
         "intervals.", call. = FALSE)
  if (is.null(H))
    stop("vcov: the fit has no Hessian of the marginal -2 log L and, with Dirichlet-process random effects, ",
         "no empirical information; confint(method = \"profile\") gives intervals.", call. = FALSE)
  object$hessian <- H
  vcov.trustfit(object, ...)
}

## The Hessian of -2 log L of a fit or, where it has none, twice its empirical information when
## that is positive definite; NULL otherwise.
.emCurvature <- function(fit) {
  if (!is.null(fit$hessian)) return(fit$hessian)
  I <- fit$information
  if (is.null(I) || !nrow(I)) return(if (!is.null(I)) 2 * I)
  ev <- eigen(I, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) > 1e-10 * max(ev)) 2 * I
}

## Standard errors of the estimates, NA where the covariance has no positive diagonal; from a
## sparse Hessian the diagonal of the covariance alone.
.structuralSE <- function(object) {
  est <- object$argument
  se <- setNames(rep(NA_real_, length(est)), names(est))
  if (.isSparse(object$hessian)) {
    d <- tryCatch(.sparseVarianceDiagonal(object$hessian, names(est)), error = function(e) NULL)
    if (!is.null(d)) se[names(d)] <- sqrt(ifelse(d < 0, NA_real_, d))
    return(se)
  }
  V <- tryCatch(vcov(object), error = function(e) NULL)
  if (is.null(V) || is.null(rownames(V))) return(se)
  cn <- intersect(rownames(V), names(est))
  d <- diag(V)[cn]
  d[d < 0] <- NA_real_
  se[cn] <- sqrt(d)
  se
}


#' Summarise a mixed-effects fit
#'
#' @param object An `em` object, see [em()].
#' @param ... Ignored.
#' @return An object of class `summary.em`: `method`, `value`, `valueNonmem`,
#'   `converged`, `iterations`, `population` (estimates with standard errors
#'   and relative standard errors in percent), `seMethod` (`"hessian"`,
#'   `"scores"` for the empirical information of the groups, see
#'   [vcov()][vcov.em], or `"none"`), `partialCurvature` (`TRUE` where FOCEI by
#'   forward over reverse reports its Hessian without the curvature of the
#'   prediction in the random effects within the volume term), `omega` and
#'   `shrinkage` per random effect.
#' @seealso [em()], [confint.em()]
#' @export
summary.em <- function(object, ...) {
  est <- object$argument
  se <- .structuralSE(object)
  pop <- data.frame(estimate = unname(est), se = unname(se), rsePct = unname(100 * se / abs(est)),
                    row.names = names(est))
  seMethod <- if (!is.null(object$hessian)) "hessian" else if (!is.null(.emCurvature(object))) "scores" else "none"
  partial <- identical(seMethod, "hessian") && identical(object$derivMode, "forward-reverse") &&
    isFALSE(object$hessianExact)
  out <- list(method = object$method, value = object$value, valueNonmem = object$valueNonmem,
              converged = object$converged, iterations = object$iterations, population = pop, seMethod = seMethod,
              partialCurvature = partial,
              omega = object$omega,
              shrinkage = object$shrinkage)
  class(out) <- c("summary.em", "list")
  out
}

#' @export
#' @rdname summary.em
#' @param x A `summary.em` object.
#' @param digits Significant digits. Default 4.
print.summary.em <- function(x, digits = 4, ...) {
  cat(sprintf("EM summary (%s): -2 log L %.6f, converged %s, %s iterations\n", x$method, x$value,
              x$converged, format(x$iterations %||% NA_integer_)))
  if (!is.na(x$valueNonmem %||% NA))
    cat(sprintf("  without the 2 pi constant of the data (NONMEM OBJ): %.6f\n", x$valueNonmem))
  cat("\nParameters:\n")
  print(signif(x$population, digits))
  if (identical(x$seMethod, "scores")) cat("Standard errors from the empirical information of the group scores\n")
  if (identical(x$seMethod, "none")) cat("No standard errors; confint(method = \"profile\") gives intervals\n")
  if (isTRUE(x$partialCurvature))
    cat("Standard errors without the curvature of the prediction in the random effects within the volume term",
        "(forward over reverse; deriv2 = TRUE gives the full one)\n")
  if (!is.null(x$omega)) {
    cat("\nRandom effects:\n")
    om <- x$omega
    om$estimate <- signif(om$estimate, digits); om$se <- signif(om$se, digits)
    print(om, row.names = FALSE)
  }
  sh <- x$shrinkage
  # several groupings: one vector per grouping
  if (is.list(sh)) sh <- unlist(unname(Map(function(s, g) setNames(s, paste0(g, ": ", names(s))), sh, names(sh))))
  if (length(sh))
    cat("\nShrinkage:", paste(sprintf("%s = %.3f", names(sh), sh), collapse = ", "), "\n")
  invisible(x)
}


#' Log-likelihood of a mixed-effects fit
#'
#' @param object An `em` object, see [em()].
#' @param ... Ignored.
#' @return The marginal log-likelihood, `-value / 2`, of class `logLik` with
#'   the number of estimated parameters as `df` and the number of
#'   observations as `nobs`, so [stats::AIC()] and [stats::BIC()] apply.
#' @export
logLik.em <- function(object, ...) {
  structure(-object$value / 2, df = length(object$argument), nobs = object$nObs,
            class = "logLik")
}


#' Predictions of a mixed-effects fit
#'
#' @description Observed values against the population prediction (`PRED`,
#' random effects 0) and the individual prediction (`IPRED`, random effects at
#' their modes) per condition, with residuals.
#'
#' @param object An `em` object, see [em()].
#' @param times Additional times of the prediction grid.
#' @param type `"all"` (default) both predictions, `"individual"` only IPRED
#'   and its residuals, `"population"` only PRED and its residuals.
#' @param ... Ignored.
#' @return A data.frame with columns `condition, group, time, name, observed,
#'   sigma, sigmaPRED, IPRED, PRED, IRES, IWRES, PRES, PWRES, source`;
#'   `sigma` is the error model at IPRED, `sigmaPRED` at PRED, `source` is
#'   `"obs"` for the rows of the data and `"grid"` for the other times.
#' @export
predict.em <- function(object, times = NULL, type = c("all", "individual", "population"), ...) {
  type <- match.arg(type)
  design <- object$design
  data <- design$data
  conds <- names(data)
  group <- setNames(rep(NA_character_, length(conds)), conds)
  for (sd in .subDesigns(design)) for (g in sd$groups) {
    cs <- sd$conditions[[g]]
    group[cs] <- ifelse(is.na(group[cs]), g, paste(group[cs], g, sep = ":"))
  }
  grids <- lapply(conds, function(cn) sort(unique(c(.emGrid(design, cn), times))))
  eta <- .etaVector(object)
  predictAt <- function(pars) setNames(.predictMany(design$prdfn, times = grids,
                                                    parsList = rep(list(pars), length(conds)),
                                                    conditions = conds, fixed = object$fixed, deriv = FALSE), conds)
  ipred <- predictAt(c(object$argument, eta))
  pred <- predictAt(c(object$argument, eta * 0))
  sigmaOf <- function(pr, cn) if (!is.null(design$errfn))
    unclass(design$errfn(out = pr[[cn]], pars = getParameters(pr[[cn]]), conditions = cn)[[cn]])
  rows <- list()
  for (cn in conds) {
    d <- data[[cn]]
    pi <- ipred[[cn]]; pp <- pred[[cn]]
    si <- sigmaOf(ipred, cn); sp <- sigmaOf(pred, cn)
    for (nm in intersect(unique(as.character(d$name)), colnames(pi))) {
      dn <- d[d$name == nm, , drop = FALSE]
      idx <- match(dn$time, pi[, "time"])
      grid_idx <- setdiff(seq_len(nrow(pi)), idx)
      # a recorded sigma wins over the error model, as in normL2()
      sig <- function(sm, i, given) {
        model <- if (!is.null(sm) && nm %in% colnames(sm)) as.numeric(sm[i, nm]) else NA_real_
        ifelse(is.na(given), model, given)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        condition = cn, group = group[[cn]], time = dn$time, name = nm, observed = dn$value,
        sigma = sig(si, idx, dn$sigma), sigmaPRED = sig(sp, idx, dn$sigma),
        IPRED = pi[idx, nm], PRED = pp[idx, nm], source = "obs", stringsAsFactors = FALSE)
      if (length(grid_idx))
        rows[[length(rows) + 1L]] <- data.frame(
          condition = cn, group = group[[cn]], time = pi[grid_idx, "time"], name = nm,
          observed = NA_real_, sigma = sig(si, grid_idx, NA_real_),
          sigmaPRED = sig(sp, grid_idx, NA_real_), IPRED = pi[grid_idx, nm],
          PRED = pp[grid_idx, nm], source = "grid", stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  out$IRES <- out$observed - out$IPRED
  out$PRES <- out$observed - out$PRED
  out$IWRES <- out$IRES / out$sigma
  out$PWRES <- out$PRES / out$sigmaPRED
  rownames(out) <- NULL
  if (type == "individual") out <- out[, setdiff(names(out), c("PRED", "PRES", "PWRES", "sigmaPRED"))]
  if (type == "population") out <- out[, setdiff(names(out), c("IPRED", "IRES", "IWRES", "sigma"))]
  out
}


#' Residuals of a mixed-effects fit
#'
#' @param object An `em` object, see [em()].
#' @param type `"IWRES"` (default; data minus individual prediction over the
#'   error model), `"IRES"` (data minus individual prediction), `"CWRES"`
#'   (conditional weighted residuals of the FOCE linearisation at the modes:
#'   per group the residual of its data against \eqn{f(\eta^*) - G \eta^*},
#'   whitened by \eqn{G \Omega G' + \Sigma}, \eqn{G} the sensitivities in the
#'   normal random effects; without those equal to IWRES), `"PRES"` or
#'   `"PWRES"` (against the population prediction), `"NPDE"` (normalised
#'   prediction distribution errors from `nsim` simulated data sets,
#'   decorrelated per group, readings at or below their limit imputed as in
#'   Nguyen et al. 2012).
#' @param nsim,seed With `"NPDE"`: number of simulated data sets and the seed
#'   of the random number generator.
#' @param ... Ignored.
#' @return Numeric vector, one entry per data row, in the order of the rows of
#'   `predict(object)` with `source == "obs"`. NPDE have the attribute `tests`,
#'   the tests of mean 0 (t-test), variance 1 (chi-square) and normality
#'   (Shapiro-Wilk) with their p-values and the global p-value, three times the
#'   smallest.
#' @example inst/examples/NLME.R
#' @references Comets, E., Brendel, K. and Mentré, F. (2008). Computing
#'   normalised prediction distribution errors to evaluate nonlinear
#'   mixed-effect models. \emph{Computer Methods and Programs in Biomedicine}
#'   90(2), 154-166.
#'
#'   Nguyen, T. H. T., Comets, E. and Mentré, F. (2012). Extension of NPDE for
#'   evaluation of nonlinear mixed effect models in presence of data below the
#'   quantification limit with applications to HIV dynamic model.
#'   \emph{Journal of Pharmacokinetics and Pharmacodynamics} 39(5), 499-518.
#' @export
residuals.em <- function(object, type = c("IWRES", "IRES", "CWRES", "PRES", "PWRES", "NPDE"),
                         nsim = 1000L, seed = NULL, ...) {
  type <- match.arg(type)
  pr <- predict(object)
  pr <- pr[pr$source == "obs", , drop = FALSE]
  if (type == "CWRES") return(.cwres(object, pr))
  if (type != "NPDE") return(pr[[type]])
  .localSeed(seed)
  out <- .npde(object, nsim, pr)
  attr(out, "tests") <- .npdeTests(out)
  out
}

## Conditional weighted residuals per group, or per component with several groupings, over the
## rows of all its conditions, in the row order of `pr`; rows outside every group get IWRES.
.cwres <- function(object, pr) {
  design <- object$design
  out <- pr$IWRES
  if (is.null(design$groupings) && !design$normal$K) return(out)
  pars <- c(object$argument, .etaVector(object))
  Omegas <- if (is.null(design$groupings)) list(object$Omega) else object$Omega
  for (u in .volumeUnits(design, c(pars, object$fixed), object$fixed)) {
    eta <- u$eta
    Om <- matrix(0, length(eta), length(eta))
    at <- 0L
    for (pc in u$pieces) {
      k <- seq_along(pc$eta) + at
      Om[k, k] <- Omegas[[pc$k]]
      at <- at + length(pc$eta)
    }
    G <- list(); rows <- integer(0)
    for (cn in u$conditions) {
      r <- which(pr$condition == cn)
      p <- design$prdfn(times = .emGrid(design, cn), pars = pars, fixed = object$fixed,
                        conditions = cn)[[cn]]
      d <- attr(p, "deriv")
      ti <- .matchNum(pr$time[r], as.numeric(p[, "time"]))
      Gc <- matrix(0, length(r), length(eta))
      have <- intersect(eta, dimnames(d)[[3]])
      for (j in seq_along(r)) Gc[j, match(have, eta)] <- d[ti[j], pr$name[r[j]], have]
      G[[cn]] <- Gc; rows <- c(rows, r)
    }
    if (!length(rows)) next
    Gm <- do.call(rbind, G)
    res <- pr$observed[rows] - (pr$IPRED[rows] - drop(Gm %*% pars[eta]))
    C <- Gm %*% Om %*% t(Gm) + diag(pr$sigma[rows]^2, length(rows))
    out[rows] <- backsolve(chol(C), res, transpose = TRUE)
  }
  out
}


#' Confidence intervals of a mixed-effects fit
#'
#' @description `method = "wald"` from the standard errors of [vcov()][vcov.em];
#' `method = "profile"` refits with each parameter held at trial values and
#' finds where -2 log L has risen by the \eqn{\chi^2_1} quantile of `level`
#' on either side, by bisection from the Wald limits.
#'
#' @param object An `em` object, see [em()].
#' @param parm Parameter names; default all estimated parameters.
#' @param level Confidence level. Default 0.95.
#' @param method `"wald"` (default) or `"profile"`.
#' @param control With `method = "profile"`: `control` of [em()] for the
#'   refits, by default that of the fit.
#' @param ... Ignored.
#' @return A data.frame with columns `name`, `value`, `lower`, `upper`.
#' @export
confint.em <- function(object, parm = NULL, level = 0.95, method = c("wald", "profile"),
                       control = NULL, ...) {
  method <- match.arg(method)
  # the profile keeps the random effects of the fit
  control <- modifyList(control %||% object$control %||% list(), list(boundary = FALSE))
  est <- object$argument
  se <- .structuralSE(object)
  z <- stats::qnorm(1 - (1 - level) / 2)
  nm <- if (is.null(parm)) names(est) else intersect(parm, names(est))
  wald <- data.frame(name = nm, value = unname(est[nm]), lower = unname(est[nm] - z * se[nm]),
                     upper = unname(est[nm] + z * se[nm]), row.names = NULL)
  if (method == "wald") return(wald)
  crit <- object$value + stats::qchisq(level, 1)
  ofv <- function(p, v)
    .emFit(object$design, est[setdiff(names(est), p)], c(object$fixed, setNames(v, p)),
           object$method, control, FALSE)$value
  side <- function(p, dir) {
    s <- if (is.finite(se[p]) && se[p] > 0) se[p] else 0.1 * max(1, abs(est[p]))
    inner <- est[p]; outer <- est[p] + dir * z * s
    k <- 0L
    while (ofv(p, outer) < crit && k < 8L) { inner <- outer; outer <- est[p] + dir * 2^(k + 1) * z * s; k <- k + 1L }
    if (k == 8L) return(dir * Inf)
    for (i in 1:12) {
      mid <- (inner + outer) / 2
      if (ofv(p, mid) < crit) inner <- mid else outer <- mid
    }
    unname((inner + outer) / 2)
  }
  wald$lower <- vapply(nm, side, 0, dir = -1)
  wald$upper <- vapply(nm, side, 0, dir = 1)
  wald
}


#' Compare nested mixed-effects fits
#'
#' @description Likelihood ratio tests of fits ordered by their number of
#' parameters, each against the one before, with AIC and BIC.
#'
#' @param object,... `em` objects of nested models on the same data, or
#'   `alignment` objects.
#' @return A data.frame with `model` (the position of the fit among the
#'   arguments), `npar`, `value` (-2 log L), `AIC`, `BIC`, `chisq`, `df` and
#'   `p` per fit.
#' @order 1
#' @export
anova.em <- function(object, ...) .anovaTable(c(list(object), list(...)))

.anovaTable <- function(fits) {
  if (length(fits) < 2L) stop("anova: compare two or more fits.", call. = FALSE)
  npar <- vapply(fits, function(f) f$npar %||% length(f$argument), 0)
  value <- vapply(fits, `[[`, 0, "value")
  nobs <- vapply(fits, function(f) f$nObs %||% nrow(f$scaled), 0)
  ord <- order(npar)
  npar <- npar[ord]; value <- value[ord]; nobs <- nobs[ord]
  chisq <- c(NA, -diff(value)); df <- c(NA, diff(npar))
  data.frame(model = ord, npar = npar, value = value, AIC = value + 2 * npar,
             BIC = value + log(nobs) * npar, chisq = chisq, df = df,
             p = stats::pchisq(chisq, df, lower.tail = FALSE))
}


#' Simulate data from a mixed-effects fit
#'
#' @description Draws the random effects of every group from their fitted
#' distributions, predicts the data of every condition, adds the noise of the
#' error model (normal, or Student-t for [normT()]) and censors readings at
#' their lower limit.
#'
#' @param object An `em` object, see [em()].
#' @param nsim Number of data sets.
#' @param seed Seed of the random number generator, as in [stats::simulate()].
#' @param ... Ignored.
#' @return A list of `nsim` [datalist]s in the layout of the data of the fit.
#' @importFrom stats simulate
#' @export
simulate.em <- function(object, nsim = 1, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  design <- object$design
  pars <- c(object$argument, object$fixed)
  lapply(seq_len(nsim), function(i) {
    eta <- .drawEtas(design, pars)
    conds <- names(design$data)
    pred <- setNames(.predictMany(design$prdfn, times = lapply(conds, function(cn) .emGrid(design, cn)),
                                  parsList = rep(list(c(object$argument, eta)), length(conds)),
                                  conditions = conds, fixed = object$fixed, deriv = FALSE), conds)
    out <- lapply(names(design$data), function(cn) {
      d <- design$data[[cn]]
      pr <- pred[[cn]]
      sig <- d$sigma
      e <- if (!is.null(design$errfn)) design$errfn(out = pr, pars = getParameters(pr), conditions = cn)[[cn]]
      if (!is.null(e)) {
        model <- e[cbind(.matchNum(d$time, as.numeric(e[, "time"])), match(as.character(d$name), colnames(e)))]
        sig <- ifelse(is.na(sig), model, sig)
      }
      mu <- pr[cbind(.matchNum(d$time, as.numeric(pr[, "time"])), match(as.character(d$name), colnames(pr)))]
      noise <- if (is.finite(design$df)) stats::rt(nrow(d), design$df) else stats::rnorm(nrow(d))
      d$value <- mu + sig * noise
      if (!is.null(d$lloq)) d$value <- pmax(d$value, ifelse(is.na(d$lloq), -Inf, d$lloq))
      d
    })
    as.datalist(do.call(rbind, Map(function(d, cn) transform(d, condition = cn), out, names(design$data))))
  })
}

## One draw of the random effects of every group: normal blocks from Omega, centre Laplace and
## spike-and-slab and cluster coordinates from their hyperparameters, clusters by the Chinese
## restaurant with the population value seated first.
.drawEtas <- function(design, pars) {
  out <- setNames(numeric(length(.allEtas(design))), .allEtas(design))
  for (sd in .subDesigns(design)) {
    E <- sd$etaNames
    lay <- sd$normal
    if (lay$K) {
      L <- .choleskyFactor(lay, pars)
      Z <- matrix(stats::rnorm(nrow(E) * lay$K), nrow(E)) %*% t(L)
      out[as.vector(E[, lay$stems, drop = FALSE])] <- as.vector(Z)
    }
    for (b in Filter(function(b) b$dist != "normal", sd$blocks)) {
      value <- function(kind) if (is.na(b[[kind]])) pars[[b$hyper[[kind]]]] else b[[kind]]
      n <- nrow(E) * length(b$stems)
      if (b$dist == "dirichletProcess") {
        for (st in b$stems) {
          sizes <- 1L; values <- 0; table <- integer(nrow(E))
          for (g in seq_len(nrow(E))) {
            k <- sample.int(length(sizes) + 1L, 1L, prob = c(sizes, value("concentration")))
            if (k > length(sizes)) { sizes <- c(sizes, 0L); values <- c(values, stats::rnorm(1, 0, value("tau"))) }
            sizes[k] <- sizes[k] + 1L
            table[g] <- k
          }
          out[E[, st]] <- values[table]
        }
        next
      }
      out[as.vector(E[, b$stems, drop = FALSE])] <- if (b$dist == "laplace")
        stats::rexp(n, value("lambda")) * sample(c(-1, 1), n, TRUE)
      else stats::rbinom(n, 1, value("prob")) * stats::rnorm(n, 0, value("tau"))
    }
  }
  out
}
