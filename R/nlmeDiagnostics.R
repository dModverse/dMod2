## Diagnostics of mixed-effects fits from simulations: NPDE, visual predictive checks, random
## effects against covariates.

utils::globalVariables(c("lo", "mid", "hi", "percentile", "blq", "blqLo", "blqMid", "blqHi", "value"))

## Per condition the data rows in the row order of predict() with source "obs": per observable of
## the data the rows of that observable.
.emObsRows <- function(design) {
  lapply(design$data, function(d) {
    nm <- as.character(d$name)
    unlist(lapply(unique(nm), function(n) which(nm == n)), use.names = FALSE)
  })
}

## Simulated readings of the data rows, uncensored, rows in the order of predict(object) with
## source "obs", one column per simulation: per column the random effects of every group drawn
## from their distribution and the noise of the error model, as simulate.em().
.simulatedRows <- function(object, nsim) {
  design <- object$design
  conds <- names(design$data)
  grids <- lapply(conds, function(cn) .emGrid(design, cn))
  pars <- c(object$argument, object$fixed)
  rowsOf <- .emObsRows(design)
  data <- Map(function(d, o) d[o, , drop = FALSE], design$data, rowsOf)
  node <- if (!is.null(design$errfn)) .fnNode(design$errfn)
  one <- function() {
    eta <- .drawEtas(design, pars, object$probPosterior)
    preds <- setNames(.predictMany(design$prdfn, times = grids, parsList = rep(list(c(object$argument, eta)), length(conds)),
                                   conditions = conds, fixed = object$fixed, deriv = FALSE), conds)
    errs <- if (!is.null(node)) {
      inner <- lapply(conds, function(cn) getParameters(preds[[cn]]))
      fixedInner <- lapply(inner, function(p) as.parvec(p[attr(p, "fixed")], deriv = FALSE, deriv2 = FALSE))
      inner <- lapply(inner, function(p) as.parvec(p[.setdiff(names(p), names(object$fixed))], deriv = FALSE,
                                                   deriv2 = FALSE))
      .evalNode(node, .bundle(conds = conds, out = unname(preds), pars = inner, fixed = fixedInner),
                FALSE, FALSE, NULL, 1L)
    } else lapply(conds, function(cn) .errorFrame(design, preds[[cn]], cn))
    unlist(lapply(seq_along(conds), function(i) {
      d <- data[[i]]; pr <- preds[[i]]; e <- errs[[i]]
      nm <- as.character(d$name)
      mu <- pr[cbind(.matchNum(d$time, as.numeric(pr[, "time"])), match(nm, colnames(pr)))]
      sig <- d$sigma
      if (!is.null(e) && anyNA(sig))
        sig <- ifelse(is.na(sig), e[cbind(.matchNum(d$time, as.numeric(e[, "time"])), match(nm, colnames(e)))], sig)
      noise <- if (is.finite(design$df)) stats::rt(nrow(d), design$df) else stats::rnorm(nrow(d))
      mu + sig * noise
    }), use.names = FALSE)
  }
  matrix(unlist(lapply(seq_len(nsim), function(k) one())), ncol = nsim)
}

## The data rows in the order of predict() with their lower limit, -Inf for none.
.emObsLimits <- function(design) {
  rowsOf <- .emObsRows(design)
  unlist(Map(function(d, o) {
    l <- d$lloq %||% rep(-Inf, nrow(d))
    l[is.na(l)] <- -Inf
    l[o]
  }, design$data, rowsOf), use.names = FALSE)
}

## Normalized prediction distribution errors per data row (Comets et al. 2008), per unit of groups.
.npde <- function(object, nsim, pr) {
  design <- object$design
  sims <- .simulatedRows(object, nsim)
  y <- pr$observed
  lloq <- .emObsLimits(design)
  # the conditions of a group, or of a component of several groupings, share the random effects
  units <- if (is.null(design$groupings)) unname(design$conditions[design$groups])
           else lapply(design$components, `[[`, "conditions")
  units <- c(units, as.list(.setdiff(names(design$data), unlist(units))))
  out <- rep(NA_real_, length(y))
  for (u in units) {
    r <- which(pr$condition %in% u)
    if (!length(r)) next
    z <- .npdeUnit(y[r], sims[r, , drop = FALSE], lloq[r])
    if (is.null(z))
      stop("residuals: the simulated readings of ", paste(u, collapse = ", "), " have a singular covariance; ",
           "raise `nsim`.", call. = FALSE)
    out[r] <- z
  }
  out
}

## NPDE of the readings `y` of one unit against their simulations `S`, rows the readings: decorrelated
## by the Cholesky factor of the simulated covariance, the pde the mid-rank among the simulations,
## censored readings imputed by .imputeCensored(). NULL for a singular covariance.
.npdeUnit <- function(y, S, lloq) {
  L <- tryCatch(t(chol(stats::cov(t(S)))), error = function(e) NULL)
  if (is.null(L)) return(NULL)
  E <- rowMeans(S)
  ys <- forwardsolve(L, .imputeCensored(y, S, lloq) - E)
  below <- rowSums(forwardsolve(L, S - E) < drop(ys))
  stats::qnorm((below + 0.5) / (ncol(S) + 1))
}

## Readings at or below their limit drawn from their law given the other readings of the unit and
## the limits: elliptical slice sampling on the Gaussian copula of the simulations, mapped back by
## their quantiles; without correlation uniform below the censored fraction (Nguyen et al. 2012).
.imputeCensored <- function(y, S, lloq, iter = 50L) {
  cens <- which(y <= lloq)
  if (!length(cens)) return(y)
  K <- ncol(S)
  obs <- .setdiff(seq_along(y), cens)
  scores <- stats::qnorm((t(apply(S, 1, rank)) - 0.5) / K)
  Q <- solve(stats::cor(t(scores)))
  w <- stats::qnorm((rowSums(S < y) + 0.5) / (K + 1))
  up <- stats::qnorm(pmax(rowMeans(S <= lloq), 1 / (2 * K)))[cens]
  # the scores of the censored readings given the others: normal with mean m and covariance C
  C <- solve(Q[cens, cens, drop = FALSE])
  m <- -drop(C %*% Q[cens, obs, drop = FALSE] %*% w[obs])
  Lc <- t(chol(C))
  x <- stats::qnorm(stats::runif(length(cens)) * stats::pnorm(up))
  for (it in seq_len(iter)) {
    nu <- drop(Lc %*% stats::rnorm(length(cens)))
    th <- stats::runif(1, 0, 2 * pi)
    lo <- th - 2 * pi
    hi <- th
    repeat {
      prop <- m + (x - m) * cos(th) + nu * sin(th)
      if (all(prop <= up)) break
      if (th < 0) lo <- th else hi <- th
      th <- stats::runif(1, lo, hi)
    }
    x <- prop
  }
  y[cens] <- vapply(seq_along(cens), function(i)
    stats::quantile(S[cens[i], ], stats::pnorm(x[i]), names = FALSE, type = 1), 0)
  y
}

## Tests of NPDE against N(0, 1): the mean by a t-test, the variance by the chi-square test of
## variance 1, the distribution by Shapiro-Wilk, and the global p-value, the smallest times 3.
.npdeTests <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  chi <- (n - 1) * stats::var(x)
  pv <- 2 * min(stats::pchisq(chi, n - 1), stats::pchisq(chi, n - 1, lower.tail = FALSE))
  sw <- stats::shapiro.test(if (n > 5000) x[seq_len(5000)] else x)
  tt <- stats::t.test(x)
  p <- c(tt$p.value, pv, sw$p.value)
  data.frame(test = c("mean", "variance", "normality", "global"),
             statistic = c(unname(tt$statistic), chi / (n - 1), unname(sw$statistic), NA),
             p = c(p, min(1, 3 * min(p))), stringsAsFactors = FALSE)
}


#' Visual predictive check
#'
#' @description Per observable and bin of time the 5th, 50th and 95th
#' percentiles of the readings against the 95 percent prediction intervals of
#' the same percentiles of `nsim` data sets simulated from the fit
#' ([simulate()][simulate.em]). Censored readings and simulations enter at
#' their limit, and a percentile at or below the limit of the bin is not
#' drawn. Where the data have a limit of quantification a second panel shows
#' the fraction of censored readings per bin against the prediction interval
#' of that fraction.
#'
#' @param x An `em` fit, see [em()].
#' @param nsim Number of simulated data sets. Default 200.
#' @param bins Number of bins, at quantiles of the times of the readings, or
#'   the breaks of the bins. Default 8.
#' @param pc `TRUE`: prediction corrected (Bergstrand et al. 2011), every
#'   reading and simulation times the median population prediction of the bin
#'   over the population prediction of that reading. Without `logScale`, readings whose
#'   population prediction is not positive are left out, with their simulations.
#' @param logScale With `pc`, `TRUE` for observables on a log scale: the
#'   correction adds the median population prediction of the bin minus the
#'   own one instead.
#' @param strata Optional column of the condition grid of the data of the fit
#'   ([datalist]), one panel per value.
#' @param seed Seed of the random number generator. Default `NULL`.
#' @param ... Ignored.
#' @return A ggplot, with the BLQ panel a list of two ggplots where cowplot is
#'   not available; attribute `vpc` holds the table of the percentiles and
#'   intervals per bin.
#' @references Bergstrand, M., Hooker, A. C., Wallin, J. E. and Karlsson, M. O.
#'   (2011). Prediction-corrected visual predictive checks for diagnosing
#'   nonlinear mixed-effects models. \emph{The AAPS Journal} 13(2), 143-151.
#' @example inst/examples/plotVPC.R
#' @export
plotVPC <- function(x, ...) UseMethod("plotVPC", x)

#' @rdname plotVPC
#' @export
plotVPC.em <- function(x, nsim = 200L, bins = 8L, pc = FALSE, logScale = FALSE, strata = NULL, seed = NULL, ...) {
  .localSeed(seed)
  design <- x$design
  pr <- predict(x)
  pr <- pr[pr$source == "obs", , drop = FALSE]
  sims <- .simulatedRows(x, nsim)
  lloq <- .emObsLimits(design)
  if (pc && !logScale) {
    # the correction is relative to the population prediction, undefined where that is not positive
    keep <- is.finite(pr$PRED) & pr$PRED > 0
    if (!any(keep))
      stop("plotVPC: no reading has a positive population prediction; use `pc = FALSE`.", call. = FALSE)
    pr <- pr[keep, , drop = FALSE]
    sims <- sims[keep, , drop = FALSE]
    lloq <- lloq[keep]
  }
  censored <- pr$observed <= lloq
  simCensored <- sims <= lloq
  y <- pmax(pr$observed, lloq)
  sims <- pmax(sims, lloq)
  times <- pr$time
  breaks <- if (length(bins) > 1L) sort(bins)
            else unique(stats::quantile(times, seq(0, 1, length.out = bins + 1L), names = FALSE))
  pr$bin <- findInterval(times, breaks, rightmost.closed = TRUE, all.inside = TRUE)
  pr$stratum <- if (is.null(strata)) "all" else {
    grid <- attr(design$data, "condition.grid")
    if (is.null(grid) || !strata %in% names(grid))
      stop("plotVPC: `strata` is a column of the condition grid of the data.", call. = FALSE)
    s <- grid[match(pr$condition, rownames(grid)), strata]
    factor(s, if (is.numeric(s)) as.character(sort(unique(s))) else .naturalLevels(s))
  }
  key <- paste(pr$name, pr$stratum, pr$bin, sep = "\r")
  if (pc) {
    medPred <- stats::ave(pr$PRED, key, FUN = stats::median)
    shift <- if (logScale) function(v) v + (medPred - pr$PRED) else function(v) v * (medPred / pr$PRED)
    y <- shift(y)
    sims <- shift(sims)
  }
  probs <- c(0.05, 0.5, 0.95)
  rows <- split(seq_len(nrow(pr)), factor(key, unique(key)))
  table <- do.call(rbind, lapply(rows, function(r) {
    simQ <- apply(sims[r, , drop = FALSE], 2, stats::quantile, probs = probs, names = FALSE)
    simQ <- matrix(simQ, length(probs))
    limit <- if (pc) -Inf else max(lloq[r])
    blqSim <- colMeans(simCensored[r, , drop = FALSE])
    data.frame(name = pr$name[r[1]], stratum = pr$stratum[r[1]], bin = pr$bin[r[1]],
               time = stats::median(times[r]), percentile = paste0(100 * probs, "%"),
               observed = stats::quantile(y[r], probs, names = FALSE),
               lo = apply(simQ, 1, stats::quantile, 0.025, names = FALSE),
               mid = apply(simQ, 1, stats::median),
               hi = apply(simQ, 1, stats::quantile, 0.975, names = FALSE),
               limit = limit, limited = any(is.finite(lloq[r])), blq = mean(censored[r]),
               blqLo = stats::quantile(blqSim, 0.025, names = FALSE), blqMid = stats::median(blqSim),
               blqHi = stats::quantile(blqSim, 0.975, names = FALSE), stringsAsFactors = FALSE)
  }))
  rownames(table) <- NULL
  # the observables in the order of the data
  table$name <- factor(table$name, unique(as.character(pr$name)))
  shown <- table
  shown$observed[shown$observed <= shown$limit] <- NA
  facet <- if (length(unique(table$stratum)) > 1L) ggplot2::facet_grid(name ~ stratum, scales = "free_y")
           else ggplot2::facet_wrap(~ name, scales = "free_y")
  p <- ggplot2::ggplot(shown, ggplot2::aes(x = time, group = percentile)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi, fill = percentile == "50%"), alpha = 0.3) +
    ggplot2::geom_line(ggplot2::aes(y = mid), linetype = "dashed", color = "gray40") +
    ggplot2::geom_line(ggplot2::aes(y = observed), color = dMod2_colors[2], linewidth = 0.8, na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(y = observed), color = dMod2_colors[2], na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = dMod2_colors[1], `FALSE` = dMod2_colors[3]), guide = "none") +
    facet + ggplot2::labs(x = "Time", y = if (pc) "Prediction-corrected value" else "Value",
                          title = sprintf("%sVPC, %d simulations: 5, 50, 95%% of the data and their 95%% intervals",
                                          if (pc) "Prediction-corrected " else "", nsim)) +
    theme_dMod2(base_size = 11)
  attr(p, "vpc") <- table
  if (!any(is.finite(lloq))) return(p)
  bt <- table[table$percentile == "50%" & table$limited, , drop = FALSE]
  blq <- ggplot2::ggplot(bt, ggplot2::aes(x = time)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = blqLo, ymax = blqHi), fill = dMod2_colors[3], alpha = 0.3) +
    ggplot2::geom_line(ggplot2::aes(y = blqMid), linetype = "dashed", color = "gray40") +
    ggplot2::geom_line(ggplot2::aes(y = blq), color = dMod2_colors[2], linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(y = blq), color = dMod2_colors[2]) +
    facet + ggplot2::labs(x = "Time", y = "Fraction censored") + theme_dMod2(base_size = 11)
  out <- if (requireNamespace("cowplot", quietly = TRUE)) cowplot::plot_grid(p, blq, ncol = 1, rel_heights = c(2, 1))
         else list(vpc = p, blq = blq)
  attr(out, "vpc") <- table
  out
}


## The random effects of every group against covariates: per random effect and covariate a
## scatter with a loess trend for a numeric covariate of more than five values, boxes otherwise.
.randomEffectCovariates <- function(fit, covariates) {
  design <- fit$design
  grid <- if (is.data.frame(covariates)) covariates else {
    g <- attr(design$data, "condition.grid")
    miss <- .setdiff(covariates, names(g))
    if (length(miss))
      stop("plotRandomEffects: ", paste(miss, collapse = ", "), " not in the condition grid of the data; ",
           "give `covariates` as a data.frame with one row per condition.", call. = FALSE)
    g[, covariates, drop = FALSE]
  }
  if ("condition" %in% names(grid) && !all(rownames(grid) %in% names(design$data)))
    rownames(grid) <- grid$condition
  cols <- .setdiff(names(grid), "condition")
  random <- fit$random
  if (is.null(random)) stop("plotRandomEffects: the fit has no random effects.", call. = FALSE)
  subs <- .subDesigns(design)
  conditionOf <- unlist(lapply(subs, function(sd) vapply(sd$conditions, `[`, "", 1)))
  first <- conditionOf[random$group]
  long <- do.call(rbind, lapply(cols, function(cv) {
    v <- grid[first, cv]
    # a covariate of few values is a category, in the order of the values
    levels <- if (is.numeric(v)) as.character(sort(unique(v))) else .naturalLevels(v)
    data.frame(eta = factor(random$eta, unique(random$eta)), mode = random$mode, covariate = cv,
               numeric = is.numeric(v) && length(unique(v)) > 5L, value = factor(as.character(v), levels),
               x = suppressWarnings(as.numeric(as.character(v))), stringsAsFactors = FALSE)
  }))
  plots <- lapply(split(long, factor(long$covariate, cols)), function(d) {
    p <- ggplot2::ggplot(d, ggplot2::aes(y = mode)) + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")
    p <- if (d$numeric[1]) p + ggplot2::geom_point(ggplot2::aes(x = x), alpha = 0.7) +
      ggplot2::geom_smooth(ggplot2::aes(x = x), method = "loess", formula = y ~ x, se = FALSE, color = dMod2_colors[2],
                           linewidth = 0.7)
    else p + ggplot2::geom_boxplot(ggplot2::aes(x = value), outlier.shape = NA) +
      ggplot2::geom_jitter(ggplot2::aes(x = value), width = 0.15, height = 0, alpha = 0.5)
    p + ggplot2::facet_wrap(~ eta, scales = "free_y", nrow = 1) +
      ggplot2::labs(x = d$covariate[1], y = "mode of the random effect") + theme_dMod2(base_size = 11)
  })
  if (length(plots) == 1L) return(plots[[1]])
  if (requireNamespace("cowplot", quietly = TRUE)) cowplot::plot_grid(plotlist = plots, ncol = 1) else plots
}
