## NLME plotting and prediction methods for `em` fits.

utils::globalVariables(c("IPRED", "PRED", "eta"))

#' Random-effect distribution diagnostics
#'
#' @description Per random effect a histogram against the estimated
#'   `N(0, Omega_kk)` density and a QQ-plot against the estimated normal, of draws
#'   from the conditional distributions of the random effects or of their modes;
#'   with several groupings one panel per grouping and random effect.
#' @param x Object to plot.
#' @param ... Method-specific arguments.
#' @return A ggplot (or a list of ggplots if `cowplot` is unavailable).
#' @export
plotRandomEffects <- function(x, ...) UseMethod("plotRandomEffects", x)

#' Observed vs predicted scatter (DV vs IPRED and DV vs PRED)
#'
#' @description S3 plot method for an `em` fit, see [em()]. Two-panel scatter: DV vs IPRED
#'   on the left, DV vs PRED on the right, identity line shown dashed. With a
#'   multi-observable fit the observable becomes a faceting row
#'   (`name ~ panel`); otherwise one row with `facet_wrap(~ panel)`.
#' @param x An `em` fit, see [em()].
#' @param ... Ignored.
#' @return A ggplot.
#' @export
plot.em <- function(x, ...) {
  fit <- x
  pf <- predict(fit)
  pf <- pf[pf$source == "obs", , drop = FALSE]
  long <- rbind(
    data.frame(condition = pf$condition, name = pf$name,
               observed = pf$observed, predicted = pf$IPRED,
               panel = "IPRED", stringsAsFactors = FALSE),
    data.frame(condition = pf$condition, name = pf$name,
               observed = pf$observed, predicted = pf$PRED,
               panel = "PRED", stringsAsFactors = FALSE))
  multiple <- length(unique(long$name)) > 1L
  p <- ggplot2::ggplot(long, ggplot2::aes(x = predicted, y = observed,
                                          color = condition)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = "grey50") +
    ggplot2::geom_point(alpha = 0.75)
  if (multiple) {
    # Observable scales can differ by orders of magnitude, so the scales are
    # free and coord_equal is not used.
    p <- p + ggplot2::facet_grid(name ~ panel, scales = "free")
  } else {
    rng <- range(c(long$observed, long$predicted), na.rm = TRUE)
    p <- p +
      ggplot2::coord_equal(xlim = rng, ylim = rng) +
      ggplot2::facet_wrap(~ panel)
  }
  p +
    scale_color_dMod() +
    ggplot2::labs(x = "Prediction", y = "Observed",
                  title = sprintf("Observed vs predicted (method = %s)", fit$method)) +
    theme_dMod(base_size = 11) +
    ggplot2::theme(legend.position = "none")
}

#' Weighted-residual diagnostics for an EM
#'
#' @description Two-panel scatter: IWRES vs IPRED and IWRES vs TIME with a
#'   loess smoother.
#' @param parframe An `em` fit, see [em()].
#' @param ... Ignored.
#' @return A ggplot.
#' @export
plotResiduals.em <- function(parframe, ...) {
  fit <- parframe
  pf <- predict(fit)
  pf <- pf[pf$source == "obs", , drop = FALSE]
  long <- rbind(
    data.frame(x = pf$IPRED, y = pf$IWRES, name = pf$name,
               panel = "IWRES vs IPRED", stringsAsFactors = FALSE),
    data.frame(x = pf$time,  y = pf$IWRES, name = pf$name,
               panel = "IWRES vs TIME",  stringsAsFactors = FALSE))
  p <- ggplot2::ggplot(long, ggplot2::aes(x = x, y = y, color = name)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::geom_smooth(se = FALSE, method = "loess", formula = y ~ x,
                         color = dMod_colors[2], linewidth = 0.6)
  if (length(unique(long$name)) > 1L) {
    p <- p + ggplot2::facet_grid(name ~ panel, scales = "free_x")
  } else {
    p <- p + ggplot2::facet_wrap(~ panel, scales = "free_x")
  }
  p +
    scale_color_dMod() +
    ggplot2::labs(x = NULL, y = "IWRES",
                  title = "Weighted residual diagnostics") +
    theme_dMod(base_size = 11)
}

#' Data and predictions of a mixed-effects fit
#'
#' @description One panel per condition: the data as points, the individual
#'   prediction (IPRED) with a band of one standard deviation of the error
#'   model, the population prediction (PRED) dashed. Several observables give
#'   one row each.
#' @param prediction An `em` fit, see [em()].
#' @param times Times of the predicted curves, by default 200 over the range
#'   of the data.
#' @param ncol Number of columns with one observable. Default 4.
#' @param showPred Logical; draw PRED. Default `TRUE`.
#' @param showBand Logical; draw the band of the error model. Default `TRUE`.
#' @param conditionsPerPage Optional integer: the conditions are split into
#'   pages of at most this many, one plot per page.
#' @param ... Ignored.
#' @return A ggplot, or a list of ggplots when `conditionsPerPage` is set.
#' @export
plotCombined.em <- function(prediction, times = NULL, ncol = 4L, showPred = TRUE, showBand = TRUE,
                            conditionsPerPage = NULL, ...) {
  fit <- prediction
  obsTimes <- sort(unique(unlist(lapply(fit$data, `[[`, "time"))))
  if (is.null(times)) times <- seq(min(obsTimes), max(obsTimes), length.out = 200L)
  pf <- predict(fit, times = times)
  # the curves over the range of the data, the data rows all
  pf <- pf[pf$source == "obs" | (pf$time >= min(obsTimes) & pf$time <= max(obsTimes)), , drop = FALSE]
  multiple <- length(unique(pf$name)) > 1L
  conditions <- if (is.factor(pf$condition)) levels(droplevels(pf$condition)) else unique(as.character(pf$condition))

  page <- function(pp, label = NULL) {
    obs <- pp[pp$source == "obs", , drop = FALSE]
    grd <- pp[order(pp$condition, pp$name, pp$time), , drop = FALSE]
    p <- ggplot2::ggplot()
    if (showBand && any(is.finite(grd$sigma)))
      p <- p + ggplot2::geom_ribbon(data = grd, ggplot2::aes(x = time, ymin = IPRED - sigma, ymax = IPRED + sigma),
                                    fill = dMod_colors[3], alpha = 0.2, linetype = 0)
    p <- p + ggplot2::geom_line(data = grd, ggplot2::aes(x = time, y = IPRED, color = "IPRED"), linewidth = 0.7)
    if (showPred)
      p <- p + ggplot2::geom_line(data = grd, ggplot2::aes(x = time, y = PRED, color = "PRED"),
                                  linewidth = 0.5, linetype = "dashed")
    p <- p + ggplot2::geom_point(data = obs, ggplot2::aes(x = time, y = observed), size = 1.6, alpha = 0.85)
    p <- p + if (multiple) ggplot2::facet_grid(name ~ condition, scales = "free_y")
             else ggplot2::facet_wrap(~ condition, ncol = ncol, scales = "free_y")
    title <- paste(sprintf("Data and predictions (method = %s)", fit$method), label)
    p + ggplot2::scale_color_manual(values = c(IPRED = dMod_colors[3], PRED = dMod_colors[2])) +
      ggplot2::labs(x = "Time", y = "Value", color = NULL, title = trimws(title)) +
      theme_dMod(base_size = 11)
  }

  if (is.null(conditionsPerPage)) return(page(pf))
  n <- suppressWarnings(as.integer(conditionsPerPage))
  if (length(n) != 1L || is.na(n) || n < 1L)
    stop("plotCombined: `conditionsPerPage` is a positive integer or NULL.", call. = FALSE)
  pages <- split(conditions, ceiling(seq_along(conditions) / n))
  lapply(seq_along(pages), function(i) {
    sub <- pf[as.character(pf$condition) %in% pages[[i]], , drop = FALSE]
    if (is.factor(sub$condition)) sub$condition <- factor(sub$condition, levels = pages[[i]])
    page(sub, sprintf("(page %d/%d)", i, length(pages)))
  })
}

#' @rdname plotRandomEffects
#' @param x An `em` fit, see [em()].
#' @param type `"conditional"` (default): per group `draws` values from the
#'   conditional distribution of its random effects, normal around the mode with
#'   the conditional standard error; under the model they follow `N(0, Omega_kk)`
#'   whatever the shrinkage. `"modes"`: the modes, narrower than `N(0, Omega_kk)`
#'   by the shrinkage.
#' @param draws Number of draws per group with `type = "conditional"`. Default 20.
#' @param bins Number of bars of every histogram. Default 30.
#' @param ... Ignored.
#' @return A ggplot (or a list with `hist` and `qq` if cowplot is unavailable).
#' @export
plotRandomEffects.em <- function(x, type = c("conditional", "modes"), draws = 20L, bins = 30, ...) {
  type <- match.arg(type)
  values <- .randomEffectValues(x, type, draws)
  # the density of every random effect over at least +/- 3.5 standard deviations
  density <- do.call(rbind, lapply(split(values, values$eta), function(v) {
    lim <- max(abs(v$value), 3.5 * v$sd[1])
    x <- seq(-lim, lim, length.out = 200L)
    data.frame(eta = v$eta[1], x = x, density = stats::dnorm(x, 0, v$sd[1]), stringsAsFactors = FALSE)
  }))
  what <- if (type == "conditional") "Conditional draws" else "Modes"
  hist <- ggplot2::ggplot(values, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = bins, fill = "grey80",
                            color = "grey40") +
    ggplot2::geom_line(data = density, ggplot2::aes(x = x, y = density), inherit.aes = FALSE,
                       color = dMod_colors[2], linewidth = 0.8) +
    ggplot2::facet_wrap(~ eta, scales = "free") +
    ggplot2::labs(x = "eta value", y = "density", title = paste(what, "vs N(0, Omega_kk)")) +
    theme_dMod(base_size = 11)
  qq <- ggplot2::ggplot(values, ggplot2::aes(sample = value / sd)) +
    ggplot2::stat_qq() +
    ggplot2::stat_qq_line(color = dMod_colors[2]) +
    ggplot2::facet_wrap(~ eta) +
    ggplot2::labs(x = "Theoretical quantile (N(0,1))", y = "Standardised eta", title = paste(what, "QQ vs N(0,1)")) +
    theme_dMod(base_size = 11)
  if (requireNamespace("cowplot", quietly = TRUE)) cowplot::plot_grid(hist, qq, ncol = 1)
  else list(hist = hist, qq = qq)
}

## The values of every normal random effect, the modes or `draws` draws per group from the normal
## around the mode with the conditional standard error, with the standard deviation of Omega; the
## panel names the grouping where there are several and the shrinkage.
.randomEffectValues <- function(fit, type, draws) {
  several <- is.list(fit$Omega)
  part <- function(f) if (several) fit[[f]] else list(fit[[f]])
  Omegas <- part("Omega"); etaModes <- part("etaModes"); etaSE <- part("etaSE"); shrinkage <- part("shrinkage")
  keep <- !vapply(Omegas, is.null, TRUE) & !vapply(etaModes, is.null, TRUE)
  if (!any(keep))
    stop("plotRandomEffects: the fit has no normal random effects.", call. = FALSE)
  do.call(rbind, unlist(lapply(which(keep), function(g) {
    E <- etaModes[[g]]; O <- Omegas[[g]]; S <- etaSE[[g]]
    if (type == "conditional" && is.null(S))
      stop("plotRandomEffects: the fit has no conditional standard errors of its random effects; ",
           "use type = \"modes\".", call. = FALSE)
    lapply(intersect(colnames(E), rownames(O)), function(k) {
      value <- if (type == "modes") unname(E[, k]) else
        rep(unname(E[, k]), each = draws) + rep(unname(S[, k]), each = draws) * stats::rnorm(nrow(E) * draws)
      sh <- if (is.null(shrinkage[[g]])) NA_real_ else shrinkage[[g]][k]
      eta <- paste0(if (several) paste0(names(Omegas)[g], ": "), k,
                    if (is.finite(sh)) sprintf(", shrinkage %.0f%%", 100 * sh))
      data.frame(eta = eta, value = value, sd = sqrt(O[k, k]), stringsAsFactors = FALSE)
    })
  }), recursive = FALSE))
}

#' @rdname plotValues
#' @export
plotValues.em <- function(x, ..., tol = 1, showSteps = FALSE, pseudoLog = FALSE) {
  if (is.null(x$multistart))
    stop("plotValues: the fit has no multistart table; msem() returns one.", call. = FALSE)
  starts <- x$multistart[!is.na(x$multistart$value), , drop = FALSE]
  plotValues(parframe(starts, parameters = character(0)), ..., tol = tol, showSteps = showSteps,
             pseudoLog = pseudoLog)
}

#' @rdname plotTrace
#' @param x An `em` fit, see [em()], with a stage trace: every method except `"focei"`
#'   with normal random effects only.
#' @param ... Ignored.
#' @return A ggplot.
#' @export
plotTrace.em <- function(x, ...) {
  tr <- x$stageTrace
  if (is.null(tr))
    stop("plotTrace needs a fit with a stage trace; method \"focei\" with normal random effects has none.",
         call. = FALSE)
  panels <- c(OFV = "OFV", deltaPsi = "|delta psi|", deltaPar = "|delta par|", deltaStruct = "|delta theta|",
              decrement = "Newton decrement", maxSoftmax = "max softmax", nEffMin = "n_eff min")
  have <- intersect(names(panels), names(tr))
  stage <- factor(tr$level %||% tr$phase %||% rep("all", nrow(tr)))
  long <- do.call(rbind, lapply(have, function(v)
    data.frame(iter = seq_len(nrow(tr)), value = tr[[v]], panel = panels[[v]], stage = stage)))
  ggplot2::ggplot(long, ggplot2::aes(x = iter, y = value, color = stage)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::facet_wrap(~ panel, scales = "free_y") +
    scale_color_dMod() +
    ggplot2::labs(x = "EM iteration", title = "EM convergence trace") +
    theme_dMod(base_size = 11)
}
