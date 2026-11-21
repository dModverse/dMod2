## NLME plotting and prediction methods for `em` fits.

utils::globalVariables(c("IPRED", "PRED", "eta", "part", "xmin", "xmax", "ymin", "ymax", "inclusion", "label",
                         "group", "estimate", "lower", "upper", "fit", "position", "width"))

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
#' @example inst/examples/emPredictions.R
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
  long$condition <- factor(long$condition, .naturalLevels(long$condition))
  long$name <- factor(long$name, unique(as.character(pf$name)))
  multiple <- length(unique(long$name)) > 1L
  p <- ggplot2::ggplot(long, ggplot2::aes(x = predicted, y = observed,
                                          color = condition)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = "gray50") +
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
    scale_color_dMod2() +
    ggplot2::labs(x = "Prediction", y = "Observed",
                  title = sprintf("Observed vs predicted (method = %s)", fit$method)) +
    theme_dMod2(base_size = 11) +
    ggplot2::theme(legend.position = "none")
}

#' Residual diagnostics of a mixed-effects fit
#'
#' @description One row of panels per type of residual: against the population
#'   prediction (IWRES against the individual one) and against time, both with
#'   a loess smoother, and a QQ-plot against N(0, 1). See
#'   [residuals()][residuals.em].
#' @param parframe An `em` fit, see [em()].
#' @param type Any of `"IWRES"`, `"CWRES"` and `"NPDE"`, by default all three.
#' @param nsim,seed With `"NPDE"`: the number of simulated data sets, default
#'   1000, and the seed of the random number generator, default `NULL`.
#' @param ... Ignored.
#' @return A ggplot.
#' @example inst/examples/plotResiduals.em.R
#' @export
plotResiduals.em <- function(parframe, type = c("IWRES", "CWRES", "NPDE"), nsim = 1000L, seed = NULL, ...) {
  fit <- parframe
  type <- match.arg(type, several.ok = TRUE)
  pf <- predict(fit)
  pf <- pf[pf$source == "obs", , drop = FALSE]
  long <- do.call(rbind, lapply(type, function(ty) {
    r <- if (ty == "IWRES") pf$IWRES else stats::residuals(fit, type = ty, nsim = nsim, seed = seed)
    vs <- if (ty == "IWRES") "IPRED" else "PRED"
    q <- stats::qnorm(stats::ppoints(length(r)))[rank(r, ties.method = "first")]
    rbind(data.frame(type = ty, panel = paste("vs", vs), x = pf[[vs]], y = r, name = pf$name, stringsAsFactors = FALSE),
          data.frame(type = ty, panel = "vs time", x = pf$time, y = r, name = pf$name, stringsAsFactors = FALSE),
          data.frame(type = ty, panel = "QQ", x = q, y = r, name = pf$name, stringsAsFactors = FALSE))
  }))
  long$type <- factor(long$type, type)
  long$panel <- factor(long$panel, c("vs IPRED", "vs PRED", "vs time", "QQ"))
  long$name <- factor(long$name, unique(as.character(pf$name)))
  scatter <- long[long$panel != "QQ", , drop = FALSE]
  qq <- long[long$panel == "QQ", , drop = FALSE]
  ggplot2::ggplot(long, ggplot2::aes(x = x, y = y, color = name)) +
    ggplot2::geom_hline(data = scatter, ggplot2::aes(yintercept = 0), linetype = "dashed", color = "gray50") +
    ggplot2::geom_abline(data = qq, ggplot2::aes(slope = 1, intercept = 0), linetype = "dashed", color = "gray50") +
    ggplot2::geom_point(alpha = 0.7, size = 1.2) +
    ggplot2::geom_smooth(data = scatter, ggplot2::aes(group = 1), se = FALSE, method = "loess", formula = y ~ x,
                         color = dMod2_colors[2], linewidth = 0.6) +
    ggplot2::facet_grid(type ~ panel, scales = "free") +
    scale_color_dMod2(guide = if (length(unique(long$name)) > 1L) "legend" else "none") +
    ggplot2::labs(x = NULL, y = "residual", color = NULL, title = "Residual diagnostics") +
    theme_dMod2(base_size = 11)
}

#' Data and predictions of a mixed-effects fit
#'
#' @description One panel per condition: the data as points, the individual
#'   prediction (IPRED) with a band of one standard deviation of the error
#'   model, the population prediction (PRED) dashed. Several observables give
#'   one row each.
#' @param prediction An `em` fit, see [em()].
#' @param times Times added to the prediction grid of [predict()][predict.em],
#'   by default 200 over the range of the data. The curves are drawn over the
#'   range of the data only.
#' @param ncol Number of columns with one observable. Default 4.
#' @param showPred Logical; draw PRED. Default `TRUE`.
#' @param showBand Logical; draw the band of the error model. Default `TRUE`.
#' @param conditionsPerPage Optional integer: the conditions are split into
#'   pages of at most this many, one plot per page.
#' @param ... Ignored.
#' @return A ggplot, or a list of ggplots when `conditionsPerPage` is set.
#' @example inst/examples/emPredictions.R
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
  conditions <- .naturalLevels(if (is.factor(pf$condition)) droplevels(pf$condition) else pf$condition)
  pf$condition <- factor(pf$condition, conditions)
  pf$name <- factor(pf$name, unique(as.character(pf$name)))

  page <- function(pp, label = NULL) {
    obs <- pp[pp$source == "obs", , drop = FALSE]
    grd <- pp[order(pp$condition, pp$name, pp$time), , drop = FALSE]
    p <- ggplot2::ggplot()
    if (showBand && any(is.finite(grd$sigma)))
      p <- p + ggplot2::geom_ribbon(data = grd, ggplot2::aes(x = time, ymin = IPRED - sigma, ymax = IPRED + sigma),
                                    fill = dMod2_colors[3], alpha = 0.2, linetype = 0)
    p <- p + ggplot2::geom_line(data = grd, ggplot2::aes(x = time, y = IPRED, color = "IPRED"), linewidth = 0.7)
    if (showPred)
      p <- p + ggplot2::geom_line(data = grd, ggplot2::aes(x = time, y = PRED, color = "PRED"),
                                  linewidth = 0.5, linetype = "dashed")
    p <- p + ggplot2::geom_point(data = obs, ggplot2::aes(x = time, y = observed), size = 1.6, alpha = 0.85)
    p <- p + if (multiple) ggplot2::facet_grid(name ~ condition, scales = "free_y")
             else ggplot2::facet_wrap(~ condition, ncol = ncol, scales = "free_y")
    title <- paste(sprintf("Data and predictions (method = %s)", fit$method), label)
    p + ggplot2::scale_color_manual(values = c(IPRED = dMod2_colors[3], PRED = dMod2_colors[2])) +
      ggplot2::labs(x = "Time", y = "Value", color = NULL, title = trimws(title)) +
      theme_dMod2(base_size = 11)
  }

  if (is.null(conditionsPerPage)) return(page(pf))
  n <- suppressWarnings(as.integer(conditionsPerPage))
  if (length(n) != 1L || is.na(n) || n < 1L)
    stop("plotCombined: `conditionsPerPage` is a positive integer or NULL.", call. = FALSE)
  pages <- split(conditions, ceiling(seq_along(conditions) / n))
  lapply(seq_along(pages), function(i) {
    sub <- pf[as.character(pf$condition) %in% pages[[i]], , drop = FALSE]
    sub$condition <- factor(sub$condition, levels = pages[[i]])
    page(sub, sprintf("(page %d/%d)", i, length(pages)))
  })
}

#' @rdname plotRandomEffects
#' @param x An `em` fit, see [em()].
#' @param type `"conditional"` (default): per group `draws` values from the
#'   conditional distribution of the random effects of the group; under the
#'   model they follow the random-effect distribution whatever the shrinkage.
#'   A normal random effect is drawn normal around the mode with the
#'   conditional standard error, against `N(0, Omega_kk)`. A spike-and-slab
#'   random effect is 0 with the posterior probability `p0` and otherwise
#'   normal, with the reported value, or where the reported structure has the
#'   spike the mode of the slab, as mean and the variance that gives the
#'   conditional standard error; drawn against the prior, the slab density
#'   `prob N(0, tau^2)` and a cap at the height of a bar of area `1 - prob`
#'   at 0. `"modes"`: the
#'   modes, or the reported values of spike-and-slab random effects, narrower
#'   than the distribution by the shrinkage. `"inclusion"`: a heatmap of the
#'   posterior probability `1 - p0` of a deviation per group and sparse random
#'   effect, labeled with the reported value where the reported structure has
#'   a deviation.
#' @param draws Number of draws per group with `type = "conditional"`. Default 20.
#' @param bins Number of bars of every histogram. Default 30.
#' @param covariates Columns of the condition grid of the data of the fit
#'   ([datalist]), or a data.frame with one row per condition (row names or a
#'   column `condition`): instead the modes of the random effects of every
#'   group against the covariates of the first condition of the group, a scatter with a
#'   loess trend for a numeric covariate of more than five values, boxes for
#'   the others.
#' @param ... Ignored.
#' @return A ggplot, or without cowplot a list with `hist` and `qq`; the QQ
#'   plot has the normal random effects only, a fit with sparse random effects
#'   alone gives the histograms as one ggplot. With `type = "inclusion"` a
#'   ggplot. With `covariates`, a ggplot for one covariate; for several a
#'   cowplot grid, or without cowplot a list of ggplots named by covariate.
#'   Spike-and-slab random effects that reach no data of the group are left
#'   out; Dirichlet-process random effects, which couple the groups, appear in
#'   the heatmap only.
#' @example inst/examples/plotRandomEffects.R
#' @export
plotRandomEffects.em <- function(x, type = c("conditional", "modes", "inclusion"), draws = 20L, bins = 30,
                                 covariates = NULL, ...) {
  if (!is.null(covariates)) return(.randomEffectCovariates(x, covariates))
  type <- match.arg(type)
  if (type == "inclusion") return(.inclusionHeatmap(x))
  values <- rbind(.randomEffectValues(x, type, draws), .spikeSlabValues(x, type, draws))
  if (is.null(values))
    stop("plotRandomEffects: the fit has no normal or spike-and-slab random effects.", call. = FALSE)
  # the panels in the order of the random effects of the fit
  values$eta <- factor(values$eta, unique(values$eta))
  bars <- do.call(rbind, lapply(split(values, values$eta), .mixtureBars, bins = bins))
  prior <- do.call(rbind, lapply(split(values, values$eta), .mixturePrior, bars = bars))
  what <- if (type == "conditional") "Conditional draws" else "Modes"
  sparse <- any(values$spike > 0)
  hist <- ggplot2::ggplot(bars) +
    ggplot2::geom_rect(ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = part),
                       color = "gray40", linewidth = 0.2) +
    ggplot2::geom_line(data = prior[prior$part == "slab", ], ggplot2::aes(x = x, y = density), inherit.aes = FALSE,
                       color = dMod2_colors[2], linewidth = 0.8) +
    ggplot2::geom_segment(data = prior[prior$part == "spike", ],
                          ggplot2::aes(x = x - 1.5 * width, xend = x + 1.5 * width, y = density, yend = density),
                          inherit.aes = FALSE, color = dMod2_colors[2], linewidth = 0.8) +
    ggplot2::facet_wrap(~ eta, scales = "free") +
    ggplot2::scale_fill_manual(values = c(slab = "gray80", spike = dMod2_colors[3]),
                               labels = c(slab = "deviation", spike = "exactly 0"),
                               guide = if (sparse) "legend" else "none") +
    ggplot2::labs(x = "eta value", y = "density", fill = NULL,
                  title = paste(what, if (sparse) "vs the prior (red), the prior spike of area 1 - prob as a cap"
                                      else "vs N(0, Omega_kk)")) +
    theme_dMod2(base_size = 11)
  normal <- values[values$spike == 0, , drop = FALSE]
  if (!nrow(normal)) return(hist)
  qq <- ggplot2::ggplot(normal, ggplot2::aes(sample = value / sd)) +
    ggplot2::stat_qq() +
    ggplot2::stat_qq_line(color = dMod2_colors[2]) +
    ggplot2::facet_wrap(~ eta) +
    ggplot2::labs(x = "Theoretical quantile (N(0,1))", y = "Standardized eta", title = paste(what, "QQ vs N(0,1)")) +
    theme_dMod2(base_size = 11)
  if (requireNamespace("cowplot", quietly = TRUE)) cowplot::plot_grid(hist, qq, ncol = 1)
  else list(hist = hist, qq = qq)
}

## The values of every normal random effect, the modes or `draws` draws per group around the mode
## with the conditional standard error, with the standard deviation of Omega; NULL for a fit
## without normal random effects.
.randomEffectValues <- function(fit, type, draws) {
  several <- is.list(fit$Omega)
  part <- function(f) if (several) fit[[f]] else list(fit[[f]])
  Omegas <- part("Omega"); etaModes <- part("etaModes"); etaSE <- part("etaSE"); shrinkage <- part("shrinkage")
  keep <- !vapply(Omegas, is.null, TRUE) & !vapply(etaModes, is.null, TRUE)
  if (!any(keep)) return(NULL)
  do.call(rbind, unlist(lapply(which(keep), function(g) {
    E <- etaModes[[g]]; O <- Omegas[[g]]; S <- etaSE[[g]]
    if (type == "conditional" && is.null(S))
      stop("plotRandomEffects: the fit has no conditional standard errors of the random effects; ",
           "use type = \"modes\".", call. = FALSE)
    lapply(intersect(colnames(E), rownames(O)), function(k) {
      value <- if (type == "modes") unname(E[, k]) else
        rep(unname(E[, k]), each = draws) + rep(unname(S[, k]), each = draws) * stats::rnorm(nrow(E) * draws)
      sh <- if (is.null(shrinkage[[g]])) NA_real_ else shrinkage[[g]][k]
      eta <- paste0(if (several) paste0(names(Omegas)[g], ": "), k,
                    if (is.finite(sh)) sprintf(", shrinkage %.0f%%", 100 * sh))
      data.frame(eta = eta, value = value, sd = sqrt(O[k, k]), spike = 0, stringsAsFactors = FALSE)
    })
  }), recursive = FALSE))
}

## The prior of every spike-and-slab random effect of a fit, by the short name of the effect:
## probability `prob` of a deviation and slab standard deviation `tau`.
.spikeSlabPriors <- function(fit) {
  sl <- fit$slab
  if (is.null(sl)) return(NULL)
  sl <- sl[sl$distribution == "spikeSlab", , drop = FALSE]
  kind <- sub("_.*$", "", sl$name)
  rows <- lapply(unique(sl$symbols), function(s) {
    b <- sl[sl$symbols == s, , drop = FALSE]
    k <- kind[sl$symbols == s]
    data.frame(eta = strsplit(s, ", ", fixed = TRUE)[[1]], prob = b$estimate[k == "prob"],
               tau = b$estimate[k == "tau"], stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

## The values of every spike-and-slab random effect that reaches the data of the group, reported or
## `draws` per group from the conditional mixture of 0 and a normal with the conditional variance;
## NULL for a fit without spike-and-slab random effects.
.spikeSlabValues <- function(fit, type, draws) {
  pri <- .spikeSlabPriors(fit)
  r <- fit$random
  if (is.null(pri) || is.null(r$p0)) return(NULL)
  keep <- r$eta %in% pri$eta
  if (!is.null(r$data)) keep <- keep & r$data
  r <- r[keep, , drop = FALSE]
  if (!nrow(r)) return(NULL)
  pri <- pri[match(r$eta, pri$eta), ]
  stem <- colnames(fit$etaModes)[match(r$eta, sub("^eta_", "", colnames(fit$etaModes)))]
  inner <- if (is.null(fit$innerModes)) r$mode else fit$innerModes[cbind(r$group, stem)]
  # the slab mean: the reported value where the structure has the deviation, else the slab mode
  m <- ifelse(r$cluster != 0, r$mode, inner)
  p0 <- pmin(pmax(r$p0, 0), 1)
  s <- sqrt(pmax((r$se^2 - p0 * (1 - p0) * m^2) / pmax(1 - p0, .Machine$double.eps), 0))
  n <- if (type == "modes") 1L else as.integer(draws)
  value <- if (type == "modes") r$mode else {
    slab <- stats::runif(nrow(r) * n) >= rep(p0, each = n)
    ifelse(slab, rep(m, each = n) + rep(s, each = n) * stats::rnorm(nrow(r) * n), 0)
  }
  data.frame(eta = rep(r$eta, each = n), value = value, sd = rep(pri$tau, each = n),
             spike = rep(1 - pri$prob, each = n), stringsAsFactors = FALSE)
}

## The histogram of the values of one panel as density bars, `bins` over the range of the values
## and centered on multiples of the bin width so that one bin is centered at 0; in a spike-and-slab panel the exact zeros as a bar of
## their own at the bottom of that bin.
.mixtureBars <- function(v, bins) {
  w <- diff(range(v$value)) / bins
  if (!(w > 0)) w <- 7 * v$sd[1] / bins
  at <- round(v$value / w)
  zero <- v$spike[1] > 0 & v$value == 0
  height <- function(sel) {
    tab <- table(at[sel])
    setNames(as.vector(tab) / (nrow(v) * w), names(tab))
  }
  spike <- height(zero); slab <- height(!zero)
  base <- ifelse(names(slab) %in% names(spike), spike[names(slab)], 0)
  bar <- function(part, mid, lo, hi) if (length(mid))
    data.frame(eta = v$eta[1], part = part, xmin = (mid - 0.5) * w, xmax = (mid + 0.5) * w, ymin = lo, ymax = hi,
               width = w, stringsAsFactors = FALSE)
  rbind(bar("slab", as.numeric(names(slab)), unname(base), unname(base + slab)),
        bar("spike", as.numeric(names(spike)), 0, unname(spike)))
}

## The density of the prior of one panel over at least +/- 3.5 standard deviations: normal, or for a
## spike-and-slab panel the slab density prob N(0, tau^2) and the height of a spike bar of area
## 1 - prob in the width of the bins of the panel.
.mixturePrior <- function(v, bars) {
  lim <- max(abs(v$value), 3.5 * v$sd[1])
  x <- seq(-lim, lim, length.out = 200L)
  spike <- v$spike[1]
  slab <- data.frame(eta = v$eta[1], part = "slab", x = x, density = (1 - spike) * stats::dnorm(x, 0, v$sd[1]),
                     width = NA_real_, stringsAsFactors = FALSE)
  if (spike == 0) return(slab)
  w <- bars$width[bars$eta == v$eta[1]][1]
  rbind(slab, data.frame(eta = v$eta[1], part = "spike", x = 0, density = spike / w, width = w,
                         stringsAsFactors = FALSE))
}

## The posterior probability of a deviation, 1 - p0, per group and sparse random effect as tiles,
## effects from top to bottom in the order of the fit, with the reported value where the reported
## structure has a deviation.
.inclusionHeatmap <- function(fit) {
  r <- fit$random
  if (is.null(r$p0) || all(is.na(r$p0)))
    stop("plotRandomEffects: type = \"inclusion\" needs spike-and-slab or Dirichlet-process random effects.",
         call. = FALSE)
  r <- r[!is.na(r$p0), , drop = FALSE]
  if (!is.null(r$grouping)) r$group <- paste0(r$grouping, ": ", r$group)
  r$eta <- factor(r$eta, rev(unique(r$eta)))
  r$group <- factor(r$group, unique(r$group))
  r$inclusion <- 1 - r$p0
  r$label <- ifelse(r$cluster != 0, sprintf("%.2f", r$mode), "")
  if (!is.null(r$data)) r$label[!r$data] <- "no data"
  ggplot2::ggplot(r, ggplot2::aes(x = group, y = eta, fill = inclusion)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = label, color = inclusion > 0.6), size = 3) +
    ggplot2::scale_fill_gradient(low = "white", high = dMod2_colors[3], limits = c(0, 1)) +
    ggplot2::scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
    ggplot2::labs(x = NULL, y = NULL, fill = "1 - p0", title = "Posterior probability of a deviation") +
    theme_dMod2(base_size = 11)
}

#' @rdname plotValues
#' @export
plotValues.em <- function(x, ..., tol = 1, showSteps = FALSE, pseudoLog = FALSE, relative = FALSE)
  .plotWaterfall(list(.valueFrame(x, ...)), tol, showSteps, pseudoLog, relative)

#' @rdname plotTrace
#' @param x An `em` fit, see [em()], with a stage trace: every method except `"focei"`
#'   with normal random effects only. The hyperparameters of spike-and-slab and
#'   Dirichlet-process random effects that the trace records (`prob_<grouping>`,
#'   `tau_<grouping>`, `concentration_<grouping>`) get one panel each.
#' @param ... Ignored.
#' @return A ggplot.
#' @example inst/examples/plotTrace.R
#' @export
plotTrace.em <- function(x, ...) {
  tr <- x$stageTrace
  if (is.null(tr))
    stop("plotTrace needs a fit with a stage trace; method \"focei\" with normal random effects has none.",
         call. = FALSE)
  panels <- c(OFV = "OFV", deltaPsi = "|delta psi|", deltaPar = "|delta par|", deltaStruct = "|delta theta|",
              decrement = "Newton decrement", maxSoftmax = "max softmax", nEffMin = "n_eff min")
  hyper <- intersect(x$slab$name, names(tr))
  panels <- c(panels, setNames(hyper, hyper))
  have <- intersect(names(panels), names(tr))
  # the stages in the order they ran
  stage <- tr$level %||% tr$phase %||% rep("all", nrow(tr))
  stage <- factor(stage, unique(stage))
  long <- do.call(rbind, lapply(have, function(v)
    data.frame(iter = seq_len(nrow(tr)), value = tr[[v]], panel = panels[[v]], stage = stage)))
  long$panel <- factor(long$panel, unname(panels[have]))
  # an iteration without a value in a panel, e.g. the step of the first one, has no point there
  long <- long[is.finite(long$value), ]
  ggplot2::ggplot(long, ggplot2::aes(x = iter, y = value, color = stage)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::facet_wrap(~ panel, scales = "free_y") +
    scale_color_dMod2() +
    ggplot2::labs(x = "EM iteration", title = "EM convergence trace") +
    theme_dMod2(base_size = 11)
}


#' @rdname plotPars
#' @param level With an `em` fit or a list of them: the level of the Wald
#'   intervals, see [confint()][confint.em]. Default 0.95.
#' @param which With an `em` fit or a list of them: names of the parameters to
#'   show, by default all estimated ones.
#' @details An `em` fit, see [em()], or a list of them, as fits of one model by
#'   several methods: the estimates of fixed effects and hyperparameters with
#'   Wald intervals, where the fit has standard errors, one row per parameter
#'   and one color per fit. The fits of a list are labeled by the names of the
#'   list, else by their method.
#' @export
plotPars.em <- function(x, ..., level = 0.95, which = NULL) .plotEstimates(list(x), level, which)

#' @rdname plotPars
#' @export
plotPars.list <- function(x, ..., level = 0.95, which = NULL) {
  if (!length(x) || !all(vapply(x, inherits, TRUE, "em")))
    stop("plotPars: a list of fits is a list of em fits.", call. = FALSE)
  .plotEstimates(x, level, which)
}

## The estimates of the fits `fits` with Wald intervals at `level`, the rows of the parameters in
## the order of their first appearance, the fits side by side within a row.
.plotEstimates <- function(fits, level, which) {
  labels <- names(fits) %||% rep("", length(fits))
  labels[!nzchar(labels)] <- vapply(fits[!nzchar(labels)], `[[`, "", "method")
  labels <- make.unique(labels, sep = " ")
  est <- do.call(rbind, Map(function(f, lab) {
    ci <- stats::confint(f, level = level)
    data.frame(fit = lab, parameter = ci$name, estimate = ci$value, lower = ci$lower, upper = ci$upper,
               stringsAsFactors = FALSE)
  }, fits, labels))
  if (!is.null(which)) est <- est[est$parameter %in% which, , drop = FALSE]
  if (!nrow(est)) stop("plotPars: none of `which` is an estimated parameter of the fits.", call. = FALSE)
  pars <- unique(est$parameter)
  est$fit <- factor(est$fit, labels)
  # the fits side by side within the row of a parameter, the first on top
  offset <- if (length(labels) > 1L) 0.6 * (0.5 - (as.integer(est$fit) - 1) / (length(labels) - 1)) else 0
  est$position <- length(pars) + 1 - match(est$parameter, pars) + offset
  ranged <- est[is.finite(est$lower) & is.finite(est$upper), , drop = FALSE]
  ggplot2::ggplot(est, ggplot2::aes(x = estimate, y = position, color = fit)) +
    ggplot2::geom_segment(data = ranged, ggplot2::aes(x = lower, xend = upper, yend = position)) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::scale_y_continuous(breaks = rev(seq_along(pars)), labels = pars, minor_breaks = NULL) +
    scale_color_dMod2(limits = labels, guide = if (length(labels) > 1L) "legend" else "none") +
    ggplot2::labs(x = "Estimate", y = NULL, color = NULL,
                  title = sprintf("Estimates with %g%% Wald intervals", 100 * level)) +
    theme_dMod2(base_size = 11)
}
