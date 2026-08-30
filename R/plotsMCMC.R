## Plotting methods for mcmc / smc results.
## Split out of plots.R: these are methods for classes the core does not define.

#' Plot the marginal posterior distributions of an [mcmc()] result
#'
#' One panel per parameter: histogram of the post-warmup samples overlaid
#' with a kernel-density curve. The vertical dashed line marks the
#' posterior mean.
#'
#' @param x An object inheriting from `mcmcResult`.
#' @param parameters Optional character vector restricting which parameter
#'   panels are drawn.
#' @param bins Histogram bin count. Default 30.
#' @param ... Ignored.
#'
#' @return A ggplot.
#' @export
plot.mcmcresult <- function(x, parameters = NULL, bins = 30L, ...) {
  long <- .bayes_samples_long(x$samples)
  if (!is.null(parameters)) {
    bad <- setdiff(parameters, levels(long$parameter))
    if (length(bad))
      stop("Unknown parameter name(s): ", paste(bad, collapse = ", "))
    long <- long[long$parameter %in% parameters, , drop = FALSE]
    long$parameter <- droplevels(long$parameter)
  }
  means <- stats::aggregate(value ~ parameter, data = long, FUN = mean)
  ggplot2::ggplot(long, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = bins, fill = "grey80", color = "grey40") +
    ggplot2::geom_density(color = dMod_colors[3], linewidth = 0.7) +
    ggplot2::geom_vline(data = means,
                        ggplot2::aes(xintercept = value),
                        color = dMod_colors[2], linetype = "dashed",
                        linewidth = 0.6) +
    ggplot2::facet_wrap(~ parameter, scales = "free") +
    ggplot2::labs(x = NULL, y = "density",
                  title = "Marginal posterior") +
    theme_dMod(base_size = 11)
}

#' Multi-chain density overlay with R-hat annotation
#'
#' @param x A `mcmcResultMulti` object.
#' @param parameters Optional character vector restricting parameters.
#' @param bins Histogram bin count (currently unused).
#' @param ... Ignored.
#' @return A ggplot.
#' @export
plot.mcmcresultmulti <- function(x, parameters = NULL, bins = 30L, ...) {
  S <- x$samples
  ch <- factor(x$chainId)
  if (!is.null(parameters)) S <- S[, parameters, drop = FALSE]
  par_names <- colnames(S)
  long <- do.call(rbind, lapply(seq_along(par_names), function(j) {
    data.frame(parameter = par_names[j], chain = ch, value = S[, j],
               stringsAsFactors = FALSE)
  }))
  long$parameter <- factor(long$parameter, levels = par_names)

  finite_rhat <- if (!is.null(x$rHat)) x$rHat[is.finite(x$rHat)] else numeric()
  rhat_label <- if (length(finite_rhat))
    sprintf("Marginal posterior (%d chains, max R-hat = %.3f)",
            x$nChains, max(finite_rhat))
  else
    sprintf("Marginal posterior (%d chains)", x$nChains)

  ggplot2::ggplot(long, ggplot2::aes(x = value,
                                       color = .data$chain,
                                       fill  = .data$chain)) +
    ggplot2::geom_density(alpha = 0.15, linewidth = 0.6) +
    ggplot2::facet_wrap(~ parameter, scales = "free") +
    scale_color_dMod() +
    ggplot2::scale_fill_manual(values = dMod_palette(nlevels(ch))) +
    ggplot2::labs(x = NULL, y = "density", title = rhat_label) +
    theme_dMod(base_size = 11)
}

#' @rdname plotTrace
#' @description Method for [mcmc()] sequential / SMC outputs: shows the
#'   adaptive tempering schedule (beta vs level), the post-resample ESS,
#'   and per-level acceptance.
#' @export
plotTrace.mcmcresultsequential <- function(x, ...) {
  diag_df <- rbind(
    data.frame(level = seq_along(x$betaPath) - 1L, value = x$betaPath,
               panel = "beta"),
    data.frame(level = seq_along(x$ESSPath) - 1L,  value = x$ESSPath,
               panel = "ESS"))
  if (length(x$acceptRates))
    diag_df <- rbind(diag_df,
                     data.frame(level = seq_along(x$acceptRates),
                                value = x$acceptRates,
                                panel = "accept rate"))
  ggplot2::ggplot(diag_df, ggplot2::aes(x = level, y = value)) +
    ggplot2::geom_line(color = dMod_colors[3], linewidth = 0.6) +
    ggplot2::geom_point(color = dMod_colors[1], size = 1.2) +
    ggplot2::facet_wrap(~ panel, scales = "free_y", ncol = 1L) +
    ggplot2::labs(x = "SMC level",
                  title = sprintf("SMC trace (log-evidence = %s)",
                                  format(x$logEvidence, digits = 4))) +
    theme_dMod(base_size = 11)
}

#' @rdname plotTrace
#' @description Method for [mcmc()] single-chain outputs: per-parameter
#'   trace of the chain across iterations.
#' @param parameters Optional character vector restricting which parameter
#'   panels are drawn.
#' @export
plotTrace.mcmcresultsingle <- function(x, parameters = NULL, ...) {
  S <- x$samples
  if (!is.null(parameters)) S <- S[, parameters, drop = FALSE]
  long <- .bayes_samples_long(S)
  long$iter <- long$sample
  ggplot2::ggplot(long, ggplot2::aes(x = iter, y = value)) +
    ggplot2::geom_line(color = dMod_colors[3], linewidth = 0.4) +
    ggplot2::facet_wrap(~ parameter, scales = "free_y") +
    ggplot2::labs(x = "iteration", y = NULL, title = "Chain trace") +
    theme_dMod(base_size = 11)
}

#' @rdname plotTrace
#' @export
plotTrace.mcmcresultblocked <- function(x, parameters = NULL, ...) {
  plotTrace.mcmcresultsingle(x, parameters = parameters, ...)
}

#' @rdname plotPairs
#' @param parameters Character, the parameters to show. `NULL` shows all.
#' @export
plotPairs.mcmcresult <- function(x, parameters = NULL, ...) {
  S <- x$samples
  if (!is.null(parameters)) S <- S[, parameters, drop = FALSE]
  par_names <- colnames(S)
  P <- ncol(S)
  if (P < 2L) stop("plotPairs needs >= 2 parameters.")

  combos <- expand.grid(xPar = par_names, yPar = par_names,
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  combos$panel <- with(combos, ifelse(xPar == yPar, "diag",
                                       ifelse(match(yPar, par_names) >
                                              match(xPar, par_names),
                                              "lower", "upper")))
  combos <- combos[combos$panel %in% c("diag", "lower"), , drop = FALSE]
  long <- do.call(rbind, lapply(seq_len(nrow(combos)), function(idx) {
    xp <- combos$xPar[idx]; yp <- combos$yPar[idx]
    data.frame(xPar  = factor(xp, levels = par_names),
               yPar  = factor(yp, levels = par_names),
               x     = S[, xp],
               y     = if (combos$panel[idx] == "diag") NA_real_ else S[, yp],
               panel = combos$panel[idx], stringsAsFactors = FALSE)
  }))

  p <- ggplot2::ggplot()
  diag_df <- long[long$panel == "diag",  , drop = FALSE]
  off_df  <- long[long$panel == "lower", , drop = FALSE]
  if (nrow(off_df) > 0L)
    p <- p + ggplot2::geom_point(data = off_df,
                                  ggplot2::aes(x = x, y = y),
                                  alpha = 0.3, size = 0.5,
                                  color = dMod_colors[3])
  if (nrow(diag_df) > 0L)
    p <- p + ggplot2::geom_density(data = diag_df,
                                    ggplot2::aes(x = x),
                                    color = dMod_colors[2],
                                    linewidth = 0.7,
                                    inherit.aes = FALSE)
  p +
    ggplot2::facet_grid(yPar ~ xPar, scales = "free", switch = "both") +
    ggplot2::labs(x = NULL, y = NULL, title = "Posterior pair-plot") +
    theme_dMod(base_size = 10) +
    ggplot2::theme(strip.placement = "outside",
                   panel.spacing   = ggplot2::unit(0.1, "lines"))
}
