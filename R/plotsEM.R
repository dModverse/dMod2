## NLME plotting and prediction methods for `em` fits.
## Split out of plots.R: these are methods for a class the core does not define.

#' Predictions from an EM object
#'
#' @description
#' Returns a long-format `data.frame` of observed values vs population (`PRED`,
#' eta = 0) and individual (`IPRED`, eta at posterior modes) predictions per
#' condition, plus residuals. Used directly by the diagnostic plot helpers.
#'
#' @param object An [EM] (from [EM]).
#' @param times Optional numeric vector of additional times for the smooth
#'   IPRED/PRED curves. Defaults to the union of observed times.
#' @param ... Ignored.
#'
#' @return A long-format `data.frame` with columns
#'   `condition, time, name, observed, sigma, IPRED, PRED, IRES, IWRES, PRES,
#'   PWRES, source` where `source = "obs"` for rows aligned with data
#'   observations and `source = "grid"` for the dense IPRED/PRED smooth grid
#'   (observed/sigma NA for grid rows).
#' @export
predict.em <- function(object, times = NULL, ...) {
  .emRequireOmega(object, "predict")
  if (is.null(object$prdfn) || is.null(object$data) || is.null(object$omega))
    stop("predict..fitNormal: fit is missing `prdfn`, `data`, or `omega` ",
         "(was the fit built by EM()?).")

  prdfn <- object$prdfn
  data  <- object$data
  om    <- object$omega
  etaModes <- object$etaModes
  pars      <- object$argument

  subjects <- rownames(om$subjectEtas)
  N <- length(subjects)
  obs_times <- sort(unique(unlist(lapply(data, `[[`, "time"))))
  grid_times <- if (is.null(times)) sort(unique(c(0, obs_times)))
                else sort(unique(c(0, obs_times, times)))

  eta_full <- as.vector(etaModes)
  names(eta_full) <- as.vector(om$subjectEtas)
  pars_ipred <- c(pars, eta_full)
  pars_pred  <- c(pars, setNames(rep(0, length(eta_full)), names(eta_full)))

  pred_ipred <- prdfn(times = grid_times, pars = pars_ipred, conditions = subjects)
  pred_pred  <- prdfn(times = grid_times, pars = pars_pred,  conditions = subjects)

  # Restrict to observables that actually appear in the data. Without this
  # filter, internal model states emitted by the prdfn (e.g. the depot Ag in a
  # PK model, present because the observation Y() is built with
  # attach.input = TRUE) would be plotted as ghost curves on the IPRED/PRED
  # axis next to the real observable, swamping the y-axis with the dose.
  data_names <- unique(unlist(lapply(data, `[[`, "name")))
  obs_names <- intersect(setdiff(colnames(pred_ipred[[1]]), "time"),
                         data_names)

  # If an errfn is attached and the data table has NA sigmas, evaluate
  # the errfn per-condition on the prediction so IWRES diagnostics get a
  # meaningful weight. Mirrors evalConditionResidual's contract.
  err <- object$errfn
  sigma_ipred <- NULL
  if (!is.null(err)) {
    sigma_ipred <- setNames(lapply(subjects, function(s) {
      pinner <- getParameters(pred_ipred[[s]])
      err(out = pred_ipred[[s]], pars = pinner, conditions = s)[[s]]
    }), subjects)
  }

  rows <- list()
  for (s in subjects) {
    d_s <- data[[s]]
    pi  <- pred_ipred[[s]]
    pp  <- pred_pred[[s]]
    si  <- if (!is.null(sigma_ipred)) sigma_ipred[[s]] else NULL
    # Coerce prdframe class away so plain matrix indexing works.
    if (!is.null(si)) si <- unclass(si)
    for (nm in obs_names) {
      d_nm <- d_s[d_s$name == nm, , drop = FALSE]
      get_sigma <- function(row_idx) {
        if (!is.null(si) && nm %in% colnames(si) && all(row_idx <= nrow(si)))
          as.numeric(si[row_idx, nm])
        else rep(NA_real_, length(row_idx))
      }
      if (nrow(d_nm)) {
        idx <- match(d_nm$time, pi[, "time"])
        sig_obs <- if (!is.null(si)) get_sigma(idx)
                   else if (all(is.na(d_nm$sigma))) rep(NA_real_, nrow(d_nm))
                   else d_nm$sigma
        rows[[length(rows) + 1L]] <- data.frame(
          condition = s, time = d_nm$time, name = nm,
          observed  = d_nm$value, sigma = sig_obs,
          IPRED     = pi[idx, nm], PRED  = pp[idx, nm],
          source    = "obs", stringsAsFactors = FALSE
        )
      }
      grid_idx <- setdiff(seq_len(nrow(pi)), match(d_nm$time, pi[, "time"]))
      if (length(grid_idx)) {
        rows[[length(rows) + 1L]] <- data.frame(
          condition = s, time = pi[grid_idx, "time"], name = nm,
          observed  = NA_real_, sigma = get_sigma(grid_idx),
          IPRED     = pi[grid_idx, nm], PRED = pp[grid_idx, nm],
          source    = "grid", stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- do.call(rbind, rows)
  out$IRES  <- out$observed - out$IPRED
  out$PRES  <- out$observed - out$PRED
  out$IWRES <- out$IRES / out$sigma
  out$PWRES <- out$PRES / out$sigma
  rownames(out) <- NULL
  out
}

#' Observed vs predicted scatter (DV vs IPRED and DV vs PRED)
#'
#' @description S3 plot method for [EM]. Two-panel scatter: DV vs IPRED
#'   on the left, DV vs PRED on the right, identity line shown dashed. With a
#'   multi-observable fit the observable becomes a faceting row
#'   (`name ~ panel`); otherwise one row with `facet_wrap(~ panel)`.
#' @param x An [EM].
#' @param ... Ignored.
#' @return A ggplot.
#' @export
plot.em <- function(x, ...) {
  .emRequireOmega(x, "plot")
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
  multi_obs <- length(unique(long$name)) > 1L
  p <- ggplot2::ggplot(long, ggplot2::aes(x = predicted, y = observed,
                                          color = condition)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = "grey50") +
    ggplot2::geom_point(alpha = 0.75)
  if (multi_obs) {
    # Per-observable scales differ by orders of magnitude (e.g. log-cp vs
    # linear PCA); free scales are required and coord_equal is incompatible
    # with that.
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
#' @param parframe An [EM].
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

#' Per-subject individual fits for an EM
#'
#' @description Per-subject IPRED curve with an IPRED +/- sigma ribbon derived
#'   from the fit's attached errfn, optional population PRED overlay dashed,
#'   and observed values as points. When the fit has a single observable the
#'   layout is `facet_wrap(~ condition)`; with multiple observables it switches
#'   to `facet_grid(name ~ condition)` (observables in rows, subjects in
#'   columns). Use `subjectsPerPage` to split very large cohorts into several
#'   plots.
#' @param x An [EM].
#' @param times Optional grid of additional times for the smooth IPRED/PRED.
#' @param ncol Facet column count for the single-observable
#'   `facet_wrap(~ condition)` layout (default 4). Ignored in
#'   `facet_grid(name ~ condition)` mode.
#' @param showPred Logical; overlay population PRED dashed (default TRUE).
#' @param showBand Logical; draw the IPRED +/- sigma ribbon if the fit carries
#'   an errfn (default TRUE).
#' @param subjectsPerPage Optional integer. If set, subjects are split into
#'   pages of at most `subjectsPerPage` each and the function returns a list
#'   of ggplots (one per page, page index appended to the title). `NULL`
#'   (default) keeps the single-plot behaviour.
#' @param ... Ignored.
#' @return A ggplot, or a list of ggplots when `subjectsPerPage` is set.
#' @export
plotIndivs.em <- function(x, times = NULL, ncol = 4L,
                              showPred = TRUE, showBand = TRUE,
                              subjectsPerPage = NULL, ...) {
  .emRequireOmega(x, "plotIndivs")
  fit <- x
  obs_times <- sort(unique(unlist(lapply(fit$data, `[[`, "time"))))
  if (is.null(times))
    times <- seq(min(obs_times), max(obs_times), length.out = 200L)

  pf <- predict(fit, times = times)
  # Trim grid rows to the observed time range. predict..fitNormal prepends t=0 to
  # the dense grid for the ODE solver, but if no subject has an observation at
  # that time the grid IPRED/PRED there can be far outside the data range
  # (e.g. log(Cc + eps) at Cc(0)=0 sits at log(eps)), which would dominate the
  # y-axis. Observation rows are kept verbatim.
  pf <- pf[pf$source == "obs" |
             (pf$time >= min(obs_times) & pf$time <= max(obs_times)),
           , drop = FALSE]

  n_obs_names <- length(unique(pf$name))
  cond_levels <- if (is.factor(pf$condition))
    levels(droplevels(pf$condition)) else unique(as.character(pf$condition))

  build_page <- function(page_pf, page_label = NULL) {
    obs <- page_pf[page_pf$source == "obs", , drop = FALSE]
    grd <- page_pf[order(page_pf$condition, page_pf$name, page_pf$time), ,
                   drop = FALSE]
    has_band <- showBand && any(is.finite(grd$sigma))

    p <- ggplot2::ggplot()
    if (has_band) {
      p <- p + ggplot2::geom_ribbon(
        data = grd,
        ggplot2::aes(x = time,
                     ymin = IPRED - sigma, ymax = IPRED + sigma),
        fill = dMod_colors[3], alpha = 0.2, linetype = 0)
    }
    p <- p +
      ggplot2::geom_line(data = grd,
                         ggplot2::aes(x = time, y = IPRED, color = "IPRED"),
                         linewidth = 0.7)
    if (showPred) {
      p <- p + ggplot2::geom_line(
        data = grd,
        ggplot2::aes(x = time, y = PRED, color = "PRED"),
        linewidth = 0.5, linetype = "dashed")
    }
    p <- p +
      ggplot2::geom_point(data = obs,
                          ggplot2::aes(x = time, y = observed),
                          size = 1.6, alpha = 0.85)

    if (n_obs_names > 1L) {
      p <- p + ggplot2::facet_grid(name ~ condition, scales = "free_y")
    } else {
      p <- p + ggplot2::facet_wrap(~ condition, ncol = ncol, scales = "free_y")
    }

    title <- sprintf("Individual fits (method = %s)", fit$method)
    if (!is.null(page_label))
      title <- paste0(title, " ", page_label)

    p +
      ggplot2::scale_color_manual(values = c(IPRED = dMod_colors[3],
                                             PRED  = dMod_colors[2])) +
      ggplot2::labs(x = "Time", y = "Value", color = NULL, title = title) +
      theme_dMod(base_size = 11)
  }

  if (is.null(subjectsPerPage))
    return(build_page(pf))

  subjectsPerPage <- as.integer(subjectsPerPage)
  if (length(subjectsPerPage) != 1L || is.na(subjectsPerPage) ||
      subjectsPerPage < 1L)
    stop("`subjectsPerPage` must be a positive integer or NULL.")

  pages <- split(cond_levels,
                 ceiling(seq_along(cond_levels) / subjectsPerPage))
  n_pages <- length(pages)
  lapply(seq_along(pages), function(i) {
    page_conds <- pages[[i]]
    sub <- pf[as.character(pf$condition) %in% page_conds, , drop = FALSE]
    if (is.factor(sub$condition))
      sub$condition <- factor(sub$condition, levels = page_conds)
    label <- sprintf("(page %d/%d)", i, n_pages)
    build_page(sub, page_label = label)
  })
}

#' @rdname plotHistIndivs
#' @param x An [EM].
#' @param ... Ignored.
#' @return A ggplot (or a list with `hist` and `qq` if cowplot is unavailable).
#' @export
plotHistIndivs.em <- function(x, ...) {
  .emRequireOmega(x, "plotHistIndivs")
  fit <- x
  etaModes <- fit$etaModes
  om <- fit$omega
  if (is.null(etaModes) || is.null(om))
    stop("plotHistIndivs..fitNormal: fit has no etaModes or omega.")
  Omega <- if (!is.null(fit$Omega)) fit$Omega else {
    L <- om$buildL(fit$argument[om$cholPars])
    tcrossprod(L)
  }
  K <- ncol(etaModes)
  eta_long <- do.call(rbind, lapply(seq_len(K), function(k) {
    sd_k <- sqrt(Omega[k, k])
    data.frame(eta_name = om$eta[k],
               value    = etaModes[, k],
               sd_est   = sd_k,
               stringsAsFactors = FALSE)
  }))
  # Per-eta Gaussian curve on a wide grid so the full N(0, Omega_kk) shape is
  # visible; using stat_function with a single mean(sd_est) painted all panels
  # with the same density (wrong when omega varies across etas) and was
  # clipped to the data extent by facet_wrap(scales = "free"). geom_line on
  # the explicit grid both picks up the panel-specific SD and widens the
  # x-range to ~ +/- 3.5 sigma.
  dens_long <- do.call(rbind, lapply(seq_len(K), function(k) {
    sd_k  <- sqrt(Omega[k, k])
    xlim  <- max(abs(etaModes[, k]), 3.5 * sd_k)
    xseq  <- seq(-xlim, xlim, length.out = 200L)
    data.frame(eta_name = om$eta[k],
               x        = xseq,
               density  = dnorm(xseq, 0, sd_k),
               stringsAsFactors = FALSE)
  }))
  p_hist <- ggplot2::ggplot(eta_long, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                            bins = 12, fill = "grey80", color = "grey40") +
    ggplot2::geom_line(data = dens_long,
                       ggplot2::aes(x = x, y = density),
                       inherit.aes = FALSE,
                       color = dMod_colors[2], linewidth = 0.8) +
    ggplot2::facet_wrap(~ eta_name, scales = "free") +
    ggplot2::labs(x = "eta value", y = "density",
                  title = "Eta distribution vs N(0, Omega_kk)") +
    theme_dMod(base_size = 11)
  p_qq <- ggplot2::ggplot(eta_long,
                          ggplot2::aes(sample = value / sd_est)) +
    ggplot2::stat_qq() +
    ggplot2::stat_qq_line(color = dMod_colors[2]) +
    ggplot2::facet_wrap(~ eta_name) +
    ggplot2::labs(x = "Theoretical quantile (N(0,1))",
                  y = "Standardised eta",
                  title = "Eta QQ vs N(0,1)") +
    theme_dMod(base_size = 11)
  if (requireNamespace("cowplot", quietly = TRUE))
    cowplot::plot_grid(p_hist, p_qq, ncol = 1)
  else list(hist = p_hist, qq = p_qq)
}

#' @rdname plotTrace
#' @param x An [EM] with `method = "quadrature"`. Errors otherwise.
#' @param ... Ignored.
#' @return A ggplot.
#' @export
plotTrace.em <- function(x, ...) {
  .emRequireOmega(x, "plotTrace")
  fit <- x
  if (!fit$method %in% c("quadrature", "foceiQuadrature") ||
      is.null(fit$stageTrace))
    stop("plotTrace requires an .fitNormal fit with method 'quadrature' or 'foceiQuadrature'.")
  tr <- fit$stageTrace
  tr$iter <- seq_len(nrow(tr))
  long <- rbind(
    data.frame(iter = tr$iter, value = tr$OFV,        panel = "OFV",
               level = factor(tr$level)),
    data.frame(iter = tr$iter, value = tr$deltaPsi,  panel = "|delta psi|",
               level = factor(tr$level)),
    data.frame(iter = tr$iter, value = tr$maxSoftmax, panel = "max softmax",
               level = factor(tr$level)),
    data.frame(iter = tr$iter, value = tr$nEffMin,   panel = "n_eff min",
               level = factor(tr$level)))
  ggplot2::ggplot(long, ggplot2::aes(x = iter, y = value, color = level)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::facet_wrap(~ panel, scales = "free_y") +
    scale_color_dMod() +
    ggplot2::labs(x = "ECM iteration",
                  title = "ECM convergence trace") +
    theme_dMod(base_size = 11)
}
