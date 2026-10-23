## parlist: the fit results of mstrust() ---------------------------------------

#' @param x list of fit results, each as returned by [trust()]. Default `NULL`
#'   returns `NULL`.
#' @rdname parlist
#' @export
as.parlist <- function(x = NULL) {
  if (is.null(x)) return(NULL)
  structure(x, class = c("parlist", "list"))
}

#' @export
#' @rdname parlist
print.parlist <- function(x, ...) {

  if (length(x) == 0L) {
    cat("Empty parlist (no fits).\n")
    return(invisible(x))
  }

  m_stat <- .statParlist(x)

  cat("Parameter list of", length(x), if (length(x) == 1L) "fit\n" else "fits\n")
  cat("... converged:     ", sum(m_stat == "converged"), "\n")
  cat("... not converged: ", sum(m_stat == "notconverged"), "\n")
  cat("... aborted:       ", sum(m_stat == "error"), "\n")
  cat("\nUse summary() for the best/worst fit, as.parframe() for a data frame.\n")

  invisible(x)

}

#' @export
#' @param object object of class `parlist`.
#' @rdname parlist
summary.parlist <- function(object, ...) {

  status <- factor(.statParlist(object), levels = c("error", "notconverged", "converged"))
  counts <- table(status)

  frame <- as.parframe(object)
  if (nrow(frame)) {
    extremes <- frame[order(frame$value)[c(1L, nrow(frame))], ]
    rownames(extremes) <- c("best", "worst")
    cat("Best and worst fit\n")
    print(extremes)
  }

  tally <- c("Fits aborted:" = counts[["error"]], "Fits not converged:" = counts[["notconverged"]],
             "Fits converged:" = counts[["converged"]], "Fits total:" = length(status))
  cat("\nStatistics of fit outcome\n", sprintf("%-20s%d\n", names(tally), tally), sep = "")

  reasons <- .stopReasonTable(object)
  if (!is.null(reasons)) {
    cat("\n\nTermination reason\n")
    cat(sprintf("%-20s%d\n", names(reasons), as.integer(reasons)), sep = "")
  }
  invisible(object)
}



# Convergence status per fit: "error", "converged" or "notconverged"
.statParlist <- function(x) {
  vapply(x, function(fit) {
    if (is.null(fit) || inherits(fit, "try-error") || !is.null(fit$error)) return("error")
    if (isTRUE(fit$converged)) "converged" else "notconverged"
  }, "", USE.NAMES = FALSE)
}


## Termination reason per fit, as reported by trust()$stopReason. NA for fits
## that errored or lack the field. "gradient" is the only certified stop;
## "stagnation" and "iterlim" mean the run ran out of resolution or budget.
.stopReasonParlist <- function(x) {
  vapply(x, function(fit) {
    if (inherits(fit, "try-error") || any(names(fit) == "error") || is.null(fit))
      return(NA_character_)
    if (is.null(fit$stopReason)) NA_character_ else as.character(fit$stopReason)
  }, NA_character_, USE.NAMES = FALSE)
}


## Tabulate termination reasons for printing; empty when nothing reports one.
.stopReasonTable <- function(x) {
  reasons <- .stopReasonParlist(x)
  reasons <- reasons[!is.na(reasons)]
  if (!length(reasons)) return(NULL)
  sort(table(reasons), decreasing = TRUE)
}


#' Plot a Parameter List
#'
#' Plots, for the converged fits, the converged against the initial value of
#' each parameter, or with `path = TRUE` the parameter values over the
#' iterations.
#'
#' @param x object of class `parlist`, e.g. from [mstrust()].
#' @param path logical, plot the path of the parameters from the initial
#'   values to convergence. Default `FALSE`. `path = TRUE` needs fits run with
#'   `trust(..., blather = TRUE)` or `mstrust(..., blather = TRUE)`.
#' @param ... not used.
#' @return A `ggplot` object.
#' @seealso [parlist()], [plotValues()]
#'   
#' @export
plot.parlist <- function(x, path = FALSE, ...) {

  fits <- x[.statParlist(x) == "converged"]
  if (!length(fits)) stop("No converged fit to plot.", call. = FALSE)
  if (path) return(.plotParlistPath(fits))

  start <- do.call(rbind, lapply(fits, `[[`, "parinit"))
  end <- do.call(rbind, lapply(fits, `[[`, "argument"))
  frame <- data.frame(parameter = rep(colnames(start), each = nrow(start)),
                      x = as.vector(start), y = as.vector(end), run = 1)
  ggplot(frame, aes(x = x, y = y)) + geom_point() + facet_wrap(~parameter)
}

# Iterates of every fit, one line per fit and one panel per parameter
.plotParlistPath <- function(fits) {
  if (!"argpath" %in% names(fits[[1]]))
    stop("No path information in the fits. Run trust() or mstrust() with blather = TRUE.")
  pars <- names(fits[[1]][["parinit"]])
  steps <- lapply(seq_along(fits), function(i) {
    iterates <- fits[[i]][["argpath"]]
    data.frame(parameter = rep(pars, each = nrow(iterates)),
               iteration = rep(seq_len(nrow(iterates)), length(pars)),
               path = as.vector(iterates), idx = as.character(i))
  })
  ggplot(do.call(rbind, steps), aes(x = iteration, y = path, colour = idx)) +
    geom_line() + facet_wrap(~parameter)
}


#' @export
#' @rdname as.parframe
#' @param sortBy character, the column by which the returned parameter frame
#' is sorted. Default `"value"`.
as.parframe.parlist <- function(x, sortBy = "value", ...) {

  keep <- which(.statParlist(x) != "error")
  fits <- x[keep]
  field <- function(name, type, convert = identity)
    vapply(fits, function(f) if (is.null(f[[name]])) type[NA_integer_] else convert(f[[name]]), type)

  out <- data.frame(index = keep,
                    value = vapply(fits, function(f) f[["value"]], numeric(1)),
                    converged = vapply(fits, function(f) isTRUE(f[["converged"]]), logical(1)),
                    iterations = vapply(fits, function(f) as.integer(f[["iterations"]]), integer(1)))

  # Diagnostics enter only when at least one fit reports them: the termination
  # reason, the objective evaluations (all, quasi-Newton, skipped) and the
  # Hessian source handovers.
  extra <- list(stopReason = .stopReasonParlist(x)[keep])
  for (name in c("neval", "qnEval", "qnSkipped", "nSwitch"))
    extra[[name]] <- field(name, integer(1), as.integer)

  # Iterations per element of trust()'s `phases`, summed over its runs, when
  # the fits had more than one.
  phases <- unique(unlist(lapply(fits, function(f)
    names(f[["phaseSweep"]]) %||% names(f[["phaseIterations"]]))))
  if (length(phases) > 1L) for (ph in phases)
    extra[[paste0("iterations.", ph)]] <- vapply(fits, function(f) {
      it <- f[["phaseIterations"]]
      if (is.null(it)) NA_integer_ else as.integer(sum(it[names(it) == ph]))
    }, integer(1))
  extra$elapsed <- field("elapsed", numeric(1))

  extra <- Filter(function(col) !all(is.na(col)), extra)
  out[names(extra)] <- extra
  meta <- names(out)

  # Fits may differ in their parameters; a parameter a fit lacks is NA there.
  pars <- as.character(unique(unlist(lapply(fits, function(f) names(f[["argument"]])))))
  for (par in pars)
    out[[par]] <- vapply(fits, function(f) as.numeric(f[["argument"]][par]), numeric(1))

  out <- parframe(out[order(out[[sortBy]]), , drop = FALSE], parameters = pars, metanames = meta)
  attr(out, "run") <- attr(x, "run")
  out
}



#'
#' @rdname parlist
#' @export
c.parlist <- function(...) {
  # a single fit is one entry, not a list of entries
  parts <- lapply(list(...), function(a)
    if (!inherits(a, "parlist") && is.list(a) && !is.null(a$argument)) list(a) else unclass(a))
  as.parlist(do.call(c, parts))
}


## Methods for the class parframe ----


#' Coerce to a Parameter Frame
#' 
#' @param x object to be coerced, e.g. a `parlist`.
#' @param ... not used.
#' @return Object of class [parframe]. For a `parlist`, the meta-information
#' columns are `index`, `value`, `converged` and `iterations`, plus the
#' termination and evaluation counts and the wall time `elapsed` the fits
#' report. Fits of [trust()] with several `phases` add one column
#' `iterations.<phase>` per phase, the iterations summed over its runs. The
#' attribute `"run"` of the parlist is kept.
#' @seealso [parframe()], [as.parvec.parframe()], and [parlist()] for an
#'   example.
#' @export
as.parframe <- function(x, ...) {
  UseMethod("as.parframe", x)
}


#' Select a Parameter Vector from a Parameter Frame
#'
#' @param x object of class `parframe`, e.g. from [as.parframe()].
#' @param index integer, select the row with the `index`-th lowest objective
#'   value. Default 1, the best fit. Without a `value` column, the `index`-th
#'   row.
#' @param ... not used.
#'
#' @return Object of class `parvec` with the parameter columns of the selected
#'   row. A warning is issued if that fit did not converge.
#' @seealso [as.parframe()], [parvec()]
#'
#' @export
as.parvec.parframe <- function(x, index = 1, ...) {
  if (length(index) != 1L || is.na(index) || index < 1 || index > nrow(x))
    stop("'index' must be one number from 1 to ", nrow(x), ", the rows of the parframe.",
         call. = FALSE)
  metanames <- attr(x, "metanames")
  row <- if ("value" %in% metanames) order(x$value)[index] else index
  if ("converged" %in% metanames && isFALSE(x$converged[row]))
    warning("The selected parameter vector is from an unconverged fit.", call. = FALSE)
  .parvecRow(x, row)
}

# Row i of a parframe as parvec, in the row order of the frame
.parvecRow <- function(x, i) {
  as.parvec(unlist(as.data.frame(x)[i, attr(x, "parameters"), drop = FALSE]))
}




#' @export
#' @rdname plotPars
plotPars.parframe <- function(x, ..., tol = 1) {

  if (!missing(...)) x <- subset(x, ...)

  x <- x[order(x$value), ]
  x$step <- as.factor(findInterval(seq_len(nrow(x)), .stepDetect(x$value, tol)))

  parNames <- attr(x, "parameters")
  parOut <- wide2long.data.frame(out = x[, c("step", "value", parNames)], keep = 1:2)
  names(parOut) <- c("step", "value", "name", "parvalue")
  plot <- ggplot2::ggplot(parOut, aes(x = name, y = parvalue, color = step)) +
    geom_boxplot(outlier.alpha = 0) + theme_dMod() + scale_color_dMod() +
    guides(x = guide_axis(angle = 90))

  attr(plot, "data") <- parOut
  plot
}


# 0 and multiples of the decades from 10 up to the largest absolute limit, both
# signs, refined from 1 over 1-2-5 to 1..9 until three fall inside the limits.
# On a pseudo-log axis, which is linear below 1, a break at 1 would crowd 0.
.decadeBreaks <- function(lims) {
  d <- 10^seq_len(max(1, ceiling(log10(max(abs(lims), 1)))))
  for (m in list(1, c(1, 2, 5), 1:9)) {
    p <- as.vector(outer(m, d))
    b <- sort(unique(c(-p, 0, p)))
    if (sum(b >= lims[1] & b <= lims[2]) >= 3L) return(b)
  }
  pretty(lims)
}


#' @export
#' @rdname plotValues
plotValues.parframe <- function(x, ..., tol = 1, showSteps = FALSE, pseudoLog = FALSE) {

  if (!missing(...)) x <- subset(x, ...)

  pars <- x[order(x$value), ]
  pars[["index"]] <- seq_len(nrow(pars))
  pars[["delta"]] <- pars$value - pars$value[1]
  jumps <- .stepDetect(pars$value, tol)

  yScale <- if (pseudoLog)
    scale_y_continuous(transform = scales::pseudo_log_trans(base = 10), breaks = .decadeBreaks)

  stepLines <- stepLabels <- NULL
  if (showSteps) {
    stepLines <- geom_vline(xintercept = jumps, lty = 2)
    stepLabels <- annotate("text", x = jumps + 1, y = pars$value[jumps], label = jumps,
                           hjust = 0, vjust = -0.5, color = "firebrick", size = 3)
  }

  # an unconverged fit is red, outside the ramp of the iterations, with a legend
  # entry of its own
  conv <- pars$converged %in% TRUE
  unconverged <- if (!all(conv)) list(
    geom_point(data = pars[!conv, ], aes(shape = "not converged"), color = dMod_colors[2]),
    scale_shape_manual(values = c(`not converged` = 19),
                       guide = guide_legend(override.aes = list(color = dMod_colors[2]),
                                            order = 1)))

  P <- ggplot2::ggplot(pars, aes(x = index, y = value)) +
    stepLines +
    geom_point(data = pars[conv, ], aes(color = iterations)) +
    unconverged +
    stepLabels +
    yScale +
    labs(x = "Fit rank", y = "Objective value", color = "Iterations", shape = NULL) +
    scale_color_dMod_c() +
    theme_dMod()

  structure(P, data = pars, jumps = jumps)

}


# Default lines of a profile plot: the chi-square thresholds, labelled with the
# level they belong to. The names are the labels, so a caller can hand in the
# thresholds its intervals are actually read at.
.profileLines <- function() {
  c("68%" = 1, "90%" = qchisq(0.90, 1), "95%" = qchisq(0.95, 1))
}

.profileLineLabels <- function(threshold) {
  labels <- names(threshold)
  if (is.null(labels)) labels <- formatC(threshold, format = "f", digits = 2)
  labels
}


# Profile sets covered by one call. Names label the colour scale.
.profileSets <- function(profs) {
  sets <- if (inherits(profs, "parframe")) list(profs) else as.list(profs)
  if (is.null(names(sets))) names(sets) <- seq_along(sets)
  sets
}

# One data.frame per profiled parameter. Failed profiles are dropped.
.profileSplit <- function(x) {
  prof <- as.data.frame(x)
  if (is.data.frame(prof)) prof <- split(prof, prof[["whichPar"]])
  ok <- vapply(prof, is.data.frame, logical(1))
  if (!all(ok)) warning(sum(!ok), " profiles discarded.", call. = FALSE)
  prof[ok]
}

# Long format for one parameter. `delta` is measured against the profile's own
# origin, where the constraint vanishes.
.profileDeltas <- function(prof, name, set, modes, maxvalue) {
  origin <- which.min(abs(prof[["constraint"]]))
  block <- function(column, mode)
    data.frame(name    = name,
               delta   = prof[[column]] - prof[[column]][origin],
               par     = prof[[name]],
               proflist = set,
               mode    = mode,
               is.zero = seq_len(nrow(prof)) == origin)
  out <- do.call(rbind, c(list(block("value", "total")),
                          lapply(modes, function(m) block(m, m))))
  out[which(out$delta <= maxvalue), , drop = FALSE]
}

# The frame behind every profile plot, over all sets and parameters.
.profileFrame <- function(sets, maxvalue) {
  out <- do.call(rbind, lapply(names(sets), function(s) {
    modes <- attr(sets[[s]], "obj.attributes")
    parts <- .profileSplit(sets[[s]])
    do.call(rbind, lapply(names(parts), function(n)
      .profileDeltas(parts[[n]], n, s, modes, maxvalue)))
  }))
  out$proflist <- as.factor(out$proflist)
  out
}

# Reference for the `parlist` points: lowest value at any profile origin.
.profileReference <- function(sets) {
  min(vapply(sets, function(p) {
    d <- as.data.frame(p)
    d[[which.min(abs(d[["constraint"]])), 1L]]
  }, numeric(1)))
}

# `parlist` overlay. A `value` column places points at their distance to the
# reference, otherwise on the zero line.
.profilePoints <- function(parlist, sets) {
  frame <- as.data.frame(parlist)
  pars <- intersect(colnames(frame), unlist(lapply(sets, attr, "parameters")))
  height <- if (is.null(frame[["value"]])) 0 else frame[["value"]] - .profileReference(sets)
  points <- data.frame(par = as.numeric(unlist(frame[pars], use.names = FALSE)),
                       name = rep(pars, each = nrow(frame)),
                       delta = as.numeric(height))
  geom_point(data = points, aes(x = par, y = delta), color = "black",
             inherit.aes = FALSE)
}

# Shared body of the plotProfile methods.
.profilePlot <- function(data, sets, parlist, ncol, maxvalue, threshold) {

  # The optimum is always drawn, whatever the thresholds are.
  keep      <- threshold != 0
  labels    <- c("optimum", .profileLineLabels(threshold)[keep])
  threshold <- c(0, threshold[keep])

  p <- ggplot(data, aes(x = par, y = delta, group = interaction(proflist, mode),
                        color = proflist, linetype = mode)) +
    facet_wrap(~name, scales = "free_x", ncol = ncol) +
    geom_hline(yintercept = threshold, lty = 2, color = "gray") +
    geom_line() +
    geom_point(data = subset(data, is.zero)) +
    scale_y_continuous(breaks = threshold, labels = labels,
                       limits = c(NA, maxvalue)) +
    xlab("parameter value") + ylab("Confidence Level") +
    # a narrow panel drops the x labels that would overlap and keeps apart from the next
    guides(x = guide_axis(check.overlap = TRUE)) +
    theme(panel.spacing.x = unit(1.2, "lines"))

  # A single set needs no colour legend to say so.
  p <- if (nlevels(data$proflist) < 2L) p + guides(color = "none")
       else p + labs(color = "set")

  if (!is.null(parlist)) p <- p + .profilePoints(parlist, sets)

  attr(p, "data") <- data
  p
}


#' @export
#' @rdname plotProfile
plotProfile.parframe <- function(profs, ..., maxvalue = 5, parlist = NULL, ncol = NULL,
                                 threshold = .profileLines()) {
  sets <- .profileSets(profs)
  data <- .profileFrame(sets, maxvalue)
  data <- droplevels(data[.filterRows(data, ...), , drop = FALSE])
  .profilePlot(data, sets, parlist, ncol, maxvalue, threshold)
}


#' @export
#' @rdname plotProfile
plotProfile.list <- function(profs, ..., maxvalue = 5, parlist = NULL, ncol = NULL,
                             threshold = .profileLines()) {
  sets <- .profileSets(profs)
  data <- .profileFrame(sets, maxvalue)
  data <- droplevels(data[.filterRows(data, ...), , drop = FALSE])
  .profilePlot(data, sets, parlist, ncol, maxvalue, threshold)
}



#' @export
#' @rdname parframe
is.parframe <- function(x) {
  inherits(x, "parframe")
}

#' @export
#' @param i row index in any format; missing keeps all rows. Used alone, as in
#'   `x[j]`, it selects columns, as for a data frame.
#' @param j column index in any format; missing keeps all columns.
#' @param drop logical. Default `FALSE`. `TRUE` returns a plain data frame or
#'   vector instead of a `parframe`.
#' @rdname parframe
"[.parframe" <- function(x, i, j, drop = FALSE) {

  out <- as.data.frame(x)
  nIndex <- nargs() - if (missing(drop)) 0L else 1L
  if (nIndex == 2L) {
    if (!missing(i)) out <- out[i]
  } else {
    if (!missing(i)) out <- out[i, , drop = FALSE]
    if (!missing(j)) out <- out[, j, drop = drop]
  }
  if (drop) return(out)

  out <- parframe(out,
                  parameters = intersect(attr(x, "parameters"), colnames(out)),
                  metanames = intersect(attr(x, "metanames"), colnames(out)),
                  objAttributes = intersect(attr(x, "obj.attributes"), colnames(out)))
  attr(out, "run") <- attr(x, "run")
  out
}


#' @export
#' @param ... for `subset()`, logical expressions in the columns of `x`; a row
#'   is kept where all of them are `TRUE`. For `rbind()`, parameter frames with
#'   the same columns.
#' @rdname parframe
subset.parframe <- function(x, ...) {
  x[.filterRows(x, ...), ]
}

#' @details `rbind()` keeps the attributes of the first frame. Profiles of one
#' parameter computed in parts, such as its left and right side, are joined
#' into one: the repeated optimum is dropped and the points are sorted by
#' `constraint`.
#' @param deparse.level not used.
#' @export
#' @rdname parframe
rbind.parframe <- function(..., deparse.level = 1) {
  frames <- Filter(Negate(is.null), list(...))
  out <- do.call(rbind.data.frame, c(unname(lapply(frames, as.data.frame)),
                                     make.row.names = FALSE))
  if (all(c("whichPar", "constraint") %in% names(out))) {
    out <- out[order(match(out$whichPar, unique(out$whichPar)), out$constraint), ]
    out <- out[!(out$constraint == 0 & duplicated(out[c("whichPar", "constraint")])), ]
  }
  first <- frames[[1L]]
  out <- parframe(out, parameters = attr(first, "parameters"),
                  metanames = attr(first, "metanames"),
                  objAttributes = attr(first, "obj.attributes"))
  attr(out, "run") <- attr(first, "run")
  out
}

# The rows of a data frame where every expression in `...` holds, NA counting as
# FALSE. An expression sees the columns of `x` first, then the environment it was
# written in, also when passed on through the `...` of another function.
.filterRows <- function(x, ...) {
  keep <- rep(TRUE, nrow(x))
  for (q in rlang::enquos(...)) keep <- keep & rlang::eval_tidy(q, data = as.list(x))
  keep & !is.na(keep)
}

#' Rows of a Parameter Frame with Unique Values
#'
#' Keeps one row per distinct objective value.
#' @param x parameter frame
#' @param incomparables not used. Argument exists for compatibility with S3 generic.
#' @param tol tolerance to decide when values are assumed to be equal, see [plotValues()].
#' @param ... Not used. Select rows first with [subset()], as in
#'   `unique(subset(x, converged))`.
#' @return The rows of `x` with distinct values, sorted by value.
#' @export
unique.parframe <- function(x, incomparables = FALSE, tol = 1, ...) {
  if (!nrow(x)) return(x)
  x <- x[order(x$value), ]
  x[.stepDetect(x$value, tol), ]
}



## Methods for the class parvec ------------------------------------------------

#' @param x numeric vector of parameter values.
#' @param names character, the parameter names. Default `NULL` uses
#'   `names(x)`.
#' @param deriv2 3D array of second derivatives, `[parameter, theta, theta]`.
#'   Default `NULL` keeps the `"deriv2"` attribute of `x`, `FALSE` drops it.
#'
#' @export
#' @rdname parvec
as.parvec <- function(x, ...) {
  UseMethod("as.parvec", x)
}


#' @export
#' @rdname parvec
as.parvec.numeric <- function(x, names = NULL, deriv = NULL, deriv2 = NULL, ...) {

  # --- Basic setup ---
  p <- as.numeric(x)
  if (is.null(names)) names(p) <- names(x) else names(p) <- names
  pnames <- names(p)

  # --- Derivative Information ---
  if (isFALSE(deriv)) {
    full_deriv <- NULL
  } else if (is.matrix(deriv)) {
    full_deriv <- deriv
  } else { # deriv == NULL
    full_deriv <- attr(x, "deriv")
  }

  # --- Second-order derivative information ---
  # `deriv2` is a 3D array [p, theta, theta], symmetric in the last two axes.
  if (isFALSE(deriv2)) {
    full_deriv2 <- NULL
  } else if (is.array(deriv2) && length(dim(deriv2)) == 3L) {
    full_deriv2 <- deriv2
  } else { # deriv2 == NULL
    full_deriv2 <- attr(x, "deriv2")
  }

  # --- Infer fixed from missing deriv rows ---
  fixed <- NULL
  if (!is.null(full_deriv)) {
    if (nrow(full_deriv) < length(pnames)) {
      fixed <- .setdiffU(pnames, rownames(full_deriv))
    }
  }

  # --- Assemble object ---
  attr(p, "deriv") <- full_deriv
  attr(p, "deriv2") <- full_deriv2
  attr(p, "fixed") <- fixed
  class(p) <- c("parvec", "numeric")
  p
}




#' Print Parameter Vectors
#'
#' Prints a parameter vector with information about its derivatives and its
#' fixed parameters.
#'
#' @param x object of class `parvec`.
#' @param ... not used.
#' @return `x`, invisibly.
#' @export
print.parvec <- function(x, ...) {
  
  par <- unclass(x)
  nms <- names(par)
  n_width <- max(0L, nchar(nms))
  
  cat("Parameter vector:\n")
  if (length(par))
    cat(sprintf("  %s : %s\n", format(nms, width = n_width, justify = "right"),
                formatC(par, digits = 6, format = "g", flag = " ")), sep = "")
  
  # "<rows> x <cols> matrix" for a derivative, "<none>" without one
  shape <- function(a, kind)
    if (is.null(a)) "<none>" else paste(paste(dim(a), collapse = " x "), kind)
  fixed <- attr(x, "fixed")
  cat("\nAttributes:\n",
      sprintf("  deriv  : %s\n", shape(attr(x, "deriv"), "matrix")),
      sprintf("  deriv2 : %s\n", shape(attr(x, "deriv2"), "array")),
      sprintf("  fixed  : %s\n", if (length(fixed)) paste(fixed, collapse = ", ") else "<none>"),
      sep = "")

  invisible(x)
}


#' Subset a Parameter Vector
#'
#' Subsets a `parvec` and its derivatives.
#'
#' @param x object of class `parvec`.
#' @param ... subsetting indices.
#' @param drop logical, drop derivative columns that are zero after
#'   subsetting. Default `FALSE`.
#'
#' @return Object of class `parvec`.
#'
#' @export
"[.parvec" <- function(x, ..., drop = FALSE) {

  # `.subset()` subsets without dispatch and without copying the attributes:
  # `unclass(x)[...]` duplicates the values and the deriv arrays first, only to
  # drop them again.
  out <- .subset(x, ...)
  nms <- names(out)
  # Reselecting every name in order reproduces x exactly, and the composition
  # protocol does it once per condition.
  if (!drop && !is.null(nms) && identical(nms, names(x))) return(x)

  deriv <- attr(x, "deriv")
  if (!drop) {
    # The row selection and the `fixed` derivation are the whole cost here, and
    # this runs a few hundred times per objective evaluation.
    fast <- parvec_attach(out, deriv, attr(x, "deriv2"))
    if (!is.null(fast)) return(fast)
  }
  if (!is.null(deriv)) {
    rn <- rownames(deriv)
    # Row set unchanged and in order: keep the matrix, skip the subset copy.
    if (!identical(rn, nms)) {
      available <- nms[match(nms, rn, 0L) > 0L]
      deriv <- if (length(available)) deriv[available, , drop = FALSE] else NULL
    }
  }

  deriv2 <- attr(x, "deriv2")
  if (!is.null(deriv2)) {
    rn2 <- dimnames(deriv2)[[1]]
    if (!identical(rn2, nms)) {
      available2 <- nms[match(nms, rn2, 0L) > 0L]
      deriv2 <- if (length(available2)) deriv2[available2, , , drop = FALSE] else NULL
    }
  }

  if (drop && !is.null(deriv)) {
    keep.cols <- colSums(abs(deriv)) > 0
    deriv <- deriv[, keep.cols, drop = FALSE]
    if (!is.null(deriv2)) {
      keep_names <- colnames(deriv)
      deriv2 <- deriv2[, keep_names, keep_names, drop = FALSE]
    }
  }

  # Assemble in place: `as.parvec()` would re-derive what is already known and
  # copy the values once more.
  attr(out, "deriv")  <- deriv
  attr(out, "deriv2") <- deriv2
  attr(out, "fixed")  <- if (!is.null(deriv) && nrow(deriv) < length(nms))
    .setdiffU(nms, rownames(deriv))
  class(out) <- c("parvec", "numeric")
  out
}

#' Concatenate Parameter Vectors
#'
#' Concatenates `parvec` objects with their derivatives.
#'
#' @param ... objects of class `parvec`. `NULL` entries are ignored.
#'
#' @return Object of class `parvec`.
#'
#' @export
c.parvec <- function(...) {

  # One operand, alone or beside NULL or an empty one, as
  # prdframe(parameters = c(pars, fixed)) has it for every condition: the same
  # result as below, without the list machinery.
  n <- ...length()
  if ((n == 1L || (n == 2L && !length(..2))) && length(..1)) {
    q <- ..1
    nms <- names(q)
    if (anyDuplicated(nms)) stop("Duplicated parameter names.")
    D <- attr(q, "deriv")
    D2 <- attr(q, "deriv2")
    v <- as.numeric(q)
    names(v) <- nms
    if (is.matrix(D)) {
      attr(v, "deriv") <- D
    } else D <- NULL
    if (is.array(D2) && length(dim(D2)) == 3L) attr(v, "deriv2") <- D2
    if (!is.null(D) && nrow(D) < length(v)) attr(v, "fixed") <- .setdiffU(names(v), rownames(D))
    class(v) <- c("parvec", "numeric")
    return(v)
  }

  p <- Filter(Negate(is.null), list(...))
  stopifnot(length(p) > 0)
  # An empty operand adds no name, value or derivative row, but it does push
  # the single-operand case into the rbind/abind path below.
  keep <- vapply(p, function(z) length(z) > 0L, TRUE)
  if (any(keep)) p <- p[keep]

  nms  <- unlist(lapply(p, names), use.names = FALSE)
  if (anyDuplicated(nms)) stop("Duplicated parameter names.")

  # One operand: rbind and abind would copy the largest arrays in the chain to
  # reproduce their own input. `prdframe(parameters = c(pars, fixed))` hits this
  # once per condition.
  if (length(p) == 1L) {
    q <- p[[1]]
    return(as.parvec(as.numeric(q), names = nms,
                     deriv  = attr(q, "deriv")  %||% FALSE,
                     deriv2 = attr(q, "deriv2") %||% FALSE))
  }

  fast <- parvec_concat(p)
  if (!is.null(fast)) return(fast)

  vals <- unlist(lapply(p, unclass), use.names = FALSE)

  d <- lapply(p, attr, "deriv")
  has_deriv <- !vapply(d, is.null, TRUE)

  if (!any(has_deriv)) {
    return(as.parvec(vals, names = nms))
  }

  J_list <- Filter(Negate(is.null), d)
  J <- do.call(rbind, J_list)

  # Concatenate deriv2 along the first axis if any input has one.
  d2 <- lapply(p, attr, "deriv2")
  has_d2 <- !vapply(d2, is.null, TRUE)
  H <- NULL
  if (any(has_d2)) {
    # Resolve a common theta basis; use the first non-null array's outer dims.
    ref <- d2[which(has_d2)[1L]][[1L]]
    theta_names <- dimnames(ref)[[2L]]
    n_theta <- dim(ref)[2L]
    H_list <- lapply(seq_along(p), function(i) {
      di <- d2[[i]]
      di_names <- names(p[[i]])
      n_i <- length(di_names)
      if (is.null(di)) {
        # Pad with zeros for parvecs that have no deriv2 attribute.
        array(0, c(n_i, n_theta, n_theta),
              dimnames = list(di_names, theta_names, theta_names))
      } else {
        di
      }
    })
    # abind() spends two thirds of a concatenation on dimname bookkeeping, and
    # this runs once per condition. Same result, assembled in place.
    n_rows <- vapply(H_list, function(a) dim(a)[1L], 0L)
    margin <- function(m) {
      for (a in H_list) if (!is.null(dimnames(a)[[m]])) return(dimnames(a)[[m]])
      NULL
    }
    rn <- lapply(H_list, function(a) dimnames(a)[[1L]])
    # abind() pads a block without names with "" and drops the margin only when
    # no block has names.
    rn <- if (any(!vapply(rn, is.null, TRUE)))
      unlist(lapply(seq_along(rn), function(i) rn[[i]] %||% rep("", n_rows[i])),
             use.names = FALSE) else NULL
    H <- array(0, c(sum(n_rows), n_theta, n_theta),
               dimnames = list(rn, margin(2L), margin(3L)))
    off <- 0L
    for (i in seq_along(H_list)) {
      if (n_rows[i]) H[off + seq_len(n_rows[i]), , ] <- H_list[[i]]
      off <- off + n_rows[i]
    }
  }

  as.parvec(vals, names = nms, deriv = J, deriv2 = H)
}


## Printing a parameter transformation (parfn) ---------------------------------

#' Print Parameter Transformations
#'
#' @param x object of class `parfn`.
#' @param ... not used.
#' @return `x`, invisibly.
#' 
#' @export
print.parfn <- function(x, ...) {
  
  conditions <- attr(x, "conditions")
  parameters <- attr(x, "parameters")
  mappings <- attr(x, "mappings")
  
  cat("Parameter transformation:\n")
  str(args(x))
  cat("\n")
  cat("... conditions:", paste0(conditions, collapse = ", "), "\n")
  cat("... parameters:", paste0(parameters, collapse = ", "), "\n")
  invisible(x)
}

#' @export
summary.parfn <- function(object, ...) {
  
  x <- object
  
  conditions <- attr(x, "conditions")
  parameters <- attr(x, "parameters")
  mappings <- attr(x, "mappings")
  
  cat("Details:\n")
  if (!inherits(x, "composed")) {
    
    
    output <- lapply(1:length(mappings), function(C) {
      
      list(
        equations = attr(mappings[[C]], "equations"),
        parameters = attr(mappings[[C]], "parameters")
      )
      
    })
    names(output) <- conditions
    
    output
    
  } else {
    
    cat("\nObject is composed. See original objects for more details.\n")
    
  }
}




## parfn / parframe / parlist / parvec constructors ----------------------------------------

## Parameter classes --------------------------------------------------------


## Body of the parfn dispatcher. Package level for the same reason as
## .Pexpl_p2p(): one parfn per condition should cost state, not code.

#' Parameter Transformation Function
#'
#' Builds a function that transforms one parameter vector into another and
#' propagates the Jacobian. [P()] calls it internally; use it directly for a
#' transformation of your own.
#' @param p2p transformation for one condition, a function
#'   `p2p(p, fixed, deriv)` that maps the parameter vector `p` and the fixed
#'   values `fixed` to a new parameter vector. With `deriv = TRUE` the result
#'   has the attribute `deriv`, the Jacobian of the transformation.
#' @param parameters character, the parameters accepted by the function.
#'   Default `NULL`.
#' @param condition character, the condition for which the transformation is
#'   defined. Default `NULL`, every condition.
#' @return Object of class `parfn`, a function
#'   `p(..., fixed, deriv, deriv2, hessian, conditions, env, cores, sweep)`
#'   that takes the parameter vector `pars` through `...`. Its attributes are
#'   `"mappings"`, the list of `p2p` functions, `"parameters"`, the union of
#'   their parameters, and `"conditions"`.
#' @seealso [+.fn], [P()]
#' @inheritSection dModfn Calling a dMod function
#' @example inst/examples/prediction.R
#' @export
parfn <- function(p2p, parameters = NULL, condition = NULL) {

  force(condition)
  st <- .leafState(p2p, "parfn", condition)
  outfn <- .fnWrap(st)
  attr(outfn, "mappings") <- setNames(list(p2p), condition)
  attr(outfn, "parameters") <- parameters
  attr(outfn, "conditions") <- condition
  attr(outfn, "compileInfo") <- attr(p2p, "compileInfo")
  attr(outfn, "resetWarmStart") <- attr(p2p, "resetWarmStart")
  class(outfn) <- c("parfn", "fn")
  outfn

}




#' Generate a Parameter Frame
#'
#' @description A parameter frame is a data frame with one parameter set per
#' row. Its columns fall into three groups: meta-information (e.g. index,
#' value, converged), attributes of an objective function (e.g. data and prior
#' contribution) and the parameters.
#' @seealso [as.parframe()], [profile()][profile.objfn], [mstrust()]
#' @param x data frame. Default `NULL` gives an empty parameter frame.
#' @param parameters character, the names of the parameter columns. Default
#'   all columns that are neither in `metanames` nor in `objAttributes`.
#' @param metanames character, the names of the meta-information columns.
#'   Default `NULL`.
#' @param objAttributes character, the names of the objective function
#'   attributes. Default `NULL`.
#' @return Object of class `parframe`, a data frame with the attributes
#' `parameters`, `metanames` and `obj.attributes`. `is.parframe()` returns a
#' logical.
#' @details Parameter frames can be subsetted either by `[ , ]` or by `subset`. If
#' `[ , index]` is used, the names of the removed columns will also be removed from
#' the corresponding attributes, i.e. metanames, obj.attributes and parameters.
#' @examples
#' pf <- parframe(data.frame(value = c(3, 1, 2), a = c(1, 2, 3), b = c(0, 1, 0)),
#'                parameters = c("a", "b"), metanames = "value")
#' pf[pf$value < 3, ]
#' as.parvec(pf)
#' @export
parframe <- function(x = NULL, parameters = setdiff(colnames(x), c(metanames, objAttributes)),
                     metanames = NULL, objAttributes = NULL) {
  out <- if (is.null(x)) data.frame() else as.data.frame(x)
  rownames(out) <- NULL
  structure(out, parameters = parameters, metanames = metanames,
            obj.attributes = objAttributes, class = c("parframe", "data.frame"))
}

#' Parameter List
#'
#' @description A parameter list holds the results of several optimisation
#' runs, as returned by [mstrust()]. `c()` concatenates parameter lists.
#' @param ... for `parlist()` and `c()`, fit results or parameter lists; for
#'   `print()` and `summary()`, not used.
#' @return Object of class `parlist`. `print()` and `summary()` return their
#'   argument invisibly.
#' @export
#' @example inst/examples/parlist.R
#' @seealso [load.parlist()], [plot.parlist()], [as.parframe()]
parlist <- function(...) {

  mylist <- list(...)
  return(as.parlist(mylist))

}



#' Parameter Vector
#'
#' @description
#' A parameter vector is a named numeric vector with the Jacobian of the
#' transformation that generated it in the attribute `"deriv"` and, optionally,
#' second derivatives in `"deriv2"`. Parameters without a row in the Jacobian
#' are listed in the attribute `"fixed"`. `parvec()` builds one from numbers,
#' `as.parvec()` from a numeric vector; use `as.parvec()` to set `deriv2`.
#'
#' @param ... for `parvec()`, numbers to be concatenated; unnamed ones are
#'   named `par1`, `par2`, ... For `as.parvec()`, not used.
#' @param deriv Jacobian matrix with one row per parameter and one column per
#'   parameter by which the vector was generated. Default `NULL`: none for
#'   `parvec()`, the `"deriv"` attribute of `x` for `as.parvec()`. `FALSE`
#'   drops it.
#'
#' @return Object of class `parvec`, a named numeric vector with the
#'   attributes `"deriv"`, `"deriv2"` and `"fixed"`. `parvec()` without
#'   arguments returns `NULL`.
#' @seealso [getDerivs()], [as.parvec.parframe()]
#'
#' @example inst/examples/parvec.R
#' @export
parvec <- function(..., deriv = NULL) {
  x <- c(...)
  if (!length(x)) return(NULL)
  as.parvec(.fillNames(vapply(x, as.numeric, 0), "par"), deriv = deriv)
}


