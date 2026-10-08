## Methods for the class parlist -----------------------------------------------

#' @param x list of fit results, each as returned by [trust()]. Default `NULL`
#'   returns `NULL`.
#' @rdname parlist
#' @export
as.parlist <- function(x = NULL) {
  if (is.null(x)) {
    return(NULL)
  } else {
    class(x) <- c("parlist", "list")
    return(x)
  }
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
  
  x <- object
  
  # Statistics
  m_stat <- .statParlist(x)
  m_error <- sum(m_stat == "error")
  m_converged <- sum(m_stat == "converged")
  m_notConverged <- sum(m_stat == "notconverged")
  m_sumStatus <- sum(m_error, m_converged, m_notConverged)
  m_total <- length(m_stat)
  
  # Best and worst fit
  m_parframe <- as.parframe(x)
  m_order <- order(m_parframe$value)
  m_bestWorst <- m_parframe[c(m_order[1], tail(m_order, 1)),]
  rownames(m_bestWorst) <- c("best", "worst")
  cat("Results of the best and worst fit\n")
  print(m_bestWorst)
  
  cat("\nStatistics of fit outcome",
      "\nFits aborted:       ", m_error,
      "\nFits not converged: ", m_notConverged,
      "\nFits converged:     ", m_converged,
      "\nFits total:         ", m_sumStatus, " [", m_total, "]", sep = "")

  m_reasons <- .stopReasonTable(x)
  if (!is.null(m_reasons)) {
    cat("\n\nTermination reason\n")
    for (nm in names(m_reasons))
      cat(formatC(nm, width = -20), m_reasons[[nm]], "\n", sep = "")
  }
  invisible(object)
}



## Gather statistics of a fitlist: one row per fit, holding the convergence
## status ("error" / "converged" / "notconverged"). Consumed by
## summary.parlist() and as.parframe.parlist().
.statParlist <- function(x) {
  status <- do.call(rbind, lapply(x, function(fit) {
    if (inherits(fit, "try-error") || any(names(fit) == "error") || any(is.null(fit))) {
      return("error")
    } else {
      if (fit$converged) {
        return("converged")
      } else {
        return("notconverged")
      }
    }
  }))
  
  rownames(status) <- 1:length(status)
  colnames(status) <- "fit status"

  return(status)
}


## Termination reason per fit, as reported by trust()$stopReason. NA for fits
## that errored or predate the field. "gradient" is the only certified stop;
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
#' @author Malenka Mader, \email{Malenka.Mader@@fdm.uni-freiburg.de}
#'   
#' @export
plot.parlist <- function(x, path = FALSE, ...) {
  
  pl <- x
  
  index <- do.call(rbind, lapply(pl, function(l) l$converged))
  fl <- pl[index]
  if (!path) {
    initPar <- do.call(rbind, lapply(fl, function(l) l$parinit))
    convPar <- do.call(rbind, lapply(fl, function(l) l$argument))
    
    ddata <- data.frame(cbind(matrix(initPar, ncol = 1), matrix(convPar, ncol = 1) ))
    ddata <- cbind(rep(colnames(initPar), each = nrow(initPar)), ddata, 1)
    names(ddata) <- c("parameter","x","y","run")
    
    #plot initial vs converged parameter values
    ggplot(data=ddata)+facet_wrap(~ parameter)+geom_point(aes(x=x,y=y))
  } else {
    if (!any (names(fl[[1]]) == "argpath")){
      stop("No path information in the fits. Run trust() or mstrust() with blather = TRUE.")
    }
    parNames <- names(fl[[1]]$parinit)
    
    pathPar <- do.call(rbind, mapply(function(l, idx) {
      mParPath <- as.data.frame(matrix(l$argpath, ncol = 1))
      mParPath <- cbind(rep(parNames,each = nrow(l$argpath), times = 1), rep(1:nrow(l$argpath), length(parNames)), mParPath, as.character(idx))
    }, l = fl, idx = 1:length(fl), SIMPLIFY = FALSE))
    names(pathPar) <- c("parameter", "iteration", "path", "idx")
    ggplot(data=pathPar)+geom_line(aes(x=iteration,y=path,colour=idx))+facet_wrap(~ parameter)
  }
}




#' @export
#' @importFrom data.table as.data.table rbindlist
#' @rdname as.parframe
#' @param sortBy character, the column by which the returned parameter frame
#' is sorted. Default `"value"`.
as.parframe.parlist <- function(x, sortBy = "value", ...) {
  .renameArgs(list(...), c(sort.by = "sortBy"), "as.parframe")
  m_stat <- .statParlist(x)
  m_metanames <- c("index", "value", "converged", "iterations")
  m_idx <- which("error" != m_stat)
  m_parframe <- data.frame(index = m_idx,
                           value = vapply(x[m_idx], function(.x) .x$value, 1.0),
                           converged = vapply(x[m_idx], function(.x) .x$converged, TRUE),
                           iterations = vapply(x[m_idx], function(.x) as.integer(.x$iterations), 1L))

  m_reasons <- .stopReasonParlist(x)[m_idx]
  if (any(!is.na(m_reasons))) {
    m_parframe$stopReason <- m_reasons
    m_metanames <- c(m_metanames, "stopReason")
  }

  # Objective (= gradient) evaluations, those spent in the quasi-Newton phase,
  # and the Hessian source handovers; present only when the optimiser reports
  # them.
  .intcol <- function(field) vapply(x[m_idx], function(.x)
    if (is.null(.x[[field]])) NA_integer_ else as.integer(.x[[field]]), 1L)
  for (field in c("neval", "qnEval", "qnSkipped", "nSwitch")) {
    col <- .intcol(field)
    if (any(!is.na(col))) { m_parframe[[field]] <- col; m_metanames <- c(m_metanames, field) }
  }

  parameters <- lapply(x[m_idx], function(x) data.table::as.data.table(as.list(x$argument)))
  parameters <- data.table::rbindlist(parameters, use.names = TRUE)
  m_parframe <- cbind(m_parframe, parameters)
  
  # Sort by value
  m_parframe <- m_parframe[order(m_parframe[[sortBy]]),]
  
  parframe(m_parframe, parameters = names(parameters), metanames = m_metanames)
}



#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#'
#' @rdname parlist
#' @export
#' @export c.parlist
c.parlist <- function(...) {
  m_fits <- lapply(list(...), unclass)
  m_fits <- do.call(c, m_fits)
  m_parlist <- mapply(function(fit, idx) {
    if (is.list(fit)) fit$index <- idx
    return(fit)
  }, fit = m_fits, idx = seq_along(m_fits), SIMPLIFY = FALSE)
  
  return(as.parlist(m_parlist))
}





## Methods for the class parframe ----


#' Coerce to a Parameter Frame
#' 
#' @param x object to be coerced, e.g. a `parlist`.
#' @param ... not used.
#' @return Object of class [parframe]. For a `parlist`, the meta-information
#' columns are `index`, `value`, `converged` and `iterations`, plus the
#' termination and evaluation counts the fits report.
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
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#'   
#' @export
as.parvec.parframe <- function(x, index = 1, ...) {
  parframe <- x
  m_order <- 1:nrow(x)
  metanames <- attr(parframe, "metanames")
  if ("value" %in% metanames) m_order <- order(parframe$value)
  best <- as.parvec(unlist(as.data.frame(parframe)[m_order[index], attr(parframe, "parameters"), drop = FALSE]))
  if ("converged" %in% metanames && !parframe[m_order[index],]$converged) {
    warning("Parameter vector of an unconverged fit is selected.", call. = FALSE)
  }
  return(best)
}




#' @export
#' @rdname plotPars
plotPars.parframe <- function(x, tol = 1, ...){
  
  if (!missing(...)) x <- subset(x, ...)
  
  jumps <- .stepDetect(x$value, tol)
  jump.index <- approx(jumps, jumps, xout = 1:length(x$value), method = "constant", rule = 2)$y
  
  x$index <- as.factor(jump.index)
  
  myparframe <- x
  parNames <- attr(myparframe,"parameters")
  parOut <- wide2long.data.frame(out = ((myparframe[, c("index", "value", parNames)])) , keep = 1:2)
  names(parOut) <- c("index", "value", "name", "parvalue")
  plot <- ggplot2::ggplot(parOut, aes(x = name, y = parvalue, color = index)) + geom_boxplot(outlier.alpha = 0) + theme_dMod() + scale_color_dMod() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  attr(plot, "data") <- parOut
  
  return(plot)
  
}


# 0 and the decades from 10 up to the largest absolute limit, both signs. On a
# pseudo-log axis, which is linear below 1, a break at 1 would crowd 0.
.decadeBreaks <- function(lims) {
  p <- 10^seq_len(max(1, ceiling(log10(max(abs(lims), 1)))))
  sort(unique(c(-p, 0, p)))
}


# Convergence always reads the same way, a filled circle for a converged fit and
# a triangle for one that stopped otherwise. Both levels stay in the legend even
# when only one occurs, so several waterfall plots can be compared directly.
.scaleConverged <- function()
  scale_shape_manual(name = "Converged",
                     values = c(`TRUE` = 16, `FALSE` = 17),
                     limits = c("TRUE", "FALSE"),
                     labels = c("yes", "no"), drop = FALSE)


#' @export
#' @rdname plotValues
plotValues.parframe <- function(x, tol = 1, ..., showSteps = FALSE) {

  if (!missing(...)) x <- subset(x, ...)

  pars <- x[order(x$value), ]
  pars[["index"]] <- seq_len(nrow(pars))
  pars[["delta"]] <- pars$value - pars$value[1]
  jumps <- .stepDetect(pars$value, tol)

  yScale <- scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                               breaks = .decadeBreaks)

  stepLines <- stepLabels <- NULL
  if (showSteps) {
    stepLines <- geom_vline(xintercept = jumps, lty = 2)
    stepLabels <- annotate("text", x = jumps + 1, y = pars$delta[jumps], label = jumps,
                           hjust = 0, vjust = -0.5, color = "firebrick", size = 3)
  }

  P <- ggplot2::ggplot(pars, aes(x = index, y = delta, pch = converged, color = iterations)) +
    stepLines +
    geom_point() +
    stepLabels +
    yScale +
    labs(x = "Fit rank", y = expression(Delta * " objective value"),
         color = "Iterations") +
    scale_color_gradient(low = "dodgerblue", high = "orange") +
    .scaleConverged() +
    theme_dMod()

  attr(P, "data") <- pars
  attr(P, "jumps") <- jumps

  return(P)

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
  delta <- 0
  if ("value" %in% colnames(parlist)) {
    delta   <- as.numeric(parlist[, "value", drop = TRUE] - .profileReference(sets))
    parlist <- parlist[, !colnames(parlist) %in%
                         c("index", "value", "converged", "iterations")]
  }
  points <- data.frame(par   = as.numeric(as.matrix(parlist)),
                       name  = rep(colnames(parlist), each = nrow(parlist)),
                       delta = delta)
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
    xlab("parameter value") + ylab("Confidence Level")

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
  data <- droplevels(subset(.profileFrame(sets, maxvalue), ...))
  .profilePlot(data, sets, parlist, ncol, maxvalue, threshold)
}


#' @export
#' @rdname plotProfile
plotProfile.list <- function(profs, ..., maxvalue = 5, parlist = NULL, ncol = NULL,
                             threshold = .profileLines()) {
  sets <- .profileSets(profs)
  data <- droplevels(subset(.profileFrame(sets, maxvalue), ...))
  .profilePlot(data, sets, parlist, ncol, maxvalue, threshold)
}



#' @export
#' @rdname parframe
is.parframe <- function(x) {
  "parframe" %in% class(x)
}

#' @export
#' @param i row index in any format. Default `NULL` keeps all rows.
#' @param j column index in any format. Default `NULL` keeps all columns.
#' @param drop logical. Default `FALSE`. `TRUE` returns a plain data frame or
#'   vector instead of a `parframe`.
#' @rdname parframe
"[.parframe" <- function(x, i = NULL, j = NULL, drop = FALSE){
  
  metanames <- attr(x, "metanames")
  obj.attributes <- attr(x, "obj.attributes")
  parameters <- attr(x, "parameters")
  
  out <- as.data.frame(x)
  if (!is.null(i)) out <- out[i, ]
  if (!is.null(j)) out <- out[, j, drop = drop]
  
  if (drop) return(out)
  
  metanames <- intersect(metanames, colnames(out))
  obj.attributes <- intersect(obj.attributes, colnames(out))
  parameters <- intersect(parameters, colnames(out))
  
  parframe(out, parameters = parameters, metanames = metanames, objAttributes = obj.attributes)
  
}


#' @export
#' @param ... for `subset()`, a logical expression in the columns of `x`.
#' @rdname parframe
subset.parframe <- function(x, ...) {

  # Columns of `x` first, then the caller's frame.
  condition <- eval(substitute(alist(...)))[[1L]]
  x[eval(condition, as.list(x), parent.frame()), ]

}

#' Rows of a Parameter Frame with Unique Values
#'
#' Keeps one row per distinct objective value.
#' @param x parameter frame
#' @param incomparables not used. Argument exists for compatibility with S3 generic.
#' @param tol tolerance to decide when values are assumed to be equal, see [plotValues()].
#' @param ... additional arguments being passed to [plotValues()], e.g. for subsetting.
#' @return A subset of the parameter frame `x`.
#' @export
unique.parframe <- function(x, incomparables = FALSE, tol = 1, ...) {
  
  
  jumps <- attr(plotValues(x = x, tol = tol, ...), "jumps")
  x[jumps, ]
  
  
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
  n_width <- max(nchar(nms))
  
  cat("Parameter vector:\n")
  for (i in seq_along(par)) {
    val <- formatC(par[i], digits = 6, format = "g")
    if (par[i] >= 0) val <- paste0(" ", val)
    cat(sprintf("  %s : %s\n", format(nms[i], width = n_width, justify = "right"), val))
  }
  
  deriv  <- attr(x, "deriv")
  deriv2 <- attr(x, "deriv2")
  fixed  <- attr(x, "fixed")

  cat("\nAttributes:\n")
  if (!is.null(deriv)) {
    d <- dim(deriv)
    cat(sprintf("  deriv  : %d x %d matrix\n", d[1], d[2]))
  } else {
    cat("  deriv  : <none>\n")
  }

  if (!is.null(deriv2)) {
    d2 <- dim(deriv2)
    cat(sprintf("  deriv2 : %d x %d x %d array\n", d2[1], d2[2], d2[3]))
  } else {
    cat("  deriv2 : <none>\n")
  }

  if (!is.null(fixed) && length(fixed) > 0) {
    cat(sprintf("  fixed  : %s\n", paste(fixed, collapse = ", ")))
  } else {
    cat("  fixed  : <none>\n")
  }

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


## Methods for the class parfn--------------------------------------------------

#' Print Parameter Transformations
#'
#' @param x object of class `parfn`.
#' @param ... not used.
#' @return `x`, invisibly.
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
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
    
    #print(output, ...)
    output
    
  } else {
    
    cat("\nObject is composed. See original objects for more details.\n")
    
  }
}




## parfn / parframe / parlist / parvec constructors (moved from classes.R) ----------------------------------------

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
#'   `colnames(x)`.
#' @param metanames character, the names of the meta-information columns.
#'   Default `NULL`.
#' @param objAttributes character, the names of the objective function
#'   attributes. Default `NULL`.
#' @param ... `obj.attributes` is deprecated, use `objAttributes`.
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
parframe <- function(x = NULL, parameters = colnames(x), metanames = NULL,
                     objAttributes = NULL, ...) {

  .renameArgs(list(...), c(obj.attributes = "objAttributes"), "parframe",
              strict = TRUE)

  if (!is.null(x)) {
    rownames(x) <- NULL
    out <- as.data.frame(x)
  } else {
    out <- data.frame()
  }

  attr(out, "parameters") <- parameters
  attr(out, "metanames") <- metanames
  attr(out, "obj.attributes") <- objAttributes
  class(out) <- c("parframe", "data.frame")

  return(out)

}

#' Parameter List
#'
#' @description A parameter list holds the results of several optimisation
#' runs, as returned by [mstrust()]. `c()` renumbers the fit index of the
#' combined list.
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
  
  mylist <- list(...)
  if (length(mylist) > 0) {
    mynames <- paste0("par", seq_along(mylist))
    is.available <- !is.null(names(mylist))
    mynames[is.available] <- names(mylist)[is.available]
    
    out <- as.numeric(unlist(mylist))
    names(out) <- mynames
    
    return(as.parvec(out, deriv = deriv))
  } else {
    return(NULL)
  }
}


