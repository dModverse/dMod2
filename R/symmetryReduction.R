#' Remove Detected Symmetries by Reparametrisation
#'
#' Turns the non-identifiable directions of a [symmetryDetection()] result into a
#' parameter transformation that removes them. A scaling fixes one coordinate to 1.
#' A general direction is replaced by its invariants, held by new parameters
#' `q_<k>`, searched up to the degree bounds `dPoly`, `dDarboux` and `dExp`. The
#' computation runs in the Python package `symident`.
#'
#' @details A chart is returned only if it is valid on the domain declared by
#'   `positive`. A `q_<k>` that takes both signs is marked `[real-valued]` and is
#'   estimated on a linear scale; a chart valid only for its positive values has
#'   `coverage = "partial"`. `print()` lists the invariant behind each `q_<k>`.
#'
#' @param object A `symmetrydetection` result from [symmetryDetection()].
#' @param fixed Coordinates with known values; scalings they remove are dropped.
#'   `NULL` (default) uses the gauge of `object`.
#' @inheritParams symmetryDetection
#' @param reportZeroCompatibility Logical. Report which coordinates each block can
#'   drive to 0. Defaults to `FALSE`.
#' @param dPoly Degree bound of polynomial and rational invariants. Defaults
#'   to 3.
#' @param dDarboux Degree bound of the Darboux polynomials. Defaults to 2.
#' @param dExp Numerator degree bound of exponential factors; `0` skips them.
#'   Defaults to 2.
#' @param separable Logical. Solve separable blocks by quadratures. Defaults to
#'   `TRUE`.
#' @param alternatives Logical. Also return the charts of the other admissible
#'   pins in `$alternatives`. Defaults to `FALSE`.
#' @param timeout Time limit in seconds for the chart search of one block; a block not
#'   solved in time keeps its invariants (`"invariantOnly"`). Defaults to 600.
#' @param verbose Logical. Report the progress per block and stage. Defaults to
#'   `FALSE`.
#' @param ... Not used.
#'
#' @return An object of class `symmetryreduction`:
#'   \describe{
#'     \item{`blocks`}{one entry per set of coupled directions, with its
#'       `invariants`, `section`, `survivorMeaning` (the invariant each `q_<k>`
#'       holds), `carrierDomain`, `coverage` and `status` (`"fixed"`,
#'       `"reduced"`, `"invariantOnly"` or `"unresolved"`).}
#'     \item{`zeroCompatibility`}{with `reportZeroCompatibility = TRUE`, the sets of
#'       coordinates that can vanish together and the `condition` for it.}
#'     \item{`trafo`}{an [eqnvec] over all coordinates for [P()] or
#'       `symmetryDetection(trafo = )`; `NULL` if nothing was reduced.}
#'     \item{`family`}{the admissible gauges of each reduced block.}
#'     \item{`removed`, `remaining`}{labels of the directions.}
#'     \item{`dependent`}{labels of directions that depend on the others and need
#'       no reduction.}
#'     \item{`partial`}{labels of removed directions whose chart covers only part of
#'       the positive orthant.}
#'     \item{`alternatives`}{the other charts, each with `trafo`, `removed`,
#'       `remaining`, `rational` and `pins`.}
#'     \item{`coordinates`, `fixed`, `settings`, `call`}{provenance.}
#'   }
#'   `print()` shows the verdict, the trafo and the new parameters, `summary()` one
#'   line per block. Both take `width`, `summary()` also `verbose`.
#'
#' @seealso [symmetryDetection()], [P()]
#' @examplesIf requireNamespace("reticulate", quietly = TRUE) && reticulate::py_module_available("symident")
#' f <- eqnvec(m = "ktx - dm*m", p = "ktl*m - dp*p")
#' g <- eqnvec(y = "p")
#' res <- symmetryDetection(f, g, reconstruct = TRUE, verbose = FALSE)
#' red <- symmetryReduction(res)
#' red
#' red$trafo
#' symmetryDetection(f, g, trafo = red$trafo, verbose = FALSE)$identifiable
#' @export
symmetryReduction <- function(object, fixed = NULL, positive = TRUE, dPoly = 3L,
                              dDarboux = 2L, dExp = 2L, separable = TRUE,
                              reportZeroCompatibility = FALSE, alternatives = FALSE,
                              timeout = 600, verbose = FALSE, ...) {
  if (!inherits(object, "symmetrydetection"))
    stop("symmetryReduction(): `object` must be a symmetrydetection result.", call. = FALSE)
  if (length(list(...)))
    warning("symmetryReduction(): unused argument(s) ignored.", call. = FALSE)
  chr <- function(x) if (is.null(x)) NULL else as.list(as.character(x))
  args <- list(fixed = chr(fixed),
               positive = if (is.logical(positive)) isTRUE(positive) else chr(positive),
               d_poly = as.integer(dPoly), d_darboux = as.integer(dDarboux),
               d_exp = as.integer(dExp), separable = isTRUE(separable),
               zero_compatibility = isTRUE(reportZeroCompatibility),
               alternatives = isTRUE(alternatives),
               timeout = if (is.infinite(timeout)) "Inf" else timeout,
               verbose = isTRUE(verbose))
  if (is.null(fixed)) args$fixed <- NULL
  .symRedFromPy(.symCall("reduce_json", .symRaw(object), args), match.call())
}

# the result of symident as a symmetryreduction object, the raw result attached
.symRedFromPy <- function(res, call) {
  raw <- attr(res, "json")
  res <- .symCamelFields(unclass(res))
  attr(res, "json") <- NULL
  chr <- function(x) if (is.null(x)) NULL else as.character(unlist(x))
  str <- function(x) if (is.null(x)) NULL else setNames(as.character(unlist(x)), names(x))
  zero <- function(z, block = FALSE) if (length(z))
    do.call(rbind, lapply(z, function(r) {
      d <- data.frame(block = if (block) r$block else NA, coordinates = r$coordinates,
                      verdict = r$verdict, limit = if (is.null(r$limit)) NA else r$limit,
                      certain = r$certain, condition = r$condition, at = r$at,
                      stringsAsFactors = FALSE)
      if (block) d else d[-1L]
    }))
  res$blocks <- lapply(res$blocks, function(b) {
    for (k in c("labels", "support", "transversal", "invariants", "certificates",
                "section", "moduleCombos", "redundantFixed")) if (!is.null(b[[k]])) b[[k]] <- chr(b[[k]])
    for (k in c("pins", "survivorMeaning", "carrierDomain")) if (!is.null(b[[k]])) b[[k]] <- str(b[[k]])
    if (!is.null(b$zeroCompatibility)) b$zeroCompatibility <- zero(b$zeroCompatibility)
    b
  })
  res$family <- lapply(res$family, function(f) {
    f$labels <- chr(f$labels)
    if (!is.null(f$admissible)) f$admissible <- lapply(f$admissible, chr)
    if (!is.null(f$matroid)) f$matroid <- lapply(f$matroid, chr)
    f
  })
  if (!is.null(res$trafo)) res$trafo <- as.eqnvec(str(res$trafo))
  for (k in c("removed", "remaining", "dependent", "coordinates", "fixed"))
    res[k] <- list(chr(res[[k]]))
  res$alternatives <- lapply(res$alternatives, function(a) {
    if (!is.null(a$trafo)) a$trafo <- as.eqnvec(str(a$trafo))
    for (k in c("removed", "remaining")) a[k] <- list(chr(a[[k]]))
    a$pins <- str(a$pins)
    a
  })
  res$partial <- as.character(unlist(lapply(res$blocks, function(b)
    if (identical(b$coverage, "partial") && identical(b$status, "reduced")) b$labels)))
  if (!is.null(res$zeroCompatibility)) res$zeroCompatibility <- zero(res$zeroCompatibility, TRUE)
  if (is.list(res$settings$positive)) res$settings$positive <- chr(res$settings$positive)
  res$call <- call
  structure(res, class = "symmetryreduction", symident = raw)
}

#' @export
print.symmetryreduction <- function(x, width = getOption("width"), ...) {
  .symReport("reduction_lines", x, width = as.integer(width))
  invisible(x)
}

#' @export
summary.symmetryreduction <- function(object, verbose = FALSE, width = getOption("width"), ...)
  structure(list(object = object, verbose = isTRUE(verbose), width = as.integer(width)),
            class = "summary.symmetryreduction")

#' @export
print.summary.symmetryreduction <- function(x, ...) {
  .symReport("reduction_summary_lines", x$object, verbose = x$verbose, width = x$width)
  invisible(x)
}
