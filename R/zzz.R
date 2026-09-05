#' Package initialization
#'
#' Declares the Python dependencies the symbolic layer needs, `getLinVars()`
#' and the `simplify` path in `eqnClass.R`, and `steadyStates()`.
#' `reticulate::py_require()` does not start Python here: it records the
#' requirement so the first call into Python provisions an env that has the
#' listed packages. The SBML/PEtab layer adds `python-libsbml` to this list.
#'
#' @keywords internal
#' @importFrom stats setNames predict density median dnorm
#' @importFrom utils globalVariables
#' @noRd
.onLoad <- function(libname, pkgname) {
  if (requireNamespace("reticulate", quietly = TRUE))
    reticulate::py_require(c("python-libsbml", "sympy", "scipy", "numpy", "symengine"))
}

#' Package attach
#'
#' @keywords internal
#' @noRd
.onAttach <- function(libname, pkgname) {
  .announceForkGuard()
}

# mstrust() and profile() fork, and a threaded BLAS deadlocks in the child. cppDE
# pins it for the width of each fork, which makes BLAS inside a worker serial, so
# say so here, and say when there is no lever and the hang is still reachable.
.announceForkGuard <- function() {
  if (isTRUE(getOption("dMod.quiet"))) return(invisible(NULL))
  g <- tryCatch(cppDE::forkGuard(), error = function(e) NULL)
  if (is.null(g) || !isTRUE(g$guard)) return(invisible(NULL))

  if (is.na(g$api)) {
    packageStartupMessage(
      "dMod2: no BLAS thread-control entry point found. If this BLAS is threaded, ",
      "mstrust() and profile() with cores > 1 can deadlock; start R with ",
      "OMP_NUM_THREADS=1. See cppDE::forkGuard().")
  } else if (!is.na(g$threads) && g$threads > 1L) {
    packageStartupMessage(sprintf(
      "dMod2: %s (%d threads) is pinned to 1 inside the forked workers of mstrust() and profile(). See cppDE::forkGuard().",
      g$api, g$threads))
  }
  invisible(NULL)
}

# Guard for optionally-Suggested feature dependencies: fail with an actionable
# message if the package a feature needs is not installed. Lets heavy/optional
# deps (Python bridge, SBML/PEtab I/O, LP solver) live in Suggests, not Imports.
.require_ns <- function(pkg, feature) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' is required for ", feature,
         "; install it with install.packages(\"", pkg, "\").", call. = FALSE)
  invisible(TRUE)
}

# cppDE's install-time SUNDIALS/KLU flags, read from its (unexported)
# cvodeConfig environment. All-empty when unavailable.
.cppDE_config <- function() {
  empty <- list(available = FALSE, cflags = "", libs = "",
                klu_available = FALSE, klu_cflags = "", klu_libs = "")
  if (!requireNamespace("cppDE", quietly = TRUE)) return(empty)
  cfg <- get0("cvodeConfig", envir = asNamespace("cppDE"), inherits = FALSE)
  if (!is.environment(cfg)) return(empty)
  utils::modifyList(empty, as.list(cfg)[names(empty)][
    !vapply(as.list(cfg)[names(empty)], is.null, logical(1))])
}
