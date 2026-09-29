#' Install external programs and libraries
#'
#' Builds msolve, required by [symmetryDetection()] with `equilibrate = TRUE`,
#' into a per-user cache. `"sundials"` and `"suitesparse"` are passed to
#' [cppDE::install_libs()].
#'
#' msolve is built with GMP, MPFR and FLINT, statically, into
#' `tools::R_user_dir("dMod2", "cache")`. Requires a C compiler, `make`, `m4`
#' and `curl` or `wget`. Not available on Windows. Re-install dMod2 afterwards.
#'
#' @param which `"msolve"`, `"sundials"` or `"suitesparse"`.
#' @param dir Install directory for msolve.
#' @param msolve_version msolve release.
#' @param jobs Parallel make jobs.
#' @param quiet If `TRUE`, suppress the build output.
#' @param ask If `TRUE`, ask before downloading.
#' @param ... Passed to [cppDE::install_libs()].
#'
#' @return The install prefix, invisibly.
#' @export
install_libs <- function(which = c("msolve", "sundials", "suitesparse"), dir = NULL,
                         msolve_version = "0.10.1", jobs = 2L, quiet = FALSE,
                         ask = interactive(), ...) {
  which <- match.arg(which)
  if (which != "msolve")
    return(cppDE::install_libs(which = which, quiet = quiet, ask = ask, ...))

  if (.Platform$OS.type == "windows")
    stop("install_libs(\"msolve\") is not supported on Windows.", call. = FALSE)
  if (is.null(dir)) dir <- tools::R_user_dir("dMod2", "cache")
  if (all(Sys.which(c("curl", "wget")) == ""))
    stop("need curl or wget to download the sources.", call. = FALSE)
  if (any(Sys.which(c("make", "m4")) == ""))
    stop("make and m4 are required for the build.", call. = FALSE)
  script <- system.file("tools", "build-msolve.sh", package = "dMod2")
  if (!nzchar(script)) script <- file.path("inst", "tools", "build-msolve.sh")

  if (isTRUE(ask)) {
    msg <- paste0("This downloads and builds msolve ", msolve_version,
                  " with GMP, MPFR and FLINT\ninto: ", dir,
                  "\nIt requires network access and takes a few minutes. Proceed?")
    if (!isTRUE(utils::askYesNo(msg, default = TRUE))) {
      message("Aborted; nothing was downloaded or written.")
      return(invisible(NULL))
    }
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  cc <- trimws(system2(file.path(R.home("bin"), "R"), c("CMD", "config", "CC"),
                       stdout = TRUE))
  out <- system2("sh", shQuote(normalizePath(script)),
                 env = c(paste0("DMOD2_LIBS_CACHE=", shQuote(normalizePath(dir))),
                         paste0("DMOD2_MSOLVE_VERSION=", shQuote(msolve_version)),
                         paste0("DMOD2_JOBS=", as.integer(jobs)),
                         paste0("CC=", shQuote(cc))),
                 stdout = TRUE, stderr = if (quiet) FALSE else "")
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0L)
    stop("the msolve build failed; see the messages above.", call. = FALSE)
  prefix <- utils::tail(out[nzchar(out)], 1L)
  if (!file.exists(file.path(prefix, "bin", "msolve")))
    stop("the build reported success but produced no msolve binary.", call. = FALSE)
  message("msolve installed in ", file.path(prefix, "bin", "msolve"),
          ". Re-install dMod2 to use it.")
  invisible(prefix)
}


# msolve recorded by configure, "" if none
.msolvePath <- function() {
  dcf <- system.file("msolveConfig.dcf", package = "dMod2")
  if (!nzchar(dcf)) return("")
  cfg <- read.dcf(dcf, fields = c("available", "path"))
  if (!isTRUE(toupper(cfg[1, "available"]) == "TRUE")) return("")
  path <- unname(cfg[1, "path"])
  if (is.na(path) || !file.exists(path)) return("")
  path
}
