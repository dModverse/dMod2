# Locates msolve and writes inst/msolveConfig.dcf; called by configure.
# Order: DMOD2_MSOLVE, a build with DMOD2_BUILD_MSOLVE=1, the per-user cache, the PATH.
# Each candidate must solve a test system. Usage: Rscript msolve-config.R <out.dcf> [<builder>]

.msolveTested <- "^0\\.10\\."

msolveCheck <- function(bin) {
  if (!nzchar(bin) || !file.exists(bin)) return(NULL)
  td <- tempfile("msolvecheck"); dir.create(td); on.exit(unlink(td, recursive = TRUE))
  fi <- file.path(td, "in.ms"); fo <- file.path(td, "out.ms")
  writeLines(c("x, y", "65521", "x^2+y-7,", "x*y-6"), fi)
  st <- suppressWarnings(system2(bin, c("-f", fi, "-o", fo, "-P", "1", "-t", "1"),
                                 stdout = FALSE, stderr = FALSE))
  if (!identical(as.integer(st), 0L) || !file.exists(fo)) return(NULL)
  out <- paste(readLines(fo, warn = FALSE), collapse = "")
  if (!grepl("^\\[0, \\[65521,\\s*2,\\s*3,", out) || !grepl("\\[36, 0, 65514, 1\\]", out))
    return(NULL)
  ver <- tryCatch(suppressWarnings(system2(bin, "-V", stdout = TRUE, stderr = FALSE))[1],
                  error = function(e) NA_character_)
  list(path = normalizePath(bin), version = if (is.na(ver)) "" else trimws(ver))
}

msolveCacheBins <- function(cache) {
  dirs <- Sys.glob(file.path(cache, "deps-msolve-*"))
  dirs <- dirs[file.exists(file.path(dirs, ".msolve-complete"))]
  if (!length(dirs)) return(character(0))
  ver <- numeric_version(sub("^deps-msolve-", "", basename(dirs)), strict = FALSE)
  file.path(dirs[order(ver, decreasing = TRUE)], "bin", "msolve")
}

msolveResolve <- function(cache = tools::R_user_dir("dMod2", "cache"), builder = "",
                          log = function(...) message("msolve: ", ...)) {
  found <- function(chk, source) c(chk, list(available = TRUE, source = source))
  explicit <- Sys.getenv("DMOD2_MSOLVE")
  if (nzchar(explicit)) {
    chk <- msolveCheck(explicit)
    if (!is.null(chk)) return(found(chk, "DMOD2_MSOLVE"))
    log("DMOD2_MSOLVE=", explicit, " fails the check; ignored")
  }
  if (identical(Sys.getenv("DMOD2_BUILD_MSOLVE"), "1") && nzchar(builder)) {
    log("building from source into ", cache)
    dir.create(cache, recursive = TRUE, showWarnings = FALSE)
    cc <- trimws(system2(file.path(R.home("bin"), "R"), c("CMD", "config", "CC"),
                         stdout = TRUE))
    out <- system2("sh", shQuote(builder), stdout = TRUE,
                   env = c(paste0("DMOD2_LIBS_CACHE=", shQuote(cache)), paste0("CC=", shQuote(cc)),
                           paste0("DMOD2_JOBS=", Sys.getenv("DMOD2_JOBS", "2"))))
    if (!is.null(attr(out, "status"))) log("the source build failed; searching further")
  }
  for (bin in msolveCacheBins(cache)) {
    chk <- msolveCheck(bin)
    if (!is.null(chk)) return(found(chk, "per-user cache"))
  }
  for (bin in unique(c(Sys.which("msolve"), path.expand("~/.local/bin/msolve")))) {
    chk <- msolveCheck(bin)
    if (is.null(chk)) next
    if (grepl(.msolveTested, chk$version)) return(found(chk, "system"))
    log(bin, " is version ", chk$version, ", not one dMod2 is tested with; ignored")
  }
  list(available = FALSE, path = "", version = "", source = "")
}

if (!interactive() && length(commandArgs(TRUE)) >= 1L && sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  res <- msolveResolve(builder = if (length(args) >= 2L) args[2] else "")
  write.dcf(data.frame(available = res$available, path = res$path, version = res$version,
                       source = res$source), args[1])
  if (isTRUE(res$available))
    message(sprintf("msolve %s (%s): %s", res$version, res$source, res$path))
  else
    message("msolve: not found; symmetryDetection(equilibrate = TRUE) needs it for ",
            "coupled steady states, see dMod2::install_libs(\"msolve\")")
}
