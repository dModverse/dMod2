# Run by testthat before any test file. Snapshot cwd and restore it at the
# end of the run so wayward setwd() inside individual tests cannot break
# testthat's relative-path lookups for the rest of the session.
.dmod_initial_wd <- getwd()

if (requireNamespace("withr", quietly = TRUE)) {
  withr::defer(setwd(.dmod_initial_wd), testthat::teardown_env())
}

# Fixture directories outside the installed package (PEtabTests/ is
# .Rbuildignore'd) are resolved into env vars once, before any test changes wd.
# Walks up at most 8 levels from the initial wd looking for the named dir.
.dmod_find_fixture <- function(name, start = .dmod_initial_wd) {
  here <- normalizePath(start, mustWork = FALSE, winslash = "/")
  for (i in seq_len(8)) {
    cand <- file.path(here, name)
    if (dir.exists(cand)) return(normalizePath(cand, winslash = "/"))
    parent <- dirname(here)
    if (parent == here) break
    here <- parent
  }
  ""
}

local({
  for (spec in list(
    c(env = "DMOD_PETABTESTS", dir = "PEtabTests")
  )) {
    if (!nzchar(Sys.getenv(spec[["env"]], unset = ""))) {
      hit <- .dmod_find_fixture(spec[["dir"]])
      if (nzchar(hit)) do.call(Sys.setenv, setNames(list(hit), spec[["env"]]))
    }
  }
})
