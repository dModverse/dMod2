# msolve lookup of configure (inst/tools/msolve-config.R)

msolve_config <- function() {
  env <- new.env()
  sys.source(system.file("tools", "msolve-config.R", package = "dMod2"), envir = env)
  env
}

# fake msolve: answers -V and writes the -P 1 output of the test system
fake_msolve <- function(dir, version, good = TRUE) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(dir, "msolve")
  out <- if (good)
    "[0, [65521, \\n2, \\n3, \\n['x', 'y'],\\n[0, 1],\\n[1,\\n[[3,\\n[36, 0, 65514, 1]],\\n[0,\\n[1]],\\n[\\n[[2,\\n[0, 10919, 54601]]]\\n]]]]]:\\n"
  else "[1, 2, -1, []]:\\n"
  writeLines(c("#!/bin/sh",
               sprintf("if [ \"$1\" = \"-V\" ]; then echo %s; exit 0; fi", version),
               "while [ $# -gt 0 ]; do [ \"$1\" = \"-o\" ] && out=$2; shift; done",
               sprintf("printf \"%s\" > \"$out\"", out)), f)
  Sys.chmod(f, "755")
  normalizePath(f)
}

test_that("configure takes DMOD2_MSOLVE, then a complete cache build, then a tested PATH msolve", {
  skip_on_os("windows")
  cfg <- msolve_config()
  root <- withr::local_tempdir()
  cache <- file.path(root, "cache")
  onPath <- file.path(root, "bin")
  withr::local_envvar(DMOD2_MSOLVE = "", DMOD2_BUILD_MSOLVE = "",
                      PATH = paste(onPath, Sys.getenv("PATH"), sep = .Platform$path.sep))
  resolve <- function() cfg$msolveResolve(cache = cache, log = function(...) NULL)

  expect_false(resolve()$available)

  # an untested version on the PATH is refused, a tested one taken
  fake_msolve(onPath, "0.9.0")
  expect_false(resolve()$available)
  sys <- fake_msolve(onPath, "0.10.3")
  r <- resolve()
  expect_true(r$available); expect_identical(r$path, sys); expect_identical(r$source, "system")

  # a cache build counts only with its completion marker, the newest first
  old <- fake_msolve(file.path(cache, "deps-msolve-0.10.1", "bin"), "0.10.1")
  expect_identical(resolve()$source, "system")
  file.create(file.path(cache, "deps-msolve-0.10.1", ".msolve-complete"))
  expect_identical(resolve()$path, old)
  new <- fake_msolve(file.path(cache, "deps-msolve-0.10.2", "bin"), "0.10.2")
  file.create(file.path(cache, "deps-msolve-0.10.2", ".msolve-complete"))
  expect_identical(resolve()$path, new)

  # a build whose output fails the check is passed over
  fake_msolve(file.path(cache, "deps-msolve-0.10.2", "bin"), "0.10.2", good = FALSE)
  expect_identical(resolve()$path, old)

  # DMOD2_MSOLVE comes first, in any version that passes the check
  own <- fake_msolve(file.path(root, "own"), "0.11.0")
  withr::local_envvar(DMOD2_MSOLVE = own)
  r <- resolve()
  expect_identical(r$path, own); expect_identical(r$source, "DMOD2_MSOLVE")
})
