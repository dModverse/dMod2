test_that(".remoteLibs puts the paths in front of R_LIBS", {
  expect_identical(dMod2:::.remoteLibs(NULL), "")
  expect_identical(dMod2:::.remoteLibs(c("~/lib", "/opt/r")),
                   'export R_LIBS="$HOME/lib:/opt/r${R_LIBS:+:$R_LIBS}"; ')
})

test_that(".remoteBuildInfo drops local paths and flags the LAPACK solver", {
  env <- new.env()
  env$m <- structure(list(), compileInfo = list(list(
    srcfile = "m.cpp",
    compileArgs = "-DFOO -I/home/u/.cache/R/cppDE/include -fopenmp -I/usr/include/openmpi",
    linkArgs = "-lsundials_sunlinsollapackdense -L/home/u/lib -lsundials_cvodes")))
  info <- dMod2:::.remoteBuildInfo(env)
  expect_identical(info$compileArgs, "-DFOO -fopenmp")
  expect_true(info$needsCVODE)
  expect_true(info$needsLapack)
  expect_false(info$needsKLU)
})

test_that("the remote build script resolves SUNDIALS flags on the remote side", {
  script <- dMod2:::.remoteBuildScript("m.cpp", "m.so", compileArgs = "-DFOO",
                                       needsCVODE = TRUE, needsLapack = TRUE)
  expect_match(script, "\"cflags\"", fixed = TRUE)
  expect_match(script, "\"cvode_lapack_libs\", \"libs\"", fixed = TRUE)
  expect_match(script, "LAPACK dense solver", fixed = TRUE)
  expect_false(grepl("\"cflags\"", dMod2:::.remoteBuildScript("m.c", "m.so"),
                     fixed = TRUE))
})

test_that("the generated include line evaluates to the installed cppDE flags", {
  skip_on_os("windows")
  cfg <- get0("cvodeConfig", envir = asNamespace("cppDE"), inherits = FALSE)
  skip_if_not(is.environment(cfg) && isTRUE(cfg$available))
  script <- dMod2:::.remoteBuildScript("m.cpp", "m.so", needsCVODE = TRUE)
  line <- grep("^PKG_CPPFLAGS=\"\\$PKG_CPPFLAGS ", strsplit(script, "\n")[[1]],
               value = TRUE)
  expect_length(line, 1L)
  # The R this test runs in, ahead of the Rscript wrapper R CMD check puts on PATH.
  path <- paste0("PATH=", shQuote(R.home("bin")), ":$PATH; ")
  out <- system2("bash", c("-c", shQuote(paste0(path, "PKG_CPPFLAGS=; ", line,
                                                "; printf %s \"$PKG_CPPFLAGS\""))),
                 stdout = TRUE)
  expect_identical(trimws(out), trimws(cfg$cflags))
})

# ssh, scp and sbatch stand-ins on PATH for which every host is the directory
# `remote`, also as its home, so a job runs on this machine as it would on
# another. sbatch starts every array index at once.
local_fake_hosts <- function(env = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = env)
  bin <- file.path(root, "bin"); remote <- file.path(root, "remote")
  dir.create(bin); dir.create(remote); dir.create(file.path(root, "local"))
  writeLines(c("#!/bin/sh", 'while [ "${1#-}" != "$1" ]; do shift; done', "shift",
               'cd "$FAKE_REMOTE" && HOME="$FAKE_REMOTE" exec sh -c "$*"'),
             file.path(bin, "ssh"))
  writeLines(c("#!/bin/sh", "n=$(sed -n 's/^#SBATCH -a 0-//p' \"$1\")",
               'for i in $(seq 0 "$n"); do',
               '  SLURM_ARRAY_TASK_ID=$i SLURM_JOB_ID=$((i + 1)) bash "$1" > "${1%.sh}_$i.out" 2>&1 &',
               "done"),
             file.path(bin, "sbatch"))
  writeLines(c("#!/bin/sh", 'args=""',
               'for a in "$@"; do',
               '  case "$a" in -*) ;; *:*) args="$args $FAKE_REMOTE/${a#*:}" ;;',
               '  *) args="$args $a" ;; esac',
               "done", "exec cp -r $args"),
             file.path(bin, "scp"))
  Sys.chmod(file.path(bin, c("ssh", "scp", "sbatch")), "755")
  withr::local_envvar(PATH = paste(bin, R.home("bin"), Sys.getenv("PATH"), sep = ":"),
                      FAKE_REMOTE = remote, .local_envir = env)
  withr::local_dir(file.path(root, "local"), .local_envir = env)
  remote
}

test_that("runbg() runs a job on every machine and leaves nothing behind", {
  skip_on_os("windows")
  skip_on_cran()
  remote <- local_fake_hosts()
  assign("runbg_x", 21, envir = globalenv())
  withr::defer(rm("runbg_x", envir = globalenv()))

  out <- runbg({ 2 * runbg_x }, machine = c("h1", "h2"), filename = "jobw",
               input = "runbg_x", wait = TRUE)
  expect_identical(out, list(h1 = 42, h2 = 42))
  expect_length(list.files(remote), 0L)
  expect_length(list.files("."), 0L)
})

test_that("runbg() hands back check(), get() and purge(), also after recover", {
  skip_on_os("windows")
  skip_on_cran()
  remote <- local_fake_hosts()
  assign("runbg_x", 21, envir = globalenv())
  withr::defer(rm("runbg_x", envir = globalenv()))

  job <- runbg({ runbg_x + 1 }, machine = "h1", filename = "jobr", input = "runbg_x")
  for (i in 1:100) {
    if (isTRUE(suppressMessages(job$check()))) break
    Sys.sleep(0.2)
  }
  expect_true(suppressMessages(job$check()))
  expect_identical(job$get(), list(h1 = 22))

  again <- runbg(NULL, machine = "h1", filename = "jobr", recover = TRUE)
  expect_identical(again$get(), list(h1 = 22))
  again$purge()
  expect_length(list.files(remote), 0L)
  expect_length(list.files("."), 0L)
})

# Polls check() of a job until it reports every result, at most for `secs`.
wait_for <- function(job, secs = 60) {
  for (i in seq_len(secs * 5)) {
    if (isTRUE(job$check())) return(TRUE)
    Sys.sleep(0.2)
  }
  FALSE
}

test_that("distributedComputing() runs one array index per value and slots the results", {
  skip_on_os("windows")
  skip_on_cran()
  skip_if(!nzchar(Sys.which("zstd")), "zstd not available")
  remote <- local_fake_hosts()
  assign("dc_x", 10, envir = globalenv())
  withr::defer(rm("dc_x", envir = globalenv()))

  job <- distributedComputing({ dc_x + as.numeric(node_ID) + var_1 }, jobname = "dcv",
                              machine = "h1", varValues = list(c(1, 2, 3)),
                              input = "dc_x")
  expect_true(wait_for(job))
  expect_identical(job$get(), list(11, 13, 15))

  # A result that is missing leaves its own slot empty, not the last one.
  unlink(file.path(remote, "dcv_folder", "dcv_1_result.RData"))
  again <- distributedComputing(NULL, jobname = "dcv", machine = "h1",
                                varValues = list(c(1, 2, 3)), recover = TRUE)
  unlink("dcv_folder", recursive = TRUE)
  expect_identical(again$get(), list(11, NULL, 15))

  again$purge(purgeLocal = TRUE)
  expect_length(list.files(remote), 0L)
  expect_length(list.files("."), 0L)
})

test_that("distributedComputing() runs the code nRep times", {
  skip_on_os("windows")
  skip_on_cran()
  skip_if(!nzchar(Sys.which("zstd")), "zstd not available")
  local_fake_hosts()

  job <- distributedComputing({ as.numeric(node_ID) }, jobname = "dcn",
                              machine = "h1", nRep = 2, input = character(0))
  expect_true(wait_for(job))
  expect_identical(job$get(), list(0, 1))
  job$purge(purgeLocal = TRUE)
})
