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

# ssh, sbatch, squeue and scancel stand-ins on PATH for which every host is the
# directory `remote`, also as its home, so a job runs on the local machine as it would
# on another. sbatch starts every array task at once.
local_fake_hosts <- function(env = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = env)
  bin <- file.path(root, "bin"); remote <- file.path(root, "remote")
  dir.create(bin); dir.create(remote); dir.create(file.path(root, "local"))
  writeLines(c("#!/bin/sh", 'while [ "${1#-}" != "$1" ]; do shift; done', "shift",
               'cd "$FAKE_REMOTE" && HOME="$FAKE_REMOTE" exec sh -c "$*"'),
             file.path(bin, "ssh"))
  writeLines(c("#!/bin/sh", "n=$(sed -n 's/^#SBATCH --array=1-//p' \"$1\")",
               'for i in $(seq 1 "$n"); do',
               '  SLURM_ARRAY_TASK_ID=$i bash "$1" > /dev/null 2>&1 &',
               "done", "echo Submitted batch job 4242"),
             file.path(bin, "sbatch"))
  writeLines(c("#!/bin/sh", "n=$(sed -n 's/^#SBATCH --array=1-//p' job.sh)",
               'for i in $(seq 1 "$n"); do [ -f "task_$i.rds" ] || echo "running $i"; done'),
             file.path(bin, "squeue"))
  writeLines(c("#!/bin/sh", "true"), file.path(bin, "scancel"))
  Sys.chmod(file.path(bin, c("ssh", "sbatch", "squeue", "scancel")), "755")
  withr::local_envvar(PATH = paste(bin, R.home("bin"), Sys.getenv("PATH"), sep = ":"),
                      FAKE_REMOTE = remote, R_LIBS = paste(.libPaths(), collapse = ":"),
                      .local_envir = env)
  withr::local_dir(file.path(root, "local"), .local_envir = env)
  remote
}

# An ssh stand-in that fails, and an empty working directory.
local_no_ssh <- function(env = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = env)
  bin <- file.path(root, "bin")
  dir.create(bin); dir.create(file.path(root, "wd"))
  writeLines(c("#!/bin/sh", "exit 99"), file.path(bin, "ssh"))
  Sys.chmod(file.path(bin, "ssh"), "755")
  withr::local_envvar(PATH = paste(bin, Sys.getenv("PATH"), sep = ":"), .local_envir = env)
  withr::local_dir(file.path(root, "wd"), .local_envir = env)
}

# Polls check() of a job until no task runs any more, at most for `secs`.
wait_for <- function(job, secs = 60) {
  for (i in seq_len(secs * 5)) {
    utils::capture.output(ok <- job$check())
    if (ok) return(TRUE)
    Sys.sleep(0.2)
  }
  FALSE
}

test_that("tasks run on localhost without ssh, with task, taskId and the calling environment", {
  skip_on_os("windows")
  skip_on_cran()
  local_no_ssh()
  f <- function() {
    x0 <- 10
    distributedComputing(x0 * task + taskId, tasks = list(a = 1, b = 2), input = "x0",
                         jobname = "jl", wait = TRUE)
  }
  expect_identical(f(), list(a = 11, b = 22))
  expect_length(list.files("."), 0L)
  expect_length(detectFreeCores("localhost"), 1L)
})

test_that("data frame rows are tasks, a task reads its cores and seeds repeat", {
  skip_on_os("windows")
  skip_on_cran()
  local_no_ssh()
  local_mocked_bindings(detectFreeCores = function(machine) setNames(rep(8, length(machine)),
                                                                     machine))
  tasks <- data.frame(a = 1:3, b = 3:1)
  run <- function() distributedComputing(c(task$a + task$b, cores, runif(1)), tasks = tasks,
                                         machine = c(localhost = 3), seed = 7L,
                                         input = character(0), jobname = "js", wait = TRUE)
  one <- run()
  expect_length(one, 3L)
  expect_identical(vapply(one, `[`, 0, 1L), rep(4, 3))
  expect_identical(vapply(one, `[`, 0, 2L), rep(3, 3))
  expect_identical(run(), one)
  expect_length(unique(vapply(one, `[`, 0, 3L)), 3L)
})

test_that("remote machines run their tasks, and recover and purge reach them", {
  skip_on_os("windows")
  skip_on_cran()
  remote <- local_fake_hosts()
  y <- 3
  job <- distributedComputing(y * taskId, tasks = 3, machine = c("h1", "h2"), input = "y",
                              jobname = "jr")
  expect_identical(attr(job, "manifest")$taskHost, c(1L, 2L, 1L))
  expect_s3_class(job, "dModJob")
  expect_true(wait_for(job))
  expect_identical(job$get(), list(3, 6, 9))
  again <- distributedComputing(NULL, jobname = "jr", recover = TRUE)
  expect_identical(again$get(), list(3, 6, 9))
  again$purge(local = TRUE)
  expect_length(list.files(remote), 0L)
  expect_length(list.files("."), 0L)
})

test_that("the model files come from dMod2.outdir", {
  skip_on_os("windows")
  skip_on_cran()
  local_fake_hosts()
  models <- withr::local_tempdir()
  writeLines("void dmod_probe(void) {}", file.path(models, "probe.c"))
  withr::with_dir(models, system2(file.path(R.home("bin"), "R"), c("CMD", "SHLIB", "probe.c"),
                                  stdout = FALSE, stderr = FALSE))
  withr::local_options(dMod2.outdir = models)
  expect_identical(distributedComputing(file.exists("probe.so"), machine = "h1",
                                        input = character(0), jobname = "jo", wait = TRUE),
                   list(TRUE))
  expect_length(list.files("."), 0L)
})

test_that("a SLURM job fills the template of its login node and runs one array task per task", {
  skip_on_os("windows")
  skip_on_cran()
  remote <- local_fake_hosts()
  withr::local_envvar(R_USER_CONFIG_DIR = withr::local_tempdir())
  tpl <- withr::local_tempfile(lines = c(
    "#!/bin/bash", "#SBATCH --job-name={{ jobname }}", "#SBATCH --array=1-{{ ntasks }}",
    "#SBATCH --time={{walltime}}", "#SBATCH --cpus-per-task={{ cores }}",
    "#SBATCH --mem-per-cpu={{ memPerCore }}gb", "export FROM_TEMPLATE=1", "{{ run }}"))
  slurmTemplate("h1", from = tpl, edit = FALSE)
  job <- distributedComputing(c(task$p, Sys.getenv("FROM_TEMPLATE"), cores),
                              tasks = data.frame(p = c("a", "b")), machine = c(h1 = 2),
                              slurm = TRUE, walltime = "00:10:00",
                              jobname = "jsl", input = character(0))
  expect_true(wait_for(job))
  expect_identical(job$get(), list(c("a", "1", "2"), c("b", "1", "2")))
  script <- readLines(file.path(remote, "jsl_folder", "1", "job.sh"))
  expect_true(all(c("#SBATCH --array=1-2", "#SBATCH --time=00:10:00",
                    "#SBATCH --cpus-per-task=2") %in% script))
  job$purge(local = TRUE)
  expect_length(list.files(remote), 0L)
})

test_that("a machine gets at most its free cores, and every machine is named once", {
  skip_on_os("windows")
  skip_on_cran()
  local_no_ssh()
  local_mocked_bindings(detectFreeCores = function(machine) setNames(rep(2, length(machine)),
                                                                     machine))
  expect_warning(out <- distributedComputing(cores, machine = c(localhost = 5),
                                             input = character(0), jobname = "jc", wait = TRUE),
                 "2 free cores; its tasks get 2 instead of 5")
  expect_identical(out, list(2L))
  expect_identical(distributedComputing(cores, input = character(0), jobname = "jf",
                                        wait = TRUE), list(2L))
  expect_identical(dMod2:::.machineCores(c(b = 3, a = 2)), c(b = 3L, a = 2L))
  expect_identical(dMod2:::.machineCores(c("a", "b")), c(a = NA_integer_, b = NA_integer_))
  expect_error(dMod2:::.machineCores(c(a = 1.5)), "whole numbers")
  expect_error(distributedComputing(1, machine = c("a", "a")), "every machine once")
  expect_error(distributedComputing(1, machine = 3), "name the machines")
})

test_that("slurmTemplate() sets up one template per login node and keeps it", {
  withr::local_envvar(R_USER_CONFIG_DIR = withr::local_tempdir())
  expect_length(slurmTemplate(), 0L)
  helix <- slurmTemplate("helix", edit = FALSE)
  expect_identical(readLines(helix),
                   readLines(system.file("slurm", "helix.sh", package = "dMod2")))
  own <- withr::local_tempfile(lines = c("#SBATCH --array=1-{{ ntasks }}", "{{ run }}"))
  mine <- slurmTemplate("mycluster", from = own, edit = FALSE)
  expect_identical(names(slurmTemplate()), c("helix", "mycluster"))
  writeLines("edited", mine)
  slurmTemplate("mycluster", from = "helix", edit = FALSE)
  expect_identical(readLines(mine), "edited")
  slurmTemplate("mycluster", from = "helix", edit = FALSE, overwrite = TRUE)
  expect_identical(readLines(mine), readLines(helix))
  expect_error(slurmTemplate("other", edit = FALSE), "no template for 'other'; give the job")
  expect_error(slurmTemplate("other", from = "none", edit = FALSE), "'none' is neither")
  expect_error(slurmTemplate("a/b", edit = FALSE), "one login node")
})

test_that("templates are checked, and a SLURM job needs one", {
  withr::local_envvar(R_USER_CONFIG_DIR = withr::local_tempdir())
  tpl <- withr::local_tempfile(lines = c("#SBATCH --time={{ walltime }}", "{{ run }}"))
  expect_error(dMod2:::.fillJobscript(tpl, list(jobname = "j")), "\\{\\{ walltime \\}\\}")
  norun <- withr::local_tempfile(lines = "#SBATCH --job-name={{ jobname }}")
  expect_error(dMod2:::.fillJobscript(norun, list(jobname = "j")), "\\{\\{ run \\}\\}")
  expect_error(distributedComputing(1, machine = c(h = 2), slurm = TRUE),
               "no job script for 'h'; set one up once with slurmTemplate\\(\"h\"")
  expect_error(distributedComputing(1, machine = c(h = 2), slurm = TRUE, jobscript = "nope.sh"),
               "nope.sh does not exist")
  expect_error(distributedComputing(1, machine = "h", slurm = TRUE), "cores of each array task")
  expect_error(distributedComputing(1, tasks = 0), "positive whole number")
})

test_that("a failing task returns its error, and terminate() stops running tasks", {
  skip_on_os("windows")
  skip_on_cran()
  local_no_ssh()
  out <- distributedComputing(stop("no"), input = character(0), jobname = "je", wait = TRUE)
  expect_s3_class(out[[1]], "try-error")

  job <- distributedComputing(Sys.sleep(60), input = character(0), jobname = "jt")
  Sys.sleep(1)
  utils::capture.output(running <- dMod2:::.taskStatus(attr(job, "manifest")))
  expect_identical(running, "running")
  job$terminate()
  for (i in 1:50) if (dMod2:::.taskStatus(attr(job, "manifest")) != "running") break else
    Sys.sleep(0.1)
  expect_identical(dMod2:::.taskStatus(attr(job, "manifest")), "failed")
  expect_warning(utils::capture.output(finished <- job$check()), "1 of 1 tasks failed")
  expect_true(finished)
  job$purge(local = TRUE)
  expect_length(list.files("."), 0L)
})
