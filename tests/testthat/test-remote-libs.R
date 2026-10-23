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
# on another. sbatch starts every array task at once and keeps its options in
# sbatch.args; squeue lists the array tasks still running. Jobs register in a
# temporary data directory.
local_fake_hosts <- function(env = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = env)
  bin <- file.path(root, "bin"); remote <- file.path(root, "remote")
  dir.create(bin); dir.create(remote); dir.create(file.path(root, "local"))
  writeLines(c("#!/bin/sh", 'while [ "${1#-}" != "$1" ]; do shift; done', "shift",
               'cd "$FAKE_REMOTE" && HOME="$FAKE_REMOTE" exec sh -c "$*"'),
             file.path(bin, "ssh"))
  writeLines(c("#!/bin/sh", 'for a; do s=$a; done', 'echo "$@" > sbatch.args',
               "n=$(sed -n 's/^#SBATCH --array=1-//p' \"$s\")",
               'for i in $(seq 1 "$n"); do',
               '  SLURM_ARRAY_TASK_ID=$i bash "$s" > /dev/null 2>&1 &',
               '  echo $! > ".array_$i"',
               "done", "echo sbatch: job accepted", "echo Submitted batch job 4242"),
             file.path(bin, "sbatch"))
  writeLines(c("#!/bin/sh",
               'for f in .array_*; do [ -f "$f" ] || continue; p=$(cat "$f")',
               '  if [ -r /proc/$p/status ]; then',
               '    grep -q "^State:[[:space:]]*Z" /proc/$p/status || echo "running ${f#.array_}"',
               '  else kill -0 "$p" 2>/dev/null && echo "running ${f#.array_}"; fi',
               "done"),
             file.path(bin, "squeue"))
  writeLines(c("#!/bin/sh", "true"), file.path(bin, "scancel"))
  Sys.chmod(file.path(bin, c("ssh", "sbatch", "squeue", "scancel")), "755")
  withr::local_envvar(PATH = paste(bin, R.home("bin"), Sys.getenv("PATH"), sep = ":"),
                      FAKE_REMOTE = remote, R_LIBS = paste(.libPaths(), collapse = ":"),
                      R_USER_DATA_DIR = file.path(root, "data"), .local_envir = env)
  withr::local_dir(file.path(root, "local"), .local_envir = env)
  remote
}

# An ssh stand-in that fails, an empty working directory and a temporary data
# directory for the registry of jobs.
local_no_ssh <- function(env = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = env)
  bin <- file.path(root, "bin")
  dir.create(bin); dir.create(file.path(root, "wd"))
  writeLines(c("#!/bin/sh", "exit 99"), file.path(bin, "ssh"))
  Sys.chmod(file.path(bin, "ssh"), "755")
  withr::local_envvar(PATH = paste(bin, Sys.getenv("PATH"), sep = ":"),
                      R_USER_DATA_DIR = file.path(root, "data"), .local_envir = env)
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
  expect_identical(suppressMessages(f()), list(a = 11, b = 22))
  expect_length(list.files("."), 0L)
  expect_identical(nrow(dMod2Jobs()), 0L)
  expect_length(detectFreeCores("localhost"), 1L)
})

test_that("data frame rows are tasks, a task reads its cores and seeds repeat", {
  skip_on_os("windows")
  skip_on_cran()
  local_no_ssh()
  local_mocked_bindings(detectFreeCores = function(machine) setNames(rep(8, length(machine)),
                                                                     machine))
  tasks <- data.frame(a = 1:3, b = 3:1)
  run <- function() suppressMessages(distributedComputing(
    c(task$a + task$b, cores, runif(1)), tasks = tasks, on = localNode(3), seed = 7L,
    input = character(0), jobname = "js", wait = TRUE))
  one <- run()
  expect_length(one, 3L)
  expect_identical(vapply(one, `[`, 0, 1L), rep(4, 3))
  expect_identical(vapply(one, `[`, 0, 2L), rep(3, 3))
  expect_identical(run(), one)
  expect_length(unique(vapply(one, `[`, 0, 3L)), 3L)
})

test_that("remote nodes run their tasks, and recover, the registry and purge reach them", {
  skip_on_os("windows")
  skip_on_cran()
  remote <- local_fake_hosts()
  local_mocked_bindings(detectFreeCores = function(machine) setNames(rep(8, length(machine)),
                                                                     machine))
  y <- 3
  expect_message(
    job <- distributedComputing(y * taskId, tasks = 3, on = sshNodes(c("h1", "h2")),
                                input = "y", jobname = "jr"),
    "jr: 3 tasks started on h1 \\(8 cores\\), h2 \\(8 cores\\)")
  expect_identical(attr(job, "manifest")$taskHost, c(1L, 2L, 1L))
  expect_s3_class(job, "dMod2Job")
  expect_true(wait_for(job))
  expect_identical(job$get(), list(3, 6, 9))
  jobs <- dMod2Jobs()
  expect_identical(jobs$jobname, "jr")
  expect_identical(jobs$on, "h1 (8 cores), h2 (8 cores)")
  expect_identical(jobs$wd, getwd())
  elsewhere <- withr::with_dir(tempdir(),
                               distributedComputing(NULL, jobname = "jr", recover = TRUE))
  expect_identical(elsewhere$get(), list(3, 6, 9))
  again <- distributedComputing(NULL, jobname = "jr", recover = TRUE)
  again$purge(local = TRUE)
  expect_length(list.files(remote), 0L)
  expect_length(list.files("."), 0L)
  expect_identical(nrow(dMod2Jobs()), 0L)
  expect_error(distributedComputing(NULL, jobname = "jr", recover = TRUE), "no job 'jr'")
})

test_that("a job ships the model files of its objects, else those of dMod2.outdir", {
  skip_on_os("windows")
  models <- withr::local_tempdir()
  file.create(file.path(models, c("a.cpp", "a.o", "a.so", "b.cpp")))
  info <- function(f) structure(list(), compileInfo = list(list(srcfile = file.path(models, f))))
  files <- function(objs, build) basename(dMod2:::.jobFiles(objs, build))
  expect_identical(files(list(o = info("a.cpp")), "compile"), "a.cpp")
  expect_identical(files(list(o = info("a.cpp")), "link"), "a.o")
  expect_identical(files(list(o = info("a.cpp")), "copy"), "a.so")
  expect_error(dMod2:::.jobFiles(list(o = info("b.cpp")), "link"), "b.o")
  withr::local_options(dMod2.outdir = models)
  expect_setequal(files(list(x = 1), "compile"), c("a.cpp", "b.cpp"))
  expect_identical(files(list(x = 1), "copy"), "a.so")
})

test_that("a remote build compiles the sources of the objects into one shared object", {
  skip_on_os("windows")
  skip_on_cran()
  local_fake_hosts()
  local_mocked_bindings(detectFreeCores = function(machine) setNames(rep(2, length(machine)),
                                                                     machine))
  models <- withr::local_tempdir()
  writeLines("void dmod_probe(void) {}", file.path(models, "probe.c"))
  writeLines("void dmod_other(void) {}", file.path(models, "other.c"))
  withr::local_options(dMod2.outdir = models)
  model <- structure(list(), compileInfo = list(list(
    srcfile = file.path(models, "probe.c"), compileArgs = "", linkArgs = "")))
  out <- suppressMessages(distributedComputing(
    c(is.loaded("dmod_probe"), file.exists("other.c")), on = sshNodes("h1"),
    input = "model", jobname = "jb", wait = TRUE))
  expect_identical(out, list(c(TRUE, FALSE)))
  expect_length(list.files("."), 0L)

  writeLines("not C", file.path(models, "broken.c"))
  attr(model, "compileInfo")[[1]]$srcfile <- file.path(models, "broken.c")
  expect_error(suppressMessages(distributedComputing(
    1, on = sshNodes("h1"), input = "model", jobname = "jx")),
    "the start on h1 failed:[\\s\\S]*error", perl = TRUE)
})

test_that("build = 'copy' ships the shared objects of dMod2.outdir", {
  skip_on_os("windows")
  skip_on_cran()
  local_fake_hosts()
  local_mocked_bindings(detectFreeCores = function(machine) setNames(rep(2, length(machine)),
                                                                     machine))
  models <- withr::local_tempdir()
  writeLines("void dmod_probe(void) {}", file.path(models, "probe.c"))
  withr::with_dir(models, system2(file.path(R.home("bin"), "R"), c("CMD", "SHLIB", "probe.c"),
                                  stdout = FALSE, stderr = FALSE))
  withr::local_options(dMod2.outdir = models)
  expect_identical(suppressMessages(distributedComputing(
    file.exists("probe.so"), on = sshNodes("h1"), build = "copy", input = character(0),
    jobname = "jo", wait = TRUE)), list(TRUE))
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
  expect_message(job <- distributedComputing(
    c(task$p, Sys.getenv("FROM_TEMPLATE"), cores), tasks = data.frame(p = c("a", "b")),
    on = slurmCluster("h1", cores = 2, walltime = "00:10:00"), jobname = "jsl",
    input = character(0)),
    "jsl: 2 tasks submitted to h1 as SLURM job 4242, 2 array tasks of 2 cores")
  expect_identical(attr(job, "manifest")$jobid, "4242")
  expect_true(wait_for(job))
  expect_identical(job$get(), list(c("a", "1", "2"), c("b", "1", "2")))
  script <- readLines(file.path(remote, "jsl_folder", "1", "job.sh"))
  expect_true(all(c("#SBATCH --array=1-2", "#SBATCH --time=00:10:00",
                    "#SBATCH --cpus-per-task=2") %in% script))
  job$purge(local = TRUE)
  expect_length(list.files(remote), 0L)
})

test_that("an array smaller than the tasks runs them in turn, with the options of sbatch", {
  skip_on_os("windows")
  skip_on_cran()
  remote <- local_fake_hosts()
  tpl <- withr::local_tempfile(lines = c("#!/bin/bash", "#SBATCH --array=1-{{ ntasks }}",
                                         "{{ run }}"))
  on <- slurmCluster("h1", cores = 1, array = 2, partition = "p", maxRunning = 1,
                     sbatch = "--constraint=x", jobscript = tpl)
  expect_output(print(on), "SLURM cluster h1: 2 array tasks of 1 cores")
  job <- distributedComputing(c(taskId, as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))),
                              tasks = 5, on = on, jobname = "ja", input = character(0))
  expect_true(wait_for(job))
  expect_identical(job$get(), lapply(1:5, function(i) c(i, (i - 1L) %% 2L + 1L)))
  expect_identical(readLines(file.path(remote, "ja_folder", "1", "sbatch.args")),
                   "--partition=p --array=1-2%1 --constraint=x job.sh")
  job$purge(local = TRUE)
})

test_that("the targets check their nodes and cores", {
  local_mocked_bindings(detectFreeCores = function(machine) setNames(rep(2, length(machine)),
                                                                     machine))
  expect_identical(sshNodes(cores = c(b = 3, a = 2))$cores, c(b = 3L, a = 2L))
  expect_identical(sshNodes(c("a", "b"))$cores, c(a = NA_integer_, b = NA_integer_))
  expect_identical(sshNodes(c("a", "b"), cores = 4)$cores, c(a = 4L, b = 4L))
  expect_output(print(localNode(2)), "localhost \\(2 cores\\)")
  expect_error(sshNodes("a", cores = 1.5), "positive whole number")
  expect_error(sshNodes(c("a", "a")), "every node once")
  expect_error(sshNodes(c("a", "b"), cores = c(a = 1)), "names of 'cores'")
  expect_error(sshNodes(c("a", "b"), cores = c(1, 2)), "one number or a vector named")
  expect_error(slurmCluster("h"), "cores of each array task")
  expect_error(slurmCluster("h", 2, walltime = "1h"), "HH:MM:SS")
  expect_error(slurmCluster("h", 2, array = 0), "'array' must be")
  expect_error(distributedComputing(1, on = c(a = 2)), "'on' must be")
  expect_error(distributedComputing(1, tasks = 0), "positive whole number")
})

test_that("a node gets at most its free cores", {
  skip_on_os("windows")
  skip_on_cran()
  local_no_ssh()
  local_mocked_bindings(detectFreeCores = function(machine) setNames(rep(2, length(machine)),
                                                                     machine))
  expect_warning(out <- suppressMessages(distributedComputing(
    cores, on = localNode(5), input = character(0), jobname = "jc", wait = TRUE)),
    "2 free cores; its tasks get 2 instead of 5")
  expect_identical(out, list(2L))
  expect_identical(suppressMessages(distributedComputing(cores, input = character(0),
                                                         jobname = "jf", wait = TRUE)),
                   list(2L))
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
  expect_error(distributedComputing(1, on = slurmCluster("h", cores = 2)),
               "no job script for 'h'; set one up once with slurmTemplate\\(\"h\"")
  expect_error(distributedComputing(1, on = slurmCluster("h", cores = 2, jobscript = "nope.sh")),
               "nope.sh does not exist")
})

test_that("a slot left as a zombie counts as finished", {
  skip_on_os(c("windows", "mac"))
  skip_if_not(file.exists("/proc/self/status"))
  local_no_ssh()
  dir.create(file.path("jz_folder", "1"), recursive = TRUE)
  # the shell forks a child and becomes `sleep 5`, which never reaps it
  system(paste("cd jz_folder/1 && sh -c 'sleep 0 & echo $! > slot_1.pid; exec sleep 5'"),
         wait = FALSE)
  state <- function(pid) {
    f <- file.path("/proc", pid, "status")
    if (!file.exists(f)) return("")
    sub("^State:\\s*(\\S).*", "\\1", grep("^State:", readLines(f), value = TRUE))
  }
  pidFile <- file.path("jz_folder", "1", "slot_1.pid")
  for (i in 1:50) if (file.exists(pidFile) && identical(state(readLines(pidFile)), "Z")) break else
    Sys.sleep(0.1)
  skip_if_not(identical(state(readLines(pidFile)), "Z"))
  m <- list(folder = "jz_folder", hosts = "localhost", slurm = FALSE, ntasks = 1L, taskHost = 1L,
            wd = getwd())
  expect_identical(dMod2:::.taskStatus(m), "failed")
})

test_that("get() joins the results and leaves out tasks without one", {
  fit <- function(v) parlist(list(value = v, argument = c(a = v), converged = TRUE,
                                  iterations = 1L))
  failed <- structure("Error : boom\n", class = "try-error", condition = simpleError("boom"))
  expect_warning(joined <- dMod2:::.combineResults(list(fit(1), NULL, failed, fit(2))),
                 "2 of 4 tasks have no result; not done: 2; failed: 3, the first with: boom")
  expect_s3_class(joined, "parlist")
  expect_length(joined, 2L)
  expect_identical(dMod2:::.combineResults(list(data.frame(a = 1), data.frame(a = 2)))$a, c(1, 2))
  expect_identical(dMod2:::.combineResults(list(a = 1, b = "x")), list(a = 1, b = "x"))
})

test_that("a failing task returns its error, a running job keeps its name, and terminate() stops it", {
  skip_on_os("windows")
  skip_on_cran()
  local_no_ssh()
  expect_warning(out <- suppressMessages(distributedComputing(
    stop("no"), input = character(0), jobname = "je", wait = TRUE)),
    "1 of 1 tasks failed; the job is kept")
  expect_length(out, 0L)
  failed <- distributedComputing(NULL, jobname = "je", recover = TRUE)
  expect_output(failed$log(1), "no")
  expect_s3_class(failed$get(combine = FALSE)[[1]], "try-error")
  failed$purge(local = TRUE)

  job <- distributedComputing(Sys.sleep(60), input = character(0), jobname = "jt")
  Sys.sleep(1)
  utils::capture.output(running <- dMod2:::.taskStatus(attr(job, "manifest")))
  expect_identical(running, "running")
  expect_error(distributedComputing(1, input = character(0), jobname = "jt"),
               "the job 'jt' submitted from .* still runs 1 tasks")
  job$terminate()
  for (i in 1:50) if (dMod2:::.taskStatus(attr(job, "manifest")) != "running") break else
    Sys.sleep(0.1)
  expect_identical(dMod2:::.taskStatus(attr(job, "manifest")), "failed")
  expect_warning(utils::capture.output(finished <- job$check()), "1 of 1 tasks failed")
  expect_true(finished)
  job$purge(local = TRUE)
  expect_length(list.files("."), 0L)
})
