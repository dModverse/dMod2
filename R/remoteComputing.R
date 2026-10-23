#' Detect Number of Free Cores
#'
#' Estimates the free cores as the number of cores minus the 1-minute load
#' average, at least 1. Works on Linux and macOS, locally or on remote
#' machines via `ssh`. On Windows it returns 1 with a warning. Under
#' `R CMD check` with `_R_CHECK_LIMIT_CORES_` set, the result is at most 2.
#' @param machine Character vector of machines: `"localhost"` is the local machine,
#'   any other entry an `ssh` destination such as `"host"` or `"user@@host"`.
#'   `NULL` (default) for the local machine.
#' @return Numeric vector of free cores, one per machine and named by it, with
#'   attributes
#'   `"ncores"` (number of cores) and `"used"` (load average).
#' @export
#' @examples
#' detectFreeCores()
detectFreeCores <- function(machine = NULL) {
  
  .getLoadAndCores <- function(prefix = NULL) {
    cmd <- function(x) {
      if (!is.null(prefix)) x <- paste("ssh", prefix, x)
      system(x, intern = TRUE)
    }
    
    # Detect OS: use uname for remote, Sys.info() for local
    os <- if (!is.null(prefix)) cmd("uname") else Sys.info()[["sysname"]]
    
    if (grepl("Windows", os, ignore.case = TRUE)) {
      # No load average on Windows, return 1 free core as safe default
      warning("detectFreeCores: load average not available on Windows, returning 1")
      nCores <- parallel::detectCores()
      return(list(free = 1, nCores = nCores, occupied = NA_real_))
    }
    
    if (grepl("Darwin", os, ignore.case = TRUE)) {
      # macOS: sysctl provides load avg as "{ x.xx x.xx x.xx }"
      occupied <- as.numeric(strsplit(cmd("sysctl -n vm.loadavg"), " ")[[1]][2])
      nCores <- as.numeric(cmd("sysctl -n hw.ncpu"))
    } else {
      # Linux: read 1-min load average from /proc/loadavg
      occupied <- as.numeric(strsplit(cmd("cat /proc/loadavg"), " ", fixed = TRUE)[[1]][1])
      nCores <- as.numeric(cmd("nproc --all"))
    }
    
    # Floor at 1: mclapply(mc.cores = 0) fails, and under heavy load the
    # 1-min average can exceed nCores.
    list(free = max(1L, round(nCores - occupied)), nCores = nCores, occupied = occupied)
  }
  
  if (!is.null(machine)) {
    output <- lapply(machine, function(m) .getLoadAndCores(if (!.isLocalhost(m)) m))
    freeCores <- setNames(vapply(output, `[[`, numeric(1), "free"), machine)
    attr(freeCores, "ncores") <- vapply(output, `[[`, numeric(1), "nCores")
    attr(freeCores, "used") <- vapply(output, `[[`, numeric(1), "occupied")
  } else {
    res <- .getLoadAndCores()
    freeCores <- res$free
    attr(freeCores, "ncores") <- res$nCores
    attr(freeCores, "used") <- res$occupied
  }

  # Under R CMD check, _R_CHECK_LIMIT_CORES_ forbids mc.cores > 2; capping here
  # keeps every mclapply caller legal.
  chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  if (nzchar(chk) && chk != "false") {
    if (length(freeCores) > 0L) freeCores[] <- pmin(freeCores, 2L)
  }

  freeCores
}


## Remote build helpers --------------------------------------------------------
## The remote build keeps the -D macros of "compileInfo" and resolves every
## path-valued flag through the R on the remote PATH.

## Shell command that puts `libs` in front of the remote R's library path. The
## paths are expanded by the remote shell, so `~` is the remote home.
.remoteLibs <- function(libs) {
  if (!length(libs)) return("")
  paths <- sub("^~", "$HOME", libs)
  paste0('export R_LIBS="', paste(paths, collapse = ":"), '${R_LIBS:+:$R_LIBS}"; ')
}

## Scan a workspace for model objects and return the portable build settings
## needed to rebuild their C/C++ sources elsewhere.
.remoteBuildInfo <- function(envir = .GlobalEnv) {

  compileArgs <- character(0)
  needsCVODE <- FALSE
  needsKLU <- FALSE
  needsLapack <- FALSE

  for (nm in ls(envir)) {
    o <- try(get(nm, envir = envir), silent = TRUE)
    if (inherits(o, "try-error")) next
    if (isTRUE(attr(o, "sparse"))) needsKLU <- TRUE
    info <- attr(o, "compileInfo")
    if (is.null(info)) next
    for (e in info) {
      ca <- trimws(e$compileArgs %||% "")
      ## Non-empty linkArgs mean external libraries (Sundials, KLU), resolved
      ## from the remote cppDE installation.
      la <- e$linkArgs %||% ""
      if (nzchar(trimws(la))) needsCVODE <- TRUE
      if (grepl("sunlinsollapackdense", la, fixed = TRUE)) needsLapack <- TRUE
      if (isTRUE(e$sparse)) needsKLU <- TRUE
      if (nzchar(ca)) compileArgs <- c(compileArgs, strsplit(ca, "\\s+")[[1]])
    }
  }

  ## Path-valued flags point into the local library tree; the remote build
  ## resolves its own from the cppDE installed there.
  compileArgs <- unique(compileArgs[nzchar(compileArgs) &
                                      !grepl("^-(I|L|l|Wl,|isystem)", compileArgs)])
  ## cppDE tags sparse models with -DKLUBTF/-DKLUAMD; compile() adds the
  ## -DKLU switch itself, so mirror that here.
  if (any(grepl("^-DKLU", compileArgs))) needsKLU <- TRUE

  list(compileArgs = paste(compileArgs, collapse = " "),
       needsCVODE  = needsCVODE,
       needsKLU    = needsKLU,
       needsLapack = needsLapack)
}

## TRUE when naming every file on one command line would overrun the shell's
## argument limit; same test as compile() applies locally.
.remoteNeedsChunking <- function(files)
  sum(nchar(files) + 1L) > .compileCmdLimit()

## Write the job workspace through zstd, as an uncompressed one can be several
## GB; without the binary it is saved uncompressed. Returns the file written.
.saveWorkspace <- function(input, file, envir = .GlobalEnv, level = 3L) {

  zstd <- Sys.which("zstd")
  if (!nzchar(zstd)) {
    save(list = input, file = file, envir = envir, compress = FALSE)
    return(file)
  }

  out <- paste0(file, ".zst")
  con <- pipe(paste0(shQuote(zstd), " -T0 -", as.integer(level),
                     " -q -f -o ", shQuote(out)), "wb")
  ## close() pcloses, so zstd has finished before tar sees the file.
  ok <- try(save(list = input, file = con, envir = envir), silent = TRUE)
  close(con)
  if (inherits(ok, "try-error")) {
    unlink(out)
    stop("Could not write the job workspace: ", attr(ok, "condition")$message,
         call. = FALSE)
  }
  out
}

## Bash script that builds the shared object on the remote machine from the
## sources, or the object files when `link = TRUE`. Beyond the argument limit it
## uses compile()'s chunked archive build, reading its inputs from `filelist`.
.remoteBuildScript <- function(files, output, compileArgs = "",
                               needsCVODE = FALSE, needsKLU = FALSE,
                               needsLapack = FALSE,
                               link = FALSE, cxx = FALSE, cores = 1,
                               workdir = NULL, filelist = NULL, chunkSize = 100,
                               bundle = 50) {

  files <- unlist(strsplit(trimws(files), "\\s+"))
  files <- files[nzchar(files)]

  cflags <- trimws(paste("-O2 -DNDEBUG -w -fPIC", compileArgs,
                         if (needsKLU) "-DKLU"))

  ldflags <- character(0)
  if (link) {
    ## A remote compiler of another GCC generation refuses the local LTO
    ## bytecode; with LTO off the link uses the fat objects' machine code.
    ldflags <- c(ldflags, "-fno-lto", "-fno-use-linker-plugin")
    ## With only .o inputs R CMD SHLIB links via the C driver, which would not
    ## pull in the C++ runtime the cppDE-generated code needs.
    if (cxx) ldflags <- c(ldflags, "-lstdc++")
  }

  ## SUNDIALS/KLU flags come from the remote cppDE install, read by Rscript in
  ## the script. Embedded in a single-quoted shell string: no single quotes.
  cfgExpr <- function(field) paste0(
    "cfg <- get0(\"cvodeConfig\", envir = asNamespace(\"cppDE\"), inherits = FALSE); ",
    "cat(if (is.environment(cfg)) paste(unlist(mget(c(", field, "), envir = cfg, ",
    "ifnotfound = \"\")), collapse = \" \") else \"\")")

  cfgFields <- function(...) paste0("\"", c(...), "\"", collapse = ", ")

  extra_libs <- NULL
  if (needsKLU || needsCVODE) {
    fields <- cfgFields(if (needsKLU) "klu_libs",
                        if (needsLapack) "cvode_lapack_libs",
                        if (needsCVODE) "libs")
    extra_libs <- paste0("PKG_LIBS=\"$PKG_LIBS $(Rscript -e '", cfgExpr(fields), "')\"")
  }

  ## Include paths, MPI's included when the remote SUNDIALS was built with it.
  ext_cppflags <- if (needsKLU || needsCVODE)
    paste0("PKG_CPPFLAGS=\"$PKG_CPPFLAGS $(Rscript -e '",
           cfgExpr(cfgFields(if (needsCVODE) "cflags", if (needsKLU) "klu_cflags")),
           "')\"")

  ## The generated source calls SUNLinSol_LapackDense when the submitting
  ## machine had it, so the remote SUNDIALS must provide it too.
  lapack_check <- if (needsLapack) c(
    paste0("if [ -z \"$(Rscript -e '", cfgExpr(cfgFields("cvode_lapack_libs")), "')\" ]; then"),
    "  echo \"dMod: the model uses SUNDIALS' LAPACK dense solver, which the\" >&2",
    "  echo \"      SUNDIALS of this machine lacks. Build one with\" >&2",
    "  echo \"      cppDE::install_libs() and reinstall cppDE here.\" >&2",
    "  exit 1",
    "fi",
    "")

  ## Job count: a fixed `cores`, or whatever the remote machine reports.
  nproc <- if (is.null(cores)) "$NPROC" else as.character(max(1L, as.integer(cores)))
  ## nproc honours OMP_NUM_THREADS, which the job sets to 1 for itself
  detect <- if (is.null(cores))
    c("NPROC=$(env -u OMP_NUM_THREADS -u OMP_THREAD_LIMIT nproc 2>/dev/null || echo 1)",
      "if [ \"$NPROC\" -gt 16 ]; then NPROC=16; fi", "")

  ## Precompiled header; a missing .gch is harmless, the sources then include
  ## the headers themselves.
  cxxSrc <- files[grepl("\\.cpp$", files, ignore.case = TRUE)]
  pchInc <- if (!link && length(cxxSrc) >= 8L) .compilePCHIncludes(cxxSrc)
  toolchain <- if (!is.null(filelist) || !is.null(pchInc)) c(
    "CC=$(R CMD config CC);   CFLAGS=$(R CMD config CFLAGS)",
    "CXX=$(R CMD config CXX); CXXFLAGS=$(R CMD config CXXFLAGS)",
    "CPICFLAGS=$(R CMD config CPICFLAGS); CXXPICFLAGS=$(R CMD config CXXPICFLAGS)",
    ## R.home("include"), not R_HOME/include: Debian keeps the headers outside R_HOME
    "RINC=\"-I$(Rscript -e 'cat(R.home(\"include\"))')\"",
    "AR=$(R CMD config AR); RANLIB=$(R CMD config RANLIB)",
    "export CC CXX CFLAGS CXXFLAGS CPICFLAGS CXXPICFLAGS RINC", "")
  pchBlock <- if (!is.null(pchInc)) c(
    "cat > dMod_pch.hpp <<'DMOD_PCH_EOF'", pchInc, "DMOD_PCH_EOF",
    paste("$CXX $RINC $PKG_CPPFLAGS $PKG_CXXFLAGS $CXXPICFLAGS $CXXFLAGS",
          "-x c++-header dMod_pch.hpp -o dMod_pch.hpp.gch || true"),
    "PKG_CXXFLAGS=\"$PKG_CXXFLAGS -include dMod_pch.hpp\"", "")

  ## Compile off the file list in parallel, then link: R CMD SHLIB does not
  ## parallelise via MAKEFLAGS.
  chunked <- !is.null(filelist) && .remoteNeedsChunking(files)
  isCxx   <- grepl("\\.cpp$", files, ignore.case = TRUE)

  ## Only the chunked path bundles: below the argument limit SHLIB gets the
  ## sources themselves, so bundles would just be compiled twice.
  doBundle <- chunked && !link && as.integer(bundle) > 1L && sum(isCxx) > 1L

  anchorLine <- if (chunked)
    paste0("ANCHOR=", shQuote(files[if (any(isCxx)) which(isCxx)[1] else 1L]))

  ## Bundling the small generated sources into few translation units parses the
  ## shared headers once per bundle.
  bundleBlock <- if (doBundle) c(
    paste0("BUNDLE=", max(2L, as.integer(bundle))),
    "grep -i -e '\\.cpp$' \"$WORK\" > dmod_cxx.lst || true",
    "grep -v -i -e '\\.cpp$' \"$WORK\" > dmod_rest.lst || true",
    "if [ -s dmod_cxx.lst ]; then",
    "  rm -f dmod_bundle_*",
    "  split -l \"$BUNDLE\" -d -a 4 dmod_cxx.lst dmod_bundle_",
    "  for b in dmod_bundle_[0-9]*; do",
    "    sed 's|.*|#include \"&\"|' \"$b\" > \"$b.cpp\"",
    "    echo \"$b.cpp\"",
    "  done > dmod_bundled.lst",
    "  cat dmod_rest.lst dmod_bundled.lst > \"$WORK\"",
    "fi", "")

  compileBlock <- if (!is.null(filelist)) c(
    paste0("FILELIST=", shQuote(filelist)),
    anchorLine, "",
    "dmod_compile_one() {",
    "  src=\"$1\"",
    "  case \"$src\" in",
    "    *.o)   : ;;                                   # link-only: already built",
    "    *.cpp) $CXX $RINC $PKG_CPPFLAGS $PKG_CXXFLAGS $CXXPICFLAGS $CXXFLAGS \\",
    "             -c \"$src\" -o \"${src%.*}.o\" ;;",
    "    *)     $CC  $RINC $PKG_CPPFLAGS $PKG_CFLAGS  $CPICFLAGS   $CFLAGS \\",
    "             -c \"$src\" -o \"${src%.*}.o\" ;;",
    "  esac",
    "}",
    "export -f dmod_compile_one", "",
    "WORK=dmod_work.lst",
    if (chunked) "grep -v -x -F \"$ANCHOR\" \"$FILELIST\" > \"$WORK\" || true"
    else "cp \"$FILELIST\" \"$WORK\"",
    "",
    bundleBlock,
    paste0("xargs -r -P ", nproc,
           " -I '{}' bash -c 'dmod_compile_one \"$@\"' _ '{}' < \"$WORK\""), "")

  ## Past the argument limit the objects go into a static archive and SHLIB gets
  ## one anchor plus the archive. The anchor is a C++ source when there is one,
  ## so SHLIB still selects the C++ linker and runtime.
  linkBlock <- if (!chunked) {
    paste0("R CMD SHLIB ", paste(files, collapse = " "), " -o ", output)
  } else {
    c(paste0("LIB=\"$PWD/", sub("\\.so$", "", output), "_objects.a\""), "",
      "# The anchor is linked directly; in the archive too it would define its",
      "# symbols twice.",
      "rm -f \"$LIB\"",
      paste0("sed 's/\\.[^.]*$/.o/' \"$WORK\" |",
             " xargs -r -n ", max(1L, as.integer(chunkSize)), " \"$AR\" qc \"$LIB\""),
      "\"$RANLIB\" \"$LIB\"", "",
      "# R resolves the entry points by name at run time, so unreferenced archive",
      "# members have to be kept.",
      "case \"$(uname -s)\" in",
      "  Darwin) PKG_LIBS=\"-Wl,-force_load,$LIB $PKG_LIBS\" ;;",
      "  *)      PKG_LIBS=\"-Wl,--whole-archive $LIB -Wl,--no-whole-archive $PKG_LIBS\" ;;",
      "esac",
      "export PKG_LIBS", "",
      paste0("R CMD SHLIB \"$ANCHOR\" -o ", output),
      "rm -f \"$LIB\"")
  }
  buildlines <- c(compileBlock, linkBlock)

  paste(c(
    "#!/bin/bash",
    "# Generated by dMod: rebuilds the model shared object on this machine.",
    "# All path-valued compiler flags are resolved here, not on the submitting",
    "# machine, because the library trees differ.",
    "set -e",
    "# R CMD check exports R_TESTS to its test subprocesses; an Rscript started",
    "# below would source that nonexistent startup file and die with exit 1.",
    "unset R_TESTS",
    "",
    if (!is.null(workdir)) c(paste0("cd ", workdir), ""),
    "CPPDE_INC=$(Rscript -e 'cat(system.file(\"include\", package = \"cppDE\"))')",
    "if [ -z \"$CPPDE_INC\" ]; then",
    "  echo \"dMod: package 'cppDE' is not installed for the R on PATH here,\" >&2",
    "  echo \"      so the model sources cannot be rebuilt. Install it remotely\" >&2",
    "  echo \"      (remotes::install_github('dModverse/cppDE')) or submit with\" >&2",
    "  echo \"      compile = FALSE and link = FALSE.\" >&2",
    "  exit 1",
    "fi",
    "",
    lapack_check,
    "PKG_CPPFLAGS=\"-I$CPPDE_INC\"",
    ext_cppflags,
    paste0("PKG_CFLAGS=\"", cflags, "\""),
    paste0("PKG_CXXFLAGS=\"", cflags, "\""),
    "PKG_LIBS=\"$(R CMD config LAPACK_LIBS) $(R CMD config BLAS_LIBS)\"",
    extra_libs,
    if (length(ldflags))
      paste0("PKG_LIBS=\"", paste(ldflags, collapse = " "), " $PKG_LIBS\""),
    detect,
    toolchain,
    pchBlock,
    "export PKG_CPPFLAGS PKG_CFLAGS PKG_CXXFLAGS PKG_LIBS",
    "",
    buildlines
  ), collapse = "\n")
}


## ---- distributedComputing() ---------------------------------------------------
## A job lives in <jobname>_folder/<k>, one subfolder per node k: in the working
## directory for "localhost", in the home directory elsewhere. The local
## <jobname>_folder holds the manifest job.rds and the files get() fetches; a copy of
## the manifest is registered in the data directory of dMod2.

## "localhost" is the local machine, reached without ssh; every other node by ssh.
.isLocalhost <- function(machine) identical(machine, "localhost")

## `cmd` run by the POSIX shell of `machine` in its job folder `dir`. On
## "localhost" `dir` lies in `wd`, and the R of this session with its library
## paths comes first.
.onMachine <- function(machine, cmd, dir, wd, ...) {
  if (.isLocalhost(machine)) {
    env <- paste0("export PATH=", shQuote(R.home("bin")), ":\"$PATH\" R_LIBS=",
                  shQuote(paste(.libPaths(), collapse = ":")), "; ")
    return(system(paste0(env, "cd ", shQuote(file.path(wd, dir)), " && ", cmd), ...))
  }
  system(paste("ssh", machine, shQuote(paste0("cd ", dir, " && ", cmd))), ...)
}

## The files of `stage` into an emptied job folder `dir` of `machine`. Symbolic
## links in the stage are followed.
.pushStage <- function(machine, stage, dir, wd) {
  unpack <- paste0("rm -rf ", dir, " && mkdir -p ", dir, " && tar -C ", dir, " -xf -")
  if (.isLocalhost(machine))
    unpack <- paste0("(cd ", shQuote(wd), " && ", unpack, ")")
  else
    unpack <- paste("ssh", machine, shQuote(unpack))
  system(paste0("tar -C ", shQuote(stage), " -chf - . | ", unpack))
}

## The files named by the shell pattern `what` of the job folder `dir` of a
## remote `machine`, copied into the same folder below `wd`.
.pullFiles <- function(machine, dir, what, wd) {
  local <- file.path(wd, dir)
  dir.create(local, recursive = TRUE, showWarnings = FALSE)
  system(paste0("ssh ", machine, " ", shQuote(paste0("cd ", dir, " && tar -cf - ", what,
                                                     " 2>/dev/null")),
                " | tar -C ", shQuote(local), " -xf -"),
         ignore.stderr = TRUE)
}

## Tasks as a list, one element per task: a count gives NULL elements, a data
## frame its rows as lists.
.asTasks <- function(tasks) {
  if (is.numeric(tasks) && length(tasks) == 1L && !is.list(tasks)) {
    n <- as.integer(tasks)
    if (is.na(n) || n < 1L || n != tasks)
      stop("'tasks' must be a positive whole number, a list or a data frame", call. = FALSE)
    return(vector("list", n))
  }
  if (is.data.frame(tasks)) {
    if (!nrow(tasks)) stop("'tasks' has no rows", call. = FALSE)
    rows <- lapply(seq_len(nrow(tasks)), function(i) as.list(tasks[i, , drop = FALSE]))
    if (!identical(rownames(tasks), as.character(seq_len(nrow(tasks)))))
      names(rows) <- rownames(tasks)
    return(rows)
  }
  if (!length(tasks)) stop("'tasks' is empty", call. = FALSE)
  as.list(tasks)
}

## Seconds of a walltime "HH:MM:SS".
.walltimeSeconds <- function(walltime) {
  p <- suppressWarnings(as.numeric(strsplit(walltime, ":", fixed = TRUE)[[1L]]))
  if (length(p) != 3L || anyNA(p))
    stop("'walltime' must be \"HH:MM:SS\"", call. = FALSE)
  sum(p * c(3600, 60, 1))
}

## The job script of a SLURM template: every {{ name }} of `values` filled in.
## The lines before {{ run }} other than #SBATCH and comments are its
## environment, which the build on the login node uses too.
.fillJobscript <- function(template, values) {
  lines <- readLines(template, warn = FALSE)
  at <- grep("^\\s*\\{\\{\\s*run\\s*\\}\\}\\s*$", lines)
  if (length(at) != 1L)
    stop("the job script ", template, " needs one line {{ run }}", call. = FALSE)
  env <- lines[seq_len(at - 1L)]
  env <- env[!grepl("^\\s*(#|$)", env)]
  values <- Filter(Negate(is.null), values)
  fill <- function(x) {
    for (k in names(values))
      x <- gsub(paste0("\\{\\{\\s*", k, "\\s*\\}\\}"), values[[k]], x)
    x
  }
  script <- fill(lines)
  left <- unique(unlist(regmatches(script, gregexpr("\\{\\{\\s*[A-Za-z]+\\s*\\}\\}", script))))
  left <- left[!grepl("^\\{\\{\\s*run\\s*\\}\\}$", left)]
  if (length(left))
    stop("the job script ", template, " has placeholders without a value: ",
         paste(left, collapse = ", "), call. = FALSE)
  list(script = script, env = fill(env))
}

## The job script templates of slurmTemplate(), one per login node.
.slurmTemplateDir <- function() file.path(tools::R_user_dir("dMod2", "config"), "slurm")

## The templates dMod2 ships, named by cluster.
.slurmShipped <- function() {
  files <- list.files(system.file("slurm", package = "dMod2"), pattern = "\\.sh$",
                      full.names = TRUE)
  stats::setNames(files, sub("\\.sh$", "", basename(files)))
}

## The template set up for a login node.
.slurmJobscript <- function(host) {
  f <- file.path(.slurmTemplateDir(), paste0(host, ".sh"))
  if (file.exists(f)) return(f)
  stop("no job script for '", host, "'; set one up once with slurmTemplate(\"", host,
       "\", from = <job script>) or pass slurmCluster(jobscript = ), see ?slurmTemplate",
       call. = FALSE)
}

## "1 task", "2 tasks".
.count <- function(n, what) paste(n, if (n == 1L) what else paste0(what, "s"))

## A positive whole number as integer; NULL stays NULL where `null` allows it.
.wholeNumber <- function(x, name, null = TRUE) {
  if (is.null(x) && null) return(NULL)
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 1 || x != round(x))
    stop("'", name, "' must be a positive whole number", call. = FALSE)
  as.integer(x)
}


#' Where distributedComputing() Runs its Tasks
#'
#' `localNode()` is the local machine, `sshNodes()` are machines reached by
#' `ssh` and `slurmCluster()` is a SLURM cluster reached through its login
#' node. [distributedComputing()] takes one of them as `on`.
#'
#' The tasks of a job run on workers: a node is one worker, and a SLURM job has
#' `array` workers, its array tasks. The tasks are dealt to the workers in turn,
#' and a worker runs its tasks one after the other, each with `cores` cores.
#' Without SLURM the processes run detached and survive the R session. With
#' SLURM the array job is submitted on the login node through the job script
#' template [slurmTemplate()] set up for it; `partition`, `account`,
#' `maxRunning` and `sbatch` are options of `sbatch`, which take precedence over
#' the `#SBATCH` lines of the template.
#'
#' @param cores The cores of every task. For `sshNodes()` one number for all
#'   nodes or a vector named by node. `NULL` (default) uses the free cores of
#'   each node, see [detectFreeCores()]; a node gets at most its free cores. For
#'   `slurmCluster()` the cores of each array task.
#' @param nodes Character, the nodes, each once: `ssh` destinations such as
#'   `"host"` or `"user@@host"` with key authentication, or `"localhost"` for the
#'   local machine. Defaults to the names of `cores`.
#' @param login Character, the login node, an `ssh` destination.
#' @param array The number of array tasks, at most one per task. `NULL`
#'   (default) is one per task.
#' @param walltime The maximum runtime of an array task, `"HH:MM:SS"`, filled
#'   into the job script. Default `NULL`, none.
#' @param memPerCore Memory per core in GB, filled into the job script. Default 2.
#' @param partition,account Character, the partition and the account of the
#'   job. Default `NULL`, those of the job script.
#' @param maxRunning The most array tasks running at once. Default `NULL`, no
#'   limit.
#' @param sbatch Character vector of further options of `sbatch`, such as
#'   `"--constraint=ib"`.
#' @param jobscript Path to a job script template, in place of the one
#'   [slurmTemplate()] set up for `login`. Default `NULL`.
#' @return An object of class `dMod2Nodes` or `dMod2Slurm`.
#' @seealso [distributedComputing()], [slurmTemplate()]
#' @examples
#' localNode(cores = 4)
#' sshNodes(c("server1", "server2"), cores = 10)
#' sshNodes(cores = c(server1 = 24, server2 = 20))
#' slurmCluster("mycluster", cores = 8, array = 32, walltime = "12:00:00")
#' @export
localNode <- function(cores = NULL) sshNodes("localhost", cores)

#' @rdname localNode
#' @export
sshNodes <- function(nodes = names(cores), cores = NULL) {
  if (!is.character(nodes) || !length(nodes) || anyNA(nodes) || any(!nzchar(nodes)) ||
      anyDuplicated(nodes))
    stop("'nodes' must name every node once, such as c(\"server1\", \"server2\")", call. = FALSE)
  if (is.null(cores)) {
    cores <- rep(NA_integer_, length(nodes))
  } else if (is.null(names(cores))) {
    if (length(cores) != 1L)
      stop("'cores' must be one number or a vector named by node", call. = FALSE)
    cores <- rep(.wholeNumber(cores, "cores"), length(nodes))
  } else {
    if (anyDuplicated(names(cores)) || !setequal(names(cores), nodes))
      stop("the names of 'cores' must be the nodes, each once", call. = FALSE)
    cores <- vapply(cores[nodes], .wholeNumber, 0L, name = "cores", null = FALSE)
  }
  structure(list(hosts = nodes, cores = stats::setNames(as.integer(cores), nodes)),
            class = c("dMod2Nodes", "dMod2Target"))
}

#' @rdname localNode
#' @export
slurmCluster <- function(login, cores, array = NULL, walltime = NULL, memPerCore = 2,
                         partition = NULL, account = NULL, maxRunning = NULL,
                         sbatch = NULL, jobscript = NULL) {
  if (!is.character(login) || length(login) != 1L || !nzchar(login))
    stop("'login' must name the login node", call. = FALSE)
  if (missing(cores))
    stop("give the cores of each array task, such as slurmCluster(\"", login,
         "\", cores = 8)", call. = FALSE)
  if (!is.null(walltime)) .walltimeSeconds(walltime)
  structure(list(login = login, cores = .wholeNumber(cores, "cores", null = FALSE),
                 array = .wholeNumber(array, "array"), walltime = walltime,
                 memPerCore = memPerCore, partition = partition, account = account,
                 maxRunning = .wholeNumber(maxRunning, "maxRunning"),
                 sbatch = as.character(sbatch), jobscript = jobscript),
            class = c("dMod2Slurm", "dMod2Target"))
}

## One line naming the nodes of a target, or the SLURM cluster and its array tasks.
.targetLabel <- function(on) {
  if (inherits(on, "dMod2Slurm"))
    return(paste0("SLURM cluster ", on$login, ": ",
                  if (is.null(on$array)) "one array task per task" else
                    paste(on$array, "array tasks"),
                  " of ", on$cores, " cores", if (!is.null(on$walltime))
                    paste0(" for ", on$walltime)))
  paste0(on$hosts, " (", ifelse(is.na(on$cores), "free", on$cores), " cores)",
         collapse = ", ")
}

#' @export
print.dMod2Target <- function(x, ...) {
  cat(.targetLabel(x), "\n", sep = "")
  invisible(x)
}

## The R script of one task, run as `Rscript job.R <taskId> <cores>` in the job folder.
.taskScript <- function(expr, jobname, seed, shared) {
  c("library(dMod2)",
    "load('workspace.RData')",
    if (is.null(shared)) "for (so in list.files(pattern = '[.]so$')) dyn.load(so)",
    if (!is.null(shared)) c(
      paste0("dyn.load(", encodeString(paste0(shared, ".so"), quote = "'"), ")"),
      "for (nm in ls()) {",
      "  o <- get(nm)",
      "  if (inherits(o, c('obsfn', 'parfn', 'prdfn', 'objfn')))",
      paste0("    { modelname(o) <- ", encodeString(shared, quote = "'"), "; assign(nm, o) }"),
      "}"),
    "taskId <- as.integer(commandArgs(TRUE)[1L])",
    "cores <- as.integer(commandArgs(TRUE)[2L])",
    "task <- readRDS('tasks.rds')[[taskId]]",
    paste0("jobname <- ", encodeString(jobname, quote = "'")),
    paste0("set.seed(", seed, " + taskId)"),
    paste0("result <- try({\n", paste(deparse(expr, width.cutoff = 500L), collapse = "\n"),
           "\n})"),
    "saveRDS(result, paste0('task_', taskId, '.rds'))")
}

## The model files a job ships: the sources ("compile"), object files ("link") or
## shared objects ("copy") of the objects in `objs`. Objects without model files
## leave it to the files of dMod2.outdir.
.jobFiles <- function(objs, build) {
  srcs <- unique(as.character(unlist(lapply(objs, function(o)
    lapply(attr(o, "compileInfo"), `[[`, "srcfile")))))
  srcs <- srcs[file.exists(srcs)]
  files <- switch(build,
    compile = srcs,
    link    = sub("\\.(c|cpp)$", ".o", srcs),
    copy    = if (length(srcs)) .dllFiles(objs) else character(0))
  if (!length(files))
    files <- Sys.glob(file.path(.dmodOutdir(), switch(build, compile = c("*.c", "*.cpp"),
                                                      link = "*.o", copy = "*.so")))
  if (build == "link" && (!length(files) || !all(file.exists(files))))
    stop("'link' needs the object files (.o) of the models: ",
         paste(basename(files[!file.exists(files)]), collapse = ", "), call. = FALSE)
  files
}

## The SLURM jobs of a manifest for squeue and scancel: by ID, else by name.
.slurmSelect <- function(m) {
  if (is.na(m$jobid)) paste0("-n ", shQuote(m$jobname), " -u \"$USER\"") else paste("-j", m$jobid)
}

## The state of every task: "done", "running" or "failed". A task waiting for
## its worker counts as running.
.taskStatus <- function(m) {
  worker <- m$taskWorker %||% if (m$slurm) seq_len(m$ntasks) else m$taskHost
  state <- rep("failed", m$ntasks)
  for (k in seq_along(m$hosts)) {
    host <- m$hosts[k]
    dir <- file.path(m$folder, k)
    probe <- if (m$slurm)
      paste0("ls; squeue -h -r ", .slurmSelect(m), " -o 'running %K' 2>/dev/null")
    else
      # a zombie is finished: a container whose init does not reap keeps it in the
      # process table, where it still answers signal 0
      paste0("ls; for f in slot_*.pid; do [ -f \"$f\" ] || continue; p=$(cat \"$f\"); ",
             "if [ -r /proc/$p/status ]; then grep -q \"^State:[[:space:]]*Z\" /proc/$p/status ",
             "|| echo alive \"$f\"; else kill -0 \"$p\" 2>/dev/null && echo alive \"$f\"; fi; done")
    out <- suppressWarnings(.onMachine(host, paste0(probe, "; true"), dir, m$wd,
                                       intern = TRUE, ignore.stderr = TRUE))
    mine <- which(m$taskHost == k)
    if (m$slurm) {
      running <- as.integer(sub("^running ", "", grep("^running [0-9]+$", out, value = TRUE)))
      state[mine[worker[mine] %in% running]] <- "running"
    } else {
      alive <- as.integer(sub("^alive slot_([0-9]+)\\.pid$", "\\1",
                              grep("^alive ", out, value = TRUE)))
      if (k %in% alive) state[mine] <- "running"
    }
    done <- as.integer(sub("^task_([0-9]+)\\.rds$", "\\1",
                           grep("^task_[0-9]+\\.rds$", out, value = TRUE)))
    state[intersect(mine, done)] <- "done"
  }
  state
}

## "<jobname>: k of n tasks done, r running, f failed".
.stateLine <- function(m, state) {
  n <- table(factor(state, c("done", "running", "failed")))
  paste0(m$jobname, ": ", n[["done"]], " of ", m$ntasks, " tasks done",
         if (n[["running"]]) paste0(", ", n[["running"]], " running"),
         if (n[["failed"]]) paste0(", ", n[["failed"]], " failed"))
}

## The states once no task runs, polled every `interval` seconds; each change is
## reported as a message.
.waitTasks <- function(m, interval) {
  last <- ""
  repeat {
    state <- .taskStatus(m)
    line <- .stateLine(m, state)
    if (!identical(line, last)) message(line)
    last <- line
    if (!any(state == "running")) return(state)
    Sys.sleep(interval)
  }
}

## The results of the tasks joined: parlists by c(), data frames by rbind(), any
## other results as a list. Tasks without a result are left out with a warning.
.combineResults <- function(out) {
  failed <- which(vapply(out, inherits, NA, "try-error"))
  missing <- which(vapply(out, is.null, NA))
  if (length(failed) || length(missing)) {
    ids <- function(i) paste0(paste(utils::head(i, 10L), collapse = ", "),
                              if (length(i) > 10L) ", ...")
    warning(length(failed) + length(missing), " of ", length(out), " tasks have no result",
            if (length(missing)) paste0("; not done: ", ids(missing)),
            if (length(failed)) paste0("; failed: ", ids(failed), ", the first with: ",
                                       conditionMessage(attr(out[[failed[1L]]], "condition"))),
            call. = FALSE)
  }
  ok <- out[!seq_along(out) %in% c(failed, missing)]
  if (!length(ok)) return(ok)
  if (all(vapply(ok, inherits, NA, "parlist"))) return(do.call(c, unname(ok)))
  if (all(vapply(ok, is.data.frame, NA)))
    return(tryCatch(do.call(rbind, unname(ok)), error = function(e) ok))
  ok
}

## The registry of submitted jobs, one manifest per job name and working directory.
.jobRegistry <- function() file.path(tools::R_user_dir("dMod2", "data"), "jobs")

.registryFile <- function(m)
  file.path(.jobRegistry(), paste0(m$jobname, "-", digest::digest(m$wd, algo = "xxhash32",
                                                                  serialize = FALSE), ".rds"))

## The manifest into the local job folder and the registry.
.saveManifest <- function(m) {
  saveRDS(m, file.path(m$wd, m$folder, "job.rds"))
  dir.create(.jobRegistry(), recursive = TRUE, showWarnings = FALSE)
  saveRDS(m, .registryFile(m))
}

.registeredJobs <- function() {
  files <- list.files(.jobRegistry(), pattern = "\\.rds$", full.names = TRUE)
  Filter(Negate(is.null), lapply(files, function(f) tryCatch(readRDS(f), error = function(e) NULL)))
}

## The manifest of the job `jobname`: the one in the working directory, else the
## one registered under that name.
.recoverManifest <- function(jobname, wd) {
  if (is.null(jobname)) stop("'recover' needs the 'jobname'", call. = FALSE)
  f <- file.path(wd, paste0(jobname, "_folder"), "job.rds")
  if (file.exists(f)) {
    m <- readRDS(f)
    m$wd <- wd
    return(m)
  }
  found <- Filter(function(m) identical(m$jobname, jobname), .registeredJobs())
  if (length(found) == 1L) return(found[[1L]])
  if (!length(found))
    stop("no job '", jobname, "' was submitted from this working directory or is ",
         "registered, see dMod2Jobs()", call. = FALSE)
  stop("several jobs '", jobname, "' are registered, submitted from ",
       paste(vapply(found, `[[`, "", "wd"), collapse = ", "),
       "; recover it from its working directory", call. = FALSE)
}

## An error when a job of this name still runs where the new one would replace
## its folders: in this working directory, or on one of the same remote nodes.
.refuseRunning <- function(jobname, wd, hosts) {
  f <- file.path(wd, paste0(jobname, "_folder"), "job.rds")
  same <- function(m) identical(m$jobname, jobname) &&
    any(m$hosts[!vapply(m$hosts, .isLocalhost, NA)] %in% hosts)
  old <- c(if (file.exists(f)) list(.recoverManifest(jobname, wd)),
           Filter(same, .registeredJobs()))
  old <- old[!duplicated(vapply(old, `[[`, "", "wd"))]
  for (m in old) {
    running <- sum(.taskStatus(m) == "running")
    if (running)
      stop("the job '", jobname, "' submitted from ", m$wd, " still runs ", running,
           " tasks; recover = TRUE returns its functions, and terminate() stops it",
           call. = FALSE)
  }
}

## The functions of the job described by the manifest `m`.
.jobFunctions <- function(m) {
  local <- file.path(m$wd, m$folder)

  check <- function() {
    state <- .taskStatus(m)
    cat(.stateLine(m, state), "\n", sep = "")
    finished <- !any(state == "running")
    failed <- sum(state == "failed")
    if (finished && failed)
      warning(failed, " of ", m$ntasks, " tasks failed; log(i) shows the output of task i",
              call. = FALSE)
    invisible(finished)
  }

  wait <- function(interval = 10) invisible(all(.waitTasks(m, interval) == "done"))

  get <- function(combine = TRUE, all = FALSE) {
    for (k in seq_along(m$hosts)) {
      if (.isLocalhost(m$hosts[k])) next
      what <- if (all) "." else "task_*.rds"
      .pullFiles(m$hosts[k], file.path(m$folder, k), what, m$wd)
    }
    out <- vector("list", m$ntasks)
    names(out) <- m$names
    for (i in seq_len(m$ntasks)) {
      f <- file.path(local, m$taskHost[i], paste0("task_", i, ".rds"))
      if (file.exists(f)) out[i] <- list(readRDS(f))
    }
    if (combine) .combineResults(out) else out
  }

  log <- function(i, lines = 20L) {
    i <- .wholeNumber(i, "i", null = FALSE)
    if (i > m$ntasks) stop("the job has ", m$ntasks, " tasks", call. = FALSE)
    k <- m$taskHost[i]
    out <- suppressWarnings(.onMachine(
      m$hosts[k], paste0("tail -n ", as.integer(lines), " task_", i, ".out 2>/dev/null; true"),
      file.path(m$folder, k), m$wd, intern = TRUE))
    cat(out, sep = "\n")
    invisible(out)
  }

  purge <- function(local = FALSE) {
    for (k in seq_along(m$hosts))
      .onMachine(m$hosts[k], paste0("cd .. && rm -rf ", k, " && cd .. && rmdir ", m$folder,
                                    " 2>/dev/null; true"),
                 file.path(m$folder, k), m$wd, ignore.stderr = TRUE)
    if (local) {
      unlink(file.path(m$wd, m$folder), recursive = TRUE)
      unlink(.registryFile(m))
    }
    invisible(NULL)
  }

  terminate <- function() {
    for (k in seq_along(m$hosts)) {
      cmd <- if (m$slurm) paste("scancel", .slurmSelect(m))
      else paste0("for f in slot_*.pid; do [ -f \"$f\" ] || continue; p=$(cat \"$f\"); ",
                  "kill -TERM -- -$p 2>/dev/null || { pkill -TERM -P $p; kill -TERM $p; }; ",
                  "done 2>/dev/null; true")
      .onMachine(m$hosts[k], cmd, file.path(m$folder, k), m$wd, ignore.stderr = TRUE)
    }
    invisible(NULL)
  }

  structure(list(check = check, wait = wait, get = get, log = log, purge = purge,
                 terminate = terminate),
            class = "dMod2Job", manifest = m)
}

## Where a job runs, from its manifest.
.jobWhere <- function(m) {
  if (!is.null(m$on)) return(m$on)
  if (m$slurm) paste("SLURM cluster", m$hosts) else paste(m$hosts, collapse = ", ")
}

#' @export
print.dMod2Job <- function(x, ...) {
  m <- attr(x, "manifest")
  cat("Job ", m$jobname, ": ", m$ntasks, if (m$ntasks == 1L) " task" else " tasks",
      " on ", .jobWhere(m), if (!is.na(m$jobid)) paste0(", SLURM job ", m$jobid), "\n",
      "check(), wait(), get(combine = TRUE, all = FALSE), log(i), purge(local = FALSE), ",
      "terminate()\n", sep = "")
  invisible(x)
}


#' Run R Code as Tasks on Other Processes, Machines or a SLURM Cluster
#'
#' Saves the objects named in `input`, copies them with the code and the model
#' files to the nodes given by `on` and runs the code there once per task, in
#' the background. Returns functions to check, wait for, collect and delete the
#' results. Needs a Unix-like system.
#'
#' The start is reported as a message, on a SLURM cluster with the job ID. A
#' build writes its output to `build.log` in the job folder; a failed build
#' shows its end in the error.
#'
#' @details The value of `expr` is the result of a task. In the code, `task`
#' is the task's element of `tasks` (`NULL` for a count), `taskId` its
#' position, `cores` the cores the task may use and `jobname` the job's name.
#' Parallel code such as [mstrust()] or [profile()][profile.objfn] takes
#' `cores = cores`. A task that raises an error returns the error as a
#' `try-error`. Task `i` starts with the seed `seed + i`.
#'
#' The tasks are dealt in turn to the workers of `on`, the nodes or the array
#' tasks of a SLURM job, and a worker runs its tasks one after the other, see
#' [localNode()]. A job is kept in `<jobname>_folder`, in the working directory
#' for `"localhost"` and in the home directory of a remote node. A job of the
#' same name that still runs there is an error. Every job is registered in
#' `tools::R_user_dir("dMod2", "data")` until `purge(local = TRUE)`, and
#' [dMod2Jobs()] lists them.
#'
#' The model files are those of the [obsfn], [parfn], [prdfn] and objective
#' objects in `input`; objects without them take every file of
#' `getOption("dMod2.outdir")`, else of the working directory.
#' \itemize{
#'   \item `build = "compile"`: the C/C++ sources are copied and compiled once
#'         per node into one shared object, against the cppDE installed there.
#'   \item `build = "link"`: the object files (`.o`) are copied and linked once
#'         per node; needs a compiler there compatible with the local one.
#'   \item `build = "copy"`: the shared objects are copied.
#'   \item `build = "auto"` (default): `"copy"` on the local machine,
#'         `"compile"` elsewhere.
#' }
#' After a build the [obsfn], [parfn] and [prdfn] objects and the objectives
#' use the new shared object.
#'
#' @param expr R code, evaluated once per task.
#' @param tasks A positive whole number of tasks, or a list or vector with one
#'   element per task, or a data frame with one row per task. `NULL` (default)
#'   is one task per worker.
#' @param on Where the tasks run: [localNode()] (default), [sshNodes()] or
#'   [slurmCluster()].
#' @param jobname Character, the name of the job and its folder. `NULL`
#'   (default) picks a random name.
#' @param input Character vector, the objects saved and copied to the nodes,
#'   looked up from the calling environment. Defaults to all objects of the
#'   calling environment.
#' @param build How the model files reach the nodes: `"auto"`, `"compile"`,
#'   `"link"` or `"copy"`, see Details.
#' @param libs Optional character vector of library paths put in front of the
#'   library search path on the nodes, for the tasks and the build. The remote
#'   shell expands them, so `~` is the remote home.
#' @param seed Integer, task `i` starts with `set.seed(seed + i)`. `NULL`
#'   (default) draws one.
#' @param wait Logical. `TRUE` waits for the tasks, returns the results of
#'   `get()` and deletes the job with `purge(local = TRUE)`; a job with failed
#'   tasks is kept. Default `FALSE`.
#' @param recover Logical. `TRUE` returns the functions of the job `jobname`,
#'   submitted from this working directory or registered, without submitting
#'   again. Default `FALSE`.
#' @param buildCores Number of parallel compiler processes of a build. `NULL`
#'   (default) uses the cores the node reports, at most 16.
#' @param buildBundle Number of generated sources a build puts into one
#'   translation unit, only for builds with more files than the shell accepts
#'   as arguments. `1` compiles one file at a time. Default 50.
#' @return With `wait = TRUE` the results of `get()`, else an object of class
#'   `dMod2Job`, a list of functions:
#' \describe{
#'   \item{`check()`}{reports how many tasks are done, running or failed and
#'     returns `TRUE` once no task runs any more, with a warning when some
#'     failed. A task fails when its process ends without a result.}
#'   \item{`wait(interval = 10)`}{polls every `interval` seconds until no task
#'     runs, reports each change as a message and returns `TRUE` when every task
#'     is done.}
#'   \item{`get(combine = TRUE, all = FALSE)`}{returns the results. With
#'     `combine = TRUE` parlists are joined by `c()` and data frames such as
#'     parframes by `rbind()`, other results come as a list; tasks without a
#'     result are left out with a warning. `combine = FALSE` returns one entry
#'     per task, named like `tasks`, `NULL` for a task not done. Results of
#'     remote nodes are copied into `<jobname>_folder/<k>`, with `all = TRUE`
#'     together with every other file of the job folder.}
#'   \item{`log(i, lines = 20)`}{prints the last lines of the output of task
#'     `i`.}
#'   \item{`purge(local = FALSE)`}{deletes the job folders on the nodes, with
#'     `local = TRUE` also the local `<jobname>_folder` and the registry entry.}
#'   \item{`terminate()`}{stops the tasks.}
#' }
#' @seealso [localNode()], [slurmTemplate()], [dMod2Jobs()], [detectFreeCores()]
#' @example inst/examples/distributedComputing.R
#' @export
distributedComputing <- function(expr, tasks = NULL, on = localNode(), jobname = NULL,
                                 input = ls(parent.frame()),
                                 build = c("auto", "compile", "link", "copy"), libs = NULL,
                                 seed = NULL, wait = FALSE, recover = FALSE,
                                 buildCores = NULL, buildBundle = 50) {

  wd <- getwd()
  if (recover) return(.jobFunctions(.recoverManifest(jobname, wd)))

  expr <- substitute(expr)
  env <- parent.frame()
  if (!inherits(on, "dMod2Target"))
    stop("'on' must be localNode(), sshNodes() or slurmCluster()", call. = FALSE)
  build <- match.arg(build)
  slurm <- inherits(on, "dMod2Slurm")
  if (slurm) {
    hosts <- on$login
    cores <- stats::setNames(on$cores, hosts)
    jobscript <- on$jobscript %||% .slurmJobscript(hosts)
    if (!file.exists(jobscript))
      stop("the job script ", jobscript, " does not exist", call. = FALSE)
  } else {
    hosts <- on$hosts
    free <- detectFreeCores(hosts)
    short <- !is.na(on$cores) & on$cores > free
    if (any(short))
      warning(paste0(hosts[short], " has ", free[short], " free cores; its tasks get ",
                     free[short], " instead of ", on$cores[short], collapse = "; "), call. = FALSE)
    cores <- stats::setNames(as.integer(ifelse(is.na(on$cores) | short, free, on$cores)), hosts)
    on$cores <- cores
  }
  tasks <- .asTasks(tasks %||% if (slurm) on$array %||% 1L else length(hosts))
  n <- length(tasks)
  workers <- if (slurm) min(on$array %||% n, n) else length(hosts)
  if (build == "auto") build <- if (all(vapply(hosts, .isLocalhost, NA))) "copy" else "compile"
  if (is.null(jobname))
    jobname <- paste0("job_", paste(sample(c(0:9, letters), 6, replace = TRUE), collapse = ""))
  if (is.null(seed)) seed <- sample.int(1e6, 1L)
  .refuseRunning(jobname, wd, hosts)

  # The worker of every task, dealt in turn: a node, or an array task of SLURM
  worker <- (seq_len(n) - 1L) %% workers + 1L
  m <- list(jobname = jobname, wd = wd, folder = paste0(jobname, "_folder"),
            on = .targetLabel(on), cores = cores, hosts = hosts, slurm = slurm, ntasks = n,
            names = names(tasks), taskHost = if (slurm) rep(1L, n) else worker,
            taskWorker = worker, workers = workers, jobid = NA_character_,
            submitted = Sys.time())
  local <- file.path(wd, m$folder)
  unlink(local, recursive = TRUE)
  dir.create(local)

  # The stage: workspace, tasks, task script, build script and model files
  stage <- tempfile("dMod_stage_")
  dir.create(stage)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  input <- setdiff(input, ".Random.seed")
  .saveWorkspace(input, file.path(stage, "workspace.RData"), envir = env)
  saveRDS(tasks, file.path(stage, "tasks.rds"))

  objs <- mget(input, envir = env, inherits = TRUE)
  files <- .jobFiles(objs, build)
  shared <- NULL
  buildCmd <- ""
  if (build %in% c("compile", "link") && length(files)) {
    shared <- paste0(jobname, "_shared_object")
    info <- .remoteBuildInfo(list2env(objs))
    writeLines(basename(files), file.path(stage, "files.txt"))
    writeLines(.remoteBuildScript(
      files = basename(files), output = paste0(shared, ".so"),
      compileArgs = info$compileArgs, needsCVODE = info$needsCVODE, needsKLU = info$needsKLU,
      needsLapack = info$needsLapack, link = build == "link",
      cxx = any(grepl("\\.cpp$", files)), cores = buildCores, filelist = "files.txt",
      bundle = buildBundle),
      file.path(stage, "build.sh"))
    # the build writes to build.log, which shows its end only when it fails
    buildCmd <- "{ bash build.sh > build.log 2>&1 || { tail -n 40 build.log; false; }; } && "
  }
  if (length(files)) file.symlink(normalizePath(files), file.path(stage, basename(files)))
  writeLines(.taskScript(expr, jobname, seed, shared), file.path(stage, "job.R"))

  unpack <- paste0("if [ -f workspace.RData.zst ]; then zstd -dqf workspace.RData.zst && ",
                   "rm -f workspace.RData.zst; fi && ")
  setup <- paste0(.remoteLibs(libs), "export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; ")

  if (slurm) {
    filled <- .fillJobscript(jobscript, list(
      jobname = jobname, ntasks = workers, cores = cores[[1L]], memPerCore = on$memPerCore,
      walltime = on$walltime))
    # array task j runs the tasks j, j + workers, ... one after the other
    filled$script[grepl("^\\s*\\{\\{\\s*run\\s*\\}\\}\\s*$", filled$script)] <- paste0(
      setup, "for i in $(seq \"$SLURM_ARRAY_TASK_ID\" ", workers, " ", n, "); do ",
      "Rscript job.R \"$i\" ", cores[[1L]], " > \"task_$i.out\" 2>&1; done")
    writeLines(filled$script, file.path(stage, "job.sh"))
    writeLines(filled$env, file.path(stage, "env.sh"))
    dir <- file.path(m$folder, 1L)
    .pushStage(hosts, stage, dir, wd)
    .saveManifest(m)
    opts <- c(if (!is.null(on$partition)) paste0("--partition=", on$partition),
              if (!is.null(on$account)) paste0("--account=", on$account),
              if (!is.null(on$maxRunning)) paste0("--array=1-", workers, "%", on$maxRunning),
              on$sbatch)
    out <- .onMachine(hosts, paste0(unpack, ". ./env.sh && ", setup, buildCmd,
                                    paste(c("sbatch", shQuote(opts), "job.sh"), collapse = " ")),
                      dir, wd, intern = TRUE)
    submitted <- grep("Submitted batch job", out, value = TRUE)
    id <- regmatches(submitted, regexpr("[0-9]+$", submitted))
    if (!length(id))
      stop("the submission on ", hosts, " reported no job ID; if the job runs, recover = TRUE ",
           "returns its functions:\n", paste(out, collapse = "\n"), call. = FALSE)
    m$jobid <- id[1L]
    .saveManifest(m)
    message(jobname, ": ", .count(n, "task"), " submitted to ", hosts, " as SLURM job ", m$jobid,
            ", ", .count(workers, "array task"), " of ", cores[[1L]], " cores")
  } else {
    run <- function(id) paste0("Rscript job.R ", id, " ", cores[[worker[id]]],
                               " > task_", id, ".out 2>&1")
    for (k in unique(worker))
      writeLines(c("#!/bin/sh", setup, vapply(which(worker == k), run, "")),
                 file.path(stage, paste0("slot_", k, ".sh")))
    .saveManifest(m)
    for (k in unique(worker)) {
      dir <- file.path(m$folder, k)
      .pushStage(hosts[k], stage, dir, wd)
      out <- suppressWarnings(.onMachine(hosts[k], paste0(
        unpack, setup, buildCmd, "{ if command -v setsid >/dev/null; then S=setsid; else S=; fi; ",
        "nohup $S sh slot_", k, ".sh > /dev/null 2>&1 < /dev/null & echo $! > slot_", k,
        ".pid; }"), dir, wd, intern = TRUE))
      if (!is.null(attr(out, "status")))
        stop("the start on ", hosts[k], " failed:\n", paste(out, collapse = "\n"), call. = FALSE)
    }
    message(jobname, ": ", .count(n, "task"), " started on ", .targetLabel(on))
  }
  job <- .jobFunctions(m)
  if (!wait) return(job)

  .waitTasks(m, 1)
  out <- job$get(combine = FALSE)
  result <- suppressWarnings(.combineResults(out))
  failed <- sum(vapply(out, function(r) is.null(r) || inherits(r, "try-error"), NA))
  if (!failed) {
    job$purge(local = TRUE)
    return(result)
  }
  warning(failed, " of ", n, " tasks failed; the job is kept: recover = TRUE returns its ",
          "functions, and log(i) shows the output of task i", call. = FALSE)
  result
}

#' Jobs of distributedComputing()
#'
#' Lists the jobs submitted with [distributedComputing()] that were not deleted
#' with `purge(local = TRUE)`, from every working directory. Their functions
#' come with `distributedComputing(NULL, jobname = , recover = TRUE)`.
#'
#' @return A data frame with one row per job: `jobname`, `on` (where it runs),
#'   `tasks`, `jobid` (of SLURM), `submitted`, `wd` (the working directory it
#'   was submitted from) and `local` (whether its local job folder is still
#'   there).
#' @seealso [distributedComputing()]
#' @examples
#' dMod2Jobs()
#' @export
dMod2Jobs <- function() {
  jobs <- .registeredJobs()
  chr <- function(name) vapply(jobs, function(m) as.character(m[[name]] %||% NA), "")
  out <- data.frame(
    jobname = chr("jobname"),
    on = vapply(jobs, .jobWhere, ""),
    tasks = vapply(jobs, function(m) as.integer(m$ntasks), 0L),
    jobid = chr("jobid"),
    submitted = .POSIXct(vapply(jobs, function(m) as.numeric(m$submitted %||% NA), 0)),
    wd = chr("wd"),
    local = vapply(jobs, function(m) dir.exists(file.path(m$wd, m$folder)), NA))
  out[order(out$submitted, decreasing = TRUE), , drop = FALSE]
}

#' Job script templates of SLURM clusters
#'
#' Sets up the job script template of a SLURM cluster, once per login node.
#' The templates are kept in `tools::R_user_dir("dMod2", "config")` and stay
#' across sessions; [distributedComputing()] on a [slurmCluster()] takes the one
#' named like its login node.
#'
#' A template is a job script with the placeholders `{{ jobname }}`,
#' `{{ ntasks }}` (the number of array tasks), `{{ cores }}`, `{{ memPerCore }}`,
#' `{{ walltime }}` and one line `{{ run }}`, which runs the tasks of the array
#' task `SLURM_ARRAY_TASK_ID`. The lines before `{{ run }}` that are not
#' comments set up the environment, for example with `module load`; the build on
#' the login node runs in it as well. A job script of the cluster becomes a
#' template by these placeholders. dMod2 ships one for the bwForCluster Helix,
#' `"helix"`.
#'
#' @param login The login node, as in [slurmCluster()]. `NULL` (default) lists
#'   the templates set up.
#' @param from The template to start from: the name of one dMod2 ships, such
#'   as `"helix"`, or the path of a job script. Defaults to the shipped
#'   template named `login`.
#' @param edit Logical, opens the template in the editor. Default
#'   `interactive()`.
#' @param overwrite Logical. `TRUE` replaces a template set up before; with
#'   `FALSE` (default) it is kept.
#'
#' @return The path of the template, invisibly. Without `login` the paths of
#'   the templates set up, named by login node.
#'
#' @seealso [distributedComputing()], [slurmCluster()]
#' @examples
#' \dontrun{
#' # once per login node, from a job script of the cluster
#' slurmTemplate("mycluster", from = "~/mycluster.sh")
#' # the template dMod2 ships for the bwForCluster Helix
#' slurmTemplate("helix")
#' slurmTemplate()
#' }
#' @export
slurmTemplate <- function(login = NULL, from = login, edit = interactive(),
                          overwrite = FALSE) {
  dir <- .slurmTemplateDir()
  if (is.null(login)) {
    files <- list.files(dir, pattern = "\\.sh$", full.names = TRUE)
    return(stats::setNames(files, sub("\\.sh$", "", basename(files))))
  }
  if (!is.character(login) || length(login) != 1L || !nzchar(login) ||
      grepl("/", login, fixed = TRUE))
    stop("'login' must name one login node", call. = FALSE)
  dest <- file.path(dir, paste0(login, ".sh"))
  if (!file.exists(dest) || overwrite) {
    shipped <- .slurmShipped()
    src <- if (from %in% names(shipped)) shipped[[from]] else from
    if (!file.exists(src))
      stop(if (missing(from)) paste0("dMod2 ships no template for '", login, "'; give ",
                                     "the job script to start from as 'from'")
           else paste0("'", from, "' is neither a job script nor a template dMod2 ships"),
           ", see ?slurmTemplate", call. = FALSE)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(src, dest, overwrite = TRUE))
      stop("cannot write ", dest, call. = FALSE)
  }
  if (edit) utils::file.edit(dest)
  invisible(dest)
}


## .sanitizeCores ----------------------------------------------------------------

# Split a two-axis core budget. `cores` is a single number (outer axis only)
# or a named vector such as c(fits = 10, conditions = 5). The product is
# reported, not capped, over-subscribing is sometimes deliberate.
.splitCores <- function(cores, outer = "outer") {
  if (is.null(cores)) return(list(outer = 1L, conditions = NULL))
  cores <- as.integer(cores)
  if (anyNA(cores) || any(cores < 1L))
    stop("cores must be positive integers.", call. = FALSE)
  if (length(cores) == 1L) return(list(outer = cores, conditions = NULL))
  if (length(cores) > 2L)  stop("cores must have length 1 or 2.", call. = FALSE)

  nm <- names(cores)
  o <- if (!is.null(nm) && outer %in% nm) cores[[outer]] else cores[[1L]]
  k <- if (!is.null(nm) && "conditions" %in% nm) cores[["conditions"]] else cores[[2L]]

  avail <- .availableCores()
  if (o * k > avail)
    warning(sprintf("cores: %d (%s) x %d (conditions) = %d exceeds the %d available.",
                    o, outer, k, o * k, avail), call. = FALSE)
  list(outer = o, conditions = k)
}


# Cores this process may use: the affinity mask, cgroup and scheduler
# allocation (SLURM, PBS, ...), not the cores of the node.
.availableCores <- function() as.integer(parallelly::availableCores())

# Caps each element of a core vector at the cores available, names kept.
.sanitizeCores <- function(cores) {
  cores[] <- pmin(as.integer(cores), .availableCores())
  cores
}




## .parallelLapply ---------------------------------------------------------------

# Cross-platform parallel-apply. Unix forks via doParallel; Windows uses a
# PSOCK cluster with explicit library-path + variable export. Driven through
# foreach::%dopar% so the caller body is the same on both platforms.
.parallelLapply <- function(X, FUN, cores = 1L,
                            extraExports = character(0),
                            envir = parent.frame(),
                            coresConditions = NULL,
                            strategy = c("auto", "fork", "psock")) {

  strategy <- match.arg(strategy)
  cores <- max(1L, as.integer(cores))
  if (cores == 1L && !is.null(coresConditions)) {
    old <- options(dMod2.cores = as.integer(coresConditions))
    on.exit(options(old), add = TRUE)
  }
  if (cores == 1L)
    return(lapply(X, FUN))

  # A forked worker cannot split a condition axis: cppDE's batch runs serially
  # inside a fork, so an inner axis selects PSOCK.
  wants_inner <- !is.null(coresConditions) && coresConditions > 1L
  use_psock <- switch(strategy,
                      psock = TRUE,
                      fork  = FALSE,
                      auto  = Sys.info()[["sysname"]] == "Windows" || wants_inner)

  if (use_psock) {
    cluster <- parallel::makeCluster(cores)
    on.exit({
      parallel::stopCluster(cluster)
      doParallel::stopImplicitCluster()
    }, add = TRUE)
    doParallel::registerDoParallel(cl = cluster)
    parallel::clusterCall(cl = cluster, function(x) .libPaths(x), .libPaths())
    parallel::clusterEvalQ(cluster, suppressMessages(library(dMod2)))
    if (length(extraExports) > 0L)
      parallel::clusterExport(cluster, envir = envir,
                              varlist = extraExports)
  } else {
    doParallel::registerDoParallel(cores = cores)
    on.exit(doParallel::stopImplicitCluster(), add = TRUE)
  }

  i <- NULL  # silence R CMD check NSE warning
  loaded_packages <- .packages()
  kcond <- if (is.null(coresConditions)) NULL else as.integer(coresConditions)
  out <- foreach::foreach(i = seq_along(X),
                          .packages = loaded_packages,
                          .options.multicore = list(preschedule = FALSE)) %dopar% {
    if (!is.null(kcond)) options(dMod2.cores = kcond)
    FUN(X[[i]])
  }
  out
}



`%dopar%` <- foreach::`%dopar%`
