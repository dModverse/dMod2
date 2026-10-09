#' Detect Number of Free Cores
#'
#' Estimates the free cores as the number of cores minus the 1-minute load
#' average, at least 1. Works on Linux and macOS, locally or on remote
#' machines via `ssh`. On Windows it returns 1 with a warning. Under
#' `R CMD check` with `_R_CHECK_LIMIT_CORES_` set, the result is at most 2.
#' @param machine Character vector of SSH hosts, e.g. `"user@@host"`. `NULL`
#'   (default) for the local machine.
#' @return Numeric vector of free cores, one per machine, with attributes
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
    output <- lapply(machine, .getLoadAndCores)
    freeCores <- vapply(output, `[[`, numeric(1), "free")
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
  ## nproc honours OMP_NUM_THREADS, which runbg() sets to 1 for the job itself
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


#' Run an R Expression in the Background
#'
#' Saves the objects named in `input` to a workspace, copies it with `scp` to
#' each machine, runs the expression there with `Rscript` and returns
#' functions to collect the results. Requires a Unix-like system and `ssh`
#' key authentication on every machine.
#'
#' @details The value of the last expression in `...` is the result of a
#' machine. `get()` returns the list of these results.
#'
#' Build-related files are handled as follows:
#' \itemize{
#'   \item `compile = TRUE`: C/C++ source files are transferred and compiled
#'         remotely via `R CMD SHLIB`.
#'   \item `link = TRUE`: object files (`.o`) are transferred and linked
#'         remotely.
#'   \item Both `FALSE` (default): shared objects (`.so`) are copied.
#' }
#' With `compile` or `link`, the `modelname` of [obsfn], [parfn] and [prdfn]
#' objects in the workspace is set to the newly built shared object. The
#' remote build needs cppDE installed for the R on the `PATH` of the remote
#' machine.
#' @param ... R code.
#' @param machine Character vector of SSH hosts, e.g. `"localhost"`,
#'   `"user@@host"` or `c("localhost", "localhost")`. Defaults to
#'   `"localhost"`.
#' @param filename Character, the base name of the job files. `NULL` (default)
#' picks a random name.
#' @param input Character vector, the objects of the global environment that are
#' saved and copied to the remote machines. Defaults to all.
#' @param compile Logical. `TRUE` transfers the C/C++ sources (`.c`, `.cpp`) of
#' the working directory and compiles them on each machine into one shared
#' object; overrides `link`. Call it from `getOption("dMod.outdir")` when that
#' is set. Defaults to `FALSE`.
#' @param link Logical. `TRUE` transfers the object files (`.o`) of the working
#' directory and links them on each machine; an error is raised if there are
#' none. Needs a remote compiler compatible with the local one. Ignored if
#' `compile = TRUE`. Defaults to `FALSE`.
#' @param buildCores Number of parallel compiler processes of the build on each
#' remote machine. `NULL` (default) uses the cores the machine reports, at most
#' 16.
#' @param buildBundle Number of generated sources the remote build puts into one
#' translation unit, as in [distributedComputing()]. Only used for builds with
#' more files than the shell accepts as arguments. `1` compiles one file at a
#' time. Defaults to 50.
#' @param wait Logical. `TRUE` waits for the job and returns its results as
#' `get()` does; if the local result files exist already, they are loaded
#' without running the job. `FALSE` (default) returns at once; the results are
#' fetched by `get()`.
#' @param recover Logical. `TRUE` returns the functions `check()`, `get()`,
#' `purge()` and `terminate()` of an earlier job without starting it again,
#' e.g. after a crashed session. Needs the `filename` of that job. Defaults to
#' `FALSE`.
#' @param walltime Optional character, the maximum runtime as `"HH:MM:SS"`,
#' after which the job is terminated. `NULL` (default) sets no limit.
#' @param libs Optional character vector of library paths put in front of the
#' remote library search path, for the job and the remote build. The remote
#' shell expands them, so `~` is the remote home.
#' @return With `wait = TRUE` the results, else a list of functions:
#' \describe{
#'   \item{`check()`}{reports whether the results are ready and returns
#'     `TRUE` or `FALSE`.}
#'   \item{`get()`}{copies the result files to the working directory and
#'     returns the results, a list named by machine.}
#'   \item{`purge()`}{deletes the job folders on the remote machines and the
#'     job files in the local working directory.}
#'   \item{`terminate()`}{kills the processes of the job on the remote
#'     machines.}
#' }
#' @seealso [distributedComputing()]
#' @export
#' @examples
#' \dontrun{
#' job <- runbg({
#'   M <- matrix(rnorm(1e2), 10, 10)
#'   solve(M)
#' }, machine = c("user@@host", "user@@host"), filename = "job1")
#' job$check()
#' result <- job$get()
#' job$purge()
#'
#' # Recover the functions of a running job, e.g. after a crashed session
#' job <- runbg(NULL, machine = c("user@@host", "user@@host"),
#'              filename = "job1", recover = TRUE)
#' job$get()
#' }
runbg <- function(..., machine = "localhost", filename = NULL, input = ls(.GlobalEnv), compile = FALSE, link = FALSE, buildCores = NULL, buildBundle = 50, wait = FALSE, recover = FALSE, walltime = NULL, libs = NULL) {
  
  expr <- as.expression(substitute(...))
  nmachines <- length(machine)
  
  # compile takes precedence over link
  if (compile) link <- FALSE
  
  # Generate a random filename if none is provided
  if (is.null(filename))
    filename <- paste0("tmp_", paste(sample(c(0:9, letters), 5, replace = TRUE), collapse = ""))
  
  filename0 <- filename
  filename <- paste(filename, 1:nmachines, sep = "_")
  
  # Initialize output list
  out <- structure(vector("list", 4), names = c("check", "get", "purge", "terminate"))
  
  # Check whether results are ready on all machines
  out[[1]] <- function() {
    
    check.out <- sapply(1:nmachines, function(m) length(suppressWarnings(
      system(paste0("ssh ", machine[m], " ls ", filename[m], "_folder/ | grep -x ", filename[m], "_result.RData"), 
             intern = TRUE))))
    
    if (all(check.out > 0)) {
      cat("Result is ready!\n")
      return(TRUE)
    } else if (any(check.out > 0)) {
      cat("Result from machines", paste(which(check.out > 0), collapse = ", "), "are ready.")
      return(FALSE)
    } else {
      cat("Not ready!\n") 
      return(FALSE)
    }
    
  }
  
  # Fetch result files from remote machines and load into workspace
  out[[2]] <- function() {
    
    result <- structure(vector(mode = "list", length = nmachines), names = machine)
    for (m in 1:nmachines) {
      .runbgOutput <- NULL
      system(paste0("scp ", machine[m], ":", filename[m], "_folder/", filename[m], "_result.RData ./"), ignore.stdout = TRUE, ignore.stderr = TRUE)
      check <- try(load(file = paste0(filename[m], "_result.RData")), silent = TRUE) 
      if (!inherits(check, "try-error")) result[[m]] <- .runbgOutput
    }
    result
  }
  
  # Remove temporary folders and files on remote machines and locally
  out[[3]] <- function() {
    
    for (m in 1:nmachines) {
      folder_exists <- suppressWarnings(
        system(paste0("ssh ", machine[m], " '[ -d ", filename[m], "_folder ] && echo 1 || echo 0'"), 
               intern = TRUE)
      )
      if (folder_exists == "1") {
        system(paste0("ssh ", machine[m], " rm -r ", filename[m], "_folder"))
      }
      
      rout_exists <- suppressWarnings(
        system(paste0("ssh ", machine[m], " '[ -f ", filename[m], ".Rout ] && echo 1 || echo 0'"), 
               intern = TRUE)
      )
      if (rout_exists == "1") {
        system(paste0("ssh ", machine[m], " rm ", filename[m], ".Rout"))
      }
    }
    
    local_files <- list.files(pattern = paste0(filename0, ".*"))
    if (length(local_files) > 0) {
      system(paste0("rm ", filename0, "*"))
    }
  }
  
  # Kill all processes associated with this job on the remote machines
  out[[4]] <- function() {
    for (m in 1:nmachines) {
      pids <- suppressWarnings(
        system(paste0("ssh ", machine[m], 
                      " 'ps aux | grep \"", filename[m], "\" | grep -v grep'"), 
               intern = TRUE)
      )
      
      if (length(pids) > 0) {
        running_pids <- sapply(strsplit(pids, "\\s+"), function(x) {
          if (grepl("R", x[8])) x[2] else NULL
        })
        running_pids <- running_pids[!sapply(running_pids, is.null)]
        
        if (length(running_pids) > 0) {
          system(paste0("ssh ", machine[m], " 'kill ", paste(running_pids, collapse = " "), "'"))
          cat("Terminated", length(running_pids), "running processes on", machine[m], "\n")
        } else {
          cat("No running processes found on", machine[m], "\n")
        }
      }
    }
  }
  
  # Recover control functions without re-submitting the job
  if (recover) return(out)
  
  # If result files already exist locally and wait = TRUE, load them directly
  resultfile <- paste(filename, "result.RData", sep = "_")
  if (all(file.exists(resultfile)) & wait) {
    
    result <- structure(vector(mode = "list", length = nmachines), names = machine)
    for (m in 1:nmachines) {
      .runbgOutput <- NULL
      load(file = resultfile[m])
      result[[m]] <- .runbgOutput
    }
    return(result)
  }
  
  # Save current workspace to be transferred to remote machines
  save(list = input, file = paste0(filename0, ".RData"), envir = .GlobalEnv)
  
  # The remote script attaches dMod2 only, not the submitting session's
  # packages; anything else belongs in the expression.
  pack <- "library(dMod2)"
  
  output <- ".runbgOutput"
  
  # The portable, model-specific flags and the file list; the rest is resolved
  # on the remote machine by the script .remoteBuildScript() writes.
  buildinfo <- list(compileArgs = "", needsCVODE = FALSE, needsKLU = FALSE,
                    needsLapack = FALSE)
  has_cxx <- FALSE
  if (compile || link) {

    buildinfo <- .remoteBuildInfo()
    has_cxx <- length(list.files(pattern = glob2rx("*.cpp"))) > 0

    if (compile) {
      sourcefiles <- c(list.files(pattern = glob2rx("*.c")),
                       list.files(pattern = glob2rx("*.cpp")))
    } else {
      object_files <- Sys.glob("*.o")
      if (length(object_files) == 0)
        stop("No .o files found for linking! You must compile first.")
      sourcefiles <- object_files
      message("runbg(): linking pre-compiled .o files remotely. This requires the ",
              "remote toolchain to be ABI-compatible with the local one; if the ",
              "link or dyn.load() fails, resubmit with compile = TRUE.")
    }

  }
  
  if (compile || link) {
    # R code to load the newly built shared object and update modelnames of known function objects
    objfns <- 'obj.fns <- ls()[sapply(ls(), function(nm) inherits(get(nm, envir=.GlobalEnv), c("obsfn", "parfn", "prdfn", "objfn")))]'
    setmn <- sprintf('for (o in obj.fns) eval(parse(text=paste0("modelname(", o, ") <- \'%s\'")))', paste0(filename0, "_shared_object"))
    load_so <- paste0("dyn.load('", filename0, "_shared_object.so')")
  } else {
    objfns <- NULL
    setmn <- NULL
    load_so <- NULL
  }
  
  # Assemble the R script to be executed on each remote machine
  program <- lapply(1:nmachines, function(m) paste(
    c(
      pack,
      paste0("setwd('~/", filename[m], "_folder')"),
      "rm(list = ls())",
      if (!is.null(walltime)) paste0("Sys.setenv(R_TIMEOUT='", walltime, "')"),
      paste0("load('", filename0, ".RData')"),
      objfns,
      setmn,
      if (!is.null(load_so)) {
        load_so
      } else {
        c("files <- list.files(pattern = '\\\\.so$')", "for (f in files) dyn.load(f)")
      },
      paste0(".node <- ", m),
      if (!is.null(walltime)) {
        paste0(".runbgOutput <- try(tools::pskill(Sys.getpid(), tools::SIGALRM, ", walltime, "); ", as.character(expr), ")")
      } else {
        paste0(".runbgOutput <- try(", as.character(expr), ")")
      },
      paste0("save(", output ,", file = '", filename[m], "_result.RData')")
    ),
    collapse = "\n"
  ))
  
  for (m in 1:nmachines) {
    
    # Write the R script for this machine
    cat(program[[m]], file = paste0(filename[m], ".R"))
    
    # Create remote working directory
    system(paste0("ssh ", machine[m], " mkdir -p ", filename[m], "_folder/"), 
           ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    # Clear any leftover files from previous runs
    system(paste0("ssh ", machine[m], " rm -r ", filename[m], "_folder/*"), 
           ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    # Transfer workspace
    system(paste0("scp ", getwd(), "/", filename0, ".RData* ", machine[m], ":", filename[m], "_folder/"))
    
    # Transfer R script
    system(paste0("scp ", getwd(), "/", filename[m], ".R* ", machine[m], ":", filename[m], "_folder/"))
    
    # Only the payload the requested mode consumes. Globs are expanded by the
    # local shell, so the command string stays short at any file count.
    payload <- if (compile) c("*.c", "*.cpp") else if (link) "*.o" else "*.so"
    payload <- payload[vapply(payload, function(g) length(Sys.glob(g)) > 0L, logical(1))]
    if (length(payload))
      system(paste0("scp ", paste0(getwd(), "/", payload, collapse = " "), " ",
                    machine[m], ":", filename[m], "_folder/"),
             ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    
    # Write and transfer compile script per machine, then build the remote command
    compile_cmd <- ""
    if (compile || link) {
      ## The build script reads its inputs from this list, not from its own
      ## command line.
      filelist_file <- paste0(filename[m], "_files.txt")
      writeLines(sourcefiles, filelist_file)
      compile_script_content <- .remoteBuildScript(
        files       = sourcefiles,
        output      = paste0(filename0, "_shared_object.so"),
        compileArgs = buildinfo$compileArgs,
        needsCVODE  = buildinfo$needsCVODE,
        needsKLU    = buildinfo$needsKLU,
        needsLapack = buildinfo$needsLapack,
        link        = link,
        cxx         = has_cxx,
        cores       = buildCores,
        workdir     = paste0(filename[m], "_folder"),
        filelist    = filelist_file,
        bundle      = buildBundle
      )
      compile_script_file <- paste0(filename[m], "_compile.sh")
      cat(compile_script_content, file = compile_script_file)
      system(paste0("scp ", getwd(), "/", compile_script_file, " ", getwd(), "/", filelist_file,
                    " ", machine[m], ":", filename[m], "_folder/"))
      compile_cmd <- paste0("bash ", filename[m], "_folder/", compile_script_file)
    }
    
    # Execute: compile/link if requested, then run the R script
    # OMP_NUM_THREADS=1 and MKL_NUM_THREADS=1 ensure single-threaded execution per job
    system(paste0(
      "ssh ", machine[m], 
      " '", .remoteLibs(libs), "export OMP_NUM_THREADS=1 && export MKL_NUM_THREADS=1",
      if (nzchar(compile_cmd)) paste0(" && ", compile_cmd),
      " && R CMD BATCH --vanilla ", filename[m], "_folder/", filename[m], ".R'"
    ), intern = FALSE, wait = wait)
  }
  
  if (!wait) return(out)
  result <- out$get()
  out$purge()
  result
}



## ---- HPC/SLURM distributed computing --------------------------------------

#' Run R Code on a Remote HPC System with SLURM
#'
#' Saves the workspace, copies it with generated R and shell scripts to a
#' remote system via `ssh`, and submits the R code there as a SLURM array job.
#' Returns functions to check, collect and delete the results.
#'
#' @details The value of `...` on each array task is its result. `get()`
#' returns the list of results, one per task.
#'
#' Repetitions of the same code are requested with `nRep`. Code that varies
#' between tasks reads the variables `var_1`, `var_2`, ..., whose values on
#' each task are taken from the corresponding vectors of `varValues`.
#'
#' The workspace is compressed with `zstd` if it is available locally. The job
#' script loads the environment modules `compiler/gnu/13.3` and `math/R`; the
#' remote build needs cppDE installed for that R. A job is only submitted if
#' the remote build succeeds. Without `ssh` key authentication, `sshPasswd` is
#' passed to `sshpass`, which must be installed locally.
#'
#' @param ... R code to be executed remotely. Variables that change between
#'   tasks are named `var_i`, see Details. `mem_per_core`, `ssh_passwd`,
#'   `var_values`, `no_rep`, `purge_local` and `custom_folders` are deprecated
#'   names of `memPerCore`, `sshPasswd`, `varValues`, `nRep`, `purgeLocal` and
#'   `customFolders`.
#' @param jobname Character, unique name of the job. A job with the same name
#'   is overwritten.
#' @param partition SLURM partition. Defaults to `"single"`.
#' @param cores Number of cores per node. Defaults to 16.
#' @param nodes Number of nodes per task. Defaults to 1.
#' @param memPerCore Memory per core in GB. Defaults to 2.
#' @param walltime Maximum runtime per task as `"hh:mm:ss"`. Defaults to
#'   `"01:00:00"`.
#' @param sshPasswd Optional password for `sshpass`. `NULL` (default) uses
#'   `ssh` key authentication, which is recommended.
#' @param machine SSH address of the remote system, e.g. `"user@@host"`.
#'   Defaults to `"cluster"`.
#' @param varValues List of vectors, one per variable `var_i`, all of the
#'   same length: the number of array tasks. Exclusive with `nRep`.
#' @param nRep Integer, number of repetitions. Exclusive with `varValues`.
#' @param recover Logical. `TRUE` submits nothing and returns the functions
#'   for a job submitted earlier under `jobname`. Default `FALSE`, which
#'   submits the job.
#' @param purgeLocal Logical, the default of the argument of the same name of
#'   `purge()`. Defaults to `FALSE`.
#' @param compile Logical. `TRUE` transfers the C/C++ sources (`*.c`,
#'   `*.cpp`) of the working directory and compiles them on the remote system
#'   into one shared object; overrides `link`. Call it from
#'   `getOption("dMod.outdir")` when that is set. Defaults to `FALSE`.
#' @param link Logical. `TRUE` transfers the object files (`*.o`) of the working
#'   directory and links them on the remote system; an error is raised if
#'   there are none. Needs a remote compiler compatible with the local one.
#'   Ignored if `compile = TRUE`. Defaults to `FALSE`.
#' @param buildCores Number of parallel compiler processes of the remote build,
#'   which runs on the login node before submission. `NULL` (default) uses the
#'   cores the login node reports, at most 16.
#' @param buildBundle Number of generated sources the remote build puts into one
#'   translation unit. Only used for builds with more files than the shell
#'   accepts as arguments. `1` compiles one file at a time. Defaults to 50.
#' @param customFolders Named character vector with the entries `"compiled"`,
#'   `"output"` and `"tmp"`, relative paths for compiled files, results and
#'   temporary files. `NULL` (default) uses the working directory.
#' @param input Character vector, the objects of the global environment to
#'   transfer. Defaults to all.
#' @param resetSeeds Logical. `TRUE` (default) removes `.Random.seed` from the
#'   transferred workspace and seeds each task with its SLURM job ID.
#' @param returnAll Logical. `TRUE` (default) makes `get()` fetch all files the
#'   job produced except the uploaded workspace and build files; `FALSE` only
#'   the result files.
#' @param libs Optional character vector of library paths put in front of the
#'   remote library search path, for the jobs and the remote build. The remote
#'   shell expands them, so `~` is the remote home.
#'
#' @return A list of functions:
#' \describe{
#'   \item{`check()`}{reports how many results are ready and returns `TRUE`
#'     when all are.}
#'   \item{`get()`}{copies the results into `<jobname>_folder/results/` and
#'     returns them, a list with one entry per task.}
#'   \item{`purge(purgeLocal)`}{deletes the job folder on the remote system,
#'     with `purgeLocal = TRUE` also the local one.}
#' }
#'
#' @seealso [runbg()], [profileParsPerNode()]
#'
#' @examples
#' \dontrun{
#' # Multi-start fits, 20 tasks
#' job <- distributedComputing(
#'   {
#'     mstrust(obj, center = pars, fits = 48, cores = 16, sd = 4)
#'   },
#'   jobname = "fits", machine = "user@@host", partition = "single",
#'   cores = 16, walltime = "02:00:00", nRep = 20
#' )
#' job$check()
#' fits <- job$get()
#' job$purge()
#'
#' # Profiles, a range of parameters per task
#' ranges <- profileParsPerNode(bestfit, 4)
#' job <- distributedComputing(
#'   {
#'     profile(obj, bestfit, whichPar = as.numeric(var_1):as.numeric(var_2),
#'             limits = c(-5, 5), cores = 16)
#'   },
#'   jobname = "profiles", machine = "user@@host", cores = 16,
#'   walltime = "02:00:00", varValues = ranges[c("from", "to")]
#' )
#' profiles <- do.call(rbind, job$get())
#' job$purge()
#' }
#'
#' @export
distributedComputing <- function(
    ...,
    jobname,
    partition = "single",
    cores = 16,
    nodes = 1,
    memPerCore = 2,
    walltime = "01:00:00",
    sshPasswd = NULL,
    machine = "cluster",
    varValues = NULL,
    nRep = NULL,
    recover = FALSE,
    purgeLocal = FALSE,
    compile = FALSE,
    link = FALSE,
    buildCores = NULL,
    buildBundle = 50,
    customFolders = NULL,
    resetSeeds = TRUE,
    returnAll = TRUE,
    input = ls(.GlobalEnv, all.names = TRUE),
    libs = NULL
){
  # The dots hold the code, unevaluated, and the deprecated argument names.
  dots <- match.call(expand.dots = FALSE)$...
  renames <- c(mem_per_core = "memPerCore", ssh_passwd = "sshPasswd",
               var_values = "varValues", no_rep = "nRep",
               purge_local = "purgeLocal", custom_folders = "customFolders")
  old <- intersect(names(dots), names(renames))
  callerEnv <- parent.frame()
  .renameArgs(lapply(dots[old], eval, envir = callerEnv), renames,
              "distributedComputing")
  code <- if (length(old)) dots[!names(dots) %in% old] else dots
  original_wd <- getwd()
  if (is.null(customFolders)) {
    output_folder_abs <- "./"
  } else if(!is.null(customFolders) & !all(length(customFolders) == 3 & sort(names(customFolders)) == c("compiled", "output", "tmp"))) {
    warning("'customFolders' must be named vector with exact three elements:\n
            'compiled', 'output', 'tmp', containing relative paths to the resp folders\n
            input is wrong, ignored.\n")
  } else {
    compiled_folder <- customFolders["compiled"]
    output_folder <- customFolders["output"]
    tmp_folder <- customFolders["tmp"]
    
    system(paste0("cp ", compiled_folder, "* ", tmp_folder))
    
    setwd(output_folder)
    output_folder_abs <- getwd()
    
    setwd(original_wd)
    setwd(tmp_folder)
  }
  
  # Safety rule: never compile and link at the same time
  if (compile) link <- FALSE
  
  
  on.exit(setwd(original_wd))
  
  # - definitions - #
  
  # relative path to the working directory
  
  wd_path <- paste0("./",jobname, "_folder/")
  
  # number of repetitions
  if(!is.null(nRep) & is.null(varValues)) {
    num_nodes <- nRep - 1
  } else if(is.null(nRep) & !is.null(varValues)) {
    num_nodes <- length(varValues[[1]]) - 1
  } else {
    stop("I dont know what you want how often done. Please set either 'nRep' or pass 'varValues' (_not_ both!)")
  }
  
  # define the ssh command depending on 'sshpass' being used
  ssh_command <- if (is.null(sshPasswd)) "ssh "
                 else paste0("sshpass -p ", sshPasswd, " ssh ")
  
  # - output functions - #
  # Structure of the output 
  out <- structure(vector("list", 3), names = c("check", "get", "purge"))
  
  # check function
  out[[1]] <- function() {

    ready <- length(suppressWarnings(system(
      paste0(ssh_command, machine, " 'ls ", jobname,
             "_folder/ 2>/dev/null | grep -E \"result[.]RData$\"'"),
      intern = TRUE)))

    if (ready >= num_nodes + 1) {
      cat("Result is ready!\n")
      return(TRUE)
    }
    cat("Result from", ready, "out of", num_nodes + 1, "nodes are ready.\n")
    FALSE
  }
  
  # get function
  out[[2]] <- function () {
    # copy all files back
    if (returnAll == T) {
      system(
        paste0(
          "mkdir -p ", output_folder_abs, "/", jobname,"_folder/results/; ",
          ssh_command, "-n ", machine, # go to remote
          # Everything the job produced, but not what was uploaded or built from
          # it. Patterns are double-quoted so the remote shell passes them to
          # tar instead of globbing them.
          " 'ZSTD_NBTHREADS=0 tar -C ", jobname, "_folder ",
          paste0("--exclude=\"", c("*_workspace.RData", "*_files.txt",
                                   "*.c", "*.cpp", "*.o", "*.a", "*.so"),
                 "\"", collapse = " "),
          " -I zstd -cf - ./'", # compress the remaining files on remote
          " | ", # pipe to local
          "",
          "tar -C ", output_folder_abs, "/", jobname,"_folder/results/ -I zstd -xf -"
        )
      )
    } else {
      # copy only result files back
      system(
        paste0(
          "mkdir -p ", output_folder_abs, "/", jobname,"_folder/results/; ",
          ssh_command, "-n ", machine, " '",
          "ZSTD_NBTHREADS=0 find ", jobname, "_folder -type f -name \"*result.RData\" -exec tar -I zstd -cf - {} +'",
          " | tar --strip-components=1 -I zstd -x -C ", shQuote(paste0(output_folder_abs, "/", jobname,"_folder/results/"))
        )
      )
    }
    
    
    # get list of all currently available output files
    result_list <- vector("list", num_nodes + 1)
    result_dir <- paste0(output_folder_abs, "/", jobname, "_folder/results/")
    result_files <- list.files(path = result_dir, pattern = glob2rx("*result.RData"))

    ## Nothing here is the normal state after an OOM kill or a walltime overrun.
    if (!length(result_files))
      cat("\n\tNo results retrieved. The jobs are still running or ended before",
          "\n\twriting a result. Check the .err files in ", jobname,
          "_folder/ on ", machine, ".\n", sep = "")
    else if (length(result_files) != num_nodes + 1)
      cat("\n\t", length(result_files), " of ", num_nodes + 1, " results ready\n",
          sep = "")

    ## Slot by SLURM array index: with partial results the file order differs.
    node_ID <- suppressWarnings(as.integer(
      sub(".*_([0-9]+)_result\\.RData$", "\\1", result_files)))
    slot <- ifelse(is.na(node_ID) | node_ID < 0L | node_ID > num_nodes,
                   seq_along(result_files), node_ID + 1L)

    for (i in seq_along(result_files)) {
      cluster_result <- NULL
      check <- try(load(file = paste0(result_dir, result_files[i])), silent = TRUE)
      if (inherits(check, "try-error"))
        warning("Could not load '", result_files[i], "': ",
                conditionMessage(attr(check, "condition")), call. = FALSE)
      else
        result_list[[slot[i]]] <- cluster_result
    }
    result_list
  }
  
  
  
  # purge function, its default taken from the argument of distributedComputing()
  purgeLocalDefault <- purgeLocal
  out[[3]] <- function (purgeLocal = purgeLocalDefault, ...) {
    .renameArgs(list(...), c(purge_local = "purgeLocal"), "purge",
                strict = TRUE)
    # remove files remote
    system(
      paste0(
        ssh_command, machine, " rm -rf ", jobname, "_folder"
      )
    )
    # also remove local files if want so
    if (purgeLocal) {
      system(
        paste0("rm -rf ", output_folder_abs, "/", jobname,"_folder")
      )
    }
  }
  
  # if recover == T, stop here
  if(recover) {
    setwd(original_wd)
    return(out)
  } 
  
  
  
  
  
  
  
  
  
  # - calculation - #
  # create wd for this run
  system(
    paste0("rm -rf ",jobname,"_folder")
  )
  system(
    paste0("mkdir ",jobname,"_folder/")
  )
  
  
  
  # Export the workspace the job needs. The default covers hidden names too,
  # because resetSeeds acts on a transferred .Random.seed.
  input <- intersect(input, ls(.GlobalEnv, all.names = TRUE))
  ws_file <- .saveWorkspace(input, paste0(wd_path, jobname, "_workspace.RData"))
  ## The job script loads the plain .RData. Not every zstd has --rm, and a failed
  ## expansion must not reach the sbatch.
  ws_unpack <- if (grepl("\\.zst$", ws_file))
    paste0("{ zstd -dqf ", basename(ws_file), " && rm -f ", basename(ws_file),
           "; } || exit 1; ") else ""
  
  
  # WRITE R
  
  
  # The node script attaches dMod2 only, not the submitting session's
  # packages; anything else belongs in the expression.
  package_list <- "library(dMod2)"
  if (compile || link) {
    objfns <- 'obj.fns <- ls()[sapply(ls(), function(nm) inherits(get(nm, envir=.GlobalEnv), c("obsfn", "parfn", "prdfn", "objfn")))]'
    setmn <- sprintf('for (o in obj.fns) eval(parse(text=paste0("modelname(", o, ") <- \'%s\'")))\n', paste0(jobname, "_shared_object"))
    load_so <- paste0("dyn.load('",jobname,"_shared_object.so')")
  } else {
    objfns <- ""
    setmn <-""
    load_so <- ""
  }
  
  
  
  # generate parameter lists
  if (!is.null(varValues)) {
    var_list <- paste(
      lapply(
        seq(1,length(varValues)),
        function(i) {
          if (is.character(varValues[[i]])) {
            paste0("var_values_", i, "=c('", paste(varValues[[i]], collapse="','"),"')")
          } else {
            paste0("var_values_", i, "=c(", paste(varValues[[i]], collapse=","),")")
          }
          
          
        }
      ),
      collapse = "\n"
    )
    # List of all names of parameters that will be changes between runs
    var_names <- paste(lapply(seq(1,length(varValues)), function(i) paste0("var_",i)))
    
    # Variables per run
    var_per_run <- paste(
      lapply(
        seq(1, length(varValues)),
        function(i) {
          paste0("var_", i, "=var_values_",i,"[(as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) + 1)]")
        }
      ),
      collapse = "\n"
    )
  } else {
    var_list <- ""
    var_per_run <- ""
  }
  
  
  # define fixed pars
  fixedpars <- paste(
    "node_ID = Sys.getenv('SLURM_ARRAY_TASK_ID')",
    "job_ID = Sys.getenv('SLURM_JOB_ID')",
    paste0("jobname = ", "'",jobname,"'"),
    sep = "\n"
  )
  
  
  # WRITE R
  if (!length(code))
    stop("distributedComputing: no R code to run.", call. = FALSE)
  expr <- as.expression(code[[1L]])
  cat(
    paste(
      "#!/usr/bin/env Rscript",
      "",
      "# Load packages",
      package_list,
      "",
      "# Load environment",
      paste0("load('",jobname,"_workspace.RData')"),
      "",
      if (resetSeeds == TRUE & exists(".Random.seed")) {
        paste0("# remove random seeds\nrm(.Random.seed)\nset.seed(as.numeric(Sys.getenv('SLURM_JOB_ID')))")
        
      },
      "",
      "# load shared object if precompiled",
      objfns,
      setmn,
      load_so,
      "",
      "# List of variablevalues",
      var_list,
      "",
      "# Define variable values per run",
      var_per_run,
      "",
      "# Fixed parameters",
      fixedpars,
      "",
      "",
      "",
      "# Paste function call",
      paste0("cluster_result <- try(", as.character(expr),")"),
      sep = "\n",
      "save(cluster_result, file = paste0(jobname,'_', node_ID, '_result.RData'))" #'_',job_ID, 
    ),
    file = paste0(wd_path, jobname,".R")
  )
  
  
  
  
  # WRITE BASH
  cat(
    paste(
      "#!/bin/bash",
      "",
      "# Job name",
      paste0("#SBATCH --job-name=",jobname),
      "# Define format of output, deactivated",
      paste0("#SBATCH --output=",jobname,"_%j-%a.out"),
      "# Define format of errorfile, deactivated",
      paste0("#SBATCH --error=",jobname,"_%j-%a.err"),
      "# Define partition",
      paste0("#SBATCH --partition=", partition),
      "# Define number of nodes per task",
      paste0("#SBATCH --nodes=", nodes),
      "# Define number of cores per node",
      paste0("#SBATCH --ntasks-per-node=",cores),
      "# Define walltime",
      paste0("#SBATCH --time=",walltime),
      "# Define of repetition",
      paste0("#SBATCH -a 0-", num_nodes),
      "# memory per CPU core",
      paste0("#SBATCH --mem-per-cpu=", memPerCore, "gb"),
      "",
      "",
      "# Load compiler modules",
      "module load compiler/gnu/13.3",
      "# Load R modules",
      "module load math/R",
      paste0("export OMP_NUM_THREADS=","1"),
      paste0("export MKL_NUM_THREADS=", "1"),
      .remoteLibs(libs),
      "",
      "# Run R script",
      paste0("Rscript ", jobname, ".R"),
      sep = "\n" 
    ),
    file = paste0(wd_path,jobname,".sh")
  )
  
  
  # The build script ships in the job folder, as a bare `R CMD SHLIB` misses
  # cppDE's includes and BLAS/LAPACK; `&&` makes a failed build skip the sbatch.
  build_script_file <- paste0(jobname, "_build.sh")
  filelist_file <- paste0(jobname, "_files.txt")
  module_cmd <- paste0("module load compiler/gnu/13.3 2>/dev/null; module load math/R; ",
                       .remoteLibs(libs))

  ## File names go in a list in the job folder, not on the command line that
  ## both tar calls share.
  transferCmds <- function(files) {
    writeLines(files, paste0(wd_path, filelist_file))
    list(locale = paste0("tar -I 'zstd -T0' -cf - -T ", wd_path, filelist_file,
                         " ", wd_path, "*"),
         remote = paste0("tar -C ./ -I zstd -xf - ; xargs -r -n 100 mv -t ./",
                         jobname, "_folder < ./", jobname, "_folder/",
                         filelist_file, "; "))
  }

  if (compile) {
    # --- FULL RECOMPILATION (.cpp/.c -> .so) ---
    sourcefiles <- c(list.files(pattern = glob2rx("*.c")),
                     list.files(pattern = glob2rx("*.cpp")))

    tarCmds    <- transferCmds(sourcefiles)
    tar_locale <- tarCmds$locale
    tar_remote <- tarCmds$remote

    buildinfo <- .remoteBuildInfo()
    cat(.remoteBuildScript(
      files       = sourcefiles,
      output      = paste0(jobname, "_shared_object.so"),
      compileArgs = buildinfo$compileArgs,
      needsCVODE  = buildinfo$needsCVODE,
      needsKLU    = buildinfo$needsKLU,
      needsLapack = buildinfo$needsLapack,
      link        = FALSE,
      cxx         = any(grepl("\\.cpp$", sourcefiles)),
      cores       = buildCores,
      filelist    = filelist_file,
      bundle      = buildBundle
    ), file = paste0(wd_path, build_script_file))

    compile_remote <- paste0(module_cmd, "bash ", build_script_file, " && ")

  } else if (link) {
    # --- LINK ONLY (.o -> .so) ---
    object_files <- Sys.glob("*.o")
    if (length(object_files) == 0)
      stop("No .o files found for linking! You must compile first.")

    message("distributedComputing(): linking pre-compiled .o files on the cluster. ",
            "This requires the cluster toolchain to be ABI-compatible with the ",
            "local one; if the link or dyn.load() fails, resubmit with compile = TRUE.")

    tarCmds    <- transferCmds(object_files)
    tar_locale <- tarCmds$locale
    tar_remote <- tarCmds$remote

    buildinfo <- .remoteBuildInfo()
    cat(.remoteBuildScript(
      files       = object_files,
      output      = paste0(jobname, "_shared_object.so"),
      compileArgs = buildinfo$compileArgs,
      needsCVODE  = buildinfo$needsCVODE,
      needsKLU    = buildinfo$needsKLU,
      needsLapack = buildinfo$needsLapack,
      link        = TRUE,
      cxx         = length(Sys.glob("*.cpp")) > 0,
      cores       = buildCores,
      filelist    = filelist_file,
      bundle      = buildBundle
    ), file = paste0(wd_path, build_script_file))

    compile_remote <- paste0(module_cmd, "bash ", build_script_file, " && ")

  } else {
    # --- NO BUILD ACTION (shared objects already available) ---
    tarCmds    <- transferCmds(Sys.glob("*.so"))
    tar_locale <- tarCmds$locale
    tar_remote <- tarCmds$remote
    compile_remote <- ""
  }
  
  # transfer and run files
  status <- system(
    paste0(
      tar_locale, # compress all files in the local working dir
      " | ", ssh_command, machine, # pipe to ssh session on remote
      " 'if [ -d ", jobname, "_folder ]; then rm -Rf ", jobname,"_folder; fi ;", # remove folder if it exists
      " mkdir -p ", jobname,"_folder; ", # create new wd on remote
      tar_remote, # uncompress files in wd on remote, if necessary move files
      "cd ", jobname, "_folder; ", # change in said wd
      ws_unpack, # expand the workspace if it was shipped compressed
      compile_remote, # compile files if said so, if not nothing happen
      "sbatch ", jobname, ".sh'" # start bash script
    )
  )

  # ssh forwards the exit status of the remote command chain, so a failed
  # remote build shows up here instead of silently queueing a broken job.
  if (status != 0)
    warning("Remote build or job submission failed (exit status ", status,
            "). See the compiler output above; nothing was submitted if the ",
            "build is what failed.", call. = FALSE)

  setwd(original_wd)
  return(out)
}



#' Split Parameters into Ranges for Distributed Profiles
#'
#' Splits the parameter indices into consecutive ranges, one per node, for
#' profile calculation with [distributedComputing()].
#'
#' @param parameters Named vector or list of parameters; only its length is
#'   used.
#' @param parsPerNode Integer, number of parameters per node.
#' @param side `"both"` (default) for one entry per range, or `"split"` for two
#'   entries per range, labelled `"left"` and `"right"`, to compute the two
#'   sides of each profile on separate nodes.
#' @param ... `fits_per_node` is deprecated, use `parsPerNode`.
#'
#' @return A list with the integer vectors `from` and `to`, the first and last
#'   parameter index of each node, and the character vector `side`.
#' @seealso [distributedComputing()], [profile()][profile.objfn]
#' @examples
#' parameters <- setNames(1:10, letters[1:10])
#' profileParsPerNode(parameters, 4)
#' profileParsPerNode(parameters, 4, side = "split")
#'
#' @export
profileParsPerNode <- function(parameters, parsPerNode, side = c("both", "split"),
                               ...) {

  .renameArgs(list(...), c(fits_per_node = "parsPerNode"), "profileParsPerNode",
              strict = TRUE)
  side <- match.arg(side)

  # get the number of parameters
  n_pars <- length(parameters)
  
  # determine the number of nodes necessary
  no_nodes <- 1:ceiling(n_pars/parsPerNode)
  
  # generate the lists which parameters are send to which node
  pars_from <- parsPerNode
  pars_to_vec <- parsPerNode
  while (pars_from < (n_pars)) {
    pars_from <- pars_from + parsPerNode
    pars_to_vec <- c(pars_to_vec, pars_from)
  }
  pars_to_vec[length(pars_to_vec)] <- n_pars
  
  pars_from_vec <- c(1, pars_to_vec+1)
  pars_from_vec <- head(pars_from_vec, -1)
  
  # adjust for sides 
  if (side == "both") {
    side_vec <-  rep("both", length(no_nodes))
  } else {
    # split pars_to_vec and pars_from_vec by repeating each element twice
    pars_to_vec <- rep(pars_to_vec, each = 2)
    pars_from_vec <- rep(pars_from_vec, each = 2)
    side_vec <- rep(c("left", "right"), length(no_nodes))
  }
  
  out <- list(from=pars_from_vec, to=pars_to_vec, side = side_vec)
  
  return(out)
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

  avail <- suppressWarnings(parallel::detectCores())
  if (!is.na(avail) && o * k > avail)
    warning(sprintf("cores: %d (%s) x %d (conditions) = %d exceeds the %d available.",
                    o, outer, k, o * k, avail), call. = FALSE)
  list(outer = o, conditions = k)
}


.sanitizeCores <- function(cores)  {
  
  max.cores <- parallel::detectCores()
  min(max.cores, cores)
  
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
  if (!is.null(coresConditions) && cores == 1L)
    options(dMod.cores = as.integer(coresConditions))
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
                          .packages = loaded_packages) %dopar% {
    if (!is.null(kcond)) options(dMod.cores = kcond)
    FUN(X[[i]])
  }
  out
}



`%dopar%` <- foreach::`%dopar%`
