#' Detect number of free cores
#'
#' @description Estimates free cores from the 1-min load average.
#' Supports Linux, macOS, and remote machines via SSH.
#' On Windows, returns 1 with a warning (no load average available).
#' Result is floored at 1 so it can be fed straight into
#' `mclapply(mc.cores = ...)` without crashing under heavy load.
#' @param machine character vector of SSH hosts, e.g. "user@@localhost".
#' NULL (default) for the local machine.
#' @return numeric vector of free cores (>= 1) with attributes "ncores" and "used".
#' @export
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
    
    # Floor at 1: callers feed `free` straight into mclapply(mc.cores = ...)
    # which rejects 0. Under heavy load the 1-min average can exceed nCores;
    # reporting 1 means "serialise, don't die" instead of crashing.
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

  # CRAN policy: R CMD check sets _R_CHECK_LIMIT_CORES_ and parallel forbids
  # mc.cores > 2 in that mode. Cap before returning so all mclapply callsites
  # (P, normL2, compile) stay legal under check without each having to
  # know about the env var.
  chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  if (nzchar(chk) && chk != "false") {
    if (length(freeCores) > 0L) freeCores[] <- pmin(freeCores, 2L)
  }

  freeCores
}


## Remote build helpers --------------------------------------------------------
##
## Rebuilding a dMod model on a remote machine takes more than a bare
## `R CMD SHLIB`: the cppDE backend emits sources that `#include <cppde/...>`
## from the cppDE package tree, and the shared object gets linked against
## BLAS/LAPACK (plus Sundials/KLU for CVODE and sparse models). `compile()`
## assembles those flags from the `"compileInfo"` attribute that every model
## object carries, but every path in there points into the *local* library
## tree and is meaningless on the cluster.
##
## `.remoteBuildInfo()` therefore keeps only the portable part of that
## information (the `-D` macro flags plus which optional backends are in play),
## and `.remoteBuildScript()` writes a build script that resolves every
## path-valued flag on the remote machine itself, by asking the R installation
## that is on PATH there.

## Scan a workspace for model objects and return the portable build settings
## needed to rebuild their C/C++ sources elsewhere.
.remoteBuildInfo <- function(envir = .GlobalEnv) {

  compileArgs <- character(0)
  needsCVODE <- FALSE
  needsKLU <- FALSE

  for (nm in ls(envir)) {
    o <- try(get(nm, envir = envir), silent = TRUE)
    if (inherits(o, "try-error")) next
    if (isTRUE(attr(o, "sparse"))) needsKLU <- TRUE
    info <- attr(o, "compileInfo")
    if (is.null(info)) next
    for (e in info) {
      ca <- trimws(e$compileArgs %||% "")
      ## Non-empty linkArgs mean the backend pulls external libraries
      ## (Sundials, and KLU on top of it for sparse CVODE models); those have
      ## to be resolved from the remote cppDE installation.
      if (nzchar(trimws(e$linkArgs %||% ""))) needsCVODE <- TRUE
      if (isTRUE(e$sparse)) needsKLU <- TRUE
      if (nzchar(ca)) compileArgs <- c(compileArgs, strsplit(ca, "\\s+")[[1]])
    }
  }

  compileArgs <- unique(compileArgs[nzchar(compileArgs)])
  ## cppDE tags sparse models with -DKLUBTF/-DKLUAMD; compile() adds the
  ## -DKLU switch itself, so mirror that here.
  if (any(grepl("^-DKLU", compileArgs))) needsKLU <- TRUE

  list(compileArgs = paste(compileArgs, collapse = " "),
       needsCVODE  = needsCVODE,
       needsKLU    = needsKLU)
}

## TRUE when naming every file on one command line would overrun the shell's
## argument limit; same test as compile() applies locally.
.remoteNeedsChunking <- function(files)
  sum(nchar(files) + 1L) > .compileCmdLimit()

## Write the job workspace through zstd. Uncompressed it hits the disk in full
## only for tar to read it back and compress it, and with one parfn per condition
## that is several GB. zstd is already required for the transfer; the fallback
## covers a submitting machine without the binary. Returns the file written.
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

## Generate the bash script that builds the shared object on the remote machine.
## `files` are the sources, or the object files when `link = TRUE`. Beyond the
## argument limit the script switches to the chunked archive build compile()
## uses locally, reading its inputs from `filelist`.
.remoteBuildScript <- function(files, output, compileArgs = "",
                               needsCVODE = FALSE, needsKLU = FALSE,
                               link = FALSE, cxx = FALSE, cores = 1,
                               workdir = NULL, filelist = NULL, chunkSize = 100,
                               bundle = 50) {

  files <- unlist(strsplit(trimws(files), "\\s+"))
  files <- files[nzchar(files)]

  cflags <- trimws(paste("-O2 -DNDEBUG -w -fPIC", compileArgs,
                         if (needsKLU) "-DKLU"))

  ldflags <- character(0)
  if (link) {
    ## Object files carry the *local* toolchain's LTO bytecode, which a remote
    ## compiler of a different GCC generation refuses to read ("bytecode stream
    ## in file 'e.o' generated with LTO version 16.0 instead of the expected
    ## 13.1"). Fat LTO objects still contain ordinary machine code, so turning
    ## LTO and the linker plugin off lets the link fall back to it.
    ldflags <- c(ldflags, "-fno-lto", "-fno-use-linker-plugin")
    ## With only .o inputs R CMD SHLIB links via the C driver, which would not
    ## pull in the C++ runtime the cppDE-generated code needs.
    if (cxx) ldflags <- c(ldflags, "-lstdc++")
  }

  ## SUNDIALS/KLU flags come from the *remote* cppDE install, so they are
  ## resolved by Rscript inside the generated script. Embedded in a
  ## single-quoted shell string, must contain no single quotes.
  cfgExpr <- function(field) paste0(
    "cfg <- get0(\"cvodeConfig\", envir = asNamespace(\"cppDE\"), inherits = FALSE); ",
    "cat(if (is.environment(cfg)) paste(unlist(mget(c(", field, "), envir = cfg, ",
    "ifnotfound = \"\")), collapse = \" \") else \"\")")

  extra_libs <- NULL
  if (needsKLU || needsCVODE) {
    fields <- paste0(c(if (needsKLU) "\"klu_libs\"", if (needsCVODE) "\"libs\""),
                     collapse = ", ")
    extra_libs <- paste0("PKG_LIBS=\"$PKG_LIBS $(Rscript -e '", cfgExpr(fields), "')\"")
  }

  ## KLU needs its include path at compile time, not just -DKLU.
  klu_cppflags <- if (needsKLU)
    paste0("PKG_CPPFLAGS=\"$PKG_CPPFLAGS $(Rscript -e '",
           cfgExpr("\"klu_cflags\""), "')\"")

  ## Job count: a fixed `cores`, or whatever the remote machine reports.
  nproc <- if (is.null(cores)) "$NPROC" else as.character(max(1L, as.integer(cores)))
  detect <- if (is.null(cores))
    c("NPROC=$(nproc 2>/dev/null || echo 1)",
      "if [ \"$NPROC\" -gt 16 ]; then NPROC=16; fi", "")

  ## Precompiled header, decided here because the prologue check needs the
  ## sources. A missing .gch is harmless, the header then just includes what
  ## the sources include anyway.
  cxxSrc <- files[grepl("\\.cpp$", files, ignore.case = TRUE)]
  pchInc <- if (!link && length(cxxSrc) >= 8L) .compilePCHIncludes(cxxSrc)
  toolchain <- if (!is.null(filelist) || !is.null(pchInc)) c(
    "CC=$(R CMD config CC);   CFLAGS=$(R CMD config CFLAGS)",
    "CXX=$(R CMD config CXX); CXXFLAGS=$(R CMD config CXXFLAGS)",
    "CPICFLAGS=$(R CMD config CPICFLAGS); CXXPICFLAGS=$(R CMD config CXXPICFLAGS)",
    "RINC=\"-I$(R RHOME)/include\"",
    "AR=$(R CMD config AR); RANLIB=$(R CMD config RANLIB)",
    "export CC CXX CFLAGS CXXFLAGS CPICFLAGS CXXPICFLAGS RINC", "")
  pchBlock <- if (!is.null(pchInc)) c(
    "cat > dMod_pch.hpp <<'DMOD_PCH_EOF'", pchInc, "DMOD_PCH_EOF",
    paste("$CXX $RINC $PKG_CPPFLAGS $PKG_CXXFLAGS $CXXPICFLAGS $CXXFLAGS",
          "-x c++-header dMod_pch.hpp -o dMod_pch.hpp.gch || true"),
    "PKG_CXXFLAGS=\"$PKG_CXXFLAGS -include dMod_pch.hpp\"", "")

  ## Compile off the file list in parallel, then link. Asking make to
  ## parallelise R CMD SHLIB via MAKEFLAGS does not work, so the objects are
  ## built here and SHLIB only links them.
  chunked <- !is.null(filelist) && .remoteNeedsChunking(files)
  isCxx   <- grepl("\\.cpp$", files, ignore.case = TRUE)

  ## Only the chunked path bundles: below the argument limit SHLIB gets the
  ## sources themselves, so bundles would just be compiled twice.
  doBundle <- chunked && !link && as.integer(bundle) > 1L && sum(isCxx) > 1L

  anchorLine <- if (chunked)
    paste0("ANCHOR=", shQuote(files[if (any(isCxx)) which(isCxx)[1] else 1L]))

  ## The sources hold one small function each and all pull in the same headers,
  ## so bundling them into few translation units parses those headers once per
  ## bundle. Every generated function survives; only the compiler runs less.
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
    "# Generated by dMod -- rebuilds the model shared object on this machine.",
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
    "PKG_CPPFLAGS=\"-I$CPPDE_INC\"",
    klu_cppflags,
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


#' Run an R expression in the background (only on UNIX)
#' 
#' @description Generate an R code of the expression that is copied via `scp`
#' to any machine (ssh-key needed). Then collect the results.
#' @details `runbg()` generates a workspace from the `input` argument
#' and copies the workspace to the remote machines via `scp`. This will only
#' work if *an ssh-key had been generated and added to the authorized keys
#' on the remote machine*. The code snippet, i.e. the `...` argument, can
#' include several intermediate results but only the last call which is not
#' redirected into a variable is returned via the variable `.runbgOutput`,
#' see example below.
#'
#' Depending on the `compile` and `link` arguments, build-related files are
#' handled as follows:
#' \itemize{
#'   \item `compile = TRUE`: C/C++ source files are transferred and compiled
#'         remotely via `R CMD SHLIB`.
#'   \item `link = TRUE`: Pre-compiled object files (`.o`) are transferred and
#'         linked remotely.
#'   \item Both `FALSE` (default): Existing shared objects (`.so`) are copied
#'         directly.
#' }
#' When `compile` or `link` is `TRUE`, objects of class `obsfn`, `parfn`, or
#' `prdfn` in the workspace automatically get their `modelname` updated to
#' point to the newly built shared object.
#'
#' The remote build runs from a generated shell script that resolves all
#' path-valued compiler flags (cppDE include directory, BLAS/LAPACK, and
#' Sundials/KLU for CVODE and sparse models) *on the remote machine*, since
#' the local library paths do not carry over. The model-specific `-D` macros
#' are read from the `"compileInfo"` attribute of the model objects in the
#' workspace. `cppDE` therefore has to be installed for the R that is on
#' `PATH` on the remote machine. Only the files the chosen mode consumes are
#' transferred, and beyond the shell's argument limit the script switches to
#' the chunked static-archive build of [compile()].
#' @param ... Some R code.
#' @param machine Character vector, e.g. `"localhost"` or `"knecht1.fdm.uni-freiburg.de"`
#' or `c("localhost", "localhost")`.
#' @param filename Character, defining the filename of the temporary file. A random
#' file name is chosen if `NULL`.
#' @param input Character vector, the objects in the workspace that are stored
#' into an R data file and copied to the remote machine.
#' @param compile Logical. If `TRUE`, C/C++ source files (`.c`, `.cpp`) are
#' transferred to the remote machine and fully recompiled into a shared object
#' (`.so`). If set to `TRUE`, this overrides `link = TRUE`.
#' @param link Logical. If `TRUE`, only existing object files (`.o`) are
#' transferred to the remote machine and linked into a shared object (`.so`),
#' skipping compilation. If no `.o` files are found, an error is raised.
#' This option is ignored if `compile = TRUE`. Object files are toolchain
#' specific, so this only works when the remote compiler is ABI-compatible
#' with the local one; prefer `compile = TRUE` when the two machines run
#' different compiler generations.
#' @param wait Logical. Wait until executed. If `TRUE`, the code checks if the
#' result file is already present in which case it is loaded. If not present,
#' `runbg()` starts, produces the result and loads it as `.runbgOutput` directly
#' into the workspace. If `wait = FALSE`, `runbg()` starts in the background
#' and the result is only loaded into the workspace when the `get()` function
#' is called, see Value section.
#' @param recover Logical. This option is useful to recover the functions
#' `check()`, `get()`, `purge()` and `terminate()`, e.g. when a session has
#' crashed. Then, the functions are recreated without restarting the job. They
#' can then be used to get the results of a job without having to do it manually.
#' Requires the correct filename, so if the previous `runbg()` was run with
#' `filename = NULL`, you have to specify the filename manually.
#' @param walltime Optional character. Maximum runtime in the format `"HH:MM:SS"`.
#' If exceeded, the job will be terminated.
#' @return List of functions `check()`, `get()`, `purge()` and `terminate()`. 
#' `check()` checks if the result is ready.
#' `get()` copies the result file
#' to the working directory and loads it into the workspace as an object called `.runbgOutput`. 
#' This object is a list named according to the machines that contains the results returned by each
#' machine.
#' `purge()` deletes the temporary folder
#' from the working directory and the remote machines.
#' `terminate()` kills all running processes associated with this job on the remote machines.
#' @export
#' @examples
#' \dontrun{
#' out_job1 <- runbg({
#'          M <- matrix(rnorm(1e2), 10, 10)
#'          solve(M)
#'          }, machine = c("localhost", "localhost"), filename = "job1")
#' out_job1$check()          
#' out_job1$get()
#' result <- .runbgOutput
#' print(result)
#' out_job1$purge()
#' }
#' \dontrun{
#' # Recover a runbg job with the option "recover"
#' out_job1 <- runbg({
#'          M <- matrix(rnorm(1e2), 10, 10)
#'          solve(M)
#'          }, machine = c("localhost", "localhost"), filename = "job1")
#' Sys.sleep(1)
#' remove(out_job1)
#' try(out_job1$check())
#' out_job1 <- runbg({
#'   "This code is not run"
#' }, machine = c("localhost", "localhost"), filename = "job1", recover = TRUE)
#' out_job1$get()
#' result <- .runbgOutput
#' print(result)
#' out_job1$purge()
#' }
runbg <- function(..., machine = "localhost", filename = NULL, input = ls(.GlobalEnv), compile = FALSE, link = FALSE, wait = FALSE, recover = FALSE, walltime = NULL) {
  
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
    
    .GlobalEnv$.runbgOutput <- result
    
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
      load(file = resultfile[m])
      result[[m]] <- .runbgOutput
    }
    .GlobalEnv$.runbgOutput <- result
    return(out)
  }
  
  # Save current workspace to be transferred to remote machines
  save(list = input, file = paste0(filename0, ".RData"), envir = .GlobalEnv)
  
  # Collect currently loaded packages to replicate the library state remotely
  pack <- sapply(strsplit(search(), "package:", fixed = TRUE), function(v) v[2])
  pack <- pack[!is.na(pack)]
  pack <- paste(paste0("try(library(", pack, "))"), collapse = "\n")
  
  output <- ".runbgOutput"
  
  # Compiler flags mirroring compile() in compile.R. The flags themselves are
  # resolved on the remote machine (see .remoteBuildScript); here only
  # collect the portable, model-specific part and the file list. Everything is
  # written into a shell script to avoid quoting issues with nested SSH commands.
  buildinfo <- list(compileArgs = "", needsCVODE = FALSE, needsKLU = FALSE)
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
    objfns <- 'obj.fns <- ls()[sapply(ls(), function(nm) inherits(get(nm, envir=.GlobalEnv), c("obsfn", "parfn", "prdfn")))]'
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
        link        = link,
        cxx         = has_cxx,
        cores       = NULL,
        workdir     = paste0(filename[m], "_folder"),
        filelist    = filelist_file
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
      " 'export OMP_NUM_THREADS=1 && export MKL_NUM_THREADS=1",
      if (nzchar(compile_cmd)) paste0(" && ", compile_cmd),
      " && R CMD BATCH --vanilla ", filename[m], "_folder/", filename[m], ".R'"
    ), intern = FALSE, wait = wait)
  }
  
  if (wait) {
    out$get()
    out$purge()
  } else {
    return(out)
  }
  
}



## ---- HPC/SLURM distributed computing (moved from toolsSeverin.R) ----------

#' Run any R function on a remote HPC system with SLURM
#'
#' @description
#' Generates R and bash scripts, transfers them to a remote HPC system via SSH,
#' and executes the given R code in parallel using the SLURM batch manager.
#' The function handles workspace export, job submission, and result retrieval.
#'
#' @details
#' `distributedComputing()` generates R and bash scripts designed to run
#' on an HPC system managed by SLURM. The current R workspace together with the
#' scripts are exported and transferred to the remote system via SSH.  
#' If ssh-key authentication is not possible, the SSH password can be provided and
#' is used by `sshpass` (which must be installed locally).
#'
#' The code to be executed remotely is passed to the `...` argument; its final
#' output is stored in `cluster_result`, which can be loaded in the local
#' workspace via the `get()` function.
#'
#' It is possible to run multiple repetitions of the same program (via `no_rep`)
#' or to pass a list of parameter arrays through `var_values`. Parameters that
#' vary between runs must be named `var_i`, where *i* matches the index
#' of the corresponding array in `var_values`.
#'
#' The workspace is serialised through `zstd` and expanded again on the cluster
#' before the build starts, so it is written and transferred once, compressed.
#' With one parameter transformation per condition it is the bulk of the upload.
#' Without `zstd` locally it is written uncompressed; the archive is compressed
#' either way.
#'
#' @param ... R code to be remotely executed. Parameters to be changed for each run
#' must be named `var_i` (see "Details").
#' @param jobname Unique name (character) for the run. Existing runs with the same
#' name will be overwritten. Must not contain the string "Minus".
#' @param partition SLURM partition name to use. Default is `"single"`.
#' @param cores Number of cores per node. Values above 16 may limit available nodes.
#' @param nodes Number of nodes per task. Default is 1; typically should not be changed.
#' @param mem_per_core Memory per CPU core in GB. Default is 2 GB.
#' @param walltime Maximum runtime in format `"hh:mm:ss"`. Default is 1 hour.
#' @param ssh_passwd Password string for SSH authentication via `sshpass`.
#' Optional, and only used if key-based authentication is unavailable.
#' @param machine SSH address of the remote HPC system, e.g. `"user@@cluster"`.
#' @param var_values List of parameter arrays. Each array corresponds to one variable
#' `var_i`. The length of each array determines the number of SLURM array jobs.
#' Mutually exclusive with `no_rep`.
#' @param no_rep Integer number of repetitions (mutually exclusive with `var_values`).
#' @param recover Logical; if `TRUE`, no computation is performed. Instead,
#' the returned list of functions `check()`, `get()`, and `purge()`
#' can be used to interact with previously submitted jobs.
#' @param purge_local Logical; if `TRUE`, the `purge()` function also
#' deletes local result files.
#' @param compile Logical; if `TRUE`, all C/C++ source files (`*.c`, `*.cpp`)
#' are transferred to the cluster and fully recompiled into shared objects (`.so`).
#' If set to `TRUE`, this overrides `link = TRUE`. The build runs from a
#' generated shell script that resolves all path-valued compiler flags
#' (cppDE include directory, BLAS/LAPACK, and Sundials/KLU for CVODE and
#' sparse models) on the cluster, so `cppDE` has to be installed for the R
#' that `module load math/R` provides there. The job is only submitted if the
#' build succeeds. Beyond the shell's argument limit the script switches to the
#' chunked static-archive build of [compile()].
#' @param link Logical; if `TRUE`, only existing object files (`*.o`) are
#' transferred to the cluster and linked into shared objects (`.so`),
#' skipping compilation. If no `.o` files are found, an error is raised.
#' This option is ignored if `compile = TRUE`. Object files are toolchain
#' specific, so this only works when the cluster compiler is ABI-compatible
#' with the local one, in particular, `.o` files produced with link-time
#' optimisation by a newer GCC cannot be read by an older one. Prefer
#' `compile = TRUE` when the two machines run different compiler generations.
#' @param buildCores Number of compiler processes the remote build runs in
#' parallel. The build happens on the login node, before the job is queued, so
#' this is unrelated to `cores`, which sizes the SLURM allocation. `NULL`
#' (default) lets the build use what the login node reports, capped at 16.
#' @param buildBundle Number of generated sources the remote build puts into one
#' translation unit. They hold one small function each and all include the same
#' headers, so bundling parses those once per bundle instead of once per file;
#' every generated function is kept either way. Only used above the shell's
#' argument limit, where the build goes through a static archive. `1` compiles
#' one file at a time.
#' @param custom_folders Named vector with exactly three elements: `"compiled"`,
#' `"output"`, and `"tmp"`. Each value is a relative path specifying where
#' compiled files, temporary data, and output results should be stored.
#' If `NULL`, all operations occur in the current working directory.
#' @param input Character vector of object names in the global environment to
#' transfer. Defaults to the whole workspace; naming only what the expression
#' uses keeps the transferred `.RData` small.
#' @param resetSeeds Logical; if `TRUE` (default), removes `.Random.seed`
#' from the transferred workspace to ensure each node has independent random seeds.
#' @param returnAll Logical; if `TRUE` (default), retrieves everything the job
#' produced, excluding the uploaded workspace and the build artefacts, which
#' are already local. If `FALSE`, only result files (`*result.RData`) are
#' fetched.
#'
#' @return
#' A list containing three functions:
#' \itemize{
#'   \item `check()` - Checks whether all remote results are complete.
#'   \item `get()` - Downloads results and loads them into
#'         `cluster_result` in the local workspace.
#'   \item `purge()` - Deletes temporary remote files; optionally removes local ones.
#' }
#'
#' @examples
#' \dontrun{
#' outDistributedComputing <- distributedComputing(
#' {
#'   mstrust(
#'     objfun=objective_function,
#'     center=outer_pars,
#'     name = "study",
#'     rinit = 1,
#'     rmax = 10,
#'     fits = 48,
#'     cores = 16,
#'     iterlim = 700,
#'     sd = 4
#'   )
#' },
#' jobname = "my_name",
#' partition = "single",
#' cores = 16,
#' nodes = 1,
#' mem_per_core = 2,
#' walltime = "02:00:00",
#' ssh_passwd = "password",
#' machine = "cluster",
#' var_values = NULL,
#' no_rep = 20,
#' recover = F,
#' compile = F,
#' link = F
#' )
#' outDistributedComputing$check()
#' outDistributedComputing$get()
#' outDistributedComputing$purge()
#' result <- cluster_result
#' print(result)
#' 
#' 
#' # calculate profiles
#' var_list <- profileParsPerNode(best_fit, 4)
#' profile_jobname <- paste0(fit_filename,"_profiles_opt")
#' method <- "optimize"
#' profilesDistributedComputing <- distributedComputing(
#'   {
#'     profile(
#'       obj = obj,
#'       pars =  best_fit,
#'       whichPar = (as.numeric(var_1):as.numeric(var_2)),
#'       limits = c(-5, 5),
#'       cores = 16,
#'       method = method,
#'       stepControl = list(
#'         stepsize = 1e-6,
#'         min = 1e-4, 
#'         max = Inf, 
#'         atol = 1e-2,
#'         rtol = 1e-2, 
#'         limit = 100
#'       ),
#'       optControl = list(iterlim = 20)
#'     )
#'   },
#'   jobname = profile_jobname,
#'   partition = "single",
#'   cores = 16,
#'   nodes = 1,
#'   walltime = "02:00:00",
#'   ssh_passwd = "password",
#'   machine = "cluster",
#'   var_values = var_list,
#'   no_rep = NULL,
#'   recover = F,
#'   compile = F,
#'   link = F
#' )
#' profilesDistributedComputing$check()
#' profilesDistributedComputing$get()
#' profilesDistributedComputing$purge()
#' profiles  <- NULL
#' for (i in cluster_result) {
#'   profiles <- rbind(profiles, i)
#' }
#' }
#' 
#' @export
distributedComputing <- function(
    ...,
    jobname,
    partition = "single",
    cores = 16,
    nodes = 1,
    mem_per_core = 2,
    walltime = "01:00:00",
    ssh_passwd = NULL,
    machine = "cluster",
    var_values = NULL,
    no_rep = NULL,
    recover = TRUE,
    purge_local = FALSE,
    compile = FALSE,
    link = FALSE,
    buildCores = NULL,
    buildBundle = 50,
    custom_folders = NULL,
    resetSeeds = TRUE,
    returnAll = TRUE,
    input = ls(.GlobalEnv, all.names = TRUE)
){
  original_wd <- getwd()
  if (is.null(custom_folders)) {
    output_folder_abs <- "./"
  } else if(!is.null(custom_folders) & !all(length(custom_folders) == 3 & sort(names(custom_folders)) == c("compiled", "output", "tmp"))) {
    warning("'custom_folders' must be named vector with exact three elements:\n
            'compiled', 'output', 'tmp', containing relative paths to the resp folders\n
            input is wrong, ignored.\n")
  } else {
    compiled_folder <- custom_folders["compiled"]
    output_folder <- custom_folders["output"]
    tmp_folder <- custom_folders["tmp"]
    
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
  
  # relative path to the working directory, will now allways be used
  
  wd_path <- paste0("./",jobname, "_folder/")
  
  # number of repetitions
  if(!is.null(no_rep) & is.null(var_values)) {
    num_nodes <- no_rep - 1
  } else if(is.null(no_rep) & !is.null(var_values)) {
    num_nodes <- length(var_values[[1]]) - 1
  } else {
    stop("I dont know what you want how often done. Please set either 'no_rep' or pass 'var_values' (_not_ both!)")
  }
  
  # define the ssh command depending on 'sshpass' being used
  ssh_command <- if (is.null(ssh_passwd)) "ssh "
                 else paste0("sshpass -p ", ssh_passwd, " ssh ")
  
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
    .GlobalEnv$cluster_result <- result_list
  }
  
  
  
  # purge function
  out[[3]] <- function (purge_local = FALSE) {
    # remove files remote
    system(
      paste0(
        ssh_command, machine, " rm -rf ", jobname, "_folder"
      )
    )
    # also remove local files if want so
    if (purge_local) {
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
  
  
  # generate list of currently loaded packages
  package_list <- sapply(strsplit(search(), "package:", fixed = TRUE), function(v) v[2])
  package_list <- package_list[!is.na(package_list)]
  package_list <- paste(paste0("try(library(", package_list, "))"), collapse = "\n")
  if (compile || link) {
    objfns <- 'obj.fns <- ls()[sapply(ls(), function(nm) inherits(get(nm, envir=.GlobalEnv), c("obsfn", "parfn", "prdfn")))]'
    setmn <- sprintf('for (o in obj.fns) eval(parse(text=paste0("modelname(", o, ") <- \'%s\'")))\n', paste0(jobname, "_shared_object"))
    load_so <- paste0("dyn.load('",jobname,"_shared_object.so')")
  } else {
    objfns <- ""
    setmn <-""
    load_so <- ""
  }
  
  
  
  # generate parameter lists
  if (!is.null(var_values)) {
    var_list <- paste(
      lapply(
        seq(1,length(var_values)),
        function(i) {
          if (is.character(var_values[[i]])) {
            paste0("var_values_", i, "=c('", paste(var_values[[i]], collapse="','"),"')")
          } else {
            paste0("var_values_", i, "=c(", paste(var_values[[i]], collapse=","),")")
          }
          
          
        }
      ),
      collapse = "\n"
    )
    # cat(variable_list)
    
    # List of all names of parameters that will be changes between runs
    var_names <- paste(lapply(seq(1,length(var_values)), function(i) paste0("var_",i)))
    
    # Variables per run
    var_per_run <- paste(
      lapply(
        seq(1, length(var_values)),
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
  
  # cat(var_per_run)
  
  
  # define fixed pars
  fixedpars <- paste(
    "node_ID = Sys.getenv('SLURM_ARRAY_TASK_ID')",
    "job_ID = Sys.getenv('SLURM_JOB_ID')",
    paste0("jobname = ", "'",jobname,"'"),
    sep = "\n"
  )
  
  
  # WRITE R
  expr <- as.expression(substitute(...))
  cat(
    paste(
      "#!/usr/bin/env Rscript",
      "",
      "# Load packages",
      package_list,
      "try(library(tidyverse))",
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
      paste0("#SBATCH --mem-per-cpu=", mem_per_core, "gb"),
      "",
      "",
      "# Load compiler modules",
      "module load compiler/gnu/13.3",
      "# Load R modules",
      "module load math/R",
      # paste0("export OPENBLAS_NUM_THREADS=",cores),
      paste0("export OMP_NUM_THREADS=","1"), # paste0("export OMP_NUM_THREADS=",cores),
      paste0("export MKL_NUM_THREADS=", "1"), # paste0("export MKL_NUM_THREADS=",cores),
      "",
      "# Run R script",
      paste0("Rscript ", jobname, ".R"),
      sep = "\n" 
    ),
    file = paste0(wd_path,jobname,".sh")
  )
  
  
  # The remote build reads its flags from a script that is shipped inside the
  # job folder; a bare `R CMD SHLIB` misses the cppDE include path and the
  # BLAS/LAPACK libraries. The script is chained with `&&` so a failed build
  # skips the sbatch instead of queueing a job that cannot run.
  build_script_file <- paste0(jobname, "_build.sh")
  filelist_file <- paste0(jobname, "_files.txt")
  module_cmd <- "module load compiler/gnu/13.3 2>/dev/null; module load math/R; "

  ## Names travel as a list, not on the command line: tar and its remote
  ## counterpart share one `system()` string. The list lives in the job folder
  ## and therefore ships inside the same archive.
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

    # Remove any old .so files before linking
    # unlink(list.files(pattern = "(\\.so)$"))

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
  
  ##
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



#' Generate parameter list for distributed profile calculation
#' 
#' @description Generates list of `WhichPar` entries to facillitate distribute
#' profile calculation.
#' @details Lists to split the parameters for which the profiles are calculated
#' on the different nodes.
#' 
#' @param parameters list of parameters 
#' @param fits_per_node numerical, number of parameters that will be send to each node.
#' @param side determine if both sides are calculated (default) or if the profiles are split in 'left' and 'right' for calculation
#' 
#' @return List with two arrays: `from` contains the number of the starting
#' parameter, while `to` stores the respective upper end of the parameter list
#' per node.
#' @examples
#' \dontrun{
#' parameter_list <- setNames(1:10, letters[1:10])
#' var_list <- profileParsPerNode(parameter_list, 4)
#' }
#' 
#' @export
profileParsPerNode <- function(parameters, fits_per_node, side = c("both", "split")[1]) {
  # sanitize side input: must be either "left", "right" or "both"
  if (!(side %in% c("both", "split"))) {
    stop("'side' must be either 'both' or 'split'")
  }
  
  # get the number of parameters
  n_pars <- length(parameters)
  
  # Get number of fits per node
  fits_per_node <- fits_per_node
  
  # determine the number of nodes necessary
  no_nodes <- 1:ceiling(n_pars/fits_per_node)
  
  # generate the lists which parameters are send to which node
  pars_from <- fits_per_node
  pars_to_vec <- fits_per_node
  while (pars_from < (n_pars)) {
    pars_from <- pars_from + fits_per_node
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



## Use Julia to calculate steady states -----------------------------------------



## .sanitizeCores (moved from tools.R) ----------------------------------------

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
 #  
 # if (Sys.info()[['sysname']] == "Windows") cores <- 1
 # return(cores)
  
}




## .parallelLapply (moved from tools.R) --------------------------------------

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

  # A forked worker cannot carry a condition axis: cppDE's batch runs serially
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
