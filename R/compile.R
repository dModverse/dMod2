## compile.R -- C/C++ compilation and DLL (de)registration helpers

## Windows-only: temp Makevars (existing user Makevars + `lines`) for R_MAKEVARS_USER.
.compileMakevarsUser <- function(lines) {
  f <- Sys.getenv("R_MAKEVARS_USER", unset = NA)
  if (is.na(f) || !file.exists(f)) {
    cand <- path.expand(c("~/.R/Makevars.ucrt", "~/.R/Makevars.win64",
                          "~/.R/Makevars.win", "~/.R/Makevars"))
    cand <- cand[file.exists(cand)]
    f <- if (length(cand)) cand[1] else NA
  }
  prev <- if (!is.na(f) && file.exists(f)) readLines(f, warn = FALSE) else character()
  mv <- tempfile(fileext = ".mk")
  writeLines(c(prev, lines), mv)
  mv
}



## Command strings reach the shell as a single argument, so the ceiling is the
## per-argument limit (128 KiB on Linux, 32 KiB on Windows), not ARG_MAX. The
## option lowers it for tests.
.compileCmdLimit <- function() {
  lim <- suppressWarnings(as.integer(getOption("dMod.compile.cmdlimit")))
  if (length(lim) == 1L && !is.na(lim)) return(lim)
  if (.Platform$OS.type == "windows") 24000L else 96000L
}

## Chunks of at most `maxN` elements and `maxChars` quoted characters, each
## short enough for one compiler/archiver invocation.
.compileChunks <- function(x, maxChars = .compileCmdLimit(), maxN = 100L) {
  if (!length(x)) return(list())
  n <- nchar(x) + 3L                    # quotes plus separator
  grp <- integer(length(x)); g <- 1L; acc <- 0L; cnt <- 0L
  for (i in seq_along(x)) {
    if (cnt >= maxN || (cnt > 0L && acc + n[i] > maxChars)) {
      g <- g + 1L; acc <- 0L; cnt <- 0L
    }
    grp[i] <- g; acc <- acc + n[i]; cnt <- cnt + 1L
  }
  unname(split(x, factor(grp, levels = seq_len(g))))
}

## Pull every archive member into the shared object: R resolves the entry points
## by name at run time, so unreferenced members would otherwise be dropped.
## Windows needs --export-all-symbols on top: R writes the export .def with
## `nm` over the objects it is handed -- only the anchor -- and a .def file
## switches ld's auto-export off, so the members would link in unexported.
.compileWholeArchive <- function(lib) {
  if (Sys.info()[["sysname"]] == "Darwin")
    return(paste0("-Wl,-force_load,", shQuote(lib)))
  paste(c(if (.Platform$OS.type == "windows") "-Wl,--export-all-symbols",
          "-Wl,--whole-archive", shQuote(lib), "-Wl,--no-whole-archive"), collapse = " ")
}

## `R CMD config` values, read once per session with a single `--all` call:
## one subprocess per variable costs ~2.5 s on every compile().
.dmodConfig <- new.env(parent = emptyenv())

.compileConfig <- function(var) {
  Rbin <- shQuote(file.path(R.home("bin"), "R"))
  if (!length(ls(.dmodConfig, all.names = TRUE))) {
    out <- tryCatch(system(paste(Rbin, "CMD config --all"), intern = TRUE),
                    error = function(e) "", warning = function(w) "")
    for (m in regmatches(out, regexec("^([A-Za-z_0-9]+) *= ?(.*)$", out)))
      if (length(m) == 3L) assign(m[2], trimws(m[3]), envir = .dmodConfig)
    assign(".read", TRUE, envir = .dmodConfig)
  }
  if (!is.null(v <- .dmodConfig[[var]])) return(v)
  assign(var, trimws(system(paste(Rbin, "CMD config", var), intern = TRUE)),
         envir = .dmodConfig)
  .dmodConfig[[var]]
}

## Include block shared by `sources`, or NULL when any prologue holds more than
## comments and #includes -- prepending it must not change how they expand.
.compilePCHIncludes <- function(sources) {
  includes <- character(0)
  for (f in sources) {
    top <- readLines(f, n = 200L, warn = FALSE)
    inc <- grep("^\\s*#\\s*include", top)
    if (!length(inc)) return(NULL)
    top <- trimws(top[seq_len(max(inc))])
    if (!all(grepl("^$|^//|^/\\*|^\\*|^#\\s*include", top))) return(NULL)
    includes <- union(includes, grep("^#", top, value = TRUE))
  }
  includes
}

## Precompile that block once: the generated sources are a few lines of
## arithmetic around template-heavy headers, so parsing them dominates.
.compilePCH <- function(sources, cmdPrefix, outdir, verbose = FALSE) {
  includes <- .compilePCHIncludes(sources)
  if (is.null(includes)) return(NULL)
  hdr <- file.path(outdir, "dMod_pch.hpp")
  writeLines(includes, hdr)
  cmd <- paste(cmdPrefix, "-x c++-header", shQuote(hdr), "-o", shQuote(paste0(hdr, ".gch")))
  if (verbose) cat(cmd, "\n")
  if (system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE) != 0) return(NULL)
  hdr
}

## cppDE remembers which shared object exports which generated entry point.
## Loading one under a name it has already seen makes that pairing stale.
## Guarded, to stay loadable against a cppDE that predates the function.
.clearSymbols <- function() {
  f <- get0("clearNativeSymbols", envir = asNamespace("cppDE"), inherits = FALSE)
  if (is.function(f)) f()
}

## Reinstalling a backend changes the headers a source compiles against without
## changing the source, so its install stamp belongs in the object cache key.
.headerStamp <- function() {
  paste(vapply(c("cppDE", "cOde"), function(pkg) {
    d <- tryCatch(find.package(pkg), error = function(e) NA_character_)
    if (is.na(d)) return("")
    paste0(pkg, utils::packageVersion(pkg),
           file.info(file.path(d, "DESCRIPTION"))$mtime)
  }, ""), collapse = "|")
}


## How many more shared objects R can dyn.load() (R_MAX_NUM_DLLS, 614 default).
.compileDLLBudget <- function() {
  lim <- suppressWarnings(as.integer(Sys.getenv("R_MAX_NUM_DLLS", "614")))
  if (is.na(lim)) lim <- 614L
  max(0L, lim - length(getLoadedDLLs()))
}


#' Compile model-related C/C++ code
#'
#' @description
#' Compiles model objects ([parfn], [obsfn], [prdfn]) related C/C++ files into shared libraries via `R CMD SHLIB`.
#'
#' @details
#' Named arguments other than the ones below are an error: `...` carries the
#' objects to compile, so a misspelled argument name would otherwise be taken
#' for an object and ignored.
#'
#' The toolchain report prints the flag set every source gets, then each
#' further set with the number of sources that declare it. Compile flags are
#' per file, `-fopenmp` reaching only the ODE sources for instance, so a single
#' bracket would name a command that never runs.
#'
#' Per-file compile and link flags are taken from the `"compileInfo"`
#' attribute that [odemodel()], [Xs()], [Xf()], [Y()] and [P()] attach to
#' their return values. Each entry carries the source file together with the
#' `compileArgs` and `linkArgs` reported by the backend that produced it
#' (`cOde::funC`, `cppDE::cppODE`, `cppDE::cvode`, ...), so solver-specific
#' libraries reach only the files that need them. Objects without
#' `compileInfo` fall back to modelname-based file discovery in the current
#' working directory.
#'
#' @param ... One or more model objects.
#' @param output Optional name for a combined shared library. When set, all
#'   files are linked into one object and the union of their `linkArgs` is
#'   applied. A bare name places the object next to the generated sources; a
#'   name carrying a directory is taken as given.
#' @param args Additional compiler/linker flags applied to every file.
#' @param cores Parallel compilation jobs (Unix only, requires `cores > 1`).
#'   Defaults to [detectFreeCores()].
#' @param chunkSize Maximum number of files per archiver call, see Details.
#' @param verbose If `TRUE`, print compiler commands.
#'
#' @section Many conditions:
#' Sources are always compiled one file per command, but naming every object on
#' the final `R CMD SHLIB` command line overruns the shell's argument limit at
#' roughly a thousand files. Beyond that, `compile()` bundles the objects into a
#' static archive in chunks of `chunkSize` and links it whole; the result is
#' identical. This needs `output` to be set, since one shared object per source
#' would exceed `R_MAX_NUM_DLLS` first.
#'
#' Objects are reused when the source bytes and the compile command are
#' unchanged, recorded in a `.dMod_objects` index next to them: the code
#' generators rewrite every source on each run, so file times say nothing about
#' whether a recompile is needed. Where the sources share their include block,
#' it is precompiled once and prepended, which is what dominates the compile of
#' the short generated files.
#'
#' @return Invisibly `TRUE` on success.
#' @export
compile <- function(..., output = NULL, args = NULL, cores = detectFreeCores(),
                    chunkSize = 100, verbose = FALSE) {

  ## save & restore env
  old <- Sys.getenv(c("PKG_CFLAGS","PKG_CXXFLAGS","PKG_CPPFLAGS","PKG_LIBS"), unset = NA)
  on.exit({
    for (n in names(old))
      if (is.na(old[n])) Sys.unsetenv(n) else Sys.setenv(structure(old[n], names = n))
  }, add = TRUE)

  objs <- list(...)
  if (!length(objs)) stop("No objects")

  # `...` is the payload, so a misspelled argument name lands among the objects
  # and is silently ignored. A named entry that carries no sources is one.
  .hasSources <- function(o)
    !is.null(attr(o, "compileInfo")) || !is.null(attr(o, "srcfile")) ||
    inherits(o, c("obsfn", "parfn", "prdfn", "odemodel"))
  nms <- names(objs)
  if (!is.null(nms)) {
    bad <- which(nzchar(nms) & !vapply(objs, .hasSources, logical(1)))
    if (length(bad)) {
      settable <- setdiff(names(formals(compile)), "...")
      hint <- vapply(nms[bad], function(n) {
        d <- utils::adist(n, settable)[1, ]
        if (min(d) <= 3) paste0("`", n, "`, did you mean `",
                                settable[which.min(d)], "`?")
        else paste0("`", n, "`")
      }, character(1))
      stop("compile: ", paste(hint, collapse = "; "),
           " Objects to compile are passed through `...` and carry sources; ",
           "everything else must match an argument name.", call. = FALSE)
    }
  }
  obj.names <- as.character(substitute(list(...)))[-1]
  Rbin  <- shQuote(file.path(R.home("bin"), "R"))
  so    <- .Platform$dynlib.ext
  cfg   <- .compileConfig
  strip <- function(x) trimws(gsub("(^| )-std=[^ ]+", "", x))

  ## classify objects
  is_dmod <- vapply(objs, inherits, logical(1), c("obsfn","parfn","prdfn"))
  is_cpp  <- vapply(objs, function(o) !is.null(attr(o, "srcfile")), logical(1))

  ## Collect per-file build info.
  ## Primary source is `attr(o, "compileInfo")` carrying
  ## (srcfile, compileArgs, linkArgs) as reported by cOde/cppDE/CVODE.
  ## Falls back to modelname-based file discovery for objects that lack
  ## compileInfo, and to the bare `srcfile` attribute for raw cppDE objects.
  info_from_compileInfo <- unlist(
    lapply(objs, function(o) attr(o, "compileInfo")),
    recursive = FALSE
  )

  info_fallback <- list()
  for (i in seq_along(objs)) {
    o <- objs[[i]]
    if (!is.null(attr(o, "compileInfo"))) next
    if (is_dmod[i]) {
      b <- outer(modelname(o), c("","_deriv","_s","_s2","_sdcv","_dfdx","_dfdp"), paste0)
      cand <- c(paste0(b, ".c"), paste0(b, ".cpp"))
      src <- cand[file.exists(cand)]
      for (s in src)
        info_fallback[[length(info_fallback) + 1]] <-
          list(srcfile = normalizePath(s, winslash = "/", mustWork = TRUE),
               compileArgs = "", linkArgs = "")
    } else if (is_cpp[i]) {
      s <- attr(o, "srcfile")
      if (length(s) && nzchar(s) && file.exists(s))
        info_fallback[[length(info_fallback) + 1]] <- list(
          srcfile     = normalizePath(s, winslash = "/", mustWork = TRUE),
          compileArgs = attr(o, "compileArgs") %||% "",
          linkArgs    = attr(o, "linkArgs")    %||% "",
          sparse      = isTRUE(attr(o, "sparse")))
    }
  }

  info <- c(info_from_compileInfo, info_fallback)

  ## Expand entries with multiple srcfiles (e.g. cOde spills _deriv.c
  ## alongside the main .c) into one entry per file.
  info <- unlist(lapply(info, function(e) {
    if (!length(e$srcfile)) return(list())
    if (length(e$srcfile) == 1L) return(list(e))
    lapply(e$srcfile, function(s) list(srcfile = s, compileArgs = e$compileArgs,
                                       linkArgs = e$linkArgs, sparse = e$sparse))
  }), recursive = FALSE)

  info <- Filter(function(e) length(e$srcfile) == 1L && nzchar(e$srcfile) && file.exists(e$srcfile), info)
  if (!length(info)) stop("No source files found")

  ## Deduplicate by srcfile, keeping the first (non-empty) flags we saw.
  ord <- order(vapply(info, function(e) e$srcfile, character(1)))
  info <- info[ord]
  keep <- !duplicated(vapply(info, function(e) e$srcfile, character(1)))
  info <- info[keep]

  files      <- vapply(info, function(e) e$srcfile, character(1))

  ## Where the combined shared object goes. A bare name lands next to the
  ## sources, a name carrying a directory is honoured as given.
  outfile <- NULL
  if (!is.null(output)) {
    if (length(output) != 1L || !is.character(output) || !nzchar(output))
      stop("compile: `output` must be a single non-empty name.", call. = FALSE)
    stem <- sub(paste0("\\", so, "$"), "", output)
    dir  <- if (basename(stem) == stem) dirname(files[1]) else dirname(stem)
    if (!dir.exists(dir))
      stop("compile: the directory of `output` does not exist: ", dir, call. = FALSE)
    # Never link into a shared object this process already holds: unloading it
    # is not portable, Windows may keep the file handle and macOS may keep the
    # image resident, so the reload would serve the old code.
    out_base <- .uniqueLibname(basename(stem))
    outfile  <- file.path(dir, paste0(out_base, so))
  }
  roots      <- sub("\\.[^.]+$", "", basename(files))
  roots_full <- sub("\\.[^.]+$", "", files)

  ## compiler flags
  if (.Platform$OS.type == "windows") cores <- 1
  pic  <- if (.Platform$OS.type == "windows") "" else "-fPIC"
  base <- paste("-O2 -DNDEBUG -w", pic)

  ## KLU flags for sparse models, mirroring cppDE::compile(). The flag lives on
  ## the cppDE object inside an odemodel, so a bare `attr(o, "sparse")` on the
  ## dMod fn objects handed to compile() never sees it -- go through the
  ## per-file info. The `-DKLU*` fallback covers objects whose compileInfo was
  ## built before `sparse` was recorded there.
  uses_klu <- any(vapply(objs, function(o) isTRUE(attr(o, "sparse")), logical(1))) ||
    any(vapply(info, function(e) isTRUE(e$sparse) ||
                 grepl("(^|\\s)-DKLU", e$compileArgs %||% ""), logical(1)))
  klu_flag <- ""; klu_lib <- ""
  if (uses_klu) {
    cfgCppDE <- .cppDE_config()
    if (!isTRUE(cfgCppDE$klu_available))
      stop("A sparse Jacobian was requested, but cppDE was installed without the ",
           "KLU linear solver.\n  Install SuiteSparse/KLU and re-install cppDE -- ",
           "see cppDE::install_libs(\"suitesparse\").", call. = FALSE)
    klu_flag <- trimws(paste("-DKLU", cfgCppDE$klu_cflags))
    klu_lib  <- cfgCppDE$klu_libs
  }

  ## shared pieces (compiler/linker) that apply to every file
  cxx_base <- paste(base, klu_flag)
  extra_args <- paste(c(args), collapse = " ")
  if (nzchar(extra_args)) {
    base     <- paste(base,     extra_args)
    cxx_base <- paste(cxx_base, extra_args)
  }
  ## BLAS/LAPACK: on Windows `R CMD config BLAS_LIBS` returns a value with
  ## unexpanded `$(R_HOME)`/`$(R_ARCH)` references. Those go into PKG_LIBS as
  ## an env var, and make should re-expand them, but in practice the
  ## expansion is unreliable inside SHLIB-generated link commands -- the
  ## final g++ invocation comes out without any BLAS libs. We sidestep that
  ## by building an absolute -L path here and skipping `R CMD config`.
  if (.Platform$OS.type == "windows") {
    ## R.home("bin") already resolves to the arch-specific bin dir (.../bin/x64)
    ## on R >= 4.2, which is where Rblas.dll / Rlapack.dll live. Appending
    ## .Platform$r_arch again produced a bogus .../bin/x64/x64 -L path (it only
    ## ever linked by accident, via the default -L.../bin/x64 that -lR adds).
    ## dMod requires R >= 4.5.0, so R.home("bin") is always the correct dir.
    r_bin   <- R.home("bin")
    blaslapack <- paste0("-L", shQuote(r_bin), " -lRlapack -lRblas")
  } else {
    blaslapack <- paste(cfg("LAPACK_LIBS"), cfg("BLAS_LIBS"))
  }
  base_libs <- paste(klu_lib, blaslapack)
  cppflags  <- paste0("-I", system.file("include", package = "cppDE"))

  ## Compiler invocation bits cached up front so parallel forks don't each
  ## re-spawn R-CMD-config. Used by compile_one_obj() for the direct
  ## $CC/$CXX -c path.
  cc_bin      <- cfg("CC")
  cxx_bin     <- cfg("CXX")
  cflags_R    <- cfg("CFLAGS")
  cxxflags_R  <- cfg("CXXFLAGS")
  cpicflags   <- cfg("CPICFLAGS")
  cxxpicflags <- cfg("CXXPICFLAGS")
  r_inc       <- paste0("-I", shQuote(R.home("include")))

  ## Toolchain report. Compile flags are per file, so this prints the set every
  ## source gets and then each further set with the number declaring it.
  .tok <- function(x) {
    t <- unlist(strsplit(trimws(paste(x, collapse = " ")), "[[:space:]]+"))
    t[nzchar(t)]
  }
  .reportToolchain <- function(label, bin, shared, pattern) {
    ent <- Filter(function(e) grepl(pattern, e$srcfile, ignore.case = TRUE), info)
    if (!length(ent)) return(invisible(NULL))

    sh    <- unique(.tok(shared))
    extra <- vapply(ent, function(e)
      paste(setdiff(unique(.tok(e$compileArgs %||% "")), sh), collapse = " "),
      character(1))

    if (!any(nzchar(extra))) {
      cat(sprintf("using %-3s compiler: %s [%s]\n", label, strip(bin), trimws(shared)))
      return(invisible(NULL))
    }
    sets  <- setdiff(unique(extra), "")
    tags  <- c("every source", sprintf("%d of %d also", vapply(sets, function(k)
      sum(extra == k), integer(1)), length(ent)))
    width <- max(nchar(tags))
    cat(sprintf("using %-3s compiler: %s\n", label, strip(bin)))
    cat(sprintf("  %-*s : %s\n", width, tags[1], trimws(shared)))
    for (i in seq_along(sets))
      cat(sprintf("  %-*s : %s\n", width, tags[i + 1L], sets[i]))
    invisible(NULL)
  }
  .reportToolchain("C",   cc_bin,  base,     "\\.c$")
  .reportToolchain("C++", cxx_bin, cxx_base, "\\.cpp$")

  ## unload stale DLLs
  loaded <- getLoadedDLLs()
  for (i in seq_along(roots))
    if (roots[i] %in% names(loaded)) try(dyn.unload(loaded[[roots[i]]][["path"]]), silent = TRUE)
  if (!is.null(outfile)) try(dyn.unload(outfile), silent = TRUE)

  ## Compile one file with its own compile/link flags applied via PKG_*.
  ## Each invocation sets the env just before shelling out to R CMD SHLIB,
  ## so per-file linkArgs (e.g. Sundials libs for CVODE) reach only the
  ## files that need them. Works inside mclapply because each fork has its
  ## own env.
  # A backend's compileArgs can repeat a flag the base already sets (-fopenmp).
  # Compile flags only: repeated -l on a link line can be load-bearing.
  .mergeFlags <- function(...) {
    tok <- unlist(strsplit(trimws(paste(...)), "[[:space:]]+"))
    tok <- tok[nzchar(tok)]
    if (any(tok %in% c("-I", "-D", "-L", "-U", "-include", "-isystem")))
      return(paste(tok, collapse = " "))
    paste(unique(tok), collapse = " ")
  }

  compile_one <- function(entry) {
    extra_c <- entry$compileArgs %||% ""
    pkg_c  <- .mergeFlags(base,     extra_c)
    pkg_cx <- .mergeFlags(cxx_base, extra_c)
    pkg_l  <- trimws(paste(base_libs, entry$linkArgs %||% ""))
    Sys.setenv(
      PKG_CFLAGS   = pkg_c,
      PKG_CXXFLAGS = pkg_cx,
      PKG_CPPFLAGS = cppflags,
      PKG_LIBS     = pkg_l
    )
    if (.Platform$OS.type == "windows") {
      mv <- .compileMakevarsUser(c(
        paste("PKG_CFLAGS =",   pkg_c),
        paste("PKG_CXXFLAGS =", pkg_cx),
        paste("PKG_CPPFLAGS =", cppflags),
        paste("PKG_LIBS =",     pkg_l)
      ))
      old_mu <- Sys.getenv("R_MAKEVARS_USER", unset = NA)
      Sys.setenv(R_MAKEVARS_USER = mv)
      on.exit({
        if (is.na(old_mu)) Sys.unsetenv("R_MAKEVARS_USER")
        else Sys.setenv(R_MAKEVARS_USER = old_mu)
        unlink(mv)
      }, add = TRUE)
    }
    cmd <- paste(Rbin, "CMD SHLIB", shQuote(entry$srcfile))
    if (verbose) cat(cmd, "\n")
    if (system(cmd, ignore.stdout = !verbose, ignore.stderr = !verbose) != 0)
      stop("Compilation failed: ", entry$srcfile)
  }

  ## Command compiling a single source to a .o via a direct $CC/$CXX -c call.
  ## The combined-output path runs these in parallel; the subsequent
  ## R CMD SHLIB link then sees the objects are up to date and only links.
  obj_cmd <- function(entry, pch = NULL) {
    src     <- entry$srcfile
    extra_c <- entry$compileArgs %||% ""
    obj     <- sub("\\.[^.]+$", ".o", src)
    if (grepl("\\.cpp$", src, ignore.case = TRUE))
      paste(cxx_bin, r_inc, cppflags, .mergeFlags(cxx_base, extra_c),
            if (!is.null(pch)) paste("-Winvalid-pch -include", shQuote(pch)),
            cxxpicflags, cxxflags_R, "-c", shQuote(src), "-o", shQuote(obj))
    else
      paste(cc_bin, r_inc, cppflags, .mergeFlags(base, extra_c),
            cpicflags, cflags_R, "-c", shQuote(src), "-o", shQuote(obj))
  }

  compile_one_obj <- function(job) {
    if (verbose) cat(job$cmd, "\n")
    if (system(job$cmd, ignore.stdout = !verbose, ignore.stderr = !verbose) != 0)
      stop("Compilation failed: ", job$srcfile)
    job$srcfile
  }

  if (is.null(output)) {
    ## One shared object per source, all dyn.load()ed -- refuse up front rather
    ## than failing halfway through a long build.
    budget <- .compileDLLBudget()
    if (length(info) > budget)
      stop(length(info), " source files would need as many shared objects, but R can ",
           "load at most ", budget, " more (R_MAX_NUM_DLLS).\n  Pass output = <name> ",
           "to link them into a single shared object.", call. = FALSE)
    if (.Platform$OS.type == "unix" && cores > 1)
      parallel::mclapply(info, compile_one, mc.cores = cores)
    else for (e in info) compile_one(e)
    for (r in roots_full) .reloadDLL(paste0(r, so))
    .clearSymbols()
  } else {
    ## Combined output: per-file compile to .o (parallel on Unix when cores>1,
    ## serial otherwise -- including on Windows), then a single R CMD SHLIB
    ## link over the original sources. Because every .o is freshly written
    ## above, make sees them as up-to-date and only runs the link recipe;
    ## passing the source list lets SHLIB pick the C++ linker when any
    ## source is .cpp, which a .o-only invocation would miss. The pre-compile
    ## also has to run on Windows: the single-call SHLIB (compile + link in
    ## one go) was occasionally producing .dll files that LoadLibrary
    ## couldn't resolve when the source pulled in BLAS via the symbolic-
    ## mode chain wrapper -- splitting compile and link sidesteps that.
    ## A precompiled header is keyed to one flag set, so it is only built when
    ## every C++ source shares its compile arguments and there are enough of
    ## them to amortise it.
    is_cxx <- grepl("\\.cpp$", files, ignore.case = TRUE)
    cxxArgs <- unique(vapply(info[is_cxx], function(e) e$compileArgs %||% "", ""))
    pch <- NULL
    if (sum(is_cxx) >= 8L && length(cxxArgs) == 1L) {
      pchdir <- tempfile("dMod_pch"); dir.create(pchdir, showWarnings = FALSE)
      on.exit(unlink(pchdir, recursive = TRUE), add = TRUE)
      pch <- .compilePCH(files[is_cxx],
                         paste(cxx_bin, r_inc, cppflags, .mergeFlags(cxx_base, cxxArgs),
                               cxxpicflags, cxxflags_R), pchdir, verbose)
    }

    ## Reuse objects whose source bytes and compile command are unchanged: the
    ## generators rewrite every source each run, so file times say nothing.
    ## Touching a reused object keeps it newer than its source for the link.
    objs <- sub("\\.[^.]+$", ".o", files)
    keys <- structure(paste(tools::md5sum(files), vapply(info, obj_cmd, ""),
                            .headerStamp()), names = files)
    cachefile <- file.path(dirname(files[1]), ".dMod_objects")
    prev <- if (file.exists(cachefile))
      tryCatch(readRDS(cachefile), error = function(e) NULL) else NULL
    hit <- if (length(prev)) match(files, names(prev)) else rep(NA_integer_, length(files))
    fresh <- !is.na(hit) & file.exists(objs)
    fresh[fresh] <- prev[hit[fresh]] == keys[fresh]
    if (any(fresh)) {
      Sys.setFileTime(objs[fresh], Sys.time())
      cat(sprintf("reusing %d unchanged object(s)\n", sum(fresh)))
    }

    jobs <- lapply(which(!fresh), function(i)
      list(srcfile = files[i], cmd = obj_cmd(info[[i]], if (is_cxx[i]) pch)))
    res <- if (.Platform$OS.type == "unix" && cores > 1)
      parallel::mclapply(jobs, compile_one_obj, mc.cores = cores)
    else lapply(jobs, compile_one_obj)
    ## mclapply returns a failing fork as a try-error instead of raising it.
    bad <- vapply(res, inherits, logical(1), "try-error")
    if (any(bad)) stop(as.character(res[bad][[1]]), call. = FALSE)
    try(saveRDS(keys, cachefile), silent = TRUE)

    ## Link step: union of every entry's linkArgs (dedup) so Sundials-dependent
    ## files still pull their libs.
    all_link <- unique(unlist(lapply(info, function(e) strsplit(trimws(e$linkArgs %||% ""), "\\s+")[[1]])))
    all_link <- all_link[nzchar(all_link)]
    all_compile <- unique(unlist(lapply(info, function(e) strsplit(trimws(e$compileArgs %||% ""), "\\s+")[[1]])))
    all_compile <- all_compile[nzchar(all_compile)]

    output <- basename(sub(paste0("\\", so, "$"), "", output))

    ## Past the argument limit the objects go into a static archive, appended in
    ## chunks, and SHLIB gets one anchor source plus the archive. The anchor is a
    ## C++ source when there is one: SHLIB picks the C++ linker from the sources
    ## it is handed, not from the archive contents.
    link_files <- files
    if (nchar(paste(shQuote(files), collapse = " ")) > .compileCmdLimit()) {
      anchor <- if (any(is_cxx)) which(is_cxx)[1] else 1L
      lib    <- file.path(dirname(files[1]), paste0(output, "_objects.a"))
      unlink(lib)
      on.exit(try(unlink(lib), silent = TRUE), add = TRUE)
      ## Via `R CMD`, which puts the Rtools toolchain on PATH under Windows.
      ar_bin     <- cfg("AR");     if (!nzchar(ar_bin))     ar_bin     <- "ar"
      ranlib_bin <- cfg("RANLIB"); if (!nzchar(ranlib_bin)) ranlib_bin <- "ranlib"
      chunks <- .compileChunks(objs[-anchor], maxN = as.integer(chunkSize))
      for (i in seq_along(chunks)) {
        ## `q` appends without an index; ranlib writes it once at the end.
        cmd <- paste(Rbin, "CMD", ar_bin, if (i == 1L) "qc" else "q", shQuote(lib),
                     paste(shQuote(chunks[[i]]), collapse = " "))
        if (verbose) cat(cmd, "\n")
        if (system(cmd, ignore.stdout = !verbose, ignore.stderr = !verbose) != 0)
          stop("Archiving failed at chunk ", i, " of ", length(chunks), call. = FALSE)
      }
      if (system(paste(Rbin, "CMD", ranlib_bin, shQuote(lib)),
                 ignore.stdout = !verbose, ignore.stderr = !verbose) != 0)
        stop("Building the archive index failed: ", lib, call. = FALSE)
      cat(sprintf("archived %d objects into %s (%d chunks)\n",
                  length(objs) - 1L, basename(lib), length(chunks)))
      link_files <- files[anchor]
      base_libs  <- paste(.compileWholeArchive(lib), base_libs)
    }

    pkg_cflags   <- trimws(paste(base,     paste(all_compile, collapse = " ")))
    pkg_cxxflags <- trimws(paste(cxx_base, paste(all_compile, collapse = " ")))
    pkg_libs     <- trimws(paste(base_libs, paste(all_link, collapse = " ")))
    Sys.setenv(
      PKG_CFLAGS   = pkg_cflags,
      PKG_CXXFLAGS = pkg_cxxflags,
      PKG_CPPFLAGS = cppflags,
      PKG_LIBS     = pkg_libs
    )

    ## Belt-and-suspenders for Windows: env-imported PKG_LIBS has been
    ## observed to vanish from SHLIB's generated link command on some R/rtools
    ## combinations, leaving the .dll unlinked against BLAS/LAPACK. Drop a
    ## per-link Makevars(.win) alongside the source files so make picks it up
    ## even if the environment doesn't make it through. We clean it up after
    ## the link so the directory state stays hermetic.
    mv_dir  <- dirname(files[1])
    mv_name <- if (.Platform$OS.type == "windows") "Makevars.win" else "Makevars"
    mv_path <- file.path(mv_dir, mv_name)
    mv_pre  <- if (file.exists(mv_path)) readLines(mv_path, warn = FALSE) else NULL
    writeLines(c(
      paste("PKG_CFLAGS =",   pkg_cflags),
      paste("PKG_CXXFLAGS =", pkg_cxxflags),
      paste("PKG_CPPFLAGS =", cppflags),
      paste("PKG_LIBS =",     pkg_libs)
    ), mv_path)
    on.exit({
      if (is.null(mv_pre)) try(unlink(mv_path), silent = TRUE)
      else                 try(writeLines(mv_pre, mv_path), silent = TRUE)
    }, add = TRUE)

    ## Windows fallback for BLAS/LAPACK: inject PKG_* via R_MAKEVARS_USER.
    if (.Platform$OS.type == "windows") {
      mv <- .compileMakevarsUser(c(
        paste("PKG_CFLAGS =",   pkg_cflags),
        paste("PKG_CXXFLAGS =", pkg_cxxflags),
        paste("PKG_CPPFLAGS =", cppflags),
        paste("PKG_LIBS =",     pkg_libs)
      ))
      old_mu <- Sys.getenv("R_MAKEVARS_USER", unset = NA)
      Sys.setenv(R_MAKEVARS_USER = mv)
      on.exit({
        if (is.na(old_mu)) Sys.unsetenv("R_MAKEVARS_USER")
        else Sys.setenv(R_MAKEVARS_USER = old_mu)
        unlink(mv)
      }, add = TRUE)
    }

    out <- outfile
    try(dyn.unload(out), silent = TRUE)
    if (file.exists(out)) unlink(out)
    ## Link, capturing stdout+stderr via system2() pipes. Strip the compiler
    ## banner afterwards: the .o files are already fresh from compile_one_obj,
    ## so make only runs the link recipe and "using C/C++ compiler:" would be
    ## misleading. We pass the streams through system2(stdout/stderr = TRUE)
    ## rather than appending a `2>&1` token to a command string: on Windows
    ## that trailing token is not consumed by a shell but swallowed by
    ## R CMD SHLIB as the make override `PKG_LIBS=2>&1`, which beats every
    ## Makevars/R_MAKEVARS_USER assignment and strips the BLAS/LAPACK libs
    ## (breaking the symbolic-mode chain_jac link with "undefined reference to
    ## dgemm_"). system2() keeps the argument vector clean and is identical on
    ## Linux/macOS, where the shell would have consumed the redirection anyway.
    Rexe <- file.path(R.home("bin"), "R")
    shlib_args <- c("CMD", "SHLIB", shQuote(link_files), "-o", shQuote(out))
    if (verbose) cat(shQuote(Rexe), paste(shlib_args, collapse = " "), "\n")
    out_lines <- suppressWarnings(
      system2(Rexe, shlib_args, stdout = TRUE, stderr = TRUE)
    )
    status <- attr(out_lines, "status")
    if (verbose) {
      out_lines <- out_lines[!grepl("^using (C|C\\+\\+) compiler:", out_lines)]
      writeLines(out_lines)
    }
    if (!is.null(status) && status != 0L)
      stop("Compilation failed:\n", paste(out_lines, collapse = "\n"))
    if (!file.exists(out))
      stop("R CMD SHLIB returned exit 0 but did not produce ", out, ":\n",
           paste(out_lines, collapse = "\n"))
    .reloadDLL(out)
    .clearSymbols()
    ## Only arguments passed as a plain variable can have their modelname
    ## updated in the caller; an expression has nothing to assign back to.
    for (i in which(is_dmod))
      if (make.names(obj.names[i]) == obj.names[i])
        eval.parent(parse(text = paste0("modelname(", obj.names[i], ") <- '", output, "'")))
  }

  invisible(TRUE)
}




#' Determine loaded DLLs available in working directory
#' 
#' @return Character vector with the names of the loaded DLLs available in the working directory
#' @export
getLocalDLLs <- function() {
  
  all.dlls <- getLoadedDLLs()
  is.local <- sapply(all.dlls, function(x) grepl(getwd(), unclass(x)$path, fixed = TRUE))
  names(is.local)[is.local]
  
}




## Loading a shared object that is already loaded is a no-op in R, so a rebuilt
## file would keep serving the old code. Entry points resolve by name, so
## unloading first is safe even for objects still in use.
## A base name no loaded shared library carries: the desired one if it is
## free, otherwise the smallest `<name>_<i>`, i >= 2, with a warning. Mirrors
## cppDE's `unique_modelname()`, which the model constructors apply to their
## own names; the two must stay in step.
.uniqueLibname <- function(name) {
  loaded <- names(getLoadedDLLs())
  if (!name %in% loaded) return(name)
  i <- 2L
  repeat {
    cand <- paste0(name, "_", i)
    if (!cand %in% loaded) {
      warning(sprintf(
        "A shared library named '%s' is already loaded; overwriting it is not portable. Using '%s' instead.",
        name, cand), call. = FALSE)
      return(cand)
    }
    i <- i + 1L
  }
}


.reloadDLL <- function(path) {
  p <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (p %in% .loadedDLLPaths()) try(dyn.unload(p), silent = TRUE)
  dyn.load(p)
}


## Absolute paths of the shared objects loaded in this process.
.loadedDLLPaths <- function() {
  dlls <- getLoadedDLLs()
  if (!length(dlls)) return(character(0))
  paths <- vapply(dlls, function(d) unclass(d)$path, character(1))
  normalizePath(paths, winslash = "/", mustWork = FALSE)
}

## Directories to search for an object's shared libraries: where its sources
## were generated, plus the working directory. A model compiled into a temp
## folder is not findable from the modelname alone.
.dllSearchDirs <- function(objects) {
  dirs <- unlist(lapply(objects, function(o)
    vapply(attr(o, "compileInfo") %||% list(),
           function(e) dirname(e$srcfile[1]), character(1))))
  unique(c(getwd(), dirs[nzchar(dirs)]))
}


#' Load shared object for a dMod object
#'
#' Usually when restarting the R session, although all objects are saved in
#' the workspace, the dynamic libraries are not linked any more. `loadDLL`
#' is a wrapper for `dyn.load` that uses the "modelname" attribute of
#' dMod objects like prediction functions, observation functions, etc. to
#' load the corresponding shared object. Searched are the directories the
#' objects were generated in and the working directory.
#'
#' Shared objects already loaded in the current process are skipped, so
#' calling `loadDLL` repeatedly is a no-op.
#'
#' @param ... objects of class prdfn, obsfn, parfn, objfn, ...
#'
#' @return Invisibly, the character vector of files loaded by this call.
#'
#' @export
loadDLL <- function(...) {

  .so    <- .Platform$dynlib.ext
  models <- modelname(...)
  names  <- paste0(outer(models, c("", "_s", "_s2", "_sdcv", "_deriv", "_dfdx", "_dfdp"),
                         paste0), .so)
  files  <- as.vector(outer(.dllSearchDirs(list(...)), names, file.path))
  files  <- normalizePath(files[file.exists(files)], winslash = "/", mustWork = FALSE)
  files  <- setdiff(unique(files), .loadedDLLPaths())
  if (!length(files)) return(invisible(character(0)))

  for (f in files) dyn.load(f)
  .clearSymbols()
  message("The following local files were dynamically loaded: ", paste(files, collapse = ", "))
  invisible(files)
}


## compileInfo plumbing (moved from classes.R) ----------------------------------------

## ODE model class -------------------------------------------------------------------

## Collect per-sub-object build info from ODE model pieces.
## Each backend (cOde::funC, cppDE::cppODE, cppDE::cvode) may attach
## `srcfile`, `compileArgs`, and `linkArgs` to its func/extended result. For
## backends that don't (cOde), we fall back to modelname-based file discovery
## in the current working directory. The resulting list is the single
## authoritative source consulted by `compile()` when given dMod fn objects.

## Dedup keys, cached on the list: merging is pairwise, so recomputing them
## every time would make summing a few thousand conditions quadratic.
.compileInfoKeys <- function(x) {
  k <- attr(x, "srckeys")
  if (!is.null(k) && length(k) == length(x)) return(k)
  vapply(x, function(e) if (length(e$srcfile)) paste(e$srcfile, collapse = "\x1f")
                        else NA_character_, character(1))
}

## Merge two compileInfo lists, deduplicating by srcfile (per file, the first
## occurrence wins -- that keeps the originating compile/link flags). Returns
## NULL when both inputs are empty so the attribute stays absent on objects
## that never had native code to begin with.
.mergeCompileInfo <- function(a, b) {
  if (!length(a) && !length(b)) return(NULL)
  ka <- .compileInfoKeys(a)
  kb <- .compileInfoKeys(b)
  keep_a <- !is.na(ka) & !duplicated(ka)
  keep_b <- !is.na(kb) & !duplicated(kb) & !(kb %in% ka[keep_a])
  out  <- c(a[keep_a], b[keep_b])
  if (!length(out)) return(NULL)
  attr(out, "srckeys") <- c(ka[keep_a], kb[keep_b])
  out
}

.collectCompileInfo <- function(...) {
  objs <- list(...)
  objs <- objs[!vapply(objs, is.null, logical(1))]
  out <- list()
  for (o in objs) {
    src <- attr(o, "srcfile")
    if (is.null(src) || !length(src) || !nzchar(src)) {
      mname <- attr(o, "modelname")
      if (is.null(mname) && is.character(o)) mname <- unname(o[1])
      if (is.null(mname) || !nzchar(mname)) next
      b <- outer(mname, c("", "_deriv", "_s", "_s2", "_sdcv", "_dfdx", "_dfdp"), paste0)
      cand <- c(paste0(b, ".c"), paste0(b, ".cpp"))
      src <- cand[file.exists(cand)]
      if (!length(src)) next
      src <- normalizePath(src, winslash = "/", mustWork = FALSE)
    } else {
      src <- normalizePath(src, winslash = "/", mustWork = FALSE)
    }
    out[[length(out) + 1]] <- list(
      srcfile     = src,
      compileArgs = attr(o, "compileArgs") %||% "",
      linkArgs    = attr(o, "linkArgs")    %||% "",
      ## cppDE flags a sparse-Jacobian model here; the flag has to travel with
      ## the file because `compile()` only ever sees the wrapping dMod fn.
      sparse      = isTRUE(attr(o, "sparse"))
    )
  }
  out
}


## odemodel() constructor lives in R/odeClass.R.


