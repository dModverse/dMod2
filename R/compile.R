## compile.R, C/C++ compilation and DLL (de)registration helpers

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



## Longest command line: the per-argument limit on Linux (128 KiB), cmd.exe's
## 8191 characters on Windows, where `R CMD` runs through it. The option lowers
## it for tests.
.compileCmdLimit <- function() {
  lim <- suppressWarnings(as.integer(getOption("dMod.compile.cmdlimit")))
  if (length(lim) == 1L && !is.na(lim)) return(lim)
  if (.Platform$OS.type == "windows") 8000L else 96000L
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

## Pull every archive member into the shared object, as R resolves entry points
## by name at run time. On Windows R's .def lists only the anchor's symbols and
## turns off auto-export, hence --export-all-symbols.
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
## comments and #includes, prepending it must not change how they expand.
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
## Guarded for a cppDE without the function.
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


## Point every cOde model inside `x` at the shared object `output`, walking the
## evaluation tree since `mappings` may lag behind the leaves.
.retargetCode <- function(x, output) {
  st <- .fnNode(x)
  if (is.null(st)) return(invisible(NULL))
  if (identical(st$op, "leaf")) {
    e <- environment(st$kernel)
    if (is.null(e)) return(invisible(NULL))
    for (nm in c("func", "extended")) {
      o <- get0(nm, envir = e, inherits = FALSE)
      if (is.null(o) || is.null(attr(o, "modelname"))) next
      attr(o, "modelname") <- output
      assign(nm, o, envir = e)
    }
    return(invisible(NULL))
  }
  for (nm in c("p1", "p2")) {
    k <- st[[nm]]
    if (is.null(k)) next
    for (kk in if (is.list(k)) k else list(k)) .retargetCode(kk, output)
  }
  invisible(NULL)
}


#' Compile Model-Related C/C++ Code
#'
#' @description
#' Compiles the generated sources of model objects ([parfn], [obsfn],
#' [prdfn], [odemodel()]) into shared objects with `R CMD SHLIB` and loads
#' them.
#'
#' @details
#' Compile and link flags are taken per source from the `"compileInfo"`
#' attribute of the objects. Objects without it are matched to sources in the
#' working directory by model name. Unchanged sources are not recompiled.
#'
#' Without `output`, each source becomes a shared object named after it, i.e.
#' after the model name. With `output`, all sources are linked into one shared
#' object of that name, and the model name of every object passed as a plain
#' variable is set to `output`, see [modelname()]. If a shared object of that
#' name is already loaded, a suffix `_2`, `_3`, ... is appended with a warning.
#' With many sources, set `output`: R limits the number of loaded shared
#' objects.
#'
#' Named arguments other than the ones below are an error.
#'
#' @param ... one or more model objects.
#' @param output character, the name of a combined shared object. Default
#'   `NULL`. A bare name places it next to the generated sources, a name with
#'   a directory is taken as given.
#' @param args character, additional compiler and linker flags for every
#'   file. Default `NULL`.
#' @param cores integer, parallel compilation jobs (Unix only). Default
#'   [detectFreeCores()].
#' @param chunkSize integer, the maximum number of object files per archiver
#'   call when the link command would exceed the command-line limit. Default
#'   100.
#' @param verbose logical, print the compiler commands. Default `FALSE`.
#'
#' @return `TRUE`, invisibly.
#' @seealso [odemodel()], [loadDLL()], [modelname()]
#' @examples
#' \donttest{
#' g <- Y(c(y = "s * x"), parameters = "s", modelname = "compile_ex_obs",
#'        outdir = tempdir())
#' p <- P(c(s = "exp(log_s)"), modelname = "compile_ex_par",
#'        outdir = tempdir())
#' compile(g, p, output = "compile_ex", cores = 1)
#' modelname(g)
#' }
#' @export
compile <- function(..., output = NULL, args = NULL, cores = detectFreeCores(),
                    chunkSize = 100, verbose = FALSE) {

  ## save & restore env; R CMD check asks SHLIB for symbol tables, which it
  ## writes as symbols.rds into the working directory
  old <- Sys.getenv(c("PKG_CFLAGS", "PKG_CXXFLAGS", "PKG_CPPFLAGS", "PKG_LIBS",
                      "_R_SHLIB_BUILD_OBJECTS_SYMBOL_TABLES_"), unset = NA)
  on.exit({
    for (n in names(old))
      if (is.na(old[n])) Sys.unsetenv(n) else do.call(Sys.setenv, as.list(old[n]))
  }, add = TRUE)
  Sys.setenv(`_R_SHLIB_BUILD_OBJECTS_SYMBOL_TABLES_` = "FALSE")

  objs <- list(...)
  if (!length(objs)) stop("No objects")

  # `...` is the payload, so a misspelled argument name lands among the objects
  # and is silently ignored. A named entry that has no sources is one.
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
           " Objects to compile are passed through `...` and contain sources; ",
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
  # `objs` is reused for the object-file paths further down, so the function
  # objects are kept aside for the post-link retargeting.
  fn_objs <- objs[is_dmod]
  is_cpp  <- vapply(objs, function(o) !is.null(attr(o, "srcfile")), logical(1))

  ## Per-file build info from `attr(o, "compileInfo")`, else from modelname-based
  ## file discovery, or the bare `srcfile` attribute of raw cppDE objects.
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

  ## Deduplicate by srcfile, keeping the first non-empty flags.
  ord <- order(vapply(info, function(e) e$srcfile, character(1)))
  info <- info[ord]
  keep <- !duplicated(vapply(info, function(e) e$srcfile, character(1)))
  info <- info[keep]

  files      <- vapply(info, function(e) e$srcfile, character(1))

  ## Where the combined shared object goes. A bare name lands next to the
  ## sources, a name that includes a directory is honoured as given.
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

  ## KLU flags for sparse models, mirroring cppDE::compile(). The flag is read
  ## from the per-file info, as fn objects do not hold it; `-DKLU` in the
  ## compile arguments counts as well.
  uses_klu <- any(vapply(objs, function(o) isTRUE(attr(o, "sparse")), logical(1))) ||
    any(vapply(info, function(e) isTRUE(e$sparse) ||
                 grepl("(^|\\s)-DKLU", e$compileArgs %||% ""), logical(1)))
  klu_flag <- ""; klu_lib <- ""
  if (uses_klu) {
    cfgCppDE <- .cppDE_config()
    if (!isTRUE(cfgCppDE$klu_available))
      stop("A sparse Jacobian was requested, but cppDE was installed without the ",
           "KLU linear solver.\n  Install SuiteSparse/KLU and re-install cppDE, ",
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
  ## On Windows `R CMD config BLAS_LIBS` holds unexpanded `$(R_HOME)` references
  ## that SHLIB's link command does not reliably expand, so the -L path is built
  ## here. R.home("bin") is the arch-specific dir holding Rblas.dll and Rlapack.dll.
  if (.Platform$OS.type == "windows") {
    r_bin   <- R.home("bin")
    blaslapack <- paste0("-L", shQuote(r_bin), " -lRlapack -lRblas")
  } else {
    blaslapack <- paste(cfg("LAPACK_LIBS"), cfg("BLAS_LIBS"))
  }
  base_libs <- paste(klu_lib, blaslapack)
  cppflags  <- paste0("-I", system.file("include", package = "cppDE"))

  ## Compiler settings read once, so parallel forks do not each call
  ## `R CMD config`.
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
      message(sprintf("using %-3s compiler: %s [%s]", label, strip(bin), trimws(shared)))
      return(invisible(NULL))
    }
    sets  <- setdiff(unique(extra), "")
    tags  <- c("every source", sprintf("%d of %d also", vapply(sets, function(k)
      sum(extra == k), integer(1)), length(ent)))
    width <- max(nchar(tags))
    message(sprintf("using %-3s compiler: %s", label, strip(bin)))
    message(sprintf("  %-*s : %s", width, tags[1], trimws(shared)))
    for (i in seq_along(sets))
      message(sprintf("  %-*s : %s", width, tags[i + 1L], sets[i]))
    invisible(NULL)
  }
  .reportToolchain("C",   cc_bin,  base,     "\\.c$")
  .reportToolchain("C++", cxx_bin, cxx_base, "\\.cpp$")

  ## unload stale DLLs
  loaded <- getLoadedDLLs()
  for (i in seq_along(roots))
    if (roots[i] %in% names(loaded)) try(dyn.unload(loaded[[roots[i]]][["path"]]), silent = TRUE)
  if (!is.null(outfile)) try(dyn.unload(outfile), silent = TRUE)

  # A backend's compileArgs can repeat a flag the base already sets (-fopenmp).
  # Compile flags only: repeated -l on a link line can be load-bearing.
  .mergeFlags <- function(...) {
    tok <- unlist(strsplit(trimws(paste(...)), "[[:space:]]+"))
    tok <- tok[nzchar(tok)]
    if (any(tok %in% c("-I", "-D", "-L", "-U", "-include", "-isystem")))
      return(paste(tok, collapse = " "))
    paste(unique(tok), collapse = " ")
  }

  ## Compile one file with its own flags set via PKG_*, so per-file linkArgs
  ## reach only the files that need them. Each mclapply fork has its own env.
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

      ## PKG_LIBS from the env can vanish from SHLIB's link command on some
      ## R/rtools combinations, hence the Makevars.win. Serial on Windows.
      mv_path <- file.path(dirname(entry$srcfile), "Makevars.win")
      mv_pre  <- if (file.exists(mv_path)) readLines(mv_path, warn = FALSE) else NULL
      writeLines(c(
        paste("PKG_CFLAGS =",   pkg_c),
        paste("PKG_CXXFLAGS =", pkg_cx),
        paste("PKG_CPPFLAGS =", cppflags),
        paste("PKG_LIBS =",     pkg_l)
      ), mv_path)
      on.exit({
        if (is.null(mv_pre)) try(unlink(mv_path), silent = TRUE)
        else                 try(writeLines(mv_pre, mv_path), silent = TRUE)
      }, add = TRUE)
    }
    ## system2() with pipes, not a `2>&1` token: without a shell on Windows the
    ## token reaches R CMD SHLIB as the make override PKG_LIBS=2>&1.
    Rexe <- file.path(R.home("bin"), "R")
    shlib_args <- c("CMD", "SHLIB", shQuote(entry$srcfile))
    if (verbose) cat(shQuote(Rexe), paste(shlib_args, collapse = " "), "\n")
    out <- suppressWarnings(system2(Rexe, shlib_args, stdout = TRUE, stderr = TRUE))
    if (verbose && length(out)) writeLines(out)
    status <- attr(out, "status")
    if (!is.null(status) && status != 0L)
      stop("Compilation failed: ", entry$srcfile, "\n",
           paste(utils::tail(out, 12), collapse = "\n"))
  }

  ## Command compiling a single source to a .o via a direct $CC/$CXX -c call.
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
    ## system2() with pipes, not a `2>&1` token, which on Windows reaches the
    ## compiler as an input file. The pipes keep the diagnostics of a failure.
    toks <- .tok(job$cmd)
    out  <- suppressWarnings(system2(toks[1], toks[-1], stdout = TRUE, stderr = TRUE))
    if (verbose && length(out)) writeLines(out)
    status <- attr(out, "status")
    if (!is.null(status) && status != 0L)
      stop("Compilation failed: ", job$srcfile, "\n",
           paste(utils::tail(out, 12), collapse = "\n"))
    job$srcfile
  }

  if (is.null(output)) {
    ## One shared object per source, all dyn.load()ed: refuse up front when R
    ## cannot load that many.
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
    ## Combined output: each source compiled to a .o, then one R CMD SHLIB link
    ## over the sources, so SHLIB picks the C++ linker. A precompiled header
    ## needs one flag set over enough C++ sources.
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
      message(sprintf("reusing %d unchanged object(s)", sum(fresh)))
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

    ## The link takes the union of every entry's linkArgs.
    all_link <- unique(unlist(lapply(info, function(e) strsplit(trimws(e$linkArgs %||% ""), "\\s+")[[1]])))
    all_link <- all_link[nzchar(all_link)]
    all_compile <- unique(unlist(lapply(info, function(e) strsplit(trimws(e$compileArgs %||% ""), "\\s+")[[1]])))
    all_compile <- all_compile[nzchar(all_compile)]

    output <- basename(sub(paste0("\\", so, "$"), "", output))

    ## Past the argument limit the objects go into a static archive in chunks,
    ## and SHLIB gets one anchor source plus the archive. The anchor is C++ when
    ## possible, as SHLIB picks the linker from the sources only.
    link_files <- files
    if (nchar(paste(shQuote(files), collapse = " ")) > .compileCmdLimit()) {
      anchor <- if (any(is_cxx)) which(is_cxx)[1] else 1L
      lib    <- file.path(dirname(files[1]), paste0(output, "_objects.a"))
      unlink(lib)
      on.exit(try(unlink(lib), silent = TRUE), add = TRUE)
      ## Via `R CMD`, which puts the Rtools toolchain on PATH under Windows.
      ar_bin     <- cfg("AR");     if (!nzchar(ar_bin))     ar_bin     <- "ar"
      ranlib_bin <- cfg("RANLIB"); if (!nzchar(ranlib_bin)) ranlib_bin <- "ranlib"
      ar_cmd <- paste(Rbin, "CMD", ar_bin, "qc", shQuote(lib))
      chunks <- .compileChunks(objs[-anchor], maxN = as.integer(chunkSize),
                               maxChars = max(1L, .compileCmdLimit() - nchar(ar_cmd)))
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
      message(sprintf("archived %d objects into %s (%d chunks)",
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

    ## PKG_LIBS from the env can vanish from SHLIB's link command on some
    ## R/rtools combinations, so a Makevars(.win) next to the sources repeats
    ## the flags. It is restored after the link.
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
    ## Link through system2() pipes, not a `2>&1` token, which on Windows
    ## becomes the make override PKG_LIBS=2>&1. The compiler banner is dropped,
    ## as only the link recipe runs.
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
    ## cOde loads a deSolve leaf's entry points from the shared object its
    ## `modelname` names, which must now be the combined one.
    for (o in fn_objs) .retargetCode(o, output)
    ## Only arguments passed as a plain variable can have their modelname
    ## updated in the caller; an expression has nothing to assign back to.
    for (i in which(is_dmod))
      if (make.names(obj.names[i]) == obj.names[i])
        eval.parent(parse(text = paste0("modelname(", obj.names[i], ") <- '", output, "'")))
  }

  invisible(TRUE)
}




#' Loaded Shared Objects in the Working Directory
#'
#' Only the working directory is searched, not `getOption("dMod.outdir")` or
#' the directories the objects were generated in.
#'
#' @return Character vector with the names of the loaded shared objects whose
#'   path lies in the working directory.
#' @seealso [loadDLL()]
#' @export
getLocalDLLs <- function() {
  
  all.dlls <- getLoadedDLLs()
  is.local <- sapply(all.dlls, function(x) grepl(getwd(), unclass(x)$path, fixed = TRUE))
  names(is.local)[is.local]
  
}




## A base name no loaded shared library uses: the desired one if it is
## free, otherwise the smallest `<name>_<i>`, i >= 2, with a warning. Must
## match cppDE's `unique_modelname()`.
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


## Loading an already loaded shared object is a no-op in R, so a rebuilt file is
## unloaded first. Entry points resolve by name, so this is safe for objects in use.
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

# Default directory of generated sources and shared objects
.dmodOutdir <- function() getOption("dMod.outdir", getwd())

## Directories to search for an object's shared libraries: where its sources
## were generated, the default output directory and the working directory.
.dllSearchDirs <- function(objects) {
  dirs <- unlist(lapply(objects, function(o)
    vapply(attr(o, "compileInfo") %||% list(),
           function(e) dirname(e$srcfile[1]), character(1))))
  unique(c(getwd(), .dmodOutdir(), dirs[nzchar(dirs)]))
}


#' Load Shared Objects for dMod Objects
#'
#' Loads the shared objects of dMod functions by their model names, e.g.
#' after restoring a workspace in a new R session. Searched are the working
#' directory, `getOption("dMod.outdir")` and the directories the sources were
#' generated in. Shared objects already loaded are skipped.
#'
#' @param ... objects of class `prdfn`, `obsfn`, `parfn` or `objfn`.
#'
#' @return Character vector of the files loaded by this call, invisibly.
#' @seealso [compile()], [modelname()]
#' @examples
#' \donttest{
#' g <- Y(c(y = "s * x"), parameters = "s", modelname = "loadDLL_ex",
#'        compile = TRUE, outdir = tempdir())
#' loadDLL(g)
#' }
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


## compileInfo plumbing ----------------------------------------------------------------

## ODE model class -------------------------------------------------------------------

## Dedup keys, cached on the list: merging is pairwise, so recomputing them
## every time would make summing a few thousand conditions quadratic.
.compileInfoKeys <- function(x) {
  k <- attr(x, "srckeys")
  if (!is.null(k) && length(k) == length(x)) return(k)
  vapply(x, function(e) if (length(e$srcfile)) paste(e$srcfile, collapse = "\x1f")
                        else NA_character_, character(1))
}

## Merge two compileInfo lists by srcfile, the first occurrence winning. NULL
## when both are empty, so objects without native code have no attribute.
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

## Build info of ODE model pieces from their `srcfile`, `compileArgs` and
## `linkArgs`, else from modelname-based file discovery in the working directory.
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
      ## Sparse-Jacobian flag, per file, as `compile()` sees only the dMod fn.
      sparse      = isTRUE(attr(o, "sparse"))
    )
  }
  out
}


## odemodel() constructor lives in R/odeClass.R.


