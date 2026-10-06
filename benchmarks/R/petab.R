## =====================================================================
##  petab.R: find and import the PEtab problems of the collection.
## =====================================================================

##  The collection has two layouts: problems/<id>/v1/problem.yaml since 2025,
##  Benchmark-Models/<id>/<id>.yaml before. Both are read, and so is a flat
##  directory of problem folders.

petab_index <- function(root) {
  root <- normalizePath(root, mustWork = TRUE)
  base <- c(file.path(root, "problems"), file.path(root, "Benchmark-Models"), root)
  base <- base[dir.exists(base)][1L]
  dirs <- list.dirs(base, recursive = FALSE)
  yaml <- vapply(dirs, function(d) {
    for (sub in c(file.path(d, "v1"), d)) {
      y <- list.files(sub, pattern = "\\.ya?ml$", full.names = TRUE)
      own <- y[basename(y) %in% c("problem.yaml", paste0(basename(d), ".yaml"))]
      if (length(own)) return(own[1L])
      if (length(y)) return(y[1L])
    }
    NA_character_
  }, "")
  out <- data.frame(name = basename(dirs), yaml = unname(yaml), stringsAsFactors = FALSE)
  out[!is.na(out$yaml), , drop = FALSE]
}

##  The tier lists short names; a problem matches when its id starts with one.
petab_select <- function(index, names) {
  hit <- vapply(names, function(n) {
    i <- which(startsWith(index$name, paste0(n, "_")) | index$name == n)
    if (length(i)) i[1L] else NA_integer_
  }, 0L)
  if (anyNA(hit)) warning("not in the collection: ", paste(names[is.na(hit)], collapse = ", "),
                          call. = FALSE)
  out <- index[hit[!is.na(hit)], , drop = FALSE]
  out$short <- names[!is.na(hit)]
  out
}

##  One import, compiled into its own directory. `backend = "Sundials"` builds
##  the CVODES adjoint for the --asa arm.
petab_import <- function(entry, outdir, backend = "cppDE", cores = 1L, tol = 1e-6,
                         sparse = NULL) {
  od <- file.path(outdir, paste0(entry$short, "_", backend))
  dir.create(od, recursive = TRUE, showWarnings = FALSE)
  o <- list(atol = tol, rtol = tol)
  t0 <- proc.time()[["elapsed"]]
  pet <- dMod2::importPEtab(entry$yaml, backend = backend, cores = cores,
                            modelname = paste0("bm_", gsub("[^A-Za-z0-9]", "", entry$short),
                                               if (backend == "Sundials") "_asa" else ""),
                            derivMode = c("forward", "reverse"), optionsOde = o,
                            optionsSens = o, sparse = sparse, outdir = od)
  attr(pet, "compile_s") <- proc.time()[["elapsed"]] - t0
  pet
}

##  Solver tolerances of every Xs() in the problem, as one setting.
petab_tolerance <- function(pet, atol, rtol) {
  o <- list(atol = atol, rtol = rtol)
  dMod2::controls(pet$x, name = "optionsOde") <- o
  dMod2::controls(pet$x, name = "optionsSens") <- o
  invisible(pet)
}
