# Shared sympy/python helpers for the symmetry test files
# (test-symmetryDetection.R, test-symmetryReduction.R). helper-*.R is sourced by
# testthat before every test file.

.sympy_works <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(FALSE)
  isTRUE(tryCatch({
    sympy <- reticulate::import("sympy", convert = TRUE)
    nzchar(sympy[["__version__"]])
  }, error = function(e) FALSE))
}


.sd_module <- function() {
  code_dir <- system.file("code", package = "dMod2")
  sysmod <- reticulate::import("sys", convert = TRUE)
  if (!(code_dir %in% sysmod$path)) sysmod$path <- c(code_dir, sysmod$path)
  reticulate::import("symmetryDetection", convert = TRUE)
}


# exact symbolic equality of two expression strings via sympy
.symExprEqual <- function(a, b) {
  spy <- reticulate::import("sympy", convert = TRUE)
  # arguments may come from $generator, which is in R's power syntax
  ab <- gsub("\\^", "**", paste0("(", a, ") - (", b, ")"))
  as.character(spy$simplify(spy$sympify(ab))) == "0"
}


# numeric tangent of a reported direction at point `pt` over `coords`. $generator
# always holds the tangent components xi_i directly (a scaling's integer weight
# w_i is expanded to xi_i = w_i * z_i at the finalisation boundary), so evaluating
# each component at the point gives the tangent.
.symTangent <- function(d, pt, coords) {
  v <- setNames(numeric(length(coords)), coords)
  for (nm in names(d$generator))
    v[nm] <- eval(parse(text = d$generator[[nm]]), pt)
  v
}


