#!/usr/bin/env Rscript
# Re-render dev/symmetries/Symmetries.Rmd -> vignettes/Symmetries.pdf
#
# Symmetries.pdf is shipped as a *static* (pre-rendered) vignette via R.rsp::asis.
# Building it needs lualatex plus CMU / New Computer Modern fonts and a Unicode-complete
# monospace font for the printed generators, and every computational chunk drives sympy
# through reticulate; neither is available on a standard CI runner. The source bundle
# (Rmd + bibliography + CSL) lives under dev/symmetries/: that whole directory is
# .Rbuildignore'd, so the tarball carries only the rendered PDF and the .asis stub.
#
# Workflow: edit dev/symmetries/Symmetries.Rmd, install the package, run this script,
# commit the regenerated vignettes/Symmetries.pdf alongside the source change. The
# chunks call library(dMod2), so the PDF shows whatever version is installed.
#
# Do not run devtools::build_vignettes(): copy_vignettes() moves the PDF into doc/
# and deletes it from vignettes/, which for an .asis vignette destroys the only copy
# of the content. R CMD build and R CMD check leave it alone.
#
# Requirements (maintainer's machine only, never CI):
#   lualatex, fontspec, unicode-math, CMU Serif, NewCMMath, DejaVu Sans Mono,
#   python sympy/numpy/scipy through reticulate.

src <- "dev/symmetries/Symmetries.Rmd"
out <- "vignettes/Symmetries.pdf"

if (!file.exists(src))
  stop("Cannot find ", src, ": run this from the package root.")

# Without sympy every computational chunk is skipped (eval = haveSympy) and the PDF
# would ship with empty outputs, which is worse than not rendering at all.
if (!isTRUE(tryCatch(reticulate::py_module_available("sympy"),
                     error = function(e) FALSE)))
  stop("sympy is not available through reticulate: the chunks would all be skipped.")

# On Windows without pandoc on PATH: try the usual RStudio bundle locations.
if (!nzchar(Sys.which("pandoc")) && !nzchar(Sys.getenv("RSTUDIO_PANDOC"))) {
  candidates <- c(
    "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools",
    "C:/Program Files/RStudio/bin/pandoc",
    "C:/Program Files/Quarto/bin/tools"
  )
  hit <- candidates[file.exists(file.path(candidates, "pandoc.exe"))][1]
  if (!is.na(hit)) {
    Sys.setenv(RSTUDIO_PANDOC = hit)
    message("Setting RSTUDIO_PANDOC = ", hit)
  }
}

# On Windows without lualatex on PATH: try the standard TeX Live install paths.
if (!nzchar(Sys.which("lualatex"))) {
  candidates <- Sys.glob(c("C:/texlive/*/bin/windows", "C:/texlive/*/bin/win32"))
  hit <- candidates[file.exists(file.path(candidates, "lualatex.exe"))][1]
  if (!is.na(hit)) {
    Sys.setenv(PATH = paste(hit, Sys.getenv("PATH"), sep = .Platform$path.sep))
    message("Prepended TeX Live to PATH: ", hit)
  }
}

# Own environment: chunks assign into `envir`, and a name collision with this
# script's own variables would otherwise go unnoticed.
rmarkdown::render(
  input       = src,
  output_file = "Symmetries.pdf",
  output_dir  = normalizePath("vignettes", mustWork = TRUE),
  envir       = new.env(parent = globalenv()),
  quiet       = FALSE
)

if (!file.exists(out))
  stop("Render finished but ", out, " was not produced.")

message("Wrote ", out, " (", format(file.size(out) / 1024, digits = 1), " KiB)")
