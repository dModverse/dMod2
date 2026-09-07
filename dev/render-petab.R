#!/usr/bin/env Rscript
# Re-render dev/petab/PEtab.Rmd -> vignettes/PEtab.pdf
#
# PEtab.pdf is shipped as a *static* (pre-rendered) vignette via R.rsp::asis.
# Building it needs lualatex plus CMU / New Computer Modern fonts, and the
# worked example needs python-libsbml through reticulate; neither is available
# on a standard CI runner. The source bundle (Rmd + bibliography + CSL) lives
# under dev/petab/: that whole directory is .Rbuildignore'd, so the tarball
# carries only the rendered PDF and the .asis stub.
#
# Workflow: edit dev/petab/PEtab.Rmd, install the package, run this script,
# commit the regenerated vignettes/PEtab.pdf alongside the source change.
#
# Do not run devtools::build_vignettes(): copy_vignettes() moves the PDF into
# doc/ and deletes it from vignettes/, which for an .asis vignette destroys the
# only copy of the content. R CMD build and R CMD check leave it alone.
#
# Requirements (maintainer's machine only, never CI):
#   lualatex, fontspec, unicode-math, CMU Serif, NewCMMath, python-libsbml.

src <- "dev/petab/PEtab.Rmd"
out <- "vignettes/PEtab.pdf"

if (!file.exists(src))
  stop("Cannot find ", src, ": run this from the package root.")

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
  output_file = "PEtab.pdf",
  output_dir  = normalizePath("vignettes", mustWork = TRUE),
  envir       = new.env(parent = globalenv()),
  quiet       = FALSE
)

if (!file.exists(out))
  stop("Render finished but ", out, " was not produced.")

message("Wrote ", out, " (", format(file.size(out) / 1024, digits = 1), " KiB)")
