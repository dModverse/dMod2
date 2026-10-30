#!/usr/bin/env Rscript
# Re-render dev/nlme/NLME.Rmd -> vignettes/NLME.pdf
#
# NLME.pdf is shipped as a static vignette via R.rsp::asis; dev/ is .Rbuildignore'd, so the
# tarball holds only the PDF and the .asis stub. Its chunks are not evaluated. Building it needs
# pandoc and xelatex. Do not run devtools::build_vignettes(), which moves the PDF out of vignettes/.

src <- "dev/nlme/NLME.Rmd"
out <- "vignettes/NLME.pdf"
if (!file.exists(src)) stop("Cannot find ", src, ": run this from the package root.")
rmarkdown::render(input = src, output_file = "NLME.pdf",
                  output_dir = normalizePath("vignettes", mustWork = TRUE),
                  envir = new.env(parent = globalenv()), quiet = FALSE)
if (!file.exists(out)) stop("Render finished but ", out, " was not produced.")
message("Wrote ", out, " (", format(file.size(out) / 1024, digits = 1), " KiB)")
