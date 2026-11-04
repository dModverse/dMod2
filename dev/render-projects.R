#!/usr/bin/env Rscript
# Knit examples/<name>/README.Rmd to README.md with its figures in examples/<name>/figures/, and
# write vignettes/articles/<name>.Rmd, the page of the project on the site. Arguments name the
# projects, by default all of them. Run from the package root, against the installed package.
#
# A README reads what the scripts of its project saved in results/ and fits nothing; run the
# scripts first, in the order the README gives.

if (!dir.exists("examples")) stop("Run this from the package root.")

projects <- basename(list.dirs("examples", recursive = FALSE))
projects <- projects[file.exists(file.path("examples", projects, "README.Rmd"))]
named <- commandArgs(TRUE)
if (length(named)) projects <- intersect(projects, named)

# The article of a project: its README without the first heading, which becomes the title, and
# the figures copied beside the article under flat names (pkgdown escapes nested image paths).
writeArticle <- function(p) {
  readme <- readLines(file.path("examples", p, "README.md"))
  title <- sub("^#\\s+", "", readme[grep("^# ", readme)[1]])
  writeLines(c("---", sprintf("title: \"%s\"", gsub("\"", "\\\\\"", title)), "---", "",
               "```{r, echo = FALSE, results = \"asis\"}",
               sprintf("dir <- \"../../examples/%s\"", p),
               "page <- readLines(file.path(dir, \"README.md\"))",
               "page <- page[-grep(\"^# \", page)[1]]",
               "figures <- unique(unlist(regmatches(page, gregexpr(\"figures/[^\\\")[:space:]]+\", page))))",
               sprintf("flat <- paste0(\"%s-\", basename(figures))", p),
               "file.copy(file.path(dir, figures), flat, overwrite = TRUE)",
               "for (i in seq_along(figures)) page <- gsub(figures[i], flat[i], page, fixed = TRUE)",
               "page <- gsub(\"<img src=\", \"<img alt=\\\"Output of the code above\\\" src=\", page, fixed = TRUE)",
               "cat(page, sep = \"\\n\")",
               "```"),
             file.path("vignettes", "articles", paste0(p, ".Rmd")))
}

# An R script as R Markdown, for README.Rmd to show a whole script with
# `cat(knitr::knit_child(text = script("analysis.R"), quiet = TRUE), sep = "\n")`: a top-level
# comment block after a blank line that opens with a capital becomes prose, `# Title ----` a
# heading, `## --- Title ---` a subheading; other comments stay in the code. The comment block
# that opens the script is left out.
script <- function(file) {
  src <- readLines(file)
  # the comment block that opens the script, which the README says in its own words
  head <- cumprod(grepl("^#", src)) == 1
  src <- src[!head]
  src <- sub("^# (.*\\S) -{3,}$", "#' ## \\1", src, perl = TRUE)
  src <- sub("^## --- (.*) ---$", "#' ### \\1", src, perl = TRUE)
  prose <- FALSE
  for (i in seq_along(src)) {
    above <- if (i == 1) "" else src[i - 1]
    opens <- !nzchar(trimws(above)) && grepl("^# [A-Z`]", src[i])
    if (grepl("^#( |$)", src[i]) && (opens || prose)) src[i] <- sub("^#( |$)", "#'\\1", src[i])
    prose <- grepl("^#'", src[i])
  }
  knitr::spin(text = src, knit = FALSE, format = "Rmd")
}

render <- function(project) {
  dir <- normalizePath(file.path("examples", project))
  # the figures knitr made last time; other images in figures/ stay
  unlink(list.files(file.path(dir, "figures"), pattern = "-[0-9]+\\.png$", full.names = TRUE))
  knitr::opts_knit$set(root.dir = dir)
  # twice the pixels at the displayed width, shown without a caption
  knitr::opts_chunk$set(fig.path = "figures/", fig.width = 8, fig.height = 5, dev = "png", dpi = 96,
                        fig.retina = 2, comment = "#>", message = FALSE, warning = FALSE, error = FALSE)
  knitr::render_markdown()
  knitr::knit_hooks$set(plot = function(x, options)
    sprintf('\n<img src="%s" width="%s" />\n\n', x, options$out.width))
  owd <- setwd(dir)
  on.exit(setwd(owd))
  env <- new.env()
  env$script <- script
  knitr::knit("README.Rmd", output = "README.md", envir = env, quiet = TRUE)
}

for (project in projects) {
  render(project)
  writeArticle(project)
  cat("rendered", project, "\n")
}
