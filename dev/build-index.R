# The subject index of the site: every span [term]{#def-...} of the articles, alphabetically, linked to
# the place where the text defines the term. An attribute index="..." in the span sets the entry text.
# Run from the package root: Rscript dev/build-index.R; writes vignettes/articles/index-terms.Rmd.
dir <- "vignettes/articles"
files <- setdiff(list.files(dir, "\\.Rmd$", full.names = TRUE), file.path(dir, "index-terms.Rmd"))
span <- "\\[((?:[^][]|\\[[^]]*\\])+)\\]\\{#(def-[A-Za-z0-9-]+)([^}]*)\\}"

entries <- do.call(rbind, lapply(files, function(f) {
  src <- readLines(f, warn = FALSE)
  title <- sub('^title:\\s*"?(.*?)"?\\s*$', "\\1", grep("^title:", src, value = TRUE)[1])
  hits <- regmatches(src, gregexpr(span, src, perl = TRUE))
  hits <- unlist(hits)
  if (!length(hits)) return(NULL)
  m <- regmatches(hits, regexec(span, hits, perl = TRUE))
  term <- vapply(m, `[`, "", 2)
  attrs <- vapply(m, `[`, "", 4)
  given <- sub('.*index="([^"]*)".*', "\\1", attrs)
  term <- ifelse(grepl('index="', attrs), given, term)
  data.frame(term = term, anchor = vapply(m, `[`, "", 3),
             page = sub("\\.Rmd$", ".html", basename(f)), title = title)
}))

key <- tolower(gsub("[^A-Za-z0-9 ]", "", gsub("\\$[^$]*\\$", "", entries$term)))
key[!nzchar(trimws(key))] <- tolower(entries$term[!nzchar(trimws(key))])
entries <- entries[order(key, entries$title), ]
key <- sort(key)
letter <- toupper(substr(trimws(key), 1, 1))
letter[!grepl("[A-Z]", letter)] <- "Symbols"

out <- c("---", 'title: "Index"', "---", "",
         "Terms of the guide and the methods articles, each linked to the place where the text defines it.", "")
for (l in unique(letter)) {
  e <- entries[letter == l, ]
  out <- c(out, paste("##", l), "",
           sprintf("- [%s](%s#%s), %s", e$term, e$page, e$anchor, e$title), "")
}
writeLines(out, file.path(dir, "index-terms.Rmd"))
cat(nrow(entries), "entries\n")
