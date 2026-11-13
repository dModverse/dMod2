# The result saved in `file` by `script`. A missing file is computed first, by sourcing the script
# in its own directory into `envir`. The attribute "stamp" of a result holds the md5 of the files
# it was computed from, relative to the directory of the script; a file that changed since warns.
stored <- function(file, script, envir = parent.frame()) {
  if (!file.exists(file)) source(script, local = envir, chdir = TRUE)
  x <- readRDS(file)
  md5 <- attr(x, "stamp")$md5
  if (is.null(md5)) {
    warning(sprintf("%s has no stamp", basename(file)), call. = FALSE)
    return(x)
  }
  now <- tools::md5sum(file.path(dirname(script), names(md5)))
  changed <- names(md5)[is.na(now) | now != md5]
  if (length(changed))
    warning(sprintf("%s was computed from an earlier version of %s", basename(file),
                    paste(changed, collapse = ", ")), call. = FALSE)
  x
}
