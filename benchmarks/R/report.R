## =====================================================================
##  report.R: the run directory, its README and the console summary.
## =====================================================================

bench_tag <- function(opts) {
  paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), "_", opts$tier,
         if (nzchar(opts$models)) paste0("_", gsub("[^A-Za-z0-9]+", "-", opts$models)),
         if (isTRUE(as.logical(opts$asa))) "_asa")
}

bench_info <- function(opts) {
  c(sprintf("date      %s", format(Sys.time())),
    sprintf("host      %s", Sys.info()[["nodename"]]),
    sprintf("dMod2     %s", format(utils::packageVersion("dMod2"))),
    sprintf("cppDE     %s", format(utils::packageVersion("cppDE"))),
    sprintf("R         %s", R.version.string),
    sprintf("options   %s", paste(names(opts), unlist(opts), sep = "=", collapse = " ")),
    sprintf("threads   OMP_NUM_THREADS=%s", Sys.getenv("OMP_NUM_THREADS", "unset")))
}

##  One row per model: the cost of each arm in value solves and its error.
bench_table <- function(df) {
  w <- stats::reshape(df[, c("model", "n_theta", "arm", "cost", "err_rel")],
                      idvar = c("model", "n_theta"), timevar = "arm", direction = "wide")
  w <- w[order(w$n_theta), , drop = FALSE]
  keep <- c("model", "n_theta", grep("^cost\\.", names(w), value = TRUE),
            grep("^err_rel\\.(forward|reverse|refine|asa)$", names(w), value = TRUE))
  w <- w[, keep[keep %in% names(w)], drop = FALSE]
  names(w) <- sub("^cost\\.", "", sub("^err_rel\\.", "err_", names(w)))
  w$value <- NULL
  rownames(w) <- NULL
  w
}

bench_print <- function(df) {
  if (!NROW(df)) return(invisible(NULL))
  cat("\nCost of a gradient in value solves, and its error against the forward\n",
      "gradient at atol 1e-12, rtol 1e-10 (max |g - g_ref| / max |g_ref|)\n\n", sep = "")
  tab <- bench_table(df)
  num <- vapply(tab, is.numeric, TRUE)
  tab[num] <- lapply(tab[num], signif, 3)
  print(tab, row.names = FALSE)
  v <- df[df$arm == "value", c("model", "sec", "solve_frac")]
  cat("\nValue solve [ms] and the share of the call in cppDE\n")
  for (i in seq_len(nrow(v)))
    cat(sprintf("  %-10s %9.3f  %3.0f%%\n", v$model[i], 1000 * v$sec[i], 100 * v$solve_frac[i]))
  invisible(tab)
}

##  A markdown table from a data frame of printable columns.
bench_md <- function(tab) {
  c(paste0("| ", paste(names(tab), collapse = " | "), " |"),
    paste0("|", strrep("---|", ncol(tab))),
    apply(tab, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")), "")
}

bench_readme <- function(outdir, opts, df, traits, skipped, figures) {
  lines <- c("# dMod2 benchmark run", "", "```", bench_info(opts), "```", "")
  if (length(traits)) {
    m <- do.call(rbind, traits)
    tt <- data.frame(model = rownames(m), ifelse(m, "x", ""), check.names = FALSE)
    lines <- c(lines, "## Traits", "", bench_md(tt))
  }
  if (NROW(df)) {
    tab <- bench_table(df)
    num <- vapply(tab, is.numeric, TRUE)
    tab[num] <- lapply(tab[num], function(x) format(signif(x, 3)))
    lines <- c(lines, "## Cost and error", "", bench_md(tab))
  }
  if (length(skipped))
    lines <- c(lines, "## Skipped", "", paste0("- ", names(skipped), ": ", skipped), "")
  if (length(figures))
    lines <- c(lines, "## Figures", "", paste0("- `", basename(figures), "`"), "")
  writeLines(lines, file.path(outdir, "README.md"))
}
