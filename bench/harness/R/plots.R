## =====================================================================
##  plots.R: the figures of a run, as PNG beside results.csv.
## =====================================================================

bench_plots <- function(df, outdir) {
  if (!NROW(df) || !requireNamespace("ggplot2", quietly = TRUE)) return(character(0))
  df$model <- factor(df$model, levels = unique(df$model[order(df$n_theta)]))
  grad <- df[df$arm != "value", , drop = FALSE]
  files <- c(file.path(outdir, "01-cost.png"), file.path(outdir, "02-error.png"))
  p1 <- ggplot2::ggplot(grad, ggplot2::aes(x = model, y = cost, fill = arm)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = NULL, y = "time / value solve",
                  title = "Cost of a gradient, models by number of parameters") +
    ggplot2::theme_bw()
  p2 <- ggplot2::ggplot(grad, ggplot2::aes(x = model, y = err_rel, colour = arm)) +
    ggplot2::geom_point(size = 2.5, position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = NULL, y = "max |g - g_ref| / max |g_ref|",
                  title = "Gradient error against forward at atol 1e-12, rtol 1e-10") +
    ggplot2::theme_bw()
  ggplot2::ggsave(files[1], p1, width = 8, height = 4.5, dpi = 120)
  ggplot2::ggsave(files[2], p2, width = 8, height = 4.5, dpi = 120)
  files
}
