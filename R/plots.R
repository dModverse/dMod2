
# ggplot2 / dplyr NSE column references; declared so R CMD check does not
# flag them as undefined globals.
utils::globalVariables(c("predicted", "observed", "sd_est", "iter", "level"))


# PDF viewer, theme and color palettes ---

#' Open a Plot in an External PDF Viewer
#'
#' Saves a plot as PDF in [tempdir()] and opens it with an external command.
#'
#' @param plot A `ggplot2` plot object. Defaults to [ggplot2::last_plot()].
#' @param command Character, the shell command that opens the PDF file.
#'   Default `NULL` opens it with the viewer of the platform: `"open"` on
#'   macOS, `"xdg-open"` on Linux and `shell.exec()` on Windows.
#' @param ... Arguments passed to [ggplot2::ggsave()].
#' @return The path of the PDF file, invisibly.
#' @export
ggopen <- function(plot = last_plot(), command = NULL, ...) {
  filename <- tempfile(pattern = "Rplot", fileext = ".pdf")
  ggsave(filename = filename, plot = plot, ...)
  if (is.null(command) && .Platform$OS.type == "windows") {
    get("shell.exec", envir = baseenv())(filename)
    return(invisible(filename))
  }
  if (is.null(command))
    command <- if (Sys.info()[["sysname"]] == "Darwin") "open" else "xdg-open"
  system(command = paste(command, shQuote(filename)))
  invisible(filename)
}


#' Standard Plotting Theme of dMod
#'
#' @param base_size Numeric, base font size. Defaults to 12.
#' @param base_family Character, font family. Defaults to `""`.
#' @param showGrid Logical, keep the panel grid. `FALSE` (default) drops it;
#'   `TRUE` keeps the grid of [ggplot2::theme_bw()].
#' @return A `ggplot2` theme.
#' @seealso [scale_color_dMod()]
#' @export
#' @examples
#' library(ggplot2)
#' ggplot(data.frame(x = 1:10, y = (1:10)^2), aes(x, y)) +
#'   geom_line() + theme_dMod()
theme_dMod <- function(base_size = 12, base_family = "", showGrid = FALSE) {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(line = element_line(colour = "black"),
          rect = element_rect(fill = "white", colour = NA),
          text = element_text(colour = "black"),
          axis.text = element_text(size = rel(1.0), colour = "black"),
          axis.text.x = element_text(margin = margin(t = 4, r = 4, b = 0, l = 4, unit = "mm")),
          axis.text.y = element_text(margin = margin(t = 4, r = 4, b = 4, l = 0, unit = "mm")),
          axis.ticks = element_line(colour = "black"),
          axis.ticks.length = unit(-2, "mm"),
          legend.key = element_rect(colour = NA),
          panel.border = element_rect(colour = "black"),
          strip.background = element_rect(fill = "white", colour = NA),
          strip.text = element_text(size = rel(1.0))) +
    (if (showGrid) NULL else theme(panel.grid = element_blank()))

}

# ---- palettes --------------------------------------------------------------
# dE: smallest pairwise CIE2000 distance over normal and simulated dichromatic
# vision (colorspace, farver::compare_colour()); 10 or above is colorblind-safe.

#' Seed Colors of the dMod Palette
#'
#' The ten qualitative house colors. They are not colorblind-safe, see
#' [dMod_palettes()] for alternatives.
#'
#' @format A character vector of hex color codes.
#' @seealso [dMod_colors_cb], [dMod_palette()]
#' @export
dMod_colors <- c("#000000", "#C5000B", "#0084D1", "#579D1C", "#FF950E",
                 "#4B1F6F", "#CC79A7", "#006400", "#F0E442", "#8B4513")

#' Colorblind-Safe Seed Colors
#'
#' Okabe-Ito palette, extended to twenty colorblind-safe colors. The first
#' eight are the Okabe-Ito colors.
#'
#' @format A character vector of 20 hex color codes.
#' @seealso [dMod_colors], [dMod_palette()]
#' @export
dMod_colors_cb <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#4F3408", "#072F62",
  "#645B68", "#4E644A", "#A6D38F", "#2600DA", "#3F3644",
  "#02ECC9", "#7B75FD", "#954A0E", "#323C2D", "#6D7F80")

# Qualitative palettes. Five colorblind-safe, five not; `dMod_palette()` and
# `scale_color_dMod()` pick from here by name.
.dMod_qual <- list(
  okabe    = list(dE = 11.1, colors = dMod_colors_cb),
  muted    = list(dE = 11.6, colors = c(
    "#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677",
    "#882255", "#AA4499", "#2D330A", "#9B7BFB", "#FC8FB1", "#1E352C", "#A884B7")),
  medium   = list(dE = 12.1, colors = c(
    "#6699CC", "#004488", "#EECC66", "#994455", "#997700", "#EE99AA", "#2D330A",
    "#541C39", "#5EE3F9", "#807987", "#6761DF", "#A79C65")),
  dark     = list(dE = 11.9, colors = c(
    "#000000", "#004488", "#994455", "#997700", "#8770FD", "#3F2D11", "#8B8491",
    "#2B303D", "#616689", "#525901", "#758E65", "#6D38DB")),
  contrast = list(dE = 21.4, colors = c("#004488", "#DDAA33", "#BB5566")),
  dMod     = list(dE =  2.1, colors = dMod_colors),
  dark2    = list(dE =  2.3, colors = c(
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D",
    "#666666")),
  set1     = list(dE =  3.2, colors = c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628",
    "#F781BF", "#999999")),
  set2     = list(dE =  1.6, colors = c(
    "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
    "#B3B3B3")),
  paired   = list(dE =  1.3, colors = c(
    "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F",
    "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")))

#' Continuous dMod Color Ramps
#'
#' `dMod_gradient` is the sequential house ramp, `dMod_divergent` the diverging
#' one; `dMod_gradient_cb` and `dMod_divergent_cb` are their colorblind-safe
#' counterparts and the defaults of [scale_color_dMod_c()] and
#' [scale_color_dMod_div()]. See [dMod_palettes()] for the full set.
#'
#' @format Character vectors of hex color codes.
#' @export
dMod_gradient <- c("#4B1F6F", "#0084D1", "#579D1C", "#FF950E", "#F0E442")

#' @export
#' @rdname dMod_gradient
dMod_gradient_cb <- c("#000000", "#4B1F6F", "#0072B2", "#009E73",
                      "#E69F00", "#F0E442")

#' @export
#' @rdname dMod_gradient
dMod_divergent <- c("#0084D1", "#FFFFFF", "#C5000B")

#' @export
#' @rdname dMod_gradient
dMod_divergent_cb <- c("#072F62", "#56B4E9", "#FFFFFF", "#E69F00", "#4F3408")

.dMod_seq <- list(
  okabe   = list(dE = 12.0, colors = dMod_gradient_cb),
  viridis = list(dE = 13.8, colors = c("#440154","#3B528B","#21908C","#5DC863","#FDE725")),
  grey    = list(dE = 13.2, colors = c("#111111", "#DDDDDD")),
  heat    = list(dE = 10.7, colors = c("#FFFFFF", "#FF950E", "#C5000B")),
  dMod    = list(dE =  5.4, colors = dMod_gradient))

.dMod_div <- list(
  okabe = list(dE = 25.6, colors = dMod_divergent_cb),
  puor  = list(dE = 15.7, colors = c("#5E3C99","#B2ABD2","#F7F7F7","#FDB863","#E66101")),
  rdbu  = list(dE = 15.2, colors = c("#CA0020","#F4A582","#F7F7F7","#92C5DE","#0571B0")),
  dMod  = list(dE = 12.5, colors = dMod_divergent),
  brbg  = list(dE =  7.9, colors = c("#A6611A","#DFC27D","#F5F5F5","#80CDC1","#018571")))

.dMod_registry <- list(qualitative = .dMod_qual, sequential = .dMod_seq,
                       diverging = .dMod_div)

.dMod_pick <- function(palette, type) {
  reg <- .dMod_registry[[type]]
  if (!is.character(palette) || length(palette) != 1L || !palette %in% names(reg))
    stop("dMod: no ", type, " palette \"", paste(palette, collapse = ", "),
         "\". Available: ", paste(names(reg), collapse = ", "),
         ". See dMod_palettes().", call. = FALSE)
  reg[[palette]]
}

#' Overview of the dMod Palettes
#'
#' One row per palette. Pass any of the names as `palette` to [dMod_palette()],
#' [scale_color_dMod()], [scale_color_dMod_c()] or [scale_color_dMod_div()].
#'
#' @param type `NULL` (default) for all palettes, or any of `"qualitative"`,
#'   `"sequential"` and `"diverging"`.
#' @return A data frame with columns `name`, `type`, `n` (number of seed
#'   colors), `dE` (smallest pairwise CIE2000 color distance under normal and
#'   simulated dichromatic vision) and `colorblind` (`dE >= 10`).
#' @export
#' @examples
#' dMod_palettes()
#' dMod_palettes("qualitative")
dMod_palettes <- function(type = NULL) {
  types <- if (is.null(type)) names(.dMod_registry) else match.arg(
    type, names(.dMod_registry), several.ok = TRUE)
  out <- do.call(rbind, lapply(types, function(ty) {
    reg <- .dMod_registry[[ty]]
    data.frame(name = names(reg), type = ty,
               n = vapply(reg, function(p) length(p$colors), integer(1)),
               colorblind = vapply(reg, function(p) p$dE >= 10, logical(1)),
               dE = vapply(reg, function(p) p$dE, numeric(1)),
               row.names = NULL, stringsAsFactors = FALSE)
  }))
  out[order(out$type, !out$colorblind, -out$dE), ]
}

#' Generate Distinct Colors from a dMod Palette
#'
#' Returns the first `n` colors of the chosen palette. Beyond its seeds the
#' colors come from `Polychrome::createPalette()` if Polychrome is installed and
#' from [grDevices::hcl.colors()] otherwise. These extra colors are not
#' colorblind-safe; a colorblind-safe palette warns when it runs out of seeds.
#'
#' @param n Integer, number of colors.
#' @param palette Name of a qualitative palette, see [dMod_palettes()].
#'   Defaults to `"okabe"`.
#' @return Character vector of length `n` with hex color codes.
#' @export
#' @examples
#' dMod_palette(3)
#' dMod_palette(3, palette = "muted")
dMod_palette <- function(n, palette = "okabe") {
  n <- as.integer(n)
  if (n <= 0L) return(character(0))
  pal <- .dMod_pick(palette, "qualitative")
  seeds <- pal$colors
  if (n <= length(seeds)) return(unname(seeds[seq_len(n)]))
  if (pal$dE >= 10)
    warning("dMod_palette(): palette \"", palette, "\" holds ", length(seeds),
            " colorblind-safe colors, ", n, " were requested. The extra ",
            n - length(seeds), " are not: encode them with linetype, shape ",
            "or facets instead.", call. = FALSE)
  if (requireNamespace("Polychrome", quietly = TRUE)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
          rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(123L)
    # createPalette() shifts the seeds, so keep them verbatim and borrow only the
    # tail. `range` caps the luminance so tail colors stay visible as lines on a
    # white panel.
    extended <- unname(Polychrome::createPalette(n, seedcolors = seeds,
                                                 range = c(25, 70)))
    return(c(seeds, extended[(length(seeds) + 1L):n]))
  }
  c(seeds, grDevices::hcl.colors(n - length(seeds), "Dark 3"))
}

#' Discrete dMod Color Scales
#'
#' @param palette Name of a qualitative palette, see [dMod_palettes()].
#'   Defaults to the colorblind-safe `"okabe"`.
#' @param ... Arguments passed to [ggplot2::discrete_scale()].
#' @return A discrete `ggplot2` scale.
#' @seealso [scale_color_dMod_c()], [dMod_palette()]
#' @export
#' @examples
#' library(ggplot2)
#' times <- seq(0, 2*pi, 0.1)
#' values <- sin(times)
#' data <- data.frame(
#'    time = times,
#'    value = c(values, 1.2*values, 1.4*values, 1.6*values),
#'    group = rep(c("C1", "C2", "C3", "C4"), each = length(times))
#' )
#' ggplot(data, aes(time, value, colour = group)) + geom_line() +
#'    theme_dMod() + scale_color_dMod()
scale_color_dMod <- function(..., palette = "okabe") {
  ggplot2::discrete_scale(aesthetics = "colour",
                          palette = function(n) dMod_palette(n, palette), ...)
}

#' @export
#' @rdname scale_color_dMod
scale_fill_dMod <- function(..., palette = "okabe") {
  ggplot2::discrete_scale(aesthetics = "fill",
                          palette = function(n) dMod_palette(n, palette), ...)
}

.dMod_ramp <- function(cols, direction) if (direction < 0) rev(cols) else cols

.dMod_mid_rescaler <- function(mid) {
  function(x, to = c(0, 1), from = range(x, na.rm = TRUE))
    scales::rescale_mid(x, to, from, mid)
}

#' Continuous dMod Color Scales
#'
#' Sequential (`_c`) and diverging (`_div`) counterparts of the discrete
#' [scale_color_dMod()]. The diverging scales map `mid` to white, also for
#' asymmetric limits.
#'
#' @param mid Numeric, the value mapped to white. Defaults to 0.
#' @param direction 1 (default) or -1; -1 reverses the ramp.
#' @param palette Name of a sequential palette for `_c`, a diverging one for
#'   `_div`, see [dMod_palettes()]. Defaults to the colorblind-safe `"okabe"`.
#' @param ... Arguments passed to [ggplot2::scale_colour_gradientn()] or
#'   [ggplot2::scale_fill_gradientn()].
#' @return A continuous `ggplot2` scale.
#' @seealso [scale_color_dMod()]
#' @export
#' @examples
#' library(ggplot2)
#' d <- expand.grid(time = seq(0, 10, 0.1), eps = seq(0, 1, 0.1))
#' ggplot(d, aes(time, exp(-eps * time), group = eps, colour = eps)) +
#'    geom_line() + theme_dMod() + scale_color_dMod_c()
scale_color_dMod_c <- function(..., direction = 1, palette = "okabe") {
  ggplot2::scale_colour_gradientn(
    colours = .dMod_ramp(.dMod_pick(palette, "sequential")$colors, direction), ...)
}

#' @export
#' @rdname scale_color_dMod_c
scale_fill_dMod_c <- function(..., direction = 1, palette = "okabe") {
  ggplot2::scale_fill_gradientn(
    colours = .dMod_ramp(.dMod_pick(palette, "sequential")$colors, direction), ...)
}

#' @export
#' @rdname scale_color_dMod_c
scale_color_dMod_div <- function(mid = 0, ..., direction = 1, palette = "okabe") {
  ggplot2::scale_colour_gradientn(
    colours = .dMod_ramp(.dMod_pick(palette, "diverging")$colors, direction),
    rescaler = .dMod_mid_rescaler(mid), ...)
}

#' @export
#' @rdname scale_color_dMod_c
scale_fill_dMod_div <- function(mid = 0, ..., direction = 1, palette = "okabe") {
  ggplot2::scale_fill_gradientn(
    colours = .dMod_ramp(.dMod_pick(palette, "diverging")$colors, direction),
    rescaler = .dMod_mid_rescaler(mid), ...)
}


ggplot <- function(...) ggplot2::ggplot(...) + scale_color_dMod() + theme_dMod()


# Transformation of plotted values ---------------------------------------

#' Coordinate Transformation for Data Frames
#'
#' Applies a symbolically defined transformation to the `value` column of a
#' data frame. A `sigma` column is transformed by first-order error
#' propagation.
#' @param data Data frame with columns `name` (character) and `value`
#'   (numeric), optionally `sigma` (numeric).
#' @param transformations Character, one transformation for all names, or a
#'   named list of characters, one per entry of the `name` column. Each
#'   transformation is an expression in a single symbol, e.g. `"log(x)"`.
#' @return The data frame with transformed `value` and `sigma`.
#' @export
#' 
#' @examples
#' d <- data.frame(name = rep(c("u", "v"), each = 3), time = 0:2, value = 1:6,
#'                 sigma = 0.2)
#' coordTransform(d, "log(y)")
#' coordTransform(d, list(u = "y^2", v = "exp(-y)"))
coordTransform <- function(data, transformations) {

  mynames <- unique(as.character(data$name))
  if (!is.list(transformations))
    transformations <- as.list(structure(rep(transformations, length(mynames)), names = mynames))

  # in place, so the rows keep their order
  for (n in intersect(mynames, names(transformations))) {
    rows <- which(data$name == n)
    subdata <- data[rows, , drop = FALSE]
    mytrafo <- str2expression(replaceSymbols(getSymbols(transformations[[n]])[1], "value",
                                             transformations[[n]]))
    # Errors propagate linearly; without one the transformation need not be differentiable.
    if (any(!is.na(subdata$sigma)))
      data$sigma[rows] <- abs(eval(D(mytrafo, "value"), subdata)) * subdata$sigma
    data$value[rows] <- eval(mytrafo, subdata)
  }
  data
}


# Generics of the prediction and data plots ---------------------------------


#' Plot a List of Model Predictions
#'
#' @param prediction A [prdlist], the output of a prediction function such as
#'   those generated by [Xs()].
#' @param ... Logical expressions to subset the plotted data.
#' @param scales The `scales` argument of [ggplot2::facet_wrap()] or
#'   [ggplot2::facet_grid()]: `"free"` (default), `"fixed"`, `"free_x"` or
#'   `"free_y"`.
#' @param facet `"wrap"` (default) or `"grid"`.
#' @param transform Named list of transformations of the states, see
#'   [coordTransform()].
#' @details The plotted data frame has columns `time`, `value`, `name` and
#'   `condition`.
#'
#' @return A `ggplot` object.
#' @seealso [plotCombined()], [plotData()]
#' @import ggplot2
#' @example inst/examples/plotting.R
#' @export
plotPrediction <- function(prediction,...) {
  UseMethod("plotPrediction", prediction)
}


#' Plot Model Predictions and Data in One Plot
#'
#' @param prediction A [prdlist], the output of a prediction function such as
#'   those generated by [Xs()].
#' @param data A [datalist], i.e. a named list of data frames with columns
#'   `name`, `time`, `value` and `sigma`.
#' @param ... Logical expressions to subset the plotted data.
#' @param scales The `scales` argument of [ggplot2::facet_wrap()] or
#'   [ggplot2::facet_grid()]: `"free"` (default), `"fixed"`, `"free_x"` or
#'   `"free_y"`.
#' @param facet `"wrap"` (default), `"grid"` or `"wrap_plain"`. `"wrap_plain"`
#'   gives one panel per combination of name and condition.
#' @param transform Named list of transformations of the states, see
#'   [coordTransform()].
#' @param aesthetics Named list of aesthetic mappings given as character, e.g.
#'   `list(linetype = "name")`. They can refer to columns of the condition grid.
#' @details The plotted data frame has columns `time`, `value`, `sigma`, `name`
#'   and `condition`.
#'
#' @return A `ggplot` object.
#' @seealso [plotPrediction()], [plotData()]
#' @example inst/examples/plotting.R
#' @export
plotCombined <- function(prediction,...) {
  UseMethod("plotCombined", prediction)
}


#' Plot a List of Data Points
#'
#' @param data A [datalist], i.e. a named list of data frames with columns
#'   `name`, `time`, `value` and `sigma`, or such a data frame with a
#'   `condition` column.
#' @param ... Logical expressions to subset the plotted data.
#' @param scales The `scales` argument of [ggplot2::facet_wrap()] or
#'   [ggplot2::facet_grid()]: `"free"` (default), `"fixed"`, `"free_x"` or
#'   `"free_y"`.
#' @param facet `"wrap"` (default), `"grid"` or `"wrap_plain"`. `"wrap_plain"`
#'   gives one panel per combination of name and condition.
#' @param transform Named list of transformations of the states, see
#'   [coordTransform()].
#' @details The plotted data frame has columns `time`, `value`, `sigma`, `name`
#'   and `condition`.
#'
#' @return A `ggplot` object.
#' @seealso [plotPrediction()], [plotCombined()]
#' @example inst/examples/plotting.R
#' @export
plotData  <- function(data,...) {
  UseMethod("plotData", data)
}

#' @export
#' @rdname plotData
plotData.data.frame <- function(data, ...) plotData(as.datalist(data), ...)

#' Profile Likelihood Plot
#'
#' @param profs A [parframe] of profiles as returned by
#'   [profile()][profile.objfn], or a list of those.
#' @param ... Logical expressions to subset the plotted data.
#' @param maxvalue Numeric, the value of the objective difference at which
#'   profiles are cut off. Defaults to 5.
#' @param parlist Matrix or data frame with parameter columns, drawn as points.
#'   With a `value` column, its differences to the lowest profile value are
#'   used.
#' @param ncol Number of columns of the plot grid. `NULL` (default) lets
#'   [ggplot2::facet_wrap()] choose.
#' @param threshold Numeric, the horizontal lines and y-axis breaks above the
#'   optimum, which is always drawn. Defaults to the chi-square thresholds for
#'   68%, 90% and 95%. Names are used as axis labels. Pass the value of
#'   [profileThreshold()] to match [confint.parframe()].
#' @return A `ggplot` object with the plotted data frame as attribute `"data"`.
#' @seealso [profile()][profile.objfn], [plotPaths()]
#' @export
#' @examples
#' pars <- c(a = 1, b = 0.5)
#' obj <- constraintL2(mu = pars, sigma = 0.1)
#' profs <- profile(obj, pars, whichPar = c("a", "b"), limits = c(-1, 1),
#'                  cores = 1)
#' plotProfile(profs)
plotProfile <- function(profs,...) {
  UseMethod("plotProfile", profs)
}


#' Profile Likelihood: Plot of the Parameter Paths
#'
#' @param profs A [parframe] of profiles as returned by
#'   [profile()][profile.objfn], or a list of those.
#' @param ... Logical expressions to subset the plotted data.
#' @param whichPar Character or index vector, the profiled parameters whose
#'   paths are drawn. `NULL` (default) takes all.
#' @param sort Logical. `TRUE` sorts each pair of parameters so that a pair
#'   appears once; `FALSE` (default) keeps both orders.
#' @param relative Logical. `TRUE` (default) shifts each path to start at the
#'   optimum.
#' @param scales Character, `"fixed"` (default) or `"free"`.
#' @return A `ggplot` object with the plotted data frame as attribute `"data"`.
#' @seealso [profile()][profile.objfn], [plotProfile()], [plotPathsMulti()]
#' @export
#' @examples
#' pars <- c(a = 1, b = 0.5)
#' obj <- constraintL2(mu = pars, sigma = 0.1)
#' profs <- profile(obj, pars, whichPar = c("a", "b"), limits = c(-1, 1),
#'                  cores = 1)
#' plotPaths(profs)
#' plotPaths(profs, whichPar = "a")
plotPaths <- function(profs, ..., whichPar = NULL, sort = FALSE, relative = TRUE, scales = "fixed") {

  paths <- droplevels(dplyr::filter(.pathFrame(profs, whichPar, relative, sort), ...))
  axisLabel <- function(k) {
    label <- paste("parameter", k)
    if (relative) as.expression(bquote(paste(Delta, .(label)))) else label
  }

  p <- ggplot2::ggplot(paths, aes(x, y, group = interaction(name, proflist), colour = name,
                                  linetype = proflist)) +
    geom_path() +
    facet_wrap(vars(combination), scales = scales) +
    labs(x = axisLabel(1), y = axisLabel(2)) +
    scale_color_dMod(name = "profiled\nparameter") +
    scale_linetype_discrete("profile\nlist") +
    theme_dMod()
  attr(p, "data") <- paths
  p
}


# Parameter paths of profiles in long format: per profile set and profiled parameter n,
# one block per pair (first, partner) of whichPar and the parameters. With `relative`
# the paths start at zero at the optimum; `absoluteX` keeps the profiled parameter absolute.
.pathFrame <- function(profs, whichPar = NULL, relative = TRUE, sort = FALSE, absoluteX = FALSE) {
  sets <- .profileSets(profs)
  do.call(rbind, lapply(names(sets), function(s) {
    parts <- .profileSplit(sets[[s]])
    parameters <- attr(sets[[s]], "parameters")
    wp <- if (is.null(whichPar)) names(parts) else whichPar
    if (is.numeric(wp)) wp <- names(parts)[wp]
    lost <- setdiff(wp, names(parts))
    if (length(lost)) stop("Parameters not profiled: ", paste(lost, collapse = ", "), call. = FALSE)
    cmb <- expand.grid.alt(wp, parameters)
    if (sort) cmb <- t(apply(cmb, 1, sort))
    cmb <- cmb[cmb[, 1] != cmb[, 2], , drop = FALSE]
    cmb <- cmb[!duplicated(paste(cmb[, 1], cmb[, 2])), , drop = FALSE]
    do.call(rbind, lapply(wp, function(n) {
      P <- as.matrix(parts[[n]][, parameters, drop = FALSE])
      origin <- P[which.min(abs(parts[[n]]$constraint)), ]
      X <- P
      if (relative) P <- sweep(P, 2, origin)
      if (!absoluteX) X <- P
      data.frame(chisquare = rep(parts[[n]]$value, nrow(cmb)), name = n, proflist = s,
                 combination = rep(paste(cmb[, 1], cmb[, 2], sep = " - \n"), each = nrow(P)),
                 partner = rep(cmb[, 2], each = nrow(P)),
                 x = as.vector(X[, cmb[, 1]]), y = as.vector(P[, cmb[, 2]]))
    }))
  }))
}


#' Plot Fluxes Given a List of Flux Equations
#'
#' Evaluates flux expressions along a model prediction and draws them stacked
#' per condition.
#'
#' @param pouter Named numeric vector of outer parameters.
#' @param x A prediction function, called as `x(times, pouter, deriv = FALSE, ...)`.
#' @param times Numeric vector of time points for the prediction.
#' @param fluxEquations Character vector or list of flux expressions in the
#'   states and inner parameters, e.g. the `rates` of an [eqnlist]. Names are
#'   shown in the legend; without names the expressions are.
#' @param legendTitle Character, the legend title. Defaults to `"Fluxes:"`.
#' @param ... Further arguments passed to `x`, such as `fixed` or `conditions`.
#'
#' @return A `ggplot` object with the flux data frame as attribute `"out"`.
#' @seealso [subset.eqnlist()] to select reactions.
#' @examples
#' times <- 0:5
#' grid <- data.frame(name = "A", time = times, row.names = paste0("A", times))
#' x <- Xd(grid)
#' pars <- structure(exp(-times / 2), names = getParameters(x))
#' plotFluxes(pars, x, seq(0, 5, 0.1),
#'            c(production = "0.2", degradation = "0.5*A"))
#' @export
plotFluxes <- function(pouter, x, times, fluxEquations, legendTitle = "Fluxes:", ...){

  formulas <- vapply(fluxEquations, as.character, "")
  labels <- names(fluxEquations)
  if (is.null(labels)) labels <- formulas
  labels[!nzchar(labels)] <- formulas[!nzchar(labels)]
  calls <- lapply(formulas, str2lang)

  # functions of the expression language that base R lacks
  lang <- new.env(parent = baseenv())
  lang$Heaviside <- function(x) (sign(x) + 1) / 2
  lang$exp10 <- function(x) 10^x
  lang$piecewise <- function(...) {
    args <- list(...)
    n <- length(args) %/% 2L
    otherwise <- if (length(args) > 2L * n) args[[length(args)]] else NA_real_
    Reduce(function(k, rest) ifelse(args[[2L * k]], args[[2L * k - 1L]], rest),
           seq_len(n), otherwise, right = TRUE)
  }

  prediction <- .fillNames(x(times, pouter, deriv = FALSE, ...), "C")
  out <- do.call(rbind, unname(Map(function(frame, condition) {
    n <- nrow(frame)
    columns <- matrix(frame, n, dimnames = dimnames(frame))
    scope <- c(as.list(as.data.frame(columns)), as.list(unclass(attr(frame, "parameters"))))
    values <- vapply(calls, function(e) rep_len(as.numeric(eval(e, scope, lang)), n), numeric(n))
    data.frame(time = frame[, "time"], name = factor(rep(labels, each = n), levels = labels),
               value = as.vector(values), condition = condition)
  }, prediction, names(prediction))))

  p <- ggplot2::ggplot(out, aes(time, value, group = name, fill = name)) +
    geom_area(position = "stack", alpha = 0.3, colour = "darkgrey", linewidth = 0.4) +
    facet_wrap(vars(condition)) +
    scale_fill_dMod(name = legendTitle) +
    labs(x = "time", y = "flux contribution") +
    theme_dMod()
  attr(p, "out") <- out
  p
}


.stepDetect <- function(x, tol) {
  
  jumps <- 1
  while (TRUE) {
    i <- which(x - x[1] > tol)[1]
    if (is.na(i)) break
    jumps <- c(jumps, tail(jumps, 1) - 1 + i)
    x <- x[-seq(1, i - 1, 1)]
  }
  
  return(jumps)
  
  
}

#' Plot Objective Values of a Collection of Fits
#'
#' Draws the waterfall plot of a fit collection: the objective values in
#' ascending order against their rank. A converged fit is colored by its number
#' of iterations, an unconverged one is red.
#'
#' @param x Data frame with columns `value`, `converged` and `iterations`,
#'   e.g. a [parframe].
#' @param ... Logical expressions to subset `x`.
#' @param tol Maximal difference between neighboring objective values that
#'   still counts as one step. Defaults to 1.
#' @param showSteps Logical. `TRUE` marks the detected steps by dashed
#'   vertical lines labelled by their index. Defaults to `FALSE`.
#' @param pseudoLog Logical. `TRUE` draws the values on a pseudo-log10 scale,
#'   linear below 1 in absolute value and defined for negative values. Defaults
#'   to `FALSE`, a linear scale.
#' @return A `ggplot` object with the plotted data frame as attribute
#'   `"data"`; its column `delta` is the distance to the best fit.
#' @seealso [mstrust()], [as.parframe()], [plotPars()]
#' @export
#' @examples
#' fits <- parframe(data.frame(value = c(10, 10.2, 15, 40), converged = TRUE,
#'                             iterations = 20, a = c(1, 1.1, 2, 3)),
#'                  parameters = "a")
#' plotValues(fits)
plotValues <- function(x, ...) UseMethod("plotValues")


#' Plot Parameter Values of a Collection of Fits
#'
#' Draws one box per parameter, colored by the step of the objective value the
#' fits belong to.
#'
#' @param x A [parframe], e.g. from `as.parframe(mstrust(...))`.
#' @param tol Maximal difference between neighboring objective values that
#'   still counts as one step. Defaults to 1.
#' @param ... Logical expressions to subset `x`.
#' @return A `ggplot` object with the plotted data frame as attribute `"data"`.
#' @seealso [plotValues()], [mstrust()]
#' @export
#' @examples
#' fits <- parframe(data.frame(value = c(10, 10.2, 15, 40), converged = TRUE,
#'                             iterations = 20, a = c(1, 1.1, 2, 3),
#'                             b = c(0, 0.1, 1, 2)),
#'                  parameters = c("a", "b"))
#' plotPars(fits)
plotPars <- function(x, ...) UseMethod("plotPars")


#' Plot Residuals of a Collection of Fits
#'
#' Sums the squared weighted residuals of each fit over all variables not
#' named in `split` and plots the sums.
#'
#' @param parframe A [parframe], e.g. from `as.parframe(mstrust(...))`. Its
#'   column `index` labels the fits; without it the row numbers are used.
#' @param x Prediction function returning a [prdlist] with names matching
#'   `data`.
#' @param data A [datalist], i.e. a named list of data frames with columns
#'   `name`, `time`, `value` and `sigma`.
#' @param split Character vector of the variables to keep, from `"time"`,
#'   `"name"`, `"condition"` and `"index"`. `split[1]` is the x-axis,
#'   `split[2]` the color grouping (defaults to `split[1]`), further entries
#'   are facets. Defaults to `"condition"`.
#' @param errmodel Optional error model, an [obsfn]. With it the sums include
#'   `log(sigma^2)` of each data point.
#' @param ... Further arguments passed to `x`.
#'
#' @return A `ggplot` object with the summed residuals as attribute `"out"`.
#' @seealso [res()], [plotValues()]
#'
#' @examples
#' times <- 0:5
#' grid <- data.frame(name = "A", time = times, row.names = paste0("A", times))
#' x <- Xd(grid, condition = "C1")
#' data <- as.datalist(data.frame(name = "A", time = times,
#'                                value = exp(-times / 2), sigma = 0.1,
#'                                condition = "C1"))
#' fits <- parframe(data.frame(value = c(1, 2), converged = TRUE,
#'                             iterations = 10,
#'                             rbind(exp(-times / 2), exp(-times / 3)) |>
#'                               `colnames<-`(getParameters(x))),
#'                  parameters = getParameters(x))
#' plotResiduals(fits, x, data, c("time", "index"))
#'
#' @export
#' @importFrom dplyr group_by summarise across
#' @importFrom rlang data_sym syms
plotResiduals <- function(parframe, ...) UseMethod("plotResiduals")

#' @export
#' @rdname plotResiduals
plotResiduals.default <- function(parframe, x, data, split = "condition",
                                  errmodel = NULL, ...) {

  times <- sort(unique(c(0, unlist(lapply(data, `[[`, "time")))))
  fits <- if ("index" %in% colnames(parframe)) parframe$index else seq_len(nrow(parframe))

  # residuals of one fit, all conditions, labelled by fit and condition
  residualsOf <- function(j) {
    pred <- x(times, .parvecRow(parframe, j), deriv = FALSE, ...)
    data.table::rbindlist(lapply(names(pred), function(C) {
      err <- if (length(errmodel))
        errmodel(out = pred[[C]], pars = getParameters(pred[[C]]), conditions = C)[[C]]
      data.frame(index = as.character(fits[j]), res(data[[C]], pred[[C]], err), condition = C)
    }))
  }
  long <- data.table::rbindlist(lapply(seq_len(nrow(parframe)), residualsOf))

  terms <- long$weighted.residual^2
  if (length(errmodel)) terms <- terms + log(long$sigma^2)
  sums <- data.table::data.table(long[, split, with = FALSE], res = terms)
  out <- as.data.frame(sums[, list(res = sum(res)), keyby = split])

  colourBy <- as.name(split[min(2L, length(split))])
  p <- ggplot(out, aes(!!as.name(split[1]), res, colour = !!colourBy, group = !!colourBy)) +
    theme_dMod() +
    geom_point() +
    geom_line()
  if (length(split) > 2L) p <- p + facet_wrap(split[-(1:2)])
  attr(p, "out") <- out
  p
}


# Plot generics for the layer branches -------------------------------------


#' Convergence Trace of an Iterative Fit
#'
#' Generic for plotting the quantities an iterative fitting method records,
#' one panel each. dMod2 defines no method for it; packages and extensions
#' that return iterative fit objects do.
#' @param x Object to plot.
#' @param ... Method-specific arguments.
#' @return A `ggplot` object, as returned by the method.
#' @export
plotTrace <- function(x, ...) UseMethod("plotTrace", x)




# Long-format samples for ggplot.
.bayes_samples_long <- function(samples, par_names = NULL) {
  if (!is.matrix(samples)) stop("samples must be a matrix.")
  if (is.null(colnames(samples))) {
    if (is.null(par_names))
      par_names <- paste0("p", seq_len(ncol(samples)))
    colnames(samples) <- par_names
  }
  N <- nrow(samples)
  data.frame(
    sample    = rep(seq_len(N), times = ncol(samples)),
    parameter = factor(rep(colnames(samples), each = N),
                       levels = colnames(samples)),
    value     = as.vector(samples),
    stringsAsFactors = FALSE)
}


#' Pair Plot of Parameter Samples
#'
#' Generic for a corner plot of sampled parameters: scatter plots of pairs in
#' the lower triangle, marginal densities on the diagonal. dMod2 defines no
#' method for it; packages and extensions that return samples do.
#'
#' @param x Object holding parameter samples.
#' @param ... Method-specific arguments.
#' @return A `ggplot` object, as returned by the method.
#' @export
plotPairs <- function(x, ...) UseMethod("plotPairs", x)




## ---- profile / parameter-path plotting ---------
#' Plot an Array of Trajectories Along the Profile of a Parameter
#'
#' Predicts the model at parameter sets taken along one profile, starting at
#' the optimum, and colors the trajectories by the value of the profiled
#' parameter.
#'
#' @param par Character, the profiled parameter.
#' @param profs A [parframe] of profiles as returned by
#'   [profile()][profile.objfn].
#' @param prd A prediction function, e.g. `g*x*p`, called as
#'   `prd(times, pars, deriv = FALSE)`.
#' @param times Numeric vector of time points for the prediction.
#' @param ... Logical expressions in the columns `time`, `name`, `value`,
#'   `condition` and those of `covtable` to subset the plotted data.
#' @param direction `"up"` (default) or `"down"`, the side of the profile that
#'   is traced from the optimum.
#' @param covtable Condition table, e.g. from [covariates()], merged with the
#'   predictions by condition so that `...` can refer to its columns. Default
#'   `NULL`.
#' @param nSim Number of trajectories, evenly spaced along the profile side
#'   and starting at the optimum; at most the number of profile points on that
#'   side. Defaults to 4.
#'
#' @return A `ggplot` object.
#' @seealso [plotProfile()], [plotPathsMulti()]
#' @examples
#' times <- 0:5
#' grid <- data.frame(name = "A", time = times, row.names = paste0("A", times))
#' x <- Xd(grid)
#' pars <- structure(exp(-times / 2), names = getParameters(x))
#' obj <- constraintL2(mu = pars, sigma = 0.1)
#' profs <- profile(obj, pars, whichPar = "A2", limits = c(-1, 1), cores = 1)
#' plotArray("A2", profs, x, seq(0, 5, 0.1))
#' @export
plotArray <- function(par, profs, prd, times, ..., direction = c("up", "down"),
                      covtable = NULL, nSim = 4) {

  direction <- match.arg(direction)
  prof <- as.data.frame(profs)[profs$whichPar == par, , drop = FALSE]
  if (!nrow(prof)) stop("plotArray: '", par, "' is not profiled.", call. = FALSE)
  origin <- which.min(abs(prof$constraint))
  side <- if (direction == "up") origin:nrow(prof) else origin:1
  side <- side[unique(round(seq(1, length(side), length.out = min(nSim, length(side)))))]
  parameters <- attr(profs, "parameters")

  # long trajectories of profile point k, NULL with a warning if the prediction fails
  trajectories <- function(k) {
    pars <- unlist(prof[side[k], parameters])
    pred <- tryCatch(prd(times, as.parvec(pars), deriv = FALSE), error = function(e) {
      warning("plotArray: parameter set ", k, " failed: ", conditionMessage(e), call. = FALSE)
      NULL
    })
    if (is.null(pred)) return(NULL)
    plain <- lapply(.fillNames(pred, ""), function(m) matrix(m, nrow(m), dimnames = dimnames(m)))
    cbind(wide2long(plain), set = k, ParValue = pars[[par]])
  }
  out <- as.data.frame(data.table::rbindlist(lapply(seq_along(side), trajectories)))

  if (length(covtable)) {
    if (!"condition" %in% names(covtable)) covtable <- cbind(condition = rownames(covtable), covtable)
    out <- merge(out, covtable, by = "condition")
  }
  out <- dplyr::filter(out, ...)

  ggplot2::ggplot(out, aes(time, value, group = set, colour = ParValue)) +
    geom_line() +
    facet_grid(rows = vars(name), cols = vars(condition), scales = "free_y") +
    scale_color_viridis_c(name = par) +
    labs(x = "time", y = "value") +
    theme_dMod()
}


.findEmptyCorner <- function(x, y) {
  xmid <- (min(x, na.rm = TRUE) + max(x, na.rm = TRUE)) / 2
  ymid <- (min(y, na.rm = TRUE) + max(y, na.rm = TRUE)) / 2
  
  corners <- list(
    bottom_left  = c(0.05, 0.05),
    bottom_right = c(0.95, 0.05),
    top_left     = c(0.05, 0.95),
    top_right    = c(0.95, 0.95)
  )
  
  counts <- c(
    bottom_left  = sum(x <= xmid & y <= ymid, na.rm = TRUE),
    bottom_right = sum(x >  xmid & y <= ymid, na.rm = TRUE),
    top_left     = sum(x <= xmid & y >  ymid, na.rm = TRUE),
    top_right    = sum(x >  xmid & y >  ymid, na.rm = TRUE)
  )
  
  corners[[which.min(counts)]]
}

#' Profile Likelihood: Plot All Parameter Paths of One Profile Together
#'
#' Draws, for each profiled parameter, the changes of all other parameters
#' along its profile in one panel. The `nPars` parameters with the largest
#' change are colored and named, the others are gray. Requires the package
#' cowplot.
#'
#' @param profs A [parframe] of profiles as returned by
#'   [profile()][profile.objfn].
#' @param whichPar Character vector, the profiled parameters to draw.
#' @param nPars Integer, number of colored and named paths. Defaults to 5.
#' @param normalizePaths Logical. `TRUE` scales each path to a maximum
#'   absolute value of 1. Defaults to `FALSE`.
#'
#' @return A `ggplot` object; for several `whichPar` the grid of
#'   [cowplot::plot_grid()].
#' @seealso [plotPaths()], [plotProfilesAndPaths()]
#' @examplesIf requireNamespace("cowplot", quietly = TRUE)
#' pars <- c(a = 1, b = 0.5, c = 2)
#' obj <- constraintL2(mu = pars, sigma = 0.1)
#' profs <- profile(obj, pars, whichPar = c("a", "b"), limits = c(-1, 1),
#'                  cores = 1)
#' plotPathsMulti(profs, c("a", "b"), nPars = 2)
#' @export
#' @import data.table
plotPathsMulti <- function(profs, whichPar, nPars = 5, normalizePaths = FALSE) {
  .require_ns("cowplot", "plotPathsMulti()")
  plots <- lapply(whichPar, .pathPanel, profs = profs, nPars = nPars,
                  normalizePaths = normalizePaths)
  if (length(plots) == 1L) return(plots[[1]])
  cowplot::plot_grid(plotlist = plots)
}

# All paths of the profile of whichPar in one panel, the nPars largest colored
.pathPanel <- function(whichPar, profs, nPars, normalizePaths) {
  data <- data.table::as.data.table(.pathFrame(profs, whichPar, absoluteX = TRUE))
  if (normalizePaths) data[, y := if (any(y != 0)) y / max(abs(y)) else y, by = partner]
  data[, max.dev := max(abs(y)), by = partner]
  data.table::setorder(data, -max.dev)
  top <- head(unique(data$partner), nPars)
  data[, label := ifelse(partner %in% top, partner, "Others")]
  colors <- c(setNames(dMod_palette(nPars + 1L)[-1L][seq_along(top)], top), Others = "gray")
  corner <- .findEmptyCorner(data$x, data$y)

  p <- ggplot2::ggplot(data, aes(x = x, y = y, color = label, group = partner)) +
    geom_line() +
    xlab(whichPar) + ylab(expression(Delta ~ p[j])) +
    scale_color_manual(name = NULL, values = colors, breaks = c(top, "Others")) + theme_dMod() +
    theme(legend.position = "inside",
          legend.position.inside = corner,
          legend.justification = corner,
          legend.background = element_rect(fill = alpha("white", 0.85), colour = "black",
                                           linewidth = 0.3),
          legend.key.size = unit(0.4, "cm"),
          legend.margin = margin(2, 4, 2, 4),
          legend.text = element_text(size = 7))
  attr(p, "data") <- as.data.frame(data)
  p
}


#' Profile Likelihood: Plot Profiles Along with Their Parameter Paths
#'
#' For each profiled parameter, draws the profile above the paths of
#' [plotPathsMulti()]. Requires the package cowplot.
#'
#' @param profs A [parframe] of profiles as returned by
#'   [profile()][profile.objfn].
#' @param whichPar Character vector, the profiled parameters to draw.
#' @param nPars Integer, number of colored and named paths. Defaults to 5.
#' @param ncol Number of columns of the plot grid. Defaults to 3.
#' @param normalizePaths Logical. `TRUE` scales each path to a maximum
#'   absolute value of 1. Defaults to `FALSE`.
#' @param modes Character vector, the contributions to the objective drawn in
#'   the profile plot besides the total. Defaults to `c("data", "prior")`.
#' @param ... Further arguments passed to [cowplot::plot_grid()].
#'
#' @return A `ggplot` object, the grid of [cowplot::plot_grid()].
#' @seealso [plotProfile()], [plotPathsMulti()]
#' @examplesIf requireNamespace("cowplot", quietly = TRUE)
#' pars <- c(a = 1, b = 0.5, c = 2)
#' obj <- constraintL2(mu = pars, sigma = 0.1)
#' profs <- profile(obj, pars, whichPar = c("a", "b"), limits = c(-1, 1),
#'                  cores = 1)
#' plotProfilesAndPaths(profs, c("a", "b"), nPars = 2, ncol = 2)
#'
#' @export
plotProfilesAndPaths <- function(profs, whichPar, nPars = 5, ncol = 3,
                                 normalizePaths = FALSE,
                                 modes = c("data", "prior"), ...) {
  dots <- list(...)
  .require_ns("cowplot", "plotProfilesAndPaths()")

  # Save original obj.attributes before any subsetting drops them
  orig_oa <- attr(profs, "obj.attributes")
  filtered_oa <- if (!is.null(orig_oa)) intersect(orig_oa, modes) else NULL
  
  profs <- profs[profs$whichPar %in% whichPar, ]
  
  cleanProfilePlot <- function(prof_sub) {
    # Remove columns for unwanted modes, so plotProfile cannot plot them
    cols_to_drop <- setdiff(orig_oa, modes)
    prof_sub <- prof_sub[, !(colnames(prof_sub) %in% cols_to_drop), drop = FALSE]
    attr(prof_sub, "obj.attributes") <- filtered_oa
    p <- plotProfile(prof_sub)
    
    # the total and the requested contributions
    pdata <- attr(p, "data")
    pdata <- pdata[pdata$mode %in% c("total", modes), , drop = FALSE]
    
    threshold <- c(1, 2.7, 3.84)
    p_new <- ggplot(pdata, aes(x = par, y = delta, group = interaction(proflist, mode), 
                               color = proflist, linetype = mode)) +
      facet_wrap(~name, scales = "free_x") +
      geom_hline(yintercept = threshold, lty = 2, color = "gray") +
      geom_line() +
      geom_point(data = subset(pdata, is.zero)) +
      ylab(expression(paste("CL /", Delta * chi^2))) +
      scale_y_continuous(breaks = c(1, 2.7, 3.84), 
                         labels = c("68% / 1   ", "90% / 2.71", "95% / 3.84"),
                         limits = c(NA, 5)) +
      xlab("parameter value") +
      labs(title = NULL, x = NULL, linetype = "contrib") +
      theme(
        strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      ) +
      guides(color = "none", fill = "none") +
      theme(legend.position = "none")
  }
  
  stacked_list <- vector("list", length(whichPar))
  
  for (z in seq_along(whichPar)) {
    prof_sub <- profs[profs$whichPar == whichPar[z], ]
    
    p_prof_noleg <- cleanProfilePlot(prof_sub)
    
    p_paths <- plotPathsMulti(prof_sub, whichPar[z], nPars, normalizePaths = normalizePaths) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    aligned_pair <- cowplot::align_plots(p_prof_noleg, p_paths, align = "v", axis = "tb")
    stacked_list[[z]] <- cowplot::plot_grid(aligned_pair[[1]], aligned_pair[[2]],
                                            ncol = 1, rel_heights = c(1, 0.7), align = "v", axis = "tb")
  }
  
  body <- do.call(cowplot::plot_grid,
                  c(list(plotlist = stacked_list, ncol = ncol), dots))
  
  return(body)
}
