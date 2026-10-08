
# ggplot2 / dplyr NSE column references; declared so R CMD check does not
# flag them as undefined globals.
utils::globalVariables(c("predicted", "observed", "sd_est", "iter", "level"))


# Custom interface to ggplot2 ---

#' Open a Plot in an External PDF Viewer
#'
#' Saves a plot as PDF in [tempdir()] and opens it with an external command.
#'
#' @param plot A `ggplot2` plot object. Defaults to [ggplot2::last_plot()].
#' @param command Character, the shell command that opens the PDF file.
#'   Defaults to `"xdg-open"`, which exists on Linux; use `"open"` on macOS.
#' @param ... Arguments passed to [ggplot2::ggsave()].
#' @return The path of the PDF file, invisibly.
#' @export
ggopen <- function(plot = last_plot(), command = "xdg-open", ...) {
  filename <- tempfile(pattern = "Rplot", fileext = ".pdf")
  ggsave(filename = filename, plot = plot, ...)
  system(command = paste(command, filename))
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
  colors <- list(
    medium = c(gray = '#737373', red = '#F15A60', green = '#7AC36A', blue = '#5A9BD4', orange = '#FAA75B', purple = '#9E67AB', maroon = '#CE7058', magenta = '#D77FB4'),
    dark = c(black = '#010202', red = '#EE2E2F', green = '#008C48', blue = '#185AA9', orange = '#F47D23', purple = '#662C91', maroon = '#A21D21', magenta = '#B43894'),
    light = c(gray = '#CCCCCC', red = '#F2AFAD', green = '#D9E4AA', blue = '#B8D2EC', orange = '#F3D1B0', purple = '#D5B2D4', maroon = '#DDB9A9', magenta = '#EBC0DA')
  )
  gray <- colors$medium["gray"]
  black <- colors$dark["black"]

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
#
# dE is the smallest pairwise CIE2000 distance over normal and simulated
# dichromatic vision; 10 or above counts as colorblind-safe. Re-derive with
# colorspace::deutan/protan/tritan and farver::compare_colour().

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
    # createPalette() shifts the seeds slightly to maximize distinctness across
    # the whole set; keep the original seeds verbatim and only borrow the tail.
    # `range` caps the luminance: unconstrained, the tail wanders up to ~224 and
    # those colors are invisible as lines on a white panel.
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


# Other ---------------------------------------------

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
#' mydata1 <- data.frame(name = c("A", "B"), time = 0:5, value = 0:5, sigma = .1)
#' coordTransform(mydata1, "log(value)")
#' coordTransform(mydata1, list(A = "exp(value)", B = "sqrt(value)"))
coordTransform <- function(data, transformations) {
  
  mynames <- unique(as.character(data$name))
  
  # Replicate transformation if not a list
  if (!is.list(transformations))
    transformations <- as.list(structure(rep(transformations, length(mynames)), names = mynames))
  
  out <- do.call(rbind, lapply(mynames, function(n) {
    
    subdata <- subset(data, name == n)
    
    if (n %in% names(transformations)) {
      
      mysymbol <- getSymbols(transformations[[n]])[1]
      mytrafo <- replaceSymbols(mysymbol, "value", transformations[[n]])
      mytrafo <- parse(text = mytrafo)
      
      if ("sigma" %in% colnames(subdata))
        subdata$sigma <- abs(with(subdata, eval(D(mytrafo, "value")))) * subdata$sigma
      subdata$value <- with(subdata, eval(mytrafo))
      
    }
    
    return(subdata)
    
  }))
  
  
  return(out)
  
  
}


# Method dispatch for plotX functions -------------


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
#' @importFrom graphics par
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
plotData.data.frame <- function(data, ...) {
  plotData.datalist(as.datalist(data), ...)
}

#' Profile Likelihood Plot
#'
#' @param profs A [parframe] of profiles as returned by [profile()], or a
#'   list of those.
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
#' @seealso [profile()], [plotPaths()]
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
#' @param profs A [parframe] of profiles as returned by [profile()], or a
#'   list of those.
#' @param ... Logical expressions to subset the plotted data.
#' @param whichPar Character or index vector, the profiled parameters whose
#'   paths are drawn. `NULL` (default) takes all.
#' @param sort Logical. `TRUE` sorts each pair of parameters so that a pair
#'   appears once; `FALSE` (default) keeps both orders.
#' @param relative Logical. `TRUE` (default) shifts each path to start at the
#'   optimum.
#' @param scales Character, `"fixed"` (default) or `"free"`.
#' @return A `ggplot` object with the plotted data frame as attribute `"data"`.
#' @seealso [profile()], [plotProfile()], [plotPathsMulti()]
#' @export
#' @examples
#' pars <- c(a = 1, b = 0.5)
#' obj <- constraintL2(mu = pars, sigma = 0.1)
#' profs <- profile(obj, pars, whichPar = c("a", "b"), limits = c(-1, 1),
#'                  cores = 1)
#' plotPaths(profs)
#' plotPaths(profs, whichPar = "a")
plotPaths <- function(profs, ..., whichPar = NULL, sort = FALSE, relative = TRUE, scales = "fixed") {
  
  if ("parframe" %in% class(profs)) 
    arglist <- list(profs)
  else
    arglist <- as.list(profs)
  
  
  if (is.null(names(arglist))) {
    profnames <- 1:length(arglist)
  } else {
    profnames <- names(arglist)
  }
  
  
  data <- do.call(rbind, lapply(1:length(arglist), function(i) {
    # choose a proflist
    proflist <- as.data.frame(arglist[[i]])
    parameters <- attr(arglist[[i]], "parameters")
    
    if (is.data.frame(proflist)) {
      whichPars <- unique(proflist$whichPar)
      proflist <- lapply(whichPars, function(n) {
        with(proflist, proflist[whichPar == n, ])
      })
      names(proflist) <- whichPars
    }
    
    if (is.null(whichPar)) whichPar <- names(proflist)
    if (is.numeric(whichPar)) whichPar <- names(proflist)[whichPar]
    
    subdata <- do.call(rbind, lapply(whichPar, function(n) {
      # matirx
      paths <- as.matrix(proflist[[n]][, parameters])
      values <- proflist[[n]][, "value"]
      origin <- which.min(abs(proflist[[n]][, "constraint"]))
      if (relative) 
        for(j in 1:ncol(paths)) paths[, j] <- as.numeric(paths[, j]) - as.numeric(paths[origin, j])
      
      combinations <- expand.grid.alt(whichPar, colnames(paths))
      if (sort) combinations <- apply(combinations, 1, sort) else combinations <- apply(combinations, 1, identity)
      combinations <- submatrix(combinations, cols = -which(combinations[1,] == combinations[2,]))
      combinations <- submatrix(combinations, cols = !duplicated(paste(combinations[1,], combinations[2,])))
      
      
      path.data <- do.call(rbind, lapply(1:dim(combinations)[2], function(j) {
        data.frame(chisquare = values, 
                   name = n,
                   proflist = profnames[i],
                   combination = paste(combinations[,j], collapse = " - \n"),
                   x = paths[, combinations[1,j]],
                   y = paths[, combinations[2,j]])
      }))
      
      return(path.data)
      
    }))
    
    return(subdata)
    
  }))
  
  data$proflist <- as.factor(data$proflist)
  
  
  if (relative)
    axis.labels <- c(expression(paste(Delta, "parameter 1")), expression(paste(Delta, "parameter 2")))  
  else
    axis.labels <- c("parameter 1", "parameter 2")
  
  
  data <- droplevels(subset(data, ...))
  data$y <- as.numeric(data$y)
  data$x <- as.numeric(data$x)
  
  suppressMessages(
    p <- ggplot(data, aes(x = x, y = y, group = interaction(name, proflist), color = name, lty = proflist)) + 
      facet_wrap(~combination, scales = scales) + 
      geom_path() + #geom_point(aes=aes(size=1), alpha=1/3) +
      xlab(axis.labels[1]) + ylab(axis.labels[2]) +
      scale_linetype_discrete(name = "profile\nlist") +
      scale_color_dMod(name = "profiled\nparameter")
  )
  
  attr(p, "data") <- data
  return(p)
  
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
#' @param nameFlux Character, the legend title. Defaults to `"Fluxes:"`.
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
plotFluxes <- function(pouter, x, times, fluxEquations, nameFlux = "Fluxes:", ...){

  if (is.null(names(fluxEquations))) names(fluxEquations) <- fluxEquations

  # Flux values evaluated in R.
  exprs <- lapply(fluxEquations, function(e) str2lang(as.character(e)))
  fluxEnv <- list2env(list(
    Heaviside = function(x) ifelse(x < 0, 0, ifelse(x == 0, 0.5, 1)),
    exp10 = function(x) 10^x,
    piecewise = function(...) {
      a <- list(...); n <- length(a)
      val <- if (n %% 2L) a[[n]] else NA_real_
      for (i in rev(seq_len(n %/% 2L))) val <- ifelse(a[[2L * i]], a[[2L * i - 1L]], val)
      val
    }), parent = baseenv())
  flux <- function(values, n)
    do.call(cbind, lapply(exprs, function(e) rep_len(as.numeric(eval(e, values, fluxEnv)), n)))
  prediction.all <- x(times, pouter, deriv = FALSE, ...)
  names.prediction.all <- names(prediction.all)
  if (is.null(names.prediction.all)) names.prediction.all <- paste0("C", 1:length(prediction.all))

  out <- lapply(1:length(prediction.all), function(cond) {
    prediction <- prediction.all[[cond]]
    pinner <- attr(prediction, "parameters")
    values <- c(as.list(as.data.frame(unclass(prediction))), as.list(unclass(pinner)))
    fluxes <- cbind(time = prediction[, "time"], flux(values, nrow(prediction)))
    return(fluxes)
  }); names(out) <- names.prediction.all
  out <- wide2long(out)

  cbPalette <- c("#999999", "#E69F00", "#F0E442", "#56B4E9", "#009E73", "#0072B2",
                 "#D55E00", "#CC79A7","#CC6666", "#9999CC", "#66CC99","red", "blue", "green","black")

  P <- ggplot(out, aes(x = time, y = value, group = name, fill = name, log = "y")) +
    facet_wrap(~condition) + scale_fill_manual(values = cbPalette, name = nameFlux) +
    geom_density(stat = "identity", position = "stack", alpha = 0.3, color = "darkgrey", linewidth = 0.4) +
    xlab("time") + ylab("flux contribution")

  attr(P, "out") <- out

  return(P)

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
#' Draws the waterfall plot of a fit collection: the objective values minus the
#' best one, in ascending order against their rank, on a pseudo-log10 scale
#' that is linear below 1. A converged fit is drawn as a circle, an
#' unconverged one as a triangle.
#'
#' @param x Data frame with columns `value`, `converged` and `iterations`,
#'   e.g. a [parframe].
#' @param ... Logical expressions to subset `x`.
#' @param tol Maximal difference between neighboring objective values that
#'   still counts as one step. Defaults to 1.
#' @param showSteps Logical. `TRUE` marks the detected steps by dashed
#'   vertical lines labelled by their index. Defaults to `FALSE`.
#' @return A `ggplot` object with the plotted data frame as attribute
#'   `"data"`; its column `delta` is the distance to the best fit.
#' @seealso [mstrust()], [as.parframe()], [plotPars()]
#' @export
#' @examples
#' fits <- parframe(data.frame(value = c(10, 10.2, 15, 40), converged = TRUE,
#'                             iterations = 20, a = c(1, 1.1, 2, 3)),
#'                  parameters = "a")
#' plotValues(fits)
plotValues <- function(x,...) {
  UseMethod("plotValues", x)
}


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
plotPars <- function(x,...) {
  UseMethod("plotPars", x)
}


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

  timesD <- sort(unique(c(0, unlist(lapply(data, function(d) d$time)))))
  
  if (!("index" %in% colnames(parframe))) {
    parframe$index <- seq_len(nrow(parframe))
  }
  
  # --- Compute residuals for all fits and conditions ---
  out <- do.call(rbind, lapply(seq_len(nrow(parframe)), function(j) {
    pred <- x(timesD, as.parvec(parframe, j), deriv = FALSE, ...)
    
    out_con <- do.call(rbind, lapply(names(pred), function(con) {
      err <- NULL
      if (!is.null(errmodel)) {
        err <- errmodel(out = pred[[con]], pars = getParameters(pred[[con]]), conditions = con)
      }
      out <- res(data[[con]], pred[[con]], err[[con]])
      cbind(out, condition = con)
    }))
    
    cbind(index = as.character(parframe[j, "index"]), out_con)
  }))
  
  # --- Summarize residuals ---
  out <- dplyr::group_by(out, across(all_of(split)))
  
  if (!is.null(errmodel)) {
    out <- dplyr::summarise(out, res = sum(weighted.residual^2 + log(sigma^2)), .groups = "drop")
  } else {
    out <- dplyr::summarise(out, res = sum(weighted.residual^2), .groups = "drop")
  }
  
  out <- as.data.frame(out)
  
  # --- Build aesthetics ---
  groupvar <- if (length(split) > 1) split[2] else split[1]
  
  p <- ggplot(out, aes(x = !!rlang::data_sym(split[1]), 
                       y = res, 
                       color = !!rlang::data_sym(groupvar), 
                       group = !!rlang::data_sym(groupvar))) + 
    theme_dMod() + 
    geom_point() + 
    geom_line()
  
  if (length(split) > 2) {
    facet_vars <- rlang::syms(split[3:length(split)])
    p <- p + facet_wrap(vars(!!!facet_vars))
  }
  
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




## ---- profile / parameter-path plotting (moved from toolsSvenja.R) ---------
#' Plot an Array of Trajectories Along the Profile of a Parameter
#'
#' Predicts the model at parameter sets taken along one profile, starting at
#' the optimum, and colors the trajectories by the value of the profiled
#' parameter. Requires the package purrr.
#'
#' @param par Character, the profiled parameter.
#' @param profs A [parframe] of profiles as returned by [profile()].
#' @param prd A prediction function, e.g. `g*x*p`, called as
#'   `prd(times, pars, deriv = FALSE)`.
#' @param times Numeric vector of time points for the prediction.
#' @param direction `"up"` (default) or `"down"`, the side of the profile that
#'   is traced from the optimum.
#' @param covtable Condition table, e.g. from [covariates()], merged with the
#'   predictions by condition so that `...` can refer to its columns, or
#'   `NULL`. Has no default.
#' @param ... Logical expressions to subset the plotted data; used only with
#'   a `covtable`.
#' @param nsimus Number of trajectories. Defaults to 4.
#'
#' @return A `ggplot` object.
#' @author Svenja Kemmer, \email{svenja.kemmer@@fdm.uni-freiburg.de}
#' @seealso [plotProfile()], [plotPathsMulti()]
#' @examplesIf requireNamespace("purrr", quietly = TRUE)
#' times <- 0:5
#' grid <- data.frame(name = "A", time = times, row.names = paste0("A", times))
#' x <- Xd(grid)
#' pars <- structure(exp(-times / 2), names = getParameters(x))
#' obj <- constraintL2(mu = pars, sigma = 0.1)
#' profs <- profile(obj, pars, whichPar = "A2", limits = c(-1, 1), cores = 1)
#' plotArray("A2", profs, x, seq(0, 5, 0.1), covtable = NULL)
#' @export
#' @import data.table
plotArray <- function (par, profs, prd, times, direction = c("up", "down"), covtable, ..., nsimus = 4) {

  direction <- match.arg(direction)

  # select subframe from profiles
  mysub <- profs %>% as.data.table() %>% .[whichPar == par, ]
  mysub[, ID := 1:nrow(mysub)]
  
  # get ID of bestfit (constraint is 0 for bestfit)
  bestID <- mysub[constraint == 0.00]$ID
  if(direction == "up") mysubF <- mysub[ID >= bestID]  
  if(direction == "down") mysubF <- mysub[ID <= bestID]
  
  # select rows according to simulation number
  partable <- mysubF[seq(1, nrow(mysubF), (round(nrow(mysubF)/nsimus)))]
  
  # remove non_parameter names
  no_pars <- c("value", "constraint", "stepsize", "gamma", "whichPar", "data", "condition_obj", "AIC", "BIC", "prior", "ID", "chisquare")
  no_pars <- intersect(no_pars, names(partable))
  partable %>% .[, (no_pars) := NULL]
  
  # make predictions
  predictionDT <- .predictArray(prd, times, pars = partable, whichpar = par)
  out_plot <- copy(predictionDT)
  
  # use covtable for subsetting of the plot
  if(!is.null(covtable)) {
    if(!"condition" %in% names(covtable)){
      covtable <- as.data.table(covtable, keep.rownames = "condition")
    } else covtable <- as.data.table(covtable)
    out_plot <- merge(out_plot, covtable, by = "condition")
    out_plot <- out_plot[...]
  }
  
  # plot
  P <- ggplot(out_plot , aes(x = time, y = value, group = ParValue, color = ParValue)) +
    facet_grid(name~condition, scales = "free_y") +
    geom_line(linewidth = 1) + 
    theme_dMod(base_size = 18) + scale_color_viridis_c() +
    theme(legend.position = "top", legend.key.size = unit(0.6,"cm")) + 
    theme(axis.line = element_line(colour = "black"), 
          panel.grid.major = element_line(colour = "grey97"), 
          panel.grid.minor = element_line(colour = "grey97"), 
          panel.background = element_blank()) +
    xlab("time") +
    ylab(paste0("value"))
  
  return(P)
}

.predictArray <- function (prd, times, pars = partable, whichpar = par, keep_names = NULL, FLAGverbose = FALSE, FLAGverbose2 = FALSE, FLAGbrowser = FALSE, ...) {
  .require_ns("purrr", ".predictArray()")
  if (FLAGverbose2) cat("Simulating", "\n")
  out <- lapply(1:nrow(pars), function(i) {
    if (FLAGverbose) cat("Parameter set", i, "\n")
    if (FLAGbrowser) browser()
    mypar <- pars[i,] %>% as.numeric()
    parval <- round(pars[i,][[whichpar]], digits = 2)
    names(mypar) <- names(pars)
    mypar <- as.parvec(mypar)
    prediction <- try(prd(times, mypar, deriv = FALSE, ...))
    if (inherits(prediction, "try-error")) {
      warning("parameter set ", i, " failed\n")
      return(NULL)
    }
    prediction <- purrr::imap(prediction, function(.x,.y){
      .x <- data.table(.x)
      if (!is.null(keep_names))
        .x[, (setdiff(names(.x), c(keep_names, "time"))) := NULL]
      .x[, `:=`(condition = .y, ParValue = parval)]
      .x
    })
    melt(rbindlist(prediction), variable.name = "name", value.name = "value", id.vars = c("time", "condition", "ParValue"))
  })
  if (FLAGverbose2) cat("postprocessing", "\n")
  out <- rbindlist(out[!is.null(out)])
  out
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

#' @importFrom ggplot2 ggplot
#' @noRd
PlotPaths <- function(profs=myprofiles, ..., whichPar, sort = FALSE, relative = TRUE, scales = "fixed", multi = TRUE, n_pars = 5, normalizePaths = FALSE) {
  
  if ("parframe" %in% class(profs)) {
    arglist <- list(profs)
  } else {
    arglist <- as.list(profs)
  }
  
  if (is.null(names(arglist))) {
    profnames <- 1:length(arglist)
  } else {
    profnames <- names(arglist)
  }
  
  
  data <- do.call(rbind, lapply(1:length(arglist), function(i) {
    # choose a proflist
    proflist <- as.data.frame(arglist[[i]])
    parameters <- attr(arglist[[i]], "parameters")
    
    if (is.data.frame(proflist)) {
      whichPars <- unique(proflist$whichPar)
      proflist <- lapply(whichPars, function(n) {
        with(proflist, proflist[whichPar == n, ])
      })
      names(proflist) <- whichPars
    }
    
    if (is.null(whichPar)) whichPar <- names(proflist)
    if (is.numeric(whichPar)) whichPar <- names(proflist)[whichPar]
    
    subdata <- do.call(rbind, lapply(whichPar, function(n) {
      # matrix
      paths <- as.matrix(proflist[[n]][, parameters])
      values <- proflist[[n]][, "value"]
      origin <- which.min(abs(proflist[[n]][, "constraint"]))
      
      # Save absolute values of profiled parameter before relativizing
      abs_profiled <- as.numeric(paths[, n])
      
      if (relative) 
        for(j in 1:ncol(paths)) paths[, j] <- as.numeric(paths[, j]) - as.numeric(paths[origin, j])
      
      # Restore absolute values for the profiled parameter (x-axis always absolute)
      paths[, n] <- abs_profiled
      
      combinations <- expand.grid.alt(whichPar, colnames(paths))
      if (sort) combinations <- apply(combinations, 1, sort) else combinations <- apply(combinations, 1, identity)
      combinations <- submatrix(combinations, cols = -which(combinations[1,] == combinations[2,]))
      combinations <- submatrix(combinations, cols = !duplicated(paste(combinations[1,], combinations[2,])))
      
      
      path.data <- do.call(rbind, lapply(1:dim(combinations)[2], function(j) {
        data.frame(chisquare = values, 
                   name = n,
                   proflist = profnames[i],
                   combination = paste(combinations[,j], collapse = " - \n"),
                   x = paths[, combinations[1,j]],
                   y = paths[, combinations[2,j]])
      }))
      
      if(multi) path.data <- path.data %>% as.data.table %>% .[, partner := tstrsplit(as.character(combination), "\n", fixed=TRUE, keep = 2)]
      
      
      return(path.data)
      
    }))
    
    return(subdata)
    
  }))
  
  data$proflist <- as.factor(data$proflist)
  
  if (relative){
    axis.labels <- c("parameter 1", expression(Delta ~ p[j]))
  } else {
    axis.labels <- c("parameter 1", "parameter 2")
  }
  
  data <- droplevels(subset(data, ...))
  removeBecauseNonsense <- c("value", "constraint", "stepsize", "chisquare", "data", "prior", "gamma", "whichPar")
  data <- data[!(partner %in% removeBecauseNonsense)]
  data$y <- as.numeric(data$y)
  data$x <- as.numeric(data$x)
  
  if (normalizePaths == TRUE) {
    data[, y := (ifelse(max(abs(y)) == 0, 0, y / abs(max(abs(y))))), by = combination] # if path is y, just return 0
    # data[, y := (2 * (y - min(y)) / (max(y) - min(y))) - 1, by = combination]
    removedCombinations <- unique(data[!is.finite(y), combination])
    data <- data[is.finite(y)]
    
    if(length(removedCombinations)>0) {warning(paste0("The following combinations have been removed due to failed paths:\n\t",paste(str_remove_all(removedCombinations, "\n"), collapse = "\n\t")))}
  }
  
  
  if(multi){
    
    # determine strength of change
    data[, max.dev := max(c(abs(max(as.numeric(y))), abs(min(as.numeric(y) )))), by = "partner"]
    setorder(data, name, -max.dev)
    
    # create new column "label" only use to assign ploting colors
    data[,label := ifelse(max.dev %in% unique(max.dev)[1:n_pars], partner, "Others")]
    
    # Define the plotting colors
    species_colors <- c(
      setNames(dMod_palette(n_pars + 1L)[-1L], unique(data$partner)[1:n_pars]),
      "Others" = "gray"
    )
    
    # Automatically find the corner with the least data density
    legend_corner <- .findEmptyCorner(data$x, data$y)
    
    suppressMessages(
      p <- ggplot2::ggplot(data, aes(x = x, y = y, color = label, group = partner)) + 
        geom_line() +
        xlab(whichPar) + ylab(expression(Delta ~ p[j])) +
        scale_linetype_discrete(name = "profile\nlist") +
        scale_color_manual(values = species_colors) + theme_dMod() +
        theme(legend.position = "inside",
              legend.position.inside = legend_corner,
              legend.justification = legend_corner,
              legend.title = element_blank(),
              legend.background = element_rect(fill = alpha("white", 0.85), colour = "black", linewidth = 0.3),
              legend.key.size = unit(0.4, "cm"),
              legend.margin = margin(2, 4, 2, 4),
              legend.text = element_text(size = 7))
    )
  } else {
    suppressMessages(
      p <- ggplot2::ggplot(data, aes(x = x, y = y, group = interaction(name, proflist), color = name, lty = proflist)) + 
        facet_wrap(~combination, scales = scales) + 
        geom_path() + #geom_point(aes=aes(size=1), alpha=1/3) +
        xlab(axis.labels[1]) + ylab(axis.labels[2]) +
        scale_linetype_discrete(name = "profile\nlist") +
        scale_color_dMod(name = "profiled\nparameter")
    )
  }
  
  attr(p, "data") <- data
  return(p)
  
}

#' Profile Likelihood: Plot All Parameter Paths of One Profile Together
#'
#' Draws, for each profiled parameter, the changes of all other parameters
#' along its profile in one panel. The `npars` parameters with the largest
#' change are colored and named, the others are gray. Requires the package
#' cowplot.
#'
#' @param profs A [parframe] of profiles as returned by [profile()].
#' @param whichpars Character vector, the profiled parameters to draw.
#' @param npars Integer, number of colored and named paths. Defaults to 5.
#' @param normalizePaths Logical. `TRUE` scales each path to a maximum
#'   absolute value of 1. Defaults to `FALSE`.
#'
#' @return A `ggplot` object; for several `whichpars` the grid of
#'   [cowplot::plot_grid()].
#' @author Svenja Kemmer, \email{svenja.kemmer@@fdm.uni-freiburg.de}
#' @seealso [plotPaths()], [plotProfilesAndPaths()]
#' @examplesIf requireNamespace("cowplot", quietly = TRUE)
#' pars <- c(a = 1, b = 0.5, c = 2)
#' obj <- constraintL2(mu = pars, sigma = 0.1)
#' profs <- profile(obj, pars, whichPar = c("a", "b"), limits = c(-1, 1),
#'                  cores = 1)
#' plotPathsMulti(profs, c("a", "b"), npars = 2)
#' @export
#' @import data.table
plotPathsMulti <- function(profs, whichpars, npars = 5, normalizePaths = FALSE) {
  .require_ns("cowplot", "plotPathsMulti()")
  if(length(whichpars) == 1){
    p <- PlotPaths(profs=profs, whichPar = whichpars, n_pars = npars, normalizePaths = normalizePaths)
    return(p)
  } else {
    PlotList <- NULL
    for(i in 1:length(whichpars)){
      par <- whichpars[i]
      p <- PlotPaths(profs=profs, whichPar = par, n_pars = npars, normalizePaths = normalizePaths)
      PlotList[[i]] <- p
    }
    pl <- cowplot::plot_grid(plotlist = PlotList)
    return(pl)
  }
}


#' Profile Likelihood: Plot Profiles Along with Their Parameter Paths
#'
#' For each profiled parameter, draws the profile above the paths of
#' [plotPathsMulti()]. Requires the package cowplot.
#'
#' @param profs A [parframe] of profiles as returned by [profile()].
#' @param whichpars Character vector, the profiled parameters to draw.
#' @param npars Integer, number of colored and named paths. Defaults to 5.
#' @param ncols Number of columns of the plot grid. Defaults to 3.
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
#' plotProfilesAndPaths(profs, c("a", "b"), npars = 2, ncols = 2)
#'
#' @export
plotProfilesAndPaths <- function(profs, whichpars, npars = 5, ncols = 3, normalizePaths = FALSE, modes = c("data", "prior"), ...) {
  .require_ns("cowplot", "plotProfilesAndPaths()")

  # Save original obj.attributes before any subsetting drops them
  orig_oa <- attr(profs, "obj.attributes")
  filtered_oa <- if (!is.null(orig_oa)) intersect(orig_oa, modes) else NULL
  
  profs <- profs[profs$whichPar %in% whichpars]
  
  cleanProfilePlot <- function(prof_sub) {
    # Remove columns for unwanted modes, so plotProfile cannot plot them
    cols_to_drop <- setdiff(orig_oa, modes)
    prof_sub <- prof_sub[, !(colnames(prof_sub) %in% cols_to_drop), drop = FALSE]
    attr(prof_sub, "obj.attributes") <- filtered_oa
    p <- plotProfile(prof_sub)
    
    # plotProfile always adds "total"; filter the underlying data to requested modes
    pdata <- attr(p, "data")
    pdata <- pdata[pdata$mode %in% modes, , drop = FALSE]
    
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
  
  stacked_list <- vector("list", length(whichpars))
  
  for (z in seq_along(whichpars)) {
    prof_sub <- profs[profs$whichPar == whichpars[z]]
    
    p_prof_noleg <- cleanProfilePlot(prof_sub)
    
    p_paths <- plotPathsMulti(prof_sub, whichpars[z], npars, normalizePaths = normalizePaths) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    aligned_pair <- cowplot::align_plots(p_prof_noleg, p_paths, align = "v", axis = "tb")
    stacked_list[[z]] <- cowplot::plot_grid(aligned_pair[[1]], aligned_pair[[2]],
                                            ncol = 1, rel_heights = c(1, 0.7), align = "v", axis = "tb")
  }
  
  body <- cowplot::plot_grid(plotlist = stacked_list, ncol = ncols, ...)
  
  return(body)
}
