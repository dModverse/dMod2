
# ggplot2 / dplyr NSE column references; declared so R CMD check does not
# flag them as undefined globals.
utils::globalVariables(c("predicted", "observed", "sd_est", "iter", "level"))


# Custom interface to ggplot2 ---

#' Open last plot in external pdf viewer
#' 
#' @description Convenience function to show last plot in an external viewer.
#' @param plot `ggplot2` plot object.
#' @param command character, indicatig which pdf viewer is started.
#' @param ... arguments going to `ggsave`.
#' @export
ggopen <- function(plot = last_plot(), command = "xdg-open", ...) {
  filename <- tempfile(pattern = "Rplot", fileext = ".pdf")
  ggsave(filename = filename, plot = plot, ...)
  system(command = paste(command, filename))
}


#' Standard plotting theme of dMod
#' 
#' @param base_size numeric, font-size
#' @param base_family character, font-name
#' @param showGrid logical, keep the panel grid. `FALSE` (the default) drops it;
#'   `TRUE` leaves [ggplot2::theme_bw()]'s grid in place.
#' @export
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
# Every palette below carries a measured number: the smallest pairwise CIE2000
# distance within it, taken as the WORST case over normal, deuteranopic,
# protanopic and tritanopic vision (for ramps, between positions at least a
# quarter of the domain apart). A palette counts as colorblind-safe at 10 or
# above, below roughly 3 two colors are indistinguishable side by side, and a
# thin line needs more headroom than a filled patch.
#
# The numbers are baked in rather than computed at load time so that the
# simulation packages stay out of the dependency list. Re-derive them with
# colorspace::deutan/protan/tritan and farver::compare_colour(method = "cie2000").

#' Seed colors of the dMod palette
#'
#' The ten qualitative house colors. Not colorblind-safe, its brown and red
#' collapse under protanopia, see [dMod_palettes()] for the alternatives.
#'
#' @export
dMod_colors <- c("#000000", "#C5000B", "#0084D1", "#579D1C", "#FF950E",
                 "#4B1F6F", "#CC79A7", "#006400", "#F0E442", "#8B4513")

#' Colorblind-safe seed colors
#'
#' The Okabe-Ito palette, extended to twenty. The first eight are Okabe-Ito
#' itself; the rest were chosen greedily to maximise the worst-case pairwise
#' distance under simulated dichromacy, subject to every color also keeping a
#' distance of 25 to the white background. The worst case stays at Okabe-Ito's
#' own 11.1 for all n up to twenty, so the twelve extra colors cost nothing.
#'
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

#' Continuous dMod colour ramps
#'
#' `dMod_gradient` is the sequential house ramp, `dMod_divergent` the diverging
#' one; `dMod_gradient_cb` and `dMod_divergent_cb` are their colorblind-safe
#' counterparts and the defaults of [scale_color_dMod_c()] and
#' [scale_color_dMod_div()]. See [dMod_palettes()] for the full set.
#'
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

#' Overview of the dMod palettes
#'
#' One row per palette, with the measured worst-case CIE2000 distance under
#' simulated dichromacy and the resulting colorblind verdict. Pass any of the
#' names as `palette =` to [dMod_palette()], [scale_color_dMod()],
#' [scale_color_dMod_c()] or [scale_color_dMod_div()].
#'
#' @param type restrict to `"qualitative"`, `"sequential"` or `"diverging"`.
#' @return A data frame with columns `name`, `type`, `n`, `colorblind`, `dE`.
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

#' Generate `n` distinct colors from a dMod palette
#'
#' Returns the first `n` colors of the chosen palette. For `n` beyond it,
#' `Polychrome::createPalette()` extends the seeds deterministically; if
#' Polychrome is not installed the overflow comes from [grDevices::hcl.colors()]
#' instead. Neither extension inherits the palette's colorblind property, which
#' is why a colorblind-safe choice warns once it runs out of seeds.
#'
#' @param n integer, number of colors to produce.
#' @param palette name of a qualitative palette, see [dMod_palettes()].
#' @return Character vector of length `n` with hex color codes.
#' @export
dMod_palette <- function(n, palette = "okabe") {
  n <- as.integer(n)
  if (n <= 0L) return(character(0))
  pal <- .dMod_pick(palette, "qualitative")
  seeds <- pal$colors
  if (n <= length(seeds)) return(unname(seeds[seq_len(n)]))
  if (pal$dE >= 10)
    warning("dMod_palette(): palette \"", palette, "\" holds ", length(seeds),
            " colorblind-safe colors, ", n, " were requested. The extra ",
            n - length(seeds), " are not -- encode them with linetype, shape ",
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

#' Discrete dMod colour scales
#'
#' @param palette name of a qualitative palette, see [dMod_palettes()]. The
#'   default is colorblind-safe.
#' @param ... arguments forwarded to [ggplot2::discrete_scale()].
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

#' Continuous dMod colour scales
#'
#' Sequential (`_c`) and diverging (`_div`) counterparts of the discrete
#' [scale_color_dMod()]. The diverging scales pin white to `mid` even for
#' asymmetric limits, so a signed quantity keeps its zero at the neutral colour.
#'
#' The suffix is `_div`, not `_d`: ggplot2 spells *discrete* `_d` (as in
#' [ggplot2::scale_colour_viridis_d()]), and here the discrete scale is the
#' unsuffixed [scale_color_dMod()].
#'
#' @param mid numeric, the value white is pinned to.
#' @param direction 1 or -1; -1 reverses the ramp.
#' @param palette name of a sequential palette for `_c`, a diverging one for
#'   `_div`, see [dMod_palettes()]. Both defaults are colorblind-safe.
#' @param ... arguments forwarded to [ggplot2::scale_colour_gradientn()] or
#'   [ggplot2::scale_fill_gradientn()].
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

#' Coordinate transformation for data frames
#' 
#' Applies a symbolically defined transformation to the `value`
#' column of a data frame. Additionally, if a `sigma` column is
#' present, those values are transformed according to Gaussian error
#' propagation.
#' @param data data frame with at least columns "name" (character) and
#' "value" (numeric). Can optionally contain a column "sigma" (numeric).
#' @param transformations character (the transformation) or named list of
#' characters. In this case, the list names must be a subset of those 
#' contained in the "name" column.
#' @return The data frame with the transformed values and sigma uncertainties.
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


#' Plot a list of model predictions
#' 
#' @param prediction Named list of matrices or data.frames, usually the output of a prediction function
#' as generated by [Xs].
#' @param ... Further arguments going to `dplyr::filter`. 
#' @param scales The scales argument of `facet_wrap` or `facet_grid`, i.e. `"free"`, `"fixed"`, 
#' `"free_x"` or `"free_y"`
#' @param facet Either `"wrap"` or `"grid"`
#' @param transform list of transformation for the states, see [coordTransform].
#' @details The data.frame being plotted has columns `time`, `value`, `name` and `condition`.
#'  
#' 
#' @return A plot object of class `ggplot`.
#' @import ggplot2
#' @example inst/examples/plotting.R
#' @export
plotPrediction <- function(prediction,...) {
  UseMethod("plotPrediction", prediction)
}


#' Plot a list of model predictions and a list of data points in a combined plot
#' 
#' @param prediction Named list of matrices or data.frames, usually the output of a prediction function
#' as generated by [Xs].
#' @param data Named list of data.frames as being used in [res], i.e. with columns `name`, `time`, 
#' `value` and `sigma`.
#' @param ... Further arguments going to `dplyr::filter`. 
#' @param scales The scales argument of `facet_wrap` or `facet_grid`, i.e. `"free"`, `"fixed"`, 
#' `"free_x"` or `"free_y"`
#' @param facet `"wrap"` or `"grid"`. Try `"wrap_plain"` for high amounts of conditions and low amounts of observables.
#' @param transform list of transformation for the states, see [coordTransform].
#' @param aesthetics Named list of aesthetic mappings, specified as character, e.g. `list(linetype = "name")`. 
#' Can refer to variables in the condition.grid
#' @details The data.frame being plotted has columns `time`, `value`, `sigma`,
#' `name` and `condition`.
#'  
#' 
#' @return A plot object of class `ggplot`.
#' @example inst/examples/plotting.R
#' @importFrom graphics par
#' @export
plotCombined <- function(prediction,...) {
  UseMethod("plotCombined", prediction)
}


#' Plot a list data points
#' 
#' @param data Named list of data.frames as being used in [res], i.e. with columns `name`, `time`, 
#' `value` and `sigma`.
#' @param ... Further arguments going to `subset`. 
#' @param scales The scales argument of `facet_wrap` or `facet_grid`, i.e. `"free"`, `"fixed"`, 
#' `"free_x"` or `"free_y"`
#' @param facet Either `"wrap"` or `"grid"`
#' @param transform list of transformation for the states, see [coordTransform].
#' @details The data.frame being plotted has columns `time`, `value`, `sigma`,
#' `name` and `condition`.
#'  
#' 
#' @return A plot object of class `ggplot`.
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

#' Profile likelihood plot
#' 
#' @param profs Lists of profiles as being returned by [profile].
#' @param ... logical going to subset before plotting.
#' @param maxvalue Numeric, the value where profiles are cut off.
#' @param parlist Matrix or data.frame with columns for the parameters to be added to the plot as points.
#' If a "value" column is contained, deltas are calculated with respect to lowest chisquare of profiles.
#' @param ncol Number of columns in the resulting plot grid.
#' @param threshold Numeric, the horizontal lines and the y-axis breaks on top
#' of the optimum the profiles start from, which is always drawn. Defaults to the
#' chi-square thresholds for 68%, 90% and 95%. Names are used as axis labels.
#' Pass the value [profileThreshold] returns for the calibration the intervals
#' are read at, so the line and [confint.parframe] agree.
#' @return A plot object of class `ggplot`.
#' @details See [profile] for examples.
#' @export
plotProfile <- function(profs,...) {
  UseMethod("plotProfile", profs)
}


#' Profile likelihood: plot of the parameter paths.
#' 
#' @param profs profile or list of profiles as being returned by [profile]
#' @param ... arguments going to subset
#' @param whichPar Character or index vector, indicating the parameters that are taken as possible reference (x-axis)
#' @param sort Logical. If paths from different parameter profiles are plotted together, possible
#' combinations are either sorted or all combinations are taken as they are.
#' @param relative logical indicating whether the origin should be shifted.
#' @param scales character, either `"free"` or `"fixed"`.
#' @return A plot object of class `ggplot`.
#' @details See [profile] for examples.
#' @export
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


#' Plot Fluxes given a list of flux Equations
#'
#' @param pouter parameters
#' @param x The model prediction function `x(times, pouter, fixed, ...)`
#' @param fluxEquations list of chars containing expressions for the fluxes,
#' if names are given, they are shown in the legend. Easy to obtain via [subset.eqnlist], see Examples.
#' @param nameFlux character, name of the legend.
#' @param times Numeric vector of time points for the model prediction
#' @param ... Further arguments going to x, such as `fixed` or `conditions`
#'
#'
#' @return A plot object of class `ggplot`.
#' @examples
#' \dontrun{
#'
#' plotFluxes(bestfit, x, times, subset(f, "B"%in%Product)$rates, nameFlux = "B production")
#' }
#' @export
plotFluxes <- function(pouter, x, times, fluxEquations, nameFlux = "Fluxes:", ...){

  if (is.null(names(fluxEquations))) names(fluxEquations) <- fluxEquations

  flux <- funCpp(fluxEquations, convenient = FALSE)$func
  prediction.all <- x(times, pouter, deriv = FALSE, ...)
  names.prediction.all <- names(prediction.all)
  if (is.null(names.prediction.all)) names.prediction.all <- paste0("C", 1:length(prediction.all))

  out <- lapply(1:length(prediction.all), function(cond) {
    prediction <- prediction.all[[cond]]
    pinner <- attr(prediction, "parameters")
    pinner.matrix <- matrix(pinner, nrow = length(pinner), ncol = nrow(prediction),
                            dimnames = list(names(pinner), NULL))
    fluxes <- cbind(time = prediction[, "time"], flux(cbind(prediction, t(pinner.matrix))))
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

#' Plotting objective values of a collection of fits
#' 
#' @param x data.frame with columns "value", "converged" and "iterations", e.g. 
#' a [parframe].
#' @param ... arguments for subsetting of x
#' @param tol maximal allowed difference between neighboring objective values
#' to be recognized as one.
#' @param showSteps logical, if `TRUE`, the detected steps are indicated by
#' dashed vertical lines and labelled by their index. Defaults to `FALSE`.
#' @export
plotValues <- function(x,...) {
  UseMethod("plotValues", x)
}


#' Plot parameter values for a fitlist
#' 
#' @param x parameter frame as obtained by as.parframe(mstrust)
#' @param tol maximal allowed difference between neighboring objective values
#' to be recognized as one.
#' @param ... arguments for subsetting of x
#' @export
plotPars <- function(x,...) {
  UseMethod("plotPars", x)
}


#' Plot residuals for a fitlist
#'
#' @description
#' Creates a plot of residuals from model fits, with flexible options for
#' grouping and faceting. Residuals can be summarized across different
#' dimensions (time, condition, observable, fit index).
#'
#' @param parframe Object of class \code{parframe}, e.g. returned by \link{mstrust}.
#' @param x Prediction function returning named list of data.frames with names 
#'   matching \code{data}.
#' @param data A \code{datalist} object, i.e. named list of data.frames with 
#'   columns \code{name}, \code{time}, \code{value}, and \code{sigma}.
#' @param split Character vector specifying how to summarize and display residuals.
#'   \itemize{
#'     \item \code{split[1]}: Variable for x-axis
#'     \item \code{split[2]}: Variable for grouping (color/line), defaults to \code{split[1]}
#'     \item \code{split[3+]}: Additional variables for \code{facet_wrap()}
#'   }
#' @param errmodel Optional error model function of type \code{prdfn}. If provided,
#'   residuals include the log-likelihood contribution from sigma.
#' @param ... Additional arguments passed to the prediction function \code{x}.
#'
#' @return A \code{ggplot} object with the summarized residual data frame
#'   attached as attribute \code{"out"}.
#'
#' @examples
#' \dontrun{
#' # Time on x-axis, faceted by condition and name
#' plotResiduals(myfitlist, g * x * p, data, 
#'               c("time", "index", "condition", "name"), 
#'               conditions = myconditions[1:4])
#'
#' # Condition on x-axis, residuals summed over time
#' plotResiduals(myfitlist, g * x * p, data, c("condition", "name", "index"))
#' }
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


#' Convergence trace of an iterative fit
#'
#' @description Generic. Methods plot the quantities their fitting method
#'   iterates on, one panel each.
#' @param x Object to plot.
#' @param ... Method-specific arguments.
#' @return A ggplot.
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


#' Pair (corner-style) plot for `mcmc()` outputs
#'
#' Lower triangle: 2D scatter of post-warmup samples. Diagonal: 1D
#' marginal density.
#'
#' @param x An `mcmcResult` (or subclass) object.
#' @param ... Method-specific arguments.
#' @return A ggplot.
#' @export
plotPairs <- function(x, ...) UseMethod("plotPairs", x)




## ---- profile / parameter-path plotting (moved from toolsSvenja.R) ---------
#' Plot an array of trajectories along the profile of a parameter
#' 
#' @param par Character of parameter name for which the array should be generated.
#' @param profs Lists of profiles as being returned by [profile]. 
#' @param prd Named list of matrices or data.frames, usually the output of a prediction function
#' as generated by [Xs].
#' @param times Numeric vector of time points for the model prediction.
#' @param direction Character "up" or "down" indicating the direction the value should be traced along the profile starting at the bestfit value.
#' @param covtable Optional covariate table or condition.grid necessary if subsetting is required.
#' @param ... Further arguments for subsetting the plot.
#' @param nsimus Number of trajectories/ simulation to be calculated.
#' 
#' @return A plot object of class `ggplot`.
#' @author Svenja Kemmer, \email{svenja.kemmer@@fdm.uni-freiburg.de}
#' @examples
#' \dontrun{
#'  plotArray("myparameter", myprofiles, g*x*p, seq(0, 250, 1), 
#'     "up", condition.grid, name == "ProteinA" & condition == "c1") 
#' }
#' @export
#' @import data.table
plotArray <- function (par, profs, prd, times, direction = c("up", "down"), covtable, ..., nsimus = 4) {
  
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
    geom_line(size = 1) + 
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

#' @keywords internal
#' @importFrom ggplot2 ggplot
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

#' Profile likelihood: plot all parameter paths belonging to one profile in one plot
#' 
#' @param profs Lists of profiles as being returned by [profile]. 
#' @param whichpars Character vector of parameter names for which the profile paths should be generated.
#' @param npars Numeric vector of number of colored and named parameter paths.
#' @param normalizePaths Logical indicating whether the paths should be normalized to absolute values of 1. Default `FALSE`; `TRUE` only useful in corner cases when you know why to do so.
#' 
#' @return A plot object of class `ggplot` for length(whichpars) = 1 and otherwise an object of class `cowplot`.
#' @author Svenja Kemmer, \email{svenja.kemmer@@fdm.uni-freiburg.de}
#' @examples
#' \dontrun{
#'  plotPathsMulti(myprofiles, c("mypar1", "mypar2"), npars = 5) 
#' }
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


#' Profile likelihood: plot profiles along with their parameter paths
#' 
#' Generates combined plots of profile likelihoods and their parameter paths without a shared legend.
#' 
#' @param profs List of profiles as returned by [profile()].
#' @param whichpars Character vector of parameter names for which the profile paths should be generated.
#' @param npars Numeric indicating number of colored and named parameter paths.
#' @param ncols Number of columns in the resulting plot grid.
#' @param normalizePaths Logical indicating whether the paths should be normalized to absolute values of 1.
#'                       Default `FALSE`.
#' @param modes Character vector of profile modes to display in the profile plot.
#'              Default `c("data", "prior")`. Use e.g. `"data"` to show only the data contribution.
#' @param ... Additional arguments passed to `cowplot::plot_grid()`.
#' 
#' @return A combined `ggplot` object containing the profiles and paths (no shared legend).
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
