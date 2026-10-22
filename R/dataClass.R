
## Class "datalist" and its constructor ------------------------------------------


#' @param x object of class `data.frame` or `list` of data frames. Data frames
#' must provide the columns "name", "time" and "value". The columns "sigma"
#' and "lloq" are optional and default to `NA` and `-Inf`.
#' @return Object of class `datalist`. `is.datalist()` returns a logical.
#' @export
#' @example inst/examples/datalist.R
#' @rdname datalist
as.datalist <- function(x, ...) {
  UseMethod("as.datalist", x)
}

#' @export
#' @param splitBy character, the columns whose combined values identify a
#' condition. Default `NULL` uses all columns except "name", "time",
#' "value", "sigma" and "lloq".
#' @param keepCovariates character, additional columns kept in the
#' condition grid. Default `NULL`.
#' @rdname datalist
as.datalist.data.frame <- function(x, splitBy = NULL, keepCovariates = NULL, ...) {

  x <- sanitizeData(as.data.frame(x))
  dataframe <- x[["data"]]
  standard.names <- x[["columns"]]
  all.names <- colnames(dataframe)

  if (is.null(splitBy)) splitBy <- setdiff(all.names, standard.names)
  # Data without a splitting column is one condition, named as datalist() names it.
  splits <- if (length(splitBy))
    do.call(paste, c(lapply(splitBy, function(n) dataframe[[n]]), list(sep = "_")))
  else rep("1", nrow(dataframe))
  if (length(splitBy) && nrow(unique(dataframe[splitBy])) > length(unique(splits)))
    stop("the values of ", paste(splitBy, collapse = ", "), " pasted with \"_\" give the ",
         "same condition name for different conditions.", call. = FALSE)


  # condition grid
  conditionframe <- dataframe[!duplicated(splits), union(splitBy, keepCovariates), drop = FALSE]
  if (!ncol(conditionframe)) conditionframe <- data.frame(condition = unique(splits))
  rownames(conditionframe) <- splits[!duplicated(splits)]


  # data list output
  dataframe <- dataframe[, standard.names, drop = FALSE]
  out <- lapply(unique(splits), function(s) dataframe[splits == s, , drop = FALSE])

  names(out) <- as.character(unique(splits))

  out <- as.datalist(out)
  attr(out, "condition.grid") <- conditionframe
  return(out)

}

#' @export
#' @param names character, the condition names. Default `NULL` takes them
#' from `names(x)`.
#' @param condition.grid data frame with one row per condition. Default is the
#' `condition.grid` attribute of `x`; if that is `NULL`, a grid with the
#' single column `condition`.
#' @rdname datalist
as.datalist.list <- function(x, names = NULL, ..., condition.grid = attr(x, "condition.grid")) {

  if (!all(vapply(x, is.data.frame, logical(1)))) stop("list of data.frame expected")
  frames <- lapply(x, function(d) sanitizeData(as.data.frame(d))[["data"]])
  conditions <- if (is.null(names)) base::names(x) else names
  if (length(conditions) != length(frames)) stop("names argument has wrong length")

  # A list without a grid gets one with the condition names as its only column.
  grid <- condition.grid %||% data.frame(condition = conditions, row.names = conditions)
  structure(stats::setNames(frames, conditions), class = c("datalist", "list"),
            condition.grid = grid)
}


## Methods for class datalist ---------------------------------------

#' @param value character, the new condition names. The rows of the
#' condition grid are renamed accordingly.
#' @export
#' @rdname datalist
"names<-.datalist" <- function(x, value) {
  grid <- attr(x, "condition.grid")
  rownames(grid) <- value
  if (!is.null(grid[["condition"]])) grid[["condition"]] <- value
  frames <- unclass(x)
  names(frames) <- value
  as.datalist(frames, condition.grid = grid)
}

#' @export
#' @rdname datalist
is.datalist <- function(x) {
  inherits(x, "datalist")
}

#' @export
#' @rdname datalist
c.datalist <- function(...) {
  parts <- list(...)
  grid <- Reduce(combine, lapply(parts, attr, "condition.grid"))
  frames <- do.call(c, lapply(parts, function(d) `attributes<-`(d, list(names = names(d)))))
  repeated <- unique(names(frames)[duplicated(names(frames))])
  if (length(repeated))
    stop("conditions are duplicated: ", paste(repeated, collapse = ", "), call. = FALSE)
  structure(frames, condition.grid = grid, class = c("datalist", "list"))
}


#' @export
print.datalist <- function(x, ...) {
  datalist <- x
  for(n in names(datalist)) {
    cat(n, ":\n", sep = "")
    print(datalist[[n]])
  }
}

# Keeps, in every condition, the rows where `subset` holds; the grid stays as it is.
#' @export
subset.datalist <- function(x, subset, ...) {
  e <- substitute(subset)
  env <- parent.frame()
  out <- unclass(x)
  out[] <- lapply(out, function(d) {
    r <- eval(e, d, env)
    d[r & !is.na(r), , drop = FALSE]
  })
  class(out) <- class(x)
  out
}

#' @export
"[.datalist" <- function(x, ...) {
  condition.grid <- attr(x, "condition.grid")
  out <- unclass(x)[...]
  attr(out, "condition.grid") <- condition.grid
  n <- names(out)
  if (!is.null(n)) {
    attr(out, "condition.grid") <- condition.grid[n, , drop = FALSE]
  }
  class(out) <- c("datalist", "list")
  return(out)
}

#' Plot Observed Data
#'
#' @description
#' Plots data with error bars and marks values below the limit of
#' quantification (BLoQ).
#'
#' @param data object of class `datalist`, or a data frame with the columns
#'   `name`, `time`, `value` and `sigma`.
#' @param ... filter expressions passed to [dplyr::filter()]. They can refer
#'   to the columns of the data and of the condition grid.
#' @param scales the `scales` argument of [ggplot2::facet_wrap()] or
#'   [ggplot2::facet_grid()]: `"free"` (default), `"fixed"`, `"free_x"` or
#'   `"free_y"`.
#' @param facet `"wrap"` (default): one panel per name, colour by condition.
#'   `"grid"`: names as rows, conditions as columns. `"wrap_plain"`: one panel
#'   per combination of name and condition.
#' @param transform list of transformations for the states, see
#'   [coordTransform()]. Default `NULL`.
#'
#' @return A `ggplot` object. Its attribute `"data"` holds the plotted data
#'   frame.
#' @seealso [plotCombined()], [plotPrediction()]
#'
#' @examples
#' data <- datalist(
#'   C1 = data.frame(name = "A", time = 0:5, value = 0:5, sigma = 0.1),
#'   C2 = data.frame(name = "A", time = 0:5, value = sin(0:5), sigma = 0.1)
#' )
#' plotData(data, time < 4)
#' plotData(data, facet = "grid")
#'
#' @export
#' @rdname plotData
#' @importFrom dplyr filter
plotData.datalist <- function(data, ..., scales = "free",
                              facet = c("wrap", "grid", "wrap_plain"),
                              transform = NULL) {

  facet <- match.arg(facet)

  # The condition grid joins the data by condition, so filters can use covariates.
  grid <- covariates(data)
  covs <- grid[setdiff(colnames(grid), "condition")]
  grid <- data.frame(condition = rownames(grid), covs, check.names = FALSE)
  frame <- base::merge(lbind(data), grid, by = "condition", all.x = TRUE)
  frame <- as.data.frame(dplyr::filter(frame, ...))
  frame$bloq <- ifelse(frame$value <= frame$lloq, "yes", "no")
  if (!is.null(transform)) frame <- coordTransform(frame, transform)

  layout <- switch(facet,
    wrap       = list(aes(group = condition, color = condition),
                      facet_wrap(~name, scales = scales)),
    grid       = facet_grid(name ~ condition, scales = scales),
    wrap_plain = facet_wrap(~name * condition, scales = scales))
  noBloq <- if (all(frame$bloq == "no")) guides(shape = "none")

  p <- ggplot(frame, aes(x = time, y = value, ymin = value - sigma, ymax = value + sigma,
                         pch = bloq)) +
    layout +
    geom_point() +
    geom_errorbar(width = 0) +
    scale_shape_manual(name = "BLoQ", values = c(yes = 4, no = 19)) +
    noBloq
  structure(p, data = frame)
}


#' @describeIn plotData S3 plot method for `datalist` objects.
#' @param x object of class `datalist`.
#' @export
plot.datalist <- function(x, ..., scales = "free",
                          facet = c("wrap", "grid", "wrap_plain"),
                          transform = NULL) {
  plotData.datalist(data = x, ..., scales = scales, facet = facet, transform = transform)
}


#' Coerce to a Data Frame
#'
#' @param x any R object
#' @return a data frame
#' @rdname as.data.frame.dMod
#' @export
as.data.frame.datalist <- function(x, ...) {

  data <- x
  condition.grid <- attr(x, "condition.grid")

  data <- lbind(data)
  if (!is.null(condition.grid)) {
    for (C in colnames(condition.grid)) {
      if (nrow(data) > 0)
        data[, C] <- condition.grid[as.character(data$condition), C]
      else
        data[, C] <- vector(mode = mode(condition.grid[[C]]))
    }
  }

  return(data)

}


#' Access the Covariates in the Data
#'
#' @param x a [datalist] or a `data.frame` with the columns
#' `c("name", "time", "value", "sigma", "lloq")`.
#' @param ... not used.
#'
#' @return A data frame with one row per condition: the `condition.grid` of a
#' datalist, or the unique combinations of the non-standard columns of a data
#' frame.
#' @seealso [as.datalist()]
#' @export
covariates <- function(x, ...) {
  UseMethod("covariates", x)
}

#' @export
#' @rdname covariates
covariates.datalist <- function(x, ...) {
  attr(x, "condition.grid", exact = TRUE)
}

#' @export
#' @rdname covariates
covariates.data.frame <- function(x, ...) {

  exclude <- c("name", "time", "value", "sigma", "lloq")
  out <- unique(x[, setdiff(names(x), exclude), drop = FALSE])
  if (!ncol(out)) stop("The data frame has no condition or covariate column.", call. = FALSE)

  if (!"condition" %in% colnames(out)) {
    rownames(out) <- do.call(paste_, out)
    return(out)
  }
  if (anyDuplicated(out[["condition"]]))
    stop("Unique entries of condition column do not correspond to unique covariate combinations.")
  rownames(out) <- out[["condition"]]
  if (ncol(out) > 1L) out <- out[, setdiff(names(out), "condition"), drop = FALSE]
  out

}


## datalist constructor ----------------------------------------

## Data classes ----------------------------------------------------------------

#' Generate a Datalist Object
#'
#' @description A datalist stores time-course data as a named list of data
#' frames, one per condition.
#' @details The standard columns of the data frames are "name" (observable
#' name), "time", "value", "sigma" (uncertainty, can be `NA`) and "lloq"
#' (lower limit of quantification, `-Inf` by default; `NA` means no limit).
#'
#' The attribute `condition.grid` holds one row per condition with further
#' information such as doses; see [covariates()]. Renaming a datalist renames
#' the rows of its condition grid.
#' @param ... for `datalist()` and `c()`, named data frames or datalists; for
#' the `as.datalist()` methods, not used.
#' @return Object of class `datalist`.
#' @seealso [plotData()], [+.datalist], [normL2()]
#' @export
datalist <- function(...) {
  mylist <- .fillNames(list(...), "")
  as.datalist(mylist, names(mylist))
}


