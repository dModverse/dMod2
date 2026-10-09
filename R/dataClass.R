
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

  .renameArgs(list(...), c(split.by = "splitBy",
                           keep.covariates = "keepCovariates"), "as.datalist")

  # Sanitize data and get names
  x <- sanitizeData(x)
  dataframe <- x[["data"]]
  standard.names <- x[["columns"]]
  all.names <- colnames(dataframe)
  
  # Get splitting information
  if (is.null(splitBy)) splitBy <- setdiff(all.names, standard.names)
  conditions <- lapply(splitBy, function(n) dataframe[, n])
  splits <- do.call(paste, c(conditions, list(sep = "_")))


  # condition grid
  conditionframe <- dataframe[!duplicated(splits), union(splitBy, keepCovariates), drop = FALSE]
  rownames(conditionframe) <- splits[!duplicated(splits)]


  # data list output
  dataframe <- cbind(data.frame(condition = splits), dataframe[, standard.names])
  out <- lapply(unique(splits), function(s) dataframe[dataframe[, 1] == s, -1])

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

  mylist <- x


  ## Check properties
  if (is.null(names)) mynames <- names(mylist) else mynames <- names
  is.data.frame <- sapply(mylist, class) == "data.frame"
  if (!all(is.data.frame)) stop("list of data.frame expected")

  # Sanitize data in list
  mylist <- lapply(mylist, function(x) sanitizeData(x)[["data"]])

  if (length(mynames) != length(mylist)) stop("names argument has wrong length")

  ## Prepare output
  names(mylist) <- mynames
  class(mylist) <- c("datalist", "list")


  if (is.null(condition.grid)) {
    condition.grid <- data.frame(condition = mynames, row.names = mynames)
  }
  attr(mylist, "condition.grid") <- condition.grid

  return(mylist)

}


## Methods for class datalist ---------------------------------------

#' @param value character, the new condition names. The rows of the
#' condition grid are renamed accordingly.
#' @export
#' @rdname datalist
"names<-.datalist" <- function(x, value) {
  x <- unclass(x)
  x <- base::`names<-`(x, value)
  attr(x, "condition.grid") <- base::`rownames<-`(attr(x, "condition.grid"), value)
  return(as.datalist(x))
}

#' @export
#' @rdname datalist
is.datalist <- function(x) {
  inherits(x, "datalist")
}

#' @export
#' @rdname datalist
c.datalist <- function(...) {
  dlist <- lapply(list(...), unclass)
  
  condition.grids <- lapply(dlist, function(i) attr(i, "condition.grid"))
  mycg <- Reduce(dMod2::combine, condition.grids)
  
  dlist <- Reduce(c, lapply(dlist, function(i) {`attr<-`(i, "condition.grid", NULL)}))
  attr(dlist, "condition.grid") <-  mycg
  class(dlist) <- "datalist"
  return(dlist)
}


#' @export
print.datalist <- function(x, ...) {
  datalist <- x
  for(n in names(datalist)) {
    cat(n, ":\n", sep = "")
    print(datalist[[n]])
  }
}

# Subset of all datalist entries
#' @export
subset.datalist <- function(x, ...){
  datalist <- lapply(x, function(i) subset(i, ...))
  return(as.datalist(datalist))
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
  
  # --- Prepare covariate table ---
  covtable <- covariates(data)
  covtable <- cbind(condition = rownames(covtable), covtable)
  covtable <- covtable[!duplicated(names(covtable))]
  
  # --- Prepare data ---
  data <- lbind(data)
  data <- base::merge(data, covtable, by = "condition", all.x = TRUE)
  data <- as.data.frame(dplyr::filter(data, ...))
  data$bloq <- ifelse(data$value <= data$lloq, "yes", "no")
  
  if (!is.null(transform)) data <- coordTransform(data, transform)
  
  # --- Construct plot ---
  p <- ggplot(data, aes(x = time, y = value, 
                        ymin = value - sigma, ymax = value + sigma, 
                        pch = bloq))
  
  if (facet == "wrap") {
    p <- p + 
      aes(group = condition, color = condition) +
      facet_wrap(~name, scales = scales)
  } else if (facet == "grid") {
    p <- p + facet_grid(name ~ condition, scales = scales)
  } else {
    p <- p + facet_wrap(~name * condition, scales = scales)
  }
  
  p <- p + 
    geom_point() + 
    geom_errorbar(width = 0) +
    scale_shape_manual(name = "BLoQ", values = c(yes = 4, no = 19))
  
  if (all(data$bloq == "no")) {
    p <- p + guides(shape = "none")
  }
  
  attr(p, "data") <- data
  p
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

  attr(x, "condition.grid")

}

#' @export
#' @rdname covariates
covariates.data.frame <- function(x, ...) {

  exclude <- c("name", "time", "value", "sigma", "lloq")
  contains.condition <- "condition" %in% colnames(x)
  out <- unique(x[, setdiff(names(x), exclude)])

  if (contains.condition) {
    if (any(duplicated(out[["condition"]]))) {
      stop("Unique entries of condition column do not correspond to unique covariate combinations.")
    }
    rownames(out) <- out[["condition"]]
    out <- out[ , setdiff(names(out), "condition")]
  } else {
    rownames(out) <- do.call(paste_, out)
  }

  return(out)

}


## datalist constructor ----------------------------------------

## Data classes ----------------------------------------------------------------

#' Generate a Datalist Object
#'
#' @description A datalist stores time-course data as a named list of data
#' frames, one per condition.
#' @details The standard columns of the data frames are "name" (observable
#' name), "time", "value", "sigma" (uncertainty, can be `NA`) and "lloq"
#' (lower limit of quantification, `-Inf` by default).
#'
#' The attribute `condition.grid` holds one row per condition with further
#' information such as doses; see [covariates()]. Renaming a datalist renames
#' the rows of its condition grid.
#' @param ... for `datalist()` and `c()`, named data frames or datalists; for
#' the `as.datalist()` methods, not used; `split.by` and `keep.covariates`
#' are deprecated names of `splitBy` and `keepCovariates`.
#' @return Object of class `datalist`.
#' @seealso [plotData()], [+.datalist], [normL2()]
#' @export
datalist <- function(...) {
  mylist <- list(...)
  mynames <- names(mylist)
  if (is.null(mynames)) mynames <- as.character(1:length(mylist))
  as.datalist(mylist, mynames)
}


