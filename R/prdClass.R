
## Methods for class prdlist ------------------------------------------------



#' @export
#' @rdname prdlist
as.prdlist <- function(x, ...) {
  UseMethod("as.prdlist", x)
}

#' @export
#' @param x list of prediction frames. Default `NULL` gives an empty list.
#' @param names character, the list names, e.g. the condition names. Default
#'   `NULL` keeps `names(x)`.
#' @rdname prdlist
as.prdlist.list <- function(x = NULL, names = NULL, ...) {

  if (is.null(x)) x <- list()
  if (is.null(names)) mynames <- names(x) else mynames <- names 

  # if (length(mynames) != length(x)) stop("names argument has wrong length")

  ## Prepare output
  names(x) <- mynames
  class(x) <- c("prdlist", "list")

  return(x)

}


#' @export
c.prdlist <- function(...) {
  
  mylist <- list(...)
  mylist <- lapply(mylist, unclass)
  newlist <- do.call(c, mylist)
  
  as.prdlist(newlist)
  
}

#' @export
"[.prdlist" <- function(x, ...) {
  out <- unclass(x)[...]
  class(out) <- c("prdlist", "list")
  return(out)
}







#' @export
print.prdlist <- function(x, ...) {
  
  mynames <- names(x)
  if (is.null(mynames)) mynames <- rep("NULL", length(x))
  
  for (i in 1:length(x)) {
    cat(mynames[i], ":\n", sep = "")
    print(x[[i]])
  }
  
}


#' @export
#' @param data object of class `datalist`. Default `NULL`. Its condition grid
#'   adds covariate columns to the output.
#' @param errfn object of class `obsfn`, the error model that predicts the
#'   column `sigma`. Default `NULL` sets `sigma` to `NaN`.
#' @param ... not used.
#' @rdname as.data.frame.dMod
as.data.frame.prdlist <- function(x, ..., data = NULL, errfn = NULL) {
  
  prediction <- x
  sigma <- NULL
  condition.grid <- attr(data, "condition.grid")
  
  if (!is.null(errfn)) {
    sigma <- as.prdlist(
      lapply(1:length(prediction), 
             function(i) errfn(prediction[[i]], 
                               getParameters(prediction[[i]]), 
                               conditions = names(prediction)[i])[[1]]),
      names = names(prediction)
    )
    sigma <- wide2long(sigma)
  }
  
  prediction <- wide2long(prediction)
  prediction$sigma <- NaN
  if (!is.null(sigma)) {
    # Keyed, not positional: the error model may order or cover its observables
    # differently from the prediction, and a row-wise copy would then hand a
    # sigma to the wrong observable.
    i <- match(paste(prediction$condition, prediction$time, prediction$name),
               paste(sigma$condition, sigma$time, sigma$name))
    prediction$sigma[!is.na(i)] <- sigma$value[i[!is.na(i)]]
  }
  
  if (!is.null(condition.grid)) {
    for (C in colnames(condition.grid)) {
      rows <- ifelse(is.na(prediction$condition), 1, as.character(prediction$condition))
      prediction[, C] <- condition.grid[rows, C]
    }
    n1 <- nrow(prediction)
  }
  
  
  return(prediction)
  
  
} 

#' @export
#' @param x object of class `prdlist` or `prdframe`.
#' @rdname plotCombined
plot.prdlist <- function(x, data = NULL, ..., scales = "free", facet = c("wrap", "grid", "wrap_plain"), transform = NULL) {
  
  prediction <- x
  
  if (is.null(names(prediction))) names(prediction) <- paste0("C", 1:length(prediction))
  if (!is.null(data) && is.null(names(data))) names(data) <- paste0("C", 1:length(data))
  
  plotCombined.prdlist(prediction = prediction, data = data, ..., scales = scales, facet = facet, transform = transform)
  
}


#' @export
#' @rdname plotCombined
#' @importFrom dplyr filter
#' @importFrom rlang parse_expr
plotCombined.prdlist <- function(prediction, data = NULL, ..., 
                                 scales = "free", 
                                 facet = c("wrap", "grid", "wrap_plain"), 
                                 transform = NULL, aesthetics = NULL) {
  
  facet <- match.arg(facet)
  
  make_aes <- function(mapping) {
    mapping <- lapply(mapping, function(col) {
      if (!is.null(col)) rlang::parse_expr(col)
    })
    do.call(aes, Filter(Negate(is.null), mapping))
  }
  
  mynames <- c("time", "name", "value", "sigma", "condition")
  covtable <- NULL
  
  # --- Prepare data ---
  if (!is.null(data)) {
    covtable <- covariates(data)
    covtable <- cbind(condition = rownames(covtable), covtable)
    covtable <- covtable[!duplicated(names(covtable))]
    
    data <- lbind(data)
    data <- base::merge(data, covtable, by = "condition", all.x = TRUE)
    data <- as.data.frame(dplyr::filter(data, ...))
    data$bloq <- ifelse(data$value <= data$lloq, "yes", "no")
    
    if (!is.null(transform)) data <- coordTransform(data, transform)
  }
  
  # --- Prepare prediction ---
  if (!is.null(prediction)) {
    prediction <- cbind(wide2long(prediction), sigma = NA)
    if (!is.null(covtable)) {
      prediction <- base::merge(prediction, covtable, by = "condition", all.x = TRUE)
    }
    prediction <- as.data.frame(dplyr::filter(prediction, ...))
    
    if (!is.null(transform)) prediction <- coordTransform(prediction, transform)
  }
  
  # --- Combine into single data frame ---
  keep_cols <- unique(c(mynames, names(covtable)))
  total <- rbind(
    if (!is.null(prediction)) prediction[, keep_cols] else NULL,
    if (!is.null(data)) data[, keep_cols] else NULL
  )
  
  # --- Build aesthetics ---
  aes_base <- list(x = "time", y = "value", 
                   ymin = "value - sigma", ymax = "value + sigma")
  if (facet == "wrap") {
    aes_base$group <- "condition"
    aes_base$color <- "condition"
  }
  aesthetics <- c(aes_base[setdiff(names(aes_base), names(aesthetics))], aesthetics)
  
  # --- Construct plot ---
  p <- ggplot(total, make_aes(aesthetics))
  
  p <- p + switch(facet,
                  wrap       = facet_wrap(~name, scales = scales),
                  grid       = facet_grid(name ~ condition, scales = scales),
                  wrap_plain = facet_wrap(~name * condition, scales = scales)
  )
  
  if (!is.null(prediction)) {
    p <- p + geom_line(data = prediction)
  }
  
  if (!is.null(data)) {
    p <- p + 
      geom_point(data = data, aes(pch = bloq)) + 
      geom_errorbar(data = data, width = 0) +
      scale_shape_manual(name = "BLoQ", values = c(yes = 4, no = 19))
    
    if (all(data$bloq == "no")) {
      p <- p + guides(shape = "none")
    }
  }
  
  attr(p, "data") <- list(data = data, prediction = prediction)
  p
}


#' @param errfn object of class `obsfn`, an error model. Default `NULL`. If
#'   given, the predicted `sigma` is drawn as a band around the prediction.
#' @export
#' @rdname plotPrediction
#' @importFrom dplyr filter
plotPrediction.prdlist <- function(prediction, ..., errfn = NULL, 
                                   scales = "free", 
                                   facet = c("wrap", "grid"), 
                                   transform = NULL) {
  
  facet <- match.arg(facet)
  
  prediction <- as.data.frame(prediction, errfn = errfn)
  prediction <- dplyr::filter(prediction, ...)
  
  if (!is.null(transform)) prediction <- coordTransform(prediction, transform)
  
  # --- Construct plot ---
  p <- ggplot(prediction, aes(x = time, y = value))
  
  if (facet == "wrap") {
    p <- p + 
      aes(group = condition, color = condition) +
      facet_wrap(~name, scales = scales)
  } else {
    p <- p + facet_grid(name ~ condition, scales = scales)
  }
  
  if (!is.null(errfn)) {
    p <- p + geom_ribbon(
      aes(ymin = value - sigma, ymax = value + sigma, fill = condition), 
      lty = 0, alpha = 0.3
    )
  }
  
  p <- p + geom_line()
  
  attr(p, "data") <- prediction
  p
}



## Methods for class prdframe ----------------------------
#' @export
#' @rdname plotCombined
plot.prdframe <- function(x, data = NULL, ..., scales = "free", facet = c("wrap", "grid", "wrap_plain"), transform = NULL) {
  
  prediction <- x
  
  prediction <- list("C1" = prediction)
  if (!is.null(data) && is.data.frame(data))
    data <- list("C1" = data)
  
  
  plotCombined.prdlist(prediction = prediction, data = data, ..., scales = scales, facet = facet, transform = transform)
  
}

#' @export
print.prdframe <- function(x, ...) {

  d1 <- attr(x, "deriv")
  d2 <- attr(x, "deriv2")

  derivs <- if (!is.null(d1)) {
    sprintf("yes [%s]", paste(dim(d1), collapse = " x "))
  } else "no"
  derivs2 <- if (!is.null(d2)) {
    sprintf("yes [%s]", paste(dim(d2), collapse = " x "))
  } else "no"

  attr(x, "deriv")      <- NULL
  attr(x, "deriv2")     <- NULL
  attr(x, "parameters") <- NULL

  print(unclass(x))
  cat("\n")
  cat("The prediction contains 1st-order derivatives: ", derivs,  "\n", sep = "")
  cat("The prediction contains 2nd-order derivatives: ", derivs2, "\n", sep = "")

}


## Methods for class prdfn ----------------------------------

#' @export
print.prdfn <- function(x, ...) {
  
  conditions <- attr(x, "conditions")
  parameters <- attr(x, "parameters")
  mappings <- attr(x, "mappings")
  
  cat("Prediction function:\n")
  str(args(x))
  cat("\n")
  cat("... conditions:", paste0(conditions, collapse = ", "), "\n")
  cat("... parameters:", paste0(parameters, collapse = ", "), "\n")
 
}

#' @export
summary.prdfn <- function(object,...) {
  
  x <- object
  
  conditions <- attr(x, "conditions")
  parameters <- attr(x, "parameters")
  mappings <- attr(x, "mappings")
  
  cat("Details:\n")
  if (!inherits(x, "composed")) {
    
    output <- lapply(1:length(mappings), function(C) {
      
      list(
        equations = attr(mappings[[C]], "equations"),
        events = .kernelSetting(mappings[[C]], "events"),
        forcings = .kernelSetting(mappings[[C]], "forcings"),
        parameters = attr(mappings[[C]], "parameters")
      )
      
    })
    names(output) <- conditions
    
    #print(output, ...)
    output
    
  } else {
    
    cat("\nObject is composed. See original objects for more details.\n")
    
  }
}

#' @export
print.obsfn <- function(x, ...) {
  
  conditions <- attr(x, "conditions")
  parameters <- attr(x, "parameters")
  mappings <- attr(x, "mappings")
  
  cat("Observation function:\n")
  str(args(x))
  cat("\n")
  cat("... conditions:", paste0(conditions, collapse = ", "), "\n")
  cat("... parameters:", paste0(parameters, collapse = ", "), "\n")
 
}

#' @export
summary.obsfn <- function(object, ...) {
  
  x <- object
  
  conditions <- attr(x, "conditions")
  parameters <- attr(x, "parameters")
  mappings <- attr(x, "mappings")
  
  cat("Details:\n")
  if (!inherits(x, "composed")) {
    
    output <- lapply(1:length(mappings), function(C) {
      
      list(
        equations = attr(mappings[[C]], "equations"),
        states = attr(mappings[[C]], "states"),
        parameters = attr(mappings[[C]], "parameters")
      )
      
    })
    names(output) <- conditions
    
    #print(output, ...)
    output
    
  } else {
    
    cat("\nObject is composed. See original objects for more details.\n")
    
  }
}


# One prediction per parameter row. Rows x conditions go out as a single
# request when the chain supports it, a parframe from mstrust() has many
# rows and the conditions repeat across them.
.predictRows <- function(x, times, pars, dots) {
  rows <- seq_len(nrow(pars))
  parsList <- lapply(rows, function(i) as.parvec(pars, i))
  serial <- function() lapply(parsList, function(pp)
    do.call(x, c(list(times, pp), list(deriv = FALSE), dots)))

  conds <- if (!is.null(dots$conditions)) dots$conditions else attr(x, "conditions")
  extra <- setdiff(names(dots), "conditions")
  if (length(rows) < 2L || is.null(conds) || length(extra) ||
      is.null(.fnNode(x))) return(serial())

  nc <- length(conds)
  flat <- .predictMany(x, times,
                       parsList = rep(parsList, each = nc),
                       conditions = rep(conds, times = length(rows)),
                       deriv = FALSE)
  lapply(rows, function(i)
    as.prdlist(setNames(flat[seq_len(nc) + (i - 1L) * nc], conds)))
}



#' Model Predictions
#' 
#' Evaluates a prediction function for every row of a parameter frame.
#'
#' @param object object of class `prdfn`.
#' @param ... further arguments going to the prediction function, e.g.
#'   `conditions`.
#' @param times numeric vector of time points.
#' @param pars a [parframe], e.g. `as.parframe(mstrust(...))` or the result
#'   of [profile()].
#' @param data object of class `datalist`. Default `NULL`. Its condition grid
#'   adds covariate columns to the output.
#' @return A data frame in long format with one block per row of `pars`. The
#'   non-parameter columns of `pars` are added with a leading dot, e.g.
#'   `.value`. The data, as a data frame, is returned as the attribute
#'   `"data"`.
#' @export
predict.prdfn <- function(object, ..., times, pars, data = NULL) {
  
  
  x <- object
  arglist <- list(...)
  if (any(names(arglist) == "conditions")) {
    C <- arglist[["conditions"]]
    if (!is.null(data)) {
      data <- data[C]
    }
  }
  if (is.null(data)) data <- data.frame()
  condition.grid.data <- attr(data, "condition.grid")
  
  preds <- .predictRows(x, times, pars, arglist)

  prediction <- do.call(combine, lapply(1:nrow(pars), function(i) {
    
    prediction <- preds[[i]]
    
    if (is.null(names(prediction))) {
      conditions <- 1
    } else {
      conditions <- names(prediction)
    }
    
    condition.grid <- data.frame(row.names = conditions)
    
    # Augment by parframe metanames and obj.attributes
    mygrid <- pars[i, !colnames(pars) %in% attr(pars, "parameters")]
    mynames <- colnames(mygrid)
    if (length(mynames) > 0) {
      mynames <- paste0(".", mynames)
      colnames(mygrid) <- mynames
      condition.grid <- cbind(condition.grid, mygrid)
    }
    
    # Augment by condition.grid of data
    if (!is.null(condition.grid.data) && ncol(condition.grid.data) > 1) 
      condition.grid <- cbind(condition.grid.data[conditions,], condition.grid)
    
    # Write condition.grid into data
    attr(data, "condition.grid") <- condition.grid

    # Return
    as.data.frame(prediction, data = data)
    
  }))
  
  n <- nrow(prediction)
  
  if (length(data) > 0) {
    attr(data, "condition.grid") <- condition.grid.data
    data <- as.data.frame(data)
    tmp <- combine(prediction, data)
    data <- tmp[-(1:n),]
  }
  
  attr(prediction, "data") <- data
  return(prediction)  
  
  
  
}



## prdfn / obsfn / prdframe / prdlist constructors (moved from classes.R) ----------------------------------------

## Prediction classes ----------------------------------------------------

#' Prediction Function
#'
#' @description A prediction function is a function 
#' `x(..., fixed, deriv, deriv2, hessian, conditions, env, cores, sweep)`.
#' Prediction functions are generated by [Xs()], [Xf()] or [Xd()].
#'
#' @param P2X Transformation function as produced by [Xs()], [Xf()] or [Xd()].
#' @param parameters character, the parameter names. Default `NULL`.
#' @param condition character, the condition name. Default `NULL`, every
#'   condition.
#'
#' @details
#' Prediction functions for different conditions are added by `+`, see
#' [+.fn]. They are composed with observation functions ([obsfn()]) and
#' parameter transformations ([parfn()]) by `*`, see [*.fn].
#'
#' @return Object of class `"prdfn"`, i.e. a function
#' `x(..., fixed, deriv, deriv2, hessian, conditions, env, cores, sweep)` returning a [prdlist].
#' The arguments `times` and `pars` should be passed via `...`, in this order.
#'
#' @example inst/examples/prediction.R
#' @export
prdfn <- function(P2X, parameters = NULL, condition = NULL) {

  st <- .leafState(P2X, "prdfn", condition)
  outfn <- .fnWrap(st)
  attr(outfn, "mappings") <- setNames(list(P2X), condition)
  attr(outfn, "parameters") <- parameters
  attr(outfn, "conditions") <- condition
  attr(outfn, "compileInfo") <- attr(P2X, "compileInfo")
  class(outfn) <- c("prdfn", "fn")
  outfn

}

#' Observation Function
#'
#' @description An observation function is concatenated with a prediction
#' function by `*` to yield a new prediction function, see [prdfn()].
#' Observation functions are generated by [Y()].
#' @param X2Y the low-level observation function, e.g. generated by [Y()].
#' @param parameters character, the parameter names. Default `NULL`.
#' @param condition character, the condition name. Default `NULL`, every
#'   condition.
#' @details Observation functions for different conditions are added by `+`,
#' see [+.fn]. They are composed with observation functions, prediction
#' functions and parameter transformations by `*`, see [*.fn].
#' @return Object of class `obsfn`, i.e. a function `x(..., fixed, deriv, deriv2, hessian, conditions, env, cores, sweep)`
#' which returns a [prdlist]. The arguments `out` (prediction) and `pars` (parameter values)
#' should be passed via the `...` argument.
#' @example inst/examples/prediction.R
#' @export
obsfn <- function(X2Y, parameters = NULL, condition = NULL) {

  st <- .leafState(X2Y, "obsfn", condition)
  outfn <- .fnWrap(st)
  attr(outfn, "mappings") <- setNames(list(X2Y), condition)
  attr(outfn, "parameters") <- parameters
  attr(outfn, "conditions") <- condition
  attr(outfn, "compileInfo") <- attr(X2Y, "compileInfo")
  class(outfn) <- c("obsfn", "fn")
  outfn

}


#' Prediction Frame
#'
#' @description
#' A prediction frame stores model predictions as a matrix with the column
#' `"time"` and one column per variable. Its attributes are `"deriv"`, the
#' first-order sensitivities with respect to the outer parameters (see [P()]),
#' `"deriv2"`, the second-order sensitivities, and `"parameters"`, the inner
#' parameters used for the prediction.
#'
#' Prediction lists ([prdlist]) returned by [Xs()], [Xd()] or [Xf()] consist
#' of prediction frames. A custom prediction function (see `P2X` in [prdfn()])
#' returns one.
#'
#' @param prediction numeric matrix of model predictions. Default `NULL` gives
#'   an empty matrix.
#' @param deriv 3D numeric array `[time, variable, parameter]` of first-order
#'   sensitivities. Default `NULL`.
#' @param deriv2 4D numeric array `[time, variable, parameter, parameter]` of
#'   second-order sensitivities, symmetric in the last two axes. Default
#'   `NULL`.
#' @param parameters named numeric vector of the inner parameters. Default
#'   `NULL`.
#'
#' @return Object of class `prdframe`, a matrix with the attributes above.
#' @seealso [getDerivs()], [getDerivs2()]
#'
#' @export
prdframe <- function(prediction = NULL,
                     deriv = NULL,
                     deriv2 = NULL,
                     parameters = NULL) {

  out <- if (!is.null(prediction)) as.matrix(prediction) else matrix(, 0, 0)

  attr(out, "deriv") <- deriv
  attr(out, "deriv2") <- deriv2
  attr(out, "parameters") <- parameters
  class(out) <- c("prdframe", "matrix")

  return(out)
}


#' Prediction List
#'
#' @description A prediction list holds model predictions, one [prdframe] per
#' condition or parameter set. Prediction functions return prediction lists.
#' @param ... for `prdlist()`, objects of class [prdframe]; unnamed entries
#'   are named `"1"`, `"2"`, ... For `as.prdlist()`, not used.
#' @return Object of class `prdlist`.
#' @seealso [prdfn()], [plotCombined()]
#' @export
prdlist <- function(...) {
  mylist <- list(...)
  mynames <- names(mylist)
  if (is.null(mynames)) mynames <- as.character(1:length(mylist))
  as.prdlist(mylist, mynames)
}



