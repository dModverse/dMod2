## Accessors for dMod objects -------------------------------------------------
## get*/controls/modelname layer, split out of classes.R.

## General purpose functions for different dMod classes ------------------------------

#' List, get and set controls for different functions
#'
#' @description Applies to objects of class `objfn`,
#' `parfn`, `prdfn` and `obsfn`. Allows to manipulate
#' different arguments that have been set when creating the
#' objects. On a sum of objectives, a control is read from and written to
#' every summand that has it, e.g. the `multipleShootingControl` of a
#' [normL2()] term. The same holds for an objective scaled by `%.*%` and
#' for an objective composed with a parameter transformation, which reach the
#' controls of the objective inside.
#' @details If called without further arguments, `controls(x)` lists the
#' available controls within an object. Calling `controls()` with `name`
#' and `condition` returns the control value. The value can be overwritten. If
#' a list or data.frame is returned, elements of those can be manipulated by the
#' `$`- or `[]`-operator.
#'
#' A control lives in the function that was built with it, which every sum or
#' composition containing that function shares. A change made on `x` is
#' therefore seen by every object built from `x`, and a change made on a
#' composition reaches the factor it came from.
#'
#' On an `fn`, the control is looked up in the functions `x` is made of, so
#' `g * x * p` reaches the controls of `g`, `x` and `p`. `condition` restricts
#' the lookup to the functions that answer for that condition; a function built
#' without a condition answers for all of them, so a control set on it for one
#' condition applies to every condition. With `condition = NULL`, the getter
#' returns the value of the first function that holds the control and the
#' setter writes every function that holds it. A name that several factors
#' share, `attach.input` for instance, is set on all of them; set it on the
#' factor itself to reach that one alone. Setting a control that no function
#' holds is an error, as is a condition `x` does not know.
#'
#' @param x function
#' @param ... arguments going to the appropriate S3 methods
#' @return Either a print-out or the values of the control. Listing the
#' controls of an `fn` prints them per function and returns the names
#' invisibly, as a list with one entry per function.
#' @examples
#' \dontrun{
#'   ## parfn with condition
#'   p <- P(eqnvec(x = "-a*x"), method = "implicit", condition = "C1")
#'   controls(p)
#'   controls(p, "C1", "keep.root")
#'   controls(p, "C1", "keep.root") <- FALSE
#'   
#'   ## obsfn with NULL condition
#'   g <- Y(g = eqnvec(y = "s*x"), f = NULL, states = "x", parameters = "s")
#'   controls(g)
#'   controls(g, NULL, "attach.input")
#'   controls(g, NULL, "attach.input") <- FALSE
#' }
#' @export
controls <- function(x, ...) {
  UseMethod("controls", x)
}



.lscontrolsObjfn <- function(x) {

  unique(unlist(lapply(.controlTargets(x), function(t) names(environment(t)$controls))))

}

# The objectives that hold a control: the objective itself, or every
# objective it is built from whose controls carry `name` (any control when
# `name` is NULL). A sum records its summands in `terms`; an objective scaled
# by %.*% or composed with a parfn records the objective it wraps in `wrapped`.
# Either is the closure that gets called, so a change there reaches the whole.
#
# `wrapped` is not `terms` on purpose: .objTerms() reads `terms` as the
# summands of a sum, and a scaled objective is not a sum of its inner one.
.controlTargets <- function(x, name = NULL) {
  has <- function(f) {
    ctl <- if (is.function(f) && !is.primitive(f)) environment(f)$controls
    is.list(ctl) && (is.null(name) || name %in% names(ctl))
  }
  if (has(x)) return(list(x))
  inner <- c(attr(x, "terms", exact = TRUE), attr(x, "wrapped", exact = TRUE))
  unlist(lapply(inner, .controlTargets, name = name), recursive = FALSE)
}

# The controls environment of a kernel, NULL for a kernel without one. `$` on
# an environment does not inherit, so a `controls` further up is not taken.
.kernelControls <- function(k) {
  e <- if (is.function(k) && !is.primitive(k)) environment(k)
  if (is.null(e) || !is.list(e$controls)) return(NULL)
  e
}

# A setting of a kernel that it may keep both as a control and as an
# attribute, `forcings` or `events`. The attribute is a copy made when the
# kernel was built, the control is what the kernel runs with and what
# controls<- changes, so the control wins where there is one.
.kernelSetting <- function(k, what) {
  e <- .kernelControls(k)
  if (!is.null(e) && what %in% names(e$controls)) e$controls[[what]]
  else attr(k, what, exact = TRUE)
}

# The leaves of `x` that hold controls, each as list(env, kind, condition),
# once per kernel. Restricted to the leaves answering `condition` and to those
# holding `name`, where given.
.fnControlLeaves <- function(x, condition = NULL, name = NULL) {
  leaves <- lapply(.fnLeaves(x, condition), function(l) {
    l$env <- .kernelControls(l$kernel)
    l
  })
  leaves <- Filter(function(l) !is.null(l$env) &&
                     (is.null(name) || name %in% names(l$env$controls)), leaves)
  # One kernel can sit in several places; format() names an environment by its
  # address, which keeps this linear in the number of leaves.
  keys <- vapply(leaves, function(l) format(l$env), "")
  leaves[!duplicated(keys)]
}

# The condition as a name, or NULL for all of them. A number picks a
# condition by position, as `mappings[[i]]` used to.
.controlCondition <- function(x, condition) {
  if (is.null(condition)) return(NULL)
  conds <- attr(x, "conditions")
  if (is.numeric(condition)) {
    n <- max(1L, length(conds))
    if (any(condition < 1 | condition > n))
      stop("controls: condition index ", paste(condition, collapse = ", "),
           " is out of range, the object has ", n, " condition(s).", call. = FALSE)
    return(if (is.null(conds)) NULL else conds[condition])
  }
  # An object without conditions answers for any condition.
  bad <- setdiff(condition, conds)
  if (!is.null(conds) && length(bad))
    stop("controls: unknown condition ", paste0("'", bad, "'", collapse = ", "),
         ". Available: ", paste(conds, collapse = ", "), ".", call. = FALSE)
  condition
}

.lscontrolsFn <- function(x, condition = NULL) {

  leaves <- .fnControlLeaves(x, condition)
  composed <- inherits(x, "composed")
  label <- function(l) {
    conds <- l$condition
    if (!is.null(conds) && !is.null(condition)) conds <- intersect(conds, condition)
    if (length(conds) > 3L)
      conds <- c(conds[1:3], paste("and", length(conds) - 3L, "more"))
    lab <- paste(conds, collapse = ", ")
    if (composed) trimws(paste(l$kind, lab)) else lab
  }
  out <- lapply(leaves, function(l) names(l$env$controls))
  names(out) <- vapply(leaves, label, "")
  for (i in seq_along(out)) {
    cat(names(out)[i], ":\n", sep = "")
    print(out[[i]])
  }
  invisible(out)

}

#' @export
#' @rdname controls
#' @param name character, the name of the control
controls.objfn <- function(x, name = NULL, ...) {

  if (is.null(name)) return(.lscontrolsObjfn(x))
  tg <- .controlTargets(x, name)
  if (!length(tg)) return(NULL)
  environment(tg[[1L]])$controls[[name]]
}

#' @export
#' @rdname controls
#' @param condition character, the condition name
controls.fn <- function(x, condition = NULL, name = NULL, ...) {

  condition <- .controlCondition(x, condition)
  if (is.null(name)) return(.lscontrolsFn(x, condition))

  tg <- .fnControlLeaves(x, condition, name)
  if (!length(tg)) return(NULL)
  tg[[1L]]$env$controls[[name]]

}


#' @export
#' @rdname controls
"controls<-" <- function(x, ..., value) {
  UseMethod("controls<-", x)
}


#' @export
#' @param value the new value
#' @rdname controls
"controls<-.objfn" <- function(x, name, ..., value) {
  tg <- .controlTargets(x, name)
  if (!length(tg))
    stop("controls<-: the objective has no control '", name, "'. Available: ",
         paste(.lscontrolsObjfn(x), collapse = ", "), ".", call. = FALSE)
  # [<- keeps an entry set to NULL, where [[<- would drop it
  for (t in tg) environment(t)$controls[name] <- list(value)
  return(x)
}

#' @export
#' @rdname controls
"controls<-.fn" <- function(x, condition = NULL, name, ..., value) {
  condition <- .controlCondition(x, condition)
  tg <- .fnControlLeaves(x, condition, name)
  if (!length(tg)) {
    avail <- unique(unlist(lapply(.fnControlLeaves(x, condition),
                                  function(l) names(l$env$controls))))
    stop("controls<-: no function in the object has a control '", name, "'",
         if (!is.null(condition))
           paste0(" for condition ", paste(condition, collapse = ", ")),
         ". Available: ",
         if (length(avail)) paste(avail, collapse = ", ") else "none", ".",
         call. = FALSE)
  }
  # [<- keeps an entry set to NULL, where [[<- would drop it
  for (l in tg) l$env$controls[name] <- list(value)
  return(x)
}


#' Extract the first derivatives of an object
#'
#' Generic function to extract first-order derivatives
#' from various model-related objects such as `parvec`, `prdframe`, or lists thereof.
#'
#' The output format depends on the class of the input object.
#'
#' @param x Object from which the first derivatives should be extracted.
#'   Supported classes are `parvec`, `prdframe`, `prdlist`, and `list`.
#' @param ... Additional arguments passed to specific methods (currently unused).
#'
#' @return The structure of the returned object depends on the class of `x`:
#' \itemize{
#'   \item `parvec`: a matrix containing first-order parameter derivatives.
#'   \item `prdframe`: a `prdframe` containing time and first-order sensitivities
#'     of each model variable with respect to all parameters.
#'   \item `prdlist`: a `prdlist` whose elements are first-derivative `prdframe`s.
#'   \item `list`: a list of derivative objects, depending on the elements.
#'   \item `objlist`: directly returns the stored gradient (named numeric vector).
#' }
#'
#' @examples
#' \dontrun{
#' # Extract sensitivities from a model prediction frame:
#' d1 <- getDerivs(myprdframe)
#'
#' # Extract parameter derivatives from a parameter vector:
#' getDerivs(myparvec)
#' }
#'
#' @export
getDerivs <- function(x, ...) {
  UseMethod("getDerivs", x)
}

#' @export
#' @rdname getDerivs
getDerivs.parvec <- function(x, ...) {

  derivs <- attr(x, "deriv")
  if (is.null(derivs))
    stop("Object does not contain first-order derivatives.")
  
  return(derivs)
}

#' @export
#' @rdname getDerivs
getDerivs.prdframe <- function(x, ...) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  times  <- x[, "time", drop = FALSE]
  derivs <- attr(x, "deriv")
  if (is.null(derivs))
    stop("Object does not contain first-order derivatives.")
  
  dn <- dimnames(derivs)
  n  <- dim(derivs)[1]
  v  <- dim(derivs)[2]
  d  <- dim(derivs)[3]

  varnames <- dn[[2]] %||% paste0("var", seq_len(v))
  parnames <- dn[[3]] %||% paste0("par", seq_len(d))

  derivswide <- times

  for (i in seq_len(v)) {
    m <- matrix(derivs[, i, ], nrow = n, ncol = d)
    colnames(m) <- paste0("\u2202", varnames[i], "/\u2202", parnames)
    derivswide <- cbind(derivswide, m)
  }
  
  prdframe(
    prediction = derivswide,
    parameters = attr(x, "parameters")
  )
}




#' @export
#' @rdname getDerivs
getDerivs.prdlist <- function(x, ...) {

  as.prdlist(
    lapply(x, function(myx) {
      getDerivs(myx, ...)
    }),
    names = names(x)
  )

}

#' @export
#' @rdname getDerivs
getDerivs.list <- function(x, ...) {

  lapply(x, function(myx) getDerivs(myx))

}


#' @export
#' @rdname getDerivs
getDerivs.objlist <- function(x, ...) {

  x$gradient

}


#' Extract second-order derivatives from an object
#'
#' Generic accessor for the `deriv2` attribute (or `hessian` field, in the
#' case of `objlist`) attached to dMod objects.
#'
#' @param x Object from which the second derivatives should be extracted.
#'   Supported classes are `parvec`, `prdframe`, `prdlist`, `list`, and
#'   `objlist`.
#' @param ... Additional arguments passed to specific methods (currently unused).
#'
#' @return The structure of the returned object depends on the class of `x`:
#' \itemize{
#'   \item `parvec`: a 3D array `[p, theta, theta]` of second derivatives.
#'   \item `prdframe`: a 4D array `[time, variable, theta, theta]`.
#'   \item `prdlist`: a list of `prdframe` second-derivative arrays.
#'   \item `objlist`: the stored `hessian` matrix.
#'   \item `list`: a list of derivative objects, depending on the elements.
#' }
#'
#' @examples
#' \dontrun{
#' d2 <- getDerivs2(myprdframe)
#' getDerivs2(myparvec)
#' }
#'
#' @export
getDerivs2 <- function(x, ...) {
  UseMethod("getDerivs2", x)
}

#' @export
#' @rdname getDerivs2
getDerivs2.parvec <- function(x, ...) {

  derivs2 <- attr(x, "deriv2")
  if (is.null(derivs2))
    stop("Object does not contain second-order derivatives.")

  return(derivs2)
}

#' @export
#' @rdname getDerivs2
getDerivs2.prdframe <- function(x, ...) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  times   <- x[, "time", drop = FALSE]
  derivs2 <- attr(x, "deriv2")
  if (is.null(derivs2))
    stop("Object does not contain second-order derivatives.")

  dn <- dimnames(derivs2)
  n  <- dim(derivs2)[1]
  v  <- dim(derivs2)[2]
  d  <- dim(derivs2)[3]

  varnames  <- dn[[2]] %||% paste0("var", seq_len(v))
  parnames1 <- dn[[3]] %||% paste0("par", seq_len(d))
  parnames2 <- dn[[4]] %||% parnames1

  # Schwarz' theorem: H is symmetric in the last two axes. Emit only the
  # independent entries (upper triangle, k1 <= k2). Diagonal entries get
  # `\u2202\u00b2f/\u2202p\u00b2`, off-diagonal `\u2202\u00b2f/\u2202p1\u2202p2`.
  pair_idx <- which(upper.tri(matrix(0, d, d), diag = TRUE), arr.ind = TRUE)

  derivswide <- times
  for (i in seq_len(v)) {
    m <- matrix(0, nrow = n, ncol = nrow(pair_idx))
    cols <- character(nrow(pair_idx))
    for (j in seq_len(nrow(pair_idx))) {
      k1 <- pair_idx[j, 1L]; k2 <- pair_idx[j, 2L]
      m[, j] <- derivs2[, i, k1, k2]
      cols[j] <- if (k1 == k2)
        paste0("\u2202\u00b2", varnames[i], "/\u2202", parnames1[k1], "\u00b2")
      else
        paste0("\u2202\u00b2", varnames[i], "/\u2202", parnames1[k1], "\u2202", parnames2[k2])
    }
    colnames(m) <- cols
    derivswide <- cbind(derivswide, m)
  }

  prdframe(
    prediction = derivswide,
    parameters = attr(x, "parameters")
  )
}

#' @export
#' @rdname getDerivs2
getDerivs2.prdlist <- function(x, ...) {

  as.prdlist(
    lapply(x, function(myx) {
      getDerivs2(myx, ...)
    }),
    names = names(x)
  )

}

#' @export
#' @rdname getDerivs2
getDerivs2.list <- function(x, ...) {

  lapply(x, function(myx) getDerivs2(myx))

}

#' @export
#' @rdname getDerivs2
getDerivs2.objlist <- function(x, ...) {

  x$hessian

}


#' Extract the parameters of an object
#'
#' @param x object from which the parameters are extracted
#' @param ... further objects; when supplied, parameters of all objects
#'   are unioned (each dispatched separately).
#' @param conditions character vector specifying the conditions to
#'   which `getParameters` is restricted (only honored by methods that
#'   carry per-condition parameter mappings).
#' @return The parameters in a format that depends on the class of `x`.
#' @export
getParameters <- function(x, ..., conditions = NULL) {
  if (...length() > 0L) {
    return(Reduce("union", lapply(list(x, ...), getParameters, conditions = conditions)))
  }
  UseMethod("getParameters")
}



#' @export
#' @rdname getParameters
getParameters.odemodel <- function(x, ..., conditions = NULL) {

  parameters <- c(
    attr(x$func, "variables"),
    attr(x$func, "parameters")
  )

  return(parameters)

}


#' @export
#' @rdname getParameters
getParameters.fn <- function(x, ..., conditions = NULL) {

  if (is.null(conditions)) {
    parameters <- attr(x, "parameters")
  } else {
    mappings <- attr(x, "mappings")
    mappings <- mappings[intersect(names(mappings), conditions)]
    parameters <- Reduce("union",
                         lapply(mappings, function(m) attr(m, "parameters"))
    )
  }

  return(parameters)

}
#' @export
#' @rdname getParameters
getParameters.parvec <- function(x, ..., conditions = NULL) {

  names(x)

}

#' @export
#' @rdname getParameters
getParameters.prdframe <- function(x, ..., conditions = NULL) {

  attr(x, "parameters")

}

#' @export
#' @rdname getParameters
getParameters.prdlist <- function(x, ..., conditions = NULL) {

  select <- 1:length(x)
  if (!is.null(conditions)) select <- intersect(names(x), conditions)
  lapply(x[select], function(myx) getParameters(myx))

}

#' @export
#' @rdname getParameters
getParameters.eqnlist <- function(x, ..., conditions = NULL) {
  comp_exprs <- character(0)
  if (!is.null(x$compartments)) {
    comp_exprs <- c(
      vapply(x$compartments, function(c) c$volume, character(1)),
      vapply(x$compartments, function(c) if (is.null(c$rule)) "" else c$rule, character(1))
    )
    comp_exprs <- comp_exprs[nzchar(comp_exprs)]
  }
  unique(c(getSymbols(x$states), getSymbols(x$rates), getSymbols(comp_exprs)))
}

#' @export
#' @rdname getParameters
getParameters.eventlist <- function(x, ..., conditions = NULL) {
  idx <- match(c("time", "value", "root"), names(x))
  idx[!is.na(idx)]
  Reduce(union, lapply(x[idx], getSymbols))
}

#' @export
#' @rdname getParameters
getParameters.eqnvec <- function(x, ..., conditions = NULL) {
  getSymbols(x)
}

#' Extract the conditions of an object
#'
#' @param x object from which the conditions should be extracted
#' @param ... additional arguments (not used right now)
#' @return The conditions in a format that depends on the class of `x`.
#' @export
getConditions <- function(x, ...) {
  UseMethod("getConditions", x)
}


#' @export
#' @rdname getConditions
getConditions.list <- function(x, ...) {

  names(x)

}


#' @export
#' @rdname getConditions
getConditions.fn <- function(x, ...) {

  attr(x, "conditions")

}

#' Get and set modelname
#'
#' @description The modelname attribute refers to the name of a C file associated with
#' a dMod function object like prediction-, parameter transformation- or
#' objective functions.
#'
#' @param x object of type `prdfn`, `parfn`, `objfn`, or a character naming
#'   such an object in the calling environment.
#' @param ... further objects; when supplied, model names of all objects
#'   are unioned (each dispatched separately).
#' @param conditions character vector of conditions
#' @return character vector of model names, corresponding to C files
#' in the local directory.
#'
#' @export
modelname <- function(x = NULL, ..., conditions = NULL) {
  if (...length() > 0L) {
    return(Reduce("union", lapply(list(x, ...), modelname, conditions = conditions)))
  }
  UseMethod("modelname")
}

#' @export
#' @rdname modelname
modelname.NULL <- function(x = NULL, ..., conditions = NULL) NULL

#' @export
#' @rdname modelname
modelname.character <- function(x = NULL, ..., conditions = NULL) {

  modelname(get(x), conditions = conditions)

}

#' @export
#' @rdname modelname
modelname.objfn <- function(x = NULL, ..., conditions = NULL) {

  attr(x, "modelname")

}

#' @export
#' @rdname modelname
modelname.fn <- function(x = NULL, ..., conditions = NULL) {

  mappings <- attr(x, "mappings")
  select <- 1:length(mappings)
  if (!is.null(conditions)) select <- intersect(names(mappings), conditions)
  modelnames <- Reduce("union",
                       lapply(mappings[select], function(m) attr(m, "modelname"))
  )

  return(modelnames)

}



#' @export
#' @rdname modelname
#' @param value character, the new modelname (does not change the C file)
"modelname<-" <- function(x, ..., value) {
  UseMethod("modelname<-", x)
}

#' @export
#' @rdname modelname
"modelname<-.fn" <- function(x, conditions = NULL, ..., value) {
  
  mappings <- attr(x, "mappings")
  if (!is.null(mappings)) {
    select <- seq_along(mappings)
    if (!is.null(conditions)) select <- intersect(names(mappings), conditions)
    if (length(value) == 1) value <- rep(value, length.out = length(select))
    
    for (i in select) {
      m <- mappings[[i]]
      
      if ("composed" %in% class(m)) {
        modelname(m) <- value[i %% length(value) + 1]  # recursive
      } else {
        attr(m, "modelname") <- value[i %% length(value) + 1]
        # A deSolve leaf keeps the compiled model in its closure, and cOde
        # reads the shared object to load from there while the entry point
        # names come from the object's own value. The rename therefore has to
        # reach the closure, including when the leaf sits inside a composition.
        e <- if (is.function(m)) environment(m) else NULL
        if (!is.null(e)) {
          if (!is.null(e[["func"]]))
            attr(e[["func"]], "modelname") <- value[i %% length(value) + 1]
          if (!is.null(e[["extended"]]))
            attr(e[["extended"]], "modelname") <- value[i %% length(value) + 1]
          # A cppDE leaf caches a prepared batch handle bound to the shared
          # object the entry point was resolved from. The rename is exactly
          # what invalidates that binding, so drop the cache with it.
          if (is.environment(e[["bcache"]]))
            rm(list = ls(e[["bcache"]], all.names = TRUE), envir = e[["bcache"]])
        }
      }
      mappings[[i]] <- m
    }
    
    attr(x, "mappings") <- mappings
    
  } else {
    attr(x, "modelname") <- value[1]
  }
  
  x
}


#' @export
#' @rdname modelname
"modelname<-.objfn" <- function(x, conditions = NULL, ..., value) {
  attr(x, "modelname") <- value
  return(x)
}





#' Extract the equations of an object
#'
#' @param x object from which the equations should be extracted
#' @param conditions character or numeric vector specifying the conditions to
#' which `getEquations` is restricted. If `conditions` has length one,
#' the result is not returned as a list.
#' @return The equations as list of `eqnvec` objects.
#' @export
getEquations <- function(x, conditions = NULL) {

    UseMethod("getEquations", x)

}



#' @export
#' @rdname getEquations
getEquations.odemodel <- function(x, conditions = NULL) {

  attr(x$func, "equations")

}



#' @export
#' @rdname getEquations
getEquations.prdfn <- function(x, conditions = NULL) {

  mappings <- attr(x, "mappings")

  if (is.null(conditions)) {
    equations <- lapply(mappings, function(m) attr(m, "equations"))
    return(equations)
  }

  if (!is.null(conditions)) {
    mappings <- mappings[conditions]
    equations <- lapply(mappings, function(m) attr(m, "equations"))
    if (length(equations) == 1) {
      return(equations[[1]])
    } else {
      return(equations)
    }
  }

}


#' @export
#' @rdname getEquations
getEquations.fn <- function(x, conditions = NULL) {

  mappings <- attr(x, "mappings")

  if (is.null(conditions)) {
    equations <- lapply(mappings, function(m) attr(m, "equations"))
    return(equations)
  }

  if (!is.null(conditions)) {
    mappings <- mappings[conditions]
    equations <- lapply(mappings, function(m) attr(m, "equations"))
    if (length(equations) == 1) {
      return(equations[[1]])
    } else {
      return(equations)
    }
  }

}

#' Extract the observables of an object
#'
#' @param x object from which the equations should be extracted
#' @param ... not used
#' @return The equations as a character.
#' @export
getObservables <- function(x, ...) {
  UseMethod("getObservables", x)
}

