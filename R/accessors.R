## Accessors for dMod objects -------------------------------------------------
## get*/controls/modelname layer, split out of classes.R.

## General purpose functions for different dMod classes ------------------------------

#' List, Get and Set Controls
#'
#' @description Reads and changes the settings an object of class `objfn`,
#' `parfn`, `prdfn` or `obsfn` was built with. `controls(x)` lists the
#' available controls; with `name`, the value of that control is returned or
#' replaced. Elements of a list or data frame value can be changed with `$` or
#' `[]`.
#'
#' @details A control belongs to the function built with it and is shared by
#' every sum or composition that contains this function. Changing it on one of
#' them changes it for all of them.
#'
#' On an `fn`, the lookup covers the functions `x` is made of, e.g. the
#' factors of `g * x * p`. `condition` restricts it to the functions defined
#' for that condition; a function built without a condition counts for every
#' condition. With `condition = NULL`, the getter returns the value of the
#' first function that holds the control and the setter changes every function
#' that holds it. On an `objfn`, including sums, scaled objectives and
#' objectives composed with a parameter transformation, the same holds for the
#' objective terms. Setting a control that no function holds is an error, as
#' is an unknown condition.
#'
#' The second positional argument is `name` for an `objfn` and `condition`
#' for an `fn`; name both arguments to avoid confusion.
#'
#' @param x object of class `objfn` or `fn`.
#' @param ... not used.
#' @return `controls(x)` prints the available controls and returns their names
#' invisibly; for an `fn`, as a list with one entry per function. With `name`,
#' the value of the control, or `NULL` if no function holds it. The setter
#' returns `x`.
#' @seealso [normL2()], [P()], [Y()]
#' @examples
#' ## parfn with condition
#' p <- P(eqnvec(x = "-a*x"), method = "implicit", condition = "C1")
#' controls(p)
#' controls(p, condition = "C1", name = "keep.root")
#' controls(p, condition = "C1", name = "keep.root") <- FALSE
#'
#' ## obsfn without condition
#' g <- Y(g = eqnvec(y = "s*x"), f = NULL, states = "x", parameters = "s")
#' controls(g)
#' controls(g, name = "attach.input")
#' controls(g, name = "attach.input") <- TRUE
#' @export
controls <- function(x, ...) {
  UseMethod("controls", x)
}



.lscontrolsObjfn <- function(x) {

  unique(unlist(lapply(.controlTargets(x), function(t) names(environment(t)$controls))))

}

# The objectives that hold a control: the objective itself, or every
# objective it is built from whose controls include `name` (any control when
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

# A control under its deprecated name answers under the current one.
.controlAlias <- function(name) {
  if (!identical(name, "optionsOde")) return(name)
  warning("controls: 'optionsOde' is deprecated, use 'options'.", call. = FALSE)
  "options"
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
#' @param name character, the name of the control. Default `NULL` lists the
#'   controls.
controls.objfn <- function(x, name = NULL, ...) {

  if (is.null(name)) return(.lscontrolsObjfn(x))
  tg <- .controlTargets(x, name)
  if (!length(tg)) return(NULL)
  environment(tg[[1L]])$controls[[name]]
}

#' @export
#' @rdname controls
#' @param condition character, the condition name. Default `NULL`, every
#'   condition.
controls.fn <- function(x, condition = NULL, name = NULL, ...) {

  condition <- .controlCondition(x, condition)
  name <- .controlAlias(name)
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
#' @param value the new value of the control.
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
  name <- .controlAlias(name)
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


#' Extract the First Derivatives of an Object
#'
#' @param x object of class `parvec`, `prdframe`, `prdlist`, `objlist`, or a
#'   list of such objects.
#' @param ... not used.
#'
#' @return Depends on the class of `x`:
#' * `parvec`: the Jacobian matrix.
#' * `prdframe`: a `prdframe` with the column `time` and one column per pair
#'   of variable and parameter.
#' * `prdlist`: a `prdlist` of such `prdframe`s.
#' * `objlist`: the gradient, a named numeric vector.
#' * `list`: a list with the result for each element.
#'
#' An error is raised if `x` has no first derivatives.
#' @seealso [getDerivs2()], [parvec()]
#'
#' @examples
#' J <- matrix(c(1, 0, 0, 2), 2, dimnames = list(c("a", "b"), c("x", "y")))
#' p <- as.parvec(c(a = 1, b = 2), deriv = J)
#' getDerivs(p)
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


#' Extract the Second Derivatives of an Object
#'
#' @param x object of class `parvec`, `prdframe`, `prdlist`, `objlist`, or a
#'   list of such objects.
#' @param ... not used.
#'
#' @return Depends on the class of `x`:
#' * `parvec`: a 3D array `[parameter, theta, theta]`.
#' * `prdframe`: a `prdframe` with the column `time` and one column per
#'   variable and pair of parameters, for the upper triangle including the
#'   diagonal.
#' * `prdlist`: a `prdlist` of such `prdframe`s.
#' * `objlist`: the Hessian matrix.
#' * `list`: a list with the result for each element.
#'
#' An error is raised if `x` has no second derivatives.
#' @seealso [getDerivs()]
#'
#' @examples
#' H <- array(0, dim = c(1, 2, 2), dimnames = list("a", c("x", "y"), c("x", "y")))
#' H["a", "x", "y"] <- H["a", "y", "x"] <- 1
#' p <- as.parvec(c(a = 1), deriv2 = H)
#' getDerivs2(p)
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


#' Extract the Parameters of an Object
#'
#' @param x object from which the parameters are extracted.
#' @param ... further objects; their parameters are added to the result.
#' @param conditions character, restrict the result to these conditions.
#'   Default `NULL`. Used by the methods for `fn` and `prdlist`.
#' @return Character vector of parameter names. For an `odemodel`, the states
#'   and parameters; for a `prdlist`, a list with one entry per condition.
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

#' Extract the Conditions of an Object
#'
#' @param x object of class `fn` or a named list such as a `datalist`.
#' @param ... not used.
#' @return Character vector of condition names, `NULL` for a function defined
#'   for every condition.
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

#' Get and Set Modelname
#'
#' @description The model name is the base name of the generated sources and
#' shared objects of a function such as a prediction function, parameter
#' transformation or objective function.
#'
#' @param x object of class `prdfn`, `parfn`, `obsfn` or `objfn`, or a
#'   character naming such an object in the calling environment.
#' @param ... further objects; their model names are added to the result.
#' @param conditions character, restrict the result to these conditions.
#'   Default `NULL`.
#' @return Character vector of model names. The setter returns `x`.
#' @seealso [compile()], [loadDLL()]
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
#' @param value character, the new model name. Files are not renamed.
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





#' Extract the Equations of an Object
#'
#' @param x object of class `odemodel` or `fn`.
#' @param conditions character or numeric, restrict the result to these
#'   conditions. Default `NULL`.
#' @return A list of `eqnvec` objects, one per condition. If `conditions` has
#'   length one, the `eqnvec` itself. For an `odemodel`, its `eqnvec`.
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

#' Extract the Observables of an Object
#'
#' Generic without methods in dMod2; other packages can provide them.
#'
#' @param x object from which the observables are extracted.
#' @param ... not used.
#' @return The observables, as a character vector.
#' @export
getObservables <- function(x, ...) {
  UseMethod("getObservables", x)
}

