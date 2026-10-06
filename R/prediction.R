

#' Prediction function of an ODE model with sensitivities
#'
#' @description Turns an [odemodel()] into a prediction function that
#' integrates the states and, with `deriv = TRUE`, their sensitivities.
#'
#' @param odemodel An [odemodel()]. Its backend decides which method runs.
#' @param forcings `data.frame` with columns `name` (character or factor),
#'   `time` and `value`, one row per data point of a forcing declared in
#'   `odemodel(forcings = )`. See section Forcings.
#' @param events `deSolve` backend only. An [eventlist], or a `data.frame`
#'   [as.eventlist()] accepts. Applied to the states only, the sensitivities
#'   are not corrected, so define events in [odemodel()] unless the prediction
#'   is used for simulation alone. The other backends take events only from
#'   [odemodel()].
#' @param names Character vector, the states and forcings to return. `NULL`
#'   returns all states followed by the forcings.
#' @param condition `NULL` for a prediction valid in every condition, or the
#'   name of the condition it belongs to.
#' @param optionsOde Named list of solver options for the solves without
#'   sensitivities. See section Solver options.
#' @param optionsSens Named list of solver options for the solves with
#'   sensitivities. See section Solver options.
#' @param optionsReverse Named list controlling the backward pass of
#'   `sweep = "reverse"`. Needs `backend = "cppDE"` and a model built with
#'   `derivMode = "reverse"`. `NULL`, the default, sweeps the step grid of the
#'   value pass. Entries:
#'   * `refine`: `TRUE` holds each step of the sweep to the error test of the
#'     CVODES backward problem under the tolerances of the solve, see
#'     [cppDE::adjointControl()].
#'   * `gradtol`: under `refine`, the absolute tolerance on each step's share
#'     of the gradient.
#' @param fcontrol `deSolve` backend only. List with the interpolation
#'   settings of the forcings, `method`, `rule`, `f` and `ties` as in
#'   [stats::approxfun()], passed to [deSolve::ode()] as its `fcontrol`.
#' @param ... Not used.
#'
#' @section Solver options:
#' On the `cppDE` and `Sundials` backends `optionsOde` and `optionsSens` take
#' the entries below and pass them to [cppDE::solveODE()]. An entry not given
#' keeps its default; an unknown entry is ignored with a warning.
#'
#' * `atol`, `rtol`: absolute and relative error tolerance, default `1e-6`.
#' * `maxsteps`: largest number of steps of one solve, default `1e6`.
#' * `maxattemps`: consecutive rejected steps before the solve fails,
#'   default `50`.
#' * `hini`: first step size, default `0`, which estimates it.
#' * `roottol`: tolerance of root-triggered events and of
#'   `rootfunc = "equilibrate"`, default `1e-6`.
#' * `maxroot`: number of times a root event may fire, default `1`.
#' * `onFailure`: `"stop"` (default), `"warn"` or `"silent"` on a failed
#'   solve. With `"warn"` and `"silent"` the result ends at the time reached,
#'   so an objective sees fewer data points rather than a failure.
#' * `traceFile`: CSV file the per-step trace is written to, default `NULL`;
#'   needs `odemodel(..., stepTrace = TRUE)`.
#' * `sensErrCon`: `optionsSens` only. `FALSE` takes the sensitivities out of
#'   the error test, which is cheaper and gives less accurate sensitivities;
#'   default `TRUE`.
#'
#' `optionsOde` applies to `deriv = FALSE` and to the value and backward pass
#' of `sweep = "reverse"`, `optionsSens` to `deriv = TRUE` and to
#' `sweep = "reverse"` with `deriv2 = TRUE`. A fit with the default forward
#' sweep therefore runs under `optionsSens`. The integration
#' method is fixed when the model is built, `odemodel(..., method = )`.
#'
#' On the `deSolve` backend both lists are arguments of [deSolve::ode()],
#' among them `method` (default `"lsoda"` for `optionsOde` and `"lsodes"` for
#' `optionsSens`), `atol`, `rtol`, `maxsteps` and `hini`. A list given there
#' replaces the default, so it names `method` if it needs one.
#'
#' @section Forcings:
#' On the `cppDE` and `Sundials` backends each forcing is the monotone cubic
#' Hermite interpolant (PCHIP) of its data points, see
#' [cppDE::forcingValues()]. Outside the points it holds the value of the
#' first and the last one; a single point gives a constant. The integrator
#' does not stop at the points, so a jump is best given by an event or by
#' two close points.
#'
#' On the `deSolve` backend the forcings are interpolated linearly and held
#' constant outside their data, unless `fcontrol` says otherwise.
#'
#' The forcings are returned after the states, without derivatives, so
#' observables of [Y()] can contain them. Different forcings per condition
#' come from one `Xs()` per condition, joined by `+`.
#'
#' @section Controls:
#' `forcings`, `names` and the solver options are kept, as given, in the
#' controls of the returned function and read at every solve, so [controls()]
#' can change them later; on the `deSolve` backend `events` and `fcontrol` as
#' well. On the `cppDE` and `Sundials` backends a replaced `optionsOde` or
#' `optionsSens` is merged over the defaults, so it only needs the entries it
#' changes.
#'
#' @return A [prdfn], called as `x(times, pars, fixed = NULL, deriv = TRUE,
#'   deriv2 = FALSE)`. It returns a [prdlist] with one [prdframe] per
#'   condition: the states over `times`, with the sensitivities in
#'   `attr(, "deriv")`, `[time, state, parameter]`, when `deriv = TRUE`, and
#'   the second-order sensitivities in `attr(, "deriv2")` when `deriv2 = TRUE`
#'   (`cppDE` backend, model built with `deriv2 = TRUE`). If `pars` carries a
#'   `deriv` attribute, as the output of a [parfn] does, the sensitivities are
#'   taken with respect to the outer parameters.
#'
#' @seealso [Xf()] for predictions without sensitivities, [cppDE::solveODE()].
#' @example inst/examples/Xs.R
#' @export
Xs <- function(odemodel, ...) {
  UseMethod("Xs", odemodel)
}

#' @export
#' @rdname Xs
Xs.deSolve <- function(odemodel, forcings = NULL, events = NULL, names = NULL, condition = NULL,
                       optionsOde = list(method = "lsoda"), optionsSens = list(method = "lsodes"),
                       fcontrol = NULL, ...) {
  
  func <- odemodel$func
  extended <- odemodel$extended
  if (is.null(extended)) warning("Element 'extended' empty. ODE model does not contain sensitivities.")
  
  myforcings <- forcings
  myevents <- events
  myfcontrol <- fcontrol
  
  if (!is.null(attr(func, "events")) & !is.null(myevents))
    warning("Events already defined in odemodel. Additional events in Xs() will be ignored.")
  if (is.null(attr(func, "events")) & !is.null(myevents))
    message("Events should be defined in odemodel(). If defined in Xs(), events will be applied, but sensitivities will not be reset accordingly.")
  
  # Variable and parameter names
  variables <- attr(func, "variables")
  parameters <- attr(func, "parameters")
  forcnames <- attr(func, "forcings")
  
  # Variable and parameter names of sensitivities
  sensvar <- attr(extended, "variables")[!attr(extended, "variables") %in% variables]
  senssplit <- strsplit(sensvar, ".", fixed = TRUE)
  senssplit.1 <- unlist(lapply(senssplit, function(v) v[1]))
  senssplit.2 <- unlist(lapply(senssplit, function(v) paste(v[-1], collapse = ".")))
  svariables <- intersect(senssplit.2, variables)
  sparameters <- setdiff(senssplit.2, variables)
  senspars <- c(svariables, sparameters)
  
  # Initial values for sensitivities
  yiniSens <- as.numeric(senssplit.1 == senssplit.2)
  names(yiniSens) <- sensvar
  
  # Only a subset of all variables/forcings is returned
  if (is.null(names)) names <- c(variables, forcnames)
  
  # Controls to be modified from outside
  controls <- list(
    forcings = myforcings,
    events = myevents,
    names = names,
    optionsOde = optionsOde,
    optionsSens = optionsSens,
    fcontrol = myfcontrol
  )
  
  P2X <- function(times, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE) {

    if (deriv2)
      stop("Xs.deSolve: second-order sensitivities require backend = 'cppDE'.")

    fixedNames <- names(fixed)
    params <- c(unclass(pars), unclass(fixed))
    yini <- params[variables]
    mypars <- params[parameters]
    
    forcings <- controls$forcings
    events <- controls$events
    optionsOde <- controls$optionsOde
    optionsSens <- controls$optionsSens
    fcontrol <- controls$fcontrol
    names <- controls$names
    
    # Add event time points (required by integrator) 
    times <- sort(union(unique(events$time), times))
    
    # Sort event time points
    if (!is.null(events)) events <- events[order(events$time), ]
    
    # Filter sensGrid by fixed parameters
    if (length(fixedNames)) senspars <- setdiff(senspars, fixedNames)
    
    myderivs <- NULL
    if (!deriv) {
      
      # Evaluate model without sensitivities
      if (!is.null(forcings)) forc <- setForcings(func, forcings) else forc <- NULL
      out <- suppressWarnings(do.call(odeC, c(list(y = unclass(yini), times = times, func = func, 
                                                   parms = mypars, forcings = forc, 
                                                   events = list(data = events), fcontrol = fcontrol), optionsOde)))
      out <- submatrix(out, cols = c("time", names))
      
    } else {
      
      # Evaluate extended model
      if (!is.null(forcings)) forc <- setForcings(extended, forcings) else forc <- NULL
      outSens <- suppressWarnings(do.call(odeC, c(list(y = c(unclass(yini), yiniSens), times = times, 
                                                       func = extended, parms = mypars, 
                                                       forcings = forc, fcontrol = fcontrol,
                                                       events = list(data = events)), optionsSens)))
      
      out <- submatrix(outSens, cols = c("time", names))

      # Forcings live in `names` so the value matrix can return them, but the
      # extended ODE only has sensitivities for state variables. Restrict
      # the deriv axis to states that are actually requested.
      svars <- intersect(names, variables)

      # Apply parameter transformation to the derivatives.
      # deSolve's outSens columns are laid out so that array() fills into
      # [time, variable, sensitivity] naturally under column-major.
      sensNames <- as.vector(outer(svars, senspars, paste, sep = "."))
      mysensitivities <- array(outSens[, sensNames],
                               dim = c(nrow(outSens), length(svars), length(senspars)))

      dP <- attr(pars, "deriv")
      if (!is.null(dP)) {
        dPsub <- dP[senspars, , drop = FALSE]
        if(any(rownames(dP) %in% senspars)) {
          myderivs <- mysensitivities %bmm% dPsub
          dimnames(myderivs) <- list(NULL, svars, colnames(dPsub))
        } else {
          myderivs <- NULL
        }
      } else {
        myderivs <- mysensitivities
        dimnames(myderivs) <- list(NULL, svars, senspars)
      }
      
    }
    
    prdframe(out, deriv = myderivs, parameters = c(pars, fixed))
    
  }
  
  attr(P2X, "parameters") <- c(variables, parameters)
  attr(P2X, "equations") <- as.eqnvec(attr(func, "equations"))
  attr(P2X, "forcings") <- forcings
  attr(P2X, "events") <- events
  attr(P2X, "modelname") <- func[1]
  attr(P2X, "compileInfo") <- attr(odemodel, "compileInfo")

  prdfn(P2X, c(variables, parameters), condition)

}

# A prepared cppDE batch handle names the shared object it resolved its entry
# point from, so it goes stale as soon as that object is gone: after a rename
# via `modelname<-`, or when a saved workspace is reopened somewhere the object
# was never built -- a cluster node, say. cppDE calls through
# `.Call(name, PACKAGE = dll)`, so a stale handle fails at the call rather than
# before it. A handle that resolved no symbol has nothing to go stale.
.batchHandleLive <- function(handle) {
  if (is.null(handle)) return(FALSE)
  dll <- handle$sym$dll
  if (is.null(dll)) return(TRUE)
  dll %in% names(getLoadedDLLs())
}

# The forcings as the cppDE solver takes them, one data.frame of time and value
# per forcing, from the data.frame the user gives.
.cppdeForcings <- function(forcings) {
  if (is.null(forcings)) return(NULL)
  if (!inherits(forcings, "data.frame"))
    stop("'forcings' must be a data.frame, data.table, or tibble")
  required_cols <- c("name", "time", "value")
  if (!all(required_cols %in% names(forcings)))
    stop("'forcings' must contain columns: ", paste(required_cols, collapse = ", "))
  if (is.factor(forcings$name)) forcings$name <- as.character(forcings$name)
  if (!is.character(forcings$name)) stop("'name' must be character or factor")
  if (!is.numeric(forcings$time)) stop("'time' must be numeric")
  if (!is.numeric(forcings$value)) stop("'value' must be numeric")
  if (anyNA(forcings[, required_cols])) stop("'forcings' contains NA values")
  split(forcings[, c("time", "value")], forcings$name)
}

# Solver options as given, merged over the defaults, warning about any name
# the solver does not know.
.cppdeOptions <- function(options, defaults, label) {
  options <- as.list(options)
  bad <- setdiff(names(options), names(defaults))
  if (length(bad) > 0)
    warning(sprintf("%s: Ignoring unknown option(s): %s", label, paste(bad, collapse = ", ")))
  modifyList(defaults, options)
}

# A cppDE leaf keeps its forcings and solver options in `controls` the way the
# user gave them, so a change made by controls<- reads like the same argument
# given to the constructor. The solver needs them derived: split by forcing,
# merged over the defaults. `derive` is redone only when its source changed,
# which identical() tells from a pointer comparison while the source stays
# the same object, so a solve pays nothing for it and a change is validated,
# and warned about, once.
.derivedControl <- function(derive) {
  cache <- new.env(parent = emptyenv())
  function(src) {
    if (!isTRUE(cache$set) || !identical(cache$src, src)) {
      cache$val <- derive(src)
      cache$src <- src
      cache$set <- TRUE
    }
    cache$val
  }
}

#' @export
#' @rdname Xs
Xs.cppDE <- function(odemodel, forcings = NULL, events = NULL, names = NULL, condition = NULL,
                      optionsOde = list(), optionsSens = list(),
                      optionsReverse = NULL, ...) {

  # Derived at every solve from what `controls` holds; derived here as well,
  # which checks the arguments and warns about unknown options when given.
  forcsOf <- .derivedControl(.cppdeForcings)
  forcsOf(forcings)

  if (!is.null(events)) {
    stop("Events should be passed to odemodel() when using backend = 'cppDE'")
  }

  optionsDefault <- list(atol = 1e-6, rtol = 1e-6, maxattemps = 50L, maxsteps = 1e6L,
                         hini = 0, roottol = 1e-6, maxroot = 1L,
                         onFailure = "stop", traceFile = NULL)
  sensDefault <- c(optionsDefault, list(sensErrCon = TRUE))
  odeOf  <- .derivedControl(function(o) .cppdeOptions(o, optionsDefault, "optionsOde"))
  sensOf <- .derivedControl(function(o) .cppdeOptions(o, sensDefault, "optionsSens"))
  odeOf(optionsOde); sensOf(optionsSens)

  func <- odemodel$func
  extended <- odemodel$extended
  extended2 <- odemodel$extended2
  reversed <- odemodel$reversed
  reversed2 <- odemodel$reversed2
  if (is.null(extended)) warning("Element 'extended' empty. ODE model does not contain sensitivities.")

  # Extract metadata
  paramNames <- c(attr(func, "variables"), attr(func, "parameters"))
  dim_names <- attr(func, "dimNames")
  inner_names <- c(attr(func, "variables"), attr(func, "parameters"))

  # Only a subset of all variables is returned
  forcNames <- attr(func, "forcings")
  outNames <- c(dim_names$variable, forcNames)
  if (is.null(names)) names <- outNames else names <- intersect(outNames, names)
  if (length(names) == 0) stop(paste("Valid names are:", paste(outNames, collapse = ", ")))

  # Controls to be modified from outside, as the user gave them. Every entry is
  # read when the leaf runs.
  controls <- list(
    forcings = forcings,
    names = names,
    optionsOde = optionsOde,
    optionsSens = optionsSens,
    optionsReverse = optionsReverse
  )

  has_deriv2 <- !is.null(extended2)
  has_reverse <- !is.null(reversed)
  has_reverse2 <- !is.null(reversed2)
  # CVODES holds its checkpoints inside the solver and runs the backward solve
  # under its own step-size control, so both the shared store and the checked
  # sweep are arrangements of the native backend alone. Checked when the
  # options are read, which catches a setting made later through controls<-.
  has_store <- has_reverse && !identical(attr(reversed, "backend"), "cvode")
  reverseOpts <- function() {
    o <- controls$optionsReverse
    if (!is.null(o) && !has_store)
      stop("Xs: 'optionsReverse' checks the native backward pass. It needs ",
           "odemodel(backend = \"cppDE\", derivMode = c(\"forward\", \"reverse\")); ",
           "the Sundials backend runs its own step-size control.", call. = FALSE)
    o
  }
  reverseOpts()

  # Checkpoints of the value pass, for the backward pass that replays the same
  # trajectory. Matched on the point they were taken at, which is what cppDE
  # fingerprints the store on, so one can never answer for another parameter.
  #
  # A store answers any number of backward sweeps, so a lookup leaves it in
  # place: a second sweep over the same trajectory replays it instead of
  # integrating again. The cache keeps the most recently used ones, at least
  # twice the widest batch the value pass has taken, and evicts the oldest.
  scache <- new.env(parent = emptyenv())
  scache$items <- list()
  scache$width <- 0L
  storeCap <- function() max(64L, 2L * scache$width)
  storeFind <- function(times, params) {
    it <- scache$items
    for (k in seq_along(it))
      if (identical(it[[k]]$times, times) && identical(it[[k]]$params, params))
        return(k)
    NA_integer_
  }
  storePut <- function(times, params, store) {
    if (is.null(store)) return(invisible(NULL))
    k <- storeFind(times, params)
    if (!is.na(k)) scache$items[[k]] <- NULL
    scache$items[[length(scache$items) + 1L]] <-
      list(times = times, params = params, store = store)
    excess <- length(scache$items) - storeCap()
    if (excess > 0L) scache$items <- scache$items[-seq_len(excess)]
    invisible(NULL)
  }
  storeGet <- function(times, params) {
    k <- storeFind(times, params)
    if (is.na(k)) return(NULL)
    e <- scache$items[[k]]
    # most recently used goes to the back
    scache$items[[k]] <- NULL
    scache$items[[length(scache$items) + 1L]] <- e
    e$store
  }
  # The checked sweep, under optionsReverse$refine.
  sweepCtl <- function() {
    o <- reverseOpts()
    if (!has_store || !isTRUE(o$refine)) return(NULL)
    cppDE::adjointControl(refine = TRUE, gradtol = o$gradtol)
  }

  # Marshalling shared by the single and the batched entry: Phi'(theta) as the
  # tangent (and Phi''(theta) as the Hessian) on the inner parameter rows.
  prep1 <- function(pars, fixed, deriv, deriv2) {
    out <- list(params = c(unclass(pars), unclass(fixed)),
                tangent = NULL, hessian = NULL)
    if (!deriv) return(out)
    deriv_in <- attr(pars, "deriv")
    if (is.null(deriv_in)) return(out)
    phi_rows <- inner_names
    # Positions, not names: a character row index costs a match() per use, and
    # this runs once per condition per evaluation.
    ridx <- match(phi_rows, rownames(deriv_in))
    hit  <- !is.na(ridx)
    src  <- ridx[hit]
    out$tangent <- matrix(0, length(phi_rows), ncol(deriv_in),
                          dimnames = list(phi_rows, colnames(deriv_in)))
    out$tangent[hit, ] <- deriv_in[src, , drop = FALSE]
    if (deriv2) {
      d2 <- attr(pars, "deriv2")
      if (!is.null(d2)) {
        out$hessian <- array(0, c(length(phi_rows), dim(d2)[2], dim(d2)[3]),
                             dimnames = c(list(phi_rows), dimnames(d2)[2:3]))
        s2 <- match(phi_rows, dimnames(d2)[[1]])
        h2 <- !is.na(s2)
        out$hessian[h2, , ] <- d2[s2[h2], , , drop = FALSE]
      }
    }
    out
  }

  # Subsetting to `names` copies the whole sensitivity array. When `names` is
  # the full variable set, the default, the copy buys nothing. Forcings come
  # after the states, as values without derivatives.
  assemble1 <- function(res, pars, fixed, deriv, deriv2) {
    nms <- intersect(controls$names, colnames(res$variable))
    fnms <- intersect(controls$names, forcNames)
    keepAll <- identical(nms, colnames(res$variable))
    out <- cbind(res$time,
                 if (keepAll) res$variable else submatrix(res$variable, cols = nms),
                 if (length(fnms))
                   cppDE::forcingValues(res$time, forcsOf(controls$forcings)[fnms]))
    colnames(out)[1] <- "time"
    dX <- dX2 <- NULL
    if (deriv) {
      dX <- if (keepAll) res$tangent else res$tangent[, nms, , drop = FALSE]
      if (deriv2 && !is.null(res$hessian))
        dX2 <- if (keepAll) res$hessian else res$hessian[, nms, , , drop = FALSE]
    }
    prdframe(out, deriv = dX, deriv2 = dX2, parameters = c(pars, fixed))
  }

  # Every solve reads its options here, so this is also where a reverse option
  # set after construction meets the check the constructor makes. `sec`, the
  # sensErrCon of a solve with tangents.
  solveOpts <- function(deriv, sec = NULL) {
    reverseOpts()
    o <- if (deriv) sensOf(controls$optionsSens) else odeOf(controls$optionsOde)
    out <- list(abstol = o$atol, reltol = o$rtol, maxattemps = o$maxattemps,
                maxsteps = o$maxsteps, hini = o$hini, roottol = o$roottol,
                maxroot = o$maxroot, onFailure = o$onFailure, traceFile = o$traceFile)
    if (!is.null(sec)) out$sensErrCon <- sec
    out
  }
  # sensErrCon of a forward solve.
  sensSec <- function() isTRUE(sensOf(controls$optionsSens)$sensErrCon)

  # Second order forward exists only where it was built, and only on cppDE:
  # Sundials provides first order both ways, deSolve forward alone. Saying which
  # rebuild would answer beats handing solveODE a NULL model.
  pickModel <- function(deriv, deriv2) {
    if (!deriv) return(func)
    if (!deriv2) return(extended)
    if (is.null(extended2))
      stop("an exact Hessian forward needs the forward-over-forward object; ",
           "rebuild via odemodel(..., deriv2 = TRUE), which needs ",
           "backend = \"cppDE\". Backwards it is ",
           "derivMode = \"forward-reverse\" instead.", call. = FALSE)
    extended2
  }

  P2X <- function(times, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE,
                  keepStore = FALSE) {

    if (deriv2 && !has_deriv2)
      stop("Xs.cppDE: model was compiled without deriv2; rebuild via odemodel(..., deriv2 = TRUE).")
    if (deriv2 && !deriv) deriv <- TRUE
    forcs <- forcsOf(controls$forcings)

    # The values of a reverse evaluation come off the reverse object itself, so
    # the checkpoints the backward pass needs are already there and the states
    # are integrated once instead of twice.
    if (keepStore && !deriv && has_store) {
      params <- c(unclass(pars), unclass(fixed))
      res <- do.call(cppDE::solveODE, c(
        list(reversed, times, params, fixed = NULL,
             forcings = forcs, keepStore = TRUE),
        solveOpts(FALSE)))
      storePut(times, params, res$store)
      return(assemble1(res, pars, fixed, FALSE, FALSE))
    }

    prep <- prep1(pars, fixed, deriv, deriv2)
    res <- do.call(cppDE::solveODE, c(
      list(pickModel(deriv, deriv2), times, prep$params,
           tangent = prep$tangent, hessian = prep$hessian, fixed = NULL,
           forcings = forcs), solveOpts(deriv, if (deriv) sensSec())))
    assemble1(res, pars, fixed, deriv, deriv2)

  }

  # Handle cache for the prepared-batch path. An optimiser re-solves the same
  # conditions thousands of times with new numbers; preparing them once keeps
  # cppDE's argument marshalling out of the parallel region's critical path.
  bcache <- new.env(parent = emptyenv())
  rcache <- new.env(parent = emptyenv())
  vcache <- new.env(parent = emptyenv())

  # A prepared handle for `model`, reused while `sig` holds; NULL where cppDE's
  # solveBatch() lacks an argument in `need`.
  batchHandle <- function(cache, sig, model, conds, opts, need = character(0)) {
    ns <- asNamespace("cppDE")
    prepFn  <- get0("prepareBatch", envir = ns, inherits = FALSE)
    solveFn <- get0("solveBatch",   envir = ns, inherits = FALSE)
    if (is.null(prepFn) || is.null(solveFn) ||
        !all(need %in% names(formals(solveFn)))) return(NULL)
    if (!identical(cache$sig, sig) || !.batchHandleLive(cache$handle)) {
      cache$handle <- do.call(prepFn, c(list(model, conditions = conds()), opts))
      cache$sig <- sig
    }
    cache$handle
  }

  # One batched cppDE call over all requests. Falls back to a loop when cppDE
  # predates the entry point. A requested trace goes through: the batch writes
  # one file per condition, using `traceFile` as a template.
  P2Xbatch <- function(times, parsList, fixedList, deriv, deriv2, cores,
                       keepStore = FALSE) {

    if (deriv2 && !has_deriv2)
      stop("Xs.cppDE: model was compiled without deriv2; rebuild via odemodel(..., deriv2 = TRUE).")
    if (deriv2 && !deriv) deriv <- TRUE
    forcs <- forcsOf(controls$forcings)

    n <- length(parsList)
    timesL <- if (is.list(times)) times else rep(list(times), n)

    if (keepStore && !deriv && has_store) {
      scache$width <- max(scache$width, n)
      o      <- solveOpts(FALSE)
      paramL <- lapply(seq_len(n), function(i)
        c(unclass(parsList[[i]]), unclass(fixedList[[i]])))
      conds  <- function() lapply(seq_len(n), function(i) list(
        times = timesL[[i]], parms = paramL[[i]],
        forcings = forcs, keepStore = TRUE))
      batch <- get0("solveODEBatch", envir = asNamespace("cppDE"),
                    inherits = FALSE)
      ob <- o[setdiff(names(o), c("traceFile", "onFailure"))]
      h <- if (is.null(batch)) NULL else
        batchHandle(rcache, list(times = timesL, forcings = forcs, opts = ob),
                    reversed, conds, ob)
      res <- if (is.null(batch))
        lapply(seq_len(n), function(i) do.call(cppDE::solveODE, c(
          list(reversed, timesL[[i]], paramL[[i]], fixed = NULL,
               forcings = forcs, keepStore = TRUE), o)))
      else if (!is.null(h))
        cppDE::solveBatch(h, parms = paramL, cores = cores,
                          traceFile = o$traceFile, onFailure = o$onFailure)
      else
        do.call(batch, c(list(reversed, conditions = conds(), cores = cores), o))
      for (i in seq_len(n))
        storePut(timesL[[i]], paramL[[i]], res[[i]]$store)
      return(lapply(seq_len(n), function(i)
        assemble1(res[[i]], parsList[[i]], fixedList[[i]], FALSE, FALSE)))
    }
    preps  <- lapply(seq_len(n), function(i)
      prep1(parsList[[i]], fixedList[[i]], deriv, deriv2))
    model <- pickModel(deriv, deriv2)
    o     <- solveOpts(deriv, if (deriv) sensSec())

    ns      <- asNamespace("cppDE")
    batch   <- get0("solveODEBatch", envir = ns, inherits = FALSE)
    prepFn  <- get0("prepareBatch",  envir = ns, inherits = FALSE)
    solveFn <- get0("solveBatch",    envir = ns, inherits = FALSE)
    obatch  <- o[setdiff(names(o), c("traceFile", "onFailure"))]

    mkConds <- function() lapply(seq_len(n), function(i)
      list(times = timesL[[i]], parms = preps[[i]]$params,
           tangent = preps[[i]]$tangent, hessian = preps[[i]]$hessian,
           forcings = forcs))

    res <- if (is.null(batch)) {
      lapply(seq_len(n), function(i) do.call(cppDE::solveODE, c(
        list(model, timesL[[i]], preps[[i]]$params,
             tangent = preps[[i]]$tangent, hessian = preps[[i]]$hessian,
             fixed = NULL, forcings = forcs), o)))
    } else if (!is.null(prepFn) && !is.null(solveFn)) {
      # The handle bakes in everything but the numbers, so it is only valid
      # while shapes and labels stay put. It names its entry point rather than
      # holding an address, so it survives a reload of the shared object -- but
      # not a rename or a move to another machine, which `sig` cannot see and
      # `.batchHandleLive()` therefore checks separately.
      sig <- list(model = as.character(model), times = timesL,
                  deriv = deriv, deriv2 = deriv2,
                  sens = lapply(preps, function(pr) dimnames(pr$tangent)),
                  forcings = forcs, opts = obatch)
      if (!identical(bcache$sig, sig) || !.batchHandleLive(bcache$handle)) {
        bcache$handle <- do.call(prepFn, c(list(model, conditions = mkConds()),
                                           obatch))
        bcache$sig <- sig
      }
      solveFn(bcache$handle,
              parms    = lapply(preps, `[[`, "params"),
              tangent  = lapply(preps, `[[`, "tangent"),
              hessian  = lapply(preps, `[[`, "hessian"),
              cores = cores, traceFile = o$traceFile, onFailure = o$onFailure)
    } else {
      do.call(batch, c(list(model, conditions = mkConds(), cores = cores), o))
    }

    lapply(seq_len(n), function(i)
      assemble1(res[[i]], parsList[[i]], fixedList[[i]], deriv, deriv2))

  }
  attr(P2X, "batchfn") <- P2Xbatch

  # ---- Reverse mode -------------------------------------------------------
  #
  # The cotangent is that of the prediction as this leaf returns it: one row
  # per output row, one column per name in `controls$names`. What comes back is
  # the cotangent of `pars`, which is where the chain above this leaf continues.
  #
  # The cotangent is widened to the model's own state set first, because the
  # solver answers on all of them and `names` may be a subset.
  #
  # A cotangent only exists after the chain above has been walked, so this is a
  # second call over the same trajectory. It integrates nothing where the value
  # pass left its checkpoints behind, and replays the recorded steps instead.
  #
  # `seeds = TRUE` reads the trailing axis as independent first-order seeds
  # rather than as directions, and answers all of them in one backward sweep:
  # cppDE integrates one adjoint per seed column on the same trajectory.
  P2Xvjp <- function(times, pars, fixed = NULL, cotangent, seeds = FALSE) {
    forcs <- forcsOf(controls$forcings)
    w <- .asCtOut(cotangent)
    K <- .ctK(w)
    if (seeds) {
      .requireReverse(has_reverse, has_reverse2, 1L)
      params <- c(unclass(pars), unclass(fixed))
      res <- do.call(cppDE::solveODE, c(
        list(reversed, times, params, fixed = NULL,
             forcings = forcs,
             cotangent = .widenSeeds(w, dim_names$variable, controls$names),
             store = if (has_store) storeGet(times, params)), solveOpts(FALSE)))
      .requireAdjoint(res)
      return(.pickCotangent(.adjointSeeds(res), names(pars)))
    }
    states <- dim_names$variable
    .requireReverse(has_reverse, has_reverse2, K)
    ct <- .widenCotangent(w, states, controls$names)
    # Second order sweeps over tangents, so it takes the directions the value
    # pass propagated, the same tangent prep1 builds for a forward solve. A
    # store cannot help it: a checkpoint's tangents do not outlive the solve
    # that took them, so the sweep integrates its own.
    pr <- prep1(pars, fixed, K > 1L, FALSE)
    o <- solveOpts(K > 1L)
    # The store and the control follow the value pass, which second order does
    # not share: cppDE refuses a store under forward-reverse.
    st <- if (K > 1L) NULL else storeGet(times, pr$params)
    call <- c(list(reversed, times, pr$params, fixed = NULL,
                   forcings = forcs, cotangent = ct$cotangent),
              if (K > 1L) NULL else
                list(store = st,
                     adjoint = sweepCtl()))
    if (K > 1L) {
      if (is.null(pr$tangent))
        stop("a second-order cotangent needs the tangents the value pass ",
             "propagated, and this input has none.", call. = FALSE)
      call[[1L]] <- reversed2
      call$tangent <- pr$tangent
      call$curvature <- ct$curvature
    }
    res <- do.call(cppDE::solveODE, c(call, o))
    .requireAdjoint(res)
    .pickCotangent(.adjointCt(res, K), names(pars))
  }
  attr(P2X, "vjpfn") <- P2Xvjp

  # Every condition's backward solve in one call, the way P2Xbatch does the
  # forward ones. Falls back to a loop where cppDE predates the entry point.
  P2Xvjpbatch <- function(times, parsList, fixedList, cotangentList, conditions,
                          cores, seeds = FALSE) {
    forcs <- forcsOf(controls$forcings)
    wList <- lapply(cotangentList, .asCtOut)
    K <- max(vapply(wList, .ctK, 1L))
    n <- length(parsList)
    timesL <- if (is.list(times)) times else rep(list(times), n)
    states <- dim_names$variable
    batch <- get0("solveODEBatch", envir = asNamespace("cppDE"), inherits = FALSE)
    condOf <- function(i)
      if (is.null(conditions)) NULL else conditions[[i]]

    # Independent first-order seeds, each request with its own number of them.
    if (seeds) {
      .requireReverse(has_reverse, has_reverse2, 1L)
      o <- solveOpts(FALSE)
      conds <- lapply(seq_len(n), function(i) {
        params <- c(unclass(parsList[[i]]), unclass(fixedList[[i]]))
        cd <- list(times = timesL[[i]], parms = params,
                   forcings = forcs,
                   cotangent = .widenSeeds(wList[[i]], states, controls$names))
        if (has_store) cd$store <- storeGet(timesL[[i]], params)
        cd
      })
      res <- if (is.null(batch))
        lapply(conds, function(a)
          do.call(cppDE::solveODE,
                  c(list(reversed, a$times, a$parms, fixed = NULL,
                         forcings = a$forcings, cotangent = a$cotangent,
                         store = a$store), o)))
      else
        do.call(batch, c(list(reversed, conditions = conds, cores = cores), o))
      return(lapply(seq_len(n), function(i) {
        .requireAdjoint(res[[i]], condOf(i))
        .pickCotangent(.adjointSeeds(res[[i]]), names(parsList[[i]]))
      }))
    }

    .requireReverse(has_reverse, has_reverse2, K)
    o <- solveOpts(K > 1L)
    conds <- lapply(seq_len(n), function(i) {
      pr <- prep1(parsList[[i]], fixedList[[i]], K > 1L, FALSE)
      ct <- .widenCotangent(wList[[i]], states, controls$names)
      cd <- list(times = timesL[[i]],
                 parms = pr$params,
                 forcings = forcs,
                 cotangent = ct$cotangent)
      if (K > 1L) {
        if (is.null(pr$tangent))
          stop("a second-order cotangent needs the tangents the value pass ",
               "propagated, and condition ", i, " has none.", call. = FALSE)
        cd$tangent <- pr$tangent
        cd$curvature <- ct$curvature
        return(cd)
      }
      # Only the store is first order only: cppDE refuses one under
      # forward-reverse, because a checkpoint's tangents do not outlive the
      # solve that took them.
      cd$store <- storeGet(timesL[[i]], pr$params)
      cd$adjoint <- sweepCtl()
      cd
    })

    model <- if (K > 1L) reversed2 else reversed
    # First order without a control reuses a prepared handle.
    h <- NULL
    if (!is.null(batch) && K == 1L &&
        all(vapply(conds, function(cd) is.null(cd$adjoint), TRUE))) {
      ob <- o[setdiff(names(o), c("traceFile", "onFailure"))]
      sig <- list(times = timesL, forcings = forcs, opts = ob,
                  ct = lapply(conds, function(cd) dim(cd$cotangent)))
      h <- batchHandle(vcache, sig, model, function() conds, ob, need = "store")
    }
    res <- if (!is.null(h))
      cppDE::solveBatch(h, parms = lapply(conds, `[[`, "parms"),
                        cotangent = lapply(conds, `[[`, "cotangent"),
                        store = lapply(conds, function(cd) cd$store),
                        cores = cores, traceFile = o$traceFile,
                        onFailure = o$onFailure)
    else if (is.null(batch))
      lapply(seq_len(n), function(i) {
        a <- conds[[i]]
        do.call(cppDE::solveODE,
                c(list(model, a$times, a$parms, fixed = NULL,
                       forcings = forcs, cotangent = a$cotangent),
                  if (K > 1L) list(tangent = a$tangent, curvature = a$curvature)
                  else list(store = a$store, adjoint = a$adjoint),
                  o))
      })
    else
      do.call(batch, c(list(model, conditions = conds, cores = cores), o))

    lapply(seq_len(n), function(i) {
      .requireAdjoint(res[[i]], condOf(i))
      .pickCotangent(.adjointCt(res[[i]], .ctK(wList[[i]])),
                     names(parsList[[i]]))
    })
  }
  attr(P2X, "vjpbatchfn") <- P2Xvjpbatch
  # Says the value pass of a reverse evaluation should run on the reverse
  # object, so the backward pass finds the checkpoints already taken.
  attr(P2X, "keepstore") <- has_store

  attr(P2X, "parameters") <- paramNames
  attr(P2X, "equations") <- as.eqnvec(attr(func, "equations"))
  attr(P2X, "forcings") <- forcings
  attr(P2X, "events") <- events
  attr(P2X, "modelname") <- func[1]
  attr(P2X, "compileInfo") <- attr(odemodel, "compileInfo")

  prdfn(P2X, paramNames, condition)
}


#' Prediction function of an ODE model without sensitivities
#'
#' @description Turns an [odemodel()] into a prediction function that
#' integrates the states alone. States missing from `pars` start at 0, so
#' further quantities, such as fluxes, can be integrated by adding them to the
#' equations.
#'
#' @param odemodel An [odemodel()]. Its backend decides which method runs.
#' @param forcings,events,condition As in [Xs()].
#' @param optionsOde Named list of solver options, as `optionsOde` of [Xs()].
#' @param fcontrol `deSolve` backend only, as in [Xs()].
#' @param ... Not used.
#'
#' @return A [prdfn], called as `x(times, pars, fixed = NULL)`, returning a
#'   [prdlist] without derivatives. `forcings` and `optionsOde` can be changed
#'   with [controls()].
#' @seealso [Xs()]
#' @export
Xf <- function(odemodel, ...) {
  UseMethod("Xf", odemodel)
}

#' @export
#' @rdname Xf
Xf.deSolve <- function(odemodel, forcings = NULL, events = NULL, condition = NULL,
                       optionsOde = list(method = "lsoda"), fcontrol = NULL, ...) {

  func <- odemodel$func

  myforcings <- forcings
  myevents <- events
  myfcontrol <- fcontrol

  variables <- attr(func, "variables")
  parameters <- attr(func, "parameters")
  yini <- rep(0, length(variables))
  names(yini) <- variables

  controls <- list(
    forcings = myforcings,
    events = myevents,
    optionsOde = optionsOde,
    fcontrol = myfcontrol
  )

  P2X <- function(times, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE) {

    if (deriv2)
      stop("Xf: second-order sensitivities are not implemented (use Xs() for deriv2).")

    events <- controls$events
    forcings <- controls$forcings
    optionsOde <- controls$optionsOde
    fcontrol <- controls$fcontrol

    # Xf has no sensitivities, so fixed/free collapses to a single pars vector.
    pars <- c(unclass(pars), unclass(fixed))

    times <- sort(union(unique(events$time), times))

    yini[names(pars[names(pars) %in% variables])] <- pars[names(pars) %in% variables]
    mypars <- pars[parameters]

    if (!is.null(forcings)) forc <- setForcings(func, forcings) else forc <- NULL
    out <- suppressWarnings(do.call(odeC, c(list(y = yini, times = times, func = func,
                                                 parms = mypars, forcings = forc,
                                                 events = list(data = events),
                                                 fcontrol = fcontrol), optionsOde)))

    prdframe(out, deriv = NULL, parameters = pars)
  }

  attr(P2X, "parameters") <- c(variables, parameters)
  attr(P2X, "equations") <- as.eqnvec(attr(func, "equations"))
  attr(P2X, "forcings") <- forcings
  attr(P2X, "events") <- events
  attr(P2X, "modelname") <- func[1]
  attr(P2X, "compileInfo") <- attr(odemodel, "compileInfo")

  prdfn(P2X, c(variables, parameters), condition)
}

#' @export
#' @rdname Xf
# Xf is the no-derivative prediction, so it has no vjp either -- not an
# omission, the point of it. A chain that needs a gradient uses Xs().
Xf.cppDE <- function(odemodel, forcings = NULL, events = NULL, condition = NULL,
                      optionsOde = list(), ...) {

  # `controls` keeps the arguments as given; the solver's forms are derived
  # from them at every solve, and here, which checks them once when given.
  forcsOf <- .derivedControl(.cppdeForcings)
  forcsOf(forcings)

  if (!is.null(events))
    stop("Events must be passed to odemodel() for backend = 'cppDE' / 'Sundials'.")

  optionsDefault <- list(atol = 1e-6, rtol = 1e-6, maxattemps = 50L, maxsteps = 1e6L,
                         hini = 0, roottol = 1e-6, maxroot = 1L,
                         onFailure = "stop", traceFile = NULL)
  odeOf <- .derivedControl(function(o) .cppdeOptions(o, optionsDefault, "optionsOde"))
  odeOf(optionsOde)

  func <- odemodel$func
  paramNames <- c(attr(func, "variables"), attr(func, "parameters"))

  controls <- list(forcings = forcings, optionsOde = optionsOde)

  P2X <- function(times, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE) {

    if (deriv2)
      stop("Xf: second-order sensitivities are not implemented (use Xs() for deriv2).")

    params <- c(unclass(pars), unclass(fixed))
    unset <- setdiff(attr(func, "variables"), names(params))
    params[unset] <- 0
    forcings <- forcsOf(controls$forcings)
    optionsOde <- odeOf(controls$optionsOde)

    out <- cppDE::solveODE(func, times, params,
                            tangent = NULL, hessian = NULL, fixed = NULL,
                            forcings = forcings,
                            abstol = optionsOde$atol, reltol = optionsOde$rtol,
                            maxattemps = optionsOde$maxattemps,
                            maxsteps = optionsOde$maxsteps,
                            hini = optionsOde$hini,
                            roottol = optionsOde$roottol,
                            maxroot = optionsOde$maxroot,
                            onFailure = optionsOde$onFailure,
                            traceFile = optionsOde$traceFile)

    out <- cbind(out$time, out$variable,
                 if (length(forcings)) cppDE::forcingValues(out$time, forcings))
    colnames(out)[1] <- "time"

    prdframe(out, deriv = NULL, parameters = c(pars, fixed))
  }

  attr(P2X, "parameters") <- paramNames
  attr(P2X, "equations") <- as.eqnvec(attr(func, "equations"))
  attr(P2X, "forcings") <- forcings
  attr(P2X, "events") <- events
  attr(P2X, "modelname") <- func[1]
  attr(P2X, "compileInfo") <- attr(odemodel, "compileInfo")

  prdfn(P2X, paramNames, condition)
}


#' Model prediction function from data.frame
#' 
#' @param data data.frame with columns "name", "time", and row names that 
#' are taken as parameter names. The data frame can contain a column "value"
#' to initialize the parameters.
#' @param condition either NULL (generic prediction for any condition) or a character, denoting
#' the condition for which the function makes a prediction.
#' @return A [prdfn], called as `x(times, pars, fixed = NULL, deriv = TRUE)`,
#'   that interpolates the parameter values linearly in time. Its parameters
#'   are the row names of `data`. Second-order derivatives are not available.
#' @examples
#' # Generate a data.frame and corresponding prediction function
#' timesD <- seq(0, 2*pi, 0.5)
#' mydata <- data.frame(name = "A", time = timesD, value = sin(timesD), 
#'                      row.names = paste0("par", 1:length(timesD)))
#' x <- Xd(mydata)
#' 
#' # Evaluate the prediction function at different time points
#' times <- seq(0, 2*pi, 0.01)
#' pouter <- structure(mydata$value, names = rownames(mydata))
#' prediction <- x(times, pouter)
#' plot(prediction)
#' @export
Xd <- function(data, condition = NULL) {
  
  states <- unique(as.character(data$name))
  
  
  # List of prediction functions with sensitivities
  predL <- lapply(states, function(s) {
    subdata <- subset(data, as.character(name) == s)
    
    M <- diag(1, nrow(subdata), nrow(subdata))
    parameters.specific <- rownames(subdata)
    if(is.null(parameters.specific)) parameters.specific <- paste("par", s, 1:nrow(subdata), sep = "_")
    sensnames <- paste(s, parameters.specific, sep = ".")
    
    # return function
    out <- function(times, pars) {
      value <- approx(x = subdata$time, y = pars[parameters.specific], xout = times, rule = 2)$y
      grad <- do.call(cbind, lapply(1:nrow(subdata), function(i) {
        approx(x = subdata$time, y = M[, i], xout = times, rule = 2)$y
      }))
      colnames(grad) <- sensnames
      attr(value, "sensitivities") <- grad
      attr(value, "sensnames") <- sensnames
      return(value)
    }
    
    attr(out, "parameters") <- parameters.specific
    
    return(out)
    
  }); names(predL) <- states
  
  # Collect parameters
  parameters <- unlist(lapply(predL, function(p) attr(p, "parameters")))
  
  # Initialize parameters if available
  pouter <- NULL
  if(any(colnames(data) == "value")) 
    pouter <- structure(data$value[match(parameters, rownames(data))], names = parameters)
  
  sensGrid <- expand.grid(states, parameters, stringsAsFactors=FALSE)
  sensNames <- paste(sensGrid[,1], sensGrid[,2], sep=".")  
  
  
  controls <- list()  
  
  P2X <- function(times, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE){

    if (deriv2)
      stop("Xd: second-order sensitivities are not implemented for data-driven prediction.")
    # `fixed` is accepted for prdfn-wrapper symmetry; Xd is purely
    # data-grid-driven, so any fixed parameters are merged into `pars`
    # for the lookup.
    if (!is.null(fixed)) pars <- c(unclass(pars), unclass(fixed))


    predictions <- lapply(states, function(s) predL[[s]](times, pars)); names(predictions) <- states
    
    out <- cbind(times, do.call(cbind, predictions))
    colnames(out) <- c("time", states)
    
    myderivs <- NULL
    if (deriv) {

      # Fill in sensitivities, column layout is state-fastest, matching the
      # expand.grid(states, parameters) ordering used to build sensNames.
      outSens <- matrix(0, nrow = length(times), ncol = length(sensNames),
                        dimnames = list(NULL, sensNames))
      for (s in states) {
        mysens   <- attr(predictions[[s]], "sensitivities")
        mynames  <- attr(predictions[[s]], "sensnames")
        outSens[, mynames] <- mysens
      }

      # Reshape to 3D [time, state, param] (batch-first), matching Xs.
      myderivs <- array(outSens,
                        dim = c(length(times), length(states), length(parameters)),
                        dimnames = list(NULL, states, parameters))

      # Chain rule via upstream parameter transformation.
      dP <- attr(pars, "deriv")
      if (!is.null(dP)) {
        dPsub <- dP[parameters, , drop = FALSE]
        myderivs <- myderivs %bmm% dPsub
        dimnames(myderivs) <- list(NULL, states, colnames(dPsub))
      }
    }

    prdframe(out, deriv = myderivs, parameters = pars)
    
  }
  
  # The interpolation is linear in the parameters it reads, and its Jacobian is
  # the same `grad` the forward path builds; contracting it with the cotangent
  # rather than with dP is the whole difference.
  attr(P2X, "vjpfn") <- function(times, pars, fixed = NULL, cotangent) {
    w <- .asCtOut(cotangent)
    p <- if (is.null(fixed)) pars else c(unclass(pars), unclass(fixed))
    K <- .ctK(w)
    out <- .ctZero(names(pars), K)
    for (s in states) {
      if (!(s %in% dimnames(w)[[2L]])) next
      pr  <- predL[[s]](times, p)
      sens <- attr(pr, "sensitivities")
      nms  <- attr(pr, "parameters") %||% attr(predL[[s]], "parameters")
      ws <- matrix(w[, s, ], nrow = dim(w)[1L], ncol = K)
      contrib <- crossprod(sens, ws)
      rownames(contrib) <- nms
      hit <- intersect(nms, rownames(out))
      if (length(hit))
        out[hit, ] <- out[hit, , drop = FALSE] + contrib[hit, , drop = FALSE]
    }
    out
  }

  attr(P2X, "parameters") <- structure(parameters, names = NULL)
  attr(P2X, "pouter") <- pouter
  
  prdfn(P2X, attr(P2X, "parameters"), condition)
  
}


#' Observation functions
#'
#' @description 
#' Creates an object of type [obsfn] that evaluates an observation function
#' and, if requested, its first and second derivatives based on the output of a model 
#' prediction function, see [prdfn], as e.g. produced by [Xs].
#' 
#' @param g Named character vector or [eqnvec] defining the observation
#' function, or a list of those named by condition, which builds one
#' observation function per condition in a single call.
#' @param f Equations, anything [as.eqnvec()] accepts, or a prediction
#'   function. The states and parameters are read off `f`; for a prediction
#'   function the forcings count as states.
#' @param states Character vector of states, added to those read off `f`.
#' @param parameters Character vector of parameters, added to the symbols of
#'   `g` and `f` that are neither states nor `time`.
#' @param condition Either `NULL` (generic prediction for any condition) or a character 
#' string specifying the condition for which the function generates predictions.
#' @param attach.input Logical, append the incoming states to the output.
#'   Defaults to `FALSE`: passing the state sensitivities through costs a
#'   full array copy per condition and objective functions do not read them.
#'   Set `TRUE` to plot states alongside observables. Can be changed with
#'   [controls()].
#' @param compile Logical. `TRUE` compiles the generated code now, `FALSE`
#'   leaves it to [compile()]. The function is evaluable only compiled.
#' @param modelname Character, the base name of the generated C++ file and
#'   shared object, default `"obsfn"`, followed by `_<condition>`.
#' @param verbose Logical, print compiler output to the R console.
#' @param cores Number of parallel jobs used to generate the sources when
#' `g` is a list; `NULL` auto-detects. Ignored for a single observation
#' function, which is one source either way.
#' @param derivMode Which derivative products to build, any of `"forward"`
#'   (default), `"reverse"` and `"forward-reverse"`.
#'   * `"forward"`: the Jacobian by forward-mode AD. Without it the function
#'     returns no `deriv` attribute.
#'   * `"reverse"`: the vector-Jacobian product the reverse sweep contracts
#'     against, a second instantiation of the expression body. It is what
#'     `obj(..., sweep = "reverse")` needs from an observation function.
#'   * `"forward-reverse"`: its derivative along a tangent, what the reverse
#'     sweep needs with `deriv2 = TRUE`.
#'
#'   Either way the observation function is evaluable only after compilation.
#' @param deriv Logical. If `TRUE` (default), attach the first-order
#'   sensitivity `attr(., "deriv")` of shape `[time, observable, theta]`.
#' @param deriv2 Logical. If `TRUE`, attach a second-order derivative
#'   `attr(., "deriv2")` array of shape `[time, observable, theta, theta]`.
#'   Requires `deriv = TRUE`. Default `FALSE`.
#' @param outdir Character. Directory for the generated C++ source and the
#'   compiled shared object. Defaults to the working directory.
#'
#' @return
#' An [obsfn], called as `g(out, pars, fixed = NULL, deriv = TRUE,
#' deriv2 = FALSE, conditions, env = NULL)` with a [prdlist] `out`, or
#' composed as `g * x`. It returns a [prdlist] of the observables and their
#' derivatives.
#' 
#' @example inst/examples/prediction.R
#' 
#' @importFrom cppDE cppFUN
#' @importFrom abind abind
#' @export
Y <- function(g, f = NULL, states = NULL, parameters = NULL,
              condition = NULL, attach.input = FALSE,
              compile = FALSE, modelname = NULL, verbose = FALSE,
              cores = NULL, deriv = TRUE, deriv2 = FALSE,
              derivMode = "forward", outdir = getwd()) {

  derivMode <- .matchDerivMode(derivMode, c("forward", "reverse", "forward-reverse"))

  # A named list of observable sets builds one obsfn per condition, generated
  # in parallel and compiled once, the way `P()` handles a trafo list.
  if (is.list(g) && !inherits(g, "eqnvec")) {
    if (is.null(names(g)) || any(!nzchar(names(g))))
      stop("Y(): a list of observables needs the conditions as its names.",
           call. = FALSE)
    cores <- if (Sys.info()[["sysname"]] == "Windows") 1L
             else if (is.null(cores)) detectFreeCores()
             else min(detectFreeCores(), cores)
    result <- Reduce("+", parallel::mclapply(seq_along(g), function(i)
      Y(g[[i]], f = f, states = states, parameters = parameters,
        condition = names(g)[i], attach.input = attach.input,
        compile = FALSE, modelname = modelname, verbose = verbose,
        deriv = deriv, deriv2 = deriv2, derivMode = derivMode, outdir = outdir),
      mc.cores = cores))
    if (compile)
      compile(result, cores = cores, output = modelname, verbose = verbose)
    return(result)
  }

  emit_d1 <- isTRUE(deriv)
  emit_d2 <- isTRUE(deriv2)
  if (emit_d2 && !emit_d1)
    stop("Y(deriv2 = TRUE) requires deriv = TRUE.", call. = FALSE)

  if (is.null(f) && is.null(states) && is.null(parameters))
    stop("Not all three arguments f, states and parameters can be NULL")

  # Define model name with condition suffix
  if (is.null(modelname)) modelname <- "obsfn"
  if (!is.null(condition)) modelname <- paste(modelname, sanitizeConditions(condition), sep = "_")

  # Identify symbols in g
  symbols <- getSymbols(unclass(g))

  # Infer states and parameters
  if (is.null(f)) {
    states <- union(states, "time")
    parameters <- union(parameters, setdiff(symbols, states))
  } else if (inherits(f, "fn")) {
    # From the leaves rather than the mappings: a composed mapping holds a
    # copy taken when it was composed, the leaf the forcings it runs with.
    myforcings <- Reduce(union, lapply(.fnLeaves(f), function(l)
      as.character(.kernelSetting(l$kernel, "forcings")$name)))
    mystates <- unique(c(do.call(c, lapply(getEquations(f), names)), "time"))
    if (length(intersect(myforcings, mystates)) > 0)
      stop("Forcings and states overlap in different conditions.")

    mystates <- c(mystates, myforcings)
    myparameters <- setdiff(union(getParameters(f), getSymbols(unclass(g))),
                            c(mystates, myforcings))
    states <- union(mystates, states)
    parameters <- union(myparameters, parameters)
  } else {
    f <- as.eqnvec(f)
    mystates <- union(names(f), "time")
    myparameters <- getSymbols(c(unclass(g), unclass(f)), exclude = mystates)
    states <- union(mystates, states)
    parameters <- union(myparameters, parameters)
  }

  observables <- names(g)
  obsParams <- intersect(symbols, parameters)
  obsStates <- setdiff(symbols, parameters)

  # Compile evaluator for g (value, Jacobian, Hessian, AD chain)
  gEval <- suppressWarnings(
    cppDE::cppFUN(
      unclass(g),
      variables  = obsStates,
      parameters = obsParams,
      compile    = compile,
      modelname  = modelname,
      outdir     = outdir,
      verbose    = verbose,
      convenient = FALSE,
      derivMode  = derivMode,
      deriv      = emit_d1,
      deriv2     = emit_d2
    )
  )

  gfun       <- gEval$func
  gjac       <- gEval$jac
  ghess      <- gEval$hess
  gevaluate  <- gEval$evaluate
  use_ad     <- "forward" %in% derivMode

  controls <- list(attach.input = attach.input)

  # Core observation mapping function
  # `.ad_out` lets the batched entry hand in a precomputed AD result; the rest
  # of the assembly is identical.
  X2Y <- function(out, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE,
                  .ad_out = NULL, .fixedObs = NULL) {

    if (deriv2 && !emit_d2)
      stop("Y() was built with deriv2 = FALSE; rebuild Y() with deriv2 = TRUE.", call. = FALSE)
    if (!emit_d1) deriv <- FALSE
    if (deriv2 && !deriv) deriv <- TRUE

    attach.input <- controls$attach.input
    fixedObsParams <- if (!is.null(.fixedObs)) .fixedObs else
      intersect(union(attr(pars, "fixed"), names(fixed)), obsParams)
    params <- c(unclass(pars), unclass(fixed))

    if (use_ad && deriv && !is.null(gevaluate)) {
      # AD path: evaluate() returns y and its tangent, chain-ruled via dX/dP.
      dX_full <- attr(out, "deriv")
      dX2_full <- attr(out, "deriv2")
      # If no obsStates appear in dX's state dim, it contains no upstream state
      # sensitivity for this observation function. Suppress it to avoid
      # spurious theta mismatches against dP (e.g. Xt() returns a deriv array
      # with unrelated layout).
      dX <- dX_full
      if (!is.null(dX) && !any(match(obsStates, dimnames(dX)[[2]], 0L) > 0L))
        dX <- NULL
      dX2 <- dX2_full
      if (!is.null(dX2) && !any(match(obsStates, dimnames(dX2)[[2]], 0L) > 0L))
        dX2 <- NULL
      ad_out <- if (!is.null(.ad_out)) .ad_out else
        gevaluate(out[, obsStates, drop = FALSE], params[obsParams],
                  tangentX = dX, tangentP = attr(pars, "deriv"),
                  hessianX = dX2, hessianP = attr(pars, "deriv2"),
                  deriv2 = deriv2,
                  attach.input = attach.input, fixed = fixedObsParams)
      # Values: evaluate() returns observables (and pass-through extras when
      # attach.input = TRUE) under attach.input semantics matching gfun.
      gAll <- ad_out$y
      # A NaN stays in the prediction: a ratio of states that all start at 0
      # is undefined at t0 only, where no data may sit (Laske_PLOSComputBiol2019).
      # normL2 stops on a NaN at a data point.
      gVal <- gAll[, observables, drop = FALSE]
      values <- cbind(time = out[, "time"], gVal)
      if (attach.input) values <- cbind(values, submatrix(out, cols = -1))
      myderivs <- ad_out$tangent
      myderivs2 <- if (deriv2) ad_out$hessian else NULL
      # Append pass-through state sensitivities for states that are attached
      # but not consumed by the observables; the AD path only emits sensitivities
      # for obsStates and would otherwise leave those rows missing.
      if (attach.input && !is.null(myderivs) && !is.null(dX_full)) {
        theta <- dimnames(myderivs)[[3]]
        outer_theta <- theta %||% dimnames(dX_full)[[3]]
        missing <- setdiff(outer_theta, dimnames(dX_full)[[3]])
        if (length(missing)) dX_full <- abind::abind(dX_full, array(0, c(dim(dX_full)[1], dim(dX_full)[2], length(missing)), dimnames = list(NULL, NULL, missing)), along = 3)
        already <- intersect(dimnames(myderivs)[[2]], dimnames(dX_full)[[2]])
        add_states <- setdiff(dimnames(dX_full)[[2]], already)
        if (length(add_states))
          myderivs <- abind::abind(myderivs, dX_full[, add_states, outer_theta, drop = FALSE], along = 2)
      }
      if (attach.input && !is.null(myderivs2) && !is.null(dX2_full)) {
        theta <- dimnames(myderivs2)[[3]]
        outer_theta <- theta %||% dimnames(dX2_full)[[3]]
        missing <- setdiff(outer_theta, dimnames(dX2_full)[[3]])
        if (length(missing)) {
          # Pad dim 3 (theta1) with zero blocks; keep dim 4 matching existing dX2_full.
          dX2_full <- abind::abind(dX2_full,
                                   array(0, c(dim(dX2_full)[1], dim(dX2_full)[2], length(missing), dim(dX2_full)[4]),
                                         dimnames = list(NULL, NULL, missing, dimnames(dX2_full)[[4]])),
                                   along = 3)
          # Pad dim 4 (theta2) with zero blocks; dim 3 now includes the missing entries.
          dX2_full <- abind::abind(dX2_full,
                                   array(0, c(dim(dX2_full)[1], dim(dX2_full)[2], dim(dX2_full)[3], length(missing)),
                                         dimnames = list(NULL, NULL, dimnames(dX2_full)[[3]], missing)),
                                   along = 4)
        }
        already <- intersect(dimnames(myderivs2)[[2]], dimnames(dX2_full)[[2]])
        add_states <- setdiff(dimnames(dX2_full)[[2]], already)
        if (length(add_states))
          myderivs2 <- abind::abind(myderivs2, dX2_full[, add_states, outer_theta, outer_theta, drop = FALSE], along = 2)
      }
    } else {
      # Values only (reverse-only build).
      gVal <- gfun(out[, obsStates, drop = FALSE], params[obsParams], attach.input, fixedObsParams)[, observables, drop = FALSE]
      values <- cbind(time = out[, "time"], gVal)
      if (attach.input) values <- cbind(values, submatrix(out, cols = -1))
      myderivs <- NULL
      myderivs2 <- NULL
    }

    prdframe(prediction = values, deriv = myderivs, deriv2 = myderivs2, parameters = c(pars, fixed))
  }

  # One evaluateBatch over all requests; the surrounding assembly stays in R.
  # The value-only path still loops.
  X2Ybatch <- function(outList, parsList, fixedList, deriv, deriv2, cores) {
    n <- length(parsList)
    loop <- function() lapply(seq_len(n), function(i)
      X2Y(outList[[i]], parsList[[i]], fixedList[[i]], deriv, deriv2))

    eb <- gEval$evaluateBatch
    if (is.null(eb) || !use_ad || !deriv || !emit_d1 ||
        (deriv2 && !emit_d2) || is.null(gevaluate)) return(loop())

    # `obsStates`/`obsParams` are constants and the incoming name sets are the
    # same for most conditions; reuse the last result instead of redoing
    # the set algebra n times.
    sRef <- NULL; sKeep <- TRUE; fRef <- NULL; fVal <- NULL
    s2Ref <- NULL; s2Keep <- TRUE
    sets <- lapply(seq_len(n), function(i) {
      pars <- parsList[[i]]; out <- outList[[i]]
      fnm <- names(fixedList[[i]])
      params <- c(unclass(pars), unclass(fixedList[[i]]))
      dX <- attr(out, "deriv")
      if (!is.null(dX)) {
        dn <- dimnames(dX)[[2]]
        if (!identical(dn, sRef)) {
          sRef <<- dn; sKeep <<- length(intersect(obsStates, dn)) > 0
        }
        if (!sKeep) dX <- NULL
      }
      dX2 <- if (deriv2) attr(out, "deriv2") else NULL
      if (!is.null(dX2)) {
        dn2 <- dimnames(dX2)[[2]]
        if (!identical(dn2, s2Ref)) {
          s2Ref <<- dn2; s2Keep <<- length(intersect(obsStates, dn2)) > 0
        }
        if (!s2Keep) dX2 <- NULL
      }
      key <- c(attr(pars, "fixed"), fnm)
      if (!identical(key, fRef)) {
        fRef <<- key
        fVal <<- intersect(union(attr(pars, "fixed"), fnm), obsParams)
      }
      list(vars = out[, obsStates, drop = FALSE], params = params[obsParams],
           tangentX = dX, tangentP = attr(pars, "deriv"), hessianX = dX2,
           hessianP = if (deriv2) attr(pars, "deriv2") else NULL,
           attach.input = controls$attach.input, fixed = fVal)
    })
    ad <- eb(sets, cores = cores, deriv2 = deriv2)
    lapply(seq_len(n), function(i)
      X2Y(outList[[i]], parsList[[i]], fixedList[[i]], deriv, deriv2,
          .ad_out = ad[[i]], .fixedObs = sets[[i]]$fixed))
  }
  # ---- Reverse mode -------------------------------------------------------
  #
  # The cotangent is that of the observables this leaf returns; back come the
  # cotangents of the two things it read, the prediction's states and its own
  # parameters. Two contractions where the forward path does two matrix
  # products, and neither of them is n_theta wide.
  #
  # attach.input passes states through untouched, so their cotangent goes
  # straight back onto the prediction.
  X2Yvjp <- function(out, pars, fixed = NULL, cotangent) {
    w <- .asCtOut(cotangent)
    if (is.null(gEval$vjp))
      stop("Y(): the reverse mode needs a vector-Jacobian product; rebuild ",
           "with Y(..., derivMode = c(\"forward\", \"reverse\"), ",
           "compile = TRUE).", call. = FALSE)
    params <- c(unclass(pars), unclass(fixed))
    fixedObsParams <- .intersectU(c(attr(pars, "fixed"), names(fixed)), obsParams)

    K <- .ctK(w)
    wm <- .ctSlice(w)
    W <- matrix(0, nrow(out), length(observables),
                dimnames = list(NULL, observables))
    iw <- match(observables, colnames(wm))
    if (any(!is.na(iw))) W[, !is.na(iw)] <- wm[, iw[!is.na(iw)], drop = FALSE]

    if (K == 1L) {
      r <- gEval$vjp(out[, obsStates, drop = FALSE], params[obsParams], W)
    } else {
      # The directions of the value pass, on the two inputs this node
      # reads: the prediction's own tangents and the parameters it passes
      # through. The cotangent brings its own beside them.
      nd <- K - 1L
      dX <- attr(out, "deriv")
      VX <- array(0, c(nrow(out), length(obsStates), nd))
      if (!is.null(dX)) {
        take <- intersect(obsStates, dimnames(dX)[[2L]])
        if (length(take))
          VX[, match(take, obsStates), ] <- dX[, take, seq_len(nd), drop = FALSE]
      }
      VP <- matrix(0, length(obsParams), nd, dimnames = list(obsParams, NULL))
      dP <- attr(pars, "deriv")
      if (!is.null(dP)) {
        take <- intersect(rownames(dP), obsParams)
        if (length(take)) VP[take, ] <- dP[take, seq_len(nd), drop = FALSE]
      }
      DW <- array(0, c(nrow(out), length(observables), 1L, nd))
      hitw <- intersect(dimnames(w)[[2L]], observables)
      if (length(hitw))
        DW[, match(hitw, observables), 1L, ] <- w[, hitw, -1L, drop = FALSE]
      r <- gEval$vjp(out[, obsStates, drop = FALSE], params[obsParams], W,
                     tangentX = VX, tangentP = VP, curvature = DW)
    }

    w_out <- array(0, c(nrow(out), ncol(out), K),
                   dimnames = list(NULL, colnames(out), NULL))
    hit <- .intersectU(obsStates, colnames(out))
    if (length(hit)) {
      w_out[, hit, 1L] <- r$cotangentX[, hit, 1L]
      if (K > 1L)
        w_out[, hit, -1L] <- r$curvatureX[, match(hit, obsStates), 1L, , drop = FALSE]
    }
    # Everything attach.input passed through keeps whatever the caller put on
    # it, the observables aside.
    if (controls$attach.input) {
      through <- intersect(setdiff(colnames(wm), c("time", observables)), colnames(out))
      if (length(through))
        w_out[, through, ] <- w_out[, through, , drop = FALSE] +
                              w[, through, , drop = FALSE]
    }

    up <- matrix(r$cotangentP[, 1L], ncol = 1L,
                 dimnames = list(rownames(r$cotangentP), NULL))
    if (K > 1L) {
      up <- cbind(up, matrix(r$curvatureP[, 1L, ], nrow = nrow(up)))
      rownames(up) <- rownames(r$cotangentP)
    }
    w_pars <- .pickCotangent(up, names(pars))
    w_pars[match(names(pars), fixedObsParams, 0L) > 0L, ] <- 0
    .ct(out = w_out, pars = w_pars)
  }

  attr(X2Y, "vjpfn") <- X2Yvjp
  attr(X2Y, "batchfn") <- X2Ybatch

  attr(X2Y, "equations")  <- as.eqnvec(g)
  attr(X2Y, "parameters") <- parameters
  attr(X2Y, "states")     <- states
  attr(X2Y, "modelname")  <- modelname
  attr(X2Y, "compileInfo") <- .collectCompileInfo(gfun, gjac, ghess, gevaluate)

  obsfn(X2Y, parameters, condition)
}

 
#' Generate a prediction function that returns times
#'
#' Function to deal with non-ODE models within the framework of dMod. See example.
#'
#' @param condition  either NULL (generic prediction for any condition) or a character, denoting
#' the condition for which the function makes a prediction.
#' @return Object of class [prdfn].
#' @examples
#' \donttest{
#' x <- Xt()
#' g <- Y(c(y = "a*time^2+b"), f = NULL, parameters = c("a", "b"),
#'        compile = TRUE, modelname = "Xt_example_obs", outdir = tempdir())
#'
#' times <- seq(-1, 1, by = .05)
#' pars <- c(a = .1, b = 1)
#'
#' plot((g*x)(times, pars))
#' }
#' @export
Xt <- function(condition = NULL) {
  P2X <- function(times, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE) {
    n_times <- length(times)
    par_names <- names(pars)
    n_pars <- length(par_names)

    out <- matrix(times, ncol = 1, dimnames = list(NULL, "time"))

    # time has no parameter dependence, both derivative arrays are zero
    # in batch-first [time, observable, ...] layout matching Xs.
    sens  <- array(0, dim = c(n_times, 1, n_pars),
                   dimnames = list(NULL, "time", par_names))
    sens2 <- if (deriv2)
      array(0, dim = c(n_times, 1, n_pars, n_pars),
            dimnames = list(NULL, "time", par_names, par_names))
    else NULL

    # Fixed parameters have no sensitivity but must stay visible downstream
    # (an error model reads them off the prediction), as in Xs.
    prdframe(out, deriv = sens, deriv2 = sens2, parameters = c(pars, fixed))
  }
  # Time depends on nothing, so its cotangent is nothing. The pass-through of
  # the parameters is the caller's business and happens above this leaf.
  attr(P2X, "vjpfn") <- function(times, pars, fixed = NULL, cotangent)
    .ctZero(names(pars), .ctK(.asCtOut(cotangent)))

  attr(P2X, "parameters") <- NULL
  attr(P2X, "equations") <- NULL
  attr(P2X, "forcings") <- NULL
  attr(P2X, "events") <- NULL
  prdfn(P2X, NULL, condition)
}


#' An identity function which vanishes upon concatenation of fns
#'
#' @return fn of class idfn
#' @export
#'
#' @examples
#' x <- Xt()
#' id <- Id()
#'
#' (id*x)(1:10, pars = c(a = 1))
#' (x*id)(1:10, pars = c(a = 1))
#' str(id*x)
#' str(x*id)
Id <- function() {
  outfn <- function() return(NULL)
  class(outfn) <- c("idfn", "fn")
  return(outfn)
}
