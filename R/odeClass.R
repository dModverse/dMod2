## Methods of class odemodel

#' @export
print.odemodel <- function(x, ...) {

  func      <- x$func
  extended  <- x$extended
  extended2 <- x$extended2

  ## attr(func, "backend") is cppDE's own marker, not dMod's `backend` argument.
  isCVODE <- inherits(x, "cppDE") && identical(attr(func, "backend"), "cvode")
  isCppDE <- inherits(x, "cppDE") && !isCVODE
  backend <- if (inherits(x, "deSolve")) "deSolve (cOde)"
             else if (isCVODE)           "Sundials (CVODE)"
             else if (isCppDE)           "cppDE"
             else                        "unknown"
  stepper <- attr(func, "method")

  suppressWarnings({

  cat("dMod odemodel\n", sep = "")
  cat("  Backend: ", backend,
      if (!is.null(stepper)) paste0(" [", stepper, "]") else "", "\n", sep = "")
  cat("  Model:   ", as.character(func), "\n", sep = "")
  if (!is.null(extended)) {
    cat("  Sens1:   ", as.character(extended), "\n", sep = "")
  } else {
    cat("  Sens1:   not compiled (deriv = FALSE)\n", sep = "")
  }
  if (isCppDE) {
    if (!is.null(extended2)) {
      cat("  Sens2:   ", as.character(extended2), "\n", sep = "")
    } else {
      cat("  Sens2:   not compiled (deriv2 = FALSE)\n", sep = "")
    }
  }

  cat("\nEquations:\n", sep = "")
  print(as.eqnvec(attr(func, "equations")))
  cat("\nStates:\n", sep = "")
  print(sort(attr(func, "variables")))
  cat("\nParameters:\n", sep = "")
  print(sort(attr(func, "parameters")))
  forcs <- attr(func, "forcings")
  if (length(forcs) > 0) {
    cat("\nForcings:\n", sep = "")
    print(sort(forcs))
  }

  })

  invisible(x)
}


#' Generate model objects for use in Xs (models with sensitivities)
#'
#' Creates and compiles model objects for systems of ordinary differential equations (ODEs)
#' with optional first- and second-order sensitivities. Depending on the selected backend,
#' the function interfaces either to [cOde::funC()] (for `backend = "deSolve"`)
#' or to [cppDE::cppODE()] / [cppDE::cvode()] (for `backend = "cppDE"` or
#' `backend = "Sundials"`).
#'
#' @param f Something that can be converted to [eqnvec], e.g. a named character vector
#'   specifying the right-hand sides of the ODE system.
#' @param deriv Logical. If `TRUE`, generate first-order sensitivities.
#'   Defaults to `TRUE`.
#' @param deriv2 Logical. If `TRUE`, also generate second-order sensitivities
#'   (requires `backend = "cppDE"`). Implies `deriv = TRUE`. Defaults to
#'   `FALSE`.
#' @param forcings Character vector with the names of external forcings.
#' @param events An [eventlist] (or `data.frame` coercible via [as.eventlist]).
#'   Must be defined here -- not on [Xs()] -- so that the sensitivity equations
#'   are extended consistently.
#' @param fixed Character vector with the names of parameters (initial values and dynamic)
#'   for which no sensitivities are required (this speeds up integration).
#' @param modelname Character. The base name of the generated C/C++ file.
#' @param backend Character string selecting the code-generation and integration
#'   backend. One of `"cppDE"`, `"Sundials"` or `"deSolve"`. The stepper itself
#'   is chosen separately (`method` for `"cppDE"`, `optionsOde$method` for
#'   `"deSolve"`).
#' @param verbose Logical. If `TRUE`, print compiler output to the R console.
#' @param outdir Character. Directory for the generated C/C++ sources and the
#'   compiled shared object. Defaults to the working directory. Only honoured
#'   for `backend = "cppDE"` / `"Sundials"`; the `deSolve` backend always
#'   writes to the working directory, so a non-default `outdir` errors there.
#' @param ... Additional arguments passed to [cppDE::cppODE()], [cppDE::cvode()]
#'   or [cOde::funC()].
#'
#' @return list with \code{func} (ODE object) and \code{extended} (ODE+Sensitivities object).
#'   Carries a \code{"compileInfo"} attribute listing source files and per-file
#'   compile/link flags collected from \code{func} and \code{extended}. This is
#'   consumed by [compile()] when the model is later compiled via a prediction
#'   function, so backend-specific linker requirements (e.g. Sundials libraries
#'   for \code{backend = "Sundials"}) are applied to the right files only.
#'
#' @seealso [cOde::funC()], [cppDE::cppODE()], [cppDE::cvode()]
#'
#' @example inst/examples/odemodel.R
#' @export
odemodel <- function(f, deriv = TRUE, deriv2 = FALSE, forcings=NULL, events = NULL,
                     fixed = NULL, modelname = "odemodel", backend = c("cppDE", "Sundials", "deSolve"),
                     verbose = FALSE, outdir = getwd(), ...) {

  f <- as.eqnvec(f)
  backend <- match.arg(backend)

  if (deriv2 && !deriv) {
    warning("`deriv2 = TRUE` implies `deriv = TRUE`. Setting deriv = TRUE.",
            call. = FALSE)
    deriv <- TRUE
  }

  dots <- list(...)

  if (deriv2 && backend == "deSolve")
    stop("Second-order sensitivities require backend = 'cppDE'.")
  if (deriv2 && backend == "Sundials")
    stop("Second-order sensitivities are not available with CVODE; use backend = 'cppDE'.")

  # cOde::funC has no outdir and always writes to the working directory, so a
  # non-default outdir would be silently ignored under the deSolve backend.
  if (backend == "deSolve" &&
      !identical(normalizePath(outdir,   mustWork = FALSE),
                 normalizePath(getwd(), mustWork = FALSE)))
    stop("`outdir` is only supported for backend = 'cppDE'/'Sundials'; ",
         "the deSolve backend always writes to the working directory. ",
         "setwd() to the target directory instead.", call. = FALSE)

  pick <- function(fn, args) {
    fm <- names(formals(fn))
    if ("..." %in% fm) args else args[intersect(names(args), fm)]
  }

  if (backend == "deSolve") {

    estimate   <- dots$estimate;   dots$estimate   <- NULL
    outputs    <- dots$outputs;    dots$outputs    <- NULL
    gridpoints <- dots$gridpoints; dots$gridpoints <- NULL

    if (is.null(gridpoints)) gridpoints <- 2
    ## `solver` is cOde::funC's own argument, not dMod's `backend`.
    func <- do.call(cOde::funC,
                    c(list(f, forcings = forcings, events = events, outputs = outputs,
                           fixed = fixed, modelname = modelname, solver = "deSolve",
                           nGridpoints = gridpoints),
                      pick(cOde::funC, dots)))
    extended <- NULL
    if (deriv) {
      modelname_s <- paste0(modelname, "_s")
      mystates <- attr(func, "variables")
      myparameters <- attr(func, "parameters")

      if (is.null(estimate) & !is.null(fixed)) {
        mystates <- setdiff(mystates, fixed)
        myparameters <- setdiff(myparameters, fixed)
      }

      if (!is.null(estimate)) {
        mystates <- intersect(mystates, estimate)
        myparameters <- intersect(myparameters, estimate)
      }

      s <- sensitivitiesSymb(f,
                             states = mystates,
                             parameters = myparameters,
                             inputs = attr(func, "forcings"),
                             events = attr(func, "events"),
                             reduce = TRUE)
      fs <- c(f, s)
      outputs <- c(attr(s, "outputs"), attr(func, "outputs"))

      events.sens <- attr(s, "events")
      events.func <- attr(func, "events")
      events <- NULL
      if (!is.null(events.func)) {
        if (is.data.frame(events.sens)) {
          events <- rbind(
            as.eventlist(events.sens),
            as.eventlist(events.func),
            stringsAsFactors = FALSE)
        } else {
          events <- do.call(rbind, lapply(1:nrow(events.func), function(i) {
            rbind(
              as.eventlist(events.sens[[i]]),
              as.eventlist(events.func[i,]),
              stringsAsFactors = FALSE)
          }))
        }

      }

      extended <- do.call(cOde::funC,
                          c(list(fs, forcings = forcings, modelname = modelname_s,
                                 solver = "deSolve", nGridpoints = gridpoints,
                                 events = events, outputs = outputs),
                            pick(cOde::funC, dots)))
    }
    out <- list(func = func, extended = extended)
    class(out) <- c("deSolve", "odemodel")
  }
  else {
    if (backend == "cppDE") {
      dots_func <- pick(cppDE::cppODE, dots[setdiff(names(dots), "nStack")])
      dots_ext  <- pick(cppDE::cppODE, dots)
      func <- do.call(cppDE::cppODE,
                      c(list(f, events = events, fixed = fixed, forcings = forcings,
                             modelname = modelname,
                             outdir = outdir, deriv = FALSE, verbose = verbose),
                        dots_func))
      extended <- NULL
      extended2 <- NULL
      if (deriv) {
        extended <- do.call(cppDE::cppODE,
                            c(list(f, events = events, fixed = fixed, forcings = forcings,
                                   modelname = paste0(modelname, "_s"), outdir = outdir,
                                   deriv = TRUE, deriv2 = FALSE, verbose = verbose),
                              dots_ext))
        if (deriv2) {
          extended2 <- do.call(cppDE::cppODE,
                               c(list(f, events = events, fixed = fixed, forcings = forcings,
                                      modelname = paste0(modelname, "_s2"), outdir = outdir,
                                      deriv = TRUE, deriv2 = TRUE, verbose = verbose),
                                 dots_ext))
        }
      }
      out <- list(func = func, extended = extended, extended2 = extended2)
      class(out) <- c("cppDE", "odemodel")
    } else if (backend == "Sundials") {
      dots_func <- pick(cppDE::cvode, dots[setdiff(names(dots), "nStack")])
      dots_ext  <- pick(cppDE::cvode, dots)
      func <- do.call(cppDE::cvode,
                      c(list(f, events = events, fixed = fixed, forcings = forcings,
                             modelname = modelname,
                             outdir = outdir, deriv = FALSE, verbose = verbose),
                        dots_func))
      extended <- NULL
      if (deriv) {
        extended <- do.call(cppDE::cvode,
                            c(list(f, events = events, fixed = fixed, forcings = forcings,
                                   modelname = paste0(modelname, "_s"), outdir = outdir,
                                   deriv = TRUE, verbose = verbose),
                              dots_ext))
      }
      out <- list(func = func, extended = extended)
      class(out) <- c("cppDE", "odemodel")
      }
  }
  attr(out, "compileInfo") <- .collectCompileInfo(out$func, out$extended, out$extended2)
  return(out)
}
