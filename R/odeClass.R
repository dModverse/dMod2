## Methods of class odemodel

#' @export
print.odemodel <- function(x, ...) {

  func      <- x$func
  extended  <- x$extended
  extended2 <- x$extended2
  reversed  <- x$reversed
  reversed2 <- x$reversed2

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
    if (!is.null(reversed))
      cat("  Adjoint: ", as.character(reversed), "\n", sep = "")
    if (!is.null(reversed2))
      cat("  Adjoint2:", as.character(reversed2), "\n", sep = "")
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


#' Compiled ODE model
#'
#' Generates and compiles the code of an ODE system and, as asked for, of its
#' first- and second-order sensitivities and its adjoint. The result is turned
#' into a prediction function by [Xs()] or [Xf()].
#'
#' @param f Named character vector, [eqnvec] or [eqnlist] with the right-hand
#'   sides; the names are the states.
#' @param deriv Logical. Build first-order sensitivities. `FALSE` builds the
#'   states alone, for [Xf()], and leaves `derivMode` without effect.
#' @param deriv2 Logical. Build second-order sensitivities, the same as
#'   `"forward-forward"` in `derivMode`. Needs `backend = "cppDE"` and implies
#'   `deriv = TRUE`.
#' @param derivMode Which derivative directions to build. More than one may be
#'   named; the default is `"forward"` alone.
#'   * `"forward"`: sensitivity equations integrated alongside the states.
#'   * `"reverse"`: the adjoint, whose cost does not grow with the number of
#'     parameters; what `obj(..., sweep = "reverse")` needs. On
#'     `backend = "cppDE"` it is the discrete adjoint of the integrator, on
#'     `backend = "Sundials"` CVODES adjoint sensitivity analysis, which
#'     refuses events. Not available on `backend = "deSolve"`.
#'   * `"forward-forward"`: second-order sensitivities as nested duals.
#'   * `"forward-reverse"`: the adjoint over tangents, which returns the
#'     Hessian of the objective in one backward sweep; what
#'     `obj(..., sweep = "reverse", deriv2 = TRUE)` needs. `cppDE` only.
#'
#'   Each direction is a compilation of its own.
#' @param forcings Character vector, the names of the forcings in `f`. Their
#'   data are given to [Xs()] or [Xf()].
#' @param events An [eventlist], or a `data.frame` [as.eventlist()] accepts.
#'   Defined here, they also act on the sensitivities and, on `cppDE`, the
#'   adjoint.
#' @param fixed Character vector, the initial values and parameters without
#'   sensitivities.
#' @param modelname Character, the base name of the generated files and
#'   symbols.
#' @param backend `"cppDE"` (default), `"Sundials"` or `"deSolve"`: the
#'   integrators of \pkg{cppDE}, CVODE(S) through \pkg{cppDE}, or \pkg{deSolve}
#'   through \pkg{cOde}.
#' @param verbose Logical. Print the compiler output.
#' @param outdir Directory for the generated sources and shared objects.
#'   `backend = "deSolve"` writes to the working directory only.
#' @param ... Passed to [cppDE::cppODE()], [cppDE::cvode()] or
#'   [cOde::funC()], according to `backend`. Among them:
#'   * `method`: the integration method, `"bdf"` (default), `"adams"`,
#'     `"rb4"` or `"tsit5"` on `cppDE`; `"bdf"` or `"adams"` on `Sundials`.
#'     On `deSolve` it is `optionsOde$method` of [Xs()] instead.
#'   * `rootfunc`: `"equilibrate"` or expressions whose roots end the
#'     integration (`cppDE`, `Sundials`).
#'   * `includeTimeZero`: integrate from 0 rather than from the first
#'     requested time, default `TRUE`.
#'   * `sparse`: sparse or dense linear algebra, `NULL` chooses from the
#'     Jacobian (`cppDE`, `Sundials`).
#'   * `stepTrace`: record per-step diagnostics, see `traceFile` in [Xs()]
#'     (`cppDE`, `Sundials`).
#'   * `compile`: `FALSE` writes the sources only, [compile()] builds them
#'     later.
#'   * `estimate`, `outputs`, `gridpoints`: the parameters with
#'     sensitivities, additional outputs and spline nodes of the forcings
#'     (`deSolve`).
#'
#' @return An object of class `odemodel`: a list with the compiled objects
#'   `func` (states), `extended` (first-order sensitivities), and on
#'   `cppDE` / `Sundials` `extended2` (second order) and `reversed`,
#'   `reversed2` (adjoints), each `NULL` where not built. The attribute
#'   `"compileInfo"` holds the sources and compiler flags [compile()] uses.
#'
#' @seealso [Xs()], [Xf()], [cppDE::cppODE()], [cppDE::cvode()],
#'   [cOde::funC()]
#'
#' @example inst/examples/odemodel.R
#' @export
odemodel <- function(f, deriv = TRUE, deriv2 = FALSE, derivMode = "forward",
                     forcings=NULL, events = NULL,
                     fixed = NULL, modelname = "odemodel", backend = c("cppDE", "Sundials", "deSolve"),
                     verbose = FALSE, outdir = getwd(), ...) {

  f <- as.eqnvec(f)
  backend <- match.arg(backend)
  derivMode <- .matchDerivMode(derivMode, c("forward", "reverse",
                                            "forward-forward", "forward-reverse"))
  reverse  <- "reverse" %in% derivMode
  reverse2 <- "forward-reverse" %in% derivMode
  # The two spellings of forward over forward name the same build product.
  if ("forward-forward" %in% derivMode) deriv2 <- TRUE
  if (reverse && backend == "deSolve")
    stop("derivMode = \"reverse\" needs backend = 'cppDE' or 'Sundials'.",
         call. = FALSE)
  if (reverse2 && backend != "cppDE")
    stop("derivMode = \"forward-reverse\" needs backend = 'cppDE'.",
         call. = FALSE)

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
      # The reverse object. A fourth compilation beside func, extended and
      # extended2, not a flag on any of them: the direction decides what the
      # generated code is, and cannot be chosen after the fact.
      reversed <- NULL
      if (reverse) {
        reversed <- do.call(cppDE::cppODE,
                            c(list(f, events = events, fixed = fixed, forcings = forcings,
                                   modelname = paste0(modelname, "_r"), outdir = outdir,
                                   derivMode = "reverse", verbose = verbose),
                              dots_func))
      }
      # Forward over reverse: the same sweep over tangents. It integrates under
      # sensitivities like `extended` and sweeps like `reversed`, so it takes
      # the sensitivity dots rather than the value ones.
      reversed2 <- NULL
      if (reverse2) {
        reversed2 <- do.call(cppDE::cppODE,
                             c(list(f, events = events, fixed = fixed, forcings = forcings,
                                    modelname = paste0(modelname, "_r2"), outdir = outdir,
                                    derivMode = "forward-reverse", verbose = verbose),
                               dots_ext))
      }
      out <- list(func = func, extended = extended, extended2 = extended2,
                  reversed = reversed, reversed2 = reversed2)
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
      # CVODES adjoint sensitivity analysis. A separate compilation, as on the
      # native backend, and with the same interface: solveODE(..., cotangent = W)
      # returns $cotangent. It refuses events, which cppDE::cvode() reports.
      reversed <- NULL
      if (reverse) {
        reversed <- do.call(cppDE::cvode,
                            c(list(f, events = events, fixed = fixed, forcings = forcings,
                                   modelname = paste0(modelname, "_r"), outdir = outdir,
                                   deriv = FALSE, derivMode = "reverse", verbose = verbose),
                              dots_func))
      }
      out <- list(func = func, extended = extended, reversed = reversed)
      class(out) <- c("cppDE", "odemodel")
      }
  }
  attr(out, "compileInfo") <- .collectCompileInfo(out$func, out$extended,
                                                  out$extended2, out$reversed,
                                                  out$reversed2)
  return(out)
}
