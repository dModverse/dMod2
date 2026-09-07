#' Non-linear optimisation via a trust-region method
#'
#' \code{trust} minimises (or maximises) a smooth objective function for which
#' value, gradient and Hessian are available. It solves a Moré-Sorensen
#' trust-region subproblem exactly, via the full eigendecomposition of the
#' Hessian.
#'
#' @section Box bounds:
#' \code{boundary = "reflective"} (default) uses the Coleman-Li interior
#' trust-region-reflective scheme: the subproblem is solved in the scaled frame
#' \code{D g}, \code{D H D + C} with \code{D = diag(|v|^(1/2))}, \code{|v_i|}
#' being the distance to the bound the gradient pushes toward. A step leaving
#' the box is truncated, reflected off the blocking faces, or replaced by the
#' scaled steepest-descent step, whichever the model rates best. Iterates stay
#' strictly inside the box and never land exactly on a bound, use
#' \code{atBound} to test bound activity.
#'
#' \code{boundary = "clip"} is the historical scheme (active-set reduction plus
#' componentwise clipping), kept for reproducing earlier fits.
#'
#' @section Convergence:
#' The run stops as soon as any of the following holds. The gradient tests are
#' evaluated at the current iterate; the value, model and step tests count only
#' on an accepted step.
#' \itemize{
#'   \item \code{gtol}: \code{max(|v * g|) <= gtol}. This is the first-order
#'     optimality measure of the box problem and vanishes both at an interior
#'     stationary point and at a bound the gradient pushes into. Without bounds
#'     it is \code{max(|g|)}.
#'   \item \code{ftol}: \code{|f - f_try| < ftol}.
#'   \item \code{mtol}: predicted reduction
#'     \code{|g^T p + 0.5 * p^T H p| < mtol}, the quantity `blather` reports as
#'     \code{preddiff}.
#'   \item \code{xtol}: step norm below \code{xtol} (disabled when \code{0}).
#'   \item stagnation: five consecutive rejected steps that left the objective
#'     flat within \code{ftol}.
#'   \item \code{rmin}: trust radius below \code{rmin}. This one reports
#'     \code{converged = FALSE} (disabled when \code{0}).
#' }
#' \code{stopReason} in the result names which test fired: \code{"gradient"},
#' \code{"fvalue"}, \code{"preddiff"}, \code{"step"}, \code{"stagnation"},
#' \code{"radius"}, \code{"objfun"} or \code{"iterlim"}. Stagnation reports
#' \code{converged = TRUE}; use \code{rmin} for a hard failure instead.
#'
#' Only \code{gtol} tests first-order optimality; the other four read the last
#' accepted step and a misleading Hessian can trip them far from a stationary
#' point. \code{hessianFallback} uses that: the first of the four hands the
#' Hessian source over instead of ending the run. In a quasi-Newton
#' phase \code{"fvalue"}, \code{"preddiff"} and \code{"step"} do not end the
#' run at all; it ends on \code{"gradient"}, \code{"stagnation"} or the radius.
#'
#' Whatever ends the run, \code{argument}, \code{value}, \code{gradient},
#' \code{hessian} and \code{atBound} describe one and the same iterate: the
#' best one visited, which under the default \code{stepControl$nonmonotone = 0}
#' is the last accepted one.
#'
#' @section Choosing tolerances:
#' The optimiser can only resolve what the objective delivers. For an ODE model
#' integrated at relative tolerance \code{rtol}, the value carries a relative
#' error of roughly \code{rtol} and the forward sensitivities about an order
#' more, giving the gradient a noise floor. Below it, further iterations chase
#' integration error. As a starting point:
#' \itemize{
#'   \item \code{ftol} at or above \code{rtol * |f|}.
#'   \item \code{gtol} at the gradient's noise floor, which for a poorly
#'     scaled problem is easier to read off a first run than to predict.
#'   \item Tighten the integrator before tightening either.
#' }
#'
#' A run ending on \code{"stagnation"} rather than \code{"gradient"} indicates
#' the gradient tolerances sit below the noise floor; \code{max(abs(gradient))}
#' of the result shows where it actually lies.
#'
#' @param objfun R function whose first argument is a numeric vector of
#'   parameters. Must return a list with components \code{value},
#'   \code{gradient}, \code{hessian}. Extra arguments accepted by
#'   \code{objfun} can be supplied via \code{...}.
#' @param parinit Named numeric starting vector. Must be finite. Values
#'   outside \code{[parlower, parupper]} are clipped with a warning; with
#'   \code{boundary = "reflective"} the result is additionally nudged just
#'   inside the box.
#' @param rinit Initial trust-region radius.
#' @param rmax Maximum allowed trust-region radius.
#' @param parscale Optional named or unnamed numeric of length
#'   \code{length(parinit)} for parameter rescaling. The subproblem
#'   operates on \code{g / parscale} and
#'   \code{H / outer(parscale, parscale)}; the trust radius and \code{xtol}
#'   are measured in that frame.
#' @param iterlim Maximum number of outer trust-region iterations.
#' @param hessianMethod Source of the model Hessian the run starts on:
#'   \code{"gn"} (default) uses the objective's Gauss-Newton Hessian,
#'   \code{"bfgs"} and \code{"sr1"} maintain a quasi-Newton update seeded from
#'   it. A run is a method and, optionally, the \code{hessianFallback} it hands
#'   over to. Quasi-Newton methods require
#'   \code{stepControl$boundary = "reflective"}. Also selects the control
#'   defaults below.
#' @param hessianFallback Source the run hands over to at the first value,
#'   model, step or stagnation stop, instead of ending there. \code{"none"}
#'   (default) never hands over. A handover restarts the trust region at
#'   \code{rinit} and discards the quasi-Newton pairs, keeping the working
#'   Hessian. Handing back to \code{"gn"} costs one extra objective
#'   evaluation.
#' @param fallbackLimit Number of handovers allowed. \code{1} (default) hands
#'   over once and runs the fallback to the end; a higher value alternates back
#'   and forth, each time at the next soft stop. \code{0} disables the fallback.
#'   The result reports the count as \code{nSwitch}.
#' @param tolControl List of convergence thresholds, merged into the defaults.
#'   \describe{
#'     \item{\code{ftol}}{Change in objective value. Default \code{1e-6}.}
#'     \item{\code{mtol}}{Predicted model decrease. Default \code{1e-6}.}
#'     \item{\code{gtol}}{First-order optimality measure. Default \code{1e-6}.
#'       Deliberately absolute: scaling a gradient by the objective value is not
#'       meaningful for a general objective, and for a sum of squares it is
#'       wrong, since \code{|f|} grows quadratically in the residuals while
#'       \code{|g|} grows linearly. Pass \code{gtol = x * |f|} for that.}
#'     \item{\code{xtol}}{Step norm. Default \code{0}, which disables it.}
#'     \item{\code{rmin}}{Lower limit on the trust radius; falling below it
#'       stops the run with \code{converged = FALSE}. Default \code{0},
#'       which disables it.}
#'   }
#' @param qnControl List of quasi-Newton settings, merged into the defaults and
#'   read only by \code{"bfgs"}, \code{"sr1"} and a quasi-Newton fallback
#'   phase.
#'   \describe{
#'     \item{\code{hessianInit}}{Seed of the approximation: \code{"gn"}
#'       (default) the objective's Hessian at \code{parinit}, \code{"identity"}
#'       the identity, which also stops \code{trust} asking for a Hessian at
#'       all. Read only when the run starts quasi-Newton, so it is inert for
#'       \code{"gn"}.}
#'     \item{\code{qnMemory}}{Number of \code{(s, y)} pairs kept. \code{0}
#'       (default) accumulates every update onto the seed. A positive value
#'       rebuilds the approximation each iteration from the last
#'       \code{qnMemory} pairs and the Shanno-Phua scaling of the newest with
#'       positive curvature, discarding the seed. Costs
#'       \code{O(qnMemory * K^2)} per iteration and saves no memory, since the
#'       subproblem needs the explicit matrix.}
#'     \item{\code{qnCautious}}{Cautious-update threshold: a BFGS pair counts
#'       only when \code{s^T y > qnCautious * |s| * |y|}. Default \code{1e-8},
#'       \code{0} disables the test, \code{"sr1"} is exempt. \code{qnSkipped}
#'       counts the rejected pairs.}
#'     \item{\code{qnRejected}}{Whether \code{"sr1"} also takes its
#'       \code{(s, y)} pair from a rejected trial point. Default \code{TRUE}
#'       for \code{"sr1"}, \code{FALSE} otherwise; \code{"gn"} and
#'       \code{"bfgs"} are unaffected either way.}
#'   }
#' @param stepControl List of step and acceptance settings, merged into the
#'   defaults.
#'   \describe{
#'     \item{\code{boundary}}{Box-bound handling, \code{"reflective"} (default)
#'       or \code{"clip"}. See Details.}
#'     \item{\code{theta.max}}{Largest fraction of the distance to a bound that
#'       a step may use, keeping iterates strictly interior. Default
#'       \code{0.99995}, used only by \code{boundary = "reflective"}.}
#'     \item{\code{nonmonotone}}{Zhang-Hager relaxation \code{eta} in
#'       \code{[0, 1)}. The step ratio is measured against a weighted average of
#'       past objective values instead of the current one, so a step may be
#'       accepted that increases the objective. Default \code{0}, the monotone
#'       rule; \code{0.85} is the value Zhang and Hager report as generally
#'       best. The result reports the best iterate visited, which is the last
#'       one when \code{eta = 0}.}
#'   }
#' @param minimize If \code{TRUE} (default) minimise; if \code{FALSE}
#'   maximise.
#' @param blather If \code{TRUE} return the per-iteration trace
#'   (\code{argpath}, \code{argtry}, \code{steptype}, \code{stepback},
#'   \code{accept}, \code{r}, \code{rho}, \code{valpath}, \code{valtry},
#'   \code{preddiff}, \code{stepnorm}, \code{hessianSource}).
#' @param parupper,parlower Named or scalar numeric bounds. If unnamed,
#'   the first element broadcasts to all parameters; if named, the
#'   entries slot by name into a length-K vector defaulting to
#'   \code{+/- Inf}.
#' @param printIter If \code{TRUE} print iteration count and objective
#'   value to the console at each function evaluation.
#' @param traceFile Optional path. If non-\code{NULL}, CSV-log per
#'   evaluation \code{iter, value, p1, p2, ...}.
#' @param ... Additional named arguments forwarded to \code{objfun}.
#'
#' @return A list with components \code{argument}, \code{value},
#'   \code{gradient}, \code{hessian}, \code{iterations}, \code{neval},
#'   \code{qnEval} (evaluations spent in the quasi-Newton phase),
#'   \code{qnSkipped} (pairs the cautious test rejected),
#'   \code{nSwitch} (Hessian source handovers),
#'   \code{evalBySource} (named integer, evaluations per source),
#'   \code{converged}, \code{atBound} (named logical, which parameters are
#'   held by a bound) and \code{stopReason}. When \code{blather = TRUE} the
#'   list also contains \code{argpath}, \code{argtry}, \code{steptype},
#'   \code{stepback}, \code{accept}, \code{r}, \code{rho}, \code{valpath},
#'   \code{valtry}, \code{preddiff}, \code{stepnorm}.
#'
#' @references Coleman, T. F. and Li, Y. (1996). An interior trust region
#'   approach for nonlinear minimization subject to bounds.
#'   \emph{SIAM Journal on Optimization} 6(2), 418-445.
#'
#'   Fröhlich, F. and Sorger, P. K. (2022). Fides: Reliable trust-region
#'   optimization for parameter estimation of ODE models.
#'   \emph{PLoS Computational Biology} 18(7), e1010322.
#'
#' @export
trust <- function(objfun, parinit, rinit = 0.1, rmax = 10,
                  iterlim   = 100L,
                  hessianMethod   = c("gn", "bfgs", "sr1"),
                  hessianFallback = c("none", "bfgs", "sr1", "gn"),
                  fallbackLimit   = 1L,
                  parscale  = NULL,
                  parupper  = NULL,
                  parlower  = NULL,
                  tolControl  = NULL,
                  qnControl   = NULL,
                  stepControl = NULL,
                  minimize  = TRUE,
                  blather   = FALSE,
                  printIter = FALSE,
                  traceFile = NULL,
                  ...) {
  hessianMethod   <- match.arg(hessianMethod)
  hessianFallback <- match.arg(hessianFallback)
  dots <- list(...)
  .trustRejectMoved(names(dots), "trust")

  ctl  <- .trustDefaults(hessianMethod, hessianFallback)
  tol  <- .mergeControl(ctl$tolControl,  tolControl,  "tolControl")
  qn   <- .mergeControl(ctl$qnControl,   qnControl,   "qnControl")
  step <- .mergeControl(ctl$stepControl, stepControl, "stepControl")

  boundary    <- match.arg(step$boundary, c("reflective", "clip"))
  hessianInit <- match.arg(qn$hessianInit, c("gn", "identity"))

  # The kernel passes hessian = FALSE in the quasi-Newton phase, so the wrapper
  # forwards it; dMod objectives skip J^T J then, others ignore it via `...`.
  fn <- function(x, hessian = TRUE)
    do.call(objfun, c(list(x, hessian = hessian), dots))
  trust_impl(fn, parinit, rinit, rmax, parscale, as.integer(iterlim),
             tol$ftol, tol$mtol, tol$gtol, tol$xtol, tol$rmin, step$theta.max,
             boundary, hessianMethod, hessianFallback,
             as.integer(fallbackLimit), hessianInit,
             as.integer(qn$qnMemory), qn$qnCautious, isTRUE(qn$qnRejected),
             step$nonmonotone, minimize, blather,
             parupper, parlower, printIter, traceFile)
}


# Defaults per Hessian source, and the single source of truth for what each
# control group accepts. Members of one group are settable together or not at
# all; see the roxygen of `trust` for what they mean.
.trustDefaults <- function(hessianMethod = "gn", hessianFallback = "none") {
  qn <- list(hessianInit = "gn", qnMemory = 0L, qnCautious = 1e-8,
             qnRejected = FALSE)
  # Nocedal and Wright, sec. 6.2: SR1 also updates from a rejected step. The
  # flag is inert outside an sr1 phase, so a fallback to sr1 arms it too.
  if ("sr1" %in% c(hessianMethod, hessianFallback)) qn$qnRejected <- TRUE
  list(
    tolControl  = list(ftol = 1e-6, mtol = 1e-6, gtol = 1e-6,
                       xtol = 0, rmin = 0),
    qnControl   = qn,
    stepControl = list(boundary = "reflective", theta.max = 0.99995,
                       nonmonotone = 0))
}


# Arguments that used to be flat, and the control member they became. Reported
# by name so a stale call says where its argument went instead of routing it to
# objfun.
.trustMoved <- c(
  ftol = "tolControl$ftol",   mtol  = "tolControl$mtol",
  gtol = "tolControl$gtol",   xtol  = "tolControl$xtol",
  rmin = "tolControl$rmin",   fterm = "tolControl$ftol",
  mterm = "tolControl$mtol",
  hessianInit = "qnControl$hessianInit", qnMemory   = "qnControl$qnMemory",
  qnCautious  = "qnControl$qnCautious",  qnRejected = "qnControl$qnRejected",
  boundary    = "stepControl$boundary",  theta.max  = "stepControl$theta.max",
  nonmonotone = "stepControl$nonmonotone")

.trustRejectMoved <- function(nms, who) {
  hit <- intersect(nms, names(.trustMoved))
  if (length(hit))
    stop(who, ": ", paste(hit, collapse = ", "), " moved into a control list. ",
         "Use ", paste(.trustMoved[hit], collapse = ", "), ".", call. = FALSE)
  invisible(NULL)
}


# Merge a user control list into the method's defaults. Names are checked
# against the defaults, so a typo errors here instead of being ignored.
.mergeControl <- function(defaults, control = NULL, label = "control") {
  if (length(control) == 0L) return(defaults)
  nms <- names(control)
  if (is.null(nms) || !all(nzchar(nms)))
    stop(label, ": all entries must be named.", call. = FALSE)
  unknown <- setdiff(nms, names(defaults))
  if (length(unknown))
    stop(label, ": unknown entries ", paste(unknown, collapse = ", "),
         ". Settable are ", paste(names(defaults), collapse = ", "), ".",
         call. = FALSE)
  modifyList(defaults, as.list(control))
}


#' @export
#' @rdname trust
#' @param ftol,mtol,gtol,xtol,rmin,theta.max,boundary \code{trustL1} only;
#'   its tolerances stay flat, it does not take the control lists of
#'   \code{trust}.
#' @param fterm,mterm Deprecated aliases for \code{ftol} and \code{mtol}.
#' @param mu Named numeric vector of reference values for the L1-penalised
#'   parameters. Names must be a subset of \code{names(parinit)}; only the
#'   named parameters receive a penalty. Each must lie strictly inside
#'   \code{[parlower, parupper]}. Defaults to a zero vector covering all of
#'   \code{parinit}.
#' @param one.sided Logical. If \code{TRUE}, the penalty is one-sided and
#'   acts as a lower wall at \code{mu}: \code{lambda * max(0, mu - p)}.
#'   Otherwise it is the symmetric \code{lambda * |p - mu|}.
#' @param lambda Strength of the L1 penalty. Either a scalar (broadcast to
#'   all entries of \code{mu}) or a named numeric aligned with \code{mu}.
trustL1 <- function(objfun, parinit, mu = 0 * parinit, one.sided = FALSE, lambda = 1,
                    rinit = 0.1, rmax = 10,
                    parscale  = NULL,
                    iterlim   = 100L,
                    ftol      = 1e-6,
                    mtol      = 1e-6,
                    gtol      = 1e-6,
                      xtol      = 0,
                    rmin      = 0,
                    theta.max = 0.99995,
                    boundary  = c("reflective", "clip"),
                    minimize  = TRUE,
                    blather   = FALSE,
                    parupper  = NULL,
                    parlower  = NULL,
                    printIter = FALSE,
                    traceFile = NULL,
                    fterm, mterm,
                    ...) {
  if (!missing(fterm)) ftol <- fterm
  if (!missing(mterm)) mtol  <- mterm
  boundary <- match.arg(boundary)

  sanePars <- sanitizePars(parinit, list(...)$fixed)
  parinit  <- sanePars$pars

  if (is.null(names(parinit)))
    stop("trustL1: parinit must be a named numeric vector")
  if (length(mu) > 0L && is.null(names(mu)))
    stop("trustL1: mu must be a named numeric vector")

  unknown <- setdiff(names(mu), names(parinit))
  if (length(unknown) > 0L)
    stop("trustL1: mu has names not present in parinit: ",
         paste(unknown, collapse = ", "))

  if (length(lambda) == 1L) {
    lambda <- structure(rep(as.numeric(lambda), length(mu)),
                        names = names(mu))
  } else {
    if (is.null(names(lambda)))
      stop("trustL1: lambda must be scalar or a named numeric vector")
    if (!setequal(names(lambda), names(mu)))
      stop("trustL1: names(lambda) must equal names(mu)")
    lambda <- lambda[names(mu)]
  }

  dots <- list(...)
  fn <- if (length(dots) > 0L) {
    function(x) do.call(objfun, c(list(x), dots))
  } else {
    objfun
  }

  mu <- structure(as.numeric(mu), names = names(mu))
  lambda <- structure(as.numeric(lambda), names = names(lambda))

  trustL1_impl(fn, parinit, mu, lambda,
               as.logical(one.sided)[1L], rinit, rmax,
               parscale, as.integer(iterlim),
               ftol, mtol, gtol, xtol, rmin, theta.max,
               boundary, minimize, blather,
               parupper, parlower, printIter, traceFile)
}


# Merge a user control list into `defaults` for a trust()/trustL1() call. Names
# are checked
# against the optimiser's formals, so every optimiser argument is settable and a
# typo errors here instead of silently reaching objfun via `...`.
.trustControl <- function(defaults, control = NULL, optimizer = trust,
                          label = "control") {
  if (length(control) == 0L) return(defaults)
  nms <- names(control)
  if (is.null(nms) || !all(nzchar(nms)))
    stop(label, ": all entries must be named.", call. = FALSE)
  settable <- setdiff(names(formals(optimizer)), c("objfun", "parinit", "..."))
  unknown <- setdiff(nms, settable)
  if (length(unknown))
    stop(label, ": unknown entries ", paste(unknown, collapse = ", "),
         ". Settable are ", paste(settable, collapse = ", "), ".", call. = FALSE)
  modifyList(defaults, as.list(control))
}
