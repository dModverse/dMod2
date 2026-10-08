#' Non-Linear Optimisation by a Trust-Region Method
#'
#' \code{trust} minimises (or maximises) a smooth objective function for which
#' value, gradient and Hessian are available. It solves a Moré-Sorensen
#' trust-region subproblem exactly, via the full eigendecomposition of the
#' Hessian.
#'
#' @section Box bounds:
#' \code{boundary = "reflective"} (default) uses the interior
#' trust-region-reflective scheme of Coleman and Li. Iterates stay strictly
#' inside the box and never land exactly on a bound; \code{atBound} in the
#' result reports which parameters a bound holds.
#'
#' \code{boundary = "clip"}: active-set reduction with componentwise clipping.
#'
#' @section Convergence:
#' The run stops as soon as one of the following holds. The gradient test is
#' evaluated at the current iterate; the value, model and step tests count only
#' on an accepted step.
#' \itemize{
#'   \item \code{gtol}: the first-order optimality measure
#'     \eqn{\max_i |v_i g_i|} is at most \code{gtol}, with \eqn{v_i} the
#'     distance to the bound the gradient points to, 1 if that bound is
#'     infinite. Without bounds it is \eqn{\max_i |g_i|}.
#'   \item \code{ftol}: \eqn{|f - f_{\mathrm{try}}| <} \code{ftol}.
#'   \item \code{mtol}: the predicted reduction
#'     \eqn{|g^\top p + \frac{1}{2} p^\top H p| <} \code{mtol}.
#'   \item \code{xtol}: the step norm is below \code{xtol}.
#'   \item stagnation: five consecutive rejected steps left the objective
#'     flat within \code{ftol}.
#'   \item \code{rmin}: the trust radius is below \code{rmin}. Reports
#'     \code{converged = FALSE}.
#' }
#' \code{stopReason} in the result names the test that fired:
#' \code{"gradient"}, \code{"fvalue"}, \code{"preddiff"}, \code{"step"},
#' \code{"stagnation"}, \code{"radius"}, \code{"objfun"} (three failed
#' evaluations in a row) or \code{"iterlim"}. Stagnation reports
#' \code{converged = TRUE}.
#'
#' Only \code{gtol} tests first-order optimality. With a
#' \code{hessianFallback}, the first \code{"fvalue"}, \code{"preddiff"},
#' \code{"step"} or \code{"stagnation"} stop hands over to the fallback source
#' instead of ending the run. In a quasi-Newton phase \code{"fvalue"},
#' \code{"preddiff"} and \code{"step"} do not end the run.
#'
#' \code{argument}, \code{value}, \code{gradient}, \code{hessian} and
#' \code{atBound} in the result belong to the best iterate visited, which is
#' the last accepted one when \code{stepControl$nonmonotone = 0}.
#'
#' @section Multiple shooting:
#' An objective built by \code{\link{normL2}(..., multipleShootingControl = )},
#' alone or in a sum, is optimised by multiple shooting with the same controls.
#' Not available then: \code{hessianMethod = "exact"}, a
#' \code{hessianFallback}, \code{qnControl$hessianReseed},
#' \code{qnControl$qnMemory}, \code{stepControl$nonmonotone},
#' \code{boundary = "clip"}, \code{minimize = FALSE} and \code{traceFile}.
#' The run reports \code{converged = TRUE} only when the scaled gaps between
#' the segments are below \code{tolControl$ctol}. The method is described in
#' \code{vignette("Optimisation", package = "dMod2")}.
#'
#' @param objfun Function whose first argument is a numeric vector of
#'   parameters, returning a list with components \code{value},
#'   \code{gradient} and \code{hessian}. It is called with the logical
#'   arguments \code{deriv} and \code{hessian}, and with \code{deriv2 = TRUE}
#'   when an exact Hessian is requested, so it must accept them or
#'   \code{...}.
#' @param parinit Named numeric starting vector. Must be finite. Values
#'   outside \code{[parlower, parupper]} are clipped with a warning; with
#'   \code{boundary = "reflective"} they are moved just inside the box.
#' @param rinit Initial trust-region radius. Default \code{0.1}.
#' @param rmax Maximum trust-region radius. Default \code{10}.
#' @param iterlim Maximum number of trust-region iterations. Default
#'   \code{100}.
#' @param hessianMethod Source of the model Hessian at the start of the run:
#'   \code{"gn"} (default) the Gauss-Newton Hessian of the objective,
#'   \code{"bfgs"} or \code{"sr1"} a quasi-Newton update seeded as
#'   \code{qnControl$hessianInit} says, \code{"exact"} the true Hessian at
#'   every iterate, which needs an objective that takes \code{deriv2}.
#'   Quasi-Newton methods require \code{stepControl$boundary = "reflective"}.
#'   Also selects the control defaults below.
#' @param hessianFallback Hessian source the run hands over to at the first
#'   value, model, step or stagnation stop: \code{"none"} (default),
#'   \code{"bfgs"}, \code{"sr1"} or \code{"gn"}. A handover restarts the trust
#'   region at \code{rinit}, discards the quasi-Newton pairs and keeps the
#'   working Hessian.
#' @param fallbackLimit Number of handovers allowed. Default \code{1}: hand
#'   over once and run the fallback to the end. A higher value alternates
#'   between the two sources, \code{0} disables the fallback. The result
#'   reports the count as \code{nSwitch}.
#' @param parscale Deprecated. Optional numeric of length
#'   \code{length(parinit)}, named or unnamed. The trust region works on the
#'   parameters multiplied by \code{parscale}, so the trust radius and
#'   \code{xtol} are measured there. Default \code{NULL}, no scaling. Fit on
#'   log scale instead.
#' @param parupper,parlower Upper and lower bounds. Unnamed: the first element
#'   applies to all parameters. Named: entries are matched by name, missing
#'   parameters are unbounded. Default \code{NULL}, unbounded.
#' @param tolControl List of convergence thresholds, merged into the defaults.
#'   \describe{
#'     \item{\code{ftol}}{Change in objective value. Default \code{1e-6}. For
#'       an ODE model integrated at relative tolerance \code{rtol}, use a value
#'       at or above \code{rtol * abs(f)}.}
#'     \item{\code{mtol}}{Predicted model decrease. Default \code{1e-6}.}
#'     \item{\code{gtol}}{First-order optimality measure, absolute. Default
#'       \code{1e-6}. For a relative test pass \code{gtol = x * abs(f)}.}
#'     \item{\code{xtol}}{Step norm. Default \code{0}, disabled.}
#'     \item{\code{rmin}}{Lower limit on the trust radius; falling below it
#'       stops the run with \code{converged = FALSE}. Default \code{0},
#'       disabled.}
#'     \item{\code{ctol}}{Multiple shooting only: scaled gap size below which
#'       the trajectory counts as continuous; no other test stops the run
#'       before. Default \code{1e-6}.}
#'   }
#' @param qnControl List of quasi-Newton settings, merged into the defaults and
#'   read only in a \code{"bfgs"} or \code{"sr1"} phase.
#'   \describe{
#'     \item{\code{hessianInit}}{Seed of the approximation: \code{"gn"}
#'       (default) the Gauss-Newton Hessian at \code{parinit};
#'       \code{"identity"} the identity matrix, after which no Hessian is
#'       requested from the objective; \code{"exact"} the true Hessian at
#'       \code{parinit}, then gradients only. \code{"exact"} needs
#'       \code{deriv2}, see \code{derivMode} in [odemodel()]. Used only when
#'       \code{hessianMethod} is \code{"bfgs"} or \code{"sr1"}.}
#'     \item{\code{hessianReseed}}{\code{"never"} (default): a stalled
#'       quasi-Newton phase stops. \code{"stall"}: it continues from a fresh
#'       Hessian of the same source at the current iterate, provided the
#'       iterate moved since the last reseed. The result reports the count as
#'       \code{nReseed}.}
#'     \item{\code{qnMemory}}{Number of \code{(s, y)} pairs kept. Default
#'       \code{0}: every update accumulates onto the seed. A positive value
#'       rebuilds the approximation in each iteration from the last
#'       \code{qnMemory} pairs, without the seed.}
#'     \item{\code{qnCautious}}{A BFGS pair is used only when
#'       \eqn{s^\top y > c\,|s|\,|y|} with \eqn{c} = \code{qnCautious}. Default
#'       \code{1e-8}, \code{0} disables the test; \code{"sr1"} is exempt. The
#'       result reports the rejected pairs as \code{qnSkipped}.}
#'     \item{\code{qnRejected}}{Whether \code{"sr1"} also updates from a
#'       rejected trial point. Default \code{TRUE} when \code{hessianMethod}
#'       or \code{hessianFallback} is \code{"sr1"}, \code{FALSE} otherwise.}
#'   }
#' @param stepControl List of step and acceptance settings, merged into the
#'   defaults.
#'   \describe{
#'     \item{\code{boundary}}{Box-bound handling, \code{"reflective"} (default)
#'       or \code{"clip"}. See section Box bounds.}
#'     \item{\code{thetaMax}}{Largest fraction of the distance to a bound a
#'       step may use. Default \code{0.99995}. \code{boundary = "reflective"}
#'       only.}
#'     \item{\code{nonmonotone}}{Zhang-Hager relaxation in \code{[0, 1)}: a
#'       step is rated against a weighted average of past objective values and
#'       may increase the objective. Default \code{0}, monotone.}
#'     \item{\code{acceptance}}{Multiple shooting only: \code{"filter"}
#'       (default), the filter of Fletcher and Leyffer on objective and gap
#'       size; \code{"merit"}, the l2 exact penalty; \code{"natural"}, Bock's
#'       natural level function, which needs \code{hessianMethod = "gn"}, a
#'       fixed sigma and no error model.}
#'     \item{\code{twoPhase}}{Multiple shooting only: if \code{TRUE}, the
#'       parameters are first fitted with the node values held fixed, then
#'       the nodes are laid out afresh and released. Default \code{FALSE}.}
#'     \item{\code{restore}}{Multiple shooting only: if \code{TRUE}, a run whose
#'       steps keep being rejected with open gaps enters a restoration phase,
#'       at most five times. Default \code{FALSE}.}
#'     \item{\code{regularise}}{Multiple shooting only: directions of the
#'       condensed problem with a singular value below \code{regularise} times
#'       the largest are held fixed in a step. Default \code{0}, off. With
#'       \code{acceptance = "natural"} directions below \code{1e-7} are always
#'       held fixed.}
#'     \item{\code{soc}}{Multiple shooting only: if \code{TRUE}, a rejected
#'       step whose gaps did not close as predicted is retried once with a
#'       second-order correction. Default \code{TRUE}.}
#'     \item{\code{anneal}}{Multiple shooting only: if \code{TRUE}, the run
#'       first minimises the objective plus a rising penalty on the gaps, then
#'       enforces continuity. Needs forward sensitivities and is skipped
#'       without them. Default \code{TRUE}.}
#'   }
#' @param minimize If \code{TRUE} (default) minimise, if \code{FALSE}
#'   maximise.
#' @param blather If \code{TRUE}, the result also contains the per-iteration
#'   trace listed under Value. Default \code{FALSE}.
#' @param printIter If \code{TRUE}, print iteration count and objective value
#'   at each evaluation. Default \code{FALSE}.
#' @param traceFile Path of a CSV file written with one row per evaluation:
#'   iteration, objective value and parameters. Default \code{NULL}, no file.
#' @param ... Further named arguments passed to \code{objfun}.
#'
#' @return A list of class \code{trustfit} with components
#'   \describe{
#'     \item{\code{argument}}{Named numeric, the best iterate.}
#'     \item{\code{value}, \code{gradient}, \code{hessian}}{Objective value,
#'       gradient and model Hessian at \code{argument}.}
#'     \item{\code{iterations}}{Number of iterations.}
#'     \item{\code{neval}}{Number of objective evaluations.}
#'     \item{\code{qnEval}}{Evaluations spent in a quasi-Newton phase.}
#'     \item{\code{qnSkipped}}{Pairs rejected by the cautious test.}
#'     \item{\code{nSwitch}}{Number of Hessian source handovers.}
#'     \item{\code{nReseed}}{Number of reseeds by
#'       \code{qnControl$hessianReseed = "stall"}.}
#'     \item{\code{evalBySource}}{Named integer, evaluations per Hessian
#'       source.}
#'     \item{\code{converged}}{Logical.}
#'     \item{\code{atBound}}{Named logical, the parameters held by a bound.}
#'     \item{\code{stopReason}}{The test that ended the run, see section
#'       Convergence.}
#'   }
#'   With \code{blather = TRUE} the list also contains the matrices
#'   \code{argpath} (accepted iterates) and \code{argtry} (trial points), one
#'   row per iteration, and the per-iteration vectors \code{steptype},
#'   \code{stepback}, \code{accept}, \code{r} (radius), \code{rho} (ratio of
#'   actual to predicted reduction), \code{valpath}, \code{valtry},
#'   \code{preddiff}, \code{stepnorm} and, when a quasi-Newton source is
#'   involved, \code{hessianSource}.
#'
#'   \code{boundary = "clip"} returns only \code{argument}, \code{value},
#'   \code{gradient}, \code{hessian}, \code{iterations}, \code{converged},
#'   \code{atBound} and \code{stopReason}, plus the trace.
#'
#'   A multiple-shooting run returns \code{argument}, \code{value},
#'   \code{gradient}, \code{hessian}, \code{iterations}, \code{neval},
#'   \code{qnSkipped}, \code{converged}, \code{stopReason}, \code{atBound} and
#'   a list \code{multipleShooting} with the final \code{nodes}, the
#'   \code{gaps} per condition, their scaled size \code{violation}, the counts
#'   \code{nRelaxed} (iterations with a node step below the full one),
#'   \code{nSOC} (accepted second-order corrections), \code{nSplit} (segments
#'   cut), \code{nRestore} (restoration phases) and \code{nStages} (annealing
#'   stages), the \code{cuts}, the \code{acceptance} and \code{sweep} used,
#'   \code{mu} (merit) or \code{filterSize} (filter), and \code{trace} when
#'   \code{blather = TRUE}.
#'
#' @seealso [mstrust()] for multi-start fits, [profile()][profile.objfn] for
#'   profile likelihoods, [normL2()] for the objective,
#'   \code{vignette("Optimisation", package = "dMod2")} for the method.
#'
#' @examples
#' rosenbrock <- function(x, ...) {
#'   a <- x[["a"]]
#'   b <- x[["b"]]
#'   list(value = 100 * (b - a^2)^2 + (1 - a)^2,
#'        gradient = c(a = -400 * a * (b - a^2) - 2 * (1 - a),
#'                     b = 200 * (b - a^2)),
#'        hessian = matrix(c(1200 * a^2 - 400 * b + 2, -400 * a, -400 * a, 200),
#'                         2, 2, dimnames = list(c("a", "b"), c("a", "b"))))
#' }
#'
#' fit <- trust(rosenbrock, c(a = -1.2, b = 1), iterlim = 200)
#' fit$argument
#' fit$stopReason
#'
#' # With an upper bound on a
#' fit <- trust(rosenbrock, c(a = -1.2, b = 1), parupper = c(a = 0.5),
#'              iterlim = 200)
#' fit$atBound
#'
#' # Quasi-Newton with a Gauss-Newton fallback
#' fit <- trust(rosenbrock, c(a = -1.2, b = 1), hessianMethod = "bfgs",
#'              hessianFallback = "gn", iterlim = 200)
#' fit$evalBySource
#'
#' @references Coleman, T. F. and Li, Y. (1996). An interior trust region
#'   approach for nonlinear minimization subject to bounds.
#'   \emph{SIAM Journal on Optimization} 6(2), 418-445.
#'
#'   Fröhlich, F. and Sorger, P. K. (2022). Fides: Reliable trust-region
#'   optimization for parameter estimation of ODE models.
#'   \emph{PLoS Computational Biology} 18(7), e1010322.
#'
#'   Bock, H. G. (1983). Recent advances in parameter identification
#'   techniques for ODE. In \emph{Numerical Treatment of Inverse Problems in
#'   Differential and Integral Equations}, 95-121. Birkhäuser.
#'
#'   Fletcher, R. and Leyffer, S. (2002). Nonlinear programming without a
#'   penalty function. \emph{Mathematical Programming} 91, 239-269.
#'
#'   Horbelt, W., Timmer, J. and Voss, H. U. (2002). Parameter estimation in
#'   nonlinear delayed feedback systems from noisy data.
#'   \emph{Physics Letters A} 299, 513-521.
#'
#'   Voss, H. U., Timmer, J. and Kurths, J. (2004). Nonlinear dynamical system
#'   identification from uncertain and indirect measurements.
#'   \emph{International Journal of Bifurcation and Chaos} 14(6), 1905-1933.
#'
#'   Peifer, M. and Timmer, J. (2007). Parameter estimation in ordinary
#'   differential equations for biochemical processes using the method of
#'   multiple shooting. \emph{IET Systems Biology} 1(2), 78-88.
#'
#' @export
trust <- function(objfun, parinit, rinit = 0.1, rmax = 10,
                  iterlim   = 100L,
                  hessianMethod   = c("gn", "bfgs", "sr1", "exact"),
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
  if (!is.null(parscale))
    warning("trust: 'parscale' is deprecated; fit on log scale instead.",
            call. = FALSE)
  stepControl <- .stepControlAlias(stepControl, "trust")

  ctl  <- .trustDefaults(hessianMethod, hessianFallback)
  tol  <- .mergeControl(ctl$tolControl,  tolControl,  "tolControl")
  qn   <- .mergeControl(ctl$qnControl,   qnControl,   "qnControl")
  step <- .mergeControl(ctl$stepControl, stepControl, "stepControl")

  boundary    <- match.arg(step$boundary, c("reflective", "clip"))
  hessianInit   <- match.arg(qn$hessianInit, c("gn", "identity", "exact"))
  hessianReseed <- match.arg(qn$hessianReseed %||% "never", c("never", "stall"))

  # The objective decides: a normL2 with a multipleShootingControl, alone or
  # in a sum, is optimised by multiple shooting, with the same controls.
  if (length(.shootingTerms(objfun)))
    return(.trustFit(.trustShooting(objfun, parinit, rinit, rmax, iterlim,
                          hessianMethod, hessianFallback, fallbackLimit,
                          parscale, parupper, parlower, tol, qn, step,
                          minimize, blather, printIter, traceFile, dots)))

  # The kernel names what it wants, 0 value, 1 gradient, 2 Gauss-Newton, 3
  # exact, and the translation into an objective's own arguments happens here,
  # where its formals are visible. One that only knows `hessian` gets the
  # logical it always got. `sweep` is not chosen here; it stays the caller's.
  #
  # The exact request is checked before the run rather than inside fn, where a
  # stop() is caught by the kernel's evaluation handler and reported as
  # "parinit not feasible".
  .fml   <- names(formals(objfun))
  .wants_exact <- identical(hessianMethod, "exact") ||
                  identical(hessianInit, "exact")
  if (.wants_exact && !any(c("deriv2", "...") %in% .fml))
    stop("trust: hessianMethod or qnControl$hessianInit is \"exact\", but ",
         "this objective has no `deriv2` argument to ask an exact Hessian ",
         "with.", call. = FALSE)

  fn <- function(x, want = 2L) {
    want <- as.integer(want)
    args <- list(x, deriv = want >= 1L, hessian = want >= 2L)
    if (want >= 3L) args$deriv2 <- TRUE
    do.call(objfun, c(args, dots))
  }
  .trustFit(trust_impl(fn, parinit, rinit, rmax, parscale, as.integer(iterlim),
             tol$ftol, tol$mtol, tol$gtol, tol$xtol, tol$rmin, step$thetaMax,
             boundary, hessianMethod, hessianFallback,
             as.integer(fallbackLimit), hessianInit, hessianReseed,
             as.integer(qn$qnMemory), qn$qnCautious, isTRUE(qn$qnRejected),
             step$nonmonotone, minimize, blather,
             parupper, parlower, printIter, traceFile))
}

# A fit is a list of class `trustfit`, so vcov() dispatches on it.
.trustFit <- function(fit) {
  if (!is.list(fit) || inherits(fit, "trustfit")) return(fit)
  class(fit) <- c("trustfit", "list")
  fit
}

#' @export
print.trustfit <- function(x, ...) {
  print(unclass(x), ...)
  invisible(x)
}


# Defaults per Hessian source, and the single source of truth for what each
# control group accepts. Members of one group are settable together or not at
# all; see the roxygen of `trust` for what they mean.
.trustDefaults <- function(hessianMethod = "gn", hessianFallback = "none") {
  qn <- list(hessianInit = "gn", hessianReseed = "never",
             qnMemory = 0L, qnCautious = 1e-8, qnRejected = FALSE)
  # Nocedal and Wright, sec. 6.2: SR1 also updates from a rejected step. The
  # flag is inert outside an sr1 phase, so a fallback to sr1 arms it too.
  if ("sr1" %in% c(hessianMethod, hessianFallback)) qn$qnRejected <- TRUE
  list(
    tolControl  = list(ftol = 1e-6, mtol = 1e-6, gtol = 1e-6,
                       xtol = 0, rmin = 0, ctol = 1e-6),
    qnControl   = qn,
    stepControl = list(boundary = "reflective", thetaMax = 0.99995,
                       nonmonotone = 0, acceptance = "filter", soc = TRUE,
                       anneal = TRUE, twoPhase = FALSE, regularise = 0,
                       restore = FALSE))
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
  hessianReseed = "qnControl$hessianReseed",
  qnCautious  = "qnControl$qnCautious",  qnRejected = "qnControl$qnRejected",
  boundary    = "stepControl$boundary",  thetaMax   = "stepControl$thetaMax",
  theta.max   = "stepControl$thetaMax",  nonmonotone = "stepControl$nonmonotone")

# `theta.max`, the former name of stepControl$thetaMax, with a warning.
.stepControlAlias <- function(stepControl, who) {
  if (!"theta.max" %in% names(stepControl)) return(stepControl)
  if ("thetaMax" %in% names(stepControl))
    stop(who, ": give stepControl$thetaMax only; 'theta.max' is its ",
         "deprecated name.", call. = FALSE)
  warning(who, ": stepControl$theta.max is deprecated, use thetaMax.",
          call. = FALSE)
  names(stepControl)[names(stepControl) == "theta.max"] <- "thetaMax"
  stepControl
}

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


# Merge a user control list into `defaults` for a trust() call. Names are checked
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
