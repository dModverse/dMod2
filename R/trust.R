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
#' Only \code{gtol} tests first-order optimality. While \code{switches}
#' remain, see \code{phases}, an \code{"fvalue"}, \code{"preddiff"},
#' \code{"step"} or \code{"stagnation"} stop hands over to the next phase
#' instead of ending the run. In a quasi-Newton phase \code{"fvalue"},
#' \code{"preddiff"} and \code{"step"} do not end the run. A phase that
#' reaches its own \code{iterlim} ends like a soft stop, with
#' \code{stopReason = "iterlim"} and \code{converged = FALSE} when no
#' handover is left. The total \code{iterlim} ends the run in any phase, with
#' \code{stopReason = "iterlim"} and \code{converged = FALSE}.
#'
#' \code{argument}, \code{value}, \code{gradient}, \code{hessian} and
#' \code{atBound} in the result belong to the best iterate visited, which is
#' the last accepted one when \code{stepControl$nonmonotone = 0}.
#'
#' @section Multiple shooting:
#' An objective built by \code{\link{normL2}(..., multipleShootingControl = )},
#' alone or in a sum, is optimised by multiple shooting with the same controls.
#' Not available then: \code{hessianMethod = "exact"}, more than one phase,
#' \code{switches}, \code{qnControl$hessianReseed},
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
#' @param iterlim Maximum number of trust-region iterations over all phases.
#'   Default \code{1000}.
#' @param hessianMethod Source of the model Hessian of a run without
#'   \code{phases}: \code{"gn"} (default) the Gauss-Newton Hessian of the
#'   objective, \code{"bfgs"} or \code{"sr1"} a quasi-Newton update seeded as
#'   \code{qnControl$hessianInit} says, \code{"exact"} the true Hessian at
#'   every iterate, which needs an objective that takes \code{deriv2}.
#'   Quasi-Newton methods require \code{stepControl$boundary = "reflective"}.
#'   Also selects the control defaults below. With \code{phases}, each phase
#'   names its own and this argument is an error.
#' @param phases List of phases, each a list of settings, run one after the
#'   other and from the first again after the last. The run hands over to the
#'   next phase at a value, model, step or stagnation stop and when a phase
#'   uses up its own \code{iterlim}. A handover restarts the trust region at
#'   \code{rinit} and discards the quasi-Newton pairs. Default \code{NULL}: one
#'   phase from \code{hessianMethod} and \code{qnControl}, evaluated with the
#'   \code{sweep} in \code{...}. Names of the list label the phases in the
#'   result, the Hessian source does otherwise. Not available with multiple
#'   shooting and \code{boundary = "clip"}.
#'   \describe{
#'     \item{\code{hessianMethod}}{Hessian source of the phase, as
#'       \code{hessianMethod}. Required.}
#'     \item{\code{iterlim}}{Maximum number of iterations of one run of the
#'       phase. Default \code{NULL}: only the total \code{iterlim}.}
#'     \item{\code{sweep}}{\code{"forward"} or \code{"reverse"}, the
#'       \code{sweep} the objective is evaluated with in the phase; a
#'       \code{sweep} in \code{...} is an error with \code{phases}. Default
#'       \code{NULL}: for \code{"bfgs"} and \code{"sr1"} reverse if the
#'       objective can be evaluated backwards at \code{parinit}, as with
#'       \code{derivMode = c("forward", "reverse")} throughout the chain,
#'       forward otherwise. A Gauss-Newton Hessian is always computed forward.
#'       A handover between phases of different sweep evaluates the current
#'       iterate once more. The result reports the sweeps as
#'       \code{phaseSweep}.}
#'     \item{\code{qnControl}}{Quasi-Newton settings of the phase, merged
#'       over \code{qnControl}; see there. \code{hessianInit} defaults to
#'       \code{"gn"} in every phase but the first: a quasi-Newton phase starts
#'       from the Gauss-Newton Hessian at the iterate of the handover.}
#'   }
#' @param switches Number of handovers allowed. Default \code{NULL}: one
#'   fewer than there are phases, so each runs once. A higher value cycles
#'   through the phases again, \code{0} runs the first phase only. The result
#'   reports the count as \code{nSwitch}.
#' @param parupper,parlower Upper and lower bounds. Unnamed: the first element
#'   applies to all parameters. Named: entries are matched by name, missing
#'   parameters are unbounded. Default \code{NULL}, unbounded.
#' @param tolControl List of convergence thresholds, merged into the defaults.
#'   \describe{
#'     \item{\code{ftol}}{Change in objective value. Default \code{1e-6}.}
#'     \item{\code{mtol}}{Predicted model decrease. Default \code{1e-6}.}
#'     \item{\code{gtol}}{First-order optimality measure. Default \code{1e-6}.}
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
#'       rejected trial point. Default \code{TRUE} for a phase whose source is
#'       \code{"sr1"}, \code{FALSE} otherwise.}
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
#'     \item{\code{nSwitch}}{Number of handovers between phases.}
#'     \item{\code{phaseIterations}}{Named integer, the iterations of each
#'       phase run, in run order and named after the phase; sums to
#'       \code{iterations}.}
#'     \item{\code{phaseSweep}}{Named character, the sweep of each element
#'       of \code{phases}, present when \code{phases} is set.}
#'     \item{\code{nReseed}}{Number of reseeds by
#'       \code{qnControl$hessianReseed = "stall"}, which allows each run of a
#'       phase the same number.}
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
#' # Gauss-Newton, then BFGS with a budget of its own, then Gauss-Newton again
#' fit <- trust(rosenbrock, c(a = -1.2, b = 1), iterlim = 200,
#'              phases = list(list(hessianMethod = "gn"),
#'                            list(hessianMethod = "bfgs", iterlim = 50)),
#'              switches = 2)
#' fit$phaseIterations
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
                  iterlim   = 1000L,
                  hessianMethod   = c("gn", "bfgs", "sr1", "exact"),
                  phases    = NULL,
                  switches  = NULL,
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
  dots <- list(...)
  .trustRejectMoved(names(dots), "trust")
  if (!is.null(phases) && !missing(hessianMethod))
    stop("trust: with `phases`, each phase names its own hessianMethod.",
         call. = FALSE)
  if (!is.null(phases) && !is.null(dots$sweep))
    stop("trust: with `phases`, each phase sets its own sweep.", call. = FALSE)
  hessianMethod <- if (is.null(phases)) match.arg(hessianMethod) else
    .phaseMethod(phases[[1L]]$hessianMethod, 1L)

  ctl  <- .trustDefaults(hessianMethod)
  tol  <- .mergeControl(ctl$tolControl,  tolControl,  "tolControl")
  qn   <- .mergeControl(ctl$qnControl,   qnControl,   "qnControl")
  step <- .mergeControl(ctl$stepControl, stepControl, "stepControl")

  boundary    <- match.arg(step$boundary, c("reflective", "clip"))
  ph <- .trustPhases(phases, switches, hessianMethod, qnControl)

  # The objective decides: a normL2 with a multipleShootingControl, alone or
  # in a sum, is optimised by multiple shooting, with the same controls.
  if (length(.shootingTerms(objfun)) > 0L) {
    if (length(ph$phases) > 1L || ph$switches > 0L)
      stop("trust: a multiple-shooting objective runs a single phase; ",
           "`phases` with more than one element and `switches` are not ",
           "supported.", call. = FALSE)
    return(.trustFit(.trustShooting(objfun, parinit, rinit, rmax, iterlim,
                          hessianMethod, parupper, parlower, tol,
                          ph$phases[[1L]]$qnControl, step,
                          minimize, blather, printIter, traceFile, dots)))
  }

  # `want` (0 value, 1 gradient, 2 Gauss-Newton, 3 exact) becomes the objective's
  # own arguments here, where its formals are visible. The exact request is checked
  # up front, since a stop() inside fn reads as "parinit not feasible".
  .fml   <- names(formals(objfun))
  .wants_exact <- any(vapply(ph$phases, function(p)
    identical(p$hessianMethod, "exact") || identical(p$qnControl$hessianInit, "exact"),
    TRUE))
  if (.wants_exact && !any(c("deriv2", "...") %in% .fml))
    stop("trust: hessianMethod or qnControl$hessianInit is \"exact\", but ",
         "this objective has no `deriv2` argument to ask an exact Hessian ",
         "with.", call. = FALSE)

  # Without `phases` the objective takes `...` as given. With `phases` each
  # phase passes its own sweep in place of one in `...`, and a request for a
  # Gauss-Newton Hessian runs forward, since it needs the forward sensitivities.
  sweeps <- if (is.null(phases)) dots$sweep %||% "forward" else {
    probe <- local({
      v <- NULL
      function() { if (is.null(v)) v <<- .reverseProbe(objfun, parinit, dots); v }
    })
    vapply(ph$phases, function(p) .phaseSweep(p$sweep, p$hessianMethod, objfun, probe), "")
  }
  phaseDots <- lapply(sweeps, function(s) {
    if (is.null(phases)) return(dots)
    d <- dots
    d$sweep <- if (identical(s, "reverse")) "reverse" else NULL
    d
  })

  fn <- function(x, want = 2L, phase = 0L) {
    want <- as.integer(want)
    args <- list(x, deriv = want >= 1L, hessian = want >= 2L)
    if (want >= 3L) args$deriv2 <- TRUE
    d <- phaseDots[[phase + 1L]]
    if (!is.null(phases) && want == 2L) d$sweep <- NULL
    .tryObjective(do.call(objfun, c(args, d)))
  }
  qfield <- function(f) lapply(ph$phases, function(p) p$qnControl[[f]])
  fit <- trust_impl(fn, parinit, rinit, rmax, NULL, as.integer(iterlim),
             tol$ftol, tol$mtol, tol$gtol, tol$xtol, tol$rmin, step$thetaMax,
             boundary,
             vapply(ph$phases, `[[`, "", "hessianMethod"),
             vapply(ph$phases, `[[`, 1L, "iterlim"),
             sweeps == "reverse",
             vapply(qfield("hessianInit"), function(v)
               match.arg(v, c("gn", "identity", "exact")), ""),
             vapply(qfield("hessianReseed"), function(v)
               match.arg(v %||% "never", c("never", "stall")), ""),
             vapply(qfield("qnMemory"), as.integer, 1L),
             vapply(qfield("qnCautious"), as.numeric, 1),
             vapply(qfield("qnRejected"), isTRUE, TRUE),
             ph$switches,
             step$nonmonotone, minimize, blather,
             parupper, parlower, printIter, traceFile)
  if (!is.null(fit$phaseIterations)) {
    n <- length(fit$phaseIterations)
    names(fit$phaseIterations) <- ph$labels[(seq_len(n) - 1L) %% length(ph$labels) + 1L]
  }
  if (!is.null(phases)) fit$phaseSweep <- stats::setNames(sweeps, ph$labels)
  .trustFit(fit)
}


.phaseMethod <- function(method, k) {
  if (is.null(method))
    stop("trust: phases[[", k, "]] names no hessianMethod.", call. = FALSE)
  match.arg(method, c("gn", "bfgs", "sr1", "exact"))
}

# The phases of a run, each with its Hessian source, budget (NA: none of its
# own), sweep and quasi-Newton settings merged over the defaults of its source
# and the top-level qnControl. Without `phases`, one phase from the top level.
.trustPhases <- function(phases, switches, hessianMethod, qnControl) {
  if (is.null(phases))
    phases <- list(list(hessianMethod = hessianMethod))
  else if (!is.list(phases) || !length(phases) ||
           !all(vapply(phases, is.list, TRUE)))
    stop("trust: `phases` must be a non-empty list of lists.", call. = FALSE)

  defaults <- list(hessianMethod = NULL, iterlim = NULL, sweep = NULL,
                   qnControl = NULL)
  out <- lapply(seq_along(phases), function(k) {
    p <- .mergeControl(defaults, phases[[k]], paste0("phases[[", k, "]]"))
    p$hessianMethod <- .phaseMethod(p$hessianMethod, k)
    lim <- p$iterlim %||% Inf
    if (length(lim) != 1L || is.na(lim) || lim < 1)
      stop("trust: phases[[", k, "]]$iterlim must be a positive number.",
           call. = FALSE)
    p$iterlim <- if (is.finite(lim)) as.integer(lim) else NA_integer_
    q <- .mergeControl(.trustDefaults(p$hessianMethod)$qnControl, qnControl,
                       "qnControl")
    # A later phase starts from the Gauss-Newton Hessian at the handover.
    if (k > 1L) q$hessianInit <- "gn"
    p$qnControl <- .mergeControl(q, p$qnControl,
                                 paste0("phases[[", k, "]]$qnControl"))
    p
  })

  switches <- switches %||% (length(out) - 1L)
  if (length(switches) != 1L || is.na(switches) || switches < 0)
    stop("trust: `switches` must be a non-negative integer.", call. = FALSE)

  labels <- names(phases) %||% rep("", length(out))
  labels <- ifelse(nzchar(labels), labels, vapply(out, `[[`, "", "hessianMethod"))
  list(phases = out, switches = as.integer(switches),
       labels = make.unique(labels, sep = "."))
}

# The sweep of one phase. NULL takes the reverse sweep for a quasi-Newton phase
# when the objective can do it at parinit, forward otherwise; an explicit
# "reverse" has to work. `probe()` is TRUE or the message of the failed probe.
.phaseSweep <- function(sweep, method, objfun, probe) {
  if (!is.null(sweep)) sweep <- match.arg(sweep, c("forward", "reverse"))
  if (identical(method, "gn")) {
    if (identical(sweep, "reverse"))
      stop("trust: a Gauss-Newton phase needs the forward sensitivities; ",
           "sweep = \"reverse\" is for \"bfgs\", \"sr1\" and \"exact\".",
           call. = FALSE)
    return("forward")
  }
  if (identical(sweep, "forward")) return("forward")
  if (is.null(sweep) && identical(method, "exact")) return("forward")
  if (!"sweep" %in% names(formals(objfun))) {
    if (is.null(sweep)) return("forward")
    stop("trust: a phase asks for sweep = \"reverse\", but the objective has ",
         "no `sweep` argument.", call. = FALSE)
  }
  ok <- probe()
  if (isTRUE(ok)) return("reverse")
  if (is.null(sweep)) return("forward")
  stop("trust: a phase asks for sweep = \"reverse\", but the objective cannot ",
       "be evaluated backwards: ", ok, call. = FALSE)
}

# TRUE if the objective evaluates value and gradient backwards at parinit, the
# error message otherwise.
.reverseProbe <- function(objfun, parinit, dots) {
  dots$sweep <- "reverse"
  tryCatch({
    suppressWarnings(do.call(objfun, c(list(parinit, deriv = TRUE,
                                            hessian = FALSE), dots)))
    TRUE
  }, error = function(e) conditionMessage(e))
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
.trustDefaults <- function(hessianMethod = "gn") {
  qn <- list(hessianInit = "gn", hessianReseed = "never",
             qnMemory = 0L, qnCautious = 1e-8, qnRejected = FALSE)
  # Nocedal and Wright, sec. 6.2: SR1 also updates from a rejected step.
  if (identical(hessianMethod, "sr1")) qn$qnRejected <- TRUE
  list(
    tolControl  = list(ftol = 1e-6, mtol = 1e-6, gtol = 1e-6,
                       xtol = 0, rmin = 0, ctol = 1e-6),
    qnControl   = qn,
    stepControl = list(boundary = "reflective", thetaMax = 0.99995,
                       nonmonotone = 0, acceptance = "filter", soc = TRUE,
                       anneal = TRUE, twoPhase = FALSE, regularise = 0,
                       restore = FALSE))
}


# Members of a control list, refused by name when given flat instead of being
# routed to objfun.
.trustMoved <- c(
  ftol = "tolControl$ftol",   mtol  = "tolControl$mtol",
  gtol = "tolControl$gtol",   xtol  = "tolControl$xtol",
  rmin = "tolControl$rmin",
  hessianInit = "qnControl$hessianInit", qnMemory   = "qnControl$qnMemory",
  hessianReseed = "qnControl$hessianReseed",
  qnCautious  = "qnControl$qnCautious",  qnRejected = "qnControl$qnRejected",
  boundary    = "stepControl$boundary",  thetaMax   = "stepControl$thetaMax",
  nonmonotone = "stepControl$nonmonotone")

.trustRejectMoved <- function(nms, who) {
  gone <- intersect(nms, c("fallbackControl", "hessianFallback", "fallbackLimit"))
  if (length(gone))
    stop(who, ": ", paste(gone, collapse = ", "), " is no argument; a run over ",
         "several Hessian sources takes `phases` and `switches`, see ?trust.",
         call. = FALSE)
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
#' @param gate Character, names in \code{mu} whose penalty is the gate
#'   \code{lambda * (p - mu)} with \code{p >= mu}: the parameter can sit exactly
#'   on \code{mu} and leaves it only upward. Starts below \code{mu} are moved
#'   onto it.
#' @param fuse List of fusion blocks. A block penalises
#'   \code{lambda * sum_{i<j} w_ij |p_i - p_j|} over its parameters, plus
#'   \code{lambda * w_ia |p_i - a|} toward an optional anchor \code{a}.
#'   Parameters with equal values move as one and split again when the data
#'   pull them apart. A block is a character vector of parameter names or a
#'   list with \code{pars}, and optionally \code{lambda} (default: the scalar
#'   \code{lambda}), \code{anchor} (default none) and \code{weights}, an
#'   \code{(m+1) x (m+1)} matrix whose last row and column refer to the
#'   anchor (default all 1). Fused parameters must not appear in \code{mu}; if
#'   \code{fuse} is given and \code{mu} is not, \code{mu} is empty.
trustL1 <- function(objfun, parinit, mu = 0 * parinit, one.sided = FALSE, lambda = 1,
                    rinit = 0.1, rmax = 10,
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
                    gate = NULL, fuse = NULL,
                    ...) {
  boundary <- match.arg(boundary)

  sanePars <- sanitizePars(parinit, list(...)$fixed)
  parinit  <- sanePars$pars

  if (is.null(names(parinit)))
    stop("trustL1: parinit must be a named numeric vector")
  if (missing(mu) && !is.null(fuse)) mu <- structure(numeric(0), names = character(0))
  if (length(mu) > 0L && is.null(names(mu)))
    stop("trustL1: mu must be a named numeric vector")

  unknown <- setdiff(names(mu), names(parinit))
  if (length(unknown) > 0L)
    stop("trustL1: mu has names not present in parinit: ",
         paste(unknown, collapse = ", "))

  lambda0 <- if (length(lambda) == 1L) as.numeric(lambda) else NULL
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

  # The kernel passes what it wants as a second argument; trustL1 always
  # takes value, gradient and Hessian, so the request is dropped here.
  dots <- list(...)
  fn <- function(x, want = 2L) .tryObjective(do.call(objfun, c(list(x), dots)))

  mu <- structure(as.numeric(mu), names = names(mu))
  lambda <- structure(as.numeric(lambda), names = names(lambda))

  if (!is.null(gate)) {
    gate <- as.character(gate)
    if (!all(gate %in% names(mu)))
      stop("trustL1: every gate must be named in mu")
    low <- gate[parinit[gate] < mu[gate]]
    parinit[low] <- mu[low]
  }
  if (!is.null(fuse)) fuse <- .fuseBlocks(fuse, names(parinit), lambda0)

  # Maximising f minus the penalty is minimising -f plus the penalty; the kernel
  # would negate the penalty with f.
  flip <- function(out) {
    if (inherits(out, "try-error")) return(out)
    for (k in c("value", "gradient", "hessian")) if (!is.null(out[[k]])) out[[k]] <- -out[[k]]
    out
  }
  if (!minimize) {
    fnMax <- fn
    fn <- function(x, want = 2L) flip(fnMax(x, want))
  }

  out <- trustL1_impl(fn, parinit, mu, lambda,
                      as.logical(one.sided)[1L], rinit, rmax,
                      NULL, as.integer(iterlim),
                      ftol, mtol, gtol, xtol, rmin, theta.max,
                      boundary, TRUE, blather,
                      parupper, parlower, printIter, traceFile,
                      gate = gate, fuse = fuse)
  if (minimize) out else flip(out)
}

# Normalise trustL1 fuse blocks to list(idx, lambda, anchor, w) with idx
# positions in `parnames`.
.fuseBlocks <- function(fuse, parnames, lambda0) {
  if (!is.list(fuse)) fuse <- list(fuse)
  lapply(fuse, function(b) {
    if (is.character(b)) b <- list(pars = b)
    pars <- b$pars
    if (length(pars) < 1L || anyDuplicated(pars))
      stop("trustL1: a fuse block needs unique parameter names")
    idx <- match(pars, parnames)
    if (anyNA(idx))
      stop("trustL1: fuse names not present in parinit: ",
           paste(pars[is.na(idx)], collapse = ", "))
    lam <- if (is.null(b$lambda)) lambda0 else b$lambda
    if (length(lam) != 1L || !is.finite(lam) || lam < 0)
      stop("trustL1: a fuse block needs one non-negative lambda")
    m <- length(pars)
    w <- if (is.null(b$weights)) matrix(1, m + 1L, m + 1L) else b$weights
    if (!is.matrix(w) || any(dim(w) != m + 1L))
      stop("trustL1: fuse weights must be a (m+1) x (m+1) matrix")
    list(idx = as.integer(idx), lambda = as.numeric(lam),
         anchor = if (is.null(b$anchor)) NA_real_ else as.numeric(b$anchor),
         w = matrix(as.numeric(w), m + 1L))
  })
}


# Merge user control into `defaults` for trust()/trustL1(). Names are checked
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


## An R error reaches the kernel as a try-error with its message; thrown
## through Rcpp, the message would be lost.
.tryObjective <- function(expr)
  tryCatch(expr, error = function(e) structure(conditionMessage(e), class = "try-error"))
