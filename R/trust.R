#' Non-linear optimisation via a trust-region method
#'
#' \code{trust} minimises (or maximises) a smooth objective function for which
#' value, gradient and Hessian are available. \code{trustL1} additionally adds
#' an L1 (lasso-style) penalty \code{lambda * sum(|p - mu|)} on a user-selected
#' subset of parameters, with kink-aware step handling. Both routines solve a
#' Moré-Sorensen trust-region subproblem exactly, via the full eigendecomposition
#' of the Hessian.
#'
#' @section Box bounds:
#' \code{boundary = "reflective"} (default) uses the Coleman-Li interior
#' trust-region-reflective scheme: the subproblem is solved in the scaled frame
#' \code{D g}, \code{D H D + C} with \code{D = diag(|v|^(1/2))}, \code{|v_i|}
#' being the distance to the bound the gradient pushes toward. A step leaving
#' the box is truncated, reflected off the blocking faces, or replaced by the
#' scaled steepest-descent step, whichever the model rates best. Iterates stay
#' strictly inside the box and never land exactly on a bound -- use
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
#' @param ftol Convergence threshold on the change in objective value.
#' @param mtol Convergence threshold on the predicted model decrease.
#' @param gtol Convergence threshold on the first-order optimality measure.
#'   Deliberately absolute: no scaling of a gradient by the objective value is
#'   meaningful for a general objective, and for a sum of squares it is actively
#'   wrong -- \code{|f|} grows quadratically in the residuals while \code{|g|}
#'   grows linearly. Pass \code{gtol = x * |f|} explicitly if you want that.
#' @param xtol Convergence threshold on the step norm. \code{0} disables it.
#' @param rmin Lower limit on the trust radius. Falling below it stops the run
#'   with \code{converged = FALSE}. \code{0} disables it.
#' @param theta.max Largest fraction of the distance to a bound that a step may
#'   use, keeping iterates strictly interior. Only used by
#'   \code{boundary = "reflective"}.
#' @param boundary Box-bound handling, \code{"reflective"} (default) or
#'   \code{"clip"}. See Details.
#' @param minimize If \code{TRUE} (default) minimise; if \code{FALSE}
#'   maximise.
#' @param blather If \code{TRUE} return the per-iteration trace
#'   (\code{argpath}, \code{argtry}, \code{steptype}, \code{stepback},
#'   \code{accept}, \code{r}, \code{rho}, \code{valpath}, \code{valtry},
#'   \code{preddiff}, \code{stepnorm}).
#' @param parupper,parlower Named or scalar numeric bounds. If unnamed,
#'   the first element broadcasts to all parameters; if named, the
#'   entries slot by name into a length-K vector defaulting to
#'   \code{+/- Inf}.
#' @param printIter If \code{TRUE} print iteration count and objective
#'   value to the console at each function evaluation.
#' @param traceFile Optional path. If non-\code{NULL}, CSV-log per
#'   evaluation \code{iter, value, p1, p2, ...}.
#' @param fterm,mterm Deprecated aliases for \code{ftol} and \code{mtol}.
#' @param ... Additional named arguments forwarded to \code{objfun}.
#'
#' @return A list with components \code{argument}, \code{value},
#'   \code{gradient}, \code{hessian}, \code{iterations},
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

  dots <- list(...)
  fn <- if (length(dots) > 0L) {
    function(x) do.call(objfun, c(list(x), dots))
  } else {
    objfun
  }
  trust_impl(fn, parinit, rinit, rmax, parscale, as.integer(iterlim),
             ftol, mtol, gtol, xtol, rmin, theta.max,
             boundary, minimize, blather,
             parupper, parlower, printIter, traceFile)
}


#' @export
#' @rdname trust
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
# are checked against the optimiser's formals, so every optimiser argument is
# settable and a typo errors here instead of silently reaching objfun via `...`.
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
