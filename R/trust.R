#' Non-Linear Optimization by a Trust-Region Method
#'
#' \code{trust} minimizes (or maximizes) a smooth objective function for which
#' value, gradient and Hessian are available. It solves a Moré-Sorensen
#' trust-region subproblem exactly, via the full eigendecomposition of the
#' Hessian.
#'
#' @section Default method:
#' Without \code{hessianMethod}, \code{phases}, \code{switches} and a
#' \code{sweep} in \code{...}, \code{trust} runs the hybrid method of Fides
#' (Fröhlich and Sorger 2022): Gauss-Newton until the trust radius stayed the
#' same for 20 consecutive iterations, then BFGS for the rest of the run. The
#' BFGS approximation starts from the Gauss-Newton Hessian at \code{parinit},
#' made positive definite, takes every accepted step of the Gauss-Newton phase
#' as an update, and the trust radius stays at the switch. These are the
#' phases
#' \preformatted{list(gn = list(hessianMethod = "gn", radiusStall = 20, onStop = "end"),
#'      bfgs = list(hessianMethod = "bfgs", sweep = "forward", keepRadius = TRUE,
#'                  qnControl = list(hessianInit = "tracked")))}
#' The run is Gauss-Newton alone, as with \code{hessianMethod = "gn"}, where
#' phases are not available (multiple shooting, \code{boundary = "clip"}) and
#' on the sparse path, see \code{sparse}.
#'
#' @section Box bounds:
#' \code{boundary = "reflective"} (default) uses the interior
#' trust-region-reflective scheme of Coleman and Li. Iterates stay strictly
#' inside the box and never land exactly on a bound; \code{atBound} in the
#' result reports which parameters a bound holds.
#'
#' \code{boundary = "clip"}: active-set reduction with componentwise clipping,
#' with a single \code{"gn"} phase. The run stops when a trial point, accepted
#' or not, changes the value by less than \code{ftol} or has a predicted
#' reduction below \code{mtol}, both reported as \code{stopReason = "fvalue"}.
#' \code{gtol}, \code{xtol}, \code{rmin} and stagnation are not tested; three
#' failed evaluations in a row end the run with \code{stopReason = "iterlim"}.
#'
#' @section Convergence:
#' With \code{boundary = "reflective"} the run stops as soon as one of the
#' following holds. The gradient test is
#' evaluated at the current iterate; the value and step tests count only on an
#' accepted step, the model test also on a rejected Newton step from inside the
#' trust region.
#' \itemize{
#'   \item \code{gtol}: the first-order optimality measure
#'     \eqn{\max_i |v_i g_i|} is at most \code{gtol}, with \eqn{v_i} the
#'     distance to the bound in the descent direction \eqn{-g_i}, 1 if that
#'     bound is infinite. Without bounds it is \eqn{\max_i |g_i|}.
#'   \item \code{ftol}: \eqn{|f - f_{\mathrm{try}}| <} \code{ftol}.
#'   \item \code{mtol}: the predicted reduction
#'     \eqn{|g^\top p + \frac{1}{2} p^\top H p| <} \code{mtol}.
#'   \item \code{xtol}: the step norm is below \code{xtol}.
#'   \item stagnation: five consecutive rejected steps left the objective
#'     flat within \code{ftol}.
#'   \item \code{rmin}: the trust radius is below \code{rmin}, or a step
#'     no longer changes the parameters in floating point. Reports
#'     \code{converged = FALSE}.
#' }
#' A \code{"gradient"}, \code{"fvalue"} or \code{"preddiff"} stop that ends the run
#' also needs the full Newton step of the model at the current iterate to
#' predict a decrease of at most \code{mtol},
#' \eqn{\frac{1}{2} g^\top H^{+} g \le} \code{mtol}, in the variables scaled
#' for the bounds. \eqn{H^{+}} is the pseudoinverse of the positive part of
#' \eqn{H} in the eigenbasis; with \code{sparse = TRUE} the decrease comes
#' from the factorization of \eqn{H} where \eqn{H} is positive definite and
#' from the dense eigendecomposition elsewhere. \code{trustL1} applies the
#' same condition to the model of the coordinates not held at a kink. Without
#' the condition the run goes on, within \code{iterlim}; \code{mtol = 0}
#' drops the condition. A run that goes on this way and then ends on the
#' radius without an accepted step that changes the value by \code{ftol} or
#' more ends on the stop held back, with \code{converged = TRUE}: the decrease
#' the model predicts is below the resolution of the objective.
#'
#' \code{stopReason} in the result names the test that fired:
#' \code{"gradient"}, \code{"fvalue"}, \code{"preddiff"}, \code{"step"},
#' \code{"stagnation"}, \code{"radius"}, \code{"objfun"} (three failed
#' evaluations in a row, ten with multiple shooting) or \code{"iterlim"}.
#' Stagnation reports \code{converged = TRUE} where the full Newton step
#' predicts a decrease of at most \code{mtol} or a stop was held back as
#' above, and \code{converged = FALSE} otherwise.
#'
#' Only \code{gtol} tests first-order optimality. While \code{switches}
#' remain, see \code{phases}, an \code{"fvalue"}, \code{"preddiff"},
#' \code{"step"} or \code{"stagnation"} stop hands over to the next phase
#' instead of ending the run. In a quasi-Newton phase \code{"fvalue"},
#' \code{"preddiff"} and \code{"step"} do not end the run. A phase that
#' reaches the \code{iterlim} of that phase ends like a soft stop, with
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
#' alone or in a sum, is optimized by multiple shooting with the same controls.
#' Not available then: \code{hessianMethod = "exact"} and the structured
#' methods, more than one phase,
#' \code{switches}, \code{qnControl$hessianInit = "exact"},
#' \code{qnControl$hessianReseed}, \code{qnControl$qnMemory},
#' \code{stepControl$nonmonotone},
#' \code{boundary = "clip"}, \code{minimize = FALSE} and \code{traceFile}.
#' The run reports \code{converged = TRUE} only when the scaled gaps between
#' the segments are below \code{tolControl$ctol}. The stops follow section
#' Convergence, with the full step of the condensed model as the Newton step;
#' stagnation also needs the gaps flat within \code{ctol}. The method is
#' described in
#' the methods article
#' [Multiple shooting](https://dmodverse.github.io/dMod2/articles/methods-multiple-shooting.html).
#'
#' @param objfun Function whose first argument is a numeric vector of
#'   parameters, returning a list with components \code{value},
#'   \code{gradient} and \code{hessian}. It is called with the logical
#'   arguments \code{deriv} and \code{hessian}, and with \code{deriv2 = TRUE}
#'   when an exact Hessian is requested, so it must accept them or
#'   \code{...}. \code{trustL1} calls it with the parameters and \code{...}
#'   only.
#' @param parinit Named numeric starting vector. Must be finite. Values
#'   outside \code{[parlower, parupper]} are clipped with a warning; with
#'   \code{boundary = "reflective"} they are moved just inside the box.
#' @param rinit Initial trust-region radius. Default \code{0.1}.
#' @param rmax Maximum trust-region radius. Default \code{10}.
#' @param iterlim Maximum number of trust-region iterations over all phases.
#'   Default \code{1000}.
#' @param hessianMethod Source of the model Hessian of a run without
#'   \code{phases}: \code{"gn"} the Gauss-Newton Hessian of the
#'   objective, \code{"bfgs"} or \code{"sr1"} a quasi-Newton update seeded as
#'   \code{qnControl$hessianInit} says, \code{"exact"} the true Hessian at
#'   every iterate, which needs an objective that takes \code{deriv2}.
#'   The structured methods take the Gauss-Newton Hessian at every iterate and
#'   add a correction learned from the steps, starting at zero:
#'   \code{"sbfgs"} by a BFGS update, \code{"ssr1"} by an SR1 update, which may
#'   be indefinite; \code{"gnsbfgs"} and \code{"gnssr1"} switch at every
#'   accepted step between the Gauss-Newton Hessian alone and the structured
#'   BFGS or SR1 model. They count as quasi-Newton phases for the stops and
#'   \code{qnEval}. The correction is learned from the secant of Dennis,
#'   Martínez and Tapia on the data rows that an objective of class
#'   \code{objfn} returns with \code{hessian = "structured"}, see [objlist()];
#'   without them from the gradient and the Gauss-Newton Hessian. The formulas
#'   are in the methods article.
#'   All but \code{"gn"} require \code{boundary = "reflective"}.
#'   Also selects the control defaults below. With \code{phases}, each phase
#'   names a method, and this argument is an error. Not given: the hybrid of
#'   section Default method.
#' @param phases List of phases, each a list of settings, run one after the
#'   other and from the first again after the last. The run hands over to the
#'   next phase at a value, model, step or stagnation stop, when a phase
#'   uses up the \code{iterlim} of that phase and when the radius stayed the same
#'   for \code{radiusStall} iterations. A handover restarts the trust region at
#'   \code{rinit}, unless \code{keepRadius}, and discards the quasi-Newton
#'   pairs. Default \code{NULL}: one phase from \code{hessianMethod} and
#'   \code{qnControl}, evaluated with the \code{sweep} in \code{...}, or the
#'   phases of section Default method. Names of the list label the phases in the
#'   result, the Hessian source does otherwise. Not available with multiple
#'   shooting and \code{boundary = "clip"}.
#'   \describe{
#'     \item{\code{hessianMethod}}{Hessian source of the phase, as
#'       \code{hessianMethod}. Required.}
#'     \item{\code{iterlim}}{Maximum number of iterations of one run of the
#'       phase. Default \code{NULL}: only the total \code{iterlim}.}
#'     \item{\code{radiusStall}}{Number of consecutive iterations that leave
#'       the trust radius unchanged after which the phase hands over, the
#'       switch of the hybrid method of Fides (Fröhlich and Sorger 2022), which
#'       takes 20 for a Gauss-Newton phase followed by BFGS. The count starts
#'       anew with every change of the radius and at a handover; without a
#'       handover left it ends nothing. Default \code{NULL}: no such handover.}
#'     \item{\code{sweep}}{\code{"forward"} or \code{"reverse"}, the
#'       \code{sweep} the objective is evaluated with in the phase; a
#'       \code{sweep} in \code{...} is an error with \code{phases}. Default
#'       \code{NULL}: for \code{"bfgs"} and \code{"sr1"} reverse if the
#'       objective can be evaluated backwards at \code{parinit}, as with
#'       \code{derivMode = c("forward", "reverse")} throughout the chain,
#'       forward otherwise. A Gauss-Newton Hessian is always computed forward,
#'       also in the structured methods.
#'       A handover between phases of different sweep evaluates the current
#'       iterate once more. The result reports the sweeps as
#'       \code{phaseSweep}.}
#'     \item{\code{keepRadius}}{\code{TRUE} starts the phase at the trust
#'       radius of the handover instead of \code{rinit}. Default
#'       \code{FALSE}.}
#'     \item{\code{onStop}}{What a value, model, step or stagnation stop does
#'       in the phase: \code{"handover"} (default) hands over while
#'       \code{switches} remain, \code{"end"} ends the run. A phase that ends
#'       on stops hands over at \code{iterlim} and \code{radiusStall} only.}
#'     \item{\code{qnControl}}{Quasi-Newton settings of the phase, merged
#'       over \code{qnControl}; see there. \code{hessianInit} defaults to
#'       \code{"gn"} in every phase but the first: a quasi-Newton phase starts
#'       from the Gauss-Newton Hessian at the iterate of the handover.
#'       \code{"tracked"}, for a \code{"bfgs"} phase after the first with
#'       \code{qnMemory = 0}, starts it from a BFGS matrix kept over the whole
#'       run instead: the Hessian at \code{parinit} of the first phase, made
#'       positive definite as a BFGS seed is, updated with the accepted steps of
#'       the phases before, under the \code{qnCautious} test of the phase, and
#'       updated on in the phase itself. Not available on the sparse path.}
#'   }
#' @param switches Number of handovers allowed. Default \code{NULL}: one
#'   fewer than there are phases, so each runs once. A higher value cycles
#'   through the phases again, \code{0} runs the first phase only. The result
#'   reports the count as \code{nSwitch}.
#' @param parupper,parlower Upper and lower bounds. Unnamed: the first element
#'   applies to all parameters. Named: entries are matched by name, missing
#'   parameters are unbounded. Default \code{NULL}, unbounded.
#' @param boundary Box-bound handling, \code{"reflective"} (default) or
#'   \code{"clip"}. See section Box bounds.
#' @param tolControl List of convergence thresholds, merged into the defaults.
#'   \describe{
#'     \item{\code{ftol}}{Change in objective value. Default \code{1e-6}.}
#'     \item{\code{mtol}}{Predicted model decrease, of the step taken and,
#'       at a convergence stop, of the full Newton step (section Convergence).
#'       \code{0} switches off both the model test and the condition on the
#'       full Newton step. Default \code{1e-6}.}
#'     \item{\code{gtol}}{First-order optimality measure. Default \code{1e-6}.}
#'     \item{\code{xtol}}{Step norm. Default \code{0}, disabled.}
#'     \item{\code{rmin}}{Lower limit on the trust radius; falling below it
#'       stops the run with \code{converged = FALSE}. Default \code{0},
#'       disabled.}
#'     \item{\code{ctol}}{Multiple shooting only: scaled gap size below which
#'       the trajectory counts as continuous; no other test stops the run
#'       before. Default \code{1e-6}.}
#'   }
#'   \code{trustL1} reads \code{ftol}, \code{mtol}, \code{gtol}, \code{xtol}
#'   and \code{rmin}.
#' @param qnControl List of quasi-Newton settings, merged into the defaults and
#'   read in a quasi-Newton or structured phase; \code{hessianInit} also in
#'   a first phase of source \code{"gn"}. A structured phase takes
#'   \code{hessianInit = "gn"} and \code{qnMemory = 0} only.
#'   \describe{
#'     \item{\code{hessianInit}}{Seed of the approximation: \code{"gn"}
#'       (default) the Gauss-Newton Hessian at \code{parinit};
#'       \code{"identity"} the identity matrix, after which the objective is
#'       asked for gradients only until a reseed; \code{"exact"} the true
#'       Hessian at \code{parinit}, then gradients only. \code{"exact"} needs
#'       \code{deriv2 = TRUE} in a call, see \code{derivMode} in [odemodel()]. In a first phase
#'       of source \code{"gn"}, \code{"exact"} makes the first step use the
#'       exact Hessian at \code{parinit}.}
#'     \item{\code{hessianReseed}}{\code{"never"} (default): a stalled
#'       quasi-Newton phase, or one whose trust radius collapses, stops.
#'       \code{"stall"}: it continues from a fresh Hessian at the current
#'       iterate, the exact one with \code{hessianInit = "exact"} and the
#'       Gauss-Newton one otherwise, provided the iterate moved since the last
#'       reseed. The result reports the count as \code{nReseed}.}
#'     \item{\code{qnMemory}}{Number of \code{(s, y)} pairs kept. Default
#'       \code{0}: every update accumulates onto the seed. A positive value
#'       rebuilds the approximation in each iteration from the last
#'       \code{qnMemory} pairs, without the seed. On the sparse path, see
#'       \code{sparse}, BFGS keeps the last \code{qnMemory} pairs, 10 where it
#'       is \code{0}, on top of the seed.}
#'     \item{\code{qnCautious}}{A BFGS pair is used only when
#'       \eqn{s^\top y > c\,|s|\,|y|} with \eqn{c} = \code{qnCautious}. Default
#'       \code{1e-8}, \code{0} disables the test; also read by \code{"sbfgs"},
#'       the SR1 methods are exempt. The result reports the rejected pairs as
#'       \code{qnSkipped}.}
#'     \item{\code{qnRejected}}{Whether an SR1 method also updates from a
#'       rejected trial point. Default \code{TRUE} for a phase whose source is
#'       \code{"sr1"}, \code{"ssr1"} or \code{"gnssr1"}, \code{FALSE}
#'       otherwise.}
#'     \item{\code{qnScaling}}{\code{"none"} (default) or \code{"residual"}:
#'       \code{"sbfgs"} and \code{"ssr1"} scale the correction with the
#'       residual norm, the norm of the derivatives of the value in the
#'       predictions and sigmas of the data, or the square root of the value
#'       for an objective without data rows, so it vanishes with the residuals
#'       (Huschens). Read by these two methods only.}
#'     \item{\code{qnHybridTol}}{Threshold of the switch of \code{"gnsbfgs"}
#'       and \code{"gnssr1"}: below it the model is the Gauss-Newton Hessian.
#'       Default \code{1e-6}.}
#'     \item{\code{qnHybridShift}}{\code{TRUE} adds the residual norm times
#'       the identity to the Gauss-Newton Hessian on the Gauss-Newton branch
#'       of \code{"gnsbfgs"} and \code{"gnssr1"}, as Zhou and Chen; default
#'       \code{FALSE}.}
#'   }
#' @param stepControl List of step and acceptance settings, merged into the
#'   defaults.
#'   \describe{
#'     \item{\code{thetaMax}}{Smallest fraction of the distance to a bound a
#'       step may use. The fraction is \code{max(thetaMax, 1 - m)}, below 1,
#'       with \code{m} the first-order optimality measure. Default
#'       \code{0.99995}. \code{boundary = "reflective"} only.}
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
#'     \item{\code{regularize}}{Multiple shooting only: directions of the
#'       condensed problem with a singular value below \code{regularize} times
#'       the largest are held fixed in a step. Default \code{0}, off. With
#'       \code{acceptance = "natural"} directions below \code{1e-7} are always
#'       held fixed.}
#'     \item{\code{soc}}{Multiple shooting only: if \code{TRUE}, a rejected
#'       step whose gaps did not close as predicted is retried once with a
#'       second-order correction. Default \code{TRUE}.}
#'     \item{\code{anneal}}{Multiple shooting only: if \code{TRUE}, the run
#'       first minimizes the objective plus a rising penalty on the gaps, then
#'       enforces continuity. Needs forward sensitivities and is skipped
#'       without them. Default \code{TRUE}.}
#'   }
#'   \code{trustL1} takes \code{thetaMax} only.
#' @param minimize If \code{TRUE} (default) minimize, if \code{FALSE}
#'   maximize.
#' @param blather If \code{TRUE}, the result also contains the per-iteration
#'   trace listed under Value. Default \code{FALSE}.
#' @param printIter If \code{TRUE}, print the number of the evaluation and the
#'   objective value at each evaluation. Default \code{FALSE}.
#' @param traceFile Path of a CSV file written with one row per evaluation:
#'   the number of the evaluation (column \code{Iteration}), objective value
#'   and parameters. Default \code{NULL}, no file.
#' @param sparse \code{"auto"} (default), \code{TRUE} or \code{FALSE}.
#'   \code{TRUE} calls \code{objfun} with \code{hessian = "sparse"} and solves
#'   every subproblem on sparse Cholesky factors of the returned \code{Matrix},
#'   a low-rank part \code{lowrank} (see [objlist()]) included; on the Hessian
#'   of the objective the steps equal those of the dense path to rounding. It
#'   needs \code{hessianMethod}
#'   \code{"gn"} or \code{"exact"}, or \code{"bfgs"} with
#'   \code{qnControl$hessianInit = "gn"}, in every phase, and
#'   \code{boundary = "reflective"}. A BFGS phase there is limited-memory BFGS:
#'   the Hessian of the objective at the start of the phase plus the compact
#'   correction of the last \code{qnMemory} pairs (Byrd, Nocedal and Schnabel
#'   1994), a low-rank part of rank twice the number of pairs.
#'   \code{"auto"} takes the sparse path where every phase takes the Hessian of
#'   the objective, outside multiple shooting, for at least \code{parameters} entries of \code{parinit} and a
#'   step cheaper than the dense one at \code{parinit}; the settings come from
#'   \code{getOption("dMod2.sparseAuto")}, see [dMod2-package].
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
#'     \item{\code{qnSkipped}}{Pairs rejected by the cautious test, and in a
#'       structured phase pairs that left the correction unchanged.}
#'     \item{\code{nSwitch}}{Number of handovers between phases.}
#'     \item{\code{phaseIterations}}{Named integer, the iterations of each
#'       phase run, in run order and named after the phase. The sum is
#'       \code{iterations} without the iteration that ends the run on failed
#'       evaluations and without the iterations that reseed instead of
#'       taking a step.}
#'     \item{\code{phaseSweep}}{Named character, the sweep of each element
#'       of \code{phases}, present when \code{phases} is set and in the
#'       default hybrid.}
#'     \item{\code{nReseed}}{Number of reseeds by
#'       \code{qnControl$hessianReseed = "stall"}, which allows each run of a
#'       phase the same number.}
#'     \item{\code{hybridGN}}{Accepted steps of a \code{"gnsbfgs"} or
#'       \code{"gnssr1"} phase after which the model is the Gauss-Newton
#'       Hessian.}
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
#'   \code{preddiff}, \code{stepnorm} and, with
#'   \code{boundary = "reflective"}, \code{hessianSource}, which is
#'   \code{"gn"} for an iteration of a hybrid method on the Gauss-Newton
#'   Hessian.
#'
#'   \code{boundary = "clip"} returns only \code{argument}, \code{value},
#'   \code{gradient}, \code{hessian}, \code{iterations}, \code{converged},
#'   \code{atBound} and \code{stopReason}, plus the trace.
#'
#'   On the sparse path, see \code{sparse}, \code{hessian} is a sparse
#'   \code{Matrix}, and \code{lowrank} where the objective returns one or a
#'   BFGS phase adds the correction: the Hessian is then \code{hessian} plus
#'   that part.
#'
#'   \code{trustL1} returns a list without class with \code{argument},
#'   \code{value}, \code{gradient}, \code{hessian}, \code{iterations},
#'   \code{converged}, \code{atBound} and \code{stopReason}, plus the trace
#'   without \code{hessianSource}. \code{value} and \code{gradient} include the
#'   penalty, whose gradient counts zero on a kink; \code{hessian} is that of
#'   \code{objfun}.
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
#'   the guide chapter
#'   [Optimization](https://dmodverse.github.io/dMod2/articles/guide-optimization.html), the methods
#'   article [Trust-region optimization](https://dmodverse.github.io/dMod2/articles/methods-trust.html).
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
#' # Gauss-Newton, then BFGS with a separate budget, then Gauss-Newton again
#' fit <- trust(rosenbrock, c(a = -1.2, b = 1), iterlim = 200,
#'              phases = list(list(hessianMethod = "gn"),
#'                            list(hessianMethod = "bfgs", iterlim = 50)),
#'              switches = 2)
#' fit$phaseIterations
#' fit$evalBySource
#'
#' @references Byrd, R. H., Nocedal, J. and Schnabel, R. B. (1994).
#'   Representations of quasi-Newton matrices and their use in limited memory
#'   methods. \emph{Mathematical Programming} 63, 129-156.
#'
#'   Coleman, T. F. and Li, Y. (1996). An interior trust region
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
#'   Dennis, J. E., Martínez, H. J. and Tapia, R. A. (1989). Convergence theory
#'   for the structured BFGS secant method with an application to nonlinear
#'   least squares. \emph{Journal of Optimization Theory and Applications}
#'   61(2), 161-178.
#'
#'   Huschens, J. (1994). On the use of product structure in secant methods
#'   for nonlinear least squares problems. \emph{SIAM Journal on Optimization}
#'   4(1), 108-129.
#'
#'   Zhou, W. and Chen, X. (2010). Global convergence of a new hybrid
#'   Gauss-Newton structured BFGS method for nonlinear least squares problems.
#'   \emph{SIAM Journal on Optimization} 20(5), 2422-2441.
#'
#'   Peifer, M. and Timmer, J. (2007). Parameter estimation in ordinary
#'   differential equations for biochemical processes using the method of
#'   multiple shooting. \emph{IET Systems Biology} 1(2), 78-88.
#'
#' @export
trust <- function(objfun, parinit, rinit = 0.1, rmax = 10,
                  iterlim   = 1000L,
                  hessianMethod   = c("gn", "bfgs", "sr1", "exact", "sbfgs", "ssr1",
                                      "gnsbfgs", "gnssr1"),
                  phases    = NULL,
                  switches  = NULL,
                  parupper  = NULL,
                  parlower  = NULL,
                  boundary  = c("reflective", "clip"),
                  tolControl  = NULL,
                  qnControl   = NULL,
                  stepControl = NULL,
                  minimize  = TRUE,
                  blather   = FALSE,
                  printIter = FALSE,
                  traceFile = NULL,
                  sparse    = "auto",
                  ...) {
  dots <- list(...)
  .checkSparse(sparse, "trust")
  if (!is.null(phases) && !missing(hessianMethod))
    stop("trust: with `phases`, every phase names the hessianMethod.",
         call. = FALSE)
  if (!is.null(phases) && !is.null(dots$sweep))
    stop("trust: with `phases`, every phase sets the sweep.", call. = FALSE)
  hybrid <- is.null(phases) && missing(hessianMethod) && is.null(switches) && is.null(dots$sweep)
  hessianMethod <- if (is.null(phases)) match.arg(hessianMethod) else
    .phaseMethod(phases[[1L]]$hessianMethod, 1L)

  if ("boundary" %in% names(stepControl))
    stop("trust: boundary is an argument of trust(), not a stepControl entry.",
         call. = FALSE)
  ctl  <- .trustDefaults(hessianMethod)
  tol  <- .mergeControl(ctl$tolControl,  tolControl,  "tolControl")
  qn   <- .mergeControl(ctl$qnControl,   qnControl,   "qnControl")
  step <- .mergeControl(ctl$stepControl, stepControl, "stepControl")

  boundary <- match.arg(boundary)
  step$boundary <- boundary
  ph <- .trustPhases(phases, switches, hessianMethod, qnControl)

  # the sparse model takes the Hessian of the objective in every phase; sparse = TRUE also
  # limited-memory BFGS on it
  fromObjective <- all(vapply(ph$phases, function(p) p$hessianMethod %in% c("gn", "exact"), TRUE))
  lbfgs <- all(vapply(ph$phases, function(p) p$hessianMethod %in% c("gn", "exact") ||
    (p$hessianMethod == "bfgs" && identical(p$qnControl$hessianInit, "gn")), TRUE))
  first <- NULL
  if (identical(sparse, "auto")) {
    choice <- if (fromObjective && boundary == "reflective" && !length(.shootingTerms(objfun)))
      .sparseChoice(objfun, parinit, dots) else list(sparse = FALSE)
    sparse <- choice$sparse
    first <- choice$first
  }
  # the default is the hybrid of Fides where the run takes the dense path and phases
  if (hybrid && !isTRUE(sparse) && boundary == "reflective" && !length(.shootingTerms(objfun))) {
    phases <- .hybridPhases
    ph <- .trustPhases(phases, NULL, hessianMethod, qnControl)
  }
  if (isTRUE(sparse) && (!lbfgs || boundary != "reflective"))
    stop("trust(sparse = TRUE) takes the Hessian of the objective, hessianMethod \"gn\" or \"exact\", ",
         "or limited-memory BFGS on it, \"bfgs\" with qnControl$hessianInit = \"gn\", ",
         "with boundary = \"reflective\".", call. = FALSE)

  # The objective decides: a normL2 with a multipleShootingControl, alone or
  # in a sum, is optimized by multiple shooting, with the same controls.
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
  # own arguments here, where the formals of the objective are visible. The exact request is checked
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
  # phase passes a separate sweep in place of one in `...`, and a request for a
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

  pnames <- names(parinit)
  # a run with a structured phase asks an objfn for the data rows with every Gauss-Newton Hessian
  structured <- !sparse && inherits(objfun, "objfn") &&
    any(vapply(ph$phases, `[[`, "", "hessianMethod") %in% .structuredMethods)
  fn <- function(x, want = 2L, phase = 0L) {
    want <- as.integer(want)
    # the evaluation of the sparse choice at parinit is the first one of the kernel
    if (!is.null(first) && want == 2L && identical(as.vector(x), as.vector(parinit))) {
      out <- first
      first <<- NULL
      return(.kernelHessian(out, pnames))
    }
    args <- list(x, deriv = want >= 1L, hessian = if (want >= 2L) (if (sparse) "sparse" else TRUE) else FALSE)
    if (want == 2L && structured) args$hessian <- "structured"
    if (want >= 3L) args$deriv2 <- TRUE
    d <- phaseDots[[phase + 1L]]
    if (!is.null(phases) && want == 2L) d$sweep <- NULL
    out <- .tryObjective(do.call(objfun, c(args, d)))
    if (is.list(out) && !is.null(out$structure)) out$structure <- .structureOver(out$structure, pnames)
    if (!sparse || want < 2L || inherits(out, "try-error")) return(out)
    if (!is.null(out$hessian) && !.isSparse(out$hessian)) {
      dense <- .denseTerms(objfun, x, d)
      return(structure(paste0("trust(sparse = TRUE): the objective returns a dense Hessian",
                              if (length(dense) && inherits(objfun, "objfn"))
                                paste0(", from term ", paste(dense, collapse = ", "), " of the sum"),
                              "; use sparse = FALSE."), class = "try-error"))
    }
    .kernelHessian(out, pnames)
  }
  qfield <- function(f) lapply(ph$phases, function(p) p$qnControl[[f]])
  fit <- trust_impl(fn, parinit, rinit, rmax, NULL, as.integer(iterlim),
             tol$ftol, tol$mtol, tol$gtol, tol$xtol, tol$rmin, step$thetaMax,
             boundary,
             vapply(ph$phases, `[[`, "", "hessianMethod"),
             vapply(ph$phases, `[[`, 1L, "iterlim"),
             vapply(ph$phases, `[[`, 1L, "radiusStall"),
             sweeps == "reverse",
             vapply(ph$phases, `[[`, TRUE, "keepRadius"),
             vapply(ph$phases, function(p) p$onStop == "end", TRUE),
             vapply(qfield("hessianInit"), function(v)
               match.arg(v, c("gn", "identity", "exact", "tracked")), ""),
             vapply(qfield("hessianReseed"), function(v)
               match.arg(v %||% "never", c("never", "stall")), ""),
             vapply(qfield("qnMemory"), as.integer, 1L),
             vapply(qfield("qnCautious"), as.numeric, 1),
             vapply(qfield("qnRejected"), isTRUE, TRUE),
             vapply(qfield("qnScaling"), function(v)
               match.arg(v, c("none", "residual")) == "residual", TRUE),
             vapply(qfield("qnHybridTol"), as.numeric, 1),
             vapply(qfield("qnHybridShift"), isTRUE, TRUE),
             ph$switches,
             step$nonmonotone, isTRUE(sparse), minimize, blather,
             parupper, parlower, printIter, traceFile)
  if (!is.null(fit$phaseIterations)) {
    n <- length(fit$phaseIterations)
    names(fit$phaseIterations) <- ph$labels[(seq_len(n) - 1L) %% length(ph$labels) + 1L]
  }
  if (!is.null(phases)) fit$phaseSweep <- stats::setNames(sweeps, ph$labels)
  .trustFit(fit)
}

# An objlist with the sparse Hessian over `pnames` as the kernel reads it, `kernel`: the upper
# triangle in compressed columns and the low-rank part, list(p, i, x, U, d). A dense or missing
# Hessian leaves it without, which the kernel takes as no Hessian.
.kernelHessian <- function(out, pnames) {
  if (!is.list(out) || !.isSparse(out$hessian)) return(out)
  out$hessian <- .sparseOver(out$hessian, pnames)
  T <- methods::as(out$hessian, "TsparseMatrix")
  up <- T@i <= T@j
  U <- Matrix::sparseMatrix(i = T@i[up] + 1L, j = T@j[up] + 1L, x = T@x[up], dims = rep(length(pnames), 2L))
  lr <- .lowrankOver(out$lowrank, pnames)
  out$lowrank <- lr
  out$kernel <- list(p = U@p, i = U@i, x = U@x, U = if (is.null(lr)) matrix(0, length(pnames), 0L) else lr$U,
                     d = if (is.null(lr)) numeric(0) else as.numeric(lr$d))
  out
}


.phaseMethod <- function(method, k) {
  if (is.null(method))
    stop("trust: phases[[", k, "]] names no hessianMethod.", call. = FALSE)
  match.arg(method, .hessianMethods)
}

.hessianMethods <- c("gn", "bfgs", "sr1", "exact", "sbfgs", "ssr1", "gnsbfgs", "gnssr1")
# The default run of trust(): the hybrid of Fides (Fröhlich and Sorger 2022).
.hybridPhases <- list(
  gn = list(hessianMethod = "gn", radiusStall = 20L, onStop = "end"),
  bfgs = list(hessianMethod = "bfgs", sweep = "forward", keepRadius = TRUE,
              qnControl = list(hessianInit = "tracked")))
# Methods that keep the Gauss-Newton Hessian of the objective and update a correction.
.structuredMethods <- c("sbfgs", "ssr1", "gnsbfgs", "gnssr1")

# The phases of a run, each with the Hessian source, budget (NA: no phase
# budget), sweep and quasi-Newton settings merged over the defaults of the source
# and the top-level qnControl. Without `phases`, one phase from the top level.
.trustPhases <- function(phases, switches, hessianMethod, qnControl) {
  if (is.null(phases))
    phases <- list(list(hessianMethod = hessianMethod))
  else if (!is.list(phases) || !length(phases) ||
           !all(vapply(phases, is.list, TRUE)))
    stop("trust: `phases` must be a non-empty list of lists.", call. = FALSE)

  defaults <- list(hessianMethod = NULL, iterlim = NULL, radiusStall = NULL,
                   sweep = NULL, keepRadius = FALSE, onStop = "handover", qnControl = NULL)
  out <- lapply(seq_along(phases), function(k) {
    p <- .mergeControl(defaults, phases[[k]], paste0("phases[[", k, "]]"))
    p$hessianMethod <- .phaseMethod(p$hessianMethod, k)
    lim <- p$iterlim %||% Inf
    if (length(lim) != 1L || is.na(lim) || lim < 1)
      stop("trust: phases[[", k, "]]$iterlim must be a positive number.",
           call. = FALSE)
    p$iterlim <- if (is.finite(lim)) as.integer(lim) else NA_integer_
    stall <- p$radiusStall %||% Inf
    if (length(stall) != 1L || is.na(stall) || stall < 1)
      stop("trust: phases[[", k, "]]$radiusStall must be a positive number.",
           call. = FALSE)
    p$radiusStall <- if (is.finite(stall)) as.integer(stall) else NA_integer_
    if (!isTRUE(p$keepRadius) && !isFALSE(p$keepRadius))
      stop("trust: phases[[", k, "]]$keepRadius must be TRUE or FALSE.", call. = FALSE)
    p$onStop <- match.arg(p$onStop, c("handover", "end"))
    q <- .mergeControl(.trustDefaults(p$hessianMethod)$qnControl, qnControl,
                       "qnControl")
    # A later phase starts from the Gauss-Newton Hessian at the handover.
    if (k > 1L) q$hessianInit <- "gn"
    p$qnControl <- .mergeControl(q, p$qnControl,
                                 paste0("phases[[", k, "]]$qnControl"))
    if (p$hessianMethod %in% .structuredMethods &&
        (!identical(p$qnControl$hessianInit, "gn") || !identical(as.integer(p$qnControl$qnMemory), 0L)))
      stop("trust: hessianMethod \"", p$hessianMethod, "\" starts from the Gauss-Newton ",
           "Hessian and keeps no pairs; it takes qnControl$hessianInit = \"gn\" and ",
           "qnMemory = 0.", call. = FALSE)
    if (identical(p$qnControl$hessianInit, "tracked") &&
        (k == 1L || p$hessianMethod != "bfgs" || !identical(as.integer(p$qnControl$qnMemory), 0L)))
      stop("trust: qnControl$hessianInit = \"tracked\" takes a \"bfgs\" phase after the first, ",
           "with qnMemory = 0.", call. = FALSE)
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
  if (method %in% c("gn", .structuredMethods)) {
    if (identical(sweep, "reverse"))
      stop("trust: a phase on the Gauss-Newton Hessian needs the forward sensitivities; ",
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
             qnMemory = 0L, qnCautious = 1e-8, qnRejected = FALSE,
             qnScaling = "none", qnHybridTol = 1e-6, qnHybridShift = FALSE)
  # Nocedal and Wright, sec. 6.2: SR1 also updates from a rejected step.
  if (hessianMethod %in% c("sr1", "ssr1", "gnssr1")) qn$qnRejected <- TRUE
  list(
    tolControl  = list(ftol = 1e-6, mtol = 1e-6, gtol = 1e-6,
                       xtol = 0, rmin = 0, ctol = 1e-6),
    qnControl   = qn,
    stepControl = list(thetaMax = 0.99995,
                       nonmonotone = 0, acceptance = "filter", soc = TRUE,
                       anneal = TRUE, twoPhase = FALSE, regularize = 0,
                       restore = FALSE))
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
#' @param mu Named numeric vector of reference values for the L1-penalized
#'   parameters. Names must be a subset of \code{names(parinit)}; only the
#'   named parameters receive a penalty. With \code{boundary = "reflective"}
#'   each must lie strictly inside \code{[parlower, parupper]}. Defaults to a
#'   zero vector covering all of
#'   \code{parinit}.
#' @param one.sided Logical. If \code{TRUE}, the penalty is one-sided,
#'   \code{lambda * max(0, mu - p)}: only values below \code{mu} are penalized.
#'   Otherwise it is the symmetric \code{lambda * |p - mu|}.
#' @param lambda Strength of the L1 penalty. Either a scalar (broadcast to
#'   all entries of \code{mu}) or a named numeric aligned with \code{mu}.
#' @param gate Character, names in \code{mu} whose penalty is the gate
#'   \code{lambda * (p - mu)} with \code{p >= mu}: the parameter can sit exactly
#'   on \code{mu} and leaves it only upward. Starts below \code{mu} are moved
#'   onto it. Needs \code{boundary = "reflective"}.
#' @param fuse List of fusion blocks. A block penalizes
#'   \code{lambda * sum_{i<j} w_ij |p_i - p_j|} over the parameters of the block, plus
#'   \code{lambda * w_ia |p_i - a|} toward an optional anchor \code{a}.
#'   Parameters with equal values move as one and split again when the data
#'   pull them apart. A block is a character vector of parameter names or a
#'   list with \code{pars}, and optionally \code{lambda} (default: the scalar
#'   \code{lambda}), \code{anchor} (default none) and \code{weights}, an
#'   \code{(m+1) x (m+1)} matrix whose last row and column refer to the
#'   anchor (default all 1). Fused parameters must not appear in \code{mu}; if
#'   \code{fuse} is given and \code{mu} is not, \code{mu} is empty. Needs
#'   \code{boundary = "reflective"}.
trustL1 <- function(objfun, parinit, mu = 0 * parinit, one.sided = FALSE, lambda = 1,
                    rinit = 0.1, rmax = 10,
                    iterlim   = 1000L,
                    parupper  = NULL,
                    parlower  = NULL,
                    boundary  = c("reflective", "clip"),
                    tolControl  = NULL,
                    stepControl = NULL,
                    minimize  = TRUE,
                    blather   = FALSE,
                    printIter = FALSE,
                    traceFile = NULL,
                    gate = NULL, fuse = NULL,
                    ...) {
  if ("hessianMethod" %in% names(list(...)))
    stop("trustL1: the model Hessian is the one of the objective; hessianMethod is an ",
         "argument of trust().", call. = FALSE)
  boundary <- match.arg(boundary)
  tolDefaults <- .trustDefaults("gn")$tolControl[c("ftol", "mtol", "gtol", "xtol", "rmin")]
  tol  <- .mergeControl(tolDefaults, tolControl, "tolControl")
  step <- .mergeControl(list(thetaMax = 0.99995), stepControl, "stepControl")

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

  # Maximizing f minus the penalty is minimizing -f plus the penalty; the kernel
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
                      tol$ftol, tol$mtol, tol$gtol, tol$xtol, tol$rmin, step$thetaMax,
                      boundary, TRUE, blather,
                      parupper, parlower, printIter, traceFile,
                      gate = gate, fuse = fuse)
  if (minimize) out else flip(out)
}

# Normalize trustL1 fuse blocks to list(idx, lambda, anchor, w) with idx
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
# against the optimizer's formals, so every optimizer argument is settable and a
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


## An R error reaches the kernel as a try-error with the message; thrown
## through Rcpp, the message would be lost.
.tryObjective <- function(expr)
  tryCatch(expr, error = function(e) structure(conditionMessage(e), class = "try-error"))
