# Entry points of mixed-effects fits: em() for one fit, msem() for a multistart. Both build the
# estimation design from the objective of the data and the random effects of the parameter
# transformation, and maximize the marginal likelihood over the random effects.


#' Mixed-effects fit by marginal maximum likelihood
#'
#' @description
#' Fits a mixed-effects model: the random effects that [insert()] with effect
#' formulas put into the parameter transformation are integrated out, and the
#' fixed parameters and the hyperparameters of the random effects are
#' estimated by maximizing the marginal likelihood.
#'
#' A group, a level of the grouping of a random term such as a patient or a
#' cell line, may own several conditions, which share the random effects of the
#' group.
#' Normal, spike-and-slab and Dirichlet-process random effects of one
#' grouping may be combined. Random terms over several groupings, nested or
#' crossed, are normal and integrated per connected component of the
#' conditions.
#'
#' Methods:
#' \describe{
#'   \item{`"focei"`}{Laplace approximation at the conditional modes with the
#'     exact derivative of the volume term. The default, and the method for
#'     spike-and-slab and Dirichlet-process random effects (at most 14
#'     groups) and for several groupings.}
#'   \item{`"quadrature"`}{adaptive sparse-grid quadrature around the modes,
#'     started from the FOCEI fit; section Quadrature.}
#'   \item{`"laplaceEM"`}{ECM on the Laplace approximation with the
#'     covariance-corrected update of `Omega`.}
#'   \item{`"saem"`}{stochastic-approximation EM, also for Dirichlet-process
#'     random effects; section SAEM.}
#' }
#'
#' @param obj Objective of the data, [normL2()] or [normT()] of the model
#'   `g * x * p`, without terms for the random effects.
#' @param random The random effects, [randomEffects()] of the parameter
#'   transformation or of `p`.
#' @param init Named numeric start with the parameters of `obj` that are
#'   neither random effects nor created by the effect formulas (those start
#'   at 0). It may also set those and the hyperparameters: `omega_<a>_<a>`,
#'   the log of the diagonal entry of the Cholesky factor of `Omega` (the log
#'   standard deviation of an uncorrelated random effect), and `omega_<a>_<b>`,
#'   the entry in row `a` and column `b` below the diagonal, with `a` and `b`
#'   the stems of the random effects without `eta_`, as in `omega_cl_cl`;
#'   `prob_*` and `tau_*` (spike and slab), `concentration_*` and `tau_*`
#'   (Dirichlet process), named `<kind>_<grouping>` and numbered
#'   `_1`, `_2`, ... where a grouping has several blocks. An `em` fit as
#'   `init` gives the start: the estimates, `fixed` values and random effects
#'   of that fit, without those it set to 0, and with sparse random effects
#'   the inner modes.
#' @param fixed Named numeric of parameters held fixed.
#' @param method `"focei"` (default), `"quadrature"`, `"laplaceEM"` or
#'   `"saem"`.
#' @param control List of options; other entries are an error.
#'   \describe{
#'     \item{`sd`}{start of the standard deviations of the random effects,
#'       default 0.3.}
#'     \item{`innerControl`}{the trust region of the modes; it stops where the
#'       Newton step is below `steptol`, default 1e-8.}
#'     \item{`trustControl`}{the outer trust region; it stops where the
#'       predicted change of -2 log L is below `mtol`, default 1e-5.}
#'     \item{`secondOrderCorrection`}{the derivative of the volume term,
#'       default `TRUE`.}
#'     \item{`derivMode`}{how FOCEI gets the second derivatives of the volume
#'       term: `"forward-reverse"` (default), sweeps for models built with
#'       `derivMode = c("forward", "reverse", "forward-reverse")`, or
#'       `"forward-forward"`, the second-order sensitivities of models built
#'       with `deriv2 = TRUE`. Without the mode asked for the fit takes the
#'       other, with a message.}
#'     \item{`boundary`}{the largest rise of -2 log L at which a variance is
#'       set to 0, default 0.1, `FALSE` for none. With one grouping and several
#'       blocks, a normal random effect that forms a block alone and whose
#'       standard deviation ends below 0.1 is dropped and the fit repeated from the
#'       estimates; where -2 log L rises by at most `boundary`, that fit is
#'       returned and names the hyperparameter in `boundary`.}
#'     \item{`quadrature`, `saem`, `sparseFit`}{lists of the sections
#'       Quadrature, SAEM and Sparse random effects.}
#'   }
#' @param cores Number of cores over the groups and conditions inside the
#'   fit; with 1 (default) no threads are started.
#' @param verbose Logical; print the progress of the iterative methods.
#' @param sparse `"auto"` (default), `TRUE` or `FALSE`: the outer trust region
#'   of FOCEI over normal random effects on the sparse Gauss-Newton Hessian.
#'   `"auto"` decides as [trust()] does.
#' @return An object of class `em` with, among others:
#'   \describe{
#'     \item{`argument`, `value`, `gradient`, `hessian`}{the estimates, -2 log
#'       L of the marginal likelihood with gradient and Hessian (a sparse
#'       `Matrix` on the sparse path, `NULL` with spike-and-slab or
#'       Dirichlet-process random effects).}
#'     \item{`information`}{with spike-and-slab random effects the
#'       empirical information that gives the standard errors: `G / (G - 1)`
#'       times the sum of the outer products of the centered scores of the `G`
#'       groups, plus half the Hessian of the conditions without random
#'       effects.}
#'     \item{`converged`, `iterations`}{}
#'     \item{`Omega`, `omega`}{the covariance of the normal random effects,
#'       their standard deviations and correlations with standard errors.}
#'     \item{`random`}{per group the modes with conditional standard errors;
#'       for spike-and-slab and Dirichlet-process random effects also the
#'       cluster, 0 for the population value, and the posterior probability
#'       `p0` of a deviation of 0; a column `data` where a sparse random effect
#'       reaches no data of the group; such an effect rests at 0 with the prior.}
#'     \item{`slab`}{the hyperparameters of the spike-and-slab and
#'       Dirichlet-process blocks.}
#'     \item{`valueSE`}{the standard error of `value` where it comes from
#'       importance sampling.}
#'     \item{`coupling`}{`"expectationConsistent"` where a group needed it,
#'       else `NULL`.}
#'     \item{`boundary`}{the hyperparameters set to 0, see `control`.}
#'     \item{`etaModes`, `etaSE`, `shrinkage`}{}
#'     \item{`certificate`}{with method `"quadrature"` per group the
#'       importance-sampling check at the estimate, section Quadrature.}
#'   }
#'   Methods: `print`, `summary`, `coef`, `vcov`, [confint][confint.em],
#'   [logLik][logLik.em], [residuals][residuals.em], [predict][predict.em],
#'   `plot`.
#' @section Quadrature:
#' Method `"quadrature"` alternates between placing the nodes around the
#' modes at the current parameters and one trust step over every parameter on
#' these frozen nodes. The nodes of the mu-referenced fixed effects are frozen
#' in `phi = X beta + eta`. Entries of `control$quadrature`:
#' \describe{
#'   \item{`tol`}{allowed change of -2 log L by the next level, default 0.1,
#'     split over the groups; each group takes the lowest level from the current
#'     stage of the group up within `tol / G` and never a lower one than before.
#'     `Inf` keeps the
#'     levels of the stages.}
#'   \item{`levels`}{the stages, default `K + 1` with `K` the number of random
#'     effects of a group; level `K + s` integrates polynomials of total
#'     degree `2 s + 1` exactly.}
#'   \item{`maxLevel`}{highest level, default `K + 8`.}
#'   \item{`tolDecrement`}{a stage ends when a Newton step on the rebuilt
#'     nodes would lower -2 log L by less than this, default 1e-4, or than the
#'     nodes resolve.}
#'   \item{`maxEcmPerStage`, `cm1Control`, `maxCm1Iter`}{rebuilds per stage
#'     (default 10), [trust()] settings and iterations of a step (default 30).}
#'   \item{`pruneTol`}{nodes whose weight lies this much in log below the
#'     largest are dropped, default `Inf`.}
#'   \item{`certificate`}{`FALSE`, or a list of `points` (default 256) and
#'     `seed` (default 1): importance sampling at the estimate from a t5
#'     proposal at the mode of each group. A group is flagged where the two
#'     -2 log L differ by more than three standard errors plus the error the
#'     level of the group allows, or the Pareto shape of the weights exceeds 0.7.}
#'   \item{`epsEcm`, `maxEcm`}{for `"laplaceEM"`: the extrapolated distance to
#'     the fixed point at which it stops, default 1e-5, and the iteration limit.}
#' }
#' @section SAEM:
#' Metropolis draws of the random effects by random-walk and independence
#' proposals at the conditional modes, Robbins-Monro averages of their
#' sufficient statistics, the mu-referenced fixed effects and `Omega` in
#' closed form, the other parameters by Newton steps on averaged quadratic
#' models. Entries of `control$saem`:
#' \describe{
#'   \item{`nBurnin`, `nEM`}{at most that many burn-in and averaging
#'     iterations, default 200 and 300.}
#'   \item{`nMcmc`, `nLaplace`}{random-walk and independence steps per
#'     iteration, default 3 each.}
#'   \item{`laplaceEvery`, `laplaceScale`}{iterations between refreshes of the
#'     modes, default 5, and the factor of the proposal covariance, default
#'     1.5.}
#'   \item{`window`, `mcTol`}{the burn-in ends when no parameter has a
#'     significant slope over the last `window` iterations (default 25), the
#'     averaging when the expected rise of -2 log L from the Monte Carlo error
#'     is below `mcTol` (default 0.05).}
#'   \item{`particles`}{with Dirichlet-process random effects the particles of
#'     the sequential importance sampling of -2 log L, default 500.}
#'   \item{`cm1Control`}{[trust()] settings of the Newton steps of the other
#'     parameters, default one iteration per SAEM iteration.}
#' }
#' With Dirichlet-process random effects a Gibbs sweep draws the clusters of
#' the groups (Neal 2000, algorithm 8), the averaging runs `nEM` iterations,
#' and the reported clusters are the most frequent partition of the late
#' draws.
#' @section Sparse random effects:
#' With spike-and-slab or Dirichlet-process random effects the
#' parameters and hyperparameters move in one trust region. The integral of a
#' group is exact where the data couple none of the factors and the sum over
#' the subsets of the spike-and-slab random effects of the group stays within a budget of
#' nodes; else it is the expectation consistent (EC) approximation, and a
#' trial point at which EC does not converge is a failed evaluation. Entries of `control$sparseFit`:
#' \describe{
#'   \item{`maxOuter`}{iterations of the trust region, default 200.}
#'   \item{`trustControl`}{the [trust()] settings; `ftol`, `mtol` and `gtol`
#'     of the `tolControl` there also end the search of the modes (else 1e-10).}
#' }
#'
#' @seealso [msem()], [randomEffects()], [insert()]; the guide chapters
#'   [Nonlinear mixed-effects models](https://dmodverse.github.io/dMod2/articles/guide-nlme.html) and
#'   [Sparse deviations](https://dmodverse.github.io/dMod2/articles/guide-sparse-deviations.html), the
#'   methods articles
#'   [Marginal likelihood of mixed-effects models](https://dmodverse.github.io/dMod2/articles/methods-nlme.html)
#'   and [Sparsity priors](https://dmodverse.github.io/dMod2/articles/methods-sparsity-priors.html).
#' @example inst/examples/em.R
#' @export
em <- function(obj, random, init = NULL, fixed = NULL,
               method = c("focei", "quadrature", "laplaceEM", "saem"),
               control = list(), cores = 1L, verbose = TRUE, sparse = "auto") {
  method <- match.arg(method)
  .checkEmControl(control)
  control$sparse <- .checkSparse(sparse, "em")
  cont <- .emContinue(obj, random, init, fixed)
  random <- cont$random; init <- cont$init; fixed <- cont$fixed
  design <- .emDesign(obj, random)
  start <- .emStart(design, init, fixed, sd = control$sd %||% 0.3)
  control$modes <- cont$modes
  .withCores(cores, .emFit(design, start, fixed, method, control, verbose))
}

## Refuses entries of `control` and of `control$sparseFit` that em() does not read.
.checkEmControl <- function(control) {
  known <- list(control = c("sd", "innerControl", "trustControl", "secondOrderCorrection", "derivMode", "quadrature",
                            "saem", "sparseFit", "boundary"),
                sparseFit = c("maxOuter", "trustControl"))
  for (k in names(known)) {
    x <- if (k == "control") control else control[[k]]
    bad <- setdiff(names(x), known[[k]])
    if (length(bad))
      stop("em: unknown ", if (k == "control") "entries of `control`: " else "entries of `control$sparseFit`: ",
           paste(bad, collapse = ", "), ".", call. = FALSE)
  }
  dm <- control$derivMode
  if (!is.null(dm) && !(length(dm) == 1L && dm %in% c("forward-reverse", "forward-forward")))
    stop("em: `control$derivMode` is \"forward-reverse\" or \"forward-forward\".", call. = FALSE)
  bd <- control$boundary
  if (!is.null(bd) && !isFALSE(bd) && !(is.numeric(bd) && length(bd) == 1L && is.finite(bd) && bd >= 0))
    stop("em: `control$boundary` is FALSE or one number of at least 0.", call. = FALSE)
}

## `expr` evaluated with `cores` cores over the conditions of an evaluation.
.withCores <- function(cores, expr) {
  old <- options(dMod2.cores = as.integer(cores))
  on.exit(options(old))
  expr
}


#' Multistart mixed-effects fit
#'
#' @description Runs [em()] from `fits` starts, the first at `init`, the others
#' perturbed by normal noise of standard deviation `sd` (on the log scale for
#' concentrations and slab standard deviations, on the logit scale
#' for slab probabilities), and returns the best fit.
#'
#' @inheritParams em
#' @param fits Number of starts. Default 20.
#' @param cores Number of fits run in parallel, or `c(fits = , conditions = )`
#'   as in [mstrust()]: fits in parallel and cores inside each fit.
#' @param sd Standard deviation of the perturbation of the starts: one
#'   number for every parameter, or a vector named by the parameters it
#'   perturbs, the others keeping their start.
#' @param retry Logical. If `TRUE` (default), a fit that fails with an error,
#'   for instance at a start where a steady state of [Pimpl()] has no stable
#'   root, is repeated from a new perturbed start, up to `nTries` attempts.
#' @param nTries Maximum number of attempts per start, including the first.
#'   Default `10`.
#' @param verbose Logical. If `TRUE`, print the table of the starts; the fits
#'   print no progress. Default `FALSE`.
#' @param keepFits Logical. If `TRUE`, the result has `fits`, the fit of every
#'   start, or the error message of a failed one, without the model. Default `FALSE`.
#' @param which Positions of the starts to fit, among the `fits` drawn. `NULL`
#'   (default) fits all. Tasks of [distributedComputing()] that set the same
#'   seed and fit disjoint parts fit the starts of one call, and `get()` joins
#'   their results into one.
#' @example inst/examples/msem.R
#' @return The fit with the smallest `value`, of class `em`, with
#'   `multistart`, a data frame with a row per start:
#'   \describe{
#'     \item{`start`}{the position of the start among the `fits` drawn.}
#'     \item{`value`, `converged`, `iterations`}{of the fit.}
#'     \item{`gradient`}{the largest absolute entry of the gradient at the end
#'       of the fit.}
#'     \item{`tries`}{the attempts it took.}
#'     \item{`stop`}{the test that ended the trust region, `"error"` where
#'       every attempt failed.}
#'     \item{`failures`}{why evaluations in the trust region failed, and how
#'       often.}
#'     \item{`error`, `boundary`}{the error of a failed start, the variances
#'       the fit set to 0.}
#'   }
#'   [plotValues()] draws it as a waterfall.
#' @seealso [em()]
#' @export
msem <- function(obj, random, init = NULL, fixed = NULL,
                 method = c("focei", "quadrature", "laplaceEM", "saem"),
                 control = list(), fits = 20L, cores = 1L, sd = 0.5, verbose = FALSE, sparse = "auto",
                 retry = TRUE, nTries = 10L, keepFits = FALSE, which = NULL) {
  method <- match.arg(method)
  .checkEmControl(control)
  which <- .msemWhich(which, fits)
  control$sparse <- .checkSparse(sparse, "msem")
  cont <- .emContinue(obj, random, init, fixed)
  random <- cont$random; init <- cont$init; fixed <- cont$fixed
  design <- .emDesign(obj, random)
  center <- .emStart(design, init, fixed, sd = control$sd %||% 0.3)
  control$modes <- cont$modes
  sd <- .perturbSd(sd, center)
  # every start and every seed is drawn here, so set.seed() reproduces a run whatever `cores`;
  # an inner axis above 1 runs the fits in PSOCK workers, as mstrust()
  starts <- c(list(center), lapply(seq_len(fits - 1L), function(i) .perturbStart(design, center, sd)))
  seeds <- sample.int(.Machine$integer.max, fits)
  # the starts and seeds of the retries come after, so those of the first attempts stay as they were
  maxTries <- if (isTRUE(retry)) max(1L, as.integer(nTries)) else 1L
  restarts <- lapply(seq_len(fits), function(i)
    lapply(seq_len(maxTries - 1L), function(k) .perturbStart(design, center, sd)))
  reseeds <- matrix(sample.int(.Machine$integer.max, fits * (maxTries - 1L)), fits)
  split <- .splitCores(cores, "fits")
  inner <- if (!is.null(split$conditions)) .sanitizeCores(split$conditions)
  # a worker returns the fit without the parts of the design, which the parent has; a fit with a
  # variance set to 0 keeps the design without that random effect
  shared <- c("prdfn", "data", "errfn", "design")
  # the message on the derivative mode travels with the fit and is given once
  one <- function(i) {
    suppressMessages(loadDLL(obj))
    note <- NULL
    for (a in seq_len(maxTries)) {
      set.seed(if (a == 1L) seeds[i] else reseeds[i, a - 1L])
      start <- if (a == 1L) starts[[i]] else restarts[[i]][[a - 1L]]
      fit <- tryCatch(withCallingHandlers({
        fit <- .emFit(design, start, fixed, method, control, FALSE)
        if (is.null(fit$boundary)) fit[shared] <- NULL
        fit
      }, dMod2DerivModeMessage = function(m) {
        note <<- conditionMessage(m)
        invokeRestart("muffleMessage")
      }), error = function(e) conditionMessage(e))
      if (!is.character(fit)) break
    }
    if (is.character(fit)) return(structure(fit, tries = a))
    fit$derivModeNote <- note
    fit$tries <- a
    fit
  }
  runs <- .parallelLapply(which, one, cores = .sanitizeCores(min(length(which), split$outer)),
                          coresConditions = inner)
  runs <- lapply(runs, function(r) if (is.null(r)) "the worker returned no result" else r)
  ok <- !vapply(runs, is.character, TRUE)
  notes <- unique(unlist(lapply(runs[ok], `[[`, "derivModeNote")))
  for (n in notes) .derivModeMessage(sub("\n$", "", n))
  runs[ok] <- lapply(runs[ok], function(r) { r$derivModeNote <- NULL; r })
  tries <- vapply(runs, function(r) as.integer(if (is.character(r)) attr(r, "tries") %||% 1L else r$tries %||% 1L), 0L)
  runs[ok] <- lapply(runs[ok], function(r) { r$tries <- NULL; r })
  table <- data.frame(start = which,
                      value = vapply(runs, function(r) if (is.character(r)) NA_real_ else r$value, 0),
                      converged = vapply(runs, function(r) !is.character(r) && isTRUE(r$converged), TRUE),
                      iterations = vapply(runs, function(r) if (is.character(r)) NA_integer_
                                          else as.integer(r$iterations %||% NA_integer_), 0L),
                      gradient = vapply(runs, function(r) if (is.character(r) || !length(r$gradient)) NA_real_
                                        else max(abs(r$gradient)), 0),
                      tries = tries,
                      stop = vapply(runs, function(r) if (is.character(r)) "error"
                                    else as.character(r$stopReason %||% NA_character_), ""),
                      failures = vapply(runs, function(r) if (is.character(r) || !length(r$failures)) ""
                                        else paste0(names(r$failures), " (", r$failures, ")", collapse = "; "), ""),
                      error = vapply(runs, function(r) if (is.character(r)) r else NA_character_, ""),
                      boundary = vapply(runs, function(r) if (is.character(r)) "" else paste(r$boundary, collapse = ", "), ""))
  if (verbose) print(table)
  if (!any(ok)) stop("msem: every fit failed; the first error: ", table$error[1], call. = FALSE)
  best <- runs[[which(ok)[which.min(table$value[ok])]]]
  if (is.null(best$boundary)) best[shared] <- list(design$prdfn, design$data, design$errfn, design)
  best$multistart <- table
  if (isTRUE(keepFits)) best$fits <- lapply(runs, function(r) { if (is.list(r)) r[shared] <- NULL; r })
  best
}


## The starts `which` of msem(), sorted.
.msemWhich <- function(which, fits) {
  if (is.null(which)) return(seq_len(fits))
  if (!is.numeric(which) || !length(which) || anyNA(which) || any(which != round(which)) ||
      any(which < 1 | which > fits) || anyDuplicated(which))
    stop("msem: `which` holds distinct positions among the ", fits, " starts.", call. = FALSE)
  sort(as.integer(which))
}

## The results of msem() for disjoint starts of one call joined: the best fit with the table of
## every start, by start, and the fits of every start where each part kept them.
.joinMsem <- function(parts) {
  table <- do.call(rbind, lapply(parts, `[[`, "multistart"))
  by <- order(table$start)
  best <- parts[[which.min(vapply(parts, function(p) p$value, 0))]]
  best$multistart <- `rownames<-`(table[by, , drop = FALSE], NULL)
  if (all(vapply(parts, function(p) !is.null(p$fits), NA)))
    best$fits <- unlist(lapply(parts, `[[`, "fits"), recursive = FALSE)[by]
  best
}



## A fit as `init`: the estimates of the fit as the start, the fixed values, the random effects without
## those it set to 0 at the boundary, whose etas are held at 0.
.emContinue <- function(obj, random, init, fixed) {
  if (!inherits(init, "em")) return(list(random = random, init = init, fixed = fixed))
  for (h in init$boundary) {
    design <- .emDesign(obj, random)
    b <- Filter(function(b) identical(paste("omega", b$short, b$short, sep = "_"), h), design$blocks)
    if (length(b) != 1L) next
    etas <- as.vector(design$etaNames[, b[[1]]$stems])
    fixed <- c(fixed, setNames(numeric(length(etas)), etas)[setdiff(etas, names(fixed))])
    random <- .dropRandomTerm(random, b[[1]]$stems)
  }
  fixed <- c(fixed, init$fixed[setdiff(names(init$fixed), names(fixed))])
  list(random = random, init = init$argument[setdiff(names(init$argument), names(fixed))], fixed = fixed,
       modes = init$innerModes)
}

## One fit of `design` from `start` by `method`, with `control`; -2 log L of the data at the
## start with the random effects at 0 must be finite.
.emFit <- function(design, start, fixed, method, control, verbose) {
  etas <- .allEtas(design)
  v <- design$objfn(pars = c(start, setNames(numeric(length(etas)), etas)), fixed = fixed, deriv = FALSE)$value
  if (!is.finite(v))
    stop("em: -2 log L of the data is not finite at the start (", format(v), "); an error model of 0 or ",
         "below at a data point gives that.", call. = FALSE)
  fit <- .emFitMethod(design, start, fixed, method, control, verbose)
  fit$control <- control
  fit$fixed <- fixed
  .emBoundary(design, fit, fixed, method, control, verbose)
}

## The fit with the normal random effects of single blocks whose standard deviation ended below 0.1
## set to 0 one at a time, smallest first, where the fit without one is at most control$boundary
## worse; every refit starts from the estimates and checks the random effects of the refit again.
.emBoundary <- function(design, fit, fixed, method, control, verbose) {
  tol <- control$boundary %||% 0.1
  if (isFALSE(tol) || !is.null(design$groupings) || is.null(design$normal)) return(fit)
  single <- Filter(function(b) b$dist == "normal" && length(b$stems) == 1L, design$blocks)
  if (!length(single)) return(fit)
  par <- c(fit$argument, fixed)
  hyper <- vapply(single, function(b) paste("omega", b$short, b$short, sep = "_"), "")
  stems <- vapply(single, `[[`, "", "stems")
  free <- hyper %in% names(fit$argument)
  logsd <- par[hyper]
  small <- which(free & is.finite(logsd) & logsd < log(0.1))
  if (length(design$blocks) <= 1L || !length(small)) return(fit)
  for (k in small[order(logsd[small])]) {
    random <- .dropRandomTerm(design$randomEffects, stems[k])
    etas <- as.vector(design$etaNames[, stems[k]])
    reduced <- tryCatch(suppressWarnings(.emFit(.emDesign(design$objfn, random),
                                                fit$argument[setdiff(names(fit$argument), hyper[k])],
                                                c(fixed, setNames(numeric(length(etas)), etas)),
                                                method, control, verbose)),
                        error = function(e) NULL)
    if (is.null(reduced) || !(reduced$value <= fit$value + tol)) next
    reduced$boundary <- c(hyper[k], reduced$boundary)
    return(reduced)
  }
  fit
}

.emFitMethod <- function(design, start, fixed, method, control, verbose) {
  normal <- all(vapply(design$blocks, function(b) b$dist == "normal", TRUE))
  if (!length(start) && normal) {
    # nothing to estimate: FOCEI evaluates the marginal at the fixed point, a structural
    # parameter moved into the start for an outer trust without steps
    if (method != "focei")
      stop("em: every parameter is fixed; method = \"focei\" evaluates the marginal there.", call. = FALSE)
    p <- setdiff(names(fixed), .hyperNames(design))[1]
    control$trustControl <- modifyList(control$trustControl %||% list(), list(iterlim = 0L))
    fit <- .emFitMethod(design, fixed[p], fixed[setdiff(names(fixed), p)], method, control, verbose)
    fit$fixed <- fixed
    fit$argument <- fit$argument[0]
    fit$converged <- TRUE
    return(fit)
  }
  if (!is.null(design$groupings)) {
    if (method != "focei" || any(vapply(design$blocks, function(b) b$dist != "normal", TRUE)))
      stop("em: random effects over several groupings are normal and fitted by method = \"focei\".",
           call. = FALSE)
    joint <- Reduce(`+`, c(list(design$objfn), lapply(design$groupings, .normalPrior)))
    return(.emResult(design, .fitFoceiJoint(design, joint, start, fixed, control), method, fixed))
  }
  if (any(vapply(design$blocks, function(b) b$dist != "normal", TRUE))) {
    if (method == "saem") return(.emResult(design, .fitSaemDirichlet(design, start, fixed, control, verbose),
                                           method, fixed))
    if (method != "focei")
      stop("em: spike-and-slab and Dirichlet-process random effects are fitted by method = \"focei\" ",
           "(Dirichlet-process ones also by \"saem\").", call. = FALSE)
    return(.emResult(design, .fitSparse(design, start, fixed, control, verbose), method, fixed))
  }
  joint <- design$objfn + .normalPrior(design)
  fit <- switch(method,
    focei = .fitFocei(design, joint, start, fixed, control),
    quadrature = {
      focei <- .fitFocei(design, joint, start, fixed, control)
      q <- .fitQuadrature(design, focei$argument, fixed, control, verbose, etaInit = focei$etaModes)
      q$trustIterations <- focei$iterations
      q
    },
    laplaceEM = .fitLaplaceEM(design, start, fixed, control, verbose),
    saem = .fitSaem(design, start, fixed, control, verbose))
  .emResult(design, fit, method, fixed)
}




## trust() for an inner step of an EM loop, whose iteration limit is part of the loop: that one
## warning is muffled, others pass.
.innerTrust <- function(args)
  withCallingHandlers(suppressMessages(do.call(trust, args)), warning = function(w)
    if (grepl("Maximum number of iterations", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning"))
