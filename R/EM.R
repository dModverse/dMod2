# Entry points of mixed-effects fits: em() for one fit, msem() for a multistart. Both build the
# estimation design from the objective of the data and the random effects of the parameter
# transformation, and maximise the marginal likelihood over the random effects.


#' Mixed-effects fit by marginal maximum likelihood
#'
#' @description
#' Fits a mixed-effects model: the random effects that [insert()] with effect
#' formulas put into the parameter transformation are integrated out, and the
#' fixed parameters, the covariance `Omega` of the normal random effects, the
#' strengths of the Laplace ones, the probability and slab standard
#' deviation of the spike-and-slab ones and the concentration and slab
#' standard deviation of the Dirichlet-process ones are estimated by
#' maximising the marginal likelihood. A group (a level of the grouping of a
#' random term, such as a patient or a cell line) may own several
#' conditions; its random effects are shared by all of them. Normal, Laplace,
#' spike-and-slab and Dirichlet-process random effects of one grouping may be
#' combined; random terms over several groupings, nested or crossed
#' (patients and centres), are normal and integrated per connected component
#' of the conditions.
#'
#' Methods:
#' \describe{
#'   \item{`"focei"`}{Laplace approximation around the conditional modes with
#'     the exact derivative of the volume term (FOCEI). The default.}
#'   \item{`"quadrature"`}{adaptive sparse-grid quadrature on nested
#'     Genz-Keister grids around the modes, started from the FOCEI fit, the
#'     level of each group set by a tolerance on -2 log L.}
#'   \item{`"laplaceEM"`}{ECM on the Laplace approximation with the
#'     covariance-corrected Omega update.}
#'   \item{`"saem"`}{stochastic-approximation EM: Metropolis draws of the
#'     random effects by random-walk and independence proposals at the
#'     conditional modes, Robbins-Monro averages of their sufficient
#'     statistics, the mu-referenced fixed effects and `Omega` in closed form;
#'     the other parameters take a trust step per burn-in iteration and then
#'     the Newton step on the Robbins-Monro average of the quadratic models of
#'     the complete-data objective. `control$saem`
#'     takes `nBurnin` and `nEM` (at most that many burn-in and averaging
#'     iterations, default 200 and 300), `nMcmc` and `nLaplace` (random-walk
#'     and independence steps per iteration, default 3 each), `laplaceEvery`
#'     (iterations between refreshes of the conditional modes, default 5),
#'     `laplaceScale` (factor of the proposal covariance, default 1.5),
#'     `window` (default 25) and `mcTol` (default 0.05): the burn-in ends when
#'     no parameter has a significant slope over the last `window` iterations,
#'     the averaging, after at least `2 * window` iterations, when the expected
#'     increase of -2 log L from the Monte Carlo error of the estimate, the mean
#'     of the iterates after the first `window`, is below `mcTol`. With
#'     Dirichlet-process random effects a Gibbs sweep draws the clusters of
#'     the groups (Neal 2000, algorithm 8) and Metropolis steps their values;
#'     the concentration and the slab follow from averaged numbers of
#'     clusters and squared values, the averaging runs `nEM` iterations, the
#'     reported clusters are the most frequent partition of the late draws,
#'     and -2 log L comes from sequential importance sampling over the
#'     partitions with `particles` (default 500) particles, its standard
#'     error in `valueSE`.}
#' }
#' Laplace and spike-and-slab random effects and several groupings are fitted
#' by `"focei"`, Dirichlet-process ones by `"focei"` (at most 14 groups) or
#' `"saem"`.
#'
#' @param obj Objective of the data, [normL2()] or [normT()] of the model
#'   `g * x * p`, without terms for the random effects.
#' @param random The random effects, [randomEffects()] of the parameter
#'   transformation or of `p`.
#' @param init Named numeric start. It needs the parameters of `obj` that
#'   are neither random effects nor created by the effect formulas (those
#'   start at 0); it may set those and the hyperparameters `omega_*` (log
#'   standard deviations and Cholesky entries) and, named
#'   `<kind>_<grouping>` and numbered `_1`, `_2`, ... when a grouping has
#'   several such blocks, `lambda_*` (Laplace), `prob_*` and `tau_*`
#'   (spike and slab), `concentration_*` and `tau_*` (Dirichlet process).
#'   An `em` fit as `init` gives its estimates as the start and its `fixed`
#'   values, and its random effects without those it set to 0 at the boundary.
#' @param fixed Named numeric of parameters held fixed.
#' @param method `"focei"` (default), `"quadrature"`, `"laplaceEM"` or
#'   `"saem"`.
#' @param control List of options: `sd` (start of the random-effect standard
#'   deviations, default 0.3), `innerControl` and `trustControl` (the inner
#'   and outer trust region of FOCEI; the inner one stops where the Newton
#'   step of the random effects is below `steptol`, default 1e-8, the outer one
#'   where the predicted change of -2 log L is below `mtol`, default 1e-5), `secondOrderCorrection` (default
#'   `TRUE`), `derivMode`, how FOCEI gets the second derivatives of its
#'   volume term: `"forward-reverse"` (default), directional sweeps forward
#'   along the random effects and backward with seeds, for models built with
#'   `derivMode = c("forward", "reverse", "forward-reverse")`: the outer model
#'   Hessian Gauss-Newton, the volume curvature learned by the SR1 correction;
#'   the reported Hessian forward over forward where
#'   the model also has `deriv2`, else the exact Schur complement of the joint
#'   Hessian (up to 500 parameters) with the volume curvature of the variance
#'   parameters, exact for models linear in the random effects; or `"forward-forward"`, the
#'   second-order sensitivities of models built with `deriv2 = TRUE`, with the
#'   exact volume Hessian. In both modes an SR1 correction from the changes of
#'   the outer gradient between the last three accepted points adds the curvature
#'   Gauss-Newton leaves out, negative curvature included, so the outer model
#'   may be indefinite. Without the mode asked for the fit takes the other
#'   with a message, also where the default falls back to `"forward-forward"`;
#'   [msem()] gives it once. `quadrature` (section Quadrature of [em()], also read by
#'   `"laplaceEM"`) and `saem` lists for those methods; `laplace`, the list of
#'   the fits with Laplace, spike-and-slab or Dirichlet-process random
#'   effects, whose parameters and hyperparameters move in one trust region:
#'   `maxOuter` (default 200), its iterations; `trustControl`, its [trust()]
#'   settings; `coupling`, `"expectationConsistent"` (default) or
#'   `"meanField"`, the integral where the data couple sparse random effects
#'   of one group that the priors keep apart. The coupling holds for the
#'   whole fit: where expectation consistent does not converge at the start,
#'   the fit is by mean field, and `fit$coupling` is `"meanFieldFallback"`; a
#'   trial point at which it does not converge is a failed evaluation of the
#'   trust region. `boundary`, the largest rise of -2 log L at which a
#'   variance is set to 0, default 0.1, `FALSE` for none: with one grouping
#'   and several blocks, a normal random effect of its own block whose
#'   standard deviation ends below 0.1 is dropped, its etas fixed at 0, and
#'   the other parameters are fitted again from the estimates; where that
#'   raises -2 log L by at most `boundary`, the fit without it is returned and
#'   names its hyperparameter in `boundary`. Other entries are an error.
#' @param cores Number of cores over the groups and conditions inside the
#'   fit; with 1 (default) no threads are started.
#' @param verbose Logical; print the progress of the iterative methods.
#' @param sparse `"auto"` (default), `TRUE` or `FALSE`: the outer trust region of
#'   FOCEI over normal random effects on the sparse Gauss-Newton Hessian, the
#'   Schur complement and the volume term of every group local to its
#'   parameters. `"auto"` decides as [trust()] does, from the number of outer
#'   parameters and the fill of the outer Hessian at the start.
#' @return An object of class `em`: `argument`, `value` (-2 log L of the
#'   marginal likelihood), `gradient`, `hessian` (a sparse `Matrix` on the
#'   sparse path, `NULL` with Laplace, spike-and-slab or Dirichlet-process
#'   random effects), `information` (with
#'   Laplace or spike-and-slab random effects the empirical information that
#'   gives the standard errors: `G / (G - 1)` times the sum over the `G` groups
#'   of the outer products of their centred scores, the gradients of their
#'   log L, plus half the Hessian of -2 log L of the conditions without random
#'   effects), `converged`, `iterations`,
#'   `Omega`, `omega` (standard deviations and correlations with standard
#'   errors), `random` (the modes of the random effects per group with
#'   conditional standard errors; with spike-and-slab or Dirichlet-process
#'   random effects also the selected cluster, 0 for the fixed effect, and
#'   the posterior probability of a deviation of 0), `strengths` and `slab`
#'   (the hyperparameters of the Laplace blocks and of the spike-and-slab and
#'   Dirichlet-process blocks), `valueSE` (the standard error of `value` when
#'   it comes from importance sampling), `coupling` (where the data couple
#'   sparse random effects: `"expectationConsistent"`, `"meanField"` or
#'   `"meanFieldFallback"`), `boundary` (the hyperparameters of the random
#'   effects set to 0, see `control`), `etaModes`, `etaSE`, `shrinkage` and
#'   the design; with method `"quadrature"` also `certificate`, per group its
#'   level, -2 log L by quadrature and by importance sampling with standard
#'   error, the Pareto shape of the weights, the estimated error of its level
#'   and a flag. Methods: `print`, `summary`, `coef`, `vcov`, [confint][confint.em],
#'   [logLik][logLik.em], [residuals][residuals.em], [predict][predict.em],
#'   `plot`.
#' @section Quadrature:
#' Method `"quadrature"` alternates between placing the nodes around the
#' modes at the current parameters and one trust step over every parameter
#' on these frozen nodes, with the exact Hessian where the model has second
#' derivatives. The nodes of the mu-referenced fixed effects (random
#' intercepts whose formula covers every condition with data, covariates
#' constant within a group) are frozen in `phi = X beta + eta`, so these
#' fixed effects need no sensitivities between two rebuilds. Entries of
#' `control$quadrature`:
#' \describe{
#'   \item{`tol`}{allowed change of -2 log L by the next level, default 0.1,
#'     split evenly over the groups: each group takes the lowest level from
#'     its stage up whose next level changes its -2 log L by at most
#'     `tol / G`, and never a lower level than before. `Inf` keeps the levels
#'     of the stages.}
#'   \item{`levels`}{the stages, default `K + 1` with `K` the number of random
#'     effects of a group; level `K + s` integrates polynomials of total
#'     degree `2 s + 1` exactly.}
#'   \item{`maxLevel`}{highest level, default `K + 8`.}
#'   \item{`tolDecrement`}{a stage ends when a Newton step on the rebuilt
#'     nodes would lower -2 log L by less than `tolDecrement`, default 1e-4,
#'     or than the nodes resolve, the change of -2 log L by the next level
#'     summed over the groups; or by less than the last step raised it on
#'     rebuilt nodes of the same levels, at its iterate of lowest -2 log L.
#'     Each counts as converged.}
#'   \item{`maxEcmPerStage`, `cm1Control`, `maxCm1Iter`}{rebuilds per stage
#'     (default 10), [trust()] settings and iterations of a step (default 30).}
#'   \item{`pruneTol`}{nodes whose weight lies this much in log below the
#'     largest are dropped, default `Inf`.}
#'   \item{`certificate`}{`FALSE`, or a list of `points` (default 256) and
#'     `seed` (default 1) of the importance-sampling check at the estimate: a
#'     multivariate t5 proposal at the mode and conditional covariance of each
#'     group, randomly shifted lattice points. A group is flagged when the two
#'     -2 log L differ by more than three standard errors plus the error its
#'     level allows (the larger of `tol / G` and its estimated level error),
#'     or when the Pareto shape of the weights exceeds 0.7.}
#' }
#' Method `"laplaceEM"` reads `epsEcm` (default 1e-5), the extrapolated
#' distance of the parameters to the fixed point at which it stops, `maxEcm`,
#' `cm1Control` and `maxCm1Iter`.
#'
#' @seealso [msem()], [randomEffects()], [insert()]
#' @example inst/examples/NLME.R
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
  .withCores(cores, .emFit(design, start, fixed, method, control, verbose))
}

## Refuses entries of `control` and of `control$laplace` that em() does not read.
.checkEmControl <- function(control) {
  known <- list(control = c("sd", "innerControl", "trustControl", "secondOrderCorrection", "derivMode", "quadrature",
                            "saem", "laplace", "boundary"),
                laplace = c("maxOuter", "trustControl", "coupling"))
  for (k in names(known)) {
    x <- if (k == "control") control else control[[k]]
    bad <- setdiff(names(x), known[[k]])
    if (length(bad))
      stop("em: unknown ", if (k == "control") "entries of `control`: " else "entries of `control$laplace`: ",
           paste(bad, collapse = ", "), ".", call. = FALSE)
  }
  dm <- control$derivMode
  if (!is.null(dm) && !(length(dm) == 1L && dm %in% c("forward-reverse", "forward-forward")))
    stop("em: `control$derivMode` is \"forward-reverse\" or \"forward-forward\".", call. = FALSE)
  bd <- control$boundary
  if (!is.null(bd) && !isFALSE(bd) && !(is.numeric(bd) && length(bd) == 1L && is.finite(bd) && bd >= 0))
    stop("em: `control$boundary` is FALSE or one number of at least 0.", call. = FALSE)
  cp <- control$laplace$coupling
  if (!is.null(cp) && !(length(cp) == 1L && cp %in% c("expectationConsistent", "meanField")))
    stop("em: `control$laplace$coupling` is \"expectationConsistent\" or \"meanField\".", call. = FALSE)
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
#' strengths, concentrations and slab standard deviations, on the logit scale
#' for slab probabilities), and returns the best fit.
#'
#' @inheritParams em
#' @param fits Number of starts. Default 20.
#' @param cores Number of fits run in parallel, or `c(fits = , conditions = )`
#'   as in [mstrust()]: fits in parallel and cores inside each fit.
#' @param sd Standard deviation of the perturbation of the starts.
#' @param budget The fit from `init` runs first; in every other start the
#'   outer trust region of FOCEI, also the one before quadrature, and that of
#'   the fits with sparse random effects end after `budget` times the
#'   iterations they took in the first fit. Quadrature, Laplace-EM and SAEM
#'   keep their own iteration limits. `Inf` lets every start run to
#'   convergence.
#' @example inst/examples/NLME.R
#' @return The fit with the smallest `value`, of class `em`, with `multistart`,
#'   a data.frame of every start with its `value`, `converged`, `iterations`,
#'   `error` and `boundary`, the variances its fit set to 0; a start that the
#'   budget ends has `converged = FALSE`.
#'   [plotValues()] draws it as a waterfall.
#' @seealso [em()]
#' @export
msem <- function(obj, random, init = NULL, fixed = NULL,
                 method = c("focei", "quadrature", "laplaceEM", "saem"),
                 control = list(), fits = 20L, cores = 1L, sd = 0.5, budget = 5, verbose = FALSE, sparse = "auto") {
  method <- match.arg(method)
  .checkEmControl(control)
  if (length(budget) != 1L || is.na(budget) || budget <= 0)
    stop("msem: `budget` is a positive number, Inf for none.", call. = FALSE)
  control$sparse <- .checkSparse(sparse, "msem")
  cont <- .emContinue(obj, random, init, fixed)
  random <- cont$random; init <- cont$init; fixed <- cont$fixed
  design <- .emDesign(obj, random)
  center <- .emStart(design, init, fixed, sd = control$sd %||% 0.3)
  # every start and every seed is drawn here, so set.seed() reproduces a run whatever `cores`;
  # an inner axis above 1 runs the fits in PSOCK workers, as mstrust()
  starts <- c(list(center), lapply(seq_len(fits - 1L), function(i) .perturbStart(design, center, sd)))
  seeds <- sample.int(.Machine$integer.max, fits)
  split <- .splitCores(cores, "fits")
  inner <- if (!is.null(split$conditions)) .sanitizeCores(split$conditions)
  # a worker returns its fit without the parts of the design, which the parent has; a fit with a
  # variance set to 0 keeps the design without that random effect
  shared <- c("prdfn", "data", "errfn", "design")
  # the message on the derivative mode travels with the fit and is given once
  one <- function(i, cap = NULL) {
    suppressMessages(loadDLL(obj))
    set.seed(seeds[i])
    note <- NULL
    tryCatch(withCallingHandlers({
      fit <- .emFit(design, starts[[i]], fixed, method, .capIterations(control, cap), FALSE)
      if (is.null(fit$boundary)) fit[shared] <- NULL
      fit$derivModeNote <- note
      fit
    }, dMod2DerivModeMessage = function(m) {
      note <<- conditionMessage(m)
      invokeRestart("muffleMessage")
    }), error = function(e) conditionMessage(e))
  }
  # the budget of the other starts comes from the first, so it does not depend on `cores`
  first <- .parallelLapply(1L, one, coresConditions = inner)
  cap <- .msemCap(first[[1]], method, budget)
  runs <- c(first, if (fits > 1L) .parallelLapply(seq_len(fits)[-1], function(i) one(i, cap),
                                                  cores = .sanitizeCores(min(fits - 1L, split$outer)),
                                                  coresConditions = inner))
  runs <- lapply(runs, function(r) if (is.null(r)) "the worker returned no result" else r)
  ok <- !vapply(runs, is.character, TRUE)
  notes <- unique(unlist(lapply(runs[ok], `[[`, "derivModeNote")))
  for (n in notes) .derivModeMessage(sub("\n$", "", n))
  runs[ok] <- lapply(runs[ok], function(r) { r$derivModeNote <- NULL; r })
  table <- data.frame(start = seq_len(fits),
                      value = vapply(runs, function(r) if (is.character(r)) NA_real_ else r$value, 0),
                      converged = vapply(runs, function(r) !is.character(r) && isTRUE(r$converged), TRUE),
                      iterations = vapply(runs, function(r) if (is.character(r)) NA_integer_
                                          else as.integer(r$iterations %||% NA_integer_), 0L),
                      error = vapply(runs, function(r) if (is.character(r)) r else NA_character_, ""),
                      boundary = vapply(runs, function(r) if (is.character(r)) "" else paste(r$boundary, collapse = ", "), ""))
  if (verbose) print(table)
  if (!any(ok)) stop("msem: every fit failed; the first error: ", table$error[1], call. = FALSE)
  best <- runs[[which(ok)[which.min(table$value[ok])]]]
  if (is.null(best$boundary)) best[shared] <- list(design$prdfn, design$data, design$errfn, design)
  best$multistart <- table
  best
}


## The iteration cap of the other starts: `budget` times the outer trust iterations of the first
## fit `first`, those of the FOCEI stage before quadrature; NULL for the methods it does not cap.
.msemCap <- function(first, method, budget) {
  if (is.character(first) || !method %in% c("focei", "quadrature") || !is.finite(budget)) return(NULL)
  done <- suppressWarnings(as.numeric(first$trustIterations %||% first$iterations)[1])
  if (!is.finite(done)) return(NULL)
  as.integer(ceiling(budget * max(done, 1)))
}


## A fit as `init`: its estimates as the start and its fixed values, its random effects without
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
  list(random = random, init = init$argument[setdiff(names(init$argument), names(fixed))], fixed = fixed)
}

## One fit of `design` from `start` by `method`, with its `control`; -2 log L of the data at the
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
## worse; every refit starts from the estimates and checks its own random effects again.
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
      stop("em: Laplace, spike-and-slab and Dirichlet-process random effects are fitted by method = \"focei\" ",
           "(Dirichlet-process ones also by \"saem\").", call. = FALSE)
    return(.emResult(design, .fitLaplace(design, start, fixed, control, verbose), method, fixed))
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


## `control` with the outer iterations of FOCEI and of the fits with sparse random effects at most
## `cap` (their defaults 200).
.capIterations <- function(control, cap) {
  if (is.null(cap)) return(control)
  control$trustControl$iterlim <- min(control$trustControl$iterlim %||% 200L, cap)
  control$laplace$maxOuter <- min(control$laplace$maxOuter %||% 200L, cap)
  control
}


## trust() for an inner step of an EM loop, whose iteration limit is part of the loop: that one
## warning is muffled, others pass.
.innerTrust <- function(args)
  withCallingHandlers(suppressMessages(do.call(trust, args)), warning = function(w)
    if (grepl("Maximum number of iterations", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning"))
