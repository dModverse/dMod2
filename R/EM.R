## Unified mixed-effects entry points. EM()/msEM()/emInit() dispatch on the
## random-effect prior attached to the objective: a Gaussian random effect from
## omega() (constraintL2 with an Omega) or an L1/Laplace random effect from
## penaltyL1() (constraintL1). Both are marginal maximum likelihood over that
## prior; only the prior and the way its parameter is estimated differ.

## Which random-effect prior does `obj` carry? Errors if neither or both.
.emPrior <- function(obj) {
  pen <- attr(obj, "penaltySpec", exact = TRUE)
  om  <- attr(obj, "omegaSpec",   exact = TRUE)
  if (!is.null(pen) && !is.null(om))
    stop("EM: `obj` carries both a penaltyL1() and an omega() random effect; ",
         "use exactly one.", call. = FALSE)
  if (!is.null(pen)) return("penaltyL1")
  if (!is.null(om))  return("omega")
  stop("EM: `obj` carries no random-effect specification. Add either ",
       "+ constraintL1(penaltyL1(...)) or an omega() random effect.",
       call. = FALSE)
}


#' Mixed-effects fit by marginal maximum likelihood
#'
#' @description
#' Fits an ODE mixed-effects (population) model by marginal maximum likelihood,
#' dispatching on the random-effect prior attached to `obj`:
#' \describe{
#'   \item{[omega]}{a Gaussian random effect with covariance `Omega`. Backends
#'     `"focei"`, `"quadrature"`, `"foceiQuadrature"`, `"laplaceEM"`, `"saem"`.}
#'   \item{[penaltyL1]}{an L1/Laplace random effect whose single strength
#'     `lambda` is estimated (the sparsity / clustering prior). Backends
#'     `"focei"` and `"saem"`.}
#' }
#' The outer marginal is FOCE(I) / EM / SAEM; the inner conditional mode uses
#' [trust] (Gaussian) or [trustL1] (Laplace).
#'
#' @param obj Composed objective (`normL2(...) + constraintL2(...)` for a
#'   Gaussian random effect, or `normL2(...) + constraintL1(penaltyL1(...))` for
#'   an L1/Laplace one) carrying the model pieces and the random-effect spec.
#' @param init Named numeric start; assemble with [emInit].
#' @param fixed Optional named numeric of parameters held fixed.
#' @param method Marginal-likelihood backend. `NULL` (default) uses `"focei"`.
#'   Allowed values depend on the prior (see Description).
#' @param control List of backend tuning options; see [EM] (Gaussian) and
#'   the penalty backends for the recognised fields.
#' @param verbose Logical; print progress.
#' @return An object of class `c("em", "list")` carrying `prior`
#'   (`"omega"`/`"penaltyL1"`), the estimated parameters (`argument`), objective
#'   `value` (-2 log L), the per-subject conditional modes (`etaModes`), and the
#'   prior-specific readout (`Omega` for a Gaussian fit; `lambda` / `sparsity` /
#'   `clusters` for an L1 fit). Methods: [print][print.em], `coef`, and (Gaussian
#'   only) `confint`, `summary`, `predict`, `plot`.
#' @seealso [msEM], [emInit], [omega], [penaltyL1], [sparsify]
#' @example inst/examples/EMlaplace.R
#' @example inst/examples/NLME.R
#' @export
EM <- function(obj, init, fixed = NULL, method = NULL, control = list(),
               verbose = TRUE) {
  ## global method check first (NULL -> "focei"); the per-prior branch narrows it.
  method <- match.arg(method, c("focei", "quadrature", "foceiQuadrature",
                                "laplaceEM", "saem"))
  if (identical(.emPrior(obj), "penaltyL1")) {
    method <- match.arg(method, c("focei", "saem"))
    .fitLaplace(obj, init, fixed = fixed, method = method, control = control,
           verbose = verbose)
  } else {
    .fitNormal(obj, init, fixed = fixed, method = method, control = control,
            verbose = verbose)
  }
}


#' Multi-start mixed-effects fit
#'
#' @description
#' The multi-start pendant of [EM]: runs [EM] from several perturbed starts and
#' returns the best fit, dispatching on the random-effect prior of `obj`.
#'
#' @param obj,fixed,method,control,verbose As in [EM].
#' @param center Named numeric start (see [emInit]); the multi-start centre.
#' @param fits Number of starts (default 20).
#' @param cores Number of parallel workers.
#' @param ... Prior-specific options forwarded to the backend (e.g. `sd`,
#'   `samplefun`).
#' @return For a Gaussian (`omega`) prior, a [parlist] of the fits; for an L1
#'   (`penaltyL1`) prior, the best [EM] fit carrying `$msTable`.
#' @seealso [EM], [emInit]
#' @export
msEM <- function(obj, center, fixed = NULL, method = NULL, control = list(),
                 fits = 20L, cores = 1L, verbose = FALSE, ...) {
  method <- match.arg(method, c("focei", "quadrature", "foceiQuadrature", "saem"))
  if (identical(.emPrior(obj), "penaltyL1")) {
    method <- match.arg(method, c("focei", "saem"))
    .msfitLaplace(obj, center, fixed = fixed, method = method, control = control,
             fits = fits, cores = cores, verbose = verbose, ...)
  } else {
    method <- match.arg(method, c("focei", "quadrature", "foceiQuadrature"))
    .msfitNormal(obj, center, fixed = fixed, method = method, control = control,
              fits = fits, cores = cores, verbose = verbose, ...)
  }
}


#' Validated start for a mixed-effects fit
#'
#' @description
#' Builds a validated named start vector from a structural start and a
#' random-effect specification, dispatching on the spec: an [omega] gives a
#' Gaussian start (with Cholesky parameters), a [penaltyL1] gives an L1 start
#' (with the penalty strength `lambda`).
#'
#' @param structural Fully named numeric of structural (and error-model) starts.
#' @param prior An [omega] or [penaltyL1] specification.
#' @param ... Forwarded to the spec-specific builder: `sd` (Gaussian, the initial
#'   random-effect SD) or `lambda` (L1, the initial penalty strength).
#' @return A named numeric start for [EM] / [msEM].
#' @seealso [EM], [omega], [penaltyL1]
#' @export
emInit <- function(structural, prior, ...) {
  if (inherits(prior, "penaltyspec")) return(.initLaplace(structural, prior, ...))
  if (inherits(prior, "omegaspec"))   return(.initNormal(structural, prior, ...))
  stop("emInit: `prior` must be a penaltyL1() or omega() specification.",
       call. = FALSE)
}


#' Structural model selection over the random-effect support
#'
#' @description
#' Cheap marginal-likelihood selection of which parameters are individual-specific
#' (and, for a clustered penalty, how the individuals group), dispatching on the
#' [penaltyL1] method attached to `obj`:
#' \describe{
#'   \item{`"lasso"` / `"fused"`}{Forward nested-support selection: one
#'     mixed-effects fit ranks the candidates by deviation magnitude, then an O(K)
#'     nested chain of supports is compared by marginal \eqn{-2\log L}. Returns the
#'     selected `support`.}
#'   \item{`"clustered"`}{Grouping selection: one all-separate fit provides the
#'     FOCE linearisation, then the clustered marginal scores the cluster path and
#'     selects the grouping per candidate parameter (jointly by default). Returns
#'     the recovered `perParam` grouping and a [plot][plot.sparsify] method.}
#' }
#' The marginal integration supplies the Occam factor that a classical penalty
#' scan needs a lambda grid plus an information criterion to approximate.
#'
#' @param obj Composed objective for the full candidate set, carrying a
#'   [penaltyL1] specification (via [constraintL1]).
#' @param center Named numeric start (see [emInit]).
#' @param fixed Optional named numeric held fixed.
#' @param method Backend for the support path (`"focei"` default, or `"saem"`);
#'   ignored for the clustered path.
#' @param control List forwarded to the backend.
#' @param fits Number of structural multistarts (default `12` for the support
#'   path, `1` for the clustered path).
#' @param cores Number of parallel workers.
#' @param sd Perturbation SD for the multistart.
#' @param verbose Logical; print the selection.
#' @param ... Clustered-path options (`quadNodes`, `quadRule`, `level`, `lamGrid`,
#'   `refine`, `joint`); see the clustered backend.
#' @return An object of class `c("sparsify", "list")`; `type = "support"` for the
#'   nested-support path (`support`, `ranking`, `chain`, `fitSelected`) or
#'   `type = "grouping"` for the clustered path (`perParam`, `structural`,
#'   `etaModes`), with a [print][print.sparsify] method (and, for a grouping,
#'   [plot][plot.sparsify]).
#' @seealso [EM], [penaltyL1], [constraintL1]
#' @example inst/examples/EMclustered.R
#' @export
sparsify <- function(obj, center, fixed = NULL, method = NULL, control = list(),
                     fits = NULL, cores = 1L, sd = 0.5, verbose = TRUE, ...) {
  pen <- attr(obj, "penaltySpec", exact = TRUE)
  if (is.null(pen))
    stop("sparsify: `obj` carries no penaltyL1() specification. Add ",
         "+ constraintL1(penaltyL1(..., subjects = ...)).", call. = FALSE)
  clustered <- any(vapply(pen$blocks, `[[`, "", "kind") == "clustered")
  if (clustered)
    return(.selectCluster(obj, center, fixed = fixed, control = control,
                            fits = if (is.null(fits)) 1L else fits,
                            cores = cores, sd = sd, verbose = verbose, ...))
  m <- if (is.null(method)) "focei" else match.arg(method, c("focei", "saem"))
  .selectLaplace(obj, center, fixed = fixed, method = m, control = control,
            fits = if (is.null(fits)) 12L else fits,
            cores = cores, sd = sd, verbose = verbose)
}
