## Restricted maximum likelihood for the error model -------------------------
##
## The error parameters solve the maximum-likelihood equation with every data
## point charged its own leverage h_ii:
##   sum_i [ 1 - h_ii - r_i^2/sigma_i^2 ] d log sigma_i^2 / d phi_k = 0
## Only first-order sensitivities enter, and the fixed point is exact although
## each step freezes h.


# One evaluation of the whole chain per L2 term, returned per term and
# condition as the objframe that res() produces: sigma, residuals,
# d pred / d par and d sigma / d par, all on the outer parameter scale.
.remlFrames <- function(objfun, pars, fixed = NULL, cores = 1L) {

  spec <- .remlTerms(objfun)
  frames <- list()
  labels <- character(0)

  for (t in seq_along(spec)) {
    term  <- spec[[t]]
    times <- term$timesD
    if (is.null(times))
      times <- sort(unique(c(0, unlist(lapply(term$data, `[[`, "time")))))

    prediction <- term$prdfn(times, pars, fixed = fixed, deriv = TRUE,
                             cores = cores)
    conditions <- intersect(names(prediction), names(term$data))

    for (cn in conditions) {
      errout <- NULL
      if (!is.null(term$errfn)) {
        pinner     <- getParameters(prediction[[cn]])
        fixedinner <- pinner[attr(pinner, "fixed")]
        errout <- term$errfn(
          out   = prediction[[cn]],
          pars  = as.parvec(pinner[setdiff(names(pinner), names(fixed))]),
          fixed = as.parvec(fixedinner, deriv = FALSE, deriv2 = FALSE),
          conditions = cn)[[cn]]
      }
      frames <- c(frames, list(res(term$data[[cn]], prediction[[cn]], errout)))
      labels <- c(labels, if (length(spec) > 1L) paste0(t, ":", cn) else cn)
    }
  }

  if (!length(frames))
    stop("reml: the objective carries no data to evaluate.", call. = FALSE)
  structure(frames, names = labels)
}


# The L2 terms of an objective. normL2() attaches one, "+" concatenates them.
# The single prdfn/data/errfn handles are first-wins on composition, so they
# only stand in when no term list is present.
.remlTerms <- function(objfun) {

  spec <- attr(objfun, "l2spec", exact = TRUE)
  if (!is.null(spec) && length(spec)) return(spec)

  x    <- attr(objfun, "prdfn")
  data <- attr(objfun, "data")
  if (is.null(x) || is.null(data))
    stop("objfun does not carry a prediction function and data. ",
         "REML needs an objective built by normL2().", call. = FALSE)
  list(list(data = data, prdfn = x, errfn = attr(objfun, "errfn"),
            timesD = attr(objfun, "timesD")))
}


# Stack the per-row derivative blocks of every frame into one matrix over
# `cols`, filling with zero where a term does not carry a parameter.
.remlBlocks <- function(frames, which, cols) {

  n <- vapply(frames, nrow, 0L)
  out <- matrix(0, sum(n), length(cols), dimnames = list(NULL, cols))
  seen <- character(0)
  off <- 0L
  for (i in seq_along(frames)) {
    d <- attr(frames[[i]], which)
    idx <- off + seq_len(n[i])
    off <- off + n[i]
    if (is.null(d)) next
    keep <- intersect(cols, colnames(d))
    seen <- union(seen, keep)
    if (length(keep)) out[idx, keep] <- d[, keep, drop = FALSE]
  }
  attr(out, "seen") <- seen
  out
}


# Stack a per-row quantity over conditions in the order .remlFrames returns.
.remlStack <- function(frames, what) {
  unlist(lapply(frames, function(z) z[[what]]), use.names = FALSE)
}


#' Leverage of the data points in the mean model
#'
#' @description
#' Hat values of the weighted mean model,
#' \eqn{h = \mathrm{diag}(W^{1/2} J (J^\top W J)^{-1} J^\top W^{1/2})} with
#' \eqn{J = \partial\mu/\partial\theta} the prediction sensitivities and
#' \eqn{W = \mathrm{diag}(1/\sigma_i^2)}. They say how much of the parameter
#' budget each data point, and by aggregation each observable, spends.
#'
#' Their sum is the numerical rank of the weighted sensitivity matrix, not the
#' nominal parameter count: a non-identifiable direction costs nothing.
#'
#' @param objfun objective function built by [normL2], carrying its prediction
#'   function, data and error model.
#' @param pars named numeric parameter vector, usually a fit.
#' @param meanpars character, the mean parameters spanning \eqn{J}. Defaults to
#'   `names(pars)` minus the error-model parameters.
#' @param fixed named numeric passed on to the prediction.
#' @param rank.tol relative threshold on the singular values of the weighted
#'   sensitivity matrix below which a direction counts as absent.
#' @param cores passed on to the prediction.
#'
#' @return A `data.frame` with one row per data point, columns `condition`,
#'   `time`, `name`, `sigma`, `residual` and `leverage`, and attributes
#'   `"rank"` (the effective number of mean parameters) and `"dof"` (a table of
#'   \eqn{n_g - \sum_{i \in g} h_{ii}} per observable).
#'
#' @seealso [reml] uses these to correct the error model.
#' @importFrom stats setNames
#' @export
remlLeverage <- function(objfun, pars, meanpars = NULL, fixed = NULL,
                         rank.tol = 1e-8,
                         cores = getOption("dMod.cores", 1L)) {

  frames <- .remlFrames(objfun, pars, fixed = fixed, cores = cores)
  if (is.null(meanpars)) meanpars <- setdiff(names(pars), .remlErrpars(objfun, pars))

  if (any(vapply(frames, function(z) is.null(attr(z, "deriv")), TRUE)))
    stop("remlLeverage: the prediction returned no sensitivities.", call. = FALSE)

  J <- .remlBlocks(frames, "deriv", meanpars)
  missing <- setdiff(meanpars, attr(J, "seen"))
  if (length(missing))
    stop("remlLeverage: meanpars not among the prediction parameters: ",
         paste(missing, collapse = ", "), call. = FALSE)

  sigma <- .remlStack(frames, "sigma")
  J <- J / sigma

  s <- svd(J)
  k <- sum(s$d > rank.tol * max(s$d))
  h <- rowSums(s$u[, seq_len(k), drop = FALSE]^2)
  # log det(J' W J) over the directions that are present, the penalty of the
  # restricted likelihood. J' W J has the squared singular values as eigenvalues.
  logdet <- 2 * sum(log(s$d[seq_len(k)]))

  out <- data.frame(
    condition = rep(names(frames), vapply(frames, nrow, 0L)),
    time      = .remlStack(frames, "time"),
    name      = .remlStack(frames, "name"),
    sigma     = sigma,
    residual  = .remlStack(frames, "residual"),
    leverage  = h,
    stringsAsFactors = FALSE)

  n.g   <- table(out$name)
  tr.g  <- tapply(h, out$name, sum)
  attr(out, "rank")   <- k
  attr(out, "logdet") <- logdet
  attr(out, "dof")    <- as.numeric(n.g) - as.numeric(tr.g[names(n.g)])
  names(attr(out, "dof")) <- names(n.g)
  out
}


# Error-model parameters that reach the estimated set: the free symbols of the
# error equations of every term, minus the observables they define.
.remlErrpars <- function(objfun, pars) {
  errs <- lapply(.remlTerms(objfun), `[[`, "errfn")
  errs <- errs[!vapply(errs, is.null, TRUE)]
  if (!length(errs)) return(character(0))
  free <- unlist(lapply(errs, function(err) {
    eqs <- unlist(getEquations(err))
    setdiff(getSymbols(eqs), names(eqs))
  }))
  intersect(unique(free), names(pars))
}


#' Estimate the error model by restricted maximum likelihood
#'
#' @description
#' Alternates between fitting the mean parameters at a fixed error model and
#' updating the error parameters from the REML stationarity condition, until
#' the error parameters stop moving. Every data point is charged its own
#' leverage instead of an equal share of the parameter budget.
#'
#' The objective is left untouched, it stays the plain \eqn{-2\log L}. Every L2
#' term of a composed objective is used, so a split
#' `normL2(d1, ...) + normL2(d2, ...)` is handled as one dataset.
#'
#' @details
#' Only first-order sensitivities are used. The condition solved is the data
#' term's: a term acting on the error parameters from outside it, a prior for
#' instance, is not part of it, and `reml()` warns when it finds one. Put such
#' terms on the mean parameters.
#'
#' @param objfun objective function built by [normL2], with an error model.
#' @param pars named numeric starting vector covering mean and error parameters.
#' @param errpars character, the error-model parameters. Defaults to the free
#'   symbols of the error equations that appear in `pars`.
#' @param fixed named numeric held fixed throughout.
#' @param iterlim maximum number of outer rounds.
#' @param tol convergence threshold on the largest change of an error parameter
#'   between two rounds.
#' @param rank.tol passed to [remlLeverage].
#' @param optimizer optimiser for the mean step, called as
#'   `optimizer(objfun, parinit, fixed = , ...)`.
#' @param control named list of arguments for `optimizer`.
#' @param cores passed on to the prediction and the objective.
#' @param ... further arguments for `objfun`.
#'
#' @return A list with components `argument` (the full parameter vector),
#'   `value` (the objective at `argument`), `errpars`, `leverage` (the frame
#'   returned by [remlLeverage]), `dof` (effective degrees of freedom per
#'   observable), `iterations`, `converged` and `fit` (the last mean fit).
#'
#' @seealso [remlLeverage], [normL2]
#' @importFrom stats optim
#' @export
reml <- function(objfun, pars, errpars = NULL, fixed = NULL,
                 iterlim = 25L, tol = 1e-6, rank.tol = 1e-8,
                 optimizer = trust, control = NULL,
                 cores = getOption("dMod.cores", 1L), ...) {

  if (is.null(attr(objfun, "errfn")))
    stop("reml: the objective has no error model to estimate.", call. = FALSE)
  if (is.null(names(pars)))
    stop("reml: pars must be a named numeric vector.", call. = FALSE)

  if (is.null(errpars)) errpars <- .remlErrpars(objfun, pars)
  unknown <- setdiff(errpars, names(pars))
  if (length(unknown))
    stop("reml: errpars not in pars: ", paste(unknown, collapse = ", "),
         call. = FALSE)
  if (!length(errpars))
    stop("reml: no error-model parameter found in pars. Pass errpars.",
         call. = FALSE)

  pars     <- structure(as.numeric(pars), names = names(pars))
  meanpars <- setdiff(names(pars), errpars)
  ctrl     <- .trustControl(list(rinit = 0.1, rmax = 10), control, optimizer,
                            "control")

  fit <- NULL
  converged <- FALSE
  iter <- 0L

  for (iter in seq_len(iterlim)) {

    fit <- do.call(optimizer,
                   c(list(objfun, pars[meanpars],
                          fixed = c(fixed, pars[errpars]), cores = cores),
                     ctrl, list(...)))
    pars[meanpars] <- fit$argument[meanpars]

    lev <- remlLeverage(objfun, pars, meanpars = meanpars, fixed = fixed,
                        rank.tol = rank.tol, cores = cores)

    # The sigma step solves the stationarity condition of the data term. Anything
    # else acting on the error parameters, a prior for instance, is not in it.
    if (iter == 1L) .remlCheckErrGrad(objfun, pars, errpars, lev, fixed, cores)

    phi.old <- pars[errpars]
    pars[errpars] <- .remlUpdateSigma(objfun, pars, errpars, lev$leverage,
                                      fixed = fixed, cores = cores)

    if (max(abs(pars[errpars] - phi.old)) < tol) {
      converged <- TRUE
      break
    }
  }

  lev <- remlLeverage(objfun, pars, meanpars = meanpars, fixed = fixed,
                      rank.tol = rank.tol, cores = cores)

  value.plain <- objfun(pars, fixed = fixed, cores = cores)$value

  list(argument    = pars,
       value       = value.plain + attr(lev, "logdet"),
       value.plain = value.plain,
       logdet      = attr(lev, "logdet"),
       errpars     = errpars,
       leverage    = lev,
       dof         = attr(lev, "dof"),
       rank        = attr(lev, "rank"),
       iterations  = iter,
       converged   = converged,
       fit         = fit)
}


# Minimise sum_i [(1 - h_i) log sigma_i^2 + r_i^2/sigma_i^2] over the error
# parameters at frozen residuals and leverages, which is the REML condition.
.remlUpdateSigma <- function(objfun, pars, errpars, h, fixed = NULL, cores = 1L) {

  cache <- new.env(parent = emptyenv())
  cache$phi <- NULL

  evaluate <- function(phi) {
    if (!is.null(cache$phi) && identical(cache$phi, phi)) return(cache$out)
    p <- pars
    p[errpars] <- phi
    frames <- .remlFrames(objfun, p, fixed = fixed, cores = cores)
    sigma  <- .remlStack(frames, "sigma")
    resid  <- .remlStack(frames, "residual")

    # The search is unconstrained, so it will offer a sigma the error model
    # cannot produce. Refuse the point instead of taking the log of it.
    if (!all(is.finite(sigma)) || any(sigma <= 0) || !all(is.finite(resid))) {
      out <- list(value = Inf, gradient = rep(0, length(errpars)))
      cache$phi <- phi
      cache$out <- out
      return(out)
    }
    dsigma <- .remlBlocks(frames, "deriv.err", errpars)
    missing <- setdiff(errpars, attr(dsigma, "seen"))
    if (length(missing))
      stop("reml: errpars do not reach any error model: ",
           paste(missing, collapse = ", "), call. = FALSE)

    value <- sum((1 - h) * 2 * log(sigma) + (resid / sigma)^2)
    w     <- 2 * (1 - h) / sigma - 2 * resid^2 / sigma^3
    out   <- list(value = value, gradient = colSums(w * dsigma))
    cache$phi <- phi
    cache$out <- out
    out
  }

  opt <- optim(pars[errpars],
               fn     = function(phi) evaluate(phi)$value,
               gr     = function(phi) evaluate(phi)$gradient,
               method = "BFGS")
  structure(opt$par, names = errpars)
}


# Compare the objective's gradient in the error parameters against the one the
# data term alone produces. A difference means a further term acts on them,
# which the sigma step does not see.
.remlCheckErrGrad <- function(objfun, pars, errpars, lev, fixed, cores) {

  frames <- .remlFrames(objfun, pars, fixed = fixed, cores = cores)
  dsigma <- .remlBlocks(frames, "deriv.err", errpars)
  if (!length(attr(dsigma, "seen"))) return(invisible(NULL))

  sigma <- lev$sigma
  resid <- lev$residual
  g.data <- colSums((2 / sigma - 2 * resid^2 / sigma^3) * dsigma)
  g.obj  <- objfun(pars, fixed = fixed, cores = cores)$gradient[errpars]

  scale <- max(1, max(abs(g.obj)))
  if (max(abs(g.obj - g.data)) > 1e-6 * scale)
    warning("reml: the objective acts on ", paste(errpars, collapse = ", "),
            " beyond the data term, a prior for instance. The sigma step solves ",
            "the data term's condition only and ignores that contribution. ",
            "Put such terms on the mean parameters instead.", call. = FALSE)
  invisible(NULL)
}
