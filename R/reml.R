## Restricted maximum likelihood for the error model -------------------------

## Error parameters solve sum_i [1 - h_ii - r_i^2/sigma_i^2] d log sigma_i^2 / d phi_k = 0
## with leverages h_ii. Only first-order sensitivities enter; the fixed point is
## exact although each step freezes h.


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

    prediction <- if (is.list(times))
      term$prdfn(unname(times), pars, fixed = fixed, deriv = TRUE,
                 conditions = names(times), cores = cores) else
      term$prdfn(times, pars, fixed = fixed, deriv = TRUE, cores = cores)
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
    stop("reml: the objective holds no data to evaluate.", call. = FALSE)
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
    stop("objfun does not hold a prediction function and data. ",
         "REML needs an objective built by normL2().", call. = FALSE)
  list(list(data = data, prdfn = x, errfn = attr(objfun, "errfn"),
            timesD = attr(objfun, "timesD")))
}


# Stack the per-row derivative blocks of every frame into one matrix over
# `cols`, filling with zero where a term lacks a parameter.
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


#' Leverage of the Data Points in the Mean Model
#'
#' @description
#' Hat values of the weighted mean model,
#' \eqn{h = \mathrm{diag}(W^{1/2} J (J^\top W J)^{-1} J^\top W^{1/2})} with
#' \eqn{J = \partial\mu/\partial\theta} the prediction sensitivities and
#' \eqn{W = \mathrm{diag}(1/\sigma_i^2)}: the leverage of each data point.
#' Their sum is the numerical rank of the weighted sensitivity matrix, so a
#' non-identifiable direction contributes nothing.
#'
#' @param objfun Objective function built by [normL2()], with prediction
#'   function, data and error model.
#' @param pars Named numeric parameter vector, usually a fit.
#' @param meanpars Character, the mean parameters spanning \eqn{J}. Default
#'   `NULL`: `names(pars)` without the error-model parameters.
#' @param fixed Named numeric, passed to the prediction. Default `NULL`.
#' @param rankTol Relative threshold on the singular values of the weighted
#'   sensitivity matrix below which a direction counts as absent. Default
#'   `1e-8`.
#' @param cores Passed to the prediction. Default `getOption("dMod2.cores", 1)`.
#'
#' @return A data frame with one row per data point and columns `condition`,
#'   `time`, `name`, `sigma`, `residual` and `leverage`. Attributes: `"rank"`,
#'   the effective number of mean parameters; `"logdet"`, the log determinant
#'   of \eqn{J^\top W J} over the directions present; `"dof"`, per observable
#'   the number of data points minus the summed leverage.
#'
#' @seealso [reml()]
#' @inherit reml examples
#' @importFrom stats setNames
#' @export
remlLeverage <- function(objfun, pars, meanpars = NULL, fixed = NULL,
                         rankTol = 1e-8,
                         cores = getOption("dMod2.cores", 1L)) {

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
  k <- sum(s$d > rankTol * max(s$d))
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


#' Estimate the Error Model by Restricted Maximum Likelihood
#'
#' @description
#' Alternates between fitting the mean parameters at a fixed error model and
#' updating the error parameters from the REML stationarity condition, until
#' the error parameters stop moving. Every data point enters with the leverage
#' of that point, see [remlLeverage()]. Every L2 term of a composed objective is
#' used, so `normL2(d1, ...) + normL2(d2, ...)` is handled as one dataset.
#'
#' @details
#' Only first-order sensitivities are used. A term acting on the error
#' parameters from outside the data term, such as a prior, is not part of the
#' REML condition; `reml()` warns when it finds one.
#'
#' @param objfun Objective function built by [normL2()], with an error model.
#' @param pars Named numeric starting vector of mean and error parameters.
#' @param errpars Character, the error-model parameters. Default `NULL`: the
#'   symbols of the error equations that appear in `pars`.
#' @param fixed Named numeric, held fixed throughout. Default `NULL`.
#' @param iterlim Maximum number of rounds. Default `25`.
#' @param tol Convergence threshold on the largest change of an error parameter
#'   between two rounds. Default `1e-6`.
#' @param rankTol Passed to [remlLeverage()]. Default `1e-8`.
#' @param optimizer Optimizer for the mean step, called as
#'   `optimizer(objfun, parinit, fixed = , cores = , ...)`. Default [trust()].
#' @param control Named list of arguments for `optimizer`, merged into
#'   `list(rinit = 0.1, rmax = 10)`. Default `NULL`.
#' @param cores Passed to the prediction and the objective. Default
#'   `getOption("dMod2.cores", 1)`.
#' @param ... Further arguments for `objfun` in the mean step, passed through
#'   `optimizer`. `value.plain` is evaluated without them.
#'
#' @return A list with components
#'   \describe{
#'     \item{`argument`}{The full parameter vector.}
#'     \item{`value`}{The REML criterion, `value.plain + logdet`.}
#'     \item{`value.plain`}{The objective at `argument`.}
#'     \item{`logdet`}{Log determinant of \eqn{J^\top W J}, see
#'       [remlLeverage()].}
#'     \item{`errpars`}{The error-model parameters.}
#'     \item{`leverage`}{The data frame returned by [remlLeverage()].}
#'     \item{`dof`}{Effective degrees of freedom per observable.}
#'     \item{`rank`}{Effective number of mean parameters.}
#'     \item{`iterations`}{Number of rounds.}
#'     \item{`converged`}{Logical.}
#'     \item{`fit`}{The last mean fit.}
#'   }
#'
#' @seealso [remlLeverage()], [normL2()]
#'
#' @examples
#' \donttest{
#' old <- options(dMod2.outdir = tempdir())
#' f <- addReaction(eqnlist(), from = "A", to = "", rate = "k*A")
#' x <- Xs(odemodel(f, modelname = "reml_x", compile = FALSE))
#' g <- Y(c(y = "A"), f = x, attachInput = FALSE, modelname = "reml_g",
#'        compile = FALSE)
#' e <- Y(c(y = "sigma_y"), f = g, attachInput = FALSE, condition = "C1",
#'        modelname = "reml_e", compile = FALSE)
#' p <- P(eqnvec(A = "exp(logA)", k = "exp(logk)", sigma_y = "sigma_y"),
#'        condition = "C1", modelname = "reml_p", compile = FALSE)
#' compile(x, g, e, p, output = "reml_example", cores = 1)
#'
#' set.seed(1)
#' times <- seq(0, 5, by = 0.25)
#' data <- datalist(C1 = data.frame(
#'   name = "y", time = times, sigma = NA_real_,
#'   value = 2 * exp(-0.5 * times) + rnorm(length(times), sd = 0.1)))
#' obj <- normL2(data, g * x * p, errmodel = e)
#'
#' fit <- reml(obj, c(logA = 0, logk = 0, sigma_y = 1))
#' fit$argument
#' fit$dof
#' head(remlLeverage(obj, fit$argument))
#' options(old)
#' }
#' @importFrom stats optim
#' @export
reml <- function(objfun, pars, errpars = NULL, fixed = NULL,
                 iterlim = 25L, tol = 1e-6, rankTol = 1e-8,
                 optimizer = trust, control = NULL,
                 cores = getOption("dMod2.cores", 1L), ...) {

  dots <- list(...)
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
                     ctrl, dots))
    pars[meanpars] <- fit$argument[meanpars]

    lev <- remlLeverage(objfun, pars, meanpars = meanpars, fixed = fixed,
                        rankTol = rankTol, cores = cores)

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
                      rankTol = rankTol, cores = cores)

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


# Minimize sum_i [(1 - h_i) log sigma_i^2 + r_i^2/sigma_i^2] over the error
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


## Restricted likelihood of mixed-effects and alignment fits ---------------------------------

## Laplace over `mean` under a flat prior, log|S / 2| - p log(2 pi), S the Schur complement over the etas of the
## information G = sum_r w_r a_r a_r' + C + P at `full`; with `deriv` the gradient in `theta`, sum_r q_r dw_r +
## tr(Gamma dP), q_r = a_r' Gamma a_r, Gamma = G^-1 minus the inverse eta blocks.
.remlTerm <- function(design, full, fixed, mean, theta = character(0), deriv = TRUE, cores = 1L) {
  units <- if (length(.allEtas(design))) .volumeUnits(design, full, fixed) else list()
  V <- c(mean, unlist(lapply(units, `[[`, "eta")))
  conds <- names(design$data)
  preds <- .predictMany(design$prdfn, times = lapply(conds, function(cn) .emGrid(design, cn)),
                        parsList = rep(list(full), length(conds)), conditions = conds, fixed = fixed,
                        deriv = TRUE, cores = cores)
  names(preds) <- conds
  frames <- .volumeFrameList(design, conds, preds, fixed, cores, deriv2 = FALSE)
  rows <- lapply(conds, function(cn) {
    A <- .frameParameters(frames[[cn]])
    rw <- .remlRows(design, cn, frames[[cn]], intersect(V, A), intersect(theta, A))
    lr <- .rowLossDerivs(rw$p, rw$s, rw$y, rw$cens, design$df)
    rw$w <- lr[, 3L, 1L]; rw$dw <- lr[, 3L, 2L]
    # Student-t readings by their expected curvature, which is positive
    if (is.finite(design$df)) {
      k <- !rw$cens
      rw$w[k] <- 2 * (design$df + 1) / ((design$df + 3) * rw$s[k]^2)
      rw$dw[k] <- -2 * rw$w[k] / rw$s[k]
    }
    rw
  })
  # the information over the mean parameters and the etas, sparse
  idx <- function(v) match(v, V)
  trip <- lapply(rows, function(rw) {
    if (!length(rw$V)) return(NULL)
    a <- matrix(rw$a, length(rw$w))
    B <- crossprod(a, rw$w * a)
    k <- idx(rw$V)
    list(i = rep(k, length(k)), j = rep(k, each = length(k)), x = as.vector(B))
  })
  for (u in units) {
    k <- idx(u$eta)
    trip[[length(trip) + 1L]] <- list(i = rep(k, length(k)), j = rep(k, each = length(k)), x = as.vector(u$P))
  }
  G <- Matrix::sparseMatrix(i = unlist(lapply(trip, `[[`, "i")), j = unlist(lapply(trip, `[[`, "j")),
                            x = unlist(lapply(trip, `[[`, "x")), dims = c(length(V), length(V)),
                            dimnames = list(V, V))
  others <- Filter(function(t) is.null(attr(t, "l2spec", exact = TRUE)), .objTerms(design$objfn))
  lr <- NULL
  if (length(others)) {
    oc <- Reduce(`+`, others)(pars = full, fixed = fixed, deriv = TRUE, hessian = "sparse")
    G <- G + .sparseOver(oc$hessian, V)
    lr <- .lowrankOver(oc$lowrank, V)
  }
  G <- Matrix::forceSymmetric(G)
  # G is the sparse part plus the low-rank one; a sparse part that is singular alone is factorized with it
  f <- .lowrankCholesky(G, lr)
  if (is.null(f) && !is.null(lr)) {
    G <- Matrix::forceSymmetric(.lowrankDense(G, lr))
    lr <- NULL
    f <- .lowrankCholesky(G)
  }
  if (is.null(f)) return(list(value = Inf))
  D <- lapply(units, function(u) {
    Ge <- as.matrix(G[u$eta, u$eta, drop = FALSE])
    if (is.null(lr)) Ge else Ge + lr$U[u$eta, , drop = FALSE] %*% (lr$d * t(lr$U[u$eta, , drop = FALSE]))
  })
  value <- f$logdet - sum(vapply(D, function(d) as.numeric(determinant(d)$modulus), 0)) - length(mean) * log(4 * pi)
  if (!deriv) return(list(value = value))
  # q_r and the eta blocks of G^-1 from L^-1 P on the rows and the unit vectors of the etas
  n <- vapply(rows, function(rw) length(rw$w), 0L)
  rhs <- Matrix::sparseMatrix(
    i = unlist(lapply(rows, function(rw) rep(idx(rw$V), each = length(rw$w)))),
    j = unlist(Map(function(rw, off) rep(off + seq_along(rw$w), length(rw$V)), rows, cumsum(c(0L, n))[seq_along(n)])),
    x = unlist(lapply(rows, function(rw) as.vector(rw$a))), dims = c(length(V), sum(n)))
  E <- Matrix::sparseMatrix(i = idx(V[-seq_along(mean)]), j = seq_along(V[-seq_along(mean)]), x = 1,
                            dims = c(length(V), length(V) - length(mean)))
  X <- .lowrankZ(f, cbind(rhs, E))
  q <- .lowrankQuadratic(f, X[, seq_len(sum(n)), drop = FALSE])
  Xe <- X[, sum(n) + seq_len(ncol(E)), drop = FALSE]
  grad <- setNames(numeric(length(theta)), theta)
  off <- 0L
  for (r in seq_along(rows)) {
    rw <- rows[[r]]
    at <- off + seq_along(rw$w)
    off <- off + length(rw$w)
    if (!length(rw$A)) next
    qr <- q[at]
    # the eta part of q_r, which the subtraction of log|G_ee| takes out
    u <- which(vapply(units, function(un) conds[r] %in% un$conditions, TRUE))
    if (length(u)) {
      e <- units[[u[1]]]$eta
      ae <- matrix(rw$a, length(rw$w))[, match(e, rw$V), drop = FALSE]
      ae[is.na(ae)] <- 0
      qr <- qr - rowSums((ae %*% solve(D[[u[1]]])) * ae)
    }
    grad[rw$A] <- grad[rw$A] + drop(crossprod(matrix(rw$sA, length(rw$w)), qr * rw$dw))
  }
  for (k in seq_along(units)) {
    cps <- intersect(names(units[[k]]$dP), theta)
    if (!length(cps)) next
    cols <- match(units[[k]]$eta, V[-seq_along(mean)])
    Ginv <- .lowrankQuadratic(f, Xe[, cols, drop = FALSE], diagonal = FALSE) - solve(D[[k]])
    for (cp in cps) grad[cp] <- grad[cp] + sum(Ginv * units[[k]]$dP[[cp]])
  }
  list(value = value, gradient = grad)
}

## The data rows of condition `cn` for .remlTerm(): prediction p, sigma s, data y, censoring, the
## derivatives `a` of p in the parameters `V` and `sA` of s in `theta`; a recorded sigma has none.
.remlRows <- function(design, cn, frames, V, theta) {
  d <- .conditionData(design, cn)
  pred <- frames$pred
  tp <- .matchNum(d$time, as.numeric(pred[, "time"]))
  ip <- tp + (match(d$name, dimnames(attr(pred, "deriv"))[[2]]) - 1L) * nrow(pred)
  rows <- list(p = pred[cbind(tp, match(d$name, colnames(pred)))], s = d$sigma, y = pmax(d$value, d$lloq),
               cens = d$value <= d$lloq, a = .rowSlices(attr(pred, "deriv"), ip, V),
               sA = .rowSlices(NULL, seq_len(nrow(d)), theta), V = V, A = theta)
  err <- frames$err
  if (is.null(err)) return(rows)
  te <- .matchNum(d$time, as.numeric(err[, "time"]))
  model <- is.na(d$sigma)
  rows$s[model] <- err[cbind(te, match(d$name, colnames(err)))][model]
  e1 <- attr(err, "deriv")
  if (is.null(e1)) return(rows)
  ie <- te + (match(d$name, dimnames(e1)[[2]]) - 1L) * nrow(err)
  none <- is.na(ie) | !model
  ie[none] <- 1L
  rows$sA <- .rowSlices(e1, ie, theta)
  rows$sA[none, ] <- 0
  rows
}

## The restricted likelihood term as an objective of the error parameters and hyperparameters
## `theta`, the other parameters of `full` and the modes held at the point `anchor` of the previous
## round, the linearization of Lindstrom and Bates. The Hessian of the term is left 0.
.remlAnchor <- function(design, anchor, mean, theta) {
  last <- NULL
  fn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, hessian = NULL, conditions = NULL,
                 env = NULL, cores = getOption("dMod2.cores", 1L)) {
    arglist <- list(...)
    pars <- arglist[match.fnargs(arglist, "pars")][[1]]
    th <- intersect(theta, names(pars))
    full <- anchor
    full[th] <- pars[th]
    key <- list(full[th], isTRUE(deriv))
    if (!identical(last$key, key))
      last <<- list(key = key, rt = .remlTerm(design, full, fixed, mean, th, deriv = isTRUE(deriv), cores = cores))
    rt <- last$rt
    if (!isTRUE(deriv) || !is.finite(rt$value)) return(objlist(value = rt$value, gradient = NULL, hessian = NULL))
    gradient <- setNames(numeric(length(pars)), names(pars))
    gradient[th] <- rt$gradient[th]
    cv <- .resolveCurvature(deriv, deriv2, hessian, "forward")
    H <- if (cv$hessian) {
      if (cv$sparse) .sparseOver(NULL, names(pars))
      else matrix(0, length(pars), length(pars), dimnames = list(names(pars), names(pars)))
    }
    objlist(value = rt$value, gradient = gradient, hessian = H)
  }
  class(fn) <- c("objfn", "fn")
  attr(fn, "parameters") <- theta
  fn
}

## The pieces of a design that .remlTerm() reads, for an objective without random effects.
.remlDesign <- function(obj) {
  spec <- attr(obj, "l2spec", exact = TRUE)[[1]]
  list(prdfn = attr(obj, "prdfn", exact = TRUE), errfn = attr(obj, "errfn", exact = TRUE),
       data = attr(obj, "data", exact = TRUE), objfn = obj, df = spec$df %||% Inf)
}
