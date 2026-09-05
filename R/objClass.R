## Methods for class "objfn" -----------------------------------------------



## Class "objlist" and its constructors ------------------------------------



#' Generate objective list from numeric vector
#' 
#' @param p Named numeric vector
#' @return list with entries value (\code{0}), 
#' gradient (\code{rep(0, length(p))}) and 
#' hessian (\code{matrix(0, length(p), length(p))}) of class \code{obj}.
#' @examples
#' p <- c(A = 1, B = 2)
#' as.objlist(p)
#' @export
as.objlist <- function(p) {
  
  objlist(value = 0,
          gradient = structure(rep(0, length(p)), names = names(p)),
          hessian = matrix(0, length(p), length(p), dimnames = list(names(p), names(p))))
  
}


#' Per-condition residual contribution to an L2 objective
#'
#' @description
#' Computes the negative-log-likelihood residual contribution of a single
#' condition (with optional error model). Exposed so quadrature node-loops
#' can evaluate one condition without paying the per-call cost of
#' [normL2]'s multi-condition setup.
#'
#' @param dataI datalist entry for one condition (data.frame with
#'   `name`, `time`, `value`, `sigma` columns).
#' @param predictionI prdframe for that condition (typically `prediction[[cn]]`
#'   from a prdfn call).
#' @param pars Named numeric parameter vector at which to evaluate.
#' @param errfn Optional obsfn defining a parameter-dependent error model.
#' @param fixed Optional fixed-parameter vector (passed through to `errfn`).
#' @param cn Character condition name. Required when `errfn` is set, used
#'   for errmodel condition routing.
#' @param eCondNames Optional character vector of condition names that have
#'   an errmodel mapping. NULL means `errfn` applies to all.
#' @param deriv,deriv2 Logical. Whether to return gradient/Hessian.
#' @param opt.BLOQ Character. BLOQ likelihood treatment (see [normL2]).
#'   One of `"M1"`, `"M3"` (default), `"M4NM"`, `"M4BEAL"`.
#'
#' @return An [objlist] for the single condition's contribution.
#' @export
evalConditionResidual <- function(dataI, predictionI, pars,
                                  errfn      = NULL,
                                  fixed      = NULL,
                                  cn         = NULL,
                                  eCondNames = NULL,
                                  deriv      = TRUE,
                                  deriv2     = FALSE,
                                  opt.BLOQ   = c("M3", "M1", "M4NM", "M4BEAL")) {
  opt.BLOQ <- match.arg(opt.BLOQ)
  err_cn <- NULL
  if (!is.null(errfn) && (is.null(eCondNames) || cn %in% eCondNames)) {
    if (is.null(cn))
      stop("evalConditionResidual: `cn` must be supplied when `errfn` is set.")
    pinner     <- getParameters(predictionI)
    fixedinner <- pinner[attr(pinner, "fixed")]
    pinner     <- as.parvec(pinner[setdiff(names(pinner), names(fixed))])
    fixedinner <- as.parvec(fixedinner, deriv = FALSE, deriv2 = FALSE)
    err_cn <- errfn(out = predictionI, pars = pinner,
                    fixed = fixedinner, conditions = cn)[[cn]]
  }

  key   <- if (is.null(cn)) "cond" else cn
  pred1 <- setNames(list(predictionI), key)
  err1  <- if (!is.null(err_cn)) setNames(list(err_cn), key) else NULL
  meta1 <- .build_normL2_meta(setNames(list(dataI), key), pred1, err1,
                              key, eCondNames)

  d_dn <- dimnames(attr(predictionI, "deriv"))
  par_names_global <- if (!is.null(d_dn)) d_dn[[3]] else character(0)

  kr <- normL2_kernel(
    prediction       = pred1,
    err_list_opt     = err1,
    meta_list        = meta1,
    par_names_global = par_names_global,
    deriv2_requested = isTRUE(deriv2),
    threads          = 1L,
    bloq_mode        = opt.BLOQ
  )
  if (deriv)
    objlist(value = kr$value, gradient = kr$gradient, hessian = kr$hessian)
  else
    objlist(value = kr$value, gradient = NULL, hessian = NULL)
}



# Internal: gradient and Hessian in the order of the parameter vector the
# objective was called with. The kernel names them after the union of the
# per-condition sensitivity blocks, while `trust()` and `optim()` read them
# positionally against their own start vector. A parameter the objective does
# not depend on gets a zero derivative rather than being left out.
.alignObjlist <- function(out, pnames) {
  g <- out$gradient
  if (is.null(g) || !length(pnames) || identical(names(g), pnames)) return(out)
  hit <- intersect(pnames, names(g))
  gradient <- setNames(numeric(length(pnames)), pnames)
  gradient[hit] <- g[hit]
  out$gradient <- gradient
  # A skipped Hessian stays NULL; only realign one that was actually built.
  if (!is.null(out$hessian)) {
    hessian <- matrix(0, length(pnames), length(pnames),
                      dimnames = list(pnames, pnames))
    if (length(hit)) hessian[hit, hit] <- out$hessian[hit, hit]
    out$hessian <- hessian
  }
  out
}

#' L2 norm between data and model prediction
#'
#' @description
#' Creates an objective function for parameter estimation based on the
#' (negative log-likelihood) L2 norm between observed data and model predictions.
#' The returned objective function can be used with optimizers such as
#' [mstrust] and supports aggregation over multiple experimental conditions.
#'
#' @param data Object of class [datalist].
#' @param x Object of class [prdfn].
#' @param errmodel Optional object of class [obsfn]. The error model may be
#'   defined only for a subset of conditions.
#' @param times Optional numeric vector of additional time points at which the
#'   prediction function is evaluated. If NULL, time points are taken from the
#'   data. Event times should be included here if the prediction model uses events.
#' @param t0 Numeric. Start of the time grid, i.e. the time at which initial
#'   values take effect. Defaults to 0.
#' @param attr.name Character string. The objective value is additionally
#'   returned as an attribute of this name, and the sum of squares behind it
#'   under `chi2`. Adding objectives pools the terms sharing an `attr.name`,
#'   so `chi2` stays one number; where two `attr.name`s meet it splits into
#'   `chi2_<attr.name>` per contribution.
#' @param cores Deprecated and ignored. Pass `cores` to the objective call, or
#'   set `options(dMod.cores = )`.
#' @param opt.BLOQ Character. NONMEM-style treatment of below-LOQ rows
#'   (those with `value <= lloq` in the data). One of `"M1"` (drop BLOQ rows
#'   from the objective), `"M3"` (censored log-likelihood, default), `"M4NM"`
#'   or `"M4BEAL"` (truncated variants; require non-negative LOQ).
#'
#' @return
#' An object of class `objfn`, i.e. a function
#' \code{obj(pars, fixed, deriv, env, cores)} returning an [objlist].
#'
#' @details
#' Combine objectives with `+` (see [sumobjfn]). `cores` is a call-time
#' argument: it sets the thread count for both the batched ODE integration and
#' the C++ residual kernel. It defaults to `getOption("dMod.cores", 1L)`.
#'
#' @example inst/examples/normL2.R
#' @export
normL2 <- function(data, x, errmodel = NULL, times = NULL, t0 = 0,
                   attr.name = "data",
                   cores = 1L,
                   opt.BLOQ = c("M3", "M1", "M4NM", "M4BEAL")) {

  if (!missing(cores))
    warning("normL2: 'cores' at construction time is deprecated and ignored. ",
            "Pass cores = to the objective call, or set options(dMod.cores = ).",
            call. = FALSE)
  opt.BLOQ <- match.arg(opt.BLOQ)

  # `t0` anchors the time grid: the prediction starts there, so that is where
  # initial values take effect.
  timesD <- sort(unique(c(t0, unlist(lapply(data, `[[`, "time")), times)))

  x.cond <- names(attr(x, "mappings"))
  d.cond <- names(data)
  stopifnot(all(d.cond %in% x.cond))

  e.cond <- if (!is.null(errmodel)) names(attr(errmodel, "mappings")) else NULL
  conditions.obj <- intersect(x.cond, d.cond)

  # Force early binding
  force(errmodel); force(conditions.obj); force(timesD)

  # Lazy meta cache for the C++ kernel path. Built on first call; rebuilt
  # if the deriv column set changes (e.g. when `fixed` toggles between
  # calls - uncommon, but cheap to detect via length+name compare).
  .meta_cache <- new.env(parent = emptyenv())
  .meta_cache$meta_list        <- NULL
  .meta_cache$par_names_global <- NULL
  .meta_cache$signature        <- NULL  # used to invalidate on shape change

  # `.prediction` lets a caller that already batched the predictions hand them
  # in; see .objEvalMany().
  myfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, hessian = TRUE,
                   conditions = NULL, env = NULL,
                   cores = getOption("dMod.cores", 1L), .prediction = NULL) {
    pars <- ..1
    if (is.null(env)) env <- new.env()
    conditions <- if (is.null(conditions)) conditions.obj else
      intersect(conditions.obj, conditions)
    if (!length(conditions)) return(NULL)

    # The Hessian is only meaningful with deriv; when it is not wanted, the
    # J^T J contraction is skipped and no second-order sensitivities are needed.
    build_hessian <- isTRUE(deriv) && isTRUE(hessian)
    deriv2 <- isTRUE(deriv2) && build_hessian

    prediction <- if (!is.null(.prediction)) .prediction else
      x(times = timesD, pars = pars, fixed = fixed,
        deriv = deriv, deriv2 = deriv2, conditions = conditions,
        cores = cores)
    if (!is.null(.prediction) && !identical(names(prediction), conditions))
      prediction <- prediction[conditions]

    # Build errmodel output per condition (if any). One batched evaluation
    # rather than one public obsfn call per condition: the shim, the bundle
    # and the prdlist wrapping were paid 32 times for 32 scalar kernel calls,
    # and the leaf's batch entry never saw more than one request.
    err_list <- NULL
    if (!is.null(errmodel)) {
      cn_eval <- if (is.null(e.cond)) conditions else intersect(conditions, e.cond)
      split <- lapply(cn_eval, function(cn) {
        pinner     <- getParameters(prediction[[cn]])
        # The `fixed` marker is derived from the sensitivity rows, so it is
        # empty under deriv = FALSE. Fall back to the outer fixed names, or
        # the error model loses those parameters entirely.
        fixedinner <- pinner[union(attr(pinner, "fixed"),
                                   intersect(names(pinner), names(fixed)))]
        list(pars  = as.parvec(pinner[setdiff(names(pinner), names(fixed))]),
             fixed = as.parvec(fixedinner, deriv = FALSE, deriv2 = FALSE))
      })
      est <- .fnNode(errmodel)
      got <- if (!is.null(est) && length(cn_eval) > 1L) {
        b <- .bundle(conds = cn_eval,
                     out   = lapply(cn_eval, function(cn) prediction[[cn]]),
                     pars  = lapply(split, `[[`, "pars"),
                     fixed = lapply(split, `[[`, "fixed"),
                     shared = FALSE)
        .evalNode(est, b, deriv, deriv2, NULL, cores)
      } else {
        lapply(seq_along(cn_eval), function(j)
          errmodel(out = prediction[[cn_eval[j]]], pars = split[[j]]$pars,
                   fixed = split[[j]]$fixed, conditions = cn_eval[j])[[cn_eval[j]]])
      }
      # keep the NULL holes: .build_normL2_meta and the kernel index positionally
      err_list <- vector("list", length(conditions))
      err_list[match(cn_eval, conditions)] <- got
    }

    # Determine current deriv signature (per-condition local par names);
    # empty for value-only (deriv = FALSE) evaluations.
    cur_sig <- lapply(prediction, function(pr) dimnames(attr(pr, "deriv"))[[3]])
    if (is.null(.meta_cache$meta_list) ||
        !identical(.meta_cache$signature, cur_sig)) {
      .meta_cache$par_names_global <- unique(unlist(cur_sig))
      .meta_cache$meta_list <- .build_normL2_meta(
        data, prediction, err_list, conditions, e.cond)
      .meta_cache$signature <- cur_sig
    }
    par_names_global <- .meta_cache$par_names_global
    if (is.null(par_names_global)) par_names_global <- character(0)

    kr <- normL2_kernel(
      prediction       = prediction,
      err_list_opt     = err_list,
      meta_list        = .meta_cache$meta_list,
      par_names_global = par_names_global,
      deriv2_requested = isTRUE(deriv2),
      threads          = as.integer(cores),
      bloq_mode        = opt.BLOQ,
      build_hessian    = build_hessian
    )
    out <- if (deriv)
      .alignObjlist(objlist(value = kr$value, gradient = kr$gradient,
                            hessian = kr$hessian), names(pars))
    else
      objlist(value = kr$value, gradient = NULL, hessian = NULL)
    attr(out, attr.name) <- out$value
    # The sum of squares alone, tagged with the contribution it belongs to.
    # Summing objectives pools terms sharing an `attr.name` into one `chi2` and
    # splits the rest into `chi2_<attr.name>`.
    attr(out, "chi2") <- setNames(kr$chi2, attr.name)
    env$prediction <- prediction
    attr(out, "env") <- env
    out
  }

  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- d.cond
  # Union of prediction-fn and errmodel parameters so the errmodel's sigma
  # parameters survive when the inner solver reads `full_pars`.
  err_pars <- if (!is.null(errmodel)) attr(errmodel, "parameters") else character(0)
  attr(myfn, "parameters") <- union(attr(x, "parameters"), err_pars)
  attr(myfn, "modelname") <- modelname(x, errmodel)
  # Carried like `*` and `+` carry it, so loadDLL() and compile() can reach the
  # shared objects of a composed objective.
  attr(myfn, "compileInfo") <- .mergeCompileInfo(attr(x, "compileInfo"),
                                                 attr(errmodel, "compileInfo"))
  # Reconstruction handles: recover the model pieces from a composed objective
  # instead of re-demanding them as arguments. Setting an attribute to NULL is
  # a no-op, so "errfn" is simply absent when there is no error model.
  # One entry per L2 term. The single handles below are first-wins on
  # composition, so a split objective would otherwise expose only its first
  # term, and reml() needs every one of them.
  attr(myfn, "l2spec") <- list(list(data = data, prdfn = x, errfn = errmodel,
                                    timesD = timesD))
  attr(myfn, "prdfn") <- x
  attr(myfn, "data")  <- data
  attr(myfn, "timesD") <- timesD
  attr(myfn, "errfn") <- errmodel
  myfn
}


# Evaluate a set of single-condition objectives, one parameter vector each, with
# the predictions gathered into one batched request. A loop over objectives
# built on a shared prdfn is one ODE solve per objective per pass; the
# condition axis is exactly what the batch parallelises.
#
# Falls back to the loop whenever the batch cannot reproduce it: a missing
# reconstruction handle, objectives from different prdfns, a condition owned by
# more than one objective, or a failing batch. The loop is the reference, so the
# fallback is always correct, only slower.
.objEvalMany <- function(objList, parsList, deriv = TRUE, deriv2 = FALSE,
                         cores = getOption("dMod.cores", 1L)) {
  n <- length(parsList)
  objs <- if (is.function(objList)) rep(list(objList), n) else objList
  stopifnot(length(objs) == n)

  serial <- function() lapply(seq_len(n), function(j)
    objs[[j]](parsList[[j]], deriv = deriv, deriv2 = deriv2, cores = cores))
  if (n < 2L) return(serial())

  prd   <- attr(objs[[1L]], "prdfn", exact = TRUE)
  conds <- vapply(objs, function(o) {
    cn <- attr(o, "conditions")
    if (length(cn) == 1L) cn else NA_character_
  }, "")
  times <- lapply(objs, attr, "timesD", exact = TRUE)
  if (is.null(prd) || is.null(.fnNode(prd)) ||
      anyNA(conds) || anyDuplicated(conds) > 0L ||
      any(vapply(times, is.null, TRUE)) ||
      !all(vapply(objs, function(o)
        identical(attr(o, "prdfn", exact = TRUE), prd) &&
        any(c(".prediction", "...") %in% names(formals(o))), TRUE)))
    return(serial())

  preds <- tryCatch(
    .predictMany(prd, times = times, parsList = parsList, conditions = conds,
                 deriv = deriv, deriv2 = deriv2, cores = cores),
    error = function(e) NULL)
  if (is.null(preds)) return(serial())

  lapply(seq_len(n), function(j)
    objs[[j]](parsList[[j]], deriv = deriv, deriv2 = deriv2,
              .prediction = setNames(preds[j], conds[j])))
}


# Build per-condition metadata for the C++ normL2_kernel. Indexes data rows
# into prediction/errmodel matrices, encodes ALOQ/BLOQ partition, and stores
# the LOQ-substituted y values (matching res()'s `pmax(value, lloq)`).
.build_normL2_meta <- function(data, prediction, err_list, conditions, e_cond) {
  err_list_named <- if (!is.null(err_list)) {
    setNames(err_list, conditions)
  } else {
    NULL
  }
  lapply(conditions, function(cn) {
    dataI <- data[[cn]]
    dataI$name <- as.character(dataI$name)
    prdfI <- prediction[[cn]]
    pcols <- colnames(prdfI)
    d_dn  <- dimnames(attr(prdfI, "deriv"))

    has_deriv <- !is.null(d_dn)
    t_idx_in_pred  <- match(dataI$time, prdfI[, "time"])
    o_idx_in_pred  <- match(dataI$name, pcols)
    o_idx_in_deriv <- if (has_deriv) match(dataI$name, d_dn[[2]])
                      else rep(0L, nrow(dataI))

    if (anyNA(t_idx_in_pred) || anyNA(o_idx_in_pred) ||
        (has_deriv && anyNA(o_idx_in_deriv))) {
      stop(".build_normL2_meta: data point not found in prediction for condition '",
           cn, "'.", call. = FALSE)
    }

    sig <- if (!is.null(dataI$sigma)) dataI$sigma else rep(NA_real_, nrow(dataI))
    sigma_is_na <- is.na(sig)
    sigma_fixed <- ifelse(sigma_is_na, 0, sig)

    t_idx_in_err <- rep(0L, nrow(dataI))
    o_idx_in_err <- rep(0L, nrow(dataI))
    o_idx_in_err_deriv <- rep(0L, nrow(dataI))
    if (any(sigma_is_na) && !is.null(err_list_named)) {
      erm <- err_list_named[[cn]]
      if (!is.null(erm)) {
        t_idx_in_err <- match(dataI$time, erm[, "time"])
        o_idx_in_err <- match(dataI$name, colnames(erm))
        e_dn <- dimnames(attr(erm, "deriv"))
        if (!is.null(e_dn)) {
          o_idx_in_err_deriv <- match(dataI$name, e_dn[[2]])
          o_idx_in_err_deriv[is.na(o_idx_in_err_deriv)] <- 0L
        }
      }
    }

    lloq <- if (!is.null(dataI$lloq)) dataI$lloq else rep(-Inf, nrow(dataI))
    val  <- pmax(dataI$value, lloq)
    bloq_mask <- as.integer(val <= lloq)

    list(
      t_idx_in_pred       = as.integer(t_idx_in_pred),
      o_idx_in_pred       = as.integer(o_idx_in_pred),
      o_idx_in_deriv      = as.integer(o_idx_in_deriv),
      t_idx_in_err        = as.integer(t_idx_in_err),
      o_idx_in_err        = as.integer(o_idx_in_err),
      o_idx_in_err_deriv  = as.integer(o_idx_in_err_deriv),
      sigma_is_na         = as.integer(sigma_is_na),
      sigma_fixed         = as.numeric(sigma_fixed),
      y_data              = as.numeric(val),
      lloq                = as.numeric(lloq),
      bloq_mask           = bloq_mask
    )
  })
}



#' Soft L2 constraint on parameters
#'
#' @param mu Named numeric vector of prior means. Its names select the
#'   constrained parameters.
#' @param sigma Named numeric or character vector. Character entries indicate
#'   log-scale sigma parameters to be estimated.
#' @param attr.name Character. Name of the attribute storing the constraint value.
#' @param condition Optional character vector of conditions.
#'
#' @details
#' Computes, depending on which path is selected,
#' \deqn{(p-\mu)^2 / \sigma^2}
#' or, if sigma is estimated,
#' \deqn{(p-\mu)^2 / \sigma^2 + 2\log(\sigma)},
#' with sigma internally transformed via \code{exp()}. This is the penalty
#' form and drops the Gaussian's normalisation; [constraintL1] and its siblings
#' are the full `-2 log` density.
#'
#' @return Object of class \code{objfn}.
#' @export
constraintL2 <- function(mu, sigma = 1, attr.name = "prior", condition = NULL) {

  est <- is.character(sigma)
  if (length(sigma) == 1) sigma <- setNames(rep(sigma, length(mu)), names(mu))
  if (is.null(names(sigma))) names(sigma) <- names(mu)
  sigma <- sigma[names(mu)]

  myfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, hessian = TRUE,
                   conditions = condition, env = NULL,
                   cores = getOption("dMod.cores", 1L)) {

    p <- list(...)[[match.fnargs(list(...), "pars")]]
    build_hessian <- isTRUE(deriv) && isTRUE(hessian)
    dP  <- if (deriv) attr(p, "deriv", exact = TRUE) else NULL
    dP2 <- if (build_hessian && deriv2) attr(p, "deriv2", exact = TRUE) else NULL

    sigma_pars <- if (est) sigma[names(mu)] else rep("", length(mu))
    sigma_vec  <- if (est) rep(0.0, length(mu)) else as.numeric(sigma[names(mu)])
    kr <- constraintL2_scalar_kernel(
      pars = p,
      dP_opt = if (!is.null(dP)) dP else NULL,
      dP2_opt = if (!is.null(dP2)) dP2 else NULL,
      inner_par_names = names(p),
      fixed_opt = fixed,
      mu_names = names(mu),
      mu = as.numeric(mu),
      sigma = sigma_vec,
      sigma_pars = as.character(sigma_pars),
      est = est,
      deriv = deriv,
      build_hessian = build_hessian
    )

    out <- objlist(value = kr$value, gradient = kr$gradient, hessian = kr$hessian)
    attr(out, attr.name) <- out$value
    attr(out, "env") <- env
    out
  }

  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- condition
  attr(myfn, "parameters") <- names(mu)
  myfn
}



# Internal: broadcast a constraint's shape argument onto the constrained
# parameters and reject the estimated-parameter spelling that only
# constraintL2() supports.
.constraintArg <- function(x, parnames, what, fn) {
  if (is.character(x))
    stop("`", what, "` must be numeric in ", fn, "(). An estimated scale is ",
         "only available in constraintL2().")
  if (anyDuplicated(parnames))
    stop("Duplicated parameter name in ", fn, "(): ",
         paste(unique(parnames[duplicated(parnames)]), collapse = ", "), ".")
  if (length(x) == 1L) x <- setNames(rep(x, length(parnames)), parnames)
  if (is.null(names(x))) names(x) <- parnames
  missing <- setdiff(parnames, names(x))
  if (length(missing))
    stop("`", what, "` has no entry for ", paste(missing, collapse = ", "), ".")
  as.numeric(x[parnames])
}


# Internal: shared assembly for the soft constraints below. `term(x, k)` takes
# the current values of the constrained parameters present in the call and the
# indices they occupy, and returns the -2 log density plus its first two
# derivatives, one entry per parameter. A parameter passed in `fixed`
# contributes to the value but not to gradient or Hessian.
.constraintTerms <- function(parnames, term, attr.name, condition) {

  myfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, hessian = TRUE,
                   conditions = condition, env = NULL,
                   cores = getOption("dMod.cores", 1L)) {

    p    <- list(...)[[match.fnargs(list(...), "pars")]]
    build_hessian <- isTRUE(deriv) && isTRUE(hessian)
    dP   <- if (deriv) attr(p, "deriv", exact = TRUE) else NULL
    dP2  <- if (build_hessian && deriv2) attr(p, "deriv2", exact = TRUE) else NULL

    allp <- c(p, fixed)
    np   <- length(p)
    gr   <- setNames(numeric(np), names(p))
    hs   <- if (build_hessian) matrix(0, np, np, dimnames = list(names(p), names(p))) else NULL

    k   <- which(parnames %in% names(allp))
    td  <- term(as.numeric(allp[parnames[k]]), k)
    val <- sum(td$value)

    if (deriv) {
      # Only the free parameters carry derivatives; the rest came in `fixed`.
      pos <- match(parnames[k], names(p))
      free <- !is.na(pos)
      gr[pos[free]] <- td$d1[free]
      if (build_hessian) hs[cbind(pos[free], pos[free])] <- td$d2[free]
    }

    # Chain rule through an upstream parfn, so `constraint * P()` is exact.
    if (deriv && !is.null(dP)) {
      gi <- gr
      gr <- drop(gi %*% dP); names(gr) <- colnames(dP)
      if (build_hessian) {
        hs <- t(dP) %*% hs %*% dP
        dimnames(hs) <- list(colnames(dP), colnames(dP))

        if (!is.null(dP2)) {
          common <- intersect(names(gi), dimnames(dP2)[[1]])
          if (length(common) > 0L) {
            theta   <- colnames(dP)
            dP2_sub <- dP2[common, theta, theta, drop = FALSE]
            flat    <- matrix(dP2_sub, nrow = length(common),
                              ncol = length(theta)^2)
            hs <- hs + matrix(crossprod(flat, gi[common]),
                              length(theta), length(theta))
          }
        }
      }
    }

    out <- objlist(value = unname(val),
                   gradient = if (deriv) gr else NULL,
                   hessian  = if (build_hessian) hs else NULL)
    attr(out, attr.name) <- out$value
    attr(out, "env") <- env
    out
  }

  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- condition
  attr(myfn, "parameters") <- parnames
  myfn
}


#' Soft L1 constraint on parameters
#'
#' @param mu Named numeric vector of prior locations. Its names select the
#'   constrained parameters.
#' @param sigma Numeric, scalar or named and aligned with `mu`. The Laplace
#'   scale, i.e. the reciprocal penalty strength.
#' @param attr.name Character. Name of the attribute storing the constraint value.
#' @param condition Optional character vector of conditions.
#'
#' @details
#' Computes the Laplace prior on the `-2 log` scale,
#' \deqn{2\log(2\sigma) + 2|p-\mu|/\sigma,}
#' the L1 counterpart of [constraintL2]. Note the difference in normalisation:
#' this constructor and its siblings are the full `-2 log` density, while
#' [constraintL2] is the penalty form and drops the Gaussian's constant.
#' The term is not differentiable at
#' \eqn{p = \mu}, where the gradient is reported as 0; [trustL1] penalises the
#' same shape inside the optimiser and handles that kink exactly, so prefer it
#' for fitting and use this constructor to score a posterior.
#'
#' @return Object of class \code{objfn}.
#' @seealso [constraintL2], [trustL1]
#' @examples
#' prior <- constraintL1(mu = c(k1 = 0), sigma = 2)
#' prior(pars = c(k1 = 1))$value
#' @export
constraintL1 <- function(mu, sigma = 1, attr.name = "prior", condition = NULL) {

  parnames <- names(mu)
  mu    <- .constraintArg(mu,    parnames, "mu",    "constraintL1")
  sigma <- .constraintArg(sigma, parnames, "sigma", "constraintL1")

  .constraintTerms(parnames, function(x, k) {
    d <- x - mu[k]; s <- sigma[k]
    list(value = 2 * log(2 * s) + 2 * abs(d) / s,
         d1 = 2 * sign(d) / s, d2 = numeric(length(x)))
  }, attr.name, condition)
}


#' Soft Cauchy constraint on parameters
#'
#' @param mu Named numeric vector of prior locations. Its names select the
#'   constrained parameters.
#' @param sigma Numeric, scalar or named and aligned with `mu`. The Cauchy
#'   scale.
#' @param attr.name Character. Name of the attribute storing the constraint value.
#' @param condition Optional character vector of conditions.
#'
#' @details
#' Computes the Cauchy prior on the `-2 log` scale,
#' \deqn{2\log(\pi\sigma) + 2\log(1 + ((p-\mu)/\sigma)^2).}
#' The heavy tails make it the tolerant alternative to [constraintL2]: a
#' parameter far from `mu` is pulled far more weakly than a Gaussian would
#' pull it.
#'
#' @return Object of class \code{objfn}.
#' @seealso [constraintL2], [constraintL1]
#' @examples
#' prior <- constraintCauchy(mu = c(k1 = 0), sigma = 2)
#' prior(pars = c(k1 = 1))$value
#' @export
constraintCauchy <- function(mu, sigma = 1, attr.name = "prior", condition = NULL) {

  parnames <- names(mu)
  mu    <- .constraintArg(mu,    parnames, "mu",    "constraintCauchy")
  sigma <- .constraintArg(sigma, parnames, "sigma", "constraintCauchy")

  .constraintTerms(parnames, function(x, k) {
    s <- sigma[k]; t <- (x - mu[k]) / s
    list(value = 2 * log(pi * s) + 2 * log1p(t^2),
         d1 = 4 * t / (s * (1 + t^2)),
         d2 = 4 * (1 - t^2) / (s^2 * (1 + t^2)^2))
  }, attr.name, condition)
}


#' Soft gamma constraint on positive parameters
#'
#' @param shape Named numeric vector of shape parameters. Its names select the
#'   constrained parameters.
#' @param scale Numeric, scalar or named and aligned with `shape`. The gamma
#'   scale, not the rate.
#' @param attr.name Character. Name of the attribute storing the constraint value.
#' @param condition Optional character vector of conditions.
#'
#' @details
#' Computes the gamma prior on the `-2 log` scale,
#' \deqn{-2(a-1)\log p + 2p/s + 2\log\Gamma(a) + 2a\log s}
#' with shape \eqn{a} and scale \eqn{s}. The value is `Inf` for a
#' non-positive parameter.
#'
#' @return Object of class \code{objfn}.
#' @seealso [constraintExponential], [constraintChisq]
#' @examples
#' prior <- constraintGamma(shape = c(k1 = 3), scale = 5)
#' prior(pars = c(k1 = 5))$value
#' @export
constraintGamma <- function(shape, scale = 1, attr.name = "prior", condition = NULL) {

  parnames <- names(shape)
  shape <- .constraintArg(shape, parnames, "shape", "constraintGamma")
  scale <- .constraintArg(scale, parnames, "scale", "constraintGamma")

  .constraintTerms(parnames, function(x, k) {
    a <- shape[k]; s <- scale[k]; ok <- x > 0
    xo <- x[ok]; ao <- a[ok]; so <- s[ok]
    value <- rep(Inf, length(x)); d1 <- d2 <- numeric(length(x))
    value[ok] <- -2 * (ao - 1) * log(xo) + 2 * xo / so +
                 2 * lgamma(ao) + 2 * ao * log(so)
    d1[ok] <- -2 * (ao - 1) / xo + 2 / so
    d2[ok] <- 2 * (ao - 1) / xo^2
    list(value = value, d1 = d1, d2 = d2)
  }, attr.name, condition)
}


#' Soft exponential constraint on positive parameters
#'
#' @param scale Named numeric vector of scale parameters, not rates. Its names
#'   select the constrained parameters.
#' @param attr.name Character. Name of the attribute storing the constraint value.
#' @param condition Optional character vector of conditions.
#'
#' @details
#' Computes the exponential prior on the `-2 log` scale,
#' \deqn{2p/s + 2\log s.}
#' It is the gamma prior at shape 1 and pulls a parameter towards 0 with a
#' constant force, which makes it the smooth one-sided counterpart of an L1
#' penalty at `mu = 0`. The value is `Inf` for a negative parameter.
#'
#' @return Object of class \code{objfn}.
#' @seealso [constraintGamma], [constraintL1]
#' @examples
#' prior <- constraintExponential(scale = c(k1 = 3))
#' prior(pars = c(k1 = 5))$value
#' @export
constraintExponential <- function(scale, attr.name = "prior", condition = NULL) {

  parnames <- names(scale)
  scale <- .constraintArg(scale, parnames, "scale", "constraintExponential")

  .constraintTerms(parnames, function(x, k) {
    s <- scale[k]; ok <- x >= 0
    value <- rep(Inf, length(x)); d1 <- d2 <- numeric(length(x))
    value[ok] <- 2 * x[ok] / s[ok] + 2 * log(s[ok])
    d1[ok] <- 2 / s[ok]
    list(value = value, d1 = d1, d2 = d2)
  }, attr.name, condition)
}


#' Soft chi-squared constraint on positive parameters
#'
#' @param df Named numeric vector of degrees of freedom. Its names select the
#'   constrained parameters.
#' @param attr.name Character. Name of the attribute storing the constraint value.
#' @param condition Optional character vector of conditions.
#'
#' @details
#' Computes the chi-squared prior on the `-2 log` scale,
#' \deqn{-(k-2)\log p + p + k\log 2 + 2\log\Gamma(k/2)}
#' with `k` degrees of freedom. The value is `Inf` for a non-positive
#' parameter.
#'
#' @return Object of class \code{objfn}.
#' @seealso [constraintGamma]
#' @examples
#' prior <- constraintChisq(df = c(k1 = 4))
#' prior(pars = c(k1 = 5))$value
#' @export
constraintChisq <- function(df, attr.name = "prior", condition = NULL) {

  parnames <- names(df)
  df <- .constraintArg(df, parnames, "df", "constraintChisq")

  .constraintTerms(parnames, function(x, k) {
    v <- df[k]; ok <- x > 0
    xo <- x[ok]; vo <- v[ok]
    value <- rep(Inf, length(x)); d1 <- d2 <- numeric(length(x))
    value[ok] <- -(vo - 2) * log(xo) + xo + vo * log(2) + 2 * lgamma(vo / 2)
    d1[ok] <- -(vo - 2) / xo + 1
    d2[ok] <- (vo - 2) / xo^2
    list(value = value, d1 = d1, d2 = d2)
  }, attr.name, condition)
}


#' Soft Rayleigh constraint on positive parameters
#'
#' @param sigma Named numeric vector of scale parameters. Its names select the
#'   constrained parameters.
#' @param attr.name Character. Name of the attribute storing the constraint value.
#' @param condition Optional character vector of conditions.
#'
#' @details
#' Computes the Rayleigh prior on the `-2 log` scale,
#' \deqn{-2\log p + 4\log\sigma + p^2/\sigma^2.}
#' Unlike the exponential it vanishes at 0, so it keeps a parameter away from
#' both 0 and large values. The value is `Inf` for a non-positive parameter.
#'
#' @return Object of class \code{objfn}.
#' @seealso [constraintGamma], [constraintExponential]
#' @examples
#' prior <- constraintRayleigh(sigma = c(k1 = 3))
#' prior(pars = c(k1 = 5))$value
#' @export
constraintRayleigh <- function(sigma, attr.name = "prior", condition = NULL) {

  parnames <- names(sigma)
  sigma <- .constraintArg(sigma, parnames, "sigma", "constraintRayleigh")

  .constraintTerms(parnames, function(x, k) {
    s <- sigma[k]; ok <- x > 0
    xo <- x[ok]; so <- s[ok]
    value <- rep(Inf, length(x)); d1 <- d2 <- numeric(length(x))
    value[ok] <- -2 * log(xo) + 4 * log(so) + xo^2 / so^2
    d1[ok] <- -2 / xo + 2 * xo / so^2
    d2[ok] <- 2 / xo^2 + 2 / so^2
    list(value = value, d1 = d1, d2 = d2)
  }, attr.name, condition)
}

#' L2 objective function for validation data point
#' 
#' @param name character, the name of the prediction, e.g. a state name.
#' @param time numeric, the time-point associated to the prediction
#' @param value character, the name of the parameter which contains the
#' prediction value.
#' @param sigma numeric, the uncertainty of the introduced test data point
#' @param attr.name character. The constraint value is additionally returned in
#' an attribute of this name. Being a squared standardised residual already, it
#' is also the term's `chi2` contribution, which pools with [normL2]'s under the
#' rules described there.
#' @param condition character, the condition for which the prediction is made.
#' @return List of class \code{objlist}, i.e. objective value, gradient and Hessian as list.
#' @seealso [normL2], [constraintL2]
#' @details Computes the constraint value 
#' \deqn{\left(\frac{x(t)-\mu}{\sigma}\right)^2}{(pred-p[names(mu)])^2/sigma^2}
#' and its derivatives with respect to p.
#' @examples
#' prediction <- list(a = matrix(c(0, 1), nrow = 1, dimnames = list(NULL, c("time", "A"))))
#' derivs <- matrix(c(0, 1, 0.1), nrow = 1, dimnames = list(NULL, c("time", "A.A", "A.k1")))
#' attr(prediction$a, "deriv") <- derivs
#' p0 <- c(A = 1, k1 = 2)
#' 
#' vali <- datapointL2(name = "A", time = 0, value = "newpoint", sigma = 1, condition = "a")
#' vali(pars = c(p0, newpoint = 1), env = .GlobalEnv)
#' @export
datapointL2 <- function(name, time, value, sigma = 1, attr.name = "validation", condition) {

  controls <- list(
    mu        = structure(name, names = value)[1], # only one data point is allowed
    time      = time[1],
    sigma     = sigma[1],
    attr.name = attr.name
  )

  myfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, hessian = TRUE,
                   conditions = NULL, env = NULL,
                   cores = getOption("dMod.cores", 1L)) {
    build_hessian <- isTRUE(deriv) && isTRUE(hessian)
    mu        <- controls$mu
    t         <- controls$time
    sigma     <- controls$sigma
    attr.name <- controls$attr.name

    arglist <- list(...)
    arglist <- arglist[match.fnargs(arglist, "pars")]
    pouter  <- arglist[[1]]
    if (is.null(env)) {
      stop("No prediction available. Use the argument env to pass an environment that contains the prediction.")
    }
    prediction <- as.list(env)$prediction

    if (!is.null(conditions) && !condition %in% conditions)
      return()
    if (is.null(conditions) && !condition %in% names(prediction))
      stop("datapointL2 requests unavailable condition. Call the objective function explicitly stating the conditions argument.")

    prdf <- prediction[[condition]]
    if (!any(prdf[, "time"] == t))
      stop("datapointL2() requests time point for which no prediction is available. Please add missing time point by the times argument in normL2()")

    dpred_attr  <- if (deriv) attr(prdf, "deriv") else NULL
    d2pred_attr <- if (build_hessian && deriv2) attr(prdf, "deriv2") else NULL
    kr <- datapointL2_kernel(
      pouter           = pouter,
      fixed_opt        = fixed,
      prdf             = prdf,
      dpred_attr_opt   = dpred_attr,
      d2pred_attr_opt  = d2pred_attr,
      obs_name         = as.character(mu),
      t                = as.numeric(t),
      sigma            = as.numeric(sigma),
      value_par        = names(mu)[1],
      deriv            = deriv,
      build_hessian    = build_hessian
    )

    out <- objlist(value = kr$value, gradient = kr$gradient, hessian = kr$hessian)
    attr(out, attr.name)    <- out$value
    # The value is the squared standardised residual itself, no normaliser.
    attr(out, "chi2")       <- setNames(out$value, attr.name)
    attr(out, "prediction") <- kr$prediction
    attr(out, "env")        <- env
    class(out) <- NULL
    out
  }
  class(myfn)             <- c("objfn", "fn")
  attr(myfn, "conditions") <- condition
  attr(myfn, "parameters") <- value[1]
  myfn
}


#' Add two lists element by element
#' 
#' @param out1 List of numerics or matrices
#' @param out2 List with the same structure as out1 (there will be no warning when mismatching)
#' @details If out1 has names, out2 is assumed to share these names. Each element of the list out1
#' is inspected. If it has a \code{names} attributed, it is used to do a matching between out1 and out2.
#' The same holds for the attributed \code{dimnames}. In all other cases, the "+" operator is applied
#' the corresponding elements of out1 and out2 as they are.
#' @return List of length of out1. 
#' @aliases sumobjlist
#' @export
#' 
"+.objlist" <- function(out1, out2) {

  if (is.null(out1)) return(out2)
  if (is.null(out2)) return(out1)

  gn1 <- names(out1$gradient)
  gn2 <- names(out2$gradient)

  # Layout of the sum: the operand spanning the other, else their union.
  pars <- if (all(gn2 %in% gn1)) gn1
          else if (all(gn1 %in% gn2)) gn2
          else union(gn1, gn2)

  addVector <- function(target, x) {
    i <- intersect(names(target), names(x))
    target[i] <- target[i] + x[i]
    target
  }
  addMatrix <- function(target, x) {
    i <- intersect(rownames(target), rownames(x))
    target[i, i] <- target[i, i] + x[i, i]
    target
  }

  what <- intersect(c("value", "gradient", "hessian"), c(names(out1), names(out2)))
  out12 <- lapply(what, function(w) switch(w,
    value    = out1$value + out2$value,
    gradient = addVector(addVector(setNames(numeric(length(pars)), pars),
                                   out1$gradient), out2$gradient),
    # A summand may carry a NULL hessian (built with hessian = FALSE); the sum
    # is NULL only when both are, otherwise the present ones add.
    hessian  = if (is.null(out1$hessian) && is.null(out2$hessian)) NULL else {
      H <- matrix(0, length(pars), length(pars), dimnames = list(pars, pars))
      if (!is.null(out1$hessian)) H <- addMatrix(H, out1$hessian)
      if (!is.null(out2$hessian)) H <- addMatrix(H, out2$hessian)
      H
    }))
  names(out12) <- what

  # Numeric attributes are summed, an absent one counting as zero. The chi2
  # contributions are kept apart from that: they are pooled by the `attr.name`
  # they belong to, not by the attribute they happen to sit under.
  numeric_attrs <- function(x) {
    a <- attributes(x)
    a <- a[vapply(a, is.numeric, logical(1))]
    a[!grepl("^chi2($|_)", names(a))]
  }
  a1 <- numeric_attrs(out1)
  a2 <- numeric_attrs(out2)
  for (n in union(names(a1), names(a2)))
    attr(out12, n) <- (if (is.null(a1[[n]])) 0 else a1[[n]]) +
                      (if (is.null(a2[[n]])) 0 else a2[[n]])

  chi2 <- c(.chi2Contributions(out1), .chi2Contributions(out2))
  if (length(chi2)) {
    chi2 <- vapply(split(unname(chi2), names(chi2)), sum, 0)
    if (length(chi2) == 1L) attr(out12, "chi2") <- chi2
    else for (n in names(chi2)) attr(out12, paste0("chi2_", n)) <- unname(chi2[n])
  }

  class(out12) <- "objlist"

  out12
}


# The chi2 contributions of an objlist as one named vector, however they are
# currently stored: a single one sits under `chi2`, several under `chi2_<name>`.
.chi2Contributions <- function(x) {
  a <- attributes(x)
  a <- a[grepl("^chi2($|_)", names(a))]
  if (!length(a)) return(numeric(0))
  out <- unlist(a, use.names = FALSE)
  names(out) <- unlist(lapply(names(a), function(n)
    if (identical(n, "chi2")) names(a[[n]]) else sub("^chi2_", "", n)))
  out
}


#' @export
print.objlist <- function(x, n1 = 20, n2 = 6, ...) {
  cat("value\n", "==================\n",x$value, "\n")
  # An objlist from a `deriv = FALSE` call carries the value alone.
  if (length(x$gradient)) {
    n1 <- min(n1, length(x$gradient))
    cat("gradient[1:",n1,"] (full length = ",length(x$gradient),")\n", "==================\n", sep = "")
    print(x$gradient[1:n1])
    cat("\n")
  }
  if (length(x$hessian)) {
    n2 <- min(n2, nrow(x$hessian))
    cat("hessian[1:",n2,",1:",n2,"]","\n", "==================\n", sep = "")
    print(x$hessian[1:n2,1:n2, drop = FALSE])
    cat("\n")
  }
  cat("\n")
  cat("attributes\n", "==================\n")
  # str() would prefix every line with the storage mode, and the chi2 tag with
  # a line of its own. The attributes carry objective contributions, so what
  # matters is the name and the number.
  a <- attributes(x)[setdiff(names(attributes(x)), c("names", "class"))]
  if (length(a)) {
    w <- max(nchar(names(a)))
    for (n in names(a)) cat(sprintf(" %-*s  %s\n", w, n, .objlistAttrText(a[[n]])))
  }
  invisible(x)
}


# One line for one objlist attribute, no storage mode and no nesting.
.objlistAttrText <- function(v) {
  if (is.environment(v)) return("<environment>")
  if (is.numeric(v) || is.character(v))
    return(paste(format(unname(v), trim = TRUE), collapse = ", "))
  paste0("<", paste(class(v), collapse = "/"), ">")
}



#' @export
print.objfn <- function(x, ...) {

  parameters <- attr(x, "parameters")

  cat("Objective function:\n")
  str(args(x))
  cat("\n")
  cat("... parameters:", paste0(parameters, collapse = ", "), "\n")

}


#' @export
summary.objfn <- function(object, ...) {

  x <- object

  parameters <- attr(x, "parameters")
  conditions <- attr(x, "conditions")
  modelnames <- attr(x, "modelname")

  cat("Details:\n")
  cat("... class:      ", paste0(class(x), collapse = ", "), "\n")
  cat("... parameters: ", paste0(parameters, collapse = ", "), "\n")
  if (!is.null(conditions))
    cat("... conditions: ", paste0(conditions, collapse = ", "), "\n")
  if (!is.null(modelnames))
    cat("... modelname:  ", paste0(modelnames, collapse = ", "), "\n")

  ctrls <- try(controls(x), silent = TRUE)
  if (!inherits(ctrls, "try-error") && length(ctrls))
    cat("... controls:   ", paste0(ctrls, collapse = ", "), "\n")

  invisible(list(class = class(x), parameters = parameters,
                 conditions = conditions, modelname = modelnames))

}



## res (moved from data.R) ---------------------------------------------------

#' Compute residuals between data and model prediction
#'
#' Matches data to predictions by time and observable, computes (weighted)
#' residuals, and propagates parameter derivatives. Values below `lloq` are
#' censored via `pmax(value, lloq)`.
#'
#' @md
#' @param data Data frame with columns `time`, `name`, `value`, `sigma`, `lloq`.
#'   Rows with `sigma = NA` are filled from `err`.
#' @param out Prediction matrix (first column = time, remaining = observables).
#'   Optional `"deriv"` attribute: `[name, param, time]` array.
#' @param err Optional error-model matrix (same layout as `out`).
#'   Optional `"deriv"` attribute: `[name, param, time]` array.
#'
#' @details
#' The returned `"deriv"` and `"deriv.err"` matrices have shape
#' \eqn{n \times p}{n x p} (residuals x parameters), extracted from the
#' `[name, param, time]` arrays on `out` and `err`.
#'
#' @return An [objframe()] with columns `time`, `name`, `value`, `prediction`,
#'   `sigma`, `residual`, `weighted.residual`, `bloq`, `weighted.0` and
#'   attributes `"deriv"` and `"deriv.err"`.
#'
#' @seealso [objframe()]
#' @export
res <- function(data, out, err = NULL) {
  
  data$name <- as.character(data$name)
  n <- nrow(data)
  times <- sort(unique(data$time))
  names <- unique(data$name)
  
  ti <- .matchNum(times, out[, 1])[.matchNum(data$time, times)]
  ni <- match(names, colnames(out))[match(data$name, names)]
  if (anyNA(ni))
    stop("Observable not found: ",
         paste(setdiff(names, colnames(out)), collapse = ", "))
  if (anyNA(ti)) stop("Some data$time not found in out[,1]")
  
  pred <- out[cbind(ti, ni)]
  
  deriv <- NULL
  if (!is.null(d <- attr(out, "deriv"))) {
    oi <- match(data$name, dimnames(d)[[2]])
    np <- dim(d)[3]
    deriv <- matrix(
      d[cbind(rep(ti, np), rep(oi, np), rep(seq_len(np), each = n))],
      n, np, dimnames = list(NULL, dimnames(d)[[3]]))
  }

  deriv2 <- NULL
  if (!is.null(d2 <- attr(out, "deriv2"))) {
    oi2 <- match(data$name, dimnames(d2)[[2]])
    np2 <- dim(d2)[3]
    # Build [n*np*np x 4] index matrix; outermost loop = k, then j, then i.
    idx <- cbind(
      rep(ti,  np2 * np2),
      rep(oi2, np2 * np2),
      rep(rep(seq_len(np2), each = n), np2),
      rep(seq_len(np2), each = n * np2)
    )
    deriv2 <- array(d2[idx], c(n, np2, np2),
                    dimnames = list(NULL, dimnames(d2)[[3]], dimnames(d2)[[4]]))
  }

  sig  <- data$sigma
  sNA  <- is.na(sig)
  derr <- NULL
  derr2 <- NULL

  if (any(sNA)) {
    if (is.null(err)) stop("NA sigmas but no errmodel")
    ti_e <- .matchNum(times, err[, 1])[.matchNum(data$time, times)]
    ni_e <- match(names, colnames(err))[match(data$name, names)]
    sig[sNA] <- err[cbind(ti_e, ni_e)][sNA]

    if (!is.null(de <- attr(err, "deriv"))) {
      oi <- match(data$name, dimnames(de)[[2]])
      np <- dim(de)[3]
      ns <- sum(sNA)
      derr <- matrix(0, n, np, dimnames = list(NULL, dimnames(de)[[3]]))
      derr[sNA, ] <- matrix(
        de[cbind(rep(ti_e[sNA], np), rep(oi[sNA], np), rep(seq_len(np), each = ns))],
        ns, np)
    }

    if (!is.null(de2 <- attr(err, "deriv2"))) {
      oi <- match(data$name, dimnames(de2)[[2]])
      np2 <- dim(de2)[3]
      ns <- sum(sNA)
      derr2 <- array(0, c(n, np2, np2),
                     dimnames = list(NULL, dimnames(de2)[[3]], dimnames(de2)[[4]]))
      idx <- cbind(
        rep(ti_e[sNA],  np2 * np2),
        rep(oi[sNA],    np2 * np2),
        rep(rep(seq_len(np2), each = ns), np2),
        rep(seq_len(np2), each = ns * np2)
      )
      derr2[sNA, , ] <- array(de2[idx], c(ns, np2, np2))
    }
  }

  val  <- pmax(data$value, data$lloq)
  resi <- pred - val
  inv  <- 1 / sig

  objframe(
    data.table::data.table(
      time = data$time, name = data$name, value = val,
      prediction = pred, sigma = sig, residual = resi,
      weighted.residual = resi * inv,
      bloq = val <= data$lloq, weighted.0 = pred * inv),
    deriv = deriv, deriv.err = derr,
    deriv2 = deriv2, deriv2.err = derr2)
}


## objlist / objframe constructors (moved from classes.R) ----------------------------------------

## Objective classes ---------------------------------------------------------


#' Generate objective list
#'
#' @description An objective list contains an objective value, a gradient, and a Hessian matrix.
#'
#' Objective lists can contain additional numeric attributes that are preserved or
#' combined with the corresponding attributes of another objective list when
#' both are added by the "+" operator, see [sumobjlist].
#'
#' Objective lists are returned by objective functions as being generated
#' by [normL2], [constraintL2] and [datapointL2].
#' @param value numeric of length 1
#' @param gradient named numeric
#' @param hessian matrix with rownames and colnames according to gradient names
#' @return Object of class `objlist`
#' @export
#' 
#' @examples
#' # objlist(1, c(a = 1, b = 2),
#' #         matrix(2, nrow = 2, ncol = 2,
#' #                dimnames = list(c("a", "b"), c("a", "b"))))
objlist <- function(value, gradient, hessian) {

  out <- list(value = value, gradient = gradient, hessian = hessian)
  class(out) <- c("objlist", "list")
  return(out)

}


#' Objective frame
#'
#' @description
#' An objective frame stores residuals and their derivatives with respect to parameters.
#' It is typically created by [res] and used internally in objective functions.
#'
#' @param mydata data.table produced by [res]
#' @param deriv numeric matrix of first-order derivatives of residuals (Jacobian)
#' @param deriv.err numeric matrix of first-order derivatives of the error model
#' @param deriv2 numeric 3D array `[n_residuals, p, p]` of second-order derivatives
#'   of residuals with respect to parameters. Optional.
#' @param deriv2.err numeric 3D array `[n_residuals, p, p]` of second-order
#'   derivatives of the error model. Optional.
#'
#' @return
#' An object of class `"objframe"` (data.table) with attributes `"deriv"` and `"deriv.err"`.
#' These arrays have the same parameter axes as those returned by [prdframe] and [res].
#' When `deriv2`/`deriv2.err` are supplied, the corresponding 3D arrays are
#' attached as `"deriv2"` / `"deriv2.err"`.
#'
#' @export
objframe <- function(mydata, deriv = NULL, deriv.err = NULL,
                     deriv2 = NULL, deriv2.err = NULL) {

  required <- c("time", "name", "value", "prediction",
                "sigma", "residual", "weighted.residual",
                "bloq", "weighted.0")
  if (!all(required %in% names(mydata)))
    stop("mydata does not have all required columns.")

  out <- data.table::as.data.table(mydata)[, ..required]
  data.table::setattr(out, "deriv",      deriv)
  data.table::setattr(out, "deriv.err",  deriv.err)
  data.table::setattr(out, "deriv2",     deriv2)
  data.table::setattr(out, "deriv2.err", deriv2.err)
  data.table::setattr(out, "class", c("objframe", "data.table", "data.frame"))
  out
}




