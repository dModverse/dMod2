## normL2 backwards --------------------------------------------------------
##
## The forward objective asks the prediction chain for dX/dtheta and contracts
## it with the residuals. The reverse one asks for values only, works out what
## the objective's derivative in those values is, and pushes that back through
## the chain. The gradient is then one sweep wide instead of n_theta.
##
## The seed is not the residual vector. The error model depends on theta too,
## so every data row seeds two things, the prediction and sigma, and the kernel
## hands both back: see `want_seed` in src/residual_kernel.cpp, where they fall
## out of the same per-row coefficients the gradient is built from.
##
## Copyright (C) 2026 Simon Beyer

.normL2_reverse <- function(pars, fixed, deriv, conditions, env, cores,
                            x, errmodel, data, timesD, e.cond, opt.BLOQ,
                            attr.name, hessian = FALSE, meta_cache = NULL) {

  # The objective's Hessian splits along a line the residual kernel already
  # draws: J' H_rho J from the forward tangents, which is what the kernel
  # computes when it is handed no second derivatives, plus the prediction's own
  # curvature weighted by the seed, which is what a dual sweep with a constant
  # seed returns. The seed and that weight are the same number by construction,
  # see src/residual_kernel.h.
  n_dir <- if (isTRUE(hessian)) length(pars) else 0L
  if (n_dir > 0L) {
    nm <- names(pars)
    attr(pars, "deriv") <- diag(length(nm))
    dimnames(attr(pars, "deriv")) <- list(nm, nm)
  }

  # --- forward, keeping the tape, and the tangents when second order needs
  #     them: they are the directions every node's vjp is differentiated along
  b <- .bundle_from_call(conditions, times = timesD, out = NULL,
                         pars = pars, fixed = fixed)
  fw <- .fwdMany(x, b, env, cores, deriv = n_dir > 0L)
  prediction <- as.prdlist(fw$values)
  prediction <- prediction[conditions]

  # --- the error model, values only ----------------------------------------
  err_list <- NULL
  err_pars <- err_fixed <- NULL
  cn_eval <- character(0)
  if (!is.null(errmodel)) {
    cn_eval <- if (is.null(e.cond)) conditions else .intersectU(conditions, e.cond)
    fixed_names <- names(fixed)
    split <- lapply(cn_eval, function(cn) {
      pinner <- getParameters(prediction[[cn]])
      own <- attr(pinner, "fixed")
      fixedinner <- pinner[c(own, .setdiffU(.intersectU(names(pinner), fixed_names), own))]
      list(pars  = as.parvec(pinner[.setdiffU(names(pinner), fixed_names)]),
           fixed = as.parvec(fixedinner, deriv = FALSE, deriv2 = FALSE))
    })
    err_pars  <- lapply(split, `[[`, "pars")
    err_fixed <- lapply(split, `[[`, "fixed")
    est <- .fnNode(errmodel)
    got <- if (!is.null(est) && length(cn_eval) > 1L) {
      eb <- .bundle(conds = cn_eval,
                    out   = lapply(cn_eval, function(cn) prediction[[cn]]),
                    pars  = err_pars, fixed = err_fixed, shared = FALSE)
      .evalNode(est, eb, FALSE, FALSE, NULL, cores)
    } else {
      lapply(seq_along(cn_eval), function(j)
        errmodel(out = prediction[[cn_eval[j]]], pars = err_pars[[j]],
                 fixed = err_fixed[[j]], deriv = FALSE,
                 conditions = cn_eval[j])[[cn_eval[j]]])
    }
    err_list <- vector("list", length(conditions))
    err_list[match(cn_eval, conditions)] <- got
  }

  # The data-to-prediction indices hold while names, tangents and row counts do.
  sig <- list(lapply(prediction, function(pr) dimnames(attr(pr, "deriv"))[[3]]),
              vapply(prediction, NROW, integer(1)),
              vapply(err_list, NROW, integer(1)))
  hit <- !is.null(meta_cache) && identical(meta_cache$sig, sig)
  meta_list <- if (hit) meta_cache$meta_list else
    .build_normL2_meta(data, prediction, err_list, conditions, e.cond)
  if (!hit && !is.null(meta_cache)) {
    meta_cache$meta_list <- meta_list
    meta_cache$sig <- sig
  }

  kr <- normL2_kernel(
    prediction       = prediction,
    err_list_opt     = err_list,
    meta_list        = meta_list,
    par_names_global = if (n_dir > 0L) names(pars) else character(0),
    deriv2_requested = FALSE,
    threads          = as.integer(cores),
    bloq_mode        = opt.BLOQ,
    build_hessian    = n_dir > 0L,
    want_seed        = isTRUE(deriv)
  )

  if (!isTRUE(deriv)) {
    out <- objlist(value = kr$value, gradient = NULL, hessian = NULL)
    attr(out, attr.name) <- out$value
    attr(out, "chi2") <- setNames(kr$chi2, attr.name)
    env$prediction <- prediction
    attr(out, "env") <- env
    return(out)
  }

  # --- the seed as cotangents on the prediction, and through the error model
  w_chain <- .normL2SeedCt(meta_list, prediction, err_list, kr, errmodel,
                           err_idx = match(cn_eval, conditions),
                           err_split = Map(function(a, b) list(pars = a, fixed = b),
                                           err_pars, err_fixed),
                           K = n_dir + 1L)

  # --- the chain backwards --------------------------------------------------
  u <- .bwdNode(fw$tape, w_chain, env, cores)
  grad <- Reduce(.addNamed,
                 lapply(u, function(z) if (is.null(z)) NULL else z$pars))
  grad <- .pickCotangent(grad, names(pars))
  gradient <- setNames(grad[, 1L], rownames(grad))

  # The two terms, on the same rows and in the same order: the kernel's is in
  # par_names_global order, the sweep's directions are the identity's columns.
  hess <- NULL
  if (n_dir > 0L) {
    hess <- kr$hessian + grad[, -1L, drop = FALSE]
    dimnames(hess) <- list(names(pars), names(pars))
  }

  out <- .alignObjlist(objlist(value = kr$value, gradient = gradient,
                               hessian = hess),
                       names(pars))
  attr(out, attr.name) <- out$value
  attr(out, "chi2") <- setNames(kr$chi2, attr.name)
  # Which direction answered. A caller used to read that off an absent Hessian,
  # which stops being a signal the moment the reverse mode can return one.
  attr(out, "sweep") <- if (n_dir > 0L) "forward-reverse" else "reverse"
  env$prediction <- prediction
  attr(out, "env") <- env
  out
}

# The objective's seed as one cotangent per condition, positionally aligned
# with `prediction`: on the prediction itself, and through the error model onto
# the prediction and the inner parameters it read. `err_idx` are the positions
# that have an error model and `err_split` the (pars, fixed) the forward pass
# handed it there. Shared by normL2 and the multiple-shooting objective.
#
# The kernel orders its rows ALOQ first, then BLOQ, so the scatter follows the
# same permutation. A row with a fixed sigma has no sigma derivative at all, and
# its sigma seed is dropped rather than multiplied by a zero.
.normL2SeedCt <- function(meta_list, prediction, err_list, kr, errmodel,
                          err_idx = NULL, err_split = NULL, K = 1L) {
  n <- length(meta_list)
  w_pred <- vector("list", n)
  w_err  <- vector("list", n)
  for (ci in seq_len(n)) {
    m  <- meta_list[[ci]]
    pr <- prediction[[ci]]
    ord <- c(which(m$bloq_mask == 0L), which(m$bloq_mask == 1L))

    W <- matrix(0, nrow(pr), ncol(pr), dimnames = list(NULL, colnames(pr)))
    pos <- m$t_idx_in_pred[ord] + (m$o_idx_in_pred[ord] - 1L) * nrow(pr)
    W[] <- .scatterAdd(length(W), pos, kr$seed$pred[[ci]])
    w_pred[[ci]] <- W

    erm <- if (is.null(err_list)) NULL else err_list[[ci]]
    if (!is.null(erm) && any(m$sigma_is_na == 1L)) {
      E <- matrix(0, nrow(erm), ncol(erm), dimnames = list(NULL, colnames(erm)))
      j <- which(m$sigma_is_na[ord] == 1L)
      r <- ord[j]
      pe <- m$t_idx_in_err[r] + (m$o_idx_in_err[r] - 1L) * nrow(erm)
      E[] <- .scatterAdd(length(E), pe, kr$seed$sigma[[ci]][j])
      w_err[[ci]] <- E
    }
  }

  # The error model reads the prediction's values and its parameters, so its
  # cotangent lands on both, and the prediction's half adds to the seed above.
  w_chain <- lapply(seq_len(n), function(ci)
    .ct(out = .ctWiden(.ct(out = .dropTime(w_pred[[ci]]))$out, K)))
  if (!is.null(errmodel) && length(err_idx)) {
    evjp <- attr(.fnLeafKernel(errmodel), "vjpfn")
    if (is.null(evjp))
      stop("normL2: the error model has no vjp entry; rebuild it with ",
           "Y(..., compile = TRUE).", call. = FALSE)
    for (j in seq_along(err_idx)) {
      ci <- err_idx[j]
      if (is.null(w_err[[ci]])) next
      u <- evjp(out = prediction[[ci]], pars = err_split[[j]]$pars,
                fixed = err_split[[j]]$fixed,
                cotangent = .ctWiden(.ct(out = .dropTime(w_err[[ci]]))$out, K))
      w_chain[[ci]] <- .addCt(w_chain[[ci]], .ct(out = .dropTime(u$out),
                                                 pars = u$pars))
    }
  }
  w_chain
}

# A zero vector of length n with v added at the linear positions idx, repeated
# positions summed.
.scatterAdd <- function(n, idx, v) {
  out <- numeric(n)
  if (!length(idx)) return(out)
  if (!anyDuplicated(idx)) {
    out[idx] <- v
    return(out)
  }
  s <- rowsum(v, idx)
  out[as.integer(rownames(s))] <- s[, 1L]
  out
}

# The raw kernel behind a single-leaf fn, which is where the vjp attribute
# sits. An error model built by Y() is exactly that.
.fnLeafKernel <- function(f) {
  st <- .fnNode(f)
  if (is.null(st) || !identical(st$op, "leaf")) return(f)
  st$kernel
}
