## normL2 backwards --------------------------------------------------------
## Copyright (C) 2026 Simon Beyer

## Values forward, the objective's derivative in them pulled back through the
## chain: one sweep instead of n_theta. Every data row seeds both prediction and
## sigma, see `want_seed` in src/residual_kernel.cpp.

.normL2_reverse <- function(pars, fixed, deriv, conditions, env, cores,
                            x, errmodel, data, timesD, e.cond, opt.BLOQ,
                            attr.name, hessian = FALSE, meta_cache = NULL) {

  # Hessian: J' H_rho J from the forward tangents (the kernel without second
  # derivatives) plus the prediction's curvature weighted by the seed (a dual
  # sweep with constant seed), see src/residual_kernel.h.
  n_dir <- if (isTRUE(hessian)) length(pars) else 0L
  if (n_dir > 0L) {
    nm <- names(pars)
    attr(pars, "deriv") <- diag(length(nm))
    dimnames(attr(pars, "deriv")) <- list(nm, nm)
  }

  # --- forward with tape, plus tangents as the vjp directions for second order
  b <- .bundle_from_call(conditions, times = timesD, out = NULL,
                         pars = pars, fixed = fixed)
  fw <- .fwdMany(x, b, env, cores, deriv = n_dir > 0L)
  prediction <- as.prdlist(fw$values)
  prediction <- prediction[conditions]

  # --- the error model on the predictions and their tangents, with a tape for its vjp
  err_list <- NULL
  err_pars <- err_fixed <- NULL
  err_tape <- NULL
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
    eb <- .bundle(conds = cn_eval,
                  out   = lapply(cn_eval, function(cn) prediction[[cn]]),
                  pars  = err_pars, fixed = err_fixed, shared = FALSE)
    ef <- .fwdMany(errmodel, eb, env, cores, deriv = n_dir > 0L)
    err_tape <- ef$tape
    err_list <- vector("list", length(conditions))
    err_list[match(cn_eval, conditions)] <- ef$values
  }

  # The data-to-prediction indices hold while names, tangents and row counts do.
  sig <- list(lapply(prediction, function(pr) dimnames(attr(pr, "deriv"))[[3]]),
              vapply(prediction, NROW, integer(1)),
              vapply(err_list, NROW, integer(1)),
              lapply(err_list, function(e) dimnames(attr(e, "deriv"))[2:3]))
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
  w_chain <- .normL2SeedCt(meta_list, prediction, err_list, kr, err_tape,
                           err_idx = match(cn_eval, conditions), env = env, cores = cores,
                           K = n_dir + 1L)

  # --- the chain backwards --------------------------------------------------
  u <- .bwdNode(fw$tape, w_chain, env, cores)
  grad <- .sumCotangents(lapply(u, function(z) if (is.null(z)) NULL else z$pars), names(pars))
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
  # Which direction answered; reverse mode can return a Hessian too.
  attr(out, "sweep") <- if (n_dir > 0L) "forward-reverse" else "reverse"
  env$prediction <- prediction
  attr(out, "env") <- env
  out
}

# The seed as one cotangent per condition aligned with `prediction`, direct and
# through the error model, whose forward tape `err_tape` covers the conditions at
# positions `err_idx`. Rows follow the kernel's ALOQ-then-BLOQ order; fixed-sigma
# rows get no sigma seed.
.normL2SeedCt <- function(meta_list, prediction, err_list, kr, err_tape,
                          err_idx = NULL, env = NULL, cores = 1L, K = 1L) {
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
  if (!is.null(err_tape) && length(err_idx)) {
    w_e <- lapply(err_idx, function(ci)
      if (is.null(w_err[[ci]])) .ct() else .ct(out = .ctWiden(.ct(out = .dropTime(w_err[[ci]]))$out, K)))
    u <- .bwdNode(err_tape, w_e, env, cores)
    for (j in seq_along(err_idx)) {
      if (j > length(u) || .ct_null(u[[j]])) next
      ci <- err_idx[j]
      w_chain[[ci]] <- .addCt(w_chain[[ci]], .ct(out = .dropTime(u[[j]]$out), pars = u[[j]]$pars))
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

