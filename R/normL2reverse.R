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
                            attr.name, hessian = FALSE) {

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
    got <- lapply(seq_along(cn_eval), function(j)
      errmodel(out = prediction[[cn_eval[j]]], pars = err_pars[[j]],
               fixed = err_fixed[[j]], deriv = FALSE,
               conditions = cn_eval[j])[[cn_eval[j]]])
    err_list <- vector("list", length(conditions))
    err_list[match(cn_eval, conditions)] <- got
  }

  meta_list <- .build_normL2_meta(data, prediction, err_list, conditions, e.cond)

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

  # --- scatter the seed back onto the two things it belongs to -------------
  #
  # The kernel orders its rows ALOQ first, then BLOQ, so the scatter follows
  # the same permutation. A row with a fixed sigma has no sigma derivative at
  # all, and its sigma seed is dropped rather than multiplied by a zero.
  w_pred <- vector("list", length(conditions))
  w_err  <- vector("list", length(conditions))
  for (ci in seq_along(conditions)) {
    m  <- meta_list[[ci]]
    pr <- prediction[[conditions[ci]]]
    ord <- c(which(m$bloq_mask == 0L), which(m$bloq_mask == 1L))

    W <- matrix(0, nrow(pr), ncol(pr), dimnames = list(NULL, colnames(pr)))
    sp <- kr$seed$pred[[ci]]
    for (j in seq_along(ord)) {
      r <- ord[j]
      W[m$t_idx_in_pred[r], m$o_idx_in_pred[r]] <-
        W[m$t_idx_in_pred[r], m$o_idx_in_pred[r]] + sp[j]
    }
    w_pred[[ci]] <- W

    erm <- err_list[[ci]]
    if (!is.null(erm) && any(m$sigma_is_na == 1L)) {
      E <- matrix(0, nrow(erm), ncol(erm), dimnames = list(NULL, colnames(erm)))
      ss <- kr$seed$sigma[[ci]]
      for (j in seq_along(ord)) {
        r <- ord[j]
        if (m$sigma_is_na[r] != 1L) next
        E[m$t_idx_in_err[r], m$o_idx_in_err[r]] <-
          E[m$t_idx_in_err[r], m$o_idx_in_err[r]] + ss[j]
      }
      w_err[[ci]] <- E
    }
  }

  # --- the error model backwards -------------------------------------------
  #
  # It reads the prediction's values and its parameters, so its cotangent lands
  # on both, and the prediction's half adds to the seed above.
  w_chain <- lapply(seq_along(conditions), function(ci)
    .ct(out = .ctWiden(.ct(out = .dropTime(w_pred[[ci]]))$out, n_dir + 1L)))
  if (!is.null(errmodel) && length(cn_eval)) {
    evjp <- attr(.fnLeafKernel(errmodel), "vjpfn")
    if (is.null(evjp))
      stop("normL2: the error model has no vjp entry; rebuild it with ",
           "Y(..., compile = TRUE).", call. = FALSE)
    for (j in seq_along(cn_eval)) {
      ci <- match(cn_eval[j], conditions)
      if (is.null(w_err[[ci]])) next
      u <- evjp(out = prediction[[cn_eval[j]]], pars = err_pars[[j]],
                fixed = err_fixed[[j]],
                w = .ctWiden(.ct(out = .dropTime(w_err[[ci]]))$out, n_dir + 1L))
      w_chain[[ci]] <- .addCt(w_chain[[ci]], .ct(out = .dropTime(u$out),
                                                 pars = u$pars))
    }
  }

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

# The raw kernel behind a single-leaf fn, which is where the vjp attribute
# sits. An error model built by Y() is exactly that.
.fnLeafKernel <- function(f) {
  st <- .fnNode(f)
  if (is.null(st) || !identical(st$op, "leaf")) return(f)
  st$kernel
}
