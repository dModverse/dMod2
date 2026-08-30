## Stage-2 d log|H_GN| / d theta correction. Called as an Rcpp::Function from
## the C++ FOCEI kernel once per outer iter (after the inner trust has
## converged at the modes), then added into the outer gradient.
##
## Math identities: envelope theorem for d/d theta at eta = eta*(theta),
## implicit chain via Newton hessian for d eta*/d theta, and the sigma-driven
## contribution when the errfn depends on theta or eta.
.normalFoceiCorrection <- function(full_pars, joint_hessian, fixed,
                                    outer_names, H_inv_list,
                                    prdfn, errfn, omega,
                                    subjects, subject_etas, K, N,
                                    data_per_subject,
                                    conv2 = 2, cores = 1L) {
  Q <- length(outer_names)
  correction0 <- setNames(numeric(Q), outer_names)

  # Shared across subjects (prediction-independent): d Omega^{-1} / d chol,
  # in closed form. buildL() writes one entry per Cholesky parameter -- exp(v)
  # on the diagonal, v off it -- so dL/dc has a single nonzero and
  #   d Omega^-1 / dc = -Omega^-1 (dL L' + L dL') Omega^-1.
  chol_in_outer <- intersect(outer_names, omega$cholPars)
  dOmega_inv_dchol <- list()
  if (length(chol_in_outer) > 0L) {
    L0        <- omega$buildL(full_pars[omega$cholPars])
    Omega_inv <- chol2inv(chol(tcrossprod(L0)))
    for (cp in chol_in_outer) {
      k <- omega$cholLoc[cp, 1L]; l <- omega$cholLoc[cp, 2L]
      dL <- matrix(0, omega$K, omega$K)
      dL[k, l] <- if (isTRUE(omega$isDiag[[cp]])) L0[k, l] else 1
      dOmega <- dL %*% t(L0) + L0 %*% t(dL)
      dOmega_inv_dchol[[cp]] <- -Omega_inv %*% dOmega %*% Omega_inv
      dimnames(dOmega_inv_dchol[[cp]]) <- list(omega$eta, omega$eta)
    }
  }

  # Per-subject Stage-2 contribution: predict the subject, build its error model,
  # accumulate its correction into a local vector, then Reduce() in subject order
  # (bit-identical to a serial running sum). This is the dominant eager-correction
  # cost -- one deriv2 ODE solve per subject -- so fork it over `cores`.
  # Predict each subject only on its own time grid (up to that subject's last
  # observation): integrating short-follow-up subjects out to the global
  # horizon would waste ODE work, and the correction reads only the subject's
  # own data times.
  t_pred <- lapply(subjects, function(cn)
    sort(unique(c(0, data_per_subject[[cn]]$time))))
  # One batched call instead of forking a deriv2 solve per subject: the fork
  # would serialise the batch anyway (cppDE bails out inside a forked child),
  # and the prdframes are expensive to ship back through a pipe.
  preds <- .predictMany(prdfn, times = t_pred,
                        parsList = rep(list(full_pars), N),
                        conditions = subjects, fixed = fixed,
                        deriv = TRUE, deriv2 = TRUE, cores = cores)

  one_subject <- function(i) {
    cn       <- subjects[i]
    t_pred_i <- t_pred[[i]]
    pred     <- structure(list(preds[[i]]), names = cn)
    err_pred <- NULL
    if (!is.null(errfn)) {
      pred_e     <- pred[[cn]]
      pinner     <- getParameters(pred_e)
      fixedinner <- pinner[attr(pinner, "fixed")]
      pinner     <- as.parvec(pinner[setdiff(names(pinner), names(fixed))])
      fixedinner <- as.parvec(fixedinner, deriv = FALSE, deriv2 = FALSE)
      err_pred <- list()
      err_pred[[cn]] <- errfn(out = pred_e, pars = pinner, fixed = fixedinner,
                              conditions = cn)[[cn]]
    }
    correction <- correction0
    eta_i_nm   <- subject_etas[i, ]
    d_i        <- data_per_subject[[cn]]
    pred_i     <- pred[[cn]]
    d_full     <- attr(pred_i, "deriv")
    d2_full    <- attr(pred_i, "deriv2")
    if (is.null(d_full) || is.null(d2_full))
      stop("`prdfn` did not produce 'deriv'/'deriv2' attributes for ",
           "condition '", cn, "'. Did you build it with deriv2-capable ",
           "Xs/Y/P?")

    avail_pars <- dimnames(d_full)[[3]]
    th_avail   <- intersect(outer_names, avail_pars)
    eta_avail  <- intersect(eta_i_nm,    avail_pars)
    if (!length(eta_avail)) return(correction)

    times_i  <- d_i$time
    names_i  <- as.character(d_i$name)
    sigma_i  <- d_i$sigma
    ti_e <- oi_e <- NULL
    if (!is.null(err_pred)) {
      err_i  <- err_pred[[cn]]
      pt     <- as.numeric(err_i[, "time"])
      ti_e   <- .matchNum(times_i, pt)
      ni_e   <- match(names_i, colnames(err_i))
      if (anyNA(ti_e) || anyNA(ni_e))
        stop("compute_correction: cannot align data point to errfn grid.")
      sigma_from_err <- err_i[cbind(ti_e, ni_e)]
      sigma_i <- ifelse(is.na(sigma_i), sigma_from_err, sigma_i)
      d_err  <- attr(err_i, "deriv")
      if (!is.null(d_err)) {
        oi_e <- match(names_i, dimnames(d_err)[[2]])
        if (anyNA(oi_e))
          stop("compute_correction: data name not in errfn deriv axis.")
      }
    }
    time_idx <- match(times_i, t_pred_i)
    if (anyNA(time_idx))
      stop("compute_correction: data times missing from prediction grid.")

    Tn <- length(times_i)
    Kn <- length(eta_avail)
    Qn <- length(th_avail)

    G    <- matrix(0,    Tn, Kn, dimnames = list(NULL, eta_avail))
    Wmix <- array (0, c(Tn, Kn, Qn))
    Veta <- array (0, c(Tn, Kn, Kn))
    for (jr in seq_len(Tn)) {
      ti <- time_idx[jr]; nm <- names_i[jr]
      G[jr, ]      <- d_full [ti, nm, eta_avail]
      Wmix[jr, , ] <- d2_full[ti, nm, eta_avail, th_avail]
      Veta[jr, , ] <- d2_full[ti, nm, eta_avail, eta_avail]
    }
    Sinv <- 1 / sigma_i^2

    H_inv_full <- H_inv_list[[i]]
    H_inv      <- H_inv_full[eta_avail, eta_avail, drop = FALSE]

    HG   <- G %*% H_inv
    qf_t <- rowSums(HG * G)

    if (Qn > 0L) {
      explicit <- numeric(Qn)
      for (k in seq_len(Qn)) {
        Mk <- crossprod(G * Sinv, Wmix[, , k, drop = TRUE])
        Mk <- Mk + t(Mk)
        explicit[k] <- conv2 * sum(H_inv * Mk)
      }
      correction[th_avail] <- correction[th_avail] + explicit
    }

    if (length(chol_in_outer) > 0L) {
      eta_pos  <- match(eta_avail, eta_i_nm)
      eta_base <- omega$eta[eta_pos]
      for (cp in chol_in_outer) {
        dOmega_inv_sub <-
          dOmega_inv_dchol[[cp]][eta_base, eta_base, drop = FALSE]
        correction[cp] <- correction[cp] +
          conv2 * sum(H_inv * dOmega_inv_sub)
      }
    }

    if (!is.null(err_pred)) {
      err_i <- err_pred[[cn]]
      d_err <- attr(err_i, "deriv")
      if (!is.null(d_err)) {
        err_par_set  <- dimnames(d_err)[[3]]
        theta_in_err <- intersect(outer_names, err_par_set)
        for (q in theta_in_err) {
          dsigma_q_t <- d_err[cbind(ti_e, oi_e, match(q, err_par_set))]
          correction[q] <- correction[q] +
            (-4) * sum(dsigma_q_t * qf_t / sigma_i^3)
        }
      }
    }

    th_for_implicit <- intersect(outer_names, colnames(joint_hessian))
    if (length(th_for_implicit) > 0L) {
      cross_block <- joint_hessian[eta_avail, th_for_implicit, drop = FALSE]
      H_newton    <- joint_hessian[eta_avail, eta_avail, drop = FALSE]
      eigN          <- eigen(H_newton, symmetric = TRUE)
      trN           <- sum(eigN$values)
      epsN          <- 1e-10 * abs(trN) / Kn
      lamN          <- pmax(eigN$values, epsN)
      H_inv_Newton  <- eigN$vectors %*% (t(eigN$vectors) / lamN)
      dimnames(H_inv_Newton) <- dimnames(H_newton)

      deta_dth <- -H_inv_Newton %*% cross_block

      dlogH_deta <- numeric(Kn)
      for (l in seq_len(Kn)) {
        Ml <- crossprod(G * Sinv, Veta[, , l, drop = TRUE])
        Ml <- Ml + t(Ml)
        dlogH_deta[l] <- conv2 * sum(H_inv * Ml)
      }
      if (!is.null(err_pred)) {
        err_i <- err_pred[[cn]]
        d_err <- attr(err_i, "deriv")
        if (!is.null(d_err)) {
          err_par_set <- dimnames(d_err)[[3]]
          for (l in seq_len(Kn)) {
            et <- eta_avail[l]
            if (et %in% err_par_set) {
              dsigma_et_t <- d_err[cbind(ti_e, oi_e,
                                          match(et, err_par_set))]
              dlogH_deta[l] <- dlogH_deta[l] +
                (-4) * sum(dsigma_et_t * qf_t / sigma_i^3)
            }
          }
        }
      }
      implicit <- as.numeric(dlogH_deta %*% deta_dth)
      correction[th_for_implicit] <-
        correction[th_for_implicit] + implicit
    }
    correction
  }

  parts <- lapply(seq_len(N), one_subject)
  failed <- vapply(parts,
                   function(p) inherits(p, "try-error") || is.null(p), TRUE)
  if (any(failed))
    stop(".normalFoceiCorrection: subject worker(s) failed: ",
         paste(subjects[failed], collapse = ", "), call. = FALSE)
  Reduce(`+`, parts, correction0)
}



# Batched predictor for the lock-step inner trust. Returns one prediction and
# one error-model frame per request; conditions may repeat across rounds.
.foceiPredictMany <- function(prdfn, errfn) {
  function(timesList, parsList, conditions, fixed) {
    cores <- getOption("dMod.cores", 1L)
    conditions <- as.character(conditions)
    preds <- .predictMany(prdfn, times = timesList, parsList = parsList,
                          conditions = conditions, fixed = fixed,
                          deriv = TRUE, cores = cores)
    n <- length(preds)
    pinner <- lapply(preds, function(q) attr(q, "parameters"))
    est <- .fnNode(errfn)
    errs <- if (!is.null(est)) {
      b <- .bundle(conds = conditions, out = as.list(preds), pars = pinner,
                   fixed = vector("list", n))
      .evalNode(est, b, TRUE, FALSE, NULL, cores)
    } else {
      lapply(seq_len(n), function(i)
        errfn(out = preds[[i]], pars = pinner[[i]],
              conditions = conditions[i])[[1L]])
    }
    list(pred = unname(as.list(preds)), err = unname(as.list(errs)))
  }
}


#' Build a quadrature-based NLME marginal-likelihood objective function
#'
#' Returns a callable `em(pars)` that integrates out the per-subject random
#' effects using a sparse-grid Gauss-Hermite rule. The ECM E-step refreshes
#' the integration nodes between outer iterations via
#' `attr(em, "rebuildQuadrature")`. Most users want [EM], which builds
#' the right `em` internally and runs a solver; use `emObjfn` directly only
#' if you need to evaluate the objective by hand.
#'
#' For the Laplace approximation of the marginal likelihood, use [EM]
#' with `method = "focei"`.
#'
#' @param obj An \code{objfn} of the form
#'   `normL2(data, g*x*p, errmodel = err) + constraintL2(mu = 0, Omega = om)`.
#'   The prediction function, data, error model, and omega spec are recovered
#'   from it automatically.
#' @param control Named list with `level` (Smolyak depth, default 4),
#'   `cores` (default 1), and `pruneTol` (node weight-pruning threshold,
#'   default `Inf` = off).
#'
#' @return A callable of class `c("emobjfn", "objfn", "fn")`. Calling it on
#'   `pars` returns an [objlist] with an `emDiag` attribute carrying
#'   quadrature diagnostics.
#'
#' @seealso [EM], [omega]
#' @examples
#' \dontrun{
#' em <- emObjfn(obj, control = list(level = 5L))
#' attr(em, "rebuildQuadrature")(init)  # E-step: build the grid at `init`
#' em(init)                             # marginal objlist (+ emDiag attribute)
#' }
#' @export
emObjfn <- function(obj, control = list()) {
  rec <- .normalReconstruct(obj)
  .normalQuadratureObjfn(obj, rec$omega,
                      prdfn    = rec$prdfn,
                      data     = rec$data,
                      errfn    = rec$errfn,
                      level    = control$level %||% 4L,
                      cores    = control$cores %||% 1L,
                      pruneTol = control$pruneTol %||% Inf)
}



## Internal quadrature-method emObjfn constructor.
##
## Closure state separates the frozen E-step (nodes_per_subject, etaModes,
## chol_value, current_level) from the trust-varying structural pars. Callers
## refresh the frozen state via attr(em, "rebuildQuadrature")(psiFull, level)
## between outer iterations and run trust(em, init = psi_structural) with the
## integration grid held fixed.
.normalQuadratureObjfn <- function(obj, omega, prdfn, data, errfn,
                                level, cores, pruneTol = Inf) {
  if (!inherits(obj, "objfn"))
    stop("`obj` must be an objfn.")
  if (!inherits(omega, "omegaspec"))
    stop("`omega` must be built by omega().")
  if (is.null(omega$subjectEtas))
    stop("`omega` must have subject expansion (call omega(..., subjects = ...)).")
  if (is.null(prdfn))
    stop("`prdfn` (the prdfn used to build `obj`) is required for ",
         "method = \"quadrature\".")
  if (is.null(data))
    stop("`data` (the datalist used for obj) is required for ",
         "method = \"quadrature\".")

  K            <- omega$K
  subject_etas <- omega$subjectEtas
  subjects     <- rownames(subject_etas)
  N            <- length(subjects)
  all_eta_nm   <- as.vector(subject_etas)
  chol_pars    <- omega$cholPars
  joint_pars   <- attr(obj, "parameters")

  # Frozen E-step state (mutable via <<- inside rebuildQuadrature only).
  nodes_per_subject <- NULL
  eta_modes_state   <- NULL
  H_i_list_state    <- NULL
  chol_value_state  <- NULL
  current_level     <- as.integer(level)
  prune_tol         <- pruneTol
  fast_meta_cache   <- NULL
  n_floored_state   <- integer(N)
  converged_state   <- logical(N)
  iter_state        <- integer(N)

  # E-step: find every subject's posterior mode and eigen-floored Gauss-Newton
  # Hessian in one call to the compiled focei_inner_trust kernel (built lazily
  # from the model on first use), then place the adapted Smolyak nodes on each
  # subject's local Gaussian. This kernel is the sole mode-finder: it bypasses
  # the per-subject R trust / normL2 glue and returns the same Gauss-Newton
  # curvature the R path used, so it is a pure speedup.
  rebuildQuadrature <- function(psiFull, level_new = NULL,
                                fixed = NULL, eta_init = NULL) {
    if (!is.null(level_new)) current_level <<- as.integer(level_new)
    if (!all(chol_pars %in% names(psiFull)))
      stop("rebuildQuadrature: `psiFull` is missing omega$cholPars.")
    outer_with_chol <- psiFull[intersect(names(psiFull),
                                         setdiff(joint_pars, all_eta_nm))]

    if (is.null(eta_init))
      eta_init <- if (!is.null(eta_modes_state)) eta_modes_state
                  else matrix(0, N, K, dimnames = dimnames(subject_etas))

    if (is.null(fast_meta_cache))
      fast_meta_cache <<- .buildBayesSubjectMeta(omega, outer_with_chol,
                                                 prdfn, data, errfn)
    sm <- fast_meta_cache$subjectMeta
    pf <- setNames(numeric(length(sm$pars_full_names)), sm$pars_full_names)
    pf[names(outer_with_chol)] <- outer_with_chol
    for (i in seq_len(N)) pf[subject_etas[i, ]] <- eta_init[i, ]
    L_om <- omega$buildL(psiFull[chol_pars])
    fr <- focei_inner_trust(
      model_cb = prdfn, err_cb = fast_meta_cache$errfn, pars_full = pf,
      eta_warmstart = eta_init, subject_meta = sm,
      Omega_inv_mat = chol2inv(t(L_om)),
      Omega_log_det = 2 * sum(log(diag(L_om))),
      fixed = fixed, control = fast_meta_cache$innerControl)

    new_modes <- matrix(0, N, K, dimnames = dimnames(subject_etas))
    new_H     <- vector("list", N)
    new_nodes <- vector("list", N)
    for (i in seq_len(N)) {
      new_modes[i, ] <- fr$eta_modes[i, ]
      new_H[[i]]     <- fr$H_GN[[i]]
      new_nodes[[i]] <- makeSubjectNodes(fr$eta_modes[i, ], fr$H_GN[[i]],
                                         current_level, pruneTol = prune_tol)
    }

    nodes_per_subject <<- new_nodes
    eta_modes_state   <<- new_modes
    H_i_list_state    <<- new_H
    chol_value_state  <<- psiFull[chol_pars]
    # nFloored is not exposed by the kernel (it floors internally); the mode
    # convergence flag and iteration count are.
    n_floored_state   <<- integer(N)
    converged_state   <<- as.logical(fr$converged)
    iter_state        <<- as.integer(fr$iterations)

    invisible(list(etaModes   = new_modes,
                   HiList     = new_H,
                   level      = current_level,
                   nFloored   = n_floored_state,
                   converged  = converged_state,
                   iterations = iter_state))
  }

  myfn <- function(..., fixed = NULL, deriv = TRUE, env = NULL,
                   cores = getOption("dMod.cores", 1L)) {
    p <- list(...)[[match.fnargs(list(...), "pars")]]
    if (is.null(env)) env <- new.env()

    if (is.null(nodes_per_subject))
      stop("emObjfn: no quadrature grid built. Call ",
           "`attr(em, 'rebuildQuadrature')(psiFull, level)` first.")

    # Build psiFull from p + frozen chol_value_state (CM-1 holds chol fixed).
    structural_names <- setdiff(names(p), c(chol_pars, all_eta_nm))
    chol_in_p <- intersect(chol_pars, names(p))
    psiFull <- c(p[structural_names],
                  if (length(chol_in_p)) p[chol_in_p] else chol_value_state)

    # Outer-active params = those varied by the trust caller = names(p).
    outer_active <- structural_names

    # Serial over subjects: the quadrature nodes go out as one batch inside
    # each subject, and a fork around that would serialise it (cppDE bails out
    # in a forked child) while shipping prdframes back through a pipe.
    per_subj <- lapply(seq_len(N), function(i) .normalEcmSubject(
      i, psiFull, eta_modes_state, omega, nodes_per_subject[[i]],
      xPred = prdfn, datalist = data, errfn = errfn,
      fixed = fixed, outerActiveNames = outer_active,
      mode = if (deriv) "with_grad" else "moments_only",
      cores = cores))

    tot_logL <- 0
    tot_gr   <- setNames(numeric(length(p)), names(p))
    tot_he   <- matrix(0, length(p), length(p),
                       dimnames = list(names(p), names(p)))
    MHatList      <- vector("list", N)
    mHatList      <- vector("list", N)
    max_softmax_per <- numeric(N)
    n_eff_per       <- numeric(N)

    for (i in seq_len(N)) {
      o <- per_subj[[i]]
      tot_logL          <- tot_logL + o$logLhat
      MHatList[[i]]   <- o$M_hat
      mHatList[[i]]   <- o$m_hat
      max_softmax_per[i] <- o$maxSoftmax
      n_eff_per[i]      <- o$n_eff
      if (deriv && !is.null(o$gradient)) {
        nm <- intersect(outer_active, names(tot_gr))
        tot_gr[nm] <- tot_gr[nm] + o$gradient[nm]
        tot_he[nm, nm] <- tot_he[nm, nm] + o$hessian[nm, nm]
      }
    }

    OFV <- -2 * tot_logL
    out <- objlist(value = OFV, gradient = tot_gr, hessian = tot_he)
    emDiag <- list(method          = "quadrature",
                    level           = current_level,
                    etaModes       = eta_modes_state,
                    HiList        = H_i_list_state,
                    mHatList      = mHatList,
                    MHatList      = MHatList,
                    maxSoftmax     = max_softmax_per,
                    n_eff           = n_eff_per,
                    nFloored       = n_floored_state,
                    innerConverged = converged_state,
                    innerIter      = iter_state)
    attr(out, "emDiag") <- emDiag
    attr(out, "env")     <- env
    out
  }

  class(myfn) <- c("emobjfn", "objfn", "fn")
  attr(myfn, "method")     <- "quadrature"
  attr(myfn, "conditions") <- attr(obj, "conditions")
  attr(myfn, "parameters") <- setdiff(joint_pars, all_eta_nm)
  attr(myfn, "omega")  <- omega
  attr(myfn, "obj")      <- obj
  attr(myfn, "prdfn")      <- prdfn
  attr(myfn, "data")       <- data
  attr(myfn, "errfn")   <- errfn
  attr(myfn, "rebuildQuadrature") <- rebuildQuadrature
  attr(myfn, "control")    <- list(level = level, cores = cores)
  myfn
}



# .fitNormal S3 constructor. Bundles solver output with the prdfn / data /
# omega references that predict.em and the diagnostic plot
# helpers (in plots.R) consume.
.normalFitMake <- function(argument, value, gradient, hessian, Omega, etaModes,
                         converged, iterations, emDiag, method,
                         foceiStart = NULL, stageTrace = NULL,
                         prdfn = NULL, data = NULL, omega = NULL,
                         errfn = NULL) {
  etaInfo  <- .normalEtaInfo(emDiag, Omega, etaModes)
  # OFV convention (identical across the FOCEI and quadrature paths): value is
  # the plain-ML marginal -2 log L with every 2*pi constant retained. This
  # equals nlmixr2's -2LL and NONMEM's
  # "OFV with constant". NONMEM's default OBJ drops the data-side Sum log(2*pi),
  # so value_nonmem re-derives that raw-.lst-comparable number.
  n_obs <- .normalNObs(data)
  value_nonmem <- if (!is.na(n_obs) && !is.null(value))
                    value - n_obs * log(2 * pi) else NA_real_
  out <- list(argument     = argument,
              value        = value,
              ofvType      = "-2LL",
              nObs         = n_obs,
              value_nonmem = value_nonmem,
              gradient     = gradient,
              hessian      = hessian,
              Omega        = Omega,
              etaModes     = etaModes,
              etaSE        = etaInfo$etaSE,
              shrinkage    = etaInfo$shrinkage,
              converged    = converged,
              iterations   = iterations,
              emDiag       = emDiag,
              method       = method,
              prior        = "omega",
              foceiStart   = foceiStart,
              stageTrace   = stageTrace,
              prdfn        = prdfn,
              data         = data,
              omega        = omega,
              errfn        = errfn)
  class(out) <- c("em", "list")
  out
}

# Total number of observations (likelihood-contributing data rows) across a
# datalist, used to convert the dMod plain-ML OFV into a NONMEM-OBJ-comparable
# value. Returns NA when the data are unavailable.
.normalNObs <- function(data) {
  if (is.null(data)) return(NA_integer_)
  n <- tryCatch(sum(vapply(data, function(d) nrow(d), integer(1))),
                error = function(e) NA_integer_)
  as.integer(n)
}

# Posterior-mode standard errors and shrinkage diagnostics for the per-subject
# random effects. Caller passes the full emDiag (which carries HInvList from
# either the R or C++ Laplace path); returns NULLs when the inverse Hessian
# list is unavailable (e.g. quadrature method).
.normalEtaInfo <- function(emDiag, Omega, etaModes) {
  out <- list(etaSE = NULL, shrinkage = NULL)
  if (is.null(emDiag) || is.null(emDiag$HInvList) || is.null(etaModes))
    return(out)
  H_inv_list <- emDiag$HInvList
  N <- nrow(etaModes); K <- ncol(etaModes)
  if (length(H_inv_list) != N) return(out)

  diag_template <- rep(NA_real_, K)
  diags <- vapply(H_inv_list, function(Hi) {
    if (is.null(Hi) || !is.matrix(Hi) || any(dim(Hi) != c(K, K)))
      return(diag_template)
    di <- diag(Hi)
    di[di < 0] <- NA_real_
    di
  }, numeric(K))
  dim(diags) <- c(K, N)
  etaSE <- sqrt(t(diags))
  dimnames(etaSE) <- dimnames(etaModes)
  out$etaSE <- etaSE

  if (!is.null(Omega) && all(dim(Omega) == c(K, K))) {
    Omega_sd <- sqrt(pmax(diag(Omega), 0))
    sd_mat <- matrix(Omega_sd, N, K, byrow = TRUE)
    shrink <- 1 - etaSE / sd_mat
    shrink[, Omega_sd <= 0] <- NA_real_
    dimnames(shrink) <- dimnames(etaModes)
    out$shrinkage <- shrink
  }
  out
}


# Synthetic err callable for the case where the user supplied no errfn
# and instead set `data$sigma` directly. Returns a prdlist whose matrix has
# the per-observation sigmas (from the data) padded onto the prediction's
# time grid. No `deriv` attribute, so the C++ kernel treats `dsigma/deta = 0`.
.normalStaticErr <- function(data) {
  data_per_subject <- lapply(seq_along(data), function(i) {
    d <- data[[i]]
    d$name <- as.character(d$name)
    d
  })
  names(data_per_subject) <- names(data)

  function(out, pars = NULL, conditions = NULL, ...) {
    s <- if (!is.null(conditions)) conditions[1] else names(data_per_subject)[1]
    out_mat <- if (inherits(out, "prdlist")) out[[s]] else out
    pred_times <- as.numeric(out_mat[, "time"])
    obs_names  <- setdiff(colnames(out_mat), "time")
    d <- data_per_subject[[s]]
    if (is.null(d))
      stop(".normalStaticErr: condition '", s,
           "' missing from data; supply an errfn.")
    if (anyNA(d$sigma))
      stop(".normalStaticErr: condition '", s,
           "' has NA in data$sigma; supply an errfn.")
    err_mat <- matrix(NA_real_, length(pred_times), length(obs_names) + 1L)
    colnames(err_mat) <- c("time", obs_names)
    err_mat[, "time"] <- pred_times
    for (o in obs_names) {
      drow <- d[d$name == o, , drop = FALSE]
      if (!nrow(drow)) next
      idx <- match(pred_times, drow$time)
      err_mat[, o] <- ifelse(is.na(idx), median(drow$sigma), drow$sigma[idx])
    }
    setNames(list(structure(err_mat, class = c("prdframe", "matrix"))), s)
  }
}


# Builds the long-format subject_meta consumed by focei_inner_trust / the
# C++ kernel. One row in t_idx_in_pred / o_idx_in_pred / y_data / ... per
# observed data point of the subject, covering all observables. eta_idx_in_*
# arrays are length K (one entry per random effect). Values 0 mark "no
# contribution" (e.g. an eta does not appear in the err prdfn's deriv).
.normalFastMeta <- function(prdfn, errfn, data, subjects,
                           eta_names_list, pars_full_names, pars_probe) {
  N <- length(subjects)
  fast_meta <- vector("list", N)
  for (i in seq_len(N)) {
    s      <- subjects[i]
    data_i <- data[[s]]
    data_i$name <- as.character(data_i$name)
    data_i <- data_i[order(data_i$time, data_i$name), , drop = FALSE]
    times_union <- sort(unique(data_i$time))

    # Probe prdfn + err once to learn deriv array structure.
    pred_probe <- prdfn(times = times_union, pars = pars_probe,
                        deriv = TRUE, conditions = s)
    prdf  <- pred_probe[[1]]
    pcols <- colnames(prdf)
    d_dn  <- dimnames(attr(prdf, "deriv"))
    eta_names_i      <- eta_names_list[[i]]
    eta_idx_in_deriv <- match(eta_names_i, d_dn[[3]])
    if (anyNA(eta_idx_in_deriv))
      stop(".normalFastMeta: eta names not present in prdfn deriv array for ",
           "subject ", s, ".")

    pinner    <- attr(prdf, "parameters")
    err_probe <- errfn(out = prdf, pars = pinner, conditions = s)
    erm       <- err_probe[[1]]
    ecols     <- colnames(erm)
    e_attr    <- attr(erm, "deriv")
    if (!is.null(e_attr)) {
      e_dn               <- dimnames(e_attr)
      eta_idx_in_err_deriv <- match(eta_names_i, e_dn[[3]])
      eta_idx_in_err_deriv[is.na(eta_idx_in_err_deriv)] <- 0L
    } else {
      e_dn               <- NULL
      eta_idx_in_err_deriv <- rep(0L, length(eta_names_i))
    }

    # Long-format per-row indices.
    t_idx_in_pred  <- match(data_i$time, prdf[, "time"])
    o_idx_in_pred  <- match(data_i$name, pcols)
    if (anyNA(o_idx_in_pred))
      stop(".normalFastMeta: observable(s) not found in prdfn output for ",
           "subject ", s, ": ",
           paste(setdiff(unique(data_i$name), pcols), collapse = ", "))
    o_idx_in_deriv <- match(data_i$name, d_dn[[2]])
    t_idx_in_err   <- match(data_i$time, erm[, "time"])
    o_idx_in_err   <- match(data_i$name, ecols)
    if (!is.null(e_dn)) {
      o_idx_in_err_deriv <- match(data_i$name, e_dn[[2]])
      o_idx_in_err_deriv[is.na(o_idx_in_err_deriv)] <- 0L
    } else {
      o_idx_in_err_deriv <- rep(0L, nrow(data_i))
    }

    fast_meta[[i]] <- list(
      times                = as.numeric(times_union),
      eta_idx_in_pars      = as.integer(match(eta_names_i, pars_full_names)),
      t_idx_in_pred        = as.integer(t_idx_in_pred),
      y_data               = as.numeric(data_i$value),
      o_idx_in_pred        = as.integer(o_idx_in_pred),
      eta_idx_in_deriv     = as.integer(eta_idx_in_deriv),
      o_idx_in_deriv       = as.integer(o_idx_in_deriv),
      t_idx_in_err         = as.integer(t_idx_in_err),
      o_idx_in_err         = as.integer(o_idx_in_err),
      eta_idx_in_err_deriv = as.integer(eta_idx_in_err_deriv),
      o_idx_in_err_deriv   = as.integer(o_idx_in_err_deriv),
      condition            = s,
      eta_names            = eta_names_i)
  }
  fast_meta
}


# Internal: run the C++ FOCEI kernel and package its output as .fitNormal.
# Used by .fitNormal(method = "focei"). Always uses the fast-inner C++ path
# with eager Stage-2 correction, calling .normalFoceiCorrection as an
# Rcpp::Function once per outer iter.
.normalFocei <- function(obj, omega, init, prdfn, data, errfn,
                         fixed = NULL,
                         innerControl = list(), trustControl = list(),
                         cores = 1L, secondOrderCorrection = FALSE,
                         methodLabel = "focei") {
  if (is.null(prdfn))
    stop(".normalFocei: `prdfn` (the prdfn) is required.")
  if (is.null(data))
    stop(".normalFocei: `data` (the datalist) is required.")
  if (is.null(omega$subjectEtas))
    stop(".normalFocei: omega has no subject expansion. Call ",
         "omega(..., subjects = ...).")

  # When the user did not supply an errfn but recorded sigma in the data
  # itself, wrap the per-row sigmas into a synthetic obsfn-like callable.
  # The fast-inner kernel treats this as "sigma constant in eta" (no err
  # deriv attribute -> Js = 0).
  if (is.null(errfn)) errfn <- .normalStaticErr(data)

  K        <- omega$K
  subjects <- rownames(omega$subjectEtas)
  N        <- length(subjects)
  outer_names <- names(init)
  eta_names_all <- as.vector(omega$subjectEtas)
  pars_full_names <- c(outer_names, eta_names_all)
  eta_idx_global <- matrix(0L, N, K)
  eta_names_list <- vector("list", N)
  for (i in seq_len(N)) {
    nms <- omega$subjectEtas[i, ]
    eta_names_list[[i]] <- nms
    eta_idx_global[i, ] <- match(nms, pars_full_names)
  }
  outer_idx_full  <- match(outer_names, pars_full_names)
  other_etas_init <- setNames(rep(0, length(eta_names_all)), eta_names_all)
  subject_meta <- list(
    subjects          = subjects,
    eta_idx_global    = eta_idx_global,
    eta_names         = eta_names_list,
    K                 = K,
    outer_names       = outer_names,
    outer_idx_in_full = outer_idx_full,
    other_etas_init   = other_etas_init,
    pars_full_names   = pars_full_names
  )
  ic <- modifyList(list(rinit = 1, rmax = 10, iterlim = 30,
                        ftol = 1e-7, mtol = 1e-7,
                        eigen_floor_relative = 1e-10),
                   innerControl)
  oc <- modifyList(list(rinit = 1, rmax = 10, iterlim = 200,
                        ftol = 1e-7, mtol = 1e-7),
                   trustControl)
  control_cpp <- list(inner = ic, outer = oc)

  pars_probe <- setNames(numeric(length(pars_full_names)), pars_full_names)
  pars_probe[outer_names] <- init
  fast_meta <- .normalFastMeta(prdfn, errfn, data, subjects,
                              eta_names_list, pars_full_names, pars_probe)
  om_meta <- list(
    chol_pars = omega$cholPars,
    chol_loc  = matrix(as.integer(omega$cholLoc), ncol = 2,
                       dimnames = NULL),
    is_diag   = as.logical(omega$isDiag))
  subject_meta$fast_meta  <- fast_meta
  subject_meta$omega_meta <- om_meta

  data_per_subject <- lapply(subjects, function(s) data[[s]])
  names(data_per_subject) <- subjects
  correction_cb <- function(full_pars, joint_hessian, H_inv_list) {
    .normalFoceiCorrection(
      full_pars = full_pars, joint_hessian = joint_hessian,
      fixed = fixed, outer_names = outer_names, H_inv_list = H_inv_list,
      prdfn = prdfn, errfn = errfn, omega = omega,
      subjects = subjects, subject_etas = omega$subjectEtas,
      K = K, N = N,
      data_per_subject = data_per_subject,
      cores = cores)
  }

  # Lock-step inner trust: one batched prediction per round instead of a
  # single-condition callback per subject. Bit-identical to the serial path,
  # so it is the default; set the option or control entry to FALSE to compare.
  lockstep <- !isFALSE(control_cpp$lockstep) &&
    !isFALSE(getOption("dMod.focei.lockstep", TRUE))

  # Stage-2 volume correction. The only part of FOCEI that needs second
  # derivatives of the prediction -- the inner mode, the Gauss-Newton Hessian
  # and the Schur block are all products of first derivatives. Off by default,
  # so a chain built without deriv2 fits out of the box.
  correction <- isTRUE(secondOrderCorrection)

  fit <- focei_run(model_cb = prdfn, err_cb = errfn,
                        joint_cb = obj, init = init,
                        subject_meta = subject_meta,
                        fixed = fixed, control = control_cpp,
                        correction_mode = if (correction) "eager" else "none",
                        correction_cb = if (correction) correction_cb else NULL,
                        predict_cb = if (lockstep)
                          .foceiPredictMany(prdfn, errfn) else NULL)

  L_omega <- if (all(omega$cholPars %in% names(fit$argument)))
    omega$buildL(fit$argument[omega$cholPars]) else NULL
  Omega <- if (!is.null(L_omega)) tcrossprod(L_omega) else NULL
  etaModes <- fit$etaModes
  rownames(etaModes) <- subjects
  colnames(etaModes) <- omega$eta
  emDiag <- list(etaStar = etaModes, logdet = fit$log_det_H,
                 sum_logdetH = fit$sum_logdetH, trace = fit$trace,
                 backend = "cpp", HInvList = fit$H_inv)
  .normalFitMake(argument = fit$argument, value = fit$value,
               gradient = fit$gradient, hessian = fit$hessian,
               Omega = Omega, etaModes = etaModes,
               converged = isTRUE(fit$converged),
               iterations = fit$iterations, emDiag = emDiag,
               method = methodLabel,
               prdfn = prdfn, data = data,
               omega = omega, errfn = errfn)
}


# Internal: run ECM (E-step / CM-1 / CM-2) polish on a quadrature-method
# emObjfn, package as .fitNormal.
.normalQuadrature <- function(em, init, fixed = NULL, foceiStart = NULL,
                              epsQuadLevels = NULL, epsEcm = 1e-4,
                              epsOfvRel = 1e-5, maxEcmPerStage = 5L,
                              maxCm1Iter = 30L, cm1Control = list(),
                              methodLabel = "quadrature", verbose = TRUE) {
  om <- attr(em, "omega")
  K  <- om$K
  chol_pars <- om$cholPars
  if (is.null(epsQuadLevels)) epsQuadLevels <- K + 1:3
  cm1 <- .trustControl(list(rinit = 1, rmax = 10, ftol = 1e-6, mtol = 1e-6),
                       cm1Control, label = "cm1Control")
  psi <- init
  if (!all(chol_pars %in% names(psi)))
    stop(".fitNormal: `init` is missing omega$cholPars (",
         paste(setdiff(chol_pars, names(psi)), collapse = ", "), ").")
  structural_names <- setdiff(names(psi), chol_pars)
  rebuild <- attr(em, "rebuildQuadrature")
  stage_rows <- list(); prev_psi <- psi; prev_ofv <- NA_real_

  for (stage in seq_along(epsQuadLevels)) {
    level <- epsQuadLevels[stage]
    if (verbose) message(sprintf("EM(%s): stage %d / %d (level = %d)",
                                 methodLabel, stage,
                                 length(epsQuadLevels), level))
    e_info <- rebuild(psi, level_new = level, fixed = fixed,
                      eta_init = if (!is.null(foceiStart))
                                   foceiStart$emDiag$etaStar else NULL)
    for (ecmIter in seq_len(maxEcmPerStage)) {
      if (ecmIter > 1L)
        e_info <- rebuild(psi, level_new = level, fixed = fixed,
                          eta_init = e_info$etaModes)
      cm1_fit <- suppressMessages(do.call(trust, modifyList(cm1, list(
        objfun = em, parinit = psi[structural_names], fixed = fixed,
        iterlim = maxCm1Iter))))
      psi[structural_names] <- cm1_fit$argument
      out_after_cm1 <- em(psi[structural_names], fixed = fixed, deriv = FALSE)
      diag_after    <- attr(out_after_cm1, "emDiag")
      psi[chol_pars] <- updateOmegaChol(diag_after$MHatList, om)
      ofv         <- out_after_cm1$value
      deltaPsi    <- max(abs(psi - prev_psi))
      deltaOfvRel <- if (is.na(prev_ofv) || abs(prev_ofv) < .Machine$double.eps)
                       Inf else abs(ofv - prev_ofv) / abs(prev_ofv)
      stage_rows[[length(stage_rows) + 1L]] <- data.frame(
        stage = stage, ecmIter = ecmIter, level = level,
        OFV = ofv, deltaPsi = deltaPsi, deltaOfvRel = deltaOfvRel,
        maxSoftmax = max(diag_after$maxSoftmax),
        nEffMin = min(diag_after$n_eff),
        cm1TrustIter = cm1_fit$iterations)
      if (verbose) message(sprintf(
        "  ecm %d : OFV=%.6f  |dpsi|=%.2e  |dOFV/OFV|=%.2e  max_smax=%.3f  nEffMin=%.1f",
        ecmIter, ofv, deltaPsi, deltaOfvRel,
        max(diag_after$maxSoftmax), min(diag_after$n_eff)))
      prev_psi <- psi; prev_ofv <- ofv
      if (deltaPsi < epsEcm || deltaOfvRel < epsOfvRel) break
    }
  }

  rebuild(psi, fixed = fixed,
          eta_init = if (length(stage_rows))
                       attr(em(psi[structural_names], deriv = FALSE),
                            "emDiag")$etaModes else NULL)
  final_out  <- em(psi[structural_names], fixed = fixed)
  final_diag <- attr(final_out, "emDiag")
  L_omega    <- om$buildL(psi[chol_pars])
  Omega      <- tcrossprod(L_omega)
  conv       <- length(stage_rows) > 0L && {
    last <- tail(stage_rows, 1)[[1]]
    last$deltaPsi < epsEcm || last$deltaOfvRel < epsOfvRel
  }
  .normalFitMake(argument   = psi,
               value      = final_out$value,
               gradient   = final_out$gradient,
               hessian    = final_out$hessian,
               Omega      = Omega,
               etaModes   = final_diag$etaModes,
               converged  = conv,
               iterations = length(stage_rows),
               emDiag     = final_diag,
               method     = methodLabel,
               foceiStart = foceiStart,
               stageTrace = do.call(rbind, stage_rows),
               prdfn      = attr(em, "prdfn"),
               data       = attr(em, "data"),
               omega      = attr(em, "omega"),
               errfn      = attr(em, "errfn"))
}


# Deterministic Laplace-EM (ECM) estimator. Reuses the quadrature machinery at
# the single-node level L = K, where the adaptive Smolyak rule collapses to the
# Laplace approximation. Per ECM iteration:
#   E-step  : rebuild per-subject modes eta_i* and eigen-floored H_i (freezes
#             the current Omega through psi[chol_pars]).
#   CM-1    : trust() on the (Laplace) marginal `em` over the structural pars.
#   CM-2    : closed-form Omega from the *covariance-corrected* second moments
#             M_i = eta_i* eta_i*^T + H_i^{-1}, via updateOmegaChol(). The
#             H_i^{-1} term is exactly what separates Laplace-EM from ITS
#             (which drops it and systematically under-estimates Omega).
# Same CM-1/CM-2 skeleton and .fitNormal packaging as .normalQuadrature; only the
# fixed level and the covariance-corrected CM-2 moments differ.
.normalLaplaceEM <- function(em, init, fixed = NULL,
                          epsEcm = 1e-4, epsOfvRel = 1e-5,
                          maxEcm = 200L, maxCm1Iter = 30L,
                          cm1Control = list(), verbose = TRUE) {
  om <- attr(em, "omega")
  K  <- om$K
  N  <- nrow(om$subjectEtas)
  chol_pars <- om$cholPars
  level <- K                       # single-node Smolyak == Laplace
  cm1 <- .trustControl(list(rinit = 1, rmax = 10, ftol = 1e-6, mtol = 1e-6),
                       cm1Control, label = "cm1Control")
  psi <- init
  if (!all(chol_pars %in% names(psi)))
    stop(".fitNormal: `init` is missing omega$cholPars (",
         paste(setdiff(chol_pars, names(psi)), collapse = ", "), ").")
  structural_names <- setdiff(names(psi), chol_pars)
  rebuild <- attr(em, "rebuildQuadrature")

  rows <- list(); prev_psi <- psi; prev_ofv <- NA_real_; conv <- FALSE
  e_info <- rebuild(psi, level_new = level, fixed = fixed)
  for (it in seq_len(maxEcm)) {
    if (it > 1L)
      e_info <- rebuild(psi, level_new = level, fixed = fixed,
                        eta_init = e_info$etaModes)

    # CM-1: structural pars via the Laplace marginal, Omega frozen.
    cm1_fit <- suppressMessages(do.call(trust, modifyList(cm1, list(
      objfun = em, parinit = psi[structural_names], fixed = fixed,
      iterlim = maxCm1Iter))))
    psi[structural_names] <- cm1_fit$argument

    # CM-2: closed-form Omega from covariance-corrected posterior moments.
    M_list <- lapply(seq_len(N), function(i)
      tcrossprod(e_info$etaModes[i, ]) + solve(e_info$HiList[[i]]))
    psi[chol_pars] <- updateOmegaChol(M_list, om)

    out <- em(psi[structural_names], fixed = fixed, deriv = FALSE)
    ofv <- out$value
    deltaPsi    <- max(abs(psi - prev_psi))
    deltaOfvRel <- if (is.na(prev_ofv) || abs(prev_ofv) < .Machine$double.eps)
                     Inf else abs(ofv - prev_ofv) / abs(prev_ofv)
    rows[[length(rows) + 1L]] <- data.frame(
      ecmIter = it, level = level, OFV = ofv,
      deltaPsi = deltaPsi, deltaOfvRel = deltaOfvRel,
      cm1TrustIter = cm1_fit$iterations)
    if (verbose) message(sprintf(
      "EM(laplaceEM) ecm %d : OFV=%.6f  |dpsi|=%.2e  |dOFV/OFV|=%.2e",
      it, ofv, deltaPsi, deltaOfvRel))
    prev_psi <- psi; prev_ofv <- ofv
    if (deltaPsi < epsEcm || deltaOfvRel < epsOfvRel) { conv <- TRUE; break }
  }

  # Report an accurate marginal -2 log L at the converged point: the single-node
  # (level-K) Laplace value used to drive the ECM is under-resolved and can sit
  # a few units above the true marginal, so evaluate the final value / gradient /
  # Hessian at level K+2 (which matches the FOCEI marginal). The *estimate* is
  # still the Laplace-EM fixed point; only the reported OFV is refined.
  reportLevel <- as.integer(K + 2L)
  final_e   <- rebuild(psi, level_new = reportLevel, fixed = fixed,
                       eta_init = e_info$etaModes)
  final_out <- em(psi[structural_names], fixed = fixed)
  final_diag <- attr(final_out, "emDiag")
  Omega     <- tcrossprod(om$buildL(psi[chol_pars]))
  .normalFitMake(argument   = psi,
               value      = final_out$value,
               gradient   = final_out$gradient,
               hessian    = final_out$hessian,
               Omega      = Omega,
               etaModes   = final_e$etaModes,
               converged  = conv,
               iterations = length(rows),
               emDiag     = final_diag,
               method     = "laplaceEM",
               foceiStart = NULL,
               stageTrace = do.call(rbind, rows),
               prdfn      = attr(em, "prdfn"),
               data       = attr(em, "data"),
               omega      = attr(em, "omega"),
               errfn      = attr(em, "errfn"))
}


# Stochastic-approximation EM (SAEM). Reuses the Bayesian per-subject eta
# sampler (.buildBayesSubjectMeta + .makeSubjectEtaObj + .run_single_chain) for
# a stochastic E-step, and updateOmegaChol for the closed-form Omega M-step.
# Per iteration k:
#   E-step (stochastic): draw eta_i ~ p(eta_i | y_i, theta, Omega) with a short
#     random-walk-Metropolis chain (nMcmc steps), one draw per subject.
#   SA    : Robbins-Monro update of the sufficient statistic
#     S <- S + gamma_k (1/N sum_i eta_i eta_i^T - S), gamma_k = 1 during the
#     exploration phase (k <= nBurnin) then 1/(k - nBurnin) for convergence.
#   CM-2  : Omega <- S via updateOmegaChol() (closed form).
#   CM-1  : structural theta via a trust() on the complete-data objective at the
#     drawn etas, RM-damped in the convergence phase.
# The final value/gradient/Hessian/etaModes are reported from an accurate
# level-(K+2) quadrature marginal at the converged point (FOCEI-comparable),
# exactly as .normalLaplaceEM does.
.normalSaem <- function(obj, omega, init, prdfn, data, errfn, fixed = NULL,
                     nBurnin = 200L, nEM = 200L, nMcmc = 10L,
                     cm1Control = list(), cores = 1L, verbose = TRUE) {
  if (is.null(errfn)) errfn <- .normalStaticErr(data)
  meta_pkg <- .buildBayesSubjectMeta(omega, init, prdfn, data, errfn)
  meta <- meta_pkg$subjectMeta
  N <- meta_pkg$N; K <- meta_pkg$K
  chol_pars <- omega$cholPars
  structural_names <- setdiff(names(init), chol_pars)
  cm1 <- .trustControl(list(rinit = 1, rmax = 10, iterlim = 30,
                            ftol = 1e-6, mtol = 1e-6),
                       cm1Control, label = "cm1Control")

  parsFull <- setNames(numeric(length(meta$pars_full_names)),
                       meta$pars_full_names)
  parsFull[names(init)] <- init
  subjEtaObjList <- lapply(seq_len(N), function(i)
    .makeSubjectEtaObj(i, meta, prdfn, errfn, parsFull, meta$eta_idx_global[i, ]))

  # CM-1 objective: structural pars at the frozen (sampled) etas + current chol.
  cm1_obj <- function(theta_struct, ...) {
    pf <- parsFull; pf[structural_names] <- theta_struct
    out <- obj(pars = pf, fixed = fixed, deriv = TRUE)
    objlist(value    = out$value,
            gradient = out$gradient[structural_names],
            hessian  = out$hessian[structural_names, structural_names,
                                   drop = FALSE])
  }

  eta_step <- rep(0.4, N)
  S_omega  <- NULL
  rows <- list(); prev_struct <- parsFull[structural_names]
  total <- as.integer(nBurnin + nEM)
  for (k in seq_len(total)) {
    gamma <- if (k <= nBurnin) 1.0 else 1.0 / (k - nBurnin)

    # E-step: one random-walk-Metropolis draw of eta_i per subject.
    L_om          <- omega$buildL(parsFull[chol_pars])
    Omega_inv     <- chol2inv(t(L_om))
    Omega_log_det <- 2 * sum(log(diag(L_om)))
    draw_one <- function(i) {
      bake <- .bake_objfun(subjEtaObjList[[i]], dots = list(),
                           extra = list(.pars_full = parsFull,
                                        .Omega_inv = Omega_inv,
                                        .Omega_log_det = Omega_log_det))
      raw <- tryCatch(.run_single_chain(
        bake, parsFull[meta$eta_names[[i]]], n = nMcmc, warmup = 0L,
        moveType = "mh",
        moveControl = list(stepsize = eta_step[i], proposalCov = "identity"),
        metricControl = list(), bounds = list(upper = rep(Inf, K),
                                              lower = rep(-Inf, K)),
        parscale = rep(1, K), dG_cb = NULL), error = function(e) NULL)
      if (is.null(raw)) return(list(eta = parsFull[meta$eta_names[[i]]],
                                    accept = 0))
      list(eta = raw$samples[nMcmc, ], accept = mean(raw$accept))
    }
    res <- if (cores > 1L)
             parallel::mclapply(seq_len(N), draw_one, mc.cores = cores)
           else lapply(seq_len(N), draw_one)

    Msum <- matrix(0, K, K)
    for (i in seq_len(N)) {
      ei <- as.numeric(res[[i]]$eta)
      parsFull[meta$eta_names[[i]]] <- ei
      Msum <- Msum + tcrossprod(ei)
      # Smooth Robbins-Monro stepsize adaptation toward a 0.4 acceptance target.
      # (A crude divide/multiply rule collapses the step under the noisy
      # few-sample acceptance estimate; the smooth log-scale rule is stable.)
      a <- res[[i]]$accept
      if (is.finite(a))
        eta_step[i] <- min(max(eta_step[i] * exp(0.3 * (a - 0.4)), 0.02), 5)
    }

    # Robbins-Monro sufficient statistic, then closed-form Omega (CM-2).
    # updateOmegaChol averages its list, so a single-element list already
    # yields S = S_omega without building N copies.
    S_new   <- Msum / N
    S_omega <- if (is.null(S_omega)) S_new else S_omega + gamma * (S_new - S_omega)
    parsFull[chol_pars] <- updateOmegaChol(list(S_omega), omega)

    # CM-1: structural M-step at the drawn etas, RM-damped in convergence phase.
    cm1_fit <- suppressMessages(do.call(trust, modifyList(cm1, list(
      objfun = cm1_obj, parinit = parsFull[structural_names]))))
    theta_hat <- cm1_fit$argument
    parsFull[structural_names] <-
      if (k <= nBurnin) theta_hat
      else parsFull[structural_names] + gamma * (theta_hat - parsFull[structural_names])

    deltaStruct <- max(abs(parsFull[structural_names] - prev_struct))
    rows[[length(rows) + 1L]] <- data.frame(
      iter = k, phase = if (k <= nBurnin) "burnin" else "converge",
      gamma = gamma, deltaStruct = deltaStruct)
    if (verbose && (k %% 50L == 0L || k == total))
      message(sprintf("EM(saem) iter %d/%d (%s) |dtheta|=%.2e",
                      k, total, if (k <= nBurnin) "burnin" else "converge",
                      deltaStruct))
    prev_struct <- parsFull[structural_names]
  }

  # Accurate marginal at the converged point (level K+2), FOCEI-comparable.
  em <- emObjfn(obj, control = list(cores = cores))
  psi <- c(parsFull[structural_names], parsFull[chol_pars])
  attr(em, "rebuildQuadrature")(psi, level_new = as.integer(K + 2L), fixed = fixed)
  final_out  <- em(psi[structural_names], fixed = fixed)
  final_diag <- attr(final_out, "emDiag")
  Omega      <- tcrossprod(omega$buildL(psi[chol_pars]))
  # Convergence heuristic: small structural drift over the last convergence iters.
  tail_rows <- tail(do.call(rbind, rows), min(20L, nEM))
  conv <- nEM > 0L && mean(tail_rows$deltaStruct) < 1e-3

  .normalFitMake(argument   = psi,
               value      = final_out$value,
               gradient   = final_out$gradient,
               hessian    = final_out$hessian,
               Omega      = Omega,
               etaModes   = final_diag$etaModes,
               converged  = conv,
               iterations = total,
               emDiag     = final_diag,
               method     = "saem",
               foceiStart = NULL,
               stageTrace = do.call(rbind, rows),
               prdfn      = prdfn,
               data       = data,
               omega      = omega,
               errfn      = errfn)
}


# Recover the model pieces (prdfn, data, errfn, omegaSpec) from a composed
# objective. normL2() and constraintL2_mvn() stamp these as attributes at
# construction and +.objfn coalesces them, so a well-formed NLME objective
# `normL2(data, g*x*p, errmodel = e) + constraintL2(mu = 0, Omega = om)`
# self-describes and callers never re-pass the pieces.
.normalReconstruct <- function(obj) {
  if (!inherits(obj, "objfn"))
    stop(".fitNormal: `obj` must be an objfn.", call. = FALSE)
  prdfn <- attr(obj, "prdfn", exact = TRUE)
  data  <- attr(obj, "data",  exact = TRUE)
  errfn <- attr(obj, "errfn", exact = TRUE)
  omega <- attr(obj, "omegaSpec", exact = TRUE)
  if (is.null(prdfn) || is.null(data))
    stop(".fitNormal: could not recover the prediction function and data from ",
         "`obj`. Build it as ",
         "normL2(data, g*x*p, errmodel = e) + constraintL2(mu = 0, Omega = om).",
         call. = FALSE)
  if (is.null(omega))
    stop(".fitNormal: `obj` carries no random-effects prior. Add ",
         "+ constraintL2(mu = 0, Omega = omega(..., subjects = ...)).",
         call. = FALSE)
  if (is.null(omega$subjectEtas))
    stop(".fitNormal: the omega in `obj` has no subject expansion. Build it with ",
         "omega(..., subjects = ...).", call. = FALSE)
  list(prdfn = prdfn, data = data, errfn = errfn, omega = omega)
}

# Warn on unrecognised control keys instead of silently ignoring them.
.normalCheckControlKeys <- function(user, known, what) {
  if (length(user) == 0L || is.null(names(user))) return(invisible())
  unknown <- setdiff(names(user), known)
  if (length(unknown))
    warning(sprintf(".fitNormal: unrecognised %s control key(s): %s. Recognised: %s.",
                    what, paste(unknown, collapse = ", "),
                    paste(known, collapse = ", ")), call. = FALSE)
  invisible()
}

.normalValidateControl <- function(control) {
  .normalCheckControlKeys(control, c("focei", "quadrature", "saem"), "top-level")
  .normalCheckControlKeys(control$saem %||% list(),
                    c("nBurnin", "nEM", "nMcmc", "cm1Control", "cores"),
                    "saem")
  fc <- control$focei %||% list()
  .normalCheckControlKeys(fc, c("innerControl", "trustControl", "cores"), "focei")
  .normalCheckControlKeys(fc$innerControl %||% list(),
                    c("rinit", "rmax", "iterlim", "ftol", "mtol",
                      "fterm", "mterm", "eigen_floor_relative"),
                    "focei$innerControl")
  .normalCheckControlKeys(fc$trustControl %||% list(),
                    c("rinit", "rmax", "iterlim", "ftol", "mtol",
                      "fterm", "mterm"),
                    "focei$trustControl")
  .normalCheckControlKeys(control$quadrature %||% list(),
                    c("level", "cores", "epsQuadLevels", "epsEcm", "epsOfvRel",
                      "maxEcmPerStage", "maxEcm", "maxCm1Iter", "cm1Control",
                      "pruneTol"),
                    "quadrature")
  invisible()
}

# Fail early (with an actionable message) if `init` is not a complete start.
.normalValidateInit <- function(init, omega) {
  if (is.null(names(init)) || any(!nzchar(names(init))))
    stop(".fitNormal: `init` must be a fully named numeric vector.", call. = FALSE)
  miss <- setdiff(omega$cholPars, names(init))
  if (length(miss))
    stop(".fitNormal: `init` is missing omega Cholesky parameter(s): ",
         paste(miss, collapse = ", "),
         ". Use emInit(structural, omega) to assemble a complete start.",
         call. = FALSE)
  bad <- names(init)[!is.finite(init)]
  if (length(bad))
    stop(".fitNormal: `init` must be finite; got NA/NaN/Inf for: ",
         paste(bad, collapse = ", "), call. = FALSE)
  invisible()
}


#' Assemble a complete NLME starting parameter vector
#'
#' Merges structural (and residual-error) starting values with default
#' starting values for every omega Cholesky parameter, producing a fully named
#' vector suitable as the `init` argument of [EM]. Diagonal (log-Cholesky)
#' entries default to `log(sd)`; off-diagonal correlations default to 0.
#'
#' @param structural Named numeric of structural and error-model starting
#'   values (everything except the omega Cholesky parameters).
#' @param omega An [omega] spec.
#' @param sd Numeric starting standard deviation for each random effect
#'   (diagonal Cholesky entries are set to `log(sd)`). Default `0.3`.
#' @return A named numeric vector containing `structural` followed by every
#'   `omega$cholPars`.
#' @seealso [EM], [omega()]
#' @noRd
.initNormal <- function(structural, omega, sd = 0.3) {
  if (!inherits(omega, "omegaspec"))
    stop("`omega` must be built by omega().", call. = FALSE)
  if (is.null(names(structural)) || any(!nzchar(names(structural))))
    stop("`structural` must be a fully named numeric vector.", call. = FALSE)
  chol <- setNames(rep(0, length(omega$cholPars)), omega$cholPars)
  chol[omega$isDiag] <- log(sd)
  dup <- intersect(names(structural), names(chol))
  if (length(dup))
    warning(".initNormal: `structural` already contains omega Cholesky ",
            "parameter(s): ", paste(dup, collapse = ", "),
            "; the supplied values are kept.", call. = FALSE)
  c(structural, chol[setdiff(names(chol), names(structural))])
}


#' Fit a nonlinear mixed-effects model
#'
#' Runs the selected marginal-likelihood estimator on a composed NLME
#' objective. The prediction function, data, error model, and random-effects
#' spec are recovered from `obj` (stamped at construction by [normL2] and
#' [constraintL2]), so they are never passed a second time. Returns an
#' `EM` S3 object consumable by [summary.em], [predict.em],
#' [plot.em], [plotIndivs] etc.
#'
#' @param obj An \code{objfn} of the form
#'   `normL2(data, g*x*p, errmodel = err) + constraintL2(mu = 0, Omega = om)`.
#'   Its model pieces are extracted automatically.
#' @param init Named numeric starting parameter vector. Must contain all
#'   structural parameters and all `omega$cholPars`; build it with [emInit].
#' @param fixed Optional named-numeric of fixed parameters.
#' @param method Estimator. \code{"focei"} runs FOCEI (Laplace + trust with
#'   the analytical \eqn{\partial \log |H_i| / \partial \theta} correction);
#'   \code{"quadrature"} runs adaptive sparse-grid Gauss-Hermite + ECM with a
#'   cold start; \code{"foceiQuadrature"} runs FOCEI first and uses the
#'   converged structural pars and modes as a warmstart for the quadrature
#'   polish; \code{"laplaceEM"} runs a deterministic Laplace-EM (ECM at the
#'   single-node Laplace level, with a closed-form covariance-corrected
#'   \eqn{\Omega} update via [updateOmegaChol]). `laplaceEM` reads its knobs
#'   from `control$quadrature` (`epsEcm`, `epsOfvRel`, `maxEcm`, `maxCm1Iter`,
#'   `cm1Control`). \code{"saem"} runs a stochastic-approximation EM (MCMC
#'   E-step of the random effects + Robbins-Monro averaging + closed-form
#'   \eqn{\Omega}); it reads `control$saem` (`nBurnin`, `nEM`, `nMcmc`,
#'   `cm1Control`, `cores`).
#' @param control Nested list of method-specific knobs. Entries:
#'   \describe{
#'     \item{`$focei`}{Recognised keys: `innerControl`, `trustControl`,
#'       `cores` (subject-level fork parallelism for the Stage-2 correction;
#'       default 1 = serial. Composes under [msEM]'s fit-level `cores` --
#'       keep the product below your core count).}
#'     \item{`$quadrature`}{Recognised keys: `level`, `cores`, `epsQuadLevels`,
#'       `epsEcm`, `epsOfvRel`, `maxEcmPerStage`, `maxCm1Iter`, `cm1Control`,
#'       `pruneTol` (node weight-pruning threshold, default `Inf` = off). Also
#'       consumed by `method = "laplaceEM"` (`maxEcm`).}
#'     \item{`$saem`}{Recognised keys: `nBurnin`, `nEM`, `nMcmc`, `cm1Control`,
#'       `cores`.}
#'   }
#'   Unrecognised keys raise a warning. `cm1Control` accepts any [trust]
#'   argument; `innerControl`/`trustControl` steer the C++ FOCEI loops and take
#'   `rinit`, `rmax`, `iterlim`, `ftol`, `mtol` (plus `eigen_floor_relative`
#'   for the inner one).
#' @param verbose Logical. If TRUE prints solver progress.
#'
#' @return An `EM` S3 list with fields `argument`, `value` (plain-ML
#'   `-2 log L`, the same additive-constant convention for both estimators;
#'   `ofvType = "-2LL"`), `value_nonmem` (`= value - nObs * log(2*pi)`, the
#'   NONMEM-`OBJ`-comparable value), `nObs`, `gradient`, `hessian`, `omega`,
#'   `etaModes`, `converged`, `iterations`, `emDiag`, `method`, `foceiStart`,
#'   `stageTrace`, `prdfn`, `data`, `omega`, `errfn`.
#'
#' @seealso [summary.em], [emObjfn], [omega], [emInit], [predict.em]
#' @noRd
.fitNormal <- function(obj, init,
                    fixed   = NULL,
                    method  = c("focei", "quadrature", "foceiQuadrature",
                                "laplaceEM", "saem"),
                    control = list(),
                    verbose = TRUE) {
  method <- match.arg(method)
  .normalValidateControl(control)
  rec   <- .normalReconstruct(obj)
  omega <- rec$omega
  .normalValidateInit(init, omega)

  fc <- control$focei      %||% list()
  qc <- control$quadrature %||% list()

  if (method == "saem") {
    sc <- control$saem %||% list()
    return(.normalSaem(obj, omega, init,
                    prdfn = rec$prdfn, data = rec$data, errfn = rec$errfn,
                    fixed   = fixed,
                    nBurnin = sc$nBurnin %||% 200L,
                    nEM     = sc$nEM     %||% 200L,
                    nMcmc   = sc$nMcmc   %||% 10L,
                    cm1Control = sc$cm1Control %||% list(),
                    cores   = sc$cores   %||% 1L,
                    verbose = verbose))
  }

  if (method == "laplaceEM") {
    em <- emObjfn(obj, control = qc)
    return(.normalLaplaceEM(em, init, fixed = fixed,
                         epsEcm     = qc$epsEcm     %||% 1e-4,
                         epsOfvRel  = qc$epsOfvRel  %||% 1e-5,
                         maxEcm     = qc$maxEcm     %||% 200L,
                         maxCm1Iter = qc$maxCm1Iter %||% 30L,
                         cm1Control = qc$cm1Control %||% list(),
                         verbose    = verbose))
  }

  if (method == "focei") {
    return(.normalFocei(obj, omega, init,
                        prdfn = rec$prdfn, data = rec$data, errfn = rec$errfn,
                        fixed = fixed,
                        innerControl = fc$innerControl %||% list(),
                        trustControl = fc$trustControl %||% list(),
                        cores        = fc$cores %||% 1L,
                        secondOrderCorrection = isTRUE(fc$secondOrderCorrection),
                        methodLabel  = "focei"))
  }

  if (method == "quadrature") {
    em <- emObjfn(obj, control = qc)
    return(.normalQuadrature(em, init, fixed = fixed,
                             foceiStart     = NULL,
                             epsQuadLevels  = qc$epsQuadLevels,
                             epsEcm         = qc$epsEcm         %||% 1e-4,
                             epsOfvRel      = qc$epsOfvRel      %||% 1e-5,
                             maxEcmPerStage = qc$maxEcmPerStage %||% 5L,
                             maxCm1Iter     = qc$maxCm1Iter     %||% 30L,
                             cm1Control     = qc$cm1Control     %||% list(),
                             methodLabel    = "quadrature",
                             verbose        = verbose))
  }

  # foceiQuadrature: FOCEI warmstart + quadrature polish.
  if (verbose) message(".fitNormal: running FOCEI warmstart ...")
  foceiStart <- .normalFocei(obj, omega, init,
                             prdfn = rec$prdfn, data = rec$data, errfn = rec$errfn,
                             fixed = fixed,
                             innerControl = fc$innerControl %||% list(),
                             trustControl = fc$trustControl %||% list(),
                             cores        = fc$cores %||% 1L,
                             secondOrderCorrection = isTRUE(fc$secondOrderCorrection),
                             methodLabel  = "focei")
  if (verbose) message(sprintf("  warmstart OFV = %.6f", foceiStart$value))
  em_qd <- emObjfn(obj, control = qc)
  .normalQuadrature(em_qd, foceiStart$argument, fixed = fixed,
                    foceiStart     = foceiStart,
                    epsQuadLevels  = qc$epsQuadLevels,
                    epsEcm         = qc$epsEcm         %||% 1e-4,
                    epsOfvRel      = qc$epsOfvRel      %||% 1e-5,
                    maxEcmPerStage = qc$maxEcmPerStage %||% 5L,
                    maxCm1Iter     = qc$maxCm1Iter     %||% 30L,
                    cm1Control     = qc$cm1Control     %||% list(),
                    methodLabel    = "foceiQuadrature",
                    verbose        = verbose)
}


#' Print a mixed-effects fit
#'
#' @param x An `EM` object (see [EM]).
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.em <- function(x, ...) {
  if (identical(x$prior, "penaltyL1")) .print_EM_penalty(x, ...)
  else .print_EM_omega(x, ...)
  invisible(x)
}

## Methods below are Gaussian-random-effect (omega) specific; a penaltyL1 fit
## does not carry an Omega / Fisher information, so guard with a clear message.
.emRequireOmega <- function(x, fn) {
  if (!identical(x$prior, "omega"))
    stop(fn, "() is implemented for omega() fits; a penaltyL1() fit exposes ",
         "print(), coef() and [sparsify].", call. = FALSE)
}

## Omega (Gaussian random-effect) branch of print.em.
.print_EM_omega <- function(x, ...) {
  cat("EM (prior = omega, method = ", x$method, ")\n", sep = "")
  cat(sprintf("  OFV (-2 log L): %.6f\n", x$value))
  if (!is.null(x$value_nonmem) && !is.na(x$value_nonmem))
    cat(sprintf("  OFV (NONMEM OBJ): %.6f\n", x$value_nonmem))
  cat(sprintf("  converged    : %s   iterations: %s\n",
              x$converged, format(x$iterations %||% NA_integer_)))
  cat("  argument     :\n")
  print(x$argument)
  if (!is.null(x$Omega)) {
    cat("  Omega        :\n")
    print(round(x$Omega, 4))
  }
  if (!is.null(x$etaModes) && !is.null(x$etaSE)) {
    cat("  eta (mode +/- SE, shrinkage):\n")
    print(.normalEtaTable(x$etaModes, x$etaSE, x$shrinkage))
  }
  invisible(x)
}

# Build a numeric data.frame with one column per quantity (mode, SE,
# optional shrinkage) per eta. Relies on print.data.frame for alignment.
.normalEtaTable <- function(etaModes, etaSE, shrinkage) {
  K <- ncol(etaModes)
  eta_names <- colnames(etaModes)
  cols <- c(rbind(eta_names, paste0("SE.", eta_names)))
  vals <- cbind(etaModes, etaSE)[, c(rbind(seq_len(K), seq_len(K) + K)),
                                 drop = FALSE]
  colnames(vals) <- cols

  if (!is.null(shrinkage)) {
    shr <- shrinkage
    colnames(shr) <- paste0("shrink.", eta_names)
    vals <- cbind(vals, shr)
  }
  round(as.data.frame(vals), 3)
}


#' Structural parameter estimates of an EM
#'
#' @param object An `EM` object (see [EM]).
#' @param ... Ignored.
#' @return The named numeric vector of fitted structural (and error-model)
#'   parameters.
#' @export
coef.em <- function(object, ...) object$argument


#' Wald confidence intervals for an EM
#'
#' Symmetric confidence intervals from the observed-information matrix
#' (`vcov(object)`, the inverse of half the outer Hessian). For
#' likelihood-based (profile) intervals use [profile()] / [confint.parframe].
#'
#' @param object An `EM` object (see [EM]).
#' @param parm Optional character vector of parameter names to report. Defaults
#'   to all structural parameters.
#' @param level Confidence level. Default 0.95.
#' @param ... Ignored.
#' @return A data.frame with columns `name`, `value`, `lower`, `upper`.
#' @export
confint.em <- function(object, parm = NULL, level = 0.95, ...) {
  .emRequireOmega(object, "confint")
  est <- object$argument
  se  <- .normalStructuralSE(object)
  z   <- stats::qnorm(1 - (1 - level) / 2)
  nm  <- if (is.null(parm)) names(est) else intersect(parm, names(est))
  data.frame(name  = nm,
             value = unname(est[nm]),
             lower = unname(est[nm] - z * se[nm]),
             upper = unname(est[nm] + z * se[nm]),
             row.names = NULL)
}

# Structural-parameter standard errors from the outer observed information.
# Reuses vcov() (statistics.R), which returns solve(0.5 * hessian) under dMod's
# -2 log L convention. Robust to partial Hessians (name-matched) and to
# non-positive diagonals (returned as NA).
.normalStructuralSE <- function(object) {
  est <- object$argument
  se  <- setNames(rep(NA_real_, length(est)), names(est))
  V   <- tryCatch(vcov(object), error = function(e) NULL)
  if (!is.null(V) && !is.null(rownames(V))) {
    cn <- intersect(rownames(V), names(est))
    d  <- diag(V)[cn]
    d[d < 0] <- NA_real_
    se[cn] <- sqrt(d)
  }
  se
}


#' Summarise an EM
#'
#' Reports the population-parameter estimates with standard errors and relative
#' standard errors (RSE%), the random-effect covariance as standard deviations
#' plus a correlation matrix, mean eta shrinkage, and the objective value.
#'
#' @param object An `EM` object (see [EM]).
#' @param ... Ignored.
#' @return An object of class `summary.em` (a list with `method`, `value`,
#'   `converged`, `iterations`, `population` data.frame, `omegaSD`, `omegaCor`,
#'   and `shrinkage`), with a `print` method.
#' @seealso [EM], [coef.em], [confint.em]
#' @export
summary.em <- function(object, ...) {
  .emRequireOmega(object, "summary")
  est <- object$argument
  se  <- .normalStructuralSE(object)
  pop <- data.frame(estimate = unname(est),
                    se       = unname(se[names(est)]),
                    rse.pct  = unname(100 * se[names(est)] / abs(est)),
                    row.names = names(est))
  Om        <- object$Omega
  omega_sd  <- if (!is.null(Om)) sqrt(diag(Om)) else NULL
  omega_cor <- if (!is.null(Om) && nrow(Om) > 1L) stats::cov2cor(Om) else NULL
  shr       <- object$shrinkage
  out <- list(method       = object$method,
              value        = object$value,
              value_nonmem = object$value_nonmem,
              ofvType      = object$ofvType %||% "-2LL",
              converged    = object$converged,
              iterations   = object$iterations,
              population   = pop,
              omegaSD      = omega_sd,
              omegaCor     = omega_cor,
              shrinkage    = if (!is.null(shr)) colMeans(shr, na.rm = TRUE) else NULL)
  class(out) <- c("summary.em", "list")
  out
}

#' Print a summary.em
#' @param x A `summary.em` object.
#' @param digits Number of significant digits. Default 4.
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.summary.em <- function(x, digits = 4, ...) {
  cat(".fitNormal summary (method = ", x$method, ")\n", sep = "")
  cat(sprintf("  OFV (%s): %.6f   converged: %s   iterations: %s\n",
              x$ofvType %||% "-2LL", x$value, x$converged,
              format(x$iterations %||% NA_integer_)))
  if (!is.null(x$value_nonmem) && !is.na(x$value_nonmem))
    cat(sprintf("  OFV (NONMEM OBJ, = value - N_obs*log(2*pi)): %.6f\n",
                x$value_nonmem))
  cat("\n  Population parameters:\n")
  print(signif(x$population, digits))
  if (!is.null(x$omegaSD)) {
    cat("\n  Random effects (Omega):\n")
    cat("    SD:  ", paste(sprintf("%s=%.4g", names(x$omegaSD), x$omegaSD),
                           collapse = "  "), "\n", sep = "")
    if (!is.null(x$omegaCor)) {
      cat("    correlation:\n")
      print(round(x$omegaCor, digits))
    }
  }
  if (!is.null(x$shrinkage)) {
    cat("\n  Mean eta shrinkage:\n    ",
        paste(sprintf("%s=%.3f", names(x$shrinkage), x$shrinkage),
              collapse = "  "), "\n", sep = "")
  }
  invisible(x)
}



#' Per-subject adaptive-quadrature marginal-likelihood evaluator
#'
#' @description
#' Evaluates the per-subject marginal likelihood
#' \eqn{\hat L_i = \int p(y_i \mid \eta) \, N(\eta \mid 0, \Omega)\, d\eta}
#' via sparse-grid Gauss-Hermite quadrature at a precomputed set of nodes
#' (output of [makeSubjectNodes]). Returns the log-likelihood, posterior
#' first and second moments of \eqn{\eta_i}, plus (optionally) the gradient
#' and Hessian of \eqn{-2 \log \hat L_i} with respect to a chosen subset of
#' outer (population) parameters.
#'
#' This helper bypasses [normL2] / [constraintL2] entirely: the data
#' likelihood contribution comes from the lifted [evalConditionResidual]
#' (one prediction call per node, single-condition), the MVN prior on the
#' subject's \eqn{\eta_b} is added in closed form via `omega$buildL`.
#' Avoids the multi-condition parameter rebinding and full-population MVN
#' contribution that calling `obj(..., conditions = subjects[i])` per node
#' would incur.
#'
#' @param subjIdx Integer in `1..N` selecting the subject.
#' @param psiFull Named numeric, the full outer parameter vector at which to
#'   evaluate (structural + chol_pars; chol_pars are frozen during CM-1).
#' @param etaModes N x K matrix of all subjects' eta values; the row at
#'   `subjIdx` is ignored. Other rows are passed through to the prediction
#'   call (required for parameter completeness in joint models).
#' @param omega An [omega] spec object with subject expansion.
#' @param nodesSubj Output of `makeSubjectNodes(eta_hat_i, H_i, level)` for
#'   the active subject.
#' @param xPred The prediction function (e.g. `g * x * p`).
#' @param datalist The [datalist] used for the joint objective.
#' @param errfn Optional obsfn (passed through to
#'   [evalConditionResidual]).
#' @param fixed Optional fixed-parameter vector.
#' @param outerActiveNames Character vector of parameter names to track
#'   gradient/Hessian for. Defaults to `names(psiFull)`. For CM-1 with frozen
#'   chol_pars, pass just the structural names.
#' @param mode `"moments_only"` (no gradient/Hessian, cheap) or
#'   `"with_grad"` (gradient + Hessian of `-2 log L_i`).
#'
#' @return A list with components:
#' \describe{
#'   \item{`logLhat`}{`log L_i` (scalar).}
#'   \item{`m_hat`}{Length-K named numeric, posterior mean of \eqn{\eta_i}.}
#'   \item{`M_hat`}{K x K named matrix, posterior 2nd moment
#'     \eqn{\hat E[\eta_i \eta_i^T \mid y_i]}.}
#'   \item{`maxSoftmax`}{Max |softmax weight| across nodes; diagnostic of
#'     grid concentration.}
#'   \item{`n_eff`}{Effective node count, 1/sum(softmax^2); diagnostic.}
#'   \item{`gradient`, `hessian`}{When `mode = "with_grad"`: gradient and
#'     Hessian of `-2 log L_i` w.r.t. `outerActiveNames`. NULL otherwise.}
#' }
#'
#' @seealso [makeSubjectNodes], [evalConditionResidual]
#' @noRd
.normalEcmSubject <- function(subjIdx, psiFull, etaModes,
                                 omega, nodesSubj,
                                 xPred, datalist,
                                 errfn           = NULL,
                                 fixed              = NULL,
                                 outerActiveNames = NULL,
                                 mode               = c("moments_only", "with_grad"),
                                 cores              = 1L) {
  mode <- match.arg(mode)
  with_grad <- (mode == "with_grad")

  K            <- omega$K
  subject_etas <- omega$subjectEtas
  subjects     <- rownames(subject_etas)
  cn           <- subjects[subjIdx]
  eta_i_names  <- subject_etas[subjIdx, ]
  chol_pars    <- omega$cholPars

  if (is.null(outerActiveNames)) outerActiveNames <- names(psiFull)
  if (!all(chol_pars %in% names(psiFull)))
    stop(".normalEcmSubject: `psiFull` must contain all omega$cholPars.")

  # Build closed-form prior log-normalisation: log N(eta_b|0,Omega)
  # = -K/2 log(2 pi) - log|L_omega| - 0.5 * z^T z, z = L_omega^{-1} eta_b.
  L_omega         <- omega$buildL(psiFull[chol_pars])
  log_det_L_omega <- sum(log(diag(L_omega)))
  log_norm_prior  <- -K / 2 * log(2 * pi) - log_det_L_omega

  # Other subjects' eta values, named with their per-subject eta names.
  if (nrow(etaModes) != length(subjects))
    stop(".normalEcmSubject: `etaModes` must have one row per subject.")
  other_eta_nm   <- as.vector(subject_etas[-subjIdx, , drop = FALSE])
  other_eta_vals <- as.vector(etaModes[-subjIdx, , drop = FALSE])
  names(other_eta_vals) <- other_eta_nm

  dataI  <- datalist[[cn]]
  times_i <- dataI$time

  B           <- nrow(nodesSubj$etaNodes)
  Q           <- length(outerActiveNames)
  log_int     <- numeric(B)
  sign_int    <- nodesSubj$weightSigns
  per_node_gr <- if (with_grad) matrix(0, B, Q,
                                       dimnames = list(NULL, outerActiveNames))
                 else NULL
  per_node_he <- if (with_grad) array(0, c(B, Q, Q),
                                       dimnames = list(NULL, outerActiveNames,
                                                       outerActiveNames))
                 else NULL

  # Constant parameter template; only the subject-i eta block varies per node,
  # so build the full vector once and overwrite that block each iteration.
  full_pars <- c(psiFull, setNames(numeric(K), eta_i_names), other_eta_vals)

  # Closed-form log-priors for ALL nodes at once: only eta_b varies and L_omega
  # is constant, so one matrix forwardsolve replaces B per-node solves.
  Z_prior       <- forwardsolve(L_omega, t(nodesSubj$etaNodes))   # K x B
  log_prior_all <- log_norm_prior - 0.5 * colSums(Z_prior^2)

  # Each node is its own nonlinear IVP, so they cannot be fused into one solve
  # -- but they are independent, so they go out as one batch against the same
  # condition. Bit-identical to the per-node form.
  pars_nodes <- lapply(seq_len(B), function(b) {
    fp <- full_pars
    fp[eta_i_names] <- nodesSubj$etaNodes[b, ]
    fp
  })
  preds_nodes <- .predictMany(xPred, times = times_i, parsList = pars_nodes,
                              conditions = rep(cn, B), fixed = fixed,
                              deriv = with_grad, cores = cores)

  for (b in seq_len(B)) {
    eta_b <- nodesSubj$etaNodes[b, ]
    full_pars[eta_i_names] <- eta_b

    pred_b <- structure(list(preds_nodes[[b]]), names = cn)
    res_b  <- evalConditionResidual(
      dataI        = dataI,
      predictionI  = pred_b[[cn]],
      pars          = full_pars,
      errfn      = errfn,
      fixed         = fixed,
      cn            = cn,
      eCondNames  = NULL,
      deriv         = with_grad,
      deriv2        = FALSE)

    # log integrand = log|W_b| + log p(y|eta_b) + log p(eta_b|Omega); the prior
    # term is precomputed for all nodes above (log_prior_all).
    # log p(y|eta_b) = -0.5 * res_b$value (the value carries -2 log p form).
    log_int[b] <- nodesSubj$logAbsWeights[b] - 0.5 * res_b$value + log_prior_all[b]

    if (with_grad) {
      gr <- res_b$gradient[outerActiveNames]
      gr[is.na(gr)] <- 0
      per_node_gr[b, ] <- -0.5 * gr
      he <- res_b$hessian[outerActiveNames, outerActiveNames, drop = FALSE]
      he[is.na(he)] <- 0
      per_node_he[b, , ] <- -0.5 * he
    }
  }

  # Signed log-sum-exp: log(sum_b sign_b * exp(log_int[b])).
  M       <- max(log_int)
  shifted <- sign_int * exp(log_int - M)
  s       <- sum(shifted)
  if (!is.finite(s) || s <= 0) {
    warning(sprintf(paste0(".normalEcmSubject(subjIdx=%d): signed-LSE ",
                           "sum = %.3e; grid is under-resolved (raise `level`)."),
                    subjIdx, s))
    logLhat <- -Inf
    softmax <- numeric(B)
  } else {
    logLhat <- M + log(s)
    softmax <- sign_int * exp(log_int - logLhat)
  }

  m_hat <- as.numeric(softmax %*% nodesSubj$etaNodes)
  M_hat <- crossprod(nodesSubj$etaNodes,
                     softmax * nodesSubj$etaNodes)
  names(m_hat) <- omega$eta
  dimnames(M_hat) <- list(omega$eta, omega$eta)

  out <- list(logLhat     = logLhat,
              m_hat       = m_hat,
              M_hat       = M_hat,
              maxSoftmax = max(abs(softmax)),
              n_eff       = if (any(softmax != 0)) 1 / sum(softmax^2) else 0)

  if (with_grad) {
    # Gradient of log L_i: sum_b softmax_b * (d log integrand_b / d theta).
    gr_lse <- as.numeric(softmax %*% per_node_gr)
    names(gr_lse) <- outerActiveNames

    # Hessian of LSE = sum softmax_b H_b + Cov_softmax(grad_b).
    # Cov = E[gg^T] - E[g] E[g]^T (with signed softmax this is an algebraic identity).
    H_avg <- matrix(0, Q, Q, dimnames = list(outerActiveNames, outerActiveNames))
    for (b in seq_len(B)) H_avg <- H_avg + softmax[b] * per_node_he[b, , ]
    Eggt   <- crossprod(per_node_gr, softmax * per_node_gr)
    cov_g  <- Eggt - tcrossprod(gr_lse)
    H_lse  <- H_avg + cov_g

    out$gradient <- -2 * gr_lse
    out$hessian  <- -2 * H_lse
  }
  out
}


# Drop heavy state (emDiag, prdfn, data, omega, errfn, foceiStart,
# stageTrace) from an .fitNormal so a parlist of .msfitNormal results stays small.
# Keeps everything as.parframe.parlist + summary.parlist + downstream
# diagnostics consume: argument, value, gradient, hessian, omega, etaModes,
# etaSE, shrinkage, converged, iterations, method.
.stripNormalFit <- function(fit) {
  keep <- c("argument", "value", "gradient", "hessian",
            "omega", "etaModes", "etaSE", "shrinkage",
            "converged", "iterations", "method", "prior")
  out <- fit[intersect(keep, names(fit))]
  class(out) <- class(fit)
  out
}


#' Multi-start nonlinear mixed-effects fit
#'
#' Runs [EM] from many starting points in parallel, returning a
#' [parlist] of fits sorted-ready for [as.parframe()] / [summary()]. `center`
#' is either a named-numeric (perturbed by `samplefun` per fit) or a
#' [parframe] (each row used as a starting point). Use this to characterise
#' the multi-modality of the marginal likelihood and pick the best optimum.
#'
#' @param obj An `objfn` passed straight to [EM] (typically
#'   `normL2(data, g*x*p, errmodel = e) + constraintL2(mu = 0, Omega = om)`).
#'   The model pieces are recovered from it automatically.
#' @param center Named numeric or [parframe]. If numeric, the population
#'   parameter vector around which random starts are sampled (structural pars
#'   plus the omega Cholesky parameters). If a parframe, each row is used as a
#'   fixed starting point and `fits` is overridden by `nrow(center)`.
#' @param fixed Optional named-numeric of fixed parameters.
#' @param method Estimator. Passed through to [EM]: `"focei"`,
#'   `"quadrature"`, or `"foceiQuadrature"`.
#' @param control Nested control list, passed through to [EM].
#' @param fits Integer, number of random starts. Ignored when `center` is a
#'   parframe.
#' @param cores Integer, number of parallel workers. On Unix uses
#'   `parallel::mclapply` (fork); on Windows a PSOCK cluster + `foreach`.
#'   Outer parallelism multiplies with the OpenMP threads each objective
#'   function uses (its `cores` argument); keep the product below your core
#'   count.
#' @param samplefun Name of a random-number generator (default `"rnorm"`)
#'   used to perturb `center`. Extra args in `...` whose names match
#'   `formals(samplefun)` are forwarded.
#' @param start1stfromCenter Logical. If `TRUE`, the first fit starts at
#'   `center` itself (no perturbation). Ignored when `center` is a parframe.
#' @param keepFull Logical. If `FALSE` (the default), each returned fit is
#'   stripped of heavy state (`emDiag`, `prdfn`, `data`, `omega`,
#'   `errfn`, `foceiStart`, `stageTrace`) so the result stays small.
#'   Set `TRUE` if you need to call [predict.em()] / [plot.em()]
#'   etc. on individual fits.
#' @param studyname Optional character. If `output = TRUE`, fits are written
#'   to `<resultPath>/<studyname>/trial-N-<timestamp>/interRes/`. Defaults to
#'   `"msem"`.
#' @param resultPath Character, base directory for the on-disk dump.
#' @param output Logical. If `TRUE`, each fit is saved as it completes
#'   (crash-resilient) and the full parlist is written at the end.
#' @param verbose Logical. If `TRUE`, prints per-fit progress and forwards
#'   `verbose = TRUE` into [EM].
#' @param retry Logical. If `TRUE` (default), a fit that throws is
#'   re-attempted with a freshly sampled `parinit` up to `nTries` times.
#'   Ignored when `center` is a parframe.
#' @param nTries Maximum number of attempts per fit slot, including the
#'   first. Default `10L`.
#' @param ... Forwarded to `samplefun` (e.g. `sd = 0.3`).
#'
#' @return A [parlist] of length `fits`, each element an `EM` (stripped
#'   per `keepFull`) with `parinit` and `index` attached. Pass to
#'   [as.parframe()] for a sorted-by-`value` table. Failed fits are stored
#'   as `list(error = ..., value = NA, converged = FALSE, ...)`.
#'
#' @seealso [EM], [mstrust()], [msParframe()], [parlist].
#' @noRd
.msfitNormal <- function(obj, center,
                      fixed    = NULL,
                      method   = c("focei", "quadrature", "foceiQuadrature"),
                      control  = list(),
                      fits     = 20,
                      cores    = 1,
                      samplefun = "rnorm",
                      start1stfromCenter = FALSE,
                      keepFull   = FALSE,
                      studyname  = NULL,
                      resultPath = ".",
                      output     = FALSE,
                      verbose    = FALSE,
                      retry      = TRUE,
                      nTries     = 10L,
                      ...) {

  method <- match.arg(method)
  cores  <- .sanitizeCores(cores)

  # Build the per-fit starting points. parframe input fixes the number of
  # fits and pulls starts directly from rows (sorted-by-value via
  # as.parvec.parframe). Numeric input perturbs by samplefun().
  varargslist <- list(...)
  argssample  <- NULL  # populated for non-parframe centers; reused on retry
  if (is.parframe(center)) {
    fits <- nrow(center)
    parInitList <- lapply(seq_len(fits), function(i) as.parvec(center, i))
  } else {
    if (is.null(names(center)) || any(!nzchar(names(center))))
      stop("`center` must be a fully named numeric vector or a parframe.")
    namessample <- intersect(names(formals(samplefun)), names(varargslist))
    argssample  <- varargslist[namessample]
    argssample$n <- length(center)
    parInitList <- lapply(seq_len(fits), function(i) {
      if (i == 1L && start1stfromCenter) {
        center
      } else {
        perturb <- do.call(samplefun, argssample)
        out <- center + perturb
        names(out) <- names(center)
        out
      }
    })
  }
  resample_init <- function() {
    perturb <- do.call(samplefun, argssample)
    out <- center + perturb
    names(out) <- names(center)
    out
  }
  cores <- min(fits, cores)

  # Nested-parallelism guard. .msfitNormal forks `cores` FIT workers; each inner
  # .fitNormal may itself fork `subject_cores` SUBJECT workers via
  # control$<method>$cores (default 1 = fits-only). The two levels compose
  # (e.g. 4 fits x 4 subjects = 16 processes) but their product must stay below
  # the physical core count or the machine oversubscribes. Warn, do not clamp:
  # the user opted into the nesting explicitly.
  inner_cores <- max(1L,
                     control$focei$cores      %||% 1L,
                     control$quadrature$cores %||% 1L,
                     control$saem$cores       %||% 1L)
  n_detected <- parallel::detectCores()
  if (is.finite(n_detected) && cores * inner_cores > n_detected)
    warning(sprintf(paste0(".msfitNormal: %d fit-workers x %d subject-cores = %d ",
                           "processes exceeds %d detected cores; expect ",
                           "oversubscription. Lower `cores` or ",
                           "control$<method>$cores."),
                    cores, inner_cores, cores * inner_cores, n_detected),
            call. = FALSE)

  # Optional on-disk dump (crash-resilient): one .Rda per fit + a final
  # parameterList.Rda.
  interResultFolder <- NULL
  resultFolder      <- NULL
  if (output) {
    if (is.null(studyname)) studyname <- ".msfitNormal"
    m_timeStamp <- format(Sys.time(), "%d-%m-%Y-%H%M%S")
    resultFolderBase <- file.path(resultPath, studyname)
    n_existing <- length(dir(resultFolderBase, pattern = "trial*"))
    m_trial <- paste0("trial-", n_existing + 1L)
    resultFolder <- file.path(resultFolderBase,
                              paste0(m_trial, "-", m_timeStamp))
    interResultFolder <- file.path(resultFolder, "interRes")
    dir.create(interResultFolder, showWarnings = FALSE, recursive = TRUE)
  }

  digits <- if (fits >= 10L) floor(log10(fits)) + 1L else 1L

  doOne <- function(i) {
    init_i <- parInitList[[i]]
    # Retry: on a try-error, draw a fresh parinit and re-run .fitNormal, up to
    # nTries times. parframe-supplied centers skip retry (rows are taken as
    # given). Pequil/Pimpl warm-start caches are flushed between attempts so
    # the retry does not inherit a dead basin from the previous start.
    max_tries <- if (isTRUE(retry) && !is.parframe(center)) as.integer(nTries) else 1L
    t0 <- Sys.time()
    fit <- NULL
    for (try_i in seq_len(max_tries)) {
      fit <- try(suppressMessages(
        .fitNormal(obj, init_i,
                fixed = fixed, method = method, control = control,
                verbose = verbose)),
        silent = !verbose)
      if (!inherits(fit, "try-error") || try_i == max_tries) break
      init_i <- resample_init()
      try(resetWarmStarts(obj, verbose = FALSE), silent = TRUE)
    }
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    if (inherits(fit, "try-error")) {
      fit <- list(error      = as.character(fit),
                  value      = NA_real_,
                  converged  = FALSE,
                  iterations = NA_integer_,
                  prior      = "omega",
                  method     = method)
      class(fit) <- c("em", "list")
    } else if (!keepFull) {
      fit <- .stripNormalFit(fit)
    }
    fit$parinit <- init_i
    fit$index   <- i
    fit$elapsed <- elapsed

    if (verbose) {
      msg <- if (!is.null(fit$error)) {
        sprintf("[.msfitNormal %s/%d] FAILED (%.1fs): %s",
                formatC(i, width = digits, flag = "0"), fits,
                elapsed, fit$error)
      } else {
        sprintf("[.msfitNormal %s/%d] OFV=%.6f  conv=%s  iter=%s  (%.1fs)",
                formatC(i, width = digits, flag = "0"), fits,
                fit$value, fit$converged,
                format(fit$iterations %||% NA_integer_), elapsed)
      }
      message(msg)
    }

    if (output) {
      saveRDS(fit, file = file.path(interResultFolder,
                                    sprintf("fit-%d.Rda", i)))
    }
    fit
  }

  # Parallel dispatch. Fork on Unix; PSOCK + foreach on Windows, and also when
  # a condition axis is asked for -- a forked worker cannot carry one. Falls
  # back to serial when cores == 1 or fits == 1.
  .cc <- .splitCores(cores, "fits")
  coresConditions <- .cc$conditions; cores <- .cc$outer
  if (!is.null(coresConditions) && cores == 1L)
    options(dMod.cores = coresConditions)
  if (cores > 1L) {
    if (Sys.info()[['sysname']] == "Windows" ||
        (!is.null(coresConditions) && coresConditions > 1L)) {
      cluster <- parallel::makeCluster(cores)
      on.exit(parallel::stopCluster(cluster), add = TRUE)
      doParallel::registerDoParallel(cluster)
      parallel::clusterCall(cl = cluster,
                            function(x) .libPaths(x), .libPaths())
      parallel::clusterExport(
        cluster, envir = environment(),
        varlist = c("obj", "omega", "parInitList", "prdfn", "data",
                    "errfn", "fixed", "method", "control", "verbose",
                    "keepFull", "output", "interResultFolder", "fits",
                    "digits", "retry", "nTries", "center", "samplefun",
                    "argssample", "resample_init", "coresConditions"))
      `%mydo%` <- foreach::`%dopar%`
      i <- NULL
      results <- foreach::foreach(
        i = seq_len(fits),
        .packages = .packages(),
        .inorder = TRUE,
        .options.multicore = list(preschedule = FALSE)) %mydo% {
          if (!is.null(coresConditions)) options(dMod.cores = coresConditions)
          doOne(i)
        }
    } else {
      results <- parallel::mclapply(seq_len(fits), doOne,
                                    mc.cores = cores,
                                    mc.preschedule = FALSE)
    }
  } else {
    results <- lapply(seq_len(fits), doOne)
  }

  if (output) {
    saveRDS(results, file = file.path(resultFolder, "parameterList.Rda"))
  }

  as.parlist(results)
}

