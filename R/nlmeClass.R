## Cholesky parameter naming convention (load-bearing):
##   diagonal:    paste0(prefix, "_", short, "_", short),  L_kk = exp(par)
##   off-diag:    paste0(prefix, "_", short_k, "_", short_l) for k > l,  L_kl = par
## where short = sub("^eta_", "", eta), so eta_Cl -> Cl, eta_V -> V, etc.



#' Specify random-effect covariance structure for NLME
#'
#' @description
#' Builds an `omegaSpec` object describing a Cholesky-parametrised random-effect
#' covariance matrix `Omega = L L^T`, with `L` lower-triangular, log-diagonal,
#' and arbitrary sparsity below the diagonal. The result is consumed by
#' [constraintL2] (when its `Omega` argument is set) and by [EM].
#'
#' @param eta Character vector of base eta names, one per random effect, e.g.
#'   `c("eta_Cl", "eta_V", "eta_Ka")`. Determines the dimension `K`.
#' @param structure One of `"diag"` (default, only diagonal Cholesky entries
#'   estimated, K free parameters) or `"full"` (full lower-triangular Cholesky,
#'   K(K+1)/2 free parameters).
#' @param correlate Optional list of length-2 character vectors selecting
#'   off-diagonal Cholesky entries to be estimated. Alternative to
#'   `structure = "full"`. Each entry is a pair of eta names whose correlation
#'   is freed.
#' @param subjects Optional character vector of subject identifiers. When given,
#'   the per-subject random-effect parameter names are constructed as
#'   `paste0(eta[k], "_", subjects[i])` and stored in `subjectEtas`. Required
#'   when the spec is fed into [EM] or `mcmc()`; optional standalone.
#' @param prefix Character string used as a name stem for Cholesky parameters.
#'   Default `"omega"`. With `eta = c("eta_Cl", "eta_V")` this yields
#'   `omega_Cl_Cl`, `omega_V_V`, `omega_V_Cl`.
#'
#' @return An `omegaSpec` (S3 list) with components:
#' \describe{
#'   \item{`K`}{Dimension.}
#'   \item{`eta`}{Base eta names.}
#'   \item{`cholPars`}{Character vector of Cholesky parameter names in
#'     column-major order.}
#'   \item{`nameMatrix`}{`K x K` character matrix; entry `[k, l]` is the
#'     parameter name of `L[k, l]` or `NA` if structurally zero.}
#'   \item{`isDiag`}{Named logical vector: TRUE for log-diagonal Cholesky
#'     parameters, FALSE for free off-diagonal entries.}
#'   \item{`cholLoc`}{`length(cholPars) x 2` integer matrix giving the
#'     `(row, col)` location of each Cholesky parameter in `L`.}
#'   \item{`subjectEtas`}{`N x K` character matrix of per-subject eta names,
#'     or `NULL` if `subjects` was not supplied.}
#'   \item{`buildL(cholVec)`}{Function returning the lower-triangular `L`
#'     given a named numeric vector containing at least `cholPars`.}
#' }
#'
#' @details
#' The Cholesky diagonal is on the log scale (enforces positive definiteness
#' of Omega); off-diagonal entries are unrestricted reals. Matches
#' `lme4`/`TMB` conventions.
#'
#' @examples
#' om1 <- omega(eta = c("eta_Cl", "eta_V"))
#' om1$cholPars
#' om1$buildL(c(omega_Cl_Cl = log(0.3), omega_V_V = log(0.2)))
#'
#' om2 <- omega(eta = c("eta_Cl", "eta_V", "eta_Ka"), structure = "full",
#'              subjects = paste0("subj", 1:3))
#' om2$subjectEtas
#' om2$cholPars
#'
#' @seealso [constraintL2], [EM]
#' @export
omega <- function(eta,
                  structure = c("diag", "full"),
                  correlate = NULL,
                  subjects  = NULL,
                  prefix    = "omega") {

  structure <- match.arg(structure)
  if (!is.character(eta) || length(eta) < 1L)
    stop("`eta` must be a non-empty character vector of random-effect names.")
  if (anyDuplicated(eta))
    stop("`eta` names must be unique.")

  K <- length(eta)
  short <- sub("^eta_", "", eta)

  name_matrix <- matrix(NA_character_, K, K, dimnames = list(eta, eta))
  for (k in seq_len(K)) {
    name_matrix[k, k] <- paste(prefix, short[k], short[k], sep = "_")
  }

  if (structure == "full") {
    if (K >= 2L) for (k in 2:K) for (l in 1:(k - 1)) {
      name_matrix[k, l] <- paste(prefix, short[k], short[l], sep = "_")
    }
  } else if (!is.null(correlate)) {
    if (!is.list(correlate))
      stop("`correlate` must be a list of length-2 character vectors.")
    for (pair in correlate) {
      if (length(pair) != 2L || !is.character(pair))
        stop("Each entry of `correlate` must be a length-2 character vector.")
      i1 <- match(pair[1], eta)
      i2 <- match(pair[2], eta)
      if (is.na(i1) || is.na(i2))
        stop("`correlate` refers to unknown eta names: ",
             paste(pair[is.na(c(i1, i2))], collapse = ", "))
      kk <- max(i1, i2); ll <- min(i1, i2)
      if (kk == ll) next
      name_matrix[kk, ll] <- paste(prefix, short[kk], short[ll], sep = "_")
    }
  }

  chol_pars <- character(0)
  chol_loc  <- matrix(integer(0), nrow = 0, ncol = 2,
                      dimnames = list(NULL, c("row", "col")))
  for (l in seq_len(K)) for (k in l:K) {
    nm <- name_matrix[k, l]
    if (is.na(nm)) next
    chol_pars <- c(chol_pars, nm)
    chol_loc  <- rbind(chol_loc, c(k, l))
  }
  rownames(chol_loc) <- chol_pars

  if (anyDuplicated(chol_pars))
    stop("Generated duplicate Cholesky parameter names. Check your `prefix`.")

  is_diag <- chol_loc[, "row"] == chol_loc[, "col"]
  names(is_diag) <- chol_pars

  subject_etas <- NULL
  if (!is.null(subjects)) {
    if (!is.character(subjects) || anyDuplicated(subjects))
      stop("`subjects` must be a non-empty character vector with unique entries.")
    subject_etas <- outer(subjects, eta, function(s, e) paste(e, s, sep = "_"))
    dimnames(subject_etas) <- list(subjects, eta)
  }

  build_L <- function(cholVec) {
    if (!all(chol_pars %in% names(cholVec)))
      stop("`cholVec` is missing Cholesky parameters: ",
           paste(setdiff(chol_pars, names(cholVec)), collapse = ", "))
    L <- matrix(0, K, K, dimnames = list(eta, eta))
    for (m in seq_along(chol_pars)) {
      k <- chol_loc[m, 1L]; l <- chol_loc[m, 2L]
      v <- cholVec[[chol_pars[m]]]
      L[k, l] <- if (is_diag[m]) exp(v) else v
    }
    L
  }

  out <- list(
    K           = K,
    eta         = eta,
    short       = short,
    structure   = structure,
    cholPars    = chol_pars,
    nameMatrix  = name_matrix,
    isDiag      = is_diag,
    cholLoc     = chol_loc,
    subjectEtas = subject_etas,
    buildL      = build_L,
    prefix      = prefix
  )
  class(out) <- c("omegaspec", "list")
  out
}



#' Print method for omega
#'
#' @param x An `omegaSpec` object.
#' @param ... Ignored.
#' @export
print.omegaspec <- function(x, ...) {
  cat("omega\n")
  cat(sprintf("  K           : %d\n", x$K))
  cat(sprintf("  eta         : %s\n", paste(x$eta, collapse = ", ")))
  cat(sprintf("  structure   : %s\n", x$structure))
  cat(sprintf("  cholPars    : %d (%s)\n",
              length(x$cholPars), paste(x$cholPars, collapse = ", ")))
  if (!is.null(x$subjectEtas)) {
    cat(sprintf("  subjects    : %d (%s)\n",
                nrow(x$subjectEtas),
                paste(rownames(x$subjectEtas), collapse = ", ")))
    cat(sprintf("  subjectEtas : %d unique random-effect parameters\n",
                length(unique(as.vector(x$subjectEtas)))))
  } else {
    cat("  subjects    : none (call omega(..., subjects = ...) before .fitNormal)\n")
  }
  invisible(x)
}



#' ECM closed-form update for the omega Cholesky parameters
#'
#' @description
#' Given a list of per-subject posterior second moments
#' `MHat_i = Ehat[eta_i %*% t(eta_i) | y_i]` (the output of the ECM E-step
#' implemented in `.normalEcmSubject`), returns the cholPars vector that
#' maximises the EM Q-function
#' \deqn{Q(\Omega) = -\frac{1}{2} N \, (\log|\Omega| + \mathrm{tr}(\Omega^{-1} S)),
#'   \quad S = \frac{1}{N}\sum_i \hat M_i,}
#' subject to the sparsity / log-Cholesky parametrisation encoded in
#' `omegaSpec`. For pure-diagonal and full structures the solution is closed
#' form; for selective `correlate` it is the minimiser of an analytic gradient
#' obtained via [trust] with finite-difference Hessian (convex small-dim
#' problem, converges in a handful of iterations).
#'
#' @param MHatList Length-N list of K x K positive-semidefinite matrices
#'   (subject posterior second moments).
#' @param omega An [omega] spec object.
#' @return Named numeric vector matching `omega$cholPars` (log-diagonal
#'   on diagonals, free real on off-diagonals).
#' @seealso [omega]
#' @export
updateOmegaChol <- function(MHatList, omega) {
  if (!inherits(omega, "omegaspec"))
    stop("`omega` must be an omegaSpec object (built by omega()).")
  K <- omega$K
  N <- length(MHatList)
  if (N < 1L) stop("updateOmegaChol: MHatList is empty.")

  S <- Reduce(`+`, MHatList) / N
  if (!is.matrix(S) || nrow(S) != K || ncol(S) != K)
    stop("updateOmegaChol: MHat entries must be K x K matrices.")
  S <- (S + t(S)) / 2

  if (rcond(S) < 1e-10) {
    warning("updateOmegaChol: sample second-moment is rank-deficient; ",
            "adding 1e-8 ridge.")
    S <- S + 1e-8 * diag(K)
  }

  chol_pars <- omega$cholPars
  chol_loc  <- omega$cholLoc
  is_diag   <- omega$isDiag

  # Pure diagonal (no off-diag cholPars): omega_kk = log(sqrt(S_kk)).
  if (all(is_diag)) {
    out <- setNames(numeric(length(chol_pars)), chol_pars)
    for (m in seq_along(chol_pars)) {
      k <- chol_loc[m, 1L]
      out[m] <- 0.5 * log(max(S[k, k], .Machine$double.eps))
    }
    return(out)
  }

  # Full lower-triangular Cholesky.
  if (omega$structure == "full") {
    L_full <- t(chol(S))
    out <- setNames(numeric(length(chol_pars)), chol_pars)
    for (m in seq_along(chol_pars)) {
      k <- chol_loc[m, 1L]; l <- chol_loc[m, 2L]
      out[m] <- if (is_diag[m]) log(L_full[k, k]) else L_full[k, l]
    }
    return(out)
  }

  # Selective correlate: minimise Q-function over the constrained cholPars.
  Q_obj <- function(chol_in) {
    v <- setNames(as.numeric(chol_in), chol_pars)
    L <- omega$buildL(v)
    Linv      <- forwardsolve(L, diag(K))
    Omega_inv <- crossprod(Linv)
    logdetO   <- 2 * sum(log(diag(L)))
    val       <- 0.5 * N * (logdetO + sum(Omega_inv * S))

    Mmat <- Omega_inv %*% S %*% Omega_inv
    ML   <- Mmat %*% L
    gr   <- setNames(numeric(length(v)), chol_pars)
    for (m in seq_along(chol_pars)) {
      k <- chol_loc[m, 1L]; l <- chol_loc[m, 2L]
      if (is_diag[m]) {
        gr[m] <- N * (1 - ML[k, k] * L[k, k])
      } else {
        gr[m] <- -N * ML[k, l]
      }
    }
    list(value = val, gradient = gr)
  }

  fd_hess <- function(chol_in, h = 1e-5) {
    M_par <- length(chol_in)
    H <- matrix(0, M_par, M_par, dimnames = list(chol_pars, chol_pars))
    for (i in seq_len(M_par)) {
      xp <- chol_in; xp[i] <- xp[i] + h
      xm <- chol_in; xm[i] <- xm[i] - h
      H[, i] <- (Q_obj(xp)$gradient - Q_obj(xm)$gradient) / (2 * h)
    }
    (H + t(H)) / 2
  }

  Q_objlist <- function(chol_in, ...) {
    qo <- Q_obj(chol_in)
    objlist(value = qo$value, gradient = qo$gradient, hessian = fd_hess(chol_in))
  }

  L0 <- t(chol(S))
  start <- setNames(numeric(length(chol_pars)), chol_pars)
  for (m in seq_along(chol_pars)) {
    k <- chol_loc[m, 1L]; l <- chol_loc[m, 2L]
    start[m] <- if (is_diag[m]) log(L0[k, k]) else L0[k, l]
  }
  fit <- suppressMessages(trust(Q_objlist, parinit = start,
                                rinit = 1, rmax = 10,
                                iterlim = 50,
                                ftol = 1e-10, mtol = 1e-10))
  setNames(as.numeric(fit$argument), chol_pars)
}



#' Parameter names of an object
#'
#' Generic for extracting the parameter names carried by an object. See
#' [parnames.omegaspec] for the random-effects spec method.
#'
#' @param x An object.
#' @param ... Passed to methods.
#' @return Character vector of parameter names.
#' @export
parnames <- function(x, ...) UseMethod("parnames")

#' All random-effect parameter names of an omega
#'
#' @description
#' Returns the union of `cholPars` and all subject-level eta names in a single
#' character vector. Useful for constructing initial parameter vectors and for
#' partitioning the outer/inner parameter space in [EM].
#'
#' @param x An `omegaSpec` object.
#' @param what One of `"all"` (default, both eta and Cholesky parameters),
#'   `"eta"` (subject-level random effects only), or `"chol"` (Cholesky
#'   parameters only).
#' @param ... Ignored.
#' @return Character vector of parameter names.
#' @export
parnames.omegaspec <- function(x, what = c("all", "eta", "chol"), ...) {
  what <- match.arg(what)
  eta_names <- if (is.null(x$subjectEtas)) character(0) else as.vector(x$subjectEtas)
  switch(what,
         all  = c(x$cholPars, eta_names),
         eta  = eta_names,
         chol = x$cholPars)
}


# Session cache for the subject-independent Smolyak grid: the z-nodes and
# weights depend only on (K, level), so build once and reuse across the
# per-subject E-step. Also precomputes the log-weight pieces reused per subject.
.sparseGridCache <- new.env(parent = emptyenv())

.getSparseGrid <- function(K, level) {
  key <- paste0(K, "_", level)
  hit <- .sparseGridCache[[key]]
  if (!is.null(hit)) return(hit)
  raw <- sparseGridGH(as.integer(K), as.integer(level))
  grid <- list(zNodes  = raw$nodes,
               logAbsW = log(abs(raw$weights)),
               signs   = sign(raw$weights),
               z2Sum   = rowSums(raw$nodes^2),
               B       = nrow(raw$nodes))
  .sparseGridCache[[key]] <- grid
  grid
}


#' Build per-subject quadrature nodes for adaptive NLME marginal integration
#'
#' Given the posterior mode `etaHat` and negative-Hessian `Hi` of the joint
#' log-density at the mode (= inverse covariance of the local Gaussian
#' approximation), constructs sparse-grid Gauss-Hermite quadrature nodes in
#' eta-space along with the log-magnitudes of their augmented weights (with
#' signs tracked separately).
#'
#' The augmented weight `W_b` for the integral `int g(eta) deta` approximated
#' as `sum_b W_b * g(eta_b)` includes: the raw Smolyak-Gauss-Hermite weight
#' `w_b_GH`, the change-of-variable Jacobian `2^(K/2) / |det L_H|`, and the
#' `exp(z_b' z_b)` factor that un-does the implicit `exp(-z'z)` weight of the
#' physicists' GH rule. Combined: `log W_b = (K/2)*log(2) - log|det L_H| +
#' log|w_b_GH| + z_b' z_b`, with sign carried in `weightSigns`.
#'
#' @param etaHat Length-K numeric, posterior mode of the random effects.
#' @param Hi K x K positive-definite matrix, the joint's negative-Hessian
#'   at `etaHat`.
#' @param level Integer >= K, Smolyak depth. K+1 to K+3 is the useful range
#'   for typical NLME problems; higher levels exponentially raise node count.
#' @param pruneTol Non-negative scalar (default `Inf`, i.e. no pruning). When
#'   finite, nodes whose augmented log-weight sits more than `pruneTol` below
#'   the per-subject maximum are dropped, cutting the per-node ODE work at a
#'   bounded (`exp(-pruneTol)`-order) relative cost in the marginal likelihood.
#'   Values of 25-40 remove negligible-weight nodes with no material accuracy
#'   loss; validate against the un-pruned value for a new model.
#'
#' @return A list with components `etaNodes` (B x K matrix, batch-first),
#'   `logAbsWeights` (length-B numeric), `weightSigns` (length-B numeric in
#'   `{-1, +1}`), and the echoed `K` and `level`.
#'
#' @seealso [sparseGridGH()]
#' @export
makeSubjectNodes <- function(etaHat, Hi, level, pruneTol = Inf) {
  K <- length(etaHat)
  if (!is.matrix(Hi) || nrow(Hi) != K || ncol(Hi) != K)
    stop("`Hi` must be a K x K matrix matching length(etaHat).")
  if (!is.numeric(level) || length(level) != 1L || level < 1L)
    stop("`level` must be a positive integer scalar.")

  grid    <- .getSparseGrid(K, as.integer(level))
  z_nodes <- grid$zNodes
  B       <- grid$B

  R         <- chol(Hi)
  log_det_L <- sum(log(abs(diag(R))))

  z_solved  <- backsolve(R, t(z_nodes))
  etaNodes <- t(sqrt(2) * z_solved) +
               matrix(etaHat, B, K, byrow = TRUE)
  if (!is.null(names(etaHat))) colnames(etaNodes) <- names(etaHat)

  log_abs <- (K / 2) * log(2) - log_det_L + grid$logAbsW + grid$z2Sum
  signs   <- grid$signs

  # Optional weight-magnitude pruning (default off, pruneTol = Inf). A node's
  # contribution to the marginal integral is bounded by exp(logAbsWeights_b)
  # times the (mode-dominated) integrand peak, so nodes whose log-weight sits
  # `pruneTol` below the per-subject maximum contribute at most exp(-pruneTol)
  # of the peak and can be dropped with bounded relative error. This directly
  # shrinks the per-node ODE work in .normalEcmSubject.
  if (is.finite(pruneTol)) {
    keep <- log_abs >= (max(log_abs) - pruneTol)
    if (!all(keep)) {
      etaNodes <- etaNodes[keep, , drop = FALSE]
      log_abs  <- log_abs[keep]
      signs    <- signs[keep]
    }
  }

  list(etaNodes      = etaNodes,
       logAbsWeights = log_abs,
       weightSigns   = signs,
       K             = K,
       level         = as.integer(level))
}


## Per-subject metadata shared by the quadrature, SAEM and Bayesian
## paths. Built from the omega spec and the model pieces, so it belongs
## with the mixed-effects layer even though the sampler also consumes it.
.buildBayesSubjectMeta <- function(omegaSpec, initFull, prdfn, data,
                                   errfn        = NULL,
                                   innerControl = list(),
                                   trustControl = list()) {
  if (is.null(errfn)) errfn <- .normalStaticErr(data)

  K        <- omegaSpec$K
  subjects <- rownames(omegaSpec$subjectEtas)
  if (is.null(subjects))
    stop(".buildBayesSubjectMeta: omegaSpec needs subject expansion ",
         "(omega(..., subjects = ...)).")
  N <- length(subjects)
  outer_names   <- names(initFull)
  eta_names_all <- as.vector(omegaSpec$subjectEtas)
  pars_full_names <- c(outer_names, eta_names_all)
  eta_idx_global  <- matrix(0L, N, K)
  eta_names_list  <- vector("list", N)
  for (i in seq_len(N)) {
    nms <- omegaSpec$subjectEtas[i, ]
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
    pars_full_names   = pars_full_names)

  ic <- modifyList(list(rinit = 1, rmax = 10, iterlim = 30,
                        ftol = 1e-7, mtol = 1e-7,
                        eigen_floor_relative = 1e-10),
                   innerControl)
  oc <- modifyList(list(rinit = 1, rmax = 10, iterlim = 200,
                        ftol = 1e-7, mtol = 1e-7),
                   trustControl)

  pars_probe <- setNames(numeric(length(pars_full_names)), pars_full_names)
  pars_probe[outer_names] <- initFull
  fast_meta <- .normalFastMeta(prdfn, errfn, data, subjects,
                              eta_names_list, pars_full_names, pars_probe)
  om_meta <- list(
    chol_pars = omegaSpec$cholPars,
    chol_loc  = matrix(as.integer(omegaSpec$cholLoc), ncol = 2, dimnames = NULL),
    is_diag   = as.logical(omegaSpec$isDiag))
  subject_meta$fast_meta  <- fast_meta
  subject_meta$omega_meta <- om_meta

  list(subjectMeta  = subject_meta,
       innerControl = ic,
       outerControl = oc,
       errfn        = errfn,
       N            = N,
       K            = K,
       outerNames   = outer_names)
}

.makeSubjectEtaObj <- function(i, subjectMeta, modelCb, errCb,
                                parsFullTemplate, etaIndices) {
  meta_i  <- subjectMeta$fast_meta[[i]]
  eta_idx <- as.integer(etaIndices)
  fn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE,
                 conditions = NULL, env = NULL, .pars_full = NULL,
                 .Omega_inv = NULL, .Omega_log_det = NULL) {
    eta_i <- ..1
    p_full <- if (is.null(.pars_full)) parsFullTemplate else .pars_full
    p_full[eta_idx] <- as.numeric(eta_i)
    kr <- focei_eval_one_subject(
      model_cb       = modelCb,
      err_cb         = errCb,
      pars_full      = p_full,
      fixed          = NULL,
      meta_i         = meta_i,
      eta_block      = as.numeric(eta_i),
      Omega_inv      = .Omega_inv,
      Omega_log_det  = .Omega_log_det)
    out <- structure(list(value = kr$value, gradient = kr$gradient,
                          hessian = kr$hessian),
                     class = c("objlist", "list"))
    attr(out, "env") <- env
    out
  }
  class(fn) <- c("objfn", "fn")
  attr(fn, "parameters") <- meta_i$eta_names
  fn
}

# Bake dots + per-call extras into a closure callable as
# objfun(pars, deriv, deriv2). `extra` threads per-block locals
# (.pars_full / .Omega_inv for the joint NLME path) that must be
# re-evaluated by the caller's R closure each iteration.
.bake_objfun <- function(raw_objfun, dots = list(), extra = list()) {
  force(raw_objfun)
  force(dots)
  function(pars, deriv = TRUE, deriv2 = FALSE) {
    args <- c(list(pars = pars), dots, extra,
              list(deriv = deriv, deriv2 = deriv2))
    do.call(raw_objfun, args)
  }
}
