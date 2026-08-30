## .fitLaplace: regularised NLME fits with a marginally-estimated penalty strength.
##
## .fitLaplace reframes the L1-penalised "which parameters are individual-specific"
## fit as a nonlinear mixed-effects model with a Laplace random-effect density
## (notes/laplace_nlme_theory.Rmd): the per-subject deviations are random
## effects, the single penalty strength `lambda` is a parameter of that density
## estimated by marginal maximum likelihood (no scan, no per-lambda multistart),
## the outer marginal is evaluated with FOCE(I) and the inner conditional mode
## with trustL1.
##
## The penalty structure is declared with penaltyL1() (penaltyClass.R)
## and stamped onto the objective by constraintL1(); the fit is assembled exactly
## like an omega() nlme fit, i.e. (via the exported EM()/emInit() dispatchers)
##   obj  <- normL2(dlist, g*x*p, errmodel = e) + constraintL1(pen)
##   init <- emInit(structural, pen)
##   fit  <- EM(obj, init, method = "focei")


## Recover the model pieces + penalty spec from a composed .fitLaplace objective.
## Sibling of .normalReconstruct() (nlmeNormal.R); requires a penaltySpec (not omegaSpec).
.laplaceReconstruct <- function(obj) {
  if (!inherits(obj, "objfn"))
    stop(".fitLaplace: `obj` must be an objfn.", call. = FALSE)
  prdfn   <- attr(obj, "prdfn",      exact = TRUE)
  data    <- attr(obj, "data",       exact = TRUE)
  errfn   <- attr(obj, "errfn",      exact = TRUE)
  penalty <- attr(obj, "penaltySpec", exact = TRUE)
  if (is.null(prdfn) || is.null(data))
    stop(".fitLaplace: could not recover the prediction function and data from ",
         "`obj`. Build it as ",
         "normL2(data, g*x*p, errmodel = e) + constraintL1(pen).",
         call. = FALSE)
  if (is.null(penalty))
    stop(".fitLaplace: `obj` carries no penalty. Add ",
         "+ constraintL1(penaltyL1(..., subjects = ...)).", call. = FALSE)
  if (is.null(penalty$subjectEtas))
    stop(".fitLaplace: the penalty in `obj` has no subject expansion. Build it with ",
         "penaltyL1(..., subjects = ...).", call. = FALSE)
  list(prdfn = prdfn, data = data, errfn = errfn, penalty = penalty)
}


#' Assemble a starting parameter vector for a regularised NLME fit
#'
#' @description
#' The L1 counterpart of [emInit]. Combines the named structural (and
#' error-model) starting values with the single penalty-strength parameter
#' named by the `penaltySpec`. The per-subject deviations are inner parameters
#' (solved by [trustL1] each outer iteration) and are not part of the outer
#' start vector, exactly as the omega Cholesky parameters are the only spread
#' parameters in an [emInit] vector.
#'
#' @param structural Named numeric of structural and error-model starting
#'   values, e.g. `c(tka = 0, tv = 3, tcl = 0.7, log_sigma_add = log(0.2))`.
#' @param penalty A `penaltySpec` (from [penaltyL1], possibly
#'   combined with `+`).
#' @param lambda Starting value(s) for the penalty strength(s) (natural, positive
#'   scale). A single number is shared by every strength name of the spec (one for
#'   a global penalty, several for a per-parameter one); a named vector sets each
#'   strength individually. Default `1`.
#' @return A named numeric suitable as the `init` of [EM].
#' @seealso [EM], [penaltyL1], [emInit]
#' @noRd
.initLaplace <- function(structural, penalty, lambda = 1) {
  if (!inherits(penalty, "penaltyspec"))
    stop("`penalty` must be built by penaltyL1().", call. = FALSE)
  if (is.null(names(structural)) || any(!nzchar(names(structural))))
    stop("`structural` must be a fully named numeric vector.", call. = FALSE)
  lam_names <- penalty$lambdaName
  if (!is.numeric(lambda) || any(lambda <= 0))
    stop("`lambda` must be positive.", call. = FALSE)
  if (length(lambda) == 1L && is.null(names(lambda))) {
    lam <- setNames(rep(as.numeric(lambda), length(lam_names)), lam_names)
  } else {
    if (is.null(names(lambda)) || !setequal(names(lambda), lam_names))
      stop("A vector-valued `lambda` must be named to match the spec's ",
           "strength names (", paste(lam_names, collapse = ", "), ").",
           call. = FALSE)
    lam <- setNames(as.numeric(lambda[lam_names]), lam_names)
  }
  dup <- intersect(names(structural), names(lam))
  if (length(dup))
    warning(".initLaplace: `structural` already contains penalty parameter(s) '",
            paste(dup, collapse = ", "), "'; the supplied value is kept.",
            call. = FALSE)
  c(structural, lam[setdiff(names(lam), names(structural))])
}


## Per-subject penalty targets for the factorising path (single/fused). Returns
## a list, one named numeric per subject, aligned to that subject's eta names
## (value = the reference the deviation is pulled toward; 0 for fused).
.laplaceFactorisingTargets <- function(penalty) {
  ## Defensive: clustered is dispatched to the coupled MAP-EM path
  ## (.clusterMap) before this is reached; this guard should never fire.
  kinds <- vapply(penalty$blocks, `[[`, "", "kind")
  if (any(kinds == "clustered"))
    stop(".fitLaplace: the factorising-target helper does not handle 'clustered' ",
         "(the coupled path is .clusterMap; see R/L1Clustering.R).",
         call. = FALSE)
  par_target <- list()
  for (b in penalty$blocks) for (par in b$params)
    par_target[[par]] <- if (identical(b$kind, "fused")) 0 else as.numeric(b$target[[par]])

  se       <- penalty$subjectEtas
  subjects <- rownames(se)
  setNames(lapply(subjects, function(s) {
    cols <- colnames(se)                         # eta_<par>
    pars <- sub(paste0("^", penalty$prefix, "_"), "", cols)
    setNames(as.numeric(unlist(par_target[pars])), as.character(se[s, ]))
  }), subjects)
}


## d diag(H)/d{theta,eta} of one subject's data Hessian, in closed form.
##
## H is the Gauss-Newton block, H = c * G' diag(1/sigma^2) G with
## G[t,k] = d yhat_t / d eta_k, so its derivative needs only second
## derivatives of the prediction -- no third ones, and no differences:
##
##   d H[k,k]/d p = c * ( 2 sum_t G[t,k] W[t,k,p] / sigma_t^2
##                        - 2 sum_t G[t,k]^2 dsigma_p,t / sigma_t^3 )
##
## with W = d2 yhat / d eta_k d p. `c` is recovered from the objective's own
## Hessian rather than assumed, which also checks the reconstruction.
##
## Returns list(dth, deta, Hd) with Hd the reconstructed diagonal, or NULL when
## the chain cannot supply deriv2.
.laplaceHdDerivs <- function(obj_s, full_mode, eta_names_i, struct_names,
                             cores = getOption("dMod.cores", 1L)) {
  prd    <- attr(obj_s, "prdfn",  exact = TRUE)
  errfn  <- attr(obj_s, "errfn",  exact = TRUE)
  dat    <- attr(obj_s, "data",   exact = TRUE)
  timesD <- attr(obj_s, "timesD", exact = TRUE)
  cn     <- attr(obj_s, "conditions")
  if (is.null(prd) || is.null(dat) || is.null(timesD) || length(cn) != 1L)
    return(NULL)

  pred <- tryCatch(
    prd(timesD, full_mode, deriv = TRUE, deriv2 = TRUE, conditions = cn,
        cores = cores)[[cn]],
    error = function(e) NULL)
  if (is.null(pred)) return(NULL)
  d1 <- attr(pred, "deriv"); d2 <- attr(pred, "deriv2")
  if (is.null(d1) || is.null(d2)) return(NULL)

  di      <- dat[[cn]]
  times_i <- di$time
  names_i <- as.character(di$name)
  tix <- .matchNum(times_i, as.numeric(pred[, "time"]))
  if (anyNA(tix)) return(NULL)

  avail <- dimnames(d1)[[3]]
  eta_a <- intersect(eta_names_i, avail)
  th_a  <- intersect(struct_names, avail)
  if (!length(eta_a)) return(NULL)

  ## sigma per data point, and its derivatives when an error model supplies them
  sigma_i <- di$sigma
  dsig <- NULL
  if (!is.null(errfn)) {
    pinner <- getParameters(pred)
    ei <- tryCatch(errfn(out = pred, pars = pinner, conditions = cn)[[cn]],
                   error = function(e) NULL)
    if (is.null(ei)) return(NULL)
    te <- .matchNum(times_i, as.numeric(ei[, "time"]))
    ne <- match(names_i, colnames(ei))
    if (anyNA(te) || anyNA(ne)) return(NULL)
    sigma_i <- ifelse(is.na(sigma_i), ei[cbind(te, ne)], sigma_i)
    de <- attr(ei, "deriv")
    if (!is.null(de)) {
      oe <- match(names_i, dimnames(de)[[2]])
      eset <- dimnames(de)[[3]]
      dsig <- list(idx = cbind(te, oe), set = eset, arr = de)
    }
  }
  if (anyNA(sigma_i)) return(NULL)

  Tn <- length(times_i); Kn <- length(eta_a)
  G <- matrix(0, Tn, Kn, dimnames = list(NULL, eta_a))
  for (jr in seq_len(Tn)) G[jr, ] <- d1[tix[jr], names_i[jr], eta_a]
  Sinv <- 1 / sigma_i^2

  ## scale factor from the objective's own Hessian diagonal
  Hd_ref <- diag(obj_s(full_mode, deriv = TRUE)$hessian[eta_a, eta_a, drop = FALSE])
  base   <- colSums(G^2 * Sinv)
  ok     <- abs(base) > 0
  if (!any(ok)) return(NULL)
  cfac <- median(Hd_ref[ok] / base[ok])
  if (!is.finite(cfac) ||
      max(abs(cfac * base - Hd_ref)) > 1e-3 * max(1, max(abs(Hd_ref))))
    return(NULL)   # not the plain GN form; fall back to the caller's path

  dfor <- function(cols) {
    if (!length(cols)) return(matrix(0, Kn, 0, dimnames = list(eta_a, NULL)))
    out <- matrix(0, Kn, length(cols), dimnames = list(eta_a, cols))
    for (jc in seq_along(cols)) {
      pcol <- cols[jc]
      W <- matrix(0, Tn, Kn)
      for (jr in seq_len(Tn)) W[jr, ] <- d2[tix[jr], names_i[jr], eta_a, pcol]
      term <- 2 * colSums(G * W * Sinv)
      if (!is.null(dsig) && pcol %in% dsig$set) {
        ds <- dsig$arr[cbind(dsig$idx, match(pcol, dsig$set))]
        term <- term - 2 * colSums(G^2 * (ds / sigma_i^3))
      }
      out[, jc] <- cfac * term
    }
    out
  }

  list(dth = dfor(th_a), deta = dfor(eta_a), Hd = cfac * base, eta = eta_a,
       theta = th_a)
}


## Build a data-only per-subject objfun eta -> objlist. `resObj_s` MUST be a
## normL2 carrying this subject's data alone: each subject needs its own marginal,
## and a shared all-subject objfn returns the sum. The outer (structural + error)
## parameters are held at `outer_struct`; all etas start at 0 and only subject i's
## block varies. Value/gradient/hessian follow dMod's -2 log L convention, so they
## feed .laplaceSubjectMarginal() directly.
.laplaceSubjectObjfun <- function(resObj_s, outer_struct, all_eta_names,
                                  eta_names_i) {
  base_full <- c(outer_struct,
                 setNames(rep(0, length(all_eta_names)), all_eta_names))
  pars_at <- function(v) {                 # v positional, in eta_names_i order
    full <- base_full
    full[eta_names_i] <- v
    full
  }
  at <- function(p) resObj_s(pars_at(p[eta_names_i]), deriv = TRUE)
  f <- function(p, ...) {
    o <- at(p)
    list(value    = o$value,
         gradient = o$gradient[eta_names_i],
         hessian  = o$hessian[eta_names_i, eta_names_i, drop = FALSE])
  }
  # The outer gradient wants the structural rows this drops, at the very point
  # the mode evaluation already visits. Hand the whole objlist back so
  # .laplaceOuterEval does not solve it a second time.
  attr(f, "full") <- at
  # The lock-step driver addresses coordinates by position and batches the
  # objectives itself, so it needs the pieces rather than the closure.
  attr(f, "pars_at")  <- pars_at
  attr(f, "resobj")   <- resObj_s
  attr(f, "etaNames") <- eta_names_i
  f
}


## Evaluate the outer marginal objective sum_i -2 log L_i at theta = (structural,
## error, lambda). With grad = TRUE also returns the outer gradient. Two modes:
##   secondOrderCorrection = FALSE (FOCE): d(-2logL)/d theta = data gradient at the mode plus
##     the first-moment (Fisher) cross term Hcross (E[eta]-etahat); lambda closed
##     form -2 sum_k (1/lambda - E|eta_k|). Misses the volume + mode-shift terms.
##   secondOrderCorrection = TRUE (FOCEI): the EXACT total derivative of the Laplace-approx
##     marginal value, d M/d theta = explicit(eta frozen) + implicit mode shift
##     (d M/d etahat)(d etahat/d theta), including lambda. The volume term's
##     d Hd/d{theta,eta} comes from .laplaceHdDerivs() in closed form, so the
##     chain has to be deriv2-capable.
## With hess = TRUE (requires grad) it also returns the structural Gauss-Newton
## Hessian sum_i o_full$hessian[struct, struct] (the data Hessian block at each
## subject's mode) as $hessian; the outer trust() CM-step consumes it. This is a
## GN approximation of the marginal Hessian (drops the d eta_hat / d theta coupling
## and the log-volume curvature), exactly as Gaussian FOCE's outer step does.
## With collect = TRUE also returns the per-subject solves for the readout.
.laplaceOuterEval <- function(theta, rec, resObjList, targets_by_subject,
                          control = list(), collect = FALSE, warm = NULL,
                          grad = FALSE, hess = FALSE, correction = FALSE) {
  penalty      <- rec$penalty
  lam_name     <- penalty$lambdaName                 # >= 1 distinct strength names
  lambda_vec   <- setNames(as.numeric(theta[lam_name]), lam_name)
  outer_struct <- theta[setdiff(names(theta), lam_name)]
  struct_names <- names(outer_struct)
  se       <- penalty$subjectEtas
  subjects <- rownames(se)
  all_eta  <- as.vector(se)
  ## strength name for each candidate coordinate, in the eta-column (= short) order
  lamByCoord <- unname(penalty$lambdaByParam[penalty$short])

  inner <- .trustControl(
    list(rinit = control$rinit %||% 1, rmax = control$rmax %||% 10,
         iterlim = control$iterlim %||% 50L, ftol = 1e-8, mtol = 1e-8),
    control$inner, optimizer = trustL1, label = "control$inner")

  total   <- 0
  gStruct <- setNames(numeric(length(struct_names)), struct_names)
  gLambda <- setNames(numeric(length(lam_name)), lam_name)
  hStruct <- if (grad && hess)
    matrix(0, length(struct_names), length(struct_names),
           dimnames = list(struct_names, struct_names)) else NULL
  perSub  <- if (collect) vector("list", length(subjects)) else NULL

  etaNames <- lapply(subjects, function(s) as.character(se[s, ]))
  objfuns  <- lapply(seq_along(subjects), function(ii)
    .laplaceSubjectObjfun(resObjList[[subjects[ii]]], outer_struct, all_eta,
                          etaNames[[ii]]))
  ## warm-start the inner conditional-mode solve from the last accepted mode:
  ## between outer steps eta_hat barely moves, so this cuts inner iterations.
  eta0List <- lapply(seq_along(subjects), function(ii) {
    s <- subjects[ii]; en <- etaNames[[ii]]
    e0 <- if (!is.null(warm) && !is.null(warm[[s]])) warm[[s]]
          else setNames(rep(0, length(en)), en)
    e0[en]
  })
  tgList <- lapply(seq_along(subjects), function(ii)
    targets_by_subject[[subjects[ii]]][etaNames[[ii]]])
  ## per-coordinate strength (each candidate coordinate uses its parameter's
  ## lambda); a single global lambda stays a scalar (bit-identical backward
  ## compat), else a per-coordinate named vector for trustL1 / normalLaplace.
  lamList <- lapply(etaNames, function(en)
    if (length(lam_name) == 1L) lambda_vec[[1L]]
    else setNames(unname(lambda_vec[lamByCoord]), en))
  smList <- .laplaceSubjectMarginals(objfuns, eta0List, tgList, lamList, inner,
                                     lockstep = !isFALSE(control$lockstep) &&
                                       !isFALSE(getOption("dMod.laplace.lockstep", TRUE)))

  for (ii in seq_along(subjects)) {
    s           <- subjects[ii]
    eta_names_i <- etaNames[[ii]]
    lambda_k    <- lamList[[ii]]
    tg          <- tgList[[ii]]
    sm          <- smList[[ii]]
    if (!is.null(warm)) warm[[s]] <- sm$etahat
    total <- total + sm$value

    if (grad) {
      ## only the correction branch still needs the point itself
      o_full    <- sm$full
      full_mode <- if (is.null(o_full) || correction) {
        f <- c(outer_struct, setNames(rep(0, length(all_eta)), all_eta))
        f[eta_names_i] <- sm$etahat
        f
      }
      if (is.null(o_full)) o_full <- resObjList[[s]](full_mode, deriv = TRUE)
      Hcross  <- o_full$hessian[struct_names, eta_names_i, drop = FALSE] # d2D/dtheta deta
      Hdd     <- o_full$hessian[eta_names_i, eta_names_i, drop = FALSE]  # d2D/deta deta
      gd      <- sm$grad; Hd <- sm$Hd; nl <- sm$nl
      K       <- length(eta_names_i)

      if (!correction) {
        ## Plain FOCE: mode-value gradient plus the first-moment (Fisher) cross
        ## term E_post[d di$value/d theta] ~= d value/d theta + Hcross (E[eta]-etahat).
        ## Misses the normal-Laplace volume term and the mode shift (the FOCEI
        ## branch below); the error-parameter (interaction) gap is several percent.
        delta   <- nl$Eeta - (sm$etahat - tg)
        gStruct <- gStruct + o_full$gradient[struct_names] +
                   as.numeric(Hcross %*% delta)
        ## each coordinate's lambda score -2 (1/lambda_k - E|eta_k|) accrues to
        ## that coordinate's strength name.
        contribL <- (-2) * (1 / lambda_k - nl$Eabs)
        gLambda  <- gLambda +
          vapply(lam_name, function(j) sum(contribL[lamByCoord == j]), 0.0)
      } else {
        ## FOCEI: the EXACT total derivative of the Laplace-approx marginal value
        ##   M = value(etahat) - (1/2) sum g_k^2/Hd_k - 2 sum log I_k,
        ## with a_k = Hd_k/2, m_k = etahat_k - g_k/Hd_k (both functions of the mode
        ## etahat(theta)). d M/d theta = explicit(eta frozen) + implicit
        ## (d M/d etahat) (d etahat/d theta). The volume term enters through
        ## d Hd_k/d{theta,eta} and d g_k/d{theta,eta}; the mode shift closes the
        ## FOCE gap (structural + lambda).
        dI_dm <- nl$dlogI_dm; dI_da <- nl$dlogI_da; dI_dl <- nl$dlogI_dlambda
        ## d Hd/d{theta,eta} in closed form -- see .laplaceHdDerivs(). The
        ## Gauss-Newton Hessian is a product of first derivatives, so its
        ## derivative needs only second derivatives of the prediction.
        .hd <- .laplaceHdDerivs(resObjList[[s]], full_mode, eta_names_i,
                                struct_names)
        if (is.null(.hd))
          stop("secondOrderCorrection = TRUE needs d Hd/d{theta,eta}, computed ",
               "in closed form from second derivatives of the prediction. ",
               "Build the chain with deriv2 = TRUE (odemodel/Y/P) and use ",
               "plain Gaussian residuals, or leave secondOrderCorrection off.",
               call. = FALSE)
        dHd_dth <- matrix(0, K, length(struct_names),
                          dimnames = list(eta_names_i, struct_names))
        dHd_dth[.hd$eta, .hd$theta] <- .hd$dth
        dHd_deta <- matrix(0, K, K, dimnames = list(eta_names_i, eta_names_i))
        dHd_deta[.hd$eta, .hd$eta] <- .hd$deta

        ## Active set: only coords off their kink (etahat != target) shift with
        ## theta; pinned (shared) coords have d etahat/d theta = 0.
        free <- which(abs(sm$etahat - tg) > 1e-6 * pmax(abs(tg), 1))

        ## d M / d etahat_l  (the mode-value gradient the shift multiplies)
        dM_deta <- numeric(K)
        for (l in seq_len(K)) {
          dgk    <- Hdd[, l]                        # d g_k / d eta_l
          dHk    <- dHd_deta[, l]                   # d Hd_k / d eta_l
          e_l    <- numeric(K); e_l[l] <- 1
          dmk    <- e_l - dgk / Hd + (gd / Hd^2) * dHk
          dak    <- 0.5 * dHk
          dM_deta[l] <- gd[l] -
            sum((gd / Hd) * dgk - 0.5 * (gd^2 / Hd^2) * dHk) -
            2 * sum(dI_da * dak + dI_dm * dmk)
        }

        ## d etahat_free / d theta from inner stationarity
        ## g_l(etahat, theta) + lambda_{p(l)}*sign(etahat_l - t_l) = 0.
        if (length(free)) {
          Hff_inv  <- solve(Hdd[free, free, drop = FALSE])
          dEta_dth <- -Hff_inv %*% t(Hcross[, free, drop = FALSE])   # |free| x Q
          sgn_free <- sign(sm$etahat - tg)[free]
        }

        ## explicit d M / d theta_q | eta frozen, plus the implicit mode shift
        for (iq in seq_along(struct_names)) {
          q     <- struct_names[iq]
          dg_dq <- Hcross[q, ]                       # d g_k / d theta_q
          dH_dq <- dHd_dth[, iq]                     # d Hd_k / d theta_q
          dm_dq <- -dg_dq / Hd + (gd / Hd^2) * dH_dq
          da_dq <- 0.5 * dH_dq
          expl  <- o_full$gradient[[q]] -
            sum((gd / Hd) * dg_dq - 0.5 * (gd^2 / Hd^2) * dH_dq) -
            2 * sum(dI_da * da_dq + dI_dm * dm_dq)
          impl  <- if (length(free)) sum(dM_deta[free] * dEta_dth[, iq]) else 0
          gStruct[[q]] <- gStruct[[q]] + expl + impl
        }
        ## each strength name j: explicit -2 sum_{k in j} dlogI_k/dlambda_j, plus
        ## its own mode shift (only coords using lambda_j perturb with it).
        for (j in lam_name) {
          in_j <- lamByCoord == j
          expl_j <- (-2) * sum(dI_dl[in_j])
          impl_j <- if (length(free)) {
            rhs_j <- sgn_free * (lamByCoord[free] == j)
            sum(dM_deta[free] * as.numeric(-Hff_inv %*% rhs_j))
          } else 0
          gLambda[[j]] <- gLambda[[j]] + expl_j + impl_j
        }
      }

      if (hess)
        hStruct <- hStruct +
          o_full$hessian[struct_names, struct_names, drop = FALSE]
    }
    if (collect) perSub[[ii]] <- sm
  }

  ## Reference conditions: data conditions carrying no random effect (the baseline
  ## line in the 2-line reference encoding). No eta, no marginalisation, no lambda
  ## term - they contribute their plain data -2 log L and structural derivatives,
  ## pinning the structural center (mu) that the penalised line deviates from.
  ref_conditions <- setdiff(names(rec$data), subjects)
  if (length(ref_conditions)) {
    full_ref <- c(outer_struct, setNames(rep(0, length(all_eta)), all_eta))
    refs <- .objEvalMany(resObjList[ref_conditions],
                         rep(list(full_ref), length(ref_conditions)),
                         deriv = grad)
    names(refs) <- ref_conditions
  }
  for (rc in ref_conditions) {
    o_ref    <- refs[[rc]]
    total    <- total + o_ref$value
    if (grad) {
      ## the reference prediction depends on only a subset of the structural
      ## names (fixed non-support etas never enter it); index that subset.
      sn <- intersect(struct_names, names(o_ref$gradient))
      gStruct[sn] <- gStruct[sn] + o_ref$gradient[sn]
      if (hess)
        hStruct[sn, sn] <- hStruct[sn, sn] +
          o_ref$hessian[sn, sn, drop = FALSE]
    }
  }

  gradient <- NULL
  if (grad) {
    gradient <- setNames(numeric(length(theta)), names(theta))
    gradient[struct_names] <- gStruct
    gradient[lam_name]     <- gLambda
  }
  if (collect) return(list(value = total, perSubject = perSub,
                           subjects = subjects, lambda = lambda_vec,
                           gradient = gradient))
  if (grad) {
    res <- list(value = total, gradient = gradient)
    if (hess) res$hessian <- hStruct
    return(res)
  }
  total
}


## SAEM cross-check backend for the factorising (penalty/fused) path.
##
## A self-contained stochastic-approximation EM that reuses .fitLaplace's per-subject
## data objfns (resObjList). Per outer iteration:
##   E-step  a short per-subject random-walk Metropolis draw of eta_i from the
##           L1-penalised posterior log p(eta) = -1/2 value_data(eta) -
##           lambda * sum_k |eta_k - t_k|. The kink needs no subgradient: MH uses
##           only the (unnormalised) log-density, not its gradient.
##   SA      a Robbins-Monro update of the sufficient statistic
##           S_abs = sum_i sum_k |eta_ik - t_k|, damped by gamma_k (1 in burn-in,
##           1/(k - nBurnin) in convergence).
##   CM-2    the closed-form lambda M-step lambda = Kn / S_abs (same as the FOCEI
##           path), the Laplace analogue of updateOmegaChol.
##   CM-1    a trust() structural M-step at the drawn etas (data + reference
##           conditions), RM-damped in the convergence phase.
## Reports a FOCE-comparable readout (argument, lambda, etaModes = posterior
## means, value = the FOCE marginal at the SAEM point) so method = "saem" and
## "focei" land side by side. See notes/laplace_nlme_theory.Rmd (SAEM estimator).
.laplaceSaem <- function(rec, resObjList, targets_by_subject, free, fixed,
                        lam_name, control, verbose) {
  penalty  <- rec$penalty
  se       <- penalty$subjectEtas
  subjects <- rownames(se)
  all_eta  <- as.vector(se)
  ## per-strength deviation counts + coordinate->strength map (eta-column order)
  lamByCoord <- unname(penalty$lambdaByParam[penalty$short])
  parsPerLam <- table(penalty$lambdaByParam)
  countByLam <- setNames(nrow(se) * as.integer(parsPerLam[lam_name]), lam_name)
  lam_free_names <- intersect(lam_name, names(free))
  lambdaMax      <- if (is.null(control$lambdaMax)) 1e4 else control$lambdaMax
  struct_names   <- setdiff(names(free), lam_name)
  ref_conditions <- setdiff(names(rec$data), subjects)
  eta_names_by   <- setNames(lapply(subjects, function(s) as.character(se[s, ])),
                             subjects)
  target_by      <- setNames(lapply(subjects, function(s)
    targets_by_subject[[s]][eta_names_by[[s]]]), subjects)

  sc <- modifyList(list(nBurnin = 150L, nEM = 150L, nMcmc = 10L, stepsize = 0.4,
                        cm1 = list()),
                   if (is.null(control$saem)) list() else control$saem)
  cm1     <- .trustControl(list(rinit = 1, rmax = 10, iterlim = 30L,
                                ftol = 1e-6, mtol = 1e-6),
                           sc$cm1, label = "control$saem$cm1")
  nBurnin <- as.integer(sc$nBurnin); nEM <- as.integer(sc$nEM)
  nMcmc   <- as.integer(sc$nMcmc)

  struct_cur <- free[struct_names]
  lambda_cur <- setNames(as.numeric(c(free, fixed)[lam_name]), lam_name)
  base0      <- setNames(rep(0, length(all_eta)), all_eta)
  etas       <- setNames(lapply(subjects, function(s)
    setNames(rep(0, length(eta_names_by[[s]])), eta_names_by[[s]])), subjects)
  step       <- setNames(rep(sc$stepsize, length(subjects)), subjects)

  ## data -2 log L for subject s at (structural, this subject's etas)
  data_value <- function(outer_struct, s, eta_s) {
    full <- c(outer_struct, base0); full[eta_names_by[[s]]] <- eta_s
    resObjList[[s]](full, deriv = FALSE)$value
  }
  ## structural M-step objective at frozen etas (data + reference lines)
  cm1_obj <- function(theta_struct) {
    outer_struct <- c(setNames(as.numeric(theta_struct), struct_names), fixed)
    total <- 0
    g <- setNames(numeric(length(struct_names)), struct_names)
    H <- matrix(0, length(struct_names), length(struct_names),
                dimnames = list(struct_names, struct_names))
    for (s in subjects) {
      full <- c(outer_struct, base0); full[eta_names_by[[s]]] <- etas[[s]]
      o <- resObjList[[s]](full, deriv = TRUE)
      total <- total + o$value
      g <- g + o$gradient[struct_names]
      H <- H + o$hessian[struct_names, struct_names, drop = FALSE]
    }
    for (rc in ref_conditions) {
      full <- c(outer_struct, base0)
      o  <- resObjList[[rc]](full, deriv = TRUE)
      total <- total + o$value
      sn <- intersect(struct_names, names(o$gradient))
      g[sn] <- g[sn] + o$gradient[sn]
      H[sn, sn] <- H[sn, sn] + o$hessian[sn, sn, drop = FALSE]
    }
    objlist(value = total, gradient = g, hessian = H)
  }

  S_abs <- NULL
  etaMeanAcc <- setNames(lapply(subjects, function(s)
    setNames(numeric(length(eta_names_by[[s]])), eta_names_by[[s]])), subjects)
  nAcc <- 0L
  rows <- list(); prev_struct <- struct_cur; total_it <- nBurnin + nEM
  for (k in seq_len(total_it)) {
    gamma <- if (k <= nBurnin) 1.0 else 1.0 / (k - nBurnin)
    outer_struct <- c(struct_cur, fixed)

    ## E-step: per-subject random-walk Metropolis on the L1-penalised posterior.
    for (s in subjects) {
      t_s   <- target_by[[s]]
      lam_s <- unname(lambda_cur[lamByCoord])   # per-coordinate strength
      logp <- function(eta)
        -0.5 * data_value(outer_struct, s, eta) - sum(lam_s * abs(eta - t_s))
      cur <- etas[[s]]; lp_cur <- logp(cur); acc <- 0L
      for (mm in seq_len(nMcmc)) {
        prop    <- cur + rnorm(length(cur), 0, step[[s]])
        lp_prop <- logp(prop)
        if (is.finite(lp_prop) && log(runif(1)) < lp_prop - lp_cur) {
          cur <- prop; lp_cur <- lp_prop; acc <- acc + 1L
        }
      }
      etas[[s]] <- cur
      a <- acc / nMcmc
      step[[s]] <- min(max(step[[s]] * exp(0.3 * (a - 0.4)), 0.02), 5)
    }

    ## SA sufficient statistic (per strength) + closed-form lambda (CM-2).
    absMat    <- matrix(vapply(subjects,
      function(s) abs(etas[[s]] - target_by[[s]]), numeric(penalty$K)),
      nrow = penalty$K)                                                # K x nSub
    s_abs_new <- vapply(lam_name,
      function(j) sum(absMat[lamByCoord == j, , drop = FALSE]), 0.0)
    S_abs     <- if (is.null(S_abs)) s_abs_new else S_abs + gamma * (s_abs_new - S_abs)
    if (length(lam_free_names))
      lambda_cur[lam_free_names] <-
        pmin(countByLam / S_abs, lambdaMax)[lam_free_names]

    ## CM-1: structural M-step at the drawn etas, RM-damped in convergence.
    fit <- suppressMessages(do.call(trust, modifyList(cm1, list(
      objfun = function(th, ...) cm1_obj(th), parinit = struct_cur))))
    struct_new <- fit$argument
    struct_cur <- if (k <= nBurnin) struct_new
                  else struct_cur + gamma * (struct_new - struct_cur)

    if (k > nBurnin) {
      nAcc <- nAcc + 1L
      for (s in subjects) etaMeanAcc[[s]] <- etaMeanAcc[[s]] + etas[[s]]
    }
    deltaStruct <- max(abs(struct_cur - prev_struct))
    rows[[length(rows) + 1L]] <- data.frame(
      iter = k, phase = if (k <= nBurnin) "burnin" else "converge",
      gamma = gamma, lambda = lambda_cur[[1]], deltaStruct = deltaStruct)
    if (verbose && (k %% 50L == 0L || k == total_it))
      cat(sprintf("  EM(saem) iter %d/%d (%s) lambda=%s |dtheta|=%.2e\n",
                  k, total_it, if (k <= nBurnin) "burnin" else "converge",
                  paste(formatC(lambda_cur, format = "g", digits = 4),
                        collapse = "/"), deltaStruct))
    prev_struct <- struct_cur
  }

  ## posterior-mean etaModes (convergence phase); FOCE-marginal value at the point.
  etaMode <- if (nAcc > 0L)
    lapply(subjects, function(s) etaMeanAcc[[s]] / nAcc) else etas
  etaModes <- do.call(rbind, lapply(etaMode, unname))
  rownames(etaModes) <- subjects; colnames(etaModes) <- penalty$eta

  theta_hat <- c(struct_cur, fixed); theta_hat[lam_name] <- lambda_cur[lam_name]
  final <- .laplaceOuterEval(theta_hat, rec, resObjList, targets_by_subject,
                         list(iterlim = 200L), collect = TRUE)
  tail_rows <- utils::tail(do.call(rbind, rows), min(20L, nEM))
  conv <- nEM > 0L && mean(tail_rows$deltaStruct) < 1e-3

  out <- list(
    argument   = theta_hat,
    value      = final$value,
    lambda     = lambda_cur,
    lambda_mle = lambda_cur,
    etaModes   = etaModes,
    sparsity   = abs(etaModes) < 1e-6,
    converged  = conv,
    iterations = total_it,
    stageTrace = do.call(rbind, rows),
    method     = "saem",
    prior      = "penaltyL1",
    penalty    = penalty,
    prdfn      = rec$prdfn,
    data       = rec$data,
    errfn      = rec$errfn)
  class(out) <- c("em", "list")
  out
}


#' Regularised NLME fit with a marginally-estimated penalty strength
#'
#' @description
#' Fits an ODE mixed-effects model in which the per-subject deviations of the
#' candidate parameters carry a Laplace (L1) random-effect density, and the
#' single penalty strength `lambda` is estimated by marginal maximum likelihood
#' (no scan, no per-lambda multistart). The outer marginal is evaluated with
#' FOCEI (normal-Laplace factorising path), the inner conditional mode with
#' [trustL1]. See `notes/laplace_nlme_theory.Rmd`.
#'
#' The objective is assembled like an [EM] objective, with [constraintL1]
#' in place of `constraintL2(Omega = )`:
#' `obj <- normL2(data, g*x*p, errmodel = e) + constraintL1(pen)`.
#'
#' @param obj Composed objective carrying the model pieces and the `penaltySpec`.
#' @param init Named numeric start (structural + error + `lambda`); see [emInit].
#' @param fixed Optional named numeric of parameters held fixed.
#' @param method Marginal-likelihood backend for the factorising
#'   `penalty`/`fused` path. `"focei"` (default) is the deterministic ECM on the
#'   normal-Laplace FOCE marginal; `"saem"` is a stochastic-approximation EM with
#'   a per-subject random-walk Metropolis E-step, provided as an independent
#'   cross-check.
#' @param control List of tuning options: inner per-subject trust-region
#'   (`rinit`, `rmax`, `iterlim`, or `inner` = a list of any [trustL1] argument);
#'   outer CM-1 trust-region (`cm1` = a list of any [trust] argument); ECM controls (`maxOuter`,
#'   `epsPar`, `epsOfvRel`); `lockstep` (default `TRUE`: solve the inner modes of
#'   all subjects in unison so each round's objectives go out as one batched
#'   call; `FALSE`, or `options(dMod.laplace.lockstep = FALSE)`, restores the
#'   per-subject loop, which it reproduces exactly); `warm` (opt-in inner-mode warm-starting); `correction`
#'   `secondOrderCorrection` (logical, default `FALSE`: the normal-Laplace FOCEI
#'   volume term in the CM-1 gradient; needs a chain built with `deriv2 = TRUE`,
#'   otherwise the plain FOCE gradient is used); and `saem` (a list of `nBurnin`, `nEM`, `nMcmc`, `stepsize`, `cm1`
#'   for `method = "saem"`).
#' @param verbose Logical; print progress.
#' @return An object of class `c("em", "list")` with the estimated
#'   parameters (`argument`), objective `value` (-2 log L), the estimated
#'   `lambda`, the per-subject conditional modes `etaModes`, the sparsity
#'   pattern, and the recovered model pieces.
#' @seealso [penaltyL1], [constraintL1], [emInit], [sparsify],
#'   [msEM], [EM]
#' @noRd
.fitLaplace <- function(obj, init, fixed = NULL, method = c("focei", "saem"),
                   control = list(), verbose = TRUE) {
  method <- match.arg(method)
  rec      <- .laplaceReconstruct(obj)
  penalty  <- rec$penalty
  lam_name <- penalty$lambdaName                     # >= 1 strength names

  miss <- setdiff(lam_name, names(init))
  if (length(miss))
    stop(".fitLaplace: `init` is missing the penalty parameter(s) '",
         paste(miss, collapse = ", "), "'. Assemble it with emInit().",
         call. = FALSE)

  ## One normL2 per data condition, each restricted to that condition's data.
  ## Each subject's marginal needs its own value; a single all-subject objfn
  ## returns the sum, which would pollute every marginal (and inflate sigma).
  ## Conditions WITH a random effect (in the
  ## penalty's subjectEtas) get the normal-Laplace marginal; conditions WITHOUT
  ## one (reference lines in the 2-line reference encoding) contribute their plain
  ## data likelihood to the structural fit (see .laplaceOuterEval).
  all_conditions <- names(rec$data)
  resObjList <- setNames(lapply(all_conditions, function(s) {
    data_s <- rec$data[s]
    class(data_s) <- class(rec$data)
    normL2(data_s, rec$prdfn, errmodel = rec$errfn)
  }), all_conditions)

  ## Clustered penalty: the coupled (complete-graph fusion) path. Dispatch to the
  ## joint MAP-EM (R/L1Clustering.R) before .laplaceFactorisingTargets, which handles
  ## only the factorising single/fused penalties. (The exact marginal-ML clustered
  ## FOCEI value/gradient is a later step; the MAP-EM recovers the grouping + a MAP
  ## lambda.)
  if (any(vapply(penalty$blocks, `[[`, "", "kind") == "clustered")) {
    free0 <- init[setdiff(names(init), names(fixed))]
    if (method == "saem")
      return(.clusterSaem(rec, resObjList, free0, fixed, lam_name, control,
                                verbose))
    return(.clusterMap(rec, resObjList, free0, fixed, lam_name, control,
                             verbose))
  }
  targets_by_subject <- .laplaceFactorisingTargets(penalty)

  free <- init[setdiff(names(init), names(fixed))]
  lam_free_names <- intersect(lam_name, names(free))   # strengths to estimate
  init_all       <- c(free, fixed)
  lambda_cur     <- setNames(as.numeric(init_all[lam_name]), lam_name)
  struct_names   <- setdiff(names(free), lam_name)
  if (length(struct_names) == 0L)
    stop(".fitLaplace: no free structural parameters to optimise.", call. = FALSE)

  ## SAEM cross-check backend: stochastic E-step + closed-form lambda M-step.
  if (method == "saem")
    return(.laplaceSaem(rec, resObjList, targets_by_subject, free, fixed,
                       lam_name, control, verbose))

  ## per-subject inner-mode cache. Warm-starting speeds the inner trustL1 solve
  ## but can cache modes from rejected trust trial points; opt in with control$warm.
  useWarm <- if (is.null(control$warm)) FALSE else isTRUE(control$warm)
  warm    <- if (useWarm) new.env(parent = emptyenv()) else NULL

  ## Outer solver = ECM (theory eq 9 + the FOCE marginal). Per outer iteration:
  ##   CM-1  structural + error pars via dMod trust() on the smooth FOCE marginal
  ##         (analytic Fisher gradient + Gauss-Newton Hessian), lambda frozen.
  ##   CM-2  lambda in closed form, lambda* = K n / sum_ik E|eta_ik|, from the
  ##         normal-Laplace posterior moments at the new structural point.
  ## The inner conditional-mode solve (per-subject trustL1 + normal-Laplace) is the
  ## shared E-step. NB the outer trust() here is a *different* trust-region from
  ## the one inside trustL1: this one works on the smooth marginal over (mu, sigma)
  ## with no L1 kink, trustL1's works on the per-subject L1-penalised mode.
  cm1 <- .trustControl(list(rinit = 1, rmax = 10, iterlim = 30L,
                            ftol = 1e-6, mtol = 1e-6),
                       control$cm1, label = "control$cm1")
  maxOuter  <- if (is.null(control$maxOuter))  50L  else as.integer(control$maxOuter)
  epsPar    <- if (is.null(control$epsPar))    1e-4 else control$epsPar
  epsOfvRel <- if (is.null(control$epsOfvRel)) 1e-5 else control$epsOfvRel
  ## FOCEI normal-Laplace volume correction in the CM-1 gradient. Needs second
  ## derivatives of the prediction, so it is off by default and the plain FOCE
  ## gradient is what a first-order chain gets. Same control name as the normal
  ## path.
  useCorrection <- isTRUE(control$secondOrderCorrection)
  ## cap the penalty strength (a floor on the scale rho = 1/lambda). A shared
  ## parameter drives its deviations to 0, so its M-step lambda = n/sum E|eta|
  ## runs to infinity; the cap lets that strength settle (== "shared") so the
  ## ECM converges. Well above any strength a genuinely individual parameter takes.
  lambdaMax <- if (is.null(control$lambdaMax)) 1e4 else control$lambdaMax
  ## per-strength deviation counts for the closed-form M-step: strength j is
  ## shared by (#subjects * #params using j) deviations.
  nSub       <- nrow(penalty$subjectEtas)
  parsPerLam <- table(penalty$lambdaByParam)
  countByLam <- setNames(nSub * as.integer(parsPerLam[lam_name]), lam_name)
  lamByCoord <- unname(penalty$lambdaByParam[penalty$short])  # aligned to eta order

  mkTheta <- function(pstruct) {
    th <- c(setNames(as.numeric(pstruct), struct_names), fixed)
    th[lam_name] <- lambda_cur[lam_name]
    th
  }
  ## CM-1 objective: FOCE marginal over the free structural pars, lambda frozen.
  emReg <- function(pstruct, ...) {
    r <- .laplaceOuterEval(mkTheta(pstruct), rec, resObjList, targets_by_subject,
                       control, warm = warm, grad = TRUE, hess = TRUE,
                       correction = useCorrection)
    list(value    = r$value,
         gradient = r$gradient[struct_names],
         hessian  = r$hessian[struct_names, struct_names, drop = FALSE])
  }

  ## a strength whose deviations collapse (parameter shared) is frozen at the
  ## cap: its marginal is monotone in lambda (boundary), so left free it would
  ## climb forever and the ECM would never settle. epsShare = the mean|eta|
  ## below which a parameter counts as shared.
  epsShare   <- if (is.null(control$epsShare)) 0.02 else control$epsShare
  shared_lams <- character(0)
  psi_struct <- free[struct_names]
  prev_par   <- c(psi_struct, lambda_cur)
  prev_ofv   <- NA_real_
  rows <- list(); conv <- FALSE; nOuter <- 0L
  for (it in seq_len(maxOuter)) {
    nOuter <- it
    ## CM-1: structural block via the outer trust(), lambda held.
    cm1_fit <- suppressMessages(do.call(trust, modifyList(cm1, list(
      objfun = emReg, parinit = psi_struct))))
    psi_struct <- cm1_fit$argument

    ## moments at the new structural point (lambda still frozen): drive CM-2 + OFV.
    coll <- .laplaceOuterEval(mkTheta(psi_struct), rec, resObjList, targets_by_subject,
                          control, collect = TRUE, warm = warm)
    ofv <- coll$value
    ## CM-2: closed-form per-strength lambda M-step. Strength j:
    ## lambda_j = (#deviations using j) / sum_{k using j} E|eta_k|. Only the free
    ## strengths are updated.
    if (length(lam_free_names)) {
      EabsMat <- matrix(vapply(coll$perSubject, function(s) s$nl$Eabs,
                               numeric(penalty$K)), nrow = penalty$K)  # K x nSub
      sumEabsByLam <- vapply(lam_name,
        function(j) sum(EabsMat[lamByCoord == j, , drop = FALSE]), 0.0)
      meanAbsByLam <- sumEabsByLam / countByLam
      ## freeze collapsed (shared) strengths at the cap; update the rest.
      newly <- setdiff(lam_free_names[meanAbsByLam[lam_free_names] < epsShare],
                       shared_lams)
      shared_lams <- c(shared_lams, newly)
      lambda_cur[shared_lams] <- lambdaMax
      active <- setdiff(lam_free_names, shared_lams)
      if (length(active))
        lambda_cur[active] <- pmin((countByLam / sumEabsByLam)[active], lambdaMax)
    }

    cur_par     <- c(psi_struct, lambda_cur)
    deltaPar    <- max(abs(cur_par - prev_par))
    deltaOfvRel <- if (is.na(prev_ofv) || abs(prev_ofv) < .Machine$double.eps)
                     Inf else abs(ofv - prev_ofv) / abs(prev_ofv)
    rows[[length(rows) + 1L]] <- data.frame(
      outerIter = it, OFV = ofv, lambda = lambda_cur[[1]],
      deltaPar = deltaPar, deltaOfvRel = deltaOfvRel,
      cm1TrustIter = cm1_fit$iterations)
    if (verbose)
      cat(sprintf(
        "  EM(%s) outer %d : -2logL=%.4f  lambda=%s  |dpar|=%.2e\n",
        method, it, ofv,
        paste(formatC(lambda_cur, format = "g", digits = 4), collapse = "/"),
        deltaPar))
    prev_par <- cur_par; prev_ofv <- ofv
    if (deltaPar < epsPar || deltaOfvRel < epsOfvRel) { conv <- TRUE; break }
  }

  theta_hat <- mkTheta(psi_struct)
  final <- .laplaceOuterEval(theta_hat, rec, resObjList, targets_by_subject,
                         control, collect = TRUE, warm = warm)

  ## readout: per-subject modes, sparsity pattern, lambda-update check
  etaModes <- do.call(rbind, lapply(final$perSubject, function(s) unname(s$etahat)))
  rownames(etaModes) <- final$subjects
  colnames(etaModes) <- penalty$eta
  Eabs <- vapply(final$perSubject, function(s) s$nl$Eabs, numeric(penalty$K))
  Eabs <- matrix(Eabs, nrow = penalty$K)               # K x nSubjects
  ## closed-form per-strength MLE check: lambda_j* = (#deviations j)/sum_j E|eta|
  lambda_mle <- setNames(
    countByLam / vapply(lam_name,
      function(j) sum(Eabs[lamByCoord == j, , drop = FALSE]), 0.0), lam_name)

  out <- list(
    argument   = theta_hat,
    value      = final$value,
    lambda     = theta_hat[lam_name],
    lambda_mle = lambda_mle,
    etaModes   = etaModes,
    sparsity   = abs(etaModes) < 1e-6,
    converged  = conv,
    iterations = nOuter,
    stageTrace = do.call(rbind, rows),
    method     = method,
    prior      = "penaltyL1",
    penalty    = penalty,
    prdfn      = rec$prdfn,
    data       = rec$data,
    errfn      = rec$errfn
  )
  class(out) <- c("em", "list")
  out
}


## Penalty (L1/Laplace random-effect) branch of print.em.
.print_EM_penalty <- function(x, ...) {
  cat("EM (prior = penaltyL1, ", x$method, ")\n", sep = "")
  if (isTRUE(x$clustered)) {
    ## clustered readout
    lbl <- if (identical(x$method, "saem")) "SAEM iterations" else "MAP-EM iterations"
    cat(sprintf("  objective  : %.4f  (%s)\n", x$value, x$valueType))
    cat(sprintf("  lambda     : %.5g\n", x$lambda))
    cat(sprintf("  converged  : %s (%d %s)\n", x$converged, x$iterations, lbl))
    struct <- x$argument[setdiff(names(x$argument), x$penalty$lambdaName)]
    cat("  structural :\n"); print(round(struct, 4))
    if (!is.null(x$clusters)) {
      cat("  clustering (per parameter):\n")
      for (nm in names(x$clusters)) {
        grp <- vapply(x$clusters[[nm]], function(g) paste(g, collapse = "="), "")
        cat(sprintf("    %-14s %s\n", sub("^eta_", "", nm), paste(grp, collapse = " | ")))
      }
    } else {
      cat("  posterior mean deviations (E[eta] reveals the grouping;\n",
          "    the hard grouping comes from sparsify()):\n", sep = "")
      print(round(x$etaModes, 4))
    }
    return(invisible(x))
  }
  cat(sprintf("  -2 log L   : %.4f\n", x$value))
  if (length(x$lambda) == 1L) {
    cat(sprintf("  lambda     : %.5g  (closed-form MLE check: %.5g)\n",
                x$lambda, x$lambda_mle))
  } else {
    cat("  lambda (per parameter):\n")
    lamtab <- rbind(estimate = x$lambda, `MLE check` = x$lambda_mle[names(x$lambda)])
    print(signif(lamtab, 4))
  }
  cat(sprintf("  converged  : %s (%d ECM iterations)\n",
              x$converged, x$iterations))
  struct <- x$argument[setdiff(names(x$argument), x$penalty$lambdaName)]
  cat("  structural :\n")
  print(round(struct, 4))
  cat(sprintf("  shared (eta_hat == 0) fraction: %.2f\n", mean(x$sparsity)))
  if (!is.null(x$msTable))
    cat(sprintf("  multistart : %d fits, best of %d converged\n",
                nrow(x$msTable), sum(x$msTable$converged)))
  invisible(x)
}


#' Multi-start regularised NLME fit
#'
#' @description
#' The multi-start pendant of [EM], mirroring [msEM]. Runs [EM] from
#' several perturbed starting points and returns the best (lowest `-2 log L`) fit,
#' with the full ranked table attached as `$msTable`. Only the free *structural*
#' start is perturbed; the penalty strength `lambda` is re-estimated from its
#' closed-form M-step inside every fit, so its starting value does not matter.
#' Marginal-likelihood surfaces of the sharp Laplace kind can carry local optima,
#' so a multi-start is the safe way to compare models (supports) on an equal
#' footing (see `notes/laplace_nlme_theory.Rmd`, the support-selection agreement).
#'
#' @param obj,fixed,method,control As in [EM].
#' @param center Named numeric start (see [emInit]); the multi-start centre.
#' @param fits Number of starts (default 20).
#' @param cores Number of parallel workers (forked via `mclapply`).
#' @param sd Standard deviation of the perturbation applied to the free
#'   structural parameters (Gaussian by default).
#' @param samplefun Name of the perturbation sampler (default `"rnorm"`).
#' @param start1stfromCenter If `TRUE` (default) the first fit starts exactly at
#'   `center`; the rest are perturbed.
#' @param verbose Logical; print a one-line summary.
#' @return The best [EM] object, additionally carrying `$msTable` (fits ranked
#'   by value) and `$nfits`; class `c("msem", "em", "list")`.
#' @seealso [EM], [emInit], [msEM]
#' @noRd
.msfitLaplace <- function(obj, center, fixed = NULL, method = c("focei", "saem"),
                     control = list(), fits = 20L, cores = 1L,
                     sd = 0.5, samplefun = "rnorm",
                     start1stfromCenter = TRUE, verbose = FALSE) {
  method <- match.arg(method)
  cores  <- .sanitizeCores(cores)
  if (is.null(names(center)) || any(!nzchar(names(center))))
    stop("`center` must be a fully named numeric vector (see .initLaplace).")

  rec      <- .laplaceReconstruct(obj)
  lam_name <- rec$penalty$lambdaName
  ## perturb the free structural start only; lambda is re-estimated each fit.
  free_struct <- setdiff(names(center), c(names(fixed), lam_name))
  fits <- as.integer(fits)

  parInitList <- lapply(seq_len(fits), function(i) {
    ini <- center
    if (!(i == 1L && start1stfromCenter) && length(free_struct)) {
      pert <- do.call(samplefun, list(n = length(free_struct), sd = sd))
      ini[free_struct] <- ini[free_struct] + pert
    }
    ini
  })
  doOne <- function(i) {
    f <- try(suppressMessages(
      .fitLaplace(obj, parInitList[[i]], fixed = fixed, method = method,
             control = control, verbose = FALSE)), silent = TRUE)
    if (inherits(f, "try-error")) NULL else f
  }
  .cc <- .splitCores(cores, "fits")
  coresConditions <- .cc$conditions; cores <- .cc$outer
  results <- if (cores > 1L || !is.null(coresConditions))
    .parallelLapply(seq_len(fits), doOne, cores = cores,
                    coresConditions = coresConditions)
  else lapply(seq_len(fits), doOne)

  ok <- Filter(Negate(is.null), results)
  if (!length(ok)) stop(".msfitLaplace: all fits failed.")
  values <- vapply(ok, function(f) f$value, 0.0)
  ord    <- order(values)
  best   <- ok[[ord[1]]]
  tab <- data.frame(
    fit       = seq_along(ok),
    value     = values,
    lambda    = vapply(ok, function(f) f$lambda[[1]], 0.0),
    converged = vapply(ok, function(f) isTRUE(f$converged), NA))[ord, ]
  rownames(tab) <- NULL
  best$msTable <- tab
  best$nfits   <- length(ok)
  class(best)  <- c("msem", class(best))
  if (verbose)
    cat(sprintf(".msfitLaplace: %d/%d fits ok; best -2logL = %.4f (lambda = %s)\n",
                length(ok), fits, best$value,
                paste(formatC(best$lambda, format = "g", digits = 4),
                      collapse = "/")))
  best
}


#' Forward nested-support model selection by marginal likelihood
#'
#' @description
#' Packages the marginal-likelihood workflow of [EM] into a cheap model
#' selector. It replaces the classical penalty-strength scan (a multistart
#' [trustL1] fit at every point of a lambda grid, then a post-hoc information
#' criterion) with **one** mixed-effects fit plus an O(K) nested-support
#' comparison:
#'
#' \enumerate{
#'   \item Fit the full candidate model once with [msEM]; this estimates the
#'     single penalty strength `lambda` and ranks the candidate parameters by the
#'     magnitude of their per-subject deviations \eqn{|\hat d|}.
#'   \item Walk the forward-nested chain of supports
#'     \eqn{S_0 \subset S_1 \subset \dots \subset S_K} induced by that ranking
#'     (\eqn{S_0} = every candidate shared, \eqn{S_K} = the full candidate set),
#'     fitting each restricted model by marginal maximum likelihood (the excluded
#'     deviations are held at zero).
#'   \item Select the support with the smallest marginal \eqn{-2\log L}. Because
#'     the marginal integrates the deviations out, a spurious individual-specific
#'     parameter lowers the data \eqn{-2\log L} but *raises* the marginal one (the
#'     Occam factor), so the minimum is the parsimonious model.
#' }
#'
#' Asymptotically in the number of subjects the selected support equals the one a
#' [trustL1] lambda-scan with BIC would choose (see
#' `notes/laplace_nlme_theory.Rmd`), at O(K) fits instead of
#' O(grid times multistart).
#'
#' @param obj Composed objective for the **full** candidate set, i.e.
#'   `normL2(data, g*x*p, errmodel = e) + constraintL1(penalty(candidates,
#'   subjects = ...))`. Its prediction function must carry the deviation
#'   parameters `eta_<par>_<subject>` for every candidate; restricted supports are
#'   obtained by holding the excluded deviations at zero.
#' @param center Named numeric start (structural + error + `lambda`); see
#'   [emInit]. Forward selection multistarts the structural fit, so this is the
#'   perturbation centre (contrast the fixed `init` of [sparsify]).
#' @param fixed Optional named numeric held fixed across all fits (typically the
#'   non-candidate structural parameters).
#' @param method,control As in [EM].
#' @param fits,cores,sd Multistart settings forwarded to [msEM] for every
#'   support fit.
#' @param verbose Logical; print the ranking and the support chain as it is built.
#' @return An object of class `c("sparsify", "list")` with `type = "support"`,
#'   the selected `support` (character vector of individual-specific parameters),
#'   the `ranking` of the candidates by \eqn{|\hat d|}, the `chain` (a data frame
#'   with one row per nested support: `support`, `size`, `lambda`, `m2ll_marg`,
#'   `selected`), the winning fit `fitSelected`, and the full-support fit `fitFull`.
#' @seealso [sparsify], [EM]
#' @noRd
.selectLaplace <- function(obj, center, fixed = NULL, method = c("focei", "saem"),
                      control = list(), fits = 12L, cores = 1L, sd = 0.5,
                      verbose = TRUE) {
  method   <- match.arg(method)
  rec      <- .laplaceReconstruct(obj)
  penFull  <- rec$penalty
  lam_name <- penFull$lambdaName
  if (length(lam_name) > 1L)
    stop("sparsify requires a single global `lambda` (support enumeration). ",
         "Per-parameter strengths give soft sparsity directly in one .fitLaplace ",
         "(see penaltyL1(..., perParam = TRUE)); do not combine the two.",
         call. = FALSE)
  cands    <- penFull$params
  K        <- length(cands)
  subj     <- rownames(penFull$subjectEtas)
  if (is.null(names(center)) || any(!nzchar(names(center))))
    stop("`center` must be a fully named numeric vector (see .initLaplace).")
  if (!lam_name %in% names(center))
    stop("sparsify: `center` is missing the penalty parameter '", lam_name,
         "'. Assemble it with emInit().", call. = FALSE)

  struct_center <- center[setdiff(names(center), lam_name)]
  all_cand_etas <- as.vector(penFull$subjectEtas)   # eta_<par>_<subject>, all cands

  ## ---- 1. rank the candidates by |d_hat| from one full-support fit ----------
  if (verbose) cat("sparsify: full-support ranking fit (.msfitLaplace) ...\n")
  fitFull <- .msfitLaplace(obj, center, fixed = fixed, method = method,
                      control = control, fits = fits, cores = cores, sd = sd,
                      verbose = FALSE)
  ## per-parameter score = sum over subjects of |mode|; columns are in cands order
  score <- colSums(abs(fitFull$etaModes))
  ord   <- order(score, decreasing = TRUE)
  ranked <- cands[ord]
  if (verbose) {
    cat("  ranking (sum_i |d_hat|):\n")
    print(round(stats::setNames(score[ord], ranked), 4))
  }

  ## ---- 2/3. forward-nested supports, marginal -2logL, pick the minimum ------
  obj_data <- normL2(rec$data, rec$prdfn, errmodel = rec$errfn)

  evalSupport <- function(j) {
    S <- ranked[seq_len(j)]                          # j = 0 -> character(0)
    if (j == 0L) {
      ## null model: every candidate deviation held at zero, no penalty. Marginal
      ## -2logL of a no-random-effect model is just its data -2logL at the MLE.
      fixed0 <- c(fixed, stats::setNames(rep(0, length(all_cand_etas)),
                                         all_cand_etas))
      best <- NULL
      for (r in seq_len(fits)) {
        p0 <- struct_center
        if (r > 1L)
          p0 <- p0 + do.call("rnorm", list(n = length(p0), sd = sd))
        f <- try(suppressWarnings(trust(
          obj_data, parinit = p0, fixed = fixed0, rinit = 1, rmax = 10,
          iterlim = 200L, ftol = 1e-8, mtol = 1e-8)), silent = TRUE)
        if (!inherits(f, "try-error") && isTRUE(f$converged) &&
            (is.null(best) || f$value < best$value)) best <- f
      }
      if (is.null(best))
        return(list(fit = NULL, value = NA_real_, lambda = NA_real_))
      return(list(fit = best, value = best$value, lambda = NA_real_))
    }
    if (j == K)                                      # full support == fitFull
      return(list(fit = fitFull, value = fitFull$value, lambda = fitFull$lambda))
    excl_cols <- paste0(penFull$prefix, "_", setdiff(cands, S))
    etas_out  <- as.vector(penFull$subjectEtas[, excl_cols, drop = FALSE])
    pen_j <- penaltyL1(S, subjects = subj)
    obj_j <- normL2(rec$data, rec$prdfn, errmodel = rec$errfn) + constraintL1(pen_j)
    fixed_j <- c(fixed, stats::setNames(rep(0, length(etas_out)), etas_out))
    init_j  <- .initLaplace(struct_center, pen_j, lambda = 1)
    f <- try(suppressMessages(.msfitLaplace(
      obj_j, init_j, fixed = fixed_j, method = method, control = control,
      fits = fits, cores = cores, sd = sd, verbose = FALSE)), silent = TRUE)
    if (inherits(f, "try-error"))
      return(list(fit = NULL, value = NA_real_, lambda = NA_real_))
    list(fit = f, value = f$value, lambda = f$lambda)
  }

  results <- vector("list", K + 1L)
  rows    <- vector("list", K + 1L)
  for (j in 0:K) {
    res <- evalSupport(j)
    results[[j + 1L]] <- res
    S <- ranked[seq_len(j)]
    rows[[j + 1L]] <- data.frame(
      support   = if (j == 0L) "(none)" else paste(S, collapse = ","),
      size      = j,
      lambda    = res$lambda,
      m2ll_marg = res$value,
      stringsAsFactors = FALSE)
    if (verbose)
      cat(sprintf("  support[%d] {%s}: -2logL_marg=%s  lambda=%s\n",
                  j, if (j == 0L) "" else paste(S, collapse = ","),
                  format(res$value, digits = 6),
                  format(res$lambda, digits = 4)))
  }
  chain <- do.call(rbind, rows)

  if (all(is.na(chain$m2ll_marg)))
    stop("sparsify: every support fit failed to produce a marginal value.",
         call. = FALSE)
  jstar <- which.min(chain$m2ll_marg)                # which.min() skips NA
  Sstar <- ranked[seq_len(chain$size[jstar])]
  chain$selected <- seq_len(nrow(chain)) == jstar

  out <- list(
    type        = "support",
    support     = Sstar,
    ranking     = stats::setNames(score[ord], ranked),
    chain       = chain,
    fitSelected = results[[jstar]]$fit,
    fitFull     = fitFull,
    penalty     = penFull,
    method      = method)
  class(out) <- c("sparsify", "list")
  if (verbose)
    cat(sprintf("sparsify: S_hat = {%s}  (marginal -2logL = %.4f)\n",
                paste(Sstar, collapse = ","), chain$m2ll_marg[jstar]))
  out
}

## Support-selection branch of print.sparsify.
.print_sparsify_support <- function(x, ...) {
  cat("sparsify (support; ", x$method, ")\n", sep = "")
  cat("  candidate ranking (sum_i |d_hat|):\n")
  print(round(x$ranking, 4))
  cat("  nested support chain:\n")
  tab <- x$chain
  tab$m2ll_marg <- round(tab$m2ll_marg, 3)
  tab$lambda    <- round(tab$lambda, 4)
  tab$selected  <- ifelse(tab$selected, "<= S_hat", "")
  print(tab, row.names = FALSE)
  cat(sprintf("  S_hat = {%s}\n", paste(x$support, collapse = ", ")))
  invisible(x)
}


## ---- normal-Laplace marginal primitives ------------------------------------
## FOCEI marginal for the .fitLaplace factorising path (penalty + fused).
##
## For one candidate coordinate the linearised per-subject marginal contribution
## is the convolution of a data Gaussian (precision `a`, mode `m`) with a
## Laplace(lambda) prior:
##   I(a, m, lambda) = \int (lambda/2) e^{-lambda|eta|} e^{-(a/2)(eta - m)^2} d eta,
## the normal-Laplace density (Reed 2006; theory eq 13-14 in
## notes/laplace_nlme_theory.Rmd). We work with log I for numerical stability
## and return, in one pass, the outer-gradient partials
##   d logI / d{m, a, lambda}
## and the posterior moments E[eta], E[|eta|], E[(eta - m)^2], which follow from
## the identities
##   E[|eta|]        = 1/lambda - d logI / d lambda
##   E[eta] - m      = (1/a)      d logI / d m
##   E[(eta - m)^2]  = -2         d logI / d a.
## E[|eta|] is exactly the quantity the closed-form lambda update (eq 9) needs.

## Whether the .fitLaplace / clustered marginal kernels use the C++ port (default)
## or the pure-R reference. Toggle with options(laplaceUseCpp = FALSE); the two
## paths agree to ~1e-10 (tests/testthat/test-.fitLaplace.R).
.laplaceUseCpp <- function() isTRUE(getOption("laplaceUseCpp", TRUE))

.TWO_OVER_SQRTPI <- 2 / sqrt(pi)

## log of the scaled complementary error function erfcx(x) = exp(x^2) erfc(x),
## stable for all real x (erfc underflows for large x, erfcx does not).
.logErfcx <- function(x) log(2) + x^2 + stats::pnorm(-x * sqrt(2), log.p = TRUE)

## log(exp(la) + exp(lb)), overflow-safe.
.logSumExp2 <- function(la, lb) pmax(la, lb) + log1p(exp(-abs(la - lb)))

## d/dx log erfcx(x) = 2x - (2/sqrt(pi)) / erfcx(x) = 2x - (2/sqrt(pi)) exp(-log erfcx).
.dLogErfcx <- function(x, logErfcx_x) 2 * x - .TWO_OVER_SQRTPI * exp(-logErfcx_x)

#' Normal-Laplace 1-D marginal, its log-derivatives and posterior moments
#'
#' @description
#' Evaluates the convolution of a data Gaussian (precision `a`, mode `m`) with a
#' Laplace(`lambda`) prior in a numerically stable, vectorised way (arguments
#' recycle to a common length). Internal helper for the FOCEI regularisation
#' backend.
#'
#' @param a Positive data precision(s).
#' @param m Data-alone mode(s).
#' @param lambda Positive Laplace rate(s) (the single penalty strength).
#' @return A list with `logI`, `dlogI_dm`, `dlogI_da`, `dlogI_dlambda`, and the
#'   posterior moments `Eeta`, `Eabs` (= `E|eta|`), `Ecen2` (= `E[(eta-m)^2]`), and
#'   `Ppos` (posterior P(eta > 0)); each a numeric of the recycled length.
#' @keywords internal
normalLaplace <- function(a, m, lambda) {
  n  <- max(length(a), length(m), length(lambda))
  a  <- rep_len(a, n); m <- rep_len(m, n); lambda <- rep_len(lambda, n)

  cc  <- sqrt(a / 2)
  u   <- (lambda / a - m) * cc          # eta > 0 branch
  v   <- (lambda / a + m) * cc          # eta < 0 branch
  lEu <- .logErfcx(u)
  lEv <- .logErfcx(v)
  lse <- .logSumExp2(lEu, lEv)

  base <- log(lambda) - log(2) + 0.5 * (log(pi) - log(2 * a)) - a * m^2 / 2
  logI <- base + lse

  w_pos <- exp(lEu - lse)               # P(eta > 0)
  w_neg <- exp(lEv - lse)               # P(eta < 0)

  dEu <- .dLogErfcx(u, lEu)
  dEv <- .dLogErfcx(v, lEv)

  ## partials of u, v w.r.t. m, lambda, a
  dc_da <- cc / (2 * a)
  du_dm <- -cc;                                dv_dm <- cc
  du_dl <- cc / a;                             dv_dl <- cc / a
  du_da <- (-lambda / a^2) * cc + (lambda / a - m) * dc_da
  dv_da <- (-lambda / a^2) * cc + (lambda / a + m) * dc_da

  ## base partials
  dB_dm <- -a * m
  dB_dl <- 1 / lambda
  dB_da <- -1 / (2 * a) - m^2 / 2

  ## chain through the log-sum-exp: weight each branch by its posterior prob
  dlogI_dm <- dB_dm + w_pos * dEu * du_dm + w_neg * dEv * dv_dm
  dlogI_dl <- dB_dl + w_pos * dEu * du_dl + w_neg * dEv * dv_dl
  dlogI_da <- dB_da + w_pos * dEu * du_da + w_neg * dEv * dv_da

  ## posterior moments from the derivative identities
  Eeta  <- m + dlogI_dm / a
  Eabs  <- 1 / lambda - dlogI_dl
  Ecen2 <- -2 * dlogI_da

  list(logI = logI, dlogI_dm = dlogI_dm, dlogI_da = dlogI_da,
       dlogI_dlambda = dlogI_dl, Eeta = Eeta, Eabs = Eabs, Ecen2 = Ecen2,
       Ppos = w_pos)
}


## Per-subject FOCEI-normal-Laplace marginal for the factorising path.
##
## `objfun(eta)` must return the subject's data objlist in dMod's -2 log L
## convention: value = -2 log p(y_i | eta), gradient = d value / d eta,
## hessian = d^2 value / d eta^2 (~ 2 * J^T Sigma^-1 J, Gauss-Newton). The
## penalised conditional mode is found by trustL1; the linearised marginal is
## then the product of per-coordinate normal-Laplace integrals (separable case,
## a_k = diag). Returns the per-subject -2 log L_i contribution plus the mode and
## the ingredients (a, m, moments) the outer correction and the lambda update
## consume.
##
## Convention (derived in notes/laplace_nlme_theory.Rmd, eq 11-14):
##   Hd = diag(hessian),  a_k = Hd_k / 2  (theory precision J^T Sigma^-1 J),
##   m_k = etahat_k - grad_k / Hd_k       (data-alone mode, separable),
##   -2 log L_i = value(etahat) - (1/2) sum grad_k^2 / Hd_k - 2 sum log I_k,
##   with I_k = normalLaplace(a_k, m_k - target_k, lambda).
.laplaceSubjectMarginal <- function(objfun, eta0, targets, lambda,
                                control = list(rinit = 1, rmax = 10,
                                               iterlim = 100L,
                                               ftol = 1e-8, mtol = 1e-8)) {
  nm <- names(eta0)
  fit <- do.call(trustL1, modifyList(control, list(
    objfun = objfun, parinit = eta0, mu = targets, lambda = lambda)))
  etahat <- fit$argument[nm]

  at <- attr(objfun, "full")
  full <- if (is.null(at)) NULL else at(etahat)
  di <- if (is.null(full)) objfun(etahat) else
    list(value = full$value, gradient = full$gradient[nm],
         hessian = full$hessian[nm, nm, drop = FALSE])
  .laplaceMarginalFrom(di, full, etahat, targets, lambda,
                       isTRUE(fit$converged))
}


## The marginal, given the objective at the mode. Split out so the lock-step
## driver below and the single solve above share one copy of the arithmetic.
.laplaceMarginalFrom <- function(di, full, etahat, targets, lambda, converged) {
  nm <- names(etahat)
  gd <- di$gradient[nm]
  Hd <- diag(di$hessian[nm, nm, drop = FALSE])

  a       <- 0.5 * Hd
  m       <- etahat - gd / Hd
  m_shift <- m - targets[nm]
  nl      <- if (.laplaceUseCpp()) normalLaplaceCpp(a, m_shift, lambda)
             else normalLaplace(a, m_shift, lambda)

  correction <- -0.5 * sum(gd^2 / Hd)
  m2ll       <- di$value + correction - 2 * sum(nl$logI)

  list(value = m2ll, etahat = etahat, a = a, m = m, m_shift = m_shift,
       grad = gd, Hd = Hd, nl = nl, full = full, converged = converged)
}


## Every subject's conditional mode, stepping in unison. Each subject runs the
## same reflective trustL1 as the single solve; the driver only collects the
## trial points of the subjects that have not stopped and evaluates them in one
## batched call, which is where the condition axis parallelises. Falls back to
## the per-subject loop for anything the lock-step driver does not carry.
.laplaceSubjectMarginals <- function(objfuns, eta0List, targetsList, lambdaList,
                                     control = list(), lockstep = TRUE,
                                     cores = getOption("dMod.cores", 1L)) {
  N <- length(objfuns)
  serial <- function() lapply(seq_len(N), function(i)
    .laplaceSubjectMarginal(objfuns[[i]], eta0List[[i]], targetsList[[i]],
                            lambdaList[[i]], control))

  K <- length(eta0List[[1L]])
  parsAt   <- lapply(objfuns, attr, "pars_at")
  resobjs  <- lapply(objfuns, attr, "resobj")
  etaNames <- lapply(objfuns, attr, "etaNames")
  if (!isTRUE(lockstep) || N < 2L || K < 1L ||
      any(vapply(eta0List, length, 0L) != K) ||
      any(vapply(parsAt, is.null, TRUE)) ||
      !identical(control$boundary %||% "reflective", "reflective") ||
      !isTRUE(control$minimize %||% TRUE) ||
      isTRUE(control$blather) || isTRUE(control$printIter) ||
      !is.null(control$traceFile) || !is.null(control$parscale) ||
      !is.null(control$parupper) || !is.null(control$parlower))
    return(serial())

  rows <- function(l, f) matrix(vapply(l, f, numeric(K)), N, K, byrow = TRUE)
  P0  <- rows(eta0List,    as.numeric)
  MU  <- rows(targetsList, as.numeric)
  LAM <- rows(lambdaList,  function(z) rep_len(as.numeric(z), K))
  colnames(P0) <- paste0("eta", seq_len(K))

  ## One batched pass, with the per-subject loop as the fallback so a batch that
  ## cannot run keeps each subject's own failure isolated.
  evalMany <- function(idx, fulls) {
    got <- tryCatch(.objEvalMany(resobjs[idx], fulls, deriv = TRUE, cores = cores),
                    error = function(e) NULL)
    if (is.null(got)) got <- lapply(seq_along(idx), function(j)
      tryCatch(resobjs[[idx[j]]](fulls[[j]], deriv = TRUE), error = function(e) NULL))
    got
  }
  objfun_many <- function(parsList, which) {
    got <- evalMany(which, lapply(seq_along(which), function(j)
      parsAt[[which[j]]](as.numeric(parsList[[j]]))))
    lapply(seq_along(which), function(j) {
      o <- got[[j]]
      if (is.null(o)) return(NULL)
      en <- etaNames[[which[j]]]
      list(value = o$value, gradient = o$gradient[en],
           hessian = o$hessian[en, en, drop = FALSE])
    })
  }

  dflt <- function(nm) eval(formals(trustL1)[[nm]])
  pick <- function(nm) control[[nm]] %||% dflt(nm)
  fits <- trustL1_lockstep_impl(
    objfun_many, P0, MU, LAM, isTRUE(control$one.sided),
    pick("rinit"), pick("rmax"), NULL, as.integer(pick("iterlim")),
    pick("ftol"), pick("mtol"), pick("gtol"), pick("xtol"), pick("rmin"),
    pick("theta.max"), TRUE, FALSE, NULL, NULL)

  etahats <- lapply(seq_len(N), function(i)
    setNames(as.numeric(fits[[i]]$argument), etaNames[[i]]))
  os <- evalMany(seq_len(N), lapply(seq_len(N), function(i)
    parsAt[[i]](as.numeric(fits[[i]]$argument))))
  lapply(seq_len(N), function(i) {
    nm <- etaNames[[i]]
    full <- os[[i]]
    if (is.null(full))
      return(.laplaceSubjectMarginal(objfuns[[i]], eta0List[[i]],
                                     targetsList[[i]], lambdaList[[i]], control))
    di <- list(value = full$value, gradient = full$gradient[nm],
               hessian = full$hessian[nm, nm, drop = FALSE])
    .laplaceMarginalFrom(di, full, etahats[[i]], targetsList[[i]],
                         lambdaList[[i]], isTRUE(fits[[i]]$converged))
  })
}
