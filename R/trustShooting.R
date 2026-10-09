## Trust-region multiple shooting ------------------------------------------------
## Bock's generalised Gauss-Newton method, condensed onto theta, on a trust region
## over parameters and nodes; acceptance by filter, l2 penalty or natural level.

## Copyright (C) 2026 Simon Beyer


# trust() on an objective holding a normL2 with a multipleShootingControl.
# Settings a shooting run cannot honour are refused. Every other summand enters
# the condensed problem with its own value, gradient and Gauss-Newton Hessian.
.trustShooting <- function(objfun, parinit, rinit, rmax, iterlim,
                           hessianMethod, hessianFallback, fallbackLimit,
                           parscale, parupper, parlower, tol, qn, step,
                           minimize, blather, printIter, traceFile, dots) {

  hm <- hessianMethod
  bad <- c(
    if (hm == "exact") "hessianMethod = \"exact\"",
    if (hessianFallback != "none") "a hessianFallback",
    if (!identical(qn$hessianInit, "gn") && !identical(qn$hessianInit, "identity"))
      "qnControl$hessianInit = \"exact\"",
    if (!identical(qn$hessianReseed %||% "never", "never")) "qnControl$hessianReseed",
    if (!identical(as.integer(qn$qnMemory), 0L)) "qnControl$qnMemory",
    if (!identical(step$boundary, "reflective")) "stepControl$boundary = \"clip\"",
    if (!identical(step$nonmonotone, 0) && !identical(step$nonmonotone, 0L)) "stepControl$nonmonotone",
    if (!isTRUE(minimize)) "minimize = FALSE",
    if (!is.null(traceFile)) "traceFile")
  if (length(bad))
    stop("trust: a multiple-shooting objective does not support ",
         paste(bad, collapse = ", "), ".", call. = FALSE)
  acc <- match.arg(step$acceptance, c("filter", "merit", "natural"))
  if (acc == "natural" && hm != "gn")
    stop("trust: stepControl$acceptance = \"natural\" is Bock's damped ",
         "Gauss-Newton method and needs hessianMethod = \"gn\".", call. = FALSE)

  # What normL2 would honour and a multiple-shooting run cannot; anything else
  # it ignores, as normL2 itself does.
  refused <- intersect(names(dots), c("deriv2", "hessian", "conditions", "env",
                                      ".prediction"))
  if (length(refused))
    stop("trust: a multiple-shooting objective takes fixed, cores and sweep, ",
         "not ", paste(refused, collapse = ", "), ".", call. = FALSE)
  fixed <- dots$fixed
  cores <- if (is.null(dots$cores)) getOption("dMod.cores", 1L) else dots$cores
  sweep <- match.arg(dots$sweep %||% "forward", c("forward", "reverse"))
  if (hm == "gn" && sweep == "reverse")
    stop("trust: a Gauss-Newton block needs forward sensitivities; with ",
         "sweep = \"reverse\" a multiple-shooting run needs ",
         "hessianMethod = \"bfgs\" or \"sr1\".", call. = FALSE)

  terms <- .objTerms(objfun)
  isS <- vapply(terms, function(t) !is.null(.shootControlOf(t)), TRUE)
  # laid out from the control as it is now, so controls<- takes effect
  sobj <- .shootObjOf(objfun)
  others <- terms[!isS]
  spec <- attr(sobj, "spec")

  tnames <- names(parinit)
  if (is.null(tnames) || !all(nzchar(tnames)))
    stop("trust: 'parinit' must be named.", call. = FALSE)
  theta <- setNames(as.numeric(parinit), tnames)
  if (!all(is.finite(theta))) stop("trust: 'parinit' must be finite.", call. = FALSE)
  K <- length(theta)
  lb <- .shootBound(parlower, tnames, -Inf)
  ub <- .shootBound(parupper, tnames, Inf)
  ps <- if (is.null(parscale)) rep(1, K) else as.numeric(parscale)
  if (length(ps) != K || !all(ps > 0))
    stop("trust: parscale must be positive, one per parameter.", call. = FALSE)
  if (any(theta < lb | theta > ub)) {
    warning("trust: 'parinit' outside the bounds was clipped.", call. = FALSE)
    theta <- pmin(pmax(theta, lb), ub)
  }
  theta <- .shootInterior(theta, lb, ub)

  wantH <- hm == "gn"
  adaptive <- isTRUE(attr(sobj, "adaptive"))
  # With deriv = FALSE a solve computes no sensitivities: what a trial point
  # needs to be judged, the data term and the gaps, at a fraction of the cost.
  ev <- function(th, nd, hessian = wantH, sw = sweep, quiet = FALSE, deriv = TRUE,
                 residuals = acc == "natural") {
    call <- function() sobj(th, nd, fixed = fixed, deriv = deriv,
                            hessian = hessian && deriv, sweep = sw, cores = cores,
                            growth = adaptive && deriv, residuals = residuals)
    r <- if (quiet)
      tryCatch(withCallingHandlers(call(), warning = function(w) invokeRestart("muffleWarning")),
               error = function(e) e)
    else tryCatch(call(), error = function(e) e)
    if (inherits(r, "error")) return(r)
    if (!is.finite(r$value)) return(simpleError("objective not finite"))
    if (nBreak) r$gaps <- zeroBreaks(r$gaps)
    fin <- vapply(r$segments, function(sg)
      all(is.finite(sg$grad)) && (is.null(sg$jac) || all(is.finite(sg$jac))) &&
        (is.null(sg$hess) || all(is.finite(sg$hess))), TRUE)
    if (!all(fin))
      return(simpleError(paste0("derivatives not finite in segment ",
                                which(!fin)[1L])))
    r
  }
  # The summands on theta alone, forward and with their Gauss-Newton Hessian:
  # a prior costs a line of algebra, another data term what it always costs.
  zeroObj <- list(value = 0, gradient = setNames(numeric(K), tnames),
                  hessian = matrix(0, K, K, dimnames = list(tnames, tnames)))
  pv <- function(th) {
    if (!length(others)) return(zeroObj)
    out <- zeroObj
    for (t in others) {
      o <- t(th, fixed = fixed, deriv = TRUE, hessian = TRUE, cores = cores)
      o <- .alignObjlist(o, tnames)
      out$value <- out$value + o$value
      out$gradient <- out$gradient + o$gradient
      if (!is.null(o$hessian)) out$hessian <- out$hessian + o$hessian
    }
    out
  }

  init <- attr(sobj, "init")
  initNodes <- function(th) {
    nd <- if (is.list(init)) init
          else .shootNodes(sobj, th, fixed = fixed,
                           observed = !identical(init, "simulation"), cores = cores,
                           smooth = identical(init, "spline"))
    # Unobserved states follow a run synchronised to the data.
    if (is.character(init) && init %in% c("data", "spline"))
      nd <- .shootSweep(sobj, th, nd, fixed, cores)
    nd
  }
  nodes <- initNodes(theta)
  scale <- .shootScale(attr(sobj, "scale"), nodes, spec, sobj, theta, fixed, cores)
  neval <- 0L

  # --- continuity breaks ------------------------------------------------------
  # Voss, Timmer and Kurths (2004): the node after a break joins the parameters
  # and its gap is no constraint. Extended variables: theta, then break nodes.
  breakTimes <- attr(sobj, "breaks")
  nBreak <- length(unlist(breakTimes))
  if (nBreak && hm != "gn")
    stop("trust: continuity breaks need hessianMethod = \"gn\".", call. = FALSE)
  bkOf <- function() {
    brk <- list(); bcol <- list(); k <- K
    nx <- length(spec$states)
    for (cn in spec$conditions) {
      b <- spec$tau[[cn]][-1L] %in% breakTimes[[cn]]
      brk[[cn]] <- b
      bcol[[cn]] <- list()
      for (i in which(b)) {
        bcol[[cn]][[as.character(i)]] <- k + seq_len(nx)
        k <- k + nx
      }
    }
    list(brk = brk, bcol = bcol)
  }
  bxNames <- if (nBreak) unlist(lapply(spec$conditions, function(cn)
    unlist(lapply(breakTimes[[cn]], function(tb) paste0(spec$states, "@break_", cn, "_", tb)))))
  tx <- c(tnames, bxNames)
  Kx <- length(tx)
  thx <- function(th) setNames(c(th, numeric(Kx - K)), tx)
  lbx <- c(lb, rep(-Inf, Kx - K))
  ubx <- c(ub, rep(Inf, Kx - K))
  psxOf <- function() c(ps, rep(1 / scale[spec$states], nBreak))
  zeroBreaks <- function(gaps) {
    if (!nBreak) return(gaps)
    b <- bkOf()$brk
    for (cn in names(gaps)) if (any(b[[cn]])) gaps[[cn]][b[[cn]], ] <- 0
    gaps
  }
  # the condensed problem in the extended variables, the other summands added
  condense <- function(E, PRg = PR$gradient) {
    bk <- bkOf()
    C <- .shootCondense(E, blocks(E), spec, tx, scale, brk = bk$brk, bcol = bk$bcol)
    C$g0 <- C$g0 + c(PRg, numeric(Kx - K))
    C$H[seq_len(K), seq_len(K)] <- C$H[seq_len(K), seq_len(K)] + PR$hessian
    C
  }

  # --- two phases -------------------------------------------------------------
  # Horbelt, Timmer and Voss (2002): theta is fitted first with the nodes held,
  # then the nodes are laid out afresh and released. A failing phase is dropped.
  start0 <- NULL
  if (isTRUE(step$twoPhase)) {
    ph <- .shootPhaseOne(sobj, theta, nodes, fixed, cores, pv, lb, ub, ps,
                         rinit, rmax, iterlim)
    neval <- neval + ph$neval
    if (!is.null(ph$theta)) {
      nd <- tryCatch(initNodes(ph$theta), error = function(e) NULL)
      E1 <- if (!is.null(nd)) ev(ph$theta, nd, quiet = TRUE)
      if (!is.null(nd) && !inherits(E1, "error")) {
        start0 <- list(theta = theta, nodes = nodes, scale = scale)
        theta <- ph$theta
        nodes <- nd
        scale <- .shootScale(attr(sobj, "scale"), nodes, spec, sobj, theta,
                             fixed, cores)
      }
      if (printIter)
        cat(sprintf("phase 1  f = %.8g  %d iterations%s\n", ph$value, ph$iterations,
                    if (is.null(nd) || inherits(E1, "error")) ", dropped" else ""))
    }
  }

  # --- adaptive nodes -------------------------------------------------------
  # A segment whose propagation matrix exceeds `growth` or whose solve fails is
  # cut; the new node takes its own state, or the start value and the data.
  kappa <- attr(sobj, "growth")
  misfit <- attr(sobj, "misfit")
  minLen <- attr(sobj, "minLength")
  minPts <- attr(sobj, "minPoints") %||% 1L
  nSplit <- 0L
  cuts <- c(growth = 0L, misfit = 0L, failure = 0L)
  eventsOf <- function(cn) {
    p0 <- .shootInner(spec, theta, fixed, cores)[[cn]]
    .shootEventTimes(spec, p0)
  }
  splittable <- function(k) {
    sg <- spec$segs[[k]]
    !is.na(.shootSplitTime(sg, minLen[[sg$cond]], minPoints = minPts))
  }
  # The spectral radius, not a norm: it is invariant under chart and state
  # scaling, and is the rate the linearised flow expands at.
  growthOf <- function(E) vapply(E$segments, function(sg) {
    G <- sg$growth
    if (is.null(G)) return(0)
    if (!all(is.finite(G))) return(Inf)
    max(Mod(eigen(G, only.values = TRUE)$values))
  }, 0)
  # segments to cut, and whether a cut anchors at the data (misfit) or keeps
  # the trajectory (growth). A misfit counts where it is local, above twice the
  # median of the condition: one everywhere is the parameters'.
  violators <- function(E) {
    rms <- vapply(E$segments, function(sg)
      if (is.null(sg$chi2) || !isTRUE(sg$ndata > 0)) NA_real_ else sqrt(sg$chi2 / sg$ndata), 0)
    med <- stats::ave(rms, vapply(spec$segs, `[[`, "", "cond"),
                      FUN = function(v) stats::median(v, na.rm = TRUE))
    byGrowth <- growthOf(E) > kappa
    byData <- !is.na(rms) & rms > misfit & rms > 2 * med
    ks <- which(byGrowth | byData)
    ok <- vapply(ks, splittable, TRUE)
    list(k = ks[ok], data = byData[ks[ok]] & !byGrowth[ks[ok]],
         why = ifelse(byGrowth[ks[ok]], "growth", "misfit"))
  }
  anchorData <- function(v, k, tn) {
    if (identical(attr(sobj, "init"), "simulation")) return(v)
    sg <- spec$segs[[k]]
    obs <- .shootObservedStates(spec)
    dc <- attr(sobj, "args")$data[[sg$cond]]
    for (st in names(obs)) {
      di <- dc[as.character(dc$name) == obs[[st]], , drop = FALSE]
      if (!nrow(di)) next
      val <- .shootObservedCurve(di, identical(attr(sobj, "init"), "spline"))(tn)
      if (spec$charts[[st]] == "log10") val <- log10(max(val, 1e-12))
      v[st] <- val
    }
    v
  }
  fallbackAt <- function(k, tn) {
    sg <- spec$segs[[k]]
    v <- if (sg$j > 1L) nodes[[sg$cond]][sg$j - 1L, spec$states]
         else {
           x0 <- unclass(.shootInner(spec, theta, fixed, cores)[[sg$cond]])[spec$states]
           .shootToChart(pmax(x0, ifelse(spec$charts == "log10", 1e-12, -Inf)),
                         spec$charts)$value
         }
    anchorData(v, k, tn)
  }
  splitSegments <- function(ks, fromTrajectory, atData = rep(FALSE, length(ks)),
                            why = rep("failure", length(ks))) {
    tau <- spec$tau
    fill <- setNames(lapply(spec$conditions, function(cn) list()), spec$conditions)
    added <- 0L
    byWhy <- cuts
    evts <- lapply(setNames(nm = spec$conditions), eventsOf)
    tns <- vapply(ks, function(k) {
      sg <- spec$segs[[k]]
      .shootSplitTime(sg, minLen[[sg$cond]], evts[[sg$cond]], minPts)
    }, 0)
    vals <- vector("list", length(ks))
    if (fromTrajectory && any(!is.na(tns)))
      vals[!is.na(tns)] <- .shootValuesAt(sobj, theta, nodes, fixed, ks[!is.na(tns)],
                                          tns[!is.na(tns)], cores)
    for (i in seq_along(ks)) {
      k <- ks[i]
      sg <- spec$segs[[k]]
      tn <- tns[i]
      if (is.na(tn)) next
      v <- vals[[i]]
      if (is.null(v)) v <- fallbackAt(k, tn)
      else if (atData[i]) v <- anchorData(v, k, tn)
      tau[[sg$cond]] <- sort(c(tau[[sg$cond]], tn))
      fill[[sg$cond]][[as.character(tn)]] <- v
      added <- added + 1L
      byWhy[why[i]] <- byWhy[why[i]] + 1L
    }
    if (!added) return(FALSE)
    oldKeys <- vapply(spec$segs, function(sg) paste(sg$cond, sg$start, sg$end), "")
    rl <- .shootRelayout(sobj, tau, nodes, fill)
    sobj <<- rl$sobj
    spec <<- attr(sobj, "spec")
    nodes <<- rl$nodes
    nSplit <<- nSplit + added
    cuts <<- byWhy
    lastKeys <<- oldKeys
    if (printIter)
      cat(sprintf("cut      %d segments (growth %d, misfit %d, failure %d)\n",
                  length(spec$segs), byWhy[["growth"]], byWhy[["misfit"]],
                  byWhy[["failure"]]))
    TRUE
  }
  lastKeys <- NULL
  # Cuts at an accepted iterate. A new layout that does not evaluate, a cut
  # that lands where the solver cannot start, is undone: the step stands on
  # the old layout. NULL when nothing was cut.
  splitAt <- function(vi) {
    keep <- list(sobj = sobj, spec = spec, nodes = nodes, nSplit = nSplit,
                 lastKeys = lastKeys, cuts = cuts)
    if (!length(vi$k) || !splitSegments(vi$k, TRUE, vi$data, vi$why)) return(NULL)
    E2 <- ev(theta, nodes, quiet = TRUE)
    neval <<- neval + 1L
    if (inherits(E2, "error")) {
      sobj <<- keep$sobj; spec <<- keep$spec; nodes <<- keep$nodes
      nSplit <<- keep$nSplit; lastKeys <<- keep$lastKeys; cuts <<- keep$cuts
      return(NULL)
    }
    E2
  }

  layoutNodes <- function() {
    E <- NULL
    for (round in seq_len(200L)) {
      E <- ev(theta, nodes, quiet = adaptive)
      neval <<- neval + 1L
      if (!adaptive) break
      if (inherits(E, "shootingSegmentError")) {
        if (!splitSegments(E$segment, FALSE)) break
        next
      }
      if (inherits(E, "error")) break
      vi <- violators(E)
      if (!length(vi$k) || !splitSegments(vi$k, TRUE, vi$data, vi$why)) break
    }
    E
  }
  sobj0 <- sobj
  spec0 <- spec
  E <- layoutNodes()
  # the parameters of a first phase can sit where a segment cut on the way
  # cannot be solved; the run then starts from where it was asked to
  if (inherits(E, "error") && !is.null(start0)) {
    theta <- start0$theta; nodes <- start0$nodes; scale <- start0$scale
    sobj <- sobj0; spec <- spec0; nSplit <- 0L; lastKeys <- NULL
    cuts[] <- 0L
    E <- layoutNodes()
  }
  if (inherits(E, "error"))
    stop("trust: parinit not feasible: ", conditionMessage(E), call. = FALSE)
  PR <- pv(theta)

  # --- annealing of the continuity, then the exact method -----------------
  nStages <- 0L
  if (isTRUE(step$anneal)) {
    evGN <- function(th, nd) ev(th, nd, hessian = TRUE, sw = "forward", quiet = TRUE)
    evV <- function(th, nd) ev(th, nd, quiet = TRUE, deriv = FALSE)
    an <- .shootAnneal(evGN, spec, theta, nodes, tnames, ps, scale, pv, lb, ub,
                       rinit, rmax, printIter, evV = evV, brk = bkOf()$brk)
    if (!is.null(an)) {
      theta <- an$theta; nodes <- an$nodes
      neval <- neval + an$neval
      nStages <- an$stages
      E <- ev(theta, nodes, quiet = TRUE)
      neval <- neval + 1L
      if (inherits(E, "error"))
        stop("trust: the annealed point failed to evaluate: ", conditionMessage(E),
             call. = FALSE)
      PR <- pv(theta)
    }
  }

  # Model Hessian per segment: its Gauss-Newton block, or a quasi-Newton one.
  # Kept by segment across a new layout; a segment that is new gets a seed.
  seedBlocks <- function(E, which) {
    if (qn$hessianInit == "gn") {
      Eg <- ev(theta, nodes, hessian = TRUE, sw = "forward")
      if (inherits(Eg, "error"))
        stop("trust: the Gauss-Newton seed failed: ", conditionMessage(Eg), call. = FALSE)
      neval <<- neval + 1L
      lapply(which, function(k) .shootAlignMat(Eg$segments[[k]]$hess, E$segments[[k]]$vars))
    } else lapply(which, function(k) {
      v <- E$segments[[k]]$vars
      Mc <- spec$segs[[k]]$M
      matrix(diag(ifelse(v %in% tnames, 1 / Mc, 1), length(v)), length(v),
             dimnames = list(v, v))
    })
  }
  segKeys <- function() vapply(spec$segs, function(sg) paste(sg$cond, sg$start, sg$end), "")
  B <- NULL
  Bkeys <- NULL
  if (hm != "gn") {
    B <- seedBlocks(E, seq_along(spec$segs))
    Bkeys <- segKeys()
  }
  remapBlocks <- function(E) {
    if (hm == "gn") return(invisible(NULL))
    keys <- segKeys()
    hit <- match(keys, Bkeys)
    Bn <- vector("list", length(keys))
    Bn[!is.na(hit)] <- B[hit[!is.na(hit)]]
    fresh <- which(is.na(hit))
    if (length(fresh)) Bn[fresh] <- seedBlocks(E, fresh)
    # the same segment under its new index: node names contain the index, the
    # order of the variables does not change
    for (k in which(!is.na(hit))) {
      v <- E$segments[[k]]$vars
      if (nrow(Bn[[k]]) == length(v)) dimnames(Bn[[k]]) <- list(v, v)
      else Bn[k] <- seedBlocks(E, k)
    }
    B <<- Bn
    Bkeys <<- keys
  }
  blocks <- function(E) if (hm == "gn") lapply(E$segments, `[[`, "hess") else B

  f <- E$value + PR$value
  h <- .shootViolation(E$gaps, scale)
  r <- rinit
  mu <- 0
  filt <- matrix(numeric(0), 0, 2)
  hmax <- max(1.25 * h, 1)
  n_fail <- 0L; n_stall <- 0L; nRelaxed <- 0L; nSOC <- 0L; qnSkipped <- 0L
  nRejRow <- 0L; nRestore <- 0L

  # Filter restoration: close the linearised gaps at theta, damped if needed;
  # failing that, put every node on theta's trajectory. NULL when neither works.
  restoreGaps <- function(C) {
    for (a in c(1, 0.5, 0.25)) {
      ndR <- .shootStepNodes(nodes, C, a, numeric(Kx))
      ER <- ev(theta, ndR, quiet = TRUE)
      neval <<- neval + 1L
      if (!inherits(ER, "error") && .shootViolation(ER$gaps, scale) < 0.5 * h)
        return(list(E = ER, nodes = ndR))
    }
    if (nBreak) return(NULL)
    ndS <- tryCatch(.shootNodes(sobj, theta, fixed = fixed, observed = FALSE,
                                cores = cores), error = function(e) NULL)
    if (is.null(ndS)) return(NULL)
    ER <- ev(theta, ndS, quiet = TRUE)
    neval <<- neval + 1L
    if (inherits(ER, "error") || !(.shootViolation(ER$gaps, scale) < h)) return(NULL)
    list(E = ER, nodes = ndS)
  }
  converged <- FALSE; stopReason <- "iterlim"
  iter <- 0L
  trace <- list()
  C <- NULL

  # --- Bock's damped Gauss-Newton method on the natural level function -------
  # Peifer and Timmer (2007), sec. 3.4: predictor-corrector damping; a step is
  # taken when the simplified step at its end has shrunk.
  handedOver <- FALSE
  if (acc == "natural") {
    if (!is.null(attr(sobj, "args")$errmodel) ||
        !all(vapply(spec$segs, function(sg) all(is.finite(sg$data$sigma)), TRUE)))
      stop("trust: stepControl$acceptance = \"natural\" needs a fixed sigma in the ",
           "data and no error model.", call. = FALSE)
    # Peifer and Timmer's constants; the corrector may damp below tmin down
    # to lmin, Deuflhard's setting for highly nonlinear problems
    tmin <- 0.01; tup <- 0.5; eta0 <- 1; eta2 <- 1.8; lmin <- 1e-6
    reg <- step$regularise %||% 0
    E <- ev(theta, nodes, quiet = TRUE)
    neval <- neval + 1L
    if (inherits(E, "error"))
      stop("trust: parinit not feasible: ", conditionMessage(E), call. = FALSE)
    f <- E$value + PR$value
    h <- .shootViolation(E$gaps, scale)
    omega <- NA_real_
    for (iter in seq_len(iterlim)) {
      C <- condense(E)
      psx <- psxOf()
      opt <- .shootOptimality(thx(theta), C$g0, lbx, ubx, psx)
      if (opt$measure <= tol$gtol && h <= tol$ctol) {
        converged <- TRUE; stopReason <- "gradient"; break
      }
      FS <- .shootFullStep(C, thx(theta), lbx, ubx, psx, reg)
      nDu <- FS$norm
      if (!(nDu > 1e-12)) {
        converged <- h <= tol$ctol; stopReason <- "step"; break
      }
      mu0 <- if (is.na(omega)) tmin else eta0 / (omega * nDu)
      lam <- if (mu0 > tup) 1 else max(mu0, tmin)
      taken <- FALSE
      for (j in seq_len(30L)) {
        thetaT <- theta + lam * FS$dth[seq_len(K)]
        nodesT <- .shootStepNodes(nodes, C, lam, lam * FS$dth)
        EV <- ev(thetaT, nodesT, quiet = TRUE, deriv = FALSE)
        neval <- neval + 1L
        if (inherits(EV, "error")) {
          if (lam <= lmin) break
          lam <- max(lam / 4, lmin)
          next
        }
        PRT <- pv(thetaT)
        # the simplified step at the trial: Jacobians of the iterate,
        # residuals and gaps of the trial
        Eb <- E
        Eb$gaps <- EV$gaps
        for (k in seq_along(Eb$segments))
          Eb$segments[[k]]$grad <- setNames(
            2 * as.numeric(crossprod(E$segments[[k]]$J, EV$segments[[k]]$res)),
            E$segments[[k]]$vars)
        Cb <- condense(Eb, PRT$gradient)
        vb <- .shootFullStep(C, thx(theta), lbx, ubx, psx, reg, g = Cb$g0 + Cb$gz,
                             Vz = Cb$Vz, eig = FS$eig, truncate = FALSE)$vec
        nb <- sqrt(sum(vb^2))
        omegaT <- 2 * sqrt(sum((vb - (1 - lam) * FS$vec)^2)) / (lam * nDu)^2
        descent <- omegaT * lam * nDu <= eta2 && nb < nDu
        if (printIter)
          cat(sprintf("natural %4d  lambda = %.3g  |du| = %.4g  |du bar| = %.4g  omega = %.3g  %s\n",
                      iter, lam, nDu, nb, omegaT, if (descent) "accept" else "reject"))
        if (descent) { taken <- TRUE; omega <- omegaT; break }
        if (lam <= lmin) break
        lam <- max(min(eta0 / (omegaT * nDu), lam / 2), lmin)
      }
      if (blather)
        trace[[length(trace) + 1L]] <- list(theta = theta, value = f, violation = h,
          lambda = lam, omega = omega, stepnorm = nDu, accept = taken,
          segments = length(spec$segs))
      if (!taken) { stopReason <- "damping"; break }
      ET <- ev(thetaT, nodesT, quiet = TRUE)
      neval <- neval + 1L
      if (inherits(ET, "error")) { stopReason <- "objfun"; break }
      df <- abs(f - (ET$value + PRT$value))
      theta <- thetaT; nodes <- nodesT; E <- ET; PR <- PRT
      f <- E$value + PR$value
      h <- .shootViolation(E$gaps, scale)
      if (adaptive) {
        E2 <- splitAt(violators(E))
        if (!is.null(E2)) {
          E <- E2
          f <- E$value + PR$value
          h <- .shootViolation(E$gaps, scale)
          omega <- NA_real_
        }
      }
      if (lam == 1 && h <= tol$ctol && df < tol$ftol) {
        converged <- TRUE; stopReason <- "fvalue"; break
      }
    }
    C <- NULL
    # Where the damping finds no step, the linearisation holds nowhere along
    # the full step; the filter's trust region takes over from here.
    if (stopReason == "damping") {
      handedOver <- TRUE
      acc <- "filter"
      stopReason <- "iterlim"
      f <- E$value + PR$value
      h <- .shootViolation(E$gaps, scale)
      hmax <- max(1.25 * h, 1)
      if (printIter) cat("natural damping failed; the filter takes over\n")
    }
  }

  # SR1 learns from rejected trials too, and needs their gradients
  trialDeriv <- hm != "gn" && isTRUE(qn$qnRejected)
  lastAccepted <- TRUE
  if (acc != "natural") for (iter in seq_len(iterlim)) {
    if (is.null(C)) C <- condense(E)
    psx <- psxOf()
    opt <- .shootOptimality(thx(theta), C$g0, lbx, ubx, psx)
    if (opt$measure <= tol$gtol && h <= tol$ctol) {
      converged <- TRUE; stopReason <- "gradient"; break
    }

    znorm <- C$znorm
    st <- .shootStep(thx(theta), C, r, lbx, ubx, psx, step$thetaMax, opt$measure,
                     step$regularise %||% 0)
    alpha <- st$alpha
    if (alpha < 1 && znorm > 0) nRelaxed <- nRelaxed + 1L
    dthx <- setNames(st$step, tx)
    dth <- dthx[seq_len(K)]
    dmf <- st$dmf
    thetaT <- theta + dth
    nodesT <- .shootStepNodes(nodes, C, alpha, dthx)
    # the gaps the linearisation predicts for this very step; (1 - alpha) h
    # unless a link was decoupled
    hl <- .shootLinViolation(E, spec, dthx, nodes, nodesT, scale, tx, bkOf()$brk)

    # after an acceptance, one solve with sensitivities serves both the
    # judgement and the next model; after a rejection values suffice
    withD <- trialDeriv || lastAccepted
    ET <- ev(thetaT, nodesT, quiet = TRUE, deriv = withD)
    neval <- neval + 1L
    ok <- !inherits(ET, "error")
    socUsed <- FALSE
    if (ok) {
      PRT <- pv(thetaT)
      fT <- ET$value + PRT$value
      hT <- .shootViolation(ET$gaps, scale)
    }

    # --- acceptance ---------------------------------------------------------
    decide <- function(fT, hT) {
      if (acc == "merit") {
        if (h - hl > 0) {
          req <- dmf / (0.9 * (h - hl))
          if (req > mu) mu <<- 1.1 * req
        }
        pred <- -dmf + mu * (h - hl)
        rho <- if (pred > 0) ((f + mu * h) - (fT + mu * hT)) / pred else -Inf
        list(accept = rho >= 0.25, rho = rho, good = rho > 0.75, htype = FALSE)
      } else {
        dq <- -dmf
        ftype <- dq > 0 && dq >= 1e-4 * h^2
        okF <- .shootFilterOk(filt, h, f, hT, fT, hmax)
        if (ftype) {
          rho <- (f - fT) / dq
          list(accept = okF && rho >= 0.25, rho = rho, good = rho > 0.75, htype = FALSE)
        } else {
          list(accept = okF, rho = NA_real_, good = hT <= 0.5 * h, htype = TRUE)
        }
      }
    }

    d <- if (ok) decide(fT, hT) else list(accept = FALSE, rho = NA_real_, good = FALSE, htype = FALSE)

    # second-order correction: close the trial gaps with the Jacobians at hand
    if (ok && !d$accept && step$soc && alpha > 0 && hT > hl + tol$ctol) {
      nodesS <- .shootSOCNodes(nodesT, C, ET$gaps, spec, bkOf()$brk)
      ES <- ev(thetaT, nodesS, quiet = TRUE, deriv = withD)
      neval <- neval + 1L
      if (!inherits(ES, "error")) {
        fS <- ES$value + PRT$value
        hS <- .shootViolation(ES$gaps, scale)
        dS <- decide(fS, hS)
        if (dS$accept) {
          ET <- ES; fT <- fS; hT <- hS; nodesT <- nodesS; d <- dS
          socUsed <- TRUE; nSOC <- nSOC + 1L
        }
      }
    }

    # an accepted trial judged on values alone gets its sensitivities now; a
    # solve that fails with them is a failed trial like any other
    if (ok && d$accept && !withD) {
      ETd <- ev(thetaT, nodesT, quiet = TRUE)
      neval <- neval + 1L
      if (inherits(ETd, "error")) {
        ok <- FALSE
        ET <- ETd
        d <- list(accept = FALSE, rho = NA_real_, good = FALSE, htype = FALSE)
      } else ET <- ETd
    }

    lastAccepted <- ok && d$accept

    if (printIter)
      cat(sprintf("iter %4d  f = %.8g  h = %.3g  r = %.3g  alpha = %.3g  %s\n",
                  iter, if (ok) fT else NA_real_, if (ok) hT else NA_real_, r,
                  alpha, if (d$accept) "accept" else "reject"))
    if (blather)
      trace[[length(trace) + 1L]] <- list(theta = theta, value = f, violation = h,
        r = r, alpha = alpha, rho = d$rho, accept = d$accept, valtry = if (ok) fT else NA_real_,
        htry = if (ok) hT else NA_real_, htype = d$htype, soc = socUsed,
        steptype = st$steptype, segments = length(spec$segs),
        error = if (ok) NA_character_ else conditionMessage(ET))

    # --- quasi-Newton pairs, from the accepted or (SR1) the rejected trial ----
    if (hm != "gn" && ok && (d$accept || isTRUE(qn$qnRejected))) {
      up <- .shootQNUpdate(B, E, ET, C, alpha, dth, nodes, nodesT, spec, hm, qn$qnCautious)
      B <- up$B
      qnSkipped <- qnSkipped + up$skipped
      if (!d$accept) C <- NULL        # the model changed without a move
    }

    stepnorm <- st$stepnorm
    # a segment the trial could not solve is too long for this step: cut it
    # at the current iterate, where its trajectory is known, and try again
    # from there
    if (!ok && adaptive && inherits(ET, "shootingSegmentError") &&
        splittable(ET$segment) && splitSegments(ET$segment, TRUE)) {
      E2 <- ev(theta, nodes, quiet = TRUE)
      neval <- neval + 1L
      if (!inherits(E2, "error")) {
        E <- E2
        f <- E$value + PR$value
        h <- .shootViolation(E$gaps, scale)
        remapBlocks(E)
        C <- NULL
        next
      }
    }
    if (!ok) {
      n_fail <- n_fail + 1L
      r <- 0.25 * r
      # a step too long for the solver is part of the search here, not a
      # broken objective
      if (n_fail >= 10L) { stopReason <- "objfun"; break }
      next
    }
    n_fail <- 0L

    if (!d$accept) {
      r <- 0.25 * min(r, stepnorm)
      if (abs(fT - f) < tol$ftol && abs(hT - h) <= tol$ctol) n_stall <- n_stall + 1L
      else n_stall <- 0L
      nRejRow <- nRejRow + 1L
      if (isTRUE(step$restore) && h > tol$ctol && nRestore < 5L &&
          (nRejRow >= 4L || n_stall >= 3L)) {
        R <- restoreGaps(if (is.null(C)) condense(E) else C)
        if (!is.null(R)) {
          if (acc == "filter") filt <- .shootFilterAdd(filt, h, f)
          nodes <- R$nodes; E <- R$E
          f <- E$value + PR$value
          h <- .shootViolation(E$gaps, scale)
          remapBlocks(E)
          C <- NULL
          r <- rinit
          nRejRow <- 0L; n_stall <- 0L
          nRestore <- nRestore + 1L
          lastAccepted <- TRUE
          if (printIter)
            cat(sprintf("restore  f = %.8g  h = %.3g\n", f, h))
          next
        }
      }
      if (n_stall >= 5L) { converged <- TRUE; stopReason <- "stagnation"; break }
      if (tol$rmin > 0 && r < tol$rmin) { stopReason <- "radius"; break }
      next
    }

    # accepted
    n_stall <- 0L
    nRejRow <- 0L
    if (acc == "filter" && d$htype) filt <- .shootFilterAdd(filt, h, f)
    if (d$good && stepnorm >= 0.9 * r) r <- min(2 * r, rmax)
    df <- abs(f - fT)
    theta <- thetaT; nodes <- nodesT; E <- ET; PR <- PRT
    f <- fT; h <- hT
    C <- NULL

    # segments grown too long at the new iterate are cut there
    if (adaptive) {
      E2 <- splitAt(violators(E))
      if (!is.null(E2)) {
        E <- E2
        f <- E$value + PR$value
        h <- .shootViolation(E$gaps, scale)
        remapBlocks(E)
      }
    }

    if (hm == "gn" && h <= tol$ctol) {
      if (df < tol$ftol) { converged <- TRUE; stopReason <- "fvalue"; break }
      if (abs(dmf) < tol$mtol) { converged <- TRUE; stopReason <- "preddiff"; break }
      if (tol$xtol > 0 && stepnorm < tol$xtol) { converged <- TRUE; stopReason <- "step"; break }
    }
  }

  if (is.null(C)) C <- condense(E)
  opt <- .shootOptimality(thx(theta), C$g0, lbx, ubx, psxOf())
  # Converged only with a continuous trajectory: open gaps let the nodes absorb
  # misfit the parameters do not explain.
  if (h > tol$ctol) converged <- FALSE

  # trust()'s own fields first, so mstrust() and the parframe tools read a
  # multiple-shooting fit like any other; what only shooting has sits apart.
  msinfo <- list(nodes = nodes, gaps = E$gaps, violation = h,
                   nRelaxed = nRelaxed, nSOC = nSOC, nSplit = nSplit, cuts = cuts,
                   nRestore = nRestore,
                   nStages = nStages,
                   acceptance = if (handedOver) "natural, then filter" else acc,
                   sweep = sweep)
  if (acc == "merit") msinfo$mu <- mu
  else if (acc == "filter") msinfo$filterSize <- nrow(filt)
  if (blather) msinfo$trace <- trace
  list(argument = theta, value = f,
       gradient = setNames(as.numeric(C$g0[seq_len(K)]), tnames),
       hessian = C$H[seq_len(K), seq_len(K), drop = FALSE],
       iterations = iter, neval = neval,
       qnSkipped = qnSkipped,
       converged = converged, stopReason = stopReason,
       atBound = opt$atBound[seq_len(K)],
       multipleShooting = msinfo)
}


# The first of two phases: theta alone by the trust region of trust(), on the
# data term of the segments with every node held where it is, plus the other
# summands. NULL for theta when the phase fails.
.shootPhaseOne <- function(sobj, theta, nodes, fixed, cores, pv, lb, ub, ps,
                           rinit, rmax, iterlim) {
  tnames <- names(theta)
  K <- length(tnames)
  neval <- 0L
  f1 <- function(pars, ...) {
    neval <<- neval + 1L
    E <- sobj(pars, nodes, fixed = fixed, deriv = TRUE, hessian = TRUE,
              sweep = "forward", cores = cores)
    g <- setNames(numeric(K), tnames)
    H <- matrix(0, K, K, dimnames = list(tnames, tnames))
    for (sg in E$segments) {
      v <- intersect(sg$vars, tnames)
      g[v] <- g[v] + sg$grad[v]
      H[v, v] <- H[v, v] + sg$hess[v, v]
    }
    P <- pv(pars)
    list(value = E$value + P$value, gradient = g + P$gradient,
         hessian = H + P$hessian)
  }
  class(f1) <- c("objfn", "fn")
  r <- tryCatch(suppressWarnings(trust(f1, theta, rinit = rinit, rmax = rmax,
                                       iterlim = iterlim, parscale = ps,
                                       parupper = ub, parlower = lb)),
                error = function(e) NULL)
  if (is.null(r) || !all(is.finite(r$argument)))
    return(list(theta = NULL, value = NA_real_, iterations = 0L, neval = neval))
  list(theta = setNames(as.numeric(r$argument), tnames), value = r$value,
       iterations = r$iterations, neval = neval)
}


## ---- Condensing ----------------------------------------------------------------
## Forwards through each condition's segments: S_{j+1} = A_k L_k and
## z_{j+1} = c_j + A_k zhat_k, with A_k the end-state Jacobian of segment k.

.shootCondense <- function(E, H, spec, tnames, scale, zbound = 10, sbound = 1e8,
                           brk = NULL, bcol = NULL) {
  K <- length(tnames)
  g0 <- numeric(K); gz <- numeric(K)
  Hh <- matrix(0, K, K, dimnames = list(tnames, tnames))
  q0 <- 0; q2 <- 0
  nseg <- length(spec$segs)
  L <- vector("list", nseg); zh <- vector("list", nseg)
  Snode <- list(); znode <- list(); Anode <- list()
  znorm2 <- 0
  nDecoupled <- 0L
  VS <- list(); Vz <- list()
  for (cn in spec$conditions) {
    ks <- which(vapply(spec$segs, `[[`, "", "cond") == cn)
    M <- length(ks)
    Sc <- vector("list", M - 1L)
    zc <- matrix(0, M - 1L, length(spec$states), dimnames = list(NULL, spec$states))
    Ac <- vector("list", M)
    S_cur <- NULL; z_cur <- NULL
    for (i in seq_len(M)) {
      k <- ks[i]
      sg <- E$segments[[k]]
      vars <- sg$vars
      Lk <- matrix(0, length(vars), K)
      tp <- match(vars, tnames)
      hit <- which(!is.na(tp))
      if (length(hit)) Lk[cbind(hit, tp[hit])] <- 1
      zk <- numeric(length(vars))
      if (i > 1L) {
        np <- match(spec$segs[[k]]$nodeNames, vars)
        Lk[np, ] <- S_cur
        zk[np] <- z_cur
      }
      Hk <- H[[k]]
      gk <- sg$grad
      HL <- Hk %*% Lk
      Hz <- as.numeric(Hk %*% zk)
      g0 <- g0 + as.numeric(crossprod(Lk, gk))
      Hh <- Hh + crossprod(Lk, HL)
      gz <- gz + as.numeric(crossprod(Lk, Hz))
      q0 <- q0 + sum(gk * zk)
      q2 <- q2 + sum(zk * Hz)
      L[[k]] <- Lk; zh[[k]] <- zk
      if (i < M && isTRUE(brk[[cn]][i])) {
        # a continuity break: the next node is a variable of its own, its
        # step the columns bcol of the extended parameters, and no gap
        Ac[[i]] <- sg$jac
        S_cur <- matrix(0, length(spec$states), K)
        S_cur[cbind(seq_along(spec$states), bcol[[cn]][[as.character(i)]])] <- 1
        z_cur <- numeric(length(spec$states))
        Sc[[i]] <- S_cur
        zc[i, ] <- z_cur
        VS[[length(VS) + 1L]] <- S_cur / scale[spec$states]
        Vz[[length(Vz) + 1L]] <- z_cur
      } else if (i < M) {
        A <- sg$jac
        Ac[[i]] <- A
        S_cur <- A %*% Lk
        z_cur <- as.numeric(E$gaps[[cn]][i, spec$states]) + as.numeric(A %*% zk)
        # A propagated gap or sensitivity past its bound is outside the
        # linearisation: the link to the node before is dropped, an inexact
        # step the acceptance judges on the true gaps.
        if (i > 1L && (max(abs(z_cur / scale[spec$states])) > zbound ||
                       max(abs(S_cur / scale[spec$states])) > sbound)) {
          np <- match(spec$segs[[k]]$nodeNames, vars)
          Ad <- A
          Ad[, np] <- 0
          S_cur <- Ad %*% Lk
          z_cur <- as.numeric(E$gaps[[cn]][i, spec$states]) + as.numeric(Ad %*% zk)
          nDecoupled <- nDecoupled + 1L
        }
        Sc[[i]] <- S_cur
        zc[i, ] <- z_cur
        vS <- S_cur / scale[spec$states]
        vz <- z_cur / scale[spec$states]
        VS[[length(VS) + 1L]] <- vS
        Vz[[length(Vz) + 1L]] <- vz
        znorm2 <- znorm2 + sum(vz^2)
      }
    }
    Snode[[cn]] <- Sc; znode[[cn]] <- zc; Anode[[cn]] <- Ac
  }
  Hh <- 0.5 * (Hh + t(Hh))
  list(g0 = setNames(g0, tnames), gz = setNames(gz, tnames), H = Hh,
       q0 = q0, q2 = q2, L = L, zh = zh, S = Snode, z = znode, A = Anode,
       znorm = sqrt(znorm2), decoupled = nDecoupled,
       VS = if (length(VS)) do.call(rbind, VS) else matrix(0, 0, K),
       Vz = unlist(Vz))
}

# The node step of a trial: alpha z + S Delta theta.
.shootStepNodes <- function(nodes, C, alpha, dth) {
  out <- nodes
  for (cn in names(C$S)) {
    Sc <- C$S[[cn]]
    for (i in seq_along(Sc))
      out[[cn]][i, colnames(C$z[[cn]])] <- nodes[[cn]][i, colnames(C$z[[cn]])] +
        alpha * C$z[[cn]][i, ] + as.numeric(Sc[[i]] %*% dth)
  }
  out
}

# Second-order correction: the trial gaps pushed forward through the end-state
# Jacobians of the current iterate, with Delta theta held.
.shootSOCNodes <- function(nodesT, C, gapsT, spec, brk = NULL) {
  out <- nodesT
  for (cn in names(C$S)) {
    ks <- which(vapply(spec$segs, `[[`, "", "cond") == cn)
    Ac <- C$A[[cn]]
    zp <- NULL
    for (i in seq_len(length(ks) - 1L)) {
      k <- ks[i]
      # nothing is closed across a break, and nothing passed over it
      if (isTRUE(brk[[cn]][i])) zp <- numeric(length(spec$states))
      else if (i == 1L) zp <- as.numeric(gapsT[[cn]][i, spec$states])
      else {
        np <- match(spec$segs[[k]]$nodeNames, colnames(Ac[[i]]))
        zp <- as.numeric(gapsT[[cn]][i, spec$states]) +
              as.numeric(Ac[[i]][, np, drop = FALSE] %*% zp)
      }
      out[[cn]][i, spec$states] <- out[[cn]][i, spec$states] + zp
    }
  }
  out
}


# Size of the gaps the linearised constraints leave after a step,
# c_j + A_j d_j - Delta s_{j+1}, scaled per state, with A_j the full end-state
# Jacobian: the exact linearisation, also where condensing decoupled a link.
.shootLinViolation <- function(E, spec, dth, nodes, nodesT, scale, tnames, brk = NULL) {
  s2 <- 0
  st <- spec$states
  for (cn in spec$conditions) {
    ks <- which(vapply(spec$segs, `[[`, "", "cond") == cn)
    M <- length(ks)
    for (i in seq_len(M - 1L)) {
      if (isTRUE(brk[[cn]][i])) next
      k <- ks[i]
      sg <- E$segments[[k]]
      d <- numeric(length(sg$vars))
      tp <- match(sg$vars, tnames)
      d[!is.na(tp)] <- dth[tp[!is.na(tp)]]
      if (i > 1L) {
        np <- match(spec$segs[[k]]$nodeNames, sg$vars)
        d[np] <- nodesT[[cn]][i - 1L, st] - nodes[[cn]][i - 1L, st]
      }
      r <- as.numeric(E$gaps[[cn]][i, st]) + as.numeric(sg$jac %*% d) -
           (nodesT[[cn]][i, st] - nodes[[cn]][i, st])
      s2 <- s2 + sum((r / scale[st])^2)
    }
  }
  sqrt(s2)
}


## ---- Partitioned quasi-Newton ------------------------------------------------
## Multipliers by lambda_{j-1} = grad_{s_j} m_j + G_j' lambda_j; segment j's pair
## is its local step and the change of the gradient of f_j + lambda_j' end_j.

.shootQNUpdate <- function(B, E, ET, C, alpha, dth, nodes, nodesT, spec, hm, cautious) {
  skipped <- 0L
  nseg <- length(spec$segs)
  d <- lapply(seq_len(nseg), function(k)
    as.numeric(C$L[[k]] %*% dth) + alpha * C$zh[[k]])
  # the node part of the actual step, which a correction may have moved
  for (k in seq_len(nseg)) {
    sg <- spec$segs[[k]]
    if (sg$j == 1L) next
    np <- match(sg$nodeNames, E$segments[[k]]$vars)
    d[[k]][np] <- nodesT[[sg$cond]][sg$j - 1L, spec$states] -
                  nodes[[sg$cond]][sg$j - 1L, spec$states]
  }
  lam <- vector("list", nseg)
  for (cn in spec$conditions) {
    ks <- which(vapply(spec$segs, `[[`, "", "cond") == cn)
    M <- length(ks)
    nxt <- NULL
    for (i in rev(seq_len(M))) {
      k <- ks[i]
      if (i < M) lam[[k]] <- nxt
      if (i == 1L) break
      sg <- E$segments[[k]]
      np <- match(spec$segs[[k]]$nodeNames, sg$vars)
      gm <- sg$grad + as.numeric(B[[k]] %*% d[[k]])
      v <- gm[np]
      if (!is.null(lam[[k]])) v <- v + as.numeric(crossprod(sg$jac[, np, drop = FALSE], lam[[k]]))
      nxt <- v
    }
  }
  for (k in seq_len(nseg)) {
    sg <- E$segments[[k]]; sgT <- ET$segments[[k]]
    y <- sgT$grad[sg$vars] - sg$grad
    if (!is.null(lam[[k]]))
      y <- y + as.numeric(crossprod(sgT$jac[, sg$vars, drop = FALSE] - sg$jac, lam[[k]]))
    s <- d[[k]]
    if (hm == "bfgs") {
      sy <- sum(s * y)
      if (!(sy > cautious * sqrt(sum(s^2) * sum(y^2)))) { skipped <- skipped + 1L; next }
    }
    B[[k]] <- qn_update_impl(B[[k]], s, as.numeric(y), hm)$B
  }
  list(B = B, skipped = skipped)
}


## ---- Acceptance --------------------------------------------------------------

.shootViolation <- function(gaps, scale) {
  s <- 0
  for (G in gaps) if (length(G)) s <- s + sum(sweep(G, 2L, scale[colnames(G)], "/")^2)
  sqrt(s)
}

# Fletcher-Leyffer: a trial is acceptable if against every filter entry, and
# the current point, it improves the gaps by a margin or the data term by one.
.shootFilterOk <- function(filt, h, f, hT, fT, hmax, beta = 0.99, gamma = 1e-5) {
  if (!(hT <= hmax)) return(FALSE)
  pts <- rbind(filt, c(h, f))
  all(hT <= beta * pts[, 1L] | fT <= pts[, 2L] - gamma * hT)
}

.shootFilterAdd <- function(filt, h, f) {
  keep <- !(filt[, 1L] >= h & filt[, 2L] >= f)
  rbind(filt[keep, , drop = FALSE], c(h, f))
}


## ---- Small helpers -------------------------------------------------------------

# The full generalised Gauss-Newton step (alpha = 1) on the pseudo-inverse,
# directions below reg dropped, and its scaled full-space vector. With `g`, `Vz`
# of another linearisation, the simplified step of the natural level function.
.shootFullStep <- function(C, theta, lb, ub, ps, reg, g = C$g0 + C$gz, Vz = C$Vz,
                           eig = NULL, truncate = TRUE) {
  if (is.null(eig)) {
    Hs <- C$H / tcrossprod(ps)
    eig <- eigen(0.5 * (Hs + t(Hs)), symmetric = TRUE)
  }
  lmax <- max(abs(eig$values), .Machine$double.eps)
  keep <- eig$values > max(reg, 1e-7)^2 * lmax
  V <- eig$vectors[, keep, drop = FALSE]
  dth <- -as.numeric(V %*% (crossprod(V, g / ps) / eig$values[keep])) / ps
  if (truncate) {
    # strictly inside the box, the node part still closing the gaps
    frac <- 1
    for (i in which(dth != 0)) {
      b <- if (dth[i] > 0) ub[i] else lb[i]
      if (is.finite(b)) frac <- min(frac, 0.99 * (b - theta[i]) / dth[i])
    }
    dth <- dth * frac
  }
  vec <- c(ps * dth, as.numeric(C$VS %*% dth) + Vz)
  list(dth = dth, vec = vec, norm = sqrt(sum(vec^2)), eig = eig)
}

# One step in the region ||M Delta theta + alpha q||^2 <= r^2, M = [D; V S],
# q = [0; V z], D scaled Coleman-Li. M is factored by QR, never squared, as S
# can span many decades; alpha is the node fraction of the full step that fits.
.shootStep <- function(theta, C, r, lb, ub, ps, thetamax, optMeasure, reg = 0) {
  K <- length(theta)
  gz0 <- C$g0 / ps
  bnd <- ifelse(gz0 < 0, ub, lb)
  fin <- is.finite(bnd)
  absv <- ifelse(fin, abs(theta - bnd) * ps, 1)
  jv <- as.numeric(fin)
  Hm <- C$H + diag(abs(gz0) * jv * ps^2 / absv, K)
  M <- rbind(diag(ps / sqrt(absv), K), C$VS)
  q <- c(numeric(K), C$Vz)
  qm <- qr(M, LAPACK = TRUE)
  # M = Q Rf with Rf = R P' for the column pivoting P, so Rf^-1 = D^-1 Q_theta
  # from the first K rows of Q, without a triangular solve
  Q <- qr.Q(qm)
  Rinv <- Q[seq_len(K), , drop = FALSE] / (ps / sqrt(absv))
  qq <- qr.qty(qm, q)
  Qtq <- qq[seq_len(K)]
  w1 <- as.numeric(Rinv %*% Qtq)                         # argmin ||M w - q||
  # the residual from the rest of Q'q; base R's qr.resid() takes no LAPACK QR
  res2 <- sum(qq[-seq_len(K)]^2)

  # length of the full Gauss-Newton step, alpha = 1, in the full space; the
  # model Hessian may be singular or indefinite, so its spectrum is floored
  g1 <- C$g0 + C$gz
  ev <- eigen(0.5 * (Hm + t(Hm)), symmetric = TRUE)
  lmax <- max(abs(ev$values), .Machine$double.eps)
  dGN <- -as.numeric(ev$vectors %*% (crossprod(ev$vectors, g1) /
                                        pmax(ev$values, 1e-10 * lmax)))
  # Only the node part of that step measures how far the gaps may close; a
  # sloppy parameter direction lengthens it without moving a node.
  nGN <- sqrt(sum(((M %*% dGN + q)[-seq_len(K)])^2))
  alpha <- if (C$znorm == 0) 0 else min(1, r / max(nGN, sqrt(res2)))

  # Delta theta = Rf^-1 u - alpha w1, with ||u|| <= rho
  ga <- C$g0 + alpha * C$gz
  w0 <- alpha * w1
  rho <- sqrt(max(0, r^2 - alpha^2 * res2))
  gu <- crossprod(Rinv, ga - Hm %*% w0)
  Hu <- crossprod(Rinv, Hm %*% Rinv)
  # Peifer and Timmer (2007), sec. 4: a direction with singular value below
  # reg times the largest loses its gradient and is damped as if fixed.
  if (reg > 0) {
    eh <- eigen(0.5 * (Hu + t(Hu)), symmetric = TRUE)
    lmax <- max(abs(eh$values), .Machine$double.eps)
    sing <- abs(eh$values) < reg^2 * lmax
    if (any(sing)) {
      Vs <- eh$vectors[, sing, drop = FALSE]
      gu <- gu - Vs %*% crossprod(Vs, gu)
      Hu <- Hu + lmax * tcrossprod(Vs)
    }
  }
  st <- trust_step_impl(numeric(K), as.numeric(gu), 0.5 * (Hu + t(Hu)), rho,
                        rep(-Inf, K), rep(Inf, K), rep(1, K), thetamax)
  dth <- as.numeric(Rinv %*% st$step) - w0

  # stay strictly inside the box
  frac <- 1
  for (i in which(dth != 0)) {
    b <- if (dth[i] > 0) ub[i] else lb[i]
    if (is.finite(b)) frac <- min(frac, (b - theta[i]) / dth[i])
  }
  if (frac < 1) dth <- dth * max(0, min(max(thetamax, 1 - optMeasure), 1 - 1e-12) * frac)

  dmf <- alpha * C$q0 + 0.5 * alpha^2 * C$q2 + sum(ga * dth) +
         0.5 * sum(dth * (Hm %*% dth))
  norm <- sqrt(sum((M %*% dth + alpha * q)^2))
  list(step = dth, alpha = alpha, dmf = dmf, stepnorm = norm,
       steptype = if (frac < 1) "truncated" else st$steptype)
}

.shootBound <- function(b, nms, default) {
  out <- setNames(rep(default, length(nms)), nms)
  if (is.null(b) || !length(b)) return(out)
  if (!is.null(names(b))) {
    hit <- intersect(names(b), nms)
    out[hit] <- b[hit]
  } else out[] <- b[1L]
  out
}

.shootInterior <- function(x, lb, ub) {
  eps <- 100 * .Machine$double.eps * pmax(1, abs(x))
  both <- is.finite(lb) & is.finite(ub) & (ub - lb <= 2 * eps)
  x[both] <- 0.5 * (lb[both] + ub[both])
  hi <- is.finite(ub) & (ub - x < eps) & !both
  lo <- is.finite(lb) & (x - lb < eps) & !both
  x[hi] <- ub[hi] - eps[hi]
  x[lo] <- lb[lo] + eps[lo]
  x
}

# The optimality measure of the box problem, as trust() computes it.
.shootOptimality <- function(x, g, lb, ub, ps) {
  gz <- g / ps
  b <- ifelse(gz < 0, ub, lb)
  fin <- is.finite(b)
  absv <- ifelse(fin, abs(x * ps - b * ps), 1)
  m <- max(abs(absv * gz))
  btol <- 1e-6
  list(measure = m,
       atBound = setNames(fin & abs(absv * gz) <= btol & abs(gz) > btol, names(x)))
}

# The scale of each state's gap: 1 on a log10 or angle chart; on a linear
# chart the largest range within one segment, from the data for an observed
# state, else from a simulation of the whole time axis at the start.
.shootScale <- function(scale, nodes, spec, sobj = NULL, theta = NULL,
                        fixed = NULL, cores = 1L) {
  states <- spec$states
  out <- setNames(rep(1, length(states)), states)
  lin <- states[spec$charts == "linear"]
  rng <- function(v) {
    v <- v[is.finite(v)]
    if (length(v) > 1L) diff(range(v)) else 0
  }
  obs <- .shootObservedStates(spec)
  sim <- NULL
  if (length(setdiff(lin, names(obs))) && !is.null(sobj))
    sim <- .shootSimulate(sobj, theta, fixed, cores)
  for (st in lin) {
    rg <- 0
    if (st %in% names(obs) && !is.null(sobj)) {
      dat <- attr(sobj, "args")$data
      rg <- rng(unlist(lapply(spec$conditions, function(cn) {
        d <- dat[[cn]]
        d$value[as.character(d$name) == obs[[st]]]
      })))
    } else if (!is.null(sim)) {
      loc <- unlist(lapply(spec$segs, function(sg) {
        m <- sim[[sg$cond]]
        if (is.null(m)) return(NULL)
        rng(m[m[, "time"] >= sg$start & m[, "time"] <= sg$end, st])
      }))
      rg <- if (length(loc)) max(loc) else 0
    }
    if (!(rg > 0)) rg <- rng(unlist(lapply(nodes, function(N) N[, st])))
    if (is.finite(rg) && rg > 0) out[st] <- rg
  }
  if (!is.null(scale)) {
    bad <- setdiff(names(scale), states)
    if (length(bad)) stop("normL2: multipleShootingControl$scale names ", paste(bad, collapse = ", "),
                          ", which are not states.", call. = FALSE)
    out[names(scale)] <- scale
  }
  out
}


## ---- Annealing of the continuity ---------------------------------------------
## Precision annealing (Ye and Abarbanel): f + w ||V c||^2 over parameters and
## nodes, uncondensed and sparse, with w rising per stage before the exact method.

# Gauss-Newton model of f + w ||V c||^2 in the scaled variables
# y = (ps * theta, s / scale), as a sparse matrix and a gradient.
.shootAnnealModel <- function(E, spec, tnames, ps, scale, w, PR, brk = NULL) {
  K <- length(tnames)
  states <- spec$states
  nx <- length(states)
  # full index of every node variable
  off <- K
  nodeIdx <- list()
  for (cn in spec$conditions) {
    ks <- which(vapply(spec$segs, `[[`, "", "cond") == cn)
    for (i in seq_along(ks)[-1L]) {
      nodeIdx[[spec$segs[[ks[i]]]$key]] <- off + seq_len(nx)
      off <- off + nx
    }
  }
  n <- off
  dsc <- c(ps, rep(1 / scale[states], (n - K) / nx))   # y = dsc * x
  # blocks are collected and bound once, so assembly is linear in the segments
  blk <- vector("list", 2L * length(spec$segs) + 1L)
  nb <- 0L
  g <- numeric(n)
  add <- function(idx, Hl, gl) {
    nb <<- nb + 1L
    blk[[nb]] <<- list(idx = idx, x = Hl)
    g[idx] <<- g[idx] + gl
  }
  locIdx <- function(k) {
    sg <- E$segments[[k]]
    idx <- match(sg$vars, tnames)
    if (spec$segs[[k]]$j > 1L) {
      np <- match(spec$segs[[k]]$nodeNames, sg$vars)
      idx[np] <- nodeIdx[[spec$segs[[k]]$key]]
    }
    idx
  }
  cont <- 0
  for (cn in spec$conditions) {
    ks <- which(vapply(spec$segs, `[[`, "", "cond") == cn)
    for (i in seq_along(ks)) {
      k <- ks[i]
      sg <- E$segments[[k]]
      idx <- locIdx(k)
      keep <- !is.na(idx)
      add(idx[keep], sg$hess[keep, keep, drop = FALSE], sg$grad[keep])
      if (i < length(ks) && !isTRUE(brk[[cn]][i])) {
        # c_j = end_j(theta, s_j) - s_{j+1}, weighted by sqrt(w) V
        cj <- as.numeric(E$gaps[[cn]][i, states]) / scale[states]
        Cj <- cbind(sg$jac[, keep, drop = FALSE] / scale[states], -diag(1 / scale[states], nx))
        cidx <- c(idx[keep], nodeIdx[[spec$segs[[ks[i + 1L]]]$key]])
        add(cidx, 2 * w * crossprod(Cj), 2 * w * as.numeric(crossprod(Cj, cj)))
        cont <- cont + sum(cj^2)
      }
    }
  }
  # the summands on theta alone
  add(seq_len(K), PR$hessian, PR$gradient)
  blk <- blk[seq_len(nb)]
  ii <- unlist(lapply(blk, function(b) rep(b$idx, times = length(b$idx))), use.names = FALSE)
  jj <- unlist(lapply(blk, function(b) rep(b$idx, each = length(b$idx))), use.names = FALSE)
  xx <- unlist(lapply(blk, function(b) as.numeric(b$x)), use.names = FALSE)
  # to the scaled variables, H_y = D^-1 H D^-1 and g_y = D^-1 g; every block is
  # symmetric, so the upper triangle is the matrix
  up <- ii <= jj
  Hy <- Matrix::sparseMatrix(i = ii[up], j = jj[up], x = xx[up] / (dsc[ii[up]] * dsc[jj[up]]),
                             dims = c(n, n), symmetric = TRUE)
  list(H = Hy, g = g / dsc, dsc = dsc, nodeIdx = nodeIdx, n = n,
       value = E$value + PR$value + w * cont, cont = cont)
}

# Moré-Sorensen on a sparse positive semi-definite model: the step of
# min g'p + p'Hp/2 over ||p|| <= r, by Newton on 1/||p(lambda)|| = 1/r with
# sparse Cholesky factors of H + lambda I.
.shootSparseTR <- function(H, g, r, maxit = 20L) {
  n <- length(g)
  dmax <- max(abs(Matrix::diag(H)), .Machine$double.eps)
  lam <- 1e-12 * dmax
  L <- Matrix::Cholesky(H, perm = TRUE, LDL = FALSE, Imult = lam)
  p <- -as.numeric(Matrix::solve(L, g, system = "A"))
  np <- sqrt(sum(p^2))
  if (np <= r) return(list(step = p, lambda = lam))
  for (it in seq_len(maxit)) {
    q2 <- sum(p * as.numeric(Matrix::solve(L, p, system = "A")))
    # Newton on 1/||p|| approaches from below and does not overshoot; the
    # guard only keeps lambda positive
    lnew <- lam + (np / r - 1) * np^2 / q2
    lam <- if (lnew > 0) lnew else lam / 10
    L <- Matrix::update(L, H, mult = lam)
    p <- -as.numeric(Matrix::solve(L, g, system = "A"))
    np <- sqrt(sum(p^2))
    if (abs(np - r) <= 1e-3 * r) break
  }
  if (np > r) p <- p * (r / np)
  list(step = p, lambda = lam)
}

# Annealing stages, each minimising f + w ||V c||^2 until a step gains less
# than `ftol`, until the RMS gap per node is below `hnode`, w reaches `wmax`, or
# the gaps stall. Returns the point the exact method starts from.
.shootAnneal <- function(evGN, spec, theta, nodes, tnames, ps, scale, pv, lb, ub,
                         r, rmax, printIter = FALSE, w0 = 1, factor = 10,
                         wmax = 1e10, hnode = 1e-3, itStage = 30L, ftol = 1e-3,
                         evV = evGN, brk = NULL) {
  K <- length(tnames)
  states <- spec$states
  nx <- length(states)
  # node rows behind the full indices
  rows <- list()
  for (cn in spec$conditions) {
    ks <- which(vapply(spec$segs, `[[`, "", "cond") == cn)
    for (i in seq_along(ks)[-1L])
      rows[[length(rows) + 1L]] <- list(cond = cn, row = i - 1L)
  }
  nnode <- max(1L, length(rows))
  moveNodes <- function(nodes, dx) {
    for (m in seq_along(rows)) {
      idx <- K + (m - 1L) * nx + seq_len(nx)
      nodes[[rows[[m]]$cond]][rows[[m]]$row, states] <-
        nodes[[rows[[m]]$cond]][rows[[m]]$row, states] + dx[idx]
    }
    nodes
  }
  E <- evGN(theta, nodes)
  if (inherits(E, "error")) return(NULL)
  neval <- 1L
  PR <- pv(theta)
  w <- w0
  stages <- 0L
  r0 <- r
  lastAccepted <- TRUE
  repeat {
    stages <- stages + 1L
    # a new weight is a new problem; the radius the last one ended on says
    # nothing about it
    r <- r0
    Mdl <- .shootAnnealModel(E, spec, tnames, ps, scale, w, PR, brk)
    cont0 <- Mdl$cont
    for (it in seq_len(itStage)) {
      st <- .shootSparseTR(Mdl$H, Mdl$g, r)
      dy <- st$step
      pred <- -(sum(Mdl$g * dy) + 0.5 * sum(dy * as.numeric(Mdl$H %*% dy)))
      dx <- dy / Mdl$dsc
      dth <- dx[seq_len(K)]
      # stay inside the box on theta
      frac <- 1
      for (i in which(dth != 0)) {
        b <- if (dth[i] > 0) ub[i] else lb[i]
        if (is.finite(b)) frac <- min(frac, 0.99995 * (b - theta[i]) / dth[i])
      }
      dx <- dx * frac; dy <- dy * frac; pred <- pred * frac
      thetaT <- theta + dx[seq_len(K)]
      nodesT <- moveNodes(nodes, dx)
      # after a rejection judged on values alone, the model built only for a
      # step taken; after an acceptance with the model at once
      EV <- if (lastAccepted) evGN(thetaT, nodesT) else evV(thetaT, nodesT)
      neval <- neval + 1L
      if (inherits(EV, "error") || !(pred > 0)) {
        r <- 0.25 * r; lastAccepted <- FALSE; next
      }
      PRT <- pv(thetaT)
      contT <- .shootViolation(EV$gaps, scale)^2
      valT <- EV$value + PRT$value + w * contT
      rho <- (Mdl$value - valT) / pred
      ET <- EV
      if (rho >= 0.25 && !lastAccepted) {
        ET <- evGN(thetaT, nodesT)
        neval <- neval + 1L
        if (inherits(ET, "error")) rho <- -Inf
      }
      lastAccepted <- rho >= 0.25
      if (printIter)
        cat(sprintf("anneal w = %.0e  F = %.8g  gaps = %.3g  r = %.3g  %s\n", w,
                    valT, sqrt(contT), r, if (rho >= 0.25) "accept" else "reject"))
      if (rho >= 0.25) {
        MT <- .shootAnnealModel(ET, spec, tnames, ps, scale, w, PRT, brk)
        dF <- Mdl$value - MT$value
        theta <- thetaT; nodes <- nodesT; E <- ET; PR <- PRT; Mdl <- MT
        if (rho > 0.75 && sqrt(sum(dy^2)) >= 0.9 * r) r <- min(2 * r, rmax)
        if (dF <= ftol * max(1, abs(Mdl$value))) break
      } else r <- 0.25 * min(r, sqrt(sum(dy^2)))
    }
    if (sqrt(Mdl$cont / nnode) <= hnode || w >= wmax) break
    # with the penalty dominant and the gaps barely shrunk, a larger weight
    # only stiffens the problem
    if (w * Mdl$cont > max(1, abs(Mdl$value - w * Mdl$cont)) && Mdl$cont > 0.81 * cont0) break
    w <- w * factor
  }
  list(theta = theta, nodes = nodes, neval = neval, stages = stages, r = r, w = w)
}
