# FOCEI: -2 log L ~ J(theta, eta*) + sum_g log|H_g / 2| around the conditional modes of every
# group, J the joint objective, H_g the Gauss-Newton Hessian of J in the etas of g. The gradient has the
# exact derivative of the volume term; the guide chapter Nonlinear mixed-effects models has the derivation.


## Derivatives up to order 4 of the -2 log likelihood of a row in the row prediction p and sigma s,
## normal or Student-t with df degrees of freedom, M3 below the limit: [n, 5, 5], entry
## [, i + 1, j + 1] = d^(i+j) / dp^i ds^j, the value itself left at 0. Normal rows above the
## limit in closed form, the others from the Taylor series in (dp, ds).
.rowLossDerivs <- function(p, s, y, censored, df = Inf) {
  if (is.finite(df)) return(.rowLossTaylor(p, s, y, censored, df))
  out <- .rowLossNormal(p, s, y)
  if (any(censored))
    out[censored, , ] <- .rowLossTaylor(p[censored], s[censored], y[censored], censored[censored], df)
  out
}

## (p - y)^2 / s^2 + 2 log s: the first term is u^2, 2 u, 2 in p times d^j s^-2 / ds^j.
.rowLossNormal <- function(p, s, y) {
  n <- length(p)
  u <- p - y
  out <- array(0, c(n, 5L, 5L))
  inP <- list(u^2, 2 * u, rep(2, n))
  for (i in 0:2) for (j in 0:(4 - i)) {
    if (i + j == 0L) next
    out[, i + 1L, j + 1L] <- inP[[i + 1L]] * (-1)^j * factorial(j + 1L) * s^(-2 - j)
  }
  for (j in 1:4) out[, 1L, j + 1L] <- out[, 1L, j + 1L] + 2 * (-1)^(j + 1L) * factorial(j - 1L) * s^(-j)
  out
}

.rowLossTaylor <- function(p, s, y, censored, df = Inf) {
  n <- length(p)
  mul <- function(A, B) {
    C <- array(0, c(n, 5L, 5L))
    for (i in 0:4) for (j in 0:(4 - i)) for (a in 0:i) for (b in 0:j)
      C[, i + 1L, j + 1L] <- C[, i + 1L, j + 1L] + A[, a + 1L, b + 1L] * B[, i - a + 1L, j - b + 1L]
    C
  }
  w <- (p - y) / s
  # w(p + dp, s + ds) - w = (u + dp) / (s + ds) - u / s, u = p - y
  D <- array(0, c(n, 5L, 5L))
  for (k in 0:4) {
    if (k > 0L) D[, 1L, k + 1L] <- (p - y) * (-1)^k / s^(k + 1L)
    if (k < 4L) D[, 2L, k + 1L] <- (-1)^k / s^(k + 1L)
  }
  # derivatives in w of the row loss above the limit: w^2, or (nu + 1) log(1 + w^2 / nu)
  q <- df + w^2
  m <- if (is.finite(df))
    cbind(2 * (df + 1) * w / q, 2 * (df + 1) * (df - w^2) / q^2,
          -4 * (df + 1) * w * (3 * df - w^2) / q^3,
          -12 * (df + 1) * (df^2 - 6 * df * w^2 + w^4) / q^4)
  else cbind(2 * w, 2, 0, 0)
  # below the limit -2 log F(-w): G = f(w) / F(-w), G' = G (h' + G) with h = log f
  if (any(censored)) {
    wc <- w[censored]
    if (is.finite(df)) {
      G0 <- exp(stats::dt(wc, df, log = TRUE) - stats::pt(-wc, df, log.p = TRUE))
      h <- -m[censored, 1:3, drop = FALSE] / 2
    } else {
      G0 <- exp(stats::dnorm(wc, log = TRUE) - stats::pnorm(-wc, log.p = TRUE))
      h <- cbind(-wc, -1, 0)
    }
    G1 <- G0 * (h[, 1] + G0)
    G2 <- G1 * (h[, 1] + G0) + G0 * (h[, 2] + G1)
    G3 <- G2 * (h[, 1] + G0) + 2 * G1 * (h[, 2] + G1) + G0 * (h[, 3] + G2)
    m[censored, ] <- 2 * cbind(G0, G1, G2, G3)
  }
  out <- array(0, c(n, 5L, 5L))
  Dk <- D
  for (k in 1:4) {
    out <- out + m[, k] / factorial(k) * Dk
    if (k < 4L) Dk <- mul(Dk, D)
  }
  # 2 log s above the limit
  a <- !censored
  for (k in 1:4) out[a, 1L, k + 1L] <- out[a, 1L, k + 1L] + 2 * (-1)^(k + 1L) / (k * s[a]^k)
  for (i in 0:4) for (j in 0:(4 - i)) out[, i + 1L, j + 1L] <- out[, i + 1L, j + 1L] * factorial(i) * factorial(j)
  out
}


## The error model of data with recorded sigmas: per data row the sigma on the prediction grid,
## without derivatives.
.dataSigmaErrfn <- function(data) {
  function(out, pars = NULL, conditions = NULL, ...) {
    cn <- conditions[1]
    frame <- if (inherits(out, "prdlist")) out[[cn]] else out
    d <- data[[cn]]
    if (anyNA(d$sigma))
      stop("em: the data of condition ", cn, " have no sigma; give an error model.", call. = FALSE)
    times <- as.numeric(frame[, "time"])
    obs <- setdiff(colnames(frame), "time")
    err <- matrix(NA_real_, length(times), length(obs) + 1L, dimnames = list(NULL, c("time", obs)))
    err[, "time"] <- times
    for (o in intersect(obs, as.character(d$name))) {
      rows <- d[as.character(d$name) == o, , drop = FALSE]
      i <- match(times, rows$time)
      err[, o] <- ifelse(is.na(i), stats::median(rows$sigma), rows$sigma[i])
    }
    setNames(list(structure(err, class = c("prdframe", "matrix"))), cn)
  }
}

.errorModel <- function(design) design$errfn %||% .dataSigmaErrfn(design$data)

## A^-1 B for a symmetric eta block A, solved on the unit diagonal of A: where the variance of a random
## effect runs to 0, the prior precision outgrows the other entries by many orders of magnitude.
.scaledSolve <- function(A, B = NULL) {
  if (is.null(B)) {
    B <- diag(nrow(A))
    dimnames(B) <- list(colnames(A), rownames(A))
  }
  d <- 1 / sqrt(abs(diag(A)))
  d[!is.finite(d)] <- 1
  d * solve(d * A * rep(d, each = nrow(A)), d * as.matrix(B))
}

## The error-model frame of condition `cn` at the prediction `pred`; the frame of the recorded
## sigmas where the error model does not cover the condition.
.errorFrame <- function(design, pred, cn) {
  out <- if (!is.null(design$errfn)) design$errfn(out = pred, pars = attr(pred, "parameters"), conditions = cn)
  out[[cn]] %||% .dataSigmaErrfn(design$data)(out = pred, conditions = cn)[[1]]
}

## The prediction grid of condition `cn`: the grid of the data objective, so that every stage
## predicts the data rows as the objective does.
.emGrid <- function(design, cn) {
  grid <- attr(design$objfn, "timesD", exact = TRUE)
  if (is.null(grid)) sort(unique(design$data[[cn]]$time)) else .gridOf(grid, cn)
}


## The data of a condition with the limit per row: censored rows have value <= lloq.
.conditionData <- function(design, cn) {
  d <- design$data[[cn]]
  d$name <- as.character(d$name)
  lloq <- d$lloq %||% rep(-Inf, nrow(d))
  lloq[is.na(lloq)] <- -Inf
  d$lloq <- lloq
  d
}


## Kernel description of one condition block: the rows, above the limit first, and their
## indices into the prediction, the derivatives and the error model at `pars`. A row with a
## recorded sigma uses it, the others the error model, as in normL2().
.conditionBlock <- function(design, cn, pars, fixed, frame = NULL) {
  d <- .conditionData(design, cn)
  cens <- d$value <= d$lloq
  d <- d[order(cens, d$time, d$name), , drop = FALSE]
  times <- .emGrid(design, cn)
  pred <- frame$pred %||% design$prdfn(times = times, pars = pars, fixed = fixed, deriv = TRUE, conditions = cn)[[1]]
  err <- frame$err %||% .errorFrame(design, pred, cn)
  dn <- dimnames(attr(pred, "deriv"))
  en <- dimnames(attr(err, "deriv"))
  obs <- match(d$name, colnames(pred))
  if (anyNA(obs))
    stop("em: observables of condition ", cn, " missing in the prediction: ",
         paste(setdiff(d$name, colnames(pred)), collapse = ", "), call. = FALSE)
  recorded <- !is.na(d$sigma)
  oe <- match(d$name, colnames(err))
  if (any(is.na(oe) & !recorded))
    stop("em: the error model of condition ", cn, " has no sigma for ",
         paste(unique(d$name[is.na(oe) & !recorded]), collapse = ", "), call. = FALSE)
  oe[recorded & is.na(oe)] <- 1L
  errDeriv <- if (!is.null(en)) match(d$name, en[[2]]) else rep(0L, nrow(d))
  errDeriv[is.na(errDeriv) | recorded] <- 0L
  list(times = as.numeric(times), condition = cn,
       t_idx_in_pred = .matchNum(d$time, as.numeric(pred[, "time"])), o_idx_in_pred = obs,
       o_idx_in_deriv = match(d$name, dn[[2]]),
       t_idx_in_err = .matchNum(d$time, as.numeric(err[, "time"])), o_idx_in_err = oe,
       o_idx_in_err_deriv = as.integer(errDeriv), sigma_data = as.numeric(d$sigma),
       y_data = as.numeric(pmax(d$value, d$lloq)), lloq = as.numeric(d$lloq),
       n_aloq = sum(d$value > d$lloq), t_df = design$df,
       parameters = unique(c(dn[[3]], if (any(!recorded)) en[[3]])))
}


## Kernel description of the groups: per group the etas and condition blocks, the layout of
## the outer and full parameter vectors and the Cholesky layout.
.foceiMeta <- function(design, start, fixed) {
  lay <- design$normal
  E <- design$etaNames[, lay$stems, drop = FALSE]
  outer <- names(start) %||% character(0)
  full <- c(outer, as.vector(E))
  probe <- setNames(numeric(length(full)), full)
  probe[outer] <- start
  groups <- design$groups
  # the frames of all conditions in one evaluation in C++ where the model allows it, laid out as the
  # frames the kernel reads
  conds <- unique(unlist(design$conditions[groups], use.names = FALSE))
  frames <- .chainFrames(design$chain, conds, lapply(conds, function(cn) .emGrid(design, cn)), list(probe), fixed,
                         FALSE, getOption("dMod2.cores", 1L))
  meta <- lapply(groups, function(g) {
    blocks <- lapply(design$conditions[[g]], function(cn) .conditionBlock(design, cn, probe, fixed, frames[[cn]]))
    used <- unique(unlist(lapply(blocks, `[[`, "parameters")))
    absent <- setdiff(E[g, ], used)
    if (length(absent))
      stop("em: the random effects ", paste(absent, collapse = ", "), " of group ", g,
           " enter none of the conditions of that group.", call. = FALSE)
    list(eta_idx_in_pars = match(E[g, ], full), eta_names = unname(E[g, ]),
         blocks = lapply(blocks, function(b) { b$parameters <- NULL; b }))
  })
  list(subjects = groups, conditions = names(design$data), K = lay$K,
       eta_idx_global = matrix(match(E, full), nrow(E)), eta_names = lapply(groups, function(g) unname(E[g, ])),
       outer_names = outer, outer_idx_in_full = match(outer, full),
       other_etas_init = setNames(rep(0, length(E)), as.vector(E)), pars_full_names = full,
       fast_meta = meta,
       omega_meta = list(chol_pars = lay$cholPars, chol_loc = matrix(as.integer(lay$loc), ncol = 2),
                         is_diag = unname(lay$isDiag)))
}


## The kernel treats censored rows by M3 only; without the kernel description `meta` the rows are
## read from the data.
.requireM3 <- function(design, meta, method) {
  censored <- if (is.null(meta))
    any(vapply(names(design$data), function(cn) { d <- .conditionData(design, cn); any(d$value <= d$lloq) }, TRUE))
  else any(unlist(lapply(meta$fast_meta, function(m)
    vapply(m$blocks, function(b) b$n_aloq < length(b$y_data), TRUE))))
  if (censored && !identical(design$bloq, "M3"))
    stop(method, " treats censored data by M3 only; build the objective with optBLOQ = \"M3\".",
         call. = FALSE)
}


## Batched predictor of the lock-step inner trust: one prediction and one error-model frame per
## request, conditions may repeat; the frame of the recorded sigmas where the error model does
## not cover a condition.
.batchPredictor <- function(design) {
  errfn <- .errorModel(design)
  function(timesList, parsList, conditions, fixed) {
    conditions <- as.character(conditions)
    # the native chain where it covers the conditions; the recorded sigmas where there is no error model
    fr <- .chainFrames(design$chain, conditions, timesList, parsList, fixed, FALSE, getOption("dMod2.cores", 1L))
    if (!is.null(fr)) {
      pred <- unname(lapply(fr, `[[`, "pred"))
      err <- unname(lapply(seq_along(fr), function(i) fr[[i]]$err %||%
        .dataSigmaErrfn(design$data)(out = fr[[i]]$pred, conditions = conditions[i])[[1L]]))
      return(list(pred = pred, err = err))
    }
    preds <- .predictMany(design$prdfn, times = timesList, parsList = parsList,
                          conditions = conditions, fixed = fixed, deriv = TRUE,
                          cores = getOption("dMod2.cores", 1L))
    n <- length(preds)
    inner <- lapply(preds, function(q) attr(q, "parameters"))
    node <- .fnNode(errfn)
    errs <- if (!is.null(node)) {
      .evalNode(node, .bundle(conds = conditions, out = as.list(preds), pars = inner,
                              fixed = vector("list", n)), TRUE, FALSE, NULL,
                getOption("dMod2.cores", 1L))
    } else lapply(seq_len(n), function(i)
      errfn(out = preds[[i]], pars = inner[[i]], conditions = conditions[i])[[1L]])
    errs <- as.list(errs)
    for (i in which(vapply(errs, is.null, TRUE)))
      errs[[i]] <- .dataSigmaErrfn(design$data)(out = preds[[i]], conditions = conditions[i])[[1L]]
    list(pred = unname(as.list(preds)), err = unname(errs))
  }
}


## The batched predictor `f` for the FOCEI kernel. Where the batch fails, every request is
## predicted alone, and one that fails gets frames the kernel takes as a failed trial of the
## subject; the message is kept in `failure`.
.guardPredictor <- function(f, failure = new.env(parent = emptyenv()))
  function(timesList, parsList, conditions, fixed) {
    failed <- function(e) {
      failure$message <- conditionMessage(e)
      list(pred = list("failed"), err = list("failed"))
    }
    tryCatch(f(timesList, parsList, conditions, fixed), error = function(e) {
      failed(e)
      one <- lapply(seq_along(conditions), function(i)
        tryCatch(f(timesList[i], parsList[i], conditions[i], fixed), error = failed))
      list(pred = lapply(one, function(x) x$pred[[1L]]), err = lapply(one, function(x) x$err[[1L]]))
    })
  }

## The prediction and error-model frames of a condition with first and second derivatives, the
## first only without `deriv2`; `err` is NULL for sigmas recorded in the data.
.volumeFrames <- function(design, cn, pred, fixed, deriv2 = TRUE) {
  if (deriv2 && (is.null(attr(pred, "deriv")) || is.null(attr(pred, "deriv2"))))
    stop("em: the prediction of condition ", cn, " has no second derivatives; build the ",
         "model with deriv2 = TRUE.", call. = FALSE)
  if (is.null(design$errfn)) return(list(pred = pred, err = NULL))
  inner <- getParameters(pred)
  fixedInner <- as.parvec(inner[attr(inner, "fixed")], deriv = FALSE, deriv2 = FALSE)
  inner <- as.parvec(inner[setdiff(names(inner), names(fixed))])
  err <- design$errfn(out = pred, pars = inner, fixed = fixedInner, deriv2 = deriv2, conditions = cn)[[cn]]
  list(pred = pred, err = err)
}

## The volume frames of the conditions `conds` with their predictions `preds`, the error model
## evaluated for all of them in one call.
.volumeFrameList <- function(design, conds, preds, fixed, cores = 1L, deriv2 = TRUE) {
  for (cn in conds)
    if (deriv2 && (is.null(attr(preds[[cn]], "deriv")) || is.null(attr(preds[[cn]], "deriv2"))))
      stop("em: the prediction of condition ", cn, " has no second derivatives; build the ",
           "model with deriv2 = TRUE.", call. = FALSE)
  node <- if (!is.null(design$errfn)) .fnNode(design$errfn)
  if (is.null(node))
    return(setNames(lapply(conds, function(cn) .volumeFrames(design, cn, preds[[cn]], fixed, deriv2)), conds))
  inner <- lapply(conds, function(cn) getParameters(preds[[cn]]))
  fixedInner <- lapply(inner, function(x) as.parvec(x[attr(x, "fixed")], deriv = FALSE, deriv2 = FALSE))
  inner <- lapply(inner, function(x) as.parvec(x[setdiff(names(x), names(fixed))]))
  errs <- .evalNode(node, .bundle(conds = conds, out = unname(preds[conds]), pars = inner,
                                  fixed = fixedInner), TRUE, deriv2, NULL, cores)
  setNames(lapply(seq_along(conds), function(i) list(pred = preds[[conds[i]]], err = errs[[i]])), conds)
}

## Parameter names the frames have derivatives in.
.frameParameters <- function(frames)
  unique(c(dimnames(attr(frames$pred, "deriv"))[[3]],
           if (!is.null(frames$err)) dimnames(attr(frames$err, "deriv"))[[3]]))

## arr[idx[r], a] or arr[idx[r], a, b] per row r of a derivative array whose time and
## observable dimensions are flattened into `idx`; names absent from `arr` give 0.
.rowSlices <- function(arr, idx, a, b = NULL) {
  out <- array(0, c(length(idx), length(a), if (!is.null(b)) length(b)))
  if (is.null(arr) || !length(a) || (!is.null(b) && !length(b))) return(out)
  dn <- dimnames(arr)
  flat <- array(arr, c(prod(dim(arr)[1:2]), dim(arr)[-(1:2)]))
  ia <- match(a, dn[[3]]); ka <- which(!is.na(ia))
  if (is.null(b)) {
    out[, ka] <- flat[idx, ia[ka], drop = FALSE]
    return(out)
  }
  ib <- match(b, dn[[4]]); kb <- which(!is.na(ib))
  out[, ka, kb] <- flat[idx, ia[ka], ib[kb], drop = FALSE]
  out
}

## Rows of one condition for the volume term: prediction p, sigma s, data y, censoring, and
## the first and second derivatives of p and s in the etas `eta` and the parameters `A`; a
## recorded sigma has no derivatives.
.volumeRows <- function(design, cn, frames, eta, A, hessian) {
  d <- .conditionData(design, cn)
  pred <- frames$pred
  tp <- .matchNum(d$time, as.numeric(pred[, "time"]))
  ip <- tp + (match(d$name, dimnames(attr(pred, "deriv"))[[2]]) - 1L) * nrow(pred)
  d1 <- attr(pred, "deriv"); d2 <- attr(pred, "deriv2")
  rows <- list(p = pred[cbind(tp, match(d$name, colnames(pred)))], s = d$sigma,
               y = pmax(d$value, d$lloq), cens = d$value <= d$lloq,
               a = .rowSlices(d1, ip, eta), pA = .rowSlices(d1, ip, A),
               Wp = .rowSlices(d2, ip, eta, A), Vp = .rowSlices(d2, ip, eta, eta),
               pAA = if (hessian) .rowSlices(d2, ip, A, A))
  n <- seq_len(nrow(d))
  rows$b <- .rowSlices(NULL, n, eta); rows$sA <- .rowSlices(NULL, n, A)
  rows$Ws <- .rowSlices(NULL, n, eta, A); rows$Vs <- .rowSlices(NULL, n, eta, eta)
  rows$sAA <- if (hessian) .rowSlices(NULL, n, A, A)
  err <- frames$err
  if (is.null(err)) return(rows)
  te <- .matchNum(d$time, as.numeric(err[, "time"]))
  model <- is.na(d$sigma)
  rows$s[model] <- err[cbind(te, match(d$name, colnames(err)))][model]
  e1 <- attr(err, "deriv"); e2 <- attr(err, "deriv2")
  if (is.null(e1)) return(rows)
  if (is.null(e2) && any(eta %in% dimnames(e1)[[3]]))
    stop("FOCEI with a random effect in the error model needs the second derivatives of the error model: ",
         "build the error model with deriv2 = TRUE.", call. = FALSE)
  # rows whose observable has no error-model derivative keep zero derivatives
  ie <- te + (match(d$name, dimnames(e1)[[2]]) - 1L) * nrow(err)
  none <- is.na(ie) | !model
  ie[none] <- 1L
  rows$b <- .rowSlices(e1, ie, eta); rows$sA <- .rowSlices(e1, ie, A)
  rows$Ws <- .rowSlices(e2, ie, eta, A); rows$Vs <- .rowSlices(e2, ie, eta, eta)
  rows$sAA <- if (hessian) .rowSlices(e2, ie, A, A)
  rows$b[none, ] <- 0; rows$sA[none, ] <- 0
  rows$Ws[none, , ] <- 0; rows$Vs[none, , ] <- 0
  if (hessian) rows$sAA[none, , ] <- 0
  rows
}

## Rows of several conditions stacked, arrays bound along their first dimension.
.stackRows <- function(parts) {
  out <- parts[[1]]
  if (length(parts) == 1L) return(out)
  for (nm in names(out)) {
    x <- lapply(parts, `[[`, nm)
    if (is.null(x[[1]])) next
    if (is.null(dim(x[[1]]))) { out[[nm]] <- unlist(x); next }
    dims <- dim(x[[1]])
    M <- do.call(rbind, lapply(x, function(a) matrix(a, dim(a)[1])))
    out[[nm]] <- array(M, c(nrow(M), dims[-1]))
  }
  out
}


## The units of the volume term, one per group, or per connected component with several
## groupings: the unit conditions and etas, the prior precision P = 2 Omega^-1 of these etas, the
## derivatives dP by the Cholesky parameters and d2P(a, b) the second ones.
.volumeUnits <- function(design, full, fixed) {
  subs <- design$groupings %||% list(design)
  derivs <- lapply(subs, function(sd) {
    lay <- sd$normal
    cd <- .choleskyDerivatives(lay, .choleskyFactor(lay, c(full, fixed)))
    list(lay = lay, P = 2 * cd$Oinv, d1 = setNames(lapply(cd$dOinv, `*`, 2), lay$cholPars),
         d2 = 2 * cd$d2Oinv)
  })
  # per grouping the groups of every condition, which may be several
  groupsOf <- lapply(subs, function(sd) split(rep(names(sd$conditions), lengths(sd$conditions)),
                                              unlist(sd$conditions, use.names = FALSE)))
  pieces <- function(cond) unlist(lapply(seq_along(subs), function(k) {
    sd <- subs[[k]]
    gs <- names(sd$conditions)[names(sd$conditions) %in% unlist(groupsOf[[k]][cond], use.names = FALSE)]
    lapply(gs, function(g) list(k = k, eta = unname(sd$etaNames[g, derivs[[k]]$lay$stems])))
  }), recursive = FALSE)
  units <- if (is.null(design$groupings)) lapply(design$groups, function(g) list(conditions = design$conditions[[g]]))
           else design$components
  lapply(units, function(u) {
    pc <- pieces(u$conditions)
    eta <- unlist(lapply(pc, `[[`, "eta"))
    K <- length(eta)
    at <- split(seq_len(K), rep(seq_along(pc), lengths(lapply(pc, `[[`, "eta"))))
    names <- unique(unlist(lapply(derivs, function(d) d$lay$cholPars)))
    P <- matrix(0, K, K)
    for (j in seq_along(pc)) P[at[[j]], at[[j]]] <- derivs[[pc[[j]]$k]]$P
    dP <- setNames(lapply(names, function(cp) {
      M <- matrix(0, K, K)
      for (j in seq_along(pc)) {
        d <- derivs[[pc[[j]]$k]]$d1[[cp]]
        if (!is.null(d)) M[at[[j]], at[[j]]] <- d
      }
      M
    }), names)
    d2P <- function(a, b) {
      M <- matrix(0, K, K)
      for (j in seq_along(pc)) {
        dk <- derivs[[pc[[j]]$k]]
        ia <- match(a, dk$lay$cholPars); ib <- match(b, dk$lay$cholPars)
        if (!is.na(ia) && !is.na(ib)) M[at[[j]], at[[j]]] <- dk$d2[ia, ib, , ]
      }
      M
    }
    list(conditions = u$conditions, eta = eta, P = P, dP = dP, d2P = d2P, pieces = pc)
  })
}

## The exact blocks H_ee and H_eA of the joint -2 log density of unit `u` in the unit etas and the
## parameters `A`, from the unit rows `rw`, their loss derivatives `lr` and the normal prior.
.jointBlocks <- function(rw, lr, u, full, A, cholIn) {
  tn <- length(rw$p); Kn <- length(u$eta)
  l_ <- function(i, j) lr[, i + 1L, j + 1L]
  Lpp <- l_(2, 0); Lps <- l_(1, 1); Lss <- l_(0, 2)
  cross <- function(xp, xs, yp, ys)
    crossprod(xp, Lpp * yp) + crossprod(xp, Lps * ys) + crossprod(xs, Lps * yp) + crossprod(xs, Lss * ys)
  curv <- function(Wp, Ws, n)
    matrix(colSums(l_(1, 0) * matrix(Wp, tn)) + colSums(l_(0, 1) * matrix(Ws, tn)), Kn, n)
  Hee <- cross(rw$a, rw$b, rw$a, rw$b) + curv(rw$Vp, rw$Vs, Kn) + u$P
  HeA <- cross(rw$a, rw$b, rw$pA, rw$sA) + curv(rw$Wp, rw$Ws, length(A))
  for (cp in intersect(cholIn, names(u$dP))) HeA[, match(cp, A)] <- HeA[, match(cp, A)] + u$dP[[cp]] %*% full[u$eta]
  dimnames(Hee) <- list(u$eta, u$eta); dimnames(HeA) <- list(u$eta, A)
  list(Hee = Hee, HeA = HeA)
}

## Gradient and Gauss-Newton Hessian of the volume term sum_u log|H_u(theta, eta_u*(theta))|, the
## mode shift d eta*/d theta = -H_ee^-1 H_eA from the exact joint blocks of the rows (`etaShift`
## per unit); the Hessian leaves out third-order sensitivities and d2 eta*/d theta2. With
## hessian = "sparse" the Hessian is a sparse Matrix and the shift of a unit has the unit columns A only.
.foceiVolume <- function(design, full, fixed, outer, H_inv_list, hessian = FALSE, cores = 1L) {
  sparse <- identical(hessian, "sparse")
  hessian <- sparse || isTRUE(hessian)
  Q <- length(outer)
  grad <- setNames(numeric(Q), outer)
  hess <- if (hessian && !sparse) matrix(0, Q, Q, dimnames = list(outer, outer))
  blocks <- list()
  units <- .volumeUnits(design, full, fixed)
  etaShift <- vector("list", length(units))
  cholIn <- intersect(outer, .hyperNames(design))

  conds <- unique(unlist(lapply(units, `[[`, "conditions"), use.names = FALSE))
  times <- lapply(conds, function(cn) .emGrid(design, cn))
  frames <- .chainFrames(design$chain, conds, times, list(full), fixed, TRUE, cores)
  if (!is.null(frames) && any(vapply(frames, function(f) is.null(f$err), TRUE))) frames <- NULL
  if (is.null(frames)) {
    preds <- .predictMany(design$prdfn, times = times, parsList = rep(list(full), length(conds)),
                          conditions = conds, fixed = fixed, deriv = TRUE, deriv2 = TRUE, cores = cores)
    names(preds) <- conds
    frames <- .volumeFrameList(design, conds, preds, fixed, cores)
  }

  # the rows of every unit in the outer parameters it depends on, directly, through the error
  # model or the prior, and the row loss of all of them in one call
  rows <- lapply(units, function(u) {
    fr <- frames[u$conditions]
    A <- unique(c(intersect(outer, unique(unlist(lapply(fr, .frameParameters)))), cholIn))
    list(A = A, rw = .stackRows(Map(function(cn, f) .volumeRows(design, cn, f, u$eta, A, hessian),
                                    u$conditions, fr)))
  })
  stacked <- function(nm) unlist(lapply(rows, function(r) r$rw[[nm]]))
  loss <- .rowLossDerivs(stacked("p"), stacked("s"), stacked("y"), stacked("cens"), design$df)
  unitOf <- rep(seq_along(rows), vapply(rows, function(r) length(r$rw$p), 0L))

  for (i in seq_along(units)) {
    u <- units[[i]]
    eta <- u$eta
    Kn <- length(eta)
    dOinv <- u$dP
    d2Oinv <- u$d2P
    A <- rows[[i]]$A
    Qa <- length(A)
    rw <- rows[[i]]$rw
    tn <- length(rw$p)
    a <- rw$a; b <- rw$b

    # the row loss and the derivatives in (p, s)
    lr <- loss[unitOf == i, , , drop = FALSE]
    l_ <- function(i, j) lr[, i + 1L, j + 1L]
    Lpp <- l_(2, 0); Lps <- l_(1, 1); Lss <- l_(0, 2)

    # d eta*/d theta from the exact joint blocks at the modes
    jb <- .jointBlocks(rw, lr, u, full, A, cholIn)
    if (!all(is.finite(jb$HeA)) || !all(is.finite(jb$Hee)))
      stop("FOCEI: the second derivatives of the model are not finite at these parameters.", call. = FALSE)
    thj <- colSums(abs(jb$HeA)) > 0
    Em <- matrix(0, Kn, Qa, dimnames = list(eta, A))
    if (any(thj)) {
      eg <- eigen(jb$Hee, symmetric = TRUE)
      lam <- pmax(eg$values, 1e-10 * abs(sum(eg$values)) / Kn)
      Em[, thj] <- -(eg$vectors %*% (t(eg$vectors) / lam)) %*% jb$HeA[, thj, drop = FALSE]
    }
    if (sparse) etaShift[[i]] <- Em
    else {
      etaShift[[i]] <- matrix(0, Kn, Q, dimnames = list(eta, outer))
      etaShift[[i]][, A] <- Em
    }

    # total first derivatives of p, s and of a, b
    Dp <- rw$pA + a %*% Em; Ds <- rw$sA + b %*% Em
    Da <- rw$Wp; Db <- rw$Ws
    for (l in seq_len(Kn)) {
      Da <- Da + array(outer(matrix(rw$Vp[, , l], tn, Kn), Em[l, ]), c(tn, Kn, Qa))
      Db <- Db + array(outer(matrix(rw$Vs[, , l], tn, Kn), Em[l, ]), c(tn, Kn, Qa))
    }

    DLpp <- l_(3, 0) * Dp + l_(2, 1) * Ds
    DLps <- l_(2, 1) * Dp + l_(1, 2) * Ds
    DLss <- l_(1, 2) * Dp + l_(0, 3) * Ds

    Hinv <- H_inv_list[[i]][eta, eta, drop = FALSE]
    aH <- a %*% Hinv; bH <- b %*% Hinv
    qaa <- rowSums(aH * a); qab <- rowSums(aH * b); qbb <- rowSums(bH * b)
    HDa <- array(0, c(tn, Kn, Qa)); HDb <- HDa
    for (al in seq_len(Qa)) {
      HDa[, , al] <- matrix(Da[, , al], tn, Kn) %*% Hinv
      HDb[, , al] <- matrix(Db[, , al], tn, Kn) %*% Hinv
    }
    # s_xy[r, alpha] = Dx_alpha' H^-1 y
    s_aa <- s_ab <- s_ba <- s_bb <- matrix(0, tn, Qa)
    for (al in seq_len(Qa)) {
      s_aa[, al] <- rowSums(matrix(HDa[, , al], tn, Kn) * a)
      s_ab[, al] <- rowSums(matrix(HDa[, , al], tn, Kn) * b)
      s_ba[, al] <- rowSums(matrix(HDb[, , al], tn, Kn) * a)
      s_bb[, al] <- rowSums(matrix(HDb[, , al], tn, Kn) * b)
    }

    gr <- colSums(2 * (Lpp * s_aa + Lps * (s_ab + s_ba) + Lss * s_bb) +
                  qaa * DLpp + 2 * qab * DLps + qbb * DLss)
    for (cp in intersect(cholIn, names(dOinv))) gr[match(cp, A)] <- gr[match(cp, A)] + sum(Hinv * dOinv[[cp]])
    grad[A] <- grad[A] + gr
    if (!hessian) next

    # D_a H, K x K per outer parameter
    V <- Lpp * a + Lps * b; Wv <- Lps * a + Lss * b
    DH <- array(0, c(Kn, Kn, Qa))
    for (al in seq_len(Qa)) {
      M <- crossprod(matrix(Da[, , al], tn, Kn), V) + crossprod(matrix(Db[, , al], tn, Kn), Wv)
      DH[, , al] <- M + t(M) + crossprod(a, DLpp[, al] * a) + crossprod(a, DLps[, al] * b) +
        crossprod(b, DLps[, al] * a) + crossprod(b, DLss[, al] * b)
    }
    for (cp in intersect(cholIn, names(dOinv))) DH[, , match(cp, A)] <- DH[, , match(cp, A)] + dOinv[[cp]]
    P <- array(0, c(Kn, Kn, Qa))
    for (al in seq_len(Qa)) P[, , al] <- Hinv %*% DH[, , al]
    Pm <- matrix(P, Kn * Kn)
    Pt <- matrix(aperm(P, c(2, 1, 3)), Kn * Kn)
    own <- intersect(cholIn, names(dOinv))

    # tr(H^-1 D_ab H) in Gauss-Newton form
    DaM <- matrix(Da, tn * Kn, Qa); DbM <- matrix(Db, tn * Kn, Qa)
    HDaM <- matrix(HDa, tn * Kn, Qa); HDbM <- matrix(HDb, tn * Kn, Qa)
    T1 <- crossprod(DaM, Lpp * HDaM) + crossprod(DaM, Lps * HDbM) +
      crossprod(DbM, Lps * HDaM) + crossprod(DbM, Lss * HDbM)
    X <- crossprod(s_aa, DLpp) + crossprod(s_ab + s_ba, DLps) + crossprod(s_bb, DLss)
    cpp <- qaa * l_(4, 0) + 2 * qab * l_(3, 1) + qbb * l_(2, 2)
    cps <- qaa * l_(3, 1) + 2 * qab * l_(2, 2) + qbb * l_(1, 3)
    css <- qaa * l_(2, 2) + 2 * qab * l_(1, 3) + qbb * l_(0, 4)
    T3 <- crossprod(Dp, cpp * Dp) + crossprod(Dp, cps * Ds) + crossprod(Ds, cps * Dp) + crossprod(Ds, css * Ds)
    # the third derivatives of the row loss times the second total derivatives of p and s
    c1 <- qaa * l_(3, 0) + 2 * qab * l_(2, 1) + qbb * l_(1, 2)
    c2 <- qaa * l_(2, 1) + 2 * qab * l_(1, 2) + qbb * l_(0, 3)
    T4 <- matrix(0, Qa, Qa)
    for (z in list(list(c1, rw$pAA, rw$Wp, rw$Vp), list(c2, rw$sAA, rw$Ws, rw$Vs))) {
      cz <- z[[1]]
      Mw <- matrix(colSums(cz * matrix(z[[3]], tn)), Kn, Qa)
      MV <- matrix(colSums(cz * matrix(z[[4]], tn)), Kn, Kn)
      T4 <- T4 + matrix(colSums(cz * matrix(z[[2]], tn)), Qa, Qa) +
        crossprod(Mw, Em) + crossprod(Em, Mw) + crossprod(Em, MV %*% Em)
    }
    vol <- 2 * T1 + 2 * (X + t(X)) + T3 + T4
    for (al in own) for (be in own)
      vol[match(al, A), match(be, A)] <- vol[match(al, A), match(be, A)] + sum(Hinv * d2Oinv(al, be))

    # - tr(H^-1 D_a H H^-1 D_b H)
    block <- vol - crossprod(Pt, Pm)
    if (sparse) blocks[[length(blocks) + 1L]] <- structure(block, dimnames = list(A, A))
    else hess[A, A] <- hess[A, A] + block
  }
  if (sparse) hess <- .sparseFromBlocks(blocks, outer)
  fx <- intersect(names(fixed), outer)
  grad[fx] <- 0
  if (hessian && length(fx)) hess[fx, ] <- hess[, fx] <- 0
  list(gradient = grad, hessian = hess, etaShift = etaShift)
}


## FOCEI fit of the normal random effects of `design` with the joint objective `joint` (data
## plus normal prior) from `start`. The conditional covariance of the etas of a group is
## (H_g / 2)^-1, H_g the inner Hessian of the -2 log density.
.fitFocei <- function(design, joint, start, fixed = NULL, control = list()) {
  if (.foceiSparseChoice(design, joint, start, fixed, control$sparse %||% "auto"))
    return(.fitFoceiSparse(design, joint, start, fixed, control))
  ic <- modifyList(list(rinit = 1, rmax = 10, iterlim = 30, steptol = 1e-8, mtol = 1e-12,
                        eigen_floor_relative = 1e-10), control$innerControl %||% list())
  oc <- modifyList(list(rinit = 1, rmax = 10, iterlim = 200, ftol = 1e-7, mtol = 1e-5),
                   control$trustControl %||% list())
  cores <- getOption("dMod2.cores", 1L)
  meta <- .foceiMeta(design, start, fixed)
  .requireM3(design, meta, "FOCEI")
  outer <- names(start)
  mode <- .secondOrderMode(design, joint, start, fixed, control)
  correction <- mode != "none"
  # the outer model: the Gauss-Newton Schur complement, forward over forward plus the volume Hessian,
  # plus the SR1 correction from the full outer gradients of the accepted points; with the volume
  # gradient, and d eta*/d theta, which starts the modes of the next trial to first order
  learn <- .outerSecant(outer)
  last <- NULL
  volume <- function(full_pars, joint_hessian, H_inv_list, data_gradient) {
    ff <- mode == "forward-forward"
    vol <- if (ff) .foceiVolume(design, full_pars, fixed, outer, H_inv_list, hessian = TRUE, cores = cores)
           else .volumeReverse(design, full_pars, fixed, outer, H_inv_list, cores = cores)
    curvature <- if (ff) vol$hessian[outer, outer] else 0
    if (ff) last <<- list(full = full_pars, hessian = vol$hessian)
    base <- .gnSchur(joint_hessian, outer, meta$eta_names, H_inv_list) + curvature
    learn$add(full_pars[outer], data_gradient[outer] + vol$gradient[outer], function(x) base %*% x)
    list(gradient = vol$gradient, hessian = curvature + learn$dense(), eta_shift = vol$etaShift)
  }
  failure <- new.env(parent = emptyenv())
  fit <- tryCatch(focei_run(model_cb = design$prdfn, err_cb = .errorModel(design), joint_cb = joint,
                            init = start, subject_meta = meta, fixed = fixed,
                            control = list(inner = ic, outer = oc),
                            correction_mode = if (correction) "eager" else "none",
                            correction_cb = if (correction) volume else NULL,
                            predict_cb = .guardPredictor(.batchPredictor(design), failure)),
                  error = function(e) stop(conditionMessage(e), if (!is.null(failure$message))
                    paste0("; the last prediction that failed: ", failure$message), call. = FALSE))

  E <- design$etaNames[, design$normal$stems, drop = FALSE]
  etaModes <- fit$etaModes
  dimnames(etaModes) <- list(design$groups, design$normal$stems)
  hessian <- fit$hessian
  full <- c(fit$argument, setNames(as.vector(etaModes), as.vector(E)))
  exact <- mode == "forward-forward" || (mode == "forward-reverse" && .jointDeriv2(joint, full, fixed))
  if (mode == "forward-reverse" && !exact)
    hessian <- .reportedReverse(design, joint, full, fixed, outer, fit$H_inv, cores)
  if (exact) {
    # exact Schur complement of the joint Hessian at the modes plus the volume Hessian
    jh <- joint(pars = full, fixed = fixed, deriv = TRUE, deriv2 = TRUE)$hessian
    on <- rownames(hessian)
    ee <- intersect(as.vector(E), rownames(jh))
    hessian <- jh[on, on, drop = FALSE] -
      crossprod(jh[ee, on, drop = FALSE], .scaledSolve(jh[ee, ee, drop = FALSE], jh[ee, on, drop = FALSE]))
    # the volume Hessian of the last correction where the fit ends on the correction point
    vol <- if (!is.null(last) && identical(names(last$full), names(full)) &&
               identical(unname(unclass(last$full)), unname(full))) last$hessian
           else .foceiVolume(design, full, fixed, outer, fit$H_inv, hessian = TRUE, cores = cores)$hessian
    hessian <- hessian + vol[on, on]
  }
  list(argument = fit$argument, value = fit$value, gradient = fit$gradient, hessian = hessian,
       etaModes = etaModes, etaCov = setNames(lapply(fit$H_inv, function(H) 2 * H), design$groups),
       converged = isTRUE(fit$converged), iterations = fit$iterations, trace = fit$trace,
       derivMode = mode, hessianExact = exact)
}

## How FOCEI gets the second derivatives: `control$derivMode`, forward over reverse by default;
## "none" with control$secondOrderCorrection = FALSE. Where the chain lacks the mode asked for,
## the other one, with a message where it was asked for explicitly; an error where it has neither.
.secondOrderMode <- function(design, joint, start, fixed, control, units = NULL) {
  if (isFALSE(control$secondOrderCorrection)) return("none")
  etas <- .allEtas(design)
  full <- c(start, setNames(numeric(length(etas)), etas))
  # without a probe where the native chain covers every condition: forward over forward with
  # second derivatives, forward over reverse with sweeps
  want <- control$derivMode %||% "forward-reverse"
  native <- .evalCpp() && !is.null(design$chain) && {
    conds <- names(design$data)
    ff <- want == "forward-forward"
    !is.null(.chainLookup(design$chain, conds, rep(list(names(full)), length(conds)), names(fixed), if (!ff) names(full),
                          TRUE, ff, lapply(conds, function(cn) .emGrid(design, cn)),
                          tangents = if (ff) 0L else 1L, seeds = if (ff) 0L else 1L))
  }
  if (native) return(want)
  units <- units %||% .volumeUnits(design, full, fixed)
  reverse <- function() !inherits(tryCatch(
    .sweepUnits(design, full, fixed, units, lapply(units, function(u)
      list(matrix(diag(length(u$eta)), length(u$eta), dimnames = list(u$eta, u$eta))))),
    error = function(e) e), "error")
  forward <- function() !inherits(tryCatch(joint(pars = full, fixed = fixed, deriv = TRUE, deriv2 = TRUE),
                                           error = function(e) e), "error")
  other <- setdiff(c("forward-reverse", "forward-forward"), want)
  if (if (want == "forward-reverse") reverse() else forward()) return(want)
  if (if (want == "forward-reverse") forward() else reverse()) {
    .derivModeMessage(if (!is.null(control$derivMode))
      paste0("FOCEI: the chain has no ", want, " second derivatives; it takes ", other, ".")
      else paste0("em: second derivatives forward over forward; models built with derivMode = ",
                  "c(\"forward\", \"reverse\", \"forward-reverse\") give them forward over reverse, which is faster."))
    return(other)
  }
  stop("FOCEI needs second derivatives of the prediction and the error model: build the models ",
       "with derivMode = c(\"forward\", \"reverse\", \"forward-reverse\") or with deriv2 = TRUE, ",
       "or set control = list(secondOrderCorrection = FALSE), which biases the estimates.", call. = FALSE)
}

## A message on the derivative mode of a fit, of a class that msem() collects to give it once.
.derivModeMessage <- function(text) {
  cond <- simpleMessage(paste0(text, "\n"))
  class(cond) <- c("dMod2DerivModeMessage", class(cond))
  message(cond)
}

## The Gauss-Newton Schur complement of the outer parameters in the joint Hessian `H` over the
## etas of every group, with their inverse blocks `H_inv_list` from the inner problem.
.gnSchur <- function(H, outer, etaNames, H_inv_list) {
  S <- H[outer, outer, drop = FALSE]
  for (i in seq_along(H_inv_list)) {
    e <- etaNames[[i]]
    Hoe <- H[outer, e, drop = FALSE]
    S <- S - Hoe %*% H_inv_list[[i]][e, e, drop = FALSE] %*% t(Hoe)
  }
  S
}

## The SR1 correction of an outer model from the steps s between the points `add()` is called at
## and the curvature the model at the newer point misses along them, y - B s with y the change of
## the full outer gradient and `times(x)` = B x; the last `memory` pairs from zero in their order,
## as U diag(d) U' (`factors()`) or a matrix (`dense()`). A pair with |r's| <= 1e-8 |r| |s| is skipped.
.outerSecant <- function(outer, memory = 2L) {
  S <- M <- matrix(0, length(outer), 0)
  last <- NULL
  keep <- function(A, a) cbind(A, a)[, max(1, ncol(A) + 2 - memory):(ncol(A) + 1), drop = FALSE]
  add <- function(theta, gradient, times) {
    theta <- theta[outer]; gradient <- gradient[outer]
    if (!is.null(last) && any(theta != last$theta)) {
      s <- theta - last$theta
      S <<- keep(S, s)
      M <<- keep(M, gradient - last$gradient - as.numeric(times(s)))
    }
    last <<- list(theta = theta, gradient = gradient)
  }
  factors <- function() {
    U <- matrix(0, length(outer), 0); d <- numeric(0)
    for (k in seq_len(ncol(S))) {
      s <- S[, k]
      r <- M[, k] - drop(U %*% (d * crossprod(U, s)))
      rs <- sum(r * s)
      if (!is.finite(rs) || abs(rs) <= 1e-8 * sqrt(sum(r^2) * sum(s^2))) next
      U <- cbind(U, r); d <- c(d, 1 / rs)
    }
    list(U = U, d = d)
  }
  dense <- function() {
    f <- factors()
    matrix(f$U %*% (f$d * t(f$U)), length(outer), length(outer), dimnames = list(outer, outer))
  }
  list(add = add, factors = factors, dense = dense)
}

## Whether the chain has the second-order sensitivities of `deriv2 = TRUE`.
.jointDeriv2 <- function(joint, full, fixed)
  !inherits(tryCatch(joint(pars = full, fixed = fixed, deriv = TRUE, deriv2 = TRUE), error = function(e) e), "error")

## The Hessian a FOCEI fit by forward over reverse reports without second-order sensitivities: the
## exact Schur complement of the joint Hessian up to 500 outer parameters, the Gauss-Newton one
## above, plus the volume curvature from dH/dtheta (`.volumeReverse(curvature = TRUE)`); `sparse`
## a sparse Matrix.
.reportedReverse <- function(design, joint, full, fixed, outer, H_inv, cores = 1L, sparse = FALSE) {
  vol <- .volumeReverse(design, full, fixed, outer, H_inv, cores = cores, local = sparse, curvature = TRUE)$hessian
  if (length(outer) <= 500L) {
    S <- .schurReverse(design, full, fixed, outer, cores)
    return(if (sparse) Matrix::Matrix(S, sparse = TRUE) + vol else S + vol)
  }
  units <- .volumeUnits(design, full, fixed)
  H <- joint(pars = full, fixed = fixed, deriv = TRUE, hessian = if (sparse) "sparse")$hessian
  S <- H[outer, outer, drop = FALSE]
  for (i in seq_along(units)) {
    e <- units[[i]]$eta
    Hoe <- H[outer, e, drop = FALSE]
    S <- S - Hoe %*% H_inv[[i]][e, e, drop = FALSE] %*% Matrix::t(Hoe)
  }
  S + vol
}


## Whether FOCEI takes the sparse outer path: `sparse` TRUE, or "auto" with at least the
## parameters of getOption("dMod2.sparseAuto") and an outer Hessian, the joint Gauss-Newton
## Hessian with the couplings of the Schur complements, that fills at most `fill`.
.foceiSparseChoice <- function(design, joint, start, fixed, sparse = "auto") {
  if (!identical(sparse, "auto")) return(isTRUE(sparse))
  th <- .sparseThresholds()
  outer <- names(start)
  if (length(outer) < th$parameters) return(FALSE)
  etas <- .allEtas(design)
  full <- c(start, setNames(numeric(length(etas)), etas))
  jt <- tryCatch(joint(pars = full, fixed = fixed, deriv = TRUE, hessian = "sparse"), error = function(e) NULL)
  if (is.null(jt) || !.isSparse(jt$hessian)) return(FALSE)
  H <- .sparseOver(jt$hessian, c(outer, etas)) != 0
  B <- H[outer, etas, drop = FALSE]
  pattern <- H[outer, outer, drop = FALSE] | (B %*% Matrix::t(B)) != 0
  Matrix::nnzero(pattern) / length(outer)^2 <= th$fill
}

## FOCEI of one grouping with the outer trust region on sparse Hessians: the modes by the inner
## trust of the kernel, the outer Hessian the sparse joint Gauss-Newton Hessian minus the Schur
## complement of every group plus the volume Hessian, each local to the parameters of the group.
## The joint objective of a trial point comes with derivatives, so the value comes from the same
## solve as on the dense path; the volume term is evaluated where a step is accepted.
.fitFoceiSparse <- function(design, joint, start, fixed = NULL, control = list()) {
  ic <- modifyList(list(rinit = 1, rmax = 10, iterlim = 30, steptol = 1e-8, mtol = 1e-12,
                        eigen_floor_relative = 1e-10), control$innerControl %||% list())
  oc <- modifyList(list(rinit = 1, rmax = 10, iterlim = 200, ftol = 1e-7, mtol = 1e-5),
                   control$trustControl %||% list())
  cores <- getOption("dMod2.cores", 1L)
  meta <- .foceiMeta(design, start, fixed)
  .requireM3(design, meta, "FOCEI")
  lay <- design$normal
  outer <- names(start)
  E <- design$etaNames[, lay$stems, drop = FALSE]
  etas <- as.vector(E)
  byGroup <- as.vector(t(E))
  N <- nrow(E); K <- lay$K
  predictor <- .guardPredictor(.batchPredictor(design))
  errfn <- .errorModel(design)
  mode <- .secondOrderMode(design, joint, start, fixed, control)
  correction <- mode != "none"

  # the modes and the value at theta, the inner trust started from `warm`
  modesAt <- function(theta, warm) {
    L <- .choleskyFactor(lay, c(theta, fixed))
    full <- setNames(numeric(length(meta$pars_full_names)), meta$pars_full_names)
    full[outer] <- theta
    inner <- focei_inner_trust(design$prdfn, errfn, full, warm, meta, chol2inv(t(L)),
                               2 * sum(log(diag(L))), fixed, control = ic, predict_cb = predictor)
    full[etas] <- as.vector(inner$eta_modes)
    jt <- joint(pars = full, fixed = fixed, deriv = TRUE, hessian = "sparse", conditions = meta$conditions)
    list(theta = theta, full = full, inner = inner, joint = jt,
         value = jt$value + sum(inner$log_det_H) - N * K * log(2))
  }
  # the gradient and the sparse Hessian of the outer objective, and d eta*/d theta per group
  derivativesAt <- function(ev, deriv2 = FALSE) {
    jt <- if (deriv2) joint(pars = ev$full, fixed = fixed, deriv = TRUE, deriv2 = TRUE, hessian = "sparse",
                            conditions = meta$conditions)
          else ev$joint
    H <- .sparseOver(jt$hessian, c(outer, etas))
    Hoe <- H[outer, byGroup, drop = FALSE]
    inv <- if (deriv2) Matrix::solve(H[byGroup, byGroup], Matrix::t(Hoe))
           else Matrix::bdiag(lapply(ev$inner$H_inv, unname)) %*% Matrix::t(Hoe)
    ev$hessian <- Matrix::forceSymmetric(H[outer, outer] - Hoe %*% inv)
    # the gradient at the exact modes to first order: the eta gradient the inner trust leaves
    # enters through the coupling
    left <- Matrix::bdiag(lapply(ev$inner$H_inv, unname)) %*% jt$gradient[byGroup]
    ev$gradient <- jt$gradient[outer] - as.numeric(Hoe %*% left)
    if (!correction) return(ev)
    # forward over reverse: the exact volume gradient on the Gauss-Newton Schur model, whose
    # curvature the SR1 correction learns
    if (mode == "forward-reverse" && !deriv2) {
      vol <- .volumeReverse(design, ev$full, fixed, outer, ev$inner$H_inv, cores = cores, local = TRUE)
      ev$gradient <- ev$gradient + vol$gradient[outer]
      ev$shift <- vol$etaShift
      return(ev)
    }
    vol <- .foceiVolume(design, ev$full, fixed, outer, ev$inner$H_inv, hessian = "sparse", cores = cores)
    ev$gradient <- ev$gradient + vol$gradient[outer]
    ev$hessian <- Matrix::forceSymmetric(ev$hessian + vol$hessian)
    ev$shift <- vol$etaShift
    ev
  }
  # the modes of `ev` moved to first order along `step`, for a group whose predicted move stays within
  # the region the inner search starts with, as in the dense kernel
  shifted <- function(ev, step) {
    w <- ev$inner$eta_modes
    for (i in seq_along(ev$shift)) {
      S <- ev$shift[[i]]
      if (!length(S)) next
      d <- drop(S %*% step[colnames(S)])
      if (is.finite(sum(d^2)) && sqrt(sum(d^2)) <= ic$rinit) w[i, ] <- w[i, ] + d
    }
    w
  }

  # the SR1 correction of the outer model from the full outer gradients of the accepted points
  learn <- .outerSecant(outer)
  cur <- derivativesAt(modesAt(start, matrix(0, N, K)))
  if (correction) learn$add(cur$theta, cur$gradient, function(x) cur$hessian %*% x)
  r <- oc$rinit
  converged <- FALSE
  it <- stall <- 0L
  stopReason <- "iterlim"
  trace <- list(value = cur$value, radius = r, rho = 0, accept = 1L)
  while (it < oc$iterlim) {
    it <- it + 1L
    lr <- if (correction) learn$factors()
    st <- .sparseTrustStep(cur$hessian, cur$gradient, r, lowrank = lr)
    step <- setNames(st$p, outer)
    Bp <- as.numeric(cur$hessian %*% st$p)
    if (length(lr$d)) Bp <- Bp + drop(lr$U %*% (lr$d * crossprod(lr$U, st$p)))
    pred <- -(sum(cur$gradient * st$p) + 0.5 * sum(st$p * Bp))
    trial <- tryCatch(modesAt(cur$theta + step, shifted(cur, step)), error = function(e) NULL)
    ft <- if (is.null(trial)) Inf else trial$value
    rho <- if (pred > 0 && is.finite(ft)) (cur$value - ft) / pred else -Inf
    # what the model can gain at all, in the largest region the fit allows, as in the dense kernel
    reach <- if (r < oc$rmax) {
      sm <- .sparseTrustStep(cur$hessian, cur$gradient, oc$rmax, lowrank = lr)
      Bm <- as.numeric(cur$hessian %*% sm$p)
      if (length(lr$d)) Bm <- Bm + drop(lr$U %*% (lr$d * crossprod(lr$U, sm$p)))
      -(sum(cur$gradient * sm$p) + 0.5 * sum(sm$p * Bm))
    } else pred
    done <- is.finite(ft) && abs(reach) < oc$mtol && (abs(cur$value - ft) < oc$ftol || abs(pred) < oc$mtol)
    ok <- if (done) ft < cur$value else rho >= 0.25
    if (!done && rho < 0.25) r <- min(r, sqrt(sum(st$p^2))) / 4
    else if (!done && rho > 0.75 && st$boundary) r <- min(2 * r, oc$rmax)
    if (ok) {
      nxt <- tryCatch(derivativesAt(trial), error = function(e) NULL)
      if (is.null(nxt)) r <- r / 4 else cur <- nxt
      if (!is.null(nxt) && correction) learn$add(cur$theta, cur$gradient, function(x) cur$hessian %*% x)
    }
    for (f in names(trace)) trace[[f]] <- c(trace[[f]], switch(f, value = cur$value, radius = r, rho = rho,
                                                               accept = as.integer(ok)))
    # as trust(): five consecutive rejected trials that leave the objective flat within ftol, the noise
    # of the evaluation, end the fit with nothing left to resolve
    stall <- if (ok || !is.finite(ft) || abs(cur$value - ft) >= oc$ftol) 0L else stall + 1L
    if (done || stall >= 5L) { converged <- TRUE; stopReason <- if (done) "model" else "stall"; break }
  }
  # the exact Schur complement of the joint Hessian at the modes plus the volume Hessian; by
  # forward over reverse with the volume curvature from dH/dtheta where there is no deriv2
  exact <- mode == "forward-forward" || (mode == "forward-reverse" && .jointDeriv2(joint, cur$full, fixed))
  final <- if (exact) derivativesAt(cur, deriv2 = TRUE) else cur
  if (mode == "forward-reverse" && !exact)
    final$hessian <- .reportedReverse(design, joint, cur$full, fixed, outer, cur$inner$H_inv, cores, sparse = TRUE)
  etaModes <- cur$inner$eta_modes
  dimnames(etaModes) <- list(design$groups, lay$stems)
  list(argument = cur$theta, value = cur$value, gradient = cur$gradient, hessian = final$hessian,
       etaModes = etaModes, etaCov = setNames(lapply(cur$inner$H_inv, function(H) 2 * H), design$groups),
       converged = converged, iterations = it, stopReason = stopReason, trace = trace, derivMode = mode,
       hessianExact = exact)
}


## FOCEI over the connected components of several groupings: per component the joint mode by
## trust, the Laplace approximation, the envelope gradient with the exact volume derivative and
## the exact Schur complement as the Hessian.
.fitFoceiJoint <- function(design, joint, start, fixed = NULL, control = list()) {
  .requireM3(design, NULL, "FOCEI")
  sparse <- .foceiSparseChoice(design, joint, start, fixed, control$sparse %||% "auto")
  curvature <- if (sparse) "sparse" else NULL
  units <- .volumeUnits(design, c(start, setNames(numeric(length(.allEtas(design))), .allEtas(design))), fixed)
  etas <- .allEtas(design)
  outer <- names(start)
  mode <- .secondOrderMode(design, joint, start, fixed, control)
  learn <- .outerSecant(outer)
  # a component has the modes where the component Newton step is below `steptol`: the volume term is linear in
  # their error, so a small change of the value does not tell they are found
  ic <- modifyList(list(rinit = 1, rmax = 10, iterlim = 100L, mtol = 1e-12, steptol = 1e-8),
                   control$innerControl %||% list())
  warm <- setNames(numeric(length(etas)), etas)
  # the objectives per condition of the R evaluation, built where the chain does not cover the model
  delayedAssign("objs", .conditionObjectives(design))
  # The modes of all components in lock step (focei_units_lockstep()): per round one batched
  # evaluation of the conditions of the components still moving, then per component a trust step on
  # the component data and the normal prior of the component etas, e' P e / 2, accepted by a separate ratio.
  modes <- function(theta) {
    held <- c(theta, fixed)
    P <- lapply(.volumeUnits(design, c(theta, fixed, warm), fixed), function(u) {
      m <- as.matrix(u$P); storage.mode(m) <- "double"; m
    })
    # the data terms of the components `idx` at the etas `etas`, summed per component in C++ where
    # the model allows it; NULL for a component whose prediction fails, the others evaluated apart
    dataUnits <- function(idx, etas) {
      cs <- unlist(lapply(units[idx], `[[`, "conditions"), use.names = FALSE)
      owner <- rep(seq_along(idx), vapply(units[idx], function(u) length(u$conditions), 0L))
      parsL <- lapply(owner, function(j) setNames(etas[[j]], units[[idx[j]]]$eta))
      loc <- tryCatch(.chainLocal(design$chain, cs, parsL, rep(list(held), length(cs)), TRUE, FALSE,
                                  getOption("dMod2.cores", 1L), design$bloq, design$df), error = function(e) NULL)
      d <- if (!is.null(loc)) chain_unit_sums(loc, owner - 1L, lapply(units[idx], `[[`, "eta"))
      else tryCatch({
        m <- .objEvalMany(objs[cs], parsL, deriv = TRUE, fixed = held)
        lapply(seq_along(idx), function(j) {
          eta <- units[[idx[j]]]$eta
          mine <- m[owner == j]
          list(value = sum(vapply(mine, `[[`, 0, "value")),
               gradient = Reduce(`+`, lapply(mine, function(x) .named(x$gradient, eta))),
               hessian = Reduce(`+`, lapply(mine, function(x) .namedHessian(x$hessian, eta))))
        })
      }, error = function(e) NULL)
      if (is.null(d)) return(if (length(idx) == 1L) list(NULL)
                             else lapply(seq_along(idx), function(j) dataUnits(idx[j], etas[j])[[1L]]))
      lapply(d, function(x) if (is.finite(x$value)) x)
    }
    w <- focei_units_lockstep(dataUnits, lapply(units, function(u) unname(warm[u$eta])), P, ic)
    for (i in seq_along(units)) warm[units[[i]]$eta] <<- w[[i]]
    warm
  }
  outerObjective <- function(theta, deriv2 = FALSE, volumeHessian = FALSE, ...) {
    theta <- theta[outer]
    full <- c(theta, modes(theta))
    gn <- joint(pars = full, fixed = fixed, deriv = TRUE, hessian = curvature)
    block <- function(u) as.matrix(gn$hessian[u$eta, u$eta, drop = FALSE])
    Hinv <- lapply(units, function(u) .scaledSolve(block(u)))
    logdet <- sum(vapply(units, function(u) as.numeric(determinant(block(u) / 2)$modulus), 0))
    cores <- getOption("dMod2.cores", 1L)
    schur <- function(h) {
      ee <- intersect(etas, rownames(h))
      h[outer, outer, drop = FALSE] -
        Matrix::crossprod(h[ee, outer, drop = FALSE], Matrix::solve(h[ee, ee, drop = FALSE], h[ee, outer, drop = FALSE]))
    }
    if (mode == "forward-forward") {
      jh <- joint(pars = full, fixed = fixed, deriv = TRUE, deriv2 = TRUE, hessian = curvature)$hessian
      vol <- .foceiVolume(design, full, fixed, outer, Hinv, hessian = if (!volumeHessian) FALSE else curvature %||% TRUE,
                          cores = cores)
      H <- schur(jh)
      if (volumeHessian) H <- H + vol$hessian[outer, outer]
      volGradient <- vol$gradient[outer]
    } else {
      # forward over reverse: the exact volume gradient on the Gauss-Newton Schur model, whose
      # curvature the SR1 correction of the dense path learns; the exact Schur complement where it
      # is reported
      vol <- if (mode != "none") .volumeReverse(design, full, fixed, outer, Hinv, cores = cores, local = sparse)
      volGradient <- if (mode == "none") 0 else vol$gradient[outer]
      H <- schur(gn$hessian)
      if (mode != "none") {
        learn$add(theta, gn$gradient[outer] + volGradient, function(x) H %*% x)
        if (!sparse) H <- H + learn$dense()
      }
      if (volumeHessian && mode != "none") {
        H <- if (.jointDeriv2(joint, full, fixed))
          schur(joint(pars = full, fixed = fixed, deriv = TRUE, deriv2 = TRUE, hessian = curvature)$hessian) +
            .foceiVolume(design, full, fixed, outer, Hinv, hessian = curvature %||% TRUE, cores = cores)$hessian[outer, outer]
        else .reportedReverse(design, joint, full, fixed, outer, Hinv, cores, sparse = sparse)
      }
    }
    if (sparse) H <- Matrix::forceSymmetric(Matrix::Matrix(H, sparse = TRUE))
    out <- objlist(value = gn$value + logdet, gradient = gn$gradient[outer] + volGradient, hessian = H)
    attr(out, "Hinv") <- Hinv
    out
  }
  oc <- modifyList(list(rinit = 1, rmax = 10, iterlim = 200L, ftol = 1e-7, mtol = 1e-5), control$trustControl %||% list())
  args <- list(objfun = outerObjective, parinit = start, rinit = oc$rinit, rmax = oc$rmax, iterlim = oc$iterlim,
               tolControl = oc[intersect(names(oc), c("ftol", "mtol", "gtol"))],
               sparse = sparse)
  # iterlim 0 evaluates the marginal at the start, which is no failure to converge
  fit <- if (oc$iterlim > 0) do.call(trust, args) else .innerTrust(args)
  final <- outerObjective(fit$argument, volumeHessian = TRUE)
  cov <- unlist(lapply(seq_along(units), function(i) {
    V <- 2 * attr(final, "Hinv")[[i]]
    dimnames(V) <- list(units[[i]]$eta, units[[i]]$eta)
    list(V)
  }), recursive = FALSE)
  list(argument = fit$argument[outer], value = final$value, gradient = final$gradient,
       hessian = final$hessian, etas = warm, etaCovUnits = cov, units = units,
       converged = isTRUE(fit$converged), iterations = fit$iterations, stopReason = fit$stopReason,
       derivMode = mode,
       hessianExact = mode == "forward-forward" || .jointDeriv2(joint, c(fit$argument[outer], warm), fixed))
}

## A gradient restricted to `names`, 0 where absent.
.named <- function(g, names) { out <- setNames(numeric(length(names)), names); h <- intersect(names, names(g)); out[h] <- g[h]; out }

.namedHessian <- function(H, names) {
  out <- matrix(0, length(names), length(names), dimnames = list(names, names))
  h <- intersect(names, rownames(H))
  if (length(h)) out[h, h] <- as.matrix(H[h, h, drop = FALSE])
  out
}
