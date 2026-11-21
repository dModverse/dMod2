# FOCEI: -2 log L ~ J(theta, eta*) + sum_g log|H_g / 2| around the conditional modes of every
# group, J the joint objective, H_g the Gauss-Newton Hessian of J in the etas of g. The gradient has the
# exact derivative of the volume term; the guide chapter Nonlinear mixed-effects models has the derivation.


## Derivatives up to order 4 of the row -2 log likelihood in prediction p and sigma s, normal or Student-t,
## M3 below the limit: [n, 5, 5] with [, i + 1, j + 1] = d^(i+j) / dp^i ds^j and the value left at 0; normal
## rows above the limit in closed form, the others from the Taylor series in (dp, ds).
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
    obs <- .setdiff(colnames(frame), "time")
    err <- matrix(NA_real_, length(times), length(obs) + 1L, dimnames = list(NULL, c("time", obs)))
    err[, "time"] <- times
    for (o in .intersect(obs, as.character(d$name))) {
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
         paste(.setdiff(d$name, colnames(pred)), collapse = ", "), call. = FALSE)
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
## the outer and full parameter vectors and the Cholesky layout. The etas are those of `stems`; with
## `reach` (groups x stems) a group has those it marks as coordinates, the others rest at 0.
.foceiMeta <- function(design, start, fixed, stems = design$normal$stems, reach = NULL) {
  lay <- design$normal
  E <- design$etaNames[, stems, drop = FALSE]
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
    on <- if (is.null(reach)) rep(TRUE, ncol(E)) else reach[g, stems]
    absent <- .setdiff(E[g, on], used)
    if (length(absent))
      stop("em: the random effects ", paste(absent, collapse = ", "), " of group ", g,
           " enter none of the conditions of that group.", call. = FALSE)
    out <- list(eta_idx_in_pars = match(E[g, on], full), eta_names = unname(E[g, on]),
                blocks = lapply(blocks, function(b) { b$parameters <- NULL; b }))
    if (!all(on)) out$coords <- which(on)
    out
  })
  list(subjects = groups, conditions = names(design$data), K = ncol(E),
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
  inner <- as.parvec(inner[.setdiff(names(inner), names(fixed))])
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
  inner <- lapply(inner, function(x) as.parvec(x[.setdiff(names(x), names(fixed))]))
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
  for (cp in .intersect(cholIn, names(u$dP))) HeA[, match(cp, A)] <- HeA[, match(cp, A)] + u$dP[[cp]] %*% full[u$eta]
  dimnames(Hee) <- list(u$eta, u$eta); dimnames(HeA) <- list(u$eta, A)
  list(Hee = Hee, HeA = HeA)
}

## Gradient and Hessian of the volume term sum_u log|H_u(theta, eta_u*(theta))| without third-order
## sensitivities and d2 eta*/d theta2, and the mode shift -H_ee^-1 H_eA per unit (`etaShift`); with
## hessian = "sparse" the row terms local to their condition and a wide unit of rank 2K + K^2 as `lowrank`.
.foceiVolume <- function(design, full, fixed, outer, H_inv_list, hessian = FALSE, cores = 1L) {
  sparse <- identical(hessian, "sparse")
  hessian <- sparse || isTRUE(hessian)
  Q <- length(outer)
  grad <- setNames(numeric(Q), outer)
  blocks <- list()
  lowU <- list(); lowD <- list()
  units <- .volumeUnits(design, full, fixed)
  etaShift <- vector("list", length(units))
  cholIn <- .intersect(outer, .hyperNames(design))

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

  for (i in seq_along(units)) {
    u <- units[[i]]
    eta <- u$eta
    Kn <- length(eta)
    Hinv <- H_inv_list[[i]][eta, eta, drop = FALSE]
    # the rows of every condition of the unit in the outer parameters of the condition
    parts <- lapply(u$conditions, function(cn) {
      Ac <- .intersect(outer, .frameParameters(frames[[cn]]))
      rw <- .volumeRows(design, cn, frames[[cn]], eta, Ac, hessian)
      rw$A <- Ac
      rw$loss <- .rowLossDerivs(rw$p, rw$s, rw$y, rw$cens, design$df)
      rw
    })
    A <- unique(c(unlist(lapply(parts, `[[`, "A"), use.names = FALSE), cholIn))
    Qa <- length(A)

    # the exact joint blocks H_ee and H_eA, and d eta*/d theta
    Hee <- u$P
    HeA <- matrix(0, Kn, Qa, dimnames = list(eta, A))
    for (rw in parts) {
      jb <- .jointBlocks(rw, rw$loss, u, full, rw$A, character(0))
      Hee <- Hee + jb$Hee - u$P
      HeA[, rw$A] <- HeA[, rw$A] + jb$HeA
    }
    for (cp in .intersect(cholIn, names(u$dP))) HeA[, cp] <- HeA[, cp] + u$dP[[cp]] %*% full[eta]
    if (!all(is.finite(HeA)) || !all(is.finite(Hee)))
      stop("FOCEI: the second derivatives of the model are not finite at these parameters.", call. = FALSE)
    thj <- colSums(abs(HeA)) > 0
    Em <- matrix(0, Kn, Qa, dimnames = list(eta, A))
    if (any(thj)) {
      eg <- eigen(Hee, symmetric = TRUE)
      lam <- pmax(eg$values, 1e-10 * abs(sum(eg$values)) / Kn)
      Em[, thj] <- -(eg$vectors %*% (t(eg$vectors) / lam)) %*% HeA[, thj, drop = FALSE]
    }
    if (sparse) etaShift[[i]] <- Em
    else {
      etaShift[[i]] <- matrix(0, Kn, Q, dimnames = list(eta, outer))
      etaShift[[i]][, A] <- Em
    }

    # per condition the local terms, their couplings N to the shift and the K x K terms Q
    gE <- numeric(Kn)
    N <- matrix(0, Qa, Kn, dimnames = list(A, NULL))
    Qm <- matrix(0, Kn, Kn)
    D <- matrix(0, Kn^2, Qa, dimnames = list(NULL, A))
    Gam <- matrix(0, Kn^2, Kn)
    for (rw in parts) {
      t <- .volumeTerms(rw, Hinv, hessian)
      Ac <- rw$A
      grad[Ac] <- grad[Ac] + t$g
      gE <- gE + t$gE
      if (!hessian) next
      if (length(Ac)) blocks[[length(blocks) + 1L]] <- structure(t$local, dimnames = list(Ac, Ac))
      N[Ac, ] <- N[Ac, ] + t$N
      Qm <- Qm + t$Q
      D[, Ac] <- D[, Ac] + t$D
      Gam <- Gam + t$Gam
    }
    grad[A] <- grad[A] + drop(crossprod(Em, gE))
    for (cp in .intersect(cholIn, names(u$dP))) grad[cp] <- grad[cp] + sum(Hinv * u$dP[[cp]])
    if (!hessian) next

    # dH by every parameter, scaled by H^-1/2 on both sides
    for (cp in .intersect(cholIn, names(u$dP))) D[, cp] <- D[, cp] + as.vector(u$dP[[cp]])
    D <- D + Gam %*% Em
    eh <- eigen((Hinv + t(Hinv)) / 2, symmetric = TRUE)
    Sh <- eh$vectors %*% (sqrt(pmax(eh$values, 0)) * t(eh$vectors))
    W <- kronecker(Sh, Sh) %*% D
    # N Em + Em' N' + Em' Q Em = [N, Em'] B [N, Em']' with B = [0, I; I, Q]
    B <- rbind(cbind(matrix(0, Kn, Kn), diag(Kn)), cbind(diag(Kn), (Qm + t(Qm)) / 2))
    eb <- eigen(B, symmetric = TRUE)
    U1 <- cbind(N, t(Em)) %*% eb$vectors
    U <- cbind(U1, t(W))
    d <- c(eb$values, rep(-1, nrow(W)))
    own <- .intersect(cholIn, names(u$dP))
    if (length(own)) {
      C2 <- matrix(0, length(own), length(own), dimnames = list(own, own))
      for (al in own) for (be in own) C2[al, be] <- sum(Hinv * u$d2P(al, be))
      blocks[[length(blocks) + 1L]] <- C2
    }
    # a unit over few parameters as one dense block, else the parts of rank 2K + K^2 apart
    if (!sparse || Qa^2 <= length(d) * Q) {
      blocks[[length(blocks) + 1L]] <- structure(U %*% (d * t(U)), dimnames = list(A, A))
    } else {
      rownames(U) <- A
      lowU[[length(lowU) + 1L]] <- U
      lowD[[length(lowD) + 1L]] <- d
    }
  }
  fx <- .intersect(names(fixed), outer)
  grad[fx] <- 0
  if (!hessian) return(list(gradient = grad, hessian = NULL, etaShift = etaShift))
  hess <- .sparseFromBlocks(blocks, outer)
  lowrank <- if (length(lowU)) .lowrankOver(Reduce(.lowrankBind, Map(function(U, d) list(U = U, d = d), lowU, lowD)), outer)
  if (length(fx)) {
    hess[fx, ] <- hess[, fx] <- 0
    if (!is.null(lowrank)) lowrank$U[fx, ] <- 0
  }
  if (!sparse) {
    hess <- as.matrix(.lowrankDense(hess, lowrank))
    lowrank <- NULL
  }
  list(gradient = grad, hessian = hess, lowrank = lowrank, etaShift = etaShift)
}

## The volume terms of the rows `rw` of one condition over the parameters rw$A of the condition, F = (Da, Db) = S + C Em
## and Y = (Dp, Ds) = T + c Em: the gradient `g` and the coefficient `gE` of it on Em; with `hessian` the
## local block, the couplings `N` and `Q` to Em, and dH as `D` and `Gam` (on Em).
.volumeTerms <- function(rw, Hinv, hessian)
  volume_terms(rw$loss, rw$a, rw$b, rw$Wp, rw$Ws, rw$Vp, rw$Vs, rw$pA, rw$sA, rw$pAA, rw$sAA, Hinv,
               length(rw$p), length(rw$A), hessian)


## FOCEI fit of the normal random effects of `design` with the joint objective `joint` (data
## plus normal prior) from `start`. The conditional covariance of the etas of a group is
## (H_g / 2)^-1, H_g the inner Hessian of the -2 log density.
.fitFocei <- function(design, joint, start, fixed = NULL, control = list()) {
  if (.foceiSparseChoice(design, joint, start, fixed, control$sparse %||% "auto"))
    return(.fitFoceiSparse(design, joint, start, fixed, control))
  ic <- modifyList(list(rinit = 1, rmax = 10, iterlim = 30, steptol = 1e-8, mtol = 1e-12,
                        eigenFloorRelative = 1e-10), control$innerControl %||% list())
  oc <- .foceiOuterControl(control$trustControl, oneGrouping = TRUE)
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
  # the volume derivatives contract dH with the weights of the integrator, the Schur complement takes H_NN^-1
  volume <- function(full_pars, joint_hessian, H_inv_list, weight_list, data_gradient) {
    ff <- mode == "forward-forward"
    vol <- if (ff) .foceiVolume(design, full_pars, fixed, outer, weight_list, hessian = TRUE, cores = cores)
           else .volumeReverse(design, full_pars, fixed, outer, weight_list, cores = cores)
    curvature <- if (ff) vol$hessian[outer, outer] else 0
    if (ff) last <<- list(full = full_pars, hessian = vol$hessian)
    base <- .gnSchur(joint_hessian, outer, meta$eta_names, H_inv_list) + curvature
    learn$add(full_pars[outer], data_gradient[outer] + vol$gradient[outer], function(x) base %*% x)
    list(gradient = vol$gradient, hessian = curvature + learn$dense(), eta_shift = vol$etaShift)
  }
  failure <- new.env(parent = emptyenv())
  run <- function(inner, start, warm, settings, newton) tryCatch(
    focei_run(model_cb = design$prdfn, err_cb = .errorModel(design), joint_cb = joint,
              init = start, subject_meta = meta, fixed = fixed,
              control = list(inner = inner, outer = settings, newtonValue = newton),
              correction_mode = if (correction) "eager" else "none",
              correction_cb = if (correction) volume else NULL,
              predict_cb = .guardPredictor(.batchPredictor(design), failure), eta_init = warm),
    error = function(e) stop(conditionMessage(e), if (!is.null(failure$message))
      paste0("; the last prediction that failed: ", failure$message), call. = FALSE))
  fit <- .sqpThenNested(run, .sqpInner(control, ic), ic, start, oc)

  E <- design$etaNames[, design$normal$stems, drop = FALSE]
  etaModes <- fit$etaModes
  dimnames(etaModes) <- list(design$groups, design$normal$stems)
  hessian <- fit$hessian
  full <- c(fit$argument, setNames(as.vector(etaModes), as.vector(E)))
  exact <- mode == "forward-forward" || (mode == "forward-reverse" && .jointDeriv2(joint, full, fixed))
  if (mode == "forward-reverse" && !exact)
    hessian <- .reportedReverse(design, joint, full, fixed, outer, fit$H_inv, cores)$hessian
  if (exact) {
    # exact Schur complement of the joint Hessian at the modes plus the volume Hessian
    jh <- joint(pars = full, fixed = fixed, deriv = TRUE, deriv2 = TRUE)$hessian
    on <- rownames(hessian)
    ee <- .intersect(as.vector(E), rownames(jh))
    hessian <- jh[on, on, drop = FALSE] -
      crossprod(jh[ee, on, drop = FALSE], .scaledSolve(jh[ee, ee, drop = FALSE], jh[ee, on, drop = FALSE]))
    # the volume Hessian of the last correction where the fit ends on the correction point
    vol <- if (!is.null(last) && identical(names(last$full), names(full)) &&
               identical(unname(unclass(last$full)), unname(full))) last$hessian
           else .foceiVolume(design, full, fixed, outer, fit$weight, hessian = TRUE, cores = cores)$hessian
    hessian <- hessian + vol[on, on]
  }
  list(argument = fit$argument, value = fit$value, gradient = fit$gradient, hessian = hessian,
       etaModes = etaModes, etaCov = setNames(lapply(fit$H_inv, function(H) 2 * H), design$groups),
       converged = isTRUE(fit$converged), iterations = fit$iterations, innerIterations = fit$innerIterations,
       sqpIterations = fit$sqpIterations, trace = fit$trace, derivMode = mode, hessianExact = exact)
}

## The inner settings of the inexact phase of control$inner = "sqp", `ic` with iterlim
## `ic$sqpIterlim` (default 2); NULL with control$inner "nested".
.sqpInner <- function(control, ic) {
  if (!identical(control[["inner"]], "sqp")) return(NULL)
  modifyList(ic, list(iterlim = as.integer(ic$sqpIterlim %||% 2L)))
}

## An outer search by `run(inner, start, warm, settings, newton)` with settings `oc`: with `sqp` first on the
## inexact modes of `sqp` to the tolerances of .sqpTolerance(), then from there on the modes of `ic`, so the result is that
## of the nested search; the iterations add up, `sqpIterations` those of the first.
.sqpThenNested <- function(run, sqp, ic, start, oc) {
  if (is.null(sqp) || oc$iterlim <= 0) {
    fit <- run(ic, start, NULL, oc, FALSE)
    fit$innerIterations <- as.integer(fit$innerIterations)
    return(fit)
  }
  first <- run(sqp, start, NULL, .sqpTolerance(oc), TRUE)
  fit <- run(ic, first$argument, first$etaModes, modifyList(oc, list(iterlim = oc$iterlim - first$iterations)),
             FALSE)
  fit$iterations <- first$iterations + fit$iterations
  fit$sqpIterations <- first$iterations
  fit$innerIterations <- as.integer(first$innerIterations + fit$innerIterations)
  if (!is.null(first$trace)) fit$trace <- Map(c, first$trace, lapply(fit$trace, `[`, -1L))
  fit
}

## The outer tolerances of the search on inexact modes, the square roots of ftol and mtol of `oc`, flat
## or in `tolControl`: the value there is off by the modes, and the nested search ends the fit.
.sqpTolerance <- function(oc) {
  loose <- function(x) lapply(x[.intersect(names(x), c("ftol", "mtol"))], sqrt)
  oc <- modifyList(oc, loose(oc))
  if (!is.null(oc$tolControl)) oc$tolControl <- modifyList(oc$tolControl, loose(oc$tolControl))
  oc
}

## The outer trust settings of FOCEI from `control$trustControl`. Flat ftol, mtol, gtol and xtol join
## `tolControl`. Over several groupings the outer search is trust(), and every other argument of
## trust() passes; one grouping runs a separate outer loop, which reads rinit, rmax and iterlim.
.foceiOuterControl <- function(tc, oneGrouping) {
  tc <- tc %||% list()
  nms <- names(tc)
  if (length(tc) && (is.null(nms) || !all(nzchar(nms))))
    stop("em: all entries of control$trustControl must be named.", call. = FALSE)
  flat <- c("ftol", "mtol", "gtol", "xtol")
  rest <- tc[.setdiff(nms, c(flat, "tolControl"))]
  tol <- modifyList(list(ftol = 1e-7, mtol = 1e-5), c(tc$tolControl %||% list(), tc[.intersect(nms, flat)]))
  if ("sparse" %in% names(rest))
    stop("em: sparse is an argument of em(), not an entry of control$trustControl.", call. = FALSE)
  if (oneGrouping) {
    bad <- .setdiff(names(rest), c("rinit", "rmax", "iterlim"))
    if (length(bad))
      stop("em: control$trustControl: ", paste(bad, collapse = ", "), " reach trust() with random ",
           "effects over several groupings only. FOCEI of one grouping runs a separate outer trust ",
           "region on rinit, rmax, iterlim, ftol and mtol.", call. = FALSE)
    return(c(modifyList(list(rinit = 1, rmax = 10, iterlim = 200L), rest), tol))
  }
  out <- .trustControl(list(rinit = 1, rmax = 10, iterlim = 200L), rest, label = "em: control$trustControl")
  if (is.null(out$phases) && is.null(out$hessianMethod)) out$hessianMethod <- "gn"
  out$tolControl <- tol
  out
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
  other <- .setdiff(c("forward-reverse", "forward-forward"), want)
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

## The SR1 correction of an outer model from the steps s between the points of `add()` and r = y - B s, y the
## change of the outer gradient and `times(x)` = B x: the last `memory` pairs as U diag(d) U' (`factors()`) or
## a matrix (`dense()`); a pair with |r's| <= 1e-8 |r| |s| is skipped.
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
    U <- matrix(0, length(outer), 0, dimnames = list(outer, NULL)); d <- numeric(0)
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
## exact Schur complement up to 500 outer parameters, the Gauss-Newton one above, plus the volume
## curvature from dH/dtheta; a list of `hessian` (sparse with `sparse`) and `lowrank` over `outer`.
.reportedReverse <- function(design, joint, full, fixed, outer, H_inv, cores = 1L, sparse = FALSE) {
  vol <- .volumeReverse(design, full, fixed, outer, H_inv, cores = cores, local = sparse, curvature = TRUE)$hessian
  if (length(outer) <= 500L) {
    S <- .schurReverse(design, full, fixed, outer, cores)
    return(list(hessian = if (sparse) Matrix::Matrix(S, sparse = TRUE) + vol else S + vol))
  }
  units <- .volumeUnits(design, full, fixed)
  jt <- joint(pars = full, fixed = fixed, deriv = TRUE, hessian = if (sparse) "sparse")
  if (sparse) {
    parts <- .hessianParts(jt, rownames(jt$hessian), outer)
    sc <- .schurParts(parts$hessian, outer, lapply(units, `[[`, "eta"),
                      lapply(seq_along(units), function(i) H_inv[[i]][units[[i]]$eta, units[[i]]$eta, drop = FALSE]),
                      parts$lowrank)
    return(list(hessian = sc$hessian + vol, lowrank = sc$lowrank))
  }
  H <- jt$hessian
  S <- H[outer, outer, drop = FALSE]
  for (i in seq_along(units)) {
    e <- units[[i]]$eta
    Hoe <- H[outer, e, drop = FALSE]
    S <- S - Hoe %*% H_inv[[i]][e, e, drop = FALSE] %*% Matrix::t(Hoe)
  }
  list(hessian = S + vol)
}


## The Schur complement H_oo - sum_g H_og H_g^-1 H_go over the etas `groups[[g]]`, Hinv[[g]] = H_g^-1,
## plus the low-rank part `lowrank` of H, as `hessian` and `lowrank`; a group whose n_g coupled outer
## parameters have n_g^2 > K_g P as the part (H_og R_g') diag(-1) (H_og R_g')', R_g' R_g = H_g^-1.
.schurParts <- function(H, outer, groups, Hinv, lowrank = NULL) {
  P <- length(outer)
  fold <- logical(length(groups))
  U <- vector("list", length(groups))
  for (g in seq_along(groups)) {
    Hoe <- H[outer, groups[[g]], drop = FALSE]
    n <- sum(Matrix::rowSums(abs(Hoe)) > 0)
    R <- if (n^2 > length(groups[[g]]) * P)
      tryCatch(chol((as.matrix(Hinv[[g]]) + t(as.matrix(Hinv[[g]]))) / 2), error = function(e) NULL)
    if (is.null(R)) { fold[g] <- TRUE; next }
    U[[g]] <- as.matrix(Hoe %*% t(R))
  }
  S <- H[outer, outer, drop = FALSE]
  if (any(fold)) {
    e <- unlist(groups[fold], use.names = FALSE)
    Hoe <- H[outer, e, drop = FALSE]
    S <- S - Hoe %*% Matrix::bdiag(lapply(Hinv[fold], function(h) unname(as.matrix(h)))) %*% Matrix::t(Hoe)
  }
  U <- do.call(cbind, U[!fold])
  if (!is.null(U)) rownames(U) <- outer
  list(hessian = Matrix::forceSymmetric(methods::as(S, "CsparseMatrix")),
       lowrank = .lowrankJoin(.lowrankOver(lowrank, outer), if (!is.null(U)) list(U = U, d = rep(-1, ncol(U)))))
}

## Whether FOCEI takes the sparse outer path: `sparse` TRUE, or "auto" with at least the parameters
## of getOption("dMod2.sparseAuto") and a step on the outer Hessian of .schurParts() that costs less
## than a dense one (.sparseCost()).
.foceiSparseChoice <- function(design, joint, start, fixed, sparse = "auto") {
  if (!identical(sparse, "auto")) return(isTRUE(sparse))
  th <- .sparseThresholds()
  outer <- names(start)
  if (length(outer) < th$parameters) return(FALSE)
  etas <- .allEtas(design)
  full <- c(start, setNames(numeric(length(etas)), etas))
  jt <- tryCatch(joint(pars = full, fixed = fixed, deriv = TRUE, hessian = "sparse"), error = function(e) NULL)
  if (is.null(jt) || !.isSparse(jt$hessian)) return(FALSE)
  H <- .hessianParts(jt, c(outer, etas), outer)$hessian != 0
  E <- design$etaNames[, design$normal$stems, drop = FALSE]
  P <- length(outer)
  wide <- vapply(seq_len(nrow(E)), function(g) sum(Matrix::rowSums(H[outer, E[g, ], drop = FALSE]) > 0)^2 > ncol(E) * P, TRUE)
  B <- H[outer, .intersect(as.vector(t(E[!wide, , drop = FALSE])), etas), drop = FALSE]
  pattern <- methods::as(H[outer, outer, drop = FALSE] | (B %*% Matrix::t(B)) != 0, "CsparseMatrix")
  if (!is.null(th$fill)) return(Matrix::nnzero(pattern) / P^2 <= th$fill)
  cost <- .sparseCost(pattern, sum(wide) * ncol(E) + length(jt$lowrank$d))
  cost[["sparse"]] < cost[["dense"]]
}

## FOCEI of one grouping on sparse Hessians by the outer loop of the dense path (focei_run_sparse()):
## the joint Gauss-Newton Hessian minus the Schur complements plus the volume Hessian, local to the
## groups or low-rank, with the SR1 correction.
.fitFoceiSparse <- function(design, joint, start, fixed = NULL, control = list()) {
  ic <- modifyList(list(rinit = 1, rmax = 10, iterlim = 30, steptol = 1e-8, mtol = 1e-12,
                        eigenFloorRelative = 1e-10), control$innerControl %||% list())
  oc <- .foceiOuterControl(control$trustControl, oneGrouping = TRUE)
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

  # the modes and the value at theta, the inner trust of the settings `settings` started from `warm`;
  # with `newton` the value at the exact modes to second order in the eta gradient the inner trust leaves
  modesAt <- function(theta, warm, settings = ic, newton = FALSE) {
    L <- .choleskyFactor(lay, c(theta, fixed))
    full <- setNames(numeric(length(meta$pars_full_names)), meta$pars_full_names)
    full[outer] <- theta
    inner <- nlme_group_kernel(design$prdfn, errfn, full, warm, meta, chol2inv(t(L)),
                               2 * sum(log(diag(L))), fixed, control = settings, predict_cb = predictor)
    full[etas] <- as.vector(inner$eta_modes)
    jt <- joint(pars = full, fixed = fixed, deriv = TRUE, hessian = "sparse", conditions = meta$conditions)
    decrement <- if (newton) sum(vapply(seq_len(N), function(g) {
      r <- jt$gradient[E[g, ]]
      sum(r * (inner$H_inv[[g]] %*% r))
    }, 0)) else 0
    list(theta = theta, full = full, inner = inner, joint = jt, iterations = inner$iterations,
         value = jt$value + sum(inner$term) + inner$constant - decrement / 2)
  }
  # the gradient and the sparse Hessian of the outer objective, and d eta*/d theta per group
  derivativesAt <- function(ev, deriv2 = FALSE) {
    jt <- if (deriv2) joint(pars = ev$full, fixed = fixed, deriv = TRUE, deriv2 = TRUE, hessian = "sparse",
                            conditions = meta$conditions)
          else ev$joint
    # a low-rank part on the outer parameters passes the Schur complement over the etas as it is
    parts <- .hessianParts(jt, c(outer, etas), outer)
    H <- parts$hessian
    Hoe <- H[outer, byGroup, drop = FALSE]
    groups <- lapply(seq_len(N), function(g) E[g, ])
    Hinv <- if (deriv2) lapply(groups, function(e) solve(as.matrix(H[e, e, drop = FALSE])))
            else ev$inner$H_inv
    sc <- .schurParts(H, outer, groups, Hinv, parts$lowrank)
    ev$hessian <- sc$hessian
    ev$lowrank <- sc$lowrank
    # the gradient at the exact modes to first order: the eta gradient the inner trust leaves
    # enters through the coupling
    left <- Matrix::bdiag(lapply(ev$inner$H_inv, unname)) %*% jt$gradient[byGroup]
    ev$gradient <- jt$gradient[outer] - as.numeric(Hoe %*% left)
    if (!correction) return(ev)
    # forward over reverse: the exact volume gradient on the Gauss-Newton Schur model, whose
    # curvature the SR1 correction learns
    if (mode == "forward-reverse" && !deriv2) {
      vol <- .volumeReverse(design, ev$full, fixed, outer, ev$inner$weight, cores = cores, local = TRUE)
      ev$gradient <- ev$gradient + vol$gradient[outer]
      ev$shift <- vol$etaShift
      return(ev)
    }
    vol <- .foceiVolume(design, ev$full, fixed, outer, ev$inner$weight, hessian = "sparse", cores = cores)
    ev$gradient <- ev$gradient + vol$gradient[outer]
    ev$hessian <- Matrix::forceSymmetric(ev$hessian + vol$hessian)
    ev$lowrank <- .lowrankJoin(ev$lowrank, vol$lowrank)
    ev$shift <- vol$etaShift
    ev
  }
  # the loop of the dense path on the sparse model: the trial value, and at an accepted trial
  # the derivatives with the SR1 correction of the outer model as part of the low-rank part and
  # the mode shifts over the outer positions of their columns
  learn <- .outerSecant(outer)
  trial <- function(theta, warm, settings = ic, newton = FALSE) {
    ev <- modesAt(setNames(as.numeric(theta), outer), warm, settings, newton)
    ev$eta_modes <- ev$inner$eta_modes
    ev
  }
  complete <- function(ev) {
    ev <- derivativesAt(ev)
    ev$gradient <- as.numeric(ev$gradient[outer])
    lr <- .lowrankOver(ev$lowrank, outer)
    if (correction) {
      learn$add(ev$theta, setNames(ev$gradient, outer), function(x) .lowrankApply(ev$hessian, lr, x))
      lr <- .lowrankJoin(lr, learn$factors())
    }
    ev$kernel <- .kernelHessian(list(hessian = ev$hessian, lowrank = lr), outer)$kernel
    ev$shift <- lapply(ev$shift, function(S)
      if (length(S)) structure(matrix(as.numeric(S), nrow(S)), idx = match(colnames(S), outer)))
    ev
  }
  run <- function(inner, start, warm, settings, newton) {
    res <- focei_run_sparse(function(theta, warm) trial(theta, warm, inner, newton), complete,
                            setNames(as.numeric(start), outer), N, K, settings, ic$rinit, warm)
    c(res, list(etaModes = res$keep$inner$eta_modes))
  }
  res <- .sqpThenNested(run, .sqpInner(control, ic), ic, start, oc)
  cur <- res$keep
  converged <- res$converged
  it <- res$iterations
  trace <- res$trace
  # the exact Schur complement of the joint Hessian at the modes plus the volume Hessian; by
  # forward over reverse with the volume curvature from dH/dtheta where there is no deriv2
  exact <- mode == "forward-forward" || (mode == "forward-reverse" && .jointDeriv2(joint, cur$full, fixed))
  final <- if (exact) derivativesAt(cur, deriv2 = TRUE) else cur
  if (mode == "forward-reverse" && !exact) {
    rh <- .reportedReverse(design, joint, cur$full, fixed, outer, cur$inner$H_inv, cores, sparse = TRUE)
    final$hessian <- rh$hessian
    final$lowrank <- rh$lowrank
  }
  etaModes <- cur$inner$eta_modes
  dimnames(etaModes) <- list(design$groups, lay$stems)
  list(argument = cur$theta, value = cur$value, gradient = setNames(cur$gradient, outer), hessian = final$hessian,
       lowrank = final$lowrank, etaModes = etaModes, etaCov = setNames(lapply(cur$inner$H_inv, function(H) 2 * H), design$groups),
       converged = converged, iterations = it, innerIterations = res$innerIterations,
       sqpIterations = res$sqpIterations, trace = trace, derivMode = mode, hessianExact = exact)
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
  # the inner settings and, on inexact modes (control$inner = "sqp"), the modes predicted along d eta*/d theta
  # from the last evaluation; `innerIterations` counts the evaluated inner steps
  inner <- ic
  newton <- FALSE
  predict <- NULL
  innerIterations <- 0L
  # the objectives per condition of the R evaluation, built where the chain does not cover the model
  delayedAssign("objs", .conditionObjectives(design))
  # The modes of all components in lock step (focei_units_lockstep()): per round one batched
  # evaluation of the conditions of the components still moving, then per component a trust step on
  # the component data and the normal prior of the component etas, e' P e / 2, accepted by a separate ratio.
  modes <- function(theta) {
    held <- c(theta, fixed)
    if (!is.null(predict)) warm <<- predict$modes + drop(predict$shift %*% (theta - predict$theta))
    P <- lapply(.volumeUnits(design, c(theta, fixed, warm), fixed), function(u) {
      m <- as.matrix(u$P); storage.mode(m) <- "double"; m
    })
    w <- focei_units_lockstep(.unitData(design, units, objs, held), lapply(units, function(u) unname(warm[u$eta])),
                              P, inner)
    innerIterations <<- innerIterations + sum(attr(w, "iterations"))
    for (i in seq_along(units)) warm[units[[i]]$eta] <<- w[[i]]
    warm
  }
  outerObjective <- function(theta, deriv2 = FALSE, volumeHessian = FALSE, ...) {
    theta <- theta[outer]
    full <- c(theta, modes(theta))
    gn <- joint(pars = full, fixed = fixed, deriv = TRUE, hessian = curvature)
    # sparse: a low-rank part on the outer parameters passes the Schur complements as it is
    sparseH <- function(h) if (sparse) .hessianParts(h, rownames(h$hessian), outer) else list(hessian = h$hessian)
    parts <- sparseH(gn)
    gn$hessian <- parts$hessian
    lowrank <- parts$lowrank
    block <- function(u) as.matrix(gn$hessian[u$eta, u$eta, drop = FALSE])
    Hinv <- lapply(units, function(u) .scaledSolve(block(u)))
    logdet <- sum(vapply(units, function(u) as.numeric(determinant(block(u) / 2)$modulus), 0))
    cores <- getOption("dMod2.cores", 1L)
    decrement <- 0
    if (newton) {
      corrected <- .inexactModes(gn, Hinv, units, etas, outer)
      decrement <- corrected$decrement
      gn$gradient <- corrected$gradient
      predict <<- list(theta = theta, modes = full[etas], shift = corrected$shift)
    }
    schur <- function(h) {
      ee <- .intersect(etas, rownames(h))
      h[outer, outer, drop = FALSE] -
        Matrix::crossprod(h[ee, outer, drop = FALSE], Matrix::solve(h[ee, ee, drop = FALSE], h[ee, outer, drop = FALSE]))
    }
    if (mode == "forward-forward") {
      jp <- sparseH(joint(pars = full, fixed = fixed, deriv = TRUE, deriv2 = TRUE, hessian = curvature))
      jh <- jp$hessian
      lowrank <- jp$lowrank
      vol <- .foceiVolume(design, full, fixed, outer, Hinv, cores = cores,
                          hessian = if (!volumeHessian) FALSE else curvature %||% TRUE)
      H <- schur(jh)
      if (volumeHessian) {
        H <- H + vol$hessian[outer, outer]
        lowrank <- .lowrankJoin(.lowrankOver(lowrank, outer), .lowrankOver(vol$lowrank, outer))
      }
      volGradient <- vol$gradient[outer]
    } else {
      # forward over reverse: the exact volume gradient on the Gauss-Newton Schur model, whose
      # curvature the SR1 correction of the dense path learns; the exact Schur complement where it
      # is reported
      vol <- if (mode != "none") .volumeReverse(design, full, fixed, outer, Hinv, cores = cores, local = sparse)
      volGradient <- if (mode == "none") 0 else vol$gradient[outer]
      H <- schur(gn$hessian)
      if (mode != "none") {
        learn$add(theta, gn$gradient[outer] + volGradient, function(x) .lowrankApply(H, lowrank, x))
        if (!sparse) H <- H + learn$dense()
      }
      if (volumeHessian && mode != "none") {
        rh <- if (.jointDeriv2(joint, full, fixed)) {
          jp <- sparseH(joint(pars = full, fixed = fixed, deriv = TRUE, deriv2 = TRUE, hessian = curvature))
          vh <- .foceiVolume(design, full, fixed, outer, Hinv, hessian = curvature %||% TRUE, cores = cores)
          list(hessian = schur(jp$hessian) + vh$hessian[outer, outer],
               lowrank = .lowrankJoin(.lowrankOver(jp$lowrank, outer), .lowrankOver(vh$lowrank, outer)))
        } else .reportedReverse(design, joint, full, fixed, outer, Hinv, cores, sparse = sparse)
        H <- rh$hessian
        lowrank <- rh$lowrank
      }
    }
    if (sparse) H <- Matrix::forceSymmetric(Matrix::Matrix(H, sparse = TRUE))
    out <- objlist(value = gn$value + logdet - decrement / 2, gradient = gn$gradient[outer] + volGradient, hessian = H)
    out$lowrank <- .lowrankOver(lowrank, outer)
    attr(out, "Hinv") <- Hinv
    out
  }
  # one outer search; the modes of the last evaluation stay in `warm`, so `etaStart` is not needed
  run <- function(innerSettings, from, etaStart, settings, newtonValue) {
    inner <<- innerSettings
    newton <<- newtonValue
    predict <<- NULL
    args <- c(list(objfun = outerObjective, parinit = from), settings, list(sparse = sparse))
    # iterlim 0 evaluates the marginal at the start, which is no failure to converge
    if (settings$iterlim > 0) do.call(trust, args) else .innerTrust(args)
  }
  fit <- .sqpThenNested(run, .sqpInner(control, ic), ic, start, .foceiOuterControl(control$trustControl, oneGrouping = FALSE))
  final <- outerObjective(fit$argument, volumeHessian = TRUE)
  cov <- unlist(lapply(seq_along(units), function(i) {
    V <- 2 * attr(final, "Hinv")[[i]]
    dimnames(V) <- list(units[[i]]$eta, units[[i]]$eta)
    list(V)
  }), recursive = FALSE)
  list(argument = fit$argument[outer], value = final$value, gradient = final$gradient,
       hessian = final$hessian, lowrank = final$lowrank, etas = warm, etaCovUnits = cov, units = units,
       converged = isTRUE(fit$converged), iterations = fit$iterations, stopReason = fit$stopReason,
       phaseIterations = fit$phaseIterations, innerIterations = innerIterations, sqpIterations = fit$sqpIterations,
       derivMode = mode,
       hessianExact = mode == "forward-forward" || .jointDeriv2(joint, c(fit$argument[outer], warm), fixed))
}

## The data terms of the components `idx` at the etas `etas`, the other parameters at `held`, summed per
## component in C++ where the model allows it; NULL for a component whose prediction fails.
.unitData <- function(design, units, objs, held) {
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
    if (!is.null(d)) return(lapply(d, function(x) if (is.finite(x$value)) x))
    if (length(idx) == 1L) return(list(NULL))
    # the components one by one, so that the others keep their value
    lapply(seq_along(idx), function(j) dataUnits(idx[j], etas[j])[[1L]])
  }
  dataUnits
}

## Value and gradient on inexact modes at the exact modes to first order in the eta gradient r the inner search
## leaves, f - r' H_ee^-1 r / 2 and g_o - H_oe H_ee^-1 r, and d eta*/d theta = -H_ee^-1 H_eo, which predicts the
## modes of the next evaluation.
.inexactModes <- function(gn, Hinv, units, etas, outer) {
  shift <- matrix(0, length(etas), length(outer), dimnames = list(etas, outer))
  decrement <- 0
  gradient <- gn$gradient
  for (i in seq_along(units)) {
    e <- units[[i]]$eta
    Hoe <- as.matrix(gn$hessian[outer, e, drop = FALSE])
    d <- drop(Hinv[[i]] %*% gradient[e])
    decrement <- decrement + sum(gradient[e] * d)
    gradient[outer] <- gradient[outer] - drop(Hoe %*% d)
    shift[e, ] <- -Hinv[[i]] %*% t(Hoe)
  }
  list(decrement = decrement, gradient = gradient, shift = shift)
}

## A gradient restricted to `names`, 0 where absent.
.named <- function(g, names) { out <- setNames(numeric(length(names)), names); h <- .intersect(names, names(g)); out[h] <- g[h]; out }

.namedHessian <- function(H, names) {
  out <- matrix(0, length(names), length(names), dimnames = list(names, names))
  h <- .intersect(names, rownames(H))
  if (length(h)) out[h, h] <- as.matrix(H[h, h, drop = FALSE])
  out
}
