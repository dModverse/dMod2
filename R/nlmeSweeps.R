# Second-order information of the mixed-effects objectives by forward over reverse: one forward
# pass along the etas, several backward sweeps with seeds on the predicted values and sigmas of
# the data rows. The etas of different units never meet in a condition, so the k-th eta of
# every unit shares one tangent.


## The forward pass along the colored eta directions of `units`, then a direction per parameter
## in `along` and the columns of `extra` (rows named by parameters): per condition of `conds` the
## prediction and error-model values with their tangents and tapes, and per data row its
## prediction p, sigma s and their derivatives pd, sd along the directions.
.sweepForward <- function(design, full, fixed, units, cores = 1L, along = character(0), extra = NULL,
                          conds = unique(unlist(lapply(units, `[[`, "conditions"), use.names = FALSE))) {
  nm <- names(full)
  Ke <- max(1L, vapply(units, function(u) length(u$eta), 0L))
  K <- Ke + length(along) + NCOL(extra) * !is.null(extra)
  V <- matrix(0, length(nm), K, dimnames = list(nm, c(paste0("eta", seq_len(Ke)), along,
                                                      if (!is.null(extra)) paste0("extra", seq_len(ncol(extra))))))
  for (u in units) V[cbind(match(u$eta, nm), seq_along(u$eta))] <- 1
  if (length(along)) V[cbind(match(along, nm), Ke + seq_along(along))] <- 1
  if (!is.null(extra)) V[match(rownames(extra), nm), Ke + length(along) + seq_len(ncol(extra))] <- extra
  n <- length(conds)
  env <- new.env()
  pars <- as.parvec(full, deriv = V)
  b <- .bundle(conds = conds, times = lapply(conds, function(cn) .emGrid(design, cn)),
               pars = rep(list(pars), n), fixed = rep(list(fixed), n), shared = FALSE)
  fw <- .fwdMany(design$prdfn, b, env, cores, deriv = TRUE)
  preds <- unname(as.list(fw$values))
  errs <- vector("list", n)
  errTape <- NULL
  if (!is.null(design$errfn)) {
    split <- lapply(preds, function(pr) {
      pin <- getParameters(pr)
      own <- attr(pin, "fixed")
      held <- pin[c(own, .setdiffU(.intersectU(names(pin), names(fixed)), own))]
      list(pars = as.parvec(pin[.setdiffU(names(pin), names(fixed))]),
           fixed = as.parvec(held, deriv = FALSE, deriv2 = FALSE))
    })
    eb <- .bundle(conds = conds, out = preds, pars = lapply(split, `[[`, "pars"),
                  fixed = lapply(split, `[[`, "fixed"), shared = FALSE)
    ef <- .fwdMany(design$errfn, eb, env, cores, deriv = TRUE)
    errs <- unname(as.list(ef$values))
    length(errs) <- n
    errTape <- ef$tape
  }
  rows <- Map(function(cn, pr, er) .sweepRows(design, cn, pr, er, colnames(V)), conds, preds, errs)
  list(conds = conds, names = nm, K = K, Ke = Ke, env = env, tape = fw$tape, errTape = errTape,
       preds = preds, errs = errs, rows = unname(rows))
}

## The data rows of condition `cn` on the frames of a forward pass: values, their derivatives
## along the directions `dirs`, and the positions of the rows in the frames. A recorded sigma has
## no derivatives and no position in the error model.
.sweepRows <- function(design, cn, pred, err, dirs) {
  d <- .conditionData(design, cn)
  along <- function(frame, t, o) {
    D <- attr(frame, "deriv")
    out <- matrix(0, length(t), length(dirs))
    if (is.null(D)) return(out)
    j <- match(dirs, dimnames(D)[[3]])
    i <- t + (match(d$name, dimnames(D)[[2]]) - 1L) * dim(D)[1]
    ok <- which(!is.na(i))
    flat <- matrix(D, prod(dim(D)[1:2]))
    out[ok, !is.na(j)] <- flat[i[ok], j[!is.na(j)], drop = FALSE]
    out
  }
  tp <- .matchNum(d$time, as.numeric(pred[, "time"]))
  op <- match(d$name, colnames(pred))
  rows <- list(p = pred[cbind(tp, op)], pd = along(pred, tp, op), s = d$sigma,
               sd = matrix(0, nrow(d), length(dirs)), y = pmax(d$value, d$lloq), cens = d$value <= d$lloq,
               tp = tp, op = op, model = is.na(d$sigma), te = NULL, oe = NULL)
  if (!any(rows$model)) return(rows)
  if (is.null(err)) err <- .dataSigmaErrfn(design$data)(out = pred, conditions = cn)[[1]]
  te <- .matchNum(d$time, as.numeric(err[, "time"]))
  oe <- match(d$name, colnames(err))
  m <- rows$model
  rows$s[m] <- err[cbind(te, oe)][m]
  rows$sd[m, ] <- along(err, te, oe)[m, , drop = FALSE]
  rows$te <- te; rows$oe <- oe
  rows
}

## The backward sweep of a forward pass `fw` with the seeds `seeds`, per condition a list of
## `pred` and `sigma`, [rows, blocks (K + 1)]: per block a seed and its derivatives along the K
## tangents. Per condition the cotangent of the parameters it reaches, [parameters, blocks
## (K + 1)]: per block J' w and the sum over rows of w Hess(f) v_k plus J' wdot_k.
.sweepBackward <- function(fw, seeds, cores = 1L, blocks = 1L) {
  S <- blocks * (fw$K + 1L)
  n <- length(fw$conds)
  place <- function(frame, t, o, v, keep) {
    W <- array(0, c(nrow(frame), ncol(frame), S), dimnames = list(NULL, colnames(frame), NULL))
    pos <- (t + (o - 1L) * nrow(frame))[keep]
    for (s in seq_len(S))
      W[, , s] <- .scatterAdd(nrow(frame) * ncol(frame), pos, v[keep, s])
    .dropTime(W)
  }
  w <- lapply(seq_len(n), function(i) {
    r <- fw$rows[[i]]
    .ct(out = place(fw$preds[[i]], r$tp, r$op, seeds[[i]]$pred, rep(TRUE, length(r$p))))
  })
  .withSeeds(blocks, {
    if (!is.null(fw$errTape)) {
      we <- lapply(seq_len(n), function(i) {
        r <- fw$rows[[i]]
        if (!any(r$model) || is.null(fw$errs[[i]])) return(.ct())
        .ct(out = place(fw$errs[[i]], r$te, r$oe, seeds[[i]]$sigma, r$model))
      })
      ue <- .bwdNode(fw$errTape, we, fw$env, cores)
      for (i in seq_len(n)) {
        if (i > length(ue) || .ct_null(ue[[i]])) next
        w[[i]] <- .addCt(w[[i]], .ct(out = .dropTime(ue[[i]]$out), pars = ue[[i]]$pars))
      }
    }
    u <- .bwdNode(fw$tape, w, fw$env, cores)
  })
  lapply(seq_len(n), function(i) {
    ct <- if (i <= length(u)) u[[i]]$pars
    if (is.null(ct)) return(matrix(0, 0, S))
    ct <- ct[intersect(rownames(ct), fw$names), , drop = FALSE]
    ct[rowSums(ct != 0) > 0, , drop = FALSE]
  })
}

## A named vector over `names` with the rows of the cotangents `cts`, column `col`, added.
.collect <- function(names, cts, col) {
  out <- setNames(numeric(length(names)), names)
  for (ct in cts) {
    hit <- intersect(rownames(ct), names)
    out[hit] <- out[hit] + ct[hit, col]
  }
  out
}


## Per unit of `units` (its `eta` and `conditions`) by forward over reverse: `rows`, the exact
## eta rows of the Hessian of the data term; `traces`, the gradients grad tr(M H_u) of the matrices
## in `Ms[[u]]`, H_u the Gauss-Newton Hessian of the data term in the unit's etas; and `hv`, the data
## Hessian times the eta vector `vs[[u]]`. All as named over `reach`, the parameters the unit's
## conditions reach. One forward pass and one backward sweep with a seed for the rows and the
## products and one per eta direction for the traces of any number of M.
.sweepUnits <- function(design, full, fixed, units, Ms = NULL, vs = NULL, cores = 1L, fw = NULL) {
  extra <- if (!is.null(vs)) {
    e <- unlist(lapply(units, `[[`, "eta"))
    matrix(unlist(vs), dimnames = list(e, NULL))
  }
  fw <- fw %||% .sweepForward(design, full, fixed, units, cores, extra = extra)
  K <- fw$K; Ke <- fw$Ke
  unitOf <- integer(0)
  for (i in seq_along(units)) unitOf[units[[i]]$conditions] <- i
  loss <- lapply(fw$rows, function(r) .rowLossDerivs(r$p, r$s, r$y, r$cens, design$df))
  l_ <- function(lr, i, j) lr[, i + 1L, j + 1L]
  zero <- function(r) matrix(0, length(r$p), K + 1L)
  # block 1, the data Hessian along every tangent: seed the loss gradient, its derivative
  # L (pd, sd); block j + 1, for grad tr(M H_u) = sum_jk M_jk U_jk with symmetric M: seed
  # 2 L J_j, and along tangent k the seeds of J_j' dL J_k, so that U has the derivatives
  # through J and through L
  B <- 1L + if (length(Ms)) Ke else 0L
  seeds <- Map(function(i, r, lr) {
    out <- list(pred = cbind(l_(lr, 1, 0), l_(lr, 2, 0) * r$pd + l_(lr, 1, 1) * r$sd),
                sigma = cbind(l_(lr, 0, 1), l_(lr, 1, 1) * r$pd + l_(lr, 0, 2) * r$sd))
    if (B == 1L) return(out)
    Ku <- length(units[[unitOf[fw$conds[i]]]]$eta)
    a <- r$pd; b <- r$sd
    for (j in seq_len(Ke)) {
      zp <- zs <- zero(r)
      if (j <= Ku) {
        zp[, 1] <- 2 * (l_(lr, 2, 0) * a[, j] + l_(lr, 1, 1) * b[, j])
        zs[, 1] <- 2 * (l_(lr, 1, 1) * a[, j] + l_(lr, 0, 2) * b[, j])
        for (k in seq_len(Ku)) {
          aa <- a[, j] * a[, k]; ab <- a[, j] * b[, k] + b[, j] * a[, k]; bb <- b[, j] * b[, k]
          zp[, k + 1L] <- l_(lr, 3, 0) * aa + l_(lr, 2, 1) * ab + l_(lr, 1, 2) * bb
          zs[, k + 1L] <- l_(lr, 2, 1) * aa + l_(lr, 1, 2) * ab + l_(lr, 0, 3) * bb
        }
      }
      out$pred <- cbind(out$pred, zp); out$sigma <- cbind(out$sigma, zs)
    }
    out
  }, seq_along(fw$rows), fw$rows, loss)
  all <- .sweepBackward(fw, seeds, cores, blocks = B)
  block <- function(j) lapply(all, function(ct) ct[, (j - 1L) * (K + 1L) + seq_len(K + 1L), drop = FALSE])
  hessRows <- block(1L)
  U <- if (B > 1L) lapply(seq_len(Ke), function(j) block(j + 1L))
  lapply(seq_along(units), function(u) {
    eta <- units[[u]]$eta; Ku <- length(eta)
    cs <- which(unitOf[fw$conds] == u)
    reach <- unique(c(eta, unlist(lapply(c(list(hessRows), U), function(cts) lapply(cts[cs], rownames)))))
    rows <- t(vapply(seq_len(Ku), function(k) .collect(reach, hessRows[cs], k + 1L), numeric(length(reach))))
    dimnames(rows) <- list(eta, reach)
    tr <- list()
    if (length(Ms[[u]])) {
      # column j + (k - 1) Ku holds U_jk, so that U %*% c(M) is the trace gradient
      Uu <- matrix(vapply(seq_len(Ku * Ku), function(q)
        .collect(reach, U[[(q - 1L) %% Ku + 1L]][cs], (q - 1L) %/% Ku + 2L), numeric(length(reach))), length(reach))
      tr <- lapply(Ms[[u]], function(M) setNames(drop(Uu %*% c(M + t(M))) / 2, reach))
    }
    list(reach = reach, rows = rows, traces = tr,
         hv = if (!is.null(vs)) .collect(reach, hessRows[cs], K + 1L))
  })
}


## The Gauss-Newton Hessian of the data term in the etas of every unit, sum over its rows of
## J' L J with J = d(p, s) / d eta from the tangents of the forward pass `fw`.
.unitCurvature <- function(design, fw, units) {
  unitOf <- integer(0)
  for (i in seq_along(units)) unitOf[units[[i]]$conditions] <- i
  lapply(seq_along(units), function(u) {
    eta <- units[[u]]$eta; Ku <- length(eta)
    H <- matrix(0, Ku, Ku, dimnames = list(eta, eta))
    for (i in which(unitOf[fw$conds] == u)) {
      r <- fw$rows[[i]]
      lr <- .rowLossDerivs(r$p, r$s, r$y, r$cens, design$df)
      a <- r$pd[, seq_len(Ku), drop = FALSE]; b <- r$sd[, seq_len(Ku), drop = FALSE]
      Lps <- lr[, 2L, 2L]
      H <- H + crossprod(a, lr[, 3L, 1L] * a) + crossprod(a, Lps * b) + crossprod(b, Lps * a) +
        crossprod(b, lr[, 1L, 3L] * b)
    }
    H
  })
}


## The volume term sum_u log|H_u| of FOCEI by forward over reverse: its gradient in `outer` and
## the mode shift d eta*/d theta = -H_ee^-1 H_eA per unit, H_u^-1 from `H_inv_list`, from the
## eta rows of the joint Hessian and the trace gradient at M = H_u^-1 (`.sweepUnits()`). With
## `curvature` also the volume Hessian from dH_u/dtheta, the trace gradients at the unit matrices:
## -tr(H^-1 D_a H H^-1 D_b H) + tr(H^-1 d2P_ab), without the terms of the curvature of the
## prediction in the etas. With `local` the shift and the Hessian of a unit have the columns of the
## parameters it depends on only.
.volumeReverse <- function(design, full, fixed, outer, H_inv_list, cores = 1L, local = FALSE,
                           curvature = FALSE) {
  units <- .volumeUnits(design, full, fixed)
  # per unit the matrices M of its trace gradients: H^-1, then E_kl for k <= l
  pairs <- lapply(units, function(u) which(upper.tri(diag(length(u$eta)), diag = TRUE), arr.ind = TRUE))
  Ms <- lapply(seq_along(units), function(i) {
    e <- units[[i]]$eta
    out <- list(H_inv_list[[i]][e, e, drop = FALSE])
    if (curvature) for (q in seq_len(nrow(pairs[[i]]))) {
      E <- matrix(0, length(e), length(e), dimnames = list(e, e))
      E[pairs[[i]][q, 1], pairs[[i]][q, 2]] <- E[pairs[[i]][q, 2], pairs[[i]][q, 1]] <- 1
      out[[q + 1L]] <- E / (1 + (pairs[[i]][q, 1] != pairs[[i]][q, 2]))
    }
    out
  })
  sw <- .sweepUnits(design, full, fixed, units, Ms, cores = cores)
  cholIn <- intersect(outer, .hyperNames(design))
  grad <- setNames(numeric(length(outer)), outer)
  etaShift <- blocks <- vector("list", length(units))
  for (u in seq_along(units)) {
    un <- units[[u]]
    eta <- un$eta; Ku <- length(eta)
    reach <- unique(c(sw[[u]]$reach, cholIn))
    A <- intersect(outer, reach)
    rowsH <- matrix(0, Ku, length(reach), dimnames = list(eta, reach))
    rowsH[, sw[[u]]$reach] <- sw[[u]]$rows
    traceOf <- function(m) {
      tr <- setNames(numeric(length(reach)), reach)
      tr[sw[[u]]$reach] <- sw[[u]]$traces[[m]]
      for (cp in intersect(cholIn, names(un$dP))) tr[cp] <- tr[cp] + sum(Ms[[u]][[m]] * un$dP[[cp]])
      tr
    }
    trace <- traceOf(1L)
    # the normal prior: its precision in the eta block, its derivative in the Cholesky parameters
    rowsH[, eta] <- rowsH[, eta] + un$P
    for (cp in intersect(cholIn, names(un$dP))) rowsH[, cp] <- rowsH[, cp] + un$dP[[cp]] %*% full[eta]
    Hee <- rowsH[, eta, drop = FALSE]
    eg <- eigen((Hee + t(Hee)) / 2, symmetric = TRUE)
    lam <- pmax(eg$values, 1e-10 * abs(sum(eg$values)) / Ku)
    Em <- -(eg$vectors %*% (t(eg$vectors) / lam)) %*% rowsH[, A, drop = FALSE]
    dimnames(Em) <- list(eta, A)
    grad[A] <- grad[A] + trace[A] + drop(crossprod(Em, trace[eta]))
    etaShift[[u]] <- Em
    if (curvature) {
      # D_a H_kl along the outer parameters with the mode shift, then the volume Hessian
      DH <- array(0, c(Ku, Ku, length(A)))
      for (q in seq_len(nrow(pairs[[u]]))) {
        tq <- traceOf(q + 1L)
        d <- tq[A] + drop(crossprod(Em, tq[eta]))
        DH[pairs[[u]][q, 1], pairs[[u]][q, 2], ] <- DH[pairs[[u]][q, 2], pairs[[u]][q, 1], ] <- d
      }
      Hinv <- Ms[[u]][[1]]
      P <- array(0, dim(DH))
      for (a in seq_along(A)) P[, , a] <- Hinv %*% DH[, , a]
      V <- -crossprod(matrix(aperm(P, c(2, 1, 3)), Ku * Ku), matrix(P, Ku * Ku))
      own2 <- intersect(cholIn, names(un$dP))
      for (al in own2) for (be in own2)
        V[match(al, A), match(be, A)] <- V[match(al, A), match(be, A)] + sum(Hinv * un$d2P(al, be))
      # the rows and columns of the Cholesky parameters only: between structural parameters the
      # positive term of the curvature of the prediction is missing, and their volume curvature is
      # of the order of one group against the data of all
      structural <- !(A %in% cholIn)
      V[structural, structural] <- 0
      blocks[[u]] <- structure(V, dimnames = list(A, A))
    }
    if (!local) {
      etaShift[[u]] <- matrix(0, Ku, length(outer), dimnames = list(eta, outer))
      etaShift[[u]][, A] <- Em
    }
  }
  fx <- intersect(names(fixed), outer)
  grad[fx] <- 0
  out <- list(gradient = grad, etaShift = etaShift)
  if (!curvature) return(out)
  if (local) out$hessian <- .sparseFromBlocks(blocks, outer)
  else {
    H <- matrix(0, length(outer), length(outer), dimnames = list(outer, outer))
    for (b in blocks) H[rownames(b), rownames(b)] <- H[rownames(b), rownames(b)] + b
    out$hessian <- H
  }
  if (length(fx)) out$hessian[fx, ] <- out$hessian[, fx] <- 0
  out
}


## The exact Hessian of the joint objective (data and normal prior) at the modes in `full`,
## reduced to `outer` by the Schur complement of the etas of every unit: one forward pass along
## the etas and every outer parameter, one backward sweep seeded by the row loss.
.schurReverse <- function(design, full, fixed, outer, cores = 1L) {
  units <- .volumeUnits(design, full, fixed)
  fw <- .sweepForward(design, full, fixed, units, cores, along = outer, conds = names(design$data))
  Ke <- fw$Ke; Q <- length(outer)
  loss <- lapply(fw$rows, function(r) .rowLossDerivs(r$p, r$s, r$y, r$cens, design$df))
  l_ <- function(lr, i, j) lr[, i + 1L, j + 1L]
  rows <- .sweepBackward(fw, Map(function(r, lr) list(
    pred = cbind(l_(lr, 1, 0), l_(lr, 2, 0) * r$pd + l_(lr, 1, 1) * r$sd),
    sigma = cbind(l_(lr, 0, 1), l_(lr, 1, 1) * r$pd + l_(lr, 0, 2) * r$sd)), fw$rows, loss), cores)
  Hoo <- vapply(seq_len(Q), function(q) .collect(outer, rows, 1L + Ke + q), numeric(Q))
  Hoo <- matrix(Hoo, Q, Q, dimnames = list(outer, outer))
  prior <- Reduce(`+`, lapply(.subDesigns(design), function(sd)
    if (sd$normal$K) .normalPrior(sd)(pars = full, fixed = fixed, deriv = TRUE)))
  ph <- prior$hessian
  po <- intersect(outer, rownames(ph))
  Hoo[po, po] <- Hoo[po, po] + ph[po, po]
  unitOf <- integer(0)
  for (i in seq_along(units)) unitOf[units[[i]]$conditions] <- i
  S <- (Hoo + t(Hoo)) / 2
  for (u in seq_along(units)) {
    eta <- units[[u]]$eta
    cs <- which(unitOf[fw$conds] == u & !is.na(unitOf[fw$conds]))
    He <- t(vapply(seq_along(eta), function(k) .collect(c(eta, outer), rows[cs], 1L + k), numeric(length(eta) + Q)))
    dimnames(He) <- list(eta, c(eta, outer))
    pe <- intersect(c(eta, outer), colnames(ph))
    He[, pe] <- He[, pe] + ph[eta, pe, drop = FALSE]
    Hee <- He[, eta, drop = FALSE]
    S <- S - crossprod(He[, outer, drop = FALSE], solve((Hee + t(Hee)) / 2, He[, outer, drop = FALSE]))
  }
  S
}

