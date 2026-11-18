# Second-order information of the mixed-effects objectives by forward over reverse: one forward
# pass along the etas, several backward sweeps with seeds on the predicted values and sigmas of
# the data rows. The etas of different units never meet in a condition, so the k-th eta of
# every unit shares one tangent.


## The forward pass along the colored eta directions of `units`, then a direction per parameter
## in `along`: per condition of `conds` the prediction and error-model values with their tangents
## and tapes, and per data row the prediction p, sigma s and their derivatives pd, sd along the
## directions. Through the native chain where it covers the model (`native`, for a backward sweep
## with `blocks` seed blocks), then without values and tapes.
.sweepForward <- function(design, full, fixed, units, cores = 1L, along = character(0),
                          conds = unique(unlist(lapply(units, `[[`, "conditions"), use.names = FALSE)),
                          blocks = 1L, chain = TRUE) {
  nm <- names(full)
  Ke <- max(1L, vapply(units, function(u) length(u$eta), 0L))
  K <- Ke + length(along)
  V <- matrix(0, length(nm), K, dimnames = list(nm, c(paste0("eta", seq_len(Ke)), along)))
  for (u in units) V[cbind(match(u$eta, nm), seq_along(u$eta))] <- 1
  if (length(along)) V[cbind(match(along, nm), Ke + seq_along(along))] <- 1
  nat <- if (chain) .sweepChainForward(design, full, fixed, V, conds, blocks, cores)
  if (!is.null(nat))
    return(list(conds = conds, names = nm, K = K, Ke = Ke, rows = nat$rows, native = nat$native,
                redo = function() .sweepForward(design, full, fixed, units, cores, along, conds, chain = FALSE)))
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

## The forward pass along the k-th tangent of the native full pass `fw` on the conditions `conds`,
## for a backward sweep with `blocks` seed blocks: the parameter stages along that tangent anew,
## the solution and the k-th tangent taken from `fw` without a solve. NULL where `fw` did not run
## natively.
.sweepSlice <- function(design, full, fixed, fw, k, conds, blocks, cores) {
  nat <- fw$native
  if (is.null(nat)) return(NULL)
  nm <- names(full)
  V <- nat$V[, k, drop = FALSE]
  times <- lapply(conds, function(cn) .emGrid(design, cn))
  planFor <- function(S) .chainLookup(design$chain, conds, rep(list(nm), length(conds)), names(fixed), nm, TRUE,
                                      FALSE, times, tangents = 1L, seeds = S)
  plan <- planFor(blocks)
  if (is.null(plan)) return(NULL)
  cores <- as.integer(max(1L, cores))
  pre <- .chainParams(plan, nat$values, TRUE, FALSE, cores, V)
  if (is.null(pre)) return(NULL)
  raw <- lapply(match(conds, fw$conds), function(i) {
    r <- nat$raw[[i]]
    nt <- nrow(r$variable); n <- ncol(r$variable)
    r$tangent <- array(r$tangent[nt * n * (k - 1L) + seq_len(nt * n)], c(nt, n, 1L))
    r
  })
  list(conds = conds, names = nm, K = 1L, Ke = 1L,
       native = list(plan = plan, planFor = planFor, blocks = blocks, pre = pre, raw = raw, values = nat$values, V = V))
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
## `pred` and `sigma`, [rows, blocks (K + 1)], or these stacked over the conditions with `n` the
## rows per condition: per block a seed and the seed derivatives along the K tangents. Per condition
## the cotangent of the parameters it reaches, [parameters, blocks (K + 1)]: per block J' w and the
## sum over rows of w Hess(f) v_k plus J' wdot_k.
.sweepBackward <- function(fw, seeds, cores = 1L, blocks = 1L) {
  stacked <- is.matrix(seeds$pred)
  if (!is.null(fw$native)) {
    out <- .sweepChainBackward(fw, if (stacked) seeds else list(
      pred = do.call(rbind, lapply(seeds, `[[`, "pred")), sigma = do.call(rbind, lapply(seeds, `[[`, "sigma"))),
      cores, blocks)
    if (!is.null(out)) return(out)
    fw <- fw$redo()
  }
  if (stacked) seeds <- .splitSeeds(seeds)
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
    ct <- ct[.intersect(rownames(ct), fw$names), , drop = FALSE]
    ct[rowSums(ct != 0) > 0, , drop = FALSE]
  })
}

## The forward pass of the sweeps through the native chain along the columns of `V` (rows the
## names of `full`): the data rows per condition as .sweepRows() gives them, and in `native` what
## .sweepChainBackward() needs. NULL where the chain does not cover the model.
.sweepChainForward <- function(design, full, fixed, V, conds, blocks, cores) {
  if (is.null(design$chain) || !.evalCpp() || !length(conds)) return(NULL)
  nm <- names(full)
  times <- lapply(conds, function(cn) .emGrid(design, cn))
  planFor <- function(S) .chainLookup(design$chain, conds, rep(list(nm), length(conds)), names(fixed), nm, TRUE,
                                      FALSE, times, tangents = ncol(V), seeds = S)
  plan <- planFor(blocks)
  if (is.null(plan)) return(NULL)
  cores <- as.integer(max(1L, cores))
  storage.mode(V) <- "double"
  values <- list(c(as.numeric(full), as.numeric(fixed)))
  pre <- .chainParams(plan, values, TRUE, FALSE, cores, V)
  raw <- if (!is.null(pre)) .chainSolveRun(plan, pre$args, cores)
  rows <- if (!is.null(raw)) chain_rows(plan, pre, raw, cores)
  if (is.null(rows)) return(NULL)
  off <- cumsum(c(0L, rows$n))
  list(rows = lapply(seq_along(conds), function(i) {
    j <- off[i] + seq_len(rows$n[i])
    list(p = rows$p[j], pd = rows$pd[j, , drop = FALSE], s = rows$s[j], sd = rows$sd[j, , drop = FALSE],
         y = rows$y[j], cens = rows$cens[j], model = rows$model[j])
  }), native = list(plan = plan, planFor = planFor, blocks = blocks, pre = pre, raw = raw, values = values,
                    V = V))
}

## The backward sweep of a native forward pass `fw` with the seeds `seeds$pred`, `seeds$sigma` on
## the rows of all conditions, as in .sweepBackward(); NULL where the chain fails, for the R path
## to report it.
.sweepChainBackward <- function(fw, seeds, cores, blocks) {
  nat <- fw$native
  plan <- if (blocks == nat$blocks) nat$plan else nat$planFor(blocks)
  if (is.null(plan)) return(NULL)
  cores <- as.integer(max(1L, cores))
  sp <- seeds$pred
  ss <- seeds$sigma
  storage.mode(sp) <- storage.mode(ss) <- "double"
  back <- chain_back(plan, nat$pre, nat$raw, sp, ss, as.integer(blocks), cores, FALSE)
  rraw <- if (!is.null(back)) .chainSolveRun(plan, back$args, cores, back = TRUE)
  pull <- if (!is.null(rraw))
    chain_pull(plan, nat$values, nat$V, nat$pre$roots, back, rraw, as.integer(blocks), cores, FALSE)
  if (is.null(pull)) return(NULL)
  lapply(pull, function(u) {
    ct <- u$ct
    rownames(ct) <- fw$names[u$rows + 1L]
    ct
  })
}

## The data rows of a forward pass stacked over the conditions of the pass, `n` the rows per condition.
.stackRowsOf <- function(rows) {
  pick <- function(f) unlist(lapply(rows, `[[`, f), use.names = FALSE)
  list(p = pick("p"), s = pick("s"), y = pick("y"), cens = pick("cens"),
       pd = do.call(rbind, lapply(rows, `[[`, "pd")), sd = do.call(rbind, lapply(rows, `[[`, "sd")),
       n = vapply(rows, function(r) length(r$p), 0L))
}

## Stacked seeds list(pred, sigma, n) as a list per condition.
.splitSeeds <- function(seeds) {
  off <- cumsum(c(0L, seeds$n))
  lapply(seq_along(seeds$n), function(i) {
    j <- off[i] + seq_len(seeds$n[i])
    list(pred = seeds$pred[j, , drop = FALSE], sigma = seeds$sigma[j, , drop = FALSE])
  })
}

## The cotangents `cts` per condition summed per unit, `unitOf` the unit of every condition:
## per unit the names and their rows.
.unitSums <- function(cts, unitOf, names) {
  P <- length(names)
  key <- unlist(Map(function(ct, u) (u - 1L) * P + match(rownames(ct), names), cts, unitOf), use.names = FALSE)
  if (!length(key)) return(list())
  S <- rowsum(do.call(rbind, lapply(cts, unname)), key)
  k <- as.integer(rownames(S))
  lapply(split(seq_along(k), (k - 1L) %/% P + 1L), function(i)
    structure(S[i, , drop = FALSE], dimnames = list(names[(k[i] - 1L) %% P + 1L], NULL)))
}

## A named vector over `names` with the rows of the cotangents `cts`, column `col`, added.
.collect <- function(names, cts, col) {
  out <- setNames(numeric(length(names)), names)
  for (ct in cts) {
    hit <- .intersect(rownames(ct), names)
    out[hit] <- out[hit] + ct[hit, col]
  }
  out
}


## Per unit of `units` (with `eta` and `conditions`) by forward over reverse: `rows`, the exact
## eta rows of the Hessian of the data term, and `traces`, the gradients grad tr(M H_u) of the
## matrices in `Ms[[u]]`, H_u the Gauss-Newton Hessian of the data term in the unit's etas. All as
## named over `reach`, the parameters the unit's conditions reach. Either one backward sweep along
## all tangents with a seed block for the rows and one per eta direction for the traces
## (`narrow = FALSE`), or per direction k a forward pass along it alone and a backward sweep with a
## seed for row k and one per M (.sweepNarrow()), whose stores have two columns in place of K + 1;
## NULL takes the narrow one where there are matrices M and few of them.
.sweepUnits <- function(design, full, fixed, units, Ms = NULL, cores = 1L, fw = NULL, narrow = NULL) {
  Ke <- max(1L, vapply(units, function(u) length(u$eta), 0L))
  nM <- max(0L, lengths(Ms))
  narrow <- nM > 0L && (narrow %||% (2L * nM + 2L < Ke + 1L))
  B <- 1L + if (nM && !narrow) Ke else 0L
  fw <- fw %||% .sweepForward(design, full, fixed, units, cores, blocks = B)
  K <- fw$K; Ke <- fw$Ke
  unitOf <- integer(0)
  for (i in seq_along(units)) unitOf[units[[i]]$conditions] <- i
  # the rows of all conditions at once, Ku the number of etas of the unit of a row
  r <- .stackRowsOf(fw$rows)
  lr <- .rowLossDerivs(r$p, r$s, r$y, r$cens, design$df)
  l_ <- function(i, j) lr[, i + 1L, j + 1L]
  if (narrow) return(.sweepNarrow(design, full, fixed, units, Ms, cores, fw, r, l_, unitOf))
  Ku <- rep(vapply(units, function(u) length(u$eta), 0L)[unitOf[fw$conds]], r$n)
  # block 1, the data Hessian along every tangent: seed the loss gradient, the derivative
  # L (pd, sd); block j + 1, for grad tr(M H_u) = sum_jk M_jk U_jk with symmetric M: seed
  # 2 L J_j, and along tangent k the seeds of J_j' dL J_k, so that U has the derivatives
  # through J and through L
  a <- r$pd; b <- r$sd
  sp <- ss <- matrix(0, length(r$p), B * (K + 1L))
  sp[, 1L] <- l_(1, 0); sp[, 1L + seq_len(K)] <- l_(2, 0) * a + l_(1, 1) * b
  ss[, 1L] <- l_(0, 1); ss[, 1L + seq_len(K)] <- l_(1, 1) * a + l_(0, 2) * b
  for (j in seq_len(B - 1L)) {
    at <- j * (K + 1L)
    on <- Ku >= j
    sp[on, at + 1L] <- (2 * (l_(2, 0) * a[, j] + l_(1, 1) * b[, j]))[on]
    ss[on, at + 1L] <- (2 * (l_(1, 1) * a[, j] + l_(0, 2) * b[, j]))[on]
    for (k in seq_len(Ke)) {
      onk <- on & Ku >= k
      aa <- a[, j] * a[, k]; ab <- a[, j] * b[, k] + b[, j] * a[, k]; bb <- b[, j] * b[, k]
      sp[onk, at + 1L + k] <- (l_(3, 0) * aa + l_(2, 1) * ab + l_(1, 2) * bb)[onk]
      ss[onk, at + 1L + k] <- (l_(2, 1) * aa + l_(1, 2) * ab + l_(0, 3) * bb)[onk]
    }
  }
  all <- .sweepBackward(fw, list(pred = sp, sigma = ss, n = r$n), cores, blocks = B)
  sums <- .unitSums(all, unitOf[fw$conds], fw$names)
  lapply(seq_along(units), function(u) {
    eta <- units[[u]]$eta; Ku <- length(eta)
    Su <- sums[[as.character(u)]]
    reach <- unique(c(eta, rownames(Su)))
    M <- matrix(0, length(reach), B * (K + 1L), dimnames = list(reach, NULL))
    if (!is.null(Su)) M[rownames(Su), ] <- Su
    rows <- t(M[, 1L + seq_len(Ku), drop = FALSE])
    dimnames(rows) <- list(eta, reach)
    tr <- list()
    if (length(Ms[[u]])) {
      # column j + (k - 1) Ku holds U_jk, the column k + 1 of block j + 1, so that U %*% c(M) is
      # the trace gradient
      Uu <- M[, as.vector(outer(seq_len(Ku), seq_len(Ku), function(j, k) j * (K + 1L) + 1L + k)), drop = FALSE]
      TM <- Uu %*% matrix(vapply(Ms[[u]], function(Mm) c(Mm + t(Mm)) / 2, numeric(Ku * Ku)), Ku * Ku)
      tr <- lapply(seq_len(ncol(TM)), function(m) setNames(TM[, m], reach))
    }
    list(reach = reach, rows = rows, traces = tr)
  })
}

## The rows and trace gradients of .sweepUnits() direction by direction. Along eta_k, from a
## forward pass along it alone over the conditions whose unit has a k-th eta: row k of the Hessian
## from the loss gradient seeded with the loss derivative L J_k along the tangent, and with M symmetric
## sum_j M_jk U_jk from the seed 2 L J M_k with, along the tangent, J M_k' dL J_k. The rows `r` and
## loss derivatives `l_` are those of the full pass `fw`; a matrix M that is 0 in every unit has no
## seed.
.sweepNarrow <- function(design, full, fixed, units, Ms, cores, fw, r, l_, unitOf) {
  nM <- max(lengths(Ms))
  Msym <- lapply(Ms, function(ml) lapply(ml, function(Mm) (Mm + t(Mm)) / 2))
  live <- which(vapply(seq_len(nM), function(m)
    any(vapply(Msym, function(ml) length(ml) >= m && any(ml[[m]] != 0), TRUE)), TRUE))
  B <- 1L + length(live)
  Ku <- vapply(units, function(u) length(u$eta), 0L)
  rowUnit <- rep(unitOf[fw$conds], r$n)
  rowCond <- rep(fw$conds, r$n)
  a <- r$pd; b <- r$sd
  rowsOf <- lapply(units, function(u) list())
  trOf <- lapply(units, function(u) list())
  for (k in seq_len(fw$Ke)) {
    us <- which(Ku >= k)
    sel <- rowUnit %in% us
    if (!any(sel)) next
    conds <- fw$conds[fw$conds %in% rowCond[sel]]
    L <- function(i, j) l_(i, j)[sel]
    x <- a[sel, k]; y <- b[sel, k]
    sp <- ss <- matrix(0, sum(sel), 2L * B)
    sp[, 1L] <- L(1, 0); sp[, 2L] <- L(2, 0) * x + L(1, 1) * y
    ss[, 1L] <- L(0, 1); ss[, 2L] <- L(1, 1) * x + L(0, 2) * y
    for (q in seq_along(live)) {
      A <- Bt <- numeric(length(r$p))
      for (u in us) {
        if (length(Msym[[u]]) < live[q]) next
        on <- rowUnit == u
        w <- Msym[[u]][[live[q]]][, k]
        A[on] <- a[on, seq_len(Ku[u]), drop = FALSE] %*% w
        Bt[on] <- b[on, seq_len(Ku[u]), drop = FALSE] %*% w
      }
      A <- A[sel]; Bt <- Bt[sel]
      sp[, 2L * q + 1L] <- 2 * (L(2, 0) * A + L(1, 1) * Bt)
      ss[, 2L * q + 1L] <- 2 * (L(1, 1) * A + L(0, 2) * Bt)
      sp[, 2L * q + 2L] <- L(3, 0) * A * x + L(2, 1) * (A * y + Bt * x) + L(1, 2) * Bt * y
      ss[, 2L * q + 2L] <- L(2, 1) * A * x + L(1, 2) * (A * y + Bt * x) + L(0, 3) * Bt * y
    }
    fk <- .sweepSlice(design, full, fixed, fw, k, conds, B, cores)
    if (is.null(fk)) {
      one <- lapply(units, function(u) list(eta = if (length(u$eta) >= k) u$eta[k] else character(0),
                                           conditions = .intersect(u$conditions, conds)))
      fk <- .sweepForward(design, full, fixed, one, cores, conds = conds, blocks = B)
      if (sum(vapply(fk$rows, function(x) length(x$p), 0L)) != sum(sel))
        stop("the data rows of a forward pass along one direction differ from those of the full pass", call. = FALSE)
    }
    back <- .sweepBackward(fk, list(pred = sp, sigma = ss, n = r$n[fw$conds %in% conds]), cores, blocks = B)
    sums <- .unitSums(back, unitOf[fk$conds], fk$names)
    for (u in names(sums)) {
      i <- as.integer(u)
      Su <- sums[[u]]
      rowsOf[[i]][[k]] <- setNames(Su[, 2L], rownames(Su))
      trOf[[i]][[k]] <- Su[, 2L * seq_along(live) + 2L, drop = FALSE]
    }
  }
  lapply(seq_along(units), function(u) {
    eta <- units[[u]]$eta
    rk <- rowsOf[[u]]; tk <- Filter(Negate(is.null), trOf[[u]])
    reach <- unique(c(eta, unlist(lapply(rk, names)), unlist(lapply(tk, rownames))))
    rows <- matrix(0, length(eta), length(reach), dimnames = list(eta, reach))
    for (k in seq_along(rk)) if (!is.null(rk[[k]])) rows[k, names(rk[[k]])] <- rk[[k]]
    TM <- matrix(0, length(reach), nM, dimnames = list(reach, NULL))
    for (t in tk) TM[rownames(t), live] <- TM[rownames(t), live, drop = FALSE] + t
    list(reach = reach, rows = rows, traces = if (length(Ms[[u]])) lapply(seq_len(nM), function(m) setNames(TM[, m], reach)) else list())
  })
}


## The Gauss-Newton Hessian of the data term in the etas of every unit, sum over the unit rows of
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


## The volume term sum_u log|H_u| of FOCEI by forward over reverse: the gradient in `outer` and
## the mode shift d eta*/d theta = -H_ee^-1 H_eA per unit, H_u^-1 from `H_inv_list`, from the
## eta rows of the joint Hessian and the trace gradient at M = H_u^-1 (`.sweepUnits()`). With
## `curvature` also the volume Hessian from dH_u/dtheta, the trace gradients at the unit matrices:
## -tr(H^-1 D_a H H^-1 D_b H) + tr(H^-1 d2P_ab), without the terms of the curvature of the
## prediction in the etas. With `local` the shift and the Hessian of a unit have the columns of the
## parameters it depends on only.
.volumeReverse <- function(design, full, fixed, outer, H_inv_list, cores = 1L, local = FALSE,
                           curvature = FALSE) {
  units <- .volumeUnits(design, full, fixed)
  # per unit the matrices M of the unit trace gradients: H^-1, then E_kl for k <= l
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
  cholIn <- .intersect(outer, .hyperNames(design))
  grad <- setNames(numeric(length(outer)), outer)
  etaShift <- blocks <- vector("list", length(units))
  for (u in seq_along(units)) {
    un <- units[[u]]
    eta <- un$eta; Ku <- length(eta)
    reach <- unique(c(sw[[u]]$reach, cholIn))
    A <- .intersect(outer, reach)
    rowsH <- matrix(0, Ku, length(reach), dimnames = list(eta, reach))
    rowsH[, sw[[u]]$reach] <- sw[[u]]$rows
    # the trace gradients of every M of the unit as columns, with the prior's share
    nM <- length(Ms[[u]])
    TR <- matrix(0, length(reach), nM, dimnames = list(reach, NULL))
    TR[sw[[u]]$reach, ] <- do.call(cbind, sw[[u]]$traces)
    own <- .intersect(cholIn, names(un$dP))
    if (length(own)) {
      Mv <- matrix(vapply(Ms[[u]], c, numeric(Ku * Ku)), Ku * Ku)
      TR[own, ] <- TR[own, , drop = FALSE] + t(crossprod(Mv, matrix(vapply(un$dP[own], c, numeric(Ku * Ku)), Ku * Ku)))
    }
    trace <- TR[, 1L]
    # the normal prior: the precision in the eta block, the derivative in the Cholesky parameters
    rowsH[, eta] <- rowsH[, eta] + un$P
    for (cp in .intersect(cholIn, names(un$dP))) rowsH[, cp] <- rowsH[, cp] + un$dP[[cp]] %*% full[eta]
    Hee <- rowsH[, eta, drop = FALSE]
    eg <- eigen((Hee + t(Hee)) / 2, symmetric = TRUE)
    lam <- pmax(eg$values, 1e-10 * abs(sum(eg$values)) / Ku)
    Em <- -(eg$vectors %*% (t(eg$vectors) / lam)) %*% rowsH[, A, drop = FALSE]
    dimnames(Em) <- list(eta, A)
    grad[A] <- grad[A] + trace[A] + drop(crossprod(Em, trace[eta]))
    etaShift[[u]] <- Em
    if (curvature) {
      # D_a H_kl along the outer parameters with the mode shift, then the volume Hessian
      pq <- pairs[[u]]
      D <- TR[A, 1L + seq_len(nrow(pq)), drop = FALSE] + crossprod(Em, TR[eta, 1L + seq_len(nrow(pq)), drop = FALSE])
      DH <- array(0, c(Ku, Ku, length(A)))
      lower <- matrix(seq_len(Ku * Ku), Ku)
      DH[lower[pq] + rep((seq_along(A) - 1L) * Ku * Ku, each = nrow(pq))] <- t(D)
      DH[lower[pq[, 2:1, drop = FALSE]] + rep((seq_along(A) - 1L) * Ku * Ku, each = nrow(pq))] <- t(D)
      Hinv <- Ms[[u]][[1]]
      P <- array(Hinv %*% matrix(DH, Ku), dim(DH))
      V <- -crossprod(matrix(aperm(P, c(2, 1, 3)), Ku * Ku), matrix(P, Ku * Ku))
      own2 <- .intersect(cholIn, names(un$dP))
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
  fx <- .intersect(names(fixed), outer)
  grad[fx] <- 0
  out <- list(gradient = grad, etaShift = etaShift)
  if (!curvature) return(out)
  if (local) out$hessian <- .sparseFromBlocks(blocks, outer)
  else {
    H <- matrix(0, length(outer), length(outer), dimnames = list(outer, outer))
    for (b in Filter(Negate(is.null), blocks)) H[rownames(b), rownames(b)] <- H[rownames(b), rownames(b)] + b
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
  po <- .intersect(outer, rownames(ph))
  Hoo[po, po] <- Hoo[po, po] + ph[po, po]
  unitOf <- integer(0)
  for (i in seq_along(units)) unitOf[units[[i]]$conditions] <- i
  S <- (Hoo + t(Hoo)) / 2
  for (u in seq_along(units)) {
    eta <- units[[u]]$eta
    cs <- which(unitOf[fw$conds] == u & !is.na(unitOf[fw$conds]))
    He <- t(vapply(seq_along(eta), function(k) .collect(c(eta, outer), rows[cs], 1L + k), numeric(length(eta) + Q)))
    dimnames(He) <- list(eta, c(eta, outer))
    pe <- .intersect(c(eta, outer), colnames(ph))
    He[, pe] <- He[, pe] + ph[eta, pe, drop = FALSE]
    Hee <- He[, eta, drop = FALSE]
    S <- S - crossprod(He[, outer, drop = FALSE], .scaledSolve((Hee + t(Hee)) / 2, He[, outer, drop = FALSE]))
  }
  S
}

