# symmetryReduction(): constructive removal of the non-identifiable directions a
# symmetryDetection() result reports. Scaling directions are gauged exactly through
# the integer weight lattice (transversal pins + invariant monomials); curved
# directions go through module reduction and an escalating exact invariant search
# (monomial -> polynomial <= dPoly -> separable quadratures -> rational with a
# single-coordinate denominator (Laurent ansatz) -> rational via Darboux
# <= dDarboux -> exponential factors <= dExp), every failed stage leaving a
# negative certificate. All arithmetic is exact: integer lattice kernels through
# the Python module, sampling over GF(p) with CRT reconstruction, and an
# independent verify prime for the final X(I) = 0 check.

# ---- coordinates and weight extraction -------------------------------------------

# the full coordinate list of the analysis; results predating info$coordinates fall
# back to the union of the direction supports (the trafo then covers only those)
.symRedCoordinates <- function(object) {
  co <- object$info$coordinates
  if (!is.null(co) && length(co)) return(as.character(co))
  warning("symmetryReduction(): this result carries no coordinate list (older dMod); ",
          "the emitted trafo covers only coordinates appearing in some direction -- ",
          "extend it with identity entries before use in P().", call. = FALSE)
  .symSort(unique(unlist(lapply(object$symmetries, .symCoords))))
}

# The integer weight rows of the scaling directions, parsed PER DIRECTION (unlike
# .symWeights, which drops everything when one weight is symbolic). A scaling with
# symbolic weights (e.g. a Hill recast "-nhill") is routed to the curved machinery.
# Returns W (rows = clean scalings, columns = union of weight names), `rows` (their
# indices into syms) and `symbolic` (indices of scaling directions W cannot hold).
.symRedWeightRows <- function(syms) {
  isScal <- vapply(syms, function(d)
    isTRUE(d$type == "scaling") && !is.null(d$weights), logical(1))
  rows <- integer(0); symbolic <- integer(0); wlist <- list()
  for (i in which(isScal)) {
    w <- suppressWarnings(as.numeric(unlist(syms[[i]]$weights)))
    if (anyNA(w) || any(w != round(w))) symbolic <- c(symbolic, i)
    else { rows <- c(rows, i); wlist[[length(wlist) + 1L]] <-
             setNames(as.integer(w), names(syms[[i]]$weights)) }
  }
  cols <- unique(unlist(lapply(wlist, names)))
  W <- matrix(0L, length(wlist), length(cols), dimnames = list(NULL, cols))
  for (r in seq_along(wlist)) W[r, names(wlist[[r]])] <- wlist[[r]]
  list(W = W, rows = rows, symbolic = symbolic)
}

# ---- small exact-arithmetic utilities --------------------------------------------

# rank over GF(p) (screening; exact verdicts go through the integer kernel)
.symRedRankModP <- function(M, p = .symPrimes[1]) {
  if (!length(M) || nrow(M) == 0L || ncol(M) == 0L) return(0L)
  .symRrefModp(M %% p, p)$rank
}

# primitive integer kernel basis of an integer matrix (columns of the returned
# matrix), exact through the Python module (GF(p) + CRT + integer validation)
.symRedIntKernel <- function(M, sd) {
  if (ncol(M) == 0L) return(matrix(0L, 0L, 0L))
  # as.list per row so a length-1 row still reaches Python as a list, not a scalar;
  # as.numeric, not as.integer -- entries run past 2^31 and Python ints do not
  ker <- sd$exactIntKernel(lapply(seq_len(nrow(M)), function(r)
    as.list(as.numeric(M[r, ]))), ncol(M))
  if (!length(ker)) return(matrix(0L, ncol(M), 0L))
  do.call(cbind, lapply(ker, as.integer))
}

# divide each row by its gcd (primitive integer rows)
.symRedPrimitiveRows <- function(M) {
  if (!nrow(M)) return(M)
  g <- apply(abs(M), 1L, function(r) Reduce(function(a, b) {
    while (b) { t <- b; b <- a %% b; a <- t }; a }, r[r != 0], accumulate = FALSE))
  g[!is.finite(g) | g == 0] <- 1
  round(M / g)
}

# an exponent vector as a monomial string "a*b^-1" ("1" if empty); `^` display
.symRedMonoString <- function(a, vars) {
  nz <- which(a != 0L)
  if (!length(nz)) return("1")
  paste(ifelse(a[nz] == 1L, vars[nz], paste0(vars[nz], "^", a[nz])), collapse = "*")
}

# ---- scaling stage ---------------------------------------------------------------

# Apply the user's `fixed` set to one scaling block: the surviving scaling space is
# the left kernel of W[, S] (combinations c with c^T W vanishing on every fixed
# column). Exact: kernel of t(W[, S]) through the integer-kernel route. Returns the
# residual weight rows (primitive), the removal count and the redundant fixed
# columns (raise no rank, the .symCatFixing peel).
.symRedScalingFixed <- function(W, fixed, sd) {
  S <- intersect(fixed, colnames(W))
  k <- nrow(W)
  if (!length(S)) return(list(Wres = W, removed = 0L, redundant = character(0),
                              combo = diag(1L, k)))
  C <- .symRedIntKernel(t(W[, S, drop = FALSE]), sd)      # k x k' residual combos
  Wres <- if (ncol(C)) .symRedPrimitiveRows(t(C) %*% W) else W[0L, , drop = FALSE]
  red <- character(0); acc <- W[, 0L, drop = FALSE]
  for (cn in S) {
    aug <- cbind(acc, W[, cn])
    if (.symRedRankModP(aug) == .symRedRankModP(acc)) red <- c(red, cn)
    else acc <- aug
  }
  list(Wres = Wres, removed = k - ncol(C), redundant = red, combo = C)
}

# Exact certificate that T is an admissible transversal of W (k rows): A = W[, T]
# invertible with M = A^{-1} W INTEGER (so every survivor's absorbed monomial stays
# Laurent). Integer certification without sympy: round the double solve, then
# verify A %*% M == W exactly in integer arithmetic (entries far below 2^53).
.symRedCertifyTransversal <- function(W, T) {
  A <- W[, T, drop = FALSE]
  if (nrow(A) != length(T)) return(NULL)
  M <- tryCatch(solve(A, W), error = function(e) NULL)
  if (is.null(M)) return(NULL)
  Mr <- round(M)
  if (max(abs(M - Mr)) > 1e-7) return(NULL)
  if (!all(A %*% Mr == W)) return(NULL)
  storage.mode(Mr) <- "integer"
  Mr
}

# The transversal machinery of one scaling block: a certified representative T
# (deterministic: coordinates entering fewest rows first), the admissible family
# (full enumeration when small, else the matroid description "one per RREF row,
# choices independent"), and the survivor meanings z_j * prod_t z_t^{-M_tj}.
.symRedTransversal <- function(Wres, fixed) {
  k <- nrow(Wres)
  supp <- colnames(Wres)[colSums(Wres != 0L) > 0L]
  supp <- setdiff(supp, fixed)                     # fixed columns are already pinned
  out <- list(T = NULL, M = NULL, admissible = NULL, matroid = NULL,
              survivorMeaning = NULL)
  if (!k) return(out)
  Ws <- Wres[, supp, drop = FALSE]
  occ <- colSums(Ws != 0L)
  ord <- supp[order(occ, match(supp, colnames(Wres)))]

  certify <- function(T) .symRedCertifyTransversal(Wres, T)

  enumerate <- length(supp) <= 30L && choose(length(supp), k) <= 200
  if (enumerate) {
    combos <- utils::combn(supp, k, simplify = FALSE)
    adm <- Filter(function(T) .symRedRankModP(Wres[, T, drop = FALSE]) == k, combos)
    out$admissible <- adm
    for (T in adm) {                               # first certified in combn order
      M <- certify(T)
      if (!is.null(M)) { out$T <- T; out$M <- M; break }
    }
  } else {
    # greedy representative; family described through the RREF row supports
    T <- character(0)
    for (cn in ord) {
      if (length(T) == k) break
      if (.symRedRankModP(Wres[, c(T, cn), drop = FALSE]) > length(T)) T <- c(T, cn)
    }
    if (length(T) == k) { M <- certify(T); if (!is.null(M)) { out$T <- T; out$M <- M } }
    ref <- .symRrefModp(Wres[, supp, drop = FALSE] %% .symPrimes[1], .symPrimes[1])
    out$matroid <- lapply(seq_len(ref$rank), function(r)
      supp[which(ref$R[r, ] != 0)])
  }
  if (!is.null(out$T)) {
    surv <- setdiff(supp, out$T)
    out$survivorMeaning <- setNames(vapply(surv, function(j) {
      ex <- setNames(integer(length(out$T) + 1L), c(j, out$T))
      ex[j] <- 1L; ex[out$T] <- -out$M[, j]
      .symRedMonoString(ex, names(ex))
    }, character(1)), surv)
  }
  out
}

# One scaling block end to end: fixed -> residual space -> invariant lattice ->
# transversal -> pins. Status: "fixed" (nothing left), "reduced" (certified
# transversal) or "invariantOnly" (no Laurent transversal; lattice still reported).
.symRedScalingBlock <- function(W, labels, fixed, sd) {
  fx <- .symRedScalingFixed(W, fixed, sd)
  Wres <- fx$Wres
  blk <- list(labels = labels, type = "scaling", support = colnames(W)[colSums(W != 0L) > 0L],
              stage = "transversal", moduleCombos = NULL,
              removedByFixed = fx$removed, redundantFixed = fx$redundant,
              certificates = character(0), reason = NULL)
  if (!nrow(Wres)) {
    blk$status <- "fixed"; blk$stage <- "fixed"
    blk$invariants <- character(0); blk$transversal <- NULL
    blk$gaugeNote <- "removed entirely by the user's fixed coordinates"
    return(blk)
  }
  K <- .symRedIntKernel(Wres[, blk$support, drop = FALSE], sd)
  rownames(K) <- blk$support
  blk$invExps <- K
  blk$invariants <- if (ncol(K)) vapply(seq_len(ncol(K)), function(j)
    .symRedMonoString(K[, j], blk$support), character(1)) else character(0)
  tv <- .symRedTransversal(Wres, fixed)
  blk$transversal <- tv$T
  blk$admissible <- tv$admissible
  blk$matroid <- tv$matroid
  blk$survivorMeaning <- tv$survivorMeaning
  blk$Wres <- Wres
  if (!is.null(tv$T)) {
    blk$status <- "reduced"
    blk$gaugeNote <- paste0("each transversal coordinate may be pinned to any ",
                            "nonzero constant (representative uses 1); the family ",
                            "ranges over all admissible transversals")
    blk$pins <- setNames(rep("1", length(tv$T)), tv$T)
    blk$certificates <- c(blk$certificates,
      "transversal certified exactly: W[,T] invertible with integer W[,T]^-1 W")
  } else {
    blk$status <- "invariantOnly"
    blk$reason <- "no transversal with an integer (Laurent) absorption matrix"
    blk$pins <- NULL
  }
  blk
}

# connected components of the weight rows through shared nonzero columns
.symRedComponents <- function(W) {
  k <- nrow(W)
  if (!k) return(list())
  parent <- seq_len(k)
  find <- function(i) { while (parent[i] != i) i <- parent[i]; i }
  for (i in seq_len(k - 1L)) for (j in (i + 1L):k)
    if (any(W[i, ] != 0L & W[j, ] != 0L)) parent[find(j)] <- find(i)
  split(seq_len(k), vapply(seq_len(k), find, integer(1)))
}

# ---- curved stage: generator prep, module reduction, block decomposition --------

# Locals table for sympify: every identifier mapped to a plain Symbol, so model
# names never collide with sympy builtins (Ci is the cosine integral, E is Euler's
# number, ...). The lookbehind keeps the exponent of "1e-5" from matching; an
# identifier followed by "(" is a function call (exp, log, ...) and must stay
# resolvable, so it never enters the table. The first lookahead pins the match to
# the full identifier so backtracking cannot smuggle a truncated name past the
# call check.
.symRedLocals <- function(strings, spy) {
  ids <- unique(unlist(regmatches(strings,
    gregexpr("(?<![0-9.A-Za-z_])[A-Za-z_][A-Za-z0-9_.]*(?![A-Za-z0-9_.])(?!\\s*\\()",
             strings, perl = TRUE))))
  ids <- setdiff(ids, "")
  reticulate::dict(setNames(lapply(ids, function(x) spy$Symbol(x)), ids))
}

# Sympify memo: every locals table in this file maps every non-call identifier to
# a plain Symbol, so the parse of a given string is the same expression whichever
# table covers it -- a global string-keyed cache is exact. The same strings recur
# across prep, module reduction, Darboux, exp, verify and solve.
.symRedSympifyCache <- new.env(parent = emptyenv())
.symRedSympify <- function(x, spy, locals) {
  key <- as.character(x)
  hit <- .symRedSympifyCache[[key]]
  if (!is.null(hit)) return(hit)
  val <- spy$sympify(key, locals = locals)
  if (length(ls(.symRedSympifyCache, all.names = TRUE)) > 4096L)
    rm(list = ls(.symRedSympifyCache, all.names = TRUE),
       envir = .symRedSympifyCache)
  .symRedSympifyCache[[key]] <- val
  val
}

# the free symbols of a sympy expression, as a character vector
.symRedFreeSyms <- function(e) {
  syms <- reticulate::iterate(e$free_symbols, function(s) as.character(s))
  .symSort(as.character(unlist(syms)))
}

# One curved generator canonicalised for the invariant search: components sympified,
# denominators cleared across the WHOLE generator (Q*X has the same invariants as X
# off Q = 0), re-serialised as polynomial strings. `vars` carries the support plus
# every coefficient symbol (a known input like `u` may enter an invariant).
.symRedGenPrep <- function(d, spy) {
  out <- list(ok = FALSE, comps = NULL, support = NULL, vars = NULL, degree = NA_integer_)
  # $generator arrives in R's power syntax (see .symPublicSymmetry); sympy reads '^'
  # as XOR, so it goes back to '**' before parsing
  comp0 <- gsub("\\^", "**", vapply(d$generator, as.character, character(1)))
  locals <- .symRedLocals(comp0, spy)
  exprs <- tryCatch(lapply(comp0, function(x) .symRedSympify(x, spy, locals)),
                    error = function(e) NULL)
  if (is.null(exprs)) { out$reason <- "component not parseable by sympy"; return(out) }
  dens <- lapply(exprs, function(e) spy$fraction(spy$cancel(spy$together(e)))[[2]])
  L <- Reduce(function(a, b) spy$lcm(a, b), dens, spy$Integer(1L))
  exprs <- lapply(exprs, function(e) spy$expand(spy$cancel(e * L)))
  isPoly <- vapply(exprs, function(e)
    tryCatch(isTRUE(e$is_polynomial()), error = function(err) FALSE), logical(1))
  if (!all(isPoly)) { out$reason <- "component not a rational function"; return(out) }
  comps <- vapply(exprs, function(e) as.character(e), character(1))
  names(comps) <- names(d$generator)
  nz <- comps != "0"
  out$ok <- TRUE
  out$comps <- comps[nz]
  out$support <- names(comps)[nz]
  out$vars <- .symSort(unique(c(out$support,
    unlist(lapply(exprs[nz], .symRedFreeSyms)))))
  degs <- vapply(exprs[nz], function(e) {
    # total_degree may hand back a sympy Integer; the string route converts both
    dg <- tryCatch(suppressWarnings(as.integer(as.character(spy$total_degree(e)))),
                   error = function(err) NA_integer_)
    if (length(dg) != 1L) NA_integer_ else dg
  }, integer(1))
  out$degree <- if (all(is.na(degs))) NA_integer_ else max(degs, na.rm = TRUE)
  out
}

# connected components of index sets through shared elements (union-find); returns
# a list of integer index groups
.symRedGroupBy <- function(sets) {
  k <- length(sets)
  if (!k) return(list())
  parent <- seq_len(k)
  find <- function(i) { while (parent[i] != i) i <- parent[i]; i }
  for (i in seq_len(k - 1L)) for (j in (i + 1L):k)
    if (length(intersect(sets[[i]], sets[[j]]))) parent[find(j)] <- find(i)
  unname(split(seq_len(k), vapply(seq_len(k), find, integer(1))))
}

# Pivot-schedule discovery for the module reduction: RREF-style elimination on the
# value matrix over GF(p), at each step choosing the (row, column) whose elimination
# leaves the fewest nonzeros. A Markowitz-style count cannot serve here: the point of
# the reduction is to discover symbolic cancellations, which appear only as numeric
# zeros AFTER a trial elimination. The fill delta is therefore computed exactly, but
# only over the affected rows and without matrix copies; a pivot in a singleton
# column eliminates nothing and is skipped. Returns the schedule and its final fill.
.symRedReduceSchedule <- function(V, p) {
  k <- nrow(V); n <- ncol(V)
  used <- logical(k); sched <- list()
  for (step in seq_len(k)) {
    nzR <- rowSums(V != 0); nzC <- colSums(V != 0)
    best <- NULL
    for (r in which(!used)) for (cc in which(V[r, ] != 0)) {
      if (nzC[cc] < 2) next
      pinv <- .symInvmod(V[r, cc], p)
      delta <- 0
      for (j in which(V[, cc] != 0)) {
        if (j == r) next
        m <- .symRedMulmodP(V[j, cc], pinv, p)
        delta <- delta + sum((V[j, ] - .symRedMulmodP(V[r, ], m, p)) %% p != 0) - nzR[j]
      }
      if (is.null(best) || delta < best$delta) best <- list(r = r, c = cc, delta = delta)
    }
    if (is.null(best)) break
    r <- best$r; cc <- best$c
    piv <- V[r, cc]
    for (j in seq_len(k)) if (j != r && V[j, cc] != 0) {
      m <- .symRedMulmodP(V[j, cc], .symInvmod(piv, p), p)
      V[j, ] <- (V[j, ] - .symRedMulmodP(V[r, ], m, p)) %% p
    }
    used[r] <- TRUE
    sched[[length(sched) + 1L]] <- c(r, cc)
  }
  list(sched = sched, fill = sum(V != 0))
}

# vectorised a*b mod p in doubles (wrapper around the scalar-b .symMulmod)
.symRedMulmodP <- function(a, b, p) .symMulmod(a %% p, b %% p, p)

# Module reduction of a curved block: generators may be combined with FUNCTION
# coefficients (pointwise span decides identifiability -- a module over the function
# field, not a vector space), so RREF with symbolic multipliers is legitimate and
# can collapse the support (M009: X13 - r2*X12 went from 7 to 3 coordinates). The
# pivot schedule is discovered on a value matrix over GF(p) (cheap, fill-minimising,
# cross-checked on a second prime), then replayed exactly in sympy. The heuristic
# affects only how small the support gets, never correctness.
.symRedModuleReduce <- function(preps, labels, sd, spy) {
  k <- length(preps)
  cols <- .symSort(unique(unlist(lapply(preps, `[[`, "support"))))
  vars <- .symSort(unique(unlist(lapply(preps, `[[`, "vars"))))
  if (k < 2L)
    return(list(preps = preps, labelsOf = as.list(seq_len(k)), combos = character(0)))

  # value matrices at one random point per prime; schedule = the lower-fill one
  has <- lapply(seq_len(k), function(r) which(cols %in% names(preps[[r]]$comps)))
  scheds <- lapply(.symPrimes[1:2], function(p) {
    pool <- .symPool()
    pt <- setNames(pool(seq_along(vars) + 7L), vars)
    exprs <- unlist(lapply(seq_len(k), function(r)
      vapply(cols[has[[r]]], function(cc) preps[[r]]$comps[[cc]], character(1))))
    vals <- .symRedEvalBatch(exprs, pt, p, sd)
    vals[is.na(vals)] <- 0
    V <- matrix(0, k, length(cols))
    off <- 0L
    for (r in seq_len(k)) {
      nv <- length(has[[r]])
      if (nv) { V[r, has[[r]]] <- vals[off + seq_len(nv)]; off <- off + nv }
    }
    .symRedReduceSchedule(V, p)
  })
  sched <- scheds[[which.min(vapply(scheds, `[[`, numeric(1), "fill"))]]$sched

  # exact symbolic replay of the schedule, batched into one Python call (the
  # per-entry cancel() loop through reticulate dominated M009-scale blocks);
  # 'combo' tracks the combination coefficients over the original rows
  rr <- sd$moduleReduceReplay(lapply(preps, function(pr) as.list(pr$comps)),
                              as.list(cols), lapply(sched, as.list))

  # Candidates = the original generators (identity combination) plus every reduced
  # row; a greedy pass picks a rank-complete basis with the smallest supports, so a
  # row is only ever replaced when the replacement's support is genuinely smaller.
  cand <- lapply(seq_len(k), function(j)
    list(pr = preps[[j]], srcIdx = j, combo = NULL))
  for (j in seq_len(k)) {
    comps <- unlist(rr$rows[[j]])                      # nonzero entries only
    if (!length(comps)) next                           # row eliminated entirely
    pr <- .symRedGenPrep(list(generator = as.list(comps)), spy)
    if (!pr$ok || !length(pr$support)) next
    Crow <- unlist(rr$combo[[j]])
    srcIdx <- which(Crow != "0")
    if (length(srcIdx) == 1L && srcIdx == j && identical(Crow[[j]], "1"))
      next                                             # unchanged: already in
    combo <- paste0(
      paste0("(", Crow[srcIdx], ")*", labels[srcIdx], collapse = " + "),
      "  [support ", length(pr$support), "]")
    cand[[length(cand) + 1L]] <- list(pr = pr, srcIdx = srcIdx, combo = combo)
  }

  p <- .symPrimes[1]
  pool <- .symPool()
  pt <- setNames(pool(seq_along(vars) + 23L), vars)
  vrow <- function(pr) {
    out <- numeric(length(cols))
    hs <- which(cols %in% names(pr$comps))
    if (length(hs)) {
      v <- .symRedEvalBatch(vapply(cols[hs], function(cc) pr$comps[[cc]],
                                   character(1)), pt, p, sd)
      v[is.na(v)] <- 0
      out[hs] <- v
    }
    out
  }
  ord <- order(vapply(cand, function(cd) length(cd$pr$support), integer(1)))
  acc <- matrix(0, 0L, length(cols))
  outPreps <- list(); labelsOf <- list(); combos <- character(0)
  for (j in ord) {
    if (nrow(acc) == k) break
    aug <- rbind(acc, vrow(cand[[j]]$pr))
    if (.symRedRankModP(aug, p) > nrow(acc)) {
      acc <- aug
      outPreps[[length(outPreps) + 1L]] <- cand[[j]]$pr
      labelsOf[[length(labelsOf) + 1L]] <- cand[[j]]$srcIdx
      if (!is.null(cand[[j]]$combo)) combos <- c(combos, cand[[j]]$combo)
    }
  }
  list(preps = outPreps, labelsOf = labelsOf, combos = combos)
}

# ---- invariant stages ------------------------------------------------------------

# Batch evaluation of rational expression strings at one integer point mod p: one
# Python round trip for the whole vector (the scalar .symEvalModq re-parses per
# call, which dominates the sampling loops at M009 scale). NA where the
# denominator vanishes mod p or the evaluation fails.
.symRedEvalBatch <- function(exprs, pt, p, sd) {
  if (!length(exprs)) return(numeric(0))
  v <- tryCatch(sd$evalRationalModBatch(as.list(as.character(exprs)),
                                        as.list(names(pt)),
                                        list(as.list(as.numeric(pt))),
                                        as.integer(p)),
                error = function(e) NULL)
  if (is.null(v)) return(rep(NA_real_, length(exprs)))
  out <- as.numeric(unlist(v[[1]]))
  out[out < 0] <- NA_real_
  out
}

# The log-derivative rows of a block at one random point mod p: one row per
# generator, entries eta_i = xi_i(z)/z_i over the moved coordinates. NULL when an
# evaluation fails (redraw the point).
.symRedEtaRows <- function(preps, vars, pt, p, sd) {
  idx <- lapply(preps, function(pr) which(vars %in% names(pr$comps)))
  exprs <- unlist(lapply(seq_along(preps), function(g)
    vapply(vars[idx[[g]]], function(v) preps[[g]]$comps[[v]], character(1))))
  vals <- .symRedEvalBatch(exprs, pt, p, sd)
  if (anyNA(vals)) return(NULL)
  rows <- matrix(0, length(preps), length(vars))
  off <- 0L
  for (g in seq_along(preps)) {
    for (k in seq_along(idx[[g]])) {
      i <- idx[[g]][k]
      rows[g, i] <- .symMulmod(vals[off + k], .symInvmod(pt[[vars[i]]] %% p, p), p)
    }
    off <- off + length(idx[[g]])
  }
  rows
}

# Common nullspace of row-sampling matrices over the 4 primes, CRT-reconstructed to
# primitive integer columns. `rowsAt(pt, p)` delivers the rows of one point (NULL =
# redraw). Sampling saturates when the GF(p) rank is flat twice. Returns the integer
# basis (columns), or NULL when pivots disagree across every prime pair or an entry
# does not reconstruct, and always the row count for the certificate text.
.symRedModularNullspace <- function(rowsAt, n, evalVars, primes = .symPrimes) {
  basisPer <- list(); pivotsPer <- list(); nRows <- 0L
  nv <- length(evalVars)
  for (pi in seq_along(primes)) {
    p <- primes[pi]
    pool <- .symPool()
    rows <- matrix(0, 0L, n); prev <- -1L; flat <- 0L
    off <- 50L * pi; draws <- 0L
    # The rank is checked once per CHUNK of drawn rows, on the carried-forward
    # reduced rows plus the chunk (the compiled RREF never re-reduces old rows):
    # a full re-reduction of the accumulated matrix per drawn row is quadratic in
    # the sample count and dominated the larger ansatz stages. Chunk 1 for small
    # systems keeps the draw-exact saturation of the scalar path.
    chunk <- max(1L, min(16L, n %/% 8L))
    red <- matrix(0, 0L, n); pend <- matrix(0, 0L, n); r <- 0L
    repeat {
      pt <- setNames(pool(seq_len(nv) + off), evalVars); off <- off + nv
      draws <- draws + 1L
      new <- rowsAt(pt, p)
      if (is.null(new)) { if (draws > 25L) break else next }
      rows <- rbind(rows, new)
      pend <- rbind(pend, new)
      if (nrow(pend) < chunk && nrow(rows) <= 4L * n + 20L) next
      rr <- .symRrefModp(rbind(red, pend) %% p, p)
      red <- rr$R[seq_len(rr$rank), , drop = FALSE]
      pend <- matrix(0, 0L, n)
      r <- rr$rank
      if (r == n) break
      flat <- if (r == prev) flat + 1L else 0L
      prev <- r
      if (flat >= 2L || nrow(rows) > 4L * n + 20L) break
    }
    nRows <- max(nRows, nrow(rows))
    if (!nrow(rows)) return(list(basis = NULL, rows = 0L))   # sampling never succeeded
    ref <- .symRrefModp(rows %% p, p)
    free <- setdiff(0:(n - 1L), ref$piv)
    if (!length(free))
      return(list(basis = matrix(0L, n, 0L), rows = nRows))
    pivotsPer[[pi]] <- ref$piv
    basisPer[[pi]] <- .symNullspaceBasis(list(dim = n, pivots = ref$piv, R = ref$R),
                                         free, p)
  }
  agree <- vapply(pivotsPer, function(pv) identical(pv, pivotsPer[[1]]), logical(1))
  use <- which(agree)
  if (length(use) < 2L) return(list(basis = NULL, rows = nRows))
  nf <- ncol(basisPer[[use[1]]])
  residues <- vapply(use, function(jj) as.integer(basisPer[[jj]]),
                     integer(n * nf))
  rec <- symRatRecon(matrix(residues, n * nf, length(use)), as.integer(primes[use]))
  if (any(rec$den == "0")) return(list(basis = NULL, rows = nRows))
  num <- suppressWarnings(as.numeric(rec$num)); den <- suppressWarnings(as.numeric(rec$den))
  if (anyNA(num) || anyNA(den) || any(abs(num) > 2^40) || any(den > 2^40))
    return(list(basis = NULL, rows = nRows))
  B <- matrix(num / den, n, nf)
  for (j in seq_len(nf)) {                    # clear each column to primitive integers
    dj <- den[(j - 1L) * n + seq_len(n)]
    L <- Reduce(function(a, b) a * b / .symRedGcd(a, b), unique(dj), 1)
    B[, j] <- B[, j] * L
  }
  B <- round(B)
  storage.mode(B) <- "integer"
  list(basis = t(.symRedPrimitiveRows(t(B))), rows = nRows)
}

.symRedGcd <- function(a, b) { a <- abs(a); b <- abs(b)
  while (b) { t <- b; b <- a %% b; a <- t }; a }

# Stage 1: Laurent-monomial invariants. m = prod z_i^{a_i} is invariant under X iff
# sum_i a_i * xi_i/z_i vanishes identically -- linear in the exponents a, one row per
# (generator, point). Exponents run over the MOVED coordinates only: a coefficient
# symbol (xi identically 0) would only multiply in a trivial constant factor.
.symRedMonomialInvariants <- function(preps, sd) {
  vars <- .symSort(unique(unlist(lapply(preps, `[[`, "support"))))
  evalVars <- .symSort(unique(c(vars, unlist(lapply(preps, `[[`, "vars")))))
  n <- length(vars)
  ns <- .symRedModularNullspace(function(pt, p)
    .symRedEtaRows(preps, vars, pt, p, sd), n, evalVars)
  cert <- sprintf(paste0("monomial stage: exact modular nullspace over %d primes, ",
                         "%d sampling rows"), length(.symPrimes), ns$rows)
  if (is.null(ns$basis))
    return(list(ok = FALSE, invariants = character(0), exps = NULL,
                cert = paste0(cert, " -- reconstruction failed (inconclusive)")))
  if (!ncol(ns$basis))
    return(list(ok = FALSE, invariants = character(0), exps = NULL,
                cert = paste0("no Laurent-monomial invariant (", cert, ")")))
  inv <- vapply(seq_len(ncol(ns$basis)), function(j)
    .symRedMonoString(ns$basis[, j], vars), character(1))
  list(ok = TRUE, invariants = inv, exps = ns$basis, vars = vars, cert = cert)
}

# an integer coefficient vector over a monomial table as a polynomial string
.symRedPolyString <- function(cf, expts, vars) {
  nz <- which(cf != 0L)
  terms <- vapply(nz, function(k) {
    mono <- .symRedMonoString(expts[k, ], vars)
    a <- cf[k]
    if (mono == "1") as.character(a)
    else if (a == 1L) mono
    else if (a == -1L) paste0("-", mono)
    else paste0(a, "*", mono)
  }, character(1))
  out <- terms[1]
  for (t in terms[-1]) out <- paste0(out, if (startsWith(t, "-")) " - " else " + ",
                                     sub("^-", "", t))
  out
}

# Stage 2: polynomial invariants of total degree <= d. For I = sum_k c_k m_k the
# condition X(I) = 0 is LINEAR in c: X(m_k) = m_k * sum_i a^(k)_i xi_i/z_i, so a row
# per (generator, point) has entries m_k(z) * <a^(k), eta(z)> mod p. The ansatz runs
# over the block support plus the coefficient symbols, but monomials with no moved
# coordinate are dropped up front: X kills them identically, so they contribute only
# the trivial invariants (the observed 2/5/9 pattern) and excluding them is lossless
# up to additive constants.
.symRedPolyInvariants <- function(preps, dPoly, sd) {
  moved <- .symSort(unique(unlist(lapply(preps, `[[`, "support"))))
  vars <- .symSort(unique(c(moved, unlist(lapply(preps, `[[`, "vars")))))
  expts <- .symMonoTable(length(vars), as.integer(dPoly))
  keep <- rowSums(expts[, match(moved, vars), drop = FALSE]) > 0L
  expts <- expts[keep, , drop = FALSE]
  N <- nrow(expts)
  if (!N) return(list(ok = FALSE, invariants = character(0),
                      cert = "polynomial stage skipped: empty ansatz"))
  rowsAt <- function(pt, p) {
    eta <- .symRedEtaRows(preps, vars, pt, p, sd)          # nGen x nvars
    if (is.null(eta)) return(NULL)
    mv <- symMonoResidues(expts, as.integer(pt %% p), p)   # monomial values
    # exponents are tiny and eta < p < 2^31, so the inner products stay far below
    # 2^53 and the plain matrix product is exact; only the outer multiply needs
    # the split-multiply (elementwise-safe for equal-length vectors)
    inner <- (expts %*% t(eta)) %% p                       # N x nGen
    t(vapply(seq_len(nrow(eta)), function(g)
      .symRedMulmodP(mv, inner[, g], p), numeric(N)))
  }
  ns <- .symRedModularNullspace(rowsAt, N, vars)
  cert <- sprintf(paste0("polynomial stage (degree <= %d, %d monomials): exact ",
                         "modular nullspace over %d primes, %d sampling rows"),
                  dPoly, N, length(.symPrimes), ns$rows)
  if (is.null(ns$basis))
    return(list(ok = FALSE, invariants = character(0),
                cert = paste0(cert, " -- reconstruction failed (inconclusive)")))
  if (!ncol(ns$basis))
    return(list(ok = FALSE, invariants = character(0),
                cert = paste0("no polynomial invariant of total degree <= ", dPoly,
                              " beyond functions of the unmoved coordinates (",
                              cert, ")")))
  sel <- .symRedIndependentPoly(ns$basis, expts, vars, moved, preps)
  inv <- vapply(sel, function(j)
    .symRedPolyString(ns$basis[, j], expts, vars), character(1))
  list(ok = TRUE, invariants = inv,
       coefs = ns$basis[, sel, drop = FALSE], expts = expts, vars = vars,
       cert = paste0(cert, sprintf("; %d functionally independent of %d",
                                   length(sel), ncol(ns$basis))))
}

# The linear nullspace over-generates: powers and products of low-degree invariants
# reappear as independent basis vectors. Select a functionally independent subset by
# a greedy Jacobian rank test at a generic point, lowest degree / fewest terms
# first; the target count is #moved - rank(xi), the dimension of the invariant
# foliation. Selection only -- correctness rests on the verify layer.
.symRedIndependentPoly <- function(B, expts, vars, moved, preps) {
  pt <- setNames(as.numeric(.symPool()(seq_along(vars) + 5L)), vars)
  Xi <- t(vapply(preps, function(pr) vapply(moved, function(v)
    if (v %in% names(pr$comps))
      eval(parse(text = pr$comps[[v]]), as.list(pt)) else 0,
    numeric(1)), numeric(length(moved))))
  target <- length(moved) - qr(Xi, tol = 1e-9)$rank
  mv <- exp(expts %*% log(pt))                       # monomial values at pt
  # The gradient spans the MOVED coordinates only, matching what `target` counts:
  # an unmoved coordinate is itself a trivial invariant, so two invariants that
  # differ only in unmoved directions do not separate orbits. Scoring them as
  # independent fills the quota with a set that leaves moved coordinates
  # untouched -- and the solve then has fewer carriers than invariants.
  mIdx <- match(moved, vars)
  grad <- function(cf) vapply(mIdx, function(i)
    sum(cf * mv * expts[, i] / pt[i]), numeric(1))
  deg <- apply(expts, 1L, sum)
  key <- vapply(seq_len(ncol(B)), function(j)
    max(deg[B[, j] != 0L]) * 1e4 + sum(B[, j] != 0L), numeric(1))
  J <- matrix(0, 0L, length(moved)); sel <- integer(0)
  for (j in order(key)) {
    if (length(sel) >= target && target > 0L) break
    Jc <- rbind(J, grad(B[, j]))
    if (qr(Jc, tol = 1e-9)$rank > nrow(J)) { J <- Jc; sel <- c(sel, j) }
  }
  if (!length(sel)) sel <- seq_len(min(1L, ncol(B)))
  sel
}

# unknown-count budget of the rational stage's Laurent ansatz
.symRedRationalCap <- 700L

# Stage 2b: rational invariants with a single-coordinate denominator. The ansatz
# is Laurent: the polynomial monomials of total degree <= dPoly plus every such
# monomial with ONE moved coordinate at exponent -1. The invariance condition
# stays LINEAR in the coefficients -- X(m_k) = m_k * <a^(k), eta> holds for
# Laurent monomials too, and symMonoResidues takes negative exponents through the
# modular inverse -- so the same exact modular-nullspace sampling applies. This is
# the cheap stage that reaches sums like (z^2*a + z*a + b*c)/z: the factor stages
# see them only as high-degree Darboux polynomials (never a factor of anything
# visible), the exp stage only at a raised numerator cap. Denominators are pruned
# to the coordinates with z_i | xi_i for every generator: for N/z_i in lowest
# terms X(N/z_i) = 0 forces z_i | N*xi_i, and z_i is irreducible. The screen is a
# numeric divisibility test (xi_i vanishes on z_i = 0), selection-only -- the
# verify layer re-proves every invariant symbolically.
.symRedRationalInvariants <- function(preps, dPoly, sd, spy) {
  moved <- .symSort(unique(unlist(lapply(preps, `[[`, "support"))))
  vars <- .symSort(unique(c(moved, unlist(lapply(preps, `[[`, "vars")))))
  skip <- function(cert) list(ok = FALSE, invariants = character(0), cert = cert)
  if (dPoly < 1L) return(skip("rational stage skipped (dPoly = 0)"))
  pool <- .symPool()
  divisible <- vapply(moved, function(v) {
    xs <- unlist(lapply(preps, function(pr)
      if (v %in% names(pr$comps)) pr$comps[[v]]))
    if (!length(xs)) return(TRUE)                    # xi identically 0
    all(vapply(1:2, function(t) {
      pt <- setNames(pool(seq_along(vars) + 31L * t), vars)
      pt[v] <- 0
      val <- .symRedEvalBatch(xs, pt, .symPrimes[1], sd)
      !anyNA(val) && all(val == 0)
    }, logical(1)))
  }, logical(1))
  denC <- moved[divisible]
  if (!length(denC))
    return(skip(paste0("rational stage skipped: no moved coordinate divides its ",
                       "own generator component (no admissible denominator)")))
  movedIdx <- match(moved, vars)
  T0 <- .symMonoTable(length(vars), as.integer(dPoly))
  keep0 <- rowSums(T0[, movedIdx, drop = FALSE]) > 0L
  laur <- lapply(match(denC, vars), function(i) {
    Ti <- T0[T0[, i] == 0L, , drop = FALSE]
    Ti[, i] <- -1L
    Ti
  })
  expts <- rbind(T0[keep0, , drop = FALSE], do.call(rbind, laur))
  N <- nrow(expts)
  if (N > .symRedRationalCap)
    return(skip(sprintf(paste0("rational stage skipped (%d Laurent monomials ",
                               "over the %d-unknown budget)"), N,
                        .symRedRationalCap)))
  rowsAt <- function(pt, p) {
    eta <- .symRedEtaRows(preps, vars, pt, p, sd)
    if (is.null(eta)) return(NULL)
    mv <- symMonoResidues(expts, as.integer(pt %% p), p)
    # entries of the inner product are bounded by dPoly*p < 2^36, exact in doubles
    inner <- (expts %*% t(eta)) %% p
    t(vapply(seq_len(nrow(eta)), function(g)
      .symRedMulmodP(mv, inner[, g], p), numeric(N)))
  }
  ns <- .symRedModularNullspace(rowsAt, N, vars)
  cert <- sprintf(paste0("rational stage (numerator degree <= %d, denominator ",
                         "one of {%s}; %d Laurent monomials): exact modular ",
                         "nullspace over %d primes, %d sampling rows"),
                  dPoly, paste(denC, collapse = ", "), N, length(.symPrimes),
                  ns$rows)
  if (is.null(ns$basis))
    return(skip(paste0(cert, " -- reconstruction failed (inconclusive)")))
  if (!ncol(ns$basis))
    return(skip(paste0("no rational invariant with a single-coordinate ",
                       "denominator and numerator degree <= ", dPoly, " (",
                       cert, ")")))
  sel <- .symRedIndependentPoly(ns$basis, expts, vars, moved, preps)
  inv <- vapply(sel, function(j)
    .symRedPolyString(ns$basis[, j], expts, vars), character(1))
  # canonical fraction form N/z instead of the raw Laurent sum: the readable
  # shape -- and the joint carrier solve is orders of magnitude faster on it
  # (sympy's solve grinds on the sum-with-embedded-quotient form)
  locals <- .symRedLocals(inv, spy)
  inv <- vapply(inv, function(s) {
    e <- tryCatch(spy$cancel(spy$together(
      .symRedSympify(gsub("\\^", "**", s), spy, locals))),
      error = function(err) NULL)
    if (is.null(e)) s else gsub("\\*\\*", "^", as.character(e))
  }, character(1), USE.NAMES = FALSE)
  list(ok = TRUE, invariants = inv,
       cert = paste0(cert, sprintf("; %d functionally independent of %d",
                                   length(sel), ncol(ns$basis))))
}

# symbolic determinant of the extactic matrix stays tractable only for a small
# monomial basis; above the cap the stage falls back to coordinate factors
.symRedExtacticCap <- 12L

# apply the generator X = sum_i xi_i d/dz_i of one prep to a sympy expression
.symRedApplyX <- function(pr, xiExprs, f, spy, locals) {
  out <- spy$Integer(0L)
  for (v in names(xiExprs))
    out <- out + xiExprs[[v]] * spy$diff(f, .symRedSympify(v, spy, locals))
  spy$expand(out)
}

# Exact rational value of a sympy expression at an integer point, as "num/den".
# Substitution pairs use explicit Symbol objects: string keys would re-sympify and
# collide with sympy builtins (Ci, E, ...), silently skipping the substitution.
.symRedSubsRational <- function(e, pt, spy) {
  pairs <- lapply(names(pt), function(n)
    reticulate::tuple(spy$Symbol(n), spy$Integer(as.integer(pt[[n]]))))
  as.character(spy$cancel(e$subs(pairs)))
}

# Exact integer sampling rows of a cofactor matrix: one row per (generator, point),
# column j evaluating cols[[j]][[g]] (sympy expressions or strings, rational in
# `vars`; values through the exact-Fraction batch evaluator, one Python call per
# generator). Each row is lcm-cleared to integers; rows with a failed evaluation
# or entries >= 2^31 are dropped. Serves the Darboux kernel [lambda], the exp
# stage [lambda | mu] and the integrating-factor solve [lambda | mu | div X].
.symRedCofactorRows <- function(cols, vars, nGen, nMon, sd) {
  rows <- list(); perGen <- integer(nGen)
  for (g in seq_len(nGen)) {
    exprs <- vapply(cols, function(cl) as.character(cl[[g]]), character(1))
    # The sample index must not walk the prime pool outwards: a cofactor of degree
    # d turns the point magnitude into its d-th power, so a stride of 101 put every
    # row of a high-degree generator past the storage bound, where the guard below
    # dropped it -- silently, leaving a rank-deficient matrix whose kernel is not
    # the cofactor kernel at all. A stride of 2 keeps the points distinct and small.
    pts <- lapply(seq_len(nMon), function(t)
      as.list(.symPool()(seq_along(vars) + 2L * t + 7L * g)))
    v <- tryCatch(sd$evalRationalBatch(as.list(exprs), as.list(vars), pts),
                  error = function(e) NULL)
    if (is.null(v)) next
    for (t in seq_len(nMon)) {
      num <- suppressWarnings(as.numeric(unlist(v$num[[t]])))
      den <- suppressWarnings(as.numeric(unlist(v$den[[t]])))
      if (anyNA(num) || anyNA(den) || any(den == 0)) next
      L <- Reduce(function(a, b) a * b / .symRedGcd(a, b), unique(den), 1)
      row <- round(num * (L / den))
      # 2^53 is where a double stops holding integers exactly -- the real bound,
      # and the kernel takes them from here as doubles
      if (all(is.finite(row)) && max(abs(row)) < 2^53) {
        rows[[length(rows) + 1L]] <- row
        perGen[g] <- perGen[g] + 1L
      }
    }
  }
  # A generator that contributed no row is simply absent from the matrix, and its
  # kernel constraints with it: report nothing rather than a kernel of the wrong
  # system (the caller turns an empty return into an inconclusive certificate).
  if (any(perGen == 0L)) return(list())
  rows
}

# Stage 3: rational invariants via Darboux polynomials (extactic route). Every
# Darboux polynomial P (X(P) = lambda*P, lambda the polynomial cofactor) of degree
# <= d divides the extactic determinant E_d of the monomial basis B_d (Pereira), so
# spy$factor(E_d) delivers ALL candidates of that degree; candidates must divide for
# every generator of the block. Rational invariants are prod P_j^{n_j} for integer
# vectors n with sum n_j lambda_j = 0 identically -- the integer kernel of the
# cofactor matrix, sampled exactly at integer points and lcm-cleared. This subsumes
# equal-cofactor pairs and the coordinate factors z_i | xi_i. Above the basis cap
# only the coordinate factors are tried (stated in the certificate).
.symRedDarboux <- function(preps, dDarboux, sd, spy, verbose = FALSE,
                           extacticOK = TRUE) {
  moved <- .symSort(unique(unlist(lapply(preps, `[[`, "support"))))
  vars <- .symSort(unique(c(moved, unlist(lapply(preps, `[[`, "vars")))))
  locals <- .symRedLocals(c(unlist(lapply(preps, `[[`, "comps")), vars), spy)
  xiOf <- lapply(preps, function(pr)
    lapply(pr$comps, function(x) .symRedSympify(x, spy, locals)))

  # Candidate Darboux polynomials from three sources: the coordinate factors, the
  # irreducible factors of the xi components themselves (an invariant hypersurface
  # often sits inside the vanishing of xi -- and when a first integral makes the
  # extactic degenerate, these are the candidates that remain reachable), and the
  # extactic factors when the basis stays below the cap. The exact division test
  # filters, so extra candidates never cost correctness.
  cands <- lapply(moved, function(v) .symRedSympify(v, spy, locals))
  for (xi in unlist(xiOf, recursive = FALSE)) {
    fl <- tryCatch(spy$factor_list(xi), error = function(e) NULL)
    if (!is.null(fl)) cands <- c(cands, lapply(fl[[2]], function(pair) pair[[1]]))
  }
  # The extactic basis runs over the MOVED coordinates only, with every
  # coefficient symbol treated as a constant of the ground field: Pereira's
  # theorem holds over any characteristic-0 field, so completeness is "factor
  # degree <= d in the moved coordinates, arbitrary degree in the parameters" --
  # a stronger statement than the all-variables basis, at a fraction of the size
  # (C(nMoved + d, d) instead of C(nvars + d, d)); module-reduced blocks with a
  # handful of moved coordinates now pass the cap that used to skip them.
  Nbasis <- choose(length(moved) + dDarboux, dDarboux)
  # the symbolic determinant's cost is driven by the entry degrees, which grow to
  # ~ dDarboux + (N-1)*(D-1) in the last extactic row -- cap that, not just N
  Dmax <- max(1L, vapply(preps, function(pr)
    if (is.na(pr$degree)) 3L else pr$degree, integer(1)))
  entryDeg <- dDarboux + (Nbasis - 1L) * (Dmax - 1L)
  extactic <- extacticOK && dDarboux >= 1L && Nbasis <= .symRedExtacticCap &&
              entryDeg <= .symRedExtacticCap
  degenerate <- FALSE
  if (extactic) {
    expts <- .symMonoTable(length(moved), as.integer(dDarboux))
    basis <- lapply(seq_len(nrow(expts)), function(k)
      .symRedSympify(.symRedMonoString(expts[k, ], moved), spy, locals))
    N <- length(basis)
    rows <- list(basis)
    for (j in seq_len(N - 1L))
      rows[[j + 1L]] <- lapply(rows[[j]], function(f)
        .symRedApplyX(preps[[1]], xiOf[[1]], f, spy, locals))
    Emat <- spy$Matrix(lapply(rows, function(r) r))
    D <- Emat$det(method = "berkowitz")               # the one expensive step
    nonzero <- any(vapply(1:2, function(t) {          # E identically 0 = degenerate
      pt <- setNames(as.list(.symPool()(seq_along(vars) + 11L * t)), vars)
      v <- tryCatch(.symRedSubsRational(D, pt, spy), error = function(e) NA_character_)
      !is.na(v) && v != "0"
    }, logical(1)))
    if (!nonzero) degenerate <- TRUE
    else {
      fl <- spy$factor_list(D)                      # convert=TRUE: already an R list
      cands <- c(cands, lapply(fl[[2]], function(pair) pair[[1]]))
    }
  }

  # keep candidates that are Darboux for EVERY generator, with their cofactors and
  # proportionality dedup -- batched into one Python call
  dc <- sd$darbouxCofactors(lapply(cands, as.character),
                            lapply(preps, function(pr) as.list(pr$comps)))
  Ps <- lapply(dc$Ps, function(s) .symRedSympify(s, spy, locals))
  lams <- lapply(dc$lams, function(lr)
    lapply(lr, function(s) .symRedSympify(s, spy, locals)))

  coverage <- if (degenerate)
      "extactic degenerate (E identically 0): coordinate and xi factors only"
    else if (extactic)
      sprintf(paste0("extactic complete for factor degree <= %d in the moved ",
                     "coordinates"), dDarboux)
    else if (!extacticOK)
      "extactic skipped (invariant count already met): coordinate and xi factors only"
    else sprintf(paste0("extactic skipped (moved-coordinate basis %d, projected ",
                        "entry degree %d, cap %d): coordinate and xi factors only"),
                 Nbasis, entryDeg, .symRedExtacticCap)
  certBase <- sprintf("Darboux stage (%d candidate factors; %s)", length(Ps), coverage)
  # Ps/lams/coverage travel on every exit path: the exp stage consumes the factor
  # basis and its cofactors even when no purely rational combination exists here.
  fail <- function(cert) list(ok = FALSE, invariants = character(0), cert = cert,
                              Ps = Ps, lams = lams, coverage = coverage)
  if (length(Ps) < 2L)
    return(fail(paste0("no rational invariant from Darboux factors -- ", certBase)))

  # integer kernel of the cofactor matrix, sampled exactly at integer points
  nMon <- choose(length(vars) + max(1L, dDarboux), length(vars)) + 3L
  rows <- .symRedCofactorRows(lams, vars, length(preps), nMon, sd)
  if (!length(rows))
    return(fail(paste0("cofactor sampling failed (inconclusive) -- ", certBase)))
  K <- .symRedIntKernel(do.call(rbind, rows), sd)
  if (!ncol(K))
    return(fail(paste0("no rational invariant with Darboux factors of degree <= ",
                       dDarboux, " -- ", certBase)))
  Pstr <- vapply(Ps, function(P) gsub("\\*\\*", "^", as.character(P)), character(1))
  inv <- vapply(seq_len(ncol(K)), function(j) {
    nz <- which(K[, j] != 0L)
    paste(vapply(nz, function(l) {
      base <- if (grepl("[+*/ -]", Pstr[l])) paste0("(", Pstr[l], ")") else Pstr[l]
      if (K[l, j] == 1L) base else paste0(base, "^", K[l, j])
    }, character(1)), collapse = "*")
  }, character(1))
  list(ok = TRUE, invariants = inv, cert = certBase,
       Ps = Ps, lams = lams, coverage = coverage)
}

# denominator-candidate cap of the exp stage (h = 1, all P_j, then all P_j^2)
.symRedExpDenCap <- 13L
# unknown-count budget per exp-stage linear system (g- plus mu-coefficients)
.symRedExpSizeCap <- 250L

# Stage 4: exponential factors E = exp(g/h) with polynomial cofactor mu = X(g/h)
# (Prelle-Singer: every Liouvillian first integral is of Darboux form once these
# are admitted). For a fixed denominator h = prod P_j^{k_j} built from the Darboux
# factors (cofactor Lambda_h = sum k_j lambda_j) the condition
#   X(g) - g*Lambda_h = mu*h
# is LINEAR in the coefficients of g (total degree <= dExp, shared across the
# block's generators) and of the per-generator cofactors mu (structural bound
# deg mu <= dExp + D - 1 - deg h), so the same exact modular-nullspace sampling as
# the polynomial stage applies. Every candidate is then certified symbolically
# (the defining identity must cancel() to 0 for every generator), and invariants
# combine Darboux and exponential factors through the integer kernel of the
# extended cofactor matrix [lambda | mu]:
#   I = prod P_j^{n_j} * exp(sum m_k g_k/h_k),  sum n lambda + sum m mu = 0.
.symRedExpFactors <- function(preps, darb, dExp, sd, spy, verbose = FALSE) {
  moved <- .symSort(unique(unlist(lapply(preps, `[[`, "support"))))
  vars <- .symSort(unique(c(moved, unlist(lapply(preps, `[[`, "vars")))))
  K <- length(preps)
  D <- max(1L, vapply(preps, function(pr)
    if (is.na(pr$degree)) 1L else pr$degree, integer(1)))
  locals <- .symRedLocals(c(unlist(lapply(preps, `[[`, "comps")), vars), spy)
  xiOf <- lapply(preps, function(pr)
    lapply(pr$comps, function(x) .symRedSympify(x, spy, locals)))
  zero <- spy$Integer(0L)

  # denominator candidates: 1, then every Darboux factor, then their squares (so
  # the cap trims the squares first)
  hC <- list(list(h = spy$Integer(1L), cof = rep(list(zero), K), deg = 0L))
  hC2 <- list()
  for (j in seq_along(darb$Ps)) {
    P <- darb$Ps[[j]]
    dj <- tryCatch(as.integer(as.character(spy$total_degree(P))),
                   error = function(e) NA_integer_)
    if (is.na(dj) || dj < 1L) next
    hC[[length(hC) + 1L]] <- list(h = P, cof = darb$lams[[j]], deg = dj)
    hC2[[length(hC2) + 1L]] <- list(h = spy$expand(P * P),
      cof = lapply(darb$lams[[j]], function(l) spy$Integer(2L) * l), deg = 2L * dj)
  }
  hC <- c(hC, hC2)
  capped <- length(hC) > .symRedExpDenCap
  if (capped) hC <- hC[seq_len(.symRedExpDenCap)]

  # g-ansatz: moved monomials plus the constant (exp(c/h) is a genuine factor)
  exptsG <- .symMonoTable(length(vars), as.integer(dExp))
  keepG <- rowSums(exptsG[, match(moved, vars), drop = FALSE]) > 0L |
    rowSums(exptsG) == 0L
  exptsG <- exptsG[keepG, , drop = FALSE]
  Ng <- nrow(exptsG)

  factors <- list()
  seen <- character(0)
  sizeSkipped <- 0L
  for (hc in hC) {
    # deg mu <= D - 1 is the DEFINING bound of an exponential factor (Christopher/
    # Llibre), not just a structural one -- and it is lossless for the combination
    # step: higher-degree cofactors could only cancel among themselves, and the
    # corresponding g-sum lies in this ansatz with a small cofactor already. It
    # also keeps the h = 1 system from being solved by every g.
    degMu <- min(D - 1L, dExp + D - 1L - hc$deg)
    exptsMu <- if (degMu >= 0L) .symMonoTable(length(vars), degMu) else
      matrix(0L, 0L, length(vars))
    Nmu <- nrow(exptsMu)
    n <- Ng + K * Nmu
    if (n > .symRedExpSizeCap) { sizeSkipped <- sizeSkipped + 1L; next }
    hStr <- as.character(hc$h)
    lamStr <- vapply(hc$cof, as.character, character(1))
    rowsAt <- function(pt, p) {
      eta <- .symRedEtaRows(preps, vars, pt, p, sd)
      if (is.null(eta)) return(NULL)
      hl <- .symRedEvalBatch(c(hStr, lamStr), pt, p, sd)
      if (anyNA(hl)) return(NULL)
      mvG <- symMonoResidues(exptsG, as.integer(pt %% p), p)
      mvMu <- if (Nmu) symMonoResidues(exptsMu, as.integer(pt %% p), p) else numeric(0)
      inner <- (exptsG %*% t(eta)) %% p                 # Ng x K
      muNeg <- if (Nmu) (p - .symRedMulmodP(mvMu, hl[1], p)) %% p else numeric(0)
      t(vapply(seq_len(K), function(g) {
        gcols <- (.symRedMulmodP(mvG, inner[, g], p) -
                    .symRedMulmodP(mvG, hl[1L + g], p)) %% p
        mucols <- numeric(K * Nmu)
        if (Nmu) mucols[(g - 1L) * Nmu + seq_len(Nmu)] <- muNeg
        c(gcols, mucols)
      }, numeric(n)))
    }
    ns <- .symRedModularNullspace(rowsAt, n, vars)
    if (is.null(ns$basis) || !ncol(ns$basis)) next
    for (jc in seq_len(ncol(ns$basis))) {
      gvec <- ns$basis[seq_len(Ng), jc]
      if (all(gvec == 0L)) next
      gExpr <- .symRedSympify(.symRedPolyString(gvec, exptsG, vars), spy, locals)
      q <- spy$cancel(gExpr / hc$h)
      if (isTRUE(q$is_number)) next                     # exp(const)
      if (hc$deg > 0L &&
          isTRUE(tryCatch(q$is_polynomial(), error = function(e) FALSE)))
        next                                            # h | g: the h = 1 run owns it
      mu <- lapply(seq_len(K), function(g) {
        if (!Nmu) return(zero)
        mv <- ns$basis[Ng + (g - 1L) * Nmu + seq_len(Nmu), jc]
        if (all(mv == 0L)) zero else
          .symRedSympify(.symRedPolyString(mv, exptsMu, vars), spy, locals)
      })
      # exact certification of the defining identity, per generator
      okAll <- all(vapply(seq_len(K), function(g) {
        Xg <- .symRedApplyX(preps[[g]], xiOf[[g]], gExpr, spy, locals)
        identical(as.character(spy$cancel(spy$together(
          Xg - gExpr * hc$cof[[g]] - mu[[g]] * hc$h))), "0")
      }, logical(1)))
      if (!okAll) next
      key <- as.character(q)
      if (key %in% seen) next
      seen <- c(seen, key)
      factors[[length(factors) + 1L]] <- list(g = gExpr, h = hc$h, mu = mu)
    }
  }

  certBase <- sprintf(paste0("exp stage (numerator degree <= %d over %d ",
                             "denominator candidate(s)%s%s): %d exponential ",
                             "factor(s)"),
                      dExp, length(hC), if (capped) ", capped" else "",
                      if (sizeSkipped) sprintf(", %d skipped over the %d-unknown budget",
                                               sizeSkipped, .symRedExpSizeCap) else "",
                      length(factors))
  if (!length(factors))
    return(list(ok = FALSE, invariants = character(0), factors = list(),
                cert = paste0("no exponential factor with numerator degree <= ",
                              dExp, " -- ", certBase)))
  if (isTRUE(verbose)) message("  exp: ", length(factors), " factor(s) certified")

  # integer kernel of the extended cofactor matrix [lambda | mu]
  J <- length(darb$Ps)
  cols <- c(darb$lams, lapply(factors, `[[`, "mu"))
  nMon <- choose(length(vars) + max(1L, dExp + D - 1L), length(vars)) + 3L
  rows <- .symRedCofactorRows(cols, vars, K, nMon, sd)
  if (!length(rows))
    return(list(ok = FALSE, invariants = character(0), factors = factors,
                cert = paste0("cofactor sampling failed (inconclusive) -- ", certBase)))
  Kk <- .symRedIntKernel(do.call(rbind, rows), sd)
  useCol <- if (ncol(Kk)) which(vapply(seq_len(ncol(Kk)), function(j)
    any(Kk[J + seq_along(factors), j] != 0L), logical(1))) else integer(0)
  if (!length(useCol))
    return(list(ok = FALSE, invariants = character(0), factors = factors,
                cert = paste0("no invariant combining the exponential factors -- ",
                              certBase)))
  Pstr <- vapply(darb$Ps, function(P) gsub("\\*\\*", "^", as.character(P)),
                 character(1))
  inv <- vapply(useCol, function(j) {
    w <- Kk[, j]
    parts <- vapply(which(w[seq_len(J)] != 0L), function(l) {
      base <- if (grepl("[+*/ -]", Pstr[l])) paste0("(", Pstr[l], ")") else Pstr[l]
      if (w[l] == 1L) base else paste0(base, "^", w[l])
    }, character(1))
    arg <- spy$cancel(spy$together(Reduce(`+`,
      lapply(which(w[J + seq_along(factors)] != 0L), function(e2)
        spy$Integer(w[J + e2]) * factors[[e2]]$g / factors[[e2]]$h))))
    argStr <- gsub("\\*\\*", "^", as.character(arg))
    # no Darboux part => sum m*mu = 0 and the argument is itself a rational
    # invariant; unwrap only when its sign is certified (exp(q) covers both
    # signs of an indefinite q, a bare positive carrier does not)
    sArg <- if (!length(parts)) .symRedSgn(arg, spy) else 0L
    if (sArg == 1L) argStr
    else if (sArg == -1L) gsub("\\*\\*", "^", as.character(spy$cancel(-arg)))
    else paste(c(parts, paste0("exp(", argStr, ")")), collapse = "*")
  }, character(1))
  list(ok = TRUE, invariants = inv, cert = certBase, factors = factors)
}

# Stage 5 (single generator, exactly two moved coordinates): integrating factor
# in Darboux form. Singer: a Liouvillian first integral exists iff an integrating
# factor M = prod P^n * exp(sum m g/h) does, with RATIONAL exponents solving the
# inhomogeneous cofactor equation sum n lambda + sum m mu = -div X. That is one
# extra column in the existing sampling kernel: integer vectors of
# [lambda | mu | div X] with a nonzero last entry, normalised by it. The first
# integral itself comes from sympy quadrature of dI = M*(xi2 dz1 - xi1 dz2) and
# is accepted only when X(I) = 0 checks exactly.
.symRedIntFactor <- function(preps, darb, exFactors, sd, spy, verbose = FALSE) {
  moved <- .symSort(unique(unlist(lapply(preps, `[[`, "support"))))
  vars <- .symSort(unique(c(moved, unlist(lapply(preps, `[[`, "vars")))))
  locals <- .symRedLocals(c(unlist(lapply(preps, `[[`, "comps")), vars), spy)
  xi <- lapply(preps[[1]]$comps[moved], function(x) .symRedSympify(x, spy, locals))
  divX <- spy$expand(Reduce(`+`, lapply(moved, function(v)
    spy$diff(xi[[v]], spy$Symbol(v)))))
  J <- length(darb$Ps); E <- length(exFactors)
  cols <- c(darb$lams, lapply(exFactors, `[[`, "mu"), list(list(divX)))
  nMon <- choose(length(vars) + 2L, length(vars)) + 3L
  rows <- .symRedCofactorRows(cols, vars, 1L, nMon, sd)
  certBase <- sprintf("integrating-factor stage (%d Darboux + %d exponential factor(s))",
                      J, E)
  if (!length(rows))
    return(list(ok = FALSE, invariants = character(0),
                cert = paste0("cofactor sampling failed (inconclusive) -- ", certBase)))
  Kk <- .symRedIntKernel(do.call(rbind, rows), sd)
  last <- J + E + 1L
  useCol <- if (ncol(Kk)) which(Kk[last, ] != 0L) else integer(0)
  if (!length(useCol))
    return(list(ok = FALSE, invariants = character(0),
                cert = paste0("no Darboux-form integrating factor -- ", certBase)))
  Pstr <- vapply(darb$Ps, function(P) gsub("\\*\\*", "^", as.character(P)),
                 character(1))
  inv <- character(0)
  for (j in useCol) {
    w <- Kk[, j]; wL <- w[last]
    parts <- vapply(which(w[seq_len(J)] != 0L), function(l) {
      base <- if (grepl("[+*/ -]", Pstr[l])) paste0("(", Pstr[l], ")") else Pstr[l]
      paste0(base, "^(", w[l], "/", wL, ")")
    }, character(1))
    expPart <- if (E && any(w[J + seq_len(E)] != 0L)) {
      arg <- Reduce(`+`, lapply(which(w[J + seq_len(E)] != 0L), function(e2)
        spy$Rational(w[J + e2], wL) * exFactors[[e2]]$g / exFactors[[e2]]$h))
      paste0("exp(", as.character(spy$cancel(spy$together(arg))), ")")
    } else character(0)
    Mstr <- paste(c(parts, expPart), collapse = "*")
    if (!nchar(Mstr)) Mstr <- "1"
    Iv <- tryCatch(sd$integratingFactorIntegral(
      gsub("\\^", "**", Mstr), as.character(xi[[moved[1]]]),
      as.character(xi[[moved[2]]]), moved[1], moved[2]),
      error = function(e) NULL)
    if (!is.null(Iv)) {
      inv <- c(inv, gsub("\\*\\*", "^", Iv))
      break
    }
  }
  if (!length(inv))
    return(list(ok = FALSE, invariants = character(0),
                cert = paste0("integrating factor found, quadrature failed -- ",
                              certBase)))
  if (isTRUE(verbose)) message("  intfactor: first integral by quadrature")
  list(ok = TRUE, invariants = inv, cert = paste0(certBase, ": solved by quadrature"))
}

# The invariant count a block needs: #moved coordinates minus the generic rank of
# the xi-matrix -- the dimension of the invariant foliation. Numeric, at a generic
# point; selection-only (correctness rests on the verify layer).
# Stage: separable characteristics. When every component of a generator involves no
# MOVED coordinate but its own, the characteristic system dz_i/xi_i = dz_j/xi_j
# decouples and the first integrals are n-1 one-dimensional quadratures
# G_i = int dz_i / xi_i, with the invariants G_i - G_j (exponentiated, so a pair of
# logarithms comes back as the rational quotient). No ansatz and no degree cap, and
# it reaches antiderivatives -- atan, log of a factored denominator -- that no
# product of Darboux factors can express; the Laurent-monomial stage is the special
# case xi_i = w_i z_i. A separable generator's own integrals need NOT be invariants
# of a multi-generator block, so every candidate is checked against every generator
# before it is offered.
.symRedSeparable <- function(preps, sd, spy) {
  moved <- .symSort(unique(unlist(lapply(preps, `[[`, "support"))))
  locals <- .symRedLocals(c(unlist(lapply(preps, `[[`, "comps")), moved), spy)
  xiOf <- lapply(preps, function(pr)
    lapply(pr$comps, function(x) .symRedSympify(x, spy, locals)))
  sepIdx <- which(vapply(seq_along(preps), function(g) {
    nm <- names(preps[[g]]$comps)
    length(nm) >= 2L && all(vapply(nm, function(v)
      !length(setdiff(intersect(.symRedFreeSyms(xiOf[[g]][[v]]), moved), v)),
      logical(1)))
  }, logical(1)))
  if (!length(sepIdx))
    return(list(ok = FALSE, invariants = character(0),
                cert = "no generator has a separable characteristic system"))
  cand <- character(0)
  for (g in sepIdx) {
    nm <- names(preps[[g]]$comps)
    # spy$Integer(1L), not R's 1: a double reaches sympy as a float and turns the
    # whole antiderivative into 1.0*atan(1.0*x)
    G <- lapply(nm, function(v) tryCatch(
      spy$integrate(spy$Integer(1L) / xiOf[[g]][[v]],
                    .symRedSympify(v, spy, locals)),
      error = function(e) NULL))
    names(G) <- nm
    ok <- !vapply(G, is.null, logical(1))
    if (sum(ok) < 2L) next
    nm <- nm[ok]
    for (j in seq_along(nm)[-1]) {
      # G_i - G_j is the invariant. Exponentiating is only worth it when the
      # antiderivatives are logarithms, where it collapses log(a) - log(b) into
      # the rational quotient; an atan or rational difference is left as it is.
      D <- G[[nm[1]]] - G[[nm[j]]]
      e <- tryCatch(
        if (length(.symRedIter(D$atoms(spy$log), function(a) a)))
          spy$simplify(spy$exp(D)) else spy$cancel(spy$together(D)),
        error = function(err) NULL)
      if (!is.null(e)) cand <- c(cand, gsub("\\*\\*", "^", as.character(e)))
    }
  }
  cand <- unique(cand[nzchar(cand) & cand != "1"])
  cert <- sprintf(paste0("separable-characteristics stage (%d of %d generator(s) ",
                         "separable, %d quadrature candidate(s))"),
                  length(sepIdx), length(preps), length(cand))
  if (!length(cand))
    return(list(ok = FALSE, invariants = character(0),
                cert = paste0("no closed-form quadrature -- ", cert)))
  # exact X(I) = 0 against EVERY generator; the same checker the verify layer uses
  keep <- tryCatch(as.logical(unlist(sd$verifyInvariants(
    as.list(gsub("\\^", "**", cand)),
    lapply(preps, function(pr) as.list(pr$comps)))), use.names = FALSE),
    error = function(e) NULL)
  if (is.null(keep) || !any(keep))
    return(list(ok = FALSE, invariants = character(0),
                cert = paste0("quadrature integrals are not invariants of every ",
                              "generator -- ", cert)))
  list(ok = TRUE, invariants = cand[keep],
       cert = paste0(cert, sprintf(": %d verified", sum(keep))))
}


.symRedBlockCorank <- function(preps) {
  moved <- .symSort(unique(unlist(lapply(preps, `[[`, "support"))))
  vars <- .symSort(unique(c(moved, unlist(lapply(preps, `[[`, "vars")))))
  pt <- setNames(as.numeric(.symPool()(seq_along(vars) + 5L)), vars)
  Xi <- t(vapply(preps, function(pr) vapply(moved, function(v)
    if (v %in% names(pr$comps))
      eval(parse(text = pr$comps[[v]]), as.list(pt)) else 0,
    numeric(1)), numeric(length(moved))))
  length(moved) - qr(Xi, tol = 1e-9)$rank
}

# a functionally independent subset of invariant strings (any stage's output, all
# R-parseable), greedy by string length under a numeric Jacobian rank test;
# keepOrder = TRUE takes the caller's ordering instead
.symRedIndependentSet <- function(invs, vars, target, moved = vars,
                                  keepOrder = FALSE) {
  if (!length(invs) || target <= 0L) return(character(0))
  pool <- as.numeric(.symPool()(seq_along(vars) + 13L))
  # exp-carrying invariants overflow at prime-sized points; small distinct values
  # keep every gradient finite without losing genericity (rank only)
  if (any(grepl("exp(", invs, fixed = TRUE))) pool <- 1 + pool / (max(pool) + 1)
  pt <- setNames(pool, vars)
  # The point spans every variable (an unmoved coefficient still has to evaluate),
  # but the gradient spans the MOVED coordinates only -- that is what `target`
  # counts. An unmoved coordinate is itself a trivial invariant, so two invariants
  # differing only in unmoved directions do not separate orbits, and scoring them
  # as independent fills the quota with a set that leaves moved coordinates free.
  grad <- function(iv) {
    ex <- parse(text = iv)[[1]]
    vapply(moved, function(v)
      tryCatch(eval(stats::D(ex, v), as.list(pt)), error = function(e) NA_real_),
      numeric(1))
  }
  J <- matrix(0, 0L, length(moved)); sel <- character(0)
  cand <- unique(invs)
  if (!keepOrder) cand <- cand[order(nchar(cand))]
  for (iv in cand) {
    if (length(sel) >= target) break
    gr <- grad(iv)
    if (!all(is.finite(gr))) next
    Jc <- rbind(J, gr)
    if (qr(Jc, tol = 1e-9)$rank > nrow(J)) { J <- Jc; sel <- c(sel, iv) }
  }
  sel
}

# run the invariant-search stages on one (already module-reduced) sub-block of
# generators -- monomial, then polynomial <= dPoly, then Darboux <= dDarboux,
# then exponential factors with numerator degree <= dExp.
.symRedSolveBlock <- function(preps, dPoly, dDarboux, dExp, sd, spy,
                              verbose = FALSE, separable = TRUE) {
  moved <- .symSort(unique(unlist(lapply(preps, `[[`, "support"))))
  vars <- .symSort(unique(c(moved, unlist(lapply(preps, `[[`, "vars")))))
  target <- .symRedBlockCorank(preps)
  certs <- character(0); found <- character(0); stage <- "monomial"
  # an invariant free of every moved coordinate is trivially constant on the
  # orbits and must never count toward the target (a Darboux factor with
  # cofactor 0, e.g. a pure parameter sum, would otherwise fake a complete
  # invariant set and an exactness claim); monomial/polynomial stages exclude
  # such invariants by construction, the factor stages cannot
  movedRe <- paste0("(?<![0-9A-Za-z_.])(", paste(moved, collapse = "|"),
                    ")(?![0-9A-Za-z_.])")
  # sign-definite invariants first: an indefinite one cannot be carried by a
  # positive outer parameter, so a definite alternative must win the slot
  sgnCache <- new.env(parent = emptyenv())
  definite <- function(iv) {
    hit <- sgnCache[[iv]]
    if (!is.null(hit)) return(hit)
    e <- tryCatch(.symRedSympify(gsub("\\^", "**", iv), spy,
                                 .symRedLocals(iv, spy)),
                  error = function(err) NULL)
    v <- !is.null(e) && .symRedSgn(e, spy) != 0L
    sgnCache[[iv]] <- v
    v
  }
  # tiers: definite rational, indefinite rational, transcendental
  pick <- function(cand) {
    cand <- unique(cand[grepl(movedRe, cand, perl = TRUE)])
    isExp <- grepl("exp(", cand, fixed = TRUE)
    def <- !isExp & vapply(cand, definite, logical(1))
    tier <- ifelse(def, 1L, ifelse(!isExp, 2L, 3L))
    .symRedIndependentSet(cand[order(tier, nchar(cand))], vars, target, moved,
                          keepOrder = TRUE)
  }
  # the stages ACCUMULATE, and a full set escalates further while it still
  # holds an indefinite invariant: a definite replacement may live one stage up
  short <- function() length(found) < target ||
    !all(vapply(found, definite, logical(1)))
  # A later stage's reshuffle of the set is always ACCEPTED -- pick's preference
  # for short representatives is load-bearing (the monomial lattice basis can
  # come back as EGF_EGFR^2/EGFR where a later stage offers EGF_EGFR*k_bind, and
  # the joint carrier solve grinds for minutes on the former) -- but the stage
  # LABEL only advances when the stage genuinely contributed: more invariants,
  # or more definite ones.
  better <- function(got) length(got) > length(found) ||
    (length(got) == length(found) &&
       sum(vapply(got, definite, logical(1))) >
       sum(vapply(found, definite, logical(1))))

  mono <- .symRedMonomialInvariants(preps, sd)
  certs <- c(certs, mono$cert)
  if (mono$ok) found <- pick(mono$invariants)
  if (short()) {
    if (isTRUE(verbose)) message("  monomial: ", length(found), "/", target,
                                 "; escalating to degree <= ", dPoly)
    poly <- .symRedPolyInvariants(preps, dPoly, sd)
    certs <- c(certs, poly$cert)
    if (poly$ok) {
      got <- pick(c(found, poly$invariants))
      if (!identical(got, found)) {
        if (better(got)) stage <- "polynomial"
        found <- got
      }
    }
  }
  # Separability is a syntactic property of the generator plus n univariate
  # quadratures -- far cheaper than the factor stages below, so it runs before
  # them. It sits AFTER the two ansatz stages because those return the closed form
  # a reader wants for the cases they do cover (a linear invariant as k1 + u*k2,
  # not as exp(-(k1 + u*k2)/u)); what is left over is where the quadrature earns
  # its place, including antiderivatives outside the Darboux language (atan).
  if (separable && short()) {
    sep <- .symRedSeparable(preps, sd, spy)
    certs <- c(certs, sep$cert)
    if (sep$ok) {
      got <- pick(c(found, sep$invariants))
      if (!identical(got, found)) {
        if (better(got)) stage <- "separable"
        found <- got
      }
    }
  } else if (!separable)
    certs <- c(certs, "separable-characteristics stage skipped (separable = FALSE)")
  if (short()) {
    if (isTRUE(verbose)) message("  polynomial: ", length(found), "/", target,
                                 "; escalating to rational")
    rat <- .symRedRationalInvariants(preps, dPoly, sd, spy)
    certs <- c(certs, rat$cert)
    if (rat$ok) {
      got <- pick(c(found, rat$invariants))
      if (!identical(got, found)) {
        if (better(got)) stage <- "rational"
        found <- got
      }
    }
  }
  dar <- NULL
  if (short()) {
    if (isTRUE(verbose)) message("  rational: ", length(found), "/", target,
                                 "; escalating to Darboux")
    dar <- .symRedDarboux(preps, dDarboux, sd, spy, verbose,
                          extacticOK = length(found) < target)
    certs <- c(certs, dar$cert)
    if (dar$ok) {
      got <- pick(c(found, dar$invariants))
      if (!identical(got, found)) {
        if (better(got)) stage <- "darboux"
        found <- got
      }
    }
  }
  if (short() && !is.null(dar)) {
    if (dExp > 0L) {
      if (isTRUE(verbose)) message("  darboux: ", length(found), "/", target,
                                   "; escalating to exponential factors")
      ex <- .symRedExpFactors(preps, dar, dExp, sd, spy, verbose)
      certs <- c(certs, ex$cert)
      if (ex$ok) {
        got <- pick(c(found, ex$invariants))
        if (!identical(got, found)) {
        if (better(got)) stage <- "exp"
        found <- got
      }
      }
      if (length(found) < target && length(preps) == 1L &&
          length(unique(unlist(lapply(preps, `[[`, "support")))) == 2L) {
        itf <- .symRedIntFactor(preps, dar, ex$factors, sd, spy, verbose)
        certs <- c(certs, itf$cert)
        # a 2-coordinate single-generator block has target 1 and `found` is empty
        # here, so the quadrature integral (already proven X(I) = 0 in Python) is
        # appended directly -- pick()'s numeric gradient cannot evaluate the
        # atan/Abs forms this stage may produce
        if (itf$ok) { found <- c(found, itf$invariants[1]); stage <- "intfactor" }
      }
    } else
      certs <- c(certs, "exp stage skipped (dExp = 0)")
  }
  list(stage = stage, invariants = found, target = target, certificates = certs,
       reason = if (length(found) >= target) NULL
         else sprintf(paste0("found %d of %d independent invariants (poly/",
                             "rational degree <= %d, Darboux degree <= %d, exp ",
                             "numerator degree <= %d); gauge pinning needs the ",
                             "full set"),
                      length(found), target, dPoly, dDarboux, dExp))
}

# Fixed-coordinate elimination for one curved group: over the function field the
# admissible generators are the combinations whose xi vanishes on every fixed
# coordinate (a fixed coordinate cannot move). Exact symbolic Gaussian elimination
# on the fixed columns; each pivot row cannot combine away from the fixed
# coordinates and is removed by the fixing -- the curved mirror of
# .symRedScalingFixed. Kept rows are exactly zero on every fixed column (the
# elimination is symbolic), so the fixed coordinates leave their support and can
# never become carriers or gauge pins. `sources` tracks which original directions
# enter each kept combination.
.symRedFixedReduce <- function(preps, fixed, spy) {
  k <- length(preps)
  cols <- .symSort(unique(unlist(lapply(preps, `[[`, "support"))))
  fixedCols <- intersect(fixed, cols)
  if (!length(fixedCols) || !k)
    return(list(preps = preps, sources = as.list(seq_len(k)),
                removed = 0L, fixedCols = character(0)))
  locals <- .symRedLocals(unlist(lapply(preps, `[[`, "comps")), spy)
  E <- lapply(preps, function(pr) {
    row <- setNames(vector("list", length(cols)), cols)
    for (cc in cols) row[[cc]] <- .symRedSympify(
      if (cc %in% names(pr$comps)) pr$comps[[cc]] else "0", spy, locals)
    row
  })
  sources <- as.list(seq_len(k))
  zero <- function(e) identical(as.character(e), "0")
  usedPivot <- logical(k)
  for (fc in fixedCols) {
    piv <- NA_integer_
    for (r in seq_len(k)) if (!usedPivot[r] && !zero(E[[r]][[fc]])) {
      piv <- r; break
    }
    if (is.na(piv)) next
    for (j in seq_len(k)) {
      if (j == piv || zero(E[[j]][[fc]])) next
      m <- spy$cancel(E[[j]][[fc]] / E[[piv]][[fc]])
      for (c2 in cols) E[[j]][[c2]] <- spy$cancel(E[[j]][[c2]] - m * E[[piv]][[c2]])
      sources[[j]] <- sort(unique(c(sources[[j]], sources[[piv]])))
    }
    usedPivot[piv] <- TRUE
  }
  outP <- list(); outS <- list()
  for (j in which(!usedPivot)) {
    comps <- vapply(cols, function(cc) as.character(E[[j]][[cc]]), character(1))
    names(comps) <- cols
    pr <- .symRedGenPrep(list(generator = as.list(comps[comps != "0"])), spy)
    if (!pr$ok || !length(pr$support)) next          # row eliminated entirely
    outP[[length(outP) + 1L]] <- pr
    outS[[length(outS) + 1L]] <- sources[[j]]
  }
  list(preps = outP, sources = outS, removed = sum(usedPivot),
       removedIdx = which(usedPivot), fixedCols = fixedCols)
}

# The curved pipeline: prep every direction with a closed-form generator, group by
# shared support, fixed-eliminate and module-reduce each group, re-split (a
# reduced row may decouple), then run the invariant stages per sub-block.
# Support-only directions and failed preps become unresolved singleton blocks.
.symRedCurved <- function(syms, idx, labels, fixed, dPoly, dDarboux, dExp, sd,
                          spy, verbose = FALSE, separable = TRUE) {
  blocks <- list()
  emit <- function(b) blocks[[length(blocks) + 1L]] <<- b
  unresolved <- function(labs, support, reason)
    list(labels = labs, type = "curved", support = support, stage = "none",
         status = "unresolved", invariants = character(0),
         certificates = character(0), reason = reason,
         moduleCombos = NULL, transversal = NULL, survivorMeaning = NULL,
         pins = NULL)

  # a block whose directions are all scalings (demoted here because their support
  # touches a curved one) is still named for what it holds, not for the machinery
  scalLabs <- labels[vapply(syms, function(d) isTRUE(d$type == "scaling"), logical(1))]
  kindOf <- function(labs)
    if (length(labs) && all(labs %in% scalLabs)) "scaling" else "curved"

  hasGen <- vapply(idx, function(i) !is.null(syms[[i]]$generator), logical(1))
  for (i in idx[!hasGen])
    emit(unresolved(labels[i], .symCoords(syms[[i]]),
      "no closed-form generator; re-run symmetryDetection(reconstruct = TRUE)"))

  gi <- idx[hasGen]
  preps <- lapply(gi, function(i) .symRedGenPrep(syms[[i]], spy))
  bad <- !vapply(preps, `[[`, logical(1), "ok")
  for (w in which(bad))
    emit(unresolved(labels[gi[w]], .symCoords(syms[[gi[w]]]), preps[[w]]$reason))
  gi <- gi[!bad]; preps <- preps[!bad]
  if (!length(gi)) return(blocks)

  for (grp in .symRedGroupBy(lapply(preps, `[[`, "support"))) {
    glab <- labels[gi[grp]]
    fx <- .symRedFixedReduce(preps[grp], fixed, spy)
    fixCert <- if (fx$removed > 0L)
      sprintf(paste0("fixed-coordinate elimination: %d of %d direction(s) cannot ",
                     "combine away from {%s}"), fx$removed, length(grp),
              paste(fx$fixedCols, collapse = ", "))
    if (fx$removed > 0L && isTRUE(verbose))
      message("curved block {", paste(glab, collapse = ", "), "}: fixing removed ",
              fx$removed, " direction(s)")
    # pivot directions that survive nowhere get their own "fixed" block, so the
    # direction bookkeeping (removed/remaining) stays complete
    fixedLabs <- if (fx$removed > 0L)
      setdiff(glab[fx$removedIdx], glab[unlist(fx$sources)]) else character(0)
    if (length(fixedLabs) || !length(fx$preps)) {
      fl <- if (length(fx$preps)) fixedLabs else glab
      emit(list(labels = fl, type = "curved", kind = kindOf(fl),
                support = .symSort(unique(unlist(lapply(preps[grp], `[[`, "support")))),
                stage = "fixed", status = "fixed",
                removedByFixed = fx$removed, invariants = character(0),
                certificates = fixCert, reason = NULL, moduleCombos = NULL,
                transversal = NULL, survivorMeaning = NULL, pins = NULL,
                gaugeNote = "removed entirely by the user's fixed coordinates"))
    }
    if (!length(fx$preps)) next
    glabK <- vapply(fx$sources, function(s)
      if (length(s) == 1L) glab[s]
      else paste0("(", paste(glab[s], collapse = "+"), ")"), character(1))
    if (isTRUE(verbose))
      message("curved block {", paste(glab, collapse = ", "), "}: module reduction")
    mr <- .symRedModuleReduce(fx$preps, glabK, sd, spy)
    for (sub in .symRedGroupBy(lapply(mr$preps, `[[`, "support"))) {
      subLabs <- .symSort(unique(glab[unlist(fx$sources[unlist(mr$labelsOf[sub])])]))
      sol <- .symRedSolveBlock(mr$preps[sub], dPoly, dDarboux, dExp, sd, spy,
                               verbose, separable)
      emit(list(labels = subLabs, type = "curved", kind = kindOf(subLabs),
                support = .symSort(unique(unlist(lapply(mr$preps[sub], `[[`, "support")))),
                stage = sol$stage, target = sol$target,
                status = if (length(sol$invariants)) "invariantOnly" else "unresolved",
                removedByFixed = if (fx$removed > 0L && !length(fixedLabs))
                  fx$removed,
                invariants = sol$invariants,
                certificates = c(fixCert, sol$certificates),
                reason = sol$reason, moduleCombos = mr$combos,
                transversal = NULL, survivorMeaning = NULL, pins = NULL,
                preps = mr$preps[sub]))
    }
  }
  blocks
}

# ---- verification and trafo assembly ---------------------------------------------

# Exact internal verification. Scaling blocks: the invariant exponents annihilate
# the residual weights over the integers (by construction; asserted). Curved blocks:
# X(I) = 0 proven symbolically -- the Lie derivative of every invariant along every
# generator must cancel() to literally "0", a full proof, not a sample. A failing
# invariant is dropped and the failure recorded (this guards our own algebra; it
# should never fire).
.symRedVerify <- function(blocks, sd) {
  for (bi in seq_along(blocks)) {
    b <- blocks[[bi]]
    if (identical(b$type, "scaling")) {
      if (!is.null(b$Wres) && !is.null(b$invExps) && ncol(b$invExps) &&
          nrow(b$Wres)) {
        prod <- b$Wres[, rownames(b$invExps), drop = FALSE] %*% b$invExps
        if (any(prod != 0))
          stop("symmetryReduction(): internal error -- scaling invariant not in the ",
               "weight kernel.", call. = FALSE)
        blocks[[bi]]$certificates <- c(b$certificates,
          "verified: integer weight annihilation, exact")
      }
      next
    }
    if (!length(b$invariants) || is.null(b$preps)) next
    keep <- as.logical(unlist(tryCatch(sd$verifyInvariants(
      as.list(gsub("\\^", "**", b$invariants)),
      lapply(b$preps, function(pr) as.list(pr$comps))),
      error = function(e) rep(FALSE, length(b$invariants)))))
    if (any(!keep)) {
      blocks[[bi]]$certificates <- c(b$certificates, sprintf(
        "verification DROPPED %d invariant(s): X(I) != 0", sum(!keep)))
      blocks[[bi]]$invariants <- b$invariants[keep]
      if (!any(keep)) blocks[[bi]]$status <- "unresolved"
    } else {
      blocks[[bi]]$certificates <- c(b$certificates,
        "verified: X(I) = 0 exactly (sympy cancel) for every generator")
    }
  }
  blocks
}

# iterate over a python iterable, tolerating reticulate's auto-conversion to list
.symRedIter <- function(x, f)
  tryCatch(reticulate::iterate(x, f), error = function(e) lapply(x, f))

# Whether a solved entry stays inside the emitted trafo language: rational
# operations, Rational powers, exp and log -- nothing else (LambertW & friends
# mean the carrier choice was wrong, not that the block is reducible). Returns
# "no", "yes" or "root" (yes, with a non-integer power needing a branch note).
.symRedEntryClass <- function(e, spy) {
  fnames <- tryCatch(unlist(.symRedIter(e$atoms(spy$Function), function(f)
    as.character(reticulate::py_get_attr(reticulate::py_get_attr(f, "func"),
                                         "__name__")))), error = function(err) "?")
  if (length(fnames) && !all(fnames %in% c("exp", "log"))) return("no")
  root <- FALSE
  pows <- tryCatch(.symRedIter(e$atoms(spy$Pow), function(pw) pw),
                   error = function(err) list())
  for (pw in pows) {
    ex <- pw$exp
    if (!isTRUE(ex$is_Rational)) return("no")
    if (!isTRUE(ex$is_Integer)) root <- TRUE
  }
  if (root) "root" else "yes"
}

# Positivity certificate: TRUE when the expression is a ratio of polynomials
# whose coefficients share one sign -- such an entry maps ANY positive outer
# point to a positive inner value, so the chart covers the whole positive
# orthant. Sufficient, not necessary.
# Sign of a RADICAL-FREE polynomial on the positive orthant: +1, -1, or 0 for
# "not decided" (includes the zero expression). Coefficient sign-purity first,
# then AM-GM square absorption (.symRedAbsorb) for mixed patterns.
.symRedSgnPoly <- function(pp, spy) {
  ex <- tryCatch(spy$expand(pp), error = function(err) NULL)
  if (is.null(ex)) return(0L)
  # N(), not Float(): once a quadratic section is in play a coefficient can be an
  # algebraic number such as -sqrt(2), which Float() refuses outright
  if (isTRUE(ex$is_number)) {
    v <- suppressWarnings(as.numeric(as.character(spy$N(ex))))
    if (is.na(v) || v == 0) return(0L)
    return(if (v > 0) 1L else -1L)
  }
  tm <- tryCatch(.symRedIter(spy$Poly(ex)$terms(), function(t) t),
                 error = function(err) NULL)
  if (is.null(tm) || !length(tm)) return(0L)
  nv <- length(unlist(tm[[1]][[1]]))
  expts <- matrix(vapply(tm, function(t) as.numeric(unlist(t[[1]])), numeric(nv)),
                  ncol = nv, byrow = TRUE)
  cf <- vapply(tm, function(t)
    suppressWarnings(as.numeric(as.character(spy$N(t[[2]])))), numeric(1))
  if (anyNA(cf) || any(!is.finite(cf))) return(0L)
  if (all(cf > 0)) return(1L)
  if (all(cf < 0)) return(-1L)
  if (.symRedAbsorb(expts, cf)) return(1L)
  if (.symRedAbsorb(expts, -cf)) return(-1L)
  0L
}

# AM-GM absorption: -c*x^m is dominated by c_u x^u + c_v x^v when u + v = 2m.
# TRUE when every negative term is absorbed AND positive budget remains (full
# consumption certifies >= 0 only). Greedy; sufficient, not complete.
.symRedAbsorb <- function(expts, cf) {
  bud <- ifelse(cf > 0, cf, 0)
  tol <- 1e-9
  for (i in which(cf < 0)) {
    needV <- -cf[i]
    tgt <- 2 * expts[i, ]
    for (u in seq_along(cf)) {
      if (bud[u] <= tol) next
      for (v in seq_along(cf)) {
        if (v <= u || bud[v] <= tol) next
        if (any(expts[u, ] + expts[v, ] != tgt)) next
        take <- min(2 * sqrt(bud[u] * bud[v]), needV)
        bud[u] <- max(0, bud[u] - (take / 2) * sqrt(bud[u] / bud[v]))
        bud[v] <- max(0, bud[v] - (take / 2) * sqrt(bud[v] / bud[u]))
        needV <- needV - take
        if (needV <= tol * abs(cf[i])) break
      }
      if (needV <= tol * abs(cf[i])) break
    }
    if (needV > tol * abs(cf[i])) return(FALSE)
  }
  sum(bud) > tol * sum(abs(cf))
}

# Sign on the positive orthant by structural recursion: numbers, Pow via base,
# Mul via factors, rational functions via num/den, and a + c*sqrt(r) via the
# exact comparison a^2 - c^2 r (the larger magnitude carries the sign). Roots
# split outermost-first, so nesting terminates. +1, -1, or 0 = not decided.
.symRedSgn <- function(pp, spy, depth = 0L) {
  if (depth > 8L) return(0L)
  ex <- tryCatch(spy$expand(pp), error = function(err) NULL)
  if (is.null(ex)) return(0L)
  if (isTRUE(ex$is_number)) {
    v <- suppressWarnings(as.numeric(as.character(spy$N(ex))))
    if (is.na(v) || v == 0) return(0L)
    return(if (v > 0) 1L else -1L)
  }
  if (isTRUE(ex$is_Function) &&
      identical(tryCatch(as.character(reticulate::py_get_attr(
        reticulate::py_get_attr(ex, "func"), "__name__")),
        error = function(err) ""), "exp")) return(1L)
  if (isTRUE(ex$is_Pow) && isTRUE(ex$exp$is_Rational)) {
    sb <- .symRedSgn(ex$base, spy, depth + 1L)
    if (sb == 1L) return(1L)
    if (sb == -1L && isTRUE(ex$exp$is_Integer))
      return(if (as.integer(as.character(ex$exp)) %% 2L == 0L) 1L else -1L)
    return(0L)
  }
  if (isTRUE(ex$is_Mul)) {
    s <- 1L
    for (arg in .symRedIter(spy$Mul$make_args(ex), function(a) a)) {
      sa <- .symRedSgn(arg, spy, depth + 1L)
      if (sa == 0L) return(0L)
      s <- s * sa
    }
    return(as.integer(s))
  }
  # only radicals over a SYMBOLIC base matter here. A numeric one -- the sqrt(2)
  # that comes out of solving 2a^2 = c -- is a positive constant, and Poly()
  # carries it as a coefficient like any other number.
  rads <- tryCatch(Filter(function(w) isTRUE(w$exp$is_Rational) &&
                            !isTRUE(w$exp$is_Integer) && !isTRUE(w$base$is_number),
                          .symRedIter(ex$atoms(spy$Pow), function(w) w)),
                   error = function(err) list())
  if (!length(rads)) {
    fr <- tryCatch(spy$fraction(spy$cancel(spy$together(ex))),
                   error = function(err) NULL)
    if (is.null(fr)) return(0L)
    s1 <- .symRedSgnPoly(fr[[1]], spy)
    s2 <- .symRedSgnPoly(fr[[2]], spy)
    if (s1 == 0L || s2 == 0L) return(0L)
    return(as.integer(s1 * s2))
  }
  len <- vapply(rads, function(w) nchar(as.character(w)), numeric(1))
  for (s in rads[order(-len)]) {
    if (!isTRUE(spy$simplify(s$exp - spy$Rational(1L, 2L))$is_zero)) next
    sp <- tryCatch(ex$as_independent(s, as_Add = TRUE), error = function(err) NULL)
    if (is.null(sp)) next
    a <- sp[[1]]
    cc <- tryCatch(spy$cancel(sp[[2]] / s), error = function(err) NULL)
    if (is.null(cc) || isTRUE(cc$has(s))) next
    if (.symRedSgn(s$base, spy, depth + 1L) != 1L) next   # radicand not positive
    sc <- .symRedSgn(cc, spy, depth + 1L)
    if (isTRUE(a$is_zero)) return(sc)
    sa <- .symRedSgn(a, spy, depth + 1L)
    if (sa == 0L || sc == 0L) next
    if (sa == sc) return(sa)
    d <- tryCatch(spy$expand(a**2L - cc**2L * s$base), error = function(err) NULL)
    if (is.null(d)) next
    sd <- .symRedSgn(d, spy, depth + 1L)
    if (sd == 0L) next
    return(if (sd > 0L) sa else sc)              # the larger magnitude carries the sign
  }
  0L
}

# "now carries" after a shift: I_l - phi_l, composed into inner coordinates.
.symRedShiftedMeaning <- function(invariants, nms, tmpN, sh) {
  m <- setNames(invariants, nms)
  for (l in seq_along(nms)) {
    tl <- tmpN[l]
    if (!is.null(sh$shift[[tl]]))
      m[nms[l]] <- gsub("\\*\\*", "^", as.character(sh$meaning[[tl]]))
  }
  m
}

# Re-carry an invariant shifted: a sign-uncertified factor linear in one
# placeholder t_l has as its root the bound phi (possibly depending on the other
# invariants); t_l -> t_l + phi is a triangular basis change and the outer
# parameter carries I_l - phi, free on (0, Inf). Two certificates per shift:
# phi >= 0, and COVERAGE -- I_l - phi composed into inner coordinates positive
# on the WHOLE orthant, which rejects lossy offsets and picks the placeholder.
# Shifts apply one at a time with re-substitution in between: an offset may
# only become certifiable after an earlier shift has landed.
.symRedShiftFix <- function(es, tmpN, Ies, spy) {
  shift <- list()
  meanE <- setNames(Ies, tmpN)
  tsym <- setNames(lapply(tmpN, spy$Symbol), tmpN)
  # factors of num/den, one extra level (a radical factor keeps its base
  # unfactored, multiplicity 1/2), plus the factors of every radicand base --
  # scanned over the entry's raw, half-power and conjugate-rationalised forms:
  # each form exposes factors the others hide (the raw form the readable ones,
  # the other two those sympy only cancels or factors modulo s^2 = z)
  factorBases <- function(e0) {
    out <- list(); seen <- character(0)
    push <- function(f) {
      k <- as.character(f)
      if (!(k %in% seen)) { seen <<- c(seen, k); out[[length(out) + 1L]] <<- f }
    }
    fparts <- function(part) {
      fl <- tryCatch(spy$factor_list(part), error = function(err) NULL)
      if (is.null(fl)) return(list())
      lapply(.symRedIter(fl[[2]], function(x) x), function(fm)
        if (is.list(fm)) fm[[1]] else reticulate::py_get_item(fm, 0L))
    }
    forms <- list(e0, tryCatch(.symRedHalfPow(e0, spy), error = function(err) NULL),
                  tryCatch(.symRedRadNorm(e0, spy), error = function(err) NULL))
    formSeen <- character(0)
    for (e in forms) {
      if (is.null(e)) next
      ek <- as.character(e)
      if (ek %in% formSeen) next
      formSeen <- c(formSeen, ek)
      fr <- tryCatch(spy$fraction(spy$cancel(spy$together(e))),
                     error = function(err) NULL)
      if (is.null(fr)) next
      lvl1 <- c(fparts(fr[[1]]), fparts(fr[[2]]))
      for (f in lvl1) push(f)
      for (f in lvl1) for (f2 in fparts(f)) push(f2)
      rads <- tryCatch(Filter(function(w) isTRUE(w$exp$is_Rational) &&
                                !isTRUE(w$exp$is_Integer) &&
                                !isTRUE(w$base$is_number),
                              .symRedIter(e$atoms(spy$Pow), function(w) w)),
                       error = function(err) list())
      for (w in rads) for (f2 in fparts(w$base)) push(f2)
    }
    out
  }
  findOne <- function() {
    for (e in es) {
      if (.symRedPosCert(e, spy)) next
      for (f in factorBases(e)) {
        if (.symRedSgn(f, spy) != 0L) next
        fs <- .symRedFreeSyms(f)
        for (tl in intersect(tmpN[!(tmpN %in% names(shift))], fs)) {
          if (length(setdiff(fs, tl))) {
            dg <- tryCatch(as.integer(spy$Poly(f, tsym[[tl]])$degree()),
                           error = function(err) NA_integer_)
            if (is.na(dg) || dg != 1L) next
          }
          rt <- tryCatch(spy$solve(f, tsym[[tl]]), error = function(err) NULL)
          if (is.null(rt) || length(rt) != 1L) next
          s <- .symRedRadNorm(rt[[1]], spy)
          if (.symRedSgn(s, spy) != 1L) next
          sub <- lapply(setdiff(tmpN, tl), function(tj)
            reticulate::tuple(tsym[[tj]], meanE[[tj]]))
          cov <- tryCatch(.symRedRadNorm(meanE[[tl]] - s$subs(sub), spy),
                          error = function(err) NULL)
          if (is.null(cov) || .symRedSgn(cov, spy) != 1L) next
          return(list(tl = tl, s = s, cov = cov))
        }
      }
    }
    NULL
  }
  # es stays in the RAW form: the conjugate rationalisation can turn a positive
  # radical denominator (sqrt(q1) + q2) into an indefinite radical-free one
  # (q2^2 - q1) that sympy cannot cancel back modulo s^2 = q1 -- the multi-form
  # certificate (.symRedPosForm) and the multi-form factor scan above see every
  # form they need without committing the entry to one
  norm <- function(e) tryCatch(spy$cancel(spy$together(e)), error = function(err) e)
  es <- lapply(es, norm)
  for (round in seq_along(tmpN)) {
    hit <- findOne()
    if (is.null(hit)) break
    shift[[hit$tl]] <- hit$s
    meanE[[hit$tl]] <- hit$cov
    pr <- list(reticulate::tuple(tsym[[hit$tl]], tsym[[hit$tl]] + hit$s))
    es <- lapply(es, function(e) norm(e$subs(pr)))
  }
  list(es = es, shift = shift, meaning = meanE)
}

# Rationalise radical denominators by conjugate multiplication (den = a + b*s
# -> times a - b*s), repeated per root. cancel() cannot see factors that exist
# only modulo s^2 = z; this rewrite makes them cancel.
.symRedRadNorm <- function(e, spy, rounds = 3L) {
  for (i in seq_len(rounds)) {
    fr <- tryCatch(spy$fraction(spy$cancel(spy$together(e))),
                   error = function(err) NULL)
    if (is.null(fr)) return(e)
    den <- fr[[2]]
    rads <- tryCatch(Filter(function(w) isTRUE(w$exp$is_Rational) &&
                              !isTRUE(w$exp$is_Integer) &&
                              !isTRUE(w$base$is_number) &&
                              isTRUE(spy$simplify(w$exp -
                                spy$Rational(1L, 2L))$is_zero),
                            .symRedIter(den$atoms(spy$Pow), function(w) w)),
                     error = function(err) list())
    if (!length(rads)) return(spy$cancel(fr[[1]] / den))
    len <- vapply(rads, function(w) nchar(as.character(w)), numeric(1))
    s <- rads[[order(-len)[1]]]
    sp <- tryCatch(den$as_independent(s, as_Add = TRUE), error = function(err) NULL)
    if (is.null(sp)) return(e)
    conj <- sp[[1]] - sp[[2]]
    e2 <- tryCatch(spy$cancel(spy$expand(fr[[1]] * conj) /
                                spy$expand(den * conj)),
                   error = function(err) NULL)
    if (is.null(e2)) return(e)
    e <- e2
  }
  e
}

# Cancel modulo s^2 = z: sympy's cancel() cannot see a common factor that exists
# only in the extension by a square root (x - sqrt(z) against x^2 - z stays
# uncancelled). Substituting every SYMBOL radicand z -> w^2 with w > 0 turns all
# half-integer powers into polynomials in w, where cancel() does see them; w is
# substituted back afterwards, so the result stays in the original symbols.
# Non-symbol radicands are left alone (the a + c*sqrt(r) split in .symRedSgn
# covers them), as are exponents that are not half-integers.
.symRedHalfPow <- function(e, spy) {
  rads <- tryCatch(Filter(function(w) isTRUE(w$exp$is_Rational) &&
                            !isTRUE(w$exp$is_Integer) && isTRUE(w$base$is_Symbol),
                          .symRedIter(e$atoms(spy$Pow), function(w) w)),
                   error = function(err) list())
  bases <- unique(vapply(rads, function(w) as.character(w$base), character(1)))
  bases <- bases[vapply(bases, function(v) all(vapply(
    Filter(function(w) identical(as.character(w$base), v), rads),
    function(w) isTRUE((spy$Integer(2L) * w$exp)$is_Integer), logical(1))),
    logical(1))]
  if (!length(bases)) return(spy$cancel(spy$together(e)))
  ws <- lapply(bases, function(v)
    spy$Symbol(paste0("dModRedH_", v), positive = TRUE))
  fwd <- lapply(seq_along(bases), function(j)
    reticulate::tuple(spy$Symbol(bases[j]), ws[[j]] * ws[[j]]))
  back <- lapply(seq_along(bases), function(j)
    reticulate::tuple(ws[[j]], spy$sqrt(spy$Symbol(bases[j]))))
  out <- spy$cancel(spy$together(e$subs(fwd)))
  spy$cancel(out$subs(back))
}

# Positivity certificate for a solved trafo entry: is it > 0 for EVERY positive value
# of the outer parameters? Numerator and denominator are certified separately and must
# agree in sign. This is what makes an emitted chart global rather than local -- a
# fit in the reduced coordinates can then never leave the model's positive domain.
# The certificate is tried on THREE forms of the entry and the first that
# certifies is returned: the entry as it stands (a conjugate rationalisation can
# only destroy a manifestly positive form like sqrt(q1)*q4/(sqrt(q1) + q2)),
# the half-power cancellation (which sees factors sympy's cancel misses modulo
# s^2 = z), and the conjugate rationalisation (which cancels factors hiding in a
# radical DENOMINATOR). NULL when no form certifies; the certified form is what
# the caller emits -- it is the readable one, and it has no spurious 0/0 points.
.symRedPosForm <- function(e, spy) {
  forms <- list(tryCatch(spy$cancel(spy$together(e)), error = function(err) NULL),
                tryCatch(.symRedHalfPow(e, spy), error = function(err) NULL),
                tryCatch(.symRedRadNorm(e, spy), error = function(err) NULL))
  seen <- character(0)
  for (f in forms) {
    if (is.null(f)) next
    k <- as.character(f)
    if (k %in% seen) next
    seen <- c(seen, k)
    fr <- tryCatch(spy$fraction(f), error = function(err) NULL)
    if (is.null(fr)) next
    s1 <- .symRedSgn(fr[[1]], spy)
    if (s1 != 0L && s1 == .symRedSgn(fr[[2]], spy)) return(f)
  }
  NULL
}

.symRedPosCert <- function(e, spy) !is.null(.symRedPosForm(e, spy))

# Sign on a MIXED domain: the symbols in `realSyms` range over ALL of R, everything
# else over the positive orthant. A carrier whose invariant is not sign-certified
# takes both signs on the positive orthant, so its outer parameter is a real one and
# the chart has to hold there -- certifying it only for positive values of that
# carrier emits a chart that covers part of the model's domain and says nothing
# about the rest. Completes the square in each real symbol, which is exactly the
# form the sum-of-squares pin produces, and hands the remainder to .symRedSgn.
# +1 or 0 = not decided; a negative verdict is never needed here.
.symRedSgnReal <- function(pp, spy, realSyms) {
  if (!length(realSyms)) return(.symRedSgn(pp, spy))
  ex <- tryCatch(spy$expand(pp), error = function(err) NULL)
  if (is.null(ex)) return(0L)
  present <- intersect(realSyms, .symRedFreeSyms(ex))
  if (!length(present)) return(.symRedSgn(ex, spy))
  rest <- ex
  corr <- spy$Integer(0L)
  for (t in present) {
    ts <- spy$Symbol(t)
    A <- tryCatch(spy$expand(rest$coeff(ts, 2L)), error = function(err) NULL)
    B <- tryCatch(spy$expand(rest$coeff(ts, 1L)), error = function(err) NULL)
    C <- tryCatch(spy$expand(rest$coeff(ts, 0L)), error = function(err) NULL)
    if (is.null(A) || is.null(B) || is.null(C)) return(0L)
    # exactly quadratic in t, with coefficients free of every real symbol: anything
    # else is out of the completion's reach and stays undecided
    chk <- tryCatch(as.character(spy$expand(rest - (A * ts^2L + B * ts + C))),
                    error = function(err) "?")
    if (!identical(chk, "0")) return(0L)
    if (length(intersect(realSyms, c(.symRedFreeSyms(A), .symRedFreeSyms(B)))))
      return(0L)
    if (.symRedSgn(A, spy) != 1L) return(0L)
    corr <- corr + B^2L / (4L * A)
    rest <- C
  }
  red <- tryCatch(spy$cancel(spy$together(rest - corr)), error = function(err) NULL)
  if (is.null(red)) 0L else if (.symRedSgn(red, spy) == 1L) 1L else 0L
}

# .symRedPosForm over the mixed domain: same three forms, the mixed sign test.
.symRedPosFormReal <- function(e, spy, realSyms) {
  if (!length(realSyms)) return(.symRedPosForm(e, spy))
  forms <- list(tryCatch(spy$cancel(spy$together(e)), error = function(err) NULL),
                tryCatch(.symRedHalfPow(e, spy), error = function(err) NULL),
                tryCatch(.symRedRadNorm(e, spy), error = function(err) NULL))
  seen <- character(0)
  for (f in forms) {
    if (is.null(f)) next
    k <- as.character(f)
    if (k %in% seen) next
    seen <- c(seen, k)
    fr <- tryCatch(spy$fraction(f), error = function(err) NULL)
    if (is.null(fr)) next
    if (.symRedSgnReal(fr[[1]], spy, realSyms) == 1L &&
        .symRedSgnReal(fr[[2]], spy, realSyms) == 1L) return(f)
  }
  NULL
}

# Gauge-section candidates for a curved block: monomial balances "m1 = m2" over
# the support, simple ones first. A balance section can intersect every positive
# orbit (a constant pin cannot -- the curved orbit may not reach it); which one
# actually yields a positive chart is decided by the certificate above.
#
# `extra` holds the summand monomials of each invariant that is a SUM, and their
# balances go in FRONT of everything: solving such an invariant for one of its
# coordinates produces a difference, positive only where the section puts the
# summands in a fixed ratio, and generic support monomials would push these
# candidates past the scan cap before they are ever tried.
.symRedSectionCands <- function(support, extra = list()) {
  head <- list()
  for (tm in extra) {
    tm <- unique(tm)
    if (length(tm) < 2L) next
    pr <- utils::combn(tm, 2L)
    head <- c(head, lapply(seq_len(ncol(pr)), function(j) c(pr[1, j], pr[2, j])))
  }
  mons <- data.frame(m = support, d = 1L, stringsAsFactors = FALSE)
  if (length(support) > 1L) {
    pr <- utils::combn(support, 2L)
    mons <- rbind(mons, data.frame(m = paste(pr[1, ], pr[2, ], sep = "*"),
                                   d = 2L, stringsAsFactors = FALSE))
  }
  if (nrow(mons) < 2L) return(head)
  idx <- utils::combn(nrow(mons), 2L)
  out <- lapply(seq_len(ncol(idx)), function(j)
    c(mons$m[idx[1, j]], mons$m[idx[2, j]]))
  out <- out[order(mons$d[idx[1, ]] + mons$d[idx[2, ]])]
  key <- function(pr) paste(.symSort(pr), collapse = " = ")
  c(head, out[!vapply(out, key, character(1)) %in% vapply(head, key, character(1))])
}

# The summand monomials of an invariant that is a sum, as printable balance operands
# (sign dropped -- the balance is between magnitudes). Empty for a single-term
# invariant, which needs no ratio pinned.
.symRedInvSummands <- function(Ie, spy) {
  num <- tryCatch(spy$expand(spy$fraction(spy$together(Ie))[[1]]),
                  error = function(e) NULL)
  if (is.null(num)) return(character(0))
  tm <- tryCatch(.symRedIter(spy$Add$make_args(num), function(x) as.character(x)),
                 error = function(e) NULL)
  tm <- unlist(tm)
  if (is.null(tm) || length(tm) < 2L) return(character(0))
  .symSort(unique(sub("^-", "", gsub("\\*\\*", "^", tm))))
}

# Monotone-transversality pre-filter for a balance m1 = m2: the numerator of
# X(log(m1/m2)) must be sign-pure per generator (zero allowed, nonzero once) --
# the section is then crossed at most once. Returns NULL on failure, else its
# values at the two points `pts`, one pair per generator: a SET of sections
# fixes the gauge only when these rows are independent (summand balances are
# all scale-blind, so blind pairs must be pruned before the solve).
.symRedSectionMonotone <- function(pr, preps, xiOf, spy, locals, pts) {
  ms <- lapply(pr, function(m) .symRedSympify(gsub("\\^", "**", m), spy, locals))
  row <- numeric(0); moved <- FALSE
  for (g in seq_along(preps)) {
    Xm <- lapply(ms, function(m)
      .symRedApplyX(preps[[g]], xiOf[[g]], m, spy, locals))
    d <- tryCatch(spy$expand(Xm[[1]] * ms[[2]] - ms[[1]] * Xm[[2]]),
                  error = function(e) NULL)
    if (is.null(d)) return(NULL)
    if (isTRUE(d$is_zero)) { row <- c(row, 0, 0); next }
    if (isTRUE(d$is_number)) {
      v <- suppressWarnings(as.numeric(as.character(spy$N(d))))
      if (is.na(v)) return(NULL)
      moved <- TRUE; row <- c(row, v, v); next
    }
    cf <- tryCatch(unlist(.symRedIter(spy$Poly(d)$coeffs(), function(cc)
      suppressWarnings(as.numeric(as.character(spy$Float(cc)))))),
      error = function(e) NULL)
    if (is.null(cf) || !length(cf) || anyNA(cf)) return(NULL)
    if (!(all(cf > 0) || all(cf < 0))) return(NULL)
    moved <- TRUE
    dstr <- gsub("\\*\\*", "^", as.character(d))
    vals <- vapply(pts, function(pt)
      tryCatch(eval(parse(text = dstr), as.list(pt)),
               error = function(e) NA_real_), numeric(1))
    vals[!is.finite(vals)] <- sign(cf[1])         # overflow: the sign still stands
    row <- c(row, vals)
  }
  if (!moved) return(NULL)
  row
}

# Solve a curved block's invariants into trafo entries. Each invariant is set
# equal to a FRESH outer parameter q_<k> (numbered from invStart, skipping
# coordinate names): a curved carrier is an invariant value, not the coordinate
# it was solved from, and reusing the coordinate's name would invite
# substituting a pinned entry into the others. The system {I_l = c_l} is solved
# jointly in sympy; per invariant a carrier coordinate is chosen that makes the
# solve exact in the trafo language:
#   - "power": v enters the rational part only in degrees {0, k} for one k >= 1
#     (k = 1 is the plain linear case) -> N - c*D = alpha*v^k + beta, a k-th root;
#   - "exp": v sits only inside the single exponential factor, linearly in its
#     argument's numerator -> a log solve.
# Chained solutions (one carrier's entry referencing another) are resolved by
# substitution, the pins of every scaling block are substituted in, and only
# solutions built from rational operations, Rational powers, exp and log are
# emitted -- a block whose invariants cannot be solved in that language stays
# invariantOnly, its invariants still reported.
.symRedSolveInvariants <- function(b, pins, spy, coords = character(0),
                                   invStart = 0L) {
  out <- list(pins = NULL, meaning = NULL, solved = FALSE)
  # a chart certified only for positive values of a real-valued carrier still covers
  # part of the domain; it is kept as the fallback and reported as partial rather
  # than thrown away, so this change can only add coverage, never remove a chart
  best <- NULL
  if (!length(b$invariants)) return(out)
  locals <- .symRedLocals(c(b$invariants, b$support, names(pins), pins), spy)
  Ies <- lapply(b$invariants, function(iv)
    tryCatch(.symRedSympify(gsub("\\^", "**", iv), spy, locals),
             error = function(e) NULL))
  if (any(vapply(Ies, is.null, logical(1)))) return(out)
  # canonical fraction form: sympy's joint solve below can grind for minutes on
  # a sum with embedded quotients where the together'd equivalent solves at once
  Ies <- lapply(Ies, function(e)
    tryCatch(spy$cancel(spy$together(e)), error = function(err) e))
  invN <- character(length(Ies))
  k <- invStart
  for (l in seq_along(Ies)) {
    repeat {
      k <- k + 1L
      if (!(paste0("q_", k) %in% coords)) break
    }
    invN[l] <- paste0("q_", k)
  }
  # the v-degrees present in a polynomial (NULL when Poly refuses)
  degsIn <- function(p, v) {
    pp <- tryCatch(spy$Poly(p, spy$Symbol(v)), error = function(e) NULL)
    if (is.null(pp)) return(NULL)
    dg <- tryCatch(unlist(.symRedIter(pp$monoms(), function(m) as.integer(m[[1]]))),
                   error = function(e) NULL)
    if (is.null(dg)) NULL else sort(unique(dg))
  }
  # One carrier per invariant that no other invariant has claimed. Candidates are
  # tried from the END of the coordinate list first: znames orders states before
  # parameters, and states are better left to the gauge pins, which are constants.
  used <- character(0); carriers <- character(0)
  for (Ie in Ies) {
    cand <- setdiff(b$support, used)
    if (length(coords)) cand <- cand[order(-match(cand, coords, nomatch = 0L))]
    expA <- tryCatch(.symRedIter(Ie$atoms(spy$exp), function(a) a),
                     error = function(e) list())
    if (length(expA) > 1L) return(out)
    G <- if (length(expA)) {
      a <- expA[[1]]$args                # tuple may auto-convert to an R list
      if (is.list(a)) a[[1]] else reticulate::py_get_item(a, 0L)
    } else NULL
    R <- if (length(expA)) spy$cancel(Ie / expA[[1]]) else Ie
    fracR <- spy$fraction(spy$cancel(spy$together(R)))
    gSyms <- if (is.null(G)) character(0) else .symRedFreeSyms(G)
    rSyms <- .symRedFreeSyms(R)
    # pure-power test: the v-degrees of numerator and denominator inside {0, k}
    # for a single k make N - c*D = alpha*v^k + beta, solvable by a k-th root
    # (k = 1 is the plain linear case)
    powerOK <- function(frac, v) {
      dn <- degsIn(frac[[1]], v); dd <- degsIn(frac[[2]], v)
      if (is.null(dn) || is.null(dd)) return(FALSE)
      length(setdiff(unique(c(dn, dd)), 0L)) == 1L
    }
    fracG <- if (!is.null(G)) spy$fraction(spy$cancel(spy$together(G)))
    pick <- NA_character_
    for (v in cand) {
      if (v %in% gSyms) {
        # log solve: v only inside the exponential; G = log(c) is then a plain
        # rational equation, so the same pure-power test applies to G
        if (v %in% rSyms) next
        if (!powerOK(fracG, v)) next
        pick <- v; break
      }
      if (powerOK(fracR, v)) { pick <- v; break }
    }
    if (is.na(pick)) return(out)
    used <- c(used, pick); carriers <- c(carriers, pick)
  }
  # joint solve of I_l = tmp_l for the carriers; tmp_l renamed to the carrier name
  # afterwards (the equation cannot hold the same symbol in both meanings)
  tmpN <- paste0("dModRedC", seq_along(Ies))
  # carriers whose invariant is not certified sign-definite range over R, so every
  # certificate below is taken on that domain and the pin has to clear it
  realTmp <- tmpN[vapply(Ies, function(Ie) .symRedSgn(Ie, spy) != 1L, logical(1))]
  eqs <- lapply(seq_along(Ies), function(l)
    Ies[[l]] - spy$Symbol(tmpN[l]))
  getE <- function(solDict, v)
    if (is.list(solDict)) solDict[[v]] else
      reticulate::py_get_item(solDict, spy$Symbol(v), silent = TRUE)
  # The carrier system does not depend on the gauge choice, so it is eliminated
  # ONCE here: it serves the constant-pin path below and, expressed in the gauge
  # coordinates, every candidate section of the search in between. Solving the
  # full system per candidate instead re-derives this elimination each time --
  # the dominant cost of the whole reduction on a block with several invariants.
  solCarr <- tryCatch(spy$solve(eqs, lapply(carriers, function(v) spy$Symbol(v)),
                                dict = TRUE), error = function(e) NULL)

  # ---- positive gauge sections -----------------------------------------------
  # A constant gauge pin is only a LOCAL chart on a curved orbit (the orbit may
  # not reach the pinned value, and the solved entries can leave the positive
  # orthant). Before falling back to pins, search for a monomial-balance section
  # m1 = m2 whose joint solve with the invariants gives entries certified
  # positive on the whole positive orthant -- a global, log-fittable chart.
  gauge0 <- setdiff(b$support, carriers)
  if (length(gauge0) >= 1L && length(gauge0) <= 2L &&
      !any(grepl("exp(", b$invariants, fixed = TRUE))) {
    scalPinPairs <- lapply(names(pins), function(nm)
      reticulate::tuple(spy$Symbol(nm), .symRedSympify(pins[[nm]], spy, locals)))
    allUnk <- c(carriers, gauge0)
    # balances may involve unmoved coordinates appearing in the invariants (on
    # each orbit they are constants, so the section still cuts transversally);
    # a pair touching no moved coordinate cannot constrain the orbit and is
    # dropped
    invSyms <- unique(unlist(lapply(Ies, .symRedFreeSyms)))
    secVars <- .symSort(unique(c(b$support,
      if (length(coords)) intersect(invSyms, coords) else invSyms)))
    cands <- .symRedSectionCands(secVars,
      lapply(Ies, function(Ie) .symRedInvSummands(Ie, spy)))
    touches <- function(pr) any(vapply(b$support, function(v)
      grepl(paste0("(?<![0-9A-Za-z_.])", v, "(?![0-9A-Za-z_.])"),
            paste(pr, collapse = " "), perl = TRUE), logical(1)))
    cands <- Filter(touches, cands)
    # certified-monotone dials only; collection continues past `need` until the
    # transversality rows span every gauge direction
    rows <- NULL
    # rank PER POINT: over concatenated points two same-generator-blind dials
    # would look independent
    ptRank <- function(M) min(vapply(1:2, function(p)
      qr(M[, seq(p, ncol(M), by = 2L), drop = FALSE], tol = 1e-9)$rank,
      integer(1)))
    if (!is.null(b$preps) && length(cands)) {
      secLoc <- .symRedLocals(c(unlist(lapply(b$preps, `[[`, "comps")),
                                unlist(cands)), spy)
      xiOf <- lapply(b$preps, function(g)
        lapply(g$comps, function(x) .symRedSympify(x, spy, secLoc)))
      evalVars <- .symSort(unique(c(secVars,
        unlist(lapply(b$preps, `[[`, "vars")))))
      pool <- as.numeric(.symPool()(seq_len(2L * length(evalVars)) + 17L))
      pool <- 1 + pool / (max(pool) + 1)          # small values: rows stay finite
      pts <- list(setNames(pool[seq_along(evalVars)], evalVars),
                  setNames(pool[length(evalVars) + seq_along(evalVars)], evalVars))
      need <- if (length(gauge0) == 1L) 40L else 10L
      keep <- list()
      for (pr in cands[seq_len(min(length(cands), 200L))]) {
        row <- .symRedSectionMonotone(pr, b$preps, xiOf, spy, secLoc, pts)
        if (is.null(row)) next
        keep[[length(keep) + 1L]] <- pr
        rows <- rbind(rows, row / max(abs(row)))
        if (length(keep) >= 3L * need) break
        if (length(keep) >= need && ptRank(rows) >= length(gauge0)) break
      }
      cands <- keep
    }
    sets <- if (length(gauge0) == 1L) lapply(cands, list)
      else if (length(cands) > 1L) {
        cmb <- utils::combn(if (is.null(rows)) min(length(cands), 10L)
                            else length(cands), 2L)
        ok <- vapply(seq_len(ncol(cmb)), function(j)
          is.null(rows) ||
            ptRank(rows[cmb[, j], , drop = FALSE]) >= length(gauge0),
          logical(1))
        lapply(which(ok), function(j)
          list(cands[[cmb[1, j]]], cands[[cmb[2, j]]]))
      } else list()
    subsOf <- function(vals) lapply(names(vals), function(v)
      reticulate::tuple(spy$Symbol(v), vals[[v]]))
    # the branches of the carrier elimination, each as (carrier value list, the
    # substitution that puts it into a section equation); empty when the hoisted
    # solve came back empty, in which case each candidate falls back to the joint
    # solve it used before
    carrBr <- if (!is.null(solCarr)) Filter(Negate(is.null), lapply(solCarr,
      function(br) {
        cv <- setNames(lapply(carriers, function(v) getE(br, v)), carriers)
        if (any(vapply(cv, is.null, logical(1)))) NULL else
          list(vals = cv, subs = subsOf(cv))
      }))
    # the gauge coordinates by triangular elimination -- shared with the zero-limit
    # face solve, which sets up the same kind of system; see .symRedTriSolve for why
    # sympy's multivariate solve() is not what runs here
    gaugeSolve <- function(sys, vars) .symRedTriSolve(sys, vars, spy)
    # One candidate's solution branches, as named lists over allUnk, from two
    # elimination orders. Carriers first (the hoisted solCarr) is the cheap one and
    # the one that pays off when several sections are tried against the same
    # invariants. It fails whenever a carrier solves through a root, because the
    # radical then sits INSIDE the section equation and no factor of it has a degree
    # the root finder can use -- the rotation a' = b, b' = -a is the small example.
    # Sections first fixes exactly that: the section is linear in a gauge coordinate,
    # substituting it into the invariants leaves a plain power equation. Both orders
    # are tried PER SECTION (carriers first), so the best-ranked section gets its
    # fallback before the search walks on -- running one order over every section
    # first burned minutes of sympy solves before the winning pair was ever tried.
    branchesOf <- function(secEqs, order) {
      if (identical(order, "sections") || !length(carrBr))
        return(gaugeSolve(c(secEqs, eqs), allUnk))
      out <- list()
      for (cb in carrBr) {
        # the carriers eliminated, the section is an equation in the gauge alone
        sub <- tryCatch(lapply(secEqs, function(e) spy$together(e$subs(cb$subs))),
                        error = function(e) NULL)
        if (is.null(sub)) next
        for (gb in gaugeSolve(sub, gauge0)) {
          if (length(gb) < length(gauge0)) next
          gp <- subsOf(gb)
          out[[length(out) + 1L]] <- c(
            lapply(cb$vals, function(e) e$subs(gp)), gb[gauge0])
        }
      }
      out
    }
    # every equation of the block holds at `vals`, denominators included
    verified <- function(vals, secEqs) {
      pairs <- lapply(allUnk, function(v)
        reticulate::tuple(spy$Symbol(v), getE(vals, v)))
      isZero <- function(e) {
        z <- tryCatch(spy$simplify(e$subs(pairs)), error = function(err) NULL)
        if (is.null(z)) NA else isTRUE(z$is_zero)
      }
      all(vapply(c(eqs, secEqs), function(e) {
        fr <- tryCatch(spy$fraction(spy$together(e)), error = function(err) NULL)
        !is.null(fr) && identical(isZero(fr[[2]]), FALSE) &&
          identical(isZero(fr[[1]]), TRUE)
      }, logical(1)))
    }
    for (st in sets[seq_len(min(length(sets), 40L))])
    for (order in c("carriers", "sections")) {
      secEqs <- lapply(st, function(pr)
        .symRedSympify(pr[1], spy, locals) - .symRedSympify(pr[2], spy, locals))
      sol2 <- branchesOf(secEqs, order)
      if (is.null(sol2) || !length(sol2)) next
      for (bi in seq_along(sol2)) {
        br <- sol2[[bi]]
        es <- list(); okBr <- TRUE
        for (v in allUnk) {
          e0 <- getE(br, v)
          if (is.null(e0) || any(allUnk %in% .symRedFreeSyms(e0))) {
            okBr <- FALSE; break
          }
          es[[v]] <- spy$cancel(e0$subs(scalPinPairs))
        }
        if (!okBr) next
        # an entry certified NEGATIVE for every positive outer value can never be
        # rescued by a carrier shift (a shifted domain is a subset of the positive
        # one), so the sign branches of the joint solve die before the shift search
        if (any(vapply(es, function(e) .symRedSgn(e, spy) == -1L, logical(1))))
          next
        sh <- .symRedShiftFix(es, tmpN, Ies, spy)
        # a shifted carrier is free on (0, Inf) by its own certificate, so it leaves
        # the real-domain set the entries are certified against
        realBr <- setdiff(realTmp, names(sh$shift))
        ent <- character(0); hasRoot <- FALSE; partial <- FALSE
        for (v in allUnk) {
          # a fractional power is fine once it is certified positive: the branch is
          # then pinned by the certificate, not left to the reader. The CERTIFIED
          # form is the one emitted -- it is the readable one and carries no
          # spurious 0/0 point from a conjugate rationalisation.
          e <- .symRedPosFormReal(sh$es[[v]], spy, realBr)
          if (is.null(e)) {
            e <- .symRedPosForm(sh$es[[v]], spy)   # positive values of the carrier only
            if (is.null(e)) { okBr <- FALSE; break }
            partial <- TRUE
          }
          cls <- .symRedEntryClass(e, spy)
          if (identical(cls, "no")) { okBr <- FALSE; break }
          if (identical(cls, "root")) hasRoot <- TRUE
          ent[v] <- gsub("\\*\\*", "^", as.character(e))
        }
        if (!okBr || !verified(br, secEqs)) next
        for (l in seq_along(carriers))
          ent <- setNames(gsub(paste0("\\b", tmpN[l], "\\b"), invN[l], ent),
                          names(ent))
        cand <- out
        cand$pins <- ent
        cand$gauge <- gauge0
        # report the balance in lowest terms: the candidates are pairs of monomials,
        # and a shared factor makes a plain pin read as a relation between two of
        # them (EGFR*k_bind = EGFR^2*k_bind is EGFR = 1)
        cand$section <- vapply(st, function(pr) {
          r <- tryCatch(spy$fraction(spy$cancel(
                 .symRedSympify(pr[1], spy, locals) /
                 .symRedSympify(pr[2], spy, locals))), error = function(e) NULL)
          if (is.null(r)) return(paste(pr[1], "=", pr[2]))
          a <- gsub("\\*\\*", "^", as.character(r[[1]]))
          bq <- gsub("\\*\\*", "^", as.character(r[[2]]))
          if (identical(a, "1")) paste(bq, "= 1") else paste(a, "=", bq)
        }, character(1))
        cand$rootNote <- if (hasRoot)
          "a solved entry carries a square root; the branch is the certified positive one"
        cand$meaning <- .symRedShiftedMeaning(b$invariants, invN, tmpN, sh)
        cand$shifted <- invN[match(names(sh$shift), tmpN)]
        cand$invNames <- invN
        cand$solved <- TRUE
        cand$coverage <- if (partial) "partial" else "total"
        cand$carrierDomain <- setNames(
          ifelse(tmpN %in% realBr, "real", "positive"), invN)
        if (!partial) return(cand)
        if (is.null(best)) best <- cand
        next
      }
    }
  }

  # ---- the constant section, certified the same way ---------------------------
  # Pinning the gauge coordinates to 1 is a section too, and a legitimate one
  # exactly when the chart it produces is positive throughout: every entry
  # sign-pure means the pinned point lies in the positive orthant for EVERY
  # positive outer value, which exhibits an orbit reaching the pin instead of
  # assuming it. That certificate is the whole difference to the old fallback,
  # which emitted the pin unchecked and so could hand back a chart valid only on
  # part of the parameter space. It runs after the balance search because a
  # balance leaves the reachable set open on both sides, where a pin fixes one
  # point of it.
  if (!is.null(solCarr) && length(solCarr)) {
    gauge <- setdiff(b$support, carriers)
    # Pin candidates for the gauge coordinates: the constant 1 first, then -- once a
    # carrier ranges over R -- the sum-of-squares pin 1 + sum_l t_l^2. On an entry
    # affine in the pin that clears every lower bound L at once, since
    # 1 + L^2 - L = (L - 1/2)^2 + 3/4 > 0; on a higher-degree entry it pushes the pin
    # past the Cauchy bound 1 + sum|a_k/a_n| on the real roots, where the entry
    # carries the sign of its leading coefficient. Either way the pinned point exists
    # on EVERY orbit -- which a constant pin cannot promise once the carrier is real.
    pinStr <- list(setNames(rep("1", length(gauge)), gauge))
    if (length(realTmp) && length(gauge))
      pinStr[[2]] <- setNames(
        rep(paste0("1 + ", paste0(realTmp, "**2", collapse = " + ")), length(gauge)),
        gauge)
    pinLoc <- .symRedLocals(c(names(pins), unlist(pins), unlist(pinStr)), spy)
    for (gp in pinStr) {
    allPins  <- c(pins, gp)
    pinPairs <- lapply(names(allPins), function(nm)
      reticulate::tuple(spy$Symbol(nm), .symRedSympify(allPins[[nm]], spy, pinLoc)))
    # every solve branch is a candidate: the first one sympy returns is routinely
    # a negative sign branch, and only one of them pins a positive chart
    for (sol in solCarr) {
    es <- vector("list", length(carriers)); okPin <- TRUE
    for (l in seq_along(carriers)) {
      e <- getE(sol, carriers[l])
      if (is.null(e)) { okPin <- FALSE; break }
      es[[l]] <- e
    }
    # chained solutions: an entry referencing another carrier means that carrier's
    # INNER value -- substitute its solved expression until every entry references
    # tmps, pins and outer symbols only
    if (okPin) for (round in seq_along(carriers)) {
      dirty <- FALSE
      for (l in seq_along(carriers)) {
        hit <- setdiff(which(carriers %in% .symRedFreeSyms(es[[l]])), l)
        if (length(hit)) {
          dirty <- TRUE
          es[[l]] <- es[[l]]$subs(lapply(hit, function(l2)
            reticulate::tuple(spy$Symbol(carriers[l2]), es[[l2]])))
        }
      }
      if (!dirty) break
    }
    sh <- NULL
    if (okPin) {
      pinned <- list()
      for (l in seq_along(carriers)) {
        if (carriers[l] %in% .symRedFreeSyms(es[[l]])) { okPin <- FALSE; break }
        pinned[[carriers[l]]] <- spy$cancel(es[[l]]$subs(pinPairs))
      }
      # a certified-negative entry cannot be rescued by a shift: next branch
      if (okPin && any(vapply(pinned, function(e)
        .symRedSgn(e, spy) == -1L, logical(1)))) okPin <- FALSE
      if (okPin) {
        sh <- .symRedShiftFix(pinned, tmpN, Ies, spy)
        pinned <- sh$es
      }
    }
    entries <- character(0); root <- FALSE; partial <- FALSE
    realPin <- if (is.null(sh)) realTmp else setdiff(realTmp, names(sh$shift))
    if (okPin) for (l in seq_along(carriers)) {
      e <- .symRedPosFormReal(pinned[[carriers[l]]], spy, realPin)
      if (is.null(e)) {
        e <- .symRedPosForm(pinned[[carriers[l]]], spy)
        if (is.null(e)) { okPin <- FALSE; break }
        partial <- TRUE
      }
      cls <- .symRedEntryClass(e, spy)
      if (identical(cls, "no")) { okPin <- FALSE; break }
      if (identical(cls, "root")) root <- TRUE
      entries[carriers[l]] <- gsub("\\*\\*", "^", as.character(e))
    }
    if (okPin) {
      toInv <- function(v) {
        for (l in seq_along(carriers))
          v <- setNames(gsub(paste0("\\b", tmpN[l], "\\b"), invN[l], v), names(v))
        setNames(gsub("\\*\\*", "^", v), names(v))
      }
      cand <- out
      cand$pins    <- c(toInv(entries), toInv(gp))
      cand$gauge   <- gauge
      cand$meaning <- .symRedShiftedMeaning(b$invariants, invN, tmpN, sh)
      cand$shifted <- invN[match(names(sh$shift), tmpN)]
      cand$invNames <- invN
      cand$solved  <- TRUE
      cand$rootNote <- if (root)
        "a solved entry carries a square root; the branch is the certified positive one"
      cand$coverage <- if (partial) "partial" else "total"
      cand$carrierDomain <- setNames(
        ifelse(tmpN %in% realPin, "real", "positive"), invN)
      if (!partial) return(cand)
      if (is.null(best)) best <- cand
    }
    }
    }
  }
  if (!is.null(best)) return(best)

  # Nothing certified. Pinning the gauge coordinates anyway would always yield a
  # chart, but only a LOCAL one: a curved orbit need not pass through the pinned
  # value at all, and the solved entries then leave the positive orthant over part
  # of the outer parameter space -- a fit started there walks off the model's
  # domain with no sign that anything is wrong. A block for which neither section
  # certified is therefore reported with its invariants and the reason.
  out$reason <- paste0(
    "no gauge section with entries certified positive on the whole positive ",
    "orthant (tried: monomial balances",
    if (any(grepl("exp(", b$invariants, fixed = TRUE)))
      " (skipped, a transcendental invariant is in the set)"
    else if (length(gauge0) < 1L || length(gauge0) > 2L)
      sprintf(" (skipped, %d gauge coordinates and the search covers 1 or 2)",
              length(gauge0)),
    ", then the constant pin)")
  out
}


# ---- zero compatibility: which coordinate the orbit can drive to 0 ----------------

# Roots of one equation in one unknown, from the FACTORS of its numerator: only
# degree <= 2 is handed to solve(). A cubic or quartic factor has no closed form the
# callers can use -- its radicals fail the trafo language's entry class -- and asking
# solve() for one is where sympy grinds, minutes on an irreducible symbolic cubic.
.symRedRootsIn <- function(e, v, spy) {
  sv <- spy$Symbol(v)
  num <- tryCatch(spy$fraction(spy$together(e))[[1]], error = function(err) NULL)
  if (is.null(num)) return(list())
  fl <- tryCatch(spy$factor_list(num, sv), error = function(err) NULL)
  if (is.null(fl)) return(list())
  out <- list()
  for (fm in .symRedIter(fl[[2]], function(x) x)) {
    f <- if (is.list(fm)) fm[[1]] else reticulate::py_get_item(fm, 0L)
    dg <- tryCatch(as.integer(spy$Poly(f, sv)$degree()),
                   error = function(err) NA_integer_)
    if (is.na(dg) || dg < 1L || dg > 2L) next
    rs <- tryCatch(spy$solve(f, sv), error = function(err) NULL)
    if (!is.null(rs)) out <- c(out, rs)
  }
  out
}

# The first `m` k-subsets of 1:n, lexicographic. combn() materialises all C(n, k) of
# them -- 9 GiB for 15 of 30 -- where only the first few are ever tried.
.symRedFirstSubsets <- function(n, k, m) {
  if (k > n || k < 0L) return(list())
  if (k == 0L) return(list(integer(0)))
  idx <- seq_len(k); out <- list()
  repeat {
    out[[length(out) + 1L]] <- idx
    if (length(out) >= m) break
    i <- k
    while (i >= 1L && idx[i] == n - k + i) i <- i - 1L
    if (i < 1L) break
    idx[i] <- idx[i] + 1L
    if (i < k) idx[(i + 1L):k] <- idx[i] + seq_len(k - i)
  }
  out
}

# Solve `sys` for `vars` by TRIANGULAR elimination: roots of one equation in one
# unknown, substitute, recurse. sympy's multivariate solve() is what makes the
# searches built on this unaffordable -- two equations with symbolic coefficients can
# grind for many minutes there, and are immediate one unknown at a time. An unknown
# with no usable root falls through to the next elimination order; a branch is
# returned only when every unknown was eliminated. Eliminating through a denominator
# can invent branches that solve a numerator only, so callers substitute an accepted
# branch back into the original equations.
.symRedTriSolve <- function(sys, vars, spy) {
  if (!length(vars)) return(list(setNames(list(), character(0))))
  subsOf <- function(vals) lapply(names(vals), function(v)
    reticulate::tuple(spy$Symbol(v), vals[[v]]))
  for (i in seq_along(sys)) for (v in vars) {
    if (!(v %in% .symRedFreeSyms(sys[[i]]))) next
    roots <- .symRedRootsIn(sys[[i]], v, spy)
    if (!length(roots)) next
    out <- list()
    for (r in roots) {
      rest <- tryCatch(lapply(sys[-i], function(e)
        spy$together(e$subs(list(reticulate::tuple(spy$Symbol(v), r))))),
        error = function(e) NULL)
      if (is.null(rest)) next
      for (s in .symRedTriSolve(rest, setdiff(vars, v), spy)) {
        s[[v]] <- if (length(s)) spy$cancel(r$subs(subsOf(s))) else r
        out[[length(out) + 1L]] <- s
      }
    }
    if (length(out)) return(out)
  }
  list()
}

# "e > 0" as the comparison it is, positive terms left, negated ones right. Valid R
# over the model's own names, so the reader decides it with one eval(). NULL when
# every term shares a sign -- the caller then has the verdict, not a condition.
.symRedIneqStr <- function(e, spy) {
  ex <- tryCatch(spy$expand(e), error = function(err) NULL)
  if (is.null(ex)) return(NULL)
  tm <- tryCatch(unlist(.symRedIter(spy$Add$make_args(ex),
                                    function(x) as.character(x))),
                 error = function(err) NULL)
  if (!length(tm)) return(NULL)
  tm <- gsub("\\*\\*", "^", tm)
  neg <- grepl("^-", tm)
  if (all(neg) || !any(neg)) return(NULL)
  paste(paste(tm[!neg], collapse = " + "), ">",
        paste(sub("^-", "", tm[neg]), collapse = " + "))
}

# What "e > 0" reduces to on the positive orthant. sign(num/den) = sign(num*den), so
# one polynomial decides; its sign-definite factors drop out (a positive one changes
# nothing, a negative one flips the comparison). Returns sign 1/-1 when certified
# either way, 0 when identically zero, NA with `cond` the surviving inequality.
.symRedCondition <- function(e, spy) {
  na <- list(sign = NA_integer_, cond = NULL)
  fr <- tryCatch(spy$fraction(spy$cancel(spy$together(e))),
                 error = function(err) NULL)
  if (is.null(fr)) return(na)
  if (isTRUE(tryCatch(spy$expand(fr[[1]])$is_zero, error = function(err) FALSE)))
    return(list(sign = 0L, cond = NULL))
  p <- tryCatch(spy$expand(fr[[1]] * fr[[2]]), error = function(err) NULL)
  if (is.null(p)) return(na)
  fl <- tryCatch(spy$factor_list(p), error = function(err) NULL)
  if (is.null(fl)) return(na)
  s <- .symRedSgnPoly(fl[[1]], spy)
  if (s == 0L) return(na)
  keep <- list()
  for (fm in .symRedIter(fl[[2]], function(x) x)) {
    f <- if (is.list(fm)) fm[[1]] else reticulate::py_get_item(fm, 0L)
    k <- if (is.list(fm)) fm[[2]] else reticulate::py_get_item(fm, 1L)
    k <- suppressWarnings(as.integer(as.character(k)))
    if (is.na(k)) return(na)
    if (k %% 2L == 0L) next                   # an even power carries no sign
    sf <- .symRedSgn(f, spy)
    if (sf == 1L) next
    if (sf == -1L) { s <- -s; next }
    keep[[length(keep) + 1L]] <- f
  }
  if (!length(keep)) return(list(sign = as.integer(s), cond = NULL))
  q <- Reduce(function(a, bb) a * bb, keep)
  if (s < 0L) q <- q * spy$Integer(-1L)
  sq <- .symRedSgnPoly(q, spy)
  if (sq != 0L) return(list(sign = sq, cond = NULL))
  list(sign = NA_integer_, cond = .symRedIneqStr(q, spy))
}

# Which coordinates can the orbit drive to 0 with NOTHING running off to infinity?
# Several at once is fine -- a set of rates switched off is still a model, and the
# flat direction ending there means it fits exactly as well. Only divergence
# disqualifies: a zero bought by sending another coordinate to infinity is no zero of
# the model.
#
# The orbit lies in the level set of the invariants, so the face {z_Z = 0} is asked
# for the point's invariant values: with the face coordinates primed and z_Z' = 0,
# I_l(z') = I_l(z) clears to num'(z')*den(z) - num(z)*den'(z') = 0, triangular-solved
# for the rest. A solved entry certified positive is silent, an undecided one IS the
# condition, and one that cannot be positive names a coordinate that has to vanish
# along -- which grows Z and asks again, so the sets come out of the solve rather
# than from enumerating subsets. An equation that survives with no unknown left pins
# an invariant at a value no face point carries: that is divergence, and it stays
# blocked however Z grows, since the equation no longer depends on Z. An equation
# that cancels identically is the 0/0 of a rational invariant on the joint face --
# vacuous, and dropped.
#
# Two gaps, both reported rather than hidden: the level set may have components
# beyond the orbit, so a positive verdict exhibits the invariant values and not a
# path to them; and an incomplete invariant set widens it further, which leaves only
# "never" exact (column `certain`).
.symRedZeroCompat <- function(b, spy, fixed = character(0), verbose = FALSE) {
  supp <- setdiff(b$support, fixed)
  invs <- b$invariants
  if (identical(b$status, "fixed") || !length(supp) || !length(invs)) return(NULL)
  locals <- .symRedLocals(c(invs, supp), spy)
  Ies <- lapply(invs, function(iv)
    tryCatch(.symRedSympify(gsub("\\^", "**", iv), spy, locals),
             error = function(e) NULL))
  if (any(vapply(Ies, is.null, logical(1)))) return(NULL)
  fr <- lapply(Ies, function(e)
    tryCatch(spy$fraction(spy$cancel(spy$together(e))), error = function(err) NULL))
  if (any(vapply(fr, is.null, logical(1)))) return(NULL)

  sym <- function(v) spy$Symbol(v)
  zero <- spy$Integer(0L)
  pv <- setNames(paste0("dModZero", seq_along(supp)), supp)
  certain <- is.null(b$target) || length(invs) >= b$target
  compLoc <- if (!is.null(b$preps))
    .symRedLocals(unlist(lapply(b$preps, `[[`, "comps")), spy)

  # v | xi_v for every generator: {v = 0} is invariant, so the orbit only approaches
  # it, as eps -> +-Inf. A scaling has xi_v = w_v*v and is always of that kind.
  darboux <- function(v) {
    if (is.null(b$preps)) return(TRUE)
    all(vapply(b$preps, function(pr) {
      if (!(v %in% names(pr$comps))) return(TRUE)
      e <- tryCatch(.symRedSympify(pr$comps[[v]], spy, compLoc),
                    error = function(err) NULL)
      !is.null(e) &&
        isTRUE(spy$expand(e$subs(list(reticulate::tuple(sym(v), zero))))$is_zero)
    }, logical(1)))
  }
  # several directions leave more face coordinates than equations: the extra ones
  # stay free and positive, and each choice of which to solve for is tried
  faceSolve <- function(eqs, unk) {
    if (!length(unk)) return(list(setNames(list(), character(0))))
    if (length(eqs) >= length(unk)) return(.symRedTriSolve(eqs, unk, spy))
    out <- list()
    for (cols in .symRedFirstSubsets(length(unk), length(eqs), 12L))
      out <- c(out, .symRedTriSolve(eqs, unk[cols], spy))
    out
  }

  # the face system for one zero set: "blocked" (divergence), NULL (not parseable),
  # or the equations with the coordinates still free
  faceEqs <- function(Z) {
    other <- setdiff(supp, Z)
    subsV <- c(lapply(other, function(u) reticulate::tuple(sym(u), sym(pv[[u]]))),
               lapply(Z, function(u) reticulate::tuple(sym(u), zero)))
    unk <- unname(pv[other])
    eqs <- list(); dens <- list()
    for (l in seq_along(fr)) {
      E <- tryCatch(spy$expand(fr[[l]][[1]]$subs(subsV) * fr[[l]][[2]] -
                               fr[[l]][[1]] * fr[[l]][[2]]$subs(subsV)),
                    error = function(err) NULL)
      if (is.null(E)) return(NULL)
      if (isTRUE(E$is_zero)) next
      if (!any(unk %in% .symRedFreeSyms(E))) return("blocked")
      eqs[[length(eqs) + 1L]] <- E
      dens[[length(dens) + 1L]] <- fr[[l]][[2]]$subs(subsV)
    }
    list(eqs = eqs, dens = dens, unk = unk, other = other)
  }

  # one zero set, answered: "yes"/"if" with the conditions and the point it lands on,
  # "grow" with the coordinates that have to vanish along, "no", or "unknown"
  analyse <- function(Z) {
    f <- faceEqs(Z)
    if (is.null(f)) return(list(verdict = "unknown"))
    if (identical(f, "blocked")) return(list(verdict = "no"))
    sols <- tryCatch(faceSolve(f$eqs, f$unk), error = function(e) list())
    grow <- character(0); unclear <- FALSE
    for (br in sols) {
      free <- vapply(f$unk, function(u) is.null(br[[u]]), logical(1))
      vals <- setNames(lapply(f$unk, function(u)
        if (is.null(br[[u]])) sym(u) else br[[u]]), f$unk)
      pairs <- lapply(f$unk, function(u) reticulate::tuple(sym(u), vals[[u]]))
      isZero <- function(e) {
        z <- tryCatch(spy$simplify(e$subs(pairs)), error = function(err) NULL)
        !is.null(z) && isTRUE(z$is_zero)
      }
      # the equations it was not eliminated from
      if (!all(vapply(f$eqs, isZero, logical(1)))) next
      cnd <- character(0); dead <- character(0); open <- FALSE
      for (i in seq_along(f$unk)) {
        s <- .symRedCondition(vals[[f$unk[i]]], spy)
        if (identical(s$sign, 1L)) next
        # cannot be positive here: it is 0 on this face, or it would have to be
        # negative -- either way its own zero belongs in the set
        if (!is.na(s$sign)) { dead <- c(dead, f$other[i]); next }
        # a condition on a free coordinate constrains no parameter: it says the
        # face is met for SOME positive value of it
        if (is.null(s$cond) || any(vapply(f$unk[free], function(w)
              grepl(paste0("\\b", w, "\\b"), s$cond), logical(1)))) {
          open <- TRUE; next
        }
        cnd <- unique(c(cnd, s$cond))
      }
      # growth is read off the entries BEFORE the denominators are consulted: a
      # rational invariant is routinely undefined on the face a coordinate vanishes
      # on (b/a at a = 0), and that face is exactly the one to grow away from
      if (length(dead)) { grow <- unique(c(grow, dead)); next }
      if (any(vapply(f$dens, isZero, logical(1)))) next   # invariant undefined there
      if (open) { unclear <- TRUE; next }
      # a free coordinate enters the others by name: `s'` is the value the
      # remaining freedom picks for s
      at <- paste(paste0(c(Z, f$other), " = ",
                         c(rep("0", length(Z)),
                           ifelse(free, "(free)", vapply(vals[f$unk], function(e)
                             gsub("\\*\\*", "^", as.character(e)), character(1))))),
                  collapse = ", ")
      for (u in names(pv))
        at <- gsub(paste0("\\b", pv[[u]], "\\b"), paste0(u, "'"), at)
      return(list(verdict = if (length(cnd)) "if" else "yes",
                  condition = paste(cnd, collapse = " & "), at = at))
    }
    if (length(grow)) return(list(verdict = "grow", grow = grow))
    list(verdict = if (unclear || !length(sols)) "unknown" else "no")
  }

  rows <- NULL; seen <- character(0)
  for (v0 in supp) {
    if (isTRUE(verbose))
      message("zero compatibility {", paste(b$labels, collapse = ", "), "}: ", v0)
    Z <- v0
    repeat {
      a <- analyse(Z)
      if (!identical(a$verdict, "grow")) break
      Z2 <- .symSort(unique(c(Z, a$grow)))
      if (length(Z2) == length(Z)) { a$verdict <- "no"; break }
      Z <- Z2
    }
    key <- paste(Z, collapse = ", ")
    if (key %in% seen) next
    seen <- c(seen, key)
    rows <- rbind(rows, data.frame(
      coordinates = key, verdict = a$verdict,
      # the joint face is met at finite eps only if none of its coordinates is a
      # Darboux one, which no orbit crosses
      limit = if (a$verdict %in% c("yes", "if")) any(vapply(Z, darboux, logical(1)))
              else NA,
      certain = certain || identical(a$verdict, "no"),
      condition = if (is.null(a$condition)) "" else a$condition,
      at = if (is.null(a$at)) "" else a$at, stringsAsFactors = FALSE))
  }
  rows
}

# One line per zero set the orbit reaches, and under which condition. A conditional
# zero is never announced as reachable: P/pP reaches {k_d = 0} from one side of its
# steady state and {P = 0} from the other, so neither is a property of the model and
# the condition is all there is to report. Coordinates no set reaches are not
# reported at all -- a zero that cannot happen is not news, and a model of thirty
# parameters would drown the two that can in the twenty-eight that cannot. They stay
# on `$zeroCompatibility` with verdict "no", as does the point each set lands on
# (`at`), data to build a reduced model from rather than a line to read past.
.symRedCatZeroCompat <- function(x, width) {
  v <- x$zeroCompatibility
  if (is.null(v) || !nrow(v)) return(invisible(NULL))
  # an incomplete invariant set widens the level set: "never" survives that,
  # anything else is an upper bound and is not stated as more
  open <- !v$certain & v$verdict != "no"
  hit <- which(!open & v$verdict %in% c("yes", "if"))
  parts <- function(i) unique(unlist(strsplit(v$coordinates[i], ", ", fixed = TRUE)))
  undecided <- setdiff(parts(which(open | v$verdict == "unknown")), parts(hit))
  if (!length(hit) && !length(undecided)) return(invisible(NULL))
  cat("\nZero limits\n")
  w <- max(nchar(v$coordinates[hit], type = "width"), 0L)
  for (i in hit) {
    txt <- paste0(if (v$verdict[i] == "yes") "everywhere"
                  else paste0("where  ", v$condition[i]),
                  if (isTRUE(v$limit[i])) ", in the limit" else "")
    lead <- paste0("  ", formatC(v$coordinates[i], width = -w), " = 0   ")
    cat(lead, .symRedWrap(txt, strrep(" ", nchar(lead, type = "width")), width),
        "\n", sep = "")
  }
  if (length(undecided))
    cat("  undecided: ", .symRedWrap(.symSort(undecided), "    ", width), "\n",
        sep = "")
  invisible(NULL)
}

# The complete P()-ready trafo: identity over every coordinate, the scaling pins,
# and the solved curved entries. Scaling pins are constants and curved entries come
# from a joint solve referencing only outer parameters, so no resolveRecurrence
# rewriting is needed; asserted: no entry references the INNER meaning of a
# different redefined coordinate through a chain a pin would break.
.symRedAssembleTrafo <- function(blocks, coords) {
  vals <- setNames(as.character(coords), coords)
  for (b in blocks) if (!is.null(b$pins)) {
    nm <- intersect(names(b$pins), names(vals))
    vals[nm] <- b$pins[nm]
  }
  as.eqnvec(structure(as.character(vals), names = names(vals)))
}

# ---- result object and print -----------------------------------------------------

.symRedResult <- function(object, blocks, trafo, coords, fixed, settings, call) {
  status <- vapply(blocks, `[[`, character(1), "status")
  labs <- lapply(blocks, `[[`, "labels")
  van <- do.call(rbind, lapply(blocks, function(b)
    if (is.null(b$zeroCompatibility)) NULL else
      cbind(block = paste(b$labels, collapse = ", "), b$zeroCompatibility,
            stringsAsFactors = FALSE)))
  structure(list(
    method      = object$method,
    coordinates = coords,
    fixed       = fixed,
    blocks      = blocks,
    trafo       = trafo,
    family      = Filter(Negate(is.null), lapply(blocks, function(b)
      if (identical(b$status, "reduced") &&
          (!is.null(b$admissible) || !is.null(b$matroid) ||
           !is.null(b$gaugeNote)))
        list(labels = b$labels, admissible = b$admissible, matroid = b$matroid,
             gaugeNote = b$gaugeNote))),
    removed     = unlist(labs[status %in% c("fixed", "reduced")]),
    remaining   = unlist(labs[!status %in% c("fixed", "reduced")]),
    zeroCompatibility = van,
    settings    = settings,
    call        = call), class = "symmetryreduction")
}

# wrap-and-indent for the print method (display columns, not bytes)
.symRedWrap <- function(items, indent, width, sep = ", ") {
  txt <- paste(items, collapse = sep)
  paste(strwrap(txt, width = max(30L, width - nchar(indent)),
                initial = "", prefix = ""), collapse = paste0("\n", indent))
}

# name = value, aligned on the widest name, continuation lines under the value
.symRedCatPairs <- function(nms, vals, width, ind = "  ") {
  w <- max(nchar(nms, type = "width"))
  for (i in seq_along(nms)) {
    lead <- paste0(ind, formatC(nms[i], width = -w), " = ")
    cat(lead, .symRedWrap(vals[i], strrep(" ", nchar(lead, type = "width")), width),
        "\n", sep = "")
  }
}

#' @export
# print() is deliberately terse, as print.symmetrydetection() is: the verdict, the
# reparametrisation and what the outer parameters now mean -- nothing a reader has
# to skip to reach the trafo. summary() adds the block report on top.
print.symmetryreduction <- function(x, width = getOption("width"), ...) {
  if (!length(x$removed) + length(x$remaining)) {
    cat("Nothing to reduce.\n"); return(invisible(x))
  }
  cat(.symRedVerdict(x), "\n", sep = "")
  .symRedCatChart(x, width)
  .symRedCatZeroCompat(x, width)
  for (b in x$blocks) if (!b$status %in% c("reduced", "fixed"))
    cat("\nNot reduced: ", paste(b$labels, collapse = ", "),
        if (length(b$invariants))
          paste0(" | invariants: ",
                 .symRedWrap(b$invariants, "    ", width)) else "",
        "\n", sep = "")
  invisible(x)
}

#' @export
summary.symmetryreduction <- function(object, verbose = FALSE,
                                     width = getOption("width"), ...)
  .symRedReport(object, verbose, width)

.symRedVerdict <- function(x) {
  nDir <- length(x$removed) + length(x$remaining)
  sprintf("Reduced %d of %d direction%s%s.",
          length(x$removed), nDir, if (nDir == 1L) "" else "s",
          if (length(x$remaining))
            sprintf("  (%s remaining)", paste(x$remaining, collapse = ", "))
          else "")
}

# the reparametrisation itself: the non-identity entries, and the invariant each
# outer name carries -- `q_1 = k_p + k_d` reads as "the invariant k_p + k_d is the
# outer parameter q_1". Shared by print() and summary(), so neither repeats it.
.symRedCatChart <- function(x, width) {
  nonid <- x$trafo[x$trafo != names(x$trafo)]
  if (length(nonid)) {
    cat("\nTrafo (non-identity entries)\n")
    .symRedCatPairs(names(nonid), as.character(nonid), width)
  }
  meaning <- unlist(lapply(x$blocks, `[[`, "survivorMeaning"))
  if (length(meaning)) {
    # which outer names are log-fittable is the one thing a reader cannot see from
    # the expression: a carrier whose invariant takes both signs on the positive
    # orthant is a REAL parameter, and putting it on a log scale would silently
    # confine the fit to half the model's domain
    dom  <- unlist(lapply(x$blocks, `[[`, "carrierDomain"))
    real <- names(meaning) %in% names(dom)[dom == "real"]
    cat("\nInvariants:\n")
    .symRedCatPairs(names(meaning),
                    paste0(unname(meaning), ifelse(real, "        [real-valued]", "")),
                    width)
    if (any(real))
      cat("  [real-valued] takes both signs on the positive orthant -- fit it",
          "linearly, not on a log scale\n")
  }
}

# one line per block: what it is, which stage answered, how it was gauged. A block
# that did NOT reduce adds its invariants and the reason -- how many invariants
# were found under which degree caps, the only thing there is to act on. The
# per-stage certificates stay on the object (`$blocks[[i]]$certificates`); printed,
# they bury the result under its own provenance.
.symRedBlockLines <- function(b, fam, verbose, width) {
  kind <- if (is.null(b$kind)) b$type else b$kind
  if (identical(kind, "curved")) kind <- "general"   # the detection report's word
  gauge <- character(0)
  if (!is.null(b$removedByFixed) && b$removedByFixed > 0L)
    gauge <- c(gauge, sprintf("fixed removed %d direction%s%s", b$removedByFixed,
                              if (b$removedByFixed == 1L) "" else "s",
                              if (length(b$redundantFixed))
                                paste0(" (redundant: ",
                                       paste(b$redundantFixed, collapse = ", "), ")")
                              else ""))
  if (length(b$section))
    gauge <- c(gauge, paste("section", paste(b$section, collapse = ", ")))
  else if (length(b$transversal))
    gauge <- c(gauge, paste("transversal",
                            paste(paste0(b$transversal, " = ",
                                         if (!is.null(b$pins)) b$pins[b$transversal]
                                         else "1"), collapse = ", ")))
  if (!is.null(fam$admissible))
    gauge <- c(gauge, sprintf("%d admissible", length(fam$admissible)))
  else if (!is.null(fam$matroid))
    gauge <- c(gauge, sprintf("matroid of %d row%s", length(fam$matroid),
                              if (length(fam$matroid) == 1L) "" else "s"))
  full <- !b$status %in% c("reduced", "fixed") || isTRUE(verbose)
  cat("  {", paste(b$labels, collapse = ", "), "} ", kind, ", ", b$status,
      if (!is.null(b$stage) && !b$stage %in% c("transversal", "fixed", "none"))
        paste0(" [", b$stage, "]") else "",
      if (length(gauge)) paste0(" | ", paste(gauge, collapse = ", ")) else "",
      "\n", sep = "")
  ind <- "      "
  if (length(b$moduleCombos))
    cat(ind, "module reduction  ",
        .symRedWrap(b$moduleCombos, paste0(ind, "                  "), width,
                    sep = ";  "), "\n", sep = "")
  # the invariants of a REDUCED block are the "outer parameters" list above, named
  # by their carrier -- printing them again says nothing new. A block that did not
  # reduce has no carriers, so there they are the result.
  if (length(b$invariants) && full)
    cat(ind, "invariants  ",
        .symRedWrap(b$invariants, paste0(ind, "            "), width), "\n",
        sep = "")
  if (!is.null(b$reason))
    cat(ind, "reason  ",
        .symRedWrap(b$reason, paste0(ind, "        "), width, sep = " "), "\n",
        sep = "")
  if (isTRUE(verbose)) {
    if (!is.null(fam$admissible))
      cat(ind, "admissible  ", .symRedWrap(vapply(fam$admissible, function(T)
        paste0("{", paste(T, collapse = ","), "}"), character(1)),
        paste0(ind, "            "), width), "\n", sep = "")
    if (!is.null(fam$matroid))
      for (r in fam$matroid)
        cat(ind, "pick one of  ", .symRedWrap(r, paste0(ind, "             "),
                                              width), "\n", sep = "")
    if (!is.null(fam$gaugeNote))
      cat(ind, "gauge  ", .symRedWrap(fam$gaugeNote, paste0(ind, "       "),
                                      width, sep = " "), "\n", sep = "")
  }
}

# the report: the chart print() gives, plus one block per coupled set of
# directions -- kind, status, stage, how it was gauged, and for a block that did
# not reduce its invariants and the reason. verbose = TRUE adds the admissible
# gauges, the invariants of the reduced blocks and the search caps.
.symRedReport <- function(x, verbose = FALSE, width = getOption("width")) {
  bar <- strrep("-", 60)
  nDir <- length(x$removed) + length(x$remaining)
  cat(bar, "\n", sep = "")
  cat(sprintf("symmetryReduction  |  from: %s   directions: %d\n", x$method, nDir))
  cat(bar, "\n", sep = "")
  if (!nDir) { cat("Nothing to reduce.\n"); return(invisible(x)) }
  if (isTRUE(verbose) && length(x$settings))
    cat("caps: ", paste(paste0(names(x$settings), "=",
                               vapply(x$settings, function(v)
                                 paste(v, collapse = ","), character(1))),
                        collapse = ", "), "\n", sep = "")
  cat(.symRedVerdict(x), "\n", sep = "")
  .symRedCatChart(x, width)
  .symRedCatZeroCompat(x, width)
  cat("\nBlocks\n")
  famOf <- function(b) {
    hit <- Filter(function(f) identical(f$labels, b$labels), x$family)
    if (length(hit)) hit[[1]] else list()
  }
  for (b in x$blocks) .symRedBlockLines(b, famOf(b), verbose, width)
  if (length(x$remaining))
    cat("\n", .symRedWrap(paste0("Remaining: ",
      paste(x$remaining, collapse = ", "),
      " | options: a structural assumption, a separating experiment, or ",
      "prediction profiles."), "", width, sep = ""), "\n", sep = "")
  if (length(x$fixed))
    cat("\n", .symRedWrap(paste0("fixed coordinates stay identity entries; keep ",
      "passing fixed = to downstream calls."), "", width, sep = ""), "\n", sep = "")
  invisible(x)
}

# ---- the user-facing function ----------------------------------------------------

#' Constructive removal of detected symmetries
#'
#' Turns the non-identifiable directions reported by [symmetryDetection()] into a
#' parameter reparametrisation. Scaling directions are gauged exactly: the integer
#' weight lattice yields the invariant monomials, and a certified transversal (one
#' coordinate per independent weight row, pinned to 1) removes the directions while
#' every surviving coordinate keeps its name and absorbs an invariant product. When
#' a whole space of transformations exists, the admissible family is reported
#' instead of hiding an arbitrary choice. Curved (polynomial/general) directions go
#' through module reduction and an escalating exact invariant search -- monomial,
#' polynomial up to `dPoly`, rational with a single-coordinate denominator (a
#' Laurent ansatz with numerator degree up to `dPoly`), rational via Darboux
#' polynomials up to `dDarboux`, Liouvillian via exponential factors `exp(g/h)`
#' with numerator degree up to `dExp` -- and every failed stage leaves a precise
#' negative certificate. By
#' Prelle--Singer the last stage completes the search class: every elementary or
#' Liouvillian first integral is of Darboux form once exponential factors are
#' admitted, so within the degree caps nothing expressible in closed form is
#' missed. A solved curved block carries each invariant on a fresh outer parameter
#' `q_<k>`; the `survivorMeaning` entries state which invariant each one holds.
#'
#' Each carrier is emitted with a DOMAIN. An invariant certified positive on the
#' positive orthant gives a positive, log-fittable carrier; one that takes both signs
#' there gives a real one, reported as `[real-valued]` and to be fitted linearly --
#' putting it on a log scale would confine the fit to part of the model's domain. The
#' gauge is chosen to match: a constant pin only reaches the orbits that happen to
#' cross it, so when a carrier is real the gauge coordinate is pinned to
#' `1 + sum_l q_l^2` instead. On an entry affine in the pin that clears every lower
#' bound at once (`1 + L^2 - L = (L - 1/2)^2 + 3/4 > 0`), and on a higher-degree entry
#' it pushes the pin past the Cauchy bound on the real roots, where the entry carries
#' the sign of its leading coefficient. The chart then reaches every orbit that meets
#' the positive orthant (`coverage = "total"`). A block where no such gauge certified
#' keeps the positive-domain chart and is reported as `coverage = "partial"`.
#'
#' With `reportZeroCompatibility = TRUE` every block also reports its ZERO LIMITS:
#' which coordinates the orbit can drive to 0 with nothing running off to infinity.
#' That is the difference between a flat direction that runs forever and one that
#' ends in a degenerate model -- a rate switched off, its whole reaction gone --
#' which fits the data exactly as well and is rarely what a symmetry is wanted for.
#' Several coordinates may vanish together (a set of rates switched off is still a
#' model); only divergence disqualifies, a zero bought by sending another coordinate
#' to infinity being no zero of the model. The orbit lies inside the level set of the
#' block's invariants, so the face `{z_Z = 0}` is asked for the invariant values a
#' point carries: the face system is solved exactly, a solved coordinate that cannot
#' be positive there is added to `Z` and the face asked again -- so the sets come out
#' of the solve, not from enumerating subsets -- and an equation left with no unknown
#' is the divergence. Whether the face is reached at finite `eps` or only approached
#' as `eps -> +-Inf` follows from `v | xi_v` (`limit`).
#'
#' A zero is never ANNOUNCED as reachable unless it is reachable from every positive
#' point; where it is not, the report states the inequality under which it is, and
#' leaves the call to the reader. Whether `k_d` can be removed from `P <-> pP`
#' depends on which side of its steady state the parameters sit (`P*k_p > k_d*pP`),
#' and on the other side `P = 0` is the reachable one instead: that is a statement
#' about the parameters, not about the model. The condition is an R expression over
#' the model's own coordinate names, so one `eval()` settles it at whatever point
#' matters -- the reduction itself stays symbolic and takes no parameter values.
#'
#' @param object A `symmetrydetection` result, from [symmetryDetection()].
#' @param fixed Character vector of coordinates the user pins at known values
#'   beforehand (same semantics as `summary(object, fixed = )`): scaling directions
#'   removed by the fixing drop out of the reparametrisation, and the fixed
#'   coordinates never enter a transversal. Unknown names warn and are ignored.
#' @param reportZeroCompatibility Logical, off by default: work out which zeros the
#'   symmetry is compatible with (see below), one exact face solve per coordinate of
#'   each block. It costs a quarter of the reduction on a six-reaction cascade and is
#'   asked for, not paid for by everyone.
#' @param dPoly Total-degree cap of the polynomial invariant search (stage 2),
#'   and the numerator-degree cap of the rational stage that follows it (a
#'   Laurent ansatz allowing one moved coordinate at exponent -1 -- the cheap
#'   route to invariants like `(z^2*a + z*b + c*d)/z` that the factor stages
#'   only reach at much higher caps).
#' @param dDarboux Degree cap of the Darboux factors in the rational invariant
#'   search (stage 3). The cofactor degree is bounded structurally (at most the
#'   generator degree minus one), not by this cap.
#' @param dExp Numerator degree cap of the exponential-factor search (stage 4):
#'   factors `exp(g/h)` with `deg(g) <= dExp` and `h` a product of Darboux
#'   factors. `dExp = 0` skips the stage. Invariants from this stage are
#'   transcendental; the emitted trafo entries may contain `log`/`exp` (fine for
#'   [P()], but not accepted by `symmetryDetection(trafo = )`, which requires
#'   rational entries).
#' @param separable Logical: run the separable-characteristics stage. When every
#'   component of a generator involves no moved coordinate but its own, the
#'   characteristic system decouples and the invariants follow from one-dimensional
#'   quadratures -- cheaper than every other stage and the only one reaching
#'   antiderivatives outside the Darboux language (`atan`). `FALSE` forces the
#'   ansatz and factor stages to carry such a block on their own.
#' @param verbose Logical: print progress per block and stage.
#' @param ... Reserved.
#'
#' @return An object of class `symmetryreduction`:
#'   \describe{
#'     \item{`blocks`}{one entry per coupled block of directions: `labels` (the
#'       `X` labels as printed by `print(object)`), `type` (`kind` names what the
#'       block holds when a scaling was demoted into a curved block), `support`, `stage`
#'       reached, `invariants` (exact strings), `certificates` (positive and
#'       negative), `transversal`/`pins`, `survivorMeaning` (which invariant each
#'       surviving coordinate or fresh `q_<k>` parameter carries), `carrierDomain`
#'       (`"positive"` or `"real"` per carrier) and `coverage` (`"total"` when the
#'       chart reaches every orbit that meets the positive orthant, `"partial"`
#'       otherwise), `moduleCombos`, `status`
#'       (`"fixed"`, `"reduced"`, `"invariantOnly"` or `"unresolved"`), `reason`,
#'       `zeroCompatibility`.}
#'     \item{`zeroCompatibility`}{the zero limits of every block, one row per set of
#'       coordinates that vanish together (`coordinates`, comma-separated):
#'       `verdict` (`"yes"` reached from every positive point, `"if"` reached exactly
#'       where `condition` holds, `"no"` not reached, `"unknown"` not decided),
#'       `limit` (`TRUE` when the face is invariant, so the
#'       zero is only approached as `eps -> +-Inf`), `certain` (`FALSE` when the
#'       block's invariant set is incomplete, which leaves anything but `"no"` an
#'       upper bound), `condition` (an R expression over the model's own coordinate
#'       names -- `eval(parse(text = condition), as.list(pars))` decides it at a
#'       point) and `at` (the degenerate point the
#'       orbit lands on; a primed name is a coordinate the remaining freedom leaves
#'       open). `NULL` unless `reportZeroCompatibility = TRUE`.}
#'     \item{`trafo`}{a complete [eqnvec] over all coordinates (identity plus the
#'       transversal pins and solved entries), ready for [P()] or
#'       `symmetryDetection(trafo = )`; `NULL` when nothing was reducible.}
#'     \item{`family`}{the general form: admissible transversals (or their matroid
#'       description) per reduced block, plus the gauge note that pins may be any
#'       nonzero constant.}
#'     \item{`removed`, `remaining`}{direction labels by outcome.}
#'     \item{`coordinates`, `fixed`, `settings`, `call`}{provenance.}
#'   }
#'   `print()` is terse: the verdict, the non-identity trafo entries, the invariant
#'   each outer parameter carries and where each zero limit is reached. `summary()`
#'   adds one line per block --
#'   kind, status, stage, how it was gauged -- plus, for a block that did not
#'   reduce, its invariants and the reason (how many invariants were found under
#'   which degree caps). `summary(verbose = TRUE)` adds the admissible gauges and
#'   the invariants of the reduced blocks. The per-stage certificates are kept in
#'   `$blocks[[i]]$certificates` and are not printed.
#'
#' @seealso [symmetryDetection()]
#' @example inst/examples/symmetryReduction.R
#' @export
symmetryReduction <- function(object, fixed = NULL, dPoly = 3L, dDarboux = 2L,
                           dExp = 2L, separable = TRUE,
                           reportZeroCompatibility = FALSE, verbose = FALSE, ...) {
  if (!inherits(object, "symmetrydetection"))
    stop("symmetryReduction(): `object` must be a symmetrydetection result.",
         call. = FALSE)
  if (length(list(...)))
    warning("symmetryReduction(): unused argument(s) ignored.", call. = FALSE)
  .require_ns("reticulate", "symmetryReduction()")
  .symCall <- match.call()
  separable <- isTRUE(separable)
  settings <- list(dPoly = as.integer(dPoly), dDarboux = as.integer(dDarboux),
                   dExp = as.integer(dExp), separable = separable,
                   primes = .symPrimes, verifyPrime = .symVerifyPrime)

  coords <- .symRedCoordinates(object)
  fixed <- unique(as.character(fixed))
  unknown <- setdiff(fixed, coords)
  if (length(unknown))
    warning("symmetryReduction(): no effect: ", paste(unknown, collapse = ", "),
            " -- not a coordinate of the analysis.", call. = FALSE)

  if (isTRUE(object$identifiable) || !length(object$symmetries))
    return(.symRedResult(object, list(), NULL, coords, fixed, settings, .symCall))

  code_dir <- system.file("code", package = "dMod2")
  sysmod <- reticulate::import("sys", convert = TRUE)
  if (!(code_dir %in% sysmod$path)) sysmod$path <- c(code_dir, sysmod$path)
  sd <- reticulate::import("symmetryDetection", convert = TRUE)
  spy <- reticulate::import("sympy", convert = TRUE)

  o <- .symOrdered(object)
  wr <- .symRedWeightRows(o$syms)

  # A scaling whose support overlaps a curved direction cannot be gauged on its
  # own: the curved invariants need not be scale-invariant, so the two reductions
  # would contradict each other (each assumes the other's coordinates fixed). Such
  # scalings join the curved block as generators (xi_i = w_i z_i, polynomial); the
  # scaling stage keeps only the components disjoint from every curved support.
  curvedIdx0 <- setdiff(seq_along(o$syms), wr$rows)
  curvedSupp <- unique(unlist(lapply(o$syms[curvedIdx0], .symCoords)))
  scalRows <- seq_len(nrow(wr$W))
  demoted <- integer(0)
  if (nrow(wr$W) && length(curvedSupp)) {
    repeat {                       # transitive closure: demotion can extend the overlap
      overlaps <- vapply(scalRows, function(r)
        any(colnames(wr$W)[wr$W[r, ] != 0L] %in% curvedSupp), logical(1))
      if (!any(overlaps)) break
      hit <- scalRows[overlaps]
      demoted <- c(demoted, hit)
      curvedSupp <- unique(c(curvedSupp,
        unlist(lapply(hit, function(r) colnames(wr$W)[wr$W[r, ] != 0L]))))
      scalRows <- setdiff(scalRows, hit)
    }
  }

  blocks <- list()
  if (length(scalRows)) {
    Wk <- wr$W[scalRows, , drop = FALSE]
    for (cp in .symRedComponents(Wk)) {
      rows <- scalRows[cp]
      Wb <- wr$W[rows, colSums(wr$W[rows, , drop = FALSE] != 0L) > 0L, drop = FALSE]
      if (isTRUE(verbose))
        message("scaling block {", paste(o$labels[wr$rows[rows]], collapse = ", "),
                "}: ", nrow(Wb), " direction(s)")
      blocks[[length(blocks) + 1L]] <-
        .symRedScalingBlock(Wb, o$labels[wr$rows[rows]], fixed, sd)
    }
  }
  curvedIdx <- sort(c(curvedIdx0, wr$rows[demoted]))
  if (length(curvedIdx))
    blocks <- c(blocks, .symRedCurved(o$syms, curvedIdx, o$labels, fixed,
                                      dPoly, dDarboux, dExp, sd, spy, verbose,
                                      separable))

  blocks <- .symRedVerify(blocks, sd)

  # solve the curved blocks' invariants into trafo entries, with every scaling pin
  # substituted in so a pinned coordinate cannot re-enter through a solution
  scalPins <- unlist(lapply(blocks, function(b)
    if (identical(b$type, "scaling")) b$pins))
  if (is.null(scalPins)) scalPins <- character(0)
  invStart <- 0L
  for (bi in seq_along(blocks)) {
    b <- blocks[[bi]]
    if (!identical(b$type, "curved") || !length(b$invariants)) next
    if (!is.null(b$target) && length(b$invariants) < b$target) next  # partial set:
    sol <- .symRedSolveInvariants(b, scalPins, spy, coords,     # gauge pin would be lossy
                                  invStart)
    if (sol$solved) {
      invStart <- invStart + length(sol$invNames)
      blocks[[bi]]$pins <- sol$pins
      blocks[[bi]]$transversal <- sol$gauge
      blocks[[bi]]$section <- sol$section
      gaugeVal <- if (length(sol$gauge)) unique(sol$pins[sol$gauge]) else character(0)
      blocks[[bi]]$coverage <- sol$coverage
      blocks[[bi]]$carrierDomain <- sol$carrierDomain
      blocks[[bi]]$gaugeNote <- paste(c(
        if (!is.null(sol$section))
          paste0("gauge section ", paste(sol$section, collapse = ", "))
        else paste0("gauge pin ", paste(sol$gauge, collapse = ", "), " = ",
                    paste(gaugeVal, collapse = ", ")),
        if (identical(sol$coverage, "partial"))
          paste0("entries certified positive for POSITIVE carrier values only -- a ",
                 "carrier that takes both signs leaves part of the positive orthant ",
                 "outside the chart")
        else "entries certified positive for every admissible outer value",
        sol$rootNote), collapse = "; ")
      blocks[[bi]]$survivorMeaning <- sol$meaning
      blocks[[bi]]$status <- "reduced"
      blocks[[bi]]$certificates <- c(b$certificates,
        "solved exactly: each invariant carried on a fresh q_<k> parameter",
        if (!is.null(sol$section))
          "section pre-certified: balance ratio strictly monotone along every orbit",
        if (identical(sol$coverage, "partial"))
          paste0("chart certified for POSITIVE carrier values only: a carrier that ",
                 "takes both signs leaves part of the positive orthant uncovered")
        else paste0("chart certified: every solved entry positive on the carrier ",
                    "domains (", paste(names(sol$carrierDomain), sol$carrierDomain,
                                       sep = " ", collapse = ", "), ")"),
        if (length(sol$shifted)) sprintf(paste0(
          "carrier offset(s) certified for %s: the shifted invariant exceeds ",
          "its offset on the whole positive orthant"),
          paste(sol$shifted, collapse = ", ")))
    } else {
      blocks[[bi]]$reason <- if (!is.null(sol$reason)) sol$reason else b$reason
      if (isTRUE(verbose))
        message("general block {", paste(b$labels, collapse = ", "), "}: ",
                blocks[[bi]]$reason, " -- reported as invariantOnly")
    }
  }

  # which coordinate the orbit can drive to zero -- read off the invariants, so it
  # runs after the solve and before the working fields are dropped
  if (isTRUE(reportZeroCompatibility)) for (bi in seq_along(blocks))
    blocks[[bi]]$zeroCompatibility <- tryCatch(
      .symRedZeroCompat(blocks[[bi]], spy, fixed, verbose),
      error = function(e) NULL)

  trafo <- .symRedAssembleTrafo(blocks, coords)
  blocks <- lapply(blocks, function(b) { b$preps <- NULL; b$Wres <- NULL
    b$invExps <- NULL; b })
  .symRedResult(object, blocks, trafo, coords, fixed, settings, .symCall)
}
