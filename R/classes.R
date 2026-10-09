## Function classes ------------------------------------------------------

#' Match Function Arguments to Choices
#'
#' Assigns unnamed entries of `arglist` to the entries of `choices` that are
#' not named in it, in order. Named entries not in `choices` are dropped.
#'
#' @param arglist list of arguments, as from `list(...)`.
#' @param choices character, the argument names to match.
#' @return Integer vector of the positions in `arglist` matching `choices`,
#'   `NA` where a choice is absent.
#' @keywords internal
#' @export
match.fnargs <- function(arglist, choices) {

  # Catch the case of names == NULL
  if (is.null(names(arglist))) names(arglist) <- rep("", length(arglist))

  # exclude named arguments which are not in choices
  arglist <- arglist[names(arglist) %in% c(choices, "")]

  # determine available arguments
  available <- choices %in% names(arglist)

  if (!all(available)) names(arglist)[names(arglist) == ""] <- choices[!available]

  if (any(duplicated(names(arglist)))) stop("duplicate arguments in prdfn/obsfn/parfn function call")

  mapping <- match(choices, names(arglist))
  return(mapping)

}


## Evaluation protocol for fn objects -------------------------------------
## The condition loop sits in the leaves, so a leaf sees all conditions at once.
## Conditions flow up through `*`: p2 gets `conditions`, p1 the names of p2's result.


# Structure descriptor ("leaf", "*", "+") in the closure's state env, or NULL
# for an fn without one. The `op` check excludes parfn's own `st`.
.fnNode <- function(f) {
  if (!is.function(f)) return(NULL)
  e <- environment(f)
  if (is.null(e)) return(NULL)
  st <- get0("st", envir = e, inherits = FALSE)
  if (is.null(st) || is.null(st$op)) return(NULL)
  st
}

# The raw kernels an fn evaluates, each as list(kernel, kind, condition),
# outermost factor first, restricted to the leaves answering for `condition`.
# This is the way from a composition to the closures that hold the controls.
.fnLeaves <- function(f, condition = NULL) {
  st <- .fnNode(f)
  if (is.null(st)) {
    m <- attr(f, "mappings")
    if (!length(m)) return(list())
    conds <- names(m)
    sel <- if (is.null(condition) || is.null(conds)) seq_along(m)
           else which(conds %in% condition)
    return(lapply(sel, function(i)
      list(kernel = m[[i]], kind = .fnKind(f),
           condition = if (is.null(conds)) NULL else conds[i])))
  }
  switch(st$op,
    leaf = {
      if (!is.null(condition) && !is.null(st$condition) &&
          !any(condition %in% st$condition))
        return(list())
      list(list(kernel = st$kernel, kind = st$kind, condition = st$condition))
    },
    "*" = c(.fnLeaves(st$p1, condition), .fnLeaves(st$p2, condition)),
    "+" = {
      idx <- if (is.null(condition)) seq_along(st$parts)
             else sort(unique(st$owner[intersect(condition, names(st$owner))]))
      unlist(lapply(st$parts[idx], .fnLeaves, condition = condition),
             recursive = FALSE)
    },
    list())
}


## ---- Condition resolution ------------------------------------------------

# Slots a leaf with condition `own` fills: all requested ones when `own` is
# NULL (one unnamed slot if none are), else those named `own`, and nothing when
# the request names other conditions only.
.resolveConditions <- function(conditions, own) {
  overlap <- test_conditions(conditions, own)
  # union() would drop repeats, and a request may name the same condition more
  # than once (quadrature nodes, parameter-frame rows).
  if (is.null(overlap))
    conditions <- if (is.null(own)) conditions
                  else if (is.null(conditions)) own
                  else union(own, conditions)
  slots <- if (is.null(own)) seq_len(max(1L, length(conditions)))
           else which(conditions %in% own)
  list(conditions = conditions,
       evaluate   = is.null(overlap) || length(overlap) > 0,
       slots      = slots)
}

# NULL holes are part of the contract: as.prdlist and do.call(c, .) rely on them.
.emptySlots <- function(conditions) {
  structure(vector("list", max(1L, length(conditions))), names = conditions)
}


## ---- Bundles -------------------------------------------------------------

# One entry per condition, or a single entry when `conds` is NULL. `shared`:
# all entries reference one request, so a leaf evaluates once and replicates.
# `times` is one vector, or a list of n for per-request grids.
.bundle <- function(conds = NULL, times = NULL, out = NULL, pars = NULL,
                    fixed = NULL, shared = FALSE) {
  list(conds = conds, times = times, out = out, pars = pars,
       fixed = fixed, shared = shared)
}

.bundle_n <- function(b) max(1L, length(b$conds))

.req_times <- function(b, i) if (is.list(b$times)) b$times[[i]] else b$times
.req_out   <- function(b, i) if (is.null(b$out))   NULL else b$out[[i]]
.req_pars  <- function(b, i) if (is.null(b$pars))  NULL else b$pars[[i]]
.req_fixed <- function(b, i) if (is.null(b$fixed)) NULL else b$fixed[[i]]

# Refcounted, so this costs pointers rather than copies.
.bundle_broadcast <- function(x, n) rep(list(x), n)

.bundle_subset <- function(b, sel) {
  if (is.null(b$conds)) return(b)
  idx <- match(sel, b$conds)
  idx <- idx[!is.na(idx)]
  .bundle(conds  = b$conds[idx],
          times  = if (is.list(b$times)) b$times[idx] else b$times,
          out    = if (is.null(b$out))   NULL else b$out[idx],
          pars   = if (is.null(b$pars))  NULL else b$pars[idx],
          fixed  = if (is.null(b$fixed)) NULL else b$fixed[idx],
          shared = b$shared)
}

.bundle_positions <- function(b, pos, conds_out) {
  if (is.null(b$conds))
    return(.bundle(conds = conds_out, times = b$times,
                   out   = if (is.null(b$out))   NULL else .bundle_broadcast(b$out[[1L]], length(pos)),
                   pars  = if (is.null(b$pars))  NULL else .bundle_broadcast(b$pars[[1L]], length(pos)),
                   fixed = if (is.null(b$fixed)) NULL else .bundle_broadcast(b$fixed[[1L]], length(pos)),
                   shared = b$shared))
  .bundle(conds  = b$conds[pos],
          times  = if (is.list(b$times)) b$times[pos] else b$times,
          out    = if (is.null(b$out))   NULL else b$out[pos],
          pars   = if (is.null(b$pars))  NULL else b$pars[pos],
          fixed  = if (is.null(b$fixed)) NULL else b$fixed[pos],
          shared = b$shared)
}

.bundle_from_call <- function(conditions, times, out, pars, fixed) {
  n <- max(1L, length(conditions))
  .bundle(conds  = conditions,
          times  = times,
          out    = if (is.null(out)) NULL else .bundle_broadcast(out, n),
          pars   = .bundle_broadcast(pars, n),
          fixed  = .bundle_broadcast(fixed, n),
          shared = TRUE)
}


## ---- Per-kind call shapes ------------------------------------------------

# inputs drives match.fnargs in the public shim; result decides prdlist wrapping.
.fnSpec <- list(
  obsfn = list(inputs = c("out", "pars"),   result = "prdlist"),
  prdfn = list(inputs = c("times", "pars"), result = "prdlist"),
  parfn = list(inputs = "pars",             result = "list"),
  objfn = list(inputs = "pars",             result = "objlist")
)

# How one element of p2's output becomes p1's (pars, fixed), per pair of kinds.
# Kept apart rather than unified: changing one needs its own test.
.handoff_prd_outerfixed <- function(v, fixed)          # obsfn * obsfn
  list(pars = attr(v, "parameters"), fixed = fixed)

.handoff_prd_innerfixed <- function(v, fixed) {        # obsfn * prdfn
  p <- attr(v, "parameters")
  f <- attr(p, "fixed")
  # Without a parameter derivative there is no transformation in between, so
  # the outer fixed parameters are the inner ones.
  if (is.null(attr(p, "deriv"))) f <- union(f, intersect(names(p), names(fixed)))
  list(pars = p, fixed = p[f])
}

.handoff_par_outerfixed <- function(v, fixed)          # obsfn * parfn
  list(pars = v, fixed = fixed)

.handoff_par_innerfixed <- function(v, fixed) {        # prdfn|parfn * parfn
  f <- attr(v, "fixed")
  list(pars = v[.setdiffU(names(v), f)], fixed = v[f])
}

.handoff_par_nofixed <- function(v, fixed)             # objfn * parfn
  list(pars = v, fixed = NULL)

.prodSpec <- list(
  "obsfn.obsfn" = list(out = "obsfn", handoff = ".handoff_prd_outerfixed", reduce = "c"),
  "obsfn.parfn" = list(out = "obsfn", handoff = ".handoff_par_outerfixed", reduce = "c"),
  "obsfn.prdfn" = list(out = "prdfn", handoff = ".handoff_prd_innerfixed", reduce = "c"),
  "prdfn.parfn" = list(out = "prdfn", handoff = ".handoff_par_innerfixed", reduce = "c"),
  "parfn.parfn" = list(out = "parfn", handoff = ".handoff_par_innerfixed", reduce = "c"),
  "objfn.parfn" = list(out = "objfn", handoff = ".handoff_par_nofixed",    reduce = "sum")
)

.fnKind <- function(f) {
  for (k in c("objfn", "obsfn", "prdfn", "parfn")) if (inherits(f, k)) return(k)
  NULL
}


## ---- Leaf evaluation -----------------------------------------------------

# Kernels expect disjoint pars / fixed.
.splitParsFixed <- function(pars, fixed) {
  if (is.null(fixed)) return(list(pars = pars, fixed = NULL))
  sub <- pars[.setdiffU(names(pars), names(fixed))]
  if (!inherits(sub, "parvec")) sub <- as.parvec(sub)   # `[.parvec` already did
  f <- as.numeric(fixed)
  names(f) <- names(fixed)
  class(f) <- c("parvec", "numeric")
  list(pars = sub, fixed = f)
}

# Report here, not three frames downstream.
.checkPrediction <- function(out, conditions) {
  # NaN passes: an observable can be undefined where no data sits (a ratio of
  # states that all start at 0), and normL2 stops on a NaN at a data point.
  if (all(is.finite(out))) return(invisible(NULL))
  bad <- (is.na(out) & !is.nan(out)) | is.infinite(out)
  if (!any(bad)) return(invisible(NULL))
  ai <- arrayInd(which(bad), dim(out))
  stop("Prediction is NA or Inf in condition ", paste0(conditions, collapse = ","),
       ".\nSubset of the prediction causing trouble:\n",
       paste0(capture.output(print(out[ai[, 1], c(1, ai[, 2])])), collapse = "\n"))
}

# `cond` is the slot's condition name; Pimpl uses it as warm-start key.
.callKernel <- function(st, b, i, cond, deriv, deriv2, keepStore = FALSE) {
  pf <- .splitParsFixed(.req_pars(b, i), .req_fixed(b, i))
  switch(st$kind,
    prdfn = if (keepStore)
              st$kernel(times = .req_times(b, i), pars = pf$pars, fixed = pf$fixed,
                        deriv = deriv, deriv2 = deriv2, keepStore = TRUE)
            else
              st$kernel(times = .req_times(b, i), pars = pf$pars, fixed = pf$fixed,
                        deriv = deriv, deriv2 = deriv2),
    obsfn = {
      o <- .req_out(b, i)
      .checkPrediction(o, cond)
      st$kernel(out = o, pars = pf$pars, fixed = pf$fixed,
                deriv = deriv, deriv2 = deriv2)
    },
    parfn = if (isTRUE(st$kernel_has_cond))
              st$kernel(pars = pf$pars, fixed = pf$fixed, deriv = deriv,
                        deriv2 = deriv2, condition = cond)
            else
              st$kernel(pars = pf$pars, fixed = pf$fixed, deriv = deriv,
                        deriv2 = deriv2))
}

# Batch entry when the leaf has one, else a loop. Not mclapply: prdframes hold
# 3-D and 4-D arrays whose trip through a fork pipe outweighs the solve.
.callKernelMany <- function(st, b, idx, conds, deriv, deriv2, cores,
                            keepStore = FALSE) {
  bf <- st$batchfn
  if (is.null(bf) || length(idx) < 2L)
    return(lapply(seq_along(idx), function(j)
      .callKernel(st, b, idx[j], conds[[j]], deriv, deriv2, keepStore)))

  split <- lapply(idx, function(i) .splitParsFixed(.req_pars(b, i), .req_fixed(b, i)))
  parsL  <- lapply(split, `[[`, "pars")
  fixedL <- lapply(split, `[[`, "fixed")

  res <- switch(st$kind,
    prdfn = if (keepStore)
              bf(times = if (is.list(b$times)) b$times[idx] else b$times,
                 parsList = parsL, fixedList = fixedL,
                 deriv = deriv, deriv2 = deriv2, cores = cores,
                 keepStore = TRUE)
            else
              bf(times = if (is.list(b$times)) b$times[idx] else b$times,
                 parsList = parsL, fixedList = fixedL,
                 deriv = deriv, deriv2 = deriv2, cores = cores),
    obsfn = {
      outL <- lapply(seq_along(idx), function(j) {
        o <- .req_out(b, idx[j]); .checkPrediction(o, conds[[j]]); o
      })
      bf(outList = outL, parsList = parsL, fixedList = fixedL,
         deriv = deriv, deriv2 = deriv2, cores = cores)
    },
    parfn = bf(parsList = parsL, fixedList = fixedL, deriv = deriv,
               deriv2 = deriv2, conditions = conds, cores = cores))

  if (isTRUE(getOption("dMod.batch.check", FALSE))) {
    ref <- lapply(seq_along(idx), function(j)
      .callKernel(st, b, idx[j], conds[[j]], deriv, deriv2, keepStore))
    cmp <- all.equal(res, ref, tolerance = 0)
    if (!isTRUE(cmp))
      stop("dMod.batch.check: batch entry of a ", st$kind,
           " leaf disagrees with the scalar kernel:\n  ",
           paste(cmp, collapse = "\n  "), call. = FALSE)
  }
  res
}

# Warm-start key for Pimpl. A leaf with its own condition keys by slot; an
# unspecific leaf answering several slots from one call keys by NULL.
.condKeys <- function(st, conds, shared) {
  if (!is.null(st$condition) || length(conds) == 1L || !shared) return(conds)
  rep(list(NULL), length(conds))
}

.evalLeaf <- function(st, b, deriv, deriv2, cores, keepStore = FALSE) {
  res <- .resolveConditions(b$conds, st$condition)
  outlist <- .emptySlots(res$conditions)
  if (!res$evaluate || length(res$slots) == 0L) return(outlist)

  slots <- res$slots
  shared <- b$shared || is.null(b$conds)
  cond_of_slot <- if (is.null(res$conditions)) list(NULL)
                  else .condKeys(st, as.list(res$conditions[slots]), shared)

  # One request behind every slot: evaluate once, replicate.
  if (shared) {
    r <- .callKernel(st, b, 1L, cond_of_slot[[1L]], deriv, deriv2, keepStore)
    for (s in slots) outlist[[s]] <- r
    return(outlist)
  }

  vals <- .callKernelMany(st, b, slots, cond_of_slot, deriv, deriv2, cores,
                          keepStore)
  for (j in seq_along(slots)) outlist[[slots[j]]] <- vals[[j]]
  outlist
}


## ---- Composition evaluation ---------------------------------------------

# p2 over every condition at once; its result names are p1's condition vector.
.evalProd <- function(st, b, deriv, deriv2, env, cores, hessian = TRUE,
                      sweep = "forward") {
  b2 <- .bundle(conds = b$conds, times = b$times, pars = b$pars,
                fixed = b$fixed, shared = b$shared,
                out = if (identical(st$p2kind, "obsfn")) b$out else NULL)
  inner <- .evalMany(st$p2, b2, deriv, deriv2, env, cores, hessian, sweep)

  conds <- names(inner)
  n <- max(1L, length(inner))
  handoff <- get(st$handoff, envir = asNamespace("dMod2"))
  hs <- lapply(seq_len(n), function(i) {
    v <- inner[[i]]
    if (is.null(v)) NULL else handoff(v, .req_fixed(b, min(i, .bundle_n(b))))
  })

  # p1 consumes p2's output when p2 yields a prdframe, the outer input when it
  # yields a parvec.
  p2_is_par <- identical(st$p2kind, "parfn")
  b1 <- .bundle(
    conds = conds,
    times = b$times,
    out   = if (p2_is_par) {
              if (is.null(b$out)) NULL else .bundle_broadcast(b$out[[1L]], n)
            } else inner,
    pars  = lapply(hs, function(h) if (is.null(h)) NULL else h$pars),
    fixed = lapply(hs, function(h) if (is.null(h)) NULL else h$fixed),
    shared = FALSE)

  res <- .evalMany(st$p1, b1, deriv, deriv2, env, cores, hessian, sweep)
  if (identical(st$reduce, "sum")) Reduce("+", res) else res
}

# One dispatch per PART: a sum of two g*x*p chains issues two batched solves.
.evalPlus <- function(st, b, deriv, deriv2, env, cores, hessian = TRUE,
                      sweep = "forward") {
  slotnames <- if (is.null(b$conds)) names(st$owner) else b$conds
  outlist <- .emptySlots(slotnames)
  own <- st$owner[slotnames]
  keep <- which(!is.na(own))
  if (!length(keep)) return(outlist)

  for (k in unique(own[keep])) {
    pos  <- keep[own[keep] == k]
    sub  <- .bundle_positions(b, pos, slotnames[pos])
    part <- .evalMany(st$parts[[k]], sub, deriv, deriv2, env, cores, hessian, sweep)
    for (j in seq_along(pos)) outlist[[pos[j]]] <- part[[j]]
  }
  outlist
}

.evalNode <- function(st, b, deriv, deriv2, env, cores, hessian = TRUE,
                      sweep = "forward") {
  switch(st$op,
    leaf = .evalLeaf(st, b, deriv, deriv2, cores),
    "*"  = .evalProd(st, b, deriv, deriv2, env, cores, hessian, sweep),
    "+"  = .evalPlus(st, b, deriv, deriv2, env, cores, hessian, sweep),
    stop(".evalNode: unknown node op '", st$op, "'.", call. = FALSE))
}

.evalMany <- function(f, b, deriv, deriv2, env, cores, hessian = TRUE,
                      sweep = "forward") {
  st <- .fnNode(f)
  if (is.null(st)) return(.evalLegacy(f, b, deriv, deriv2, env, hessian, sweep))
  .evalNode(st, b, deriv, deriv2, env, cores, hessian, sweep)
}


## ---- Reverse evaluation --------------------------------------------------
## Two phases: .fwdNode() returns values and a tape of every node's forward
## values, .bwdNode() walks the tape back with one .ct() cotangent per condition.

.fwdNode <- function(st, b, env, cores, deriv = FALSE) {
  switch(st$op,
    # A leaf keeps its checkpoints for the backward pass to replay. Under
    # second order it propagates the tangents each vjp is differentiated along.
    leaf = list(values = .evalLeaf(st, b, deriv, FALSE, cores,
                                   keepStore = isTRUE(st$keepstore)),
                tape   = list(op = "leaf", st = st, b = b)),
    "*"  = .fwdProd(st, b, env, cores, deriv),
    "+"  = .fwdPlus(st, b, env, cores, deriv),
    stop(".fwdNode: unknown node op '", st$op, "'.", call. = FALSE))
}

.fwdMany <- function(f, b, env, cores, deriv = FALSE) {
  st <- .fnNode(f)
  if (is.null(st))
    stop("reverse mode needs an fn object with a structure descriptor; this one ",
         "predates the evaluation protocol and has none.", call. = FALSE)
  .fwdNode(st, b, env, cores, deriv)
}

.fwdProd <- function(st, b, env, cores, deriv = FALSE) {
  b2 <- .bundle(conds = b$conds, times = b$times, pars = b$pars,
                fixed = b$fixed, shared = b$shared,
                out = if (identical(st$p2kind, "obsfn")) b$out else NULL)
  f2 <- .fwdNode(.fnNode(st$p2), b2, env, cores, deriv)
  inner <- f2$values

  conds <- names(inner)
  n <- max(1L, length(inner))
  handoff <- get(st$handoff, envir = asNamespace("dMod2"))
  hs <- lapply(seq_len(n), function(i) {
    v <- inner[[i]]
    if (is.null(v)) NULL else handoff(v, .req_fixed(b, min(i, .bundle_n(b))))
  })

  p2_is_par <- identical(st$p2kind, "parfn")
  b1 <- .bundle(
    conds = conds,
    times = b$times,
    out   = if (p2_is_par) {
              if (is.null(b$out)) NULL else .bundle_broadcast(b$out[[1L]], n)
            } else inner,
    pars  = lapply(hs, function(h) if (is.null(h)) NULL else h$pars),
    fixed = lapply(hs, function(h) if (is.null(h)) NULL else h$fixed),
    shared = FALSE)

  f1 <- .fwdNode(.fnNode(st$p1), b1, env, cores, deriv)
  values <- if (identical(st$reduce, "sum")) Reduce("+", f1$values) else f1$values
  list(values = values,
       tape = list(op = "*", st = st, b = b, t1 = f1$tape, t2 = f2$tape,
                   inner = inner, p2_is_par = p2_is_par, n = n))
}

.fwdPlus <- function(st, b, env, cores, deriv = FALSE) {
  slotnames <- if (is.null(b$conds)) names(st$owner) else b$conds
  outlist <- .emptySlots(slotnames)
  own <- st$owner[slotnames]
  keep <- which(!is.na(own))
  parts <- list()
  if (length(keep)) for (k in unique(own[keep])) {
    pos  <- keep[own[keep] == k]
    sub  <- .bundle_positions(b, pos, slotnames[pos])
    f    <- .fwdNode(.fnNode(st$parts[[k]]), sub, env, cores, deriv)
    for (j in seq_along(pos)) outlist[[pos[j]]] <- f$values[[j]]
    parts[[as.character(k)]] <- list(pos = pos, tape = f$tape)
  }
  list(values = outlist,
       tape = list(op = "+", st = st, b = b, parts = parts,
                   n = length(slotnames)))
}

# `seeds = TRUE` reads the trailing axis of every cotangent as independent
# first-order seeds instead of directions.
.bwdNode <- function(tape, w, env, cores, seeds = FALSE) {
  switch(tape$op,
    leaf = .bwdLeaf(tape, w, cores, seeds),
    "*"  = .bwdProd(tape, w, env, cores, seeds),
    "+"  = .bwdPlus(tape, w, env, cores, seeds),
    stop(".bwdNode: unknown node op '", tape$op, "'.", call. = FALSE))
}

# The leaf's own vjp; a leaf that answered no slot leaves a NULL hole. Slots
# behind one shared request have their cotangents summed and take one vjp, as
# the vjp is linear in its seed; separate requests go through the batch entry.
.bwdLeaf <- function(tape, w, cores, seeds = FALSE) {
  st <- tape$st; b <- tape$b
  vjp <- st$vjpfn
  if (is.null(vjp))
    stop("reverse mode: the ", st$kind, " leaf has no vjp entry. A prediction ",
         "needs odemodel(derivMode = c(\"forward\", \"reverse\")) and Xs(); ",
         "Xf() computes no derivatives in either direction, which is what it ",
         "is for. An observation or a ",
         "transformation needs derivMode = \"reverse\" and compile = TRUE.",
         call. = FALSE)

  res <- .resolveConditions(b$conds, st$condition)
  n <- max(1L, length(res$conditions))
  out <- vector("list", n)
  if (!res$evaluate || !length(res$slots)) return(out)

  shared <- b$shared || is.null(b$conds)
  live <- Filter(function(s) !.ct_null(w[[s]]), res$slots)
  if (!length(live)) return(out)

  # One vjp call at input index i, seeded with ws. A node reached only through
  # the parameters it passes on solves nothing.
  call_one <- function(i, ws, cond) {
    pf <- .splitParsFixed(.req_pars(b, i), .req_fixed(b, i))
    r <- if (!identical(st$kind, "parfn") && is.null(ws$out)) .ct() else
      switch(st$kind,
      prdfn = .ct(pars = vjp(times = .req_times(b, i), pars = pf$pars,
                             fixed = pf$fixed, cotangent = ws$out)),
      obsfn = vjp(out = .req_out(b, i), pars = pf$pars, fixed = pf$fixed,
                  cotangent = ws$out),
      parfn = .ct(pars = vjp(pars = pf$pars, fixed = pf$fixed,
                             cotangent = ws$pars, condition = cond)),
      stop(".bwdLeaf: no reverse mode for a ", st$kind, " leaf.", call. = FALSE))
    # Whatever the node passed through untouched keeps its cotangent.
    if (!identical(st$kind, "parfn")) {
      K <- max(.ctK(ws$pars), .ctK(ws$out))
      r <- .addCt(r, .ct(pars = .pickCotangent(ws$pars, names(pf$pars), K)))
    }
    r
  }

  # Seed mode: a prediction whose vjp takes `seeds` answers all in one sweep,
  # any other leaf is called once per seed and the answers are stacked.
  vjp_seeds <- identical(st$kind, "prdfn") && "seeds" %in% names(formals(vjp))
  call_seeds <- function(i, ws, cond) {
    if (vjp_seeds) {
      pf <- .splitParsFixed(.req_pars(b, i), .req_fixed(b, i))
      S  <- max(.ctK(ws$pars), .ctK(ws$out))
      r  <- if (is.null(ws$out)) .ct() else
        .ct(pars = vjp(times = .req_times(b, i), pars = pf$pars,
                       fixed = pf$fixed, cotangent = ws$out, seeds = TRUE))
      return(.addCt(r, .ct(pars = .pickCotangent(ws$pars, names(pf$pars), S))))
    }
    S <- max(.ctK(ws$pars), .ctK(ws$out))
    parts <- lapply(seq_len(S), function(k) call_one(i, .ct(
      out  = if (is.null(ws$out))  NULL else ws$out[, , k, drop = FALSE],
      pars = if (is.null(ws$pars)) NULL else ws$pars[, k, drop = FALSE]), cond))
    .stackSeeds(parts)
  }
  if (seeds) call_one_mode <- call_seeds else call_one_mode <- call_one

  if (shared) {
    ws <- Reduce(.addCt, lapply(live, function(s) w[[s]]))
    cond <- if (is.null(res$conditions)) NULL else res$conditions[live[1L]]
    out[[live[1L]]] <- call_one_mode(1L, ws, cond)
    return(out)
  }

  batchable <- identical(st$kind, "prdfn") && !is.null(st$vjpbatchfn) &&
               length(live) > 1L &&
               all(vapply(live, function(s) !is.null(w[[s]]$out), TRUE)) &&
               (!seeds || "seeds" %in% names(formals(st$vjpbatchfn)))
  if (batchable) {
    split <- lapply(live, function(s) .splitParsFixed(.req_pars(b, s),
                                                      .req_fixed(b, s)))
    vals <- st$vjpbatchfn(
      times     = if (is.list(b$times)) b$times[live] else b$times,
      parsList  = lapply(split, `[[`, "pars"),
      fixedList = lapply(split, `[[`, "fixed"),
      cotangentList = lapply(live, function(s) w[[s]]$out),
      conditions = if (is.null(res$conditions)) NULL else as.list(res$conditions[live]),
      cores     = cores,
      seeds     = seeds)
    if (isTRUE(getOption("dMod.batch.check", FALSE))) {
      ref <- lapply(seq_along(live), function(j)
        if (seeds)
          vjp(times = .req_times(b, live[j]), pars = split[[j]]$pars,
              fixed = split[[j]]$fixed, cotangent = w[[live[j]]]$out, seeds = TRUE)
        else
          vjp(times = .req_times(b, live[j]), pars = split[[j]]$pars,
              fixed = split[[j]]$fixed, cotangent = w[[live[j]]]$out))
      cmp <- all.equal(vals, ref, tolerance = 0)
      if (!isTRUE(cmp))
        stop("dMod.batch.check: the batched vjp of a ", st$kind,
             " leaf disagrees with the scalar one:\n  ",
             paste(cmp, collapse = "\n  "), call. = FALSE)
    }
    # As in call_one, the pass-through half is picked at the width the node
    # hands on, not at its own.
    for (j in seq_along(live)) {
      wj <- w[[live[j]]]
      K  <- max(.ctK(wj$pars), .ctK(wj$out))
      out[[live[j]]] <- .addCt(
        .ct(pars = vals[[j]]),
        .ct(pars = .pickCotangent(wj$pars, names(split[[j]]$pars), K)))
    }
    return(out)
  }

  for (s in live) {
    cond <- if (is.null(res$conditions)) NULL else res$conditions[s]
    out[[s]] <- call_one_mode(s, w[[s]], cond)
  }
  out
}

# Per-seed answers of a leaf stacked back onto one seed axis: the out-halves
# along their third axis, the pars-halves on the union of their rows.
.stackSeeds <- function(parts) {
  S <- length(parts)
  outs <- lapply(parts, `[[`, "out")
  out <- NULL
  if (!all(vapply(outs, is.null, TRUE))) {
    ref <- outs[[which(!vapply(outs, is.null, TRUE))[1L]]]
    d <- dim(ref)
    out <- array(0, c(d[1L], d[2L], S), dimnames = c(dimnames(ref)[1:2], list(NULL)))
    for (k in seq_len(S)) if (!is.null(outs[[k]])) out[, , k] <- outs[[k]][, , 1L]
  }
  pl <- lapply(parts, `[[`, "pars")
  rows <- unique(unlist(lapply(pl, rownames)))
  pars <- if (!length(rows)) NULL else
    do.call(cbind, lapply(pl, function(p) .pickCotangent(p, rows, 1L)))
  if (!is.null(pars)) dimnames(pars) <- list(rows, NULL)
  list(out = out, pars = pars)
}

.bwdProd <- function(tape, w, env, cores, seeds = FALSE) {
  st <- tape$st
  # A summed objective hands every term the same cotangent, because the sum's
  # derivative in each term is one.
  w1 <- if (identical(st$reduce, "sum")) rep(list(w[[1L]]), tape$n) else w
  u1 <- .bwdNode(tape$t1, w1, env, cores, seeds)

  if (tape$p2_is_par) {
    # p1 read p2's parvec as its parameters, and the outer `out` as its input.
    w2 <- lapply(u1, function(u) if (is.null(u)) NULL else .ct(pars = u$pars))
    down <- .bwdNode(tape$t2, w2, env, cores, seeds)
    outer_out <- lapply(u1, function(u) if (is.null(u)) NULL else .ct(out = u$out))
    return(.mergeCt(down, outer_out))
  }
  # p1 read p2's own output, values and parameters both.
  .bwdNode(tape$t2, u1, env, cores, seeds)
}

.bwdPlus <- function(tape, w, env, cores, seeds = FALSE) {
  out <- vector("list", tape$n)
  for (p in tape$parts) {
    u <- .bwdNode(p$tape, w[p$pos], env, cores, seeds)
    for (j in seq_along(p$pos)) out[[p$pos[j]]] <- u[[j]]
  }
  out
}

.mergeCt <- function(a, b) {
  n <- max(length(a), length(b))
  lapply(seq_len(n), function(i) .addCt(a[[i]], b[[i]]))
}

## ---- Public shim ---------------------------------------------------------

# Every fn object is a thin wrapper over its descriptor. `cores` must be a
# formal: match.fnargs drops unknown named arguments from `...`.
.fnWrap <- function(st) {
  function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, hessian = TRUE,
           conditions = st$default_conditions, env = NULL,
           cores = getOption("dMod.cores", 1L), sweep = "forward")
    .fnCall(st, list(...), fixed, deriv, deriv2, hessian, conditions, env, cores,
            sweep)
}

.fnCall <- function(st, arglist, fixed, deriv, deriv2, hessian, conditions, env,
                    cores, sweep = "forward") {
  spec <- .fnSpec[[st$kind]]
  arglist <- arglist[match.fnargs(arglist, spec$inputs)]
  names(arglist) <- spec$inputs
  b <- .bundle_from_call(conditions,
                         times = arglist$times, out = arglist$out,
                         pars  = arglist$pars,  fixed = fixed)
  out <- .evalNode(st, b, deriv, deriv2, env, cores, hessian, sweep)
  if (identical(spec$result, "prdlist")) as.prdlist(out) else out
}

# Evaluate a prediction chain for many parameter sets in one batch. One
# condition per parameter set, which may repeat; `times` is one grid, or a
# list of one per request.
.predictMany <- function(x, times, parsList, conditions, fixed = NULL,
                         deriv = TRUE, deriv2 = FALSE, env = NULL,
                         cores = getOption("dMod.cores", 1L)) {
  n <- length(parsList)
  if (length(conditions) != n)
    stop(".predictMany: one condition per parameter set is required.", call. = FALSE)
  b <- .bundle(conds = conditions, times = times, out = NULL,
               pars  = parsList,
               fixed = if (is.list(fixed)) fixed else .bundle_broadcast(fixed, n),
               shared = FALSE)
  st <- .fnNode(x)
  out <- if (is.null(st)) .evalLegacy(x, b, deriv, deriv2, env)
         else .evalNode(st, b, deriv, deriv2, env, cores)
  as.prdlist(out)
}


# Leaf descriptor. `kernel` is the raw P2X / X2Y / p2p, which stays the
# mapping every accessor already walks.
.leafState <- function(kernel, kind, condition) {
  list2env(list(op = "leaf", kind = kind, kernel = kernel,
                batchfn = attr(kernel, "batchfn"),
                vjpfn = attr(kernel, "vjpfn"),
                vjpbatchfn = attr(kernel, "vjpbatchfn"),
                keepstore = isTRUE(attr(kernel, "keepstore")),
                kernel_has_cond = "condition" %in% names(formals(kernel)),
                condition = condition, default_conditions = condition),
           parent = emptyenv())
}

# Drives an fn without a descriptor one condition at a time.
.evalLegacy <- function(f, b, deriv, deriv2, env, hessian = TRUE,
                        sweep = "forward") {
  kind <- .fnKind(f)
  # A request without conditions asks the fn for all of its own, so the result
  # keeps the names `.evalProd` reads back as p1's condition vector.
  conds <- if (is.null(b$conds)) attr(f, "conditions") else b$conds
  outlist <- .emptySlots(conds)
  nb <- .bundle_n(b)
  for (i in seq_len(max(1L, length(conds)))) {
    j <- min(i, nb)
    cond <- if (is.null(conds)) NULL else conds[i]
    r <- switch(kind,
      prdfn = f(times = .req_times(b, j), pars = .req_pars(b, j),
                fixed = .req_fixed(b, j), deriv = deriv, deriv2 = deriv2,
                conditions = cond, env = env),
      obsfn = f(out = .req_out(b, j), pars = .req_pars(b, j),
                fixed = .req_fixed(b, j), deriv = deriv, deriv2 = deriv2,
                conditions = cond, env = env),
      parfn = f(pars = .req_pars(b, j), fixed = .req_fixed(b, j),
                deriv = deriv, deriv2 = deriv2, conditions = cond, env = env),
      objfn = f(pars = .req_pars(b, j), fixed = .req_fixed(b, j),
                deriv = deriv, deriv2 = deriv2, hessian = hessian,
                conditions = cond, env = env, sweep = sweep),
      stop(".evalLegacy: cannot drive an fn of class ",
           paste(class(f), collapse = "/"), call. = FALSE))
    # an objfn returns its objlist directly, not a per-condition list
    outlist[[i]] <- if (identical(kind, "objfn")) r
                    else if (length(r) >= 1L) r[[1L]] else NULL
  }
  outlist
}


## ---- Composition bookkeeping ---------------------------------------------

# A `+` node contributes its parts, anything else is one part owning all of
# its conditions. Flattening keeps a long Reduce("+", .) one node.
.fnParts <- function(f) {
  st <- .fnNode(f)
  if (!is.null(st) && identical(st$op, "+"))
    return(list(parts = st$parts, owner = st$owner))
  conds <- attr(f, "conditions")
  list(parts = list(f), owner = setNames(rep(1L, length(conds)), conds))
}

# Overlapping conditions: the later operand wins.
.mergeOwnership <- function(x1, x2) {
  m1 <- attr(x1, "mappings"); m2 <- attr(x2, "mappings")
  if (is.null(names(m1)) || is.null(names(m2)))
    stop("General transformations (NULL names) cannot be coerced.")

  c1 <- attr(x1, "conditions"); c2 <- attr(x2, "conditions")
  overlap <- intersect(c1, c2)
  if (length(overlap) > 0) {
    warning(paste("Condition", overlap, "existed and has been overwritten."))
    m1 <- m1[!c1 %in% overlap]
    c1 <- c1[!c1 %in% overlap]
  }

  p1 <- .fnParts(x1); p2 <- .fnParts(x2)
  owner <- c(p1$owner[c1], p2$owner[c2] + length(p1$parts))
  parts <- c(p1$parts, p2$parts)

  keep  <- sort(unique(owner))               # drop parts nothing owns
  list(parts = parts[keep],
       owner = setNames(match(owner, keep), names(owner)),
       mappings = c(m1, m2),
       conditions = c(c1, c2))
}

# A condition-unspecific operand is asked with NULL, not with the composed
# condition name, otherwise getParameters/modelname return nothing.
.condFor <- function(f, cond) if (is.null(attr(f, "conditions"))) NULL else cond

# Read `what` off an operand's mapping for one condition.
.mapAttrAt <- function(f, cond, what) {
  m <- attr(f, "mappings")
  if (is.null(m) || !length(m)) return(NULL)
  sel <- if (is.null(cond) || is.null(names(m))) seq_along(m) else match(cond, names(m))
  sel <- sel[!is.na(sel)]
  for (i in sel) {
    v <- .kernelSetting(m[[i]], what)
    if (!is.null(v)) return(v)
  }
  NULL
}

# Metadata a composed mapping keeps. Without this getEquations, summary.*,
# Y(f = <composed>), compare() and petabExport all see NULL.
.composedMappingAttrs <- function(m, p1, p2, cond, p1kind, p2kind) {
  c1 <- .condFor(p1, cond); c2 <- .condFor(p2, cond)
  attr(m, "parameters") <- getParameters(p2, conditions = c2)
  attr(m, "modelname")  <- union(modelname(p1, conditions = c1),
                                 modelname(p2, conditions = c2))
  # equations: the prdfn-classed operand owns the state names, otherwise p2.
  eqsrc <- if (identical(p1kind, "prdfn")) p1 else if (identical(p2kind, "prdfn")) p2 else p2
  attr(m, "equations") <- .mapAttrAt(eqsrc, .condFor(eqsrc, cond), "equations")
  # forcings / events live on the prediction side only.
  for (what in c("forcings", "events", "states")) {
    v <- .mapAttrAt(p1, c1, what)
    if (is.null(v)) v <- .mapAttrAt(p2, c2, what)
    if (!is.null(v)) attr(m, what) <- v
  }
  m
}

# Per-condition callables for the consumers that walk `mappings`.
.composeMapping <- function(st, cond, kind) {
  force(cond)
  switch(kind,
    obsfn = function(out, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE,
                     cores = getOption("dMod.cores", 1L))
      .fnCall(st, list(out = out, pars = pars), fixed, deriv, deriv2, TRUE,
              cond, NULL, cores)[[1]],
    prdfn = function(times, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE,
                     cores = getOption("dMod.cores", 1L))
      .fnCall(st, list(times = times, pars = pars), fixed, deriv, deriv2, TRUE,
              cond, NULL, cores)[[1]],
    function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE,
             cores = getOption("dMod.cores", 1L))
      .fnCall(st, list(pars = pars), fixed, deriv, deriv2, TRUE,
              cond, NULL, cores)[[1]])
}

.composeMappings <- function(st, p1, p2, conditions, kind) {
  n <- max(1L, length(conditions))
  out <- lapply(seq_len(n), function(i) {
    cond <- if (is.null(conditions)) NULL else conditions[i]
    .composedMappingAttrs(.composeMapping(st, cond, kind), p1, p2, cond,
                          st$p1kind, st$p2kind)
  })
  setNames(out, conditions)
}

# Relabel a leaf so it answers for several conditions through one kernel. The
# descriptor is relabelled too, as `+` dispatches on it, not on `mappings`.
.fnWithConditions <- function(fn, conds) {
  st <- .fnNode(fn)
  if (is.null(st) || !identical(st$op, "leaf"))
    stop(".fnWithConditions: expects a leaf fn.", call. = FALSE)

  st2 <- list2env(as.list(st, all.names = TRUE), parent = emptyenv())
  st2$condition <- conds
  st2$default_conditions <- conds

  out <- .fnWrap(st2)
  for (a in c("parameters", "compileInfo", "resetWarmStart"))
    attr(out, a) <- attr(fn, a, exact = TRUE)
  attr(out, "mappings")   <- setNames(rep(list(st$kernel), length(conds)), conds)
  attr(out, "conditions") <- conds
  class(out) <- class(fn)
  out
}


# `forcings` on a summed fn is read opportunistically by compare().
.unionMappingAttr <- function(mappings, what) {
  vals <- lapply(mappings, attr, what)
  vals <- vals[!vapply(vals, is.null, logical(1))]
  if (!length(vals)) NULL else vals[[1]]
}


## General concatenation of functions ------------------------------------------

# The summands of an objective: the flat list `+` recorded, or the objective
# itself. A wrapper (%.*%, objfn * parfn) is one summand, `wrapped` is not read.
.objTerms <- function(f) {
  t <- attr(f, "terms", exact = TRUE)
  if (is.null(t)) list(f) else t
}

#' Direct Sum of Objective Functions
#'
#' @param e1 function of class `objfn`.
#' @param e2 function of class `objfn`.
#' @return Object of class `objfn` whose value, gradient and Hessian are the
#'   sums of those of `e1` and `e2`.
#' @seealso [normL2()], [constraintL2()], [datapointL2()], \code{\link{\%.*\%}}
#' @aliases sumobjfn
#' @example inst/examples/objective.R
#' @export
"+.objfn" <- function(e1, e2) {

  if (is.null(e1)) return(e2)

  conditions.x1 <- attr(e1, "conditions")
  conditions.x2 <- attr(e2, "conditions")
  conditions12 <- union(conditions.x1, conditions.x2)

  parameters.x1 <- attr(e1, "parameters")
  parameters.x2 <- attr(e2, "parameters")
  parameters12 <- union(parameters.x1, parameters.x2)

  modelname.x1 <- attr(e1, "modelname")
  modelname.x2 <- attr(e2, "modelname")
  modelname12 <- union(modelname.x1, modelname.x2)


  # objfn + objfn
  if (inherits(e1, "objfn") & inherits(e2, "objfn")) {

    outfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, hessian = NULL,
                      conditions = conditions12, env = NULL,
                      cores = getOption("dMod.cores", 1L),
                      sweep = "forward") {

      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("pars"))]
      pars <- arglist[[1]]

      # A term understands `sweep` only if it names it as a formal; one that
      # does not runs forward, and under an exact request still returns its
      # Hessian. Every wrapper in this package declares `sweep`.
      .call <- function(f, conds, e) {
        if (identical(sweep, "reverse") && "sweep" %in% names(formals(f)))
          f(pars = pars, fixed = fixed, deriv = deriv, deriv2 = deriv2,
            hessian = hessian, conditions = conds, env = e, cores = cores,
            sweep = sweep)
        else
          f(pars = pars, fixed = fixed, deriv = deriv, deriv2 = deriv2,
            hessian = if (identical(sweep, "reverse")) isTRUE(deriv2)
                      else hessian,
            conditions = conds, env = e, cores = cores)
      }
      # A term without conditions is evaluated once, any other on its
      # intersection with `conditions`, and not at all when that is empty.
      v1 <- v2 <- NULL
      if (is.null(conditions.x1)) {
        v1 <- .call(e1, conditions.x1, env)
      } else if (any(conditions %in% conditions.x1)) {
        v1 <- .call(e1, intersect(conditions, conditions.x1), env)
      }

      if (is.null(conditions.x2)) {
        v2 <- .call(e2, conditions.x2, env)
      } else if (any(conditions %in% conditions.x2)) {
        v2 <- .call(e2, intersect(conditions, conditions.x2), attr(v1, "env"))
      }

      # .sumobjlist adds an absent Hessian as zero, which is wrong when only
      # one term has one.
      .h <- function(v) !is.null(v) && !is.null(v$hessian)
      if (!is.null(v1) && !is.null(v2) && xor(.h(v1), .h(v2)))
        stop("a summed objective got a Hessian from ",
             if (.h(v1)) "its first" else "its second", " term and none from ",
             if (.h(v1)) "its second" else "its first",
             ". Adding them would drop that term's curvature silently.",
             call. = FALSE)

      out <- v1 + v2
      attr(out, "env") <- attr(v1, "env")
      return(out)
    }

    class(outfn) <- c("objfn", "fn")
    attr(outfn, "conditions") <- conditions12
    attr(outfn, "parameters") <- parameters12
    attr(outfn, "modelname") <- modelname12
    # Reconstruction handles, coalesced from either operand.
    for (.a in c("prdfn", "data", "errfn", "timesD")) {
      .v <- attr(e1, .a, exact = TRUE)
      if (is.null(.v)) .v <- attr(e2, .a, exact = TRUE)
      if (!is.null(.v)) attr(outfn, .a) <- .v
    }
    # l2spec is CONCATENATED: every L2 term keeps its own data, prediction and
    # error model, which is what reml() needs from a split objective.
    attr(outfn, "l2spec") <- c(attr(e1, "l2spec", exact = TRUE),
                               attr(e2, "l2spec", exact = TRUE))
    # The summands themselves, flat however the sum was nested, so a caller
    # can take a data term apart from the priors beside it.
    attr(outfn, "terms") <- c(.objTerms(e1), .objTerms(e2))
    return(outfn)

  }


}


#' Multiplication of Objective Functions with Scalars
#'
#' @description \code{x1 \%.*\% x2} multiplies an object of class `objfn` or
#' `objlist` by a scalar.
#'
#' @param x1 numeric of length one.
#' @param x2 object of class `objfn` or `objlist`.
#' @return For an `objfn`, an `objfn` whose value, gradient and Hessian are
#'   scaled by `x1`. For an `objlist`, the `objlist` with every numeric entry
#'   and attribute scaled. Otherwise `x1 * x2`.
#' @seealso [+.objfn]
#' @examples
#' obj <- constraintL2(mu = c(a = 0, b = 0), sigma = 1)
#' obj2 <- 2 %.*% obj
#' obj2(c(a = 1, b = 2))$value
#' @export
"%.*%" <- function(x1, x2) {

  if (inherits(x2, "objlist")) {

    out <- lapply(x2, function(x) {
      x1*x
    })
    # Multiply attributes
    out2.attributes <- attributes(x2)[sapply(attributes(x2), is.numeric)]
    attr.names <- names(out2.attributes)
    out.attributes <- lapply(attr.names, function(n) {
      x1*attr(x2, n)
    })
    attributes(out) <- attributes(x2)
    attributes(out)[attr.names] <- out.attributes

    return(out)


  } else if (inherits(x2, "objfn")) {

    conditions12 <- attr(x2, "conditions")
    parameters12 <- attr(x2, "parameters")
    modelname12 <- attr(x2, "modelname")
    outfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, hessian = NULL,
                      conditions = conditions12, env = NULL,
                      cores = getOption("dMod.cores", 1L),
                      sweep = "forward") {

      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("pars"))]
      pars <- arglist[[1]]

      # A scaled objective is still the same objective, so the direction goes
      # through; a term that cannot take it keeps the forward one.
      v2 <- if (identical(sweep, "reverse") && "sweep" %in% names(formals(x2)))
        x2(pars = pars, fixed = fixed, deriv = deriv, deriv2 = deriv2,
           hessian = hessian, conditions = conditions, env = env, cores = cores,
           sweep = sweep)
      else
        x2(pars = pars, fixed = fixed, deriv = deriv, deriv2 = deriv2,
           hessian = if (identical(sweep, "reverse")) isTRUE(deriv2)
                     else hessian,
           conditions = conditions, env = env, cores = cores)

      out <- x1 %.*% v2
      attr(out, "env") <- attr(v2, "env")
      return(out)
    }

    class(outfn) <- c("objfn", "fn")
    attr(outfn, "conditions") <- conditions12
    attr(outfn, "parameters") <- parameters12
    attr(outfn, "modelname") <- modelname12
    # The objective inside, for controls(). Not `terms`: a scaled objective
    # is not a sum.
    attr(outfn, "wrapped") <- list(x2)
    return(outfn)

  } else {

    x1*x2

  }

}


#' Direct Sum of Functions
#'
#' Adds prediction functions, parameter transformations or observation
#' functions of the same class.
#'
#' @param e1 function of class `obsfn`, `prdfn` or `parfn`.
#' @param e2 function of the same class as `e1`.
#' @details Each function is defined for a set of conditions. The sum is
#' defined for their union. A condition present in both is taken from `e2`,
#' with a warning.
#' @return Object of the same class as `e1` and `e2`, defined for the union of
#' conditions.
#' @aliases sumfn
#' @seealso [P()], [Y()], [Xs()], [*.fn]
#' @example inst/examples/prediction.R
#' @export
"+.fn" <- function(e1, e2) {

  if (is.null(e1)) return(e2)

  k1 <- .fnKind(e1); k2 <- .fnKind(e2)
  if (is.null(k1) || is.null(k2) || !identical(k1, k2))
    stop("\"+.fn\": cannot add ", paste(class(e1), collapse = "/"), " and ",
         paste(class(e2), collapse = "/"), ".", call. = FALSE)

  own <- .mergeOwnership(e1, e2)

  st <- list2env(list(op = "+", kind = k1, parts = own$parts, owner = own$owner,
                      default_conditions = own$conditions), parent = emptyenv())
  outfn <- .fnWrap(st)

  attr(outfn, "mappings")    <- own$mappings
  attr(outfn, "parameters")  <- union(attr(e1, "parameters"), attr(e2, "parameters"))
  attr(outfn, "compileInfo") <- .mergeCompileInfo(attr(e1, "compileInfo"),
                                                  attr(e2, "compileInfo"))
  attr(outfn, "conditions")  <- own$conditions
  attr(outfn, "forcings")    <- .unionMappingAttr(own$mappings, "forcings")

  # Keep "composed" only when a composed operand went in, so summary() keeps
  # its detail branch for a sum of leaves and drops it for a sum of chains.
  cls <- c(k1, "fn")
  if (inherits(e1, "composed") || inherits(e2, "composed")) cls <- c(cls, "composed")
  class(outfn) <- cls

  outfn

}


#' Direct Sum of Datasets
#'
#' Combines two datalists.
#'
#' @param e1 object of class `datalist`.
#' @param e2 object of class `datalist`.
#' @details A condition present in both datalists is taken from `e2`, with
#' a warning. The condition grids are combined.
#' @return Object of class `datalist` for the union of conditions.
#' @seealso [as.datalist()]
#' @aliases sumdatalist
#' @example inst/examples/sumdatalist.R
#' @export
"+.datalist" <- function(e1, e2) {

  overlap <- names(e2)[names(e2) %in% names(e1)]
  if (length(overlap) > 0) {
    warning(paste("Condition", overlap, "existed and has been overwritten."))
    e1 <- e1[!names(e1) %in% names(e2)]
  }

  conditions <- union(names(e1), names(e2))
  data <- lapply(conditions, function(C) rbind(e1[[C]], e2[[C]]))
  names(data) <- conditions

  grid1 <- attr(e1, "condition.grid")
  grid2 <- attr(e2, "condition.grid")

  grid <- combine(grid1, grid2)




  if (is.data.frame(grid)) grid <- grid[!duplicated(rownames(grid)), , drop = FALSE]

  out <- as.datalist(data)
  attr(out, "condition.grid") <- grid

  return(out)
}

out_conditions <- function(c1, c2) {

  if (!is.null(c1)) return(c1)
  if (!is.null(c2)) return(c2)
  return(NULL)

}

test_conditions <- function(c1, c2) {
  if (is.null(c1)) return(NULL)
  if (is.null(c2)) return(NULL)
  .intersectU(c1, c2)
}

#' Concatenation of Functions
#'
#' Concatenates observation functions, prediction functions and parameter
#' transformations: `(e1 * e2)(times, pars)` evaluates `e2` first and passes
#' its output to `e1`.
#'
#' @param e1 function of class `objfn`, `obsfn`, `prdfn`, `parfn` or `idfn`.
#' @param e2 function of class `obsfn`, `prdfn`, `parfn` or `idfn`.
#' @details Both functions must be defined for the same conditions, or one of
#' them for no specific condition (`conditions = NULL`).
#' @return `obsfn * obsfn` and `obsfn * parfn` return an `obsfn`,
#'   `obsfn * prdfn` and `prdfn * parfn` a `prdfn`, `parfn * parfn` a `parfn`
#'   and `objfn * parfn` an `objfn`. `Id() * f` and `f * Id()` return `f`.
#' @seealso [+.fn], [Id()]
#' @aliases prodfn
#' @example inst/examples/prediction.R
#' @export
"*.fn" <- function(e1, e2) {

  # An unspecific function combines with any other; two specific ones must
  # share their conditions, and one condition never combines with several.

  conditions.p1 <- attr(e1, "conditions")
  conditions.p2 <- attr(e2, "conditions")

  is_unspecific <- function(x) is.null(x)
  is_specific   <- function(x) !is.null(x) && length(x) == 1
  is_multiple   <- function(x) !is.null(x) && length(x) > 1

  if (!is_unspecific(conditions.p1) &&
      !is_unspecific(conditions.p2)) {

    # one specific, one multiple -> forbidden
    if ((is_specific(conditions.p1) && is_multiple(conditions.p2)) ||
        (is_specific(conditions.p2) && is_multiple(conditions.p1))) {

      stop(
        "Invalid composition of functions:\n",
        "Incompatible condition sets.\n\n",
        "Left-hand function conditions:  ",
        paste(conditions.p1, collapse = ", "), "\n",
        "Right-hand function conditions: ",
        paste(conditions.p2, collapse = ", "), "\n\n",
        "A function defined for a single condition cannot be\n",
        "combined with a function defined for multiple conditions.\n",
        "Either both functions must cover all conditions,\n",
        "or one function must be condition-unspecific."
      )
    }
  }

  if (inherits(e1, "idfn")) return(e2)
  if (inherits(e2, "idfn")) return(e1)

  key  <- paste(.fnKind(e1), .fnKind(e2), sep = ".")
  spec <- if (length(key) == 1L) .prodSpec[[key]] else NULL
  if (is.null(spec))
    stop("\"*.fn\": no composition defined for ",
         paste(class(e1), collapse = "/"), " * ",
         paste(class(e2), collapse = "/"), ".", call. = FALSE)

  conditions.out <- out_conditions(conditions.p1, conditions.p2)

  st <- list2env(list(op = "*", kind = spec$out, p1 = e1, p2 = e2,
                      p1kind = .fnKind(e1), p2kind = .fnKind(e2),
                      handoff = spec$handoff, reduce = spec$reduce,
                      default_conditions = NULL), parent = emptyenv())
  outfn <- .fnWrap(st)

  attr(outfn, "conditions")  <- conditions.out
  attr(outfn, "parameters")  <- attr(e2, "parameters")
  attr(outfn, "compileInfo") <- .mergeCompileInfo(attr(e1, "compileInfo"),
                                                  attr(e2, "compileInfo"))

  if (identical(spec$out, "objfn")) {
    # An objfn has no mappings; without these an objfn * parfn loses its
    # parameter set, its model name and the reconstruction handles.
    attr(outfn, "modelname") <- union(attr(e1, "modelname"), attr(e2, "modelname"))
    for (.a in c("data", "errfn", "timesD")) {
      .v <- attr(e1, .a, exact = TRUE)
      if (!is.null(.v)) attr(outfn, .a) <- .v
    }
    # The reconstructed prediction has to live in the outer coordinates.
    .prd <- attr(e1, "prdfn", exact = TRUE)
    if (!is.null(.prd))
      attr(outfn, "prdfn") <- tryCatch(.prd * e2, error = function(e) .prd)
    .l2 <- attr(e1, "l2spec", exact = TRUE)
    if (!is.null(.l2))
      attr(outfn, "l2spec") <- lapply(.l2, function(tm) {
        tm$prdfn <- tryCatch(tm$prdfn * e2, error = function(e) tm$prdfn)
        tm
      })
    # The objective inside, for controls(). Not `terms`: the composition
    # reads other parameters than the objective it wraps.
    attr(outfn, "wrapped") <- list(e1)
  } else {
    attr(outfn, "mappings") <- .composeMappings(st, e1, e2, conditions.out, spec$out)
  }

  class(outfn) <- c(spec$out, "fn", "composed")

  outfn

}
