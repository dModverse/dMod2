## Function classes ------------------------------------------------------

#' dMod match function arguments
#' 
#' The function is exported for dependency reasons
#' 
#' @param arglist list
#' @param choices character
#' 
#' @export
match.fnargs <- function(arglist, choices) {

  # Catch the case of names == NULL
  if (is.null(names(arglist))) names(arglist) <- rep("", length(arglist))

  # exlude named arguments which are not in choices
  arglist <- arglist[names(arglist) %in% c(choices, "")]

  # determine available arguments
  available <- choices %in% names(arglist)

  if (!all(available)) names(arglist)[names(arglist) == ""] <- choices[!available]

  if (any(duplicated(names(arglist)))) stop("duplicate arguments in prdfn/obsfn/parfn function call")

  mapping <- match(choices, names(arglist))
  return(mapping)

}


## Evaluation protocol for fn objects -------------------------------------
##
## Moves the loop over conditions from outside a chain into the leaves, so a
## leaf sees every condition at once and can batch them.
##
##   st        structure descriptor ("leaf" / "*" / "+") in a state env, the
##             pattern parfn already uses. Not an attribute: modelname<-.fn and
##             the PEtab relabeller rewrite `mappings` by hand.
##   bundle    one (times, pars, fixed, out) per condition, plus `shared`.
##   .evalMany recursive evaluator; shared bundle + one condition reduces to
##             the pre-rebuild call sequence.
##
## Conditions flow UP through `*`: "*.fn" passes `conditions` only to p2 and
## calls p1 with names(<p2 result>), which is NULL when both are unspecific.


# Structure descriptor, or NULL for objects built before the rebuild. The
# `op` check also excludes parfn's legacy `st`.
.fnNode <- function(f) {
  if (!is.function(f)) return(NULL)
  e <- environment(f)
  if (is.null(e)) return(NULL)
  st <- get0("st", envir = e, inherits = FALSE)
  if (is.null(st) || is.null(st$op)) return(NULL)
  st
}


## ---- Condition resolution ------------------------------------------------

# The truth table every leaf reproduces (own = leaf's condition):
#
#   own    conditions    -> slots filled
#   NULL   NULL             one unnamed slot
#   NULL   c("C1","C2")     both, same result replicated
#   "C1"   NULL             one slot named "C1"
#   "C1"   c("C1","C2")     "C1" only; "C2" stays NULL
#   "C1"   "C2"             nothing evaluated; one NULL slot
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

# One entry per condition, or a single entry when `conds` is NULL.
#
# shared = TRUE: pars/fixed/out are n references to one object. A leaf must
# then evaluate once and replicate -- batching would turn one solve into n.
#
# times is one vector for all requests, or a list of n for per-request grids;
# every composition call site shares one.
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

# How one element of p2's output becomes p1's (pars, fixed). Five behaviours
# across six branches; named rather than unified because at least one is
# probably wrong and changing that needs its own oracle. Lines are pre-rebuild
# classes.R.
.handoff_prd_outerfixed <- function(v, fixed)          # obsfn * obsfn  (:464)
  list(pars = attr(v, "parameters"), fixed = fixed)

.handoff_prd_innerfixed <- function(v, fixed) {        # obsfn * prdfn  (:568)
  p <- attr(v, "parameters")
  list(pars = p, fixed = p[attr(p, "fixed")])
}

.handoff_par_outerfixed <- function(v, fixed)          # obsfn * parfn  (:516)
  list(pars = v, fixed = fixed)

.handoff_par_innerfixed <- function(v, fixed) {        # prdfn|parfn * parfn (:627)
  f <- attr(v, "fixed")
  list(pars = v[setdiff(names(v), f)], fixed = v[f])
}

.handoff_par_nofixed <- function(v, fixed)             # objfn * parfn  (:730)
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
  sub <- pars[setdiff(names(pars), names(fixed))]
  if (!inherits(sub, "parvec")) sub <- as.parvec(sub)   # `[.parvec` already did
  f <- as.numeric(fixed)
  names(f) <- names(fixed)
  class(f) <- c("parvec", "numeric")
  list(pars = sub, fixed = f)
}

# Report here, not three frames downstream.
.checkPrediction <- function(out, conditions) {
  bad <- is.na(out) | is.infinite(out)
  if (!any(bad)) return(invisible(NULL))
  ai <- arrayInd(which(bad), dim(out))
  stop("Prediction is NA or Inf in condition ", paste0(conditions, collapse = ","),
       ".\nSubset of the prediction causing trouble:\n",
       paste0(capture.output(print(out[ai[, 1], c(1, ai[, 2])])), collapse = "\n"))
}

# `cond` is the slot's condition name; Pequil/Pimpl use it as warm-start key.
.callKernel <- function(st, b, i, cond, deriv, deriv2) {
  pf <- .splitParsFixed(.req_pars(b, i), .req_fixed(b, i))
  switch(st$kind,
    prdfn = st$kernel(times = .req_times(b, i), pars = pf$pars, fixed = pf$fixed,
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

# Batch entry when the leaf has one, else a loop. Not mclapply: prdframes carry
# 3-D and 4-D arrays whose trip through a fork pipe outweighs the solve.
.callKernelMany <- function(st, b, idx, conds, deriv, deriv2, cores) {
  bf <- st$batchfn
  if (is.null(bf) || length(idx) < 2L)
    return(lapply(seq_along(idx), function(j)
      .callKernel(st, b, idx[j], conds[[j]], deriv, deriv2)))

  split <- lapply(idx, function(i) .splitParsFixed(.req_pars(b, i), .req_fixed(b, i)))
  parsL  <- lapply(split, `[[`, "pars")
  fixedL <- lapply(split, `[[`, "fixed")

  res <- switch(st$kind,
    prdfn = bf(times = if (is.list(b$times)) b$times[idx] else b$times,
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
      .callKernel(st, b, idx[j], conds[[j]], deriv, deriv2))
    cmp <- all.equal(res, ref, tolerance = 0)
    if (!isTRUE(cmp))
      stop("dMod.batch.check: batch entry of a ", st$kind,
           " leaf disagrees with the scalar kernel:\n  ",
           paste(cmp, collapse = "\n  "), call. = FALSE)
  }
  res
}

# Warm-start key for Pequil/Pimpl. A leaf with its own condition keys by slot;
# an unspecific leaf answering several slots from ONE call has no single
# condition and keys by NULL, as before the rebuild.
.condKeys <- function(st, conds, shared) {
  if (!is.null(st$condition) || length(conds) == 1L || !shared) return(conds)
  rep(list(NULL), length(conds))
}

.evalLeaf <- function(st, b, deriv, deriv2, cores) {
  res <- .resolveConditions(b$conds, st$condition)
  outlist <- .emptySlots(res$conditions)
  if (!res$evaluate || length(res$slots) == 0L) return(outlist)

  slots <- res$slots
  shared <- b$shared || is.null(b$conds)
  cond_of_slot <- if (is.null(res$conditions)) list(NULL)
                  else .condKeys(st, as.list(res$conditions[slots]), shared)

  # One request behind every slot: evaluate once, replicate.
  if (shared) {
    r <- .callKernel(st, b, 1L, cond_of_slot[[1L]], deriv, deriv2)
    for (s in slots) outlist[[s]] <- r
    return(outlist)
  }

  vals <- .callKernelMany(st, b, slots, cond_of_slot, deriv, deriv2, cores)
  for (j in seq_along(slots)) outlist[[slots[j]]] <- vals[[j]]
  outlist
}


## ---- Composition evaluation ---------------------------------------------

# p2 over every condition at once; its result names are p1's condition vector.
.evalProd <- function(st, b, deriv, deriv2, env, cores) {
  b2 <- .bundle(conds = b$conds, times = b$times, pars = b$pars,
                fixed = b$fixed, shared = b$shared,
                out = if (identical(st$p2kind, "obsfn")) b$out else NULL)
  inner <- .evalMany(st$p2, b2, deriv, deriv2, env, cores)

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

  res <- .evalMany(st$p1, b1, deriv, deriv2, env, cores)
  if (identical(st$reduce, "sum")) Reduce("+", res) else res
}

# One dispatch per PART: a sum of two g*x*p chains issues two batched solves.
.evalPlus <- function(st, b, deriv, deriv2, env, cores) {
  slotnames <- if (is.null(b$conds)) names(st$owner) else b$conds
  outlist <- .emptySlots(slotnames)
  own <- st$owner[slotnames]
  keep <- which(!is.na(own))
  if (!length(keep)) return(outlist)

  for (k in unique(own[keep])) {
    pos  <- keep[own[keep] == k]
    sub  <- .bundle_positions(b, pos, slotnames[pos])
    part <- .evalMany(st$parts[[k]], sub, deriv, deriv2, env, cores)
    for (j in seq_along(pos)) outlist[[pos[j]]] <- part[[j]]
  }
  outlist
}

.evalNode <- function(st, b, deriv, deriv2, env, cores) {
  switch(st$op,
    leaf = .evalLeaf(st, b, deriv, deriv2, cores),
    "*"  = .evalProd(st, b, deriv, deriv2, env, cores),
    "+"  = .evalPlus(st, b, deriv, deriv2, env, cores),
    stop(".evalNode: unknown node op '", st$op, "'.", call. = FALSE))
}

.evalMany <- function(f, b, deriv, deriv2, env, cores) {
  st <- .fnNode(f)
  if (is.null(st)) return(.evalLegacy(f, b, deriv, deriv2, env))
  .evalNode(st, b, deriv, deriv2, env, cores)
}


## ---- Public shim ---------------------------------------------------------

# Every fn object is this: a thin wrapper over its descriptor. The signature
# is the pre-rebuild one plus `cores`, which must be a formal -- match.fnargs
# drops named arguments it does not know, so a `cores` in `...` would vanish.
.fnWrap <- function(st) {
  function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE,
           conditions = st$default_conditions, env = NULL,
           cores = getOption("dMod.cores", 1L))
    .fnCall(st, list(...), fixed, deriv, deriv2, conditions, env, cores)
}

.fnCall <- function(st, arglist, fixed, deriv, deriv2, conditions, env, cores) {
  spec <- .fnSpec[[st$kind]]
  arglist <- arglist[match.fnargs(arglist, spec$inputs)]
  names(arglist) <- spec$inputs
  b <- .bundle_from_call(conditions,
                         times = arglist$times, out = arglist$out,
                         pars  = arglist$pars,  fixed = fixed)
  out <- .evalNode(st, b, deriv, deriv2, env, cores)
  if (identical(spec$result, "prdlist")) as.prdlist(out) else out
}

# Evaluate a prediction chain for many parameter sets in one batch.
# `conditions` names one condition per parameter set and MAY repeat: quadrature
# nodes and parameter-frame rows both send the same condition several times.
# `times` is one grid, or a list of one per request.
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
                kernel_has_cond = "condition" %in% names(formals(kernel)),
                condition = condition, default_conditions = condition),
           parent = emptyenv())
}

# Pre-rebuild path: drive an fn without a descriptor one condition at a time.
.evalLegacy <- function(f, b, deriv, deriv2, env) {
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
                deriv = deriv, deriv2 = deriv2, conditions = cond, env = env),
      stop(".evalLegacy: cannot drive an fn of class ",
           paste(class(f), collapse = "/"), call. = FALSE))
    # an objfn returns its objlist directly, not a per-condition list
    outlist[[i]] <- if (identical(kind, "objfn")) r
                    else if (length(r) >= 1L) r[[1L]] else NULL
  }
  outlist
}


## ---- Composition bookkeeping ---------------------------------------------

# A `+` node contributes its parts; anything else is one part owning all of
# its conditions. Flattening keeps a+b+c one node rather than two nested ones,
# which matters for P()'s Reduce("+", .) over thousands of conditions.
.fnParts <- function(f) {
  st <- .fnNode(f)
  if (!is.null(st) && identical(st$op, "+"))
    return(list(parts = st$parts, owner = st$owner))
  conds <- attr(f, "conditions")
  list(parts = list(f), owner = setNames(rep(1L, length(conds)), conds))
}

# Overlapping conditions: later operand wins, as before the rebuild.
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
# condition name -- otherwise getParameters/modelname return nothing.
.condFor <- function(f, cond) if (is.null(attr(f, "conditions"))) NULL else cond

# Read `what` off an operand's mapping for one condition.
.mapAttrAt <- function(f, cond, what) {
  m <- attr(f, "mappings")
  if (is.null(m) || !length(m)) return(NULL)
  sel <- if (is.null(cond) || is.null(names(m))) seq_along(m) else match(cond, names(m))
  sel <- sel[!is.na(sel)]
  for (i in sel) {
    v <- attr(m[[i]], what)
    if (!is.null(v)) return(v)
  }
  NULL
}

# Metadata a composed mapping carries. Without this getEquations, summary.*,
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

# Per-condition callable, kept so the 25+ consumers that walk `mappings` and
# the external callers that invoke them keep working.
.composeMapping <- function(st, cond, kind) {
  force(cond)
  switch(kind,
    obsfn = function(out, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE,
                     cores = getOption("dMod.cores", 1L))
      .fnCall(st, list(out = out, pars = pars), fixed, deriv, deriv2,
              cond, NULL, cores)[[1]],
    prdfn = function(times, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE,
                     cores = getOption("dMod.cores", 1L))
      .fnCall(st, list(times = times, pars = pars), fixed, deriv, deriv2,
              cond, NULL, cores)[[1]],
    function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE,
             cores = getOption("dMod.cores", 1L))
      .fnCall(st, list(pars = pars), fixed, deriv, deriv2,
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

# Relabel a leaf so it answers for several conditions through one kernel.
# Rewriting only the `mappings` attribute is not enough: `+` dispatches on the
# descriptor, and a leaf that still says condition = NULL would evaluate once
# and replicate instead of answering per condition.
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

#' Direct sum of objective functions
#'
#' @param x1 function of class `objfn`
#' @param x2 function of class `objfn`
#' @details The objective functions are evaluated and their results as added. Sometimes,
#' the evaluation of an objective function depends on results that have been computed
#' internally in a preceding objective function. Therefore, environments are forwarded
#' and all evaluations take place in the same environment. The first objective function
#' in a sum of functions generates a new environment.
#' @return Object of class `objfn`.
#' @seealso [normL2], [constraintL2], [datapointL2]
#' @aliases sumobjfn
#' @example inst/examples/objective.R
#' @export
"+.objfn" <- function(x1, x2) {

  if (is.null(x1)) return(x2)

  conditions.x1 <- attr(x1, "conditions")
  conditions.x2 <- attr(x2, "conditions")
  conditions12 <- union(conditions.x1, conditions.x2)

  parameters.x1 <- attr(x1, "parameters")
  parameters.x2 <- attr(x2, "parameters")
  parameters12 <- union(parameters.x1, parameters.x2)

  modelname.x1 <- attr(x1, "modelname")
  modelname.x2 <- attr(x2, "modelname")
  modelname12 <- union(modelname.x1, modelname.x2)


  # objfn + objfn
  if (inherits(x1, "objfn") & inherits(x2, "objfn")) {

    outfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE,
                      conditions = conditions12, env = NULL,
                      cores = getOption("dMod.cores", 1L)) {

      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("pars"))]
      pars <- arglist[[1]]

      # 1. If conditions.xi is null, always evaluate xi, but only once
      # 2. If not null, evaluate at intersection with conditions
      # 3. If not null & intersection is empty, don't evaluate xi at all
      v1 <- v2 <- NULL
      if (is.null(conditions.x1)) {
        v1 <- x1(pars = pars, fixed = fixed, deriv = deriv, deriv2 = deriv2, conditions = conditions.x1, env = env, cores = cores)
      } else if (any(conditions %in% conditions.x1)) {
        v1 <- x1(pars = pars, fixed = fixed, deriv = deriv, deriv2 = deriv2, conditions = intersect(conditions, conditions.x1), env = env, cores = cores)
      }

      if (is.null(conditions.x2)) {
        v2 <- x2(pars = pars, fixed = fixed, deriv = deriv, deriv2 = deriv2, conditions = conditions.x2, env = env, cores = cores)
      } else if (any(conditions %in% conditions.x2)) {
        v2 <- x2(pars = pars, fixed = fixed, deriv = deriv, deriv2 = deriv2, conditions = intersect(conditions, conditions.x2), env = attr(v1, "env"), cores = cores)
      }

      out <- v1 + v2
      attr(out, "env") <- attr(v1, "env")
      return(out)
    }

    class(outfn) <- c("objfn", "fn")
    attr(outfn, "conditions") <- conditions12
    attr(outfn, "parameters") <- parameters12
    attr(outfn, "modelname") <- modelname12
    # Propagate the reconstruction handles so a composed objective exposes its
    # model pieces regardless of term order or nesting. Coalesce from either
    # operand.
    for (.a in c("prdfn", "data", "errfn", "timesD")) {
      .v <- attr(x1, .a, exact = TRUE)
      if (is.null(.v)) .v <- attr(x2, .a, exact = TRUE)
      if (!is.null(.v)) attr(outfn, .a) <- .v
    }
    # l2spec is CONCATENATED: every L2 term keeps its own data, prediction and
    # error model, which is what reml() needs from a split objective.
    attr(outfn, "l2spec") <- c(attr(x1, "l2spec", exact = TRUE),
                               attr(x2, "l2spec", exact = TRUE))
    return(outfn)

  }


}


#' Multiplication of objective functions with scalars
#'
#' @description The `\%.*\%` operator allows to multiply objects of class objlist or objfn with
#' a scalar.
#'
#' @param x1 object of class objfn or objlist.
#' @param x2 numeric of length one.
#' @return An objective function or objlist object.
#'
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
    outfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE,
                      conditions = conditions12, env = NULL,
                      cores = getOption("dMod.cores", 1L)) {

      arglist <- list(...)
      arglist <- arglist[match.fnargs(arglist, c("pars"))]
      pars <- arglist[[1]]

      v2 <- x2(pars = pars, fixed = fixed, deriv = deriv, deriv2 = deriv2,
               conditions = conditions, env = env, cores = cores)

      out <- x1 %.*% v2
      attr(out, "env") <- attr(v2, "env")
      return(out)
    }

    class(outfn) <- c("objfn", "fn")
    attr(outfn, "conditions") <- conditions12
    attr(outfn, "parameters") <- parameters12
    attr(outfn, "modelname") <- modelname12
    return(outfn)

  } else {

    x1*x2

  }

}


#' Direct sum of functions
#'
#' Used to add prediction function, parameter transformation functions or observation functions.
#'
#' @param x1 function of class `obsfn`, `prdfn` or `parfn`
#' @param x2 function of class `obsfn`, `prdfn` or `parfn`
#' @details Each prediction function is associated to a number of conditions. Adding functions
#' means merging or overwriting the set of conditions.
#' @return Object of the same class as `x1` and `x2` which returns results for the
#' union of conditions.
#' @aliases sumfn
#' @seealso [P], [Y], [Xs]
#' @example inst/examples/prediction.R
#' @export
"+.fn" <- function(x1, x2) {

  if (is.null(x1)) return(x2)

  k1 <- .fnKind(x1); k2 <- .fnKind(x2)
  if (is.null(k1) || is.null(k2) || !identical(k1, k2))
    stop("\"+.fn\": cannot add ", paste(class(x1), collapse = "/"), " and ",
         paste(class(x2), collapse = "/"), ".", call. = FALSE)

  own <- .mergeOwnership(x1, x2)

  st <- list2env(list(op = "+", kind = k1, parts = own$parts, owner = own$owner,
                      default_conditions = own$conditions), parent = emptyenv())
  outfn <- .fnWrap(st)

  attr(outfn, "mappings")    <- own$mappings
  attr(outfn, "parameters")  <- union(attr(x1, "parameters"), attr(x2, "parameters"))
  attr(outfn, "compileInfo") <- .mergeCompileInfo(attr(x1, "compileInfo"),
                                                  attr(x2, "compileInfo"))
  attr(outfn, "conditions")  <- own$conditions
  attr(outfn, "forcings")    <- .unionMappingAttr(own$mappings, "forcings")

  # Keep "composed" only when a composed operand went in, so summary() keeps
  # its detail branch for a sum of leaves and drops it for a sum of chains.
  cls <- c(k1, "fn")
  if (inherits(x1, "composed") || inherits(x2, "composed")) cls <- c(cls, "composed")
  class(outfn) <- cls

  outfn

}


#' Direct sum of datasets
#'
#' Used to merge datasets with overlapping conditions.
#'
#' @param data1 dataset of class `datalist`
#' @param data2 dataset of class `datalist`
#' @details Each data list contains data frames for a number of conditions.
#' The direct sum of datalist is meant as merging the two data lists and
#' returning the overarching datalist.
#' @return Object of class `datalist` for the
#' union of conditions.
#' @aliases sumdatalist
#' @example inst/examples/sumdatalist.R
#' @export
"+.datalist" <- function(data1, data2) {

  overlap <- names(data2)[names(data2) %in% names(data1)]
  if (length(overlap) > 0) {
    warning(paste("Condition", overlap, "existed and has been overwritten."))
    data1 <- data1[!names(data1) %in% names(data2)]
  }

  conditions <- union(names(data1), names(data2))
  data <- lapply(conditions, function(C) rbind(data1[[C]], data2[[C]]))
  names(data) <- conditions

  grid1 <- attr(data1, "condition.grid")
  grid2 <- attr(data2, "condition.grid")

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
  return(intersect(c1, c2))
}

#' Concatenation of functions
#'
#' Used to concatenate observation functions, prediction functions and parameter transformation functions.
#'
#' @param p1 function of class `obsfn`, `prdfn`, `parfn` or `idfn`
#' @param p2 function of class `obsfn`, `prdfn`, `parfn` or `idfn`
#' @return Object of the same class as `x1` and `x2`.
#' @aliases prodfn
#' @example inst/examples/prediction.R
#' @export
"*.fn" <- function(p1, p2) {

  # ============================================================
  # Global consistency check for condition handling
  #
  # Rules:
  # - A condition-unspecific function (conditions = NULL) may be
  #   combined with any other function.
  # - Two condition-specific functions must cover the same set
  #   of conditions.
  # - It is NOT allowed to combine a single-condition function
  #   with a multi-condition function.
  # ============================================================

  conditions.p1 <- attr(p1, "conditions")
  conditions.p2 <- attr(p2, "conditions")

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

  if (inherits(p1, "idfn")) return(p2)
  if (inherits(p2, "idfn")) return(p1)

  key  <- paste(.fnKind(p1), .fnKind(p2), sep = ".")
  spec <- if (length(key) == 1L) .prodSpec[[key]] else NULL
  if (is.null(spec))
    stop("\"*.fn\": no composition defined for ",
         paste(class(p1), collapse = "/"), " * ",
         paste(class(p2), collapse = "/"), ".", call. = FALSE)

  conditions.out <- out_conditions(conditions.p1, conditions.p2)

  st <- list2env(list(op = "*", kind = spec$out, p1 = p1, p2 = p2,
                      p1kind = .fnKind(p1), p2kind = .fnKind(p2),
                      handoff = spec$handoff, reduce = spec$reduce,
                      default_conditions = NULL), parent = emptyenv())
  outfn <- .fnWrap(st)

  attr(outfn, "conditions")  <- conditions.out
  attr(outfn, "parameters")  <- attr(p2, "parameters")
  attr(outfn, "compileInfo") <- .mergeCompileInfo(attr(p1, "compileInfo"),
                                                  attr(p2, "compileInfo"))

  if (identical(spec$out, "objfn")) {
    # An objfn carries no mappings; without these an objfn * parfn loses its
    # parameter set, its model name and the reconstruction handles.
    attr(outfn, "modelname") <- union(attr(p1, "modelname"), attr(p2, "modelname"))
    for (.a in c("data", "errfn", "timesD")) {
      .v <- attr(p1, .a, exact = TRUE)
      if (!is.null(.v)) attr(outfn, .a) <- .v
    }
    # The reconstructed prediction has to live in the outer coordinates.
    .prd <- attr(p1, "prdfn", exact = TRUE)
    if (!is.null(.prd))
      attr(outfn, "prdfn") <- tryCatch(.prd * p2, error = function(e) .prd)
    .l2 <- attr(p1, "l2spec", exact = TRUE)
    if (!is.null(.l2))
      attr(outfn, "l2spec") <- lapply(.l2, function(tm) {
        tm$prdfn <- tryCatch(tm$prdfn * p2, error = function(e) tm$prdfn)
        tm
      })
  } else {
    attr(outfn, "mappings") <- .composeMappings(st, p1, p2, conditions.out, spec$out)
  }

  class(outfn) <- c(spec$out, "fn", "composed")

  outfn

}
