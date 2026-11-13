## ---- Native chain ---------------------------------------------------------------------------
## normL2(data, g * x * p, errmodel = e) evaluated per condition in C++: the transformations by
## their own compiled functions, the batch solve of the model, the observables, the error model
## and the data term (src/chain_kernel.cpp). A plan per signature holds the indices; whatever it
## does not cover takes the R path, as does everything under options(dMod2.evaluation = "R").

## Whether objectives, FOCEI and the quadrature evaluate in C++: getOption("dMod2.evaluation"),
## "C++" (the default) or "R".
.evalCpp <- function() {
  ev <- getOption("dMod2.evaluation", "C++")
  if (!identical(ev, "C++") && !identical(ev, "R"))
    stop('options(dMod2.evaluation = ) is "C++" or "R".', call. = FALSE)
  identical(ev, "C++")
}

## A stage of the transformation for the parfn kernel `k` on the names `prev` before it:
## list(kind, fn, vjp, inNames, fnout, through, guessNames, state, order), a compiled function
## (kind 0) or a steady state (kind 1) whose native `state` solves the roots, with derivatives of
## order `order`. With `vjp` the vector-Jacobian product of the sweeps (`"_vjp_ad"`) or of the
## adjoint (`"_vjp"`). NULL where the kernel cannot take part.
.chainStage <- function(k, prev, cond, deriv2, order, vjp = NULL) {
  addr <- function(name) tryCatch(getNativeSymbolInfo(name)$address, error = function(e) NULL)
  if (isTRUE(attr(k, "pexpl", exact = TRUE))) {
    st <- environment(k)$st
    if (!isTRUE(st$use_ad) || !isTRUE(st$emit_d1) || (deriv2 && !isTRUE(st$emit_d2))) return(NULL)
    if (anyNA(match(st$parameters, prev))) return(NULL)
    fn <- addr(if (deriv2) st$ad2_symbol else st$ad_symbol)
    v <- if (!is.null(vjp)) addr(paste0(st$modelname, vjp))
    if (is.null(fn) || (!is.null(vjp) && is.null(v))) return(NULL)
    through <- if (isTRUE(environment(k)$controls$attachInput)) setdiff(prev, st$outnames) else character(0)
    return(list(kind = 0L, fn = fn, vjp = v, inNames = st$parameters, fnout = st$outnames, through = through,
                guessNames = character(0), state = NULL, order = 0L))
  }
  impl <- attr(k, "chainImpl", exact = TRUE)
  if (!is.function(impl)) return(NULL)
  im <- impl(cond)
  # a resting state absent from `prev` enters as 0 and leaves as 0
  rest0 <- intersect(setdiff(im$parameters, prev), im$resting)
  if ((order >= 1L && !im$d1) || (order >= 2L && !im$d2) ||
      anyNA(match(setdiff(im$parameters, rest0), prev))) return(NULL)
  list(kind = 1L, fn = NULL, vjp = NULL, inNames = im$parameters, fnout = c(im$dependent, im$zero, rest0),
       through = setdiff(prev, c(im$dependent, im$zero)), guessNames = im$dependent, state = im$state,
       order = as.integer(order))
}

## The stages of the parfn kernels `kernels`, in the order they apply, from the outer names
## `outer`: their indices into the vector before each, with the passed entries reduced to what a
## later stage or the names `needed` read; `o` the outer names the first one reads, `final` the
## names after the last. NULL where a kernel cannot take part.
.chainStages <- function(kernels, outer, needed, cond, deriv2, order, vjp = NULL) {
  st <- vector("list", length(kernels))
  prev <- outer
  for (s in seq_along(kernels)) {
    d <- .chainStage(kernels[[s]], prev, cond, deriv2, order, vjp)
    if (is.null(d)) return(NULL)
    st[[s]] <- d
    prev <- c(d$fnout, d$through)
  }
  need <- needed
  for (s in rev(seq_along(st))) {
    st[[s]]$through <- st[[s]]$through[st[[s]]$through %in% need]
    need <- unique(c(st[[s]]$inNames, st[[s]]$guessNames, st[[s]]$through))
  }
  o <- intersect(need, outer)
  prev <- o
  for (s in seq_along(st)) {
    d <- st[[s]]
    st[[s]] <- c(d, list(`in` = match(d$inNames, prev, 0L) - 1L, pass = match(d$through, prev) - 1L,
                         guess = match(d$guessNames, prev, 0L) - 1L, nfo = length(d$fnout)))
    prev <- c(d$fnout, d$through)
  }
  list(stages = st, o = o, final = prev)
}

## The plan for the objective of `spec` (an entry of attr(obj, "l2spec")) on `conds`: per
## condition the free outer names `free[[i]]`, the fixed ones `fixedNames` behind them, and the
## positions of the condition directions among `global`. With `times`, a list per condition, the plan of the
## prediction and error model on these times, without the data. With `tangents` the plan of the
## sweeps of forward over reverse instead: the solve on `times` along that many tangents of the
## free names, the data rows, and the backward solve with `seeds` seeds. With `adjoint` the plan of
## the adjoint gradient: values forward keeping the checkpoints, the backward solve of first order
## and the raw `_vjp` entries. NULL where a piece cannot take part.
.chainPlan <- function(spec, conds, free, fixedNames, global, deriv, deriv2 = FALSE, times = NULL,
                       tangents = 0L, seeds = 0L, adjoint = FALSE) {
  prd <- spec$prdfn
  err <- spec$errfn
  data <- spec$data
  if (is.null(.fnNode(prd)) || !length(conds)) return(NULL)
  if (is.null(global)) global <- character(0)
  addr <- function(name) tryCatch(getNativeSymbolInfo(name)$address, error = function(e) NULL)
  ad <- if (deriv2) "_eval_ad2" else "_eval_ad"
  sweep <- tangents > 0L
  if (sweep && (!deriv || deriv2 || seeds < 1L)) return(NULL)
  if (adjoint && (deriv || deriv2 || sweep)) return(NULL)
  frames <- !is.null(times) && !sweep
  vnames <- paste0("v", seq_len(tangents))
  # the derivatives a steady state gives with the root, and the vector-Jacobian products
  order <- if (deriv2 || sweep) 2L else if (deriv || adjoint) 1L else 0L
  vjp <- if (adjoint) "_vjp" else if (sweep) "_vjp_ad"
  timesOf <- function(i) if (!is.null(times)) times[[i]] else if (is.list(spec$timesD)) spec$timesD[[conds[i]]] else spec$timesD
  xk <- NULL
  leaves <- vector("list", length(conds))
  for (i in seq_along(conds)) {
    lv <- .fnLeaves(prd, conds[i])
    kinds <- vapply(lv, `[[`, "", "kind")
    if (length(kinds) < 2L || !identical(kinds[1:2], c("obsfn", "prdfn")) || any(kinds[-(1:2)] != "parfn"))
      return(NULL)
    if (is.null(xk)) xk <- lv[[2L]]$kernel
    else if (!identical(xk, lv[[2L]]$kernel)) return(NULL)
    gi <- attr(lv[[1L]]$kernel, "chainObs", exact = TRUE)
    gi <- if (is.function(gi)) gi()
    if (is.null(gi) || (deriv2 && !gi$d2)) return(NULL)
    ei <- NULL
    if (!is.null(err)) {
      le <- .fnLeaves(err, conds[i])
      if (length(le) > 1L) return(NULL)
      if (length(le)) {
        f <- attr(le[[1L]]$kernel, "chainObs", exact = TRUE)
        ei <- if (is.function(f)) f()
        if (is.null(ei) || (deriv2 && !ei$d2)) return(NULL)
      }
      if (frames && is.null(ei)) return(NULL)
    }
    needed <- c(attr(xk, "parameters"), gi$st$parameters, ei$st$parameters)
    sg <- .chainStages(rev(lapply(lv[-(1:2)], `[[`, "kernel")), c(free[[i]], fixedNames), needed, conds[i],
                       deriv2, order, vjp)
    if (is.null(sg)) return(NULL)
    leaves[[i]] <- list(sg = sg, gi = gi, ei = ei, kernels = lapply(lv[-(1:2)], `[[`, "kernel"),
                        dirs = if (sweep) vnames else if (adjoint) character(0) else sg$o[sg$o %in% free[[i]]])
  }
  solve <- attr(xk, "chainSolve", exact = TRUE)
  if (!is.function(solve)) return(NULL)
  tl <- lapply(seq_along(conds), timesOf)
  dl <- lapply(leaves, `[[`, "dirs")
  xs <- tryCatch(solve(tl, dl, deriv, deriv2, 0L, adjoint), error = function(e) NULL)
  if (is.null(xs)) return(NULL)
  xr <- if (sweep) tryCatch(solve(tl, dl, TRUE, FALSE, seeds), error = function(e) NULL)
        else if (adjoint) tryCatch(solve(tl, dl, FALSE, FALSE, max(1L, seeds), TRUE), error = function(e) NULL)
  if ((sweep || adjoint) && (is.null(xr) || !identical(xr$handle$required, xs$handle$required))) return(NULL)

  cond <- vector("list", length(conds))
  for (i in seq_along(conds)) {
    l <- leaves[[i]]
    sg <- l$sg
    gst <- l$gi$st; est <- l$ei$st
    outer <- c(free[[i]], fixedNames)
    oin <- match(sg$o, outer)
    xparm <- match(xs$handle$required, sg$final)
    # an observable reading `time` gets it from the solve, index 0 (-1 after the shift)
    yvar <- match(gst$innames, xs$variables)
    yvar[is.na(yvar) & gst$innames == "time"] <- 0L
    ypar <- match(gst$parameters, sg$final)
    if (anyNA(xparm) || anyNA(yvar) || anyNA(ypar)) return(NULL)
    yfn <- addr(paste0(gst$modelname, ad))
    if (is.null(yfn)) return(NULL)
    entry <- list(
      p_oidx = oin - 1L, p_odir = match(sg$o, l$dirs, 0L) - 1L,
      p_ovrow = ifelse(oin <= length(free[[i]]), oin - 1L, -1L),
      stages = lapply(sg$stages, function(d)
        d[c("kind", "fn", "vjp", "in", "pass", "guess", "nfo", "state", "order")]),
      nd = length(l$dirs), n_pout = length(sg$final), x_args = xs$handle$preps[[i]]$args, x_parm = xparm - 1L,
      y_fn = yfn, y_var = yvar - 1L, y_par = ypar - 1L, n_yout = length(gst$outnames), y_names = gst$outnames,
      e_fn = NULL, e_var = integer(0), e_par = integer(0), n_eout = 0L, e_vjp = NULL, e_vjp1 = NULL,
      e_names = if (!is.null(est)) est$outnames,
      dir_names = l$dirs, dir_global = match(l$dirs, global, 0L) - 1L)
    if (sweep || adjoint) {
      entry$xr_args <- xr$handle$preps[[i]]$args
      entry[[if (adjoint) "y_vjp1" else "y_vjp"]] <- addr(paste0(gst$modelname, vjp))
      if (is.null(entry[[if (adjoint) "y_vjp1" else "y_vjp"]])) return(NULL)
    }
    if (!is.null(est)) {
      evar <- match(est$innames, gst$outnames)
      evar[is.na(evar) & est$innames == "time"] <- 0L
      epar <- match(est$parameters, sg$final)
      efn <- addr(paste0(est$modelname, ad))
      evjp <- if (sweep || adjoint) addr(paste0(est$modelname, vjp))
      if (!anyNA(evar) && !anyNA(epar) && !is.null(efn) && (!(sweep || adjoint) || !is.null(evjp)))
        entry[c("e_fn", "e_var", "e_par", "n_eout", if (adjoint) "e_vjp1" else "e_vjp")] <-
          list(efn, evar - 1L, epar - 1L, length(est$outnames), evjp)
      else if (frames) return(NULL)
    }
    if (frames) {
      cond[[i]] <- entry
      next
    }
    d <- data[[conds[i]]]
    tIdx <- match(d$time, timesOf(i))
    oIdx <- match(as.character(d$name), gst$outnames)
    sig <- if (!is.null(d$sigma)) d$sigma else rep(NA_real_, nrow(d))
    na <- is.na(sig)
    eIdx <- rep(NA_integer_, nrow(d))
    if (any(na)) {
      if (is.null(entry$e_fn)) return(NULL)
      eIdx <- match(as.character(d$name), est$outnames)
      if (anyNA(eIdx[na])) return(NULL)
    }
    if (anyNA(tIdx) || anyNA(oIdx)) return(NULL)
    lloq <- if (!is.null(d$lloq)) d$lloq else rep(-Inf, nrow(d))
    val <- pmax(d$value, lloq)
    cond[[i]] <- c(entry, list(
      t_data = as.numeric(d$time), o_idx = oIdx - 1L, e_idx = ifelse(na, eIdx, 1L) - 1L, sigma_na = as.integer(na),
      sigma_fixed = ifelse(na, 0, sig), y_data = as.numeric(val), lloq = as.numeric(lloq),
      bloq = as.integer(val <= lloq)))
  }
  list(cond = cond, n_global = length(global), sym = xs$handle$sym, handle = xs$handle, run = xs$run,
       rsym = xr$handle$sym, rhandle = xr$handle, rrun = xr$run,
       kernels = c(list(xk), unlist(lapply(leaves, `[[`, "kernels"), recursive = FALSE)),
       controls = c(list(environment(xk)$controls),
                    lapply(unlist(lapply(leaves, `[[`, "kernels"), recursive = FALSE), function(k) environment(k)$controls)),
       observables = lapply(leaves, function(l) l$gi), onFailure = xs$onFailure)
}

## Whether `plan` still holds: the shared objects loaded and the controls of the leaves as they
## were.
.chainLive <- function(plan) {
  live <- function(h) is.null(h$sym) || .batchHandleLive(h)
  live(plan$handle) && (is.null(plan$rhandle) || live(plan$rhandle)) &&
    identical(lapply(plan$kernels, function(k) environment(k)$controls), plan$controls)
}

## The transformations of chain_params(), the steady states solved on the way; `roots` in the
## result for chain_pull(). NULL where a function fails.
.chainParams <- function(plan, values, deriv, deriv2, cores, tangents = NULL) {
  pre <- chain_params(plan, values, deriv, deriv2, as.integer(cores), tangents)
  for (w in pre$warnings) warning(w, call. = FALSE)
  pre
}

## Frees the checkpoints of reverse solves at once; a cppDE without releaseStore() leaves them to
## the garbage collector.
.releaseStores <- function(stores) {
  release <- get0("releaseStore", envir = asNamespace("cppDE"), inherits = FALSE)
  if (!is.null(release) && length(stores)) release(stores)
  invisible(NULL)
}

## The forward solve of `plan` on the arguments `args`, or with `back` the backward solve; NULL
## where it fails.
.chainSolveRun <- function(plan, args, cores, back = FALSE) {
  run <- if (back) plan$rrun else plan$run
  sym <- if (back) plan$rsym else plan$sym
  cores <- as.integer(max(1L, cores))
  tryCatch(if (is.function(run)) run(args, cores) else do.call(.Call, list(sym$name, args, cores, PACKAGE = sym$dll)),
           error = function(e) NULL)
}

## The data term of `plan` at `values`, the outer values per condition or one vector for all, in
## the order of the plan; NULL where the chain fails, for the R path to report it.
.chainEval <- function(plan, values, deriv, deriv2, hessian, local, cores, bloq, df) {
  pre <- .chainParams(plan, values, deriv, deriv2, cores)
  raw <- if (!is.null(pre)) .chainSolveRun(plan, pre$args, cores)
  if (is.null(raw)) return(NULL)
  chain_data(plan, pre, raw, deriv, deriv2, hessian, local, as.integer(cores), bloq, if (is.finite(df)) df else 0,
             FALSE)
}


## The data term of `ctx` per request in C++: request j on condition `conds[j]`, which may repeat,
## at the free values `parsList[[j]]` and fixed values `fixedList[[j]]` (NULL for none, the same
## names in every request), with value, chi2 and, with `deriv`, gradient and Hessian in the names of
## the free values; NULL where the chain does not cover the model or fails, for the R path.
.chainLocal <- function(ctx, conds, parsList, fixedList, deriv, deriv2, cores, bloq, df) {
  if (is.null(ctx) || !.evalCpp() || !length(conds)) return(NULL)
  fixedList <- fixedList %||% rep(list(NULL), length(conds))
  plan <- .chainLookup(ctx, conds, lapply(parsList, names), names(fixedList[[1L]]), NULL, deriv, deriv2)
  if (is.null(plan)) return(NULL)
  values <- Map(function(p, f) c(unclass(p), unclass(f)), parsList, fixedList)
  .chainEval(plan, values, deriv, deriv2, deriv, TRUE, cores, bloq, df)$local
}


## The value and adjoint gradient of the data term of `ctx` on `conditions` at `pars` and `fixed`:
## values forward, the derivatives of the value in the predictions and sigmas pulled back in one
## sweep; NULL where the chain does not cover the model or fails, for the R path to take it.
.chainAdjoint <- function(ctx, conditions, pars, fixed, cores, bloq, df) {
  if (is.null(ctx) || !.evalCpp()) return(NULL)
  nm <- names(pars)
  # under refine a second seed for chi2, as in .normL2_reverse()
  correct <- .refines(ctx$spec$prdfn, conditions)
  S <- if (correct) 2L else 1L
  plan <- .chainLookup(ctx, conditions, rep(list(nm), length(conditions)), names(fixed), nm, FALSE, FALSE,
                       seeds = S, adjoint = TRUE)
  if (is.null(plan)) return(NULL)
  cores <- as.integer(max(1L, cores))
  values <- list(c(as.numeric(pars), as.numeric(fixed)))
  pre <- .chainParams(plan, values, FALSE, FALSE, cores)
  raw <- if (!is.null(pre)) .chainSolveRun(plan, pre$args, cores)
  # the checkpoints serve this one backward solve
  on.exit(if (is.list(raw)) .releaseStores(lapply(raw, function(r) if (is.list(r)) r$store)), add = TRUE)
  dat <- if (!is.null(raw))
    chain_data(plan, pre, raw, FALSE, FALSE, FALSE, FALSE, cores, bloq, if (is.finite(df)) df else 0, TRUE)
  back <- if (!is.null(dat))
    chain_back(plan, pre, raw, cbind(dat$seed_pred, if (correct) dat$seed_chi2_pred),
               cbind(dat$seed_sigma, if (correct) dat$seed_chi2_sigma), S, cores, TRUE)
  rraw <- if (!is.null(back)) .chainSolveRun(plan, back$args, cores, back = TRUE)
  pull <- if (!is.null(rraw))
    chain_pull(plan, values, matrix(0, length(nm), 0), pre$roots, back, rraw, S, cores, TRUE)
  if (is.null(pull)) return(NULL)
  gradient <- setNames(numeric(length(nm)), nm)
  for (u in pull) gradient[u$rows + 1L] <- gradient[u$rows + 1L] + u$ct[, 1L]
  if (!correct) return(list(value = dat$value, chi2 = dat$chi2, gradient = gradient))
  corr <- c(0, 0)
  for (r in rraw) if (!is.null(r$correction)) corr <- corr + as.numeric(r$correction)
  list(value = dat$value + corr[1L], chi2 = dat$chi2 + corr[2L], gradient = gradient,
       correction = c(value = corr[1L], chi2 = corr[2L]))
}


## The plans of a data term `spec` (prdfn, errfn, data, timesD, as in attr(obj, "l2spec")): the
## entries per condition, kind and derivative order, and the plans they come from by number.
.chainContext <- function(spec) {
  ctx <- new.env(parent = emptyenv())
  ctx$spec <- spec
  ctx$entries <- new.env(parent = emptyenv())
  ctx$plans <- new.env(parent = emptyenv())
  ctx$n <- 0L
  ctx$last <- NULL
  ctx
}

## The plan of `ctx` for `conditions`, built where missing; NULL where the chain does not cover
## them. `times`, a list per condition, asks for the frames on these times instead of the data,
## with `tangents` and `seeds` for the sweeps (.chainPlan()).
.chainLookup <- function(ctx, conditions, free, fixedNames, global, deriv, deriv2, times = NULL,
                         retry = TRUE, tangents = 0L, seeds = 0L, adjoint = FALSE) {
  # a repeated request: the plan of the last one while the shared objects of that plan are loaded
  req <- list(conditions, free, fixedNames, global, deriv, deriv2, times, tangents, seeds, adjoint)
  last <- ctx$last
  if (!is.null(last) && identical(last$req, req) && all(vapply(last$plans, .chainLive, TRUE)))
    return(last$out)
  frames <- !is.null(times)
  # per condition and kind the plans for the names asked for, a few side by side, as the inner and
  # the outer evaluations of a fit ask for them; the newest first
  keys <- paste(conditions, deriv, deriv2, frames, tangents, seeds, adjoint)
  slots <- mget(keys, envir = ctx$entries, ifnotfound = list(NULL))
  fits <- function(e, i) identical(e$free, free[[i]]) && identical(e$fixed, fixedNames) &&
    identical(e$global, global) && (!frames || identical(e$times, times[[i]]))
  have <- lapply(seq_along(keys), function(i) {
    for (e in slots[[i]]) if (fits(e, i)) return(e)
    NULL
  })
  stale <- which(vapply(have, is.null, TRUE))
  if (length(stale)) {
    plan <- .chainPlan(ctx$spec, conditions[stale], free[stale], fixedNames, global, deriv, deriv2,
                       if (frames) times[stale], tangents, seeds, adjoint)
    id <- NULL
    if (!is.null(plan)) {
      ctx$n <- ctx$n + 1L
      id <- as.character(ctx$n)
      assign(id, plan, envir = ctx$plans)
    }
    for (j in seq_along(stale)) {
      i <- stale[j]
      e <- list(free = free[[i]], fixed = fixedNames, global = global, none = is.null(plan),
                times = if (frames) times[[i]], entry = if (!is.null(plan)) plan$cond[[j]], plan = id)
      slot <- c(list(e), get0(keys[i], envir = ctx$entries, inherits = FALSE))
      assign(keys[i], slot[seq_len(min(length(slot), 4L))], envir = ctx$entries)
      have[[i]] <- e
    }
  }
  if (any(vapply(have, function(e) isTRUE(e$none), TRUE))) return(NULL)
  plans <- mget(unique(vapply(have, `[[`, "", "plan")), envir = ctx$plans)
  # a plan whose shared objects or controls changed is built again, once
  if (!all(vapply(plans, .chainLive, TRUE))) {
    rm(list = ls(ctx$entries), envir = ctx$entries)
    rm(list = ls(ctx$plans), envir = ctx$plans)
    ctx$last <- NULL
    return(if (retry) .chainLookup(ctx, conditions, free, fixedNames, global, deriv, deriv2, times, FALSE,
                                   tangents, seeds, adjoint))
  }
  out <- list(cond = lapply(have, `[[`, "entry"), n_global = length(global), sym = plans[[1L]]$sym,
              rsym = plans[[1L]]$rsym, run = plans[[1L]]$run, rrun = plans[[1L]]$rrun)
  ctx$last <- list(req = req, plans = plans, out = out)
  out
}

## The prediction and error model of `conds` on `times`, per condition list(pred, err) as matrices
## with the time first and the derivatives as attributes, from the outer values `parsList` (one
## vector for all or one per condition) and `fixed`; NULL where the chain does not cover them.
.chainFrames <- function(ctx, conds, times, parsList, fixed, deriv2, cores) {
  if (is.null(ctx) || !.evalCpp()) return(NULL)
  free <- if (length(parsList) == 1L) rep(list(names(parsList[[1L]])), length(conds)) else lapply(parsList, names)
  plan <- .chainLookup(ctx, conds, free, names(fixed), NULL, TRUE, deriv2, times)
  if (is.null(plan)) return(NULL)
  values <- lapply(parsList, function(p) c(unclass(p), unclass(fixed)))
  pre <- .chainParams(plan, values, TRUE, deriv2, cores)
  raw <- if (!is.null(pre)) .chainSolveRun(plan, pre$args, cores)
  if (is.null(raw)) return(NULL)
  out <- chain_frames(plan, pre, raw, deriv2, as.integer(cores))
  if (!is.null(out)) names(out) <- conds
  out
}
