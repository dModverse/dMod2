## Functions to generate parameter transformation ----

#' Parameter transformation
#'
#' Builds a [parfn] with [Pexpl()] or [Pimpl()].
#'
#' @param trafo An [eqnvec], named character, [eqnlist], or list thereof.
#' @param parameters Outer parameter names.
#' @param condition Condition label.
#' @param compile,modelname,verbose Forwarded to [cppDE::cppFUN()].
#' @param method `"explicit"`, `"implicit"`, or `NULL` for `"implicit"` if
#'   `trafo` is an [eqnlist] and `"explicit"` otherwise.
#' @param cores Number of cores. `NULL` uses [detectFreeCores()]; 1 on Windows.
#' @param deriv,deriv2 Attach first and second order sensitivities. `deriv2`
#'   requires `deriv = TRUE`.
#' @param ... Forwarded to the backend, e.g. `outdir`.
#'
#' @return A [parfn].
#' @seealso [Pexpl()], [Pimpl()]
#' @export
P <- function(trafo = NULL, parameters = NULL, condition = NULL,
              compile = FALSE, modelname = NULL, method = NULL,
              cores = NULL, verbose = FALSE,
              deriv = TRUE, deriv2 = FALSE, ...) {

  if (is.null(trafo)) return()
  if (isTRUE(deriv2) && !isTRUE(deriv))
    stop("P(deriv2 = TRUE) requires deriv = TRUE.", call. = FALSE)

  if (!is.list(trafo) || inherits(trafo, "eqnlist") || inherits(trafo, "eqnvec")) {
    trafo_list <- list(trafo); names(trafo_list) <- condition
  } else trafo_list <- trafo

  ## detectFreeCores() has side effects (warnings, SSH on remotes); call once.
  cores <- if (Sys.info()[['sysname']] == "Windows") 1L
           else if (is.null(cores)) detectFreeCores()
           else min(detectFreeCores(), cores)

  ## Codegen-only inside mclapply: dyn.load in a forked worker is lost when
  ## the fork exits, so the actual compile happens in the parent below.
  result <- Reduce("+", mclapply(seq_along(trafo_list), function(i) {
    tr   <- trafo_list[[i]]
    cond <- names(trafo_list[i])
    m <- if (!is.null(method)) match.arg(method, c("explicit", "implicit"))
         else if (inherits(tr, "eqnlist")) "implicit" else "explicit"
    switch(m,
      explicit = Pexpl(as.eqnvec(tr), parameters = parameters, condition = cond,
                       compile = FALSE, modelname = modelname, verbose = verbose,
                       deriv = deriv, deriv2 = deriv2, ...),
      implicit = Pimpl(trafo = tr, parameters = parameters, condition = cond,
                       compile = FALSE, modelname = modelname, verbose = verbose,
                       deriv = deriv, deriv2 = deriv2, ...))
  }, mc.cores = cores))

  if (compile) compile(result, cores = cores, output = modelname, verbose = verbose)
  result
}


## Per-condition warm-start caches of a Pimpl. A condition-less Pimpl composed
## with a condition-specific Pexpl is called once per condition; each condition
## keeps its own cache, keyed by the condition parfn() passes on. A NULL or
## empty key uses a shared slot.
.warmstart_registry <- function() {
  caches <- new.env(parent = emptyenv())
  get_cache <- function(key) {
    k <- if (is.null(key) || !nzchar(key)) "__default__" else key
    cc <- get0(k, envir = caches, inherits = FALSE)
    if (is.null(cc)) { cc <- new.env(parent = emptyenv()); assign(k, cc, envir = caches) }
    cc
  }
  reset <- function() {
    nms <- ls(caches, all.names = TRUE)
    for (k in nms) {
      cc <- get(k, envir = caches, inherits = FALSE)
      inner <- ls(cc, all.names = TRUE)
      if (length(inner)) rm(list = inner, envir = cc)
    }
    invisible(nms)
  }
  list(get = get_cache, reset = reset, caches = caches)
}



## Body of the parameter transformation built by Pexpl(). Package level, so
## one parfn per condition carries the state, not another copy of this code.

.Pexpl_p2p <- function(st, pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE, .ad_out = NULL,
                       attach.input = FALSE) {

  if (deriv2 && !st$emit_d2)
    stop("Pexpl was built with deriv2 = FALSE; rebuild with deriv2 = TRUE.", call. = FALSE)
  if (!st$emit_d1) deriv <- FALSE
  if (deriv2 && !deriv) deriv <- TRUE

  p <- c(pars, fixed)
  ## A missing or unloaded build is reported by the AD entry itself.
  ad_ok <- st$use_ad && !is.null(st$evaluate)

  Jac <- NULL; Hess <- NULL

  if (ad_ok && deriv) {
    ## Dual-mode entry reads `params` and `dP` positionally against the
    ## codegen order, so reorder both to `st$parameters`.
    dP  <- attr(pars, "deriv")
    dP2 <- if (deriv2) attr(pars, "deriv2") else NULL
    if (is.null(dP)) {
      active <- setdiff(st$parameters, names(fixed))
      dP <- diag(length(active)); dimnames(dP) <- list(active, active)
    }
    out <- if (!is.null(.ad_out)) .ad_out else
      st$evaluate(NULL, p[st$parameters], tangentX = NULL, tangentP = dP,
                  hessianX = NULL, hessianP = dP2,
                  deriv2 = deriv2, attach.input = attach.input,
                  fixed = intersect(names(fixed), st$parameters))
    pinnerVal <- out$y[1, ]
    tg <- out$tangent
    if (!is.null(tg))
      Jac <- matrix(tg, dim(tg)[2], dim(tg)[3],
                    dimnames = list(dimnames(tg)[[2]], dimnames(tg)[[3]]))
    if (deriv2 && !is.null(out$hessian))
      Hess <- array(out$hessian, dim(out$hessian)[2:4],
                    dimnames = dimnames(out$hessian)[2:4])
  } else {
    ## Values only (reverse-only build). The inputs are passed through below,
    ## as on the derivative path, so both return the same names.
    pinnerVal <- st$fun(NULL, p, fixed = names(fixed))[, ]
  }

  if (any(is.nan(pinnerVal)))
    stop("Inner parameter(s) evaluate to NaN:\n\t",
         paste(names(pinnerVal)[is.nan(pinnerVal)], collapse = "\n\t"),
         ".\nLikely cause: division by zero or missing inputs.", call. = FALSE)

  ## attach.input appends every input the transformation does not map, a fixed
  ## one included, as Pimpl does.
  through <- if (attach.input) setdiff(names(p), names(pinnerVal)) else character(0)
  val <- if (length(through)) c(pinnerVal, .subset(p, through)) else pinnerVal

  ## A build without forward derivatives has none to give. Rows for the inputs
  ## alone would mark every output as fixed.
  if (!deriv || is.null(Jac))
    return(as.parvec(val, deriv = FALSE, deriv2 = FALSE))

  Jac  <- Jac[rowSums(Jac != 0) > 0, , drop = FALSE]
  Hess <- if (deriv2 && !is.null(Hess)) Hess[rownames(Jac), , , drop = FALSE]

  ## An input without a derivative row counts as fixed downstream, and the
  ## forward gradient along it would be zero while the reverse path, which
  ## hands each input its cotangent back, has it right. So every input that
  ## varies keeps its row. What the caller fixed has none and stays fixed.
  moving <- setdiff(through, names(fixed))
  if (length(moving)) {
    d <- .Pexpl_through(pars, moving, colnames(Jac), second = !is.null(Hess))
    Jac <- .Pexpl_widen(Jac, d$theta)
    Jac <- rbind(Jac, d$deriv)
    if (!is.null(Hess))
      Hess <- .Pexpl_bind3(.Pexpl_widen(Hess, d$theta), d$deriv2)
  }
  as.parvec(val, deriv = Jac, deriv2 = if (is.null(Hess)) FALSE else Hess)
}


# The derivatives of the inputs Pexpl passes through, in the basis `theta` of
# the transformation's own Jacobian. An input carries its derivatives in; where
# none come in, it is a parameter of the chain and its own direction, which the
# basis gains if the transformation does not read it. An input that comes in
# without a row was fixed further up and keeps no row.
.Pexpl_through <- function(pars, moving, theta, second = FALSE) {
  dP <- attr(pars, "deriv")
  if (is.null(dP)) {
    theta <- union(theta, moving)
    D <- matrix(0, length(moving), length(theta), dimnames = list(moving, theta))
    D[cbind(moving, moving)] <- 1
    D2 <- if (second)
      array(0, c(length(moving), length(theta), length(theta)),
            dimnames = list(moving, theta, theta))
    return(list(theta = theta, deriv = D, deriv2 = D2))
  }
  rows <- intersect(moving, rownames(dP))
  D <- dP[rows, theta, drop = FALSE]
  D2 <- NULL
  if (second) {
    D2 <- array(0, c(length(rows), length(theta), length(theta)),
                dimnames = list(rows, theta, theta))
    dP2 <- attr(pars, "deriv2")
    have <- intersect(rows, dimnames(dP2)[[1L]])
    if (length(have)) D2[have, , ] <- dP2[have, theta, theta, drop = FALSE]
  }
  list(theta = theta, deriv = D, deriv2 = D2)
}

# A Jacobian [p, theta] or a Hessian [p, theta, theta] zero-padded to a larger
# basis, the columns it already has kept in place.
.Pexpl_widen <- function(x, theta) {
  have <- dimnames(x)[[2L]]
  if (identical(have, theta)) return(x)
  if (length(dim(x)) == 2L) {
    out <- matrix(0, nrow(x), length(theta), dimnames = list(rownames(x), theta))
    out[, have] <- x
  } else {
    out <- array(0, c(dim(x)[1L], length(theta), length(theta)),
                 dimnames = list(dimnames(x)[[1L]], theta, theta))
    out[, have, have] <- x
  }
  out
}

# Stack two Hessians [p, theta, theta] of one basis along their rows.
.Pexpl_bind3 <- function(a, b) {
  na <- dim(a)[1L]; nb <- dim(b)[1L]; nt <- dim(a)[2L]
  out <- array(0, c(na + nb, nt, nt),
               dimnames = list(c(dimnames(a)[[1L]], dimnames(b)[[1L]]),
                               dimnames(a)[[2L]], dimnames(a)[[3L]]))
  if (na) out[seq_len(na), , ] <- a
  if (nb) out[na + seq_len(nb), , ] <- b
  out
}


# One evaluateBatch over all conditions; the value-only path still loops.
.Pexpl_batch <- function(st, parsList, fixedList, deriv, deriv2, cores,
                         attach.input = FALSE) {
  n <- length(parsList)
  loop <- function() lapply(seq_len(n), function(i)
    .Pexpl_p2p(st, parsList[[i]], fixedList[[i]], deriv, deriv2,
               attach.input = attach.input))

  eb <- st$evaluateBatch
  ad_ok <- st$use_ad && !is.null(st$evaluate) && is.loaded(st$ad_symbol)
  if (is.null(eb) || !ad_ok || !deriv || !st$emit_d1) return(loop())
  if (deriv2 && !isTRUE(st$emit_d2)) return(loop())

  sets <- lapply(seq_len(n), function(i) {
    pars <- parsList[[i]]; fixed <- fixedList[[i]]
    dP <- attr(pars, "deriv")
    if (is.null(dP)) {
      active <- setdiff(st$parameters, names(fixed))
      dP <- diag(length(active)); dimnames(dP) <- list(active, active)
    }
    list(vars = NULL, params = c(pars, fixed)[st$parameters],
         tangentX = NULL, tangentP = dP, hessianX = NULL,
         hessianP = if (deriv2) attr(pars, "deriv2") else NULL,
         attach.input = attach.input,
         fixed = intersect(names(fixed), st$parameters))
  })
  ad <- eb(sets, cores = cores, deriv2 = deriv2)
  lapply(seq_len(n), function(i)
    .Pexpl_p2p(st, parsList[[i]], fixedList[[i]], deriv, deriv2,
               .ad_out = ad[[i]], attach.input = attach.input))
}


# The kernel with its batch and reverse entries. All three close over this
# frame, which holds `st` and `controls` alone, and read `controls` when they
# run: controls<- rebinds it here, so a change reaches every entry at once.
.Pexpl_wrap <- function(st, controls) {
  p2p <- function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE)
    .Pexpl_p2p(st, pars, fixed, deriv, deriv2,
               attach.input = controls$attach.input)
  attr(p2p, "vjpfn") <- function(pars, fixed = NULL, cotangent, condition = NULL)
    .Pexpl_vjp(st, pars, fixed, cotangent, condition,
               attach.input = controls$attach.input)
  attr(p2p, "batchfn") <- function(parsList, fixedList, deriv, deriv2,
                                   conditions, cores)
    .Pexpl_batch(st, parsList, fixedList, deriv, deriv2, cores,
                 attach.input = controls$attach.input)
  p2p
}

# w' Jac, where the forward path forms Jac %*% dP. The transformation is one
# evaluation with no variables and one observation, so the vjp is the same call
# the observation functions make, with the cotangent on the inner parameters.
#
# attach.input passes the outer parameters through untouched, so their cotangent
# adds to whatever the transformation itself puts on them.
.Pexpl_vjp <- function(st, pars, fixed = NULL, cotangent, condition = NULL,
                       attach.input = FALSE) {
  if (is.null(st$vjp))
    stop("Pexpl(): the reverse mode needs a vector-Jacobian product; rebuild ",
         "with derivMode = c(\"forward\", \"reverse\") and compile = TRUE.",
         call. = FALSE)
  w <- .asCtPars(cotangent)
  K <- .ctK(w)
  p <- c(pars, fixed)
  outnames <- st$outnames
  W <- matrix(0, 1L, length(outnames), dimnames = list(NULL, outnames))
  hit <- intersect(rownames(w), outnames)
  if (length(hit)) W[1L, hit] <- w[hit, 1L]

  if (K == 1L) {
    r <- st$vjp(NULL, p[st$parameters], W)
    u <- matrix(r$cotangentP[, 1L], ncol = 1L,
                dimnames = list(rownames(r$cotangentP), NULL))
  } else {
    # The node has no variables, so only the parameters carry tangents in, and
    # the cotangent brings its own. Both halves of d/dv (w' J) come back in one
    # pass; nothing here forms a Hessian.
    nd <- K - 1L
    V <- matrix(0, length(st$parameters), nd,
                dimnames = list(st$parameters, NULL))
    dp <- attr(pars, "deriv")
    if (!is.null(dp)) {
      take <- intersect(rownames(dp), st$parameters)
      if (length(take)) V[take, ] <- dp[take, seq_len(nd), drop = FALSE]
    }
    DW <- array(0, c(1L, length(outnames), 1L, nd))
    if (length(hit))
      DW[1L, match(hit, outnames), 1L, ] <- w[hit, -1L, drop = FALSE]
    r <- st$vjp(NULL, p[st$parameters], W, tangentP = V, curvature = DW)
    u <- cbind(r$cotangentP[, 1L, drop = FALSE],
               matrix(r$curvatureP[, 1L, ], ncol = nd))
    rownames(u) <- rownames(r$cotangentP)
  }
  wp <- .pickCotangent(u, names(pars))

  if (attach.input) {
    through <- setdiff(rownames(w), outnames)
    keep <- intersect(through, rownames(wp))
    if (length(keep))
      wp[keep, ] <- wp[keep, , drop = FALSE] + w[keep, , drop = FALSE]
  }
  wp
}

#' Parameter transformation (explicit, algebraic)
#'
#' Builds `p_inner = f(p_outer)` from symbolic expressions via
#' [cppDE::cppFUN], with derivatives by AD. The returned [parfn] attaches the
#' Jacobian and, optionally, the Hessian. It is evaluable only after
#' compilation.
#'
#' @param trafo Named character / [eqnvec]; names are inner parameters,
#'   values are expressions in the outer parameters.
#' @param parameters Outer parameters; defaults to `getSymbols(trafo)`.
#' @param attach.input Append the outer inputs the transformation does not map
#'   to the output, fixed ones included. An input that varies keeps its
#'   derivatives, the identity when it enters the chain here, and one passed in
#'   `fixed` stays fixed. Kept as a control of the returned function and read at
#'   every call, see [controls()].
#' @param condition Condition label.
#' @param compile,modelname,verbose Forwarded to [cppDE::cppFUN].
#' @param deriv,deriv2 Attach `attr(., "deriv")` `[p, theta]` and/or
#'   `attr(., "deriv2")` `[p, theta, theta]`. `deriv2` needs `deriv = TRUE`.
#' @param derivMode Which derivative products to build, any of `"forward"`
#'   (AD, default), `"reverse"` (the vector-Jacobian product the reverse sweep
#'   contracts against) and `"forward-reverse"` (its derivative along a tangent,
#'   for the reverse sweep with `deriv2 = TRUE`).
#' @param outdir Directory for the generated source and shared object,
#'   default the working directory.
#'
#' @return A [parfn].
#' @seealso [Pimpl], [P].
#' @importFrom cppDE cppFUN
#' @export
Pexpl <- function(trafo, parameters = NULL, attach.input = FALSE, condition = NULL,
                  compile = FALSE, modelname = NULL, verbose = FALSE,
                  deriv = TRUE, deriv2 = FALSE,
                  derivMode = "forward",
                  outdir = getwd()) {

  derivMode <- .matchDerivMode(derivMode, c("forward", "reverse", "forward-reverse"))
  emit_d1   <- isTRUE(deriv)
  emit_d2   <- isTRUE(deriv2)
  if (emit_d2 && !emit_d1)
    stop("Pexpl(deriv2 = TRUE) requires deriv = TRUE.", call. = FALSE)

  if (is.null(parameters)) {
    parameters <- getSymbols(trafo)
  } else {
    identity <- setNames(parameters[!(parameters %in% names(trafo))],
                         parameters[!(parameters %in% names(trafo))])
    trafo <- c(trafo, identity)
    parameters <- getSymbols(trafo)
  }

  if (is.null(modelname)) modelname <- "expl_parfn"
  if (!is.null(condition)) modelname <- paste(modelname, sanitizeConditions(condition), sep = "_")

  PEval <- suppressWarnings(cppDE::cppFUN(
    unclass(trafo), variables = NULL, parameters = parameters, fixed = NULL,
    compile = compile, modelname = modelname, outdir = outdir,
    verbose = verbose, convenient = FALSE, derivMode = derivMode,
    deriv = emit_d1, deriv2 = emit_d2))

  fun <- PEval$func; jac <- PEval$jac; hess <- PEval$hess; evaluate <- PEval$evaluate
  use_ad     <- "forward" %in% derivMode
  ad_symbol  <- paste0(modelname, "_eval_ad")
  ad2_symbol <- paste0(modelname, "_eval_ad2")

  ## The wrapper closes over `st` and `controls` alone, not over Pexpl's frame.
  st <- list2env(list(fun = fun, jac = jac, hess = hess, evaluate = evaluate,
                      evaluateBatch = PEval$evaluateBatch, vjp = PEval$vjp,
                      outnames = names(trafo),
                      parameters = parameters,
                      use_ad = use_ad, ad_symbol = ad_symbol,
                      ad2_symbol = ad2_symbol, emit_d1 = emit_d1,
                      emit_d2 = emit_d2), parent = emptyenv())
  p2p <- .Pexpl_wrap(st, controls = list(attach.input = attach.input))

  attr(p2p, "equations")   <- as.eqnvec(trafo)
  attr(p2p, "parameters")  <- parameters
  attr(p2p, "modelname")   <- modelname
  attr(p2p, "compileInfo") <- .collectCompileInfo(fun, jac, hess, evaluate)
  parfn(p2p, parameters, condition)
}


#' Conserved-quantity coefficient matrix
#'
#' @param totals Named list of conserved-quantity expressions (from
#'   [getTotals()]).
#' @param states State names.
#' @return Matrix with one row per total and one column per participating
#'   state, or `NULL` for no totals.
#' @keywords internal
.cq_matrix <- function(totals, states) {
  if (!length(totals)) return(NULL)
  sp <- unique(unlist(lapply(totals, function(e) intersect(getSymbols(e), states)),
                      use.names = FALSE))
  C <- matrix(0, length(totals), length(sp), dimnames = list(names(totals), sp))
  e0 <- setNames(as.list(rep(0, length(sp))), sp)
  for (i in seq_along(totals)) {
    expr <- parse(text = totals[[i]])
    base <- eval(expr, envir = e0)
    for (k in intersect(getSymbols(totals[[i]]), sp)) {
      e1 <- e0; e1[[k]] <- 1
      C[i, k] <- eval(expr, envir = e1) - base
    }
  }
  C
}


#' Detect states that are structurally zero in steady state
#'
#' Iterated three-layer test on the stoichiometric matrix of an `eqnlist`
#' (matches AlyssaPetit v1.2):
#' \enumerate{
#'   \item *Neg-only column*: only outflux, must be zero at SS.
#'   \item *Pos-only column with single-state feeder*: each feeding flux
#'     vanishes at SS; if a feeder's rate involves exactly one state, that
#'     state must be zero.
#'   \item *Sink cluster (LP)*: subsets whose combined mass leaks
#'     monotonically (mass-balance LP via `lpSolve::lp`).
#' }
#' For each zero-state the column is dropped, the state symbol is substituted
#' by `"0"` in remaining rates, structurally-zero reactions are removed, and
#' detection re-runs (removals can expose new zero-states).
#'
#' @param eqnlist_obj An [eqnlist].
#' @return `list(zero_states, eqnlist)` with the reduced system.
#' @keywords internal
.zeroStatesFromSmatrix <- function(eqnlist_obj) {
  S0 <- eqnlist_obj$smatrix
  if (is.null(S0) || ncol(S0) == 0L || nrow(S0) == 0L)
    return(list(zero_states = character(0), eqnlist = eqnlist_obj))

  S <- suppressWarnings(matrix(as.numeric(S0), nrow = nrow(S0), ncol = ncol(S0),
                               dimnames = dimnames(S0)))
  S[is.na(S)] <- 0
  rates       <- eqnlist_obj$rates
  description <- eqnlist_obj$description
  reactionCompartment <- eqnlist_obj$reactionCompartment

  zero_states <- character(0)

  # SBML rate rules may carry conditionals, so the structural probes below
  # evaluate them with SBML's flat `piecewise(v1, c1, ..., otherwise)`.
  ss_env <- new.env(parent = baseenv())
  ss_env$piecewise <- function(...) {
    a <- list(...); n <- length(a); i <- 1L
    while (i + 1L <= n) {
      if (isTRUE(a[[i + 1L]])) return(a[[i]])
      i <- i + 2L
    }
    if (n %% 2L == 1L) a[[n]] else 0
  }

  .is_struct_zero <- function(expr) {
    v <- tryCatch({
      syms <- getSymbols(expr)
      if (!length(syms)) eval(parse(text = expr), ss_env)
      else eval(parse(text = replaceSymbols(syms, rep("1", length(syms)), expr)),
                ss_env)
    }, error = function(e) NA_real_)
    isTRUE(v == 0)
  }

  # "Pure influx whose rate involves one state" implies that state is zero at
  # steady state only if the rate cannot vanish at a nonzero value. A
  # conditional rate can, so probe before concluding.
  .vanishes_only_at_zero <- function(expr, st) {
    others <- setdiff(getSymbols(expr), st)
    base <- if (length(others))
              replaceSymbols(others, rep("1", length(others)), expr) else expr
    vals <- vapply(c(1, 1e3, 1e6), function(v) {
      e <- replaceSymbols(st, format(v, scientific = FALSE), base)
      out <- tryCatch(eval(parse(text = e), ss_env), error = function(err) NA_real_)
      if (length(out) == 1L) as.numeric(out) else NA_real_
    }, numeric(1))
    !any(!is.na(vals) & vals == 0)
  }

  ## Drop a-priori-zero rates: their +1 stoichiometry otherwise masks the
  ## sink-cluster LP downstream.
  drop_rate0 <- vapply(rates, .is_struct_zero, logical(1))
  if (any(drop_rate0)) {
    keep <- !drop_rate0
    S     <- S[keep, , drop = FALSE]
    rates <- rates[keep]
    description <- description[keep]
    if (!is.null(reactionCompartment))
      reactionCompartment <- reactionCompartment[keep]
  }

  ## No influx (column has no strictly positive entry) => zero at SS;
  ## also catches all-zero columns from prior reductions.
  .neg_col <- function(M) {
    for (j in seq_len(ncol(M))) if (!any(M[, j] > 0)) return(j)
    NA_integer_
  }

  .pos_col_zero_state <- function(M, rates_chr) {
    cn <- colnames(M)
    for (j in seq_len(ncol(M))) {
      col <- M[, j]
      if (any(col < 0) || !any(col > 0)) next
      for (k in which(col > 0)) {
        in_rate <- intersect(getSymbols(rates_chr[k]), cn)
        if (length(in_rate) == 1L &&
            .vanishes_only_at_zero(rates_chr[k], in_rate)) return(in_rate)
      }
    }
    NA_character_
  }

  .sink_cluster <- function(M, eps = 1e-8, Mbig = 1e4) {
    .require_ns("lpSolve", "steady-state sink-cluster detection")
    nF <- nrow(M); nS <- ncol(M)
    if (nF == 0L || nS == 0L) return(integer(0))
    c_obj <- colSums(M)
    id    <- diag(nS)
    for (i in seq_len(nS)) {
      lb <- rep(0, nS); ub <- rep(Mbig, nS); lb[i] <- 1; ub[i] <- 1
      res <- tryCatch(
        lpSolve::lp("min", c_obj,
                    rbind(M, id, id),
                    c(rep("<=", nF), rep(">=", nS), rep("<=", nS)),
                    c(rep(0,  nF), lb, ub)),
        error = function(e) NULL)
      if (is.null(res) || res$status != 0 || res$objval >= -eps) next
      # The support may add a conserved moiety to a leaking cluster. Only its
      # species that reach a leaking reaction along reactions of the support
      # are 0: educts of the leaking reactions, then the educts of every
      # reaction producing a species already found.
      w    <- res$solution
      sup  <- which(w > eps)
      leak <- which(drop(M %*% w) < -eps)
      out  <- intersect(sup, which(colSums(M[leak, , drop = FALSE] < 0) > 0))
      repeat {
        feed <- which(rowSums(M[, out, drop = FALSE] > 0) > 0)
        add  <- setdiff(intersect(sup, which(colSums(M[feed, , drop = FALSE] < 0) > 0)), out)
        if (!length(add)) break
        out <- c(out, add)
      }
      if (length(out)) return(sort(out))
    }
    integer(0)
  }

  ## Substitute a state -> 0 stoichiometrically AND kinetically; drop any
  ## reaction whose rate is structurally zero after substitution.
  .zero_out <- function(j) {
    state_name <- colnames(S)[j]
    zero_states <<- c(zero_states, state_name)
    keep_rxn <- rep(TRUE, nrow(S))
    for (k in seq_len(nrow(S))) {
      in_stoich <- S[k, j] != 0
      in_rate   <- state_name %in% getSymbols(rates[k])
      if (!in_stoich && !in_rate) next
      new_rate <- if (in_rate) replaceSymbols(state_name, "0", rates[k]) else rates[k]
      if (.is_struct_zero(new_rate)) keep_rxn[k] <- FALSE
      else rates[k] <<- new_rate
    }
    S    <<- S[keep_rxn, -j, drop = FALSE]
    rates       <<- rates[keep_rxn]
    description <<- description[keep_rxn]
    if (!is.null(reactionCompartment))
      reactionCompartment <<- reactionCompartment[keep_rxn]
  }

  repeat {
    progressed <- FALSE
    repeat {
      step <- FALSE
      jneg <- .neg_col(S)
      if (!is.na(jneg)) { .zero_out(jneg); step <- TRUE; progressed <- TRUE; next }
      pzs <- .pos_col_zero_state(S, rates)
      if (!is.na(pzs)) {
        j <- match(pzs, colnames(S))
        if (!is.na(j)) { .zero_out(j); step <- TRUE; progressed <- TRUE; next }
      }
      if (!step) break
    }
    sink <- .sink_cluster(S)
    if (!length(sink)) break
    for (j in sort(sink, decreasing = TRUE)) { .zero_out(j); progressed <- TRUE }
    if (!progressed) break
  }

  if (!length(zero_states))
    return(list(zero_states = character(0), eqnlist = eqnlist_obj))

  S_out <- S
  S_out[S_out == 0] <- NA
  storage.mode(S_out) <- storage.mode(S0)

  new_obj <- eqnlist_obj
  new_obj$smatrix     <- S_out
  new_obj$states      <- colnames(S_out)
  new_obj$rates       <- rates
  new_obj$description <- description
  if (!is.null(reactionCompartment))
    new_obj$reactionCompartment <- reactionCompartment
  if (!is.null(new_obj$compartmentOf))
    new_obj$compartmentOf <- new_obj$compartmentOf[colnames(S_out)]

  list(zero_states = zero_states, eqnlist = new_obj)
}


#' Reset warm-start caches in `Pimpl` parameter transformations
#'
#' Walks `fn` and its closure environments and clears the warm-start
#' cache on every reachable [Pimpl] parfn. Use before workflows
#' that cross basins of attraction (multistart, profile after a structural
#' change) where a stale root pins the solver in the wrong region.
#'
#' @param fn A `parfn`, `prdfn`, `obsfn`, `objfn`, or composed `fn`.
#' @param verbose Print one-line summary of cleared caches.
#' @return Invisibly, labels of the cleared caches (empty if none found).
#' @export
resetWarmStarts <- function(fn, verbose = TRUE) {
  if (!is.function(fn))
    stop("`fn` must be a function (parfn / prdfn / obsfn / objfn / composed fn).",
         call. = FALSE)

  visited_envs   <- new.env(parent = emptyenv())
  invoked_resets <- new.env(parent = emptyenv())
  labels         <- character(0)

  call_reset <- function(r) {
    key <- format(environment(r))
    if (key %in% names(invoked_resets)) return()
    assign(key, TRUE, envir = invoked_resets)
    new_labels <- tryCatch(r(),
                           error = function(e) sprintf("<reset error: %s>",
                                                       conditionMessage(e)))
    labels <<- c(labels, as.character(new_labels))
  }

  walk <- function(x) {
    env <- NULL
    if (is.function(x)) {
      r <- attr(x, "resetWarmStart")
      if (!is.null(r) && is.function(r)) call_reset(r)
      env <- environment(x)
    } else if (is.environment(x)) {
      env <- x
    } else if (is.list(x)) {
      for (v in x) walk(v)
      return()
    } else {
      return()
    }
    if (is.null(env)) return()
    key <- format(env)
    if (key %in% names(visited_envs)) return()
    assign(key, TRUE, envir = visited_envs)
    for (nm in ls(env, all.names = TRUE)) {
      val <- tryCatch(get(nm, envir = env, inherits = FALSE),
                      error = function(e) NULL)
      if (is.function(val) || is.environment(val) || is.list(val)) walk(val)
    }
  }
  walk(fn)

  if (isTRUE(verbose)) {
    if (length(labels))
      message("resetWarmStarts: cleared ", length(labels),
              " warm-start cache(s):\n  ",
              paste(labels, collapse = "\n  "))
    else
      message("resetWarmStarts: no warm-start caches found.")
  }
  invisible(labels)
}


## Solve A X = B, A square or tall with consistent rows, by QR on A scaled to
## relative coordinates: columns by cs, rows by their norm in those columns. A
## rank-deficient A uses the pseudoinverse (minimum-norm sensitivity) and warns
## with the null-space directions.
#' @keywords internal
.pimpl_solve_dfdx <- function(A, B, cs = rep(1, ncol(A))) {
  if (ncol(A) == 0L) return(B[0L, , drop = FALSE])
  A0 <- A
  A  <- sweep(A, 2L, cs, `*`)
  rs <- pmax(rowSums(abs(A)), .Machine$double.xmin)
  A  <- A / rs; B <- B / rs
  dimnames(A) <- dimnames(A0)
  sv  <- svd(A); d <- sv$d
  tol <- max(dim(A)) * d[1L] * .Machine$double.eps
  rnk <- sum(d > tol)

  if (rnk == ncol(A)) {
    X <- qr.coef(qr(A), B)
  } else {
    nd <- sv$v[, (rnk + 1L):ncol(sv$v), drop = FALSE]
    rownames(nd) <- colnames(A)
    warning(.pimpl_format_singularity(d, rnk, nd), call. = FALSE)
    inv_d <- ifelse(d > tol, 1 / d, 0)
    X <- sv$v %*% (inv_d * crossprod(sv$u, B))
  }
  X <- cs * matrix(X, ncol(A), ncol(B))
  dimnames(X) <- list(colnames(A), colnames(B))
  X
}

## Zero forcing symbols in all rates and drop their equations, then return
## an eqnvec. Forcings (a state held at 0) leave the equation set so they
## are never solved for. Shared by the steady-state preamble and by
## symmetryDetection().
#' @keywords internal
.zero_and_drop_forcings <- function(trafo, forcings) {
  if (is.null(forcings) || !length(forcings)) return(as.eqnvec(trafo))
  if (inherits(trafo, "eqnlist"))
    trafo$rates <- replaceSymbols(forcings, rep("0", length(forcings)), trafo$rates)
  else
    trafo <- replaceSymbols(forcings, rep("0", length(forcings)), trafo)
  trafo <- as.eqnvec(trafo)
  trafo[setdiff(names(trafo), forcings)]
}

#' Normalise steady-state inputs
#'
#' Coerces `trafo` to an `eqnvec`, zeroes and removes forcings, removes states
#' that are structurally zero at steady state and promotes states with a zero
#' right-hand side to parameters.
#'
#' @param trafo An [eqnlist], [eqnvec] or named character vector.
#' @param parameters Outer parameters (may be `NULL`).
#' @param forcings Forcing names.
#' @return A named list with `trafo`, `states`, `zero_states`, `dependent`,
#'   `parameters`, `parms_all`, `totals`, `C_mat` (see [.cq_matrix()]) and,
#'   for an [eqnlist], `influx`, the sum of the producing fluxes per state.
#' @keywords internal
.normalize_ss_inputs <- function(trafo, parameters, forcings) {
  zero_states <- character(0); influx <- NULL
  original_params <- character(0)
  totals <- list()

  if (inherits(trafo, "eqnlist")) {
    original_params <- setdiff(getParameters(trafo), trafo$states)
    if (!is.null(forcings))
      trafo$rates <- replaceSymbols(forcings, rep("0", length(forcings)), trafo$rates)
    zs <- .zeroStatesFromSmatrix(trafo)
    zero_states <- zs$zero_states
    trafo       <- zs$eqnlist
    if (!length(trafo$states))
      stop("All states are structurally zero in steady state; no dynamical ",
           "state remains to solve for. The network likely has irreversible ",
           "drains without matching influx (an open system with trivial ",
           "all-zero equilibrium).", call. = FALSE)
    totals <- getTotals(trafo)
    influx <- vapply(getFluxes(trafo), function(t) {
      pos <- trimws(t[!startsWith(trimws(t), "-")])
      if (length(pos)) paste(sub("^\\+", "", pos), collapse = " + ") else "0"
    }, "")
  } else if (inherits(trafo, "eqnvec") || is.character(trafo)) {
    if (!is.null(forcings))
      trafo <- replaceSymbols(forcings, rep("0", length(forcings)), trafo)
  } else stop("'trafo' must be an eqnlist, eqnvec or character vector", call. = FALSE)
  trafo <- as.eqnvec(trafo)
  if (!is.null(forcings)) trafo <- trafo[setdiff(names(trafo), forcings)]

  states <- names(trafo)
  const_states <- states[vapply(unclass(trafo), function(x)
    tryCatch(identical(eval(parse(text = x)), 0),
             error = function(e) FALSE), logical(1))]
  if (length(const_states))
    parameters <- union(parameters %||% character(0), const_states)

  dependent <- setdiff(states, parameters %||% character(0))
  if (!length(dependent))
    stop("No dependent states to solve for. All states appear in 'parameters'.",
         call. = FALSE)
  parameters <- Reduce(union, list(getSymbols(trafo[dependent], exclude = dependent),
                                   parameters %||% character(0),
                                   original_params, names(totals)))
  list(trafo = trafo, states = states, zero_states = zero_states,
       dependent = dependent, parameters = parameters,
       parms_all = setdiff(parameters, dependent),
       totals = totals, C_mat = .cq_matrix(totals, states),
       influx = if (!is.null(influx)) influx[intersect(dependent, names(influx))])
}

#' @keywords internal
.pimpl_format_singularity <- function(d, rnk, null_dirs) {
  lines <- vapply(seq_len(ncol(null_dirs)), function(j) {
    v   <- null_dirs[, j]
    sig <- abs(v) > 0.05 * max(abs(v)); if (!any(sig)) sig <- abs(v) == max(abs(v))
    paste0("  null #", j, ": ",
           paste(sprintf("%+.3f*%s", v[sig], rownames(null_dirs)[sig]), collapse = " "))
  }, character(1))
  paste0(
    sprintf("df/dx is rank-deficient (rank %d of %d; smallest sv %.2e of largest %.2e); ",
            rnk, length(d), d[length(d)], d[1L]),
    "using SVD pseudoinverse (minimum-norm sensitivity on the constraint manifold).\n",
    "Null-space direction(s) in dependent-state coordinates:\n",
    paste(lines, collapse = "\n"),
    "\nLikely cause: an unmodelled conserved quantity, redundant equations, or a ",
    "continuum of steady states. Pass an `eqnlist` so `Pimpl` can auto-detect CQs ",
    "for cleaner sensitivities."
  )
}


# Pimpl solver controls merged over their defaults. Unknown names are an error.
.pimplPTC <- function(controlsPTC) {
  def <- list(rtol = 1e-10, atol = 1e-14, flowTol = 1, maxit = NULL, dtInit = 1e-2,
              positive = TRUE, stability = TRUE, archive = 8L,
              nStarts = 20L, startRange = c(-5, 5), startScale = NULL, seed = 1L)
  given <- as.list(controlsPTC)
  if (length(given) && (is.null(names(given)) || any(!nzchar(names(given)))))
    stop("Pimpl: controlsPTC must be a named list.", call. = FALSE)
  bad <- setdiff(names(given), names(def))
  if (length(bad))
    stop("Pimpl: unknown controlsPTC entr", if (length(bad) > 1L) "ies " else "y ",
         paste(bad, collapse = ", "), ". Known: ", paste(names(def), collapse = ", "), ".",
         call. = FALSE)
  ctrl <- modifyList(def, given, keep.null = TRUE)
  if (is.null(ctrl$startScale)) ctrl$startScale <- if (ctrl$positive) "log10" else "linear"
  ctrl$startScale <- match.arg(ctrl$startScale, c("log10", "linear"))
  if (length(ctrl$startRange) != 2L || !(ctrl$startRange[1] < ctrl$startRange[2]))
    stop("Pimpl: startRange must be c(lower, upper) with lower < upper.", call. = FALSE)
  if (ctrl$positive && ctrl$startScale == "linear" && ctrl$startRange[1] < 0)
    stop("Pimpl: a linear startRange below 0 needs positive = FALSE.", call. = FALSE)
  ctrl
}

#' Implicit parameter transformation
#'
#' Solves \eqn{f(x, p) = 0} for the states \eqn{x} by pseudo-transient
#' continuation ([cppDE::ptc()]) and returns a [parfn] with sensitivities from
#' the implicit function theorem.
#'
#' @param trafo Named character, [eqnvec] or [eqnlist] defining \eqn{f}.
#' @param parameters Outer parameter names. For an [eqnlist] the totals
#'   \eqn{T} of the conserved quantities are added.
#' @param forcings Forcing names, set to 0.
#' @param condition Condition label.
#' @param keep.root If `TRUE`, roots are kept per condition as initial guesses
#'   and repeated calls are answered from memory.
#' @param flow If `TRUE`, \eqn{\dot{x} = f(x, p)}{dx/dt = f(x, p)} and a stable
#'   steady state is returned. If `FALSE`, any regular root.
#' @param compile,modelname,verbose Forwarded to [cppDE::cppFUN()].
#' @param deriv,deriv2 Attach first and second order sensitivities. `deriv2`
#'   requires `deriv = TRUE`.
#' @param controlsPTC Named list of solver controls:
#'   \describe{
#'     \item{`rtol`, `atol`}{Convergence tolerances, default `1e-10`, `1e-14`.}
#'     \item{`flowTol`}{Relative local error per pseudo-time step, default `1`.}
#'     \item{`maxit`}{Iterations per start, default
#'       `ceiling(70 * log(n + 1))` for `n` states.}
#'     \item{`dtInit`}{Initial pseudo-time step relative to the fastest rate,
#'       default `1e-2`.}
#'     \item{`positive`}{Keep states positive, default `TRUE`.}
#'     \item{`stability`}{Require a stable root if `flow = TRUE`, default
#'       `TRUE`.}
#'     \item{`archive`}{Roots kept per condition, default `8`.}
#'     \item{`nStarts`}{Random starts after the deterministic ones, default
#'       `20`.}
#'     \item{`startRange`}{Range of a random start, drawn uniformly per
#'       state on `startScale`, default `c(-5, 5)`.}
#'     \item{`startScale`}{`"log10"` or `"linear"`, default `"log10"` if
#'       `positive` and `"linear"` otherwise.}
#'     \item{`seed`}{Seed of the random starts, default `1`. Each warm-start
#'       cache (per condition, emptied by [resetWarmStarts()] before every fit
#'       of [mstrust()]) draws from its own stream, derived from `seed`, the
#'       condition and the global RNG state when the cache first needs one. The
#'       global RNG is read, not advanced.}
#'   }
#' @param outdir Directory for the generated files.
#'
#' @details Conserved quantities of an [eqnlist] enter as \eqn{C x = T}.
#' States without influx are 0. Initial guesses are kept roots, then the
#' states in `pars`, missing ones at 1, then `nStarts` random starts. If all
#' fail, it is an error. Identical parameter values are solved
#' once, also across conditions. `keep.root` and `controlsPTC` can be changed
#' with [controls()].
#'
#' @return A [parfn].
#' @seealso [Pexpl()], [P()]
#' @export
#' @import cppDE
#' @importFrom digest digest
Pimpl <- function(trafo, parameters = NULL, forcings = NULL, condition = NULL,
                  keep.root = TRUE, flow = inherits(trafo, "eqnlist"),
                  compile = FALSE, modelname = NULL, verbose = FALSE,
                  deriv = TRUE, deriv2 = FALSE, controlsPTC = list(),
                  outdir = getwd()) {

  flow    <- isTRUE(flow)
  emit_d1 <- isTRUE(deriv)
  emit_d2 <- isTRUE(deriv2)
  if (emit_d2 && !emit_d1)
    stop("Pimpl(deriv2 = TRUE) requires deriv = TRUE.", call. = FALSE)
  .pimplPTC(controlsPTC)

  norm <- .normalize_ss_inputs(trafo, parameters, forcings)
  zero_states <- norm$zero_states
  dependent   <- norm$dependent
  parameters  <- norm$parameters
  n_dep       <- length(dependent)

  if (is.null(modelname)) modelname <- "impl_parfn"
  if (!is.null(condition)) modelname <- paste(modelname, sanitizeConditions(condition), sep = "_")

  # conservation rows C x - T; species that are parameters stay in the expression
  tn <- if (!is.null(norm$C_mat)) names(norm$totals) else character(0)
  cons <- if (length(tn))
    setNames(paste0("(", unlist(norm$totals), ") - ", tn), tn) else character(0)
  C_dep <- matrix(0, length(tn), n_dep, dimnames = list(tn, dependent))
  if (length(tn)) {
    cs <- intersect(colnames(norm$C_mat), dependent)
    C_dep[, cs] <- norm$C_mat[tn, cs, drop = FALSE]
  }
  all_exprs <- c(unclass(norm$trafo[dependent]), cons)
  n_eq      <- length(all_exprs)
  parms_all <- intersect(norm$parms_all, getSymbols(all_exprs))

  PEval <- suppressWarnings(cppDE::cppFUN(
    all_exprs, variables = dependent, parameters = parms_all, fixed = NULL,
    compile = compile, modelname = modelname, outdir = outdir,
    verbose = verbose, convenient = FALSE,
    deriv = TRUE, deriv2 = emit_d2, derivMode = "forward"))

  # producing fluxes per state, for states that rest at 0 at the given values
  influx <- if (flow) norm$influx
  PIn <- if (length(influx)) suppressWarnings(cppDE::cppFUN(
    unclass(influx), variables = dependent,
    parameters = intersect(parms_all, getSymbols(unclass(influx))), fixed = NULL,
    compile = compile, modelname = paste0(modelname, "_influx"), outdir = outdir,
    verbose = verbose, convenient = FALSE, deriv = FALSE))

  rowsOf <- function(x) {
    x <- if (is.matrix(x)) x else matrix(x[dependent], 1, dimnames = list(NULL, dependent))
    x
  }
  eval_G <- function(x, pv) {
    F <- PEval$func(rowsOf(x), pv[parms_all])
    matrix(F, dim(F)[1], dim(F)[2], dimnames = list(NULL, dimnames(F)[[2]]))
  }
  eval_J <- function(x, pv) {
    J <- PEval$jac(rowsOf(x), pv[parms_all])
    matrix(c(J), dim(J)[2], dim(J)[3], dimnames = list(dimnames(J)[[2]], dimnames(J)[[3]]))
  }
  eval_H <- function(x, pv) {
    if (is.null(PEval$hess)) return(NULL)
    H4 <- PEval$hess(rowsOf(x), pv[parms_all])
    array(c(H4), dim(H4)[2:4], dimnames = dimnames(H4)[2:4])
  }

  reg    <- .warmstart_registry()
  solved <- new.env(parent = emptyenv())   # roots and sensitivities by parameter values
  # solve counts and successful starts
  stats  <- new.env(parent = emptyenv())
  statsReset <- function() {
    stats$calls <- 0L; stats$memo <- 0L; stats$solves <- 0L; stats$iter <- 0L
    stats$how <- integer(0); stats$failed <- 0L
  }
  statsReset()
  solvedKeys <- character(0)
  remember <- function(key, val) {
    assign(key, val, envir = solved)
    solvedKeys <<- c(setdiff(solvedKeys, key), key)
    if (length(solvedKeys) > 64L) {
      rm(list = solvedKeys[1L], envir = solved); solvedKeys <<- solvedKeys[-1L]
    }
  }

  # read at every call, see controls()
  controls <- list(keep.root = keep.root, controlsPTC = controlsPTC)

  # States whose producing fluxes vanish and whose right-hand side is negative
  # at two fixed positive probes rest at 0. Repeated with the states found set
  # to 0 until no further state is found.
  probes <- rbind(exp(sin(seq_len(n_dep))), exp(cos(seq_len(n_dep))))
  colnames(probes) <- dependent
  numericZeros <- function(pv) {
    if (is.null(PIn)) return(character(0))
    pin  <- intersect(parms_all, getSymbols(unclass(influx)))
    zero <- character(0)
    repeat {
      X <- probes; X[, zero] <- 0
      F <- PIn$func(X, pv[pin])
      F <- matrix(F, dim(F)[1], dim(F)[2], dimnames = list(NULL, dimnames(F)[[2]]))
      R <- eval_G(X, pv)[, colnames(F), drop = FALSE]
      new <- setdiff(colnames(F)[F[1, ] == 0 & F[2, ] == 0 & R[1, ] < 0 & R[2, ] < 0], zero)
      if (!length(new)) break
      zero <- c(zero, new)
    }
    zero
  }

  # Projection onto C x = T_eff minimising relative entropy: x exp(C' lambda),
  # lambda by Newton. Keeps states positive; coefficients of any sign.
  onManifold <- function(x, act, Teff) {
    if (!length(Teff)) return(x)
    Ca <- C_dep[, act, drop = FALSE]
    xa <- x[act]; lam <- numeric(nrow(Ca))
    res <- function(l) drop(Ca %*% (xa * exp(drop(crossprod(Ca, l))))) - Teff
    r <- res(lam); tol <- 1e-15 * pmax(abs(Teff), drop(abs(Ca) %*% xa))
    for (it in 1:100) {
      if (all(abs(r) <= tol)) break
      e <- xa * exp(drop(crossprod(Ca, lam)))
      H <- Ca %*% (e * t(Ca))
      step <- tryCatch(solve(H, r), error = function(err) NULL)
      if (is.null(step)) break
      a <- 1
      repeat {
        rn <- res(lam - a * step)
        if (all(is.finite(rn)) && sum(rn^2) < sum(r^2) || a < 1e-8) break
        a <- a / 2
      }
      lam <- lam - a * step; r <- rn
    }
    x[act] <- xa * exp(drop(crossprod(Ca, lam)))
    x
  }

  # largest real part of the eigenvalues of J on ker C
  stabilityOf <- function(x, pv, act) {
    J  <- eval_J(x, pv)[dependent, , drop = FALSE][act, act, drop = FALSE]
    Ca <- C_dep[, act, drop = FALSE]
    Nb <- if (nrow(Ca)) MASS::Null(t(Ca)) else diag(length(act))
    if (!ncol(Nb)) return(list(maxRe = -Inf))
    ev <- eigen(t(Nb) %*% J %*% Nb)
    k  <- which.max(Re(ev$values))
    list(maxRe = Re(ev$values[k]), maxIm = abs(Im(ev$values[k])), scale = max(Mod(ev$values)),
         dir = setNames(as.numeric(Nb %*% Re(ev$vectors[, k])), act))
  }

  # starts in a fixed order, then the stability check
  solveRoot <- function(pv, x_user, arch, ctrl, cache = NULL, condition = NULL) {
    if (is.null(ctrl$maxit)) ctrl$maxit <- ceiling(70 * log(n_dep + 1))
    zero <- numericZeros(pv)
    act  <- setdiff(dependent, zero)
    if (!length(act)) return(list(x = setNames(rep(0, n_dep), dependent), zero = zero, how = "all zero"))
    # totals carried by the dependent species
    Teff <- if (length(tn)) -eval_G(setNames(rep(0, n_dep), dependent), pv)[1, tn] else numeric(0)
    Ca <- C_dep[, act, drop = FALSE]
    if (length(Teff) && any(Teff > 0 & rowSums(Ca != 0) == 0))
      stop("Pimpl: a conserved total is positive but every species carrying it rests at 0.",
           call. = FALSE)
    prep <- function(x) {
      x <- x[dependent]; x[zero] <- 0
      if (ctrl$positive) x[act][!(x[act] > 0)] <- 1
      onManifold(x, act, Teff)
    }
    attempts <- list()
    if (length(arch)) {
      d  <- vapply(arch, function(a) sum((log(pmax(abs(a$pv), 1e-300)) -
                                           log(pmax(abs(pv[names(a$pv)]), 1e-300)))^2), 0)
      a  <- arch[[which.min(d)]]
      xw <- a$x
      if (!is.null(a$dxdp)) {
        xp <- xw + drop(a$dxdp %*% (pv[colnames(a$dxdp)] - a$pv[colnames(a$dxdp)]))
        if (all(is.finite(xp)) && (!ctrl$positive || all(xp[act] > 0))) xw <- xp
      }
      attempts$warm_newton <- list(x = prep(xw), dt = if (flow) 1e8 else 1e8)
      attempts$warm_flow   <- list(x = prep(a$x), dt = if (flow) ctrl$dtInit else 1)
    }
    attempts$guess <- list(x = prep(x_user), dt = if (flow) ctrl$dtInit else 1)
    if (flow) attempts$guess_fine <- c(attempts$guess, fine = TRUE)
    if (any(x_user[act] != 1)) {
      attempts$ones <- list(x = prep(setNames(rep(1, n_dep), dependent)),
                            dt = if (flow) ctrl$dtInit else 1)
      if (flow) attempts$ones_fine <- c(attempts$ones, fine = TRUE)
    }

    log <- character(0)
    run <- function(x, dt, fine = FALSE) {
      r <- cppDE::ptc(PEval, x = x[dependent], parms = pv[parms_all], solve = act,
                      rows = act, C = if (length(Teff)) Ca, total = Teff, flow = flow,
                      positive = ctrl$positive,
                      controls = list(rtol = ctrl$rtol, atol = ctrl$atol,
                                      flowTol = ctrl$flowTol / if (fine) 10 else 1,
                                      maxit = ctrl$maxit, dtInit = dt))
      list(x = r$x[act], ok = r$converged, iter = r$iterations, reason = r$message)
    }
    # multistart: states uniform over startRange on startScale, from the stream of
    # this warm-start cache; the global RNG is read, not advanced
    if (ctrl$nStarts > 0L) {
      oldSeed <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
      on.exit(if (is.null(oldSeed)) rm(".Random.seed", envir = globalenv())
              else assign(".Random.seed", oldSeed, envir = globalenv()), add = TRUE)
      if (is.null(cache)) cache <- new.env(parent = emptyenv())
      if (is.null(cache$msBase)) {
        h <- digest::digest(list(ctrl$seed, condition, oldSeed), algo = "xxhash32")
        cache$msBase <- strtoi(substr(h, 1, 7), 16L); cache$msCount <- 0L
      }
      cache$msCount <- cache$msCount + 1L
      set.seed((cache$msBase + cache$msCount) %% .Machine$integer.max)
      for (i in seq_len(ctrl$nStarts)) {
        u <- stats::runif(n_dep, ctrl$startRange[1], ctrl$startRange[2])
        attempts[[paste0("start", i)]] <- list(
          x = prep(setNames(if (ctrl$startScale == "log10") 10^u else u, dependent)),
          dt = if (flow) ctrl$dtInit else 1)
      }
    }
    for (nm in names(attempts)) {
      r <- run(attempts[[nm]]$x, attempts[[nm]]$dt, isTRUE(attempts[[nm]]$fine))
      stats$iter <- stats$iter + r$iter
      log <- c(log, paste0(nm, ": ", r$reason))
      if (!r$ok) next
      x <- setNames(numeric(n_dep), dependent); x[act] <- r$x
      if (!(flow && ctrl$stability)) return(list(x = x, zero = zero, how = nm))
      st <- stabilityOf(x, pv, act)
      if (st$maxRe <= 1e-8 * st$scale) return(list(x = x, zero = zero, how = nm))
      # unstable: restart along the unstable eigenvector, both directions
      log <- c(log, sprintf("%s: unstable root, leading eigenvalue %.2e%s, leaving it", nm, st$maxRe,
                            if (st$maxIm > 0) sprintf(" +- %.2ei (oscillatory)", st$maxIm) else ""))
      v  <- st$dir; nz <- v != 0 & x[act] > 0
      if (!any(nz)) next
      eta <- 0.1 * min(x[act][nz] / abs(v[nz]))
      for (s in c(1, -1)) {
        r2 <- run(replace(x, act, x[act] + s * eta * v), ctrl$dtInit)
        if (!r2$ok) next
        x2 <- setNames(numeric(n_dep), dependent); x2[act] <- r2$x
        if (stabilityOf(x2, pv, act)$maxRe <= 1e-8 * st$scale)
          return(list(x = x2, zero = zero, how = paste0(nm, "+escape")))
      }
    }
    stats$failed <- stats$failed + 1L
    stop("Pimpl: no ", if (flow) "stable " else "", "root found. Attempts:\n  ",
         paste(log, collapse = "\n  "), call. = FALSE)
  }

  ift <- function(root, pv, want_d2) {
    Jall <- eval_J(root, pv)
    A <- Jall[, dependent, drop = FALSE]
    n_par <- length(parms_all)
    cs <- ifelse(root[dependent] > 0, root[dependent], 1)
    dxdp <- if (n_par) .pimpl_solve_dfdx(A, -Jall[, parms_all, drop = FALSE], cs)
            else matrix(numeric(0), n_dep, 0, dimnames = list(dependent, character(0)))
    d2 <- NULL
    if (want_d2) {
      H_all <- eval_H(root, pv)
      if (is.null(H_all))
        stop("Pimpl(deriv2 = TRUE) requires hess(); rebuild with deriv2 = TRUE.", call. = FALSE)
      if (n_par) {
        f_xx <- H_all[, dependent, dependent, drop = FALSE]
        f_xp <- H_all[, dependent, parms_all, drop = FALSE]
        f_pp <- H_all[, parms_all, parms_all, drop = FALSE]
        T1 <- t(dxdp) %bmm% f_xx %bmm% dxdp
        T2 <- t(dxdp) %bmm% f_xp
        RHS <- T1 + T2 + aperm(T2, c(1L, 3L, 2L)) + f_pp
        d2 <- array(.pimpl_solve_dfdx(A, -matrix(RHS, n_eq, n_par * n_par), cs),
                    c(n_dep, n_par, n_par), dimnames = list(dependent, parms_all, parms_all))
      } else d2 <- array(0, c(n_dep, 0L, 0L))
    }
    list(dxdp = dxdp, d2 = d2)
  }

  # output Jacobian and Hessian, chained with the derivatives of the inputs
  build_derivs <- function(dxdp, d2, out, p, emptypars, fixed, dP, dP2, want_d2) {
    input_cols <- setdiff(names(p), c(dependent, names(fixed)))
    n_in <- length(input_cols)
    par_input <- intersect(parms_all, input_cols)

    jacobian <- matrix(0, length(out), n_in, dimnames = list(names(out), input_cols))
    ep <- intersect(emptypars, input_cols)
    if (length(ep)) jacobian[cbind(ep, ep)] <- 1
    cd <- intersect(colnames(dxdp), input_cols)
    if (length(cd)) jacobian[dependent, cd] <- dxdp[, cd, drop = FALSE]

    hessian <- NULL
    if (want_d2) {
      hessian <- array(0, c(length(out), n_in, n_in),
                       dimnames = list(names(out), input_cols, input_cols))
      if (length(par_input))
        hessian[dependent, par_input, par_input] <- d2[, par_input, par_input, drop = FALSE]
    }

    if (!is.null(dP)) {
      dPsub <- submatrix(dP, rows = colnames(jacobian))
      th    <- colnames(dPsub); n_th <- length(th)
      if (want_d2) {
        new_hess <- t(dPsub) %bmm% hessian %bmm% dPsub
        dimnames(new_hess) <- list(names(out), th, th)
        if (!is.null(dP2)) {
          dP2sub <- dP2[input_cols, th, th, drop = FALSE]
          new_hess <- new_hess + array(
            jacobian %*% matrix(dP2sub, n_in, n_th * n_th),
            c(length(out), n_th, n_th), dimnames = list(names(out), th, th))
        }
        hessian <- new_hess
      }
      jacobian <- jacobian %*% dPsub
    }

    keep <- rowSums(jacobian != 0) > 0
    jacobian <- jacobian[keep, , drop = FALSE]
    if (!is.null(hessian)) hessian <- hessian[keep, , , drop = FALSE]
    list(jacobian = jacobian, hessian = hessian)
  }

  p2p <- function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE, condition = NULL) {
    if (deriv2 && !emit_d2)
      stop("Pimpl was built with deriv2 = FALSE; rebuild with deriv2 = TRUE.", call. = FALSE)
    if (!emit_d1) deriv <- FALSE
    if (deriv2 && !deriv) deriv <- TRUE

    keep.root <- controls$keep.root
    ctrl  <- .pimplPTC(controls$controlsPTC)
    cache <- reg$get(condition)
    p   <- pars
    dP  <- attr(p, "deriv")
    dP2 <- if (deriv2) attr(p, "deriv2") else NULL
    if (!is.null(fixed)) {
      p <- p[!names(p) %in% names(fixed)]
      p <- c(p, fixed)
    }
    emptypars <- setdiff(names(p), c(dependent, names(fixed)))
    miss <- setdiff(dependent, names(p)); if (length(miss)) p[miss] <- 1
    pv <- setNames(as.numeric(p[parms_all]), parms_all)
    if (anyNA(pv))
      stop("Pimpl: missing value(s) for ", paste(parms_all[is.na(pv)], collapse = ", "),
           call. = FALSE)

    # memo key: parameter values and solver controls
    key <- digest::digest(list(unname(pv), ctrl), algo = "xxhash64")
    got <- if (keep.root) get0(key, envir = solved, inherits = FALSE) else NULL
    stats$calls <- stats$calls + 1L
    if (is.null(got)) {
      arch <- if (keep.root) cache$arch else NULL
      got  <- solveRoot(pv, setNames(as.numeric(p[dependent]), dependent), arch, ctrl,
                        cache = cache, condition = condition)
      got$ift <- NULL
      stats$solves <- stats$solves + 1L
      stats$how[got$how] <- (if (is.na(stats$how[got$how])) 0L else stats$how[got$how]) + 1L
    } else stats$memo <- stats$memo + 1L
    root <- got$x
    need_ift <- deriv && (is.null(got$ift) || (deriv2 && is.null(got$ift$d2)))
    if (need_ift) {
      got$ift <- tryCatch(ift(root, pv, deriv2), error = function(e) {
        warning("Pimpl: IFT-based sensitivities unavailable at the current root (",
                conditionMessage(e), "). Returning value only.", call. = FALSE)
        NULL
      })
    }
    if (keep.root) {
      remember(key, got)
      entry <- list(pv = pv, x = root, dxdp = got$ift$dxdp)
      cache$arch <- c(list(entry), Filter(function(a) !identical(a$pv, pv), cache$arch))
      if (length(cache$arch) > ctrl$archive) cache$arch <- cache$arch[seq_len(ctrl$archive)]
    }

    zero_vec <- if (length(zero_states))
      setNames(rep(0, length(zero_states)), zero_states) else NULL
    out <- c(root, zero_vec, p[setdiff(names(p), c(dependent, zero_states))])

    d <- if (deriv && !is.null(got$ift))
      build_derivs(got$ift$dxdp, got$ift$d2, out, p, emptypars, fixed, dP, dP2, deriv2)
    as.parvec(out,
              deriv  = if (deriv  && !is.null(d)) d$jacobian else NULL,
              deriv2 = if (deriv2 && !is.null(d)) d$hessian  else if (deriv2) NULL else FALSE)
  }

  attr(p2p, "vjpfn")       <- .parfnVjpFromJacobian(p2p)
  attr(p2p, "equations")   <- as.eqnvec(all_exprs)
  attr(p2p, "parameters")  <- parameters
  attr(p2p, "modelname")   <- modelname
  attr(p2p, "compileInfo") <- .collectCompileInfo(PEval$func, PEval$jac, PEval$hess,
                                                  if (!is.null(PIn)) PIn$func)
  attr(p2p, "resetWarmStart") <- local({
    reg_ref <- reg; mn <- modelname; cond <- condition
    function() {
      reg_ref$reset()
      rm(list = ls(solved, all.names = TRUE), envir = solved)
      solvedKeys <<- character(0)
      paste0("Pimpl(", mn, if (!is.null(cond)) paste0(":", cond) else "", ")")
    }
  })
  parfn(p2p, parameters, condition)
}


## Values passed through `...` resolve against the condition row, then the
## per-branch symbols, then the calling frame. The last of those is what makes
## a call from inside a function work.
.evalDots <- function(dots, row, currentTrafo, currentSymbols, callerEnv) {
  env <- list2env(as.list(row), envir = new.env(parent = callerEnv))
  env$.currentTrafo   <- currentTrafo
  env$.currentSymbols <- currentSymbols
  lapply(eval(dots), eval, envir = env)
}


#' Construct and modify parameter transformations
#'
#' Symbolic helpers used by [P()] and [Xs()] to build, substitute, and
#' branch transformation rules. The condition table from [branch()] is
#' stored on `attr(., "tree")`.
#'
#' - `define` resets the LHS of `expr` to its RHS.
#' - `insert` substitutes the RHS for the LHS wherever it occurs.
#' - `branch` duplicates a trafo across conditions, optionally applying
#'   per-condition substitutions taken from `table`.
#'
#' @param trafo Named character / [eqnvec], or a list thereof.
#' @param expr `"lhs ~ rhs"` formula string.
#' @param table Condition table (row per condition, column per parameter)
#'   carried as the `tree` attribute when branching.
#' @param conditions Condition names; default `rownames(table)`.
#' @param apply One of `"nothing"`, `"insert"`, `"define"`; how the
#'   `table` entries are folded into each branch.
#' @param conditionMatch Regex on condition names; restricts the operation.
#' @param ... Named values to substitute into `expr` symbols.
#'
#' @return Same shape as `trafo` (or a per-condition list if branched).
#' @export
#' @example inst/examples/define.R
define <- function(trafo, expr, ..., conditionMatch = NULL) {
  if (missing(trafo)) trafo <- NULL
  tree <- attr(trafo, "tree")
  if (is.list(trafo) && is.null(names(trafo)))
    stop("If trafo is a list, elements must be named.", call. = FALSE)
  if (is.list(trafo) && !all(names(trafo) %in% rownames(tree)))
    stop("List names must be a subset of rownames(attr(trafo, 'tree')).", call. = FALSE)
  mytrafo <- if (is.list(trafo)) trafo else list(trafo)

  dots      <- substitute(alist(...))
  callerEnv <- parent.frame()
  out  <- lapply(seq_along(mytrafo), function(i) {
    .currentTrafo   <- mytrafo[[i]]
    .currentSymbols <- if (is.null(.currentTrafo)) NULL else getSymbols(.currentTrafo)
    row <- if (is.list(trafo)) tree[names(mytrafo)[i], , drop = FALSE]
           else tree[1, , drop = FALSE]
    if (!is.null(conditionMatch) && !str_detect(rownames(row), conditionMatch))
      return(.currentTrafo)
    args <- .evalDots(dots, row, .currentTrafo, .currentSymbols, callerEnv)
    do.call(repar, c(list(expr = expr, trafo = .currentTrafo, reset = TRUE), args))
  })
  names(out) <- names(mytrafo)
  if (!is.list(trafo)) out <- out[[1]]
  attr(out, "tree") <- tree
  out
}


#' @export
#' @rdname define
insert <- function(trafo, expr, ..., conditionMatch = NULL) {
  if (missing(trafo)) trafo <- NULL
  tree <- attr(trafo, "tree")
  if (is.list(trafo) && is.null(names(trafo)))
    stop("If trafo is a list, elements must be named.", call. = FALSE)
  if (is.list(trafo) && !all(names(trafo) %in% rownames(tree)))
    stop("List names must be a subset of rownames(attr(trafo, 'tree')).", call. = FALSE)
  mytrafo <- if (is.list(trafo)) trafo else list(trafo)

  dots      <- substitute(alist(...))
  callerEnv <- parent.frame()
  out  <- lapply(seq_along(mytrafo), function(i) {
    .currentTrafo   <- mytrafo[[i]]
    .currentSymbols <- if (is.null(.currentTrafo)) NULL else getSymbols(.currentTrafo)
    row <- if (is.list(trafo)) tree[names(mytrafo)[i], , drop = FALSE]
           else tree[1, , drop = FALSE]
    if (!is.null(conditionMatch) && !str_detect(rownames(row), conditionMatch))
      return(.currentTrafo)
    args <- .evalDots(dots, row, .currentTrafo, .currentSymbols, callerEnv)
    if (!length(args))
      return(do.call(repar, list(expr = expr, trafo = .currentTrafo)))
    ## Logical dots gate the substitution per condition rather than naming a
    ## symbol, so they decide and are then dropped.
    isGate <- vapply(args, is.logical, logical(1))
    gate   <- unlist(args[isGate], use.names = FALSE)
    if (length(gate) && any(!gate)) return(.currentTrafo)
    do.call(repar, c(list(expr = expr, trafo = .currentTrafo), args[!isGate]))
  })
  names(out) <- names(mytrafo)
  if (!is.list(trafo)) out <- out[[1]]
  attr(out, "tree") <- tree
  out
}


#' @export
#' @rdname define
branch <- function(trafo, table = NULL,
                   conditions = rownames(table),
                   apply = c("nothing", "insert", "define")) {
  apply <- match.arg(apply)
  if (is.null(table) && is.null(conditions)) return(trafo)
  if (is.null(conditions)) conditions <- paste0("C", seq_len(nrow(table)))
  if (is.null(table))      table      <- data.frame(condition = conditions,
                                                    row.names = conditions)
  rownames(table) <- conditions

  out <- setNames(lapply(conditions, function(x) trafo), conditions)
  attr(out, "tree") <- table
  if (apply == "nothing") return(out)

  for (cn in conditions) {
    row <- table[cn, !colnames(table) %in% c("condition", "conditions"), drop = FALSE]
    for (par in colnames(row)) {
      val <- row[[par]]; if (is.na(val)) next
      single <- out[[cn]]; attr(single, "tree") <- row
      out[[cn]] <- if (apply == "insert") insert(single, paste0(par, " ~ ", val))
                   else                    define(single, paste0(par, " ~ ", val))
      attr(out[[cn]], "tree") <- NULL
    }
  }
  out
}



## Split compound identifiers into single symbols: "_" -> ":" inside names
## only, digit-leading and empty parts prefixed to stay syntactic. Inverted by
## .decolonize.
.numprefix <- "..dModnum.."

.colonize <- function(x) {
  m <- gregexpr("[A-Za-z.][A-Za-z0-9._]*", x)
  regmatches(x, m) <- lapply(regmatches(x, m), function(ids) {
    vapply(ids, function(id) {
      parts <- strsplit(id, "_", fixed = TRUE)[[1]]
      # strsplit() drops trailing empty parts; `gamma_` must not become `gamma`.
      parts <- c(parts, rep("", attr(regexpr("_*$", id), "match.length")))
      digit <- grepl("^[0-9]", parts)
      parts[digit] <- paste0(.numprefix, parts[digit])
      parts[parts == ""] <- .numprefix
      paste(parts, collapse = ":")
    }, character(1), USE.NAMES = FALSE)
  })
  x
}

.decolonize <- function(x)
  gsub(.numprefix, "", gsub(":", "_", x, fixed = TRUE), fixed = TRUE)


#' Reparameterization
#'
#' Replaces symbols on either side of `"lhs ~ rhs"`. With `reset = TRUE`
#' the LHS entry of `trafo` is overwritten by the RHS (per row when `...`
#' supplies vector replacements); otherwise the LHS is substituted into
#' `trafo` wherever it occurs. Symbols separated by `_` are recognised as
#' compound identifiers, e.g. `Delta_x` -> symbols `"Delta"` and `"x"`.
#'
#' @param expr `"lhs ~ rhs"` string (or a formula).
#' @param trafo Character / [eqnvec] / list. `NULL` builds a fresh trafo
#'   from the LHS.
#' @param ... Named character/numeric vectors; each row supplies one
#'   substitution.
#' @param reset Overwrite (`TRUE`) or substitute into (`FALSE`).
#' @return Same shape as `trafo`.
#' @export
#' @importFrom stats as.formula
#' @examples
#' innerpars   <- letters[1:3]
#' constraints <- c(a = "b + c")
#' mycondition <- "cond1"
#' trafo <- repar("x ~ x",        x = innerpars)
#' trafo <- repar("x ~ y",        trafo, x = names(constraints), y = constraints)
#' trafo <- repar("x ~ exp(x)",   trafo, x = innerpars)
#' trafo <- repar("x ~ x + Delta_x_condition",
#'                trafo, x = innerpars, condition = mycondition)
repar <- function(expr, trafo = NULL, ..., reset = FALSE) {
  if (inherits(expr, "formula")) expr <- deparse(expr)
  parsed <- as.character(stats::as.formula(.colonize(expr)))
  lhs <- parsed[2]; rhs <- parsed[3]

  args <- lapply(list(...), as.character)
  if (length(args)) {
    reps <- as.data.frame(args, stringsAsFactors = FALSE)
    apply_repl <- function(side) vapply(seq_len(nrow(reps)), function(i)
      .decolonize(replaceSymbols(colnames(reps), reps[i, ], side)),
      character(1))
    lhs <- apply_repl(lhs); rhs <- apply_repl(rhs)
  } else {
    lhs <- .decolonize(lhs)
    rhs <- .decolonize(rhs)
  }

  if (is.null(trafo))                           as.eqnvec(structure(lhs, names = lhs))
  else if (is.list(trafo)      && !reset)       lapply(trafo, function(t) replaceSymbols(lhs, rhs, t))
  else if (is.character(trafo) && !reset)       replaceSymbols(lhs, rhs, trafo)
  else if (is.list(trafo)      &&  reset)       lapply(trafo, function(t) { t[lhs] <- rhs; t })
  else { trafo[lhs] <- rhs; trafo }
}

paste_ <- function(...) paste(..., sep = "_")