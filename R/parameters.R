## Functions to generate parameter transformation ----

#' Parameter Transformation
#'
#' @description Builds a [parfn] with [Pexpl()] or [Pimpl()].
#'
#' @param trafo An [eqnvec], named character or [eqnlist], or a list of
#'   those named by condition, which builds one transformation per condition.
#'   `NULL` (default) returns `NULL`. Expressions follow [cppDE::expressions].
#' @param parameters Outer parameter names. Default `NULL`, the symbols of
#'   `trafo`.
#' @param condition Condition label of a single `trafo`, default `NULL`; a
#'   list takes the list names instead.
#' @param compile Logical. `TRUE` compiles all conditions into one shared
#'   object named `modelname`, or with `modelname = NULL` into one shared
#'   object per condition; `FALSE` (default) leaves it to [compile()].
#' @param modelname Base name of the generated files, default `NULL`, which
#'   is `"expl_parfn"` or `"impl_parfn"`, followed by `_<condition>`.
#' @param verbose Logical. Print the compiler output, default `FALSE`.
#' @param method `"explicit"`, `"implicit"`, or `NULL` (default) for
#'   `"implicit"` if `trafo` is an [eqnlist] and `"explicit"` otherwise.
#' @param cores Number of parallel jobs generating and compiling the sources of
#'   a list; on Windows the sources are generated one after the other. `NULL`
#'   (default) uses [detectFreeCores()].
#' @param deriv,deriv2 Build first and second derivatives, default `TRUE` and
#'   `FALSE`. `deriv2` requires `deriv = TRUE`.
#' @param ... Further arguments of [Pexpl()] (`attachInput`, `derivMode`,
#'   `outdir`) or [Pimpl()] (`forcings`, `keepRoot`, `flow`, `controlsPTC`,
#'   `outdir`), as `method` selects.
#'
#' @return A [parfn], called as `p(pars, ...)`, see section Calling a dMod
#'   function. The random effects of a transformation built by [insert()]
#'   stay attached, so [randomEffects()] reads them from the result.
#' @inheritSection dModfn Calling a dMod function
#' @seealso [Pexpl()], [Pimpl()], [define()]
#' @example inst/examples/P.R
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
  cores <- if (is.null(cores)) detectFreeCores() else min(detectFreeCores(), cores)

  ## Codegen-only inside mclapply: dyn.load in a forked worker is lost when
  ## the fork exits, so the actual compile happens in the parent below.
  result <- Reduce("+", parallel::mclapply(seq_along(trafo_list), function(i) {
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
  }, mc.cores = if (.Platform$OS.type == "windows") 1L else cores))

  if (compile) compile(result, cores = cores, output = modelname, verbose = verbose)
  attr(result, "randomEffects") <- attr(trafo, "randomEffects", exact = TRUE)
  result
}


## Body of the parameter transformation built by Pexpl(). Package level, so
## one parfn per condition holds the state, not another copy of this code.

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
      active <- .setdiffU(st$parameters, names(fixed))
      dP <- diag(length(active)); dimnames(dP) <- list(active, active)
    }
    out <- if (!is.null(.ad_out)) .ad_out else
      st$evaluate(NULL, p[st$parameters], tangentX = NULL, tangentP = dP,
                  hessianX = NULL, hessianP = dP2,
                  deriv2 = deriv2, attach.input = attach.input,
                  fixed = as.character(.intersectU(names(fixed), st$parameters)))
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
  through <- if (attach.input) .setdiffU(names(p), names(pinnerVal)) else character(0)
  val <- if (length(through)) c(pinnerVal, .subset(p, through)) else pinnerVal

  ## A build without forward derivatives has none to give. Rows for the inputs
  ## alone would mark every output as fixed.
  if (!deriv || is.null(Jac))
    return(as.parvec(val, deriv = FALSE, deriv2 = FALSE))

  ## An output is fixed when the output expression has no varying input, not when the
  ## derivative vanishes at this point: the Hessian may not.
  varying <- if (is.null(attr(pars, "deriv"))) .setdiffU(names(pars), names(fixed))
             else rownames(attr(pars, "deriv"))
  # kept for the last set of varying inputs, which consecutive conditions mostly share
  dependent <- if (identical(varying, st$depKey) && identical(rownames(Jac), names(st$depVal))) st$depVal else {
    st$depKey <- varying
    st$depVal <- vapply(rownames(Jac), function(o)
      is.null(st$depends[[o]]) || any(st$depends[[o]] %in% varying), logical(1))
  }
  Jac  <- Jac[dependent | rowSums(Jac != 0) > 0, , drop = FALSE]
  Hess <- if (deriv2 && !is.null(Hess)) Hess[rownames(Jac), , , drop = FALSE]

  ## Every varying input keeps a derivative row: without one it counts as fixed
  ## downstream and the forward gradient along it would be zero. What the caller
  ## fixed has no row and stays fixed.
  moving <- .setdiffU(through, names(fixed))
  if (length(moving)) {
    d <- .Pexpl_through(pars, moving, colnames(Jac), second = !is.null(Hess))
    Jac <- .Pexpl_widen(Jac, d$theta)
    Jac <- rbind(Jac, d$deriv)
    if (!is.null(Hess))
      Hess <- .Pexpl_bind3(.Pexpl_widen(Hess, d$theta), d$deriv2)
  }
  as.parvec(val, deriv = Jac, deriv2 = if (is.null(Hess)) FALSE else Hess)
}


# Derivatives of the inputs Pexpl passes through, in the basis `theta` of the incoming
# Jacobian. Without incoming derivatives an input is a direction of itself, added to
# the basis if needed; an incoming input without a row was fixed and keeps none.
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
      active <- .setdiffU(st$parameters, names(fixed))
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


# The kernel with batch and reverse entries. All three close over this
# frame, which holds `st` and `controls` alone, and read `controls` when they
# run: controls<- rebinds it here, so a change reaches every entry at once.
.Pexpl_wrap <- function(st, controls) {
  p2p <- function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE)
    .Pexpl_p2p(st, pars, fixed, deriv, deriv2,
               attach.input = controls$attachInput)
  attr(p2p, "vjpfn") <- function(pars, fixed = NULL, cotangent, condition = NULL)
    .Pexpl_vjp(st, pars, fixed, cotangent, attach.input = controls$attachInput)
  attr(p2p, "batchfn") <- function(parsList, fixedList, deriv, deriv2,
                                   conditions, cores)
    .Pexpl_batch(st, parsList, fixedList, deriv, deriv2, cores,
                 attach.input = controls$attachInput)
  # the leaves of a sum of such kernels go to .Pexpl_many() together
  attr(p2p, "pexpl") <- TRUE
  p2p
}

# The requests of several Pexpl() kernels, each with a separate compiled function, in one
# evaluateMany() of cppDE; one by one where cppDE has none or a kernel lacks the forward
# derivatives asked for.
.Pexpl_many <- function(kernels, parsList, fixedList, deriv, deriv2, cores) {
  sts <- lapply(kernels, function(k) environment(k)$st)
  attach <- vapply(kernels, function(k) isTRUE(environment(k)$controls$attachInput), TRUE)
  many <- get0("evaluateMany", envir = asNamespace("cppDE"), inherits = FALSE)
  ok <- !is.null(many) && deriv && all(vapply(sts, function(st)
    st$use_ad && st$emit_d1 && !is.null(st$evaluateBatch) && (!deriv2 || isTRUE(st$emit_d2)), TRUE))
  if (!ok) return(lapply(seq_along(kernels), function(i)
    kernels[[i]](parsList[[i]], fixedList[[i]], deriv, deriv2)))
  sets <- lapply(seq_along(sts), function(i) {
    st <- sts[[i]]; pars <- parsList[[i]]; fixed <- fixedList[[i]]
    dP <- attr(pars, "deriv")
    if (is.null(dP)) {
      active <- .setdiffU(st$parameters, names(fixed))
      dP <- diag(length(active)); dimnames(dP) <- list(active, active)
    }
    list(vars = NULL, params = c(unclass(pars), unclass(fixed))[st$parameters],
         tangentP = dP, hessianP = if (deriv2) attr(pars, "deriv2"), attach.input = attach[i],
         fixed = .intersectU(names(fixed), st$parameters))
  })
  ad <- many(lapply(sts, `[[`, "evaluateBatch"), sets, cores = cores, deriv2 = deriv2)
  lapply(seq_along(sts), function(i)
    .Pexpl_p2p(sts[[i]], parsList[[i]], fixedList[[i]], deriv, deriv2, .ad_out = ad[[i]],
               attach.input = attach[i]))
}

# w' Jac, the reverse of the forward Jac %*% dP: the observation functions' vjp
# call with no variables and the cotangent on the inner parameters. Under
# attach.input the passed-through outer parameters add their own cotangent.
.Pexpl_vjp <- function(st, pars, fixed = NULL, cotangent, attach.input = FALSE) {
  if (is.null(st$vjp))
    stop("Pexpl(): the reverse mode needs a vector-Jacobian product; rebuild ",
         "with derivMode = c(\"forward\", \"reverse\") and compile = TRUE.",
         call. = FALSE)
  w <- .asCtPars(cotangent)
  K <- .ctK(w)
  b <- .ctBlocks(K)
  p <- c(pars, fixed)
  outnames <- st$outnames
  hit <- intersect(rownames(w), outnames)

  dp <- attr(pars, "deriv")
  if (b$D == 0L || is.null(dp) || !any(dp != 0)) {
    # No tangent arrives, so every slice is pulled back as the first, all as seeds of one call.
    W <- array(0, c(1L, length(outnames), K), dimnames = list(NULL, outnames, NULL))
    if (length(hit)) W[1L, hit, ] <- w[hit, , drop = FALSE]
    r <- st$vjp(NULL, p[st$parameters], W)
    u <- matrix(r$cotangentP, ncol = K, dimnames = list(rownames(r$cotangentP), NULL))
  } else {
    # The node has no variables, so only the parameters bring tangents in, and
    # the cotangent brings separate ones. Both halves of d/dv (w' J) come back in one
    # pass; nothing here forms a Hessian.
    nd <- b$D
    W <- array(0, c(1L, length(outnames), b$S), dimnames = list(NULL, outnames, NULL))
    if (length(hit)) W[1L, hit, ] <- w[hit, b$first, drop = FALSE]
    V <- matrix(0, length(st$parameters), nd,
                dimnames = list(st$parameters, NULL))
    take <- intersect(rownames(dp), st$parameters)
    if (length(take)) V[take, ] <- dp[take, seq_len(nd), drop = FALSE]
    DW <- array(0, c(1L, length(outnames), b$S, nd))
    if (length(hit))
      for (s in seq_len(b$S)) DW[1L, match(hit, outnames), s, ] <- w[hit, b$rest[, s], drop = FALSE]
    r <- st$vjp(NULL, p[st$parameters], W, tangentP = V, curvature = DW)
    u <- matrix(0, nrow(r$cotangentP), K, dimnames = list(rownames(r$cotangentP), NULL))
    u[, b$first] <- r$cotangentP
    for (s in seq_len(b$S)) u[, b$rest[, s]] <- matrix(r$curvatureP[, s, ], nrow = nrow(u))
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

#' Explicit Parameter Transformation
#'
#' @description Builds a [parfn] that evaluates inner parameters as
#' expressions in outer parameters, with their derivatives. It is evaluable
#' only after compilation.
#'
#' @param trafo Named character or [eqnvec]; names are inner parameters,
#'   values are expressions in the outer parameters, in the language of
#'   [cppDE::expressions].
#' @param parameters Outer parameters; default `NULL`, the symbols of
#'   `trafo`. Names given here that `trafo` does not map pass through as
#'   identity.
#' @param attachInput Logical. Append the inputs that `trafo` does not map;
#'   varying ones with identity derivatives, fixed ones as fixed. Default
#'   `FALSE`. Can be changed with [controls()].
#' @param condition `NULL` (default) for a transformation valid in every
#'   condition, or the name of the condition it belongs to.
#' @param compile Logical. `TRUE` compiles now, `FALSE` (default) leaves it
#'   to [compile()].
#' @param modelname Base name of the generated files, default `NULL`, which
#'   is `"expl_parfn"`, followed by `_<condition>`.
#' @param verbose Logical. Print the compiler output, default `FALSE`.
#' @param deriv,deriv2 Build the first derivatives, `attr(, "deriv")`,
#'   `[inner, outer]`, and the second, `attr(, "deriv2")`,
#'   `[inner, outer, outer]`. Default `TRUE` and `FALSE`; `deriv2` needs
#'   `deriv = TRUE`.
#' @param derivMode Which derivative code to build, any of `"forward"`
#'   (default), `"reverse"` and `"forward-reverse"`. `"reverse"` is needed for
#'   `sweep = "reverse"`, `"forward-reverse"` for `sweep = "reverse"` with
#'   `deriv2 = TRUE`.
#' @param outdir Directory for the generated source and shared object,
#'   default `getOption("dMod2.outdir")`, else the working directory.
#'
#' @return A [parfn], called as `p(pars, ...)`, see section Calling a dMod
#'   function.
#' @inheritSection dModfn Calling a dMod function
#' @seealso [Pimpl()], [P()]
#' @example inst/examples/Pexpl.R
#' @importFrom cppDE cppFUN
#' @export
Pexpl <- function(trafo, parameters = NULL, attachInput = FALSE, condition = NULL,
                  compile = FALSE, modelname = NULL, verbose = FALSE,
                  deriv = TRUE, deriv2 = FALSE,
                  derivMode = "forward",
                  outdir = .dmodOutdir()) {

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
                      outnames = names(trafo), modelname = modelname,
                      depends = lapply(setNames(as.character(trafo), names(trafo)), getSymbols),
                      parameters = parameters,
                      use_ad = use_ad, ad_symbol = ad_symbol,
                      ad2_symbol = ad2_symbol, emit_d1 = emit_d1,
                      emit_d2 = emit_d2), parent = emptyenv())
  p2p <- .Pexpl_wrap(st, controls = list(attachInput = attachInput))

  attr(p2p, "equations")   <- as.eqnvec(trafo)
  attr(p2p, "parameters")  <- parameters
  attr(p2p, "modelname")   <- modelname
  attr(p2p, "compileInfo") <- .collectCompileInfo(fun, jac, hess, evaluate)
  parfn(p2p, parameters, condition)
}


## The conserved totals of the network without the structurally zero states; one proportional to a
## total of the full network with the zero states at 0 is written as that total, under the name of that total.
## Totals of constant states (`held`) are left out.
.restingTotals <- function(full, reduced, zero, held) {
  red <- conservedQuantities(reduced)
  red <- Filter(function(e) !all(intersect(getSymbols(e), reduced$states) %in% held), red)
  if (!length(red) || is.null(full$smatrix)) return(red)
  ref <- conservedQuantities(full)
  if (length(zero)) ref <- lapply(ref, function(e) replaceSymbols(zero, rep("0", length(zero)), e))
  states <- reduced$states
  others <- setdiff(getSymbols(c(unlist(red), unlist(ref))), states)
  values <- setNames(as.list(1 + 0.37 * seq_along(others)), others)
  coef <- function(e) .linearCoefs(e, states, values)
  refCoef <- lapply(ref, coef)
  out <- list()
  for (k in seq_along(red)) {
    v <- coef(red[[k]])
    hit <- which(vapply(refCoef, function(w) !is.null(v) && !is.null(w) && any(w != 0) &&
                          qr(cbind(v, w))$rank == 1L, TRUE))
    hit <- setdiff(hit, match(names(out), names(ref)))
    if (length(hit)) out[[names(ref)[hit[1]]]] <- ref[[hit[1]]]
    else out[[.uniqueName(names(red)[k], c(names(out), names(ref)))]] <- red[[k]]
  }
  out
}


#' Detect states that are structurally zero in steady state
#'
#' Iterated three-layer test on the stoichiometric matrix of an `eqnlist`, the
#' one the steady-state solver of [steadyStates()] applies:
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
#' detection re-runs (removals can expose new zero-states). A state that some
#' reaction changes but none is left to change, outside every conserved
#' moiety, rests: the steady state leaves the value of such a state free.
#'
#' @param eqnlist_obj An [eqnlist].
#' @return `list(zero_states, resting, eqnlist)` with the reduced system.
#' @keywords internal
.zeroStatesFromSmatrix <- function(eqnlist_obj) {
  S0 <- eqnlist_obj$smatrix
  if (is.null(S0) || ncol(S0) == 0L || nrow(S0) == 0L)
    return(list(zero_states = character(0), resting = character(0), eqnlist = eqnlist_obj))

  S <- suppressWarnings(matrix(as.numeric(S0), nrow = nrow(S0), ncol = ncol(S0),
                               dimnames = dimnames(S0)))
  S[is.na(S)] <- 0
  # states some reaction changes; the others, such as catalysts, are constant
  changed <- colnames(S)[colSums(S != 0) > 0]
  # the conserved moieties, one column of states each
  moieties <- .conservedBasis(S)
  moiety <- if (is.null(moieties)) character(0) else rownames(moieties)[rowSums(moieties != 0) > 0]
  # whether a state is the last one of some moiety not yet at 0, which then holds the moiety total
  .last_holder <- function(st) {
    if (is.null(moieties) || !(st %in% rownames(moieties))) return(FALSE)
    for (k in which(moieties[st, ] != 0)) {
      others <- setdiff(rownames(moieties)[moieties[, k] != 0], c(st, zero_states))
      if (!length(others)) return(TRUE)
    }
    FALSE
  }
  rates       <- eqnlist_obj$rates
  description <- eqnlist_obj$description
  reactionCompartment <- eqnlist_obj$reactionCompartment

  zero_states <- character(0)

  # SBML rate rules may contain conditionals, so the structural probes below
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

  ## Outflux without influx => zero at SS. A state without any reaction left is
  ## constant instead, and so is the last state of a moiety not at 0, which
  ## holds the moiety's total.
  .neg_col <- function(M) {
    for (j in seq_len(ncol(M))) {
      if (any(M[, j] > 0) || !any(M[, j] < 0) || .last_holder(colnames(M)[j])) next
      return(j)
    }
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
      # Only species of the support reaching a leaking reaction along support
      # reactions are 0: educts of the leaking reactions, then the educts of
      # every reaction producing a species already found.
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
    return(list(zero_states = character(0), resting = character(0), eqnlist = eqnlist_obj))

  idle <- colnames(S)[colSums(S != 0) == 0]
  resting <- setdiff(intersect(idle, changed), moiety)
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

  list(zero_states = zero_states, resting = resting, eqnlist = new_obj)
}


#' Normalize steady-state inputs
#'
#' Coerces `trafo` to an `eqnvec`, zeroes and removes forcings, removes states
#' that are structurally zero at steady state and promotes states with a zero
#' right-hand side to parameters; the resting ones among them default to 0.
#'
#' @param trafo An [eqnlist], [eqnvec] or named character vector.
#' @param parameters Outer parameters (may be `NULL`).
#' @param forcings Forcing names.
#' @return A named list with `trafo`, `states`, `zero_states`, `resting`,
#'   `dependent`, `parameters`, `parms_all`, `totals` and, for an [eqnlist],
#'   `influx`, the sum of the producing fluxes per state.
#' @keywords internal
.normalize_ss_inputs <- function(trafo, parameters, forcings) {
  zero_states <- resting <- character(0); influx <- NULL
  original_params <- character(0)
  totals <- list()

  if (inherits(trafo, "eqnlist")) {
    original_params <- setdiff(getParameters(trafo), trafo$states)
    if (!is.null(forcings))
      trafo$rates <- replaceSymbols(forcings, rep("0", length(forcings)), trafo$rates)
    # a state no reaction changes is a parameter, not a total
    S <- trafo$smatrix; S[is.na(S)] <- 0
    held <- colnames(S)[colSums(abs(S)) == 0]
    full <- trafo
    zs <- .zeroStatesFromSmatrix(trafo)
    zero_states <- zs$zero_states
    resting     <- zs$resting
    trafo       <- zs$eqnlist
    if (!length(trafo$states))
      stop("All states are structurally zero in steady state; no dynamical ",
           "state remains to solve for. The network likely has irreversible ",
           "drains without matching influx (an open system with trivial ",
           "all-zero equilibrium).", call. = FALSE)
    totals <- .restingTotals(full, trafo, zero_states, c(held, resting))
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
  # a state in a total takes the state value from the total
  const_states <- setdiff(const_states, unlist(lapply(totals, getSymbols)))
  if (length(const_states))
    parameters <- union(parameters %||% character(0), const_states)

  dependent <- setdiff(states, parameters %||% character(0))
  if (!length(dependent))
    stop("No dependent states to solve for. All states appear in 'parameters'.",
         call. = FALSE)
  parameters <- Reduce(union, list(getSymbols(trafo[dependent], exclude = dependent),
                                   parameters %||% character(0),
                                   original_params, names(totals)))
  list(trafo = trafo, states = states, zero_states = zero_states, resting = resting,
       dependent = dependent, parameters = parameters,
       parms_all = setdiff(parameters, dependent),
       totals = totals,
       influx = if (!is.null(influx)) influx[intersect(dependent, names(influx))])
}

# Pimpl solver controls merged over their defaults. Unknown names are an error.
.pimplPTC <- function(controlsPTC) {
  def <- list(reltol = 1e-10, abstol = 1e-14, flowTol = 1, maxit = NULL, dtInit = 1e-2,
              positive = TRUE, stability = TRUE, stabilityTol = 0,
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

#' Implicit Parameter Transformation
#'
#' @description Solves \eqn{f(x, p) = 0} for the states \eqn{x} by pseudo-transient
#' continuation ([cppDE::ptc()]) and returns a [parfn] with sensitivities from
#' the implicit function theorem.
#'
#' @param trafo Named character, [eqnvec] or [eqnlist] defining \eqn{f}, in the
#'   language of [cppDE::expressions].
#' @param parameters Outer parameter names. For an [eqnlist] the totals
#'   \eqn{T} of the conserved quantities are added.
#'   States named here are held at their value instead of solved for, and a
#'   state whose right-hand side is `"0"` becomes a parameter.
#' @param forcings Forcing names, set to 0. Default `NULL`.
#' @param condition `NULL` (default) for a transformation valid in every
#'   condition, or the name of the condition it belongs to.
#' @param keepRoot If `TRUE` (default), roots are kept and repeated calls with
#'   the same values are answered from memory.
#' @param flow If `TRUE`, \eqn{\dot{x} = f(x, p)}{dx/dt = f(x, p)} and a stable
#'   steady state is returned. If `FALSE`, any regular root. Default `TRUE`
#'   for an [eqnlist], else `FALSE`.
#' @param compile Logical. `TRUE` compiles now, `FALSE` (default) leaves it
#'   to [compile()].
#' @param modelname Base name of the generated files, default `NULL`, which
#'   is `"impl_parfn"`, followed by `_<condition>`.
#' @param verbose Logical. Print the compiler output, default `FALSE`.
#' @param deriv,deriv2 Build first and second derivatives, default `TRUE` and
#'   `FALSE`. `deriv2` requires `deriv = TRUE`.
#' @param controlsPTC Named list of solver controls; an entry not given keeps
#'   the default, an unknown entry is an error:
#'   \describe{
#'     \item{`reltol`, `abstol`}{Convergence tolerances, default `1e-10`, `1e-14`.}
#'     \item{`flowTol`}{Relative local error per pseudo-time step, default `1`.}
#'     \item{`maxit`}{Iterations per start, default
#'       `ceiling(70 * log(n + 1))` for `n` states.}
#'     \item{`dtInit`}{Initial pseudo-time step relative to the fastest rate,
#'       default `1e-2`.}
#'     \item{`positive`}{Keep states positive, default `TRUE`.}
#'     \item{`stability`}{Require a stable root if `flow = TRUE`, default
#'       `TRUE`; where no start reaches one, the evaluation fails. With
#'       `FALSE` an unstable root is returned, on which a
#'       [constraintStability()] term rises.}
#'     \item{`stabilityTol`}{With `stability`, a root counts as stable while
#'       the real part of the leading eigenvalue on the manifold of the conserved
#'       quantities is at most `stabilityTol`, in the units of the rates; default
#'       `0`. A small positive value admits a resting state that departs from
#'       rest more slowly than the experiment lasts, as near a Hopf
#'       bifurcation.}
#'     \item{`nStarts`}{Random starts after the deterministic ones, default
#'       `20`.}
#'     \item{`startRange`}{Range of a random start, drawn uniformly per
#'       state on `startScale`, default `c(-5, 5)`.}
#'     \item{`startScale`}{`"log10"` or `"linear"`, default `"log10"` if
#'       `positive` and `"linear"` otherwise.}
#'     \item{`seed`}{Seed of the random starts, default `1`. The starts
#'       depend on the seed and the values; the global RNG is not used.}
#'   }
#' @param outdir Directory for the generated files, default
#'   `getOption("dMod2.outdir")`, else the working directory.
#'
#' @details Conserved quantities of an [eqnlist] enter as \eqn{C x = T}.
#' States without influx are 0. A state whose reactions all vanish with them
#' rests: it takes the value from `pars`, 0 where `pars` has none, and is not
#' among the parameters. The root is solved from the states in `pars`,
#' missing ones at 1, then from all states at 1, then from `nStarts` random
#' starts; with `flow = TRUE` the stable root the flow leads to. If all fail, it
#' is an error. A root depends on the parameters, the states in `pars` and the
#' controls only, not on the calls before it. Identical values are solved once,
#' also across conditions, and conditions are solved in parallel with the
#' cores of the evaluation. `keepRoot` and `controlsPTC` can be changed with
#' [controls()].
#'
#' @return A [parfn], called as `p(pars, ...)`, see section Calling a dMod
#'   function.
#' @inheritSection dModfn Calling a dMod function
#' @seealso [Pexpl()], [P()]; the guide chapter
#'   [Setting up models](https://dmodverse.github.io/dMod2/articles/guide-model.html).
#' @example inst/examples/Pimpl.R
#' @export
#' @import cppDE
#' @importFrom digest digest
Pimpl <- function(trafo, parameters = NULL, forcings = NULL, condition = NULL,
                  keepRoot = TRUE, flow = inherits(trafo, "eqnlist"),
                  compile = FALSE, modelname = NULL, verbose = FALSE,
                  deriv = TRUE, deriv2 = FALSE, controlsPTC = list(),
                  outdir = .dmodOutdir()) {

  flow    <- isTRUE(flow)
  emit_d1 <- isTRUE(deriv)
  emit_d2 <- isTRUE(deriv2)
  if (emit_d2 && !emit_d1)
    stop("Pimpl(deriv2 = TRUE) requires deriv = TRUE.", call. = FALSE)
  .pimplPTC(controlsPTC)

  norm <- .normalize_ss_inputs(trafo, parameters, forcings)
  zero_states <- norm$zero_states
  resting     <- norm$resting
  dependent   <- norm$dependent
  parameters  <- setdiff(norm$parameters, resting)
  n_dep       <- length(dependent)

  if (is.null(modelname)) modelname <- "impl_parfn"
  if (!is.null(condition)) modelname <- paste(modelname, sanitizeConditions(condition), sep = "_")

  # conservation rows C x - T, C possibly depending on parameters such as volumes; species that
  # are parameters stay in the expression
  tn <- names(norm$totals)
  cons <- if (length(tn))
    setNames(paste0("(", unlist(norm$totals), ") - ", tn), tn) else character(0)
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
  pinflux <- if (length(influx)) intersect(parms_all, getSymbols(unclass(influx))) else character(0)
  PIn <- if (length(influx)) suppressWarnings(cppDE::cppFUN(
    unclass(influx), variables = dependent, parameters = pinflux, fixed = NULL,
    compile = compile, modelname = paste0(modelname, "_influx"), outdir = outdir,
    verbose = verbose, convenient = FALSE, deriv = FALSE))

  # the native state, built at the first root after compile(): the roots of repeated values
  # (keepRoot) and the statistics
  native <- NULL
  ctl <- list(keepRoot = keepRoot, controlsPTC = controlsPTC)
  live <- function() !is.null(native) && pimpl_live(native)
  pushControls <- function() pimpl_controls(native, isTRUE(ctl$keepRoot), .pimplPTC(ctl$controlsPTC))
  nativeState <- function() {
    if (live()) return(native)
    entry <- get0(".ptcEntry", envir = asNamespace("cppDE"), inherits = FALSE)
    if (!is.function(entry))
      stop("Pimpl needs the steady-state solver of a newer cppDE; reinstall cppDE.", call. = FALSE)
    addr <- function(name) {
      a <- tryCatch(getNativeSymbolInfo(name)$address, error = function(e) NULL)
      if (is.null(a)) stop("Pimpl: compile '", modelname, "' first.", call. = FALSE)
      a
    }
    mn <- attr(PEval, "modelname")
    native <<- pimpl_state(list(
      val = addr(paste0(mn, "_eval")), ad = addr(paste0(mn, "_eval_ad")),
      ad2 = if (emit_d2) addr(paste0(mn, "_eval_ad2")),
      influx = if (!is.null(PIn)) addr(paste0(attr(PIn, "modelname"), "_eval")),
      ptc = entry(), nDep = n_dep, nPar = length(parms_all), nTotals = length(tn),
      inDep = match(names(influx), dependent) - 1L, inPar = match(pinflux, parms_all) - 1L,
      depNames = dependent, flow = flow))
    pushControls()
    native
  }
  # keepRoot and controlsPTC, see controls(); checked when set
  makeActiveBinding("controls", function(value) {
    if (missing(value)) return(ctl)
    .pimplPTC(value$controlsPTC)
    ctl <<- value
    if (live()) pushControls()
    invisible(value)
  }, environment())
  makeActiveBinding("stats", function() {
    if (!live()) return(list(calls = 0L, memo = 0L, solves = 0L, iter = 0, how = integer(0), failed = 0L))
    pimpl_stats(native)
  }, environment())
  statsReset <- function() if (live()) invisible(pimpl_stats(native, TRUE))

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
      dPsub <- .submatrix(dP, rows = colnames(jacobian))
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

    # a solved state and a varying input passed on keep their rows where the Jacobian vanishes here:
    # the second derivative may not, and an output without a row counts as fixed downstream
    varying <- if (is.null(dP)) input_cols else intersect(input_cols, rownames(dP))
    keep <- rownames(jacobian) %in% c(dependent, varying) | rowSums(jacobian != 0) > 0
    jacobian <- jacobian[keep, , drop = FALSE]
    if (!is.null(hessian)) hessian <- hessian[keep, , , drop = FALSE]
    list(jacobian = jacobian, hessian = hessian)
  }

  p2p <- function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE, condition = NULL) {
    if (deriv2 && !emit_d2)
      stop("Pimpl was built with deriv2 = FALSE; rebuild with deriv2 = TRUE.", call. = FALSE)
    if (!emit_d1) deriv <- FALSE
    if (deriv2 && !deriv) deriv <- TRUE

    p   <- pars
    dP  <- attr(p, "deriv")
    dP2 <- if (deriv2) attr(p, "deriv2") else NULL
    # a resting state not given rests at 0, without derivatives
    rest0 <- setdiff(resting, c(names(p), names(fixed)))
    if (length(rest0)) fixed <- c(fixed, setNames(rep(0, length(rest0)), rest0))
    if (!is.null(fixed)) {
      p <- p[!names(p) %in% names(fixed)]
      p <- c(p, fixed)
    }
    emptypars <- setdiff(names(p), c(dependent, names(fixed)))
    miss <- setdiff(dependent, names(p)); if (length(miss)) p[miss] <- 1
    pv <- setNames(as.numeric(p[parms_all]), parms_all)
    if (anyNA(pv))
      stop("Pimpl: missing value(s) for ", paste(parms_all[is.na(pv)], collapse = ", "), call. = FALSE)
    got <- pimpl_root(nativeState(), pv, as.numeric(p[dependent]), if (deriv2) 2L else if (deriv) 1L else 0L)
    for (w in got$warnings) warning(w, call. = FALSE)
    root <- setNames(got$x, dependent)

    zero_vec <- if (length(zero_states))
      setNames(rep(0, length(zero_states)), zero_states) else NULL
    out <- c(root, zero_vec, p[setdiff(names(p), c(dependent, zero_states))])

    np <- length(parms_all)
    d <- if (deriv && !is.null(got$D))
      build_derivs(matrix(got$D, n_dep, np, dimnames = list(dependent, parms_all)),
                   if (!is.null(got$D2)) array(got$D2, c(n_dep, np, np), list(dependent, parms_all, parms_all)),
                   out, p, emptypars, fixed, dP, dP2, deriv2)
    as.parvec(out,
              deriv  = if (deriv  && !is.null(d)) d$jacobian else NULL,
              deriv2 = if (deriv2 && !is.null(d)) d$hessian  else if (deriv2) NULL else FALSE)
  }

  attr(p2p, "vjpfn")       <- .parfnVjpFromJacobian(p2p)
  # the steady state as a stage of the native chain (.chainStage()): the unknowns, the states at
  # 0, the parameters and the native state that solves the roots
  attr(p2p, "chainImpl") <- function(condition) list(
    dependent = dependent, zero = zero_states, resting = resting, parameters = parms_all, d1 = emit_d1, d2 = emit_d2,
    state = nativeState())
  # the eigenvalues of the Jacobian at the steady state on the manifold of the conserved quantities,
  # with the gradients in `parms_all` of the real parts above `onset` (constraintStability())
  attr(p2p, "spectrum") <- function(pars, onset) {
    root <- p2p(pars, deriv = FALSE)
    p <- unclass(pars)
    p[setdiff(resting, names(p))] <- 0
    pv <- setNames(as.numeric(p[parms_all]), parms_all)
    if (anyNA(pv)) stop("Pimpl: missing value(s) for ", paste(parms_all[is.na(pv)], collapse = ", "), call. = FALSE)
    sp <- pimpl_spectrum(nativeState(), pv, as.numeric(root[dependent]), onset)
    colnames(sp$grad) <- parms_all
    sp
  }
  attr(p2p, "stabilityTol") <- function() .pimplPTC(ctl$controlsPTC)$stabilityTol
  attr(p2p, "equations")   <- as.eqnvec(all_exprs)
  attr(p2p, "parameters")  <- parameters
  attr(p2p, "modelname")   <- modelname
  attr(p2p, "compileInfo") <- .collectCompileInfo(PEval$func, PEval$jac, PEval$hess,
                                                  if (!is.null(PIn)) PIn$func)
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


#' Construct and Modify Parameter Transformations
#'
#' @description Build, change and branch the symbolic transformations that
#' [P()] turns into a [parfn].
#'
#' - `define()` sets the entry named by the left-hand side of `expr` to the
#'   right-hand side, adding it if missing.
#' - `insert()` replaces the left-hand side by the right-hand side wherever it
#'   occurs in the expressions.
#' - `branch()` copies `trafo` once per condition and keeps `table` as
#'   attribute `tree`.
#'
#' @param trafo Named character or [eqnvec], or a list of those named by
#'   condition as `branch()` returns it. `NULL` starts a new one.
#' @param expr `"lhs <- rhs"` string. The symbols in it are replaced by the
#'   arguments in `...` of the same name; see [repar()]. For `insert()` on a
#'   branched transformation also one or several effect formulas, e.g.
#'   `cl ~ wt + (1 | ID)`, which replace the symbol by a linear predictor
#'   over the columns of the branch table with fixed and random effects; see
#'   [randomEffects()]. Named arguments in `...` are then placeholders of the
#'   formulas, `insert(trafo, x ~ 1 + (1 | ID), x = c("ka", "cl"))`. A
#'   substitution or a second effect formula for a symbol that has an effect
#'   formula is an error: set the scale of a symbol before adding effects.
#' @param table `data.frame` with one row per condition, row names the
#'   conditions. Kept as attribute `tree`; the columns can be read in `...`
#'   by name.
#' @param conditions Condition names; default `rownames(table)`.
#' @param apply One of `"nothing"` (default), `"insert"` and `"define"`: what
#'   is done with every non-missing entry `table[condition, column]`, with
#'   the column name as left-hand side and the entry as right-hand side.
#'   Columns named `condition` or `conditions` are skipped.
#' @param conditionMatch Regular expression; only the conditions whose names
#'   match are changed. Default `NULL`, all. Needs the tree of a
#'   transformation built by `branch()`.
#' @param ... Named vectors replacing the symbols of `expr` of the same name,
#'   one substitution per element. They are evaluated per condition, with
#'   the columns of `attr(trafo, "tree")`, `.currentTrafo` (the
#'   transformation of the condition) and `.currentSymbols` (the symbols of
#'   that transformation) in
#'   scope. An unnamed logical argument of `insert()` selects the conditions
#'   where it is `TRUE`. Next to effect formulas they are more formulas and
#'   placeholders, evaluated once in the calling environment.
#'
#' @return Same shape as `trafo`, or a list named by condition for
#'   `branch()`.
#' @seealso [repar()], [P()]
#' @export
#' @example inst/examples/define.R
define <- function(trafo, expr, ..., conditionMatch = NULL) {
  if (missing(trafo)) trafo <- NULL
  dots <- substitute(alist(...))
  callerEnv <- parent.frame()
  re <- attr(trafo, "randomEffects", exact = TRUE)
  .reparBranches(trafo, dots, callerEnv, conditionMatch, function(tr, args) {
    .guardEffectSymbols(.substitutionSides(expr, args)$lhs, re)
    do.call(repar, c(list(expr = expr, trafo = tr, reset = TRUE), args))
  })
}


#' @export
#' @rdname define
insert <- function(trafo, expr, ..., conditionMatch = NULL) {
  if (missing(trafo)) trafo <- NULL
  dots <- substitute(alist(...))
  callerEnv <- parent.frame()
  # effect formulas, one or several; named dots are their placeholders
  if (inherits(expr, "formula")) {
    more <- lapply(eval(dots), eval, envir = callerEnv)
    nm <- names(more) %||% rep("", length(more))
    if (!all(vapply(more[!nzchar(nm)], inherits, TRUE, "formula")))
      stop("insert: next to effect formulas, pass more formulas or named placeholders.", call. = FALSE)
    formulas <- .expandEffectFormulas(c(list(expr), more[!nzchar(nm)]), more[nzchar(nm)])
    return(.insertEffects(trafo, formulas, conditionMatch))
  }
  re <- attr(trafo, "randomEffects", exact = TRUE)
  .reparBranches(trafo, dots, callerEnv, conditionMatch, function(tr, args) {
    ## Logical dots gate the substitution per condition rather than naming a
    ## symbol: any FALSE leaves the branch unchanged.
    gates <- vapply(args, is.logical, logical(1))
    if (!all(unlist(args[gates]))) return(tr)
    .guardEffectSymbols(.substitutionSides(expr, args[!gates])$lhs, re)
    do.call(repar, c(list(expr = expr, trafo = tr), args[!gates]))
  })
}


# `edit(branch, args)` for each branch of `trafo` (a single trafo is one
# branch) whose condition matches `conditionMatch`, with the dots evaluated
# against that condition. The result has the shape and tree of `trafo`.
.reparBranches <- function(trafo, dots, callerEnv, conditionMatch, edit) {
  tree <- attr(trafo, "tree")
  branched <- is.list(trafo)
  if (branched && is.null(names(trafo)))
    stop("If trafo is a list, elements must be named.", call. = FALSE)
  if (branched && !all(names(trafo) %in% rownames(tree)))
    stop("List names must be a subset of rownames(attr(trafo, 'tree')).", call. = FALSE)
  if (!is.null(conditionMatch) && is.null(tree))
    stop("conditionMatch needs the condition tree of a trafo built by branch().", call. = FALSE)

  branches <- if (branched) trafo else list(trafo)
  keys <- if (branched) names(trafo) else 1L
  out <- Map(function(tr, key) {
    cond <- tree[key, , drop = FALSE]
    if (!is.null(conditionMatch) && !str_detect(rownames(cond), conditionMatch))
      return(tr)
    symbols <- if (!is.null(tr)) getSymbols(tr)
    edit(tr, .evalDots(dots, cond, tr, symbols, callerEnv))
  }, branches, keys)
  if (!branched) out <- out[[1L]]
  attr(out, "tree") <- tree
  attr(out, "randomEffects") <- attr(trafo, "randomEffects", exact = TRUE)
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
      out[[cn]] <- if (apply == "insert") insert(single, paste0(par, " <- ", val))
                   else                    define(single, paste0(par, " <- ", val))
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
#' @description Replaces the left-hand side of `"lhs <- rhs"` by the
#' right-hand side in `trafo`, after the symbols of `expr` are replaced by the
#' arguments in `...`. [define()] and [insert()] build on it.
#'
#' @param trafo Named character, [eqnvec] or a list of those. Default `NULL`
#'   starts from the identity of the left-hand sides, so the result maps each
#'   left-hand side to the matching right-hand side.
#' @param expr `"lhs <- rhs"` string. Inside a name, `_` separates symbols,
#'   so `"Delta_x"` with `x = "a"` becomes `"Delta_a"`. A `"lhs ~ rhs"`
#'   string or a formula is an error: `~` is reserved for model formulas.
#' @param ... Named character or numeric vectors of equal length, replacing
#'   the symbols of `expr` of the same name; element `i` of each gives the
#'   `i`-th substitution.
#' @param reset `FALSE` (default) replaces the left-hand side wherever it
#'   occurs in the expressions, `TRUE` sets the entry it names.
#' @return Same shape as `trafo`.
#' @seealso [define()], [insert()]
#' @export
#' @examples
#' innerpars   <- letters[1:3]
#' constraints <- c(a = "b + c")
#' mycondition <- "cond1"
#' trafo <- repar(NULL,  "x <- x",      x = innerpars)
#' trafo <- repar(trafo, "x <- y",      x = names(constraints), y = constraints)
#' trafo <- repar(trafo, "x <- exp(x)", x = innerpars)
#' trafo <- repar(trafo, "x <- x + Delta_x_condition",
#'                x = innerpars, condition = mycondition)
repar <- function(trafo = NULL, expr, ..., reset = FALSE) {
  if (.isReparExpr(trafo) && (missing(expr) || !.isReparExpr(expr)))
    stop("repar: the order is repar(trafo, expr).", call. = FALSE)
  sides <- .substitutionSides(expr, list(...))
  lhs <- sides$lhs; rhs <- sides$rhs

  if (is.null(trafo)) return(as.eqnvec(structure(rhs, names = lhs)))
  if (is.list(trafo)      && !reset)       lapply(trafo, function(t) replaceSymbols(lhs, rhs, t))
  else if (is.character(trafo) && !reset)       replaceSymbols(lhs, rhs, trafo)
  else if (is.list(trafo)      &&  reset)       lapply(trafo, function(t) { t[lhs] <- rhs; t })
  else { trafo[lhs] <- rhs; trafo }
}

# A substitution string, or a formula or `~` string given in place of one, as
# opposed to a transformation.
.isReparExpr <- function(x)
  inherits(x, "formula") ||
    (is.character(x) && length(x) == 1L && is.null(names(x)) &&
       grepl("<-|~", x))

# Left- and right-hand sides of a substitution, one pair per element of the
# placeholder vectors in `args`.
.substitutionSides <- function(expr, args) {
  parsed <- .parseSubstitution(expr)
  args <- lapply(args, as.character)
  if (!length(args)) return(list(lhs = .decolonize(parsed[1]), rhs = .decolonize(parsed[2])))
  reps <- as.data.frame(args, stringsAsFactors = FALSE)
  side <- function(x) vapply(seq_len(nrow(reps)), function(i)
    .decolonize(replaceSymbols(colnames(reps), reps[i, ], x)), character(1))
  list(lhs = side(parsed[1]), rhs = side(parsed[2]))
}

# Left- and right-hand side of a `"lhs <- rhs"` string, colonized.
.parseSubstitution <- function(expr) {
  if (inherits(expr, "formula") || (is.character(expr) && !grepl("<-", expr, fixed = TRUE) &&
                                      grepl("~", expr, fixed = TRUE)))
    stop("Substitutions are written \"lhs <- rhs\", e.g. \"x <- exp(x)\"; ",
         "'~' is reserved for model formulas.", call. = FALSE)
  if (!is.character(expr) || length(expr) != 1L)
    stop("A substitution is one string \"lhs <- rhs\".", call. = FALSE)
  e <- str2lang(.colonize(expr))
  if (!is.call(e) || !identical(e[[1]], as.name("<-")) || (is.call(e[[3]]) && identical(e[[3]][[1]], as.name("<-"))))
    stop("Substitutions are written \"lhs <- rhs\", got \"", expr, "\".", call. = FALSE)
  c(deparse1(e[[2]]), deparse1(e[[3]]))
}

paste_ <- function(...) paste(..., sep = "_")