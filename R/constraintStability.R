#' Stability of the steady state of Pimpl() as a constraint
#'
#' @description A term of the objective that keeps the steady state of a
#'   [Pimpl()] stable. For every eigenvalue \eqn{\lambda} of the Jacobian of the
#'   right-hand side at the steady state, on the manifold of the conserved
#'   quantities, it adds the square of the residual
#'   \eqn{\sqrt{w}\, r^{3/2}} for \eqn{r \le 1} and
#'   \eqn{\sqrt{w}\, (1 + 3 (r - 1) / 2)} beyond, with
#'   \eqn{r = \max(0, \mathrm{Re}\,\lambda - \lambda_0) / (\lambda_{tol} - \lambda_0)}.
#'   With \eqn{s} the largest eigenvalue modulus there, \eqn{\lambda_{tol}} the
#'   `stabilityTol` and \eqn{m} the `stabilityMargin` of the `controlsPTC` of
#'   `pSS`, the ramp begins at \eqn{\lambda_0 = \lambda_{tol} - m s}. The term is
#'   zero on a steady state whose leading real part lies more than \eqn{m s}
#'   below \eqn{\lambda_{tol}}, reaches the weight \eqn{w} for each eigenvalue, both of
#'   a complex pair, whose real part reaches \eqn{\lambda_{tol}}, and grows
#'   quadratically beyond. Value, gradient and Gauss-Newton Hessian are
#'   continuous. Steady states that agree in the parameters, resting states and
#'   start values they read count once. The gradient is exact.
#'
#' @details The term is an assumption of the model, that the system rests where
#'   the steady state is taken, not data: [profile()] with `stop = "data"` and
#'   likelihood ratios read the data term without it. Parameters of `p` that the
#'   evaluation does not give, such as random effects, enter at 0.
#'
#'   Where `pSS` requires a stable root, an evaluation beyond the stability
#'   tolerance fails. With `Pimpl(controlsPTC = list(stability = FALSE))` it returns the
#'   unstable root there, and the term goes on rising.
#'
#'   The attribute `"eigenvalues"` of the result lists the eigenvalues above the
#'   onset, one row each: the conditions of the steady state, real and imaginary
#'   part, the half-life or doubling time \eqn{\log 2 / |\mathrm{Re}\,\lambda|}, the
#'   onset, the contribution to the term, and the states with the largest
#'   entries of the right eigenvector, scaled to a largest entry of 1. The sum
#'   with other terms keeps the value only, so the attribute is read from a call
#'   of the term itself.
#'
#' @param pSS A parfn of [Pimpl()], built with `deriv2 = TRUE`.
#' @param p The parfn whose output `pSS` reads, or `NULL` where `pSS` reads the
#'   parameters of the objective.
#' @param weight \eqn{w}, the value of the term where the real part of an
#'   eigenvalue reaches the stability tolerance, in units of \eqn{-2 \log L}.
#'   Default 100.
#' @param attrName Name of the attribute that holds the value of the term,
#'   default `"stability"`.
#' @return An objfn of the parameters of `p`, or of `pSS` without `p`, to add to
#'   the objective.
#' @seealso [Pimpl()], [constraintL2()]
#' @example inst/examples/constraintStability.R
#' @export
constraintStability <- function(pSS, p = NULL, weight = 100, attrName = "stability") {
  kernels <- attr(pSS, "mappings")
  if (!length(kernels) || !all(vapply(kernels, function(k) is.function(attr(k, "spectrum")), TRUE)))
    stop("constraintStability: `pSS` is a parfn of Pimpl().", call. = FALSE)
  if (!(is.numeric(weight) && length(weight) == 1L && weight > 0))
    stop("constraintStability: `weight` is a positive number.", call. = FALSE)
  inputs <- getParameters(pSS)
  # the parameters of p that enter the inputs of the steady state
  parameters <- if (is.null(p)) inputs else {
    eqs <- Filter(Negate(is.null), getEquations(p))
    used <- unique(unlist(lapply(eqs, function(e) getSymbols(unclass(e)[intersect(names(e), inputs)]))))
    if (length(eqs)) intersect(getParameters(p), c(used, setdiff(inputs, unlist(lapply(eqs, names)))))
    else getParameters(p)
  }
  kernelAt <- function(cond) if (!is.null(cond) && !is.null(kernels[[cond]])) cond else 1L
  # what the root and the Jacobian at the root read: parameters, resting states, start values;
  # known once the steady state is compiled
  reads <- NULL

  myfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, hessian = NULL,
                   conditions = NULL, env = NULL, cores = getOption("dMod2.cores", 1L)) {
    args <- list(...)
    pars <- args[[match.fnargs(args, "pars")]]
    free <- names(pars)
    tol <- attr(kernels[[1]], "stabilityTol")()
    lamTol <- tol[["stabilityTol"]]
    margin <- tol[["stabilityMargin"]]
    if (is.null(p)) {
      inner <- list(c(unclass(pars), unclass(fixed)))
      jac <- list(matrix(diag(1, length(free)), length(free), dimnames = list(free, free)))
      conds <- list(NULL)
    } else {
      missing <- setdiff(getParameters(p), c(free, names(fixed)))
      fx <- c(fixed, setNames(rep(0, length(missing)), missing))
      out <- p(pars, fixed = fx, deriv = TRUE)
      inner <- lapply(out, unclass)
      jac <- lapply(out, function(o) attr(o, "deriv", exact = TRUE))
      conds <- if (is.null(names(out))) vector("list", length(out)) else as.list(names(out))
    }
    value <- 0
    gradient <- setNames(numeric(length(free)), free)
    H <- matrix(0, length(free), length(free), dimnames = list(free, free))
    if (is.null(reads))
      reads <<- lapply(kernels, function(k) unlist(attr(k, "chainImpl")(NULL)[c("parameters", "resting", "dependent")]))
    at <- lapply(conds, kernelAt)
    keys <- vapply(seq_along(inner), function(i) {
      x <- inner[[i]]
      paste(c(at[[i]], format(x[intersect(names(x), reads[[at[[i]]]])], digits = 15)), collapse = ",")
    }, "")
    rows <- list()
    for (i in which(!duplicated(keys))) {
      sp <- attr(kernels[[at[[i]]]], "spectrum")(inner[[i]], lamTol, margin)
      if (!length(sp$active)) next
      J <- jac[[i]]
      cols <- intersect(free, colnames(J))
      # the gradient in `free` of a row vector in the parameters of pSS
      chain <- function(v) {
        g <- setNames(numeric(length(free)), free)
        on <- intersect(names(v), rownames(J))
        if (length(on) && length(cols)) g[cols] <- drop(v[on] %*% J[on, cols, drop = FALSE])
        g
      }
      # the ramp ends at the tolerance and is `margin` times the largest eigenvalue modulus wide
      ds <- chain(sp$scaleGrad)
      span <- margin * sp$scale
      lam0 <- lamTol - span
      for (k in seq_along(sp$active)) {
        re <- sp$re[sp$active[k]]
        r <- (re - lam0) / span
        dr <- (chain(sp$grad[k, ]) + (1 - r) * margin * ds) / span
        # residual sqrt(w) r^(3/2) up to the tolerance, then linear with the same slope
        rho <- sqrt(weight) * if (r <= 1) r^1.5 else 1 + 1.5 * (r - 1)
        drho <- sqrt(weight) * 1.5 * sqrt(min(r, 1)) * dr
        value <- value + rho^2
        gradient <- gradient + 2 * rho * drho
        H <- H + 2 * tcrossprod(drho)
        u <- sp$vectors[, k]
        top <- head(order(u, decreasing = TRUE), 3L)
        rows[[length(rows) + 1L]] <- data.frame(
          condition = paste(unlist(conds[keys == keys[i]]), collapse = ", "),
          re = re, im = sp$im[sp$active[k]], halfLife = log(2) / abs(re), onset = lam0,
          penalty = rho^2,
          states = paste0(names(u)[top], " (", formatC(u[top], digits = 2, format = "f"), ")", collapse = ", "))
      }
    }
    # gradient and Hessian only where asked for, as the terms they are summed with
    cv <- .resolveCurvature(deriv, deriv2, hessian, "forward")
    res <- if (!isTRUE(deriv)) objlist(value = value, gradient = NULL, hessian = NULL)
           else if (cv$sparse) .sparseLift(value, gradient, H, free)
           else objlist(value = value, gradient = gradient, hessian = if (cv$hessian) H)
    attr(res, attrName) <- value
    attr(res, "eigenvalues") <- if (length(rows)) do.call(rbind, rows)
      else data.frame(condition = character(0), re = numeric(0), im = numeric(0), halfLife = numeric(0),
                      onset = numeric(0), penalty = numeric(0), states = character(0))
    attr(res, "env") <- env
    res
  }
  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- NULL
  attr(myfn, "parameters") <- parameters
  myfn
}
