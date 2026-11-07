#' Stability of the steady state of Pimpl() as a constraint
#'
#' @description A term of the objective that keeps the steady state of a
#'   [Pimpl()] stable. For every eigenvalue \eqn{\lambda} of the Jacobian of the
#'   right-hand side at the steady state, on the manifold of the conserved
#'   quantities, it adds \eqn{w\, r^2} with
#'   \eqn{r = \max(0, \mathrm{Re}\,\lambda - \lambda_0) / (\lambda_{tol} - \lambda_0)}:
#'   nothing while every real part lies below the onset \eqn{\lambda_0}, and the
#'   weight \eqn{w} for each eigenvalue, both of a complex pair, whose real part
#'   reaches \eqn{\lambda_{tol}}, the `stabilityTol` of `pSS`. Steady states
#'   that agree in all inputs count once. The gradient is exact, the Hessian
#'   Gauss-Newton.
#'
#' @details The term is an assumption of the model, that the system rests where
#'   the steady state is taken, not data: [profile()] with `stop = "data"` and
#'   likelihood ratios read the data term without it. Parameters of `p` that the
#'   evaluation does not give, such as random effects, enter at 0.
#'
#'   Where `pSS` requires a stable root, an evaluation beyond `stabilityTol`
#'   fails. With `Pimpl(controlsPTC = list(stability = FALSE))` it returns the
#'   unstable root there, and the term goes on rising.
#'
#' @param pSS A parfn of [Pimpl()], built with `deriv2 = TRUE`.
#' @param p The parfn whose output `pSS` reads, or `NULL` where `pSS` reads the
#'   parameters of the objective.
#' @param onset \eqn{\lambda_0}, below `stabilityTol`. `NULL` (default) for
#'   `stabilityTol - 1 / duration`: the term begins where a perturbation of the
#'   steady state decays more slowly than the data last.
#' @param duration The time the data span, for the default `onset`.
#' @param weight \eqn{w}, the value of the term where the real part of an
#'   eigenvalue reaches `stabilityTol`, in units of \eqn{-2 \log L}. Default 100.
#' @param attrName Name of the attribute that holds the value of the term,
#'   default `"stability"`.
#' @return An objfn of the parameters of `p`, or of `pSS` without `p`, to add to
#'   the objective.
#' @seealso [Pimpl()], [constraintL2()]
#' @example inst/examples/constraintStability.R
#' @export
constraintStability <- function(pSS, p = NULL, onset = NULL, duration = NULL, weight = 100,
                                attrName = "stability") {
  kernels <- attr(pSS, "mappings")
  if (!length(kernels) || !all(vapply(kernels, function(k) is.function(attr(k, "spectrum")), TRUE)))
    stop("constraintStability: `pSS` is a parfn of Pimpl().", call. = FALSE)
  if (is.null(onset) && is.null(duration))
    stop("constraintStability: give `onset`, or `duration`, the time the data span.", call. = FALSE)
  if (!is.null(duration) && !(is.numeric(duration) && length(duration) == 1L && duration > 0))
    stop("constraintStability: `duration` is a positive number.", call. = FALSE)
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
  kernelOf <- function(cond) if (!is.null(cond) && !is.null(kernels[[cond]])) kernels[[cond]] else kernels[[1]]

  myfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, hessian = NULL,
                   conditions = NULL, env = NULL, cores = getOption("dMod2.cores", 1L)) {
    args <- list(...)
    pars <- args[[match.fnargs(args, "pars")]]
    free <- names(pars)
    tol <- attr(kernels[[1]], "stabilityTol")()
    lam0 <- if (is.null(onset)) tol - 1 / duration else onset
    if (!(is.numeric(lam0) && length(lam0) == 1L && lam0 < tol))
      stop("constraintStability: `onset` must be below the stabilityTol of `pSS` (", tol, ").", call. = FALSE)
    span <- tol - lam0
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
    seen <- character(0)
    for (i in seq_along(inner)) {
      x <- inner[[i]]
      key <- paste(format(x[intersect(names(x), inputs)], digits = 15), collapse = ",")
      if (key %in% seen) next
      seen <- c(seen, key)
      sp <- attr(kernelOf(conds[[i]]), "spectrum")(x, lam0)
      if (!length(sp$active)) next
      J <- jac[[i]]
      on <- intersect(colnames(sp$grad), rownames(J))
      cols <- intersect(free, colnames(J))
      for (k in seq_along(sp$active)) {
        r <- (sp$re[sp$active[k]] - lam0) / span
        g <- setNames(numeric(length(free)), free)
        if (length(on) && length(cols))
          g[cols] <- drop(sp$grad[k, on, drop = FALSE] %*% J[on, cols, drop = FALSE]) / span
        value <- value + weight * r^2
        gradient <- gradient + 2 * weight * r * g
        H <- H + 2 * weight * tcrossprod(g)
      }
    }
    # gradient and Hessian only where asked for, as the terms they are summed with
    cv <- .resolveCurvature(deriv, deriv2, hessian, "forward")
    res <- if (!isTRUE(deriv)) objlist(value = value, gradient = NULL, hessian = NULL)
           else if (cv$sparse) .sparseLift(value, gradient, H, free)
           else objlist(value = value, gradient = gradient, hessian = if (cv$hessian) H)
    attr(res, attrName) <- value
    attr(res, "env") <- env
    res
  }
  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- NULL
  attr(myfn, "parameters") <- parameters
  myfn
}
