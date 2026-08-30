## Penalty specification for regularised NLME fits (.fitLaplace).
##
## A `penaltySpec` is the L1/Laplace counterpart of an `omegaSpec`
## (see nlmeClass.R): it declares which model parameters are candidates for
## being individual-specific and how their per-subject deviations are
## penalised, but carries a single regularisation strength `lambda` instead of
## a Cholesky-parametrised covariance. It shares the structural fields of
## `omegaSpec` (`eta`, `short`, `K`, `subjects`, `subjectEtas`) so the two can
## eventually be unified.
##
## One constructor, penaltyL1(pars, method = ...), with three methods:
##   method = "lasso"     : element-wise L1 of the named deviations toward a
##                          target (default 0). Factorises over subjects.
##   method = "fused"     : star, deviations toward one shared estimated centre
##                          (equivalently eta toward 0). Factorises.
##   method = "clustered" : complete-graph pairwise |eta_i - eta_j|; individuals
##                          crystallise into clusters. Couples subjects.
## Specs combine with `+` (one shared lambda, blocks concatenated).
##
## The eta-name convention matches omega():
##   base parameter "Ka" -> eta name "eta_Ka" -> subject eta "eta_Ka_<subject>".


## Resolve the penalty-strength names for one block's parameters. `lambda` is
## either a single name (shared by all `pars`; per-parameter names when
## `perParam = TRUE`, i.e. "<lambda>_<par>") or a named character vector giving
## one name per parameter explicitly. Returns a named character vector
## (par -> lambda name).
.penaltyLambdaMap <- function(lambda, pars, perParam = FALSE) {
  if (!is.character(lambda) || length(lambda) < 1L || any(!nzchar(lambda)))
    stop("`lambda` must be a non-empty character vector of parameter name(s).")
  if (length(lambda) > 1L || !is.null(names(lambda))) {
    if (is.null(names(lambda)) || !setequal(names(lambda), pars))
      stop("A vector-valued `lambda` must be named to match `pars`.")
    return(lambda[pars])
  }
  if (isTRUE(perParam))
    return(stats::setNames(paste(lambda, pars, sep = "_"), pars))
  stats::setNames(rep(lambda, length(pars)), pars)
}


## Internal constructor. `blocks` is a list of penalty blocks, each a list with
## fields `kind` ("single"/"fused"/"clustered"), `params` (base parameter
## names), `lambda` (named char vector par -> penalty-strength name) and, for
## "single", `target` (named numeric aligned to `params`). A single global
## strength (the default) gives one `lambda` name shared by every param; a
## per-parameter penalty gives one name per param (soft sparsity in the fit).
.newPenaltySpec <- function(blocks, subjects = NULL, prefix = "eta") {

  if (!is.list(blocks) || length(blocks) < 1L)
    stop("`blocks` must be a non-empty list of penalty blocks.")

  params_by_block <- lapply(blocks, `[[`, "params")
  pars_all <- unlist(params_by_block, use.names = FALSE)
  if (anyDuplicated(pars_all))
    stop("A parameter appears in more than one penalty block: ",
         paste(unique(pars_all[duplicated(pars_all)]), collapse = ", "),
         ". Each candidate parameter may carry only one penalty.")

  ## per-parameter strength names, then the distinct names (order preserved)
  lambda_by_param <- unlist(lapply(blocks, `[[`, "lambda"))
  lambda_by_param <- lambda_by_param[pars_all]
  if (is.null(lambda_by_param) || any(!nzchar(lambda_by_param)))
    stop("Every penalty block must carry non-empty `lambda` name(s).")
  lambda_names <- unique(unname(lambda_by_param))

  eta   <- paste0(prefix, "_", pars_all)
  short <- pars_all
  K     <- length(pars_all)

  subject_etas <- NULL
  if (!is.null(subjects)) {
    if (!is.character(subjects) || anyDuplicated(subjects) ||
        length(subjects) < 1L)
      stop("`subjects` must be a non-empty character vector with unique entries.")
    subject_etas <- outer(subjects, eta, function(s, e) paste(e, s, sep = "_"))
    dimnames(subject_etas) <- list(subjects, eta)
  }

  out <- list(
    blocks        = blocks,
    params        = pars_all,
    eta           = eta,
    short         = short,
    K             = K,
    subjects      = subjects,
    subjectEtas   = subject_etas,
    lambdaName    = lambda_names,       # distinct strength names (length 1 = global)
    lambdaByParam = lambda_by_param,    # par -> its strength name
    prefix        = prefix
  )
  class(out) <- c("penaltyspec", "list")
  out
}


## Expand a scalar-or-named `target` to a named numeric aligned to `pars`.
.expandPenaltyTarget <- function(target, pars) {
  if (length(target) == 1L && is.null(names(target)))
    return(stats::setNames(rep(as.numeric(target), length(pars)), pars))
  if (is.null(names(target))) {
    if (length(target) != length(pars))
      stop("Unnamed `target` must have length 1 or length(pars).")
    return(stats::setNames(as.numeric(target), pars))
  }
  if (!setequal(names(target), pars))
    stop("`target` names must match `pars`.")
  target[pars]
}


#' L1 penalty for a mixed-effects fit
#'
#' @description
#' Declares an L1/Laplace penalty on the per-subject deviations
#' `eta_<par>_<subject>` of the named parameters: the sparsity/clustering
#' counterpart of [omega], which declares a Gaussian random effect. In the
#' mixed-effects reading the deviations are random effects drawn from a Laplace
#' density whose single strength `lambda` is not fixed but ESTIMATED by marginal
#' maximum likelihood (one [EM] fit, no penalty-strength scan). A deviation
#' driven to zero says the subject shares the population value of that parameter,
#' so the estimated support of the deviations is the parsimonious "which
#' parameters are individual-specific" model. Consumed by [constraintL1] (which
#' stamps the spec onto a `normL2` objective) and by [EM]; combine several
#' penalty terms with `+`, and turn the fit into a hard structural model with
#' [sparsify] (the marginal likelihood supplies the Occam factor a classical
#' penalty scan needs a lambda grid plus an information criterion to mimic).
#'
#' The `method` chooses which relationship between subjects is penalised, i.e.
#' which question the fit answers. All three are members of the (generalised)
#' fused-lasso family:
#' \describe{
#'   \item{`"lasso"` (default)}{Element-wise L1 of each deviation toward a fixed
#'     `target` (default 0), \eqn{\lambda \sum_i |\eta_i - t|}. The plain lasso:
#'     it shrinks individual deviations to the target and answers "is this
#'     subject different from the reference value?". It factorises over subjects
#'     (the cheap, per-subject path) and also drives the classical two-group
#'     reference encoding, where one subject is the baseline (no deviation) and
#'     the other carries the single penalised difference. Works for any number of
#'     subjects.}
#'   \item{`"fused"`}{Star penalty pulling every deviation toward one shared,
#'     estimated centre (equivalently `eta` toward 0, with the centre absorbed
#'     into the structural mean). It answers "is this subject an outlier from the
#'     population mean?", the direct L1 analog of a Gaussian random effect. It
#'     factorises, but needs every subject to carry an anchored deviation (no
#'     reference baseline).}
#'   \item{`"clustered"`}{Complete-graph pairwise penalty
#'     \eqn{\sum_{i<j} |\eta_i - \eta_j|}. It answers the richer "how do the
#'     subjects group?": individuals crystallise into an unknown NUMBER of
#'     clusters that share a value, with all-shared (one cluster) and
#'     all-distinct (n clusters) as the two extremes, so the shared vs individual
#'     question is just its two-cluster special case. It COUPLES all
#'     subjects (the penalty is not separable), needs every subject to carry an
#'     anchored deviation, is meaningful from three subjects up, and is scored by
#'     the exact clustered marginal in [sparsify]. Higher cost than the
#'     factorising methods.}
#' }
#'
#' @param pars Character vector of base parameter names that are candidates for
#'   being individual-specific, e.g. `c("Ka", "Cl")`. Their per-subject
#'   deviations `eta_<par>_<subject>` are penalised.
#' @param method One of `"lasso"` (element-wise, default), `"fused"` (star
#'   toward a shared centre), or `"clustered"` (pairwise, groups emerge).
#' @param target Scalar or named numeric giving the reference value each
#'   deviation is pulled toward. Default `0`. Used only for `method = "lasso"`.
#' @param perParam Logical; if `TRUE`, give each parameter its own penalty
#'   strength `lambda_<par>` (estimated separately by marginal ML: the route to
#'   soft sparsity in the fit, where shared parameters acquire a large strength
#'   that shrinks their deviations, individual ones a small strength). Default
#'   `FALSE` (a single global strength `lambda`). Not supported for
#'   `method = "clustered"`.
#' @param subjects Optional character vector of subject identifiers. Required
#'   before the spec is fed to [EM]; optional standalone.
#'
#' @return A `penaltySpec` (S3 list) with components `blocks` (the penalty
#'   blocks), `params`, `eta`, `short`, `K`, `subjects`, `subjectEtas`,
#'   `lambdaName` (distinct strength names), `lambdaByParam`, `prefix`.
#' @seealso [constraintL1], [EM], [sparsify], [omega]
#' @examples
#' pen <- penaltyL1(c("Ka", "Cl"), subjects = paste0("s", 1:3))
#' pen$eta
#' gp <- penaltyL1(c("kon", "koff"), method = "clustered",
#'                 subjects = paste0("s", 1:4))
#' gp$blocks[[1]]$kind
#' @export
penaltyL1 <- function(pars, method = c("lasso", "fused", "clustered"),
                      target = 0, perParam = FALSE, subjects = NULL) {
  method <- match.arg(method)
  if (!is.character(pars) || length(pars) < 1L)
    stop("`pars` must be a non-empty character vector of parameter names.")
  if (anyDuplicated(pars))
    stop("`pars` names must be unique.")
  if (identical(method, "clustered") && isTRUE(perParam))
    stop("penaltyL1(method = \"clustered\") uses a single global `lambda`; ",
         "per-parameter strengths are not supported on the coupled path.")
  kind  <- switch(method, lasso = "single", fused = "fused",
                  clustered = "clustered")
  block <- list(kind   = kind,
                params = pars,
                lambda = .penaltyLambdaMap("lambda", pars, perParam))
  if (identical(method, "lasso"))
    block$target <- .expandPenaltyTarget(target, pars)
  .newPenaltySpec(list(block), subjects = subjects, prefix = "eta")
}


## Merge two penalty specs into one (concatenated blocks; each block carries its
## own penalty-strength name(s), so the merged spec is the union - a single
## global lambda when both share it, distinct names otherwise).
.mergePenaltySpec <- function(s1, s2) {
  if (is.null(s1)) return(s2)
  if (is.null(s2)) return(s1)
  if (!inherits(s1, "penaltyspec") || !inherits(s2, "penaltyspec"))
    stop("Both operands must be penaltySpec objects.")
  if (s1$prefix != s2$prefix)
    stop("Cannot combine penalties with different `prefix` (",
         s1$prefix, " vs ", s2$prefix, ").")
  subj1 <- s1$subjects; subj2 <- s2$subjects
  if (!is.null(subj1) && !is.null(subj2) && !setequal(subj1, subj2))
    stop("Cannot combine penalties defined on different subject sets.")
  subjects <- if (is.null(subj1)) subj2 else subj1
  .newPenaltySpec(c(s1$blocks, s2$blocks),
                  subjects = subjects, prefix = s1$prefix)
}

#' Combine two penalty specifications
#'
#' @param e1,e2 `penaltySpec` objects sharing one `lambda`, `prefix` and subject
#'   set. Their penalty blocks are concatenated; a parameter may appear in only
#'   one block.
#' @return A merged `penaltySpec`.
#' @export
"+.penaltyspec" <- function(e1, e2) .mergePenaltySpec(e1, e2)


#' @export
c.penaltyspec <- function(...) {
  specs <- list(...)
  Reduce(.mergePenaltySpec, specs)
}


#' Print method for a penalty specification
#'
#' @param x A `penaltySpec` object.
#' @param ... Ignored.
#' @export
print.penaltyspec <- function(x, ...) {
  cat("penaltySpec\n")
  cat(sprintf("  K       : %d\n", x$K))
  if (length(x$lambdaName) == 1L) {
    cat(sprintf("  lambda  : %s (single global strength)\n", x$lambdaName))
  } else {
    cat(sprintf("  lambda  : %s (per-parameter strengths)\n",
                paste(x$lambdaName, collapse = ", ")))
  }
  for (b in x$blocks) {
    meth <- switch(b$kind, single = "lasso", b$kind)   # display the method name
    tgt  <- if (identical(b$kind, "single"))
      sprintf(" -> target %s", paste(format(b$target), collapse = ", ")) else ""
    cat(sprintf("  method  : %-9s [%s]%s\n",
                meth, paste(b$params, collapse = ", "), tgt))
  }
  if (!is.null(x$subjectEtas)) {
    cat(sprintf("  subjects: %d (%s)\n",
                nrow(x$subjectEtas),
                paste(rownames(x$subjectEtas), collapse = ", ")))
  } else {
    cat("  subjects: none (set subjects = ... before EM)\n")
  }
  invisible(x)
}


#' Parameter names of a penalty specification
#'
#' @param x A `penaltySpec` object.
#' @param what One of `"all"` (subject etas plus the lambda name), `"eta"`
#'   (subject-level deviation parameters), or `"lambda"` (the strength name).
#' @param ... Ignored.
#' @return Character vector of parameter names.
#' @export
parnames.penaltyspec <- function(x, what = c("all", "eta", "lambda"), ...) {
  what <- match.arg(what)
  eta_names <- if (is.null(x$subjectEtas)) character(0) else as.vector(x$subjectEtas)
  switch(what,
         all    = c(eta_names, x$lambdaName),
         eta    = eta_names,
         lambda = x$lambdaName)
}


## Expand a penaltySpec into a flat list of penalised terms. Each term is a
## list(i, j, target, lambda): the penalty contributes lambda * |value(i) - ref|
## where ref = value(j) for a difference term (clustered, j an eta name) or
## `target` for an element-wise term (single/fused, j = NA), and `lambda` is the
## penalty-strength name for that term's parameter (one shared name for a global
## penalty, a per-parameter name otherwise). Requires subject expansion.
.penaltyTerms <- function(spec) {
  se <- spec$subjectEtas
  if (is.null(se))
    stop("penaltySpec has no subject expansion; set subjects = ... .")
  subjects <- rownames(se)
  n <- length(subjects)
  terms <- list()
  push <- function(i, j, target, lambda)
    terms[[length(terms) + 1L]] <<-
      list(i = i, j = j, target = target, lambda = lambda)

  for (b in spec$blocks) {
    for (par in b$params) {
      col    <- paste0(spec$prefix, "_", par)
      etas_k <- se[, col]                       # named by subjects
      lam_k  <- spec$lambdaByParam[[par]]       # this parameter's strength name
      if (b$kind %in% c("single", "fused")) {
        tgt <- if (identical(b$kind, "fused")) 0 else as.numeric(b$target[[par]])
        for (i in seq_len(n))
          push(etas_k[[i]], NA_character_, tgt, lam_k)
      } else if (identical(b$kind, "clustered")) {
        if (n >= 2L) for (a in seq_len(n - 1L)) for (cc in (a + 1L):n)
          push(etas_k[[a]], etas_k[[cc]], 0, lam_k)
      } else {
        stop("Unknown penalty block kind: ", b$kind)
      }
    }
  }
  terms
}


#' L1 / Laplace regularisation constraint over subject-level deviations
#'
#' @description
#' The L1 counterpart of [constraintL2] with an `Omega`. Given a [penaltyL1]
#' specification it contributes the penalty
#' \eqn{\lambda \sum_t |(D\eta)_t|}, where the terms \eqn{(D\eta)_t} are the
#' element-wise deviations toward a target (`lasso`/`fused`) or the pairwise
#' between-subject differences (`clustered`), and \eqn{\lambda} is the single
#' regularisation strength named by the spec. Stamps the `penaltySpec` on the
#' returned objfn so a composed objective
#' (`normL2(...) + constraintL1(...)`) self-describes for [EM] (see
#' `.nlmeReconstruct`).
#'
#' The Hessian contribution of the L1 term is zero (its second derivative
#' vanishes away from the kinks); the conditional-mode solve is handled by
#' [trustL1] inside [EM], not by feeding this term to [trust].
#'
#' @param penalty A `penaltySpec` (with subject expansion) from [penaltyL1],
#'   possibly combined with `+`.
#' @param attr.name Name of the numeric attribute holding the penalty value.
#'   Default `"prior_l1"` (distinct from constraintL2's `"prior"`, so the two
#'   accumulate independently under `+`).
#' @param condition Optional condition (default `NULL`, condition-unspecific).
#' @return An `objfn` carrying `attr(., "penaltySpec")`.
#' @seealso [penaltyL1], [EM], [sparsify], [constraintL2]
#' @export
constraintL1 <- function(penalty, attr.name = "prior_l1", condition = NULL) {

  if (!inherits(penalty, "penaltyspec"))
    stop("`penalty` must be a penaltySpec built by penaltyL1().")
  if (is.null(penalty$subjectEtas))
    stop("`penalty` must have subject expansion. Set subjects = ... .")

  spec         <- penalty
  lambda_names <- spec$lambdaName
  terms        <- .penaltyTerms(spec)
  parnames     <- c(as.vector(spec$subjectEtas), lambda_names)

  # Clustered blocks give one term per subject pair, so the term list is
  # O(n^2). Flatten it once and resolve names to positions per parameter
  # vector, not per term: the loop below keeps its arithmetic order exactly.
  t_i   <- vapply(terms, function(z) z$i, "")
  t_j   <- vapply(terms, function(z) z$j, NA_character_)
  t_tgt <- vapply(terms, function(z)
    if (length(z$target) == 1L) as.numeric(z$target) else NA_real_, 0)
  t_lam <- vapply(terms, function(z) z$lambda, "")
  .idx  <- new.env(parent = emptyenv())

  myfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE,
                   conditions = condition, env = NULL,
                   cores = getOption("dMod.cores", 1L)) {

    p    <- list(...)[[match.fnargs(list(...), "pars")]]
    allp <- c(p, fixed)
    np   <- names(p)

    if (!all(parnames %in% names(allp)))
      return(objlist(value = 0,
                     gradient = setNames(numeric(length(p)), np),
                     hessian  = matrix(0, length(p), length(p),
                                       dimnames = list(np, np))))

    key <- list(names(allp), np)
    if (!identical(.idx$key, key)) {
      .idx$key <- key
      .idx$ai <- match(t_i,   names(allp)); .idx$gi <- match(t_i,   np)
      .idx$aj <- match(t_j,   names(allp)); .idx$gj <- match(t_j,   np)
      .idx$al <- match(t_lam, names(allp)); .idx$gl <- match(t_lam, np)
    }
    ai <- .idx$ai; aj <- .idx$aj; al <- .idx$al
    gi <- .idx$gi; gj <- .idx$gj; gl <- .idx$gl
    av <- as.numeric(allp)

    gr  <- setNames(numeric(length(p)), np)
    val <- 0
    for (k in seq_along(t_i)) {
      lambda <- av[al[k]]                       # this term's strength
      a  <- av[ai[k]]
      b  <- if (is.na(aj[k])) t_tgt[k] else av[aj[k]]
      d  <- a - b
      ad <- abs(d)
      val <- val + lambda * ad
      s   <- if (d > 0) 1 else if (d < 0) -1 else 0
      if (!is.na(gi[k])) gr[gi[k]] <- gr[gi[k]] + lambda * s
      if (!is.na(gj[k])) gr[gj[k]] <- gr[gj[k]] - lambda * s
      if (!is.na(gl[k])) gr[gl[k]] <- gr[gl[k]] + ad
    }

    out <- objlist(value = val, gradient = gr,
                   hessian = matrix(0, length(p), length(p),
                                    dimnames = list(np, np)))
    attr(out, attr.name) <- val
    attr(out, "env") <- env
    out
  }

  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- condition
  attr(myfn, "parameters") <- parnames
  # NLME reconstruction handle: expose the penaltySpec so a composed objective
  # (normL2 + constraintL1) self-describes (see .normalReconstruct()).
  attr(myfn, "penaltySpec") <- spec
  myfn
}
