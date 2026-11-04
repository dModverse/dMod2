#' Random effects of a parameter transformation
#'
#' @description
#' Returns the random-effect structure that [insert()] with effect formulas
#' attached to a branched transformation: the random terms, their groupings,
#' distributions and correlation blocks, and the parameters the formulas
#' created. [em()] takes it as `random` and builds the prior of the random
#' effects from it. [P()] keeps the structure, so `randomEffects()` takes the
#' transformation or the parameter function built from it.
#'
#' An effect formula `cl ~ wt + (1 + dose | p | ID)` replaces the symbol `cl`
#' wherever it occurs in a branch by a linear predictor. The scale of the
#' symbol is the one the transformation gives it: after
#' `insert(trafo, "x <- exp(x)", x = "cl")` the random effects are normal on
#' the log scale. Its parts:
#' \describe{
#'   \item{left-hand side}{one or more symbols of the transformation,
#'     `ka + v ~ ...` applies the right-hand side to each.}
#'   \item{fixed effects}{the terms outside the bars, coded as by
#'     [stats::model.matrix()] (treatment contrasts) over the columns of the
#'     branch table. The intercept keeps the name of the symbol (`cl`), the
#'     other columns are `cl_<column>` (`cl_wt`, `cl_log_wt`).}
#'   \item{random terms}{`(1 | ID)`, `(1 + dose | ID)` with a random slope,
#'     `(1 | p | ID)` with an ID `p`: terms with the same ID and grouping are
#'     correlated across formulas, the columns of one term are correlated, other
#'     terms are independent. `gr(ID, dist = laplace())` instead of `ID` sets
#'     the distribution, see [normal()]. The random effect of a group is
#'     `eta_<symbol>_<group>` for the intercept and
#'     `eta_<symbol>_<column>_<group>` for a slope; a symbol with random terms
#'     in a second grouping gets its name in the stem, `eta_cl_centre_<group>`.}
#' }
#'
#' @param trafo A branched transformation from [insert()] with effect
#'   formulas, or a parameter function [P()] built from one.
#' @param x A `randomeffects` object.
#' @param ... Ignored.
#'
#' @return An object of class `randomeffects`, a list with `terms` (one row
#'   per random effect column: symbol, term of its formula, column, grouping,
#'   ID, distribution, its fixed hyperparameters and the name stem of
#'   its etas),
#'   `blocks` (the correlated groups of rows), `subjects` (the groups of each
#'   grouping and the conditions they own), `fixed` (the parameters the
#'   formulas created), `etas` (eta names, one row per condition, one column
#'   per term row) and `fixedDesign` (per symbol the model matrix of its fixed
#'   effects, one row per condition).
#' @seealso [insert()], [em()], [normal()]
#' @example inst/examples/randomEffects.R
#' @export
randomEffects <- function(trafo) {
  re <- attr(trafo, "randomEffects", exact = TRUE)
  if (is.null(re))
    stop("randomEffects: no random effects; add them with ",
         "insert(branch(trafo, table), par ~ ... + (1 | group)).", call. = FALSE)
  re
}


#' @export
#' @rdname randomEffects
print.randomeffects <- function(x, ...) {
  tm <- x$terms
  cat(sprintf("Random effects: %d %s over %s\n", nrow(tm), if (nrow(tm) == 1L) "term" else "terms",
              paste(unique(tm$group), collapse = ", ")))
  for (b in seq_along(x$blocks)) {
    r <- tm[x$blocks[[b]], , drop = FALSE]
    given <- function(v) if (is.na(v)) "estimated" else format(v)
    dist <- switch(r$dist[1], normal = "normal",
                   laplace = sprintf("laplace, lambda %s", given(r$lambda[1])),
                   spikeSlab = sprintf("spikeSlab, prob %s, tau %s", given(r$prob[1]), given(r$tau[1])),
                   dirichletProcess = sprintf("dirichletProcess, concentration %s, tau %s",
                                              given(r$concentration[1]), given(r$tau[1])))
    cat(sprintf("  block %d (%s, %s): %s\n", b, dist, r$group[1], paste(r$stem, collapse = ", ")))
  }
  if (length(x$fixed)) cat("  fixed effects:", paste(x$fixed, collapse = ", "), "\n")
  invisible(x)
}


#' Distributions of random effects
#'
#' @description The distribution of the random effects of a term, set in an
#' effect formula by `gr(group, dist = )`, e.g.
#' `cl ~ 1 + (1 | gr(cellline, dist = laplace()))`. Their parameters are
#' estimated by [em()] unless given.
#'
#' - `normal()`: multivariate normal with mean 0; its covariance is estimated
#'   per correlation block. The default.
#' - `laplace()`: Laplace with mean 0 and strength `lambda`, shared by all
#'   symbols of the term, so the deviation of a group from the fixed effects
#'   is exactly 0 unless the data need it.
#' - `spikeSlab()`: the deviation of a group is exactly 0 (spike) or drawn
#'   from \eqn{N(0, \tau^2)}{N(0, tau^2)} (slab); every group deviates with
#'   probability `prob`, independently. Up to 12 spike-and-slab random
#'   effects of a group are integrated jointly and exactly; with more, each is
#'   integrated alone and the coupling by the data is approximated as set by
#'   `control$laplace$coupling` of [em()].
#' - `dirichletProcess()`: per symbol, the deviations follow a Dirichlet
#'   process with base measure \eqn{N(0, \tau^2)}{N(0, tau^2)} and
#'   `concentration`, the fixed effect seated first at deviation 0: the
#'   groups and the fixed effect form clusters under the Ewens prior of the
#'   Chinese restaurant process, the cluster of the fixed effect at exactly
#'   0 and every other cluster one common deviation from the base measure.
#'   The prior depends on the sizes of the clusters only, not on the order of
#'   the groups. [em()] sums exactly over all partitions of at most 14
#'   groups and reports the most probable clusters; `method = "saem"` draws
#'   the clusters of any number of groups.
#'
#' [em()] estimates `prob`, `concentration` and `tau` unless given.
#'
#' @param lambda Strength of the Laplace distribution; `NULL` (default)
#'   estimates it, a positive number fixes it.
#' @param prob Probability that a group deviates, in (0, 1); `NULL`
#'   (default) estimates it.
#' @param concentration Concentration of the Ewens prior, positive; `NULL`
#'   (default) estimates it.
#' @param tau Standard deviation of the slab; `NULL` (default) estimates it.
#' @return An object of class `redist`.
#' @seealso [randomEffects()], [insert()]
#' @example inst/examples/normal.R
#' @export
normal <- function() .redist("normal")

#' @export
#' @rdname normal
laplace <- function(lambda = NULL) {
  if (!is.null(lambda) && !(is.numeric(lambda) && length(lambda) == 1L && lambda > 0))
    stop("laplace: `lambda` is NULL or one positive number.", call. = FALSE)
  .redist("laplace", lambda = lambda)
}

#' @export
#' @rdname normal
spikeSlab <- function(prob = NULL, tau = NULL) {
  if (!.positiveOrNull(prob) || (!is.null(prob) && prob >= 1))
    stop("spikeSlab: `prob` is NULL or one number in (0, 1).", call. = FALSE)
  if (!.positiveOrNull(tau)) stop("spikeSlab: `tau` is NULL or one positive number.", call. = FALSE)
  .redist("spikeSlab", prob = prob, tau = tau)
}

#' @export
#' @rdname normal
dirichletProcess <- function(concentration = NULL, tau = NULL) {
  if (!.positiveOrNull(concentration))
    stop("dirichletProcess: `concentration` is NULL or one positive number.", call. = FALSE)
  if (!.positiveOrNull(tau)) stop("dirichletProcess: `tau` is NULL or one positive number.", call. = FALSE)
  .redist("dirichletProcess", concentration = concentration, tau = tau)
}

.positiveOrNull <- function(x) is.null(x) || (is.numeric(x) && length(x) == 1L && x > 0)

.redist <- function(dist, lambda = NULL, prob = NULL, concentration = NULL, tau = NULL)
  structure(list(dist = dist, lambda = lambda %||% NA_real_, prob = prob %||% NA_real_,
                 concentration = concentration %||% NA_real_, tau = tau %||% NA_real_), class = "redist")

#' @export
print.redist <- function(x, ...) {
  arg <- function(v) if (is.na(v)) "NULL" else format(v)
  cat(switch(x$dist, normal = "normal()",
             laplace = sprintf("laplace(lambda = %s)", arg(x$lambda)),
             spikeSlab = sprintf("spikeSlab(prob = %s, tau = %s)", arg(x$prob), arg(x$tau)),
             dirichletProcess = sprintf("dirichletProcess(concentration = %s, tau = %s)", arg(x$concentration),
                                        arg(x$tau))), "\n")
  invisible(x)
}


.hasBar <- function(e) {
  if (!is.call(e)) return(FALSE)
  if (identical(e[[1]], as.name("|"))) return(TRUE)
  any(vapply(as.list(e)[-1], .hasBar, TRUE))
}

## Top-level summands of a formula side with their sign.
.effectSummands <- function(e, sign = 1) {
  if (is.call(e) && length(e) == 3L && (identical(e[[1]], as.name("+")) || identical(e[[1]], as.name("-"))))
    return(c(.effectSummands(e[[2]], sign),
             .effectSummands(e[[3]], if (identical(e[[1]], as.name("-"))) -sign else sign)))
  list(list(e = e, sign = sign))
}

## Replace the names in `e` by the values of `map`, a named list of names.
.substituteNames <- function(e, map) {
  if (is.name(e)) {
    v <- map[[as.character(e)]]
    return(if (is.null(v)) e else as.name(v))
  }
  if (!is.call(e)) return(e)
  as.call(lapply(as.list(e), .substituteNames, map = map))
}

## Effect formulas with their placeholders expanded: element i of every named
## vector in `placeholders` replaces the name it is given under.
.expandEffectFormulas <- function(formulas, placeholders) {
  if (!length(placeholders))
    return(Map(function(f, k) structure(f, origin = k), formulas, seq_along(formulas)))
  n <- unique(lengths(placeholders))
  if (length(n) != 1L)
    stop("insert: the placeholders of effect formulas have different lengths.", call. = FALSE)
  unlist(lapply(seq_along(formulas), function(k) lapply(seq_len(n), function(i) {
    f <- formulas[[k]]
    map <- lapply(placeholders, function(v) as.character(v[[i]]))
    g <- .substituteNames(f, map)
    environment(g) <- environment(f)
    structure(g, class = "formula", origin = attr(f, "origin") %||% k)
  })), recursive = FALSE)
}

## The symbols of the left-hand side `a + b`.
.effectSymbols <- function(lhs) {
  if (is.name(lhs)) return(as.character(lhs))
  if (is.call(lhs) && identical(lhs[[1]], as.name("+")) && length(lhs) == 3L)
    return(c(.effectSymbols(lhs[[2]]), .effectSymbols(lhs[[3]])))
  stop("insert: the left-hand side of an effect formula is one or more symbols of the trafo, ",
       "as in `ka + v ~ ...`, not ", deparse1(lhs), ". Set the scale in the trafo first, e.g. ",
       "insert(trafo, \"x <- exp(x)\", x = \"cl\").", call. = FALSE)
}

## One effect formula: its symbols, the fixed part as a one-sided formula and the random terms.
.parseEffectFormula <- function(f) {
  if (length(f) != 3L)
    stop("insert: an effect formula has a left-hand side, as in `cl ~ wt + (1 | ID)`.", call. = FALSE)
  symbols <- .effectSymbols(f[[2]])
  sm <- .effectSummands(f[[3]])
  isBar <- vapply(sm, function(t) {
    e <- t$e
    while (is.call(e) && identical(e[[1]], as.name("("))) e <- e[[2]]
    is.call(e) && identical(e[[1]], as.name("|"))
  }, TRUE)
  fixedTerms <- vapply(sm[!isBar], function(t) paste(if (t$sign < 0) "-" else "+", deparse1(t$e)), "")
  fixed <- stats::as.formula(paste("~", if (length(fixedTerms)) paste(fixedTerms, collapse = " ") else "1"))
  random <- lapply(sm[isBar], function(t) {
    if (t$sign < 0) stop("insert: a random term cannot be subtracted.", call. = FALSE)
    e <- t$e
    while (identical(e[[1]], as.name("("))) e <- e[[2]]
    .parseBar(e, environment(f))
  })
  lapply(symbols, function(s) list(symbol = s, fixed = fixed, random = random,
                                   origin = attr(f, "origin") %||% 1L))
}

## `expr | group`, `expr | id | group`; group a column, an interaction a:b or gr(group, dist = ).
.parseBar <- function(e, env) {
  id <- NA_character_
  left <- e[[2]]; grp <- e[[3]]
  if (is.call(left) && identical(left[[1]], as.name("|"))) {
    if (!is.name(left[[3]]))
      stop("insert: the ID of a random term is a name, as in (1 | p | ID).", call. = FALSE)
    id <- as.character(left[[3]]); left <- left[[2]]
  }
  dist <- normal()
  if (is.call(grp) && identical(grp[[1]], as.name("gr"))) {
    a <- as.list(grp)[-1]
    nm <- names(a) %||% rep("", length(a))
    if (any(!nm %in% c("", "dist")))
      stop("insert: gr() takes the grouping and `dist`, as in gr(ID, dist = laplace()).", call. = FALSE)
    if ("dist" %in% nm) dist <- .evalDistribution(a[["dist"]], env)
    grp <- a[[which(nm == "")[1]]]
  }
  vars <- all.vars(grp)
  if (!length(vars)) stop("insert: a random term needs a grouping column.", call. = FALSE)
  list(expr = stats::as.formula(call("~", left)), group = paste(vars, collapse = ":"), vars = vars,
       id = id, dist = dist)
}

## `dist =` of gr(): normal(), laplace(), spikeSlab() or dirichletProcess() of this package,
## whatever else is attached.
.evalDistribution <- function(e, env) {
  scope <- new.env(parent = env %||% globalenv())
  scope$normal <- normal; scope$laplace <- laplace; scope$spikeSlab <- spikeSlab
  scope$dirichletProcess <- dirichletProcess
  d <- eval(e, scope)
  if (!inherits(d, "redist"))
    stop("insert: gr(dist = ) is normal(), laplace(), spikeSlab() or dirichletProcess(), not ", deparse1(e), ".",
         call. = FALSE)
  d
}

.effectName <- function(x) gsub("[^[:alnum:]_]+", "_", gsub("^\\(|\\)$", "", x))

## Stops where two distinct values of `x` get the same parameter name.
.uniqueNames <- function(x, ...) {
  u <- unique(x)
  dup <- duplicated(.effectName(u))
  if (any(dup))
    stop("insert: the ", ..., " give the same name to ", paste(u[.effectName(u) %in% .effectName(u)[dup]],
         collapse = " and "), "; rename them.", call. = FALSE)
}

.formatNumber <- function(x) sprintf("%.17g", x)

## Columns of the branch table a formula reads; a missing one may be a substitution in disguise.
.checkEffectColumns <- function(pf, tree) {
  vars <- unique(c(all.vars(pf$fixed), unlist(lapply(pf$random, function(r) c(all.vars(r$expr), r$vars)))))
  miss <- setdiff(vars, names(tree))
  if (!length(miss)) return(invisible())
  stop("insert: ", paste(miss, collapse = ", "), " in the effect formula of ", pf$symbol,
       " is not a column of the branch table. To substitute a symbol write \"", pf$symbol, " <- ",
       sub("^\\+\\s*", "", deparse1(pf$fixed[[2]])), "\".", call. = FALSE)
}

## The design of one formula over the rows of `tree`: the fixed-effect matrix with the names of
## its parameters, and per random-effect column its stem, the eta and coefficient of every row.
## A stem already taken by a term of another grouping gets the grouping appended. With `sparse`
## the fixed-effect matrix is a sparse Matrix.
.effectMatrices <- function(pf, tree, taken = NULL, sparse = FALSE) {
  .checkEffectColumns(pf, tree)
  sym <- pf$symbol
  frame <- stats::model.frame(pf$fixed, tree, na.action = stats::na.fail)
  X <- if (sparse) Matrix::sparse.model.matrix(pf$fixed, frame) else stats::model.matrix(pf$fixed, frame)
  cn <- colnames(X)
  .uniqueNames(cn, "fixed effects of ", sym)
  colnames(X) <- ifelse(cn == "(Intercept)", sym, paste0(sym, "_", .effectName(cn)))
  random <- list()
  for (ti in seq_along(pf$random)) {
    rt <- pf$random[[ti]]
    g <- do.call(paste, c(lapply(tree[rt$vars], as.character), sep = "_"))
    # a condition with a missing grouping belongs to no group and keeps the fixed effects
    none <- Reduce(`|`, lapply(tree[rt$vars], is.na))
    g[none] <- NA_character_
    .uniqueNames(g[!none], "groups of ", rt$group)
    Z <- stats::model.matrix(rt$expr, stats::model.frame(rt$expr, tree, na.action = stats::na.fail))
    Z[none, ] <- 0
    zn <- colnames(Z)
    stems <- ifelse(zn == "(Intercept)", paste0("eta_", sym), paste0("eta_", sym, "_", .effectName(zn)))
    other <- Filter(function(r) r$group != rt$group, random)
    own <- c(taken$stem[taken$group != rt$group], vapply(other, `[[`, "", "stem"))
    clash <- stems %in% own
    stems[clash] <- paste0(stems[clash], "_", .effectName(rt$group))
    for (k in seq_along(zn))
      random[[length(random) + 1L]] <- list(symbol = sym, term = ti, column = zn[k], group = rt$group,
                                            id = rt$id, dist = rt$dist$dist,
                                            lambda = rt$dist$lambda, prob = rt$dist$prob,
                                            concentration = rt$dist$concentration, tau = rt$dist$tau,
                                            stem = stems[k],
                                            etas = ifelse(none, NA_character_, paste0(stems[k], "_", .effectName(g))),
                                            z = Z[, k],
                                            groups = g)
  }
  list(X = X, random = random, columns = intersect(all.vars(pf$fixed), names(tree)))
}

## Linear predictor of one formula for every row of `tree`, the fixed parameters it creates and
## one row per random-effect column.
.effectPredictor <- function(pf, tree, taken = NULL) {
  em <- .effectMatrices(pf, tree, taken)
  term <- function(coef, name) ifelse(coef == 1, name, paste0(name, " * ", .formatNumber(coef)))
  pieces <- lapply(seq_len(nrow(tree)), function(i) {
    v <- em$X[i, ]
    out <- term(v[v != 0], colnames(em$X)[v != 0])
    for (r in em$random) if (r$z[i] != 0) out <- c(out, term(r$z[i], r$etas[i]))
    out
  })
  lin <- vapply(pieces, function(p) if (length(p)) paste(p, collapse = " + ") else "0", "")
  list(expr = paste0("(", lin, ")"), fixed = colnames(em$X), rows = em$random,
       X = `rownames<-`(em$X, rownames(tree)))
}

## Apply effect formulas to the branches of `trafo` and attach the random-effect structure.
.insertEffects <- function(trafo, formulas, conditionMatch) {
  tree <- attr(trafo, "tree")
  if (is.null(tree) || !is.list(trafo))
    stop("insert: effect formulas need a branched trafo, branch(trafo, table = ), whose table ",
         "holds the covariates and groupings; a substitution is a string, \"x <- exp(x)\".", call. = FALSE)
  conds <- names(trafo)
  use <- if (is.null(conditionMatch)) rep(TRUE, length(conds)) else grepl(conditionMatch, conds)
  tr <- tree[conds[use], , drop = FALSE]
  prev <- attr(trafo, "randomEffects", exact = TRUE)
  terms <- if (!is.null(prev)) prev$terms else NULL
  etaCols <- if (!is.null(prev)) as.list(prev$etas) else list()
  groups <- if (!is.null(prev)) prev$groups else list()
  fixedPars <- if (!is.null(prev)) prev$fixed else character(0)
  fixedDesign <- if (!is.null(prev)) prev$fixedDesign else list()
  # formulas of one call share an origin with their placeholder and `+` expansions
  base <- if (!is.null(terms) && nrow(terms)) max(terms$origin) else 0L
  parsed <- unlist(lapply(formulas, .parseEffectFormula), recursive = FALSE)
  out <- trafo
  for (pf in parsed) {
    pf$origin <- base + pf$origin
    if (!any(vapply(trafo[use], function(t) pf$symbol %in% getSymbols(t), TRUE)))
      stop("insert: the symbol ", pf$symbol, " of an effect formula is not in the trafo.", call. = FALSE)
    if (pf$symbol %in% .effectSymbolsOf(prev, terms, fixedPars))
      stop("insert: ", pf$symbol, " already has an effect formula.", call. = FALSE)
    ep <- .effectPredictor(pf, tr, terms)
    for (i in seq_along(ep$expr)) {
      cn <- rownames(tr)[i]
      out[[cn]] <- replaceSymbols(pf$symbol, ep$expr[i], out[[cn]])
    }
    fixedPars <- union(fixedPars, ep$fixed)
    fixedDesign[[pf$symbol]] <- ep$X
    for (r in ep$rows) {
      if (!is.null(terms) && r$stem %in% terms$stem)
        stop("insert: the random effect ", r$stem, " is defined twice.", call. = FALSE)
      row <- data.frame(symbol = r$symbol, origin = pf$origin, term = r$term, column = r$column,
                        group = r$group, id = r$id, dist = r$dist,
                        lambda = r$lambda, prob = r$prob, concentration = r$concentration, tau = r$tau,
                        stem = r$stem, stringsAsFactors = FALSE)
      terms <- rbind(terms, row)
      col <- setNames(rep(NA_character_, length(conds)), conds)
      col[rownames(tr)] <- r$etas
      etaCols[[r$stem]] <- unname(col)
      groups[[r$stem]] <- setNames(r$groups, rownames(tr))
    }
  }
  attr(out, "tree") <- tree
  etas <- as.data.frame(etaCols, row.names = conds, stringsAsFactors = FALSE, optional = TRUE)
  attr(out, "randomEffects") <- .newRandomEffects(terms, etas, groups, fixedPars, fixedDesign)
  out
}

## Symbols that have an effect formula: the fixed parameters the formulas created.
.effectSymbolsOf <- function(re, terms = re$terms, fixed = re$fixed)
  unique(c(terms$symbol, fixed))

## Substituting a symbol that has an effect formula would rewrite its linear predictor.
.guardEffectSymbols <- function(lhs, re) {
  hit <- intersect(lhs, .effectSymbolsOf(re))
  if (length(hit))
    stop("insert/define: ", paste(hit, collapse = ", "), " has an effect formula; set the ",
         "scale of a parameter before adding effects to it.", call. = FALSE)
}

## Blocks: normal rows of one grouping with one ID, or of one term without an ID, are correlated;
## the other rows of one term of one formula share their hyperparameters. `fixedDesign` holds per
## symbol the model matrix of its fixed effects, rows the conditions.
.newRandomEffects <- function(terms, etas, groups, fixedPars, fixedDesign = list(), subjects = NULL,
                              etaNames = NULL) {
  if (is.null(terms))
    terms <- data.frame(symbol = character(0), origin = integer(0), term = integer(0),
                        column = character(0),
                        group = character(0), id = character(0), dist = character(0),
                        lambda = numeric(0), prob = numeric(0),
                        concentration = numeric(0), tau = numeric(0), stem = character(0))
  key <- ifelse(terms$dist != "normal", paste(terms$dist, terms$origin, terms$term, sep = "\r"),
         ifelse(is.na(terms$id), paste("term", terms$symbol, terms$term, sep = "\r"),
                paste("id", terms$id, terms$group, terms$dist, sep = "\r")))
  blocks <- unname(split(seq_len(nrow(terms)), factor(key, unique(key))))
  if (any(terms$dist != "normal" & !is.na(terms$id)))
    stop("insert: Laplace, spike-and-slab and Dirichlet-process random effects are not correlated; drop the ID of their ",
         "terms.", call. = FALSE)
  subjects <- subjects %||% lapply(split(as.character(names(groups)), terms$group), function(st) {
    g <- groups[[st[1]]]
    split(names(g), factor(g, unique(g)))
  })
  structure(list(terms = terms, blocks = blocks, subjects = subjects, fixed = fixedPars,
                 etas = etas, groups = groups, fixedDesign = fixedDesign, etaNames = etaNames),
            class = "randomeffects")
}
