# The estimation design of a mixed-effects fit: the model pieces of the objective, the groups of
# the random effects with their conditions and etas, the distribution blocks and the names of
# their hyperparameters. Everything downstream of em() reads the design, never the spec.


## The plans of the native chain of the data term `spec` of `obj`: those of its normL2() term, so
## that they last over the fits of one objective; a context of their own where it has none.
.objectiveChain <- function(obj, spec) {
  if (is.null(spec)) return(NULL)
  term <- Find(function(t) identical(attr(t, "l2spec", exact = TRUE)[[1]], spec), .objTerms(obj))
  ctx <- if (is.function(term)) get0(".chain", envir = environment(term), inherits = FALSE)
  if (is.environment(ctx)) ctx else .chainContext(spec)
}

## The design of `obj` and the random effects `random`: with one grouping its groups, etas and
## blocks are fields of the design, with several `groupings` holds one design per grouping and
## `components` the sets of conditions their groups connect.
.emDesign <- function(obj, random) {
  if (!inherits(obj, "objfn"))
    stop("em: `obj` is the objective of the data, normL2() or normT().", call. = FALSE)
  if (!inherits(random, "randomeffects"))
    stop("em: `random` is randomEffects(trafo) or randomEffects(p) of the parameter ",
         "transformation with effect formulas.", call. = FALSE)
  prdfn <- attr(obj, "prdfn", exact = TRUE)
  data  <- attr(obj, "data", exact = TRUE)
  if (is.null(prdfn) || is.null(data))
    stop("em: `obj` holds no prediction and data; build it as normL2(data, g * x * p).",
         call. = FALSE)
  spec <- attr(obj, "l2spec", exact = TRUE)[[1]]
  if (!nrow(random$terms))
    stop("em: the effect formulas have no random term, as in (1 | ID).", call. = FALSE)
  subs <- lapply(setNames(nm = unique(random$terms$group)), function(gr) .groupingDesign(random, data, gr))
  base <- list(prdfn = prdfn, errfn = attr(obj, "errfn", exact = TRUE), data = data, objfn = obj,
               df = spec$df %||% Inf, bloq = spec$optBLOQ %||% "M3", fixedEffects = random$fixed,
               randomEffects = random, chain = .objectiveChain(obj, spec))
  inGroup <- unique(unlist(lapply(subs, function(sd) unlist(sd$conditions))))
  free <- setdiff(names(data), inGroup)
  if (length(subs) == 1L)
    return(structure(c(base, subs[[1]], list(freeConditions = free)), class = "emdesign"))
  structure(c(base, list(groupings = subs, components = .components(subs), freeConditions = free,
                         blocks = unlist(lapply(subs, `[[`, "blocks"), recursive = FALSE))),
            class = "emdesign")
}

## The groups of one grouping: their conditions with data, etas, distribution blocks and the
## Cholesky layout of the normal ones.
.groupingDesign <- function(random, data, grouping) {
  terms <- random$terms
  conditions <- lapply(random$subjects[[grouping]], intersect, names(data))
  conditions <- conditions[lengths(conditions) > 0L]
  groups <- names(conditions)
  own <- Filter(function(b) terms$group[b[1]] == grouping, random$blocks)
  stems <- terms$stem[unlist(own)]
  # the etas of a group, given per grouping where a condition holds several groups of it
  etaNames <- if (!is.null(random$etaNames[[grouping]])) random$etaNames[[grouping]][groups, stems, drop = FALSE]
              else do.call(rbind, lapply(groups, function(g)
                unlist(random$etas[conditions[[g]][1], stems], use.names = FALSE)))
  dimnames(etaNames) <- list(groups, stems)
  blocks <- lapply(own, function(b) .designBlock(terms[b, , drop = FALSE]))
  dist <- vapply(blocks, `[[`, "", "dist")
  blocks[dist == "laplace"] <- .nameHypers(blocks[dist == "laplace"], grouping, "lambda")
  slab <- dist %in% c("spikeSlab", "dirichletProcess")
  blocks[slab] <- .nameHypers(blocks[slab], grouping, function(b)
    c(if (b$dist == "spikeSlab") "prob" else "concentration", "tau"))
  list(grouping = grouping, groups = groups, conditions = conditions, etaNames = etaNames,
       blocks = blocks, normal = .normalLayout(blocks[dist == "normal"]))
}

## The connected components of the conditions: two conditions are connected when they share a
## group of some grouping.
.components <- function(subs) {
  conds <- unique(unlist(lapply(subs, function(sd) unlist(sd$conditions))))
  root <- setNames(conds, conds)
  find <- function(x) { while (root[[x]] != x) x <- root[[x]]; x }
  for (sd in subs) for (cs in sd$conditions) {
    r <- find(cs[1])
    for (c2 in cs[-1]) { r2 <- find(c2); if (r2 != r) root[r2] <- r }
  }
  top <- vapply(conds, find, "")
  lapply(unname(split(conds, factor(top, unique(top)))), function(cs) list(conditions = cs))
}

## The random effects `re` without the term of `stem`, whose etas the caller fixes at 0.
.dropRandomTerm <- function(re, stem) {
  keep <- re$terms$stem != stem
  etaNames <- lapply(re$etaNames, function(E) E[, setdiff(colnames(E), stem), drop = FALSE])
  .newRandomEffects(re$terms[keep, , drop = FALSE], re$etas[, keep, drop = FALSE], re$groups[keep], re$fixed,
                    re$fixedDesign, subjects = re$subjects[unique(re$terms$group[keep])],
                    etaNames = if (length(re$etaNames)) etaNames)
}

## The grouping designs of a design, one with a single grouping.
.subDesigns <- function(design) unname(design$groupings) %||% list(design)

## Every eta of the design.
.allEtas <- function(design) unlist(lapply(.subDesigns(design), function(sd) as.vector(sd$etaNames)))


## One distribution block: its stems (random-effect columns) and short names, the distribution
## and its fixed hyperparameters, NA when estimated.
.designBlock <- function(rows) {
  list(dist = rows$dist[1], lambda = rows$lambda[1], prob = rows$prob[1],
       concentration = rows$concentration[1], tau = rows$tau[1], stems = rows$stem,
       short = sub("^eta_", "", rows$stem))
}

## Names of the hyperparameters of the blocks of one distribution, `hyper` named by `kinds`
## (a character vector or a function of the block): <kind>_<grouping>, numbered when the grouping
## has several such blocks.
.nameHypers <- function(blocks, grouping, kinds) {
  suffix <- if (length(blocks) == 1L) "" else paste0("_", seq_along(blocks))
  Map(function(b, s) {
    k <- if (is.function(kinds)) kinds(b) else kinds
    b$hyper <- setNames(paste0(k, "_", .effectName(grouping), s), k)
    b
  }, blocks, suffix)
}

## The estimated hyperparameters of a block, named by kind.
.estimatedHypers <- function(b) b$hyper[is.na(unlist(b[names(b$hyper)]))]


## Cholesky layout of the normal random effects: one lower-triangular factor over all normal
## stems, block diagonal by the blocks, log on the diagonal. Names omega_<a>_<a> on the
## diagonal and omega_<a>_<b> below it, a the later stem.
.normalLayout <- function(blocks) {
  stems <- unlist(lapply(blocks, `[[`, "stems"))
  short <- unlist(lapply(blocks, `[[`, "short"))
  K <- length(stems)
  loc <- matrix(integer(0), 0, 2)
  names <- character(0)
  start <- 0L
  for (b in blocks) {
    k <- start + seq_along(b$stems)
    for (l in k) for (r in k[k >= l]) {
      loc <- rbind(loc, c(r, l))
      names <- c(names, paste("omega", short[r], short[l], sep = "_"))
    }
    start <- start + length(b$stems)
  }
  dimnames(loc) <- list(names, c("row", "col"))
  list(K = K, stems = stems, short = short, cholPars = names, loc = loc,
       isDiag = setNames(loc[, 1] == loc[, 2], names))
}

## The Cholesky factor of Omega from the Cholesky parameters in `pars`.
.choleskyFactor <- function(layout, pars) {
  L <- matrix(0, layout$K, layout$K, dimnames = list(layout$stems, layout$stems))
  v <- pars[layout$cholPars]
  L[layout$loc] <- ifelse(layout$isDiag, exp(v), v)
  L
}


## The parameters of the fit: the outer parameters of the prediction without the etas and the
## hyperparameters; the error model sees the inner ones.
.structuralNames <- function(design) {
  all <- attr(design$prdfn, "parameters")
  setdiff(all, c(.allEtas(design), .hyperNames(design)))
}

.hyperNames <- function(design)
  c(unlist(lapply(.subDesigns(design), function(sd) sd$normal$cholPars)),
    unname(unlist(lapply(design$blocks, `[[`, "hyper"))))


## A complete start: the structural parameters from `init`, the fixed effects the formulas
## created at 0, the random-effect and slab standard deviations at `sd`, Laplace strengths at 1,
## split probabilities at 1/2.
.emStart <- function(design, init = NULL, fixed = NULL, sd = 0.3) {
  structural <- setdiff(.structuralNames(design), names(fixed))
  miss <- setdiff(structural, c(design$fixedEffects, names(init)))
  if (length(miss))
    stop("em: `init` needs values for ", paste(miss, collapse = ", "),
         "; the parameters the effect formulas created start at 0.", call. = FALSE)
  start <- setNames(rep(0, length(structural)), structural)
  given <- intersect(names(init), structural)
  start[given] <- init[given]
  chol <- unlist(unname(lapply(.subDesigns(design), function(sub)
    setNames(ifelse(sub$normal$isDiag, log(sd), 0), sub$normal$cholPars))))
  initial <- c(lambda = 1, prob = 0.5, concentration = 1, tau = sd)
  sparse <- unlist(unname(lapply(design$blocks, function(b) {
    h <- .estimatedHypers(b)
    setNames(initial[names(h)], h)
  })))
  hyper <- c(chol, sparse)
  given <- intersect(names(init), names(hyper))
  hyper[given] <- init[given]
  hyper <- hyper[setdiff(names(hyper), names(fixed))]
  unknown <- setdiff(names(init), c(structural, names(hyper), names(fixed)))
  if (length(unknown))
    stop("em: `init` names parameters the fit does not have: ", paste(unknown, collapse = ", "),
         ".", call. = FALSE)
  c(start, hyper)
}


## The mu-referenced fixed effects of grouping design `sub`, phi = X beta + eta: per normal random
## intercept, the first of its symbol, whose formula covers every condition with data and whose
## fixed-effect rows are constant within each group, the fixed parameters and their rows per group.
.muReference <- function(design, sub = design) {
  random <- design$randomEffects
  terms <- random$terms
  dataConditions <- names(design$data)
  out <- list()
  for (stem in sub$normal$stems) {
    r <- match(stem, terms$stem)
    sym <- terms$symbol[r]
    intercepts <- terms$stem[terms$symbol == sym & terms$column == "(Intercept)"]
    if (terms$column[r] != "(Intercept)" || stem != intercepts[1]) next
    if (anyNA(random$etas[dataConditions, stem])) next
    X <- random$fixedDesign[[sym]]
    if (is.null(X)) next
    rows <- lapply(sub$conditions, function(cs) unique(X[cs, , drop = FALSE]))
    if (any(vapply(rows, nrow, 1L) != 1L)) next
    out[[stem]] <- list(stem = stem, symbol = sym, pars = colnames(X),
                        X = matrix(unlist(rows), length(rows), ncol(X), byrow = TRUE,
                                   dimnames = list(names(rows), colnames(X))))
  }
  out
}

## Generalised least squares for the referenced fixed effects: with `phi` the posterior means per
## group (columns the normal stems, phi in the referenced ones), the free entries of `beta` that
## minimise sum_g (phi_g - D_g beta)' Omega^-1 (phi_g - D_g beta).
.glsFixedEffects <- function(ref, phi, Omega, beta, free = names(beta)) {
  stems <- colnames(phi)
  pars <- unlist(lapply(unname(ref), `[[`, "pars"))
  W <- solve(Omega)
  free <- intersect(free, pars)
  fixed <- setdiff(pars, free)
  A <- matrix(0, length(free), length(free)); b <- numeric(length(free))
  for (g in rownames(phi)) {
    D <- matrix(0, length(stems), length(pars), dimnames = list(stems, pars))
    for (e in ref) D[e$stem, e$pars] <- e$X[g, ]
    rest <- phi[g, ] - D[, fixed, drop = FALSE] %*% beta[fixed]
    Df <- D[, free, drop = FALSE]
    A <- A + crossprod(Df, W %*% Df)
    b <- b + crossprod(Df, W %*% rest)
  }
  setNames(drop(solve(A, b)), free)
}

## A start perturbed by normal noise of standard deviation `sd` on the scale where each entry is
## unbounded: log for Laplace strengths, concentrations and slab standard deviations, logit for
## slab probabilities, the identity for the structural and Cholesky parameters.
.perturbStart <- function(design, center, sd) {
  kinds <- unlist(lapply(design$blocks, function(b) setNames(names(b$hyper), unname(b$hyper))))
  kind <- kinds[names(center)]
  noise <- stats::rnorm(length(center), 0, sd)
  out <- center + noise
  pos <- !is.na(kind) & kind %in% c("lambda", "concentration", "tau")
  out[pos] <- center[pos] * exp(noise[pos])
  prob <- !is.na(kind) & kind == "prob"
  out[prob] <- stats::plogis(stats::qlogis(center[prob]) + noise[prob])
  out
}
