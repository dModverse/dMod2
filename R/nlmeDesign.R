# The estimation design of a mixed-effects fit: the model pieces of the objective, the groups of
# the random effects with their conditions and etas, the distribution blocks and the names of
# their hyperparameters. Everything downstream of em() reads the design, never the spec.


## The plans of the native chain of the data term `spec` of `obj`: those of the normL2() term of `obj`, so
## that they last over the fits of one objective; a context of their own where it has none.
.objectiveChain <- function(obj, spec) {
  if (is.null(spec)) return(NULL)
  term <- Find(function(t) identical(attr(t, "l2spec", exact = TRUE)[[1]], spec), .objTerms(obj))
  ctx <- if (is.function(term)) get0(".chain", envir = environment(term), inherits = FALSE)
  if (is.environment(ctx)) ctx else .chainContext(spec)
}

## The design of `obj` and the random effects `random`: with one grouping the groups, etas and
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
  if (.refines(prdfn))
    stop("em: the prediction has optionsReverse = list(refine = TRUE), which applies to the ",
         "reverse gradient of normL2(). The sweeps of em() run forward over reverse without ",
         "refine; build the prediction for em() without it.", call. = FALSE)
  if (!nrow(random$terms))
    stop("em: the effect formulas have no random term, as in (1 | ID).", call. = FALSE)
  subs <- lapply(setNames(nm = unique(random$terms$group)), function(gr) .groupingDesign(random, data, gr))
  base <- list(prdfn = prdfn, errfn = attr(obj, "errfn", exact = TRUE), data = data, objfn = obj,
               df = spec$df %||% Inf, bloq = spec$optBLOQ %||% "M3", fixedEffects = random$fixed,
               randomEffects = random, chain = .objectiveChain(obj, spec))
  inGroup <- unique(unlist(lapply(subs, function(sd) unlist(sd$conditions))))
  free <- setdiff(names(data), inGroup)
  if (length(subs) == 1L) {
    design <- structure(c(base, subs[[1]], list(freeConditions = free)), class = "emdesign")
    design$reach <- .etaReach(design)
    return(design)
  }
  structure(c(base, list(groupings = subs, components = .components(subs), freeConditions = free,
                         blocks = unlist(lapply(subs, `[[`, "blocks"), recursive = FALSE))),
            class = "emdesign")
}

## Per group, whether each sparse eta reaches the data of the group conditions, along the dependencies
## of the parameter transformations, steady states, the ODE with events and the observables or
## the error model to a measured name; the normal etas count as reaching. A piece without
## equations counts as passing everything on, so that an eta is marked only where no path exists.
.etaReach <- function(design) {
  E <- design$etaNames
  sparse <- vapply(colnames(E), function(st) any(vapply(design$blocks, function(b)
    b$dist != "normal" && st %in% b$stems, TRUE)), TRUE)
  reach <- matrix(TRUE, nrow(E), ncol(E), dimnames = dimnames(E))
  if (!any(sparse)) return(reach)
  for (i in seq_len(nrow(E))) {
    etas <- E[i, sparse]
    hit <- setNames(logical(length(etas)), etas)
    for (cn in design$conditions[[design$groups[i]]]) hit <- hit | .reachData(design, cn, etas)
    reach[i, sparse] <- hit
  }
  reach
}

## Which of `etas` reach the data of condition `cn`: a graph of named nodes in layers, the outer
## parameters, the outputs of every transformation, the states, the observables and the error
## model, searched in C++ (graph_reach) from every eta to the measured names. Every expression is
## read with the names that are 0 in this condition set to 0 (.zeroOut()): a state that starts at 0,
## whose right-hand side vanishes with the other zeros and that no event moves stays 0 and passes
## nothing on. A piece without equations, such as a function written by hand, reaches all.
.reachData <- function(design, cn, etas) {
  every <- setNames(rep(TRUE, length(etas)), etas)
  lv <- .fnLeaves(design$prdfn, cn)
  kinds <- vapply(lv, `[[`, "", "kind")
  at <- which(kinds == "prdfn")
  if (length(at) != 1L || any(kinds[seq_len(at - 1L)] != "obsfn") || any(kinds[-seq_len(at)] != "parfn"))
    return(every)
  P2X <- lv[[at]]$kernel
  ode <- attr(P2X, "equations", exact = TRUE)
  xd <- attr(P2X, "depends", exact = TRUE)
  if (isTRUE(attr(P2X, "timeOnly", exact = TRUE)) || !is.null(xd)) ode <- list()
  if (is.null(ode)) return(every)
  # nodes are keys "layer\rname", numbered in C++; per layer the names it holds and those at 0
  node <- function(layer, name) paste(layer, name, sep = "\r")
  held <- list()
  zeroAt <- list()
  hold <- function(layer, name) {
    L <- as.character(layer)
    held[[L]][[length(held[[L]]) + 1L]] <<- name
  }
  from <- list(); to <- list()
  edge <- function(fl, fn, tl, tn) {
    if (!length(fn) || !length(tn)) return(invisible())
    from[[length(from) + 1L]] <<- node(fl, fn); to[[length(to) + 1L]] <<- node(tl, tn)
    hold(fl, fn); hold(tl, tn)
  }
  namesAt <- function(layer) unique(unlist(held[[as.character(layer)]], use.names = FALSE))
  zerosAt <- function(layer) zeroAt[[as.character(layer)]] %||% character(0)
  # the transformations from the outer parameters inwards, layer by layer; the names an input
  # passes on unchanged keep their node in the next layer
  hold(0L, etas)
  layer <- 0L
  for (leaf in rev(lv[-seq_len(at)])) {
    k <- leaf$kernel
    env <- environment(k)
    eq <- attr(k, "equations", exact = TRUE)
    zin <- zerosAt(layer)
    zout <- character(0)
    if (!is.null(attr(k, "chainImpl", exact = TRUE))) {
      # a steady state: every root from the inputs it depends on, the zero states at 0, the other
      # inputs passed on
      dep <- .pimplDepends(env, zin)
      for (pin in names(dep)) edge(layer, rep(pin, length(dep[[pin]])), layer + 1L, dep[[pin]])
      outs <- c(env$dependent, env$zero_states)
      hold(layer + 1L, outs)
      through <- setdiff(namesAt(layer), outs)
      edge(layer, through, layer + 1L, through)
      zout <- c(env$zero_states, intersect(through, zin))
    } else if (!is.null(eq)) {
      for (o in names(eq)) {
        d <- .zeroDepends(eq[[o]], zin)
        if (d$zero) zout <- c(zout, o)
        edge(layer, d$symbols, layer + 1L, rep(o, length(d$symbols)))
        hold(layer + 1L, o)
      }
      if (isTRUE(env$controls$attachInput)) {
        through <- setdiff(namesAt(layer), names(eq))
        edge(layer, through, layer + 1L, through)
        zout <- c(zout, intersect(through, zin))
      }
    } else return(every)
    zeroAt[[as.character(layer + 1L)]] <- zout
    layer <- layer + 1L
  }
  inner <- layer
  zeroPar <- zerosAt(inner)
  # the states that stay 0: they start at 0, their right-hand side vanishes with the other zeros,
  # and no event sets them to anything else
  events <- if (length(ode)) attr(environment(P2X)$func, "events", exact = TRUE)
  if (!is.data.frame(events)) events <- NULL
  moves <- function(zs) {
    if (is.null(events)) return(character(0))
    v <- as.character(events$var)
    keep <- vapply(seq_len(nrow(events)), function(r) {
      method <- as.character(events$method[r] %||% "replace")
      method == "multiply" || .zeroDepends(as.character(events$value[r]), zs)$zero
    }, TRUE)
    unique(v[!keep])
  }
  zeroState <- intersect(names(ode), zeroPar)
  repeat {
    zs <- c(setdiff(zeroPar, names(ode)), zeroState)
    stay <- zeroState[vapply(zeroState, function(st) .zeroDepends(ode[[st]], zs)$zero, TRUE)]
    stay <- setdiff(stay, moves(zs))
    if (setequal(stay, zeroState)) break
    zeroState <- stay
  }
  zs <- c(setdiff(zeroPar, names(ode)), zeroState)
  for (st in setdiff(names(ode), zeroState)) {
    sy <- .zeroDepends(ode[[st]], zs)$symbols
    edge(inner, sy, inner, rep(st, length(sy)))
  }
  for (st in names(xd)) edge(inner, xd[[st]], inner, rep(st, length(xd[[st]])))
  if (!is.null(events))
    for (r in seq_len(nrow(events))) {
      st <- as.character(events$var[r])
      if (st %in% zeroState) next
      sy <- unique(unlist(lapply(setdiff(names(events), c("var", "method")), function(col)
        .zeroDepends(as.character(events[[col]][r]), zs)$symbols)))
      edge(inner, sy, inner, rep(st, length(sy)))
    }
  # the observables, the first reading the states and inner parameters, every further one the
  # outputs of the one before and the inner parameters
  top <- inner
  zeroTop <- zs
  for (leaf in rev(lv[seq_len(at - 1L)])) {
    obs <- attr(leaf$kernel, "equations", exact = TRUE)
    if (is.null(obs)) return(every)
    below <- namesAt(top)
    zout <- character(0)
    for (o in names(obs)) {
      d <- .zeroDepends(obs[[o]], unique(c(zeroTop, zs)))
      if (d$zero) zout <- c(zout, o)
      lower <- intersect(d$symbols, below)
      edge(top, lower, top + 1L, rep(o, length(lower)))
      edge(inner, setdiff(d$symbols, lower), top + 1L, rep(o, length(setdiff(d$symbols, lower))))
      hold(top + 1L, o)
    }
    if (isTRUE(environment(leaf$kernel)$controls$attachInput)) {
      through <- setdiff(below, names(obs))
      edge(top, through, top + 1L, through)
      zout <- c(zout, intersect(through, zeroTop))
    }
    zeroTop <- zout
    top <- top + 1L
  }
  d <- design$data[[cn]]
  measured <- unique(as.character(d$name))
  byModel <- unique(as.character(d$name[if (is.null(d$sigma)) TRUE else is.na(d$sigma)]))
  sinks <- node(top, measured)
  if (length(byModel)) {
    err <- if (!is.null(design$errfn)) .fnLeaves(design$errfn, cn)
    if (!length(err) || any(vapply(err, `[[`, "", "kind") != "obsfn")) return(every)
    # the error model reads the observables and the inner parameters
    base <- top
    for (leaf in rev(err)) {
      eq <- attr(leaf$kernel, "equations", exact = TRUE)
      if (is.null(eq)) return(every)
      below <- namesAt(base)
      for (o in names(eq)) {
        sy <- .zeroDepends(eq[[o]], unique(c(zeroTop, zs)))$symbols
        lower <- intersect(sy, below)
        edge(base, lower, base + 1L, rep(o, length(lower)))
        edge(inner, setdiff(sy, lower), base + 1L, rep(o, length(setdiff(sy, lower))))
      }
      base <- base + 1L
    }
    sinks <- c(sinks, node(base, byModel))
  }
  setNames(graph_reach(as.character(unlist(from)), as.character(unlist(to)), node(0L, etas), sinks), etas)
}

## The roots of a steady state (environment `env` of Pimpl()) each input moves: the states whose
## equation holds the input, then every state whose equation holds a moved state and every state of
## a conservation row with the input or a moved state. The other states solve a subsystem that holds
## neither the input nor a moved state; where it lacks full structural column rank, the input moves
## every root. Inputs in `zeros` are 0.
.pimplDepends <- function(env, zeros) {
  dep <- env$dependent
  rows <- env$all_exprs
  syms <- lapply(rows, function(r) .zeroDepends(r, zeros)$symbols)
  cons <- !(names(rows) %in% dep)
  inputs <- setdiff(unique(unlist(syms)), dep)
  inputs <- intersect(inputs, env$parms_all)
  setNames(lapply(inputs, function(pin) {
    moved <- character(0)
    repeat {
      hit <- vapply(syms, function(s) pin %in% s || any(moved %in% s), TRUE)
      add <- unique(c(names(rows)[hit & !cons], unlist(lapply(syms[hit & cons], intersect, dep))))
      add <- intersect(add, dep)
      if (all(add %in% moved)) break
      moved <- union(moved, add)
    }
    rest <- setdiff(dep, moved)
    keep <- (!cons & names(rows) %in% rest) | (cons & vapply(syms, function(s) !any(intersect(s, dep) %in% moved), TRUE))
    if (length(rest) && .structuralRank(syms[keep], rest) < length(rest)) dep else moved
  }), inputs)
}

## The size of a maximum matching between the rows, each the symbols it holds, and `cols`.
.structuralRank <- function(rows, cols) {
  owner <- setNames(rep(NA_integer_, length(cols)), cols)
  adj <- lapply(rows, intersect, cols)
  seen <- NULL
  augment <- function(r) {
    for (c in adj[[r]]) {
      if (seen[[c]]) next
      seen[[c]] <<- TRUE
      if (is.na(owner[[c]]) || augment(owner[[c]])) {
        owner[[c]] <<- r
        return(TRUE)
      }
    }
    FALSE
  }
  for (r in seq_along(adj)) {
    seen <- setNames(rep(FALSE, length(cols)), cols)
    augment(r)
  }
  sum(!is.na(owner))
}

## The symbols an expression holds once the names in `zeros` are 0, and whether it is then 0.
.zeroDepends <- function(expr, zeros) {
  e <- .zeroOut(str2lang(as.character(expr)), zeros)
  list(zero = is.numeric(e) && length(e) == 1L && !is.na(e) && e == 0,
       symbols = if (is.numeric(e)) character(0) else setdiff(all.vars(e), "time"))
}

## An expression with the names in `zeros` set to 0: a product with a factor 0, a quotient with
## numerator 0 and a positive power of 0 vanish, sums drop their terms 0, and an elementary function
## of constants is evaluated. Returns a number where the result is constant.
.zeroOut <- function(e, zeros) {
  if (is.numeric(e) || is.logical(e)) return(as.numeric(e))
  if (is.name(e)) return(if (as.character(e) %in% zeros) 0 else e)
  if (!is.call(e) || !is.name(e[[1]])) return(e)
  f <- as.character(e[[1]])
  a <- lapply(as.list(e)[-1], .zeroOut, zeros = zeros)
  zero <- function(x) is.numeric(x) && length(x) == 1L && !is.na(x) && x == 0
  const <- function(x) is.numeric(x) && length(x) == 1L && !is.na(x)
  if (f == "(") return(a[[1]])
  if (f == "*" && any(vapply(a, zero, TRUE))) return(0)
  if (f == "/" && zero(a[[1]]) && !zero(a[[2]])) return(0)
  if (f == "^" && zero(a[[1]]) && const(a[[2]]) && a[[2]] > 0) return(0)
  if (f == "+") {
    a <- Filter(Negate(zero), a)
    if (!length(a)) return(0)
    if (length(a) == 1L) return(a[[1]])
  }
  if (f == "-") {
    if (all(vapply(a, zero, TRUE))) return(0)
    if (length(a) == 2L && zero(a[[2]])) return(a[[1]])
    if (length(a) == 2L && zero(a[[1]])) a <- a[2]
  }
  out <- as.call(c(e[[1]], a))
  if (f %in% .elementary && all(vapply(a, const, TRUE))) {
    v <- tryCatch(suppressWarnings(eval(out, baseenv())), error = function(err) NULL)
    if (const(v)) return(v)
  }
  out
}

.elementary <- c("+", "-", "*", "/", "^", "exp", "log", "log10", "log2", "log1p", "expm1", "sqrt",
                 "abs", "sin", "cos", "tan", "sinh", "cosh", "tanh", "asin", "acos", "atan",
                 "asinh", "acosh", "atanh")

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


## One distribution block: the stems (random-effect columns) and short names, the distribution
## and the fixed hyperparameters, NA when estimated.
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
## intercept, the first of the symbol, whose formula covers every condition with data and whose
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

## Generalized least squares for the referenced fixed effects: with `phi` the posterior means per
## group (columns the normal stems, phi in the referenced ones), the free entries of `beta` that
## minimize sum_g (phi_g - D_g beta)' Omega^-1 (phi_g - D_g beta).
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

## The standard deviation of the perturbation of every parameter of `center` from msem(sd = ): one
## number for all, or a vector named by the parameters it perturbs, 0 for the others.
.perturbSd <- function(sd, center) {
  if (!is.numeric(sd) || !length(sd) || anyNA(sd) || any(sd < 0))
    stop("msem: `sd` is a non-negative number or a vector of them named by parameters.", call. = FALSE)
  if (length(sd) == 1L && is.null(names(sd))) return(setNames(rep(sd, length(center)), names(center)))
  if (is.null(names(sd)) || any(names(sd) == ""))
    stop("msem: `sd` with several entries names the parameters it perturbs.", call. = FALSE)
  unknown <- setdiff(names(sd), names(center))
  if (length(unknown))
    stop("msem: `sd` names ", paste(unknown, collapse = ", "), ", which the fit does not estimate.", call. = FALSE)
  out <- setNames(numeric(length(center)), names(center))
  out[names(sd)] <- sd
  out
}

## A start perturbed by normal noise of standard deviations `sd` on the scale where each entry is
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
