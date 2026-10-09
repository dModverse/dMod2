# scanL1(lambda = "em"): the strength of the Lq penalty estimated by an EM on
# the marginal likelihood, one strength per family (gates, reference
# parameters). The penalty lambda |b|^q is the exponential power prior
# p(b) = q lambda^(1/q) / (2 Gamma(1/q)) exp(-lambda |b|^q), one-sided for gates.

# -2 log of the normalising constant of the prior.
.l1EmLogNorm <- function(lambda, q, oneSided)
  -2 * (log(q) + log(lambda) / q - lgamma(1 / q) - if (oneSided) 0 else log(2))

# Posterior of one term from the data around the MAP value `b` (gradient `g`
# and precision `c` of -log L) and the prior: log of the integral without the
# prior constant, and E|b|^q.
.l1EmTerm <- function(b, g, c, lambda, q, oneSided) {
  s0 <- lambda^(-1 / q)
  # Data without curvature on the scale of the prior: the posterior is the prior.
  if (c * s0^2 < 1e-8) { c <- 0; g <- 0 }
  m  <- if (c > 0) b - g / c else b
  w  <- if (c > 0) min(10 / sqrt(c), 1e3 * s0 + abs(m)) else 1e3 * s0 + abs(m)
  lo <- if (oneSided) 0 else min(m - w, -20 * s0)
  hi <- max(m + w, 20 * s0)
  f  <- function(x) -(g * (x - b) + 0.5 * c * (x - b)^2) - lambda * abs(x)^q
  sd <- if (c > 0) min(1 / sqrt(c), w) else w
  cut <- c(lo, hi, 0, m + c(-1, 1) * 8 * sd, b + c(-1, 1) * 8 * sd)
  cut <- sort(unique(cut[cut >= lo & cut <= hi]))
  grid <- c(cut, seq(lo, hi, length.out = 401), m + seq(-8, 8, by = 0.25) * sd)
  shift <- max(f(grid[grid >= lo & grid <= hi]))
  piece <- function(k, a, z)
    stats::integrate(function(x) abs(x)^k * exp(f(x) - shift), a, z,
                     subdivisions = 500L, rel.tol = 1e-8, stop.on.error = FALSE)$value
  Z  <- sum(vapply(seq_len(length(cut) - 1L), function(i) piece(0, cut[i], cut[i + 1L]), 0))
  Eq <- sum(vapply(seq_len(length(cut) - 1L), function(i) piece(q, cut[i], cut[i + 1L]), 0)) / Z
  list(logZ = log(Z) + shift, Eq = Eq)
}

# MAP of the Lq penalised objective at strengths `lam` (named per term), by
# reweighted L1 for q < 1.
.l1EmMap <- function(obj, st, lam, q, gates, fixed, ctl, scale) {
  terms <- names(lam)
  run <- function(st, w) {
    args <- c(list(obj, st), ctl$trust,
              list(fixed = fixed, mu = stats::setNames(rep(0, length(terms)), terms),
                   lambda = 2 * lam * w / scale, gate = if (length(gates)) gates))
    do.call(trustL1, args)
  }
  w <- stats::setNames(rep(1, length(terms)), terms)
  if (q < 1) w <- q * pmax(abs(st[terms] / scale), ctl$eps)^(q - 1)
  fit <- run(st, w)
  if (q < 1) for (k in seq_len(ctl$nq)) {
    w <- q * pmax(abs(fit$argument[terms] / scale), ctl$eps)^(q - 1)
    nxt <- try(run(fit$argument, w), silent = TRUE)
    if (inherits(nxt, "try-error")) break
    fit <- nxt
  }
  fit
}

# E-step at a MAP fit: data gradient and Hessian, precision of every term
# conditional on the others with the unpenalised parameters profiled, and the
# posterior of every term.
.l1EmEstep <- function(obj, fit, lam, q, gates, fixed, scale) {
  terms <- names(lam)
  arg <- fit$argument
  d <- obj(arg, fixed = fixed, deriv = TRUE)
  H <- d$hessian / 2
  N <- setdiff(rownames(H), terms)
  P <- .l1EmSchur(H, terms, N)
  a <- scale[terms]
  post <- lapply(terms, function(j)
    .l1EmTerm(arg[[j]] / a[[j]], d$gradient[[j]] / 2 * a[[j]], max(P[j, j], 0) * a[[j]]^2,
              lam[[j]], q, j %in% gates))
  names(post) <- terms
  list(data = d$value, logZ = vapply(post, `[[`, 0, "logZ") + log(a),
       Eq = vapply(post, `[[`, 0, "Eq"), precision = diag(P))
}

.l1EmSchur <- function(A, I, N) {
  if (!length(N)) return(A[I, I, drop = FALSE])
  e <- eigen(A[N, N, drop = FALSE], symmetric = TRUE)
  keep <- e$values > 1e-10 * max(abs(e$values), 1e-300)
  V <- e$vectors[, keep, drop = FALSE]
  P <- A[I, I, drop = FALSE] -
    A[I, N, drop = FALSE] %*% V %*% (t(V) / e$values[keep]) %*% A[N, I, drop = FALSE]
  (P + t(P)) / 2
}

# One EM from `st`: MAP, posterior of the terms, closed-form strengths per
# family, until the strengths settle. The prior acts on every term divided by
# its `scale`. Value: approximate -2 log marginal
# likelihood plus the -2 log Gamma hyperprior of the strengths.
.l1EmRun <- function(obj, st, family, q, gates, fixed, ctl, scale) {
  em <- ctl$em
  fams <- unique(family)
  lamF <- stats::setNames(rep(em$init, length(fams)), fams)
  trace <- NULL
  for (it in seq_len(ctl$nem)) {
    lam <- stats::setNames(lamF[family], names(family))
    fit <- .l1EmMap(obj, st, lam, q, gates, fixed, ctl, scale)
    es  <- .l1EmEstep(obj, fit, lam, q, gates, fixed, scale)
    new <- vapply(fams, function(f) {
      j <- family == f
      (sum(j) / q + em$a - 1) / (sum(es$Eq[j]) + em$b)
    }, 0)
    value <- es$data - 2 * sum(es$logZ) +
      sum(vapply(names(family), function(j)
        .l1EmLogNorm(lam[[j]], q, j %in% gates), 0)) -
      2 * sum((em$a - 1) * log(lamF) - em$b * lamF)
    trace <- rbind(trace, c(iteration = it, value = value, lamF))
    st <- fit$argument
    done <- max(abs(log(new) - log(lamF))) < em$tol
    lamF <- new
    if (done) break
  }
  fit$value  <- value
  fit$lambdaEM <- stats::setNames(trace[nrow(trace), fams, drop = TRUE], fams)
  fit$emTrace  <- trace
  fit$emIterations <- it
  fit$terms <- data.frame(term = names(family), family = unname(family),
                          scale = unname(scale[names(family)]),
                          estimate = unname(fit$argument[names(family)]),
                          Eq = unname(es$Eq), se = unname(1 / sqrt(es$precision)),
                          stringsAsFactors = FALSE)
  fit
}

# The EM variant of scanL1(): one multistart of the EM, then refits without
# penalty along the terms ordered by their size at the best run.
.l1EmScan <- function(obj, start0, sparse, gates, reference, fixSel, q, ctl, fits, sd,
                      cores, wf, full, fullPars, zero, fixed, alpha) {
  family <- c(stats::setNames(rep("gate", length(gates)), gates),
              stats::setNames(rep("reference", length(reference)), reference))
  # Adaptive: a reference parameter is penalised relative to its full estimate,
  # as a gate is relative to its full rate.
  scale <- stats::setNames(rep(1, length(family)), names(family))
  if (isTRUE(ctl$em$adaptive) && length(reference))
    scale[reference] <- pmax(abs(fullPars[reference]), ctl$eps)
  best <- .l1Multistart(function(st) .l1EmRun(obj, st, family, q, gates, fixSel, ctl, scale),
                        start0, fits, sd, cores, obj, extra = list(sparse),
                        positive = gates, wf = wf, levelTol = ctl$tolHits)
  if (is.null(best)) stop("scanL1: every EM run failed.", call. = FALSE)
  # Refits along the terms ordered by their penalised size: from the MAP
  # structure, terms are removed while the test against the full model does
  # not reject, and added back while it does.
  terms <- names(family)
  u <- abs(best$argument[terms]) / scale[terms]
  ord <- terms[order(-u)]
  nOn <- sum(u > 0)
  nFull <- length(fullPars)
  refit <- function(k) {
    th <- best$argument
    th[setdiff(terms, ord[seq_len(k)])] <- 0
    th[intersect(ord[seq_len(k)], terms[th[terms] == 0])] <- 1
    st <- .l1Structure(th, gates, reference, list())
    r <- .l1Refit(st, best$argument, obj, fullPars, zero, fixed, ctl, fits, sd, cores, wf)
    if (is.null(r)) return(NULL)
    df <- nFull - r$nfree
    stat <- max(0, r$value - full$value)
    list(st = st, r = r, row = data.frame(key = st$key, value = r$value, nfree = r$nfree,
         starts = r$starts %||% NA_integer_, hits = r$hits %||% NA_integer_, stat = stat,
         df = df, p = if (df > 0) stats::pchisq(stat, df, lower.tail = FALSE) else 1,
         bic = NA_real_, size = k, stringsAsFactors = FALSE))
  }
  cache <- list()
  get <- function(k) {
    kk <- as.character(k)
    if (is.null(cache[[kk]])) cache[[kk]] <<- refit(k)
    cache[[kk]]
  }
  pass <- function(k) { g <- get(k); !is.null(g) && g$row$p >= alpha }
  k <- nOn
  if (pass(k)) {
    while (k > 0L && pass(k - 1L)) k <- k - 1L
  } else {
    while (k < length(terms) && !pass(k)) k <- k + 1L
  }
  done <- Filter(Negate(is.null), cache)
  if (!length(done)) stop("scanL1: every refit failed.", call. = FALSE)
  refitTab <- do.call(rbind, lapply(done, `[[`, "row"))
  refitTab <- refitTab[order(refitTab$size), ]
  rownames(refitTab) <- NULL
  selR <- get(k)
  key <- if (is.null(selR)) "" else selR$st$key
  structs <- stats::setNames(lapply(done, `[[`, "st"), vapply(done, function(d) d$st$key, ""))
  mapKey <- .l1Structure(best$argument, gates, reference, list())$key
  lev <- c(list(best), best$level)
  keys <- vapply(lev, function(f) .l1Structure(f$argument, gates, reference, list())$key, "")
  pathTab <- data.frame(lambda = NA_real_, value = best$value, key = mapKey,
                        removed = sum(u == 0),
                        p = refitTab$p[match(mapKey, refitTab$key)],
                        converged = isTRUE(best$converged), starts = best$starts,
                        hits = best$hits, em = best$emIterations, stringsAsFactors = FALSE)
  out <- list(path = pathTab, coefficients = NULL, arguments = list(best$argument),
              refits = refitTab,
              level = data.frame(key = keys, value = vapply(lev, `[[`, 0, "value")),
              levelFits = list(lapply(lev, `[`, c("argument", "value"))),
              full = list(value = full$value, argument = fullPars, values = full$values,
                          starts = full$starts, hits = full$hits),
              selected = key, lambdaSelected = NA_real_,
              structure = structs, select = "lrt", alpha = alpha, q = q,
              fit = if (is.null(selR)) fullPars else selR$r$argument,
              groups = list(), gates = gates, reference = reference,
              em = list(lambda = best$lambdaEM, terms = best$terms, trace = best$emTrace,
                        values = best$values, scale = scale))
  class(out) <- "scanL1"
  out
}
