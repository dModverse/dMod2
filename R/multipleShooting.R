## Multiple shooting ------------------------------------------------------------
## Segment j starts from node s_j, continuity c_j = phi(x_j(tau_{j+1})) - s_{j+1}
## in the chart phi. Nodes enter as inner initial states; nothing recompiles.

## Copyright (C) 2026 Simon Beyer


# obj(pars, nodes, ...) returning the data term, the gaps and, per segment, the
# blocks on (theta_loc, s_j) and the end state's chart value with its Jacobian.
# Segment j holds the data and fixed-time events in [tau_j, tau_{j+1}).
.shootObjective <- function(data, x, p, g = NULL, errmodel = NULL, nodes,
                            charts = NULL, times = NULL,
                            attr.name = "data", opt.BLOQ = "M3") {

  if (!inherits(x, "prdfn")) stop("multiple shooting: 'x' must be a prediction function.", call. = FALSE)
  if (!inherits(p, "parfn")) stop("multiple shooting: 'p' must be a parameter transformation.", call. = FALSE)
  if (!is.null(g) && !inherits(g, "obsfn"))
    stop("multiple shooting: 'g' must be an observation function or NULL.", call. = FALSE)

  states <- .shootStates(x)
  charts <- .shootCharts(charts, states)

  x.cond <- names(attr(x, "mappings"))
  d.cond <- names(data)
  conditions <- if (is.null(x.cond) || all(!nzchar(x.cond))) d.cond else intersect(d.cond, x.cond)
  if (!length(conditions)) stop("multiple shooting: no condition of the data is predicted by 'x'.", call. = FALSE)
  e.cond <- if (!is.null(errmodel)) names(attr(errmodel, "mappings")) else NULL

  nodesL <- .shootPerCondition(nodes, conditions, "nodes")
  span <- .shootSpan(data, times, conditions)

  # Segments of every condition, in one flat list: the order every evaluation
  # and the driver index by.
  segs <- list()
  for (cn in conditions) {
    dc <- data[[cn]]
    tstart <- span[[cn]][1L]
    tend <- span[[cn]][2L]
    tcn <- .shootTimesOf(times, cn)
    tau <- .shootSnap(sort(unique(c(tstart, as.numeric(nodesL[[cn]])))), c(dc$time, tcn))
    tau <- tau[!duplicated(tau)]
    if (any(tau < tstart))
      stop("multiple shooting: condition '", cn, "' has nodes before its start, t = ",
           tstart, ".", call. = FALSE)
    if (any(tau[-1L] >= tend))
      stop("multiple shooting: condition '", cn, "' has nodes at or after its last data ",
           "point, t = ", tend, "; every segment needs time to run.", call. = FALSE)
    M <- length(tau)
    upper <- c(tau[-1L], tend)
    for (j in seq_len(M)) {
      last <- j == M
      inseg <- dc$time >= tau[j] & (if (last) dc$time <= upper[j] else dc$time < upper[j])
      extra <- tcn[tcn >= tau[j] & tcn <= upper[j]]
      key <- paste0(cn, "#", j)
      segs[[length(segs) + 1L]] <- list(
        cond = cn, j = j, M = M, key = key,
        start = tau[j], end = upper[j], last = last,
        times = sort(unique(c(tau[j], dc$time[inseg], extra, upper[j]))),
        data = dc[inseg, , drop = FALSE],
        nodeNames = if (j > 1L) paste0(states, "@", key) else character(0))
    }
  }
  keys <- vapply(segs, `[[`, "", "key")
  dataSeg <- setNames(lapply(segs, `[[`, "data"), keys)
  ekeys <- if (is.null(e.cond)) NULL else keys[vapply(segs, function(s) s$cond %in% e.cond, TRUE)]

  spec <- list2env(list(
    x = x, p = p, g = g, errmodel = errmodel, states = states, charts = charts,
    conditions = conditions, e.cond = e.cond, segs = segs, keys = keys,
    dataSeg = dataSeg, ekeys = ekeys, opt.BLOQ = opt.BLOQ,
    tau = lapply(setNames(conditions, conditions), function(cn)
      vapply(Filter(function(s) s$cond == cn, segs), `[[`, 0, "start"))),
    parent = emptyenv())
  meta_cache <- new.env(parent = emptyenv())

  myfn <- function(pars, nodes, fixed = NULL, deriv = TRUE, hessian = deriv,
                   sweep = c("forward", "reverse"),
                   cores = getOption("dMod.cores", 1L), env = NULL,
                   growth = FALSE, residuals = FALSE) {
    sweep <- match.arg(sweep)
    if (is.null(env)) env <- new.env()
    .shootCheckNodes(nodes, spec)
    if (identical(sweep, "reverse") && isTRUE(deriv))
      .shootReverse(spec, meta_cache, pars, nodes, fixed, cores, env, attr.name)
    else
      .shootForward(spec, meta_cache, pars, nodes, fixed, isTRUE(deriv),
                 isTRUE(deriv) && isTRUE(hessian), cores, env, attr.name,
                 growth = isTRUE(growth), residuals = isTRUE(residuals))
  }

  class(myfn) <- "shootingobj"
  attr(myfn, "conditions") <- conditions
  attr(myfn, "parameters") <- attr(p, "parameters")
  attr(myfn, "states")     <- states
  attr(myfn, "charts")     <- charts
  attr(myfn, "nodes")      <- spec$tau
  attr(myfn, "modelname")  <- modelname(x, g, p, errmodel)
  attr(myfn, "compileInfo") <- Reduce(.mergeCompileInfo, lapply(
    list(g, x, p, errmodel), function(f) attr(f, "compileInfo")))
  attr(myfn, "spec") <- spec
  # what a new layout of the nodes is built from
  attr(myfn, "args") <- list(data = data, x = x, p = p, g = g,
                             errmodel = errmodel, charts = charts,
                             times = times, attr.name = attr.name,
                             opt.BLOQ = opt.BLOQ)
  myfn
}


## ---- Setup helpers -----------------------------------------------------------

# The state names of a prediction function from Xs(): the equations of its ODE.
.shootStates <- function(x) {
  maps <- attr(x, "mappings")
  eq <- if (length(maps)) attr(maps[[1L]], "equations") else NULL
  if (is.null(eq))
    stop("multiple shooting: 'x' has no equations; it has to come from Xs().", call. = FALSE)
  names(eq)
}

.shootCharts <- function(charts, states) {
  out <- setNames(rep("linear", length(states)), states)
  if (is.null(charts)) return(out)
  if (is.null(names(charts)) || !all(nzchar(names(charts))))
    stop("multiple shooting: 'charts' must be named by state.", call. = FALSE)
  bad <- setdiff(names(charts), states)
  if (length(bad))
    stop("multiple shooting: 'charts' names ", paste(bad, collapse = ", "),
         ", which are not states of 'x'.", call. = FALSE)
  if (!all(charts %in% c("log10", "linear", "angle")))
    stop("multiple shooting: a chart is \"log10\", \"linear\" or \"angle\".", call. = FALSE)
  out[names(charts)] <- charts
  out
}

# A value for all conditions, or a list named by condition.
.shootPerCondition <- function(v, conditions, what) {
  if (is.list(v)) {
    miss <- setdiff(conditions, names(v))
    if (length(miss))
      stop("multiple shooting: '", what, "' is a list without entries for ",
           paste(miss, collapse = ", "), ".", call. = FALSE)
    return(v[conditions])
  }
  setNames(rep(list(v), length(conditions)), conditions)
}

# Node times within a relative 1e-9 of a grid time become that time: the
# solver merges times that close, and a segment's grid has to keep its length.
.shootSnap <- function(tau, grid) {
  if (!length(grid)) return(tau)
  i <- findInterval(tau, sort(unique(grid)), all.inside = FALSE)
  g <- sort(unique(grid))
  for (d in c(0L, 1L)) {
    j <- pmin(pmax(i + d, 1L), length(g))
    hit <- abs(tau - g[j]) <= 1e-9 * pmax(1, abs(tau))
    tau[hit] <- g[j][hit]
  }
  tau
}

# The extra times of condition `cn`: normL2's `times` for all conditions, or
# its entry of the list.
.shootTimesOf <- function(times, cn)
  as.numeric(if (is.list(times)) times[[cn]] else times)

# Start and end of the time axis of every condition: the start of its grid in
# normL2, and its last data point or extra time.
.shootSpan <- function(data, times, conditions) {
  grid <- .normL2Grid(data, times)
  setNames(lapply(conditions, function(cn) {
    ts <- .gridOf(grid, cn)[1L]
    tcn <- .shootTimesOf(times, cn)
    c(ts, max(c(data[[cn]]$time, tcn[tcn >= ts], ts)))
  }), conditions)
}

.shootCheckNodes <- function(nodes, spec) {
  if (!is.list(nodes) || !all(spec$conditions %in% names(nodes)))
    stop("multiple shooting: 'nodes' must be a list with one matrix per condition, ",
         "as a previous fit returns them.", call. = FALSE)
  for (cn in spec$conditions) {
    nd <- nodes[[cn]]
    need <- length(spec$tau[[cn]]) - 1L
    if (!is.matrix(nd) || nrow(nd) != need || !all(spec$states %in% colnames(nd)))
      stop("multiple shooting: the nodes of condition '", cn, "' must be a matrix with ",
           need, " rows and a column per state.", call. = FALSE)
  }
  invisible(NULL)
}


## ---- Charts ------------------------------------------------------------------

# State value from a node coordinate, and its derivative.
.shootFromChart <- function(s, charts) {
  lg <- charts == "log10"
  v <- s; d <- rep(1, length(s))
  v[lg] <- 10^s[lg]
  d[lg] <- log(10) * v[lg]
  list(value = setNames(v, names(charts)), deriv = setNames(d, names(charts)))
}

# Node coordinate of a state value, and its derivative.
.shootToChart <- function(v, charts, where = "") {
  lg <- charts == "log10"
  if (any(lg & !(v > 0)))
    stop("multiple shooting: state ", paste(names(charts)[lg & !(v > 0)], collapse = ", "),
         " is not positive", where, ", so its log10 chart is undefined.",
         call. = FALSE)
  s <- v; d <- rep(1, length(v))
  s[lg] <- log10(v[lg])
  d[lg] <- 1 / (log(10) * v[lg])
  list(value = setNames(s, names(charts)), deriv = setNames(d, names(charts)))
}


## ---- Inner parameters of the segments ----------------------------------------

# Inner parameters of the segments `which`: p's own for a condition's first
# segment, node values otherwise. With derivatives each solve integrates theta
# and its node (with `x0cols` also the initial state) as directions.
.shootSegmentParsAll <- function(pin, spec, nodes, deriv, x0cols = FALSE,
                                 which = seq_along(spec$segs)) {
  states <- spec$states
  nx <- length(states)
  lg <- spec$charts[states] == "log10"
  out <- vector("list", length(which))
  segs <- spec$segs[which]
  conds <- vapply(segs, `[[`, "", "cond")
  for (cn in unique(conds)) {
    pos <- which(conds == cn)
    p0 <- pin[[cn]]
    v0 <- unclass(p0)
    attributes(v0) <- NULL
    names(v0) <- names(p0)
    si <- match(states, names(v0))
    D <- if (deriv) attr(p0, "deriv") else NULL
    if (deriv && is.null(D)) D <- matrix(0, 0, 0, dimnames = list(character(0), character(0)))
    tmpl <- NULL
    if (deriv) {
      rows <- names(v0)[names(v0) %in% union(rownames(D), states)]
      tmpl <- matrix(0, length(rows), ncol(D) + nx)
      keep <- intersect(setdiff(rownames(D), states), rows)
      if (length(keep)) tmpl[match(keep, rows), seq_len(ncol(D))] <- D[keep, , drop = FALSE]
      sr <- cbind(match(states, rows), ncol(D) + seq_len(nx))
      fixedRows <- if (length(rows) < length(v0)) setdiff(names(v0), rows) else NULL
    }
    for (i in pos) {
      sg <- segs[[i]]
      if (sg$j == 1L) {
        if (!deriv) { q <- p0; attr(q, "deriv") <- NULL; out[[i]] <- q; next }
        if (!x0cols) { out[[i]] <- p0; next }
        Dn <- matrix(0, length(rows), ncol(D) + nx,
                     dimnames = list(rows, c(colnames(D), paste0(states, "@x0"))))
        dr <- intersect(rownames(D), rows)
        Dn[match(dr, rows), seq_len(ncol(D))] <- D[dr, , drop = FALSE]
        Dn[sr] <- 1
        out[[i]] <- .shootParvec(v0, Dn, fixedRows)
        next
      }
      sv <- nodes[[cn]][sg$j - 1L, states]
      v <- v0
      v[si] <- ifelse(lg, 10^sv, sv)
      if (!deriv) { out[[i]] <- .shootParvec(v, NULL, NULL); next }
      Dn <- tmpl
      Dn[sr] <- ifelse(lg, log(10) * v[si], 1)
      dimnames(Dn) <- list(rows, c(colnames(D), sg$nodeNames))
      out[[i]] <- .shootParvec(v, Dn, fixedRows)
    }
  }
  out
}

# A parvec as as.parvec() builds it, without its checks.
.shootParvec <- function(v, D, fixed) {
  attr(v, "deriv") <- D
  attr(v, "fixed") <- fixed
  class(v) <- c("parvec", "numeric")
  v
}

# p for every condition in one call.
.shootInner <- function(spec, pars, fixed, cores) {
  out <- spec$p(pars, fixed = fixed, deriv = TRUE, conditions = spec$conditions,
                cores = cores)
  out <- out[spec$conditions]
  miss <- spec$conditions[vapply(out, is.null, TRUE)]
  if (length(miss))
    stop("multiple shooting: 'p' returns nothing for condition ", paste(miss, collapse = ", "),
         ".", call. = FALSE)
  out
}

# A segment has to run from its node to its end, which the grid encodes. A
# model that forces t = 0 or an event time into the grid moves the start.
.shootCheckGrid <- function(pr, seg, k = NA_integer_) {
  tt <- pr[, "time"]
  if (!length(tt) || tt[1L] != seg$start || tt[length(tt)] != seg$end ||
      length(tt) != length(seg$times)) {
    why <- if (length(tt) && tt[1L] < seg$start && tt[1L] == 0)
      paste0(" The ODE backend integrated from t = 0: build the model with ",
             "odemodel(..., includeTimeZero = FALSE).")
    else if (length(tt) && tt[length(tt)] < seg$end)
      paste0(" The solver stopped at t = ", tt[length(tt)], ".")
    else " The backend extended the grid; fixed-time events need cppDE >= 0.10.2."
    msg <- paste0("multiple shooting: segment ", seg$j, " of condition '", seg$cond,
                  "' should run on [", seg$start, ", ", seg$end,
                  "], the solver returned [", if (length(tt)) tt[1L] else NA, ", ",
                  if (length(tt)) tt[length(tt)] else NA, "].", why)
    # classed, so the driver can split the segment that failed
    stop(structure(class = c("shootingSegmentError", "error", "condition"),
                   list(message = msg, call = NULL, segment = k)))
  }
  invisible(NULL)
}

# The error model on the observations of every segment, as normL2 evaluates it.
.shootErrors <- function(spec, obs, fixed, deriv, cores) {
  if (is.null(spec$errmodel)) return(NULL)
  segs <- spec$segs
  idx <- which(vapply(segs, function(s) is.null(spec$e.cond) || s$cond %in% spec$e.cond, TRUE))
  if (!length(idx)) return(NULL)
  split <- lapply(idx, function(k) {
    pinner <- getParameters(obs[[k]])
    fixedinner <- pinner[union(attr(pinner, "fixed"),
                               intersect(names(pinner), names(fixed)))]
    list(pars  = as.parvec(pinner[setdiff(names(pinner), names(fixed))]),
         fixed = as.parvec(fixedinner, deriv = FALSE, deriv2 = FALSE))
  })
  conds <- vapply(segs[idx], `[[`, "", "cond")
  b <- .bundle(conds = conds, out = obs[idx],
               pars = lapply(split, `[[`, "pars"),
               fixed = lapply(split, `[[`, "fixed"), shared = FALSE)
  got <- .evalMany(spec$errmodel, b, deriv, FALSE, NULL, cores)
  err <- vector("list", length(segs))
  err[idx] <- got
  list(err = err, idx = idx, split = split)
}

.shootMeta <- function(spec, cache, obs, err) {
  sig <- lapply(obs, function(pr) dimnames(attr(pr, "deriv"))[[3]])
  shape <- c(vapply(obs, NROW, integer(1)), vapply(err, NROW, integer(1)))
  if (is.null(cache$meta) || !identical(cache$sig, sig) || !identical(cache$shape, shape)) {
    cache$meta <- .build_normL2_meta(spec$dataSeg, setNames(obs, spec$keys),
                                     err, spec$keys, spec$ekeys)
    cache$sig <- sig
    cache$shape <- shape
  }
  cache$meta
}

# The observations of every segment: g on the states, or the states themselves.
.shootObserve <- function(spec, preds, deriv, cores, env) {
  if (is.null(spec$g)) return(list(values = preds, b = NULL))
  hs <- lapply(preds, function(v) .handoff_prd_innerfixed(v, NULL))
  conds <- vapply(spec$segs, `[[`, "", "cond")
  b <- .bundle(conds = conds, out = preds, pars = lapply(hs, `[[`, "pars"),
               fixed = lapply(hs, `[[`, "fixed"), shared = FALSE)
  list(values = .evalMany(spec$g, b, deriv, FALSE, env, cores), b = b)
}

# Chart value of every segment's end state; with `deriv`, its Jacobian on the
# segment's own columns.
.shootEnds <- function(spec, preds, deriv) {
  states <- spec$states
  lapply(seq_along(spec$segs), function(k) {
    seg <- spec$segs[[k]]
    if (seg$last) return(NULL)
    pr <- preds[[k]]
    r <- nrow(pr)
    xe <- pr[r, states]
    # a state of a log10 chart that the solve drives to zero or below, as a
    # gating variable can under wild parameters, has no end value in its
    # chart: classed like a failed solve, so the driver cuts the segment
    bad <- spec$charts[states] == "log10" & !(xe > 0)
    if (any(bad))
      stop(structure(class = c("shootingSegmentError", "error", "condition"),
                     list(message = paste0(
                       "multiple shooting: state ", paste(states[bad], collapse = ", "),
                       " is not positive at the end of segment ", seg$j,
                       " of condition '", seg$cond, "', t = ", seg$end,
                       ", so its log10 chart is undefined."),
                       call = NULL, segment = k)))
    ch <- .shootToChart(xe, spec$charts)
    J <- NULL
    if (deriv) {
      dX <- attr(pr, "deriv")
      J <- dX[r, states, , drop = TRUE]
      if (!is.matrix(J))
        J <- matrix(J, nrow = length(states),
                    dimnames = list(states, dimnames(dX)[[3]]))
      J <- J * ch$deriv
    }
    list(value = ch$value, jac = J)
  })
}

# The gaps end - node per condition; on an angle chart taken modulo 2 pi into
# [-pi, pi), so that a trajectory one turn further on is continuous.
.shootGaps <- function(spec, ends, nodes) {
  out <- list()
  ang <- spec$charts[spec$states] == "angle"
  for (cn in spec$conditions) {
    ks <- which(vapply(spec$segs, `[[`, "", "cond") == cn)
    M <- length(ks)
    G <- matrix(0, M - 1L, length(spec$states),
                dimnames = list(spec$tau[[cn]][-1L], spec$states))
    for (i in seq_len(M - 1L))
      G[i, ] <- ends[[ks[i]]]$value - nodes[[cn]][i, spec$states]
    if (any(ang)) G[, ang] <- (G[, ang] + pi) %% (2 * pi) - pi
    out[[cn]] <- G
  }
  out
}


## ---- Forward evaluation ------------------------------------------------------

.shootForward <- function(spec, cache, pars, nodes, fixed, deriv, hessian, cores,
                       env, attr.name, growth = FALSE, residuals = FALSE) {
  segs <- spec$segs
  pin <- .shootInner(spec, pars, fixed, cores)
  parsL <- .shootSegmentParsAll(pin, spec, nodes, deriv, x0cols = growth)
  conds <- vapply(segs, `[[`, "", "cond")
  preds <- .predictMany(spec$x, lapply(segs, `[[`, "times"), parsL, conds,
                        deriv = deriv, env = env, cores = cores)
  preds <- unname(unclass(preds))
  for (k in seq_along(segs)) .shootCheckGrid(preds[[k]], segs[[k]], k)

  obs <- .shootObserve(spec, preds, deriv, cores, env)$values
  ee <- .shootErrors(spec, obs, fixed, deriv, cores)
  err <- if (is.null(ee)) NULL else ee$err
  meta <- .shootMeta(spec, cache, obs, err)
  pg <- if (deriv) unique(unlist(lapply(obs, function(pr) dimnames(attr(pr, "deriv"))[[3]])))
        else character(0)
  kr <- normL2_kernel(
    prediction = setNames(obs, spec$keys), err_list_opt = err, meta_list = meta,
    par_names_global = if (is.null(pg)) character(0) else pg,
    deriv2_requested = FALSE, threads = as.integer(cores),
    bloq_mode = spec$opt.BLOQ, build_hessian = hessian,
    local_blocks = deriv)

  ends <- .shootEnds(spec, preds, deriv)
  out <- list(value = kr$value, chi2 = kr$chi2,
              gaps = .shootGaps(spec, ends, nodes),
              sweep = "forward")
  if (deriv) {
    out$segments <- lapply(seq_along(segs), function(k) {
      loc <- kr$local[[k]]
      e <- ends[[k]]
      gn <- names(loc$gradient)
      jn <- if (is.null(e)) NULL else colnames(e$jac)
      vars <- if (is.null(jn) || identical(gn, jn)) gn else union(gn, jn)
      x0 <- if (growth && segs[[k]]$j == 1L) grep("@x0$", vars, value = TRUE) else character(0)
      if (length(x0)) {
        vars <- setdiff(vars, x0)
        gn <- setdiff(gn, x0)
        jn <- setdiff(jn, x0)
      }
      G <- if (is.null(e)) NULL
           else if (segs[[k]]$j > 1L) e$jac[, segs[[k]]$nodeNames, drop = FALSE]
           else if (length(x0)) .shootChartGrowth(e$jac[, x0, drop = FALSE],
                                                  pin[[segs[[k]]$cond]], spec)
      H <- if (hessian) loc$hessian else NULL
      sg <- list(vars = vars,
           grad = if (identical(names(loc$gradient), vars)) loc$gradient
                  else .shootAlignVec(loc$gradient[gn], vars),
           hess = if (!hessian) NULL
                  else if (identical(rownames(H), vars) && identical(colnames(H), vars)) H
                  else .shootAlignMat(H, vars),
           end = if (is.null(e)) NULL else e$value,
           jac = if (is.null(e)) NULL
                 else if (identical(colnames(e$jac), vars))
                   `dimnames<-`(e$jac, list(rownames(e$jac), vars))
                 else .shootAlignCols(e$jac[, jn, drop = FALSE], vars),
           growth = G, chi2 = loc$chi2, ndata = nrow(segs[[k]]$data))
      if (residuals) {
        rs <- .shootResiduals(obs[[k]], segs[[k]]$data, TRUE)
        sg$res <- rs$r
        sg$J <- if (is.null(rs$J)) matrix(0, 0L, length(vars), dimnames = list(NULL, vars))
                else .shootAlignCols(rs$J[, intersect(colnames(rs$J), vars), drop = FALSE], vars)
      }
      sg
    })
  } else {
    out$segments <- lapply(seq_along(ends), function(k) {
      e <- ends[[k]]
      sg <- list(end = if (is.null(e)) NULL else e$value)
      if (residuals) sg$res <- .shootResiduals(obs[[k]], segs[[k]]$data, FALSE)$r
      sg
    })
  }
  attr(out, attr.name) <- out$value
  env$prediction <- setNames(obs, spec$keys)
  env$states <- setNames(preds, spec$keys)
  out$env <- env
  out
}

# Weighted residuals (prediction - data) / sigma of a segment and, with
# `deriv`, their Jacobian on the prediction's directions, for Bock's natural
# level function. Fixed sigma only.
.shootResiduals <- function(pr, dat, deriv) {
  if (!NROW(dat)) return(list(r = numeric(0), J = NULL))
  ti <- match(dat$time, pr[, "time"])
  nm <- as.character(dat$name)
  r <- (pr[cbind(ti, match(nm, colnames(pr)))] - dat$value) / dat$sigma
  J <- NULL
  if (deriv) {
    D <- attr(pr, "deriv")
    ni <- match(nm, dimnames(D)[[2]])
    J <- matrix(0, length(ti), dim(D)[3], dimnames = list(NULL, dimnames(D)[[3]]))
    for (i in seq_along(ti)) J[i, ] <- D[ti[i], ni[i], ] / dat$sigma[i]
  }
  list(r = as.numeric(r), J = J)
}

.shootAlignVec <- function(v, vars) {
  out <- setNames(numeric(length(vars)), vars)
  if (length(v)) out[names(v)] <- v
  out
}
.shootAlignMat <- function(H, vars) {
  out <- matrix(0, length(vars), length(vars), dimnames = list(vars, vars))
  if (!is.null(H) && length(H)) {
    keep <- intersect(rownames(H), vars)
    out[keep, keep] <- H[keep, keep]
  }
  out
}

# The end-state Jacobian of a first segment along its initial state, mapped
# into the chart of the nodes: d end / d s_1 = d end / d x_0 * T'(s_1).
.shootChartGrowth <- function(J, pin, spec) {
  # T'(s_1) = ln(10) x_0 on a log10 chart, which is zero, not undefined, at a
  # start of zero
  x0 <- unclass(pin)[spec$states]
  d <- ifelse(spec$charts == "log10", log(10) * pmax(x0, 0), 1)
  J * rep(d, each = nrow(J))
}
.shootAlignCols <- function(J, vars) {
  out <- matrix(0, nrow(J), length(vars), dimnames = list(rownames(J), vars))
  out[, colnames(J)] <- J
  out
}


## ---- Reverse evaluation ------------------------------------------------------
## One backward sweep per segment with 1 + n_x seeds, the data term and each end
## state; p's forward Jacobian maps the inner cotangents onto theta.

.shootReverse <- function(spec, cache, pars, nodes, fixed, cores, env, attr.name) {
  segs <- spec$segs
  states <- spec$states
  nx <- length(states)
  S <- 1L + nx
  pin <- .shootInner(spec, pars, fixed, cores)
  parsV <- .shootSegmentParsAll(pin, spec, nodes, FALSE)
  conds <- vapply(segs, `[[`, "", "cond")
  n <- length(segs)

  bx <- .bundle(conds = conds, times = lapply(segs, `[[`, "times"), pars = parsV,
                fixed = rep(list(NULL), n), shared = FALSE)
  fx <- .fwdMany(spec$x, bx, env, cores)
  preds <- unname(unclass(fx$values))
  for (k in seq_len(n)) .shootCheckGrid(preds[[k]], segs[[k]], k)

  fg <- NULL
  obs <- preds
  if (!is.null(spec$g)) {
    hs <- lapply(preds, function(v) .handoff_prd_innerfixed(v, NULL))
    bg <- .bundle(conds = conds, out = preds, pars = lapply(hs, `[[`, "pars"),
                  fixed = lapply(hs, `[[`, "fixed"), shared = FALSE)
    fg <- .fwdMany(spec$g, bg, env, cores)
    obs <- unname(unclass(fg$values))
  }
  ee <- .shootErrors(spec, obs, fixed, FALSE, cores)
  err <- if (is.null(ee)) NULL else ee$err
  meta <- .shootMeta(spec, cache, obs, err)
  kr <- normL2_kernel(
    prediction = setNames(obs, spec$keys), err_list_opt = err, meta_list = meta,
    par_names_global = character(0), deriv2_requested = FALSE,
    threads = as.integer(cores), bloq_mode = spec$opt.BLOQ,
    build_hessian = FALSE, want_seed = TRUE, local_blocks = TRUE)

  # data seeds on the observations, and through the error model
  w_obs <- .normL2SeedCt(meta, obs, err, kr, spec$errmodel,
                         err_idx = ee$idx, err_split = ee$split)

  # the observation function backwards, onto the states and inner parameters
  u_g <- if (is.null(fg)) w_obs else .bwdNode(fg$tape, w_obs, env, cores)

  ends <- .shootEnds(spec, preds, FALSE)
  Wx <- lapply(seq_len(n), function(k) {
    pr <- preds[[k]]
    W <- array(0, c(nrow(pr), ncol(pr), S), dimnames = list(NULL, colnames(pr), NULL))
    u <- u_g[[k]]$out
    if (!is.null(u)) {
      hit <- intersect(dimnames(u)[[2L]], colnames(pr))
      W[, hit, 1L] <- u[, hit, 1L]
    }
    if (!segs[[k]]$last) {
      xe <- pr[nrow(pr), states]
      dphi <- .shootToChart(xe, spec$charts)$deriv
      for (i in seq_len(nx)) W[nrow(pr), states[i], 1L + i] <- dphi[i]
    }
    .ct(out = W)
  })
  u_x <- .bwdNode(fx$tape, Wx, env, cores, seeds = TRUE)

  segments <- lapply(seq_len(n), function(k) {
    seg <- segs[[k]]
    U <- u_x[[k]]$pars
    up <- u_g[[k]]$pars
    if (!is.null(up)) {
      rows <- union(rownames(U), rownames(up))
      U <- .pickCotangent(U, rows, S)
      U[rownames(up), 1L] <- U[rownames(up), 1L] + up[, 1L]
    }
    D <- attr(pin[[seg$cond]], "deriv")
    inner <- rownames(D)
    G1 <- NULL
    if (seg$j == 1L && !seg$last)
      G1 <- .shootChartGrowth(t(.pickCotangent(U, states, S)[, -1L, drop = FALSE]),
                              pin[[seg$cond]], spec)
    if (seg$j > 1L) inner <- setdiff(inner, states)
    Ui <- .pickCotangent(U, inner, S)
    G <- crossprod(D[inner, , drop = FALSE], Ui)          # theta x S
    if (seg$j > 1L) {
      s <- nodes[[seg$cond]][seg$j - 1L, states]
      Us <- .pickCotangent(U, states, S) * .shootFromChart(s, spec$charts)$deriv
      rownames(Us) <- seg$nodeNames
      G <- rbind(G, Us)
    }
    vars <- rownames(G)
    list(vars = vars,
         grad = setNames(G[, 1L], vars),
         hess = NULL,
         end = if (seg$last) NULL else ends[[k]]$value,
         jac = if (seg$last) NULL else
           matrix(t(G[, -1L, drop = FALSE]), nx, length(vars),
                  dimnames = list(states, vars)),
         growth = if (seg$last) NULL
                  else if (seg$j > 1L) t(G[seg$nodeNames, -1L, drop = FALSE])
                  else G1,
         chi2 = kr$local[[k]]$chi2, ndata = nrow(seg$data))
  })

  out <- list(value = kr$value, chi2 = kr$chi2,
              gaps = .shootGaps(spec, ends, nodes),
              segments = segments, sweep = "reverse")
  attr(out, attr.name) <- out$value
  env$prediction <- setNames(obs, spec$keys)
  env$states <- setNames(preds, spec$keys)
  out$env <- env
  out
}


# An observed state's trajectory as a function of time: the replicate means
# interpolated linearly, or with `smooth` a GCV smoothing spline (Horbelt,
# Timmer and Voss 2002).
.shootObservedCurve <- function(di, smooth = FALSE) {
  ag <- stats::aggregate(di$value, list(time = di$time), mean)
  if (nrow(ag) == 1L) return(function(t) rep(ag$x, length(t)))
  if (smooth && nrow(ag) >= 4L) {
    sp <- stats::smooth.spline(ag$time, ag$x)
    return(function(t) stats::predict(sp, t)$y)
  }
  function(t) stats::approx(ag$time, ag$x, xout = t, rule = 2)$y
}

# Start values of the nodes in chart coordinates, one matrix per condition:
# directly observed states from the data, others from a simulation at `pars`,
# or where that fails the last state reached or p's initial value.
.shootNodes <- function(objfun, pars, fixed = NULL, observed = TRUE,
                        cores = getOption("dMod.cores", 1L), smooth = FALSE) {
  spec <- attr(objfun, "spec")
  states <- spec$states
  pin <- .shootInner(spec, pars, fixed, cores)
  obsState <- .shootObservedStates(spec)
  out <- list()
  for (cn in spec$conditions) {
    tau <- spec$tau[[cn]]
    segsC <- Filter(function(s) s$cond == cn, spec$segs)
    grid <- sort(unique(unlist(lapply(segsC, `[[`, "times"))))
    p0 <- pin[[cn]]
    attr(p0, "deriv") <- NULL
    init <- unclass(p0)[states]
    X <- matrix(rep(init, each = length(tau) - 1L), length(tau) - 1L,
                dimnames = list(tau[-1L], states))
    # A node is the state before the events at its time, which its segment
    # applies. Solver and data report it after, so a node on an event reads the
    # state just before and observed states subtract the simulated jump.
    evt <- .shootEventTimes(spec, p0)
    onEvent <- tau[-1L] %in% evt
    dt <- 1e-8 * max(1, abs(tau))
    before <- tau[-1L][onEvent] - dt
    jump <- matrix(0, length(tau) - 1L, length(states), dimnames = list(NULL, states))
    sim <- try(suppressWarnings(.predictMany(
      spec$x, list(sort(unique(c(grid, before)))),
      list(as.parvec(unclass(p0)[names(p0)], deriv = FALSE)),
      cn, deriv = FALSE, cores = cores)[[1L]]), silent = TRUE)
    # only a simulation that ran to the end and stayed finite says anything
    # about the unobserved states; otherwise they keep their initial values
    if (!inherits(sim, "try-error") && NROW(sim) &&
        max(sim[, "time"]) >= max(tau) && all(is.finite(sim[, states]))) {
      st <- sim[, "time"]
      for (i in seq_len(length(tau) - 1L)) {
        X[i, ] <- sim[match(tau[i + 1L], st), states]
        if (onEvent[i]) {
          pre <- sim[match(tau[i + 1L] - dt, st), states]
          jump[i, ] <- X[i, ] - pre
          X[i, ] <- pre
        }
      }
    }
    if (observed && length(obsState)) {
      dc <- Reduce(rbind, lapply(segsC, `[[`, "data"))
      for (st in names(obsState)) {
        di <- dc[as.character(dc$name) == obsState[[st]], , drop = FALSE]
        if (nrow(di) < 1L) next
        X[, st] <- .shootObservedCurve(di, smooth)(tau[-1L]) - jump[, st]
      }
    }
    lg <- spec$charts == "log10"
    for (st in states[lg]) {
      bad <- !(X[, st] > 0)
      if (any(bad)) X[bad, st] <- min(c(X[!bad, st], 1e-8))
      X[, st] <- log10(X[, st])
    }
    out[[cn]] <- X
  }
  out
}

# Sets the unobserved states of every node by replacement synchronisation: one
# run from p's initial state with observed states reset to the data at every
# data time. A failed solve restarts from the next node as it was.
.shootSweep <- function(sobj, theta, nodes, fixed, cores) {
  spec <- attr(sobj, "spec")
  smooth <- identical(attr(sobj, "init"), "spline")
  unobs <- setdiff(spec$states, names(.shootObservedStates(spec)))
  if (!length(unobs)) return(nodes)
  pin <- .shootInner(spec, theta, fixed, cores)
  for (cn in spec$conditions) {
    tau <- spec$tau[[cn]]
    if (length(tau) < 2L) next
    p0 <- pin[[cn]]
    attr(p0, "deriv") <- NULL
    x <- unclass(p0)[spec$states]
    tr <- .shootSyncTrack(spec, p0, cn, x, tau[1L], tau[-1L], cores,
                          restart = nodes[[cn]], smooth = smooth)
    for (i in seq_along(tr)) {
      if (is.null(tr[[i]])) next
      ch <- try(.shootToChart(tr[[i]], spec$charts), silent = TRUE)
      if (!inherits(ch, "try-error")) nodes[[cn]][i, unobs] <- ch$value[unobs]
    }
  }
  nodes
}

# The data of the observed states of condition `cn` as reset values: one row
# per data time, the mean over replicates, without event times, where the data
# report the state after the event and a segment starts from the one before.
.shootResets <- function(spec, cn, events, smooth = FALSE) {
  obs <- .shootObservedStates(spec)
  segsC <- Filter(function(s) s$cond == cn, spec$segs)
  dc <- Reduce(rbind, lapply(segsC, `[[`, "data"))
  if (!length(obs) || is.null(dc) || !NROW(dc)) return(NULL)
  tt <- sort(unique(dc$time))
  tt <- tt[!(tt %in% events)]
  V <- matrix(NA_real_, length(tt), length(obs), dimnames = list(NULL, names(obs)))
  for (st in names(obs)) {
    di <- dc[as.character(dc$name) == obs[[st]], , drop = FALSE]
    if (!nrow(di)) next
    have <- tt %in% di$time
    V[have, st] <- .shootObservedCurve(di, smooth)(tt[have])
  }
  list(time = tt, value = V)
}

# States at the increasing times `at` of a run from `x` at `a`, observed states
# reset to the data on the way; each before any reset or event, NULL past a
# failure. With `restart` (node matrix), a failed run resumes at the next node.
.shootSyncTrack <- function(spec, p0, cn, x, a, at, cores, restart = NULL,
                            smooth = FALSE) {
  events <- .shootEventTimes(spec, p0)
  rs <- .shootResets(spec, cn, events, smooth)
  rt <- if (is.null(rs)) numeric(0) else rs$time[rs$time > a & rs$time <= max(at)]
  brk <- sort(unique(c(rt, at)))
  out <- vector("list", length(at))
  v <- unclass(p0)
  cur <- a
  ok <- TRUE
  for (b in brk) {
    if (ok) {
      v[spec$states] <- x
      pr <- try(suppressWarnings(.predictMany(spec$x, list(c(cur, b)),
                                              list(as.parvec(v, deriv = FALSE)), cn,
                                              deriv = FALSE, cores = cores)[[1L]]),
                silent = TRUE)
      ok <- !inherits(pr, "try-error") && NROW(pr) >= 2L &&
        pr[NROW(pr), "time"] == b && all(is.finite(pr[NROW(pr), spec$states]))
      if (ok) x <- pr[NROW(pr), spec$states]
    }
    i <- match(b, at)
    if (!is.na(i)) {
      if (ok) out[[i]] <- x
      else if (!is.null(restart)) {
        x <- .shootFromChart(restart[i, spec$states], spec$charts)$value
        ok <- TRUE
      }
    }
    if (ok && !is.null(rs)) {
      j <- match(b, rs$time)
      if (!is.na(j)) {
        r <- rs$value[j, ]
        r <- r[is.finite(r)]
        x[names(r)] <- r
      }
    }
    cur <- b
  }
  out
}

# The states on the whole time axis of every condition, simulated at `theta`
# in one piece; NULL for a condition whose solve does not reach its end.
.shootSimulate <- function(sobj, theta, fixed, cores) {
  spec <- attr(sobj, "spec")
  pin <- .shootInner(spec, theta, fixed, cores)
  out <- list()
  for (cn in spec$conditions) {
    grid <- sort(unique(unlist(lapply(Filter(function(s) s$cond == cn, spec$segs),
                                      `[[`, "times"))))
    p0 <- pin[[cn]]
    attr(p0, "deriv") <- NULL
    sim <- try(suppressWarnings(.predictMany(
      spec$x, list(grid), list(as.parvec(unclass(p0)[names(p0)], deriv = FALSE)),
      cn, deriv = FALSE, cores = cores)[[1L]]), silent = TRUE)
    if (!inherits(sim, "try-error") && NROW(sim) &&
        max(sim[, "time"]) >= max(grid) && all(is.finite(sim[, spec$states])))
      out[[cn]] <- sim
  }
  out
}

# Times of the fixed-time events of the prediction `x`, evaluated at the inner
# parameters `p0` of a condition. A time given as a parameter or an expression
# of parameters is evaluated; root events have no time and are left out.
.shootEventTimes <- function(spec, p0) {
  maps <- attr(spec$x, "mappings")
  if (!length(maps)) return(numeric(0))
  env <- environment(maps[[1L]])
  ev <- NULL
  func <- if (!is.null(env)) get0("func", envir = env, inherits = FALSE)
  if (!is.null(func)) ev <- attr(func, "events")
  if (is.null(ev) && !is.null(env)) {
    ctl <- get0("controls", envir = env, inherits = FALSE)
    if (is.list(ctl)) ev <- ctl$events
  }
  if (is.null(ev) || !NROW(ev)) return(numeric(0))
  ev <- as.data.frame(ev)
  fixedTime <- if (is.null(ev$root)) rep(TRUE, nrow(ev)) else is.na(ev$root)
  vals <- as.list(unclass(p0))
  out <- vapply(as.character(ev$time[fixedTime]), function(tx) {
    v <- suppressWarnings(as.numeric(tx))
    if (!is.na(v)) return(v)
    r <- try(eval(parse(text = tx), vals), silent = TRUE)
    if (inherits(r, "try-error") || !is.numeric(r)) NA_real_ else as.numeric(r)[1L]
  }, 0)
  unique(out[is.finite(out)])
}

# Observables whose equation is exactly a state: state -> observable name.
.shootObservedStates <- function(spec) {
  states <- spec$states
  if (is.null(spec$g)) return(setNames(as.list(states), states))
  maps <- attr(spec$g, "mappings")
  eq <- if (length(maps)) attr(maps[[1L]], "equations") else NULL
  if (is.null(eq)) return(list())
  eq <- unclass(eq)
  hit <- which(trimws(eq) %in% states)
  setNames(as.list(names(eq)[hit]), trimws(eq[hit]))
}


## ---- From a normL2 objective ------------------------------------------------

# The pieces of normL2(multipleShootingControl = ): the chain split into
# observation, prediction and transformation, and the segments laid out.
.shootFromNormL2 <- function(data, prd, errmodel, times, attr.name,
                             opt.BLOQ, shooting) {
  if (isTRUE(shooting)) shooting <- list()
  if (!is.list(shooting))
    stop("normL2: multipleShootingControl is TRUE or a list, e.g. ",
         "list(charts = c(n = \"log10\")).", call. = FALSE)
  known <- c("nodes", "charts", "scale", "init", "growth", "misfit", "minPoints",
             "breaks")
  bad <- setdiff(names(shooting), known)
  if (length(bad) || (length(shooting) && is.null(names(shooting))))
    stop("normL2: multipleShootingControl takes ", paste(known, collapse = ", "),
         if (length(bad)) paste0(", not ", paste(bad, collapse = ", ")), ".",
         call. = FALSE)
  nodes <- if (is.null(shooting$nodes)) "auto" else shooting$nodes
  transitions <- identical(nodes, "transitions")
  adaptive <- identical(nodes, "auto")
  growth <- if (is.null(shooting$growth)) 10 else shooting$growth
  if (!is.numeric(growth) || length(growth) != 1L || !(growth > 1))
    stop("normL2: multipleShootingControl$growth is a number above 1.", call. = FALSE)
  misfit <- if (is.null(shooting$misfit)) 3 else shooting$misfit
  if (!is.numeric(misfit) || length(misfit) != 1L || !(misfit > 0))
    stop("normL2: multipleShootingControl$misfit is a positive number.", call. = FALSE)
  init <- if (is.null(shooting$init)) "data" else shooting$init
  if (!is.list(init) && !(is.character(init) && length(init) == 1L &&
                          init %in% c("data", "spline", "simulation")))
    stop("normL2: multipleShootingControl$init is \"data\", \"spline\", \"simulation\" ",
         "or a list of node values.", call. = FALSE)
  if ((adaptive || transitions) && is.list(init))
    stop("normL2: node values in multipleShootingControl$init need the node ",
         "times in multipleShootingControl$nodes.", call. = FALSE)
  parts <- .shootSplitChain(prd)
  minPoints <- if (is.null(shooting$minPoints)) length(.shootStates(parts$x)) + 1L
               else shooting$minPoints
  if (!is.numeric(minPoints) || length(minPoints) != 1L || !(minPoints >= 1) ||
      minPoints != round(minPoints))
    stop("normL2: multipleShootingControl$minPoints is a whole number of at least 1.",
         call. = FALSE)
  minPoints <- as.integer(minPoints)
  if (transitions) {
    nodes <- .shootTransitionNodes(data, parts$x, times)
    # without a fast change in the data the layout is "auto" after all
    if (is.null(nodes)) adaptive <- TRUE
  }
  if (adaptive) nodes <- .shootAutoNodes(data, parts$x, times, oscillations = TRUE,
                                         minPoints = minPoints)
  # a break is a node whose continuity is not enforced; "auto" lays one there
  breaks <- NULL
  if (!is.null(shooting$breaks)) {
    x.cond <- names(attr(parts$x, "mappings"))
    conds <- if (is.null(x.cond) || all(!nzchar(x.cond))) names(data)
             else intersect(names(data), x.cond)
    breaks <- .shootPerCondition(shooting$breaks, conds, "breaks")
    nodesL <- .shootPerCondition(nodes, conds, "nodes")
    for (cn in conds) {
      b <- as.numeric(breaks[[cn]])
      if (adaptive) nodesL[[cn]] <- sort(unique(c(nodesL[[cn]], b)))
      else if (!all(b %in% nodesL[[cn]]))
        stop("normL2: multipleShootingControl$breaks must be node times.", call. = FALSE)
      breaks[[cn]] <- b
    }
    nodes <- nodesL
  }
  sobj <- .shootObjective(data, parts$x, parts$p, parts$g, errmodel,
                          nodes = nodes, charts = shooting$charts,
                          times = times, attr.name = attr.name,
                          opt.BLOQ = opt.BLOQ)
  if (is.list(init)) .shootCheckNodes(init, attr(sobj, "spec"))
  attr(sobj, "scale") <- shooting$scale
  attr(sobj, "init") <- init
  attr(sobj, "adaptive") <- adaptive
  attr(sobj, "growth") <- growth
  attr(sobj, "misfit") <- misfit
  attr(sobj, "minLength") <- .shootMinLength(data, attr(sobj, "conditions"))
  attr(sobj, "minPoints") <- minPoints
  attr(sobj, "breaks") <- breaks
  sobj
}

## ---- Adaptive nodes -----------------------------------------------------------
## The driver splits a segment on growth beyond `growth`, a failed solve, or a
## local misfit; each half keeps `minPoints` data points.

.shootAutoNodes <- function(data, x, times, n = 10L, oscillations = FALSE,
                            minPoints = 1L) {
  x.cond <- names(attr(x, "mappings"))
  conds <- if (is.null(x.cond) || all(!nzchar(x.cond))) names(data)
           else intersect(names(data), x.cond)
  span <- .shootSpan(data, times, conds)
  setNames(lapply(conds, function(cn) {
    ts <- span[[cn]][1L]; te <- span[[cn]][2L]
    m <- n
    if (oscillations)
      m <- min(max(n, .shootTurns(data[[cn]]) + 1L),
               max(n, nrow(data[[cn]]) %/% minPoints))
    ts + (te - ts) * seq_len(m - 1L) / m
  }), conds)
}

# Turning points of a condition's data, the most of any observable: extrema of
# the replicate means the data reverse from by more than `k` times their noise.
.shootTurns <- function(d, k = 5) {
  best <- 0L
  for (nm in unique(as.character(d$name))) {
    di <- d[as.character(d$name) == nm, , drop = FALSE]
    if (nrow(di) < 3L) next
    v <- stats::aggregate(di$value, list(time = di$time), mean)$x
    s <- if (!is.null(di$sigma)) stats::median(di$sigma, na.rm = TRUE) else NA_real_
    if (!is.finite(s) || !(s > 0)) s <- stats::mad(diff(v)) / sqrt(2)
    h <- k * s
    if (!is.finite(h) || !(h > 0)) next
    n <- 0L; dir <- 0L; ext <- v[1L]
    for (x in v[-1L]) {
      if (dir == 0L) {
        if (abs(x - ext) > h) { dir <- if (x > ext) 1L else -1L; ext <- x }
      } else if ((x - ext) * dir > 0) {
        ext <- x
      } else if (abs(x - ext) > h) {
        n <- n + 1L; dir <- -dir; ext <- x
      }
    }
    best <- max(best, n)
  }
  best
}

# Nodes just before fast changes of the data, where an observable's rate lies
# more than `k` robust deviations off its median, a twentieth of the median
# interval before the change. NULL when no condition shows one.
.shootTransitionNodes <- function(data, x, times, k = 8, n = 10L) {
  base <- .shootAutoNodes(data, x, times, n)
  span <- .shootSpan(data, times, names(base))
  found <- FALSE
  out <- setNames(lapply(names(base), function(cn) {
    d <- data[[cn]]
    ts <- span[[cn]][1L]
    onsets <- numeric(0); spacing <- Inf
    for (nm in unique(as.character(d$name))) {
      di <- d[as.character(d$name) == nm, , drop = FALSE]
      ag <- stats::aggregate(di$value, list(time = di$time), mean)
      if (nrow(ag) < 5L) next
      dt <- diff(ag$time)
      spacing <- min(spacing, dt)
      rate <- diff(ag$x) / dt
      # noise-free data between changes can leave no spread at all; then any
      # rate off the median is a change
      md <- stats::median(rate)
      sc <- stats::mad(rate)
      if (!(sc > 0)) sc <- 64 * .Machine$double.eps * max(1, abs(md))
      fast <- which(abs(rate - md) > k * sc)
      if (!length(fast)) next
      # consecutive fast samples are one change, which starts at its first
      first <- fast[c(TRUE, diff(fast) > 2L)]
      onsets <- c(onsets, ag$time[first])
    }
    onsets <- sort(unique(onsets))
    te <- max(d$time)
    if (!length(onsets)) return(base[[cn]])
    found <<- TRUE
    # Changes closer than a quarter of the median of the longer half of the
    # intervals are one, such as a spike's rise and fall.
    if (length(onsets) > 2L) {
      di <- diff(onsets)
      typical <- stats::median(di[di >= stats::median(di)])
      onsets <- onsets[c(TRUE, di > typical / 4)]
    }
    gap <- if (length(onsets) > 1L) stats::median(diff(onsets)) else te - ts
    lead <- max(2 * spacing, gap / 20)
    pre <- onsets - lead
    pre <- pre[pre > ts + spacing / 2 & pre < te]
    # the even nodes stay where no change is near
    far <- base[[cn]][vapply(base[[cn]], function(b) all(abs(b - onsets) > gap), TRUE)]
    nd <- sort(unique(c(pre, far)))
    # no two nodes closer than two data spacings
    keep <- c(TRUE, diff(nd) > 2 * spacing)
    nd[keep]
  }), names(base))
  if (found) out else NULL
}

# Half the finest spacing of a condition's data: no segment gets shorter.
.shootMinLength <- function(data, conds) {
  setNames(lapply(conds, function(cn) {
    d <- diff(sort(unique(data[[cn]]$time)))
    d <- d[d > 0]
    if (length(d)) min(d) / 2 else Inf
  }), conds)
}

# Where segment k is cut: at the grid time nearest its middle, or the middle
# itself when no grid time lies in the middle half. NA when the halves would
# fall below the minimal length.
.shootSplitTime <- function(seg, minLength, events = numeric(0), minPoints = 1L) {
  len <- seg$end - seg$start
  if (!(len / 2 >= minLength)) return(NA_real_)
  mid <- seg$start + len / 2
  # a node on an event time starts from the state before it, which neither
  # the data nor the trajectory at that time report; keep off them
  inner <- seg$times[seg$times > seg$start + len / 4 & seg$times < seg$end - len / 4]
  inner <- inner[!(inner %in% events)]
  cand <- if (length(inner)) inner[order(abs(inner - mid))]
          else if (mid %in% events) mid + len / 8 else mid
  # each half keeps at least minPoints data points of its own: a node is
  # there to be told by data, not to take the place of the parameters
  if (minPoints > 1L) {
    td <- seg$data$time
    ok <- vapply(cand, function(tn) sum(td >= seg$start & td < tn) >= minPoints &&
                                    sum(td >= tn & td <= seg$end) >= minPoints, TRUE)
    cand <- cand[ok]
  }
  if (length(cand)) cand[1L] else NA_real_
}

# The chart value of every state at the times `tn` inside the segments `ks`,
# from each segment's own trajectory at (theta, nodes), in one batched solve;
# NULL for a segment whose solve fails.
.shootValuesAt <- function(sobj, theta, nodes, fixed, ks, tn, cores) {
  spec <- attr(sobj, "spec")
  if (!length(ks)) return(list())
  pin <- .shootInner(spec, theta, fixed, cores)
  pk <- .shootSegmentParsAll(pin, spec, nodes, FALSE, which = ks)
  starts <- vapply(spec$segs[ks], `[[`, 0, "start")
  conds <- vapply(spec$segs[ks], `[[`, "", "cond")
  one <- function(pr, t) {
    if (inherits(pr, "try-error") || NROW(pr) < 2L || pr[NROW(pr), "time"] != t)
      return(NULL)
    xv <- pr[NROW(pr), spec$states]
    if (!all(is.finite(xv))) return(NULL)
    ch <- try(.shootToChart(xv, spec$charts), silent = TRUE)
    if (inherits(ch, "try-error")) NULL else ch$value
  }
  grids <- Map(c, starts, tn)
  prs <- try(suppressWarnings(unname(unclass(.predictMany(
    spec$x, grids, pk, conds, deriv = FALSE, cores = cores)))), silent = TRUE)
  # a failure anywhere fails the batch; then each segment on its own
  if (inherits(prs, "try-error"))
    prs <- lapply(seq_along(ks), function(i) try(suppressWarnings(.predictMany(
      spec$x, grids[i], pk[i], conds[i], deriv = FALSE, cores = cores)[[1L]]), silent = TRUE))
  Map(one, prs, tn)
}

# A new layout: the objective rebuilt on the node times `tau` (per condition,
# the start included) and the node values kept, the new ones from `fill`, a
# list per condition of chart values named by time.
.shootRelayout <- function(sobj, tau, nodes, fill) {
  a <- attr(sobj, "args")
  new <- .shootObjective(a$data, a$x, a$p, a$g, a$errmodel,
                         nodes = lapply(tau, function(t) t[-1L]),
                         charts = a$charts, times = a$times,
                         attr.name = a$attr.name, opt.BLOQ = a$opt.BLOQ)
  for (nm in c("scale", "init", "adaptive", "growth", "misfit", "minLength", "minPoints",
              "breaks"))
    attr(new, nm) <- attr(sobj, nm)
  spec <- attr(new, "spec")
  nodesN <- list()
  for (cn in spec$conditions) {
    tn <- spec$tau[[cn]][-1L]
    old <- nodes[[cn]]
    M <- matrix(NA_real_, length(tn), length(spec$states),
                dimnames = list(tn, spec$states))
    hit <- match(as.character(tn), rownames(old))
    M[!is.na(hit), ] <- old[hit[!is.na(hit)], spec$states]
    for (i in which(is.na(hit))) M[i, ] <- fill[[cn]][[as.character(tn[i])]][spec$states]
    nodesN[[cn]] <- M
  }
  list(sobj = new, nodes = nodesN)
}

# A composition flattened into its factors, left to right.
.shootFactors <- function(f) {
  st <- .fnNode(f)
  if (!is.null(st) && identical(st$op, "*"))
    return(c(.shootFactors(st$p1), .shootFactors(st$p2)))
  list(f)
}

# g * x * p: observation functions, exactly one prediction function, and the
# parameter transformations it is fed from.
.shootSplitChain <- function(prd) {
  fs <- .shootFactors(prd)
  kinds <- vapply(fs, function(f) { k <- .fnKind(f); if (is.null(k)) "" else k }, "")
  ix <- which(kinds == "prdfn")
  n <- length(fs)
  if (length(ix) != 1L || ix == n || !all(kinds[(ix + 1L):n] == "parfn") ||
      (ix > 1L && !all(kinds[seq_len(ix - 1L)] == "obsfn")))
    stop("normL2: multiple shooting needs the prediction as a chain g * x * p ",
         "of observation functions, one prediction function from Xs() and ",
         "parameter transformations.", call. = FALSE)
  list(g = if (ix > 1L) Reduce(`*`, fs[seq_len(ix - 1L)]) else NULL,
       x = fs[[ix]],
       p = Reduce(`*`, fs[(ix + 1L):n]))
}

# The multipleShootingControl of an objective, NULL for one without.
.shootControlOf <- function(f) {
  if (!is.function(f) || is.primitive(f)) return(NULL)
  ctl <- environment(f)$controls
  if (is.list(ctl)) ctl$multipleShootingControl else NULL
}

# The summands multiple shooting applies to.
.shootingTerms <- function(objfun)
  Filter(function(t) !is.null(.shootControlOf(t)), .objTerms(objfun))

# The multiple-shooting objective of a normL2, alone or as the one summand of
# a sum that holds the control, laid out from the control as it is now.
.shootObjOf <- function(objfun) {
  tm <- .shootingTerms(objfun)
  if (length(tm) != 1L)
    stop("trust: an objective may hold one multiple-shooting term, this one ",
         "holds ", length(tm), ".", call. = FALSE)
  e <- environment(tm[[1L]])
  .shootFromNormL2(e$data, e$x, e$errmodel, e$times, e$attrName,
                   e$optBLOQ, e$controls$multipleShootingControl)
}

# The same objective, optimised by single shooting: a wrapper that answers as
# it does and keeps none of its terms or controls. The objective itself, and
# whoever else holds it, keeps its control.
.singleShooting <- function(f) {
  if (!length(.shootingTerms(f))) return(f)
  g <- function(...) f(...)
  at <- attributes(f)
  attributes(g) <- at[setdiff(names(at), c("terms", "srcref"))]
  g
}
