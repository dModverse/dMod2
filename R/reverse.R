## Reverse-mode helpers -----------------------------------------------------
## Copyright (C) 2026 Simon Beyer

# A cotangent: `out` on the node's output matrix, `pars` on its parameters,
# both optional and additive. The trailing axis of extent K holds the cotangent
# in slice 1 and its directional derivatives after it, or K seeds in seed mode.
.ct <- function(out = NULL, pars = NULL) {
  if (!is.null(pars) && is.null(dim(pars)))
    pars <- matrix(pars, ncol = 1L, dimnames = list(names(pars), NULL))
  if (!is.null(out) && length(dim(out)) == 2L) {
    dn <- dimnames(out)
    dim(out) <- c(dim(out), 1L)
    dimnames(out) <- c(dn, list(NULL))
  }
  list(out = out, pars = pars)
}

# How many directions a half has. Works on both shapes, [p, K] and
# [n, m, K], because the axis is the last one either way.
.ctK <- function(x) {
  if (is.null(x)) return(1L)
  d <- dim(x)
  as.integer(d[length(d)])
}

.ctZero <- function(nms, K = 1L)
  matrix(0, length(nms), K, dimnames = list(nms, NULL))

# Inside .withSeeds(S, ...) the trailing axis of a second-order cotangent holds S blocks of
# equal width, each a seed followed by its derivatives along the tangents, so that one backward
# sweep answers S seeds; outside it there is one block.
.reverseState <- new.env(parent = emptyenv())
.withSeeds <- function(S, expr) {
  old <- .reverseState$seeds
  .reverseState$seeds <- S
  on.exit(.reverseState$seeds <- old)
  expr
}

# The blocks of a trailing axis of width K: `first`, the position of every seed, and `rest`,
# [D, S], the positions of its derivatives along the D tangents.
.ctBlocks <- function(K, S = .reverseState$seeds %||% 1L) {
  if (K %% S)
    stop("a cotangent of width ", K, " does not split into ", S, " seed blocks.", call. = FALSE)
  B <- K %/% S
  first <- (seq_len(S) - 1L) * B + 1L
  list(S = S, D = B - 1L, first = first,
       rest = matrix(rep(first, each = B - 1L) + seq_len(B - 1L), B - 1L, S))
}

# The solver's answer as a K-column cotangent: per seed the gradient and its directional
# derivatives along the chain's tangents beside it.
.adjointCt <- function(res, K = 1L) {
  b <- .ctBlocks(K)
  g <- res$cotangent
  if (is.null(dim(g))) g <- matrix(g, ncol = 1L, dimnames = list(names(g), NULL))
  u <- matrix(0, nrow(g), K, dimnames = list(rownames(g), NULL))
  u[, b$first] <- g[, seq_len(b$S)]
  if (b$D == 0L) return(u)
  if (is.null(res$curvature))
    stop("a second-order cotangent needs a model built with ",
         "derivMode = \"forward-reverse\"", call. = FALSE)
  for (s in seq_len(b$S))
    u[, b$rest[, s]] <- matrix(res$curvature[, seq_len(b$D), s], nrow = nrow(u))
  u
}

.ct_null <- function(w) is.null(w) || (is.null(w$out) && is.null(w$pars))

# The sum of row-named matrices of one width on the rows `nms`, in one pass; rows outside `nms`
# are dropped, as .pickCotangent() does.
.sumCotangents <- function(cts, nms) {
  cts <- Filter(Negate(is.null), cts)
  if (!length(cts)) return(.ctZero(nms, 1L))
  K <- ncol(cts[[1L]])
  out <- matrix(0, length(nms), K, dimnames = list(nms, NULL))
  i <- match(unlist(lapply(cts, rownames), use.names = FALSE), nms)
  v <- do.call(rbind, cts)
  hit <- !is.na(i)
  if (any(hit)) out[sort(unique(i[hit])), ] <- rowsum(v[hit, , drop = FALSE], i[hit], reorder = TRUE)
  out
}

# Adds two row-named matrices on the union of their rows, keeping the order of
# the first. Positions are resolved once with match(), as this runs per node.
.addNamed <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  if (ncol(a) != ncol(b))
    stop("two cotangent halves have ", ncol(a), " and ", ncol(b),
         " directions; a node handed on a width its neighbour does not have.",
         call. = FALSE)
  na <- rownames(a); nb <- rownames(b)
  nms <- c(na, nb[is.na(match(nb, na))])
  out <- matrix(0, length(nms), ncol(a), dimnames = list(nms, NULL))
  out[seq_along(na), ] <- a
  ib <- match(nb, nms)
  out[ib, ] <- out[ib, , drop = FALSE] + b
  out
}

.addCt <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  o <- if (is.null(a$out)) b$out else if (is.null(b$out)) a$out else a$out + b$out
  .ct(o, .addNamed(a$pars, b$pars))
}

# Which object a cotangent of width K needs from the prediction. K > 1 is
# second order and needs the forward-over-reverse object.
.requireReverse <- function(has_reverse, has_reverse2, K) {
  if (K > 1L && !has_reverse2)
    stop("the reverse mode propagates directions here, which needs the ",
         "forward-over-reverse object; rebuild via odemodel(..., derivMode = ",
         "c(\"forward\", \"forward-reverse\")).", call. = FALSE)
  if (K == 1L && !has_reverse)
    stop("the model has no reverse object; rebuild via odemodel(..., ",
         "derivMode = c(\"forward\", \"reverse\")).", call. = FALSE)
  invisible(NULL)
}

# A seeded solve without a cotangent has failed; a silent zero would drop
# its block of the gradient.
.requireAdjoint <- function(res, condition = NULL) {
  if (!is.null(res$cotangent)) return(invisible(NULL))
  stop("the backward solve returned no cotangent",
       if (is.null(condition)) "" else paste0(" for condition ", condition),
       ", though it reported success. The model was seeded, so this is a ",
       "backend fault rather than a modelling one.", call. = FALSE)
}

# A matrix on exactly the rows `nms`, zero where `w` has none. `K` is the width
# expected back: an absent half is zero at the width of its neighbour.
.pickCotangent <- function(w, nms, K = .ctK(w)) {
  out <- .ctZero(nms, K)
  if (is.null(w) || !length(nms)) return(out)
  i <- match(nms, rownames(w))
  hit <- !is.na(i)
  if (any(hit)) out[hit, ] <- w[i[hit], , drop = FALSE]
  out
}

# The rows of `w` among `nms`, without the zero rows of the others; at another width as
# .pickCotangent().
.pickRows <- function(w, nms, K = .ctK(w)) {
  if (is.null(w) || ncol(w) != K) return(.pickCotangent(w, nms, K))
  w[!is.na(match(rownames(w), nms)), , drop = FALSE]
}

# A cotangent on a subset of columns, widened to the full state set with zeros.
# Slice 1 becomes the solver's `cotangent`, the further direction slices its
# `curvature` (NULL at first order), since a seed column costs a whole sweep.
.widenCotangent <- function(w, full, subset) {
  n <- dim(w)[1L]
  b <- .ctBlocks(.ctK(w))
  W <- array(0, c(n, length(full), b$S),
             dimnames = list(NULL, full, NULL))
  hit <- intersect(subset, full)
  if (length(hit)) W[, hit, ] <- w[, hit, b$first, drop = FALSE]
  if (b$D == 0L) return(list(cotangent = W, curvature = NULL))
  # The curvature is positional, in the column order of `full`.
  cv <- array(0, c(n, length(full), b$S, b$D))
  if (length(hit))
    for (s in seq_len(b$S)) cv[, match(hit, full), s, ] <- w[, hit, b$rest[, s], drop = FALSE]
  list(cotangent = W, curvature = cv)
}

# Seed-mode counterpart: every slice of the trailing axis goes to the solver
# as a seed column of its own, none becomes a curvature.
.widenSeeds <- function(w, full, subset) {
  n <- dim(w)[1L]
  S <- .ctK(w)
  W <- array(0, c(n, length(full), S), dimnames = list(NULL, full, NULL))
  hit <- intersect(subset, full)
  if (length(hit)) W[, hit, ] <- w[, hit, , drop = FALSE]
  W
}

# The solver's answer to a multi-seed sweep, one column per seed.
.adjointSeeds <- function(res) {
  u <- res$cotangent
  if (is.null(dim(u)))
    u <- matrix(u, ncol = 1L, dimnames = list(names(u), NULL))
  dimnames(u) <- list(rownames(u), NULL)
  u
}

# A vjp called from outside the chain may get a plain matrix or named vector;
# these add the direction axis.
.asCtOut  <- function(w) if (is.null(w) || length(dim(w)) == 3L) w else .ct(out = w)$out
.asCtPars <- function(w) if (is.null(w) || !is.null(dim(w))) w else .ct(pars = w)$pars

# A cotangent widened to K directions with zeros: an objective's seed is a
# constant, so its tangents vanish.
.ctWiden <- function(w, K) {
  if (K <= 1L) return(w)
  d <- dim(w)
  out <- array(0, c(d[1L], d[2L], K), dimnames = dimnames(w))
  out[, , 1L] <- w[, , 1L]
  out
}

# The "time" column is a label, not an output, and never takes a cotangent.
.dropTime <- function(m) {
  if (is.null(m) || length(dim(m)) < 2L) return(m)
  keep <- dimnames(m)[[2L]] != "time"
  if (all(keep)) return(m)
  if (length(dim(m)) == 2L) m[, keep, drop = FALSE] else m[, keep, , drop = FALSE]
}

# Vjp of a parfn that builds its own Jacobian, as Pimpl does by the implicit
# function theorem: the transposed Jacobian applied to the cotangent. The
# Jacobian's width is the node's own parameter set, not n_theta.
.parfnVjpFromJacobian <- function(p2p) {
  function(pars, fixed = NULL, cotangent, condition = NULL) {
    w <- .asCtPars(cotangent)
    K <- .ctK(w)
    # Incoming tangents: the node's curvature is contracted along them. `deriv`
    # is then stripped so the Jacobian is not chained to the outer parameters.
    V <- if (K > 1L) attr(pars, "deriv") else NULL
    # Directions with no tangent here meet no curvature: the further slices
    # are pulled back as the first, without second derivatives.
    curved <- K > 1L && !is.null(V) && any(V != 0)
    attr(pars, "deriv") <- NULL
    attr(pars, "deriv2") <- NULL
    v <- p2p(pars, fixed = fixed, deriv = TRUE, deriv2 = curved,
             condition = condition)
    J <- attr(v, "deriv")
    # A missing Jacobian is an error, not a zero cotangent.
    if (is.null(J) || !is.matrix(J))
      stop("a transformation returned no Jacobian, so the backward pass has ",
           "nothing to contract here. A preceding warning usually names the ",
           "cause; Pimpl falls back to value only when the implicit ",
           "function theorem cannot be applied at the current root.",
           call. = FALSE)
    wv <- .pickCotangent(w, rownames(J))
    u  <- crossprod(J, wv)
    rownames(u) <- colnames(J)
    if (curved) {
      H <- attr(v, "deriv2")
      if (is.null(H))
        stop("a second-order cotangent needs this transformation's own second ",
             "derivatives, and it returned none; rebuild it with deriv2 = TRUE.",
             call. = FALSE)
      # w' H, the node's curvature seen through each seed, times the
      # tangents that arrived: the second half of d/dv (J' w).
      b <- .ctBlocks(K)
      Vk <- matrix(0, ncol(J), b$D, dimnames = list(colnames(J), NULL))
      if (!is.null(V)) {
        hit <- intersect(rownames(V), colnames(J))
        if (length(hit))
          Vk[hit, ] <- V[hit, seq_len(b$D), drop = FALSE]
      }
      for (s in seq_len(b$S)) {
        M <- apply(H, c(2L, 3L), function(col) sum(col * wv[, b$first[s]]))
        u[, b$rest[, s]] <- u[, b$rest[, s], drop = FALSE] + M %*% Vk
      }
    }
    .pickCotangent(u, names(pars))
  }
}
