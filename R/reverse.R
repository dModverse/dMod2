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

# The solver's answer as a K-column cotangent: the gradient in slice 1, its
# directional derivatives along the chain's tangents beside it.
.adjointCt <- function(res, K = 1L) {
  g <- res$cotangent[, 1L]
  u <- matrix(g, ncol = 1L, dimnames = list(names(g), NULL))
  if (K == 1L) return(u)
  if (is.null(res$curvature))
    stop("a second-order cotangent needs a model built with ",
         "derivMode = \"forward-reverse\"", call. = FALSE)
  cbind(u, matrix(res$curvature[, seq_len(K - 1L), 1L], nrow = nrow(u)))
}

.ct_null <- function(w) is.null(w) || (is.null(w$out) && is.null(w$pars))

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

# A cotangent on a subset of columns, widened to the full state set with zeros.
# Slice 1 becomes the solver's `cotangent`, the further direction slices its
# `curvature` (NULL at first order), since a seed column costs a whole sweep.
.widenCotangent <- function(w, full, subset) {
  n <- dim(w)[1L]
  K <- .ctK(w)
  W <- array(0, c(n, length(full), 1L),
             dimnames = list(NULL, full, NULL))
  hit <- intersect(subset, full)
  if (length(hit)) W[, hit, 1L] <- w[, hit, 1L, drop = FALSE]
  if (K == 1L) return(list(cotangent = W, curvature = NULL))
  # The curvature is positional, in the column order of `full`.
  cv <- array(0, c(n, length(full), 1L, K - 1L))
  if (length(hit)) cv[, match(hit, full), 1L, ] <- w[, hit, -1L, drop = FALSE]
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

# One direction slice of an out-half, as the matrix a leaf's vjp takes.
.ctSlice <- function(w, k = 1L) {
  m <- w[, , k, drop = FALSE]
  dim(m) <- dim(w)[1:2]
  dimnames(m) <- dimnames(w)[1:2]
  m
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
      # w' H, the node's curvature seen through the cotangent, times the
      # tangents that arrived: the second half of d/dv (J' w).
      M <- apply(H, c(2L, 3L), function(col) sum(col * wv[, 1L]))
      Vk <- matrix(0, ncol(J), K - 1L, dimnames = list(colnames(J), NULL))
      if (!is.null(V)) {
        hit <- intersect(rownames(V), colnames(J))
        if (length(hit))
          Vk[hit, ] <- V[hit, seq_len(K - 1L), drop = FALSE]
      }
      u[, -1L] <- u[, -1L, drop = FALSE] + M %*% Vk
    }
    .pickCotangent(u, names(pars))
  }
}
