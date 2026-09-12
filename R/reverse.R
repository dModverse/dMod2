## Reverse-mode helpers -----------------------------------------------------
##
## The forward chain threads a Jacobian: every leaf receives dP/dtheta on its
## input and returns dX/dtheta on its output, so the chain rule happens inside
## the leaves and the width of everything in flight is n_theta. The reverse
## chain threads a cotangent the other way, and its width is the number of
## seeds, which for an objective gradient is one.
##
## A cotangent travelling through the tree is the same shape whichever node it
## passes, because every node's output is the same pair:
##
##   out    the matrix a prediction or an observation carries, NULL for a parfn
##   pars   the parameters the node passes through on its `parameters` attribute
##
## Both halves are optional and both accumulate. The `pars` half is what makes
## the tree a graph rather than a chain: an observation function reads the
## prediction's parameters as well as its values, so its cotangent lands on the
## same vector the prediction's own vjp writes to.
##
## Copyright (C) 2026 Simon Beyer

# The cotangent one node hands to the node below it. Both halves carry a
# trailing direction axis of extent K: slice 1 is the cotangent itself, slices
# 2..K are its derivatives along the directions the value pass carried. First
# order is K = 1, the same storage plus a dim attribute and the same code path.
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

# How many directions a half carries. Works on both shapes, [p, K] and
# [n, m, K], because the axis is the last one either way.
.ctK <- function(x) if (is.null(x)) 1L else as.integer(utils::tail(dim(x), 1L))

.ctZero <- function(nms, K = 1L)
  matrix(0, length(nms), K, dimnames = list(nms, NULL))

# The solver's answer as a K-column cotangent: the gradient in slice 1, its
# directional derivatives beside it. adjoint2's middle axis is the sensitivity
# set the call seeded, which is exactly the directions the chain carries.
.adjointCt <- function(res, K = 1L) {
  g <- res$adjoint[, 1L]
  u <- matrix(g, ncol = 1L, dimnames = list(names(g), NULL))
  if (K == 1L) return(u)
  if (is.null(res$adjoint2))
    stop("a second-order cotangent needs a model built with ",
         "derivMode = \"forward-reverse\"", call. = FALSE)
  cbind(u, matrix(res$adjoint2[, seq_len(K - 1L), 1L], nrow = nrow(u)))
}

.ct_null <- function(w) is.null(w) || (is.null(w$out) && is.null(w$pars))

# Adds two row-named matrices on the union of their rows. Cheaper than a merge
# and it keeps the order of the first, which is the order a caller expects back.
# Positions are resolved once with match(); indexing by name resolves them again
# on every use, and this runs once per node per condition.
.addNamed <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  if (ncol(a) != ncol(b))
    stop("two cotangent halves carry ", ncol(a), " and ", ncol(b),
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

# A named vector on exactly `nms`, zero where the cotangent says nothing.
# A solve that answers a seed and returns no adjoint has failed, and an absent
# cotangent reads as a zero one everywhere above. Saying so beats a gradient
# that is quietly zero in the directions the solver dropped.
# Which object a cotangent of width K asks the prediction for, and what to say
# when it is not there. K > 1 is second order and needs the fifth compilation.
.requireReverse <- function(has_reverse, has_reverse2, K) {
  if (K > 1L && !has_reverse2)
    stop("the reverse mode carries directions here, which needs the ",
         "forward-over-reverse object; rebuild via odemodel(..., derivMode = ",
         "c(\"forward\", \"forward-reverse\")).", call. = FALSE)
  if (K == 1L && !has_reverse)
    stop("the model has no reverse object; rebuild via odemodel(..., ",
         "derivMode = c(\"forward\", \"reverse\")).", call. = FALSE)
  invisible(NULL)
}

.requireAdjoint <- function(res, condition = NULL) {
  if (!is.null(res$adjoint)) return(invisible(NULL))
  stop("the backward solve returned no adjoint",
       if (is.null(condition)) "" else paste0(" for condition ", condition),
       ", though it reported success. The model was seeded, so this is a ",
       "backend fault rather than a modelling one.", call. = FALSE)
}

# One match() rather than an intersect() and two name lookups: this is the
# hottest thing in a reverse objective that is not the solver.
# `K` is the width the caller expects back, which is not always the width `w`
# has: an absent half is zero at whatever width its neighbour carries.
.pickCotangent <- function(w, nms, K = .ctK(w)) {
  out <- .ctZero(nms, K)
  if (is.null(w) || !length(nms)) return(out)
  i <- match(nms, rownames(w))
  hit <- !is.na(i)
  if (any(hit)) out[hit, ] <- w[i[hit], , drop = FALSE]
  out
}

# A seed on a subset of columns, widened to the full set with zeros. The
# solver answers on every state; a leaf may only return some of them.
#
# This is the one place the chain's direction axis and the solver's seed axis
# meet, and they are not the same thing: a seed column is its own functional and
# costs a whole sweep, a direction rides inside the dual. So slice 1 becomes the
# seed and the rest becomes its tangents.
.widenSeed <- function(w, full, subset) {
  n <- dim(w)[1L]
  K <- .ctK(w)
  W <- array(0, c(n, length(full), 1L),
             dimnames = list(NULL, full, NULL))
  hit <- intersect(subset, full)
  if (length(hit)) W[, hit, 1L] <- w[, hit, 1L, drop = FALSE]
  if (K > 1L) {
    # The tangent block is positional: it has no dimnames of its own, and the
    # seed's column order is `full`.
    tg <- array(0, c(n, length(full), 1L, K - 1L))
    if (length(hit)) tg[, match(hit, full), 1L, ] <- w[, hit, -1L, drop = FALSE]
    attr(W, "seedTangent") <- tg
  }
  W
}

# A vjp is also called from outside the chain, where the natural shape of a
# cotangent is the matrix or the named vector it was before the direction axis.
# Normalising at the leaf keeps both callers on one path.
.asCtOut  <- function(w) if (is.null(w) || length(dim(w)) == 3L) w else .ct(out = w)$out
.asCtPars <- function(w) if (is.null(w) || !is.null(dim(w))) w else .ct(pars = w)$pars

# A cotangent widened to K directions with zeros. The objective's seed is a
# constant of the functional it seeds, so its own tangents vanish and only the
# width has to travel.
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

# The "time" column is a label, not an output: nothing differentiates it and no
# leaf seeds it. Dropping it here keeps every vjp free of the special case.
.dropTime <- function(m) {
  if (is.null(m) || length(dim(m)) < 2L) return(m)
  keep <- dimnames(m)[[2L]] != "time"
  if (all(keep)) return(m)
  if (length(dim(m)) == 2L) m[, keep, drop = FALSE] else m[, keep, , drop = FALSE]
}

# A parfn whose Jacobian is a matrix it already builds -- Pimpl solves it by the
# implicit function theorem, Pequil reads it off the endpoint sensitivity of a
# nested steady-state solve. The vjp is that matrix transposed onto w.
#
# The split is deliberate and not a shortcut. What makes the forward mode
# expensive is that its width is n_theta, and the outer chain is where n_theta
# lives; a nested transformation's width is its own parameter set, which does
# not grow when the outer parametrisation does. So the trajectory goes backwards
# and the sub-problem stays forward, and the cost of the whole is still
# independent of n_theta.
#
# `deriv` is stripped off `pars` first: in a reverse chain the forward pass ran
# without one, and a leftover would chain the Jacobian to the outer parameters
# here instead of one node further down, where it belongs.
.parfnVjpFromJacobian <- function(p2p) {
  function(pars, fixed = NULL, w, condition = NULL) {
    w <- .asCtPars(w)
    K <- .ctK(w)
    # The incoming tangents, read before the strip below takes them off: at
    # second order they are what the node's curvature is contracted along.
    V <- if (K > 1L) attr(pars, "deriv") else NULL
    attr(pars, "deriv") <- NULL
    attr(pars, "deriv2") <- NULL
    v <- p2p(pars, fixed = fixed, deriv = TRUE, deriv2 = (K > 1L),
             condition = condition)
    J <- attr(v, "deriv")
    if (is.null(J) || !is.matrix(J))
      return(.ctZero(names(pars), K))
    wv <- .pickCotangent(w, rownames(J))
    u  <- crossprod(J, wv)
    rownames(u) <- colnames(J)
    if (K > 1L) {
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
