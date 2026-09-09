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

# The cotangent one node hands to the node below it.
.ct <- function(out = NULL, pars = NULL) list(out = out, pars = pars)

.ct_null <- function(w) is.null(w) || (is.null(w$out) && is.null(w$pars))

# Adds two named vectors on the union of their names. Cheaper than a merge and
# it keeps the order of the first, which is the order a caller expects back.
.addNamed <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  nms <- union(names(a), names(b))
  out <- setNames(numeric(length(nms)), nms)
  out[names(a)] <- out[names(a)] + a
  out[names(b)] <- out[names(b)] + b
  out
}

.addCt <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  o <- if (is.null(a$out)) b$out else if (is.null(b$out)) a$out else a$out + b$out
  .ct(o, .addNamed(a$pars, b$pars))
}

# A named vector on exactly `nms`, zero where the cotangent says nothing.
.pickCotangent <- function(w, nms) {
  out <- setNames(numeric(length(nms)), nms)
  if (is.null(w) || !length(nms)) return(out)
  hit <- intersect(names(w), nms)
  if (length(hit)) out[hit] <- w[hit]
  out
}

# A seed on a subset of columns, widened to the full set with zeros. The
# solver answers on every state; a leaf may only return some of them.
.widenSeed <- function(w, full, subset) {
  n <- nrow(w)
  W <- array(0, c(n, length(full), 1L),
             dimnames = list(NULL, full, NULL))
  hit <- intersect(subset, full)
  if (length(hit)) W[, hit, 1L] <- w[, hit, drop = FALSE]
  W
}

# The "time" column is a label, not an output: nothing differentiates it and no
# leaf seeds it. Dropping it here keeps every vjp free of the special case.
.dropTime <- function(m) {
  if (is.null(m)) return(NULL)
  if (is.matrix(m) && "time" %in% colnames(m)) m[, setdiff(colnames(m), "time"), drop = FALSE]
  else m
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
    attr(pars, "deriv") <- NULL
    attr(pars, "deriv2") <- NULL
    v <- p2p(pars, fixed = fixed, deriv = TRUE, deriv2 = FALSE,
             condition = condition)
    J <- attr(v, "deriv")
    if (is.null(J) || !is.matrix(J))
      return(setNames(numeric(length(pars)), names(pars)))
    wv <- .pickCotangent(w, rownames(J))
    .pickCotangent(setNames(as.numeric(crossprod(J, wv)), colnames(J)),
                   names(pars))
  }
}
