# Fixtures of the scanL1() tests: least squares objectives without compilation.

# Weighted least squares over a model y = f(theta) with Jacobian, as an
# objective of the shape trust() and trustL1() call: value, gradient and
# Gauss-Newton Hessian over the free parameters.
.lsq_obj <- function(model, data, sigma) {
  function(pars, fixed = NULL, ...) {
    all <- c(pars, fixed)
    m   <- model(all)
    r   <- (m$y - data) / sigma
    J   <- m$J[, names(pars), drop = FALSE] / sigma
    objlist(value = sum(r^2), gradient = setNames(drop(2 * crossprod(J, r)), names(pars)),
            hessian = 2 * crossprod(J))
  }
}

# Linear model whose coefficients are a log10 magnitude times a gate.
.gated_model <- function(X) function(p) {
  k  <- paste0("k", seq_len(ncol(X)))
  a  <- 10^p[paste0("log10_", k)]; s <- p[paste0("s_", k)]
  J  <- cbind(sweep(X, 2, log(10) * a * s, `*`), sweep(X, 2, a, `*`))
  colnames(J) <- c(paste0("log10_", k), paste0("s_", k))
  list(y = drop(X %*% (a * s)), J = J)
}

# Cell types around one common level, with fold changes to a reference type.
.type_model <- function(type) function(p) {
  r <- c(R = 0, p[c("r_B", "r_C", "r_D")])
  names(r) <- c("R", "B", "C", "D")
  J <- cbind(mu = 1, sapply(c("B", "C", "D"), function(t) as.numeric(type == t)))
  colnames(J) <- c("mu", "r_B", "r_C", "r_D")
  list(y = p[["mu"]] + r[type], J = J)
}

# Exponential decay per cell type with fold changes on the log10 rate: the toy
# model of Hauber et al. (2023), solved analytically.
.decay_model <- function(type, times, types = unique(type)) function(p) {
  rn <- paste0("r_", types[-1])
  r  <- c(0, p[rn]); names(r) <- types
  k  <- 10^(p[["lk"]] + r[type])
  y  <- exp(-k * times)
  dk <- -times * y * k * log(10)
  J  <- cbind(lk = dk, sapply(types[-1], function(t) dk * (type == t)))
  colnames(J) <- c("lk", rn)
  list(y = y, J = J)
}
