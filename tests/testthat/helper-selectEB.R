# Fixtures of the selectEB() tests: least squares objectives without
# compilation, on top of .lsq_obj() from helper-scanL1.R.

# Linear model y = X beta, every coefficient a reference parameter.
.linear_model <- function(X) function(p) {
  b <- p[colnames(X)]
  list(y = drop(X %*% b), J = X)
}

# Intercept plus gated coefficients: y = b0 + X (10^log10_k * s_k).
.gated_lin_model <- function(X) function(p) {
  k <- paste0("k", seq_len(ncol(X)))
  a <- 10^p[paste0("log10_", k)]; s <- p[paste0("s_", k)]
  J <- cbind(b0 = 1, sweep(X, 2, log(10) * a * s, `*`), sweep(X, 2, a, `*`))
  colnames(J) <- c("b0", paste0("log10_", k), paste0("s_", k))
  list(y = p[["b0"]] + drop(X %*% (a * s)), J = J)
}

# Gated linear problem with `on` the present candidates; returns objective,
# gates, center and the key of the truth.
.gated_lin_problem <- function(seed, J = 8, n = 40, on = c(1, 3, 6),
                               k = c(1, 0.5, 2), sigma = 0.1) {
  set.seed(seed)
  X  <- matrix(stats::runif(n * J), n, J)
  ks <- paste0("k", seq_len(J))
  truth <- c(b0 = 0.3, stats::setNames(rep(-1, J), paste0("log10_", ks)),
             stats::setNames(rep(0, J), paste0("s_", ks)))
  truth[paste0("log10_", ks[on])] <- log10(k)
  truth[paste0("s_", ks[on])] <- 1
  m <- .gated_lin_model(X)
  y <- m(truth)$y + stats::rnorm(n, 0, sigma)
  list(obj = .lsq_obj(m, y, sigma),
       zero = stats::setNames(as.list(paste0("log10_", ks)), paste0("s_", ks)),
       center = c(b0 = 0, stats::setNames(rep(0, J), paste0("log10_", ks))),
       truth = paste0("-", paste(paste0("s_", ks[-on]), collapse = " -")))
}

# Two parallel paths A -> B with gated rates and a gated linear drift:
# y = 1 - exp(-(ka s_a + kb s_b) t) + kc s_c t. The paths can replace each other.
.paths_model <- function(t) function(p) {
  ka <- 10^p[["la"]]; kb <- 10^p[["lb"]]; kc <- 10^p[["lc"]]
  r  <- ka * p[["s_a"]] + kb * p[["s_b"]]
  e  <- exp(-r * t)
  J  <- cbind(la = t * e * log(10) * ka * p[["s_a"]], lb = t * e * log(10) * kb * p[["s_b"]],
              lc = log(10) * kc * p[["s_c"]] * t,
              s_a = t * e * ka, s_b = t * e * kb, s_c = kc * t)
  list(y = 1 - e + kc * p[["s_c"]] * t, J = J)
}

# Exponential decay per cell line with amplitude and rate, each as a log10
# fold change to the reference line R.
.lines_model <- function(line, times, lines) function(p) {
  ra <- c(0, p[paste0("ra_", lines[-1])]); names(ra) <- lines
  rk <- c(0, p[paste0("rk_", lines[-1])]); names(rk) <- lines
  A <- 10^(p[["la"]] + ra[line]); k <- 10^(p[["lk"]] + rk[line])
  e <- exp(-k * times)
  dA <- log(10) * A * e
  dk <- -log(10) * A * times * k * e
  J <- cbind(la = dA, lk = dk,
             sapply(lines[-1], function(l) dA * (line == l)),
             sapply(lines[-1], function(l) dk * (line == l)))
  colnames(J) <- c("la", "lk", paste0("ra_", lines[-1]), paste0("rk_", lines[-1]))
  list(y = A * e, J = J)
}
