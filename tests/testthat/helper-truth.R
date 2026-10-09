# Pure-R closed-form references, no compilation: tests check the compiled
# machinery against the mathematics rather than another implementation.


## ---- Closed-form solutions ----------------------------------------------

# Closed-form linear decay with its first and second derivatives in (x0, k),
# vectorised over t.
truth_decay <- function(t, x0, k) {
  E <- exp(-k * t)
  value <- x0 * E
  grad  <- cbind(x0 = E, k = -t * x0 * E)
  hess <- array(0, c(length(t), 2L, 2L),
                dimnames = list(NULL, c("x0", "k"), c("x0", "k")))
  hess[, "x0", "k"] <- -t * E
  hess[, "k", "x0"] <- -t * E
  hess[, "k", "k"]  <- t^2 * x0 * E
  list(value = value, gradient = grad, hessian = hess)
}

# Closed-form two-step linear cascade started from A only. Equal rates take
# the degenerate branch of the formula.
truth_two_step <- function(t, A0, k1, k2) {
  A <- A0 * exp(-k1 * t)
  if (isTRUE(all.equal(k1, k2))) {
    B <- A0 * k1 * t * exp(-k1 * t)
  } else {
    B <- A0 * k1 / (k2 - k1) * (exp(-k1 * t) - exp(-k2 * t))
  }
  C <- A0 - A - B
  data.frame(time = t, A = A, B = B, C = C)
}


## ---- Data simulation -----------------------------------------------------

# Simulate a noisy data frame in dMod's expected layout. truth_fn is a
# function (times, pars) -> numeric vector of length(times).
make_noisy_data <- function(truth_fn, pars, times,
                            name = "y",
                            sigma = 0.05,
                            condition = "C1",
                            lloq = NULL,
                            seed = 1L) {
  set.seed(seed)
  vals_true <- truth_fn(times, pars)
  noise <- rnorm(length(times), mean = 0, sd = sigma)
  out <- data.frame(
    name      = name,
    time      = times,
    value     = vals_true + noise,
    sigma     = sigma,
    condition = condition,
    stringsAsFactors = FALSE
  )
  if (!is.null(lloq)) out$lloq <- lloq
  out
}


## ---- numDeriv wrappers ---------------------------------------------------

# numDeriv with Richardson extrapolation; skips when numDeriv is absent.
numderiv_grad <- function(fn, x, ...) {
  if (!requireNamespace("numDeriv", quietly = TRUE))
    testthat::skip("numDeriv not installed")
  numDeriv::grad(fn, x, method = "Richardson", ...)
}

numderiv_hess <- function(fn, x, ...) {
  if (!requireNamespace("numDeriv", quietly = TRUE))
    testthat::skip("numDeriv not installed")
  numDeriv::hessian(fn, x, method = "Richardson", ...)
}


## ---- Closed-form objective values ---------------------------------------

# Reference value of nll_ALOQ, -2 log-likelihood of above-LOQ Gaussian data.
truth_nll_aloq <- function(pred, obs, sigma) {
  wr <- (pred - obs) / sigma
  sum(wr^2) + sum(log(2 * pi * sigma^2))
}

# Reference value of nll_BLOQ under M3, pred taken at the BLOQ time points.
truth_nll_bloq_m3 <- function(pred, lloq, sigma) {
  wr <- (pred - lloq) / sigma
  -2 * sum(stats::pnorm(-wr, log.p = TRUE))
}

# Reference M4NM / M4BEAL BLOQ contribution, without the stability fallbacks
# the package applies at extreme arguments.
truth_nll_bloq_m4 <- function(pred, lloq, sigma) {
  wr <- (pred - lloq) / sigma
  w0 <- pred / sigma
  -2 * sum(log(1 - stats::pnorm(wr) / stats::pnorm(w0)))
}
