## selectEB on simulated linear models: the error rate of rule = "alpha", null
## candidates as their number grows, and the estimate of g. No compilation;
## minutes on several cores.
library(dMod2)
cores <- max(1L, parallel::detectCores() - 1L)

## Least squares on y = X beta + noise, every coefficient a reference parameter
linearProblem <- function(seed, J, n, beta, sigma = 1) {
  set.seed(seed)
  X <- cbind(b0 = 1, matrix(rnorm(n * J), n, J, dimnames = list(NULL, paste0("b", seq_len(J)))))
  y <- drop(X %*% c(0.5, beta)) + rnorm(n, 0, sigma)
  obj <- function(pars, fixed = NULL, ...) {
    b <- c(pars, fixed)[colnames(X)]
    r <- (drop(X %*% b) - y) / sigma
    Jp <- X[, names(pars), drop = FALSE] / sigma
    objlist(value = sum(r^2), gradient = setNames(drop(2 * crossprod(Jp, r)), names(pars)),
            hessian = 2 * crossprod(Jp))
  }
  list(obj = obj, ref = colnames(X)[-1], center = setNames(rep(0, J + 1), colnames(X)),
       X = X[, -1])
}
quick <- list(fullFits = 1, refitFits = 0, swaps = FALSE)

## rule = "alpha" with one null candidate: share of datasets that select it ------
selected <- unlist(parallel::mclapply(1:1000, function(s) {
  f <- linearProblem(s, 1, 30, 0)
  selectEB(f$obj, f$center, reference = f$ref, rule = "alpha", fits = 1,
           control = c(quick, seed = s))$terms$selected
}, mc.cores = cores))
mean(selected)                                 # about alpha = 0.05
qbinom(c(0.005, 0.995), 1000, 0.05) / 1000     # 99 % range around it

## Only null candidates, J = 5, 20, 40: share with any candidate in the best structure
falsePos <- sapply(c(J5 = 5, J20 = 20, J40 = 40), function(J) mean(unlist(
  parallel::mclapply(1:100, function(s) {
    f <- linearProblem(s, J, 100, rep(0, J))
    any(selectEB(f$obj, f$center, reference = f$ref, fits = 2,
                 control = c(quick, seed = s))$best$z)
  }, mc.cores = cores))))
falsePos                                       # stays at or below 0.05

## 20 of 60 effects drawn from a g-prior with g = 50 -------------------------------
g0 <- 50
gRes <- do.call(rbind, parallel::mclapply(1:20, function(s) {
  set.seed(1000 + s)
  f <- linearProblem(s, 60, 200, rep(0, 60))
  on <- 1:20
  beta <- numeric(60)
  beta[on] <- drop(t(chol(g0 * solve(crossprod(f$X[, on])))) %*% rnorm(20))
  y <- drop(f$X %*% beta) + rnorm(200)
  X <- f$X
  obj <- function(pars, fixed = NULL, ...) {
    r <- drop(X %*% c(pars, fixed)[colnames(X)]) - y
    Jp <- X[, names(pars), drop = FALSE]
    objlist(value = sum(r^2), gradient = setNames(drop(2 * crossprod(Jp, r)), names(pars)),
            hessian = 2 * crossprod(Jp))
  }
  fit <- selectEB(obj, setNames(rep(0, 60), colnames(X)), reference = colnames(X),
                  fits = 4, control = c(quick, seed = s))
  zs <- abs(beta[on]) * sqrt(diag(crossprod(X[, on])))
  data.frame(g = fit$best$g, strong = mean(fit$best$z[on][zs > 5]), fpr = mean(fit$best$z[-on]))
}, mc.cores = cores))
median(gRes$g)                                 # near g0 = 50
mean(gRes$strong, na.rm = TRUE)                # strong effects found, above 0.95
mean(gRes$fpr)                                 # null effects kept out, below 0.02
