# selectEB() on simulated data: error rates, multiplicity, estimate of g and
# recovery of an ODE mechanism. Minutes on several cores, run with
# DMOD_LONGTESTS=1.

.eb_long <- function() {
  skip_if_not(identical(Sys.getenv("DMOD_LONGTESTS"), "1"), "DMOD_LONGTESTS is not 1")
  skip_on_cran()
}

.eb_linprob <- function(seed, J, n, beta, sigma = 1) {
  set.seed(seed)
  X <- cbind(b0 = 1, matrix(rnorm(n * J), n, J, dimnames = list(NULL, paste0("b", seq_len(J)))))
  y <- drop(X %*% c(0.5, beta)) + rnorm(n, 0, sigma)
  list(obj = .lsq_obj(.linear_model(X), y, sigma), ref = colnames(X)[-1],
       center = stats::setNames(rep(0, J + 1), colnames(X)), X = X[, -1])
}

.eb_quick <- list(fullFits = 1, refitFits = 0, swaps = FALSE)

test_that("rule alpha keeps a missing term at the level alpha", {
  .eb_long()
  sel <- unlist(parallel::mclapply(1:1000, function(s) {
    f <- .eb_linprob(s, 1, 30, 0)
    selectEB(f$obj, f$center, reference = f$ref, rule = "alpha", fits = 1,
             control = c(.eb_quick, seed = s))$terms$selected
  }, mc.cores = test_cores()))
  ci <- stats::qbinom(c(0.005, 0.995), 1000, 0.05) / 1000
  expect_gte(mean(sel), ci[1])
  expect_lte(mean(sel), ci[2])
})

test_that("null candidates stay out as their number grows", {
  .eb_long()
  for (J in c(5, 20, 40)) {
    fp <- unlist(parallel::mclapply(1:100, function(s) {
      f <- .eb_linprob(s, J, 100, rep(0, J))
      any(selectEB(f$obj, f$center, reference = f$ref, fits = 2,
                   control = c(.eb_quick, seed = s))$best$z)
    }, mc.cores = test_cores()))
    expect_lte(mean(fp), 0.05)
  }
})

test_that("g is estimated as the squared effect size in units of its standard error", {
  .eb_long()
  g0 <- 50
  res <- do.call(rbind, parallel::mclapply(1:20, function(s) {
    set.seed(1000 + s)
    f <- .eb_linprob(s, 60, 200, rep(0, 60))
    on <- 1:20
    beta <- numeric(60)
    beta[on] <- drop(t(chol(g0 * solve(crossprod(f$X[, on])))) %*% rnorm(20))
    y <- drop(f$X %*% beta) + rnorm(200)
    obj <- .lsq_obj(.linear_model(f$X), y, 1)
    r <- selectEB(obj, stats::setNames(rep(0, 60), colnames(f$X)), reference = colnames(f$X),
                  fits = 4, control = c(.eb_quick, seed = s))
    zs <- abs(beta[on]) * sqrt(diag(crossprod(f$X[, on])))
    data.frame(g = r$best$g, strong = mean(r$best$z[on][zs > 5]), fpr = mean(r$best$z[-on]))
  }, mc.cores = test_cores()))
  expect_gt(stats::median(res$g), g0 / 2)
  expect_lt(stats::median(res$g), 2 * g0)
  expect_gt(mean(res$strong, na.rm = TRUE), 0.95)
  expect_lt(mean(res$fpr), 0.02)
})

test_that("the reactions of an ODE network are recovered", {
  .eb_long()
  skip_if_not_installed("cppDE")
  dir <- tempfile("selectEBnet"); dir.create(dir)
  old <- setwd(dir); on.exit(setwd(old), add = TRUE)
  f <- eqnvec(A = "-k1*A - k4*A - k6*A + k7*C", B = "k1*A - k2*B - k5*B + k8*D",
              C = "k2*B + k4*A - k3*C - k7*C", D = "k3*C + k5*B + k6*A - k8*D")
  x <- Xs(odemodel(f, modelname = "ebnet_x", compile = FALSE, outdir = dir), compile = FALSE)
  g <- Y(eqnvec(oA = "A", oB = "B", oC = "C", oD = "D"), f, modelname = "ebnet_g",
         compile = FALSE, attachInput = FALSE, outdir = dir)
  ks <- paste0("k", 1:8)
  trafo <- gateL1(eqnvec(A = "1", B = "0", C = "0", D = "0",
                         stats::setNames(paste0("10^log10_", ks), ks)), ks)
  p <- P(trafo, condition = "C1", modelname = "ebnet_p", compile = FALSE, outdir = dir)
  compile(g, x, p, output = "ebnet_model")
  prd <- g * x * p
  truth <- c(stats::setNames(rep(-1, 8), paste0("log10_", ks)),
             stats::setNames(c(1, 1, 1, rep(0, 5)), paste0("s_", ks)))
  truth[c("log10_k1", "log10_k2", "log10_k3")] <- log10(c(0.5, 0.3, 0.2))
  center <- stats::setNames(rep(-0.5, 8), paste0("log10_", ks))
  key <- "-s_k4 -s_k5 -s_k6 -s_k7 -s_k8"
  hit <- unlist(parallel::mclapply(1:10, function(s) {
    set.seed(s)
    sim <- wide2long(prd(seq(0.5, 15, 0.5), truth, deriv = FALSE))
    sim$value <- sim$value + rnorm(nrow(sim), 0, 0.02)
    sim$sigma <- 0.02
    obj <- normL2(as.datalist(sim), prd)
    selectEB(obj, center, zero = attr(trafo, "gates"), fits = 6,
             control = list(seed = s))$selected == key
  }, mc.cores = test_cores()))
  expect_gte(mean(hit), 0.9)
})
