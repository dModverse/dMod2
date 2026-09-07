# Behavioral tests for reduceReplicates() and fitErrorModel().
#
# reduceReplicates aggregates multi-replicate (name, time, condition) data
# down to a per-group mean and standard error. fitErrorModel fits a
# parametric variance model to the reduced data.
#
# Tests:
#   reduceReplicates: synthetic n-replicate samples with known mean and
#     standard deviation; the reduction should produce mean(value) and
#     standard error (= sd / sqrt(n)) up to a tolerance set by N(0,1)
#     sampling noise.
#   fitErrorModel: feed the reduced data to a constant-variance error
#     model and check that the recovered sigma equals the simulated
#     sigma_true within the n-replicates sampling envelope.


## ---- reduceReplicates: mean and standard error -----------------------

test_that("reduceReplicates returns the sample mean and standard error per group", {
  set.seed(42)
  sigma_true <- 0.1
  n_rep <- 5
  # 3 time points, single condition, 5 replicates each.
  raw <- do.call(rbind, lapply(c(1, 2, 5), function(tt) {
    data.frame(name = "A", time = tt,
               value = 1.0 + rnorm(n_rep, sd = sigma_true),
               condition = "C1",
               stringsAsFactors = FALSE)
  }))
  red <- reduceReplicates(raw)

  # Each (time, condition) group is a separate row.
  expect_equal(nrow(red), 3)
  expect_setequal(red$time, c(1, 2, 5))

  for (tt in c(1, 2, 5)) {
    row <- red[red$time == tt, ]
    grp <- raw[raw$time == tt, "value"]
    # `value` column is the mean.
    expect_equal(row$value, mean(grp), tolerance = 1e-12)
    # `sigma` is the standard error (sd / sqrt(n)).
    expect_equal(row$sigma, sd(grp) / sqrt(n_rep), tolerance = 1e-10)
    expect_equal(row$n, n_rep)
  }
})


## ---- fitErrorModel: recovers sigma_true on constant-variance data ----

test_that("fitErrorModel recovers exp(s0) ~ sigma_true^2 for constant variance", {
  testthat::skip_if_not_installed("optimx")
  set.seed(123)
  sigma_true <- 0.15
  n_rep <- 20
  # Many time points, many replicates, single condition: lots of evidence
  # to pin down the variance.
  raw <- do.call(rbind, lapply(seq(0.5, 5, by = 0.5), function(tt) {
    data.frame(name = "A", time = tt,
               value = 1.0 + rnorm(n_rep, sd = sigma_true),
               condition = "C1",
               stringsAsFactors = FALSE)
  }))
  red <- reduceReplicates(raw)

  # Constant variance error model: sigma^2 = exp(s0). Optimum: s0 ~ log(sigma_true^2).
  fit <- fitErrorModel(red, factors = "condition",
                       errorModel = "exp(s0)",
                       par = c(s0 = log(0.01)),
                       plotting = FALSE, blather = TRUE)
  s0_hat <- unique(fit$s0)[1]
  sigma_hat_sq <- exp(s0_hat)
  # Variance of sample variance with N-1 dof ~= 2 * sigma^4 / (N - 1).
  # Across 10 time points * 20 replicates we expect ~5% accuracy.
  expect_lt(abs(sqrt(sigma_hat_sq) - sigma_true) / sigma_true, 0.20)
})


# ============================================================================
# as.data.frame(<prdlist>, errfn = )
# ============================================================================

test_that("as.data.frame joins the error model by observable, not by row order", {
  skip_if_not_installed("cppDE")
  skip_on_cran()
  withr::local_dir(tempdir())

  reactions <- addReaction(eqnlist(), "A", "B", "k*A")
  m <- odemodel(reactions, modelname = "adf_ode", compile = FALSE)
  x <- Xs(m)
  g <- Y(eqnvec(obsA = "A", obsB = "B"), x, attach.input = FALSE,
         modelname = "adf_obs", compile = FALSE)
  # the error model lists its observables the other way round
  e <- Y(eqnvec(obsB = "sB", obsA = "sA"), g, attach.input = FALSE,
         modelname = "adf_err", compile = FALSE)
  compile(x, g, e, cores = 1)

  pars <- c(A = 2, B = 0, k = 0.3, sA = 0.11, sB = 0.77)
  prd  <- (g * x)(seq(0, 3, length.out = 4), pars)
  df   <- as.data.frame(prd, errfn = e)

  expect_equal(unique(df$sigma[df$name == "obsA"]), 0.11)
  expect_equal(unique(df$sigma[df$name == "obsB"]), 0.77)

  # an unnamed prediction numbers its conditions rather than losing the column
  names(prd) <- NULL
  bare <- as.data.frame(prd)
  expect_true("condition" %in% names(bare))
  expect_equal(unique(bare$condition), "1")
})


test_that("wide2long forwards keep and na.rm to the matrix method", {
  m <- matrix(c(0, 1, 2, NA, 4, 5), ncol = 3,
              dimnames = list(NULL, c("time", "a", "b")))
  out <- wide2long(list(C1 = m), na.rm = TRUE)
  expect_false(anyNA(out$value))
  expect_equal(nrow(out), 3L)
  expect_equal(unique(out$condition), "C1")
})
