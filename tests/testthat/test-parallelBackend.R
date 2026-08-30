# Smoke tests for the cross-platform parallel-apply helper. Mirrors the
# pattern used by mstrust / profile() (foreach + doParallel under the
# hood) so we can be confident the Windows path works without an actual
# Windows runner.


test_that(".parallelLapply with cores=1 falls back to lapply", {
  out <- dMod2:::.parallelLapply(1:5, function(i) i * 2, cores = 1L)
  expect_identical(out, lapply(1:5, function(i) i * 2))
})


test_that(".parallelLapply with cores>1 returns correct results", {
  # Either Unix fork via doParallel or PSOCK via makeCluster -- both go
  # through the same %dopar% loop and must produce identical results to
  # the serial path.
  out_par <- dMod2:::.parallelLapply(1:8, function(i) i^2, cores = 2L)
  out_ser <- lapply(1:8, function(i) i^2)
  expect_identical(out_par, out_ser)
})


test_that(".parallelLapply preserves order across workers", {
  # R CMD check caps cores at 2 via _R_CHECK_LIMIT_CORES_, so honour that
  # ceiling here. 2 workers across 20 items still exercises the order-
  # preservation property of the foreach %dopar% backend.
  n_cores <- if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) 2L else 4L
  set.seed(7L)
  X <- as.list(rnorm(20))
  out <- dMod2:::.parallelLapply(X, function(x) x + 1, cores = n_cores)
  expect_identical(out, lapply(X, function(x) x + 1))
})


test_that(".parallelLapply propagates worker errors", {
  # foreach surfaces worker errors as a regular R error -- the inner
  # stop("kaboom") should reach the caller via tryCatch.
  err <- tryCatch(dMod2:::.parallelLapply(1:3, function(i) {
    if (i == 2L) stop("kaboom") else i
  }, cores = 2L), error = function(e) e)
  expect_s3_class(err, "error")
})


## ---- Two-axis core budget ------------------------------------------------

test_that(".splitCores reads a scalar as the outer axis only", {
  r <- dMod2:::.splitCores(4L, "fits")
  expect_identical(r$outer, 4L)
  expect_null(r$conditions)

  expect_identical(dMod2:::.splitCores(NULL)$outer, 1L)
})


test_that(".splitCores reads both axes, by name or by position", {
  r <- dMod2:::.splitCores(c(fits = 3L, conditions = 2L), "fits")
  expect_identical(r$outer, 3L)
  expect_identical(r$conditions, 2L)

  # positional: first is the outer axis
  r <- dMod2:::.splitCores(c(3L, 2L), "fits")
  expect_identical(r$outer, 3L)
  expect_identical(r$conditions, 2L)

  # the outer axis is named after the caller's axis
  r <- dMod2:::.splitCores(c(pars = 5L, conditions = 2L), "pars")
  expect_identical(r$outer, 5L)
})


test_that(".splitCores rejects nonsense and warns on over-subscription", {
  expect_error(dMod2:::.splitCores(0L), "positive")
  expect_error(dMod2:::.splitCores(c(1L, 2L, 3L)), "length 1 or 2")
  expect_warning(dMod2:::.splitCores(c(fits = 10000L, conditions = 10000L), "fits"),
                 "exceeds")
})
