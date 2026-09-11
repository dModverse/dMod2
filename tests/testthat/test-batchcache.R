# The prepared cppDE batch handle in Xs.cppDE names the shared object it
# resolved its entry point from. Renaming the model, or reopening a saved
# workspace where that object does not exist, leaves the cache pointing at a
# .so that cannot be called.

test_that("a batch handle whose shared object is gone is re-prepared", {

  skip_on_cran()

  fx <- fx_decay_multicond_compiled()
  times <- seq(0, 5, by = 1)

  invisible(fx$prd(times, fx$outerpars, deriv = TRUE))

  bcache <- environment(attr(fx$xfn, "mappings")[[1]])$bcache
  expect_false(is.null(bcache$handle))

  # What the cluster node sees: the handle survived, its shared object did not.
  bcache$handle$sym$dll <- "no_such_shared_object"

  out <- expect_silent(fx$prd(times, fx$outerpars, deriv = TRUE))
  expect_false(identical(bcache$handle$sym$dll, "no_such_shared_object"))

  ref <- exp(fx$outerpars["s_C1_log"]) * exp(-0.5 * times)
  expect_equal(unname(out$C1[, "y"]), unname(ref), tolerance = 1e-5)

})

test_that("modelname<- drops the cached batch handle", {

  skip_on_cran()

  fx <- fx_decay_multicond_compiled()
  invisible(fx$prd(seq(0, 5, by = 1), fx$outerpars, deriv = TRUE))

  bcache <- environment(attr(fx$xfn, "mappings")[[1]])$bcache
  expect_false(is.null(bcache$handle))

  x <- fx$xfn
  modelname(x) <- "renamed_shared_object"
  expect_true(is.null(bcache$handle))

  # And the next call rebuilds it against whatever is loaded now.
  out <- fx$prd(seq(0, 5, by = 1), fx$outerpars, deriv = TRUE)
  expect_false(is.null(bcache$handle))
  expect_true(all(is.finite(out$C1[, "y"])))

})
