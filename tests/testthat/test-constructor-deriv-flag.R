# Constructor-level `deriv = TRUE/FALSE` gating for P, Pexpl, Pimpl,
# Pequil, and Y. Symmetric with the existing `deriv2` flag: the
# constructor decides whether the artifact carries first-order
# sensitivities, and the runtime call errors out if it asks for
# something the construction didn't produce.

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("cppDE")
  testthat::skip_on_cran()
}

# The compiled models of this file, generated with compile = FALSE on first
# use and linked into one shared object.
ctor_models <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    dir <- file.path(tempdir(), "ctor_models")
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    withr::local_dir(dir)
    nm <- function(x) paste0(x, "_", as.integer(Sys.time()))

    pexpl <- Pexpl(c(A = "a * x", B = "b + y"), deriv = FALSE,
                   modelname = nm("test_pexpl_nod1"), derivMode = "forward",
                   verbose = FALSE)
    pequil <- Pequil(c(A = "k_in - k_out * A"),
                     parameters = c("k_in", "k_out"), deriv = FALSE,
                     modelname = nm("test_pequil_nod1"), verbose = FALSE,
                     attach.input = FALSE)
    pimpl <- Pimpl(c(x = "x - a"), parameters = "a", deriv = FALSE,
                   modelname = nm("test_pimpl_nod1"), verbose = FALSE)
    gfn <- Y(c(obs = "k * A"), states = c("A", "time"), parameters = "k",
             deriv = FALSE, modelname = nm("test_y_nod1"),
             derivMode = "forward", verbose = FALSE, attach.input = FALSE)
    pdisp <- P(c(A = "a * x"), method = "explicit", deriv = FALSE,
               modelname = nm("test_P_nod1"), verbose = FALSE)

    compile(pexpl, pequil, pimpl, gfn, pdisp, output = nm("ctor_models"),
            cores = 4L)

    cache <<- list(pexpl = pexpl, pequil = pequil, pimpl = pimpl, gfn = gfn,
                   pdisp = pdisp)
    cache
  }
})


test_that("Pexpl(deriv = FALSE) yields a parvec without deriv attribute", {
  skip_if_no_compile()
  pf <- ctor_models()$pexpl
  out <- pf(c(a = 2, b = 3, x = 4, y = 5))
  expect_null(attr(out[[1]], "deriv"))
  # Default runtime deriv = TRUE is silently capped by the constructor:
  # no error, just no deriv attribute on the result.
  out2 <- pf(c(a = 2, b = 3, x = 4, y = 5), deriv = TRUE)
  expect_null(attr(out2[[1]], "deriv"))
})


test_that("Pequil(deriv = FALSE) skips the sensitivity model", {
  skip_if_no_compile()
  pf <- ctor_models()$pequil
  out <- pf(c(k_in = 1, k_out = 0.5, A = 0.1))
  expect_null(attr(out[[1]], "deriv"))
})


test_that("Pimpl(deriv = FALSE) drops the IFT chain rule from output", {
  skip_if_no_compile()
  pf <- ctor_models()$pimpl
  out <- pf(c(a = 1.5, x = 0.5))
  expect_null(attr(out[[1]], "deriv"))
})


test_that("Y(deriv = FALSE) produces output without deriv attribute", {
  skip_if_no_compile()
  gfn <- ctor_models()$gfn

  prd <- structure(
    cbind(time = c(0, 1), A = c(1, 2)),
    parameters = structure(c(k = 1.5), fixed = NULL),
    class = c("prdframe", "matrix", "array"))

  res <- gfn(out = prd, pars = c(k = 1.5))
  expect_null(attr(res[[1]], "deriv"))
})


test_that("Constructors reject deriv = FALSE combined with deriv2 = TRUE", {
  expect_error(Pexpl(c(A = "x"), deriv = FALSE, deriv2 = TRUE, compile = FALSE),
               "requires deriv = TRUE")
  expect_error(Pimpl(c(x = "x - a"), parameters = "a", deriv = FALSE, deriv2 = TRUE),
               "requires deriv = TRUE")
  expect_error(Pequil(c(A = "k - A"), parameters = "k", deriv = FALSE, deriv2 = TRUE),
               "requires deriv = TRUE")
  expect_error(Y(c(obs = "A"), states = c("A", "time"),
                 deriv = FALSE, deriv2 = TRUE),
               "requires deriv = TRUE")
  expect_error(P(c(A = "k"), method = "explicit",
                 deriv = FALSE, deriv2 = TRUE),
               "requires deriv = TRUE")
})


test_that("P() dispatcher forwards deriv to each method", {
  skip_if_no_compile()
  pf <- ctor_models()$pdisp
  out <- pf(c(a = 2, x = 3))
  expect_null(attr(out[[1]], "deriv"))
})
