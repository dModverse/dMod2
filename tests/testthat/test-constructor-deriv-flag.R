# derivMode selects what a constructor builds; deriv, deriv2 and sweep of a call
# select what the call computes, and a call that asks for what was not built errors.

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

    pexpl <- Pexpl(c(A = "a * x", B = "b + y"), modelname = nm("test_pexpl_nod1"),
                   derivMode = "none", verbose = FALSE)
    pimpl <- Pimpl(c(x = "x - a"), parameters = "a", derivMode = "none",
                   modelname = nm("test_pimpl_nod1"), verbose = FALSE)
    gfn <- Y(c(obs = "k * A"), states = c("A", "time"), parameters = "k",
             modelname = nm("test_y_nod1"), derivMode = "none", verbose = FALSE,
             attachInput = FALSE)
    pdisp <- P(c(A = "a * x"), method = "explicit", derivMode = "none",
               modelname = nm("test_P_nod1"), verbose = FALSE)

    f <- eqnvec(A = "-k1*A", B = "k1*A - k2*B")
    mff <- odemodel(f, modelname = nm("dm_ode"), derivMode = "forward-forward",
                    compile = FALSE)
    mrev <- odemodel(f, modelname = nm("dm_rev"), derivMode = "reverse", compile = FALSE)
    gff <- Y(c(obsB = "s*B"), f, modelname = nm("dm_obs"), derivMode = "forward-forward")
    tr <- c(A = "exp(lA)", B = "0", k1 = "exp(lk1)", k2 = "exp(lk2)", s = "exp(ls)")
    pff <- P(tr, condition = "C1", modelname = nm("dm_p"), derivMode = "forward-forward")
    pfw <- P(tr, condition = "C1", modelname = nm("dm_pf"))
    xff <- Xs(mff, options = list(atol = 1e-10, rtol = 1e-10))
    xrev <- Xs(mrev)

    compile(pexpl, pimpl, gfn, pdisp, xff, xrev, gff, pff, pfw,
            output = nm("ctor_models"), cores = test_cores())

    data <- as.datalist(data.frame(name = "obsB", time = c(0.5, 1, 2, 4),
                                   value = c(0.8, 1.1, 0.9, 0.5), sigma = 0.1,
                                   condition = "C1"))
    cache <<- list(pexpl = pexpl, pimpl = pimpl, gfn = gfn, pdisp = pdisp,
                   xff = xff, xrev = xrev, gff = gff, pff = pff, pfw = pfw, data = data,
                   pars = c(lA = log(2), lk1 = log(0.6), lk2 = log(0.3), ls = log(1.5)))
    cache
  }
})


test_that("constructors reject deriv and deriv2 as build arguments", {
  expect_error(odemodel(eqnvec(A = "-k*A"), deriv = FALSE, compile = FALSE),
               "are not build arguments; derivMode selects")
  expect_error(P(c(A = "k"), method = "explicit", deriv2 = TRUE),
               "`deriv2` is not an argument; derivMode selects")
  expect_error(P(c(A = "k"), method = "implicit", deriv = FALSE),
               "are not build arguments; derivMode selects")
  expect_error(Pexpl(c(A = "x"), deriv2 = TRUE), "unused argument")
  expect_error(Y(c(obs = "A"), states = c("A", "time"), deriv2 = TRUE), "unused argument")
  # `deriv` matches derivMode partially, a logical value names the change
  expect_error(Pimpl(c(x = "x - a"), parameters = "a", deriv = FALSE),
               "`deriv` and `deriv2` are not build arguments")
  expect_error(Y(c(obs = "A"), states = c("A", "time"), deriv = TRUE),
               "`deriv` and `deriv2` are not build arguments")
})


test_that("derivMode = \"none\" combines with no other direction", {
  expect_error(Pexpl(c(A = "x"), derivMode = c("none", "forward")), "none")
  expect_error(Y(c(obs = "A"), states = c("A", "time"), derivMode = c("reverse", "none")),
               "none")
  expect_error(P(c(A = "k"), method = "explicit", derivMode = character(0)), "derivMode")
})


test_that("functions built with derivMode = \"none\" return values and refuse derivatives", {
  skip_if_no_compile()
  m <- ctor_models()
  pars <- c(a = 2, b = 3, x = 4, y = 5)
  out <- m$pexpl(pars, deriv = FALSE)
  expect_equal(unname(unclass(out[[1]])[c("A", "B")]), c(8, 8))
  expect_null(attr(out[[1]], "deriv"))
  expect_error(m$pexpl(pars), "Pexpl\\(\\) was built without derivMode = \"forward\"")

  expect_equal(unname(unclass(m$pimpl(c(a = 1.5, x = 0.5), deriv = FALSE)[[1]])["x"]), 1.5)
  expect_error(m$pimpl(c(a = 1.5, x = 0.5)), "Pimpl\\(\\) was built without derivMode = \"forward\"")

  expect_null(attr(m$pdisp(c(a = 2, x = 3), deriv = FALSE)[[1]], "deriv"))
  expect_error(m$pdisp(c(a = 2, x = 3)), "derivMode = \"forward\"")

  prd <- structure(cbind(time = c(0, 1), A = c(1, 2)),
                   parameters = structure(c(k = 1.5), fixed = NULL),
                   class = c("prdframe", "matrix", "array"))
  res <- m$gfn(out = prd, pars = c(k = 1.5), deriv = FALSE)
  expect_equal(unname(res[[1]][, "obs"]), c(1.5, 3))
  expect_error(m$gfn(out = prd, pars = c(k = 1.5)), "Y\\(\\) was built without derivMode = \"forward\"")
})


test_that("a call asking for an unbuilt direction names the missing derivMode", {
  skip_if_no_compile()
  m <- ctor_models()
  times <- c(0, 1, 2)
  inner <- c(A = 2, B = 0, k1 = 0.6, k2 = 0.3)
  expect_error(m$xrev(times, inner), "derivMode = \"forward\"")
  expect_null(attr(m$xrev(times, inner, deriv = FALSE)[[1]], "deriv"))
  expect_error(m$pfw(m$pars, deriv2 = TRUE), "derivMode = \"forward-forward\"")
  obj <- normL2(m$data, m$gff * m$xff * m$pff)
  expect_error(obj(m$pars, sweep = "reverse"), "derivMode = \"reverse\"")
})


test_that("derivMode = \"forward-forward\" gives the Hessian of deriv2 = TRUE in dMod2 0.16.0", {
  skip_if_no_compile()
  m <- ctor_models()
  obj <- normL2(m$data, m$gff * m$xff * m$pff)
  out <- obj(m$pars, deriv2 = TRUE)
  # value and exact Hessian of the same objective built with deriv2 = TRUE in dMod2 0.16.0
  ref <- matrix(c(1496.45262699787, 577.328918012139, -620.715428312784, 1496.45262699787,
                  577.328918012139, 164.606366756938, -196.700167236959, 577.328918012139,
                  -620.715428312784, -196.700167236959, 167.720578617145, -620.715428312784,
                  1496.45262699787, 577.328918012139, -620.715428312784, 1496.45262699787),
                4, 4, dimnames = rep(list(c("lA", "lk1", "lk2", "ls")), 2))
  expect_equal(out$value, 82.3524114951828, tolerance = 1e-10)
  expect_equal(out$hessian[rownames(ref), colnames(ref)], ref, tolerance = 1e-10)
})
