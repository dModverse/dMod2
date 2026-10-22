# Prediction functions without compiled code: Xd(), predict(), prdlist() and %.*%.

.xdGrid <- data.frame(name = "A", time = c(0, 1, 2), row.names = c("a0", "a1", "a2"))

# a1 = exp(t1); fixed parameters leave the Jacobian
.xdP2p <- function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE) {
  allp <- c(unclass(pars), unclass(fixed))
  out <- c(a0 = allp[["a0"]], a1 = exp(allp[["t1"]]), a2 = allp[["a2"]])
  J <- diag(c(1, exp(allp[["t1"]]), 1))
  dimnames(J) <- list(names(out), c("a0", "t1", "a2"))
  as.parvec(out, deriv = J[, setdiff(colnames(J), names(fixed)), drop = FALSE])
}

test_that("Xd applies the chain rule of the parameter transformation", {
  x <- Xd(.xdGrid, condition = "C1")
  p <- parfn(.xdP2p, c("a0", "t1", "a2"), "C1")
  th <- c(a0 = 1, t1 = log(3), a2 = 2)
  d <- attr((x * p)(c(0.5, 1), th)[[1]], "deriv")
  expect_equal(dimnames(d)[[3]], c("a0", "t1", "a2"))
  expect_equal(unname(d[, "A", "t1"]), c(0.5, 1) * 3)
  expect_equal(unname(d[, "A", "a0"]), c(0.5, 0))
})

test_that("Xd gives derivatives with respect to the free parameters only", {
  x <- Xd(.xdGrid)
  d <- attr(x(c(0.5, 1.5), c(a1 = 2, a2 = 3), fixed = c(a0 = 1))[[1]], "deriv")
  expect_equal(dimnames(d)[[3]], c("a1", "a2"))
  expect_equal(unname(d[, "A", "a1"]), c(0.5, 0.5))
})

test_that("Xd names a missing grid parameter", {
  x <- Xd(.xdGrid)
  expect_error(x(c(0, 1), c(a0 = 1, a2 = 2)), "a1")
})

test_that("predict() labels every prediction with its own parameter row", {
  x <- Xd(data.frame(name = "A", time = 0:1, row.names = c("a0", "a1")), "C1")
  pf <- parframe(data.frame(value = c(5, 1, 3), a0 = c(10, 20, 30), a1 = 0),
                 parameters = c("a0", "a1"), metanames = "value")
  r <- predict(x, times = 0, pars = pf)
  expect_equal(r$value[order(r$.value)], c(20, 30, 10))
})

test_that("predict() adds a condition grid with one covariate", {
  grid <- data.frame(name = "A", time = 0:2, row.names = c("a0", "a1", "a2"))
  x <- Xd(grid, "DEM") + Xd(grid, "APAP")
  pf <- parframe(data.frame(value = 1, a0 = 1, a1 = 2, a2 = 3),
                 parameters = c("a0", "a1", "a2"), metanames = "value")
  data <- as.datalist(data.frame(name = "A", time = 0:1, value = 1:2, sigma = 0.1,
                                 dose = c("DEM", "APAP")), splitBy = "dose")
  r <- predict(x, times = c(0, 1), pars = pf, data = data)
  expect_true("dose" %in% colnames(r))
  expect_equal(r$dose, r$condition)
})

test_that("%.*% scales an objective and keeps absent derivatives absent", {
  o1 <- constraintL2(mu = c(a = 1, b = 2), sigma = 1)
  o2 <- constraintL2(mu = c(a = 0), sigma = 1)
  p <- c(a = 0.5, b = 1)
  s <- (2 %.*% o1 + o2)(p)
  r1 <- o1(p); r2 <- o2(p)
  expect_equal(s$value, 2 * r1$value + r2$value)
  expect_equal(s$gradient[["b"]], 2 * r1$gradient[["b"]])
  expect_null((2 %.*% o1)(p, deriv = FALSE)$gradient)
  expect_equal((2 %.*% o1 + o2)(p, deriv = FALSE)$value, s$value)
  expect_no_error((2 %.*% o1 + o2)(p, hessian = FALSE))
  expect_error(c(1, 2) %.*% o1, "scalar")
})

test_that("prdlist() names unnamed entries by position", {
  expect_length(prdlist(), 0L)
  expect_equal(names(prdlist(1, b = 2)), c("1", "b"))
})

test_that("dMod functions refuse unknown, surplus and missing arguments", {
  x <- Xd(.xdGrid, "C1") + Xd(.xdGrid, "C2")
  p <- c(a0 = 1, a1 = 2, a2 = 3)
  expect_error(x(0:1, p, condition = "C1"), "condition")
  expect_error(x(0:1, p, 7), "too many")
  expect_error(x(0:1), "pars")
  expect_named(x(0:1, p, conditions = "C1"), "C1")
  obj <- constraintL2(mu = p, sigma = 1)
  expect_error(obj(p, fixd = c(a0 = 1)), "fixd")
})

test_that("an empty prdlist prints", {
  expect_no_error(capture.output(print(prdlist())))
})

test_that("modelname<- names each condition and honours conditions", {
  mk <- function(cond, mn) {
    f <- function(...) NULL
    attr(f, "modelname") <- mn
    prdfn(f, parameters = "a", condition = cond)
  }
  x <- mk("C1", "m1") + mk("C2", "m2") + mk("C3", "m3")
  y <- x
  modelname(y) <- c("n1", "n2", "n3")
  expect_equal(unname(vapply(attr(y, "mappings"), attr, "", "modelname")), c("n1", "n2", "n3"))
  y <- x
  modelname(y, conditions = "C2") <- "new"
  expect_equal(unname(vapply(attr(y, "mappings"), attr, "", "modelname")), c("m1", "new", "m3"))
})

test_that("getParameters of a function for every condition takes conditions", {
  grid <- data.frame(name = "A", time = 0:2, row.names = paste0("A", 0:2))
  x <- Xd(grid)
  expect_equal(getParameters(x, conditions = "C1"), getParameters(x))
})
