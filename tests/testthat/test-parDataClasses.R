# Parameter lists, parameter frames, parameter vectors and datalists, without compiled code.

.fitsFixture <- function() {
  f1 <- list(value = 3, argument = c(a = 1, b = 2), parinit = c(a = 0, b = 0), converged = TRUE,
             iterations = 5L, stopReason = "gradient", neval = 7L, elapsed = 0.5,
             phaseIterations = c(first = 2L, second = 3L),
             argpath = matrix(c(0, 0.5, 1, 0, 1, 2), 3, dimnames = list(NULL, c("a", "b"))))
  f2 <- list(value = 1, argument = c(a = 2, b = 3), parinit = c(a = 1, b = 1), converged = TRUE,
             iterations = 3L, stopReason = "gradient", neval = 4L,
             argpath = matrix(c(1, 2, 1, 3), 2, dimnames = list(NULL, c("a", "b"))))
  f3 <- list(value = 2, argument = c(a = 5, c = 9), parinit = c(a = 0, c = 0), converged = FALSE,
             iterations = 10L, stopReason = "iterlim")
  parlist(f1, structure("boom", class = "try-error"), f2, f3, list(error = "x"))
}

# dput with attributes in alphabetical order, at every level of a list
.canon <- function(x) {
  if (is.list(x)) x[] <- lapply(x, .canon)
  a <- attributes(x)
  if (!is.null(a)) attributes(x) <- a[order(names(a))]
  x
}
.canonDput <- function(x) dput(.canon(x))

.xdGrid <- function() {
  as.datalist(data.frame(name = rep(c("A", "B"), 4), time = rep(0:1, each = 2), value = 1:8,
                         sigma = 0.5, lloq = c(-Inf, -Inf, 3.5, -Inf),
                         dose = rep(c(1, 2), each = 4)), splitBy = "dose")
}

test_that("as.parlist and parlist", {
  expect_null(as.parlist())
  expect_null(as.parlist(NULL))
  pl <- as.parlist(list(list(value = 1)))
  expect_identical(class(pl), c("parlist", "list"))
  expect_identical(class(parlist()), c("parlist", "list"))
  expect_length(parlist(), 0L)
})

test_that("summary of a parlist prints best and worst fit and the outcome statistics", {
  expect_snapshot(summary(.fitsFixture()))
  expect_snapshot(summary(parlist(structure("boom", class = "try-error"))))
  expect_snapshot(summary(parlist()))
  fits <- .fitsFixture()
  expect_identical(withVisible(capture.output(s <- summary(fits)))$visible, TRUE)
  expect_identical(s, fits)
})

test_that("as.parframe of a parlist", {
  fits <- .fitsFixture()
  attr(fits, "run") <- "r1"
  expect_snapshot(.canonDput(as.parframe(fits)))
  expect_snapshot(.canonDput(as.parframe(fits, sortBy = "iterations")))
  expect_snapshot(.canonDput(as.parframe(parlist(list(error = "x")))))
  expect_snapshot(.canonDput(as.parframe(parlist())))
  expect_snapshot(.canonDput(as.parframe(parlist(list(value = 1, argument = c(a = 1), iterations = 2)))))
})

test_that("plot of a parlist, initial against converged values and the paths", {
  fits <- .fitsFixture()
  p <- plot(fits)
  expect_s3_class(p, "ggplot")
  expect_snapshot(.canonDput(as.data.frame(p$data)))
  b <- ggplot2::ggplot_build(p)$data[[1]]
  expect_equal(b$x, c(0, 1, 0, 1))
  expect_equal(b$y, c(1, 2, 2, 3))

  q <- plot(fits, path = TRUE)
  expect_snapshot(.canonDput(as.data.frame(q$data)))

  expect_error(plot(parlist(fits[[4]])), "No converged fit")
  expect_error(plot(parlist(fits[[1]][-10]), path = TRUE), "No path information")
})

test_that("as.parvec of a parframe warns for an unconverged fit", {
  pf <- parframe(data.frame(value = c(2, 1), converged = c(TRUE, FALSE), a = 1:2),
                 parameters = "a", metanames = c("value", "converged"))
  expect_warning(v <- as.parvec(pf), "unconverged")
  expect_equal(unname(v), 2)
  expect_no_warning(as.parvec(pf, 2))
})

test_that("plotPars and plotValues of a parframe", {
  pf <- as.parframe(.fitsFixture())
  pp <- plotPars(pf)
  expect_snapshot(.canonDput(attr(pp, "data")))
  expect_no_error(suppressWarnings(ggplot2::ggplot_build(pp)))
  pv <- plotValues(pf)
  d <- attr(pv, "data")
  expect_equal(d$value, sort(pf$value))
  expect_equal(d$delta, sort(pf$value) - min(pf$value))
  expect_equal(d$index, seq_len(nrow(pf)))
  expect_identical(attr(pv, "jumps"), dMod2:::.stepDetect(sort(pf$value), 1))
  expect_no_error(ggplot2::ggplot_build(pv))
})

test_that("the parlist overlay of a profile plot", {
  sets <- list(`1` = parframe(data.frame(value = c(2, 1, 3), constraint = c(-1, 0, 1), a = 1:3),
                              parameters = "a", metanames = c("value", "constraint")))
  pl <- parframe(data.frame(value = c(1.5, 4), a = c(1, 2), b = c(3, 4)),
                 parameters = c("a", "b"), metanames = "value")
  layer <- dMod2:::.profilePoints(pl, sets)
  expect_snapshot(.canonDput(layer$data))
  m <- cbind(a = c(1, 2), b = c(3, 4))
  expect_snapshot(.canonDput(dMod2:::.profilePoints(m, sets)$data))
})

test_that("subsetting a parvec", {
  d <- matrix(c(1, 0, 0, 0, 2, 0), 3, dimnames = list(c("a", "b", "c"), c("u", "v")))
  x <- as.parvec(c(a = 1, b = 2, c = 3), deriv = d)
  expect_snapshot(.canonDput(unclass(x[c("b", "c")])))
  expect_snapshot(.canonDput(unclass(x["a", drop = TRUE])))
  expect_identical(x[c("a", "b", "c")], x)
})

test_that("parframe() and is.parframe()", {
  expect_snapshot(.canonDput(parframe()))
  df <- data.frame(value = 1:2, a = 3:4, row.names = c("x", "y"))
  expect_snapshot(.canonDput(parframe(df, metanames = "value")))
  expect_snapshot(.canonDput(parframe(cbind(value = 1, a = 2, b = 3), parameters = "a",
                                metanames = "value", objAttributes = "b")))
  expect_snapshot(.canonDput(parframe(data.table::data.table(value = 1, a = 2))))
  expect_true(is.parframe(parframe()))
  expect_false(is.parframe(data.frame()))
})

test_that("as.datalist of a list", {
  l <- list(C1 = data.frame(name = "A", time = 0, value = 1),
            C2 = data.frame(name = "A", time = 1, value = 2, lloq = NA))
  expect_snapshot(.canonDput(as.datalist(l)))
  expect_snapshot(.canonDput(as.datalist(unname(l), names = c("u", "v"))))
  g <- data.frame(dose = 1:2, row.names = c("C1", "C2"))
  expect_identical(attr(as.datalist(l, condition.grid = g), "condition.grid"), g)
  expect_snapshot(.canonDput(as.datalist(list())))
  expect_error(as.datalist(list(C1 = 1)), "list of data.frame expected")
  expect_error(as.datalist(l, names = "u"), "wrong length")
  expect_error(as.datalist(unname(l)), "wrong length")
  expect_error(as.datalist(list(C1 = data.frame(name = "A"))), "mandatory")
  expect_true(is.datalist(as.datalist(l)))
  expect_false(is.datalist(l))
})

test_that("renaming a datalist renames the grid", {
  dl <- .xdGrid()
  names(dl) <- c("lo", "hi")
  expect_snapshot(.canonDput(dl))
  dl <- datalist(C1 = data.frame(name = "A", time = 0, value = 1))
  names(dl) <- "X"
  expect_snapshot(.canonDput(attr(dl, "condition.grid")))
})

test_that("c() of datalists combines the grids", {
  d <- .xdGrid()
  e <- datalist(C3 = data.frame(name = "A", time = 0, value = 9, sigma = 1))
  expect_snapshot(.canonDput(c(d, e)))
  expect_snapshot(.canonDput(c(d[1], d[2])))
  expect_error(c(d, d), "conditions are duplicated: 1, 2")
})

test_that("plotData of a datalist", {
  dl <- .xdGrid()
  p <- plotData(dl)
  expect_identical(attr(p, "data")$bloq, c("no", "no", "yes", "no", "no", "no", "no", "no"))
  g <- plotData(dl, dose == 1, facet = "grid", scales = "fixed")
  expect_identical(as.character(unique(attr(g, "data")$condition)), "1")
  expect_identical(attr(plot(dl), "data"), attr(p, "data"))

  tr <- plotData(dl, transform = list(A = "log(x)"))
  expect_equal(attr(tr, "data")$value[attr(tr, "data")$name == "A"], log(c(1, 3, 5, 7)))
})

test_that("covariates of a datalist are the condition grid", {
  dl <- .xdGrid()
  expect_identical(covariates(dl), attr(dl, "condition.grid"))
  expect_null(covariates(structure(list(), class = c("datalist", "list"))))
})

test_that("c() of a parvec names unnamed operands by their argument names", {
  pv <- as.parvec(c(a = 1, b = 2))
  out <- c(pv, z = 3, w = c(4, 5))
  expect_s3_class(out, "parvec")
  expect_identical(names(out), c("a", "b", "z", "w1", "w2"))
  expect_equal(unclass(out)[["z"]], 3)
  expect_identical(names(c(pv, z = 3)), c("a", "b", "z"))
  expect_identical(names(c(pv, d = c(d = 7))), c("a", "b", "d"))
})
