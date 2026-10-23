# Data, parameter and fit classes, without compiled code.

.dlGrid <- function() {
  as.datalist(data.frame(name = "A", time = rep(0:1, 2), value = 1:4, sigma = 1,
                         dose = rep(c(1, 2), each = 2)), splitBy = "dose")
}

test_that("lloq = NA means no limit of quantification", {
  d <- data.frame(name = "A", time = 0:2, value = c(1, 2, 3), sigma = 1, condition = "C1")
  x <- Xd(data.frame(name = "A", time = 0:2, row.names = c("a0", "a1", "a2")), condition = "C1")
  p <- c(a0 = 0, a1 = 0, a2 = 0)
  ref <- normL2(as.datalist(transform(d, lloq = -Inf)), x)(p)$value
  expect_equal(normL2(as.datalist(transform(d, lloq = NA)), x)(p)$value, ref)
})

test_that("renaming a datalist renames the conditions of its data frame", {
  dl <- as.datalist(data.frame(name = "A", time = 0, value = 1:2, sigma = 1, condition = c("C1", "C2")))
  names(dl) <- c("low", "high")
  expect_equal(unique(as.data.frame(dl)$condition), c("low", "high"))
})

test_that("subset of a datalist keeps the condition grid", {
  dl <- subset(.dlGrid(), time > 0)
  expect_equal(attr(dl, "condition.grid")$dose, c(1, 2))
  expect_equal(unname(vapply(dl, nrow, 1L)), c(1L, 1L))
})

test_that("reduceReplicates keeps the levels of factor columns", {
  d <- data.frame(name = "A", time = 0, value = c(1, 3, 2, 4),
                  condition = factor(c("ctrl", "ctrl", "treat", "treat")))
  r <- suppressMessages(reduceReplicates(d))
  expect_setequal(as.character(r$condition), c("ctrl", "treat"))
})

test_that("as.datalist takes data without covariates, tibbles and data.tables", {
  d <- data.frame(name = "A", time = 0:1, value = 1:2)
  dl <- as.datalist(d)
  expect_equal(names(dl), "1")
  expect_equal(nrow(dl[[1]]), 2L)
  d$condition <- c("C1", "C2")
  expect_equal(names(as.datalist(data.table::as.data.table(d))), c("C1", "C2"))
  expect_equal(names(as.datalist(dplyr::as_tibble(d))), c("C1", "C2"))
})

test_that("covariates of a data frame are a data frame for any number of covariates", {
  d <- data.frame(name = "A", time = 0, value = 1:2, condition = c("C1", "C2"), dose = c(1, 2))
  cv <- covariates(d)
  expect_s3_class(cv, "data.frame")
  expect_equal(rownames(cv), c("C1", "C2"))
  expect_equal(cv$dose, c(1, 2))
  expect_equal(nrow(covariates(d[, 1:4])), 2L)
})

test_that("wide2long refuses duplicated names", {
  expect_error(wide2long(cbind(time = 0:1, A = 1:2, A = 3:4)), "A")
})

test_that("a parframe keeps its class for one column and selects columns by one index", {
  ms <- msParframe(c(a = 0), 3)
  sub <- ms[2:3, ]
  expect_s3_class(sub, "parframe")
  expect_equal(colnames(sub), "a")
  expect_equal(attr(sub, "parameters"), "a")
  pf <- parframe(data.frame(value = 1:3, a = 4:6, b = 7:9), parameters = c("a", "b"),
                 metanames = "value")
  expect_equal(pf["a"]$a, 4:6)
  expect_equal(colnames(pf[2]), "a")
})

test_that("as.parvec of a parframe checks the index", {
  pf <- parframe(data.frame(value = 1:3, converged = TRUE, a = 4:6), parameters = "a",
                 metanames = c("value", "converged"))
  expect_equal(unname(as.parvec(pf, 2)), 5)
  expect_error(as.parvec(pf, 10), "3")
  expect_error(as.parvec(pf, 0))
})

test_that("a parvec with NA prints", {
  expect_output(print(as.parvec(c(a = NA, b = 1))), "NA")
})

test_that("parvec() names unnamed entries by position", {
  expect_equal(names(parvec(a = 1, 2)), c("a", "par2"))
  expect_equal(names(parvec(c(a = 1, b = 2))), c("a", "b"))
})

test_that("summary of a parlist works without fits and without convergence flags", {
  expect_output(summary(parlist()), "total")
  fits <- parlist(list(value = 1, argument = c(a = 1), iterations = 1),
                  list(value = 2, argument = c(a = 2), iterations = 1, converged = TRUE))
  expect_output(summary(fits), "converged")
})

test_that("c() of parlists concatenates the fits", {
  f <- list(value = 1, argument = c(a = 1), converged = TRUE, iterations = 1)
  fits <- c(parlist(f), parlist(f, f))
  expect_s3_class(fits, "parlist")
  expect_length(fits, 3L)
  expect_null(fits[[3]]$index)
})

test_that("empty datalists and the default parameters of a parframe", {
  expect_length(datalist(), 0L)
  expect_equal(names(datalist(data.frame(name = "A", time = 0, value = 1))), "1")
  pf <- parframe(data.frame(value = 1, a = 2), metanames = "value")
  expect_equal(attr(pf, "parameters"), "a")
})

test_that("conditions with different extra columns bind", {
  d <- datalist(C1 = data.frame(name = "A", time = 0, value = 1, sigma = 1, rep = 1),
                C2 = data.frame(name = "A", time = 0, value = 2, sigma = 1))
  out <- lbind(unclass(d))
  expect_equal(nrow(out), 2L)
  expect_true(is.na(out$rep[2]))
})

test_that("c() of datalists keeps the class and refuses a condition twice", {
  d <- datalist(C1 = data.frame(name = "A", time = 0, value = 1, sigma = 1))
  e <- datalist(C2 = data.frame(name = "A", time = 0, value = 2, sigma = 1))
  expect_s3_class(c(d, e), "list")
  expect_error(c(d, d), "duplicated")
})

test_that("as.datalist refuses condition names that collide", {
  df <- data.frame(name = "A", time = 0, value = 1, sigma = 1,
                   a = c("x_y", "x"), b = c("z", "y_z"))
  expect_error(as.datalist(df), "same condition name")
})

test_that("subset() of a parframe takes all expressions in the caller's frame", {
  pf <- parframe(data.frame(value = c(10, 15, 40, NA), converged = c(TRUE, TRUE, TRUE, FALSE),
                            a = 1:4), metanames = c("value", "converged"))
  expect_equal(subset(pf, converged, value < 20)$value, c(10, 15))
  cut <- 20
  f <- function(pf) { cut <- 12; subset(pf, value < cut) }
  expect_equal(f(pf)$value, 10)
  expect_equal(subset(pf, value > 12)$value, c(15, 40))
  g <- function(pf) { lim <- 12; plotValues(pf, value < lim) }
  expect_s3_class(g(pf), "ggplot")
})

test_that("a parlist binds fits with different parameters, and c() takes a single fit", {
  f1 <- list(value = 1, argument = c(a = 1), converged = TRUE, iterations = 1)
  f2 <- list(value = 2, argument = c(a = 2, b = 3), converged = TRUE, iterations = 1)
  pf <- as.parframe(parlist(f1, f2))
  expect_true(is.na(pf$b[1]))
  expect_length(c(parlist(f1, f1), f2), 3L)
})
