# normL2(times = ) as a list named by condition: each condition solves on its
# own grid from its first time, and an event before that start does not fire.

.times_fx <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    .dmod_with_fx_workdir({
      f  <- addReaction(eqnlist(), from = "A", to = "", rate = "k*A")
      ev <- addEvent(eventlist(), var = "A", time = -10, value = "dose", method = "add")
      x  <- Xs(odemodel(f, events = ev, modelname = "t0pc_x", compile = FALSE))
      g  <- Y(c(y = "A"), f = x, attach.input = FALSE, modelname = "t0pc_g", compile = FALSE)
      p1 <- P(eqnvec(A = "A0", k = "k", dose = "dose"), condition = "C1",
              modelname = "t0pc_p1", compile = FALSE)
      p2 <- P(eqnvec(A = "A0", k = "k", dose = "0"), condition = "C2",
              modelname = "t0pc_p2", compile = FALSE)
      compile(x, g, p1, p2, output = "t0pc_all", cores = 2L)
      cache <<- list(prd = g * x * (p1 + p2))
    })
    cache
  }
})

.times_data <- function() {
  d <- function(cn, v) data.frame(name = "y", time = c(0, 1, 2), value = v,
                                  sigma = 0.1, condition = cn)
  as.datalist(rbind(d("C1", c(0.6, 0.4, 0.2)), d("C2", c(0.9, 0.5, 0.3))))
}

test_that("per-condition times start each condition at its own first time", {
  skip_if_not_installed("cppDE"); skip_on_cran()
  prd  <- .times_fx()$prd
  data <- .times_data()
  pars <- c(A0 = 1, k = 0.05, dose = 2)

  early <- normL2(data, prd, times = -10)(pars, deriv = TRUE)
  split <- normL2(data, prd, times = list(C1 = -10, C2 = 0))(pars, deriv = TRUE)
  late  <- normL2(data, prd)(pars, deriv = TRUE)

  # C2 decays from -10 to 0 under `early`, from 0 under `split`
  pr <- attr(split, "env")$prediction
  expect_equal(unname(pr$C1[1, "time"]), -10)
  expect_equal(unname(pr$C2[1, "time"]), 0)
  expect_equal(unname(pr$C2[pr$C2[, "time"] == 0, "y"]), 1, tolerance = 1e-6)
  expect_equal(unname(pr$C1[pr$C1[, "time"] == 0, "y"]), 3 * exp(-0.5), tolerance = 1e-6)
  # the event before a start of 0 does not fire
  prl <- attr(late, "env")$prediction
  expect_equal(unname(prl$C1[prl$C1[, "time"] == 0, "y"]), 1, tolerance = 1e-6)

  expect_false(isTRUE(all.equal(split$value, early$value)))
  expect_false(isTRUE(all.equal(split$value, late$value)))
  # a condition the list leaves out solves on its data times alone
  expect_equal(normL2(data, prd, times = list(C1 = -10))(pars, deriv = TRUE)$value,
               split$value)

  # gradient against central differences
  fd <- vapply(names(pars), function(nm) {
    h <- 1e-6; up <- dn <- pars; up[nm] <- up[nm] + h; dn[nm] <- dn[nm] - h
    obj <- normL2(data, prd, times = list(C1 = -10, C2 = 0))
    (obj(up, deriv = FALSE)$value - obj(dn, deriv = FALSE)$value) / (2 * h)
  }, numeric(1))
  expect_equal(unname(split$gradient[names(pars)]), unname(fd), tolerance = 1e-4)
})

test_that("times as a list is named by the conditions of the data", {
  skip_if_not_installed("cppDE"); skip_on_cran()
  prd  <- .times_fx()$prd
  data <- .times_data()
  expect_error(normL2(data, prd, times = list(C1 = -10, C3 = 0)), "not C3")
  expect_error(normL2(data, prd, times = list(-10, 0)), "named by the conditions")
})
