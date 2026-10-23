# Characterisation of Xd(), Id(), the event and start checks of the deSolve
# prediction functions, the forcing check of Y(), objective sums and their
# printing, res() with a factor name column, and define()/insert().

.xdGrid <- data.frame(name = c("A", "A", "A", "B", "B"),
                       time = c(0, 1, 3, 2, 0),
                       value = c(1, 2, 4, 5, 3),
                       row.names = c("a0", "a1", "a3", "b2", "b0"))

test_that("Xd interpolates per state and clamps outside its grid", {
  x <- Xd(.xdGrid, condition = "C1")
  expect_identical(getParameters(x), c("a0", "a1", "a3", "b2", "b0"))
  pouter <- attr(attr(x, "mappings")[[1]], "pouter")
  expect_identical(pouter, c(a0 = 1, a1 = 2, a3 = 4, b2 = 5, b0 = 3))
  times <- c(-1, 0, 0.5, 1, 2, 3, 4)
  pr <- x(times, pouter)[["C1"]]
  expect_identical(colnames(pr), c("time", "A", "B"))
  expect_identical(unname(pr[, "A"]), c(1, 1, 1.5, 2, 3, 4, 4))
  expect_identical(unname(pr[, "B"]), c(3, 3, 3.5, 4, 5, 5, 5))
  d <- attr(pr, "deriv")
  expect_identical(dim(d), c(7L, 2L, 5L))
  expect_identical(dimnames(d)[[3]], names(pouter))
  expect_identical(unname(d[, "A", "a1"]), c(0, 0, 0.5, 1, 0.5, 0, 0))
  expect_identical(unname(d[, "A", "a3"]), c(0, 0, 0, 0, 0.5, 1, 1))
  expect_identical(unname(d[, "B", "b0"]), c(1, 1, 0.75, 0.5, 0, 0, 0))
  expect_identical(unname(d[, "B", "a0"]), rep(0, 7))
  expect_null(attr(x(times, pouter, deriv = FALSE)[["C1"]], "deriv"))
  expect_error(x(times, pouter, deriv2 = TRUE), "second-order")
})

test_that("Xd without value column has no pouter, one point is constant", {
  x <- Xd(data.frame(name = "A", time = 2, row.names = "p"))
  expect_null(attr(attr(x, "mappings")[[1]], "pouter"))
  pr <- x(c(0, 2, 5), c(p = 7))[[1]]
  expect_identical(unname(pr[, "A"]), c(7, 7, 7))
  expect_identical(unname(attr(pr, "deriv")[, "A", "p"]), c(1, 1, 1))
})

test_that("Xd vjp equals the forward Jacobian contracted with the cotangent", {
  x <- Xd(.xdGrid)
  p <- c(a0 = 1, a1 = 2, a3 = 4, b2 = 5)
  fixed <- c(b0 = 3)
  times <- c(0, 0.5, 2.5, 4)
  pr <- x(times, p, fixed = fixed)[[1]]
  set.seed(2)
  w <- matrix(rnorm(8), 4, 2, dimnames = list(NULL, c("A", "B")))
  ref <- apply(attr(pr, "deriv") * as.vector(w), 3, sum)
  vjp <- attr(attr(x, "mappings")[[1]], "vjpfn")
  got <- vjp(times, p, fixed, w)
  expect_identical(rownames(got), names(p))
  expect_equal(got[, 1], ref[names(p)], tolerance = 1e-14)
  wA <- w[, "A", drop = FALSE]
  gotA <- vjp(times, p, fixed, wA)
  expect_equal(gotA[c("b2"), 1], c(b2 = 0))
})

test_that("Id is neutral in a product", {
  id <- Id()
  expect_s3_class(id, c("idfn", "fn"), exact = TRUE)
  expect_null(id())
  x <- Xt()
  expect_identical(id * x, x)
  expect_identical(x * id, x)
})

test_that("Xs on deSolve reports events given outside the odemodel", {
  withr::local_dir(tempdir())
  f <- eqnlist() |> addReaction("A", "", "k*A")
  ev <- eventlist(var = "A", time = 1, value = 1, method = "add")
  m <- odemodel(f, modelname = "fnb_ev0", backend = "deSolve", compile = FALSE)
  expect_message(Xs(m, events = ev), "Events should be defined in odemodel")
  expect_silent(Xs(m))
  m2 <- odemodel(f, events = ev, modelname = "fnb_ev1", backend = "deSolve",
                 compile = FALSE)
  expect_warning(Xs(m2, events = ev), "Events already defined in odemodel")
  expect_silent(Xs(m2))
})

test_that("Xf on deSolve starts states missing from pars at zero", {
  withr::local_dir(tempdir())
  f <- eqnlist() |> addReaction("A", "B", "k*A") |> addReaction("", "C", "1")
  m <- odemodel(f, modelname = "fnb_xf", backend = "deSolve", deriv = FALSE,
                compile = FALSE)
  x <- Xf(m)
  compile(x)
  expect_identical(getParameters(x), c("A", "B", "C", "k"))
  pr <- x(c(0, 1, 2), c(A = 1, k = 0.5))[[1]]
  expect_identical(unname(pr[1, c("A", "B", "C")]), c(1, 0, 0))
  expect_equal(unname(pr[, "A"] + pr[, "B"]), c(1, 1, 1), tolerance = 1e-5)
  expect_equal(unname(pr[, "C"]), c(0, 1, 2), tolerance = 1e-5)
})

test_that("Y infers states, forcings and parameters from a prediction function", {
  withr::local_dir(tempdir())
  f <- eqnlist() |> addReaction("A", "", "k*A*F")
  m <- odemodel(f, forcings = "F", modelname = "fnb_yf", compile = FALSE)
  x <- Xs(m, forcings = data.frame(name = "F", time = 0:1, value = 1))
  g <- Y(c(y = "s*A + F + off", z = "time*k"), f = x, compile = FALSE,
         modelname = "fnb_obs")
  expect_setequal(getParameters(g), c("s", "off", "k"))
  expect_false("F" %in% getParameters(g))
  expect_false("time" %in% getParameters(g))
  xbad <- Xs(m, forcings = data.frame(name = "A", time = 0:1, value = 1))
  expect_error(Y(c(y = "A"), f = xbad, compile = FALSE, modelname = "fnb_bad"),
               "Forcings and states overlap in different conditions.", fixed = TRUE)
})

# Objective stub that records the conditions and environment it was called with.
.logObj <- function(conds, value, log) {
  f <- function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE, hessian = NULL,
                conditions = NULL, env = NULL, cores = 1L, sweep = "forward") {
    log$calls <- c(log$calls, list(list(id = value, conditions = conditions, env = env)))
    out <- objlist(value, c(a = value), NULL)
    attr(out, "env") <- if (is.null(env)) log$env else env
    out
  }
  class(f) <- c("objfn", "fn")
  attr(f, "conditions") <- conds
  attr(f, "parameters") <- "a"
  f
}

test_that("+.objfn evaluates every term on its own conditions", {
  log <- new.env(); log$env <- new.env()
  o <- .logObj(c("C1", "C2"), 1, log) + .logObj(NULL, 10, log) + .logObj("C3", 100, log)
  expect_identical(attr(o, "conditions"), c("C1", "C2", "C3"))
  r <- o(c(a = 0))
  expect_identical(r$value, 111)
  expect_identical(lapply(log$calls, `[[`, "conditions"),
                   list(c("C1", "C2"), NULL, "C3"))
  # A term without conditions gets the caller's env, a later conditioned one
  # the env the terms before it returned.
  expect_null(log$calls[[2]]$env)
  expect_identical(log$calls[[3]]$env, log$env)
  expect_identical(attr(r, "env"), log$env)

  log$calls <- NULL
  r <- o(c(a = 0), conditions = c("C2", "C9"))
  expect_identical(r$value, 11)
  expect_identical(lapply(log$calls, `[[`, "conditions"), list("C2", NULL))
  expect_null(log$calls[[1]]$env)

  log$calls <- NULL
  # The unconditioned term inside the sum over C1, C2 is evaluated as in a flat sum.
  r <- o(c(a = 0), conditions = "C9")
  expect_identical(r$value, 10)
  expect_identical(lapply(log$calls, `[[`, "conditions"), list(NULL))

  e0 <- new.env()
  log$calls <- NULL
  o(c(a = 0), env = e0)
  expect_true(all(vapply(log$calls, function(cl) identical(cl$env, e0), TRUE)))
})

test_that("+.objlist sums on the union of parameters and keeps attributes", {
  a <- objlist(1, c(a = 1, b = 2), matrix(1, 2, 2, dimnames = list(c("a", "b"), c("a", "b"))))
  b <- objlist(2, c(c = 3, a = 1), NULL)
  attr(a, "data") <- 0.5
  s <- a + b
  expect_identical(names(s), c("value", "gradient", "hessian"))
  expect_identical(s$gradient, c(a = 2, b = 2, c = 3))
  expect_identical(attr(s, "data"), 0.5)
  v <- unclass(objlist(1, NULL, NULL))[1]
  class(v) <- "objlist"
  s2 <- v + v
  expect_identical(names(s2), "value")
  expect_identical(s2$value, 2)
})

test_that("print.objlist prints value, truncated gradient and Hessian, attributes", {
  g <- setNames(seq_len(25) / 10, paste0("p", seq_len(25)))
  H <- diag(25); dimnames(H) <- list(names(g), names(g))
  o <- objlist(3.5, g, H)
  attr(o, "chi2") <- c(data = 2.25)
  attr(o, "data") <- 2.25
  expect_snapshot(print(o))
  expect_snapshot(print(o, n1 = 3, n2 = 2))
  expect_snapshot(print(objlist(1, c(a = 1), NULL)))
  expect_snapshot(print(structure(list(value = 4), class = "objlist")))
  expect_invisible(print(o))
})

test_that("res accepts a factor name column and returns character names", {
  prdf <- matrix(c(0, 1, 2, 3, 5, 7), 2, 3, dimnames = list(NULL, c("time", "A", "B")))
  dat <- data.frame(time = c(1, 0), name = factor(c("B", "A")), value = c(6, 2),
                    sigma = c(1, 2), lloq = -Inf)
  r <- res(dat, prdf)
  expect_identical(r$name, c("B", "A"))
  expect_identical(r$prediction, c(7, 2))
  expect_identical(r$weighted.residual, c(1, 0))
})

test_that("define and insert respect conditionMatch and the condition table", {
  tr <- define(NULL, "x ~ x", x = c("A", "k"))
  tab <- data.frame(dose = c(1, 10), row.names = c("low", "high"))
  trL <- branch(tr, table = tab)
  d <- define(trL, "A ~ dose", dose = dose, conditionMatch = "^h")
  expect_identical(names(d), c("low", "high"))
  expect_identical(as.character(d$low["A"]), "A")
  expect_identical(as.character(d$high["A"]), "10")
  expect_identical(attr(d, "tree"), tab)
  i <- insert(trL, "x ~ exp(x)", x = "k", conditionMatch = "low")
  expect_identical(as.character(i$low["k"]), "exp(k)")
  expect_identical(as.character(i$high["k"]), "k")
  g <- insert(trL, "x ~ 2*x", x = "k", dose > 5)
  expect_identical(as.character(g$low["k"]), "k")
  expect_identical(as.character(g$high["k"]), "2 * k")
  s <- insert(tr, "x ~ exp(x)", x = .currentSymbols)
  expect_identical(as.character(s), c("exp(A)", "exp(k)"))
  expect_identical(as.character(define(tr, "B ~ 0")), c("A", "k", "0"))
  expect_error(define(list(tr), "A ~ 1"), "must be named")
  expect_error(insert(list(foo = tr), "A ~ 1"), "subset")
  expect_error(define(tr, "A ~ 1", conditionMatch = "x"), "branch")
  f <- function(v) define(tr, "A ~ v", v = v)
  expect_identical(as.character(f(3)["A"]), "3")
})
