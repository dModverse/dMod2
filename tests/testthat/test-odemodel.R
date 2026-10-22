test_that("deSolve forcings follow linear or constant interpolation in Xs and Xf", {
  withr::local_dir(tempdir())
  f <- eqnlist() |>
    addReaction("", "A", "scale*F") |>
    addReaction("A", "", "k*A")
  model <- odemodel(f, forcings = "F", fixed = "scale", backend = "deSolve", compile = FALSE)
  compile(Xs(model))

  knots <- 0:5
  inputs <- list(wave = data.frame(name = "F", time = knots, value = cos(knots)),
                 decay = data.frame(name = "F", time = knots, value = 2^-knots),
                 flat = data.frame(name = "F", time = 0, value = 0.5))
  build <- function(ctor, ...) Reduce(`+`, Map(function(fc, cn) ctor(model, fc, condition = cn, ...),
                                               inputs, names(inputs)))
  grid <- seq(0, 5, by = 0.002)
  forcing <- function(x) {
    p <- getParameters(x)
    x(grid, setNames(seq(0.4, 1, length.out = length(p)), p))$wave[, "F"]
  }

  # Linear: the slope changes at each of the four inner knots.
  kinks <- abs(diff(diff(forcing(build(Xs))))) > 1e-8
  expect_identical(sum(kinks), 4L)

  # Constant: the forcing jumps at each of the four inner knots, with or without sensitivities.
  piecewise <- list(method = "constant", rule = 2)
  expect_identical(sum(abs(diff(forcing(build(Xs, fcontrol = piecewise)))) > 1e-8), 4L)
  expect_identical(sum(abs(diff(forcing(build(Xf, fcontrol = piecewise)))) > 1e-8), 4L)
})
test_that("Xs and Xf on deSolve check their arguments and return the requested times", {
  withr::local_dir(tempdir())
  f <- eqnlist() |> addReaction("A", "B", "k*A") |> addReaction("B", "", "d*B")
  ev <- eventlist(var = "A", time = 1.5, value = 1, method = "add")
  model <- odemodel(f, events = ev, modelname = "desolve_args", backend = "deSolve",
                    compile = FALSE)
  x <- Xs(model)
  compile(x)
  expect_error(Xs(model, names = "Z"), "Z")
  expect_error(Xs(model, optionsReverse = list(refine = TRUE)), "optionsReverse")
  expect_error(Xf(model, optionsSens = list()), "optionsSens")
  pars <- c(A = 1, B = 0, k = 1, d = 0.5)
  pr <- x(0:2, pars)[[1]]
  expect_equal(unname(pr[, "time"]), 0:2)
  expect_equal(dim(attr(pr, "deriv"))[1], 3L)
  expect_error(x(0:2, pars[-4]), "d")
  x3 <- Xs(model, options = list(maxsteps = 3))
  expect_warning(x3(0:2, pars, deriv = FALSE), "stopped")
  expect_warning(x3(0:2, pars), "stopped")
  m0 <- odemodel(f, modelname = "desolve_noderiv", backend = "deSolve", deriv = FALSE,
                 compile = FALSE)
  expect_error(Xs(m0), "sensitivities")
})

test_that("odemodel(deriv = FALSE) builds no derivative directions", {
  m <- odemodel(eqnvec(A = "-k*A"), deriv = FALSE, derivMode = c("forward", "reverse"),
                modelname = "test_om_noderiv_r", compile = FALSE)
  expect_null(m$reversed)
  m2 <- odemodel(eqnvec(A = "-k*A"), deriv = FALSE,
                 derivMode = c("forward", "forward-forward"),
                 modelname = "test_om_noderiv_ff", compile = FALSE)
  expect_null(m2$extended)
  expect_null(m2$extended2)
})
