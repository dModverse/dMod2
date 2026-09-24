# symmetryDetection() and symmetryReduction() on models with states inside
# exponentials.

symdet <- function(...) symmetryDetection(..., verbose = FALSE)

fAB <- function(pre = "kE*")
  eqnvec(A = "kA - kd*A", B = paste0("kB - ", pre, "exp(A)*B"))
gAB <- eqnvec(y = "s*B")
general <- function(res) Filter(function(d) identical(d$type, "general"), res$symmetries)


test_that("an exponential of a state blocks its scaling", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  sc <- symdet(fAB(""), gAB, method = "scaling")
  expect_length(sc$symmetries, 1L)
  expect_setequal(names(sc$symmetries[[1]]$generator), c("B", "kB", "s"))
  # an invariant exponent lets A scale with a
  sc2 <- symdet(eqnvec(A = "kA - kd*A", B = "kB - kE*exp(-A/a)*B"), gAB,
                method = "scaling")
  expect_true(any(vapply(sc2$symmetries, function(d)
    setequal(names(d$generator), c("A", "a", "kA")), logical(1))))
})


test_that("the modular engine carries a state inside exp()", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  obs <- symdet(fAB(""), gAB)
  expect_equal(c(obs$rank, obs$dim), c(5L, 6L))
  expect_setequal(obs$info$coordinates, c("A", "B", "kA", "kB", "kd", "s"))

  tr <- symdet(fAB(), gAB, reconstruct = TRUE)
  expect_equal(c(tr$rank, tr$dim), c(5L, 7L))
  d <- general(tr)
  expect_length(d, 1L)
  gen <- d[[1]]$generator
  expect_setequal(names(gen), c("A", "kA", "kE"))
  expect_true(.symExprEqual(gen[["kA"]], paste0("kd*(", gen[["A"]], ")")))
  expect_true(.symExprEqual(gen[["kE"]], paste0("-kE*(", gen[["A"]], ")")))

  # a known initial value removes the translation
  expect_length(general(symdet(fAB(), gAB, trafo = eqnvec(A = "0"))), 0L)
  # an observable inside exp() carries the same information as its argument
  expect_true(symdet(eqnvec(A = "-k*A"), eqnvec(y = "s*exp(A)"))$identifiable)
})


test_that("other bases, constant offsets and hyperbolic functions", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  g10 <- general(symdet(eqnvec(A = "kA - kd*A", B = "kB - kE*exp10(A)*B"), gAB,
                        reconstruct = TRUE))[[1]]$generator
  expect_true(.symExprEqual(g10[["kE"]], paste0("-log(10)*kE*(", g10[["A"]], ")")))

  mixed <- symdet(eqnvec(A = "kA - kd*A",
                         B = "kB - kE*exp((10 - A)/10)*B + 2^(A/3)*c"), gAB,
                  reconstruct = TRUE)
  expect_equal(c(mixed$rank, mixed$dim), c(6L, 8L))
  gm <- general(mixed)[[1]]$generator
  expect_true(.symExprEqual(gm[["c"]], paste0("-log(2)/3*c*(", gm[["A"]], ")")))
  expect_true(.symExprEqual(gm[["kE"]], paste0("kE/10*(", gm[["A"]], ")")))

  # a shift of u, Vh and a, and a scaling of u, Vh, k and a
  th <- symdet(eqnvec(x = "tanh((u - Vh)/k) - x", u = "a - u"), eqnvec(y = "x"))
  expect_equal(c(th$rank, th$dim), c(3L, 5L))
})


test_that("later events on a state inside exp() move the exponential with it", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  ev0 <- addEvent(eventlist(), var = "B", time = 0, value = "0", method = "add")
  add <- symdet(fAB(), gAB, reconstruct = TRUE,
                events = addEvent(ev0, var = "A", time = 5, value = "dA", method = "add"))
  expect_equal(c(add$rank, add$dim), c(6L, 8L))
  expect_false("dA" %in% unlist(lapply(add$symmetries, `[[`, "support")))

  rep <- symdet(fAB(), gAB, reconstruct = TRUE,
                events = addEvent(ev0, var = "A", time = 5, value = "a1",
                                  method = "replace"))
  expect_true("a1" %in% general(rep)[[1]]$support)

  expect_error(symdet(fAB(), gAB, events = addEvent(ev0, var = "A", time = 5,
                                                    value = "2", method = "multiply")),
               "not supported")
  expect_error(symdet(fAB(), gAB, equilibrate = TRUE), "equilibrate")
})


test_that("Hodgkin-Huxley: conductances relative to C, and only I + gL*EL", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  gate <- function(a, b, x) sprintf("(%s)*(1 - %s) - (%s)*%s", a, x, b, x)
  hh <- eqnvec(
    V = "(I - gNa*m^3*h*(V - ENa) - gK*n^4*(V - EK) - gL*(V - EL))/C",
    m = gate("0.1*(V + 40)/(1 - exp(-(V + 40)/10))", "4*exp(-(V + 65)/18)", "m"),
    h = gate("0.07*exp(-(V + 65)/20)", "1/(1 + exp(-(V + 35)/10))", "h"),
    n = gate("0.01*(V + 55)/(1 - exp(-(V + 55)/10))", "0.125*exp(-(V + 65)/80)", "n"))
  obs <- symdet(hh, eqnvec(y = "V"), reconstruct = TRUE)
  expect_equal(c(obs$rank, obs$dim), c(10L, 12L))
  sc <- Filter(function(d) identical(d$type, "scaling"), obs$symmetries)
  expect_setequal(names(sc[[1]]$generator), c("C", "I", "gK", "gL", "gNa"))
  gl <- general(obs)[[1]]$generator
  expect_setequal(names(gl), c("EL", "I"))
  expect_true(.symExprEqual(gl[["I"]], paste0("-gL*(", gl[["EL"]], ")")))
  expect_setequal(names(symdet(hh, eqnvec(y = "V"), method = "scaling")$symmetries[[1]]$generator),
                  c("C", "I", "gK", "gL", "gNa"))
})


test_that("the polynomial engine treats exponentials as independent variables", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  pol <- symdet(fAB(), gAB, method = "polynomial",
                polynomial = polynomialControl(ansatz = "par", pMax = 1L))
  expect_true(all(vapply(pol$symmetries, function(d) isTRUE(d$verified), logical(1))))
  expect_true(any(vapply(pol$symmetries, function(d)
    setequal(names(d$generator), c("A", "kA", "kE")), logical(1))))
})


test_that("a translation gauge reduces a direction with real coordinates", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  obs <- symdet(fAB(), gAB, reconstruct = TRUE)
  # with every coordinate positive the orbit leaves the domain: no chart
  expect_length(symmetryReduction(obs)$remaining, 1L)
  red <- symmetryReduction(obs, positive = c("B", "kB", "kd", "kE", "s"))
  expect_length(red$remaining, 0L)
  expect_true(symdet(fAB(), gAB, trafo = red$trafo)$identifiable)
})
