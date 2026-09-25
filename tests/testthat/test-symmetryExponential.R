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


test_that("the exp stage reaches a wide block at a lower numerator degree", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  skip_on_cran()

  # thermal abuse of a Li-ion cell (Hatchard et al. 2001), T measured: eight
  # directions over 25 coordinates, one of them a translation of tsei
  Rsei <- "Asei*exp(-Esei/T)*csei"
  Rne  <- "Ane*exp(-tsei/tref)*exp(-Ene/T)*cneg"
  Rpe  <- "Ape*exp(-Epe/T)*alpha*(1 - alpha)"
  Re   <- "Ae*exp(-Ee/T)*ce"
  f <- eqnvec(
    csei = paste0("-", Rsei), cneg = paste0("-", Rne), tsei = Rne, alpha = Rpe,
    ce = paste0("-", Re),
    T = sprintf("(Hsei*Wc*%s + Hne*Wc*%s + Hpe*Wp*%s + He*We*%s - hA*(T - Ta))/rhocp",
                Rsei, Rne, Rpe, Re))
  obs <- symdet(f, eqnvec(y = "T"), reconstruct = TRUE)
  expect_length(obs$symmetries, 8L)
  red <- suppressWarnings(symmetryReduction(obs))
  expect_length(red$remaining, 0L)
  inv <- unlist(lapply(red$blocks, `[[`, "invariants"))
  expect_true(any(grepl("exp(", inv, fixed = TRUE)))
  expect_true(symdet(f, eqnvec(y = "T"), trafo = red$trafo)$identifiable)
})


test_that("time in f is a clock, not a parameter", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  f <- eqnvec(x = "A*exp(-lam*time) - d*x")
  obs <- symdet(f, eqnvec(y = "s*x"), reconstruct = TRUE)
  expect_false("time" %in% obs$info$coordinates)
  expect_length(obs$symmetries, 1L)
  expect_setequal(names(obs$symmetries[[1]]$generator), c("A", "s", "x"))
  sc <- symdet(f, eqnvec(y = "s*x"), method = "scaling")
  expect_length(sc$symmetries, 1L)
  expect_setequal(names(sc$symmetries[[1]]$generator), c("A", "s", "x"))
})


test_that("log() and fractional powers of positive coordinates", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # Gompertz growth: x and K scale together, the scale s inversely
  gz <- eqnvec(x = "a*x*log(K/x)")
  gg <- eqnvec(y = "s*x")
  for (m in c("observability", "scaling")) {
    r <- symdet(gz, gg, method = m, reconstruct = TRUE)
    expect_length(r$symmetries, 1L)
    expect_identical(r$symmetries[[1]]$type, "scaling")
    w <- vapply(r$symmetries[[1]]$weights, as.numeric, numeric(1))
    expect_equal(w[c("K", "s", "x")] / w[["K"]], c(K = 1, s = -1, x = 1),
                 ignore_attr = TRUE)
  }
  obs <- symdet(gz, gg, reconstruct = TRUE)
  expect_setequal(obs$info$coordinates, c("x", "K", "a", "s"))
  red <- suppressWarnings(symmetryReduction(obs))
  expect_length(red$remaining, 0L)
  expect_true(symdet(gz, gg, trafo = red$trafo)$identifiable)
  # log() needs the positive domain
  expect_error(symdet(gz, gg, positive = FALSE), "positive")

  # sqrt: x scales by lambda^2, k by lambda, s by lambda^-2
  sq <- eqnvec(x = "-k*sqrt(x)")
  for (m in c("observability", "scaling")) {
    r <- symdet(sq, gg, method = m, reconstruct = TRUE)
    expect_length(r$symmetries, 1L)
    w <- vapply(r$symmetries[[1]]$weights, as.numeric, numeric(1))
    expect_equal(w[c("k", "s", "x")] / w[["k"]], c(k = 1, s = -2, x = 2),
                 ignore_attr = TRUE)
  }
})


test_that("events, abs() and logarithms of sums in the log chart", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  gz <- eqnvec(x = "a*x*log(K/x)")
  gg <- eqnvec(y = "s*x")
  # a known dose separates s and K; log(2) enters exactly, as a known constant
  dose <- addEvent(eventlist(), var = "x", time = 0, value = "2", method = "replace")
  expect_true(symdet(gz, gg, events = dose)$identifiable)
  # a dose in a parameter joins the chart
  dD <- addEvent(eventlist(), var = "x", time = 0, value = "D", method = "replace")
  r <- symdet(gz, gg, events = dD, reconstruct = TRUE)
  expect_setequal(names(r$symmetries[[1]]$generator), c("D", "K", "s"))
  expect_error(symdet(gz, gg, events = addEvent(eventlist(), var = "x", time = 3,
                                                value = "1", method = "add")), "add")

  # abs() of a positive state, max() refused
  r <- symdet(eqnvec(x = "-k*abs(x)*x"), gg, reconstruct = TRUE)
  expect_setequal(names(r$symmetries[[1]]$generator), c("k", "s", "x"))
  expect_error(symdet(eqnvec(x = "-k*max(x, c)"), gg), "not analytic")

  # log(c + x): c and x trade places, mapped back by the chain rule
  fc <- eqnvec(x = "-k*log(c + x)")
  gc <- eqnvec(y = "s*(c + x)")
  r <- symdet(fc, gc, reconstruct = TRUE)
  expect_length(r$symmetries, 1L)
  gen <- r$symmetries[[1]]$generator
  expect_setequal(names(gen), c("c", "x"))
  expect_true(.symExprEqual(gen[["x"]], paste0("-(", gen[["c"]], ")")))
  red <- suppressWarnings(symmetryReduction(r))
  expect_true(symdet(fc, gc, trafo = red$trafo)$identifiable)
  # two different logarithms of x
  expect_error(symdet(eqnvec(x = "r*x*log(1 + K/x)"), gg), "two different")
})
