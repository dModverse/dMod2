# The wrapper around symident: argument and field names, the model specification,
# the objects and their reports. The analysis itself is tested in symident.

test_that("argument names go to snake_case and field names back to camelCase", {
  expect_identical(dMod2:::.symSnake(c("relevanceCapDir", "laurentDegNum", "timeout",
                                       "perprimeMinPrimes", "stackedMin")),
                   c("relevance_cap_dir", "laurent_deg_num", "timeout",
                     "perprime_min_primes", "stacked_min"))
  expect_identical(dMod2:::.symCamel(c("lie_order_used", "rank", "gauge_suggestion")),
                   c("lieOrderUsed", "rank", "gaugeSuggestion"))
  x <- list(info = list(lie_order_used = 2L, model_exprs = list("k_on")),
            symmetries = list(list(complete_generator = list(k_on = "k_on"),
                                   generator = list(k_off = "-k_off"),
                                   state_components = list(x_c1 = "x_c1"))))
  y <- dMod2:::.symCamelFields(x)
  expect_identical(names(y$info), c("lieOrderUsed", "modelExprs"))
  expect_identical(names(y$symmetries[[1]]), c("completeGenerator", "generator", "stateComponents"))
  expect_identical(names(y$symmetries[[1]]$completeGenerator), "k_on")
  expect_identical(names(y$symmetries[[1]]$generator), "k_off")
  expect_identical(names(y$symmetries[[1]]$stateComponents), "x_c1")
})

test_that("reconstControl() passes homotopy to symident as a flag", {
  expect_true(reconstControl()$homotopy)
  expect_false(reconstControl(homotopy = FALSE)$homotopy)
  expect_error(reconstControl(homotopy = NA))
  expect_identical(dMod2:::.symControlSnake("homotopy"), "homotopy")
})

test_that("the optional reconstruction settings travel only when set", {
  ctrl <- reconstControl()
  for (k in c("homotopyDegree", "degreeCapMax", "exactProof", "exactWorkCap", "seed"))
    expect_true(k %in% names(ctrl) && is.null(ctrl[[k]]))
  ctrl <- reconstControl(degreeCapMax = 6, exactProof = FALSE, exactWorkCap = 1e6, seed = 3)
  expect_identical(ctrl$degreeCapMax, 6L)
  expect_identical(ctrl$seed, 3L)
  expect_false(ctrl$exactProof)
  expect_error(reconstControl(homotopyDegree = 0))
  expect_identical(dMod2:::.symControlSnake(c("degreeCapMax", "exactWorkCap", "homotopyDegree")),
                   c("degree_cap_max", "exact_work_cap", "homotopy_degree"))
})

test_that("a missing msolve is offered, and the call repeated once it is installed", {
  skip_if_not_installed("jsonlite")
  calls <- 0L
  missingMsolve <- structure(class = c("symident.sym.msolve.MsolveMissingError", "error",
                                       "condition"), list(message = "needs msolve", call = NULL))
  local_mocked_bindings(.symident = function(...) list(f = function(...) {
    calls <<- calls + 1L
    if (calls == 1L) stop(missingMsolve)
    '{"ok": true}'
  }), .package = "dMod2")
  local_mocked_bindings(.offerMsolve = function() FALSE, .package = "dMod2")
  expect_error(dMod2:::.symCall("f"), "confirm its installation")
  calls <- 0L
  local_mocked_bindings(.offerMsolve = function() TRUE, .package = "dMod2")
  expect_true(dMod2:::.symCall("f")$ok)
  expect_identical(calls, 2L)
})

test_that("the specification holds reactions, totals, conditions and events", {
  r <- eqnlist() |>
    addReaction("A", "B", "k1 * A") |>
    addReaction("B", "A", "k2 * B") |>
    customTotals(list(totC = "A+B"))
  ev <- addEvent(eventlist(), var = "A", time = 0, value = "dose", method = "replace")
  cd <- data.frame(dose = c(1, 2), row.names = c("lo", "hi"))
  s <- dMod2:::.symSpec(r, eqnvec(y = "s*A"), NULL, cd, ev)
  expect_identical(unlist(s$reactions$species), c("A", "B"))
  expect_identical(names(s$totals), "totC")
  expect_identical(names(s$f), c("A", "B"))
  expect_identical(s$g, list(y = "s*A"))
  expect_identical(unlist(s$conditions$rows), c("lo", "hi"))
  expect_identical(unlist(s$conditions$cols$dose), c(1, 2))
  expect_identical(s$events[[1]][c("var", "time", "value", "method")],
                   list(var = "A", time = "0", value = "dose", method = "replace"))
  g2 <- dMod2:::.symSpec(eqnvec(x = "-k*x"), list(eqnvec(y = "x"), eqnvec(y = "s*x")),
                         NULL, NULL, NULL)$g
  expect_length(g2, 2L)
})

skip_if_not_installed("reticulate")
skip_if_not_installed("jsonlite")
skip_if_not(reticulate::py_module_available("symident"), "symident not available")

f_ab <- eqnvec(A = "-k1*A + k2*B", B = "k1*A - k2*B")

test_that("a scaling and a general direction, with the result object", {
  r <- symmetryDetection(f_ab, eqnvec(y = "s*A"), reconstruct = TRUE, verbose = FALSE)
  expect_s3_class(r, "symmetrydetection")
  expect_false(r$identifiable)
  expect_identical(c(r$rank, r$dim), c(3L, 5L))
  expect_identical(vapply(r$symmetries, `[[`, "", "type"), c("scaling", "general"))
  expect_s3_class(r$symmetries[[1]]$generator, "eqnvec")
  expect_identical(unclass(r$symmetries[[1]]$generator), c(A = "A", B = "B", s = "-s"))
  expect_identical(r$info$lieOrderDriver, 1L)
  expect_identical(vapply(r$symmetries, `[[`, "", "level"), c("exact", "exact"))
  expect_identical(r$info$rankStatus$status, "proven")
  expect_true(all(c("A", "B", "k1", "k2", "s") %in% r$info$coordinates))
  expect_false(is.null(attr(r, "symident")))
  expect_identical(r$call[[1]], quote(symmetryDetection))
  out <- capture.output(print(r))
  expect_match(out[1], "rank 3 / 5")
  expect_true(any(grepl("Scalings:", out)))
  sm <- capture.output(summary(r))
  expect_true(any(grepl("Computation:", sm)))
  expect_true(any(grepl("Lie order", sm)))
})

test_that("an unknown input and its change along each direction", {
  r <- symmetryDetection(eqnvec(L = "-(c*T + L)/tau"), eqnvec(y = "L"), unknownInputs = "T",
                         reconstruct = TRUE, verbose = FALSE)
  expect_identical(c(r$rank, r$dim), c(1L, 3L))
  expect_identical(r$info$coordinates, c("L", "c", "tau"))
  ch <- lapply(r$symmetries, `[[`, "inputChange")
  expect_true(all(vapply(ch, function(x) inherits(x, "eqnvec") && names(x) == "T", TRUE)))
  expect_identical(names(r$info$unknownInputs$jets), "T")
  expect_true(any(grepl("input change:", capture.output(print(r)))))
})

test_that("an identifiable model, and conditions that fix a scale", {
  r <- symmetryDetection(eqnvec(x = "-k*x"), eqnvec(y = "x"), verbose = FALSE)
  expect_true(r$identifiable)
  expect_length(r$symmetries, 0L)
  expect_match(capture.output(print(r))[1], "identifiable")
  r <- symmetryDetection(eqnvec(x = "-k*x"), eqnvec(y = "s*x"), verbose = FALSE)
  expect_false(r$identifiable)
  dose <- addEvent(eventlist(), var = "x", time = 0, value = "dose", method = "replace")
  r <- symmetryDetection(eqnvec(x = "b - a*x"), eqnvec(y = "s*x"), events = dose,
                         conditions = data.frame(dose = 2, row.names = "stim"),
                         verbose = FALSE)
  expect_true(r$identifiable)
})

test_that("events split the time line and equilibrate starts at rest", {
  f <- eqnvec(R = "kpr - kdg*R + kon*u*R", u = "0")
  ev <- addEvent(eventlist(), var = "u", time = 0, value = "init_u", method = "replace") |>
    addEvent(var = "u", time = 60, value = "0", method = "replace")
  r <- symmetryDetection(f, eqnvec(y = "scale*R"), equilibrate = TRUE, events = ev,
                         conditions = data.frame(init_u = 1, row.names = "Ctrl"),
                         forcings = "u", verbose = FALSE)
  expect_identical(r$info$segments, 2L)
  expect_true(isTRUE(r$info$settings$equilibrate))
})

test_that("the gap order is raised to its rank bound, else the rank is a lower bound", {
  f <- eqnvec(x = "u*w - k*x", w = "c*(1 - u)*v", v = "1 - u", u = "0")
  ev <- addEvent(eventlist(), var = "u", time = 0, value = "0", method = "replace") |>
    addEvent(var = "u", time = 1, value = "1", method = "replace")
  tr <- eqnvec(w = "0", v = "0", u = "0")
  r <- symmetryDetection(f, eqnvec(y = "x"), trafo = tr, events = ev, verbose = FALSE)
  expect_true(r$identifiable)
  expect_true(r$info$gapCertified)
  cap <- reconstControl(gapOrderCap = 1L)
  r <- symmetryDetection(f, eqnvec(y = "x"), trafo = tr, events = ev, reconstruct = TRUE,
                         verify = FALSE, control = cap, verbose = FALSE)
  expect_false(r$info$gapCertified)
  expect_match(r$info$gapCertificate, "below its bound 3")
  expect_identical(r$info$rankStatus$status, "lower bound")
  expect_false(r$identifiable)
  expect_true(any(grepl("^gap order: ", unlist(r$info$limits))))
  expect_null(r$info$bound$probability)
  expect_true(any(grepl("gap order is not proven", capture.output(print(r)))))
  expect_warning(r <- symmetryDetection(f, eqnvec(y = "x"), trafo = tr, events = ev, reconstruct = TRUE,
                                        control = cap, verbose = FALSE),
                 "still grows")
  expect_identical(r$info$verification$method, "gap guard")
  r <- symmetryDetection(f_ab, eqnvec(y = "s*A"), reconstruct = TRUE, verbose = FALSE)
  expect_null(r$info$gapCertified)
  expect_length(r$info$branches, 0L)
})

test_that("regimes name the regimes an experiment visits", {
  expect_identical(dMod2:::.symRegimes(c("x > c", "x < c")), list("x > c", "x < c"))
  expect_identical(dMod2:::.symRegimes(list(a = "x > c", b = list(c("x > c", "y > 0")))),
                   list(a = list("x > c"), b = list(list("x > c", "y > 0"))))
  r <- symmetryDetection(eqnvec(x = "piecewise(-k1*x, x > c, -k2*x)"), eqnvec(y = "x"),
                         trafo = eqnvec(x = "2"), regimes = "x > c", verbose = FALSE)
  expect_true(all(c("c", "k2") %in% unlist(lapply(r$symmetries, `[[`, "support"))))
})

test_that("scalings only, the gauge and the control settings", {
  r <- symmetryDetection(f_ab, eqnvec(y = "s*A"), scalingsOnly = TRUE, verbose = FALSE)
  expect_identical(r$identifiable, NA)
  expect_true(all(vapply(r$symmetries, `[[`, "", "type") == "scaling"))
  r <- symmetryDetection(f_ab, eqnvec(y = "s*A"), gaugePreference = NULL,
                         control = reconstControl(degreeCap = 3L, timeout = 100),
                         verbose = FALSE)
  expect_length(r$gauge, 1L)
  expect_identical(r$info$settings$degreeCap, 3L)
  r <- symmetryDetection(f_ab, eqnvec(y = "s*A"), reconstruct = TRUE,
                         control = reconstControl(degreeCapMax = 6L, seed = 2L), verbose = FALSE)
  expect_identical(r$info$settings$control$degreeCapMax, 6L)
})

test_that("the reduction removes the directions and reports its chart", {
  r <- symmetryDetection(f_ab, eqnvec(y = "s*A"), reconstruct = TRUE, verbose = FALSE)
  red <- symmetryReduction(r, reportZeroCompatibility = TRUE)
  expect_s3_class(red, "symmetryreduction")
  expect_identical(red$removed, vapply(red$blocks, function(b) b$labels, ""))
  expect_length(red$remaining, 0L)
  expect_identical(red$partial, character(0))
  expect_s3_class(red$trafo, "eqnvec")
  expect_identical(unname(unclass(red$trafo)["A"]), "1")
  expect_true(all(c("block", "coordinates", "verdict", "condition") %in%
                    names(red$zeroCompatibility)))
  expect_true(any(grepl("Trafo", capture.output(print(red)))))
  expect_true(any(grepl("Blocks", capture.output(summary(red)))))
  fx <- symmetryReduction(r, fixed = "s")
  expect_identical(unname(unclass(fx$trafo)["s"]), "s")
  expect_error(symmetryReduction(list()), "symmetrydetection")
  expect_true(red$level %in% c("exact", "open"))
  expect_true(isTRUE(red$assumptions$domain))
  expect_true(all(vapply(red$blocks, function(b) length(b$unproven) == 0L || red$level == "open", TRUE)))
  expect_true(all(vapply(red$blocks, function(b) b$generatorLevel, "") == "exact"))
  expect_length(red$limits, 0L)
  lim <- symmetryReduction(r, workCap = 0L, check = FALSE)
  expect_true(any(grepl("work_cap=0", unlist(lim$limits))))
  expect_true(red$checks$pins && red$checks$invariants)
  expect_true(red$checks$loop$ok)
  expect_identical(red$checks$loop$rankStatus, "proven")
  expect_null(symmetryReduction(r, check = FALSE)$checks$loop)
})

test_that("the reduction takes the declared positive coordinates of the detection", {
  pos <- c("k", "e", "s", "a")
  r <- symmetryDetection(eqnvec(x = "-k*x", z = "a*x - e*z"), eqnvec(y = "s*z"), positive = pos,
                         reconstruct = TRUE, verbose = FALSE)
  expect_identical(as.character(unlist(symmetryReduction(r)$assumptions$domain)), pos)
  expect_true(isTRUE(symmetryReduction(r, positive = TRUE)$assumptions$domain))
})

test_that("the reduction refuses open directions unless allowOpen", {
  r <- symmetryDetection(f_ab, eqnvec(y = "s*A"), verbose = FALSE)
  expect_error(symmetryReduction(r), "no closed form")
  red <- symmetryReduction(r, allowOpen = TRUE)
  expect_identical(red$removed, "X\u2081")
  expect_identical(red$remaining, "X\u2082")
})
