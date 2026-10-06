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
                                   generator = list(k_off = "-k_off"))))
  y <- dMod2:::.symCamelFields(x)
  expect_identical(names(y$info), c("lieOrderUsed", "modelExprs"))
  expect_identical(names(y$symmetries[[1]]), c("completeGenerator", "generator"))
  expect_identical(names(y$symmetries[[1]]$completeGenerator), "k_on")
  expect_identical(names(y$symmetries[[1]]$generator), "k_off")
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

test_that("scalings only, the gauge and the control settings", {
  r <- symmetryDetection(f_ab, eqnvec(y = "s*A"), scalingsOnly = TRUE, verbose = FALSE)
  expect_identical(r$identifiable, NA)
  expect_true(all(vapply(r$symmetries, `[[`, "", "type") == "scaling"))
  r <- symmetryDetection(f_ab, eqnvec(y = "s*A"), gaugePreference = NULL,
                         control = reconstControl(degreeCap = 3L, timeout = 100),
                         verbose = FALSE)
  expect_length(r$gauge, 1L)
  expect_identical(r$info$settings$degreeCap, 3L)
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
})
