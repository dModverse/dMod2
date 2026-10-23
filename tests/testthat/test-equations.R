# Equation lists, equation vectors and their linear algebra, without compiled code.

.eqChain <- function() {
  eqnlist() |>
    addReaction("A", "B", "k1*A", "r1") |>
    addReaction("B", "C", "k2*B", "r2") |>
    addReaction("C", "D", "k3*C", "r3")
}

test_that("addReaction reads numeric coefficients of any form", {
  f <- eqnlist() |>
    addReaction("3/2*A", "B", "k*A") |>
    addReaction("(2)*A + C", "2*3*D", "k*A*C")
  expect_equal(unname(f$smatrix[1, "A"]), -1.5)
  expect_equal(unname(f$smatrix[2, c("A", "C", "D")]), c(-2, -1, 6))
  expect_error(addReaction(eqnlist(), "A*k", "B", "k*A"), "coefficient")
})

test_that("getReactions writes coefficients as they are", {
  f <- eqnlist() |> addReaction("11*A", "0.1*B", "k*A") |> addReaction("2.1*C", "21*D", "k2*C")
  r <- getReactions(f)
  expect_equal(r$Educt, c("11*A", "2.1*C"))
  expect_equal(r$Product, c("0.1*B", "21*D"))
  expect_equal(subset(f, grepl("k2", Rate))$rates, "k2*C")
})

test_that("subset of an eqnlist matches %in% on ordinary vectors as base R does", {
  f <- .eqChain()
  expect_equal(subset(f, Description %in% "r3")$description, "r3")
  expect_equal(subset(f, Description %in% c("r2", "r3"))$description, c("r2", "r3"))
  expect_equal(subset(f, "B" %in% Educt)$description, "r2")
})

test_that("conserved quantities have no rounding residues", {
  S <- matrix(c(0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -2, 0, 0, 0, 1, 0, 0, 0,
                1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, -1, 1, 0, -2, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 3,
                0, 0), 8, 10, dimnames = list(NULL, paste0("X", 1:10)))
  Z <- nullZ(S)
  expect_false(any(Z != 0 & abs(Z) < 1e-10))
  expect_equal(max(abs(S %*% Z)), 0)
  cq <- conservedQuantities(S)[[1]]
  expect_false(any(grepl("e-1", cq)))
})

test_that("eqnvec() names unnamed entries by position", {
  expect_equal(names(eqnvec(a = "1", "2")), c("a", "eqn2"))
  expect_equal(names(eqnvec(c(a = "1", b = "2"))), c("a", "b"))
  expect_equal(unclass(eqnvec(a = 1)), c(a = "1"))
})

test_that("reactions need a species and states cannot be named Rate or Description", {
  expect_error(addReaction(eqnlist(), "", "", "k"), "species")
  expect_error(addReaction(eqnlist(), "Rate", "B", "k*Rate"), "Rate")
})

test_that("empty equation and event lists are handled", {
  tf <- tempfile(fileext = ".csv")
  expect_error(write.eqnlist(eqnlist(), file = tf), "no reaction")
  expect_s3_class(as.eventlist(eventlist()), "eventlist")
})

test_that("is.eqnlist checks that the stoichiometric columns are the states", {
  f <- addReaction(eqnlist(), "A", "B", "k*A")
  expect_true(is.eqnlist(f))
  f$states <- rev(f$states)
  expect_false(is.eqnlist(f))
})

test_that("compartments or compartmentOf alone set the layout", {
  S <- matrix(c(-1, 1), 1, 2)
  e <- eqnlist(S, c("A", "B"), "k*A", compartments = list(cyt = "Vc"),
               reactionCompartment = "cyt")
  expect_setequal(names(e$compartments), c("cyt", "defaultComp"))
  expect_equal(unname(e$compartmentOf), c("defaultComp", "defaultComp"))
  expect_error(eqnlist(S, c("A", "B"), "k*A", compartmentOf = c(A = "cyt", B = "nuc")),
               "undefined compartments")
})

test_that("c() keeps the layout of an eqnlist without reactions and custom totals", {
  decl <- assignCompartment(eqnlist(), A = "cyt", volume = "Vc")
  r <- addReaction(eqnlist(), "A", "B", "k*A")
  expect_warning(both <- c(decl, r), "the layout wins")
  expect_equal(both$compartmentOf[["A"]], "cyt")
  f <- customTotals(addReaction(eqnlist(), "A", "B", "k1*A"), list(tot = "A + B"))
  g <- addReaction(eqnlist(), "C", "D", "k2*C")
  expect_true("tot" %in% names(getTotals(c(f, g))))
})

test_that("log10Transform keeps names that sympy knows as constants", {
  skip_if_not_installed("reticulate")
  expect_match(log10Transform(eqnvec(x = "-E^2*x"))[["x_l10"]], "E")
  expect_match(log10Transform(eqnvec(x = "I^2*x"))[["x_l10"]], "I")
})

test_that("format.eqnvec pads operators, not powers or exponents", {
  expect_equal(unname(format.eqnvec(c(a = "x**2*k-1e-3"))), "x^2 * k - 1e-3")
})

test_that("equation vectors need unique names and subset() returns an eqnlist", {
  expect_error(as.eqnvec(c("x", "y"), names = c("", "")), "name")
  expect_error(as.eqnvec(c("x", "y"), names = c("a", "a")), "duplicated")
  f <- addReaction(eqnlist(), "A", "B", "k*A")
  expect_s3_class(subset(f, "C" %in% Educt), "eqnlist")
  expect_s3_class(subset(eqnlist(), "C" %in% Educt), "eqnlist")
})
