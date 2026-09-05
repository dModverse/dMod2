test_that("insert on a symbol that is not in the trafo is a no-op", {

  trafo <- define(eqnvec(), "x~x", x = c("k_pr_R2mRNA", "alpha", "offset_a"))

  expect_identical(insert(trafo, "xxx ~ 1"), trafo)
  expect_identical(insert(trafo, "k_pr_R2mRNA_Pert ~ 1"), trafo)

})

test_that("compound identifiers may contain digit-leading parts", {

  # "21mRNA" is not a syntactic name on its own, but "k_pr_21mRNA" is.
  trafo <- define(eqnvec(), "x~x",
                  x = c("k_pr_21mRNA", "k_pr_21mRNA_Pert", "alpha"))

  expect_equal(trafo[["k_pr_21mRNA"]], "k_pr_21mRNA")
  expect_equal(insert(trafo, "k_pr_21mRNA ~ 1")[["k_pr_21mRNA"]], "1")
  expect_equal(insert(trafo, "x ~ exp10(x)", x = "k_pr_21mRNA_Pert")[["k_pr_21mRNA_Pert"]],
               "exp10(k_pr_21mRNA_Pert)")

})

test_that("substitution into compound identifiers still works", {

  trafo <- define(eqnvec(), "x~x", x = c("alpha", "k_pr_21mRNA"))
  out   <- insert(trafo, "x ~ x + Delta_x_condition", x = "alpha", condition = "C1")

  expect_equal(out[["alpha"]], "alpha+Delta_alpha_C1")
  expect_equal(out[["k_pr_21mRNA"]], "k_pr_21mRNA")

})

test_that("numeric literals survive the compound-identifier rewrite", {

  trafo <- define(eqnvec(), "x~x", x = "alpha")
  expect_equal(insert(trafo, "alpha ~ 1e-4 * alpha + 0.5")[["alpha"]],
               "1e-04 * alpha + 0.5")

})

test_that("branch tolerates grid columns that match no parameter", {

  trafo <- define(eqnvec(), "x~x", x = c("k_pr_R2mRNA_Pert", "alpha"))
  grid  <- data.frame(k_pr_21mRNA_Pert = c("1", "k_knd"),
                      row.names = c("C1", "C2"))

  out <- branch(trafo, table = grid, apply = "insert")

  expect_named(out, c("C1", "C2"))
  expect_equal(unname(out[["C1"]][["k_pr_R2mRNA_Pert"]]), "k_pr_R2mRNA_Pert")
  expect_equal(unname(out[["C2"]][["k_pr_R2mRNA_Pert"]]), "k_pr_R2mRNA_Pert")

})

test_that("dots are resolved in the calling frame", {

  # A caller inside a function must see its own variables, not only globals.
  build <- function() {
    pars <- c("alpha", "beta")
    repl <- "0.5"
    define(eqnvec(), "x~x", x = pars) |> insert("beta ~ v", v = repl)
  }

  out <- build()

  expect_equal(out[["alpha"]], "alpha")
  expect_equal(out[["beta"]], "0.5")

})

test_that("condition columns and .currentSymbols outrank the calling frame", {

  grid <- data.frame(shift = c("s1", "s2"), row.names = c("C1", "C2"))

  build <- function() {
    shift <- "caller"
    define(eqnvec(), "x~x", x = c("alpha", "beta")) |>
      branch(table = grid, apply = "nothing") |>
      insert("x ~ x_y", x = .currentSymbols, y = shift)
  }

  out <- build()

  expect_equal(out[["C1"]][["alpha"]], "alpha_s1")
  expect_equal(out[["C2"]][["beta"]], "beta_s2")

})


test_that("subset() on an eqnlist sees the calling frame", {

  eq <- eqnlist() |>
    addReaction("A", "B", "k1*A") |>
    addReaction("B", "C", "k2*B")

  # The condition may name columns of the reaction table, may use the `%in%`
  # this method overloads, and may name variables of the calling frame.
  byRate  <- function(pattern) subset(eq, grepl(pattern, Rate))
  byEduct <- function(species) subset(eq, species %in% Educt)

  expect_equal(nrow(getReactions(byRate("k1"))), 1L)
  expect_equal(nrow(getReactions(byEduct("A"))), 1L)
  expect_equal(nrow(getReactions(subset(eq, grepl("k", Rate)))), 2L)

})
