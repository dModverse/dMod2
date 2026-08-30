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
