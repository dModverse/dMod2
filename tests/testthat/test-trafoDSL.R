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

test_that("identifiers may end in an underscore", {

  # SBML ids may end in underscores, which must survive substitution.
  trafo <- define(eqnvec(), "x~x", x = c("gamma_", "I0_", "a__b"))
  trafo <- define(trafo, "S_ ~ N_ - I0_")

  expect_equal(insert(trafo, "gamma_ ~ gamma_CA")[["gamma_"]], "gamma_CA")
  expect_equal(repar(trafo, "I0_ ~ I0_NY")[["S_"]], "N_ - I0_NY")
  expect_equal(insert(trafo, "x ~ exp(x)", x = "a__b")[["a__b"]], "exp(a__b)")

})

test_that("substitution into compound identifiers still works", {

  trafo <- define(eqnvec(), "x~x", x = c("alpha", "k_pr_21mRNA"))
  out   <- insert(trafo, "x ~ x + Delta_x_condition", x = "alpha", condition = "C1")

  expect_equal(out[["alpha"]], "alpha + Delta_alpha_C1")
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

test_that("a compound replacement keeps its precedence", {

  trafo <- define(eqnvec(), "x~x", x = c("a", "b"))
  trafo <- define(trafo, "c ~ x/tot*k^tot", x = "a")
  out   <- insert(trafo, "tot ~ 500*(1 + tot)")

  expect_equal(out[["c"]], "a/(500 * (1 + tot)) * k^(500 * (1 + tot))")
  expect_equal(insert(trafo, "a ~ exp(a)")[["c"]], "exp(a)/tot * k^tot")
  expect_equal(insert(trafo, "tot ~ -1")[["c"]], "a/(-1) * k^(-1)")

})

test_that("replaceSymbols replaces variables only", {

  expect_equal(replaceSymbols("k", "p+q", "f(k = k) + l$k + exp(k)"),
               "f(k = (p+q)) + l$k + exp(p+q)")
  expect_equal(replaceSymbols(c("a", "b"), c("b", "a"), "a - b"), "b - a")
  x <- eqnvec(A = "k1*B", B = "0")
  out <- replaceSymbols("k1", "k2 + k3", x)
  expect_s3_class(out, "eqnvec")
  expect_equal(unclass(out), c(A = "(k2 + k3)*B", B = "0"))

})

test_that("getSymbols lists variables once, in order", {

  expect_equal(getSymbols(c("k1*A/(Km + A)", "0", "exp(-k2*time)")),
               c("k1", "A", "Km", "k2", "time"))
  expect_equal(getSymbols("k1*A + f(x = B)", exclude = "A"), c("k1", "B"))
  expect_equal(getSymbols("2*pi*k*time"), c("k", "time"))
  expect_null(getSymbols(NULL))

})

test_that("conditionMatch needs a trafo with a condition tree", {
  tr <- eqnvec(x = "x", y = "y")
  expect_error(define(tr, "x ~ 2", conditionMatch = "C"), "branch")
  expect_error(insert(tr, "x ~ a", a = "b", conditionMatch = "C"), "branch")
})
