## Characterisation of the equation classes and helpers (eqnClass.R, linalg.R,
## symbolics.R, utils.R, eventClass.R).

eqnNet <- function() {
  eq <- eqnlist()
  eq <- addReaction(eq, "A", "pA", "act_A * A * stimulus", "Phosphorylation of A")
  eq <- addReaction(eq, "pA", "A", "deact_A * pA", "Dephosphorylation of pA")
  eq <- addReaction(eq, "2*pA", "pA_pA", "form_complex_pA * pA^2", "Complex formation of pA")
  eq <- addReaction(eq, "B", "pB", "act_B * B * pA_pA", "Phosphorylation of B")
  addReaction(eq, "pB", "B", "deact_B * pB", "Dephosphorylation of pB")
}

test_that("eqnlist() checks dimensions and numbers descriptions", {
  S <- matrix(c(-1, 1, 0, 0, -1, 1), 2, 3, byrow = TRUE)
  expect_error(eqnlist(S, c("A", "B"), c("r1", "r2")),
               "Number of states does not coincide with number of columns of stoichiometric matrix")
  expect_error(eqnlist(S, c("A", "B", "C"), "r1"),
               "Number of rates does not coincide with number of rows of stoichiometric matrix")
  el <- eqnlist(S, c("A", "B", "C"), c("k1*A", "k2*B"))
  expect_identical(el$description, c("1", "2"))
  expect_identical(el$smatrix,
                   matrix(c(-1, 1, NA, NA, -1, 1), 2, 3, byrow = TRUE,
                          dimnames = list(NULL, c("A", "B", "C"))))
  expect_true(is.eqnlist(el))
  el2 <- eqnlist(S, c("A", "B", "C"), c("k1*A", "k2*B"), description = c("one", "two"))
  expect_identical(el2$description, c("one", "two"))
  expect_identical(class(el), c("eqnlist", "list"))
  df <- data.frame(a = c(-1, 0), b = c(1, -1))
  el3 <- eqnlist(df, rates = c("k1*a", "k2*b"))
  expect_true(is.matrix(el3$smatrix))
  expect_identical(el3$states, c("a", "b"))
  e0 <- eqnlist()
  expect_null(e0$smatrix)
  expect_identical(e0$description, character(0))
  expect_true(is.eqnlist(e0))
  # the description length is not checked
  expect_identical(eqnlist(S, c("A", "B", "C"), c("k1*A", "k2*B"), description = "x")$description, "x")
})

test_that("is.eqnlist() and is.eqnvec() recognise valid objects only", {
  el <- eqnNet()
  expect_true(is.eqnlist(el))
  expect_false(is.eqnlist(unclass(el)))
  bad <- el; bad$rates <- bad$rates[-1]
  expect_false(is.eqnlist(bad))
  bad <- el; bad$smatrix <- as.vector(bad$smatrix)
  expect_false(is.eqnlist(bad))
  bad <- el; bad$states <- rev(bad$states)
  expect_false(is.eqnlist(bad))
  bad <- el; bad$compartmentOf <- NULL
  expect_false(is.eqnlist(bad))
  bad <- el; bad$volumes <- NULL
  expect_false(is.eqnlist(bad))
  bad <- el; bad$amountStates <- "nope"
  expect_false(is.eqnlist(bad))
  bad <- el; bad$reactionCompartment <- rep("nowhere", length(el$rates))
  expect_false(is.eqnlist(bad))
  e0 <- eqnlist()
  bad <- e0; bad$rates <- "k"
  expect_false(is.eqnlist(bad))
  bad <- e0; bad$volumes <- c(A = "1")
  expect_false(is.eqnlist(bad))
  bad <- e0; bad$description <- "x"
  expect_false(is.eqnlist(bad))
  bad <- e0; bad$reactionCompartment <- "c"
  expect_false(is.eqnlist(bad))
  bad <- e0; bad$compartments <- list(c1 = list(volume = "1", rule = NULL))
  expect_false(is.eqnlist(bad))
  bad$compartmentOf <- c(A = "c1")
  expect_true(is.eqnlist(bad))
  bad$compartmentOf <- c(A = "c2")
  expect_false(is.eqnlist(bad))

  # the reaction-less branch does not test the class
  expect_true(is.eqnlist(list()))
  expect_true(is.eqnlist(unclass(eqnlist())))
  expect_true(is.eqnvec(eqnvec(a = "1")))
  expect_false(is.eqnvec(c(a = "1")))
  expect_false(is.eqnvec(structure("1", class = c("eqnvec", "character"))))
  expect_identical(is.eqnvec(NULL), FALSE)
})

test_that("c.eqnvec concatenates and refuses duplicate names", {
  out <- c(eqnvec(a = "x"), eqnvec(b = "y+1"))
  expect_identical(out, structure(c(a = "x", b = "y+1"), class = c("eqnvec", "character")))
  expect_error(c(eqnvec(a = "x"), eqnvec(a = "y")), "Names must be unique")
  out2 <- c(eqnvec(a = "x"), c(b = "z"))
  expect_s3_class(out2, "eqnvec")
  expect_error(c(eqnvec(a = "x"), "y"), "every equation needs a name")
  expect_identical(names(out2), c("a", "b"))
})

test_that("as.eqnvec.character validates names and parses", {
  expect_error(as.eqnvec(c(a = "x+")), "equations cannot be parsed: ")
  expect_error(as.eqnvec("x"), "equations need names")
  expect_error(as.eqnvec(c("x", "y"), names = "a"), "Length of names and equations do not coincide")
  expect_error(as.eqnvec(c(a = "x", a = "y")), "equation names are duplicated: a")
  expect_null(as.eqnvec.character(NULL))
  expect_identical(as.eqnvec("x", names = "b"),
                   structure(c(b = "x"), class = c("eqnvec", "character")))
})

test_that("conservedQuantities returns the kernel as named expressions", {
  cq <- conservedQuantities(eqnNet())
  expect_identical(unname(unlist(cq)), c("0.5*A+0.5*pA+1*pA_pA", "1*B+1*pB"))
  expect_identical(names(cq), c("total_1", "total_2"))
  expect_identical(conservedQuantities(addReaction(eqnlist(), "A", "B", "k*A") |> addReaction("B", "", "d*B")),
                   list())
  S2 <- rbind(c(1, -1, 0, 1)); colnames(S2) <- c("w", "x", "y", "z")
  e2 <- eqnlist(S2, rates = "k*w", description = "r1")
  expect_identical(unname(unlist(conservedQuantities(e2))), c("1*w+1*x", "1*y", "-1*w+1*z"))
  S3 <- rbind(c(1, NA, -1)); colnames(S3) <- c("a", "b", "c")
  expect_identical(unname(unlist(conservedQuantities(eqnlist(S3, rates = "k*c", description = "r1")))),
                   c("1*b", "1*a+1*c"))
  expect_error(conservedQuantities(S3), "must be an eqnlist")
})

test_that("getSymbols skips zero entries and excluded symbols", {
  expect_identical(getSymbols(c("0", "k1*A/(Km + A)", "exp(-k2*time)")),
                   c("k1", "A", "Km", "k2", "time"))
  expect_identical(getSymbols(c("a*pi", "0"), exclude = "b"), "a")
  expect_null(getSymbols(NULL))
  expect_null(getSymbols("0"))
  expect_identical(getSymbols(c(x = "f(a = b)")), "b")
})

test_that("attrs keeps or drops the selected attributes", {
  x <- structure(1:3, names = c("a", "b", "c"), note = "temporary", dim2 = 1)
  expect_identical(attrs(x), structure(1:3, names = c("a", "b", "c")))
  expect_identical(attrs(x, "note"), structure(1:3, note = "temporary"))
  expect_identical(attrs(x, "note", keep = FALSE),
                   structure(1:3, names = c("a", "b", "c"), dim2 = 1))
  m <- structure(matrix(1:4, 2, dimnames = list(c("r", "s"), NULL)), extra = TRUE)
  expect_identical(attrs(m), matrix(1:4, 2, dimnames = list(c("r", "s"), NULL)))
  df <- structure(data.frame(a = 1), tag = "x")
  expect_identical(attrs(df), data.frame(a = 1))
  expect_identical(attrs(1:2), 1:2)
  expect_identical(attrs(x, character(0), keep = FALSE), x)
})

test_that("nullZ returns the rref null space basis", {
  A <- rbind(c(1, -1, 0), c(0, 1, -1))
  expect_identical(nullZ(A), matrix(c(1, 1, 1), 3, 1))
  expect_identical(dim(nullZ(diag(3))), c(3L, 0L))
  expect_identical(nullZ(matrix(0, 1, 2)), diag(2))
  B <- rbind(c(1, 2, 3), c(2, 4, 6))
  expect_equal(B %*% nullZ(B), matrix(0, 2, 2))
})

test_that("eventlist helpers build and coerce event tables", {
  ev <- eventlist(var = "A", time = "5", value = 1, method = "add")
  expect_identical(as.list(as.data.frame(ev)),
                   list(var = "A", time = "5", value = 1, root = NA, method = "add"))
  ev <- addEvent(NULL, var = "A", time = "5", value = 1, method = "add")
  ev <- addEvent(ev, var = "A", time = "10", value = 1, method = "add")
  expect_identical(nrow(ev), 2L)
  expect_s3_class(ev, "eventlist")
  expect_identical(as.eventlist(as.data.frame(ev)), ev)
  expect_identical(as.eventlist(list()), eventlist())
  expect_error(as.eventlist(list(var = "A")), "x needs to provide var, time, value, and method.")
})
