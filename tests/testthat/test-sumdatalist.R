test_that("the sum of datalists unites the conditions, the right summand wins a shared one", {
  first <- data.frame(name = "A", time = c(0, 1), value = c(1, 2), sigma = 0.1,
                      drug = c("u", "v"), level = "1")
  second <- data.frame(name = "A", time = c(0, 1), value = c(3, 4), sigma = 0.1,
                       drug = c("v", "w"), level = "1")
  d1 <- as.datalist(first, splitBy = c("drug", "level"))
  d2 <- as.datalist(second, splitBy = c("drug", "level"))

  expect_warning(total <- d1 + d2, "existed and has been overwritten")
  expect_setequal(names(total), c("u_1", "v_1", "w_1"))
  expect_identical(nrow(attr(total, "condition.grid")), 3L)
  expect_identical(total[["v_1"]]$value, 3)
  expect_identical(total[["u_1"]]$value, 1)
})
