test_that("Direct sum of datalists", {
  
  
  # Start with two data frames
  mydata1 <- data.frame(
    name = "A",
    time = 0:1,
    value = 1:2,
    sigma = .1,
    compound = c("DEM", "APAP"),
    dose = "0.1"
  )
  
  mydata2 <- data.frame(
    name = "A",
    time = 0:1,
    value = 3:4,
    sigma = .1,
    compound = c("APAP", "DCF"),
    dose = "0.1"
  )
 
  # Create datalists from dataframes
  data1 <- as.datalist(mydata1, splitBy = c("compound", "dose")) 
  data2 <- as.datalist(mydata2, splitBy = c("compound", "dose")) 
  
  

  
  data <- suppressWarnings(data1 + data2)
  condition.grid <- attr(data, "condition.grid")
  
  expect_equal(nrow(condition.grid), 3)
  expect_equal(length(data), 3)
  expect_warning(data1 + data2, "Condition .* existed and has been overwritten.")
  expect_equal(nrow(data[["APAP_0.1"]]), 1)
  })
