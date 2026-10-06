# Compile jobs per test file: DMOD_TEST_CORES, four by default.
test_cores <- function() {
  n <- suppressWarnings(as.integer(Sys.getenv("DMOD_TEST_CORES", "4")))
  if (is.na(n) || n < 1L) 4L else n
}
