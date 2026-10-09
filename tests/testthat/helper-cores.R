# Compile jobs per test file: DMOD_TEST_CORES, two by default.
test_cores <- function() {
  n <- suppressWarnings(as.integer(Sys.getenv("DMOD_TEST_CORES", "2")))
  if (is.na(n) || n < 1L) 2L else n
}
