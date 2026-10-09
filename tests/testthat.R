library(testthat)
library(dMod2)

# Some test files setwd(tempdir()) without restoring; restoring the original
# wd on exit keeps test_check's relative-path lookups consistent.
.tt_initial_wd <- getwd()
on.exit(setwd(.tt_initial_wd), add = TRUE)

test_check("dMod2")
