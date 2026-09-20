test_that(".remoteLibs puts the paths in front of R_LIBS", {
  expect_identical(dMod2:::.remoteLibs(NULL), "")
  expect_identical(dMod2:::.remoteLibs(c("~/lib", "/opt/r")),
                   'export R_LIBS="$HOME/lib:/opt/r${R_LIBS:+:$R_LIBS}"; ')
})
