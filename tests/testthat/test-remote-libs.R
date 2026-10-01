test_that(".remoteLibs puts the paths in front of R_LIBS", {
  expect_identical(dMod2:::.remoteLibs(NULL), "")
  expect_identical(dMod2:::.remoteLibs(c("~/lib", "/opt/r")),
                   'export R_LIBS="$HOME/lib:/opt/r${R_LIBS:+:$R_LIBS}"; ')
})

test_that(".remoteBuildInfo drops local paths and flags the LAPACK solver", {
  env <- new.env()
  env$m <- structure(list(), compileInfo = list(list(
    srcfile = "m.cpp",
    compileArgs = "-DFOO -I/home/u/.cache/R/cppDE/include -fopenmp -I/usr/include/openmpi",
    linkArgs = "-lsundials_sunlinsollapackdense -L/home/u/lib -lsundials_cvodes")))
  info <- dMod2:::.remoteBuildInfo(env)
  expect_identical(info$compileArgs, "-DFOO -fopenmp")
  expect_true(info$needsCVODE)
  expect_true(info$needsLapack)
  expect_false(info$needsKLU)
})

test_that("the remote build script resolves SUNDIALS flags on the remote side", {
  script <- dMod2:::.remoteBuildScript("m.cpp", "m.so", compileArgs = "-DFOO",
                                       needsCVODE = TRUE, needsLapack = TRUE)
  expect_match(script, "\"cflags\"", fixed = TRUE)
  expect_match(script, "\"cvode_lapack_libs\", \"libs\"", fixed = TRUE)
  expect_match(script, "LAPACK dense solver", fixed = TRUE)
  expect_false(grepl("\"cflags\"", dMod2:::.remoteBuildScript("m.c", "m.so"),
                     fixed = TRUE))
})

test_that("the generated include line evaluates to the installed cppDE flags", {
  skip_on_os("windows")
  cfg <- get0("cvodeConfig", envir = asNamespace("cppDE"), inherits = FALSE)
  skip_if_not(is.environment(cfg) && isTRUE(cfg$available))
  script <- dMod2:::.remoteBuildScript("m.cpp", "m.so", needsCVODE = TRUE)
  line <- grep("^PKG_CPPFLAGS=\"\\$PKG_CPPFLAGS ", strsplit(script, "\n")[[1]],
               value = TRUE)
  expect_length(line, 1L)
  out <- system2("bash", c("-c", shQuote(paste0("PKG_CPPFLAGS=; ", line,
                                                "; printf %s \"$PKG_CPPFLAGS\""))),
                 stdout = TRUE)
  expect_identical(trimws(out), trimws(cfg$cflags))
})
