## Past the shell's argument limit compile() links through a chunked static
## archive instead of naming every object on the SHLIB command line.

test_that(".compileChunks caps both count and command length", {

  files <- sprintf("/a/b/file%04d.o", 1:250)

  chunks <- dMod2:::.compileChunks(files, maxN = 100L)
  expect_equal(lengths(chunks), c(100L, 100L, 50L))
  expect_identical(unlist(chunks, use.names = FALSE), files)

  ## An over-long name forms a chunk of its own rather than an empty one.
  long <- c(strrep("x", 500), files[1:3])
  expect_equal(lengths(dMod2:::.compileChunks(long, maxChars = 100L, maxN = 100L)),
               c(1L, 3L))

  expect_identical(dMod2:::.compileChunks(character(0)), list())

})


test_that(".mergeCompileInfo dedups by srcfile, first flags win", {

  a <- list(list(srcfile = "x.cpp", compileArgs = "A"),
            list(srcfile = "y.cpp", compileArgs = "B"))
  b <- list(list(srcfile = "y.cpp", compileArgs = "OVERRIDDEN"),
            list(srcfile = "z.cpp", compileArgs = "C"),
            list(srcfile = "z.cpp", compileArgs = "D"))

  m <- dMod2:::.mergeCompileInfo(a, b)
  expect_equal(vapply(m, function(e) e$srcfile,     ""), c("x.cpp", "y.cpp", "z.cpp"))
  expect_equal(vapply(m, function(e) e$compileArgs, ""), c("A", "B", "C"))

  ## The cached keys must stay in step with the list they annotate.
  expect_equal(attr(m, "srckeys"), c("x.cpp", "y.cpp", "z.cpp"))
  expect_equal(dMod2:::.compileInfoKeys(m), attr(m, "srckeys"))
  expect_equal(dMod2:::.mergeCompileInfo(m, m), m)

  expect_null(dMod2:::.mergeCompileInfo(NULL, NULL))
  expect_null(dMod2:::.mergeCompileInfo(list(list(srcfile = character(0))), NULL))

})


test_that("compile() links many conditions through a chunked archive", {

  withr::local_dir(tempdir())
  ## Force the archive route without generating a thousand sources.
  withr::local_options(dMod.compile.cmdlimit = 1L)

  conditions <- sprintf("chk%02d", 1:12)
  p <- Reduce("+", lapply(conditions, function(cn)
    P(eqnvec(k1 = "exp(logk1)", A0 = "exp(logA0)*scale"),
      condition = cn, compile = FALSE, modelname = paste0("chunked_p_", cn))))

  expect_length(attr(p, "compileInfo"), length(conditions))
  expect_output(compile(p, output = "chunked_all", cores = 1), "archived 11 objects")

  expect_true(file.exists(paste0("chunked_all", .Platform$dynlib.ext)))
  expect_equal(list.files(pattern = "^chunked_all.*\\.a$"), character(0))
  expect_equal(unique(modelname(p)), "chunked_all")

  ## Every condition's entry point resolves out of the single shared object.
  out <- p(c(logk1 = log(2), logA0 = log(3), scale = 5), deriv = TRUE)
  expect_equal(vapply(conditions, function(cn) unname(out[[cn]]["A0"]), 0),
               setNames(rep(15, length(conditions)), conditions))
  J <- attr(out[[conditions[length(conditions)]]], "deriv")
  expect_equal(J["A0", "scale"], 3)
  expect_equal(J["k1", "logk1"], 2)

})


test_that("compile() refuses more shared objects than R can load", {

  withr::local_dir(tempdir())
  withr::local_envvar(R_MAX_NUM_DLLS = "1")

  p <- P(eqnvec(k1 = "exp(logk1)"), condition = "dll",
         compile = FALSE, modelname = "chunked_budget")

  expect_error(compile(p), "R_MAX_NUM_DLLS")
  expect_error(compile(p), "output = <name>", fixed = TRUE)

})


test_that("compile() reuses objects whose source and command are unchanged", {

  withr::local_dir(withr::local_tempdir())

  trafo <- eqnvec(k1 = "exp(logk1)", A0 = "exp(logA0)*scale")
  conditions <- sprintf("reuse%02d", 1:9)
  mk <- function() Reduce("+", lapply(conditions, function(cn)
    P(trafo, condition = cn, compile = FALSE, modelname = paste0("reuse_p_", cn))))

  p <- mk()
  expect_no_match(capture.output(compile(p, output = "reuse_all", cores = 1)), "reusing")

  ## Codegen rewrites every source, so only the content decides.
  p <- mk()
  expect_output(compile(p, output = "reuse_all", cores = 1),
                paste("reusing", length(conditions), "unchanged"))
  out <- p(c(logk1 = log(2), logA0 = log(3), scale = 5), deriv = TRUE)
  expect_equal(as.numeric(out[[conditions[1]]]["A0"]), 15)
  expect_equal(attr(out[[conditions[9]]], "deriv")["k1", "logk1"], 2)

  ## A changed compile command has to invalidate the entry.
  expect_no_match(capture.output(compile(p, output = "reuse_all", cores = 1,
                                         args = "-DDMOD_TEST_FLAG")),
                  "reusing")

})


test_that(".compilePCHIncludes only accepts prologues it can safely prepend", {

  d <- withr::local_tempdir()
  clean <- file.path(d, "clean.cpp")
  writeLines(c("/** generated **/", "", "#include <cmath>", "#include <vector>",
               "double f() { return 0; }"), clean)
  expect_equal(dMod2:::.compilePCHIncludes(clean), c("#include <cmath>", "#include <vector>"))

  ## A macro ahead of the includes can change how they expand, so no header.
  guarded <- file.path(d, "guarded.cpp")
  writeLines(c("#define _USE_MATH_DEFINES", "#include <cmath>", "double f() { return 0; }"),
             guarded)
  expect_null(dMod2:::.compilePCHIncludes(c(clean, guarded)))
  expect_null(dMod2:::.compilePCHIncludes(file.path(d, "none.cpp") |>
                                    (\(f) { writeLines("int f();", f); f })()))

})


skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("cppDE")
  testthat::skip_on_cran()
}


## `...` carries the objects, so a misspelled argument name would be taken for
## one and dropped without a word.

test_that("compile rejects a named argument that carries no sources", {
  fake <- structure(list(), class = c("obsfn", "fn"))

  expect_error(compile(fake, outout = "all"), "did you mean `output`")
  expect_error(compile(fake, coress = 2), "did you mean `cores`")
  expect_error(compile(fake, wumpus = 1), "`wumpus`")
  # naming the objects themselves is not what the check is after
  expect_false(grepl("did you mean",
                     tryCatch({ compile(myfn = fake); "" },
                              error = function(e) conditionMessage(e))))
})


test_that("the toolchain report separates shared from per-file flags", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()

  out <- capture.output(compile(bench$gfn, bench$xfn, cores = 1))
  head  <- grep("^using C\\+\\+ compiler:", out, value = TRUE)
  every <- grep("^ +every source +:", out, value = TRUE)
  also  <- grep("^ +[0-9]+ of [0-9]+ also +:", out, value = TRUE)

  expect_length(head, 1L)

  # No source declares a flag of its own, a toolchain without OpenMP for
  # instance, and the report collapses to the one-line form that names the
  # shared set in brackets.
  if (!length(also)) {
    expect_match(head, "[", fixed = TRUE)
    expect_length(every, 0L)
    return(invisible(NULL))
  }

  # no bracket on the header once the sets differ, it would name no real command
  expect_false(grepl("\\[", head))
  expect_length(every, 1L)
  expect_false(grepl("fopenmp", every))
  expect_length(also, 1L)
  expect_true(grepl("-fopenmp", also))
})


test_that("output places the shared object where it says", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  so <- .Platform$dynlib.ext
  d <- file.path(tempdir(), paste0("cmp_", as.integer(runif(1, 1e6, 9e6))))
  dir.create(d)

  # a name carrying a directory is taken as given
  invisible(capture.output(
    compile(bench$gfn, bench$xfn, output = file.path(d, "combined"), cores = 1)))
  expect_true(file.exists(file.path(d, paste0("combined", so))))

  # a bare name lands next to the sources
  src <- attr(bench$xfn, "compileInfo")[[1]]$srcfile
  invisible(capture.output(compile(bench$gfn, bench$xfn, output = "bare", cores = 1)))
  expect_true(file.exists(file.path(dirname(src), paste0("bare", so))))

  expect_error(compile(bench$gfn, output = file.path(d, "nope", "x")),
               "does not exist")
  expect_error(compile(bench$gfn, output = character(0)), "single non-empty name")
})
