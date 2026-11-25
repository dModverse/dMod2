# Effect formulas in insert() and their random-effect structure.

.cov <- data.frame(ID = c("s1", "s2", "s3"), wt = c(62, 75, 81), dose = c(100, 100, 200),
                   sex = c("f", "m", "f"), row.names = c("C1", "C2", "C3"))
.tr <- function() branch(define(NULL, "x <- x", x = c("ka", "cl", "v")), table = .cov)

test_that("an effect formula replaces the symbol by the linear predictor on the scale of the trafo", {
  tr <- insert(.tr(), "x <- exp(x)", x = "cl")
  tr <- insert(tr, cl ~ wt + sex + (1 + dose | ID), v ~ 1 + (1 | ID))
  expect_equal(as.character(tr$C2["cl"]),
               "exp((cl + cl_wt * 75 + cl_sexm + eta_cl_s2 + eta_cl_dose_s2 * 100))")
  expect_equal(as.character(tr$C3["v"]), "(v + eta_v_s3)")
  expect_equal(as.character(tr$C3["ka"]), "ka")
  re <- randomEffects(tr)
  expect_s3_class(re, "randomeffects")
  expect_equal(re$terms$stem, c("eta_cl", "eta_cl_dose", "eta_v"))
  expect_setequal(re$fixed, c("cl", "cl_wt", "cl_sexm", "v"))
  expect_equal(re$etas["C2", "eta_cl_dose"], "eta_cl_dose_s2")
  # the columns of one term are correlated, other terms independent
  expect_equal(lapply(re$blocks, function(b) re$terms$stem[b]), list(c("eta_cl", "eta_cl_dose"), "eta_v"))
  expect_equal(re$subjects$ID, list(s1 = "C1", s2 = "C2", s3 = "C3"))
})

test_that("a condition with a missing grouping belongs to no group and keeps the fixed effect", {
  cov <- .cov; cov$ID[2] <- NA
  tr <- insert(branch(define(NULL, "x <- x", x = c("ka", "cl")), table = cov), cl ~ 1 + (1 | ID))
  expect_equal(as.character(tr$C2["cl"]), "(cl)")
  re <- randomEffects(tr)
  expect_true(is.na(re$etas["C2", "eta_cl"]))
  expect_equal(re$subjects$ID, list(s1 = "C1", s3 = "C3"))
})

test_that("placeholders and sums on the left apply one right-hand side to many symbols", {
  a <- insert(.tr(), x ~ 1 + (1 | ID), x = c("ka", "cl"))
  b <- insert(.tr(), ka + cl ~ 1 + (1 | ID))
  expect_identical(a, b)
  re <- randomEffects(a)
  expect_equal(lapply(re$blocks, function(b) re$terms$stem[b]), list("eta_ka", "eta_cl"))
  re <- randomEffects(insert(.tr(), x ~ 1 + (1 | p | ID), x = c("ka", "cl", "v")))
  expect_equal(lapply(re$blocks, function(b) re$terms$stem[b]), list(c("eta_ka", "eta_cl", "eta_v")))
  expect_equal(as.character(insert(.tr(), cl ~ log(wt))$C1["cl"]),
               paste0("(cl + cl_log_wt * ", sprintf("%.17g", log(62)), ")"))
})

test_that("an ID correlates terms across formulas and gr() sets the distribution", {
  tr <- insert(.tr(), ka ~ 1 + (1 | p | ID), cl ~ 1 + (1 | p | ID),
               v ~ 1 + (1 | gr(ID, dist = spikeSlab())))
  re <- randomEffects(tr)
  expect_equal(lapply(re$blocks, function(b) re$terms$stem[b]), list(c("eta_ka", "eta_cl"), "eta_v"))
  expect_equal(re$terms$dist, c("normal", "normal", "spikeSlab"))
  expect_true(is.na(re$terms$prob[3]))
  share <- 0.2
  re <- randomEffects(insert(.tr(), v ~ 1 + (1 | gr(ID, dist = spikeSlab(prob = share)))))
  expect_equal(re$terms$prob, 0.2)
})

test_that("P() keeps the random effects and effect formulas stack over calls", {
  tr <- insert(insert(.tr(), ka ~ 1 + (1 | ID)), cl ~ wt + (1 | ID))
  expect_equal(randomEffects(tr)$terms$stem, c("eta_ka", "eta_cl"))
  expect_equal(as.character(tr$C3["cl"]), "(cl + cl_wt * 81 + eta_cl_s3)")
  withr::local_dir(tempdir())
  p <- P(tr, method = "explicit", compile = FALSE, modelname = "re_keep")
  expect_identical(randomEffects(p), randomEffects(tr))
})

test_that("misspecified effect formulas and late substitutions stop with a message", {
  expect_error(insert(define(NULL, "x <- x", x = "cl"), cl ~ 1 + (1 | ID)), "branched")
  expect_error(insert(.tr(), log(cl) ~ 1), "scale in the trafo first")
  expect_error(insert(.tr(), ka ~ kz), "\"ka <- kz\"")
  expect_error(insert(.tr(), k ~ 1 + (1 | ID)), "not in the trafo")
  expect_error(insert(.tr(), cl ~ 1 + (1 | site)), "not a column of the branch table")
  expect_error(insert(.tr(), cl ~ 1 + (1 | p | gr(ID, dist = spikeSlab())),
                      v ~ 1 + (1 | p | gr(ID, dist = spikeSlab()))), "not correlated")
  expect_error(insert(.tr(), cl ~ 1 + (1 | gr(ID, dist = "t"))), "normal\\(\\), spikeSlab\\(\\) or dirichletProcess\\(\\)")
  expect_error(spikeSlab(prob = 1.5), "in \\(0, 1\\)")
  expect_identical(dirichletProcess(concentration = 1.5)$concentration, 1.5)
  expect_error(dirichletProcess(concentration = -1), "positive")
  expect_error(spikeSlab(tau = -1), "positive")
  expect_error(randomEffects(.tr()), "no random effects")
  tr <- insert(.tr(), cl ~ 1 + (1 | ID))
  expect_error(insert(tr, "x <- exp(x)", x = "cl"), "has an effect formula")
  expect_error(insert(tr, cl ~ 1 + (1 | ID)), "already has an effect formula")
  expect_silent(insert(tr, "x <- exp(x)", x = "v"))
})

test_that("a cluster term shares the concentration and tau over the symbols of the term", {
  re <- randomEffects(insert(.tr(), cl + v ~ 1 + (1 | gr(ID, dist = dirichletProcess()))))
  expect_equal(unique(re$terms$dist), "dirichletProcess")
  expect_length(re$blocks, 1L)
  expect_true(all(is.na(re$terms$concentration)))
})

test_that("spikeSlab() takes a slab mean and a spike, -Inf included", {
  d <- spikeSlab(prob = 0.3, mean = -1, spike = -Inf)
  expect_equal(d[c("mean", "spike")], list(mean = -1, spike = -Inf))
  expect_equal(spikeSlab()[c("mean", "spike")], list(mean = 0, spike = 0))
  expect_output(print(d), "spikeSlab\\(prob = 0.3, tau = NULL, mean = -1, spike = -Inf\\)")
  expect_error(spikeSlab(mean = Inf), "finite")
  expect_error(spikeSlab(mean = c(0, 1)), "finite")
  expect_error(spikeSlab(spike = Inf), "-Inf")
  expect_error(spikeSlab(spike = NA_real_), "-Inf")
  re <- randomEffects(insert(.tr(), v ~ 1 + (1 | gr(ID, dist = spikeSlab(mean = 2, spike = -Inf)))))
  expect_equal(re$terms[c("mean", "spike")], data.frame(mean = 2, spike = -Inf))
  expect_output(print(re), "mean 2, spike -Inf")
  expect_equal(randomEffects(insert(.tr(), v ~ 1 + (1 | ID)))$terms[c("mean", "spike")],
               data.frame(mean = 0, spike = 0))
})

test_that("gr() without a grouping puts all conditions in one group", {
  tr <- insert(.tr(), ka + cl ~ 0 + (1 | gr(dist = spikeSlab(spike = -Inf))))
  expect_equal(as.character(tr$C2[c("ka", "cl")]), c("(eta_ka_all)", "(eta_cl_all)"))
  re <- randomEffects(tr)
  expect_equal(re$subjects$all, list(all = c("C1", "C2", "C3")))
  expect_equal(unique(re$terms$group), "all")
  expect_length(re$fixed, 0)
  # one spike-and-slab block, so the symbols share prob and tau
  expect_equal(lapply(re$blocks, function(b) re$terms$stem[b]), list(c("eta_ka", "eta_cl")))
  spikeParams <- c("ka", "cl")
  expect_identical(insert(.tr(), x ~ 0 + (1 | gr(dist = spikeSlab(spike = -Inf))), x = spikeParams), tr)
  cov <- .cov; names(cov)[1] <- "all"
  expect_error(insert(branch(define(NULL, "x <- x", x = "cl"), table = cov), cl ~ 1 + (1 | all)),
               "rename the column")
})

test_that("the base of an exponential is read off the equations", {
  eqs <- list(c(a = "exp(x1)", b = "10^(x2)", c = "exp10(x3) + 2.3^x4", d = "6^(k + x5 - z)"),
              c(a = "exp(x1)", e = "x6 * exp(x6)", f = "exp(2 * x7)", g = "kb^x8", h = "2^(k - x9)"))
  ids <- paste0("x", 1:9)
  expect_equal(dMod2:::.wrapBases(eqs, ids),
               setNames(c(exp(1), 10, 10, 2.3, NA, NA, NA, NA, NA), ids))
  expect_equal(dMod2:::.wrapBases(eqs, ids, summand = TRUE),
               setNames(c(exp(1), 10, 10, 2.3, 6, NA, NA, NA, NA), ids))
  expect_equal(dMod2:::.petab_detect_scales(eqs, ids[1:4]),
               c(x1 = "log", x2 = "log10", x3 = "log10", x4 = "lin"))
})
