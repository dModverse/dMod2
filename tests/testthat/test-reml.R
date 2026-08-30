# ============================================================================
# Restricted maximum likelihood for the error model.
#
# The claims are checked against the stationarity condition itself,
#
#   sum_i [ 1 - h_ii - r_i^2/sigma_i^2 ] d log sigma_i^2 / d phi = 0,
#
# not against another dMod code path. For one constant sigma per observable
# that is sigma^2 = RSS/(n - sum h), which is checked in the same block.
# ============================================================================

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("cppDE")
  testthat::skip_on_cran()
}


# Decay chain with a constant error model. The observation carries a scale `s`,
# so the error model reports it as one of its inner parameters; `fix.s` has the
# transformation fix it, which separates the estimated parameter count from the
# inner one.
.build_reml_chain <- function(bench, mn_suffix, fix.s = FALSE) {
  .dmod_with_fx_workdir({
    gfn_s <- Y(c(y = "s*A"), f = bench$xfn, attach.input = FALSE,
               modelname = paste0("fx_decay_obs_s_", mn_suffix), compile = TRUE)
    e_const <- Y(c(y = "sigma_y"), f = gfn_s, attach.input = FALSE,
                 condition = "C1",
                 modelname = paste0("fx_decay_err_", mn_suffix),
                 compile = TRUE)
    trafo <- eqnvec(A = "A", k = "k", sigma_y = "sigma_y",
                    s = if (fix.s) "1" else "s")
    pfn <- P(trafo, condition = "C1",
             modelname = paste0("fx_decay_p_reml_", mn_suffix), compile = TRUE)
    prd <- gfn_s * bench$xfn * pfn
  })
  list(prd = prd, e = e_const)
}


test_that("leverages are hat values: in [0, 1] and summing to the rank", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.1)
  data$C1$sigma <- NA_real_

  ec <- .build_reml_chain(bench, "lev")
  obj  <- normL2(data, ec$prd, errmodel = ec$e)
  pars <- c(bench$outerpars_id, sigma_y = 0.1, s = 1)

  lev <- remlLeverage(obj, pars)
  expect_equal(nrow(lev), nrow(data$C1))
  expect_true(all(lev$leverage >= 0 & lev$leverage <= 1))
  expect_equal(sum(lev$leverage), attr(lev, "rank"), tolerance = 1e-8)
  # s and A enter only as their product, so one direction is absent
  expect_equal(attr(lev, "rank"), 2L)
  expect_equal(as.numeric(attr(lev, "dof")),
               nrow(data$C1) - sum(lev$leverage), tolerance = 1e-8)
})


test_that("reml converges to the stationary point of the restricted likelihood", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.1)
  data$C1$sigma <- NA_real_

  ec <- .build_reml_chain(bench, "fix", fix.s = TRUE)
  obj  <- normL2(data, ec$prd, errmodel = ec$e)
  pars <- c(bench$outerpars_id, sigma_y = 0.1)

  out <- reml(obj, pars, tol = 1e-8)
  expect_true(out$converged)
  expect_setequal(out$errpars, "sigma_y")

  lev <- out$leverage
  expect_equal(sum(1 - lev$leverage - (lev$residual / lev$sigma)^2), 0,
               tolerance = 1e-5)
  # same condition read as the closed form for one constant sigma
  expect_equal(sum((lev$residual / lev$sigma)^2),
               nrow(data$C1) - sum(lev$leverage), tolerance = 1e-5)
  # value.plain is the likelihood at the returned point, value adds the penalty
  expect_equal(out$value.plain,
               normL2(data, ec$prd, errmodel = ec$e)(out$argument)$value,
               tolerance = 1e-10)
  expect_equal(out$value, out$value.plain + out$logdet, tolerance = 1e-12)
  # REML trades likelihood for an unbiased sigma, so it cannot beat ML
  ml <- trust(obj, pars, rinit = 0.1, rmax = 10, iterlim = 500)$value
  expect_gte(out$value.plain, ml - 1e-6)
  expect_gt(out$argument["sigma_y"], 0)
})


test_that("confint takes the finite-sample F threshold", {
  grid <- seq(-4, 4, by = 0.005)
  prof <- parframe(
    data.frame(value = grid^2, constraint = grid, stepsize = 0.005, gamma = 1,
               whichPar = "p1", data = grid^2, p1 = grid),
    parameters = "p1",
    metanames = c("value", "constraint", "stepsize", "gamma", "whichPar"),
    obj.attributes = "data")

  ci_chisq <- confint(prof, level = 0.95)
  expect_equal(ci_chisq$upper, sqrt(qchisq(0.95, 1)), tolerance = 1e-3)
  expect_equal(ci_chisq$lower, -sqrt(qchisq(0.95, 1)), tolerance = 1e-3)

  n <- 48; p <- 6
  thr <- n * log(1 + qf(0.95, 1, n - p) / (n - p))
  ci_F <- confint(prof, level = 0.95, method = "F", n = n, p = p)
  expect_equal(ci_F$upper, sqrt(thr), tolerance = 1e-3)
  expect_gt(ci_F$upper, ci_chisq$upper)

  expect_error(confint(prof, method = "F"), "needs n and p")
  expect_error(confint(prof, method = "F", n = 4, p = 6), "must be positive")
})


test_that("reml sees every term of a split objective", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.1)
  data$C1$sigma <- NA_real_

  ec <- .build_reml_chain(bench, "split", fix.s = TRUE)
  pars <- c(bench$outerpars_id, sigma_y = 0.1)

  n <- nrow(data$C1)
  d_a <- data; d_a$C1 <- data$C1[1:6, , drop = FALSE]
  d_b <- data; d_b$C1 <- data$C1[7:n, , drop = FALSE]
  tt <- data$C1$time

  pooled <- normL2(data, ec$prd, errmodel = ec$e, times = tt)
  split  <- normL2(d_a, ec$prd, errmodel = ec$e, times = tt) +
            normL2(d_b, ec$prd, errmodel = ec$e, times = tt)

  # the split objective carries both terms, not just the first
  expect_length(attr(split, "l2spec"), 2L)
  expect_equal(nrow(remlLeverage(split, pars)), n)
  expect_equal(sum(remlLeverage(split, pars)$leverage),
               sum(remlLeverage(pooled, pars)$leverage), tolerance = 1e-8)

  out_pooled <- reml(pooled, pars, tol = 1e-8)
  out_split  <- reml(split,  pars, tol = 1e-8)
  expect_equal(out_split$argument[["sigma_y"]],
               out_pooled$argument[["sigma_y"]], tolerance = 1e-5)
})


test_that("profileThreshold matches its definition", {
  expect_equal(profileThreshold(0.95), qchisq(0.95, 1))
  n <- 48; p <- 6
  expect_equal(profileThreshold(0.95, method = "F", n = n, p = p),
               n * log(1 + qf(0.95, 1, n - p) / (n - p)))
  expect_gt(profileThreshold(0.95, method = "F", n = n, p = p),
            profileThreshold(0.95))
  expect_error(profileThreshold(method = "F"), "needs n and p")
  expect_error(profileThreshold(method = "F", n = 4, p = 6), "must be positive")
})


test_that("profile runs to the threshold it is given", {
  nm <- c("a", "b")
  obj <- function(pars, fixed = NULL, deriv = TRUE, ...) {
    pp <- c(pars, fixed)[nm]
    g <- c(a = 2 * pp[["a"]], b = 2 * pp[["b"]])
    H <- matrix(c(2, 0, 0, 2), 2, 2, dimnames = list(nm, nm))
    free <- names(pars)
    structure(list(value = as.numeric(pp[["a"]]^2 + pp[["b"]]^2),
                   gradient = g[free], hessian = H[free, free, drop = FALSE]),
              class = c("objlist", "list"))
  }
  class(obj) <- c("objfn", "fn")
  attr(obj, "modelname") <- character(0)

  d_chisq <- profileThreshold(0.95)
  d_F     <- profileThreshold(0.95, method = "F", n = 48, p = 6)

  p1 <- profile(obj, c(a = 0, b = 0), "a", limits = c(-5, 5), cores = 1)
  p2 <- profile(obj, c(a = 0, b = 0), "a", limits = c(-5, 5), cores = 1,
                delta = d_F)

  expect_gte(max(p1$value), d_chisq)
  expect_gte(max(p2$value), d_F)
  # the profile of a is exactly a^2, so it turns around just past sqrt(delta)
  expect_gt(max(abs(p2$constraint)), max(abs(p1$constraint)))
  expect_gte(max(abs(p2$constraint)), sqrt(d_F))
  expect_lt(max(abs(p2$constraint)), 2 * sqrt(d_F))
})


test_that("plotProfile draws the threshold it is given", {
  nm <- c("a", "b")
  obj <- function(pars, fixed = NULL, deriv = TRUE, ...) {
    pp <- c(pars, fixed)[nm]
    g <- c(a = 2 * pp[["a"]], b = 2 * pp[["b"]])
    H <- matrix(c(2, 0, 0, 2), 2, 2, dimnames = list(nm, nm))
    free <- names(pars)
    structure(list(value = as.numeric(pp[["a"]]^2 + pp[["b"]]^2),
                   gradient = g[free], hessian = H[free, free, drop = FALSE]),
              class = c("objlist", "list"))
  }
  class(obj) <- c("objfn", "fn")
  attr(obj, "modelname") <- character(0)

  prof <- profile(obj, c(a = 0, b = 0), "a", limits = c(-5, 5), cores = 1)

  p <- plotProfile(prof)
  expect_s3_class(p, "ggplot")
  # the level of the fit is always the first line, whatever the thresholds are
  expect_equal(unname(ggplot2::layer_scales(p)$y$breaks),
               c(0, 1, qchisq(0.90, 1), qchisq(0.95, 1)))
  expect_equal(ggplot2::layer_scales(p)$y$get_labels()[1], "optimum")

  thr <- profileThreshold(0.95, method = "F", n = 48, p = 6)
  p2 <- plotProfile(prof, threshold = c("95% / F" = thr))
  expect_equal(unname(ggplot2::layer_scales(p2)$y$breaks), c(0, thr))
  expect_equal(ggplot2::layer_scales(p2)$y$get_labels(), c("optimum", "95% / F"))
  expect_equal(unname(unique(ggplot2::layer_data(p2, 1)$yintercept)), c(0, thr))
})


test_that("the profile legend names sets and disappears for a single one", {
  nm <- c("a", "b")
  obj <- function(pars, fixed = NULL, deriv = TRUE, ...) {
    pp <- c(pars, fixed)[nm]
    g <- c(a = 2 * pp[["a"]], b = 2 * pp[["b"]])
    H <- matrix(c(2, 0, 0, 2), 2, 2, dimnames = list(nm, nm))
    free <- names(pars)
    structure(list(value = as.numeric(pp[["a"]]^2 + pp[["b"]]^2),
                   gradient = g[free], hessian = H[free, free, drop = FALSE]),
              class = c("objlist", "list"))
  }
  class(obj) <- c("objfn", "fn")
  attr(obj, "modelname") <- character(0)

  prof <- profile(obj, c(a = 0, b = 0), "a", limits = c(-3, 3), cores = 1)

  # what the reader sees, not what the object stores
  drawn <- function(p) {
    gt <- ggplot2::ggplotGrob(p)
    out <- character(0)
    walk <- function(g) {
      if (inherits(g, "titleGrob") || inherits(g, "text"))
        out <<- c(out, as.character(g$label))
      if (!is.null(g$children)) lapply(g$children, walk)
      if (!is.null(g$grobs))    lapply(g$grobs, walk)
      invisible(NULL)
    }
    walk(gt)
    unique(out)
  }

  one <- drawn(plotProfile(prof))
  expect_false("proflist" %in% one)
  expect_false("set" %in% one)
  expect_true("mode" %in% one)

  two <- drawn(plotProfile(list(A = prof, B = prof)))
  expect_false("proflist" %in% two)
  expect_true("set" %in% two)
  expect_true(all(c("A", "B") %in% two))
})
