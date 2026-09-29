# Multiple shooting: normL2(shooting = ) and trust() on the objective it returns.
#
# The oracle is single shooting. With nodes read off one continuous trajectory
# the gaps vanish, the data term is the single-shooting value, and the condensed
# gradient and Gauss-Newton Hessian are the single-shooting ones. The segment
# blocks are checked against finite differences, the reverse sweep against the
# forward one, and a converged fit against the single-shooting fit.

skip_on_cran()

.ms_dir <- function() {
  d <- file.path(tempdir(), "dmod_ms")
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  d
}
.ms_tol <- list(atol = 1e-10, rtol = 1e-10)

# FitzHugh-Nagumo, V observed, and a positive two-state chain for log charts.
# One shared object for everything the file evaluates.
.ms_models <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    d <- .ms_dir()
    owd <- setwd(d); on.exit(setwd(owd))
    fr <- c("forward", "reverse")
    f <- eqnvec(V = "c*(V - V^3/3 + R)", R = "-(V - a + b*R)/c")
    m <- odemodel(f, modelname = "ms_fhn", derivMode = fr,
                  includeTimeZero = FALSE, compile = FALSE, outdir = d)
    m0 <- odemodel(f, modelname = "ms_fhn0", compile = FALSE, outdir = d)
    x <- Xs(m, optionsOde = .ms_tol, optionsSens = .ms_tol)
    x0 <- Xs(m0, optionsOde = .ms_tol, optionsSens = .ms_tol)
    g <- Y(c(y = "V"), f = x, modelname = "ms_fhn_obs", compile = FALSE,
           derivMode = fr)
    p <- P(eqnvec(V = "V", R = "R", a = "a", b = "b", c = "exp(lc)"),
           condition = "C", modelname = "ms_fhn_p", compile = FALSE,
           derivMode = fr)

    re <- addReaction(addReaction(eqnlist(), "A", "B", "k1*A", "conversion"),
                      "B", "", "k2*B", "decay")
    mc <- odemodel(re, modelname = "ms_chain", derivMode = fr,
                   includeTimeZero = FALSE, compile = FALSE, outdir = d)
    xc <- Xs(mc, optionsOde = .ms_tol, optionsSens = .ms_tol)
    gc <- Y(c(yA = "s*A", yB = "s*B"), f = xc, modelname = "ms_chain_obs",
            compile = FALSE, derivMode = fr)
    pc <- P(eqnvec(A = "exp(lA)", B = "0", k1 = "exp(lk1)", k2 = "exp(lk2)",
                   s = "exp(ls)"),
            condition = "C", modelname = "ms_chain_p", compile = FALSE,
            derivMode = fr)
    # a second condition with an initial amount of its own
    pc2 <- P(eqnvec(A = "exp(lA2)", B = "0", k1 = "exp(lk1)", k2 = "exp(lk2)",
                    s = "exp(ls)"),
             condition = "D", modelname = "ms_chain_p2", compile = FALSE,
             derivMode = fr)
    # an error model with an estimated standard deviation
    ec <- Y(c(yA = "sd", yB = "sd"), f = gc, states = c("yA", "yB"),
            parameters = "sd", modelname = "ms_chain_err", compile = FALSE,
            derivMode = fr)
    pce <- P(eqnvec(A = "exp(lA)", B = "0", k1 = "exp(lk1)", k2 = "exp(lk2)",
                    s = "exp(ls)", sd = "exp(lsd)"),
             condition = "C", modelname = "ms_chain_pe", compile = FALSE,
             derivMode = fr)
    # a dose on a node and one inside a segment
    evd <- eventlist(var = c("A", "A"), time = c(4, 7.5), value = c("dose", "dose"),
                     method = c("add", "add"))
    md <- odemodel(re, modelname = "ms_dose", events = evd, derivMode = fr,
                   includeTimeZero = FALSE, compile = FALSE, outdir = d)
    xd <- Xs(md, optionsOde = .ms_tol, optionsSens = .ms_tol)
    gd <- Y(c(yA = "s*A", yB = "s*B"), f = xd, modelname = "ms_dose_obs",
            compile = FALSE, derivMode = fr)
    pd <- P(eqnvec(A = "exp(lA)", B = "0", k1 = "exp(lk1)", k2 = "exp(lk2)",
                   s = "exp(ls)", dose = "exp(ldose)"),
            condition = "C", modelname = "ms_dose_p", compile = FALSE,
            derivMode = fr)
    suppressMessages(compile(x, x0, g, p, xc, gc, pc, pc2, ec, pce, xd, gd, pd,
                             output = "ms_all", cores = 4))

    truth <- c(V = -1, R = 1, a = 0.2, b = 0.2, lc = log(3))
    tt <- seq(0, 20, by = 0.25)
    set.seed(1)
    y <- (g * x * p)(tt, truth, deriv = FALSE)[[1]][, "y"]
    data <- as.datalist(data.frame(name = "y", time = tt,
                                   value = y + rnorm(length(tt), 0, 0.1),
                                   sigma = 0.1, condition = "C"))

    truthc <- c(lA = 0, lk1 = log(0.5), lk2 = log(0.2), ls = 0)
    tc <- seq(0, 10, by = 0.5)
    pr <- (gc * xc * pc)(tc, truthc, deriv = FALSE)[[1]]
    datac <- as.datalist(rbind(
      data.frame(name = "yA", time = tc, value = pr[, "yA"] * exp(rnorm(length(tc), 0, 0.02)),
                 sigma = 0.02, condition = "C"),
      data.frame(name = "yB", time = tc[-1], value = pr[-1, "yB"] * exp(rnorm(length(tc) - 1, 0, 0.02)),
                 sigma = 0.02, condition = "C")))
    truth2 <- c(truthc, lA2 = log(2))
    pr2 <- (gc * xc * pc2)(tc, truth2, deriv = FALSE)[[1]]
    datac2 <- datac + as.datalist(rbind(
      data.frame(name = "yA", time = tc, value = pr2[, "yA"] * exp(rnorm(length(tc), 0, 0.02)),
                 sigma = 0.02, condition = "D"),
      data.frame(name = "yB", time = tc[-1], value = pr2[-1, "yB"] * exp(rnorm(length(tc) - 1, 0, 0.02)),
                 sigma = 0.02, condition = "D")))
    datace <- datac
    datace$C$sigma <- NA
    truthd <- c(truthc, ldose = 0)
    prd_d <- (gd * xd * pd)(tc, truthd, deriv = FALSE)[[1]]
    datad <- as.datalist(rbind(
      data.frame(name = "yA", time = tc, value = prd_d[, "yA"] * exp(rnorm(length(tc), 0, 0.02)),
                 sigma = 0.02, condition = "C"),
      data.frame(name = "yB", time = tc[-1], value = prd_d[-1, "yB"] * exp(rnorm(length(tc) - 1, 0, 0.02)),
                 sigma = 0.02, condition = "C")))
    cache <<- list(xd = xd, gd = gd, pd = pd, datad = datad, truthd = truthd,
                   x = x, x0 = x0, g = g, p = p, data = data, truth = truth,
                   xc = xc, gc = gc, pc = pc, datac = datac, truthc = truthc,
                   pc2 = pc2, datac2 = datac2, truth2 = truth2,
                   ec = ec, pce = pce, datace = datace)
    cache
  }
})

# The data term and the gap vector of one evaluation, as a function of the
# local variables of segment k, for numDeriv.
.ms_local_fn <- function(sobj, th, nd, k) {
  spec <- attr(sobj, "spec")
  seg <- spec$segs[[k]]
  function(v) {
    thv <- th; ndv <- nd
    for (nm in names(v)) {
      if (nm %in% names(thv)) thv[nm] <- v[nm]
      else ndv[[seg$cond]][seg$j - 1L, sub("@.*", "", nm)] <- v[nm]
    }
    E <- sobj(thv, ndv, deriv = FALSE)
    c(E$value, if (!seg$last) E$gaps[[seg$cond]][seg$j, spec$states])
  }
}


## ---- Construction ----

test_that("normL2 with multipleShootingControl evaluates as single shooting", {
  mo <- .ms_models()
  obj0 <- normL2(mo$data, mo$g * mo$x * mo$p)
  obj <- normL2(mo$data, mo$g * mo$x * mo$p,
                multipleShootingControl = list(nodes = seq(2, 18, 2)))
  expect_s3_class(obj, "objfn")
  expect_equal(obj(mo$truth)$value, obj0(mo$truth)$value)
  expect_equal(obj(mo$truth)$gradient, obj0(mo$truth)$gradient)
  expect_length(dMod2:::.shootingTerms(obj), 1L)
  expect_length(dMod2:::.shootingTerms(obj0), 0L)
  sm <- obj + constraintL2(mo$truth, 10)
  expect_length(attr(sm, "terms"), 2L)
  expect_length(dMod2:::.shootingTerms(sm), 1L)
})

test_that("the control is read and written through controls(), also on a sum", {
  mo <- .ms_models()
  obj <- normL2(mo$data, mo$g * mo$x * mo$p, multipleShootingControl = TRUE)
  expect_true(isTRUE(controls(obj, "multipleShootingControl")))
  expect_true("multipleShootingControl" %in% controls(obj))
  sm <- obj + constraintL2(mo$truth, 10)
  controls(sm, "multipleShootingControl") <- list(nodes = c(5, 10))
  # the summand is the objective's own closure, so the change reaches both
  expect_equal(controls(obj, "multipleShootingControl"), list(nodes = c(5, 10)))
  expect_equal(unname(attr(dMod2:::.shootObjOf(sm), "nodes")$C), c(0, 5, 10))
  controls(sm, "multipleShootingControl") <- NULL
  expect_true("multipleShootingControl" %in% controls(obj))
  expect_length(dMod2:::.shootingTerms(sm), 0L)
  expect_error(controls(sm, "nonsense") <- 1, "no control")
})

test_that("the segment layout rejects nodes outside the data", {
  mo <- .ms_models()
  expect_error(normL2(mo$data, mo$g * mo$x * mo$p, multipleShootingControl = list(nodes = c(5, 20))),
               "at or after its last data")
  expect_error(normL2(mo$data, mo$g * mo$x * mo$p, multipleShootingControl = list(nodes = 5, foo = 1)),
               "not foo")
  expect_error(normL2(mo$data, mo$g * mo$x * mo$p, multipleShootingControl = list(growth = 0.5)),
               "above 1")
  expect_error(normL2(mo$data, mo$g * mo$x * mo$p,
                      multipleShootingControl = list(init = list(C = matrix(0, 1, 2)))),
               "need the node times")
  expect_error(normL2(mo$data, mo$g * mo$x * mo$p,
                      multipleShootingControl = list(nodes = 5, charts = c(Q = "log10"))),
               "not states")
})

test_that("an end state below zero in a log10 chart is a segment failure", {
  spec <- list(states = "A", charts = c(A = "log10"),
               segs = list(list(j = 1L, cond = "C", end = 2, last = FALSE),
                           list(j = 2L, cond = "C", end = 4, last = TRUE)))
  pr <- function(a) matrix(c(0, 2, 1, a), 2, dimnames = list(NULL, c("time", "A")))
  err <- tryCatch(dMod2:::.shootEnds(spec, list(pr(-1e-9), NULL), FALSE),
                  error = function(e) e)
  expect_s3_class(err, "shootingSegmentError")
  expect_identical(err$segment, 1L)
  expect_match(conditionMessage(err), "not positive at the end of segment 1")
  ends <- dMod2:::.shootEnds(spec, list(pr(0.01), NULL), FALSE)
  expect_equal(unname(ends[[1]]$value), -2)
})

test_that("a model that forces t = 0 into the grid is refused at evaluation", {
  mo <- .ms_models()
  obj <- normL2(mo$data, mo$g * mo$x0 * mo$p, multipleShootingControl = list(nodes = seq(2, 18, 2)))
  sobj <- dMod2:::.shootObjOf(obj)
  nd <- dMod2:::.shootNodes(sobj, mo$truth)
  expect_error(sobj(mo$truth, nd), "includeTimeZero = FALSE")
})


test_that("multipleShootingControl = TRUE lays out ten segments that adapt", {
  mo <- .ms_models()
  obj <- normL2(mo$data, mo$g * mo$x * mo$p, multipleShootingControl = TRUE)
  sobj <- dMod2:::.shootObjOf(obj)
  expect_true(attr(sobj, "adaptive"))
  expect_equal(unname(attr(sobj, "nodes")$C), seq(0, 18, by = 2))
  expect_equal(attr(sobj, "growth"), 10)
  expect_equal(attr(sobj, "minLength")$C, 0.125)
})

test_that("a new layout keeps the data term and the continuity", {
  mo <- .ms_models()
  obj0 <- normL2(mo$data, mo$g * mo$x * mo$p)
  obj <- normL2(mo$data, mo$g * mo$x * mo$p, multipleShootingControl = TRUE)
  sobj <- dMod2:::.shootObjOf(obj)
  spec <- attr(sobj, "spec")
  nd <- dMod2:::.shootNodes(sobj, mo$truth, observed = FALSE)
  k <- 4L
  tn <- dMod2:::.shootSplitTime(spec$segs[[k]], 0.125)
  expect_true(tn > spec$segs[[k]]$start && tn < spec$segs[[k]]$end)
  v <- dMod2:::.shootValuesAt(sobj, mo$truth, nd, NULL, k, tn, 1L)[[1L]]
  tau <- spec$tau
  tau$C <- sort(c(tau$C, tn))
  fill <- list(C = setNames(list(v), as.character(tn)))
  rl <- dMod2:::.shootRelayout(sobj, tau, nd, fill)
  expect_equal(length(attr(rl$sobj, "spec")$segs), length(spec$segs) + 1L)
  E <- rl$sobj(mo$truth, rl$nodes, deriv = FALSE)
  expect_lt(max(abs(unlist(E$gaps))), 1e-6)
  expect_equal(E$value, obj0(mo$truth, deriv = FALSE)$value, tolerance = 1e-5)
  expect_true(attr(rl$sobj, "adaptive"))
})

test_that("the growth of every segment is its propagation matrix, forward and reverse", {
  mo <- .ms_models()
  obj <- normL2(mo$data, mo$g * mo$x * mo$p, multipleShootingControl = TRUE)
  sobj <- dMod2:::.shootObjOf(obj)
  nd <- dMod2:::.shootNodes(sobj, mo$truth)
  Ef <- sobj(mo$truth, nd, growth = TRUE)
  Er <- sobj(mo$truth, nd, sweep = "reverse")
  for (k in 1:9) {
    expect_equal(dim(Ef$segments[[k]]$growth), c(2L, 2L))
    expect_equal(unname(Er$segments[[k]]$growth), unname(Ef$segments[[k]]$growth),
                 tolerance = 1e-5, info = paste("segment", k))
  }
  # no stray directions of the first segment's start among its variables
  expect_false(any(grepl("@x0", Ef$segments[[1L]]$vars)))
})


## ---- Single shooting as the oracle ----

test_that("continuous nodes reproduce the single-shooting value, gradient and Hessian", {
  mo <- .ms_models()
  obj0 <- normL2(mo$data, mo$g * mo$x * mo$p)
  obj <- normL2(mo$data, mo$g * mo$x * mo$p, multipleShootingControl = list(nodes = seq(2, 18, 2)))
  sobj <- dMod2:::.shootObjOf(obj)
  spec <- attr(sobj, "spec")
  th <- mo$truth + c(0.05, -0.05, 0.02, 0.02, -0.02)
  nd <- dMod2:::.shootNodes(sobj, th, observed = FALSE)
  E <- sobj(th, nd)
  ss <- obj0(th)
  expect_lt(max(abs(unlist(E$gaps))), 1e-6)
  expect_equal(E$value, ss$value, tolerance = 1e-5)
  C <- dMod2:::.shootCondense(E, lapply(E$segments, `[[`, "hess"), spec,
                              names(th), c(V = 1, R = 1))
  expect_equal(unname(C$g0), unname(ss$gradient[names(th)]), tolerance = 1e-5)
  expect_equal(unname(C$H), unname(ss$hessian[names(th), names(th)]), tolerance = 1e-6)
})

test_that("segment gradients and end-state Jacobians match finite differences", {
  skip_if_not_installed("numDeriv")
  mo <- .ms_models()
  obj <- normL2(mo$data, mo$g * mo$x * mo$p, multipleShootingControl = list(nodes = seq(2, 18, 2)))
  sobj <- dMod2:::.shootObjOf(obj)
  th <- mo$truth + c(0.1, -0.1, 0.05, 0.03, -0.05)
  nd <- dMod2:::.shootNodes(sobj, mo$truth)
  E <- sobj(th, nd)
  k <- 3L
  sg <- E$segments[[k]]
  loc <- sg$vars
  nv <- nd$C[2L, sub("@.*", "", grep("@", loc, value = TRUE))]
  names(nv) <- grep("@", loc, value = TRUE)
  v0 <- c(th, nv)[loc]
  J <- numDeriv::jacobian(.ms_local_fn(sobj, th, nd, k), v0)
  # the data term of all segments depends on theta, the one of segment k alone
  # on its node
  gtot <- Reduce(`+`, lapply(E$segments, function(s) s$grad[names(th)]))
  it <- match(names(th), loc)
  inode <- setdiff(seq_along(loc), it)
  expect_equal(J[1L, it], unname(gtot), tolerance = 1e-5)
  expect_equal(J[1L, inode], unname(sg$grad[loc[inode]]), tolerance = 1e-5)
  expect_equal(J[-1L, ], unname(sg$jac[, loc]), tolerance = 1e-5)
})

test_that("log10 charts carry the chain rule on the nodes and the ends", {
  skip_if_not_installed("numDeriv")
  mo <- .ms_models()
  obj <- normL2(mo$datac, mo$gc * mo$xc * mo$pc,
                multipleShootingControl = list(nodes = c(3, 6), charts = c(A = "log10", B = "log10")))
  sobj <- dMod2:::.shootObjOf(obj)
  th <- mo$truthc + 0.1
  nd <- dMod2:::.shootNodes(sobj, mo$truthc)
  expect_true(all(nd$C[, "A"] < 0))       # log10 of values below one
  E <- sobj(th, nd)
  k <- 2L
  sg <- E$segments[[k]]
  loc <- sg$vars
  nv <- nd$C[1L, sub("@.*", "", grep("@", loc, value = TRUE))]
  names(nv) <- grep("@", loc, value = TRUE)
  v0 <- c(th, nv)[loc]
  J <- numDeriv::jacobian(.ms_local_fn(sobj, th, nd, k), v0)
  inode <- grep("@", loc)
  expect_equal(J[-1L, ], unname(sg$jac[, loc]), tolerance = 1e-5)
  expect_equal(J[1L, inode], unname(sg$grad[loc[inode]]), tolerance = 1e-5)
})

test_that("the reverse sweep reproduces the forward blocks", {
  mo <- .ms_models()
  obj <- normL2(mo$data, mo$g * mo$x * mo$p, multipleShootingControl = list(nodes = seq(2, 18, 2)))
  sobj <- dMod2:::.shootObjOf(obj)
  th <- mo$truth + c(0.1, -0.1, 0.05, 0.03, -0.05)
  nd <- dMod2:::.shootNodes(sobj, mo$truth)
  Ef <- sobj(th, nd)
  Er <- sobj(th, nd, sweep = "reverse")
  expect_equal(Er$value, Ef$value, tolerance = 1e-6)
  for (k in seq_along(Ef$segments)) {
    v <- Ef$segments[[k]]$vars
    expect_equal(Er$segments[[k]]$grad[v], Ef$segments[[k]]$grad,
                 tolerance = 1e-5, info = paste("segment", k))
    if (!is.null(Ef$segments[[k]]$jac))
      expect_equal(Er$segments[[k]]$jac[, v], Ef$segments[[k]]$jac,
                   tolerance = 1e-5, info = paste("segment", k))
  }
})


test_that("an estimated error model differentiates through the segments", {
  skip_if_not_installed("numDeriv")
  mo <- .ms_models()
  obj <- normL2(mo$datace, mo$gc * mo$xc * mo$pce, errmodel = mo$ec,
                multipleShootingControl = list(nodes = c(3, 6), charts = c(A = "log10", B = "log10")))
  sobj <- dMod2:::.shootObjOf(obj)
  th <- c(mo$truthc, lsd = log(0.03)) + 0.05
  nd <- dMod2:::.shootNodes(sobj, th)
  E <- sobj(th, nd)
  Er <- sobj(th, nd, sweep = "reverse")
  f1 <- function(v) sobj(v, nd, deriv = FALSE)$value
  gfd <- numDeriv::grad(f1, th)
  gan <- Reduce(`+`, lapply(E$segments, function(s) s$grad[names(th)]))
  expect_equal(unname(gan), gfd, tolerance = 1e-5)
  grv <- Reduce(`+`, lapply(Er$segments, function(s) s$grad[names(th)]))
  expect_equal(unname(grv), gfd, tolerance = 1e-5)
})

test_that("several conditions share theta and keep nodes of their own", {
  mo <- .ms_models()
  prd <- mo$gc * mo$xc * (mo$pc + mo$pc2)
  obj0 <- normL2(mo$datac2, prd)
  obj <- normL2(mo$datac2, prd,
                multipleShootingControl = list(nodes = list(C = c(3, 6), D = c(2, 4, 8)),
                                charts = c(A = "log10", B = "log10")))
  sobj <- dMod2:::.shootObjOf(obj)
  expect_equal(lengths(attr(sobj, "nodes")), c(C = 3L, D = 4L))
  ref <- trust(obj0, mo$truth2, iterlim = 200)
  fit <- trust(obj, mo$truth2 + 0.2, iterlim = 200)
  expect_true(fit$converged)
  expect_equal(fit$value, ref$value, tolerance = 1e-6)
  # y = s*A sees the scale and the initial amounts only as products
  ident <- function(a) c(a[["lA"]] + a[["ls"]], a[["lA2"]] + a[["ls"]],
                         a[["lk1"]], a[["lk2"]])
  expect_equal(ident(fit$argument), ident(ref$argument), tolerance = 1e-4)
  expect_equal(nrow(fit$multipleShooting$nodes$D), 3L)
})


test_that("fixed-time events fire once, on a node and inside a segment", {
  mo <- .ms_models()
  prd <- mo$gd * mo$xd * mo$pd
  obj0 <- normL2(mo$datad, prd)
  obj <- normL2(mo$datad, prd,
                multipleShootingControl = list(nodes = c(2, 4, 6, 8),
                                               charts = c(A = "log10", B = "log10")))
  sobj <- dMod2:::.shootObjOf(obj)
  th <- mo$truthd + 0.05
  nd <- dMod2:::.shootNodes(sobj, th, observed = FALSE)
  E <- sobj(th, nd)
  # continuous nodes: no gap at the node that carries the dose, and the data
  # term of the single solve
  expect_lt(max(abs(unlist(E$gaps))), 1e-6)
  expect_equal(E$value, obj0(th, deriv = FALSE)$value, tolerance = 1e-6)
  ref <- trust(obj0, mo$truthd, iterlim = 200)
  fit <- trust(obj, mo$truthd + 0.2, iterlim = 200)
  expect_equal(fit$value, ref$value, tolerance = 1e-6)
})


## ---- trust() ----

test_that("trust() converges to the single-shooting optimum", {
  mo <- .ms_models()
  obj0 <- normL2(mo$data, mo$g * mo$x * mo$p)
  obj <- normL2(mo$data, mo$g * mo$x * mo$p, multipleShootingControl = list(nodes = seq(1, 19, 1)))
  ref <- trust(obj0, mo$truth, iterlim = 200)
  start <- mo$truth + c(0.8, -0.6, 0.4, -0.3, 0.3)
  for (acc in c("filter", "merit")) {
    fit <- trust(obj, start, iterlim = 200, stepControl = list(acceptance = acc))
    expect_true(fit$converged, info = acc)
    expect_lt(fit$multipleShooting$violation, 1e-6)
    expect_equal(fit$value, ref$value, tolerance = 1e-6, info = acc)
    expect_equal(fit$argument, ref$argument, tolerance = 1e-4, info = acc)
    # at a continuous trajectory the condensed gradient is the single-shooting one
    expect_lt(max(abs(fit$gradient - obj0(fit$argument)$gradient[names(start)])), 1e-3)
  }
})

test_that("automatic nodes converge to the single-shooting optimum", {
  mo <- .ms_models()
  obj0 <- normL2(mo$data, mo$g * mo$x * mo$p)
  obj <- normL2(mo$data, mo$g * mo$x * mo$p, multipleShootingControl = TRUE)
  ref <- trust(obj0, mo$truth, iterlim = 200)
  start <- mo$truth + c(0.8, -0.6, 0.4, -0.3, 0.3)
  fit <- trust(obj, start, iterlim = 200)
  expect_true(fit$converged)
  expect_equal(fit$value, ref$value, tolerance = 1e-6)
  expect_gte(fit$multipleShooting$nSplit, 0L)
  expect_equal(nrow(fit$multipleShooting$nodes$C), length(fit$multipleShooting$gaps$C[, 1L]))
})

test_that("annealing the continuity reaches the same optimum", {
  mo <- .ms_models()
  obj0 <- normL2(mo$data, mo$g * mo$x * mo$p)
  obj <- normL2(mo$data, mo$g * mo$x * mo$p, multipleShootingControl = TRUE)
  ref <- trust(obj0, mo$truth, iterlim = 200)
  start <- mo$truth + c(0.8, -0.6, 0.4, -0.3, 0.3)
  fit <- trust(obj, start, iterlim = 200, stepControl = list(anneal = TRUE))
  expect_true(fit$converged)
  expect_gte(fit$multipleShooting$nStages, 1L)
  expect_equal(fit$value, ref$value, tolerance = 1e-6)
  expect_lt(fit$multipleShooting$violation, 1e-6)
})

test_that("the natural level function converges to the single-shooting optimum", {
  mo <- .ms_models()
  obj0 <- normL2(mo$data, mo$g * mo$x * mo$p)
  ref <- trust(obj0, mo$truth, iterlim = 200)
  start <- mo$truth + c(0.8, -0.6, 0.4, -0.3, 0.3)
  for (nodes in list(seq(1, 19, 1), "auto")) {
    obj <- normL2(mo$data, mo$g * mo$x * mo$p,
                  multipleShootingControl = list(nodes = nodes))
    fit <- trust(obj, start, iterlim = 200, stepControl = list(acceptance = "natural"))
    expect_true(fit$converged, info = paste(nodes, collapse = " "))
    expect_lt(fit$multipleShooting$violation, 1e-6)
    expect_equal(fit$value, ref$value, tolerance = 1e-6)
    expect_equal(fit$argument, ref$argument, tolerance = 1e-4)
  }
  # Bock's method is a Gauss-Newton method on residuals of a fixed sigma
  obj <- normL2(mo$data, mo$g * mo$x * mo$p, multipleShootingControl = TRUE)
  expect_error(trust(obj, start, hessianMethod = "bfgs",
                     stepControl = list(acceptance = "natural")), "hessianMethod")
})

test_that("the residuals of a segment make up its value, gradient and Hessian", {
  mo <- .ms_models()
  obj <- normL2(mo$data, mo$g * mo$x * mo$p, multipleShootingControl = list(nodes = c(5, 10, 15)))
  sobj <- dMod2:::.shootObjOf(obj)
  th <- mo$truth + 0.1
  nd <- dMod2:::.shootNodes(sobj, th, cores = 1L)
  E <- sobj(th, nd, deriv = TRUE, hessian = TRUE, residuals = TRUE)
  for (sg in E$segments) {
    expect_equal(sum(sg$res^2), sg$chi2, tolerance = 1e-10)
    expect_equal(2 * as.numeric(crossprod(sg$J, sg$res)), unname(sg$grad), tolerance = 1e-8)
    expect_equal(2 * crossprod(sg$J), sg$hess, tolerance = 1e-8, ignore_attr = TRUE)
  }
  Ev <- sobj(th, nd, deriv = FALSE, residuals = TRUE)
  expect_equal(Ev$segments[[2]]$res, E$segments[[2]]$res, tolerance = 1e-6)
})

test_that("two phases, a spline start and regularisation reach the optimum", {
  mo <- .ms_models()
  obj0 <- normL2(mo$data, mo$g * mo$x * mo$p)
  ref <- trust(obj0, mo$truth, iterlim = 200)
  start <- mo$truth + c(0.8, -0.6, 0.4, -0.3, 0.3)
  obj <- normL2(mo$data, mo$g * mo$x * mo$p, multipleShootingControl = TRUE)
  fit <- trust(obj, start, iterlim = 200, stepControl = list(twoPhase = TRUE))
  expect_true(fit$converged)
  expect_equal(fit$value, ref$value, tolerance = 1e-6)
  expect_identical(fit$multipleShooting$nRestore, 0L)
  fit <- trust(obj, start, iterlim = 200, stepControl = list(restore = TRUE))
  expect_true(fit$converged)
  expect_equal(fit$value, ref$value, tolerance = 1e-6)
  objS <- normL2(mo$data, mo$g * mo$x * mo$p,
                 multipleShootingControl = list(init = "spline"))
  fit <- trust(objS, start, iterlim = 200)
  expect_true(fit$converged)
  expect_equal(fit$value, ref$value, tolerance = 1e-6)
  for (acc in c("filter", "natural")) {
    fit <- trust(obj, start, iterlim = 200,
                 stepControl = list(acceptance = acc, regularise = 1e-6))
    expect_true(fit$converged, info = acc)
    expect_equal(fit$value, ref$value, tolerance = 1e-6, info = acc)
  }
})

test_that("a spline start reads the observed states off a smoothing spline", {
  mo <- .ms_models()
  obj <- normL2(mo$data, mo$g * mo$x * mo$p,
                multipleShootingControl = list(nodes = c(5, 10, 15), init = "spline"))
  sobj <- dMod2:::.shootObjOf(obj)
  nd <- dMod2:::.shootNodes(sobj, mo$truth, cores = 1L, smooth = TRUE)
  d <- mo$data$C
  sp <- stats::smooth.spline(d$time, d$value)
  expect_equal(unname(nd$C[, "V"]), stats::predict(sp, c(5, 10, 15))$y, tolerance = 1e-10)
  expect_error(normL2(mo$data, mo$g * mo$x * mo$p,
                      multipleShootingControl = list(init = "splines")), "spline")
})

test_that("minPoints keeps enough data on both sides of a cut", {
  seg <- list(start = 0, end = 10, times = seq(0, 10, 1),
              data = data.frame(time = c(1, 2, 3, 8, 9)))
  expect_equal(dMod2:::.shootSplitTime(seg, 0.5), 5)
  # the halves around 5 hold 3 and 2 data times; 3 each is too many
  expect_true(is.na(dMod2:::.shootSplitTime(seg, 0.5, minPoints = 3L)))
  # 2 each: a cut in the middle half that leaves both sides two data times
  tn <- dMod2:::.shootSplitTime(seg, 0.5, minPoints = 2L)
  expect_true(tn > 2.5 && tn < 7.5)
  expect_gte(sum(seg$data$time < tn), 2L); expect_gte(sum(seg$data$time >= tn), 2L)
  expect_error(normL2(.ms_models()$data, with(.ms_models(), g * x * p),
                      multipleShootingControl = list(minPoints = 0)), "minPoints")
})

test_that("a continuity break frees the node after it", {
  mo <- .ms_models()
  prd <- mo$g * mo$x * mo$p
  # the state is kicked at t = 10, which the model cannot do; data before and
  # after, noise-free
  t1 <- seq(0, 9.75, by = 0.25); t2 <- seq(0, 10, by = 0.25)
  xs <- (mo$x * mo$p)(c(0, 10), mo$truth, deriv = FALSE)[[1]]
  kick <- mo$truth
  kick["V"] <- xs[nrow(xs), "V"] + 1; kick["R"] <- xs[nrow(xs), "R"]
  d <- as.datalist(data.frame(name = "y", time = c(t1, 10 + t2),
    value = c(prd(t1, mo$truth, deriv = FALSE)[[1]][, "y"],
              prd(t2, kick, deriv = FALSE)[[1]][, "y"]),
    sigma = 0.1, condition = "C"))
  perfect <- length(c(t1, t2)) * log(2 * pi * 0.1^2)
  for (acc in c("filter", "natural")) {
    obj <- normL2(d, prd, multipleShootingControl = list(breaks = 10))
    fit <- trust(obj, mo$truth + 0.2, iterlim = 200, stepControl = list(acceptance = acc))
    expect_true(fit$converged, info = acc)
    expect_equal(fit$value, perfect, tolerance = 1e-6, info = acc)
    expect_equal(fit$argument, mo$truth, tolerance = 1e-6, info = acc)
    expect_equal(unname(fit$multipleShooting$nodes$C["10", "V"]), unname(kick["V"]),
                 tolerance = 1e-6, info = acc)
  }
  # without the break the kick cannot be fitted
  obj <- normL2(d, prd, multipleShootingControl = TRUE)
  fit <- trust(obj, mo$truth + 0.2, iterlim = 200)
  expect_gt(fit$value, perfect + 10)
  expect_error(normL2(d, prd, multipleShootingControl = list(nodes = c(5, 15), breaks = 10)),
               "node times")
  expect_error(trust(normL2(d, prd, multipleShootingControl = list(breaks = 10)), mo$truth,
                     hessianMethod = "bfgs"), "breaks")
})

test_that("auto lays out one segment per half oscillation of the data", {
  tt <- seq(0, 10, by = 0.05)
  set.seed(2)
  osc <- data.frame(name = "y", time = tt, value = sin(2 * pi * tt) + rnorm(length(tt), 0, 0.02),
                    sigma = 0.02)
  flat <- transform(osc, value = exp(-tt) + rnorm(length(tt), 0, 0.02))
  # ten periods reverse twenty times; noise below five sigma never counts
  expect_equal(dMod2:::.shootTurns(osc), 20L)
  expect_equal(dMod2:::.shootTurns(flat), 0L)
  d <- as.datalist(rbind(transform(osc, condition = "A"), transform(flat, condition = "B")))
  nd <- dMod2:::.shootAutoNodes(d, list(), NULL, oscillations = TRUE, minPoints = 3L)
  expect_equal(lengths(nd), c(A = 20L, B = 9L))
  # no more segments than minPoints data points each allow, and never fewer than ten
  nd <- dMod2:::.shootAutoNodes(d, list(), NULL, oscillations = TRUE, minPoints = 20L)
  expect_equal(lengths(nd), c(A = 9L, B = 9L))
})

test_that("minPoints counts data points and defaults to one more than the states", {
  mo <- .ms_models()
  sobj <- dMod2:::.shootObjOf(normL2(mo$data, mo$g * mo$x * mo$p,
                                     multipleShootingControl = TRUE))
  expect_equal(attr(sobj, "minPoints"), 3L)
  seg <- list(start = 0, end = 10, times = seq(0, 10, 1),
              data = data.frame(time = rep(c(1, 2, 8, 9), each = 2)))
  expect_equal(dMod2:::.shootSplitTime(seg, 0.5, minPoints = 4L), 5)
  expect_true(is.na(dMod2:::.shootSplitTime(seg, 0.5, minPoints = 5L)))
})

test_that("a misfit everywhere does not refine the layout", {
  mo <- .ms_models()
  obj <- normL2(mo$data, mo$g * mo$x * mo$p, multipleShootingControl = TRUE)
  # far from the truth every segment misses its data alike
  start <- mo$truth + c(1.5, -1.5, 1, 1, -1)
  fit <- trust(obj, start, iterlim = 0, stepControl = list(anneal = FALSE))
  expect_lt(nrow(fit$multipleShooting$nodes$C) + 1L, 20L)
  expect_named(fit$multipleShooting$cuts, c("growth", "misfit", "failure"))
})

test_that("an angle chart takes gaps modulo 2 pi", {
  spec <- list(conditions = "C", states = c("x", "phi"),
               charts = c(x = "linear", phi = "angle"),
               segs = rep(list(list(cond = "C")), 3L), tau = list(C = c(0, 1, 2)))
  ends <- list(list(value = c(x = 1, phi = 2 * pi + 0.1)),
               list(value = c(x = 2, phi = -0.2)), NULL)
  nodes <- list(C = matrix(c(1.5, 0, 2, 4 * pi), 2L, 2L, byrow = TRUE,
                           dimnames = list(NULL, c("x", "phi"))))
  G <- dMod2:::.shootGaps(spec, ends, nodes)$C
  expect_equal(unname(G[, "x"]), c(-0.5, 0))
  expect_equal(unname(G[, "phi"]), c(0.1, -0.2))
  mo <- .ms_models()
  expect_error(normL2(mo$data, mo$g * mo$x * mo$p,
                      multipleShootingControl = list(charts = c(R = "circle"))), "angle")
  sobj <- dMod2:::.shootObjOf(normL2(mo$data, mo$g * mo$x * mo$p,
                                     multipleShootingControl = list(charts = c(R = "angle"))))
  nd <- dMod2:::.shootNodes(sobj, mo$truth)
  expect_equal(unname(dMod2:::.shootScale(NULL, nd, attr(sobj, "spec"), sobj, mo$truth)["R"]), 1)
})

test_that("node times next to a grid time snap onto it", {
  g <- seq(0, 20, by = 0.1)
  expect_identical(dMod2:::.shootSnap(c(0.7, 0.33333), g), c(g[8L], 0.33333))
  mo <- .ms_models()
  obj <- normL2(mo$data, mo$g * mo$x * mo$p, times = g,
                multipleShootingControl = list(nodes = c(0.7, 5.3)))
  sobj <- dMod2:::.shootObjOf(obj)
  nd <- dMod2:::.shootNodes(sobj, mo$truth)
  expect_true(is.finite(sobj(mo$truth, nd, deriv = FALSE)$value))
})

test_that("transitions lays a node just before every fast change of the data", {
  tt <- 0:100
  # a slow rise that drops back at 30, 60 and 90
  d <- as.datalist(data.frame(name = "y", time = tt, value = (tt %% 30) / 30,
                              sigma = 1, condition = "C"))
  x <- structure(function(...) NULL, mappings = NULL)
  nd <- dMod2:::.shootTransitionNodes(d, x, numeric(0))$C
  expect_equal(nd, c(27, 57, 87))
  # no fast change: no layout of its own
  d$C$value <- sin(tt / 10)
  expect_null(dMod2:::.shootTransitionNodes(d, x, numeric(0)))
  # through normL2, which falls back to "auto" on FitzHugh-Nagumo
  mo <- .ms_models()
  obj <- normL2(mo$data, mo$g * mo$x * mo$p,
                multipleShootingControl = list(nodes = "transitions"))
  fit <- trust(obj, mo$truth + c(0.8, -0.6, 0.4, -0.3, 0.3), iterlim = 200)
  expect_true(fit$converged)
})

test_that("the sparse trust-region step solves the shifted system", {
  H <- Matrix::forceSymmetric(Matrix::Matrix(c(4, 1, 0, 1, 3, 0, 0, 0, 2), 3, sparse = TRUE))
  g <- c(1, -2, 0.5)
  big <- dMod2:::.shootSparseTR(H, g, 10)
  expect_equal(big$step, -as.numeric(solve(as.matrix(H), g)), tolerance = 1e-8)
  small <- dMod2:::.shootSparseTR(H, g, 0.1)
  expect_equal(sqrt(sum(small$step^2)), 0.1, tolerance = 1e-3)
  expect_equal(small$step, -as.numeric(solve(as.matrix(H) + diag(small$lambda, 3), g)),
               tolerance = 1e-6)
})

test_that("quasi-Newton blocks run forward and on the reverse sweep", {
  mo <- .ms_models()
  obj0 <- normL2(mo$data, mo$g * mo$x * mo$p)
  obj <- normL2(mo$data, mo$g * mo$x * mo$p, multipleShootingControl = list(nodes = seq(1, 19, 1)))
  ref <- trust(obj0, mo$truth, iterlim = 200)
  start <- mo$truth + c(0.3, -0.3, 0.2, -0.2, 0.1)
  fb <- trust(obj, start, iterlim = 300, hessianMethod = "bfgs")
  fs <- trust(obj, start, iterlim = 300, hessianMethod = "sr1", sweep = "reverse")
  expect_equal(fb$value, ref$value, tolerance = 1e-6)
  expect_equal(fs$value, ref$value, tolerance = 1e-6)
  expect_equal(fs$multipleShooting$sweep, "reverse")
})

test_that("a prior enters the condensed problem", {
  mo <- .ms_models()
  prior <- constraintL2(mo$truth, sigma = 0.5)
  obj0 <- normL2(mo$data, mo$g * mo$x * mo$p) + prior
  obj <- normL2(mo$data, mo$g * mo$x * mo$p, multipleShootingControl = list(nodes = seq(2, 18, 2))) + prior
  ref <- trust(obj0, mo$truth, iterlim = 200)
  fit <- trust(obj, mo$truth + 0.2, iterlim = 200)
  expect_equal(fit$value, ref$value, tolerance = 1e-6)
})

test_that("trust() refuses what multiple shooting cannot honour", {
  mo <- .ms_models()
  obj <- normL2(mo$data, mo$g * mo$x * mo$p, multipleShootingControl = list(nodes = seq(2, 18, 2)))
  expect_error(trust(obj, mo$truth, hessianMethod = "exact"), "does not support")
  expect_error(trust(obj, mo$truth, hessianFallback = "sr1"), "does not support")
  expect_error(trust(obj, mo$truth, stepControl = list(nonmonotone = 0.8)), "does not support")
  expect_error(trust(obj, mo$truth, sweep = "reverse"), "bfgs")
  expect_error(trust(obj, mo$truth, deriv2 = TRUE), "takes fixed, cores and sweep")
})

test_that("mstrust fits a multiple-shooting objective through trust", {
  mo <- .ms_models()
  obj <- normL2(mo$data, mo$g * mo$x * mo$p, multipleShootingControl = list(nodes = seq(2, 18, 2)))
  fits <- mstrust(obj, mo$truth, fits = 2, sd = 0.1, iterlim = 100)
  expect_length(fits, 2L)
  expect_true(all(vapply(fits, function(f) is.list(f$multipleShooting), TRUE)))
  expect_s3_class(as.parframe(fits), "parframe")
})

test_that("profile() re-optimises by single shooting", {
  mo <- .ms_models()
  obj <- normL2(mo$data, mo$g * mo$x * mo$p, multipleShootingControl = list(nodes = seq(2, 18, 2)))
  single <- dMod2:::.singleShooting(obj)
  expect_length(dMod2:::.shootingTerms(single), 0L)
  expect_s3_class(single, "objfn")
  expect_equal(single(mo$truth)$value, obj(mo$truth)$value)
  # the objective itself keeps its control
  expect_length(dMod2:::.shootingTerms(obj), 1L)
})


## ---- The pieces multiple shooting is built from ----

test_that("mstrust routes arguments by the formals of its optimiser", {
  seen <- new.env()
  myopt <- function(objfun, parinit, special = 0, ...) {
    seen$special <- special
    seen$dots <- names(list(...))
    list(argument = parinit, value = 0, converged = TRUE, iterations = 0L)
  }
  f <- function(pars, ...) list(value = 0, gradient = pars * 0,
                                hessian = diag(length(pars)))
  class(f) <- c("objfn", "fn")
  mstrust(f, c(a = 1), fits = 1, optmethod = myopt, special = 3)
  expect_equal(seen$special, 3)
  expect_false("special" %in% seen$dots)
})

test_that("the backward walk carries several first-order seeds in one sweep", {
  mo <- .ms_models()
  pin <- mo$p(mo$truth, deriv = FALSE)[["C"]]
  tt <- seq(0, 5, by = 0.5)
  b <- dMod2:::.bundle(conds = "C", times = list(tt), pars = list(pin),
                       fixed = list(NULL), shared = FALSE)
  fw <- dMod2:::.fwdMany(mo$x, b, new.env(), 1L)
  pr <- fw$values[[1L]]
  W <- array(0, c(nrow(pr), ncol(pr), 3L), dimnames = list(NULL, colnames(pr), NULL))
  W[nrow(pr), "V", 1L] <- 1
  W[nrow(pr), "R", 2L] <- 1
  W[5L, "V", 3L] <- 1
  u <- dMod2:::.bwdNode(fw$tape, list(dMod2:::.ct(out = W)), new.env(), 1L, seeds = TRUE)
  U <- u[[1L]]$pars
  expect_equal(ncol(U), 3L)
  # each seed alone, as a first-order sweep of its own, replaying the same store
  for (s in 1:3) {
    us <- dMod2:::.bwdNode(fw$tape, list(dMod2:::.ct(out = W[, , s, drop = FALSE])),
                           new.env(), 1L)
    expect_equal(U[rownames(us[[1L]]$pars), s], us[[1L]]$pars[, 1L], tolerance = 1e-7)
  }
})

test_that("trust_step_impl takes the Newton step when it fits the radius", {
  H <- matrix(c(4, 1, 1, 3), 2)
  g <- c(1, -2)
  st <- dMod2:::trust_step_impl(c(a = 0, b = 0), g, H, 10, rep(-Inf, 2),
                                rep(Inf, 2), c(1, 1), 0.99995)
  expect_equal(unname(st$step), unname(-solve(H, g)), tolerance = 1e-12)
  expect_equal(st$quad, sum(g * st$step) + 0.5 * sum(st$step * (H %*% st$step)))
  st2 <- dMod2:::trust_step_impl(c(a = 0, b = 0), g, H, 0.1, rep(-Inf, 2),
                                 rep(Inf, 2), c(1, 1), 0.99995)
  expect_equal(sqrt(sum(st2$step^2)), 0.1, tolerance = 1e-8)
})

test_that("qn_update_impl reproduces the textbook updates", {
  B <- diag(2); s <- c(1, 0.5); y <- c(2, 1.5)
  sr1 <- dMod2:::qn_update_impl(B, s, y, "sr1")$B
  w <- y - B %*% s
  expect_equal(sr1, B + w %*% t(w) / sum(w * s), ignore_attr = TRUE)
  bfgs <- dMod2:::qn_update_impl(B, s, y, "bfgs")$B
  expect_equal(bfgs %*% s, matrix(y), tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("the normL2 kernel returns every condition's own block", {
  mo <- .ms_models()
  obj <- normL2(mo$data, mo$g * mo$x * mo$p, multipleShootingControl = list(nodes = seq(2, 18, 2)))
  sobj <- dMod2:::.shootObjOf(obj)
  E <- sobj(mo$truth, dMod2:::.shootNodes(sobj, mo$truth))
  expect_true(all(vapply(E$segments, function(s) is.matrix(s$hess), TRUE)))
  expect_true(all(vapply(E$segments, function(s)
    identical(rownames(s$hess), s$vars), TRUE)))
})
