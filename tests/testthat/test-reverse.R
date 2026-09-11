# Reverse-mode derivatives through the chain: stage 7 of
# cppDE/dev/adjoint-plan.md.
#
# The oracle is the forward mode, and it is not sharp: a forward-sensitivity
# solve carries n_theta tangent columns and its error norm takes the maximum
# over all of them, so it steps finer than a value-only run does. The two modes
# therefore differentiate two discretisations that differ by O(tol), each of
# them exactly. The tolerances below are set tight enough that the gap sits far
# under them, and one test measures the gap itself rather than bounding it.
#
# Where the adjoint is checked at rounding level, on one shared step sequence,
# is cppDE's dev/cxx/test_reverse_*.cpp.

skip_on_cran()

.rev_dir <- function() {
  d <- file.path(tempdir(), "dmod_rev")
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  d
}

# Tight, so the discretisation gap between the two modes stays negligible.
.rev_opt <- list(atol = 1e-11, rtol = 1e-11)

.rev_reactions <- function()
  addReaction(addReaction(eqnlist(), "A", "B", "k1*A", "conversion"),
              "B", "", "k2*B", "decay")

# One compiled chain for the whole file: three compilations per test would
# dominate its runtime and prove nothing extra.
.rev_fx <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    d <- .rev_dir()
    owd <- setwd(d); on.exit(setwd(owd))
    re <- .rev_reactions()
    m <- odemodel(re, modelname = "rv_ode", deriv = TRUE, derivMode = c("forward", "reverse"),
                  outdir = d, compile = TRUE)
    x <- Xs(m, optionsOde = .rev_opt, optionsSens = .rev_opt)
    g <- Y(c(obsA = "s*A", obsB = "s*B"), re, compile = TRUE,
           modelname = "rv_obs", outdir = d)
    e <- Y(c(obsA = "sd_rel*obsA + sd_abs", obsB = "sd_rel*obsB + sd_abs"), g,
           states = c("obsA", "obsB"), parameters = c("sd_rel", "sd_abs"),
           compile = TRUE, modelname = "rv_err", outdir = d)
    tr <- c(A = "exp(logA)", B = "0", k1 = "exp(logk1)", k2 = "exp(logk2)",
            s = "exp(logs)")
    p <- P(tr, condition = "C1", compile = TRUE, modelname = "rv_p", outdir = d)
    pe <- P(c(tr, sd_rel = "exp(logsdrel)", sd_abs = "exp(logsdabs)"),
            condition = "C1", compile = TRUE, modelname = "rv_pe", outdir = d)
    cache <<- list(dir = d, m = m, x = x, g = g, e = e, p = p, pe = pe,
                   times = seq(0, 8, length.out = 41),
                   pars = c(logA = log(2), logk1 = log(0.6),
                            logk2 = log(0.3), logs = log(1.5)))
    cache
  }
})

# Data on the prediction's own grid, so nothing has to be interpolated onto it.
.rev_data <- function(fx, prd, pars, conds = "C1", sigma = 0.1, seed = 4L) {
  pred <- prd(fx$times, pars)
  set.seed(seed)
  idx <- c(4L, 9L, 16L, 26L, 36L)
  d <- do.call(rbind, lapply(conds, function(cn)
    do.call(rbind, lapply(c("obsA", "obsB"), function(nm) {
      out <- data.frame(name = nm, time = pred[[cn]][idx, "time"],
                        value = pred[[cn]][idx, nm] *
                                exp(rnorm(length(idx), 0, 0.05)),
                        condition = cn, stringsAsFactors = FALSE)
      if (!is.null(sigma)) out$sigma <- sigma
      out
    }))))
  as.datalist(d)
}

# The two gradients, aligned and compared on the forward one's scale.
expect_modes_agree <- function(obj, pars, tolerance = 1e-6, info = NULL) {
  f <- obj(pars, deriv = TRUE)
  r <- obj(pars, deriv = TRUE, sweep = "reverse")
  expect_equal(r$value, f$value, tolerance = 1e-8, info = info)
  expect_equal(unname(r$gradient[names(f$gradient)]), unname(f$gradient),
               tolerance = tolerance, info = info)
  invisible(list(forward = f, reverse = r))
}


test_that("the solver alone answers what the forward sensitivities answer", {
  fx <- .rev_fx()
  inner <- c(A = 2, B = 0, k1 = 0.6, k2 = 0.3)
  pred <- fx$x(fx$times, inner)[[1]]

  set.seed(1)
  w <- matrix(rnorm(nrow(pred) * 2), nrow(pred), 2,
              dimnames = list(NULL, c("A", "B")))
  ref <- apply(attr(pred, "deriv") * as.vector(w), 3, sum)

  vjp <- attr(attr(fx$x, "mappings")[[1]], "vjpfn")
  expect_false(is.null(vjp))
  got <- vjp(fx$times, inner, NULL, w)

  expect_equal(unname(got[names(ref)]), unname(ref), tolerance = 1e-7)
})

test_that("normL2 carries the whole chain backwards", {
  fx  <- .rev_fx()
  prd <- fx$g * fx$x * fx$p
  obj <- normL2(.rev_data(fx, prd, fx$pars), prd)

  both <- expect_modes_agree(obj, fx$pars)
  # The reverse objective is gradient-only by construction: there is no J to
  # contract, which is exactly what the quasi-Newton arm wants.
  expect_null(both$reverse$hessian)
  expect_false(is.null(both$forward$hessian))
})

test_that("the reverse mode reaches every condition and every branch", {
  fx <- .rev_fx()
  d  <- fx$dir
  owd <- setwd(d); on.exit(setwd(owd))
  conds <- c("C1", "C2")
  tr <- c(A = "exp(logA)", B = "0", k1 = "exp(logk1)", k2 = "exp(logk2)",
          s = "exp(logs)")
  p2 <- Reduce("+", lapply(conds, function(cn)
    P(repar(paste0("logk1 ~ logk1 + dk_", cn), tr), condition = cn,
      compile = TRUE, modelname = paste0("rv_p2_", cn), outdir = d)))

  pars <- c(fx$pars, dk_C1 = 0.1, dk_C2 = -0.15)
  prd  <- fx$g * fx$x * p2
  obj  <- normL2(.rev_data(fx, prd, pars, conds, seed = 7L), prd)

  both <- expect_modes_agree(obj, pars)
  # A condition-specific parameter must come back with a contribution from its
  # own branch and from nowhere else, which a shared-parameter-only check would
  # not notice.
  expect_true(all(abs(both$reverse$gradient[c("dk_C1", "dk_C2")]) > 1e-6))

  # Two conditions is where the batched backward solve takes over from the
  # loop, so this is where it has to answer the same thing.
  withr::with_options(list(dMod.batch.check = TRUE), {
    r2 <- obj(pars, deriv = TRUE, sweep = "reverse")
    expect_equal(r2$gradient, both$reverse$gradient)
  })
})

test_that("an estimated error model seeds the prediction a second time", {
  fx <- .rev_fx()
  prd <- fx$g * fx$x * fx$pe
  pars <- c(fx$pars, logsdrel = log(0.08), logsdabs = log(0.02))
  # sigma unknown, so the error model has to supply it and carries theta itself.
  obj <- normL2(.rev_data(fx, prd, pars, sigma = NULL), prd, fx$e)

  both <- expect_modes_agree(obj, pars)
  expect_true(all(abs(both$reverse$gradient[c("logsdrel", "logsdabs")]) > 1e-6))
})

test_that("the gap to the forward mode is the discretisation, not the adjoint", {
  fx <- .rev_fx()
  prd <- fx$g * fx$x * fx$p
  data <- .rev_data(fx, prd, fx$pars)

  rel <- vapply(10^-c(4, 10), function(tt) {
    o  <- list(atol = tt, rtol = tt)
    xx <- Xs(fx$m, optionsOde = o, optionsSens = o)
    oo <- normL2(data, fx$g * xx * fx$p)
    a  <- oo(fx$pars, deriv = TRUE)
    b  <- oo(fx$pars, deriv = TRUE, sweep = "reverse")
    max(abs(a$gradient - b$gradient[names(a$gradient)])) / max(abs(a$gradient))
  }, numeric(1))

  # Six decades of tolerance buy most of six decades of agreement. A channel
  # the reverse pass had missed would leave a floor here instead.
  expect_lt(rel[2], rel[1] * 1e-3)
})

test_that("an event with an estimated dose goes backwards too", {
  fx <- .rev_fx()
  d  <- fx$dir
  owd <- setwd(d); on.exit(setwd(owd))
  ev <- eventlist(var = "A", time = "t_dose", value = "d_amt", method = "add")
  m <- odemodel(.rev_reactions(), events = ev, modelname = "rv_ev",
                deriv = TRUE, derivMode = c("forward", "reverse"), outdir = d, compile = TRUE)
  xv <- Xs(m, optionsOde = .rev_opt, optionsSens = .rev_opt)
  pv <- P(c(A = "exp(logA)", B = "0", k1 = "exp(logk1)", k2 = "exp(logk2)",
            s = "exp(logs)", t_dose = "3", d_amt = "exp(logdose)"),
          condition = "C1", compile = TRUE, modelname = "rv_pev", outdir = d)

  pars <- c(fx$pars, logdose = log(0.8))
  prd  <- fx$g * xv * pv
  obj  <- normL2(.rev_data(fx, prd, pars), prd)

  both <- expect_modes_agree(obj, pars, tolerance = 1e-5)
  # The jump itself has to carry a derivative, or logdose comes back at zero.
  expect_gt(abs(both$reverse$gradient[["logdose"]]), 1e-6)
})

test_that("a summed objective passes the direction on", {
  fx <- .rev_fx()
  prd <- fx$g * fx$x * fx$p
  obj <- normL2(.rev_data(fx, prd, fx$pars), prd) +
         constraintL2(fx$pars * 0, sigma = 4)
  # The constraint has no reverse path of its own and keeps the forward one;
  # the sum is still the same number, and still carries no Hessian.
  both <- expect_modes_agree(obj, fx$pars)
  expect_null(both$reverse$hessian)
})

test_that("a model without a reverse object says so", {
  fx <- .rev_fx()
  d  <- fx$dir
  owd <- setwd(d); on.exit(setwd(owd))
  m <- odemodel(.rev_reactions(), modelname = "rv_noRev", deriv = TRUE,
                outdir = d, compile = TRUE)
  xf <- Xs(m)
  prd <- fx$g * xf * fx$p
  obj <- normL2(.rev_data(fx, fx$g * fx$x * fx$p, fx$pars), prd)
  expect_error(obj(fx$pars, deriv = TRUE, sweep = "reverse"),
               "no reverse object")
  expect_error(odemodel(.rev_reactions(), modelname = "rv_deSolve",
                        derivMode = c("forward", "reverse"), backend = "deSolve", outdir = d),
               "cppDE")
})


test_that("a steady-state transformation goes backwards too", {
  fx <- .rev_fx()
  d  <- fx$dir
  owd <- setwd(d); on.exit(setwd(owd))

  # A* = k_in / k_out feeding the decay chain's initial A, so the gradient has
  # to pass through the nested steady state to reach logkin and logkout.
  # The reverse path solves the nested steady state twice, once for the value
  # and once for the Jacobian, and each lands within roottol of the fixed point.
  # That gap is the sub-solve's own and has nothing to do with the adjoint, so
  # it is tightened out of the way rather than tolerated.
  pq <- Pequil(c(A = "k_in - k_out * A"), parameters = c("k_in", "k_out"),
               modelname = "rv_equil", compile = TRUE, attach.input = TRUE,
               deriv = TRUE, outdir = d, verbose = FALSE,
               controlsODE = list(abstol = 1e-12, reltol = 1e-12,
                                  roottol = 1e-12))
  pl <- P(c(k_in = "exp(logkin)", k_out = "exp(logkout)", B = "0",
            k1 = "exp(logk1)", k2 = "exp(logk2)", s = "exp(logs)"),
          condition = "C1", compile = TRUE, modelname = "rv_pq", outdir = d)

  prd  <- fx$g * fx$x * pq * pl
  pars <- c(logkin = log(1.5), logkout = log(0.75),
            logk1 = log(0.6), logk2 = log(0.3), logs = log(1.5))
  obj  <- normL2(.rev_data(fx, prd, pars, seed = 11L), prd)

  both <- expect_modes_agree(obj, pars, tolerance = 1e-5)
  # Both steady-state parameters have to arrive; a Jacobian contracted on the
  # wrong side would silently zero one of them.
  expect_true(all(abs(both$reverse$gradient[c("logkin", "logkout")]) > 1e-6))
})


test_that("censored rows go backwards on every BLOQ treatment", {
  fx  <- .rev_fx()
  prd <- fx$g * fx$x * fx$p
  d   <- as.data.frame(.rev_data(fx, prd, fx$pars))
  # A limit above the smaller observable's tail, so some rows really are
  # censored and the kernel's BLOQ branch is the one under test.
  d$lloq <- stats::quantile(d$value, 0.35)
  expect_gt(sum(d$value <= d$lloq), 1L)
  dl <- as.datalist(d)

  for (mode in c("M3", "M4NM", "M4BEAL", "M1")) {
    obj <- normL2(dl, prd, opt.BLOQ = mode)
    expect_modes_agree(obj, fx$pars, info = mode)
  }
})

test_that("a fixed sigma seeds the prediction and nothing else", {
  fx <- .rev_fx()
  prd <- fx$g * fx$x * fx$pe
  pars <- c(fx$pars, logsdrel = log(0.08), logsdabs = log(0.02))
  # sigma given in the data, so the error model is there but carries no
  # derivative for these rows: its cotangent has to be dropped rather than
  # multiplied by a zero that is never formed.
  obj <- normL2(.rev_data(fx, prd, pars, sigma = 0.1), prd, fx$e)
  both <- expect_modes_agree(obj, pars)
  expect_equal(unname(both$reverse$gradient[c("logsdrel", "logsdabs")]),
               c(0, 0), tolerance = 1e-10)
})
test_that("the Sundials backend goes backwards too", {
  skip_if_not(isTRUE(cppDE:::cvodeConfig$available),
              "CVODE backend not available")
  fx <- .rev_fx()
  d  <- .rev_dir()
  owd <- setwd(d); on.exit(setwd(owd))

  # CVODES adjoint sensitivity analysis under the same chain the native reverse
  # mode uses. It is a third discretisation: the adjoint is solved as its own
  # ODE over checkpointed forward states rather than by replaying the steps, so
  # this is a cross-check by foreign mathematics and not a repeat.
  m <- odemodel(.rev_reactions(), modelname = "rv_sun", deriv = TRUE,
                backend = "Sundials", derivMode = c("forward", "reverse"),
                outdir = d, compile = TRUE)
  expect_false(is.null(m$reversed))

  x   <- Xs(m, optionsOde = .rev_opt, optionsSens = .rev_opt)
  prd <- fx$g * x * fx$p
  obj <- normL2(.rev_data(fx, fx$g * fx$x * fx$p, fx$pars), prd)

  fo <- obj(fx$pars, deriv = TRUE, hessian = FALSE)
  rv <- obj(fx$pars, deriv = TRUE, sweep = "reverse")

  expect_null(rv$hessian)
  expect_equal(rv$value, fo$value, tolerance = 1e-8)
  expect_equal(unname(rv$gradient[names(fo$gradient)]), unname(fo$gradient),
               tolerance = 1e-5)
})

test_that("the Sundials reverse object refuses events", {
  skip_if_not(isTRUE(cppDE:::cvodeConfig$available),
              "CVODE backend not available")
  d <- .rev_dir()
  owd <- setwd(d); on.exit(setwd(owd))
  ev <- data.frame(var = "A", time = 1, value = 0.2, method = "add",
                   stringsAsFactors = FALSE)
  # The native backend replays the jump; CVODES integrates the adjoint over
  # checkpointed states and has no way to be told about one.
  expect_error(odemodel(.rev_reactions(), modelname = "rv_sun_ev",
                        backend = "Sundials", events = ev,
                        derivMode = c("forward", "reverse"), outdir = d),
               "does not support events")
})

test_that("odemodel builds the forward-reverse object and names it", {
  d <- .rev_dir()
  owd <- setwd(d); on.exit(setwd(owd))
  m <- odemodel(.rev_reactions(), modelname = "rv_fr", outdir = d,
                derivMode = c("forward", "forward-reverse"), compile = TRUE,
                nStack = 4)

  expect_null(m$extended2)
  expect_null(m$reversed)
  expect_false(is.null(m$reversed2))
  expect_identical(attr(m$reversed2, "derivMode"), "forward-reverse")

  # forward-forward is the older deriv2 = TRUE under its own name.
  m2 <- odemodel(.rev_reactions(), modelname = "rv_ffm", outdir = d,
                 derivMode = c("forward", "forward-forward"), compile = FALSE)
  expect_false(is.null(m2$extended2))
  expect_identical(attr(m2$extended2, "derivMode"), "forward-forward")

  expect_error(odemodel(.rev_reactions(), modelname = "rv_fr_bad", outdir = d,
                        backend = "Sundials", derivMode = "forward-reverse",
                        compile = FALSE),
               "forward-reverse")
})
