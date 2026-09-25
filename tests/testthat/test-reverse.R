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

# Every model the file evaluates, generated first and linked into one shared
# object: a build per model and derivative direction would dominate the runtime
# and prove nothing extra.
.rev_models <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    d <- .rev_dir()
    owd <- setwd(d); on.exit(setwd(owd))
    re <- .rev_reactions()
    fr <- c("forward", "reverse")
    tr <- c(A = "exp(logA)", B = "0", k1 = "exp(logk1)", k2 = "exp(logk2)",
            s = "exp(logs)")

    m <- odemodel(re, modelname = "rv_ode", deriv = TRUE, derivMode = fr,
                  outdir = d, compile = FALSE)
    g <- Y(c(obsA = "s*A", obsB = "s*B"), re, derivMode = fr, modelname = "rv_obs",
           outdir = d)
    e <- Y(c(obsA = "sd_rel*obsA + sd_abs", obsB = "sd_rel*obsB + sd_abs"), g,
           states = c("obsA", "obsB"), parameters = c("sd_rel", "sd_abs"),
           derivMode = fr, modelname = "rv_err", outdir = d)
    p <- P(tr, condition = "C1", derivMode = fr, modelname = "rv_p", outdir = d)
    pe <- P(c(tr, sd_rel = "exp(logsdrel)", sd_abs = "exp(logsdabs)"),
            condition = "C1", derivMode = fr, modelname = "rv_pe", outdir = d)
    # One trafo per condition, each with a parameter of its own.
    pc <- Reduce("+", lapply(c("C1", "C2"), function(cn)
      P(repar(paste0("logk1 ~ logk1 + dk_", cn), tr), condition = cn, derivMode = fr,
        modelname = paste0("rv_p2_", cn), outdir = d)))

    ev <- eventlist(var = "A", time = "t_dose", value = "d_amt", method = "add")
    mev <- odemodel(re, events = ev, modelname = "rv_ev", deriv = TRUE,
                    derivMode = fr, outdir = d, compile = FALSE)
    pev <- P(c(tr, t_dose = "3", d_amt = "exp(logdose)"), condition = "C1",
             derivMode = fr, modelname = "rv_pev", outdir = d)
    mnr <- odemodel(re, modelname = "rv_noRev", deriv = TRUE, outdir = d,
                    compile = FALSE)

    # The steady state is solved twice on the reverse path, once for the value
    # and once for the Jacobian; roottol keeps that sub-solve gap out of the way.
    pq <- Pequil(c(A = "k_in - k_out * A"), parameters = c("k_in", "k_out"),
                 modelname = "rv_equil", attach.input = TRUE, deriv = TRUE,
                 outdir = d, verbose = FALSE,
                 controlsODE = list(abstol = 1e-12, reltol = 1e-12,
                                    roottol = 1e-12))
    pl <- P(c(k_in = "exp(logkin)", k_out = "exp(logkout)", B = "0",
              k1 = "exp(logk1)", k2 = "exp(logk2)", s = "exp(logs)"),
            condition = "C1", derivMode = fr, modelname = "rv_pq", outdir = d)

    sun <- if (isTRUE(cppDE:::cvodeConfig$available))
      odemodel(re, modelname = "rv_sun", deriv = TRUE, backend = "Sundials",
               derivMode = fr, outdir = d, compile = FALSE)

    m2 <- odemodel(re, modelname = "rv2_ode", deriv = TRUE, deriv2 = TRUE,
                   outdir = d, compile = FALSE,
                   derivMode = c(fr, "forward-forward", "forward-reverse"))
    g2 <- Y(c(obsA = "s*A", obsB = "s*B"), re, deriv2 = TRUE, derivMode = c(fr, "forward-reverse"),
            modelname = "rv2_obs", outdir = d)
    q <- P(tr, condition = "C1", deriv2 = TRUE, derivMode = c(fr, "forward-reverse"),
           modelname = "rv2_p", outdir = d)
    # A second condition on the same ODE. The batched backward path only
    # engages with more than one live condition.
    q2 <- q + P(tr, condition = "C2", deriv2 = TRUE, derivMode = c(fr, "forward-reverse"),
                modelname = "rv2_p2", outdir = d)

    compile(m, g, e, p, pe, pc, mev, pev, mnr, pq, pl, sun, m2, g2, q, q2,
            output = "rv_all", cores = 4L)

    pars <- c(logA = log(2), logk1 = log(0.6), logk2 = log(0.3), logs = log(1.5))
    cache <<- list(
      first = list(dir = d, m = m,
                   x = Xs(m, optionsOde = .rev_opt, optionsSens = .rev_opt),
                   g = g, e = e, p = p, pe = pe, pc = pc, mev = mev, pev = pev,
                   mnr = mnr, pq = pq, pl = pl, sun = sun,
                   times = seq(0, 8, length.out = 41), pars = pars),
      second = list(x = Xs(m2, optionsOde = .rev_opt, optionsSens = .rev_opt),
                    g = g2, p = q, p2 = q2,
                    times = seq(0, 8, length.out = 21), pars = pars))
    cache
  }
})

.rev_fx  <- function() .rev_models()$first
.rev2_fx <- function() .rev_models()$second

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

  expect_equal(unname(got[names(ref), 1L]), unname(ref), tolerance = 1e-7)
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
  conds <- c("C1", "C2")
  pars <- c(fx$pars, dk_C1 = 0.1, dk_C2 = -0.15)
  prd  <- fx$g * fx$x * fx$pc
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
  xv <- Xs(fx$mev, optionsOde = .rev_opt, optionsSens = .rev_opt)

  pars <- c(fx$pars, logdose = log(0.8))
  prd  <- fx$g * xv * fx$pev
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
  xf <- Xs(fx$mnr)
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

  # A* = k_in / k_out feeding the decay chain's initial A, so the gradient has
  # to pass through the nested steady state to reach logkin and logkout.
  prd  <- fx$g * fx$x * fx$pq * fx$pl
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

  # CVODES adjoint sensitivity analysis under the same chain the native reverse
  # mode uses. It is a third discretisation: the adjoint is solved as its own
  # ODE over checkpointed forward states rather than by replaying the steps, so
  # this is a cross-check by foreign mathematics and not a repeat.
  m <- fx$sun
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
                        derivMode = c("forward", "reverse"), outdir = d,
                        compile = FALSE),
               "does not support events")
})

test_that("odemodel builds the forward-reverse object and names it", {
  d <- .rev_dir()
  owd <- setwd(d); on.exit(setwd(owd))
  # Code generation only: the second-order models compile the same
  # forward-reverse source.
  m <- odemodel(.rev_reactions(), modelname = "rv_fr", outdir = d,
                derivMode = c("forward", "forward-reverse"), compile = FALSE)

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

# ---------------------------------------------------------------------------
#  Second order: the exact Hessian through the whole chain.
#
#  Oracle is deriv2 = TRUE, forward over forward, on the same objective. The two
#  differentiate two discretisations, the forward one under sensitivity error
#  control and the backward one on the value run's grid, so the gap is O(tol)
#  and the same one the first-order tests measure.
# ---------------------------------------------------------------------------

# Data on whatever the chain's own columns are called, on its own grid.
.rev2_data <- function(fx, chain, nms, seed = 4L, conditions = "C1") {
  pred <- chain(fx$times, fx$pars)
  set.seed(seed)
  idx <- c(4L, 8L, 12L, 16L, 20L)
  d <- do.call(rbind, lapply(conditions, function(cc)
    do.call(rbind, lapply(nms, function(nm)
      data.frame(name = nm, time = pred[[cc]][idx, "time"],
                 value = pred[[cc]][idx, nm] * exp(rnorm(length(idx), 0, 0.05)),
                 sigma = 0.1, condition = cc, stringsAsFactors = FALSE)))))
  as.datalist(d)
}

test_that("the chain answers the Hessian forward over forward answers", {
  fx <- .rev2_fx()
  for (nm in c("x * p", "g * x * p")) {
    chain <- if (nm == "x * p") fx$x * fx$p else fx$g * fx$x * fx$p
    cols  <- if (nm == "x * p") c("A", "B") else c("obsA", "obsB")
    obj   <- normL2(.rev2_data(fx, chain, cols), chain)

    fwd <- obj(fx$pars, deriv2 = TRUE)
    rev <- obj(fx$pars, sweep = "reverse", deriv2 = TRUE)

    expect_identical(attr(rev, "sweep"), "forward-reverse", info = nm)
    expect_equal(rev$gradient, fwd$gradient, tolerance = 1e-3, info = nm)
    expect_equal(rev$hessian, fwd$hessian, tolerance = 1e-4, info = nm)
    # Symmetric on one grid, and on the objective's own parameter order.
    expect_equal(rev$hessian, t(rev$hessian), tolerance = 1e-8, info = nm)
    expect_identical(dimnames(rev$hessian),
                     list(names(fx$pars), names(fx$pars)), info = nm)
  }
})

test_that("the batched backward path carries the directions too", {
  # The batched leaf sizes its pass-through half from the cotangent it was
  # handed, not from its own. With one condition the batch entry never engages,
  # so this is the first place a K-column answer meets a one-column neighbour.
  fx    <- .rev2_fx()
  chain <- fx$x * fx$p2
  obj   <- normL2(.rev2_data(fx, chain, c("A", "B"), conditions = c("C1", "C2")),
                  chain)

  # Guard against a vacuous pass: the batched route needs a batch vjp on the
  # prediction leaf and more than one live condition. Without both, the loop
  # this test is about never runs.
  expect_false(is.null(attr(attr(fx$x, "mappings")[[1L]], "vjpbatchfn")))
  expect_length(chain(fx$times, fx$pars), 2L)

  fwd <- obj(fx$pars, deriv2 = TRUE)
  rev <- obj(fx$pars, sweep = "reverse", deriv2 = TRUE)

  expect_identical(attr(rev, "sweep"), "forward-reverse")
  expect_equal(rev$gradient, fwd$gradient, tolerance = 1e-3)
  expect_equal(rev$hessian, fwd$hessian, tolerance = 1e-4)
  expect_equal(rev$hessian, t(rev$hessian), tolerance = 1e-8)
  # The batched and the scalar route are the same arithmetic; dMod.batch.check
  # re-runs every batched leaf through the scalar kernel and compares.
  withr::local_options(dMod.batch.check = TRUE)
  expect_equal(obj(fx$pars, sweep = "reverse", deriv2 = TRUE)$hessian,
               rev$hessian, tolerance = 0)
})

test_that("a reverse evaluation says which direction answered it", {
  fx <- .rev2_fx()
  chain <- fx$x * fx$p
  obj <- normL2(.rev2_data(fx, chain, c("A", "B")), chain)

  first <- obj(fx$pars, sweep = "reverse")
  expect_identical(attr(first, "sweep"), "reverse")
  expect_null(first$hessian)
})

test_that("a summed objective keeps every term's curvature", {
  # A constraint has no reverse path of its own, so it runs forward. Under an
  # exact request it must still hand back its Hessian: .sumobjlist adds an
  # absent one as zero, so a dropped term would leave the total short of that
  # term's curvature and nothing would fail.
  fx    <- .rev2_fx()
  chain <- fx$x * fx$p
  dat   <- .rev2_data(fx, chain, c("A", "B"))
  mu    <- fx$pars * 0.9
  obj   <- normL2(dat, chain) + constraintL2(mu = mu, sigma = 0.3)

  fwd <- obj(fx$pars, deriv2 = TRUE)
  rev <- obj(fx$pars, sweep = "reverse", deriv2 = TRUE)

  expect_equal(rev$value, fwd$value, tolerance = 1e-8)
  expect_equal(rev$gradient, fwd$gradient, tolerance = 1e-3)
  expect_equal(rev$hessian, fwd$hessian, tolerance = 1e-4)

  # The constraint's own curvature is 2/sigma^2 on the diagonal and does not
  # vanish, so a total that dropped it would differ by exactly that much.
  bare <- normL2(dat, chain)(fx$pars, sweep = "reverse", deriv2 = TRUE)
  expect_gt(max(abs(rev$hessian - bare$hessian)), 1)
})

test_that("trust drives the exact Hessian, forwards and backwards", {
  fx    <- .rev2_fx()
  chain <- fx$x * fx$p
  obj   <- normL2(.rev2_data(fx, chain, c("A", "B")), chain)
  start <- fx$pars + c(0.4, -0.35, 0.3, 0.2)

  gn <- trust(obj, start, rinit = 0.1, rmax = 10, iterlim = 100L)
  expect_true(gn$converged)

  # A Newton run: the objective's own Hessian at every iterate. The subproblem
  # solver takes an indefinite matrix natively, so this needed no new algebra.
  nw <- trust(obj, start, rinit = 0.1, rmax = 10, iterlim = 100L,
              hessianMethod = "exact")
  expect_true(nw$converged)
  expect_equal(nw$value, gn$value, tolerance = 1e-6)

  # The shape this whole mode exists for: one exact Hessian at the start, then
  # a descent on reverse gradients alone.
  rv <- trust(obj, start, rinit = 0.1, rmax = 10, iterlim = 100L,
              hessianMethod = "sr1", sweep = "reverse",
              qnControl = list(hessianInit = "exact"))
  expect_true(rv$converged)
  expect_equal(rv$value, gn$value, tolerance = 1e-6)

  # A stalled quasi-Newton phase fetching a fresh curvature rather than
  # stopping. Whether it stalls at all is decided at the solver's noise floor
  # and therefore by the platform, which is the point: "stall" has to reach the
  # same place either way and may never cost more than "never" does.
  rs <- trust(obj, start, rinit = 0.1, rmax = 10, iterlim = 100L,
              hessianMethod = "sr1", sweep = "reverse",
              qnControl = list(hessianInit = "exact", hessianReseed = "stall"))
  expect_true(rs$converged)
  expect_true(rs$stopReason %in% c("gradient", "stagnation"))
  expect_equal(rs$value, gn$value, tolerance = 1e-6)
  # A reseed at a standing iterate would refetch the same matrix, so one is
  # taken only after the iterate moves and the run cannot loop on it.
  expect_lt(rs$iterations, 100L)

  # Reverse and forward Newton reach the same place from the same start.
  nr <- trust(obj, start, rinit = 0.1, rmax = 10, iterlim = 100L,
              hessianMethod = "exact", sweep = "reverse")
  expect_true(nr$converged)
  expect_equal(nr$value, nw$value, tolerance = 1e-6)
})

test_that("a reseed needs a moved iterate and cannot outlast one", {
  # Whether the run above stalls is decided at the solver's noise floor, so it
  # does on some platforms and not on others. gtol shut off forces the case.
  fx    <- .rev2_fx()
  chain <- fx$x * fx$p
  obj   <- normL2(.rev2_data(fx, chain, c("A", "B")), chain)
  start <- fx$pars + c(0.4, -0.35, 0.3, 0.2)
  qn    <- list(hessianInit = "exact", hessianReseed = "stall")
  args  <- list(obj, start, rinit = 0.1, rmax = 10, iterlim = 100L,
                hessianMethod = "sr1", sweep = "reverse",
                tolControl = list(gtol = 1e-14))

  rs <- do.call(trust, c(args, list(qnControl = qn)))
  qn$hessianReseed <- "never"
  nv <- do.call(trust, c(args, list(qnControl = qn)))

  # A stall is made of rejected steps, so refetching at the same iterate would
  # return the same matrix. Reseeds are taken, and then the run stops the way
  # "never" does instead of spending its budget on them.
  expect_gt(rs$nReseed, 0L)
  expect_true(rs$converged)
  expect_identical(rs$stopReason, "stagnation")
  expect_lt(rs$iterations, 100L)
  expect_identical(nv$nReseed, 0L)
  expect_equal(rs$value, nv$value, tolerance = 1e-6)
})

test_that("an exact Hessian asked of an objective that cannot give one says so", {
  # A configuration error, caught before the run. Inside the closure it would be
  # swallowed by the kernel's evaluation handler and come back as "parinit not
  # feasible", which names the wrong thing.
  plain <- function(pars, deriv = TRUE, hessian = NULL) {
    objlist(value = sum(pars^2), gradient = 2 * pars,
            hessian = if (isTRUE(hessian))
              diag(2, length(pars), length(pars)) else NULL)
  }
  st <- c(a = 1, b = -1)
  expect_error(trust(plain, st, hessianMethod = "exact"), "deriv2")
  expect_error(trust(plain, st, hessianMethod = "sr1",
                     qnControl = list(hessianInit = "exact")), "deriv2")
  # Without an exact request it runs as it always did.
  expect_true(trust(plain, st, iterlim = 50L)$converged)
})

test_that("a contradictory request is overridden and says so", {
  fx    <- .rev2_fx()
  chain <- fx$x * fx$p
  obj   <- normL2(.rev2_data(fx, chain, c("A", "B")), chain)

  # A Gauss-Newton Hessian needs J, which the reverse mode does not build. The
  # request is dropped to the cheaper answer and warns rather than going quiet.
  expect_warning(gn <- obj(fx$pars, sweep = "reverse", hessian = TRUE),
                 "Gauss-Newton")
  expect_null(gn$hessian)

  # An exact Hessian and no Hessian cannot both hold; hessian = FALSE wins,
  # again with a warning.
  expect_warning(no <- obj(fx$pars, deriv2 = TRUE, hessian = FALSE),
                 "overrides")
  expect_null(no$hessian)

  # Taking the default is not a contradiction and stays silent, in either
  # direction.
  expect_silent(obj(fx$pars, sweep = "reverse"))
  expect_silent(obj(fx$pars, deriv2 = TRUE, sweep = "reverse"))
})
