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
    m <- odemodel(re, modelname = "rv_ode", deriv = TRUE, reverse = TRUE,
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
                deriv = TRUE, reverse = TRUE, outdir = d, compile = TRUE)
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
               "reverse = TRUE")
  expect_error(odemodel(.rev_reactions(), modelname = "rv_deSolve",
                        reverse = TRUE, backend = "deSolve", outdir = d),
               "cppDE")
})
