# Behavioral tests for trust(): quadratic convergence, truth recovery, fixed
# parameters, bounds, Hessian sources and the trace.

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("cppDE")
  testthat::skip_on_cran()
}

fx_register(multicond = TRUE, extra = function(bench) list(
  trust_err = Y(c(y = "sqrt(s0^2 + srel^2 * y^2)"), f = bench$gfn, attachInput = FALSE,
                condition = "C1", modelname = "trust_err", compile = FALSE),
  trust_p_err = P(eqnvec(A = "exp(A_log)", k = "exp(k_log)", s0 = "exp(s0_log)", srel = "srel"),
                  condition = "C1", modelname = "trust_p_err", compile = FALSE)))


# Synthetic quadratic objective in arbitrary dimension. Minimum at `target`.
.quadratic_objfn <- function(target) {
  d <- length(target)
  function(p, ...) {
    list(
      value    = 0.5 * sum((p - target)^2),
      gradient = p - target,
      hessian  = diag(d))
  }
}


## ---- Quadratic: one-step convergence ----------------------------------

test_that("trust ends where a Newton step it rejects predicts less than mtol", {
  # the values rise once the iterate is near the minimum, so the Newton steps there are rejected
  hit <- FALSE
  obj <- function(pars, ...) {
    d <- pars[["x"]] - 1
    v <- d^2 + if (hit) 1e-5 else 0
    if (abs(d) < 1e-3) hit <<- TRUE
    objlist(v, c(x = 2 * d), matrix(2.2, 1, 1, dimnames = list("x", "x")))
  }
  fit <- trust(obj, c(x = 0), rinit = 1, rmax = 10, iterlim = 100)
  expect_true(fit$converged)
  expect_identical(fit$stopReason, "preddiff")
  expect_lt(fit$iterations, 10)
})

test_that("trust converges to the exact minimum of a quadratic in one outer step", {
  target <- c(a = 1.0, b = -0.5, c = 2.3)
  obj <- .quadratic_objfn(target)
  init <- c(a = 0, b = 0, c = 0)

  fit <- trust(obj, init, rinit = 5, rmax = 100, iterlim = 50,
               printIter = FALSE)
  expect_true(fit$converged)
  expect_equal(fit$argument, target, tolerance = 1e-8, ignore_attr = TRUE)
  # A pure quadratic with rinit large enough to admit the Newton step is
  # accepted on iteration 1; trust then runs one more iteration to verify
  # termination criteria.
  expect_lte(fit$iterations, 3L)
})

# A quadratic with minimum 0 at the origin and a flat second direction.
.flat_objfn <- function(p, ...) {
  h <- c(1, 1e-3)
  list(value = 0.5 * sum(h * p^2), gradient = h * p, hessian = diag(h))
}

test_that("a gradient stop needs the full Newton step to predict at most mtol", {
  fit <- trust(.flat_objfn, c(x = 0, y = 9), rinit = 1, rmax = 10,
               tolControl = list(ftol = 1e-4, mtol = 1e-4, gtol = 1e-2))
  expect_true(fit$converged)
  expect_identical(fit$stopReason, "gradient")
  expect_lte(fit$value, 1e-4)
})

test_that("a value stop needs the full Newton step to predict at most mtol", {
  fit <- trust(.flat_objfn, c(x = 0, y = 9), rinit = 0.05, rmax = 0.05, hessianMethod = "gn",
               tolControl = list(ftol = 1e-3, mtol = 1e-4, gtol = 1e-6))
  expect_true(fit$converged)
  expect_identical(fit$stopReason, "fvalue")
  expect_lte(fit$value, 1e-4)
})


test_that("a stop held back by the Newton decrease ends the run when the radius collapses", {
  # the gradient predicts a decrease along the flat direction that the noisy value does not follow
  h <- c(1, 1e-4)
  obj <- function(p, ...) {
    objlist(0.5 * sum(h * p^2) + 5e-5 * sin(1e6 * p[["y"]]), h * p + c(0, 0.01),
            diag(h) |> `dimnames<-`(list(names(p), names(p))))
  }
  fit <- trust(obj, c(x = 1, y = 1), rinit = 1, rmax = 10, iterlim = 200,
               tolControl = list(ftol = 1e-4, mtol = 1e-4, gtol = 1e-6, rmin = 1e-3))
  expect_true(fit$converged)
  expect_identical(fit$stopReason, "fvalue")
})


## ---- Simulated-truth recovery ----------------------------------------

test_that("trust(normL2) recovers true (A, k) from simulated decay data within noise budget", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  true_pars <- c(A = 1.0, k = 0.5)
  sigma_sim <- 0.02  # very low noise so fit is well-determined
  data <- fx_decay_data(pars = true_pars, sigma = sigma_sim,
                        times = seq(0, 8, by = 0.5), seed = 17L)

  obj <- normL2(data, bench$prd_id)
  init <- c(A = 0.6, k = 0.9)  # away from truth
  fit <- trust(obj, init, rinit = 1, rmax = 10, iterlim = 100,
               printIter = FALSE)

  expect_true(fit$converged)
  # Tolerances match the practical envelope at this noise level.
  expect_lt(abs(fit$argument[["A"]] - true_pars[["A"]]) / true_pars[["A"]],
            0.10)
  expect_lt(abs(fit$argument[["k"]] - true_pars[["k"]]) / true_pars[["k"]],
            0.10)
})


## ---- fixed = ... holds parameters constant ---------------------------

test_that("trust honors fixed = ... (held parameters unchanged in argument)", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(sigma = 0.05)
  obj   <- normL2(data, bench$prd_id)

  init  <- c(A = 0.7, k = 0.6)
  fix   <- c(k = 0.5)   # hold k, vary only A
  free  <- init[setdiff(names(init), names(fix))]

  fit <- trust(obj, free, rinit = 0.5, rmax = 5, iterlim = 50,
               fixed = fix, printIter = FALSE)
  expect_true(fit$converged)
  # `argument` contains only the free parameters; k should be absent.
  expect_false("k" %in% names(fit$argument))
  expect_true("A" %in% names(fit$argument))
})


## ---- bounds clamp the optimum ---------------------------------------

test_that("trust honors parupper on one component while leaving the other free", {
  # An upper bound on one component clamps it; the other converges freely.
  target <- c(a = 5.0, b = 1.0)
  obj <- .quadratic_objfn(target)

  fit <- trust(obj, parinit = c(a = 0, b = 0), rinit = 1, rmax = 10,
               iterlim = 50,
               parupper = c(a = 2, b = Inf),
               printIter = FALSE)
  expect_equal(unname(fit$argument[["a"]]), 2.0, tolerance = 1e-5)
  expect_equal(unname(fit$argument[["b"]]), 1.0, tolerance = 1e-5)
})


## ---- blather returns the per-iter trace ------------------------------

## ---- Near-singular Hessian: stay within trust radius ------------------

test_that("trust step stays inside the trust radius when the Hessian is near-singular", {
  # Hessian eigenvalues span many orders of magnitude, so LAPACK returns the
  # smallest as tiny negative noise; the Moré-Sorensen easy branch must still
  # keep ||p|| ~ r.
  K <- 7
  vals <- c(4.4e-3, 6.0e-4, 1.9e-5, 9.6e-7, 8.0e-10, 8.2e-19, 9.4e-21)
  set.seed(42)
  Q <- qr.Q(qr(matrix(rnorm(K * K), K, K)))    # random orthogonal basis
  H <- Q %*% diag(vals) %*% t(Q)
  H <- 0.5 * (H + t(H))
  g <- as.numeric(Q %*% c(0.5, 0.6, 0.4, 0.3, 0.2, 0.0, 0.0))
  par_target <- setNames(rep(0, K), paste0("p", seq_len(K)))
  init       <- setNames(rep(0.1, K), paste0("p", seq_len(K)))

  # Constant near-singular Hessian: many minima are acceptable, only the step
  # length is checked.
  objfun <- function(p, ...) {
    dp <- p - init
    list(value    = sum(g * dp) + 0.5 * as.numeric(t(dp) %*% H %*% dp),
         gradient = as.numeric(g + H %*% dp),
         hessian  = H)
  }

  fit <- suppressWarnings(
    trust(objfun, init, rinit = 0.1, rmax = 10, iterlim = 50,
          blather = TRUE, printIter = FALSE))
  expect_true(all(is.finite(fit$argument)))
  # Every accepted or proposed step stays within rmax and the per-iteration radius.
  expect_true(all(fit$stepnorm <= fit$r + 1e-6))
  expect_true(all(fit$stepnorm <= 10 + 1e-6))
})


test_that("trust(blather = TRUE) returns all trace fields with finite numbers", {
  target <- c(a = 1.0, b = -0.5, c = 2.3)
  obj <- .quadratic_objfn(target)
  init <- c(a = 0, b = 0, c = 0)

  fit <- trust(obj, init, rinit = 0.3, rmax = 5, iterlim = 20,
               blather = TRUE, printIter = FALSE)
  expect_true(fit$converged)
  n <- fit$iterations
  expect_equal(nrow(fit$argpath), n)
  expect_equal(ncol(fit$argpath), 3L)
  expect_equal(nrow(fit$argtry),  n)
  expect_equal(length(fit$steptype), n)
  expect_equal(length(fit$accept),   n)
  expect_equal(length(fit$r),        n)
  expect_equal(length(fit$rho),      n)
  expect_equal(length(fit$valpath),  n)
  expect_equal(length(fit$valtry),   n)
  expect_equal(length(fit$preddiff), n)
  expect_equal(length(fit$stepnorm), n)
  expect_true(all(fit$steptype %in% c("Newton", "easy-easy", "hard-easy", "hard-hard")))
  expect_true(all(is.finite(fit$valpath)))
  expect_true(all(fit$r >= 0))
})


## ---- Coleman-Li boundary handling -------------------------------------
# The reflective scheme keeps iterates strictly inside the box, so the bound
# is approached but never touched; `atBound` reports the activity instead.

# General quadratic 0.5 * (p - target)' A (p - target), A positive definite.
.quadratic_objfn_A <- function(target, A) {
  function(p, ...) {
    d <- as.numeric(p - target)
    list(value    = as.numeric(0.5 * t(d) %*% A %*% d),
         gradient = as.numeric(A %*% d),
         hessian  = A)
  }
}

# The Coleman-Li model value at step `s`: the plain quadratic model plus the
# curvature term of the metric itself, 0.5 * shat' C shat with C = |g| * jv.
.coleman_li_model <- function(theta, g, H, s, lb, ub) {
  bnd  <- ifelse(g < 0, ub, lb)
  absv <- ifelse(is.finite(bnd), abs(theta - bnd), 1)
  jv   <- ifelse(is.finite(bnd), 1, 0)
  shat <- s / sqrt(absv)
  as.numeric(t(g) %*% s + 0.5 * t(s) %*% H %*% s + 0.5 * sum(abs(g) * jv * shat^2))
}


test_that("trust hits the analytic KKT point of a bound-active quadratic", {
  # A coupled quadratic with an active upper bound: the result meets the KKT point.
  A      <- matrix(c(2, 1, 1, 2), 2, 2)
  target <- c(x = 3, y = 3)
  # Tolerances pinned: the assertions below measure KKT accuracy, which must
  # not drift with the package defaults.
  fit <- trust(.quadratic_objfn_A(target, A), c(x = 0, y = 0),
               rinit = 1, rmax = 10, iterlim = 100,
               tolControl = list(gtol = 1e-10, ftol = 0, mtol = 0),
               parupper = c(x = 1, y = Inf))

  expect_true(fit$converged)
  expect_equal(fit$stopReason, "gradient")
  expect_equal(unname(fit$argument[["x"]]), 1, tolerance = 1e-6)
  expect_equal(unname(fit$argument[["y"]]), 4, tolerance = 1e-6)
  # Strictly interior: never exactly on the bound, always below it.
  expect_lt(fit$argument[["x"]], 1)
  expect_equal(unname(fit$atBound), c(TRUE, FALSE))
  expect_equal(unname(fit$gradient), c(-3, 0), tolerance = 1e-5)
  # The extra diagonal term is what buys fast convergence at a bound; without
  # it this crawls in geometrically over dozens of iterations.
  expect_lte(fit$iterations, 15L)

  # The shipped defaults still end on the gradient test, at the accuracy
  # gtol implies.
  deflt <- trust(.quadratic_objfn_A(target, A), c(x = 0, y = 0),
                 rinit = 1, rmax = 10, iterlim = 100,
                 parupper = c(x = 1, y = Inf))
  expect_equal(deflt$stopReason, "gradient")
  expect_equal(unname(deflt$argument[["x"]]), 1, tolerance = 1e-4)
  expect_true(deflt$atBound[["x"]])
})


test_that("preddiff describes the step actually taken, bound active or not", {
  A      <- matrix(c(2, 1, 1, 2), 2, 2)
  target <- c(x = 3, y = 3)
  obj    <- .quadratic_objfn_A(target, A)
  lb     <- c(-Inf, -Inf)
  ub     <- c(1, Inf)

  fit <- trust(obj, c(x = 0, y = 0), rinit = 0.3, rmax = 10, iterlim = 100,
               parupper = c(x = 1, y = Inf), blather = TRUE)

  for (k in seq_len(fit$iterations)) {
    theta <- fit$argpath[k, ]
    s     <- fit$argtry[k, ] - theta
    o     <- obj(theta)
    expect_equal(fit$preddiff[k],
                 .coleman_li_model(theta, o$gradient, o$hessian, s, lb, ub),
                 tolerance = 1e-9)
  }
  expect_true(all(fit$stepnorm <= fit$r + 1e-9))
})


test_that("truncated and reflected steps are taken when the box blocks", {
  # Only off-diagonal coupling pushes the scaled Newton step out of the box,
  # so this is a correlated problem with the optimum far outside an asymmetric box.
  A <- matrix(c(4.789636, 1.795686,  0.196845,
                1.795686, 1.377690, -0.720839,
                0.196845, -0.720839, 1.548688), 3, 3)
  tg   <- c(a = -10.187971, b =  5.856299, c = 5.418100)
  lb   <- c(a =  -0.895015, b = -0.548275, c = -1.473029)
  ub   <- c(a =   0.530507, b =  0.821269, c = 0.988922)
  init <- c(a =  -0.539109, b = -0.150667, c = -0.316588)

  fit <- trust(.quadratic_objfn_A(tg, A), init, rinit = 1.014421, rmax = 20,
               iterlim = 100, tolControl = list(gtol = 1e-12),
               parlower = lb, parupper = ub, blather = TRUE)

  expect_true(fit$converged)
  expect_true(all(c("truncated", "reflected") %in% fit$stepback))
  # Strictly interior at every iterate, not just at the end.
  for (nm in names(init)) {
    expect_true(all(fit$argpath[, nm] > lb[[nm]]))
    expect_true(all(fit$argpath[, nm] < ub[[nm]]))
  }
  expect_true(all(fit$argument > lb) && all(fit$argument < ub))
})


test_that("a collapsing trust radius reports failure, not convergence", {
  # Gradient points uphill, so no step is ever accepted and the radius decays.
  # A flat rejected step must not count as convergence.
  misleading <- function(p, ...) {
    list(value = sum(p^2), gradient = -2 * p, hessian = diag(length(p)))
  }
  expect_warning(
    fit <- trust(misleading, c(a = 1, b = 1), rinit = 1, rmax = 10,
                 iterlim = 200, tolControl = list(rmin = 0.05)),
    "rmin")
  expect_false(fit$converged)
  expect_equal(fit$stopReason, "radius")
})


test_that("a step that no longer moves the iterate ends the run as a collapsed radius", {
  # every evaluation after the first is higher, at the start as well, so no step is accepted and
  # the values never turn flat
  calls <- 0L
  obj <- function(p, ...) {
    calls <<- calls + 1L
    h <- diag(2, length(p))
    dimnames(h) <- list(names(p), names(p))
    objlist(sum(p^2) + if (calls > 1L) 3 else 0, 2 * p, h)
  }
  expect_warning(fit <- trust(obj, c(a = 1, b = 1), rinit = 1, rmax = 10, iterlim = 500),
                 "resolution")
  expect_false(fit$converged)
  expect_identical(fit$stopReason, "radius")
  expect_lt(fit$iterations, 100L)
})

test_that("a flat objective under repeated rejection stops as stagnation", {
  # Without an rmin floor the run ends once the objective stops moving. That is
  # "stagnation", never "gradient": the optimizer cannot tell a noise floor from
  # a bad model.
  misleading <- function(p, ...) {
    list(value = sum(p^2), gradient = -2 * p, hessian = diag(length(p)))
  }
  fit <- trust(misleading, c(a = 1, b = 1), rinit = 1, rmax = 10, iterlim = 500)
  expect_equal(fit$stopReason, "stagnation")
  expect_lt(fit$iterations, 500L)
})


test_that("stagnation is convergence only where the Newton decrease is at most mtol", {
  # a value that follows no step: the gradient predicts a decrease the run never finds
  flat <- function(p, ...) list(value = 1, gradient = c(1, 1), hessian = diag(2))
  fit <- trust(flat, c(a = 0, b = 0), iterlim = 500)
  expect_identical(fit$stopReason, "stagnation")
  expect_false(fit$converged)
  # noise above ftol around the minimum of a quadratic: no step is accepted, the model is at rest
  start <- c(a = 1e-4, b = 0)
  noisy <- function(p, ...) list(value = sum(p^2) + 5e-5 * any(p != start), gradient = 2 * p,
                                 hessian = diag(2, 2))
  fit <- trust(noisy, start, rinit = 1e-6, iterlim = 500,
               tolControl = list(ftol = 1e-4, mtol = 1e-4, gtol = 1e-6))
  expect_identical(fit$stopReason, "stagnation")
  expect_true(fit$converged)
})

test_that("a rank-deficient Gauss-Newton Hessian gives a minimum-norm step", {
  # Rank-deficient Gauss-Newton Hessian with the gradient orthogonal to the null
  # space: the step must not extend along the zero-curvature direction.
  obj <- function(p, ...) {
    res <- as.numeric(p[1] + p[2] - 2)
    J   <- matrix(c(1, 1), nrow = 1)
    list(value = 0.5 * res^2, gradient = as.numeric(t(J) * res), hessian = t(J) %*% J)
  }
  fit <- trust(obj, c(x = 0, y = 0), rinit = 1, rmax = 10, iterlim = 100)

  expect_true(all(is.finite(fit$argument)))
  expect_equal(sum(fit$argument), 2, tolerance = 1e-6)
  # The minimum-norm solution, without null-space drift.
  expect_equal(unname(fit$argument[["x"]] - fit$argument[["y"]]), 0, tolerance = 1e-6)
})


test_that("curvature at roundoff level does not steer the degenerate step", {
  # A shifted Hessian makes lam_min negative by construction rather than by LAPACK
  # rounding, so the drift shows on every build.
  obj <- function(d) function(p, ...) {
    res <- as.numeric(p[1] + p[2] - 2)
    J   <- matrix(c(1, 1), nrow = 1)
    list(value = 0.5 * res^2, gradient = as.numeric(t(J) * res),
         hessian = t(J) %*% J - diag(d, 2))
  }
  drift <- function(d) {
    fit <- trust(obj(d), c(x = 0, y = 0), rinit = 1, rmax = 10, iterlim = 100)
    expect_equal(sum(fit$argument), 2, tolerance = 1e-6)   # the fit itself is unaffected
    unname(fit$argument[["x"]] - fit$argument[["y"]])
  }
  expect_equal(drift(1e-15), 0, tolerance = 1e-6)   # roundoff-level: ignored
  expect_gt(abs(drift(1e-6)), 1)                    # real: extended to the boundary
})


test_that("a null space wider than one direction still gives the minimum-norm step", {
  # A rank-1 Hessian with a two-dimensional null space: roundoff splits the zero
  # eigenvalues, so no LAPACK build returns a basis g is exactly orthogonal to.
  J <- matrix(c(1, 2, 3), nrow = 1)
  obj <- function(p, ...) {
    res <- as.numeric(J %*% p - 2)
    list(value = 0.5 * res^2, gradient = as.numeric(t(J) * res), hessian = t(J) %*% J)
  }
  fit <- trust(obj, c(x = 0, y = 0, z = 0), rinit = 1, rmax = 10, iterlim = 100)

  expect_equal(as.numeric(J %*% fit$argument), 2, tolerance = 1e-6)
  # The minimum-norm solution.
  expect_equal(unname(fit$argument), as.numeric(t(J)) * 2 / sum(J^2),
               tolerance = 1e-6)
})


test_that("bounds compose with parinit on a bound and minimize = FALSE", {
  target <- c(x = 5, y = 5)

  # parinit sitting exactly on the bound must be nudged inside, not frozen.
  on_bound <- trust(.quadratic_objfn(target), c(x = 1, y = 0),
                    rinit = 1, rmax = 10, parupper = c(x = 1, y = Inf))
  expect_true(on_bound$converged)
  expect_equal(unname(on_bound$argument[["y"]]), 5, tolerance = 1e-6)
  expect_true(on_bound$atBound[["x"]])

  # Maximizing a concave objective under the same bound.
  negobj <- function(p, ...) {
    o <- .quadratic_objfn(target)(p)
    list(value = -o$value, gradient = -o$gradient, hessian = -o$hessian)
  }
  mx <- trust(negobj, c(x = 0, y = 0), rinit = 1, rmax = 10,
              parupper = c(x = 1, y = Inf), minimize = FALSE)
  expect_true(mx$converged)
  expect_equal(unname(mx$argument[["y"]]), 5, tolerance = 1e-6)
  expect_true(mx$atBound[["x"]])
})


test_that("without bounds the two boundary schemes agree", {
  # Without bounds |v| == 1 and C == 0, so the reflective scheme reduces to the plain one.
  target <- c(a = 1.0, b = -0.5, c = 2.3)
  init   <- c(a = 0, b = 0, c = 0)
  refl <- trust(.quadratic_objfn(target), init, rinit = 0.5, rmax = 10,
                boundary = "reflective")
  clip <- trust(.quadratic_objfn(target), init, rinit = 0.5, rmax = 10,
                boundary = "clip")
  expect_equal(unname(refl$argument), unname(clip$argument), tolerance = 1e-10)
})


test_that("boundary = 'clip' keeps landing exactly on the bound", {
  fit <- trust(.quadratic_objfn(c(a = 5, b = 1)), c(a = 0, b = 0),
               rinit = 1, rmax = 10, parupper = c(a = 2, b = Inf),
               boundary = "clip")
  expect_identical(unname(fit$argument[["a"]]), 2)
  expect_error(
    trust(.quadratic_objfn(c(a = 1)), c(a = 0), rinit = 1, rmax = 10,
          boundary = "nonsense"),
    "should be one of")
})


test_that("a control entry given flat reaches the objective, which refuses it", {
  obj <- constraintL2(mu = c(a = 1), sigma = 1)
  expect_error(trust(obj, c(a = 0), gtol = 1e-8), "unknown argument gtol")
  expect_error(trust(obj, c(a = 0), qnMemory = 3L), "unknown argument qnMemory")
  expect_error(trust(.quadratic_objfn(c(a = 1)), c(a = 0), stepControl = list(boundary = "clip")),
               "argument of trust")
})


test_that("a flat, high-value start makes progress instead of stopping at once", {
  # Large residuals with weak sensitivities: |f| is huge while |g| is small, so a
  # gradient criterion scaled by |f| would stop before the first step. gtol has
  # no relative counterpart.
  obj <- function(p, ...) {
    s <- 5e-4
    r <- as.numeric(p * s - 100)
    list(value = sum(r^2), gradient = 2 * s * r, hessian = 2 * s^2 * diag(length(p)))
  }
  init <- c(a = 0, b = 0)
  f0 <- obj(init)
  expect_lt(max(abs(f0$gradient)), 1e-5 * f0$value)   # the trap would be armed

  fit <- suppressWarnings(trust(obj, init, iterlim = 50))
  expect_gt(fit$iterations, 0L)
  expect_lt(fit$value, f0$value)
})


## ---- Interchangeable Hessian source (gn / bfgs / sr1) -----------------

# A flat objective: constant value, nonzero gradient. Every trial step is
# rejected with an unchanged value, so trust stalls deterministically: one of
# the hooks a handover fires on.
.flat_objfn <- function(d) {
  function(p, ...) list(value = 1.0, gradient = rep(1.0, d), hessian = diag(d))
}

# A quadratic reported with a Hessian `stiff` times too large. The model step is
# that much too short, so the value test fires long before the minimum: the
# stopped-short case a handover has to survive.
.misscaled_objfn <- function(target, stiff = 20) {
  d <- length(target)
  function(p, ...) {
    list(
      value    = 0.5 * sum((p - target)^2),
      gradient = p - target,
      hessian  = stiff * diag(d))
  }
}

# Rosenbrock: a curved valley. A quasi-Newton model mispredicts along it, so
# steps get rejected.
.rosenbrock_objfn <- function() {
  function(p, ...) {
    x <- p[[1]]; y <- p[[2]]
    gr <- c(-2 * (1 - x) - 400 * x * (y - x^2), 200 * (y - x^2))
    hs <- matrix(c(2 - 400 * y + 1200 * x^2, -400 * x, -400 * x, 200), 2, 2)
    names(gr) <- names(p); dimnames(hs) <- list(names(p), names(p))
    list(value = (1 - x)^2 + 100 * (y - x^2)^2, gradient = gr, hessian = hs)
  }
}

test_that("quasi-Newton methods reach the same optimum as gn on a nonlinear fit", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(pars = c(A = 1.0, k = 0.5), sigma = 0.02,
                         times = seq(0, 8, by = 0.5), seed = 17L)
  obj   <- normL2(data, bench$prd_id)
  init  <- c(A = 0.6, k = 0.9)

  ref <- trust(obj, init, rinit = 1, rmax = 10, iterlim = 200)
  expect_true(ref$converged)
  runs <- list(bfgs = list(hessianMethod = "bfgs"),
               sr1  = list(hessianMethod = "sr1"),
               fb   = list(phases = list(list(hessianMethod = "gn"), list(hessianMethod = "bfgs"))))
  for (hm in names(runs)) {
    fit <- do.call(trust, c(list(obj, init, rinit = 1, rmax = 10, iterlim = 200),
                            runs[[hm]]))
    expect_true(fit$converged, info = hm)
    expect_equal(fit$argument, ref$argument, tolerance = 1e-3, ignore_attr = TRUE,
                 info = hm)
  }
})

test_that("qnEval counts objective evaluations spent in the quasi-Newton phase", {
  target <- c(a = 1.0, b = -0.5)
  obj <- .quadratic_objfn(target)
  init <- c(a = 0, b = 0)
  expect_equal(trust(obj, init)$qnEval, 0L)                     # gn: never
  expect_gt(trust(obj, init, hessianMethod = "bfgs")$qnEval, 0L)
  expect_gt(trust(obj, init, hessianMethod = "sr1")$qnEval, 0L)
})

test_that("evalBySource and nSwitch account for every evaluation", {
  target <- c(a = 1.0, b = -0.5)
  obj    <- .quadratic_objfn(target)
  init   <- c(a = 0, b = 0)

  # "exact" is a fourth source: the objective's own Hessian, told apart from the
  # Gauss-Newton pass-through by which one the run asked for. A quadratic
  # objective has no `deriv2` argument, so it is not driven here.
  for (hm in c("gn", "bfgs", "sr1")) {
    fit <- trust(obj, init, hessianMethod = hm)
    expect_identical(names(fit$evalBySource), c("gn", "bfgs", "sr1", "exact", "sbfgs",
                                                "ssr1", "gnsbfgs", "gnssr1"), info = hm)
    expect_identical(fit$evalBySource[["exact"]], 0L, info = hm)
    expect_identical(sum(fit$evalBySource), fit$neval, info = hm)
    expect_identical(fit$nSwitch, 0L, info = hm)                 # no handover
    expect_gt(fit$evalBySource[[hm]], 0L)                        # charged to itself
  }
  expect_identical(trust(obj, init)$evalBySource[["bfgs"]], 0L)  # gn: only gn

  # Two phases hand over once and both are charged separately.
  hy <- suppressWarnings(trust(.flat_objfn(2L), init,
                               phases = list(list(hessianMethod = "gn"), list(hessianMethod = "bfgs"))))
  expect_identical(hy$nSwitch, 1L)
  expect_identical(sum(hy$evalBySource), hy$neval)
  expect_gt(hy$evalBySource[["gn"]], 0L)
  expect_gt(hy$evalBySource[["bfgs"]], 0L)
})

test_that("sr1 also learns from rejected steps, and gets there sooner for it", {
  obj  <- .rosenbrock_objfn()
  init <- c(x = -1.2, y = 1)
  fit  <- function(qr)
    trust(obj, init, rinit = 0.1, rmax = 10, iterlim = 2000,
          hessianMethod = "sr1", blather = TRUE,
          qnControl = list(hessianInit = "identity", qnRejected = qr))

  off <- fit(FALSE)
  on  <- fit(TRUE)

  # There are rejected steps to learn from.
  expect_gt(sum(off$accept == 0), 0L)
  for (f in list(off, on)) {
    expect_true(f$converged)
    expect_identical(f$stopReason, "gradient")
    expect_equal(unname(f$argument), c(1, 1), tolerance = 1e-6)
  }
  expect_lt(on$neval, off$neval)
})

test_that("qnRejected leaves gn and bfgs alone", {
  # BFGS updates only on accepted steps and gn forms no pair, so both are
  # bit-identical either way.
  obj  <- .rosenbrock_objfn()
  init <- c(x = -1.2, y = 1)
  for (hm in c("gn", "bfgs")) {
    args <- list(obj, init, rinit = 0.1, rmax = 10, iterlim = 2000,
                 hessianMethod = hm)
    off <- do.call(trust, c(args, list(qnControl = list(qnRejected = FALSE))))
    on  <- do.call(trust, c(args, list(qnControl = list(qnRejected = TRUE))))
    expect_identical(on$neval, off$neval, info = hm)
    expect_identical(on$value, off$value, info = hm)
    expect_identical(on$argument, off$argument, info = hm)
  }
})

test_that("blather reports the Hessian source per iteration", {
  target <- c(a = 1.0, b = -0.5, c = 2.3)
  obj <- .quadratic_objfn(target)
  init <- c(a = 0, b = 0, c = 0)
  expect_true(all(trust(obj, init, blather = TRUE)$hessianSource == "gn"))
  expect_true(all(trust(obj, init, hessianMethod = "bfgs",
                        blather = TRUE)$hessianSource == "bfgs"))
})

test_that("a bfgs phase takes over from gn on stagnation", {
  obj  <- .flat_objfn(2L)
  init <- c(a = 0, b = 0)

  gn <- suppressWarnings(trust(obj, init, hessianMethod = "gn", blather = TRUE))
  expect_true(all(gn$hessianSource == "gn"))

  hy <- suppressWarnings(trust(obj, init, blather = TRUE,
                               phases = list(list(hessianMethod = "gn"), list(hessianMethod = "bfgs"))))
  expect_identical(hy$hessianSource[1], "gn")   # always starts on gn
  expect_true("bfgs" %in% hy$hessianSource)      # and switches
  expect_gt(hy$qnEval, 0L)
})

test_that("the next phase takes over on a value stop too, and gets past where gn stops", {
  target <- c(a = 1.0, b = -0.5, c = 2.3)
  obj    <- .misscaled_objfn(target)
  init   <- c(a = 0, b = 0, c = 0)

  gn <- trust(obj, init, rinit = 1, rmax = 10, iterlim = 500,
              hessianMethod = "gn")
  expect_true(gn$stopReason %in% c("fvalue", "preddiff"))
  expect_true(gn$converged)
  expect_gt(max(abs(gn$argument - target)), 1e-3)   # nowhere near it

  hy <- trust(obj, init, rinit = 1, rmax = 10, iterlim = 500,
              phases = list(list(hessianMethod = "gn"), list(hessianMethod = "bfgs")),
              blather = TRUE)
  expect_true("bfgs" %in% hy$hessianSource)
  expect_identical(hy$stopReason, "gradient")       # first-order, not stopped short
  expect_lt(hy$value, gn$value)
  expect_equal(hy$argument, target, tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("a soft stop past the handover ends the run", {
  # Only one switch is available, so the run must terminate like bfgs afterwards
  # rather than run to iterlim.
  obj  <- .flat_objfn(2L)
  init <- c(a = 0, b = 0)
  hy <- suppressWarnings(trust(obj, init, iterlim = 500,
                               phases = list(list(hessianMethod = "gn"), list(hessianMethod = "bfgs"))))
  # the gradient the value never follows leaves a Newton decrease above mtol: not converged
  expect_false(hy$converged)
  expect_identical(hy$stopReason, "stagnation")
  expect_lt(hy$iterations, 500L)
})

test_that("hessianInit is inert for the methods that need the objective Hessian", {
  target <- c(a = 1.0, b = -0.5, c = 2.3)
  obj <- .quadratic_objfn(target)
  init <- c(a = 0, b = 0, c = 0)
  runs <- list(gn = list(), fb = list(phases = list(list(hessianMethod = "gn"), list(hessianMethod = "bfgs"))))
  for (hm in names(runs)) {
    args <- c(list(obj, init, rinit = 1), runs[[hm]])
    ref <- do.call(trust, args)
    alt <- do.call(trust, c(args, list(qnControl = list(hessianInit = "identity"))))
    expect_equal(alt$argument, ref$argument, ignore_attr = TRUE, info = hm)
    expect_identical(alt$neval, ref$neval, info = hm)
  }
})

test_that("a BFGS phase seeded from an indefinite Hessian starts positive definite", {
  # a double well in a, coupled to a quadratic in b, from a start of negative curvature
  obj <- function(x, ...) {
    a <- x[["a"]]; b <- x[["b"]]
    list(value = a^4 / 4 - a^2 / 2 + (b - 1)^2 + 0.3 * a * b,
         gradient = c(a = a^3 - a + 0.3 * b, b = 2 * (b - 1) + 0.3 * a),
         hessian = matrix(c(3 * a^2 - 1, 0.3, 0.3, 2), 2, dimnames = list(c("a", "b"), c("a", "b"))))
  }
  init <- c(a = 0.1, b = 0)
  runs <- list(
    seed = list(hessianMethod = "bfgs", qnControl = list(hessianInit = "gn")),
    handover = list(phases = list(list(hessianMethod = "gn", iterlim = 1), list(hessianMethod = "bfgs"))))
  for (k in names(runs)) {
    fit <- do.call(trust, c(list(obj, init, iterlim = 200), runs[[k]]))
    expect_true(fit$converged, info = k)
    expect_identical(fit$stopReason, "gradient", info = k)
    expect_gt(min(eigen(fit$hessian, symmetric = TRUE, only.values = TRUE)$values), 0)
  }
})

test_that("quasi-Newton methods require the reflective boundary and a known name", {
  obj <- .quadratic_objfn(c(a = 1.0))
  init <- c(a = 0)
  expect_error(trust(obj, init, boundary = "clip",
                     hessianMethod = "bfgs"), "reflective")
  expect_error(trust(obj, init, hessianMethod = "nope"))   # match.arg rejects
})


test_that("model tests do not end a quasi-Newton run", {

  fx    <- fx_decay_multicond_compiled()
  obj   <- normL2(fx_decay_data_multi(), fx$prd)
  start <- fx$outerpars + 3

  runs <- lapply(c(gn = "gn", bfgs = "bfgs", sr1 = "sr1"), function(m)
    trust(obj, start, rinit = 0.1, rmax = 10, iterlim = 300, hessianMethod = m,
          qnControl = list(hessianInit = if (m == "gn") "gn" else "identity")))

  # The model tests may end a run only while the model is the objective's own
  # Hessian. A quasi-Newton run is left with gradient, stagnation and radius.
  model <- c("fvalue", "preddiff", "step")
  expect_true(runs$gn$stopReason %in% model)
  expect_false(runs$bfgs$stopReason %in% model)
  expect_false(runs$sr1$stopReason %in% model)

  expect_true(all(vapply(runs, `[[`, TRUE, "converged")))
  for (fit in runs) expect_equal(fit$value, runs$gn$value, tolerance = 1e-4)

})


## ---- The iterate the result reports -----------------------------------

test_that("the reported iterate is the last accepted one under monotone acceptance", {
  obj  <- .rosenbrock_objfn()
  init <- c(x = -1.2, y = 1)
  for (hm in c("gn", "bfgs", "sr1")) {
    f   <- trust(obj, init, rinit = 0.1, rmax = 10, iterlim = 2000,
                 hessianMethod = hm, blather = TRUE)
    acc <- which(f$accept == 1)
    last <- acc[length(acc)]

    expect_identical(f$value, f$valtry[last], info = hm)
    expect_equal(unname(f$argument), unname(f$argtry[last, ]),
                 tolerance = 1e-14, info = hm)
    # The reported gradient has to belong to the reported argument.
    expect_equal(unname(f$gradient), unname(obj(f$argument)$gradient),
                 tolerance = 1e-12, info = hm)
  }
})

test_that("nonmonotone accepts an increase and still reports the best iterate", {
  obj  <- .rosenbrock_objfn()
  init <- c(x = -1.2, y = 1)
  for (hm in c("gn", "bfgs")) {
    f   <- trust(obj, init, rinit = 0.1, rmax = 10, iterlim = 2000,
                 hessianMethod = hm, stepControl = list(nonmonotone = 0.85),
                 blather = TRUE)
    acc <- which(f$accept == 1)

    # The mechanism under test: a step that raises the objective is taken.
    expect_gt(sum(f$valtry[acc] > f$valpath[acc]), 0)
    expect_identical(f$value, min(c(f$valpath[1], f$valtry[acc])), info = hm)
    expect_equal(unname(f$argument), c(1, 1), tolerance = 1e-4, info = hm)
  }
})

test_that("nonmonotone = 0 is the monotone rule and a positive eta changes the run", {
  obj  <- .rosenbrock_objfn()
  init <- c(x = -1.2, y = 1)
  run  <- function(...)
    trust(obj, init, rinit = 0.1, rmax = 10, iterlim = 2000, blather = TRUE, ...)

  mono <- run()
  expect_identical(run(stepControl = list(nonmonotone = 0))$valpath, mono$valpath)

  acc <- which(mono$accept == 1)
  expect_true(all(mono$valtry[acc] < mono$valpath[acc]))
  expect_false(identical(run(stepControl = list(nonmonotone = 0.85))$valpath,
                         mono$valpath))
})

test_that("nonmonotone outside [0, 1) is rejected", {
  obj <- .quadratic_objfn(c(a = 1.0))
  expect_error(trust(obj, c(a = 0), stepControl = list(nonmonotone = 1)),
               "nonmonotone")
  expect_error(trust(obj, c(a = 0), stepControl = list(nonmonotone = -0.1)),
               "nonmonotone")
})

test_that("an iteration-limit stop also reports the last accepted iterate", {
  # The iteration limit ends the run after an accepted step just as the soft
  # stops do, and it is the exit a quasi-Newton run reaches it through.
  obj  <- .rosenbrock_objfn()
  init <- c(x = -1.2, y = 1)
  for (hm in c("gn", "bfgs")) {
    f <- suppressWarnings(trust(obj, init, rinit = 0.1, rmax = 10, iterlim = 8,
                                hessianMethod = hm, blather = TRUE))
    expect_identical(f$stopReason, "iterlim", info = hm)
    acc <- which(f$accept == 1)
    expect_identical(f$value, f$valtry[acc[length(acc)]], info = hm)
  }
})


## ---- Phases --------------------------------------------------------------

test_that("a single phase or no switch reproduces the plain run", {
  obj  <- .misscaled_objfn(c(a = 1, b = -2))
  init <- c(a = 0, b = 0)
  plain <- trust(obj, init, iterlim = 200, hessianMethod = "gn")
  expect_identical(plain$phaseIterations, c(gn = plain$iterations))
  expect_null(plain$phaseSweep)
  for (alt in list(list(phases = list(list(hessianMethod = "gn"))),
                   list(phases = list(list(hessianMethod = "gn"),
                                      list(hessianMethod = "bfgs")),
                        switches = 0L))) {
    f <- do.call(trust, c(list(obj, init, iterlim = 200), alt))
    expect_identical(f$neval, plain$neval)
    expect_identical(f$value, plain$value)
    expect_identical(f$nSwitch, 0L)
    expect_identical(f$phaseIterations, plain$phaseIterations)
  }
})

test_that("the source of a later phase is selectable and is charged for the phase evaluations", {
  obj  <- .misscaled_objfn(c(a = 1, b = -2))
  init <- c(a = 0, b = 0)
  for (fb in c("bfgs", "sr1")) {
    f <- trust(obj, init, iterlim = 200, blather = TRUE,
               phases = list(list(hessianMethod = "gn"), list(hessianMethod = fb)))
    expect_identical(f$nSwitch, 1L, info = fb)
    expect_true(fb %in% f$hessianSource, info = fb)
    expect_gt(f$evalBySource[[fb]], 0)
    expect_identical(sum(f$evalBySource), f$neval, info = fb)
    expect_identical(f$phaseSweep, c(gn = "forward", setNames("forward", fb)))
  }
})

test_that("switches beyond the list cycle back to the first phase", {
  # A flat objective stalls under either source, so every allowed handover
  # fires and the sources alternate.
  obj  <- .flat_objfn(2L)
  init <- c(a = 0, b = 0)
  ph   <- list(list(hessianMethod = "gn"), list(hessianMethod = "bfgs"))
  one <- trust(obj, init, iterlim = 300, blather = TRUE, phases = ph)
  two <- trust(obj, init, iterlim = 300, blather = TRUE, phases = ph, switches = 2L)

  expect_identical(one$nSwitch, 1L)
  expect_identical(two$nSwitch, 2L)
  expect_named(two$phaseIterations, c("gn", "bfgs", "gn"))
  # Back on gn after the second handover, and the reseed evaluation is charged.
  expect_identical(utils::tail(two$hessianSource, 1L), "gn")
  expect_identical(sum(two$evalBySource), two$neval)
})

test_that("phases are checked by name and the fallback arguments refused", {
  obj  <- .quadratic_objfn(c(a = 1.0))
  init <- c(a = 0)
  expect_error(trust(obj, init, hessianMethod = "nope"), "should be one of")
  expect_error(trust(obj, init, phases = list(list(hessianMethod = "nope"))),
               "should be one of")
  expect_error(trust(obj, init, phases = list(list(hessianMethod = "gn", budget = 3))),
               "unknown entries budget")
  expect_error(trust(obj, init, phases = list(list(iterlim = 3))), "no hessianMethod")
  expect_error(trust(obj, init, phases = list(list(hessianMethod = "gn", iterlim = 0))),
               "positive")
  expect_error(trust(obj, init, phases = list(list(hessianMethod = "gn")), switches = -1),
               "non-negative")
  expect_error(trust(obj, init, hessianMethod = "sr1",
                     phases = list(list(hessianMethod = "gn"))), "every phase names the hessianMethod")
  expect_error(trust(obj, init, phases = list(list(hessianMethod = "gn")), sweep = "reverse"),
               "every phase sets the sweep")
  expect_error(trust(obj, init, boundary = "clip",
                     phases = list(list(hessianMethod = "gn"), list(hessianMethod = "gn"))),
               "single \"gn\" phase")
})

test_that("a phase reports the sweep it used and refuses one it cannot run", {
  # A plain R function has no `sweep` argument, so every phase runs forward.
  obj  <- .misscaled_objfn(c(a = 1, b = -2))
  init <- c(a = 0, b = 0)
  f <- trust(obj, init, blather = TRUE,
             phases = list(gn = list(hessianMethod = "gn"),
                           qn = list(hessianMethod = "sr1")))
  expect_identical(f$phaseSweep, c(gn = "forward", qn = "forward"))
  expect_named(f$phaseIterations, c("gn", "qn"))
  expect_error(trust(obj, init, phases = list(list(hessianMethod = "gn"),
                                              list(hessianMethod = "sr1", sweep = "reverse"))),
               "no `sweep` argument")
  expect_error(trust(obj, init, phases = list(list(hessianMethod = "sr1"),
                                              list(hessianMethod = "gn", sweep = "reverse"))),
               "forward sensitivities")
})

test_that("a phase that uses up the phase iterations ends unconverged", {
  obj  <- .misscaled_objfn(c(a = 1, b = -2))
  init <- c(a = 0, b = 0)
  ph   <- list(list(hessianMethod = "gn"), list(hessianMethod = "bfgs", iterlim = 2L))
  f <- trust(obj, init, blather = TRUE, phases = ph)
  expect_identical(f$nSwitch, 1L)
  expect_identical(f$stopReason, "iterlim")
  expect_false(f$converged)
  expect_identical(sum(f$hessianSource == "bfgs"), 2L)
  expect_identical(f$phaseIterations[["bfgs"]], 2L)

  # With a handover left, the budget hands back instead.
  g <- trust(obj, init, blather = TRUE, phases = ph, switches = 2L)
  expect_identical(g$nSwitch, 2L)
  expect_identical(utils::tail(g$hessianSource, 1L), "gn")
})

test_that("a phase hands over after radiusStall iterations at an unchanged radius", {
  # every step ends on the boundary of the region, so the radius grows to rmax and stays there
  obj  <- function(p, ...) list(value = 0.5 * sum((p - 10)^2), gradient = p - 10, hessian = diag(2))
  init <- c(a = 0, b = 0)
  ph   <- list(list(hessianMethod = "gn", radiusStall = 3), list(hessianMethod = "bfgs"))
  f <- trust(obj, init, rinit = 0.1, rmax = 0.1, blather = TRUE, phases = ph)
  expect_identical(f$phaseIterations[["gn"]], 3L)
  expect_identical(f$hessianSource[1:4], c("gn", "gn", "gn", "bfgs"))
  expect_true(f$converged)

  # a change of the radius restarts the count
  g <- trust(obj, init, rinit = 0.1, rmax = 0.4, blather = TRUE, phases = ph)
  expect_equal(g$r[1:6], c(0.1, 0.2, 0.4, 0.4, 0.4, 0.1))
  expect_identical(g$phaseIterations[["gn"]], 5L)

  # without a handover left the count stops nothing
  h <- trust(obj, init, rinit = 0.1, rmax = 0.1, blather = TRUE, phases = ph, switches = 0)
  expect_identical(names(h$phaseIterations), "gn")
  expect_true(h$converged)
  expect_error(trust(obj, init, phases = list(list(hessianMethod = "gn", radiusStall = 0))),
               "radiusStall must be a positive number")
})

## ---- Default: the hybrid of Fides -------------------------------------------

# The misscaled quadratic with names on the Hessian, in a dense and a sparse version.
.misscaled_named <- function(target, sparse = FALSE) {
  f <- .misscaled_objfn(target)
  function(p, ..., hessian = TRUE) {
    out <- f(p)
    dimnames(out$hessian) <- list(names(p), names(p))
    if (sparse && identical(hessian, "sparse")) out$hessian <- Matrix::Matrix(out$hessian, sparse = TRUE)
    out
  }
}

test_that("the sparse choice evaluates the reflective start, which the kernel then reuses", {
  target <- c(a = 1, b = -2, c = 3)
  seen <- list()
  obj <- function(p, ..., hessian = TRUE) {
    seen[[length(seen) + 1L]] <<- p
    H <- diag(2, length(p))
    dimnames(H) <- list(names(p), names(p))
    if (identical(hessian, "sparse")) H <- Matrix::Matrix(H, sparse = TRUE)
    list(value = sum((p - target)^2), gradient = 2 * (p - target), hessian = H)
  }
  fit <- withr::with_options(list(dMod2.sparseAuto = list(parameters = 0, fill = 1)),
    trust(obj, c(a = 0, b = 0, c = 0), parlower = c(a = 0), iterlim = 50))
  expect_true(fit$converged)
  expect_gt(seen[[1]][["a"]], 0)
  expect_gt(max(abs(seen[[2]] - seen[[1]])), 1e-8)
})

# The damped BFGS update of the kernel with the cautious test.
.bfgs_r <- function(B, s, y, cautious = 1e-8) {
  if (!(sum(s * y) > cautious * sqrt(sum(s^2) * sum(y^2)))) return(B)
  Bs <- as.vector(B %*% s)
  sy <- sum(s * y); sBs <- sum(s * Bs)
  th <- if (sy < 0.2 * sBs) 0.8 * sBs / (sBs - sy) else 1
  r <- th * y + (1 - th) * Bs
  B + tcrossprod(r) / sum(s * r) - tcrossprod(Bs) / sBs
}

test_that("the default runs Gauss-Newton, then BFGS after 20 iterations at an unchanged radius", {
  obj  <- .misscaled_named(c(a = 1, b = -2, c = 3))
  init <- c(a = 0, b = 0, c = 0)
  fit <- trust(obj, init, iterlim = 200, blather = TRUE)
  expect_true(fit$converged)
  expect_identical(fit$nSwitch, 1L)
  expect_identical(names(fit$phaseIterations), c("gn", "bfgs"))
  n <- fit$phaseIterations[["gn"]]
  expect_identical(fit$hessianSource, rep(c("gn", "bfgs"), fit$phaseIterations))
  # the radius at the head of the last 20 Gauss-Newton iterations is the one at the switch
  expect_true(all(fit$r[(n - 19):(n + 1)] == fit$r[n + 1]))
  expect_false(fit$r[n - 20] == fit$r[n + 1])
  expect_equal(fit$argument, c(a = 1, b = -2, c = 3), tolerance = 1e-6)
})

test_that("the tracked BFGS matrix is the start Hessian updated with every accepted pair", {
  obj  <- .misscaled_named(c(a = 1, b = -2, c = 3))
  init <- c(a = 0, b = 0, c = 0)
  fit <- trust(obj, init, iterlim = 200, blather = TRUE)
  expect_identical(fit$nSwitch, 1L)
  B <- obj(init)$hessian
  for (i in which(fit$accept == 1)) {
    s <- fit$argtry[i, ] - fit$argpath[i, ]
    B <- .bfgs_r(B, s, obj(fit$argtry[i, ])$gradient - obj(fit$argpath[i, ])$gradient)
  }
  expect_equal(fit$hessian, B, tolerance = 1e-10, ignore_attr = TRUE)
  expect_false(isTRUE(all.equal(B, obj(init)$hessian, check.attributes = FALSE)))

  # an indefinite start Hessian enters with the absolute values of the eigenvalues
  ind <- function(p, ...) {
    out <- obj(p)
    out$hessian[] <- diag(c(20, -20, 20))
    out
  }
  expect_warning(f2 <- trust(ind, init, iterlim = 3, blather = TRUE, phases = list(
    list(hessianMethod = "gn", iterlim = 3),
    list(hessianMethod = "bfgs", qnControl = list(hessianInit = "tracked")))), "Maximum number")
  B <- diag(20, 3)
  for (i in which(f2$accept == 1)) {
    s <- f2$argtry[i, ] - f2$argpath[i, ]
    B <- .bfgs_r(B, s, ind(f2$argtry[i, ])$gradient - ind(f2$argpath[i, ])$gradient)
  }
  expect_identical(f2$nSwitch, 1L)
  expect_equal(f2$hessian, B, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("keepRadius keeps the radius at the handover, and onStop = \"end\" ends the run", {
  obj  <- .misscaled_named(c(a = 1, b = -2, c = 3))
  init <- c(a = 0, b = 0, c = 0)
  ph <- function(keep) list(list(hessianMethod = "gn", radiusStall = 20, onStop = "end"),
                            list(hessianMethod = "bfgs", keepRadius = keep,
                                 qnControl = list(hessianInit = "tracked")))
  kept <- trust(obj, init, iterlim = 200, blather = TRUE, phases = ph(TRUE))
  restart <- trust(obj, init, iterlim = 200, blather = TRUE, phases = ph(FALSE))
  n <- kept$phaseIterations[["gn"]]
  expect_identical(restart$phaseIterations[["gn"]], n)
  expect_gt(kept$r[n], 0.1)
  expect_identical(kept$r[n + 1], kept$r[n])
  expect_identical(restart$r[n + 1], 0.1)

  # a value stop of the first phase hands over by default and ends the run with "end"
  flat <- function(p, ...) list(value = 0.5 * sum(c(1, 1e-3) * p^2), gradient = c(1, 1e-3) * p,
                                hessian = diag(c(1, 1e-3)))
  run <- function(stop) trust(flat, c(x = 0, y = 9), rinit = 0.05, rmax = 0.05,
                              tolControl = list(ftol = 1e-3, mtol = 1e-4, gtol = 1e-6),
                              phases = list(list(hessianMethod = "gn", onStop = stop),
                                            list(hessianMethod = "bfgs")))
  expect_identical(run("handover")$nSwitch, 1L)
  ended <- run("end")
  expect_identical(ended$nSwitch, 0L)
  expect_identical(ended$stopReason, "fvalue")
  expect_true(ended$converged)

  expect_error(trust(obj, init, phases = list(list(hessianMethod = "bfgs",
                                                   qnControl = list(hessianInit = "tracked")))),
               "after the first")
  expect_error(trust(obj, init, phases = list(list(hessianMethod = "gn"),
                                              list(hessianMethod = "sr1",
                                                   qnControl = list(hessianInit = "tracked")))),
               "after the first")
  expect_error(trust(obj, init, phases = list(list(hessianMethod = "gn", keepRadius = NA))),
               "TRUE or FALSE")
  expect_error(trust(obj, init, phases = list(list(hessianMethod = "gn", onStop = "nope"))),
               "should be one of")
})

test_that("explicit gn, the sparse path and the clip boundary stay Gauss-Newton", {
  target <- c(a = 1, b = -2, c = 3)
  init <- c(a = 0, b = 0, c = 0)
  dense <- .misscaled_named(target)
  expect_identical(names(trust(dense, init, iterlim = 200)$phaseIterations), c("gn", "bfgs"))
  gn <- trust(dense, init, iterlim = 200, hessianMethod = "gn")
  expect_identical(gn$phaseIterations, c(gn = gn$iterations))
  expect_null(gn$phaseSweep)

  sp <- .misscaled_named(target, sparse = TRUE)
  forced <- trust(sp, init, iterlim = 200, sparse = TRUE)
  auto <- withr::with_options(list(dMod2.sparseAuto = list(parameters = 0, fill = 1)),
                              trust(sp, init, iterlim = 200))
  for (f in list(forced, auto)) {
    expect_s4_class(f$hessian, "sparseMatrix")
    expect_identical(f$phaseIterations, c(gn = f$iterations))
    expect_identical(f$iterations, gn$iterations)
  }
  expect_error(trust(sp, init, sparse = TRUE, phases = list(
    list(hessianMethod = "gn"), list(hessianMethod = "bfgs", qnControl = list(hessianInit = "tracked")))),
    "limited-memory BFGS")
  clip <- trust(dense, init, iterlim = 200, boundary = "clip")
  expect_true(clip$converged)
  expect_null(clip$phaseIterations)
})

## ---- Counters across repeated handovers -------------------------------------

# The flat objective with a `sweep` argument; `log()` lists the sweep and the
# Hessian request of every evaluation.
.logged_flat <- function(d) {
  log <- character(0)
  list(f = function(p, ..., sweep = "forward") {
         log <<- c(log, paste(sweep, isTRUE(list(...)$hessian)))
         list(value = 1.0, gradient = rep(1.0, d), hessian = diag(d))
       },
       log = function() log)
}

test_that("phase iterations add up, and the handover iteration counts on the phase it ends", {
  obj  <- .flat_objfn(2L)
  init <- c(a = 0, b = 0)
  ph   <- list(list(hessianMethod = "gn"), list(hessianMethod = "bfgs"))
  for (s in 1:3) {
    f <- suppressWarnings(trust(obj, init, phases = ph, switches = s, blather = TRUE))
    expect_identical(f$nSwitch, s)
    expect_identical(sum(f$phaseIterations), f$iterations)
    expect_identical(rle(f$hessianSource)$lengths, unname(f$phaseIterations))
    expect_identical(f$stopReason, "stagnation")
    expect_false(f$converged)
  }

  # A budget on the first phase ends it after exactly that many iterations.
  g <- suppressWarnings(trust(obj, init, phases = list(list(hessianMethod = "gn", iterlim = 3),
                                                       list(hessianMethod = "bfgs"))))
  expect_identical(g$phaseIterations[["gn"]], 3L)
  expect_identical(g$nSwitch, 1L)
  expect_false(g$converged)
})

test_that("a budget run out after the last handover ends unconverged on either source", {
  obj  <- .flat_objfn(2L)
  init <- c(a = 0, b = 0)
  ph   <- list(list(hessianMethod = "gn", iterlim = 2), list(hessianMethod = "bfgs", iterlim = 2))
  for (s in 1:3) {
    f <- trust(obj, init, phases = ph, switches = s, blather = TRUE)
    expect_identical(f$nSwitch, s)
    expect_identical(unname(f$phaseIterations), rep(2L, s + 1L))
    expect_identical(f$stopReason, "iterlim")
    expect_false(f$converged)
    expect_identical(utils::tail(f$hessianSource, 1L), c("gn", "bfgs")[s %% 2 + 1])
  }
})

test_that("the total iterlim ends the run inside a later phase", {
  obj  <- .flat_objfn(2L)
  init <- c(a = 0, b = 0)
  expect_warning(
    f <- trust(obj, init, iterlim = 7L,
               phases = list(list(hessianMethod = "gn"), list(hessianMethod = "bfgs"))),
    "Maximum number of iterations")
  expect_identical(f$phaseIterations, c(gn = 5L, bfgs = 2L))
  expect_identical(f$stopReason, "iterlim")
  expect_false(f$converged)
})

test_that("every evaluation is charged once, with the sweep of the phase", {
  lf   <- .logged_flat(2L)
  init <- c(a = 0, b = 0)
  f <- suppressWarnings(trust(lf$f, init, blather = TRUE, switches = 3L,
                              phases = list(list(hessianMethod = "bfgs", sweep = "forward"),
                                            list(hessianMethod = "sr1"))))
  expect_identical(f$phaseSweep, c(bfgs = "forward", sr1 = "reverse"))
  expect_named(f$phaseIterations, c("bfgs", "sr1", "bfgs", "sr1"))
  expect_identical(rle(f$hessianSource)$lengths, unname(f$phaseIterations))
  # One evaluation more than neval: the probe of the reverse sweep.
  expect_length(lf$log(), f$neval + 1L)
  expect_identical(sum(f$evalBySource), f$neval)
  # A Hessian is always asked for forward, and the reverse phase evaluates backwards.
  expect_false("reverse TRUE" %in% lf$log())
  expect_true("reverse FALSE" %in% lf$log()[-1L])
})

# Drops by one every eighth evaluation, wherever it is evaluated, so a phase
# alternates an accepted step with a stagnation.
.stepping_objfn <- function(d) {
  n <- 0L
  function(p, ...) {
    n <<- n + 1L
    list(value = 1000 - n %/% 8L, gradient = rep(1.0, d), hessian = diag(d))
  }
}

test_that("each quasi-Newton phase has a separate reseed allowance", {
  init <- c(a = 0, b = 0)
  qn   <- list(hessianReseed = "stall")
  one <- trust(.stepping_objfn(2L), init,
               phases = list(list(hessianMethod = "sr1", qnControl = qn)))
  two <- trust(.stepping_objfn(2L), init,
               phases = list(list(hessianMethod = "sr1", qnControl = qn),
                             list(hessianMethod = "bfgs", qnControl = qn)))
  expect_identical(one$nReseed, 8L)
  expect_identical(one$stopReason, "stagnation")
  expect_identical(two$nSwitch, 1L)
  expect_identical(two$nReseed, 16L)
  expect_identical(sum(two$evalBySource), two$neval)
})

# An objective may decline to build a Hessian, and NULL is how it says so.
# `give` decides at which points it returns one.
.declining_objfn <- function(give = function(p) FALSE) {
  function(p, hessian = TRUE, ...) {
    x <- p[[1]]; y <- p[[2]]
    gr <- c(-2 * (1 - x) - 400 * x * (y - x^2), 200 * (y - x^2))
    hs <- matrix(c(2 - 400 * y + 1200 * x^2, -400 * x, -400 * x, 200), 2, 2)
    names(gr) <- names(p); dimnames(hs) <- list(names(p), names(p))
    list(value = (1 - x)^2 + 100 * (y - x^2)^2, gradient = gr,
         hessian = if (isTRUE(hessian) && give(p)) hs else NULL)
  }
}

test_that("a quasi-Newton start survives an objective that returns no Hessian", {
  init <- c(x = -1.2, y = 1)
  expect_warning(
    fit <- trust(.declining_objfn(), init, rinit = 1, rmax = 10, iterlim = 200,
                 hessianMethod = "bfgs"),
    "no Hessian at parinit")
  expect_true(fit$converged)
  expect_lt(fit$value, 1e-8)
})

test_that("gn says which argument it cannot honor when no Hessian arrives", {
  init <- c(x = -1.2, y = 1)
  expect_error(trust(.declining_objfn(), init, rinit = 1, rmax = 10),
               "no Hessian at parinit")
  expect_error(trust(.declining_objfn(), init, rinit = 1, rmax = 10,
                     boundary = "clip"),
               "no Hessian at parinit")
})

test_that("a trial point without a Hessian is a failed evaluation, not an error", {
  # A Hessian at the start only: the kernel must treat the trial as infeasible
  # rather than let the conversion throw from outside the kernel handler.
  init <- c(x = -1.2, y = 1)
  expect_warning(
    fit <- trust(.declining_objfn(function(p) isTRUE(all.equal(unname(p), unname(init)))), init, rinit = 1,
                 rmax = 10, iterlim = 20),
    "evaluation failed")
  expect_identical(fit$stopReason, "objfun")
  expect_false(fit$converged)
})


## ---- Structured quasi-Newton methods --------------------------------

# A least-squares objective from residuals: value sum(r^2), gradient 2 J'r, Gauss-Newton Hessian
# 2 J'J, and with deriv2 the exact Hessian. `res(p)` returns r, J and the residual Hessians Hr.
.lsq_objfn <- function(res) {
  function(p, deriv = TRUE, hessian = TRUE, deriv2 = FALSE, ...) {
    e <- res(p); n <- names(p)
    H <- 2 * crossprod(e$J)
    if (deriv2) H <- H + 2 * Reduce(`+`, Map(`*`, e$r, e$Hr))
    dimnames(H) <- list(n, n)
    objlist(sum(e$r^2), setNames(as.vector(2 * crossprod(e$J, e$r)), n), H)
  }
}

# The same with the data rows of the residuals, as objlist() describes them: Jacobian J, seeds 2 r
# and the Gauss-Newton Hessian, under hessian = "structured".
.lsq_rows_objfn <- function(res) {
  base <- .lsq_objfn(res)
  function(p, deriv = TRUE, hessian = TRUE, ...) {
    out <- base(p, deriv = deriv, hessian = hessian, ...)
    if (identical(hessian, "structured")) {
      e <- res(p)
      out$structure <- list(jacobian = unname(e$J), seeds = 2 * e$r, hessian = out$hessian)
    }
    out
  }
}

# Residuals (x + 1, l x^2 + x - 1) per coordinate: minimum at 0 with value 4 for every l, where
# the second-order term S = -2 l is not small against J'J = 2.
.ds_objfn <- function(l) .lsq_objfn(function(p) {
  x <- p[[1]]; y <- p[[2]]
  Z <- matrix(0, 2, 2)
  list(r  = c(x + 1, l * x^2 + x - 1, y + 1, l * y^2 + y - 1),
       J  = rbind(c(1, 0), c(2 * l * x + 1, 0), c(0, 1), c(0, 2 * l * y + 1)),
       Hr = list(Z, matrix(c(2 * l, 0, 0, 0), 2), Z, matrix(c(0, 0, 0, 2 * l), 2)))
})

# Exponential decay fitted to exact data: residual zero at (a, k) = (1, 0.5).
.exp_objfn <- local({
  tt <- seq(0, 4, by = 0.5)
  .lsq_objfn(function(p) {
    a <- p[[1]]; k <- p[[2]]; e <- exp(-k * tt)
    list(r = a * e - exp(-0.5 * tt), J = cbind(e, -a * tt * e),
         Hr = lapply(seq_along(tt), function(i)
           tt[i] * e[i] * matrix(c(0, -1, -1, a * tt[i]), 2)))
  })
})

.structured <- c("sbfgs", "ssr1", "gnsbfgs", "gnssr1")
.tight <- list(gtol = 1e-6, ftol = 0, mtol = 0)

test_that("every Hessian method reaches the minimum of a large-residual problem", {
  obj <- .ds_objfn(-0.8)
  for (hm in c("gn", "bfgs", "sr1", "exact", .structured)) {
    fit <- trust(obj, c(x = 1, y = -1), rinit = 1, iterlim = 200, hessianMethod = hm,
                 tolControl = .tight)
    expect_true(fit$converged, info = hm)
    expect_equal(unname(fit$argument), c(0, 0), tolerance = 1e-6, info = hm)
    expect_equal(fit$value, 4, tolerance = 1e-10, info = hm)
    expect_gt(fit$evalBySource[[hm]], 0L)
  }
})

test_that("structured methods need fewer iterations than gn where the residual is large", {
  obj <- .ds_objfn(-0.8)
  it <- vapply(c("gn", .structured), function(hm)
    trust(obj, c(x = 1, y = -1), rinit = 1, iterlim = 200, hessianMethod = hm,
          tolControl = .tight)$iterations, 1L)
  expect_gt(it[["gn"]], 40L)
  expect_true(all(it[.structured] <= 15L), info = paste(names(it), it, collapse = " "))
})

test_that("the hybrid methods and residual scaling keep Gauss-Newton speed at zero residual", {
  init <- c(a = 3, k = 2)
  run <- function(hm, ...) trust(.exp_objfn, init, rinit = 1, iterlim = 200, hessianMethod = hm,
                                 tolControl = .tight, ...)
  gn <- run("gn")
  expect_true(gn$converged)
  fast <- list(run("gnsbfgs"), run("gnssr1"),
               run("sbfgs", qnControl = list(qnScaling = "residual")),
               run("ssr1", qnControl = list(qnScaling = "residual")))
  for (f in fast) {
    expect_true(f$converged)
    expect_equal(unname(f$argument), c(1, 0.5), tolerance = 1e-6)
    expect_lte(f$iterations, gn$iterations + 1L)
  }
  # without the scaling a structured correction remains and slows the end
  expect_gt(run("sbfgs")$iterations, gn$iterations + 2L)
  expect_gt(fast[[1]]$hybridGN, 0L)
})

test_that("the hybrid switch follows the ratio z's / s's against qnHybridTol", {
  # z = (g+ - g - G+ s) sqrt(f+ / f); "gnsbfgs" takes the Gauss-Newton branch where z's / s's
  # is below the tolerance, "gnssr1" where |z's| / s's is
  for (hm in c("gnsbfgs", "gnssr1")) for (tol in c(1e-6, 1e-2)) {
    fit <- trust(.exp_objfn, c(a = 3, k = 2), rinit = 1, iterlim = 200, hessianMethod = hm,
                 tolControl = .tight, qnControl = list(qnHybridTol = tol), blather = TRUE)
    acc <- which(fit$accept == 1)
    expect_gt(length(acc), 2L)
    branch <- vapply(acc, function(i) {
      x0 <- fit$argpath[i, ]; x1 <- fit$argtry[i, ]
      o0 <- .exp_objfn(x0); o1 <- .exp_objfn(x1)
      s <- x1 - x0
      z <- (o1$gradient - o0$gradient - o1$hessian %*% s) * sqrt(o1$value / o0$value)
      q <- sum(z * s) / sum(s * s)
      if (hm == "gnssr1") q <- abs(q)
      if (q < tol) "gn" else hm
    }, "")
    inner <- acc < fit$iterations
    expect_identical(fit$hessianSource[acc[inner] + 1L], branch[inner], info = paste(hm, tol))
    expect_identical(fit$hybridGN, sum(branch == "gn"), info = paste(hm, tol))
  }
  # a negative curvature term keeps "gnsbfgs" on Gauss-Newton, "gnssr1" takes it
  neg <- .ds_objfn(0.5)
  b <- trust(neg, c(x = 1, y = 1), rinit = 1, hessianMethod = "gnsbfgs", tolControl = .tight)
  s <- trust(neg, c(x = 1, y = 1), rinit = 1, hessianMethod = "gnssr1", tolControl = .tight)
  expect_gte(b$hybridGN, b$iterations - 2L)
  expect_identical(s$hybridGN, 0L)
  expect_lt(s$iterations, b$iterations)
})

test_that("qnHybridShift adds the residual norm to the Gauss-Newton branch", {
  neg <- .ds_objfn(0.5)
  fit <- trust(neg, c(x = 1, y = 1), rinit = 1, hessianMethod = "gnsbfgs", tolControl = .tight,
               qnControl = list(qnHybridShift = TRUE), blather = TRUE)
  expect_true(fit$converged)
  gn <- which(fit$hessianSource == "gn")
  expect_gt(length(gn), 2L)
  for (i in gn) {
    x0 <- fit$argpath[i, ]; s <- fit$argtry[i, ] - x0; o <- neg(x0)
    m <- sum(o$gradient * s) + 0.5 * sum(s * ((o$hessian + sqrt(o$value) * diag(2)) %*% s))
    expect_equal(fit$preddiff[i], m, tolerance = 1e-8, info = i)
  }
})

test_that("a structured SR1 pair orthogonal to the step is skipped", {
  # f = p'Qp / 2 reported with G = Q + D, so the structured secant is w = -D s, orthogonal to
  # a first step s = (-1, -1)
  Q <- c(2, 3); D <- c(1, -1)
  obj <- function(p, ...) {
    n <- names(p)
    objlist(sum(Q * p^2) / 2, setNames(Q * p, n), matrix(diag(Q + D), 2, 2, dimnames = list(n, n)))
  }
  one <- function(init) suppressWarnings(
    trust(obj, init, rinit = 2, iterlim = 1L, hessianMethod = "ssr1", blather = TRUE))
  skip <- one(c(a = 1.5, b = 2 / 3))
  expect_equal(unname(skip$argtry[1, ] - skip$argpath[1, ]), c(-1, -1), tolerance = 1e-12)
  expect_identical(skip$qnSkipped, 1L)
  expect_equal(unname(skip$hessian), diag(Q + D), tolerance = 1e-12)
  # off the orthogonal direction the pair updates G by w w' / w's
  upd <- one(c(a = 1.5, b = 1))
  s <- upd$argtry[1, ] - upd$argpath[1, ]
  w <- -D * s
  expect_identical(upd$qnSkipped, 0L)
  expect_equal(unname(upd$hessian), diag(Q + D) + tcrossprod(w) / sum(w * s), tolerance = 1e-10)
})

test_that("structured methods run as phases and take over from any source", {
  obj  <- .ds_objfn(-0.8)
  init <- c(x = 1, y = -1)
  gn <- trust(obj, init, rinit = 1, tolControl = .tight)
  for (first in c("gn", "bfgs", "exact")) for (hm in .structured) {
    fit <- trust(obj, init, rinit = 1, tolControl = .tight,
                 phases = list(list(hessianMethod = first, iterlim = 3), list(hessianMethod = hm)))
    info <- paste(first, hm)
    expect_true(fit$converged, info = info)
    expect_identical(fit$nSwitch, 1L, info = info)
    expect_identical(names(fit$phaseIterations), c(first, hm), info = info)
    expect_gt(fit$evalBySource[[hm]], 0L)
    expect_identical(sum(fit$evalBySource), fit$neval, info = info)
    expect_lt(fit$iterations, gn$iterations)
    expect_equal(unname(fit$argument), c(0, 0), tolerance = 1e-6, info = info)
  }
})

test_that("structured methods honor bounds and maximization", {
  obj <- .ds_objfn(-0.8)
  for (hm in .structured) {
    fit <- trust(obj, c(x = 1, y = -1), rinit = 1, hessianMethod = hm, tolControl = .tight,
                 parlower = c(x = 0.2))
    expect_true(fit$converged, info = hm)
    expect_equal(fit$argument[["x"]], 0.2, tolerance = 1e-6, info = hm)
    expect_true(fit$atBound[["x"]], info = hm)
    neg <- function(p, ...) { o <- obj(p, ...); objlist(-o$value, -o$gradient, -o$hessian) }
    mx <- trust(neg, c(x = 1, y = -1), rinit = 1, hessianMethod = hm, tolControl = .tight,
                minimize = FALSE)
    expect_equal(unname(mx$argument), c(0, 0), tolerance = 1e-6, info = hm)
  }
})

test_that("structured methods refuse settings they cannot honor", {
  obj <- .ds_objfn(-0.8)
  init <- c(x = 1, y = -1)
  expect_error(trust(obj, init, hessianMethod = "sbfgs", qnControl = list(qnMemory = 3)),
               "keeps no pairs")
  expect_error(trust(obj, init, hessianMethod = "ssr1",
                     qnControl = list(hessianInit = "identity")), "Gauss-Newton")
  expect_error(trust(obj, init, phases = list(list(hessianMethod = "gnsbfgs", sweep = "reverse"))),
               "forward sensitivities")
  expect_error(trust(obj, init, hessianMethod = "sbfgs", sparse = TRUE), "sparse")
  expect_error(trust(obj, init, hessianMethod = "sbfgs", boundary = "clip"), "reflective")
  expect_error(trust(obj, init, hessianMethod = "gnsbfgs", qnControl = list(qnScaling = "x")))
  expect_error(trustL1(obj, init, hessianMethod = "sbfgs"), "argument of trust")
})

test_that("structured methods fit a sum of normL2 and a prior to the gn optimum", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  data  <- fx_decay_data(pars = c(A = 1.0, k = 0.5), sigma = 0.02,
                         times = seq(0, 8, by = 0.5), seed = 17L)
  obj   <- normL2(data, bench$prd_id) + constraintL2(mu = c(A = 0.8, k = 0.4), sigma = 0.1)
  init  <- c(A = 0.6, k = 0.9)
  ref <- trust(obj, init, rinit = 1, iterlim = 200)
  for (hm in .structured) {
    fit <- trust(obj, init, rinit = 1, iterlim = 200, hessianMethod = hm)
    expect_true(fit$converged, info = hm)
    expect_equal(fit$argument, ref$argument, tolerance = 1e-4, ignore_attr = TRUE, info = hm)
  }
})

test_that("the structured secant of an objective with data rows is (J+ - J)'u+", {
  # one accepted SR1 step from A = 0: the model is G+ + w w' / w's with w the secant
  res <- function(p) {
    a <- p[[1]]; k <- p[[2]]; tt <- c(0.5, 1, 2, 3); e <- exp(-k * tt)
    list(r = a * e - c(2, 0.5, 1.5, 0.2), J = cbind(e, -a * tt * e))
  }
  obj <- .lsq_rows_objfn(res)
  class(obj) <- c("objfn", "fn")
  fit <- suppressWarnings(trust(obj, c(a = 1, k = 0.3), rinit = 0.2, iterlim = 1L,
                                hessianMethod = "ssr1", blather = TRUE))
  expect_true(fit$accept)
  x0 <- fit$argpath[1, ]; x1 <- fit$argtry[1, ]
  w <- drop(crossprod(res(x1)$J - res(x0)$J, 2 * res(x1)$r))
  o1 <- obj(x1)
  expect_equal(unname(fit$hessian), unname(o1$hessian + tcrossprod(w) / sum(w * (x1 - x0))),
               tolerance = 1e-10)
})

test_that("normL2 returns the data rows, and the secant of a sum with a prior is that of the data", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  x <- fx_extra()
  data <- fx_decay_data(pars = c(A = 1.0, k = 0.5), sigma = 0.02, times = seq(0, 8, by = 0.5),
                        seed = 17L)
  data[["C1"]]$sigma <- NA
  data[["C1"]]$value <- data[["C1"]]$value * rep_len(c(1.3, 0.7), nrow(data[["C1"]]))
  term <- normL2(data, bench$gfn * bench$xfn * x$trust_p_err, x$trust_err)
  obj <- term + constraintL2(mu = c(A_log = 0, k_log = log(0.4)), sigma = 0.5)
  init <- c(A_log = log(2), k_log = log(0.9), s0_log = log(0.05), srel = 0.3)
  rowsOf <- list()
  for (ev in c("C++", "R")) {
    withr::local_options(dMod2.evaluation = ev)
    o <- obj(init, hessian = "structured")
    d <- term(init)
    st <- o$structure
    # a row per data point and one per sigma that depends on parameters
    expect_identical(dim(st$jacobian), c(2L * nrow(data[["C1"]]), 4L), info = ev)
    expect_equal(drop(crossprod(st$jacobian, st$seeds)), d$gradient[names(init)], tolerance = 1e-10,
                 ignore_attr = TRUE, info = ev)
    expect_equal(st$hessian, d$hessian[names(init), names(init)], tolerance = 1e-10, info = ev)
    expect_null(obj(init)$structure)
    rowsOf[[ev]] <- st
  }
  expect_equal(rowsOf[["C++"]], rowsOf[["R"]], tolerance = 1e-8)
  # a scaled objective scales seeds and Hessian
  sc <- (2 %.*% term)(init, hessian = "structured")$structure
  expect_equal(sc$seeds, 2 * rowsOf[["C++"]]$seeds, tolerance = 1e-12)

  fit <- suppressWarnings(trust(obj, init, rinit = 0.1, iterlim = 1L, hessianMethod = "ssr1",
                                blather = TRUE))
  expect_true(fit$accept)
  x0 <- fit$argpath[1, ]; x1 <- fit$argtry[1, ]
  o0 <- obj(x0, hessian = "structured"); o1 <- obj(x1, hessian = "structured")
  w <- drop(crossprod(o1$structure$jacobian - o0$structure$jacobian, o1$structure$seeds))
  expect_equal(unname(fit$hessian), unname(o1$hessian + tcrossprod(w) / sum(w * (x1 - x0))),
               tolerance = 1e-8)
})
