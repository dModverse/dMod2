# Behavioral tests for trust() (trust-region optimizer, C++ implementation).
#
# Verifies:
#   * convergence to the exact minimum of a quadratic in one trust step
#   * recovery of simulated-truth parameters from a noisy decay dataset
#   * extra objfun args passed via closure binding (no more ... forwarding)
#   * parscale rescaling lands at the same minimum
#   * parupper clamps a single component
#   * blather returns the per-iter trace
#
# Trust is also exercised end-to-end via normL2 -> trust in
# test-mstrust-profile.R and the existing FOCEI tests.

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("cppDE")
  testthat::skip_on_cran()
}


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
  # Practical envelope at this noise level: ~10% relative on (A, k). The
  # exponential decay has a known A * k -> kA scale identifiability that
  # makes k a few % less constrained than A.
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


## ---- parscale invariance --------------------------------------------

test_that("trust with parscale lands at the same minimum as the unscaled run", {
  target <- c(a = 1.0, b = -0.5, c = 2.3)
  obj <- .quadratic_objfn(target)
  init <- c(a = 0, b = 0, c = 0)

  fit_plain <- trust(obj, init, rinit = 1, rmax = 100, iterlim = 50,
                     printIter = FALSE)
  fit_scaled <- trust(obj, init, parscale = c(a = 2, b = 0.5, c = 10),
                      rinit = 1, rmax = 100, iterlim = 50,
                      printIter = FALSE)
  expect_equal(unname(fit_scaled$argument[names(target)]),
               unname(fit_plain$argument[names(target)]),
               tolerance = 1e-6)
})


## ---- bounds clamp the optimum ---------------------------------------

test_that("trust honors parupper on one component while leaving the other free", {
  # 2D problem; unconstrained min at (a, b) = (5, 1). Upper bound on `a`
  # at 2 clamps that component; `b` converges freely to 1.
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
  # Construct an objective whose Hessian has eigenvalues spanning ~20 orders
  # of magnitude (4e-3 down to 1e-21). LAPACK's smallest eigenvalue then
  # comes back as tiny negative numerical noise, which previously sent the
  # subproblem into the hard-hard branch and produced steps thousands of
  # times larger than the trust radius. The proper Moré-Sorensen easy
  # branch keeps ||p|| ~ r.
  K <- 7
  # Eigenvalues spanning many orders of magnitude.
  vals <- c(4.4e-3, 6.0e-4, 1.9e-5, 9.6e-7, 8.0e-10, 8.2e-19, 9.4e-21)
  set.seed(42)
  Q <- qr.Q(qr(matrix(rnorm(K * K), K, K)))    # random orthogonal basis
  H <- Q %*% diag(vals) %*% t(Q)
  H <- 0.5 * (H + t(H))
  g <- as.numeric(Q %*% c(0.5, 0.6, 0.4, 0.3, 0.2, 0.0, 0.0))
  par_target <- setNames(rep(0, K), paste0("p", seq_len(K)))
  init       <- setNames(rep(0.1, K), paste0("p", seq_len(K)))

  # A static objective with constant H and gradient g + H * (p - init).
  # Minimum is wherever H(p) = -g; with H near-singular many minima are
  # acceptable. We only care that trust does not blow up.
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
  # Every accepted/proposed step must stay within rmax (and the per-iter
  # radius). Without the fix this assertion fails by ~5 orders of magnitude.
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
#
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
  # min 0.5 (p - 3)' A (p - 3) s.t. x <= 1, with A = [[2,1],[1,2]].
  # Active bound at x = 1 leaves d2 = -d1/2, hence y = 4 and g = (-3, 0).
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
  # A diagonal Hessian never needs stepback: the metric's curvature term C
  # already keeps the scaled Newton step inside the box. It is off-diagonal
  # coupling that pushes a coordinate out, so this is a correlated problem with
  # the optimum far outside an asymmetric box.
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
  # Under the old fterm-only test this reported converged = TRUE on the very
  # first flat rejected step.
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


test_that("a flat objective under repeated rejection stops as stagnation", {
  # Without an rmin floor the run ends once the radius has been cut far enough
  # that the objective no longer moves. That is reported as "stagnation", never
  # as "gradient" -- the distinction a caller needs, because the optimiser
  # cannot tell a genuine noise floor from a bad model.
  misleading <- function(p, ...) {
    list(value = sum(p^2), gradient = -2 * p, hessian = diag(length(p)))
  }
  fit <- trust(misleading, c(a = 1, b = 1), rinit = 1, rmax = 10, iterlim = 500)
  expect_equal(fit$stopReason, "stagnation")
  expect_lt(fit$iterations, 500L)
})


test_that("a rank-deficient Gauss-Newton Hessian gives a minimum-norm step", {
  # f = 0.5 (x + y - 2)^2. H = J'J is rank 1 with null space (1, -1), and
  # g = J'r is orthogonal to it -- the exactly-degenerate hard case. Extending
  # along the zero-curvature direction would drift along x - y for no gain.
  obj <- function(p, ...) {
    res <- as.numeric(p[1] + p[2] - 2)
    J   <- matrix(c(1, 1), nrow = 1)
    list(value = 0.5 * res^2, gradient = as.numeric(t(J) * res), hessian = t(J) %*% J)
  }
  fit <- trust(obj, c(x = 0, y = 0), rinit = 1, rmax = 10, iterlim = 100)

  expect_true(all(is.finite(fit$argument)))
  expect_equal(sum(fit$argument), 2, tolerance = 1e-6)
  # Minimum-norm solution of x + y = 2 from (0, 0) is (1, 1): no null-space drift.
  expect_equal(unname(fit$argument[["x"]] - fit$argument[["y"]]), 0, tolerance = 1e-6)
})


test_that("curvature at roundoff level does not steer the degenerate step", {
  # H = J'J - d*I has eigenvalues {-d, 2-d}, so lam_min is negative by construction
  # rather than by LAPACK rounding -- the test above only catches the drift on builds
  # whose dsyevr happens to return the singular direction as negative.
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
  # f = 0.5 (J p - 2)^2 with J = (1, 2, 3). H = J'J is rank 1, so its null space
  # is two-dimensional and comes back as a pair of eigenvalues split by roundoff.
  # No LAPACK build returns a basis for that space which g is exactly orthogonal
  # to, so the drift the two tests above only see on some builds shows up here on
  # every one of them.
  J <- matrix(c(1, 2, 3), nrow = 1)
  obj <- function(p, ...) {
    res <- as.numeric(J %*% p - 2)
    list(value = 0.5 * res^2, gradient = as.numeric(t(J) * res), hessian = t(J) %*% J)
  }
  fit <- trust(obj, c(x = 0, y = 0, z = 0), rinit = 1, rmax = 10, iterlim = 100)

  expect_equal(as.numeric(J %*% fit$argument), 2, tolerance = 1e-6)
  # Minimum-norm solution of J p = 2 from the origin is J' * 2 / ||J||^2.
  expect_equal(unname(fit$argument), as.numeric(t(J)) * 2 / sum(J^2),
               tolerance = 1e-6)
})


test_that("bounds compose with parscale, parinit on a bound, and minimize = FALSE", {
  target <- c(x = 5, y = 5)

  # parinit sitting exactly on the bound must be nudged inside, not frozen.
  on_bound <- trust(.quadratic_objfn(target), c(x = 1, y = 0),
                    rinit = 1, rmax = 10, parupper = c(x = 1, y = Inf))
  expect_true(on_bound$converged)
  expect_equal(unname(on_bound$argument[["y"]]), 5, tolerance = 1e-6)
  expect_true(on_bound$atBound[["x"]])

  # parscale must not move the optimum.
  # Tight tolerances so both runs reach the optimum rather than stopping at
  # their own frame-dependent distance from the bound.
  scaled <- trust(.quadratic_objfn(target), c(x = 0, y = 0),
                  rinit = 1, rmax = 10, parupper = c(x = 1, y = Inf),
                  parscale = c(10, 0.1), iterlim = 200,
                  tolControl = list(gtol = 1e-12))
  plain  <- trust(.quadratic_objfn(target), c(x = 0, y = 0),
                  rinit = 1, rmax = 10, parupper = c(x = 1, y = Inf),
                  iterlim = 200, tolControl = list(gtol = 1e-12))
  expect_equal(unname(scaled$argument), unname(plain$argument), tolerance = 1e-5)

  # Maximising a concave objective under the same bound.
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
  # |v| == 1 and C == 0 there, so the reflective scheme reduces to the old one.
  target <- c(a = 1.0, b = -0.5, c = 2.3)
  init   <- c(a = 0, b = 0, c = 0)
  refl <- trust(.quadratic_objfn(target), init, rinit = 0.5, rmax = 10,
                stepControl = list(boundary = "reflective"))
  clip <- trust(.quadratic_objfn(target), init, rinit = 0.5, rmax = 10,
                stepControl = list(boundary = "clip"))
  expect_equal(unname(refl$argument), unname(clip$argument), tolerance = 1e-10)
})


test_that("boundary = 'clip' keeps landing exactly on the bound", {
  fit <- trust(.quadratic_objfn(c(a = 5, b = 1)), c(a = 0, b = 0),
               rinit = 1, rmax = 10, parupper = c(a = 2, b = Inf),
               stepControl = list(boundary = "clip"))
  expect_identical(unname(fit$argument[["a"]]), 2)
  expect_error(
    trust(.quadratic_objfn(c(a = 1)), c(a = 0), rinit = 1, rmax = 10,
          stepControl = list(boundary = "nonsense")),
    "should be one of")
})


test_that("an argument that moved into a control list is rejected by name", {
  obj <- .quadratic_objfn(c(a = 1))
  expect_error(trust(obj, c(a = 0), gtol = 1e-8), "tolControl\\$gtol")
  expect_error(trust(obj, c(a = 0), qnMemory = 3L), "qnControl\\$qnMemory")
  expect_error(trust(obj, c(a = 0), boundary = "clip"), "stepControl\\$boundary")
})


test_that("a flat, high-value start makes progress instead of stopping at once", {
  # Weak sensitivities with large residuals: |f| is huge while |g| is small, as
  # when an ODE model's parameters run into a saturated regime. Any gradient
  # criterion scaled by |f| declares this converged before the first step --
  # |f| grows quadratically in the residuals, |g| only linearly -- which is why
  # gtol has no relative counterpart.
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
               fb   = list(hessianFallback = "bfgs"))
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

  for (hm in c("gn", "bfgs", "sr1")) {
    fit <- trust(obj, init, hessianMethod = hm)
    expect_identical(names(fit$evalBySource), c("gn", "bfgs", "sr1"), info = hm)
    expect_identical(sum(fit$evalBySource), fit$neval, info = hm)
    expect_identical(fit$nSwitch, 0L, info = hm)                 # no handover
    expect_gt(fit$evalBySource[[hm]], 0L)                        # charged to itself
  }
  expect_identical(trust(obj, init)$evalBySource[["bfgs"]], 0L)  # gn: only gn

  # A fallback hands over once and both phases are charged separately.
  hy <- suppressWarnings(trust(.flat_objfn(2L), init, hessianFallback = "bfgs"))
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

test_that("a bfgs fallback takes over from gn on stagnation", {
  obj  <- .flat_objfn(2L)
  init <- c(a = 0, b = 0)

  gn <- suppressWarnings(trust(obj, init, hessianMethod = "gn", blather = TRUE))
  expect_true(all(gn$hessianSource == "gn"))

  hy <- suppressWarnings(trust(obj, init, hessianFallback = "bfgs", blather = TRUE))
  expect_identical(hy$hessianSource[1], "gn")   # always starts on gn
  expect_true("bfgs" %in% hy$hessianSource)      # and switches
  expect_gt(hy$qnEval, 0L)
})

test_that("the fallback takes over on a value stop too, and gets past where gn stops", {
  target <- c(a = 1.0, b = -0.5, c = 2.3)
  obj    <- .misscaled_objfn(target)
  init   <- c(a = 0, b = 0, c = 0)

  gn <- trust(obj, init, rinit = 1, rmax = 10, iterlim = 500,
              hessianMethod = "gn")
  expect_true(gn$stopReason %in% c("fvalue", "preddiff"))
  expect_true(gn$converged)
  expect_gt(max(abs(gn$argument - target)), 1e-3)   # nowhere near it

  hy <- trust(obj, init, rinit = 1, rmax = 10, iterlim = 500,
              hessianFallback = "bfgs", blather = TRUE)
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
  hy <- suppressWarnings(trust(obj, init, hessianFallback = "bfgs",
                               iterlim = 500))
  expect_true(hy$converged)
  expect_identical(hy$stopReason, "stagnation")
  expect_lt(hy$iterations, 500L)
})

test_that("hessianInit is inert for the methods that need the objective Hessian", {
  target <- c(a = 1.0, b = -0.5, c = 2.3)
  obj <- .quadratic_objfn(target)
  init <- c(a = 0, b = 0, c = 0)
  runs <- list(gn = list(), fb = list(hessianFallback = "bfgs"))
  for (hm in names(runs)) {
    args <- c(list(obj, init, rinit = 1), runs[[hm]])
    ref <- do.call(trust, args)
    alt <- do.call(trust, c(args, list(qnControl = list(hessianInit = "identity"))))
    expect_equal(alt$argument, ref$argument, ignore_attr = TRUE, info = hm)
    expect_identical(alt$neval, ref$neval, info = hm)
  }
})

test_that("quasi-Newton methods require the reflective boundary and a known name", {
  obj <- .quadratic_objfn(c(a = 1.0))
  init <- c(a = 0)
  expect_error(trust(obj, init, stepControl = list(boundary = "clip"),
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


## ---- Selectable Hessian fallback --------------------------------------

test_that("a fallback run is the same whichever way it is spelled", {
  init <- c(a = 0, b = 0)
  for (nm in c("flat", "misscaled")) {
    obj <- if (nm == "flat") .flat_objfn(2L) else .misscaled_objfn(c(a = 1, b = -2))
    old <- trust(obj, init, iterlim = 200, hessianFallback = "bfgs",
                 blather = TRUE)
    new <- trust(obj, init, iterlim = 200, hessianFallback = "bfgs",
                 fallbackLimit = 1L, blather = TRUE)
    expect_identical(new$neval, old$neval, info = nm)
    expect_identical(new$value, old$value, info = nm)
    expect_identical(new$nSwitch, old$nSwitch, info = nm)
    expect_identical(new$hessianSource, old$hessianSource, info = nm)
    expect_identical(new$argument, old$argument, info = nm)
  }
})

test_that("no fallback reproduces the plain run", {
  obj  <- .misscaled_objfn(c(a = 1, b = -2))
  init <- c(a = 0, b = 0)
  plain <- trust(obj, init, iterlim = 200)
  for (off in list(list(hessianFallback = "none"),
                   list(hessianFallback = "bfgs", fallbackLimit = 0L))) {
    f <- do.call(trust, c(list(obj, init, iterlim = 200), off))
    expect_identical(f$neval, plain$neval)
    expect_identical(f$value, plain$value)
    expect_identical(f$nSwitch, 0L)
  }
})

test_that("the fallback source is selectable and is charged for its evaluations", {
  obj  <- .misscaled_objfn(c(a = 1, b = -2))
  init <- c(a = 0, b = 0)
  for (fb in c("bfgs", "sr1")) {
    f <- trust(obj, init, iterlim = 200, hessianFallback = fb, blather = TRUE)
    expect_identical(f$nSwitch, 1L, info = fb)
    expect_true(fb %in% f$hessianSource, info = fb)
    expect_gt(f$evalBySource[[fb]], 0)
    expect_identical(sum(f$evalBySource), f$neval, info = fb)
  }
})

test_that("fallbackLimit above one hands control back to the primary source", {
  # A flat objective stalls under either source, so every allowed handover
  # fires and the sources alternate.
  obj  <- .flat_objfn(2L)
  init <- c(a = 0, b = 0)
  one <- trust(obj, init, iterlim = 300, hessianFallback = "bfgs",
               fallbackLimit = 1L, blather = TRUE)
  two <- trust(obj, init, iterlim = 300, hessianFallback = "bfgs",
               fallbackLimit = 2L, blather = TRUE)

  expect_identical(one$nSwitch, 1L)
  expect_identical(two$nSwitch, 2L)
  # Back on gn after the second handover, and the reseed evaluation is charged.
  expect_identical(utils::tail(two$hessianSource, 1L), "gn")
  expect_identical(sum(two$evalBySource), two$neval)
})

test_that("an unknown Hessian source is rejected by name", {
  obj <- .quadratic_objfn(c(a = 1.0))
  expect_error(trust(obj, c(a = 0), hessianMethod = "nope"), "should be one of")
  expect_error(trust(obj, c(a = 0), hessianFallback = "nope"), "should be one of")
})

# An objective may decline to build a Hessian, and NULL is how it says so.
# `give` decides which evaluations carry one.
.declining_objfn <- function(give = function(n) FALSE) {
  n <- 0L
  function(p, hessian = TRUE, ...) {
    n <<- n + 1L
    x <- p[[1]]; y <- p[[2]]
    gr <- c(-2 * (1 - x) - 400 * x * (y - x^2), 200 * (y - x^2))
    hs <- matrix(c(2 - 400 * y + 1200 * x^2, -400 * x, -400 * x, 200), 2, 2)
    names(gr) <- names(p); dimnames(hs) <- list(names(p), names(p))
    list(value = (1 - x)^2 + 100 * (y - x^2)^2, gradient = gr,
         hessian = if (isTRUE(hessian) && give(n)) hs else NULL)
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

test_that("gn says which argument it cannot honour when no Hessian arrives", {
  init <- c(x = -1.2, y = 1)
  expect_error(trust(.declining_objfn(), init, rinit = 1, rmax = 10),
               "no Hessian at parinit")
  expect_error(trust(.declining_objfn(), init, rinit = 1, rmax = 10,
                     stepControl = list(boundary = "clip")),
               "no Hessian at parinit")
})

test_that("a trial point without a Hessian is a failed evaluation, not an error", {
  # Seeded, then declining: the kernel must treat the trial as infeasible
  # rather than let the conversion throw from outside its own handler.
  init <- c(x = -1.2, y = 1)
  expect_warning(
    fit <- trust(.declining_objfn(function(n) n == 1L), init, rinit = 1,
                 rmax = 10, iterlim = 20),
    "evaluation failed")
  expect_identical(fit$stopReason, "objfun")
  expect_false(fit$converged)
})
