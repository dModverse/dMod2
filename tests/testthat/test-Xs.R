# ============================================================================
# Behavioral tests for the prediction functions Xs / Xd / Xf.
#
# Closed-form analytical references throughout (linear decay, two-step
# cascade, Xd grid recovery, forced linear ODE). Hessian and second-order
# chain-rule semantics live in test-deriv2.R.
# ============================================================================

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("cppDE")
  testthat::skip_on_cran()
}

# The models of this file beyond the shared decay fixture, generated with
# compile = FALSE on first use and linked into one shared object.
xs_models <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    dir <- file.path(tempdir(), "xs_models")
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    withr::local_dir(dir)
    nm <- function(x) paste0(x, "_", as.integer(Sys.time()))

    cascade <- eqnlist() |>
      addReaction("A", "B", "k1 * A") |>
      addReaction("B", "C", "k2 * B")
    m_cascade <- odemodel(cascade, modelname = nm("xs_cascade"), compile = FALSE)
    p_cascade <- P(eqnvec(A = "A", B = "0", C = "0", k1 = "k1", k2 = "k2"),
                   condition = "C1", modelname = nm("xs_cascade_p"),
                   compile = FALSE)

    p_C2 <- P(eqnvec(A = "A_C2", k = "k_C2"), condition = "C2",
              modelname = nm("xs_p"), compile = FALSE)

    m_event <- odemodel(eqnlist() |> addReaction("A", "", "k * A", "decay"),
                        events = eventlist(var = "A", time = 5, value = "A_add",
                                           method = "add"),
                        modelname = nm("xs_event"), compile = FALSE)
    p_event <- P(eqnvec(A = "A", k = "k", A_add = "A_add"), condition = "C1",
                 modelname = nm("xs_event_p"), compile = FALSE)

    forced <- eqnlist() |>
      addReaction("",  "A", "F",     "production by forcing") |>
      addReaction("A", "",  "k * A", "decay")
    m_forc <- odemodel(forced, forcings = "F", modelname = nm("xs_forc"),
                       compile = FALSE)
    p_forc <- P(eqnvec(A = "A", k = "k"), condition = "C1",
                modelname = nm("xs_forc_p"), compile = FALSE)

    # Two separate builds of one model, compared against each other.
    f <- c(A = "-k1*A + k2*B",
           B =  "k1*A - k2*B")
    m_rep1 <- odemodel(f, modelname = nm("xs_rep_v1"), backend = "cppDE",
                       compile = FALSE)
    m_rep2 <- odemodel(f, modelname = nm("xs_rep_v2"), backend = "cppDE",
                       compile = FALSE)
    trafo <- c(A = "A", B = "B", k1 = "exp(log_k1)", k2 = "exp(log_k2)")
    p_rep <- P(trafo, modelname = nm("xs_rep_trafo"), compile = FALSE)
    p_rep_cl <- P(trafo, condition = "closed",
                  modelname = nm("xs_rep_trafo_cl"), compile = FALSE)
    p_rep_op <- P(c(A = "A", B = "B", k1 = "exp(log_k_open)", k2 = "exp(log_k2)"),
                  condition = "open", modelname = nm("xs_rep_trafo_op"),
                  compile = FALSE)

    compile(m_cascade, p_cascade, p_C2, m_event, p_event, m_forc, p_forc,
            m_rep1, m_rep2, p_rep, p_rep_cl, p_rep_op,
            output = nm("xs_models"), cores = 4L)

    cache <<- list(m_cascade = m_cascade, p_cascade = p_cascade, p_C2 = p_C2,
                   m_event = m_event, p_event = p_event, m_forc = m_forc,
                   p_forc = p_forc, m_rep1 = m_rep1, m_rep2 = m_rep2,
                   p_rep = p_rep, p_rep_cl = p_rep_cl, p_rep_op = p_rep_op)
    cache
  }
})


# ---- Xs: state trajectories (closed-form) -------------------------------

test_that("Xs on linear decay matches A0 * exp(-k * t)", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  prd_states <- bench$xfn * bench$pfn_id

  times <- c(0, 0.5, 1, 2, 5, 10)
  pars  <- c(A = 1.7, k = 0.42)
  out <- prd_states(times = times, pars = pars, deriv = FALSE)

  closed <- pars[["A"]] * exp(-pars[["k"]] * times)
  expect_equal(out$C1[, "A"], closed, tolerance = 1e-5)
})


test_that("Xs on two-step cascade matches the closed-form A(t), B(t), C(t)", {
  skip_if_no_compile()
  testthat::skip_if_not_installed("cppDE")
  mods <- xs_models()
  prd <- Xs(mods$m_cascade) * mods$p_cascade

  times <- c(0, 0.5, 1, 2, 4)
  pars <- c(A = 1.0, k1 = 0.7, k2 = 0.3)
  out <- prd(times = times, pars = pars, deriv = FALSE)

  closed <- truth_two_step(times, pars["A"], pars["k1"], pars["k2"])
  # Numerical integration against the closed form: the solver's local error
  # shows up around 1e-5 relative. testthat 3e compares element-wise, so the
  # bound must cover the worst time point rather than the vector average.
  expect_equal(out$C1[, "A"], closed$A, tolerance = 1e-4)
  expect_equal(out$C1[, "B"], closed$B, tolerance = 1e-4)
  expect_equal(out$C1[, "C"], closed$C, tolerance = 1e-4)
})


# ---- Xs: first-order sensitivities --------------------------------------

test_that("Xs sensitivities on linear decay match analytical d/d(A0,k)", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  prd_states <- bench$xfn * bench$pfn_id

  times <- c(0, 1, 2, 5)
  pars  <- c(A = 1.0, k = 0.5)
  out <- prd_states(times = times, pars = pars, deriv = TRUE)
  d <- attr(out$C1, "deriv")  # shape [time, var, par]

  expect_equal(d[, "A", "A"], exp(-pars[["k"]] * times),     tolerance = 1e-5)
  expect_equal(d[, "A", "k"], -times * pars[["A"]] * exp(-pars[["k"]] * times),
               tolerance = 1e-5)
})


test_that("Xs predictions across conditions are independent and parameter-local", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  prd_multi <- bench$xfn * (bench$pfn_id + xs_models()$p_C2)
  times <- c(0, 1, 2, 5)
  pars <- c(A = 1.0, k = 0.5, A_C2 = 2.0, k_C2 = 1.0)
  out <- prd_multi(times = times, pars = pars, deriv = TRUE)

  expect_equal(out$C1[, "A"], 1.0 * exp(-0.5 * times), tolerance = 1e-5)
  expect_equal(out$C2[, "A"], 2.0 * exp(-1.0 * times), tolerance = 1e-5)

  dC1 <- attr(out$C1, "deriv")
  dC2 <- attr(out$C2, "deriv")
  expect_setequal(dimnames(dC1)[[3]], c("A", "k"))
  expect_setequal(dimnames(dC2)[[3]], c("A_C2", "k_C2"))
})


test_that("Xs sensitivities match the analytical decay formulas at non-integer times", {
  skip_if_no_compile()
  bench <- fx_decay_compiled()
  prd_states <- bench$xfn * bench$pfn_id
  times <- c(0, 0.7, 1.4, 3.5)
  pars  <- c(A = 1.3, k = 0.65)

  out <- prd_states(times = times, pars = pars, deriv = TRUE)
  d <- attr(out$C1, "deriv")

  expect_equal(d[, "A", "A"], exp(-pars[["k"]] * times), tolerance = 1e-5)
  expect_equal(d[, "A", "k"],
               -times * pars[["A"]] * exp(-pars[["k"]] * times),
               tolerance = 1e-5)
})


# ---- Xs: events ---------------------------------------------------------

test_that("Xs with an 'add' event reproduces the analytical post-event trajectory", {
  # Linear decay with an additive jump at t0:
  #   pre   A(t)        = A0 * exp(-k * t)
  #   at t0 A(t0)       = A0 * exp(-k * t0) + Delta
  #   post  A(t)        = (A0 * exp(-k * t0) + Delta) * exp(-k * (t - t0))
  skip_if_no_compile()
  mods <- xs_models()
  prd <- Xs(mods$m_event) * mods$p_event

  A0 <- 1.0; k <- 0.4; Delta <- 0.5; t0 <- 5
  pars <- c(A = A0, k = k, A_add = Delta)
  times <- c(0, 1, 3, 4.99, 5.01, 6, 8, 10)
  out <- prd(times = times, pars = pars, deriv = FALSE)$C1

  A_t0_pre <- A0 * exp(-k * t0)
  A_t0_post <- A_t0_pre + Delta
  expected <- ifelse(times < t0,
                     A0 * exp(-k * times),
                     A_t0_post * exp(-k * (times - t0)))

  expect_equal(out[match(times, out[, "time"]), "A"], expected,
               tolerance = 1e-4)
})


# ---- Xs: forcings -------------------------------------------------------

test_that("Xs with constant forcing input matches the closed-form linear ODE solution", {
  # dA/dt = F - k*A with constant F gives A(t) = (A0 - F/k)*exp(-k*t) + F/k.
  skip_if_no_compile()
  mods <- xs_models()

  u_const <- 0.6
  forc <- data.frame(name = "F",
                     time = seq(0, 20, by = 0.5),
                     value = u_const)
  xf <- Xs(mods$m_forc, forcings = forc, condition = "C1")
  prd <- xf * mods$p_forc

  A0 <- 0.1; k <- 0.3
  pars <- c(A = A0, k = k)
  times <- c(0, 1, 2, 5, 10, 20)
  out <- prd(times = times, pars = pars, deriv = FALSE)$C1

  expected <- (A0 - u_const / k) * exp(-k * times) + u_const / k
  expect_equal(out[match(times, out[, "time"]), "A"], expected,
               tolerance = 1e-4)

  expect_lt(abs(out[match(20, out[, "time"]), "A"] - u_const / k),
            abs(A0 - u_const / k) * exp(-k * 20) * 2)
})


# ---- Xd: linear-interpolation grid prediction ---------------------------

test_that("Xd returns the grid values at the grid times", {
  grid <- data.frame(
    name = "A",
    time = c(0, 1, 2, 3),
    row.names = c("p0", "p1", "p2", "p3"))
  pars <- c(p0 = 1.0, p1 = 0.5, p2 = 0.25, p3 = 0.125)

  xfn <- Xd(grid, condition = "C1")
  out <- xfn(times = c(0, 1, 2, 3), pars = pars)$C1
  expect_equal(out[match(c(0, 1, 2, 3), out[, "time"]), "A"],
               c(1.0, 0.5, 0.25, 0.125), tolerance = 1e-12)
})


test_that("Xd linearly interpolates between grid points", {
  grid <- data.frame(
    name = "A",
    time = c(0, 1, 2, 3),
    row.names = c("p0", "p1", "p2", "p3"))
  pars <- c(p0 = 1.0, p1 = 0.5, p2 = 0.25, p3 = 0.125)

  xfn <- Xd(grid, condition = "C1")
  mid_times <- c(0.5, 1.5, 2.5)
  out <- xfn(times = mid_times, pars = pars)$C1
  expected <- c((1.0 + 0.5) / 2, (0.5 + 0.25) / 2, (0.25 + 0.125) / 2)
  expect_equal(out[match(mid_times, out[, "time"]), "A"],
               expected, tolerance = 1e-12)
})


# ---- Xf: no-sensitivity ODE prediction ----------------------------------

test_that("Xf reproduces the linear-decay closed form and emits no deriv attribute", {
  skip_if_no_compile()
  xfn <- Xf(fx_decay_compiled()$m, condition = "C1")

  times <- c(0, 1, 2, 5)
  pars <- c(A = 1.3, k = 0.42)
  out <- xfn(times = times, pars = pars)$C1

  expect_equal(out[match(times, out[, "time"]), "A"],
               pars[["A"]] * exp(-pars[["k"]] * times),
               tolerance = 1e-5)
  expect_null(attr(out, "deriv"))
})


# ============================================================================
# Xs.cppDE theta-sensitivity path: heap vs stack AD slab parity
# (Phi'(theta) as tangent; per-condition varying theta counts via two-condition setup)
# ============================================================================

test_that("Heap and stack AD slabs match on a single-condition linear model", {

  # Default heap slab vs explicit stack slab (B,log_k1,log_k2}), same
  # parameter transformation for both.
  mods <- xs_models()
  tight <- list(atol = 1e-10, rtol = 1e-10)
  x1 <- Xs(mods$m_rep1, optionsSens = tight) * mods$p_rep
  x2 <- Xs(mods$m_rep2, optionsSens = tight) * mods$p_rep

  theta <- c(A = 1.0, B = 0.2, log_k1 = log(0.5), log_k2 = log(0.3))
  times <- seq(0, 3, length.out = 7)

  pred1 <- x1(times, theta,
              conditions = NULL,
              deriv = TRUE)
  pred2 <- x2(times, theta,
              conditions = NULL,
              deriv = TRUE)

  d1 <- getDerivs(pred1)[[1]]
  d2 <- getDerivs(pred2)[[1]]

  arr1 <- attr(d1, "deriv")
  arr2 <- attr(d2, "deriv")
  expect_equal(dim(arr1), dim(arr2))
  expect_equal(dimnames(arr1)[[2]], dimnames(arr2)[[2]])
  expect_equal(dimnames(arr1)[[3]], dimnames(arr2)[[3]])
  expect_equal(as.numeric(arr1), as.numeric(arr2), tolerance = 1e-6)
})


test_that("Heap/stack parity holds with per-condition varying theta subsets", {

  # Stack upper bound: any condition may activate up to 4 thetas.
  # Condition "closed" uses log_k1; condition "open" uses log_k_open instead.
  # Global theta set has 5 elements; each condition activates 4.
  mods <- xs_models()
  p <- mods$p_rep_cl + mods$p_rep_op

  tight <- list(atol = 1e-10, rtol = 1e-10)
  x1 <- Xs(mods$m_rep1, optionsSens = tight) * p
  x2 <- Xs(mods$m_rep2, optionsSens = tight) * p

  theta <- c(A = 1.0, B = 0.2,
             log_k1 = log(0.5), log_k_open = log(0.8), log_k2 = log(0.3))
  times <- seq(0, 2, length.out = 5)

  pred1 <- x1(times, theta, deriv = TRUE)
  pred2 <- x2(times, theta, deriv = TRUE)

  for (cond in c("closed", "open")) {
    arr1 <- attr(pred1[[cond]], "deriv")
    arr2 <- attr(pred2[[cond]], "deriv")
    expect_equal(dim(arr1), dim(arr2),
                 info = paste("shape mismatch for condition", cond))
    expect_equal(dim(arr2)[3], 4L,
                 info = paste("ncol mismatch for condition", cond))
    expect_equal(dimnames(arr1)[[3]], dimnames(arr2)[[3]],
                 info = paste("theta names mismatch for condition", cond))
    expect_equal(as.numeric(arr1), as.numeric(arr2), tolerance = 1e-6,
                 info = paste("values mismatch for condition", cond))
  }
})
