# Tests for Pimpl() (implicit parameter transformation, steady states by
# pseudo-transient continuation), structural zero-state detection and the
# conserved-quantity basis of eqnlist objects.

skip_if_no_compile <- function() {
  testthat::skip_if_not_installed("cppDE")
  testthat::skip_on_cran()
}


## ---- Shared models -----------------------------------------------------

# Every model a test evaluates, generated once and compiled on first use. The
# ODE sources need -fopenmp and the algebraic ones do not, so each flag set
# gets its own shared object.
.pp_env <- new.env(parent = emptyenv())

pp_models <- function() {
  if (!is.null(.pp_env$models)) return(.pp_env$models)
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)
  stamp <- as.integer(Sys.time())
  mn <- function(x) paste0("test_pp_", x, "_", stamp)

  el_2moiety <- eqnlist() |>
    addReaction("A", "B", "k1*A") |> addReaction("B", "A", "k2*B") |>
    addReaction("C", "D", "k3*C") |> addReaction("D", "C", "k4*D")

  el_dimer <- eqnlist() |>
    addReaction("2*M", "D", "ka*M^2") |>
    addReaction("D", "2*M", "kd*D") |>
    customTotals(list(total_MD = "M + 2*D"))

  el_zero <- eqnlist() |>
    addReaction("",       "R",  "k_pr_R") |>
    addReaction("R",      "",   "k_dg_R * R") |>
    addReaction("L + R",  "LR", "k_on * L * R") |>
    addReaction("LR",     "L + R", "k_off * LR") |>
    addReaction("LR",     "",   "k_dg_LR * LR")
  zero <- Pimpl(el_zero, modelname = mn("zero"), verbose = FALSE)

  el_ERK <- eqnlist() |>
    addReaction("ERK",  "pERK", "k1 * ERK") |>
    addReaction("pERK", "ERK",  "k2 * pERK")

  # production-decay, two independent objects
  el_pd <- eqnlist() |>
    addReaction("", "A", "k_in") |> addReaction("A", "", "k_out * A")
  pd  <- Pimpl(el_pd, modelname = mn("pd"), verbose = FALSE)
  pd2 <- Pimpl(el_pd, modelname = mn("pd2"), verbose = FALSE)

  # X exists only with the plasmid; Y is produced constitutively and from X
  el_pl <- eqnlist() |>
    addReaction("", "X", "k_tx * plasmid") |> addReaction("X", "", "k_dg * X") |>
    addReaction("", "Y", "k_y + k_xy * X") |> addReaction("Y", "", "k_dy * Y")
  plasmid <- Pimpl(el_pl, modelname = mn("plasmid"), verbose = FALSE)

  # bistable: low and high stable steady state, unstable one between
  el_bi <- eqnlist() |>
    addReaction("", "X", "k0 + k1 * X^2 / (K^2 + X^2)") |>
    addReaction("X", "", "d * X")
  bistable <- Pimpl(el_bi, modelname = mn("bistable"), verbose = FALSE)

  # a linear residual read as a flow: the only root is unstable
  unstable <- Pimpl(c(x = "x - a"), parameters = "a", flow = TRUE,
                    modelname = mn("unstable"), verbose = FALSE)

  ## Pimpl and Pexpl
  lin <- Pimpl(c(x = "x - a"), parameters = "a", modelname = mn("lin"), verbose = FALSE)

  d2 <- Pimpl(c(C = "k1*(totA - C)*(totB - C) - km*C"),
              parameters = c("k1","km","totA","totB"), deriv2 = TRUE,
              modelname = mn("d2"), verbose = FALSE, controlsPTC = list(reltol = 1e-13))

  d2_2 <- Pimpl(c(x1 = "a*x1 - b*x2 - 1", x2 = "x1*x2 - c"),
                parameters = c("a","b","c"), deriv2 = TRUE,
                modelname = mn("d2_2"), verbose = FALSE, controlsPTC = list(reltol = 1e-13))

  el_AB <- eqnlist() |>
    addReaction("A", "B", "k*A") |>
    addReaction("B", "A", "km*B")
  d2_cq <- Pimpl(el_AB, parameters = c("k","km"), deriv2 = TRUE,
                 modelname = mn("d2_cq"), verbose = FALSE, controlsPTC = list(reltol = 1e-13))

  moiety2_im <- Pimpl(el_2moiety, modelname = mn("2moiety_im"))

  dimer_im <- Pimpl(el_dimer, modelname = mn("dimer_im"))

  el_recycle <- eqnlist() |>
    addReaction("G + S", "GS", "k1*G*S") |>
    addReaction("GS", "G + S", "k1r*GS") |>
    addReaction("GS", "G + P", "k1c*GS") |>
    addReaction("P", "S", "k3*P")
  recycle <- Pimpl(el_recycle, modelname = mn("recycle"))

  sing <- Pimpl(c(x1 = "x1 + x2 - s", x2 = "x1 + x2 - s"), parameters = "s",
                deriv2 = FALSE, modelname = mn("sing"), verbose = FALSE,
                controlsPTC = list(positive = FALSE))

  noroot <- Pimpl(c(x = "x*x + 1"), parameters = character(0),
                  modelname = mn("noroot"), verbose = FALSE,
                  controlsPTC = list(positive = FALSE, maxit = 50L))

  noparam <- Pimpl(c(A = "A - 1.0"), parameters = NULL, modelname = mn("noparam"))

  erk_im <- Pimpl(el_ERK, parameters = c("k1", "k2"),
                  modelname = mn("erk_im"), verbose = FALSE)

  px <- P(list(C1 = c(k_in = "s",     k_out = "1", A = "1"),
               C2 = c(k_in = "2 * s", k_out = "1", A = "1")),
          modelname = mn("px"), verbose = FALSE)
  # two conditions with the same steady-state inputs
  px_same <- P(list(C1 = c(k_in = "s", k_out = "1", A = "1"),
                    C2 = c(k_in = "s", k_out = "1", A = "1")),
               modelname = mn("px_same"), verbose = FALSE)

  compile(lin, d2, d2_2, d2_cq, moiety2_im, dimer_im, recycle, sing, noroot,
          noparam, erk_im, px, px_same, zero, pd, pd2, plasmid, bistable, unstable,
          output = mn("impl"), cores = test_cores())

  .pp_env$models <- list(
    lin = lin, d2 = d2, d2_2 = d2_2, d2_cq = d2_cq, moiety2_im = moiety2_im,
    dimer_im = dimer_im, recycle = recycle, el_recycle = el_recycle, sing = sing,
    noroot = noroot, noparam = noparam, erk_im = erk_im, px = px, px_same = px_same,
    zero = zero, pd = pd, pd2 = pd2, plasmid = plasmid, bistable = bistable,
    unstable = unstable)
}


## ---- Pimpl: linear constraint (value) ---------------------------------

test_that("Pimpl on x - a = 0 returns x = a (root-finding correctness)", {
  skip_if_no_compile()
  pf <- pp_models()$lin

  pars <- c(a = 2.0, x = 0.5)
  out <- pf(pars, deriv = TRUE)[[1]]

  # Pimpl returns inner state plus pass-through inputs; pick by name.
  expect_equal(as.numeric(out["x"]), pars[["a"]], tolerance = 1e-8)

  # The IFT sensitivity has the sign of the closed form.
  J <- attr(out, "deriv")
  expect_equal(unname(J["x", "a"]), 1, tolerance = 1e-6)
})


## ---- Pimpl: IFT second-order on a 1D quadratic SS ---------------------

test_that("Pimpl(deriv2) matches FD on a 1D mass-action steady state", {
  skip_if_no_compile()
  # A + B <-> C with the totals substituted in: one quadratic residual.
  pf <- pp_models()$d2

  p0 <- c(k1 = 1, km = 0.5, totA = 2, totB = 3, C = 0.5)
  inputs <- c("k1","km","totA","totB")
  v <- pf(p0, deriv2 = TRUE)[[1]]
  J_an <- attr(v, "deriv")
  H_an <- attr(v, "deriv2")

  # The Jacobian matches the closed form.
  expect_equal(unname(J_an["C", inputs]),
               c(0.3, -0.6, 0.6, 0.2), tolerance = 1e-6)

  # FD-check the Hessian against the analytical Jacobian
  jac_at <- function(par_in) {
    pp <- p0; pp[names(par_in)] <- par_in
    attr(pf(pp, deriv2 = FALSE)[[1]], "deriv")["C", inputs, drop = FALSE]
  }
  H_fd <- matrix(0, length(inputs), length(inputs), dimnames = list(inputs, inputs))
  h <- 1e-5
  for (j in seq_along(inputs)) {
    pj <- p0[inputs]; pj[j] <- pj[j] + h
    pn <- p0[inputs]; pn[j] <- pn[j] - h
    H_fd[, j] <- (jac_at(pj) - jac_at(pn)) / (2*h)
  }
  H_fd <- 0.5 * (H_fd + t(H_fd))
  expect_equal(unname(H_an["C", inputs, inputs]), unname(H_fd),
               tolerance = 1e-5)
})


## ---- Pimpl: IFT second-order with coupled n_dep = 2 -------------------

test_that("Pimpl(deriv2) matches FD on a 2-state coupled SS", {
  skip_if_no_compile()
  pf <- pp_models()$d2_2

  p0 <- c(a = 2, b = 0.5, c = 1, x1 = 0.5, x2 = 0.5)
  inputs <- c("a","b","c")
  v <- pf(p0, deriv2 = TRUE)[[1]]
  J_an <- attr(v, "deriv")
  H_an <- attr(v, "deriv2")

  jac_at <- function(par_in) {
    pp <- p0; pp[names(par_in)] <- par_in
    attr(pf(pp, deriv2 = FALSE)[[1]], "deriv")[c("x1","x2"), inputs, drop = FALSE]
  }
  H_fd <- array(0, c(2, length(inputs), length(inputs)),
                dimnames = list(c("x1","x2"), inputs, inputs))
  h <- 1e-5
  for (j in seq_along(inputs)) {
    pj <- p0[inputs]; pj[j] <- pj[j] + h
    pn <- p0[inputs]; pn[j] <- pn[j] - h
    H_fd[, , j] <- (jac_at(pj) - jac_at(pn)) / (2*h)
  }
  H_fd <- 0.5 * (H_fd + aperm(H_fd, c(1,3,2)))
  expect_equal(unname(H_an[c("x1","x2"), inputs, inputs]),
               unname(H_fd), tolerance = 1e-4)
})


## ---- Pimpl: IFT second-order through CQ-elim reconstruction -----------

test_that("Pimpl(deriv2) propagates Hessian through CQ-eliminated species", {
  skip_if_no_compile()
  # A <-> B with one conserved quantity: one species is eliminated, one stays dependent.
  pf <- pp_models()$d2_cq

  p0 <- c(k = 2, km = 0.5, total_1 = 3, A = 0.5, B = 0.5)
  inputs <- c("k","km","total_1")
  v <- pf(p0, deriv2 = TRUE)[[1]]
  J_an <- attr(v, "deriv")
  H_an <- attr(v, "deriv2")

  # States match the closed form.
  Bs <- p0[["k"]] * p0[["total_1"]] / (p0[["k"]] + p0[["km"]])
  expect_equal(as.numeric(v["B"]), Bs, tolerance = 1e-6)
  expect_equal(as.numeric(v["A"]), p0[["total_1"]] - Bs, tolerance = 1e-6)

  # The conservation law ties the two Hessians together.
  expect_equal(unname(H_an["A", inputs, inputs]),
               -unname(H_an["B", inputs, inputs]), tolerance = 1e-10)

  # FD-check on B's analytical Hessian
  jac_at <- function(par_in) {
    pp <- p0; pp[names(par_in)] <- par_in
    attr(pf(pp, deriv2 = FALSE)[[1]], "deriv")["B", inputs, drop = FALSE]
  }
  H_fd <- matrix(0, length(inputs), length(inputs), dimnames = list(inputs, inputs))
  h <- 1e-5
  for (j in seq_along(inputs)) {
    pj <- p0[inputs]; pj[j] <- pj[j] + h
    pn <- p0[inputs]; pn[j] <- pn[j] - h
    H_fd[, j] <- (jac_at(pj) - jac_at(pn)) / (2*h)
  }
  H_fd <- 0.5 * (H_fd + t(H_fd))
  expect_equal(unname(H_an["B", inputs, inputs]), unname(H_fd),
               tolerance = 1e-5)
})


## ---- getEquations with conservation rows ------------------------------

test_that("Pimpl solves the full system together with the conservation rows", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  # Pimpl keeps every moiety species and every rate equation and adds the
  # conservation row; nothing is reconstructed by subtraction.
  el <- eqnlist()
  el <- addReaction(el, "A", "B", "k*A")
  el <- addReaction(el, "B", "A", "km*B")

  pf <- Pimpl(el, parameters = c("k","km"),
              modelname = paste0("test_Pimpl_eq_", as.integer(Sys.time())),
              compile = FALSE)
  eqs <- getEquations(pf)[[1]]
  expect_setequal(names(eqs), c("A", "B", "total_1"))
  expect_match(eqs[["total_1"]], "total_1", fixed = TRUE)
  expect_true(grepl("k", eqs[["A"]]))
  expect_true("total_1" %in% getParameters(pf))
})

## ---- CQ spectrum: independent / non-unit-coef / overlapping -----------

# Conserved quantities enter as rows C x = T: independent CQs give one row and
# total each, non-unit stoichiometry a coefficient in C, overlapping CQs shared species.

test_that("two independent conserved moieties solve to closed form", {
  skip_if_no_compile()
  # two disjoint moieties
  m <- pp_models()

  p0 <- c(k1 = 2, k2 = 1, k3 = 0.5, k4 = 1.5, total_1 = 4, total_2 = 6,
          A = 1, B = 1, C = 1, D = 1)
  Bs <- p0[["k1"]] * p0[["total_1"]] / (p0[["k1"]] + p0[["k2"]]); As <- p0[["total_1"]] - Bs
  Ds <- p0[["k3"]] * p0[["total_2"]] / (p0[["k3"]] + p0[["k4"]]); Cs <- p0[["total_2"]] - Ds

  o <- m$moiety2_im(p0)[[1]]
  expect_equal(as.numeric(o[c("A","B","C","D")]), c(As, Bs, Cs, Ds), tolerance = 1e-8)
})

test_that("non-unit stoichiometric coefficient is handled (2*M <-> D)", {
  skip_if_no_compile()
  # Dimerisation conserves the monomer count.
  m <- pp_models()

  p0 <- c(ka = 1, kd = 2, total_MD = 5)
  r  <- p0[["ka"]] / p0[["kd"]]
  Ms <- (-1 + sqrt(1 + 8 * r * p0[["total_MD"]])) / (4 * r)
  Ds <- r * Ms^2

  oi <- m$dimer_im(p0)[[1]]
  expect_equal(as.numeric(oi["M"]), Ms, tolerance = 1e-8)
  expect_equal(as.numeric(oi["D"]), Ds, tolerance = 1e-8)
  expect_equal(as.numeric(oi["M"] + 2 * oi["D"]), p0[["total_MD"]], tolerance = 1e-12)
})

test_that("overlapping conserved quantities reconstruct consistently (recycle enzyme)", {
  skip_if_no_compile()
  # A closed catalytic cycle with two overlapping CQs; one has a negative
  # coefficient, so reconstruction divides by it and nests another eliminated species.
  m <- pp_models()

  totals <- getTotals(m$el_recycle)
  expect_length(totals, 2L)                       # two independent CQs

  p0 <- c(k1 = 2, k1r = 1, k1c = 3, k3 = 1, total_1 = 1, total_2 = 4,
          G = 0.3, GS = 0.2, S = 2, P = 1)
  o <- m$recycle(p0)[[1]]

  # all four species reconstructed and non-negative at the steady state
  expect_setequal(intersect(c("G","GS","S","P"), names(o)), c("G","GS","S","P"))
  expect_true(all(as.numeric(o[c("G","GS","S","P")]) >= -1e-6))
  # the enzyme moiety G + GS equals its total to the solver tolerance
  expect_equal(as.numeric(o["G"] + o["GS"]), p0[["total_1"]], tolerance = 1e-5)
  expect_true(is.finite(as.numeric(o["S"])))
})

test_that("a fully open network (all states drain to zero) errors clearly", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  # No influx anywhere; every species is structurally zero in steady state.
  el <- eqnlist()
  el <- addReaction(el, "E + S", "ES", "k1*E*S")
  el <- addReaction(el, "ES", "E + S", "k1r*ES")
  el <- addReaction(el, "ES", "E + P", "k1c*ES")

  expect_error(Pimpl(el, compile = FALSE), "structurally zero in steady state")
})


## ---- Pimpl: singular df/dx triggers a diagnostic error ----------------

test_that("Pimpl uses the pseudoinverse when df/dx is rank-deficient", {
  skip_if_no_compile()
  # Duplicate residuals make df/dx rank-deficient. Pimpl warns and returns the
  # minimum-norm sensitivity from the Moore-Penrose pseudoinverse.
  pf <- pp_models()$sing

  p0 <- c(s = 1, x1 = 0.5, x2 = 0.5)
  expect_warning(out <- pf(p0, deriv = TRUE)[[1]], "rank-deficient")
  J <- attr(out, "deriv")
  expect_true(is.matrix(J))
  # The sensitivity satisfies the differentiated constraint.
  expect_equal(sum(J[c("x1", "x2"), "s"]), 1, tolerance = 1e-6)
})


## ---- resetWarmStarts ---------------------------------------------------

test_that("resetWarmStarts clears Pimpl's cache by name", {
  skip_if_no_compile()
  pf <- pp_models()$lin
  pf(c(a = 1, x = 0.5))

  reset_env <- environment(attr(pf, "resetWarmStart"))$reg_ref$caches[["__default__"]]
  expect_false(is.null(reset_env$arch))

  labels <- resetWarmStarts(pf, verbose = FALSE)
  expect_length(labels, 1L)
  expect_match(labels, "^Pimpl\\(")
  expect_null(reset_env$arch)
})


test_that("resetWarmStarts walks into composed functions", {
  skip_if_no_compile()
  # resetWarmStarts() walks this frame through `wrap`, so no other
  # warm-starting object may sit in it.
  p1 <- pp_models()$pd
  p2 <- pp_models()$pd2

  p1(c(k_in = 1, k_out = 0.5))
  p2(c(k_in = 2, k_out = 0.4))

  ref1 <- environment(attr(p1, "resetWarmStart"))$reg_ref$caches[["__default__"]]
  ref2 <- environment(attr(p2, "resetWarmStart"))$reg_ref$caches[["__default__"]]
  expect_false(is.null(ref1$arch))
  expect_false(is.null(ref2$arch))

  wrap <- function(pars) list(p1 = p1(pars), p2 = p2(pars))
  labels <- resetWarmStarts(wrap, verbose = FALSE)
  expect_setequal(sub("\\(.*", "", labels), c("Pimpl", "Pimpl"))
  expect_length(labels, 2L)
  expect_null(ref1$arch)
  expect_null(ref2$arch)
})


test_that("a condition-less Pimpl keeps an independent warm start per condition", {
  skip_if_no_compile()
  m <- pp_models()
  resetWarmStarts(m$pd, verbose = FALSE)

  # each condition reaches its own closed-form steady state
  pf  <- m$pd * m$px
  out <- pf(c(s = 3), deriv = FALSE)
  expect_equal(as.numeric(out$C1["A"]), 3, tolerance = 1e-10)
  expect_equal(as.numeric(out$C2["A"]), 6, tolerance = 1e-10)

  caches <- environment(attr(m$pd, "resetWarmStart"))$reg_ref$caches
  expect_true(all(c("C1", "C2") %in% ls(caches)))
  expect_equal(as.numeric(caches[["C1"]]$arch[[1]]$x["A"]), 3, tolerance = 1e-10)
  expect_equal(as.numeric(caches[["C2"]]$arch[[1]]$x["A"]), 6, tolerance = 1e-10)

  # the next call starts each condition from its own root
  st <- environment(attr(m$pd, "mappings")[[1]])$stats
  how0 <- st$how
  pf(c(s = 3.3), deriv = FALSE)
  expect_equal(unname(st$how["warm_newton"]) - (if (is.na(how0["warm_newton"])) 0 else unname(how0["warm_newton"])), 2)
})

test_that("conditions with identical steady-state inputs are solved once", {
  skip_if_no_compile()
  m <- pp_models()
  resetWarmStarts(m$pd, verbose = FALSE)
  E <- environment(attr(m$pd, "mappings")[[1]]); E$statsReset()
  out <- (m$pd * m$px_same)(c(s = 2), deriv = TRUE)
  expect_equal(as.numeric(out$C2["A"]), 2, tolerance = 1e-10)
  expect_equal(E$stats$calls, 2L)
  expect_equal(E$stats$solves, 1L)
  expect_equal(E$stats$memo, 1L)
})


test_that("resetWarmStarts is idempotent and silent when nothing to reset", {
  f <- function(x) x + 1
  labels <- resetWarmStarts(f, verbose = FALSE)
  expect_length(labels, 0L)
})


test_that("resetWarmStarts rejects non-function inputs", {
  expect_error(resetWarmStarts(1L), "must be a function")
  expect_error(resetWarmStarts(NULL), "must be a function")
})


## ---- Pimpl hard error when residual cannot be driven below ftol -------

test_that("Pimpl throws when no attempt converges", {
  skip_if_no_compile()
  # a residual without a real root
  pf <- pp_models()$noroot

  resetWarmStarts(pf, verbose = FALSE)
  expect_error(pf(c(x = 0.1))[[1]], "no root found")
})


# ==== Structural zero states (.zeroStatesFromSmatrix): NegCol, PosCol, sink-cluster LP ====

helper <- dMod2:::.zeroStatesFromSmatrix


## ---- Layer 1: only-outflux column --------------------------------------

test_that("NegCol: a state with only outflux is detected as zero", {
  el <- eqnlist() |>
    addReaction("A", "", "k_dg * A")
  zs <- helper(el)
  expect_setequal(zs$zero_states, "A")
  expect_equal(ncol(zs$eqnlist$smatrix), 0L)
})


## ---- Layer 2: PosCol with single-state feeder --------------------------

test_that("PosCol: only-influx state whose feeder has a single reactant", {
  # An only-influx state forces its single-reactant feeder to zero, and the
  # state itself drops on the next iteration.
  el <- eqnlist() |>
    addReaction("", "X", "k_pr * Z") |>
    addReaction("Z", "",  "k_dg * Z")
  zs <- helper(el)
  expect_setequal(zs$zero_states, c("Z", "X"))
})


## ---- Layer 3: structural sink cluster ----------------------------------

test_that("Sink cluster: TGFb + R_TGFb + R_TGFb_int (combined mass leaks)", {
  skip_if_not_installed("lpSolve")
  # Every column has both signs, so layers 1 and 2 miss the cluster; its total
  # mass leaks through one degradation.
  el <- eqnlist() |>
    addReaction("L + R", "RL",  "k_on  * L * R") |>
    addReaction("RL",    "L + R", "k_off * RL") |>
    addReaction("RL",    "RLi", "k_int * RL") |>
    addReaction("RLi",   "L",   "k_dec * RLi") |>
    addReaction("RLi",   "",    "k_dg  * RLi") |>
    addReaction("",      "R",   "k_pr_R") |>
    addReaction("R",     "",    "k_dg_R * R")
  zs <- helper(el)
  expect_true(all(c("L", "RL", "RLi") %in% zs$zero_states))
  # R is not in the sink cluster: it has its own production/degradation.
  expect_false("R" %in% zs$zero_states)
})


test_that("Sink cluster: a conserved moiety listed first stays nonzero", {
  skip_if_not_installed("lpSolve")
  # A <-> Ap is conserved and comes first; the leaking cluster {L, RL} comes
  # after it. Their union also passes the LP, but only {L, RL} drains.
  el <- eqnlist() |>
    addReaction("A",     "Ap",  "k1 * A") |>
    addReaction("Ap",    "A",   "k2 * Ap") |>
    addReaction("",      "R",   "k_pr_R") |>
    addReaction("R",     "",    "k_dg_R * R") |>
    addReaction("L + R", "RL",  "k_on * L * R") |>
    addReaction("RL",    "L + R", "k_off * RL") |>
    addReaction("RL",    "",    "k_dg * RL")
  zs <- helper(el)
  expect_setequal(zs$zero_states, c("L", "RL"))
})


## ---- Idempotence on a model with no zero states ------------------------

test_that("Pure production-decay has no zero states", {
  el <- eqnlist() |>
    addReaction("", "A", "k_in") |>
    addReaction("A", "",  "k_out * A")
  zs <- helper(el)
  expect_equal(zs$zero_states, character(0))
  expect_identical(zs$eqnlist, el)
})


test_that("Pimpl reports zero states with value 0 and drops them from getParameters()", {
  skip_if_no_compile()
  # Receptor R alone has production + degradation (nonzero baseline).
  # Ligand L and complex LR form a sink cluster (LR degrades).
  pf <- pp_models()$zero

  expect_false(any(c("L", "LR") %in% getParameters(pf)))
  pars <- c(k_pr_R = 2, k_dg_R = 0.5, k_on = 1, k_off = 0.5, k_dg_LR = 0.1)
  out  <- pf(pars)[[1]]
  expect_identical(as.numeric(out[c("L", "LR")]), c(0, 0))
  expect_equal(as.numeric(out["R"]), pars[["k_pr_R"]] / pars[["k_dg_R"]], tolerance = 1e-10)
})

test_that("a state without influx at the given parameters is exactly 0", {
  skip_if_no_compile()
  pf <- pp_models()$plasmid
  p0 <- c(k_tx = 2, k_dg = 0.5, k_y = 1, k_xy = 3, k_dy = 0.25)

  out <- pf(c(p0, plasmid = 0))[[1]]
  expect_identical(as.numeric(out["X"]), 0)
  expect_equal(as.numeric(out["Y"]), 4, tolerance = 1e-10)
  # the sensitivity to the plasmid is still the one at plasmid > 0
  J <- attr(out, "deriv")
  expect_equal(unname(J["X", "plasmid"]), 4, tolerance = 1e-10)
  expect_equal(unname(J["Y", "plasmid"]), 48, tolerance = 1e-10)

  out1 <- pf(c(p0, plasmid = 1))[[1]]
  expect_equal(as.numeric(out1[c("X", "Y")]), c(4, 52), tolerance = 1e-10)
})

test_that("Pimpl returns a stable steady state of a bistable system", {
  skip_if_no_compile()
  pf <- pp_models()$bistable
  p0 <- c(k0 = 0.02, k1 = 1, K = 1, d = 0.4)
  f  <- function(x) p0[["k0"]] + p0[["k1"]] * x^2 / (p0[["K"]]^2 + x^2) - p0[["d"]] * x
  lo <- uniroot(f, c(0.01, 0.1), tol = 1e-14)$root
  mi <- uniroot(f, c(0.1, 1), tol = 1e-14)$root
  hi <- uniroot(f, c(1, 3), tol = 1e-14)$root

  root <- function(x0) {
    resetWarmStarts(pf, verbose = FALSE)
    as.numeric(pf(c(p0, X = x0), deriv = FALSE)[[1]]["X"])
  }
  expect_equal(root(0.8 * mi), lo, tolerance = 1e-8)
  expect_equal(root(1.2 * mi), hi, tolerance = 1e-8)
  # started on the unstable root, it leaves along the unstable direction
  expect_true(abs(root(mi) - lo) < 1e-8 || abs(root(mi) - hi) < 1e-8)
})

test_that("Pimpl with flow = TRUE refuses an unstable root", {
  skip_if_no_compile()
  pf <- pp_models()$unstable
  resetWarmStarts(pf, verbose = FALSE)
  expect_error(pf(c(a = 2, x = 1)), "no stable root found")
})


## ---- No-progress: hard error instead of stale initial values -----------

# ==== Edge case: Pimpl with no outer parameters ====

test_that("Pimpl with no outer parameters does not crash in build_jacobian", {
  skip_if_no_compile()
  p <- pp_models()$noparam
  out <- p(c(dummy = 1.0))
  expect_true(is.numeric(unclass(out[[1]])["A"]))
  expect_equal(unname(unclass(out[[1]])["A"]), 1.0, tolerance = 1e-3)
})


# ==== Totals as parameters: `totalXxx` from the longest common substring, else `total_<index>` ====

## ---- Smart naming: pERK + ERK -> totalERK ----------------------------

test_that(".smartTotalName picks the longest common substring", {
  expect_equal(dMod2:::.smartTotalName(c("pERK", "ERK"), character(0),
                                        character(0), 1), "totalERK")
  expect_equal(dMod2:::.smartTotalName(c("TGFb", "R1_TGFb", "R1_TGFb_int"),
                                        character(0), character(0), 1),
               "totalTGFb")
  # No common substring -> fallback to total_<index>
  expect_equal(dMod2:::.smartTotalName(c("A", "B"), character(0),
                                        character(0), 1), "total_1")
  # Single character common -> below the >=2 threshold, fallback
  expect_equal(dMod2:::.smartTotalName(c("Ax", "Ay"), character(0),
                                        character(0), 1), "total_1")
  # Collision with an existing parameter -> disambiguate
  expect_equal(dMod2:::.smartTotalName(c("pERK", "ERK"), character(0),
                                        "totalERK", 1), "totalERK_2")
})


## ---- Pimpl smart naming: ERK + pERK -> totalERK ---------------------

test_that("P(method='implicit') preserves the eqnlist smatrix for CQ detection", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  el <- eqnlist() |>
    addReaction("ERK",  "pERK", "k1 * ERK") |>
    addReaction("pERK", "ERK",  "k2 * pERK")

  pf <- P(el, method = "implicit", compile = FALSE, verbose = FALSE,
          modelname = paste0("test_P_implicit_cq_", as.integer(Sys.time())))
  expect_true("totalERK" %in% getParameters(pf))
})


test_that("Pimpl on ERK <-> pERK introduces totalERK via LCS", {
  skip_if_no_compile()
  pf <- pp_models()$erk_im
  expect_true("totalERK" %in% getParameters(pf))
  expect_false("total_1" %in% getParameters(pf))
})


test_that("Pimpl solves ERK <-> pERK to closed form in totals", {
  skip_if_no_compile()
  pf <- pp_models()$erk_im
  pars <- c(k1 = 1, k2 = 3, totalERK = 4)
  out  <- pf(pars, deriv = FALSE)[[1]]
  # states match the closed form in the total
  expect_equal(as.numeric(out[c("ERK", "pERK")]), c(3, 1), tolerance = 1e-10)
})


test_that("Pimpl runs its random starts from a fixed seed and leaves the global RNG alone", {
  skip_if_no_compile()
  pf <- pp_models()$erk_im
  on.exit(controls(pf, name = "controlsPTC") <- list(), add = TRUE)
  pars <- c(k1 = 1, k2 = 3, totalERK = 4)
  set.seed(7); before <- .Random.seed
  resetWarmStarts(pf, verbose = FALSE)
  controls(pf, name = "controlsPTC") <- list(maxit = 1L, nStarts = 3L)
  expect_error(pf(pars, deriv = FALSE), "start3: ")
  expect_identical(.Random.seed, before)
  controls(pf, name = "controlsPTC") <- list(nStarts = 3L)
  out <- pf(pars, deriv = FALSE)[[1]]
  expect_equal(as.numeric(out[c("ERK", "pERK")]), c(3, 1), tolerance = 1e-10)
  expect_identical(.Random.seed, before)

  # one stream per warm-start cache: per condition, gone after resetWarmStarts()
  env <- environment(attr(pf, "mappings")[[1]]); reg <- env$reg
  resetWarmStarts(pf, verbose = FALSE)
  controls(pf, name = "controlsPTC") <- list(maxit = 1L, nStarts = 2L)
  for (cn in c("a", "b")) expect_error(env$p2p(pars, deriv = FALSE, condition = cn), "start2: ")
  expect_false(identical(reg$get("a")$msBase, reg$get("b")$msBase))
  expect_identical(reg$get("a")$msCount, 1L)
  resetWarmStarts(pf, verbose = FALSE)
  expect_null(reg$get("a")$msBase)
})


test_that("Pimpl solves a network whose R-Smad rows barely turn over", {
  # A row that barely turns over beside a fast row of the same flux must not make
  # the diagonally scaled step matrix singular.
  skip_if_no_compile()
  fx <- readRDS(test_path("fixtures", "pimpl_smad_row.rds"))
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)
  mn <- paste0("test_pp_smadrow_", as.integer(Sys.time()))
  pf <- Pimpl(fx$reactions, forcings = fx$forcings, deriv = FALSE, modelname = mn,
              controlsPTC = list(nStarts = 0L))
  compile(pf, output = mn, cores = test_cores())
  out <- pf(fx$pv, deriv = FALSE)[[1]]
  expect_equal(out[["Smad2"]] + out[["pSmad2"]] + out[["C234"]], fx$pv[["tSmad2"]],
               tolerance = 1e-10)
})


test_that("Pimpl rejects unknown arguments and controls", {
  el <- eqnlist() |>
    addReaction("ERK",  "pERK", "k1 * ERK") |>
    addReaction("pERK", "ERK",  "k2 * pERK")
  expect_error(Pimpl(el, expressInTotals = FALSE), "unused argument")
  expect_error(Pimpl(el, controlsPTC = list(ftol = 1e-9)), "unknown controlsPTC entry ftol")
  ctl <- dMod2:::.pimplPTC
  expect_identical(ctl(list())$startScale, "log10")
  expect_identical(ctl(list(positive = FALSE))$startScale, "linear")
  expect_error(ctl(list(startScale = "linear", startRange = c(-1, 1))), "needs positive = FALSE")
  expect_error(ctl(list(startRange = c(5, -5))), "lower < upper")
})


# ==== CQ basis as an eqnlist field: getTotals(), customTotals(), mutators, pre-totals eqnlists ====

## ---- getTotals auto-detection ---------------------------------------

test_that("getTotals returns smart-named auto-detected totals", {
  el_AB <- eqnlist() |>
    addReaction("A", "B", "k1 * A") |>
    addReaction("B", "A", "k2 * B")
  tot_AB <- getTotals(el_AB)
  expect_length(tot_AB, 1L)
  expect_equal(names(tot_AB), "total_1")  # LCS of {A, B} is "" -> fallback

  el_ERK <- eqnlist() |>
    addReaction("ERK",  "pERK", "k1 * ERK") |>
    addReaction("pERK", "ERK",  "k2 * pERK")
  tot_ERK <- getTotals(el_ERK)
  expect_length(tot_ERK, 1L)
  expect_equal(names(tot_ERK), "totalERK")
})


## ---- customTotals strict validation ---------------------------------

test_that("customTotals validates structure, rank, and CQ membership", {
  el <- eqnlist() |>
    addReaction("ERK",  "pERK", "k1 * ERK") |>
    addReaction("pERK", "ERK",  "k2 * pERK")

  el2 <- customTotals(el, list(totalE = "ERK + pERK"))
  tot <- getTotals(el2)
  expect_equal(names(tot), "totalE")
  expect_equal(tot[["totalE"]], "ERK + pERK")
  expect_true(isTRUE(attr(el2$totals, "custom")))

  expect_error(customTotals(el, list(bogus = "ERK")),
               "not a conservation quantity")
  expect_error(customTotals(el, list(too = "ERK + pERK", many = "ERK + pERK")),
               "Expected 1 conservation quantit(y|ies), got 2")
  expect_error(customTotals(el, list(weird = "ERK * pERK")),
               "not linear")

  el_reset <- customTotals(el2, NULL)
  expect_null(el_reset$totals)
  el_reset2 <- customTotals(el2, list())
  expect_null(el_reset2$totals)
})


## ---- Mutator preservation of custom totals --------------------------

test_that("addReaction preserves custom totals when CQ structure survives", {
  el <- eqnlist() |>
    addReaction("ERK",  "pERK", "k1 * ERK") |>
    addReaction("pERK", "ERK",  "k2 * pERK")
  el <- customTotals(el, list(totalE = "ERK + pERK"))

  # a reaction outside the moiety keeps the totals
  el2 <- addReaction(el, "X", "", "k3 * X")
  expect_equal(names(el2$totals), "totalE")
  expect_true(isTRUE(attr(el2$totals, "custom")))
})


test_that("addReaction warns and resets when custom totals are invalidated", {
  el <- eqnlist() |>
    addReaction("ERK",  "pERK", "k1 * ERK") |>
    addReaction("pERK", "ERK",  "k2 * pERK")
  el <- customTotals(el, list(totalE = "ERK + pERK"))

  # Adding a degradation of ERK breaks the conservation
  expect_warning(
    el2 <- addReaction(el, "ERK", "", "k_dg * ERK"),
    "customTotals invalidated"
  )
  expect_null(el2$totals)
})


## ---- Backward compat with pre-totals eqnlist ------------------------

test_that("is.eqnlist accepts pre-totals eqnlists (missing $totals field)", {
  el <- eqnlist() |>
    addReaction("A", "B", "k1 * A") |>
    addReaction("B", "A", "k2 * B")
  el_old <- el; el_old$totals <- NULL
  el_old <- el_old[setdiff(names(el_old), "totals")]
  class(el_old) <- c("eqnlist", "list")
  expect_true(is.eqnlist(el_old))
  # getTotals still works (computes lazily)
  expect_length(getTotals(el_old), 1L)
})


test_that("Pimpl uses customTotals names in the parvec interface", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  el <- eqnlist() |>
    addReaction("ERK",  "pERK", "k1 * ERK") |>
    addReaction("pERK", "ERK",  "k2 * pERK")
  el <- customTotals(el, list(totalERKpool = "ERK + pERK"))

  pf <- Pimpl(el, parameters = c("k1", "k2"),
              modelname = paste0("test_Pimpl_custom_totals_",
                                 as.integer(Sys.time())),
              compile = FALSE, verbose = FALSE)
  params <- getParameters(pf)
  expect_true("totalERKpool" %in% params)
  expect_false(any(grepl("^total[^E]", params)))  # no auto-named total slipped in
})


## ---- print.eqnlist shows totals section -----------------------------

test_that("print.eqnlist shows the conserved quantities by name", {
  el <- eqnlist() |>
    addReaction("ERK",  "pERK", "k1 * ERK") |>
    addReaction("pERK", "ERK",  "k2 * pERK")
  out <- capture.output(print(el))
  expect_true(any(grepl("Conserved quantities", out)))
  expect_true(any(grepl("totalERK", out)))

  el2 <- customTotals(el, list(totalE = "ERK + pERK"))
  out2 <- capture.output(print(el2))
  expect_true(any(grepl("Conserved quantities .custom.", out2)))
  expect_true(any(grepl("totalE\\s*=\\s*ERK \\+ pERK", out2)))
})
