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
# gets a separate shared object.
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
              parameters = c("k1","km","totA","totB"), derivMode = "forward-forward",
              modelname = mn("d2"), verbose = FALSE, controlsPTC = list(reltol = 1e-13))

  d2_2 <- Pimpl(c(x1 = "a*x1 - b*x2 - 1", x2 = "x1*x2 - c"),
                parameters = c("a","b","c"), derivMode = "forward-forward",
                modelname = mn("d2_2"), verbose = FALSE, controlsPTC = list(reltol = 1e-13))

  el_AB <- eqnlist() |>
    addReaction("A", "B", "k*A") |>
    addReaction("B", "A", "km*B")
  d2_cq <- Pimpl(el_AB, parameters = c("k","km"), derivMode = "forward-forward",
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
                modelname = mn("sing"), verbose = FALSE,
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

  # a root whose first derivative vanishes at a = 1 while the second does not
  kink <- Pimpl(c(x = "(a-1)^2 + 1 - x"), parameters = "a", derivMode = "forward-forward",
                modelname = mn("kink"), verbose = FALSE,
                controlsPTC = list(positive = FALSE))

  # E is a catalyst, constant rather than zero
  el_cat <- eqnlist() |>
    addReaction("", "S", "k0") |> addReaction("S", "", "kd*S") |>
    addReaction("S + E", "P + E", "k1*E*S") |> addReaction("P", "S", "k2*P")
  catalyst <- Pimpl(el_cat, modelname = mn("catalyst"), verbose = FALSE)

  # a moiety across two compartments, one of them of symbolic volume
  el_vol <- eqnlist() |>
    assignCompartment(B = "nuc", volume = "V") |>
    addReaction("A", "B", "k1*A") |> addReaction("B", "A", "k2*B")
  vol <- Pimpl(el_vol, modelname = mn("vol"), derivMode = "forward-forward", verbose = FALSE)

  # a conserved cycle driven by a species that is structurally zero at the steady state
  el_rest <- eqnlist() |>
    addReaction("L", "L + R", "u * L") |>
    addReaction("I", "pI", "a * R * I") |> addReaction("pI", "iI", "b * pI") |>
    addReaction("iI", "I", "c * iI")
  rest <- Pimpl(el_rest, modelname = mn("rest"), verbose = FALSE)
  # the same network as g * x * rest * p, the resting state absent from p or given by it
  rest_ode <- odemodel(el_rest, modelname = mn("rest_ode"), compile = FALSE)
  rest_obs <- Y(c(yI = "I", yP = "pI"), f = el_rest, modelname = mn("rest_obs"))
  tr_rest <- c(u = "exp(lu)", a = "exp(la)", b = "exp(lb)", c = "1",
               setNames("exp(lT)", grep("^total", attr(rest, "parameters"), value = TRUE)))
  rest_p0 <- P(tr_rest, condition = "C1", modelname = mn("rest_p0"))
  rest_p1 <- P(c(tr_rest, R = "exp(lR)"), condition = "C1", modelname = mn("rest_p1"))

  # the Brusselator: x = a, y = b / a, with eigenvalues of real part (b - 1 - a^2) / 2 while complex
  bruss <- Pimpl(c(x = "a - (b + 1) * x + x^2 * y", y = "b * x - x^2 * y"), modelname = mn("bruss"),
                 derivMode = "forward-forward", verbose = FALSE)
  bruss_p <- P(c(a = "exp(la)", b = "exp(lb)"), modelname = mn("bruss_p"))
  # a real eigenvalue through 0 where two roots exchange stability, next to a mode that stays
  trans <- Pimpl(c(x = "x * (q - x)", y = "x - y"), parameters = "q", modelname = mn("trans"),
                 derivMode = "forward-forward", verbose = FALSE, controlsPTC = list(positive = FALSE))
  # the same steady state as a flow with a stable root required, passing a parameter on, which
  # differs between two conditions
  bruss_s <- Pimpl(c(x = "a - (b + 1) * x + x^2 * y", y = "b * x - x^2 * y"), parameters = c("a", "b", "s"),
                   flow = TRUE, modelname = mn("bruss_s"), derivMode = "forward-forward", verbose = FALSE)
  bruss_c <- P(list(C1 = c(a = "exp(la)", b = "exp(lb)", s = "1"), C2 = c(a = "exp(la)", b = "exp(lb)", s = "2")),
               modelname = mn("bruss_c"))
  # two states in a chain whose rates can coincide, a defective double eigenvalue, beside a fast mode
  chain <- Pimpl(c(x = "k0 - k1 * x", y = "k1 * x - k2 * y + c * z", z = "kf0 - kf * z + d * y"),
                 parameters = c("k0", "k1", "k2", "kf0", "kf", "c", "d"), modelname = mn("chain"),
                 derivMode = "forward-forward", verbose = FALSE)

  compile(lin, d2, d2_2, d2_cq, moiety2_im, dimer_im, recycle, sing, noroot,
          noparam, erk_im, px, px_same, zero, pd, pd2, plasmid, bistable, unstable,
          kink, catalyst, vol, rest, rest_ode, rest_obs, rest_p0, rest_p1, bruss, bruss_p,
          bruss_s, bruss_c, trans, chain, output = mn("impl"), cores = test_cores())

  .pp_env$models <- list(
    lin = lin, d2 = d2, d2_2 = d2_2, d2_cq = d2_cq, moiety2_im = moiety2_im,
    dimer_im = dimer_im, recycle = recycle, el_recycle = el_recycle, sing = sing,
    noroot = noroot, noparam = noparam, erk_im = erk_im, px = px, px_same = px_same,
    zero = zero, pd = pd, pd2 = pd2, plasmid = plasmid, bistable = bistable,
    unstable = unstable, kink = kink, catalyst = catalyst, vol = vol, rest = rest,
    rest_ode = rest_ode, rest_obs = rest_obs, rest_p0 = rest_p0, rest_p1 = rest_p1,
    bruss = bruss, bruss_p = bruss_p, bruss_s = bruss_s, bruss_c = bruss_c, trans = trans,
    chain = chain)
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
  # Dimerization conserves the monomer count.
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

  totals <- conservedQuantities(m$el_recycle)
  expect_length(totals, 2L)                       # two independent CQs

  p0 <- c(k1 = 2, k1r = 1, k1c = 3, k3 = 1, total_1 = 1, total_2 = 4,
          G = 0.3, GS = 0.2, S = 2, P = 1)
  o <- m$recycle(p0)[[1]]

  # all four species reconstructed and non-negative at the steady state
  expect_setequal(intersect(c("G","GS","S","P"), names(o)), c("G","GS","S","P"))
  expect_true(all(as.numeric(o[c("G","GS","S","P")]) >= -1e-6))
  # the enzyme moiety G + GS equals the moiety total to the solver tolerance
  expect_equal(as.numeric(o["G"] + o["GS"]), p0[["total_1"]], tolerance = 1e-5)
  expect_true(is.finite(as.numeric(o["S"])))
})

test_that("a network without influx keeps the conserved totals, an open one errors clearly", {
  skip_if_no_compile()
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)

  # species consumed without influx are zero, the others hold their conserved totals
  el <- eqnlist()
  el <- addReaction(el, "E + S", "ES", "k1*E*S")
  el <- addReaction(el, "ES", "E + S", "k1r*ES")
  el <- addReaction(el, "ES", "E + P", "k1c*ES")
  n <- dMod2:::.normalize_ss_inputs(el, NULL, NULL)
  expect_setequal(n$zero_states, c("S", "ES"))
  expect_setequal(n$dependent, c("E", "P"))
  expect_length(n$totals, 2L)

  # every species drains
  open <- eqnlist() |> addReaction("A", "", "k1*A") |> addReaction("A", "B", "k2*A") |>
    addReaction("B", "", "k3*B")
  expect_error(Pimpl(open, compile = FALSE), "structurally zero in steady state")
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


## ---- Kept roots --------------------------------------------------------

test_that("keepRoot answers repeated values from memory with the same root", {
  skip_if_no_compile()
  pf <- pp_models()$lin
  E <- environment(attr(pf, "mappings")[[1]])
  on.exit(controls(pf, name = "keepRoot") <- TRUE, add = TRUE)
  pars <- c(a = 1.25, x = 0.5)
  E$statsReset()
  first <- pf(pars); again <- pf(pars)
  expect_identical(again, first)
  expect_equal(c(E$stats$solves, E$stats$memo), c(1L, 1L))
  controls(pf, name = "keepRoot") <- FALSE
  E$statsReset()
  expect_identical(pf(pars), first)
  expect_equal(c(E$stats$solves, E$stats$memo), c(1L, 0L))
})


test_that("a condition-less Pimpl solves each condition from the condition values", {
  skip_if_no_compile()
  m <- pp_models()
  out <- (m$pd * m$px)(c(s = 3), deriv = FALSE)
  expect_equal(as.numeric(out$C1["A"]), 3, tolerance = 1e-10)
  expect_equal(as.numeric(out$C2["A"]), 6, tolerance = 1e-10)
})

test_that("a root does not depend on the calls before it", {
  skip_if_no_compile()
  pf <- pp_models()$bistable
  p0 <- c(k0 = 0.02, k1 = 1, K = 1, d = 0.4)
  f  <- function(x) p0[["k0"]] + p0[["k1"]] * x^2 / (p0[["K"]]^2 + x^2) - p0[["d"]] * x
  lo <- uniroot(f, c(0.01, 0.1), tol = 1e-14)$root
  # first a root on another branch, at other parameters
  on.exit(controls(pf, name = "keepRoot") <- TRUE, add = TRUE)
  pf(c(replace(p0, "d", 0.2), X = 0.1), deriv = FALSE)
  after <- pf(c(p0, X = 0.1), deriv = FALSE)[[1]][["X"]]
  controls(pf, name = "keepRoot") <- FALSE
  fresh <- pf(c(p0, X = 0.1), deriv = FALSE)[[1]][["X"]]
  expect_identical(after, fresh)
  expect_equal(fresh, lo, tolerance = 1e-10)
})

test_that("conditions with identical steady-state inputs are solved once", {
  skip_if_no_compile()
  m <- pp_models()
  E <- environment(attr(m$pd, "mappings")[[1]]); E$statsReset()
  out <- (m$pd * m$px_same)(c(s = 2.375), deriv = TRUE)
  expect_equal(as.numeric(out$C2["A"]), 2.375, tolerance = 1e-10)
  expect_equal(E$stats$calls, 2L)
  expect_equal(E$stats$solves, 1L)
  expect_equal(E$stats$memo, 1L)
})


## ---- Pimpl hard error when residual cannot be driven below ftol -------

test_that("Pimpl throws when no attempt converges", {
  skip_if_no_compile()
  # a residual without a real root
  pf <- pp_models()$noroot

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
  # An only-influx state forces the single-reactant feeder to zero; the state
  # itself is left without reactions and rests.
  el <- eqnlist() |>
    addReaction("", "X", "k_pr * Z") |>
    addReaction("Z", "",  "k_dg * Z")
  zs <- helper(el)
  expect_setequal(zs$zero_states, "Z")
  expect_identical(zs$resting, "X")
})


## ---- Layer 3: structural sink cluster ----------------------------------

test_that("Sink cluster: TGFb + R_TGFb + R_TGFb_int (combined mass leaks)", {
  skip_if_not_installed("lpSolve")
  # Every column has both signs, so layers 1 and 2 miss the cluster; the total
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
  # R is not in the sink cluster: it has a separate production/degradation.
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
  # A held species names the total after the species
  held <- addReaction(eqnlist(), "Erk + u", "pErk + u", "k*Erk*u") |>
    addReaction("pErk", "Erk", "kb*pErk")
  expect_setequal(names(conservedQuantities(held)), c("totalErk", "u"))
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


test_that("Pimpl runs the random starts from a fixed seed and leaves the global RNG alone", {
  skip_if_no_compile()
  pf <- pp_models()$erk_im
  on.exit(controls(pf, name = "controlsPTC") <- list(), add = TRUE)
  pars <- c(k1 = 1, k2 = 3, totalERK = 4)
  set.seed(7); before <- .Random.seed
  controls(pf, name = "controlsPTC") <- list(maxit = 1L, nStarts = 3L)
  expect_error(pf(pars, deriv = FALSE), "start3: ")
  expect_identical(.Random.seed, before)
  controls(pf, name = "controlsPTC") <- list(nStarts = 3L)
  out <- pf(pars, deriv = FALSE)[[1]]
  expect_equal(as.numeric(out[c("ERK", "pERK")]), c(3, 1), tolerance = 1e-10)
  expect_identical(.Random.seed, before)

  # the random starts depend on the values, not on the condition or the calls before
  env <- environment(attr(pf, "mappings")[[1]])
  controls(pf, name = "controlsPTC") <- list(maxit = 1L, nStarts = 2L)
  msg <- function(cn) tryCatch(env$p2p(pars, deriv = FALSE, condition = cn), error = conditionMessage)
  expect_match(msg("a"), "start2: ")
  expect_identical(msg("a"), msg("b"))
  expect_identical(.Random.seed, before)
})


test_that("Pimpl solves a network whose R-Smad rows barely turn over", {
  # A row that barely turns over beside a fast row of the same flux must not make
  # the diagonally scaled step matrix singular.
  skip_if_no_compile()
  fx <- readRDS(test_path("fixtures", "pimpl_smad_row.rds"))
  oldwd <- setwd(.dmod_fx_workdir()); on.exit(setwd(oldwd), add = TRUE)
  mn <- paste0("test_pp_smadrow_", as.integer(Sys.time()))
  pf <- Pimpl(fx$reactions, forcings = fx$forcings, derivMode = "none", modelname = mn,
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


# ==== CQ basis as an eqnlist field: conservedQuantities(), customTotals(), mutators, pre-totals eqnlists ====

## ---- conservedQuantities auto-detection ---------------------------------------

test_that("conservedQuantities returns smart-named auto-detected totals", {
  el_AB <- eqnlist() |>
    addReaction("A", "B", "k1 * A") |>
    addReaction("B", "A", "k2 * B")
  tot_AB <- conservedQuantities(el_AB)
  expect_length(tot_AB, 1L)
  expect_equal(names(tot_AB), "total_1")  # LCS of {A, B} is "" -> fallback

  el_ERK <- eqnlist() |>
    addReaction("ERK",  "pERK", "k1 * ERK") |>
    addReaction("pERK", "ERK",  "k2 * pERK")
  tot_ERK <- conservedQuantities(el_ERK)
  expect_length(tot_ERK, 1L)
  expect_equal(names(tot_ERK), "totalERK")
})


## ---- customTotals strict validation ---------------------------------

test_that("customTotals validates structure, rank, and CQ membership", {
  el <- eqnlist() |>
    addReaction("ERK",  "pERK", "k1 * ERK") |>
    addReaction("pERK", "ERK",  "k2 * pERK")

  el2 <- customTotals(el, list(totalE = "ERK + pERK"))
  tot <- conservedQuantities(el2)
  expect_equal(names(tot), "totalE")
  expect_equal(tot[["totalE"]], "ERK + pERK")
  expect_true(isTRUE(attr(el2$totals, "custom")))

  expect_error(customTotals(el, list(bogus = "ERK")),
               "not a conserved quantity")
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
  # conservedQuantities still works (computes lazily)
  expect_length(conservedQuantities(el_old), 1L)
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


test_that("a solved state keeps the second derivative where the first vanishes", {
  skip_if_no_compile()
  w <- pp_models()$kink(c(a = 1), deriv2 = TRUE)[[1]]
  expect_equal(unname(w[["x"]]), 1, tolerance = 1e-8)
  expect_false("x" %in% attr(w, "fixed"))
  expect_equal(attr(w, "deriv2")["x", "a", "a"], 2, tolerance = 1e-6)
})

test_that("a catalyst is a constant state, not a zero one", {
  skip_if_no_compile()
  pc <- pp_models()$catalyst
  expect_true("E" %in% attr(pc, "parameters"))
  expect_false(any(grepl("^total", attr(pc, "parameters"))))
  given <- c(k0 = 1, kd = 1, k1 = 1, k2 = 2, E = 4)
  pars <- given[intersect(names(given), attr(pc, "parameters"))]
  v <- pc(c(pars, setNames(rep(4, length(setdiff(attr(pc, "parameters"), names(given)))),
                           setdiff(attr(pc, "parameters"), names(given)))),
          deriv = FALSE)[[1]]
  expect_equal(unname(v[c("S", "P", "E")]), c(1, 2, 4), tolerance = 1e-6)
})

test_that("the conserved amount holds the volumes, also symbolic ones", {
  skip_if_no_compile()
  pv <- pp_models()$vol
  tn <- grep("^total", attr(pv, "parameters"), value = TRUE)
  expect_length(tn, 1L)
  pars <- c(k1 = 1.5, k2 = 0.5, V = 0.4, setNames(2, tn))
  o <- pv(pars, deriv = TRUE)[[1]]
  A <- 2 * 0.5 / 2; B <- 1.5 * A / (0.5 * 0.4)
  expect_equal(as.numeric(o[c("A", "B")]), c(A, B), tolerance = 1e-8)
  J <- attr(o, "deriv")
  expect_equal(unname(J[c("A", "B"), "V"]), c(0, -B / 0.4), tolerance = 1e-7)
})

test_that("a resting state is 0 unless given, and a given value drives the moiety", {
  skip_if_no_compile()
  pr <- pp_models()$rest
  expect_false("R" %in% attr(pr, "parameters"))
  tn <- grep("^total", attr(pr, "parameters"), value = TRUE)
  expect_length(tn, 1L)
  pars <- c(u = 1, a = 1, b = 1, c = 1, setNames(2, tn))
  o <- pr(pars[intersect(names(pars), attr(pr, "parameters"))], deriv = FALSE)[[1]]
  expect_equal(as.numeric(o[c("I", "pI", "iI", "R")]), c(2, 0, 0, 0))
  expect_identical(as.numeric(o[c("pI", "iI")]), c(0, 0))

  o <- pr(c(pars[intersect(names(pars), attr(pr, "parameters"))], R = 0.1))[[1]]
  expect_equal(as.numeric(o[c("I", "pI", "iI", "R")]), c(20, 2, 2, 1.2) / 12, tolerance = 1e-8)
  J <- attr(o, "deriv")
  expect_equal(unname(J[c("I", "pI", "R"), "R"]), c(-4, 2, 1.44) / 1.44, tolerance = 1e-6)
})


test_that("the native chain reads a resting state absent from the input as 0", {
  skip_if_no_compile()
  m <- pp_models()
  x <- Xs(m$rest_ode)
  data <- as.datalist(data.frame(
    name = rep(c("yI", "yP"), each = 3), time = rep(c(0, 1, 2), 2),
    value = c(1.7, 1.6, 1.7, 0.1, 0.2, 0.15), sigma = 0.1, condition = "C1"), splitBy = "condition")
  outer <- c(lu = 0, la = 0.1, lb = -0.2, lT = log(2), lR = log(0.1))
  for (p in list(m$rest_p0, m$rest_p1)) {
    obj <- normL2(data, m$rest_obs * x * m$rest * p)
    pars <- outer[getParameters(p)]
    a <- withr::with_options(list(dMod2.evaluation = "C++"), obj(pars))
    b <- withr::with_options(list(dMod2.evaluation = "R"), obj(pars))
    entries <- environment(obj)$.chain$entries
    expect_false(isTRUE(get(ls(entries)[1], envir = entries)[[1]]$none))
    expect_equal(a$value, b$value, tolerance = 1e-10)
    expect_equal(a$gradient, b$gradient, tolerance = 1e-10)
    expect_equal(a$hessian, b$hessian, tolerance = 1e-10)
  }
  expect_true(abs(a$gradient[["lR"]]) > 1)
})


test_that("the spectrum of a steady state has the exact gradients of the eigenvalue real parts", {
  skip_if_no_compile()
  m <- pp_models()
  sp <- attr(attr(m$bruss, "mappings")[[1]], "spectrum")(c(a = 1, b = 1.5), onset = -1)
  expect_equal(sp$re, c(-0.25, -0.25), tolerance = 1e-10)
  expect_equal(abs(sp$im), rep(sqrt(1 - 0.25^2), 2), tolerance = 1e-10)
  expect_equal(unname(sp$grad[, c("a", "b")]), rbind(c(-1, 0.5), c(-1, 0.5)), tolerance = 1e-8)
  # the rate scale, the largest eigenvalue modulus, with gradient
  sp <- attr(attr(m$bruss, "mappings")[[1]], "spectrum")(c(a = 1.3, b = 2.1), onset = Inf)
  expect_length(sp$active, 0L)
  expect_equal(sp$scale, 1.3, tolerance = 1e-10)
  expect_equal(unname(sp$scaleGrad[c("a", "b")]), c(1, 0), tolerance = 1e-10)
})

test_that("coinciding eigenvalues have the exact gradient of their sum", {
  skip_if_no_compile()
  m <- pp_models()
  spec <- attr(attr(m$chain, "mappings")[[1]], "spectrum")
  # a double eigenvalue with one eigenvector: each takes the mean of the derivatives
  pars <- c(k0 = 2, k1 = 1.3, k2 = 1.3, kf0 = 1, kf = 10, c = 0, d = 0)
  sp <- spec(pars, onset = -Inf)
  slow <- which(abs(sp$re + 1.3) < 1e-8)
  expect_length(slow, 2L)
  expect_equal(unname(sp$grad[slow, c("k1", "k2", "kf")]), rbind(c(-0.5, -0.5, 0), c(-0.5, -0.5, 0)),
               tolerance = 1e-10)
  # the term on a ramp that holds both, smooth in the rates
  w <- 10
  stab <- constraintStability(m$chain, weight = w, margin = 0.2)
  r <- 1 - 1.3 / (0.2 * 10)
  out <- stab(pars)
  expect_equal(out$value, 2 * w * r^3, tolerance = 1e-8)
  expect_equal(out$gradient[c("k1", "k2", "kf")], c(k1 = -3 * w * r^2 / 2, k2 = -3 * w * r^2 / 2,
                                                    kf = 2 * 3 * w * r^2 * 1.3 / (0.2 * 100)), tolerance = 1e-8)
  # with feedback through the fast mode the pair splits; the gradient is that of the value
  for (q in list(pars, replace(pars, c("c", "d"), c(0.5, 0.3)))) {
    num <- numDeriv::grad(function(v) stab(setNames(v, names(q)))$value, unname(q))
    expect_equal(unname(stab(q)$gradient[names(q)]), num, tolerance = 1e-7)
  }
})

test_that("the rate scale has the exact gradient, also where a volume moves the manifold", {
  skip_if_no_compile()
  pv <- pp_models()$vol
  tn <- grep("^total", attr(pv, "parameters"), value = TRUE)
  pars <- c(k1 = 1.5, k2 = 0.5, V = 0.4, setNames(2, tn))
  spec <- attr(attr(pv, "mappings")[[1]], "spectrum")
  sp <- spec(pars, onset = -Inf)
  # on a one-dimensional manifold the scale is the modulus of the eigenvalue
  expect_equal(sp$scale, abs(sp$re), tolerance = 1e-10)
  expect_equal(sp$scaleGrad[names(pars)], -sp$grad[1, names(pars)], tolerance = 1e-10)
  num <- numDeriv::grad(function(q) spec(setNames(q, names(pars)), onset = Inf)$scale, unname(pars))
  expect_equal(unname(sp$scaleGrad[names(pars)]), num, tolerance = 1e-7)
})

# The term of one eigenvalue from the real part `re` and the rate scale `s`, with both gradients,
# for the ramp of width `margin * s` below `tol`.
.stabilityExact <- function(re, dre, s, ds, w, tol = 0, margin = 1e-6) {
  span <- margin * s
  r <- max(0, (re - tol + span) / span)
  dr <- (dre + (1 - r) * margin * ds) / span
  rho <- sqrt(w) * if (r <= 1) r^1.5 else 1 + 1.5 * (r - 1)
  drho <- sqrt(w) * 1.5 * sqrt(min(r, 1)) * dr
  list(value = rho^2, gradient = 2 * rho * drho, hessian = 2 * tcrossprod(drho), r = r)
}

test_that("constraintStability is zero on stable steady states and rises at the Hopf boundary", {
  skip_if_no_compile()
  m <- pp_models()
  w <- 10
  stab <- constraintStability(m$bruss, m$bruss_p, weight = w, margin = 1e-6)
  expect_setequal(attr(stab, "parameters"), c("la", "lb"))
  # the term of the complex pair from the closed-form eigenvalues and rate scale
  exact <- function(la, lb) {
    a <- exp(la); b <- exp(lb)
    one <- .stabilityExact((b - 1 - a^2) / 2, c(-a^2, b / 2), a, c(a, 0), w)
    list(value = 2 * one$value, gradient = 2 * one$gradient, hessian = 2 * one$hessian, r = one$r)
  }
  # zero on stable steady states, also with slow modes
  expect_equal(stab(c(la = 0, lb = log(1.999)))$value, 0)
  expect_equal(stab(c(la = log(0.05), lb = log(1e-3)))$value, 0)
  expect_equal(nrow(attr(stab(c(la = 0, lb = log(1.999))), "eigenvalues")), 0L)
  # below the onset, on the ramp, at the tolerance and beyond; on the ramp to the accuracy of
  # the root relative to the width of the ramp
  for (re in c(-2e-6, -1e-6, -5e-7, -1e-9, 5e-7, 1e-6, 0.25)) {
    lb <- log(2 + 2 * re)
    out <- stab(c(la = 0, lb = lb))
    ref <- exact(0, lb)
    expect_equal(out$value, ref$value, tolerance = 1e-3)
    expect_equal(unname(out$gradient), ref$gradient, tolerance = 1e-3)
    expect_equal(unname(out$hessian), ref$hessian, tolerance = 1e-3)
  }
  # each eigenvalue at the tolerance adds the weight, beyond it the term grows quadratically
  span <- 1e-6
  expect_equal(stab(c(la = 0, lb = log(2)))$value, 2 * w, tolerance = 1e-3)
  expect_equal(stab(c(la = 0, lb = log(2 + 2 * span)))$value, 2 * w * 2.5^2, tolerance = 1e-3)
  pars <- c(la = 0, lb = log(2.5))
  num <- numDeriv::grad(function(q) stab(setNames(q, names(pars)))$value, unname(pars))
  expect_equal(unname(stab(pars)$gradient), num, tolerance = 1e-6)
  # value, gradient and Gauss-Newton curvature go to 0 at the onset as r^3, r^2 and r
  near <- vapply(c(1e-3, 1e-2), function(r) {
    o <- stab(c(la = 0, lb = log(2 + 2 * (r - 1) * span)))
    c(o$value, sqrt(sum(o$gradient^2)), max(abs(o$hessian)))
  }, numeric(3))
  expect_equal(near[, 1] / near[, 2], c(1e-3, 1e-2, 1e-1), tolerance = 1e-2)
  # and are continuous at the tolerance
  side <- lapply(c(1 - 1e-7, 1 + 1e-7), function(r) stab(c(la = 0, lb = log(2 + 2 * (r - 1) * span))))
  expect_equal(side[[1]]$value, side[[2]]$value, tolerance = 1e-5)
  expect_equal(side[[1]]$gradient, side[[2]]$gradient, tolerance = 1e-5)
  expect_equal(side[[1]]$hessian, side[[2]]$hessian, tolerance = 1e-5)
  # the eigenvalues on the ramp and beyond, with the conditions and the leading states
  ev <- attr(stab(c(la = 0, lb = log(2 + 1e-6))), "eigenvalues")
  expect_equal(nrow(ev), 2L)
  expect_equal(ev$re, rep(5e-7, 2), tolerance = 1e-3)
  expect_equal(sort(ev$im), c(-1, 1), tolerance = 1e-6)
  expect_equal(sum(ev$penalty), exact(0, log(2 + 1e-6))$value, tolerance = 1e-3)
  expect_match(ev$states, "^y \\(1\\.00\\), x \\(0\\.71\\)$")
  # the derivatives only where asked for, as the terms it is summed with
  expect_null(stab(pars, deriv = FALSE)$gradient)
  expect_null(stab(pars, hessian = FALSE)$hessian)
  expect_s4_class(stab(pars, hessian = "sparse")$hessian, "sparseMatrix")
  # beyond the boundary the unstable steady state, on which the term rises
  expect_equal(as.numeric(m$bruss(c(a = 1, b = 2.5))[[1]][c("x", "y")]), c(1, 2.5), tolerance = 1e-8)
  # without p on the inputs of the steady state, the same value
  expect_equal(constraintStability(m$bruss, weight = w, margin = 1e-6)(c(a = 1, b = 2 + 1e-6))$value,
               stab(c(la = 0, lb = log(2 + 1e-6)))$value, tolerance = 1e-3)
  expect_error(constraintStability(m$bruss, m$bruss_p, weight = 0), "weight")
  expect_error(constraintStability(m$bruss, m$bruss_p, margin = 0), "margin")
  expect_error(constraintStability(m$bruss, m$bruss_p, width = -1), "width")
  # an absolute width beyond margin * s: a ramp that does not move with the scale
  wide <- constraintStability(m$bruss, m$bruss_p, weight = w, width = 0.05)
  for (re in c(-0.06, -0.03, -0.01, 0.02)) {
    pars <- c(la = 0, lb = log(2 + 2 * re))
    ref <- .stabilityExact(re, c(-1, 1 + re), 1, c(0, 0), w, margin = 0.05)
    out <- wide(pars)
    expect_equal(out$value, 2 * ref$value, tolerance = 1e-6)
    expect_equal(unname(out$gradient), 2 * ref$gradient, tolerance = 1e-6)
    expect_equal(unname(out$hessian), 2 * ref$hessian, tolerance = 1e-6)
  }
})

test_that("constraintStability follows a real eigenvalue through 0 smoothly", {
  skip_if_no_compile()
  m <- pp_models()
  w <- 10
  stab <- constraintStability(m$trans, weight = w, margin = 1e-6)
  # a start that leads to the root whose eigenvalue crosses 0
  pars <- function(re) c(q = -re, x = -sign(re), y = 1)
  at <- function(re) {
    o <- stab(pars(re))
    expect_lt(abs(as.numeric(m$trans(pars(re))[[1]]["x"]) + re), 1e-10)
    c(value = o$value, gradient = o$gradient[["q"]], hessian = o$hessian["q", "q"])
  }
  # the term from the closed-form eigenvalue and rate scale
  exact <- function(re) {
    e <- .stabilityExact(re, -1, 1, 0, w)
    c(value = e$value, gradient = e$gradient, hessian = e$hessian[1, 1])
  }
  for (re in c(-1e-3, -2e-6, -1e-6, -5e-7, -1e-9, 1e-7, 1e-6, 0.1, 0.5))
    expect_equal(at(re), exact(re), tolerance = 1e-3)
  span <- 1e-6
  inner <- at((1e-3 - 1) * span) / at((1e-2 - 1) * span)
  expect_equal(unname(inner), c(1e-3, 1e-2, 1e-1), tolerance = 0.05)
  # beyond the tolerance the term is smooth enough for a numerical derivative
  for (re in c(0.1, 0.5)) {
    g <- stab(pars(re))$gradient[["q"]]
    expect_equal(g, numDeriv::grad(function(z) stab(c(q = z, x = -sign(re), y = 1))$value, -re),
                 tolerance = 1e-6)
  }
  # margin sets the width of the ramp relative to the scale
  wide <- constraintStability(m$trans, weight = w, margin = 0.1)
  for (re in c(-0.2, -0.05, 0.05)) {
    o <- wide(pars(re))
    expect_equal(o$value, .stabilityExact(re, -1, 1, 0, w, margin = 0.1)$value, tolerance = 1e-8)
    expect_equal(o$gradient[["q"]], numDeriv::grad(function(z) wide(c(q = z, x = -sign(re), y = 1))$value, -re),
                 tolerance = 1e-6)
  }
  old <- controls(m$trans, name = "controlsPTC")
  on.exit(controls(m$trans, name = "controlsPTC") <- old, add = TRUE)
  # an absolute stabilityTol moves the ramp
  controls(m$trans, name = "controlsPTC") <- list(stabilityTol = 0.5, positive = FALSE)
  for (re in c(0.2, 0.5 - 5e-7, 0.7)) {
    o <- stab(pars(re))
    expect_equal(o$value, .stabilityExact(re, -1, 1, 0, w, tol = 0.5)$value, tolerance = 1e-3)
  }
  expect_equal(stab(pars(0.7))$gradient[["q"]],
               numDeriv::grad(function(z) stab(c(q = z, x = -1, y = 1))$value, -0.7), tolerance = 1e-6)
})

test_that("constraintStability counts a steady state once over conditions that differ elsewhere", {
  skip_if_no_compile()
  m <- pp_models()
  pars <- c(la = 0, lb = log(2 - 1e-6))
  one <- constraintStability(m$bruss, m$bruss_p, weight = 10, margin = 1e-6)(pars)
  two <- constraintStability(m$bruss_s, m$bruss_c, weight = 10, margin = 1e-6)(pars)
  expect_gt(one$value, 0)
  expect_equal(two$value, one$value, tolerance = 1e-3)
  expect_equal(two$gradient, one$gradient, tolerance = 1e-3)
  expect_equal(attr(two, "eigenvalues")$condition, c("C1, C2", "C1, C2"))
})

test_that("Pimpl admits no unstable root beyond the precision of the eigenvalues", {
  skip_if_no_compile()
  m <- pp_models()
  root <- m$bruss_s(c(a = 1, b = 2 + 1e-9, s = 1))[[1]]
  expect_equal(as.numeric(root[c("x", "y")]), c(1, 2 + 1e-9), tolerance = 1e-8)
  expect_error(m$bruss_s(c(a = 1, b = 2 + 1e-6, s = 1)), "no stable root")
})
