test_that("a closed cycle started in its steady state stays there", {
  withr::local_dir(tempdir())

  cycle <- eqnlist() |>
    addReaction("U", "V", "a*U") |>
    addReaction("V", "U", "b*V") |>
    addReaction("V", "W", "c*V") |>
    addReaction("W", "U", "d*W")
  ss <- steadyStates(cycle)
  expect_false(identical(ss, 0))

  x <- Xs(odemodel(cycle, modelname = "ssCycle", compile = FALSE))
  inner <- getParameters(x)
  p <- P(repar(setNames(inner, inner), "lhs~rhs", lhs = names(ss), rhs = ss),
         condition = "steady", compile = FALSE)
  compile(x, p, output = "ssCycle_all", cores = test_cores())

  outer <- setNames(seq(0.3, 0.9, length.out = length(getParameters(p))), getParameters(p))
  traj <- wide2long((x * p)(seq(0, 10, 0.5), outer, deriv = FALSE))
  drift <- tapply(traj$value, traj$name, function(v) diff(range(v)) / max(abs(v), 1))
  expect_true(all(drift <= 1e-8))
})

test_that("additive fluxes are split instead of derailing the solver", {

  # A production with a basal and an induced arm has no rate constant of its
  # own; fluxpars must still align with the flux columns.
  withr::local_dir(tempdir())

  reactions <- eqnlist()
  reactions <- addReaction(reactions, "", "B", "k_pr_B", "B production")
  reactions <- addReaction(reactions, "B", "", "k_dg_B*B", "B degradation")
  reactions <- addReaction(reactions, "", "A", "(k_basal_A + k_ind_A*B)", "A production")
  reactions <- addReaction(reactions, "A", "", "k_dg_A*A", "A degradation")

  mysteadies <- steadyStates(reactions)

  expect_setequal(names(mysteadies), c("A", "B"))

  # both states match the closed form
  pars <- c(k_pr_B = 0.7, k_dg_B = 0.3, k_basal_A = 0.2, k_ind_A = 1.3, k_dg_A = 0.9)
  ss <- lapply(mysteadies, function(eqn) eval(parse(text = eqn), as.list(pars)))
  B_expected <- unname(pars["k_pr_B"] / pars["k_dg_B"])
  A_expected <- unname((pars["k_basal_A"] + pars["k_ind_A"] * B_expected) / pars["k_dg_A"])
  expect_equal(ss$B, B_expected)
  expect_equal(ss$A, A_expected)

})

test_that("a binding cascade is solved as positive sums", {

  withr::local_dir(tempdir())

  reactions <- eqnlist()
  reactions <- addReaction(reactions, "", "R", "k_pr_R", "R production")
  reactions <- addReaction(reactions, "R", "", "k_dg_R*R", "R degradation")
  reactions <- addReaction(reactions, "", "L", "k_pr_L", "L production")
  reactions <- addReaction(reactions, "L", "", "k_dg_L*L", "L degradation")
  reactions <- addReaction(reactions, "", "K", "k_pr_K", "K production")
  reactions <- addReaction(reactions, "K", "", "k_dg_K*K", "K degradation")
  reactions <- addReaction(reactions, "L + R", "LR", "k_on*L*R", "L binds R")
  reactions <- addReaction(reactions, "LR", "", "k_dg_LR*LR", "LR degradation")
  reactions <- addReaction(reactions, "LR + K", "LRK", "k_on2*LR*K", "K binds LR")
  reactions <- addReaction(reactions, "LRK", "", "k_dg_LRK*LRK", "LRK degradation")

  mysteadies <- steadyStates(reactions)

  # Without solveQuadratic the output is manifestly positive: no subtractions.
  expect_false(any(grepl("-", mysteadies, fixed = TRUE)))

  # And it is a steady state, at a random positive point.
  odes <- as.eqnvec(reactions)
  set.seed(4)
  symbolsOf <- function(x) unique(unlist(lapply(x, function(e) all.vars(parse(text = e)))))
  # Entries kept as their own name are free parameters, not definitions.
  defs <- mysteadies[trimws(mysteadies) != names(mysteadies)]
  free <- setdiff(symbolsOf(c(as.character(odes), mysteadies)), names(defs))
  env <- as.list(stats::runif(length(free), 0.2, 2))
  names(env) <- free

  todo <- defs
  while (length(todo)) {
    ready <- vapply(todo, function(eqn) all(all.vars(parse(text = eqn)) %in% names(env)),
                    logical(1))
    if (!any(ready)) break
    env[names(todo)[ready]] <- lapply(todo[ready], function(eqn) eval(parse(text = eqn), env))
    todo <- todo[!ready]
  }
  expect_length(todo, 0)

  expect_true(all(unlist(env[names(mysteadies)]) > 0))
  residual <- vapply(odes, function(e) eval(parse(text = e), env), numeric(1))
  expect_lt(max(abs(residual)), 1e-10)

})

test_that("compartment volume ratios reach the steady-state backend", {

  # getFluxes() scales fluxes by volume ratios the backend does not see, so the
  # ratios must be folded in before solving.
  withr::local_dir(tempdir())

  reactions <- eqnlist()
  reactions <- addReaction(reactions, "", "L", "k_pr_L", "L influx", compartment = "ext")
  reactions <- addReaction(reactions, "L", "", "k_dg_L*L", "L decay", compartment = "ext")
  reactions <- addReaction(reactions, "", "R", "k_pr_R", "R synthesis", compartment = "cell")
  reactions <- addReaction(reactions, "R", "", "k_dg_R*R", "R decay", compartment = "cell")
  reactions <- addReaction(reactions, "L + R", "LR", "k_on*L*R", "binding",
                           compartment = "cell", rateCompartment = "cell")
  reactions <- addReaction(reactions, "LR", "", "k_dg_LR*LR", "LR decay", compartment = "cell")
  reactions$compartments$cell$volume <- "Vc"
  reactions$compartments$ext$volume  <- "Ve"

  mysteadies <- steadyStates(reactions)

  # The unscaled system mentions neither volume.
  expect_true(any(grepl("Vc", mysteadies, fixed = TRUE)))
  expect_true(any(grepl("Ve", mysteadies, fixed = TRUE)))

  # The one that matters: the fixed point solves the ODE dMod integrates.
  odes <- as.eqnvec(reactions)
  symbolsOf <- function(x) unique(unlist(lapply(x, function(e) all.vars(parse(text = e)))))
  defs <- mysteadies[trimws(mysteadies) != names(mysteadies)]
  set.seed(3)
  free <- setdiff(symbolsOf(c(as.character(odes), mysteadies)), names(defs))
  env <- as.list(stats::runif(length(free), 0.3, 2))
  names(env) <- free
  for (nm in names(defs)) env[[nm]] <- eval(parse(text = defs[[nm]]), env)

  residual <- vapply(odes, function(e) eval(parse(text = e), env), numeric(1))
  expect_lt(max(abs(residual)), 1e-10)

})

test_that("a leaking cluster across compartments is 0, a moiety listed first is not", {

  # A binding across compartments leaks as a cluster only in amounts. A conserved
  # moiety listed before it stays nonzero; only states draining into a leaking
  # reaction are 0.
  withr::local_dir(tempdir())
  toy <- function(Vc, Ve) eqnlist() |>
    assignCompartment(L = "ext", volume = Ve) |>
    addReaction("A", "Ap", "k_act*R*A", compartment = "cell") |>
    addReaction("Ap", "A", "k_deact*Ap", compartment = "cell") |>
    addReaction("", "R", "k_pr_R", compartment = "cell") |>
    addReaction("R", "", "k_dg_R*R", compartment = "cell") |>
    addReaction("L + R", "LR", "k_on*L*R", compartment = "cell", rateCompartment = "cell") |>
    addReaction("LR", "L + R", "k_off*LR", compartment = "cell", rateCompartment = "cell") |>
    addReaction("LR", "", "k_dg_LR*LR", compartment = "cell") |>
    setCompartmentVolume(cell = Vc)
  symbolsOf <- function(x) unique(unlist(lapply(x, function(e) all.vars(parse(text = e)))))

  for (vols in list(c("Vc", "Ve"), c("1.4", "0.45"))) {
    lab <- paste0("volumes ", vols[1], "/", vols[2])
    r <- toy(vols[1], vols[2])
    ss <- steadyStates(r, verbose = FALSE)
    expect_false(identical(ss, 0), label = lab)
    if (identical(ss, 0)) next
    expect_setequal(names(ss)[trimws(ss) == "0"], c("L", "LR"))

    odes <- as.eqnvec(r)
    defs <- ss[trimws(ss) != names(ss)]
    set.seed(7)
    free <- setdiff(symbolsOf(c(as.character(odes), defs)), names(defs))
    env <- as.list(setNames(stats::runif(length(free), 0.3, 2), free))
    for (n in names(defs)) env[[n]] <- eval(parse(text = defs[[n]]), env)
    expect_true(all(unlist(env[c("A", "Ap", "R")]) > 0), label = lab)
    residual <- vapply(odes, function(e) eval(parse(text = e), env), numeric(1))
    expect_lt(max(abs(residual)), 1e-10, label = lab)
  }

})

test_that("a saturating binding rate is solved without a difference", {

  withr::local_dir(tempdir())

  reactions <- eqnlist()
  reactions <- addReaction(reactions, "", "Src", "k_pr_Src", "Source production")
  reactions <- addReaction(reactions, "Src", "Lig", "k_sec*Src", "Ligand secretion")
  reactions <- addReaction(reactions, "", "Rec", "k_pr_Rec", "Receptor synthesis")
  reactions <- addReaction(reactions, "Rec", "", "k_dg_Rec*Rec", "Receptor degradation")
  reactions <- addReaction(reactions, "Lig + Rec", "Cpx", "k_form*Lig*Rec/(Km + Rec)",
                           "Saturating complex formation")
  reactions <- addReaction(reactions, "Cpx", "", "k_dg_Cpx*Cpx", "Complex degradation")

  mysteadies <- steadyStates(reactions)

  # No result contains a subtraction.
  expect_false(any(grepl("-", mysteadies, fixed = TRUE)))


  # And it is still a steady state, at a random positive point.
  odes <- as.eqnvec(reactions)
  symbolsOf <- function(x) unique(unlist(lapply(x, function(e) all.vars(parse(text = e)))))
  defs <- mysteadies[trimws(mysteadies) != names(mysteadies)]
  set.seed(11)
  free <- setdiff(symbolsOf(c(as.character(odes), mysteadies)), names(defs))
  env <- as.list(stats::runif(length(free), 0.3, 2))
  names(env) <- free

  todo <- defs
  while (length(todo)) {
    ready <- vapply(todo, function(eqn) all(all.vars(parse(text = eqn)) %in% names(env)),
                    logical(1))
    if (!any(ready)) break
    env[names(todo)[ready]] <- lapply(todo[ready], function(eqn) eval(parse(text = eqn), env))
    todo <- todo[!ready]
  }
  expect_length(todo, 0)

  expect_true(all(unlist(env[names(defs)]) > 0))
  residual <- vapply(odes, function(e) eval(parse(text = e), env), numeric(1))
  expect_lt(max(abs(residual)), 1e-10)

})

test_that("a production with a basal arm is solved as positive sums", {

  # The additive flux is split, and no single rate constant absorbs the balance.
  withr::local_dir(tempdir())

  reactions <- eqnlist()
  reactions <- addReaction(reactions, "", "Src", "k_pr_Src", "Source production")
  reactions <- addReaction(reactions, "Src", "Lig", "k_sec*Src", "Ligand secretion")
  reactions <- addReaction(reactions, "", "Rec", "(k_pr_Rec + k_basal_Rec)",
                           "Receptor synthesis, basal plus induced")
  reactions <- addReaction(reactions, "Rec", "", "k_dg_Rec*Rec", "Receptor degradation")
  reactions <- addReaction(reactions, "Lig + Rec", "Cpx", "k_form*Lig*Rec/(Km + Rec)",
                           "Saturating complex formation")
  reactions <- addReaction(reactions, "Cpx", "", "k_dg_Cpx*Cpx", "Complex degradation")

  mysteadies <- steadyStates(reactions)
  expect_false(any(grepl("-", mysteadies, fixed = TRUE)))
  expect_true("k_pr_Rec" %in% names(mysteadies))

  # And it is a steady state, at a random positive point (r_Rec_1 included).
  odes <- as.eqnvec(reactions)
  symbolsOf <- function(x) unique(unlist(lapply(x, function(e) all.vars(parse(text = e)))))
  defs <- mysteadies[trimws(mysteadies) != names(mysteadies)]
  set.seed(12)
  free <- setdiff(symbolsOf(c(as.character(odes), mysteadies)), names(defs))
  env <- as.list(stats::runif(length(free), 0.3, 2))
  names(env) <- free
  for (nm in names(defs)) env[[nm]] <- eval(parse(text = defs[[nm]]), env)

  residual <- vapply(odes, function(e) eval(parse(text = e), env), numeric(1))
  expect_lt(max(abs(residual)), 1e-10)

})

test_that("a recycling loop is solved without a definition referring to itself", {

  withr::local_dir(tempdir())

  reactions <- eqnlist()
  reactions <- addReaction(reactions, "", "Lig", "k_pr_Lig", "Ligand production")
  reactions <- addReaction(reactions, "Lig", "", "k_dg_Lig*Lig", "Ligand clearance")
  reactions <- addReaction(reactions, "Lig", "Cpx", "k_form*Lig", "Complex formation")
  reactions <- addReaction(reactions, "Cpx", "Cpx_int", "k_int*Cpx", "Internalisation")
  reactions <- addReaction(reactions, "Cpx_int", "Cpx", "k_rec*Cpx_int", "Recycling")
  reactions <- addReaction(reactions, "Cpx_int", "", "k_dg_int*Cpx_int", "Internalised degradation")
  reactions <- addReaction(reactions, "Cpx", "", "k_dg_Cpx*Cpx", "Complex degradation")
  reactions <- addReaction(reactions, "Cpx", "", "k_dg_Cpx_Inh*Cpx*Inh", "Inhibitor-driven removal")
  reactions <- addReaction(reactions, "", "Inh", "k_pr_Inh*Cpx", "Inhibitor induction")
  reactions <- addReaction(reactions, "Inh", "", "k_dg_Inh*Inh", "Inhibitor decay")

  mysteadies <- steadyStates(reactions)

  defs <- mysteadies[trimws(mysteadies) != names(mysteadies)]
  # No definition may reference its own left-hand side.
  selfref <- mapply(function(nm, eqn) nm %in% all.vars(parse(text = eqn)),
                    names(defs), defs)
  expect_false(any(selfref))
  expect_false(any(grepl("-", mysteadies, fixed = TRUE)))

  # And it is a steady state, at a random positive point.
  odes <- as.eqnvec(reactions)
  symbolsOf <- function(x) unique(unlist(lapply(x, function(e) all.vars(parse(text = e)))))
  set.seed(13)
  free <- setdiff(symbolsOf(c(as.character(odes), mysteadies)), names(defs))
  env <- as.list(stats::runif(length(free), 0.3, 2))
  names(env) <- free

  todo <- defs
  while (length(todo)) {
    ready <- vapply(todo, function(eqn) all(all.vars(parse(text = eqn)) %in% names(env)),
                    logical(1))
    if (!any(ready)) break
    env[names(todo)[ready]] <- lapply(todo[ready], function(eqn) eval(parse(text = eqn), env))
    todo <- todo[!ready]
  }
  expect_length(todo, 0)

  expect_true(all(unlist(env[names(defs)]) > 0))
  residual <- vapply(odes, function(e) eval(parse(text = e), env), numeric(1))
  expect_lt(max(abs(residual)), 1e-10)

})

test_that("a symbol defined twice or in a cycle is an error", {
  reticulate::py_require(c("numpy", "sympy", "scipy"))
  ap <- reticulate::import_from_path("AlyssaPetit",
                                     path = system.file("code", package = "dMod2"))
  topo <- ap$`_topo_sort_eqs`
  expect_equal(unlist(topo(list("a = b", "b = c", "c = c"))), c("b = c", "a = b", "c = c"))
  expect_error(topo(list("a = b", "a = c")), "more than once")
  expect_error(topo(list("a = 2*a")), "itself")
  expect_error(topo(list("a = b + 1", "b = a")), "cyclically")
})

test_that("a block that no single balance resolves is solved jointly", {

  # Receptors recycling through endosomes and a ligand-free complex: every
  # balance alone leaves a difference; the balances together cancel it.
  withr::local_dir(tempdir())
  r <- eqnlist() |>
    addReaction("", "R1", "k_transl_R1") |>
    addReaction("", "R2", "k_transl_R2") |>
    addReaction("R1", "R1_endo", "k_int_R1*R1") |>
    addReaction("R1_endo", "R1", "k_rec_R1_endo*R1_endo") |>
    addReaction("R1_endo", "", "k_dg_R1_endo*R1_endo") |>
    addReaction("R2", "R2_endo", "k_int_R2*R2") |>
    addReaction("R2_endo", "R2", "k_rec_R2_endo*R2_endo") |>
    addReaction("R2_endo", "", "k_dg_R2_endo*R2_endo") |>
    addReaction("R1 + R2", "R1_R2", "k_form_R1_R2*R1*R2") |>
    addReaction("R1_R2", "R1 + R2", "k_dec_R1_R2*R1_R2") |>
    addReaction("R1_R2", "R1_R2_int", "k_int_R1_R2*R1_R2") |>
    addReaction("R1_R2_int", "R1_R2", "k_rec_R1_R2_int*R1_R2_int") |>
    addReaction("R1_R2_int", "R1_endo + R2_endo", "k_dec_R1_R2_int*R1_R2_int") |>
    addReaction("R1", "", "k_dg_R1_S7*R1*S7/(Km_S7 + R1)")

  ss <- steadyStates(r)
  expect_false(identical(ss, 0))
  defs <- ss[trimws(ss) != names(ss)]
  expect_false(any(grepl("-", defs, fixed = TRUE)))

  odes <- as.eqnvec(r)
  symbolsOf <- function(x) unique(unlist(lapply(x, function(e) all.vars(parse(text = e)))))
  set.seed(7)
  free <- setdiff(symbolsOf(c(as.character(odes), defs)), names(defs))
  env <- as.list(setNames(stats::runif(length(free), 0.3, 2), free))
  todo <- defs
  while (length(todo)) {
    ready <- vapply(todo, function(e) all(all.vars(parse(text = e)) %in% names(env)), TRUE)
    if (!any(ready)) break
    env[names(todo)[ready]] <- lapply(todo[ready], function(e) eval(parse(text = e), env))
    todo <- todo[!ready]
  }
  expect_length(todo, 0)
  expect_true(all(unlist(env[names(defs)]) > 0))
  residual <- vapply(odes, function(e) eval(parse(text = e), env), numeric(1))
  expect_lt(max(abs(residual)), 1e-10)
})

test_that("every model above is solved positively", {

  withr::local_dir(tempdir())
  vol <- eqnlist() |>
    addReaction("", "L", "k_pr_L", compartment = "ext") |>
    addReaction("L", "", "k_dg_L*L", compartment = "ext") |>
    addReaction("", "R", "k_pr_R", compartment = "cell") |>
    addReaction("R", "", "k_dg_R*R", compartment = "cell") |>
    addReaction("L + R", "LR", "k_on*L*R", compartment = "cell", rateCompartment = "cell") |>
    addReaction("LR", "", "k_dg_LR*LR", compartment = "cell")
  vol$compartments$cell$volume <- "Vc"
  vol$compartments$ext$volume  <- "Ve"
  models <- list(
    transport = eqnlist() |>
      addReaction("Tca_buffer", "Tca_cyto", "import_Tca*Tca_buffer") |>
      addReaction("Tca_cyto", "Tca_buffer", "export_Tca_baso*Tca_cyto") |>
      addReaction("Tca_cyto", "Tca_canalicular", "export_Tca_cana*Tca_cyto") |>
      addReaction("Tca_canalicular", "Tca_buffer", "transport_Tca*Tca_canalicular"),
    additive = eqnlist() |>
      addReaction("", "B", "k_pr_B") |> addReaction("B", "", "k_dg_B*B") |>
      addReaction("", "A", "(k_basal_A + k_ind_A*B)") |> addReaction("A", "", "k_dg_A*A"),
    chain = eqnlist() |>
      addReaction("", "R", "k_pr_R") |> addReaction("R", "", "k_dg_R*R") |>
      addReaction("", "L", "k_pr_L") |> addReaction("L", "", "k_dg_L*L") |>
      addReaction("", "K", "k_pr_K") |> addReaction("K", "", "k_dg_K*K") |>
      addReaction("L + R", "LR", "k_on*L*R") |> addReaction("LR", "", "k_dg_LR*LR") |>
      addReaction("LR + K", "LRK", "k_on2*LR*K") |> addReaction("LRK", "", "k_dg_LRK*LRK"),
    volumes = vol,
    basal = eqnlist() |>
      addReaction("", "Src", "k_pr_Src") |> addReaction("Src", "Lig", "k_sec*Src") |>
      addReaction("", "Rec", "(k_pr_Rec + k_basal_Rec)") |>
      addReaction("Rec", "", "k_dg_Rec*Rec") |>
      addReaction("Lig + Rec", "Cpx", "k_form*Lig*Rec/(Km + Rec)") |>
      addReaction("Cpx", "", "k_dg_Cpx*Cpx"),
    recycling = eqnlist() |>
      addReaction("", "R1", "k_transl_R1") |> addReaction("", "R2", "k_transl_R2") |>
      addReaction("R1", "R1_endo", "k_int_R1*R1") |>
      addReaction("R1_endo", "R1", "k_rec_R1_endo*R1_endo") |>
      addReaction("R1_endo", "", "k_dg_R1_endo*R1_endo") |>
      addReaction("R2", "R2_endo", "k_int_R2*R2") |>
      addReaction("R2_endo", "R2", "k_rec_R2_endo*R2_endo") |>
      addReaction("R2_endo", "", "k_dg_R2_endo*R2_endo") |>
      addReaction("R1 + R2", "R1_R2", "k_form_R1_R2*R1*R2") |>
      addReaction("R1_R2", "R1 + R2", "k_dec_R1_R2*R1_R2") |>
      addReaction("R1_R2", "R1_R2_int", "k_int_R1_R2*R1_R2") |>
      addReaction("R1_R2_int", "R1_R2", "k_rec_R1_R2_int*R1_R2_int") |>
      addReaction("R1_R2_int", "R1_endo + R2_endo", "k_dec_R1_R2_int*R1_R2_int") |>
      addReaction("R1", "", "k_dg_R1_S7*R1*S7/(Km_S7 + R1)")
  )
  symbolsOf <- function(x) unique(unlist(lapply(x, function(e) all.vars(parse(text = e)))))

  for (nm in names(models)) {
    r <- models[[nm]]
    ss <- steadyStates(r, verbose = FALSE)
    expect_false(identical(ss, 0), label = nm)
    defs <- ss[trimws(ss) != names(ss)]
    expect_false(any(grepl("-", defs, fixed = TRUE)), label = nm)
    selfref <- mapply(function(n, e) n %in% all.vars(parse(text = e)), names(defs), defs)
    expect_false(any(selfref), label = nm)

    odes <- as.eqnvec(r)
    set.seed(5)
    free <- setdiff(symbolsOf(c(as.character(odes), defs)), names(defs))
    env <- as.list(setNames(stats::runif(length(free), 0.3, 2), free))
    for (n in names(defs)) env[[n]] <- eval(parse(text = defs[[n]]), env)
    expect_true(all(unlist(env[names(defs)]) > 0), label = nm)
    residual <- vapply(odes, function(e) eval(parse(text = e), env), numeric(1))
    expect_lt(max(abs(residual)), 1e-10, label = nm)
  }
})

test_that("solveQuadratic takes the positive root", {
  withr::local_dir(tempdir())
  models <- list(
    dimer = eqnlist() |> addReaction("", "A", "k_pr_A") |> addReaction("A", "", "k_dg_A*A") |>
      addReaction("2*A", "D", "k_dim*A^2") |> addReaction("D", "", "k_dg_D*D"),
    squareOnly = eqnlist() |> addReaction("", "A", "k_pr_A") |>
      addReaction("2*A", "D", "k_dim*A^2") |> addReaction("D", "", "k_dg_D*D"),
    saturating = eqnlist() |> addReaction("", "X", "k_pr_X") |> addReaction("X", "", "k_dg_X*X") |>
      addReaction("X", "", "k_cat*E*X/(Km + X)") |>
      addReaction("", "E", "k_pr_E") |> addReaction("E", "", "k_dg_E*E")
  )
  symbolsOf <- function(x) unique(unlist(lapply(x, function(e) all.vars(parse(text = e)))))
  for (nm in names(models)) {
    r <- models[[nm]]
    ss <- steadyStates(r, solveQuadratic = TRUE, verbose = FALSE)
    defs <- ss[trimws(ss) != names(ss)]
    # the state is a root, no rate constant was spent on its balance
    expect_true(any(grepl("sqrt", defs)), label = nm)
    odes <- as.eqnvec(r)
    set.seed(6)
    free <- setdiff(symbolsOf(c(as.character(odes), defs)), names(defs))
    env <- as.list(setNames(stats::runif(length(free), 0.3, 2), free))
    for (n in names(defs)) env[[n]] <- eval(parse(text = defs[[n]]), env)
    expect_true(all(unlist(env[names(defs)]) > 0), label = nm)
    residual <- vapply(odes, function(e) eval(parse(text = e), env), numeric(1))
    expect_lt(max(abs(residual)), 1e-10, label = nm)
  }
  # a positive linear coefficient keeps the root free of subtractions
  ss <- steadyStates(models$dimer, solveQuadratic = TRUE, verbose = FALSE)
  expect_false(any(grepl("-", ss, fixed = TRUE)))
})

test_that("the first state of a given conserved quantity stays free", {
  withr::local_dir(tempdir())
  r <- eqnlist() |> addReaction("", "K", "k_pr_K") |> addReaction("K", "", "k_dg_K*K") |>
    addReaction("A", "pA", "k_phos*K*A") |> addReaction("pA", "A", "k_dephos*pA")
  rates <- c("k_phos", "k_dephos", "k_pr_K", "k_dg_K")
  for (cq in c("A + pA = tA", "pA + A = tA")) {
    ss <- steadyStates(r, givenCQs = cq, neglect = rates, verbose = FALSE)
    kept <- sub(" .*", "", cq)
    expect_identical(ss[[kept]], kept)
    expect_false(identical(ss[[setdiff(c("A", "pA"), kept)]], setdiff(c("A", "pA"), kept)))
  }
})


test_that("species named like sympy objects (Ci, E, S, Q) are plain symbols", {
  withr::local_dir(tempdir())
  # Names of sympy objects must not resolve to those objects.
  r <- eqnlist() |>
    addReaction("", "R", "ksR", "synthesis") |>
    addReaction("R", "", "kdR*R", "turnover") |>
    addReaction("R + E", "C", "kon*R*E", "binding") |>
    addReaction("C", "Ci + E", "kint*C", "internalisation") |>
    addReaction("Ci", "", "kdeg*Ci", "degradation") |>
    addReaction("C + S", "C + Q", "kp*C*S", "phosphorylation") |>
    addReaction("Q", "S", "kdp*Q", "dephosphorylation")
  ss <- steadyStates(r, verbose = FALSE)
  expect_false(identical(ss, 0))
  expect_true(length(ss) > 0L)
})


test_that("steadyStates writes nothing without file and writes the RDS with it", {
  dir <- withr::local_tempdir()
  withr::local_dir(dir)
  el <- eqnlist() |>
    addReaction("", "A", "k_in") |>
    addReaction("A", "", "k_out * A")
  before <- list.files(tempdir(), pattern = "^steadyStates")

  ss <- steadyStates(el)
  expect_equal(list.files(dir), character(0))
  expect_setequal(list.files(tempdir(), pattern = "^steadyStates"), before)

  ss2 <- steadyStates(el, file = "ss.rds")
  expect_true(file.exists("ss.rds"))
  expect_identical(readRDS("ss.rds"), ss2)
  expect_identical(ss2, ss)
})
