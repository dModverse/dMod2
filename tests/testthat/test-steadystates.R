## Context: "SteadyStates"  (context() is deprecated in testthat 3e; kept as a note)
test_that("steady_states_are_steady", {

  # Run inside tempdir so steadyStates() (writes reactions_for_Alyssa*) and
  # the Pexpl + odemodel codegen don't pollute tests/testthat/.
  withr::local_dir(tempdir())

  #-!Start example code
  #-! library(dMod2)
  #-! setwd(tempdir())

  reactions <- eqnlist()
  reactions <-   addReaction(reactions, "Tca_buffer", "Tca_cyto", "import_Tca*Tca_buffer", "Basolateral uptake")
  reactions <-   addReaction(reactions, "Tca_cyto", "Tca_buffer", "export_Tca_baso*Tca_cyto", "Basolateral efflux")
  reactions <-   addReaction(reactions, "Tca_cyto", "Tca_canalicular", "export_Tca_cana*Tca_cyto", "Canalicular efflux")
  reactions <-   addReaction(reactions, "Tca_canalicular", "Tca_buffer", "transport_Tca*Tca_canalicular", "Transport bile")

  mysteadies <- steadyStates(reactions)
  #-! print(mysteadies)

  x <- Xs(odemodel(reactions, modelname = "ssTest", compile = FALSE))
  compile(x)

  parameters <- getParameters(x)
  trafo <- `names<-`(parameters, parameters)
  trafo <- repar("inner~steadyEqn", trafo, inner = names(mysteadies), steadyEqn = mysteadies)

  pSS <- P(trafo, condition = "steady")

  set.seed(2)
  pars <- structure(runif( length(getParameters(pSS)), 0,1), names = getParameters(pSS))

  prediction <- (x*pSS)(seq(0,10, 0.1), pars, deriv = F)
  #-! plot(prediction)
  #-!End example code

  is_steady <- function(prediction) {
    mypred <- wide2long(prediction)
    steady_conds <- sapply(split(mypred, mypred[c("name", "condition")]), function(i) {
      sd(i$value) <= 1e-8 * max(abs(i$value), 1)
    })
    return(all(steady_conds))
  }

  # Define your expectations here
  expect_length(mysteadies, 3)
  expect_true(is_steady(prediction))

})

test_that("additive fluxes are split instead of derailing the solver", {

  # A production with a basal and an induced arm has no rate constant of its
  # own, which used to leave fluxpars one entry short of the flux columns --
  # every later fluxpars[column] lookup shifted.
  withr::local_dir(tempdir())

  reactions <- eqnlist()
  reactions <- addReaction(reactions, "", "B", "k_pr_B", "B production")
  reactions <- addReaction(reactions, "B", "", "k_dg_B*B", "B degradation")
  reactions <- addReaction(reactions, "", "A", "(k_basal_A + k_ind_A*B)", "A production")
  reactions <- addReaction(reactions, "A", "", "k_dg_A*A", "A degradation")

  mysteadies <- steadyStates(reactions)

  expect_setequal(names(mysteadies), c("A", "B"))

  # B = k_pr_B/k_dg_B and A = (k_basal_A + k_ind_A*B)/k_dg_A
  pars <- c(k_pr_B = 0.7, k_dg_B = 0.3, k_basal_A = 0.2, k_ind_A = 1.3, k_dg_A = 0.9)
  ss <- lapply(mysteadies, function(eqn) eval(parse(text = eqn), as.list(pars)))
  B_expected <- unname(pars["k_pr_B"] / pars["k_dg_B"])
  A_expected <- unname((pars["k_basal_A"] + pars["k_ind_A"] * B_expected) / pars["k_dg_A"])
  expect_equal(ss$B, B_expected)
  expect_equal(ss$A, A_expected)

})

test_that("a locked rate constant never costs the solution its positivity", {

  # Solving R, K and LRK up front locks k_on and k_on2, which leaves LR with no
  # usable pivot side. Dropping a locked constant from the side instead of
  # disqualifying it made LR's two-flux side look safe, and the solve returned
  # influx - other_outflux -- a rate that can go negative.
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
  free <- setdiff(symbolsOf(c(as.character(odes), mysteadies)), names(mysteadies))
  env <- as.list(stats::runif(length(free), 0.2, 2))
  names(env) <- free

  # Entries kept as their own name are free parameters, not definitions.
  todo <- mysteadies[trimws(mysteadies) != names(mysteadies)]
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

  # getFluxes() scales every flux by V_ref / V_X, but the backend sees only
  # rates and stoichiometry. Without folding the ratios in, steadyStates()
  # silently returned the fixed point of the unscaled system.
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

test_that("positive = TRUE spends the row on a rate constant, not on a difference", {

  # Rec's ODE is not affine in Rec, so the direct positive-solve pass skips it.
  # By the time this phase reaches it Lig is substituted and the Lig-driven sink
  # has become the constant secretion flux, so solving for Rec would give
  # (k_pr_Rec - k_pr_Src)/k_dg_Rec. Linear in Rec, so there is no second root to
  # pick: only the pivot onto k_pr_Rec keeps the result a sum.
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

  # The invariant: positive = TRUE (the default) must not emit a subtraction.
  expect_false(any(grepl("-", mysteadies, fixed = TRUE)))

  # The row was spent on the rate constant, so Rec stays a free parameter.
  expect_true("k_pr_Rec" %in% names(mysteadies))
  expect_false("Rec" %in% names(mysteadies))

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

test_that("a sign-indefinite fixed point is refused (1.2) or repaired (1.3)", {

  # With a basal arm the additive flux is split, so no single rate constant
  # can absorb the row. Version 1.2 substitutes solutions eagerly: by the
  # time Rec is reached, the Lig-driven sink has become a constant, every
  # candidate pivot is a difference, and the backend must refuse rather than
  # return a negative-capable trafo. Version 1.3 keeps the fluxes symbolic,
  # so the two production arms are still available as a ratio-parameter
  # pivot and the model is solved manifestly positively.
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

  # The diagnosis goes to Python stdout, which R cannot capture, so only the
  # return is asserted.
  expect_identical(steadyStates(reactions, version = "1.2"), 0)

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

test_that("a recycled pivot rate constant is resolved, not defined by itself", {

  # Cpx_int is removed first and spends its row on its two out-fluxes, so the
  # recycling column in Cpx's row is replaced by the influx it was solved
  # against -- k_int*Cpx. Cpx's own influx sum then contains k_int, which is
  # also the first pivot Cpx spends its row on, so distributing that sum used
  # to emit k_int = f(..., k_int): a trafo defined in terms of itself. The sum
  # is linear in k_int and has to be resolved against the pivot relations
  # first.
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

  # The row was spent on the recycling loop's rate constants.
  expect_true(all(c("k_int", "k_rec") %in% names(mysteadies)))

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
