# -------------------------------------------------------------------------#
# Symmetry motifs: small models on which symmetryDetection()/symmetryReduction()
# are slow or fail
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# Each motif is a function returning the arguments of symmetryDetection() and,
# where known, the number of non-identifiable directions (`nDir`). The motifs
# isolate one difficulty each: curved directions whose entries depend on many
# parameters, invariants outside the Darboux language, steady states before a
# stimulus, and parameters that act only after an event. run_motif.R runs one
# motif, run_all.sh the whole collection.
#
# [AUTHOR]
# Simon Beyer
#
# [Date]
# Sun 27 Sep 2026
# -------------------------------------------------------------------------#

.motifs <- list()

# ---- linear compartment models ---------------------------------------------------
# Input and output in compartment 1, leaks everywhere. Identifiable: the 2N-1
# coefficients of the transfer function plus the dose, so 3N-1 rates leave N curved
# directions whose entries mix every rate.

# catenary chain 1 - 2 - ... - N
.catenary <- function(N) {
  f <- character(N)
  for (i in seq_len(N)) {
    out <- c(paste0("k0", i), if (i < N) paste0("kf", i), if (i > 1) paste0("kb", i - 1))
    terms <- paste0("-(", paste(out, collapse = " + "), ")*x", i)
    if (i > 1) terms <- c(terms, paste0("kf", i - 1, "*x", i - 1))
    if (i < N) terms <- c(terms, paste0("kb", i, "*x", i + 1))
    f[i] <- paste(terms, collapse = " + ")
  }
  f <- setNames(f, paste0("x", seq_len(N)))
  ic <- setNames(c("D", rep("0", N - 1)), paste0("x", seq_len(N)))
  list(f = do.call(eqnvec, as.list(f)), g = eqnvec(y = "x1"),
       trafo = do.call(eqnvec, as.list(ic)), nDir = N)
}

# mammillary: central compartment 1 exchanging with 2..N
.mammillary <- function(N) {
  f <- character(N)
  out1 <- c("k01", paste0("k", 2:N, "1"))
  f[1] <- paste(c(paste0("-(", paste(out1, collapse = " + "), ")*x1"),
                  paste0("k1", 2:N, "*x", 2:N)), collapse = " + ")
  for (i in 2:N)
    f[i] <- paste0("k", i, "1*x1 - (k0", i, " + k1", i, ")*x", i)
  f <- setNames(f, paste0("x", seq_len(N)))
  ic <- setNames(c("D", rep("0", N - 1)), paste0("x", seq_len(N)))
  list(f = do.call(eqnvec, as.list(f)), g = eqnvec(y = "x1"),
       trafo = do.call(eqnvec, as.list(ic)), nDir = N)
}

.motifs$cat2  <- function() .catenary(2)
.motifs$cat3  <- function() .catenary(3)
.motifs$cat4  <- function() .catenary(4)
.motifs$cat5  <- function() .catenary(5)
.motifs$mam3  <- function() .mammillary(3)
.motifs$mam4  <- function() .mammillary(4)
.motifs$cat6  <- function() .catenary(6)

# ---- enzyme kinetics --------------------------------------------------------------
# full mass action, only the product seen through a scale; S0 known by dose
.motifs$enzyme <- function() {
  f <- eqnvec(E = "-k1*E*S + km1*C + k2*C", S = "-k1*E*S + km1*C",
              C = "k1*E*S - km1*C - k2*C", P = "k2*C")
  list(f = f, g = eqnvec(y = "s*P"),
       trafo = eqnvec(E = "E0", S = "S0", C = "0", P = "0"), fixed = "S0")
}

# ---- autocrine loop with a steady state before the stimulus -----------------------
# gene -> intracellular ligand pool -> secretion -> receptor complex -> gene.
# The pool is unobserved, its degradation competes with secretion.
.autocrine <- function() {
  f <- eqnvec(
    m = "ktx*(b0 + C) - dm*m",
    P = "ktl*m - (kdg + ksec)*P",
    L = "ksec*P/vol - kon*L*R + koff*C - kL*L",
    R = "ksR - kdR*R - kon*L*R + koff*C",
    C = "kon*L*R - koff*C - kint*C")
  g <- eqnvec(ym = "sm*m", yC = "sC*C", yL = "L")
  list(f = f, g = g)
}
.motifs$autocrine_free <- function() {
  a <- .autocrine()
  a$trafo <- eqnvec(m = "m0", P = "P0", L = "L0", R = "R0", C = "C0")
  a
}

# ---- rates acting only after an event ----------------------------------------------
# a knock-down switch enters at t = 10, after the transient of the first segment
.motifs$switch_late <- function() {
  f <- eqnvec(x = "k1*u - (d1 + ks*sw)*x", y2 = "k2*x - d2*y2", u = "0", sw = "0")
  ev <- addEvent(eventlist(), var = "sw", time = 10, value = "1", method = "replace")
  ev <- addEvent(ev, var = "u", time = 0, value = "dose", method = "replace")
  list(f = f, g = eqnvec(obs = "s*y2"), events = ev,
       conditions = data.frame(dose = c(1, 2), row.names = c("lo", "hi")))
}

# ---- nonlinear curved direction with a wide invariant -------------------------------
# two unobserved routes into one observed state; the routes differ only through
# products of their rates
.motifs$tworoute <- function() {
  f <- eqnvec(A = "-(ka + kb)*A", B1 = "ka*A - k1*B1", B2 = "kb*A - k2*B2",
              C = "k3*B1 + k4*B2 - kc*C")
  list(f = f, g = eqnvec(y = "C"), trafo = eqnvec(A = "A0", B1 = "0", B2 = "0", C = "0"),
       fixed = "A0")
}

# n unobserved routes: invariants k_i*r_i and the sum of the r_i. No face (a route at 0
# kills its product); the chart needs n-1 balances at once.
.routes <- function(n) {
  r <- paste0("r", seq_len(n)); k <- paste0("k", seq_len(n))
  f <- c(A = paste0("-(", paste(r, collapse = " + "), ")*A"),
         setNames(paste0(r, "*A - d", seq_len(n), "*B", seq_len(n)), paste0("B", seq_len(n))),
         C = paste0(paste0(k, "*B", seq_len(n), collapse = " + "), " - kc*C"))
  f <- do.call(eqnvec, as.list(f))
  tr <- do.call(eqnvec, as.list(c(A = "A0", setNames(rep("0", n), paste0("B", seq_len(n))),
                                  C = "0")))
  list(f = f, g = eqnvec(y = "C"), trafo = tr, fixed = "A0")
}
.motifs$route3 <- function() .routes(3)
.motifs$route4 <- function() .routes(4)

# ---- Hill-type production with an unobserved activator ------------------------------
.motifs$hill <- function() {
  f <- eqnvec(u = "-ku*u", x = "V*u^2/(K^2 + u^2) - d*x", y2 = "kt*x - d2*y2")
  list(f = f, g = eqnvec(obs = "s*y2"), trafo = eqnvec(u = "u0", x = "0", y2 = "0"))
}

# ---- autocrine loop from its steady state, stimulated by a ligand dose ---------------
# M011 in small: the fixed point as trafo (steadyStates()), a dose event at t = 0 in two
# conditions. The resting state ties the unobserved pool to every rate around it.
.autocrineReactions <- function() {
  eqnlist() |>
    addReaction("", "m", "ktx*(b0 + C)", "transcription") |>
    addReaction("m", "", "dm*m", "mRNA decay") |>
    addReaction("", "P", "ktl*m", "translation") |>
    addReaction("P", "", "kdg*P", "degradation") |>
    addReaction("P", "L", "ksec*P", "secretion") |>
    addReaction("L", "", "kL*L", "ligand decay") |>
    addReaction("", "R", "ksR", "receptor synthesis") |>
    addReaction("R", "", "kdR*R", "receptor turnover") |>
    addReaction("L + R", "C", "kon*L*R", "binding") |>
    addReaction("C", "L + R", "koff*C", "unbinding") |>
    addReaction("C", "", "kint*C", "internalisation")
}
.motifs$autocrine_ss <- function() {
  r <- .autocrineReactions()
  ss <- .motifCache("autocrine_ss_steady", function()
    steadyStates(r, verbose = FALSE))
  ev <- addEvent(eventlist(), var = "L", time = 0, value = "dose", method = "add")
  list(f = r, g = eqnvec(ym = "sm*m", yC = "sC*C", yL = "L"), trafo = ss, events = ev,
       conditions = data.frame(dose = c(0, 1), row.names = c("ctrl", "stim")))
}

# ---- receptor pre-equilibrated, knockdown switched on long before the stimulus -------
# M011's event layout: a switch at t = -1000 changes the receptor turnover in one
# condition, the ligand arrives at t = 0; the resting state is the trafo.
.motifs$receptor_kd <- function() {
  r <- eqnlist() |>
    addReaction("", "R", "ksR", "synthesis") |>
    addReaction("R", "", "(kdR + kkd*bool_kd)*R", "turnover") |>
    addReaction("L + R", "C", "kon*L*R", "binding") |>
    addReaction("C", "L + R", "koff*C", "unbinding") |>
    addReaction("C", "Ci", "kint*C", "internalisation") |>
    addReaction("Ci", "", "kdeg*Ci", "degradation") |>
    addReaction("C + S", "C + pS", "kp*C*S", "phosphorylation") |>
    addReaction("pS", "S", "kdp*pS", "dephosphorylation")
  ss <- .motifCache("receptor_kd_steady", function()
    steadyStates(r, forcings = "bool_kd", verbose = FALSE))
  ev <- addEvent(eventlist(), var = "bool_kd", time = -1000, value = "kd", method = "replace")
  ev <- addEvent(ev, var = "L", time = 0, value = "dose", method = "add")
  list(f = r, g = eqnvec(ypS = "s1*pS", yR = "s2*(R + C + Ci)"), trafo = ss, events = ev,
       conditions = data.frame(kd = c(0, 1, 0), dose = c(1, 1, 0),
                               row.names = c("ctrl", "knockdown", "unstim")))
}

# steady states are computed once per motif and kept next to the results
.motifCache <- function(key, fn) {
  dir <- file.path(Sys.getenv("DMOD_MOTIF_CACHE", tempdir()))
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(dir, paste0(key, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  val <- fn(); saveRDS(val, f); val
}
