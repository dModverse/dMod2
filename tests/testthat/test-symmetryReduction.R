# Behavioral tests for symmetryReduction() (constructive symmetry reduction).
#
# Assertions are mathematical, never string comparison: invariants are checked by
# X(I) = 0 (exact via .symExprEqual / numeric tangents), reparametrisations by
# re-running symmetryDetection with the emitted trafo and reading the rank.

redquiet <- function(...) suppressWarnings(symmetryReduction(...))
symdet2 <- function(...) symmetryDetection(..., verbose = FALSE)

# a synthetic symmetrydetection result carrying explicit generators
.mkdir <- function(gen, type = "polynomial", weights = NULL)
  structure(list(type = type, generator = gen, weights = weights,
                 support = names(gen), explicit = TRUE),
            class = "symmetrygenerator")
.mkobj <- function(dirs, coords)
  structure(list(method = "observability", identifiable = FALSE,
                 rank = 1L, dim = length(coords), symmetries = dirs,
                 info = list(engine = "modular", settings = list(),
                             coordinates = coords),
                 call = NULL), class = "symmetrydetection")

# X(I) at a random point for a generator given as named component strings
.lieAt <- function(gen, inv, pt, eps = 1e-7) {
  f <- function(z) eval(parse(text = inv), as.list(z))
  sum(vapply(names(gen), function(v) {
    zp <- pt; zp[v] <- zp[v] + eps
    eval(parse(text = gen[[v]]), as.list(pt)) * (f(zp) - f(pt)) / eps
  }, numeric(1)))
}


test_that("one scaling: transversal pin, family, end-to-end identifiable", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  f <- eqnvec(m = "ktx - dm*m", p = "ktl*m - dp*p")
  g <- eqnvec(y = "p")
  res <- symdet2(f, g, method = "observability", reconstruct = TRUE)
  red <- redquiet(res)
  b <- red$blocks[[1]]
  expect_identical(b$type, "scaling")
  expect_identical(b$status, "reduced")
  # the invariant lattice is spanned by ktx*ktl and ktx/m (weights 1,-1,1)
  expect_true(any(vapply(b$invariants, function(iv)
    .symExprEqual(gsub("\\^", "**", iv), "ktx*ktl") ||
    .symExprEqual(gsub("\\^", "**", iv), "1/(ktx*ktl)"), logical(1))))
  # every listed transversal is one coordinate of the support
  expect_true(all(unlist(red$family[[1]]$admissible) %in% b$support))
  # the emitted trafo removes the direction
  res2 <- symdet2(f, g, method = "observability", trafo = red$trafo)
  expect_true(res2$identifiable)
})


test_that("two chained scalings with a shared coordinate", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  W1 <- list(a = "1", b = "-1")            # a/b trade-off
  W2 <- list(c = "1", a = "-1")            # c/a trade-off (shares a)
  obj <- .mkobj(list(.mkdir(list(a = "a", b = "-b"), "scaling", W1),
                     .mkdir(list(c = "c", a = "-a"), "scaling", W2)),
                c("a", "b", "c", "other"))
  red <- redquiet(obj)
  b <- red$blocks[[1]]
  expect_identical(b$status, "reduced")
  expect_length(b$transversal, 2L)
  # {b, c} is not admissible alone? rank({b,c}) = 2 (rows are independent there),
  # but {a} plus nothing cannot cover both rows: every admissible pair has rank 2
  for (T in red$family[[1]]$admissible)
    expect_identical(dMod2:::.symRedRankModP(
      rbind(c(1L, -1L, 0L), c(-1L, 0L, 1L))[, match(T, c("a", "b", "c")),
                                            drop = FALSE]), 2L)
  # fixing the shared coordinate removes only one direction
  red2 <- redquiet(obj, fixed = "a")
  b2 <- red2$blocks[[1]]
  expect_identical(b2$removedByFixed, 1L)
  expect_identical(b2$status, "reduced")
})


test_that("fixed: rank verdict, redundancy, unknown names", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  W <- list(x = "1", y = "-1", z = "-1")
  obj <- .mkobj(list(.mkdir(list(x = "x", y = "-y", z = "-z"), "scaling", W)),
                c("x", "y", "z"))
  expect_warning(symmetryReduction(obj, fixed = "nothere"), "no effect")
  red <- redquiet(obj, fixed = c("x", "y"))
  b <- red$blocks[[1]]
  expect_identical(b$status, "fixed")
  expect_identical(b$removedByFixed, 1L)
  expect_identical(b$redundantFixed, "y")     # x already kills the single row
  # trafo is pure identity: fixed coordinates are the user's, not pinned here
  expect_true(all(red$trafo == names(red$trafo)))
})


test_that("curved polynomial invariant: solved and identifiable end-to-end", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  f <- eqnvec(A = "k1 + u*k2 - kdeg*A")
  g <- eqnvec(y = "A")
  res <- symdet2(f, g, method = "observability", fixed = "u", reconstruct = TRUE)
  red <- redquiet(res)
  b <- red$blocks[[1]]
  expect_identical(b$status, "reduced")
  expect_identical(b$stage, "polynomial")
  expect_true(any(vapply(b$invariants, function(iv)
    .symExprEqual(gsub("\\^", "**", iv), "k1 + u*k2"), logical(1))))
  # the monomial stage certified negative before escalating
  expect_true(any(grepl("no Laurent-monomial invariant", b$certificates)))
  res2 <- symdet2(f, g, method = "observability", fixed = "u", trafo = red$trafo)
  expect_true(res2$identifiable)
})


test_that("overlapping scaling and curved directions merge into one block", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  f <- eqnvec(A = "-k1*A + k2*B", B = "k1*A - k2*B")
  g <- eqnvec(y = "alpha*A")
  res <- symdet2(f, g, method = "observability", reconstruct = TRUE)

  # the scaling (A, B, alpha) shares B with the general direction, so the two are
  # gauged together instead of through the weight lattice. The grouping is asserted
  # with the search caps switched off; solving the joint invariant set over all five
  # coordinates is the slow free path (see the opt-in test below).
  merged <- redquiet(res, dPoly = 0L, dDarboux = 0L, dExp = 0L)
  expect_length(merged$blocks, 1L)
  expect_identical(merged$blocks[[1]]$type, "curved")
  expect_length(merged$blocks[[1]]$labels, 2L)   # both directions, one block

  # pre-gauging the readout leaves one curved direction, which reduces end-to-end
  red <- redquiet(res, fixed = "alpha")
  expect_setequal(vapply(red$blocks, `[[`, character(1), "status"),
                  c("fixed", "reduced"))
  cb <- red$blocks[[which(vapply(red$blocks, `[[`, character(1),
                                 "status") == "reduced")]]
  expect_length(cb$invariants, 2L)            # k1 + k2 and k2*(A + B)
  res2 <- symdet2(f, g, method = "observability", fixed = "alpha",
                  trafo = red$trafo)
  expect_true(res2$identifiable)
})


test_that("rational invariant via Darboux; degree cap gives a certificate", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  obj <- .mkobj(list(.mkdir(list(x = "x**2", y = "y**2"))), c("x", "y"))
  # separable = FALSE and dPoly = 0 throughout: this generator IS separable, and
  # its invariant (x - y)/(x*y) sits in the rational stage's Laurent ansatz, so
  # both earlier stages would take the block before the factor stages ever run
  # (covered elsewhere) -- here the point is what the factor stages do on their own
  # cap 0 (and exp stage off): coordinate factors only, honest negative
  # certificate -- with dExp > 0 this block now RESOLVES through exp(1/x - 1/y)
  red0 <- redquiet(obj, dPoly = 0L, dDarboux = 0L, dExp = 0L, separable = FALSE)
  b0 <- red0$blocks[[1]]
  expect_identical(b0$status, "unresolved")
  expect_true(any(grepl("coordinate and xi factors only", b0$certificates)))
  # the exp stage picks the block up where the rational stages certified
  # failure; the sign-indefinite argument stays exp-wrapped
  redE <- redquiet(obj, dPoly = 0L, dDarboux = 0L, separable = FALSE)
  expect_identical(redE$blocks[[1]]$stage, "exp")
  expect_gt(length(redE$blocks[[1]]$invariants), 0L)
  expect_true(any(grepl("exp(", redE$blocks[[1]]$invariants, fixed = TRUE)))
  # cap 1: the extactic factors x, y, x - y give (x - y)/(x*y)
  red1 <- redquiet(obj, dPoly = 0L, dDarboux = 1L, separable = FALSE)
  b1 <- red1$blocks[[1]]
  expect_identical(b1$stage, "darboux")
  expect_length(b1$invariants, 1L)
  pt <- c(x = 1.31, y = 0.77)
  expect_lt(abs(.lieAt(list(x = "x^2", y = "y^2"), b1$invariants[1], pt)), 1e-4)
})


test_that("module reduction collapses the M009-shaped pair", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  X12 <- list(C = "Ci*rr + Ci", Ci = "-(Ci*rr - Ci)",
              k2i = "-(k2*rr - k2 + k2i*rr + k2i)",
              r1 = "r1*rr + r1", r2 = "r2*rr + r2 + rr")
  X13 <- list(C = "Ci*r2*rr + Ci*r2", Ci = "-(Ci*r2*rr - Ci*r2)",
              k2i = "-(k2*r2*rr - k2*r2 + k2i*r2*rr + k2i*r2)",
              r1 = "-r1*rr",
              rr = "-(r2*rr**3 - 2*r2*rr**2 - r2*rr - rr**3 - rr**2)")
  obj <- .mkobj(list(.mkdir(X12), .mkdir(X13)),
                c("C", "Ci", "k2", "k2i", "r1", "r2", "rr"))
  red <- redquiet(obj)
  b <- red$blocks[[1]]
  # the discovered combination (-r2)*X1 + X2 lives on 3 coordinates
  expect_true(any(grepl("support 3", b$moduleCombos)))
})


test_that("multi-generator block: only common invariants survive", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  # X1 scales (a, b); X2 scales (b, c): the common monomial lattice is spanned by
  # a*b*c (exponents alpha = beta = gamma) -- a*b is X1-invariant but not
  # X2-invariant, so the joint block must NOT report it
  obj <- .mkobj(list(.mkdir(list(a = "a", b = "-b"), "general"),
                     .mkdir(list(b = "b", c = "-c"), "general")),
                c("a", "b", "c"))
  red <- redquiet(obj)
  b <- red$blocks[[1]]
  expect_length(b$invariants, 1L)
  expect_true(.symExprEqual(gsub("\\^", "**", b$invariants[1]), "a*b*c") ||
              .symExprEqual(gsub("\\^", "**", b$invariants[1]), "1/(a*b*c)"))
})


test_that("support-only directions and stripped coordinates degrade gracefully", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  d <- structure(list(type = "general", generator = NULL, weights = NULL,
                      support = c("p", "q"), explicit = FALSE),
                 class = "symmetrygenerator")
  obj <- .mkobj(list(d), c("p", "q", "r"))
  red <- redquiet(obj)
  b <- red$blocks[[1]]
  expect_identical(b$status, "unresolved")
  expect_match(b$reason, "reconstruct = TRUE")
  # no coordinate list: support-union fallback with a warning
  obj$info$coordinates <- NULL
  expect_warning(red2 <- symmetryReduction(obj), "no coordinate list")
  expect_setequal(names(red2$trafo), c("p", "q"))
})


test_that("trafo entries parse and carry integer powers only", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  f <- eqnvec(A = "-k1*A + k2*B", B = "k1*A - k2*B")
  g <- eqnvec(y = "alpha*A")
  res <- symdet2(f, g, method = "observability", reconstruct = TRUE)
  red <- redquiet(res, fixed = "alpha")
  expect_true(inherits(red$trafo, "eqnvec"))
  expect_setequal(names(red$trafo), red$coordinates)
  for (v in red$trafo) {
    expect_error(parse(text = v), NA)
    expect_false(grepl("\\^\\s*\\(?\\s*[0-9]*\\.[0-9]", v))   # no fractional powers
  }
})


test_that("identifiable results return an empty reduction", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  f <- eqnvec(x = "-k*x")
  res <- symdet2(f, eqnvec(y = "x"), method = "observability")
  red <- redquiet(res)
  expect_length(red$blocks, 0L)
  expect_null(red$trafo)
  expect_output(print(red), "Nothing to reduce")
})


test_that("print stays lean; summary carries the block report", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  f <- eqnvec(m = "ktx - dm*m", p = "ktl*m - dp*p")
  res <- symdet2(f, eqnvec(y = "p"), method = "observability", reconstruct = TRUE)
  red <- redquiet(res)
  out <- capture.output(print(red))
  expect_true(any(grepl("^Reduced 1 of 1 direction", out)))
  expect_true(any(grepl("Trafo", out)))
  expect_true(any(grepl("ktx *= ktx\\*ktl", out)))          # aligned meaning line
  # nothing of the report leaks into print()
  expect_false(any(grepl("Scaling block|transversal:|admissible|\\[", out)))
  expect_lt(length(out), 12L)

  # summary() adds one line per block; the invariants of a REDUCED block are the
  # Invariants list above, so they are not repeated, and no certificate is printed
  rep <- capture.output(summary(red))
  expect_true(any(grepl("^Blocks$", rep)))
  # [^}] not . -- under an ASCII session charset the label subscript reaches the
  # captured line as its raw bytes, which a byte-wise "." cannot span
  expect_true(any(grepl("\\{X[^}]+\\} scaling, reduced \\| transversal ktl = 1", rep)))
  expect_true(any(grepl("admissible", rep)))
  expect_true(any(grepl("ktx *= ktx\\*ktl", rep)))
  expect_false(any(grepl("invariants  |certificate|certified", rep)))
  expect_lt(length(rep), 20L)
  # verbose adds the admissible sets and the raw invariants -- still no certificates
  vrb <- capture.output(summary(red, verbose = TRUE))
  expect_true(any(grepl("admissible  \\{", vrb)))
  expect_true(any(grepl("invariants  ", vrb)))
  expect_false(any(grepl("\\[transversal certified|\\[verified", vrb)))
  expect_gt(length(vrb), length(rep))
})


test_that("EGF/MEK/ERK cascade: all four directions reduced end-to-end", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  # The receptor scaling merges with the rational confounder into one curved block
  # that is gauged FREELY, with no pre-fixing: the widest path the solver takes.
  reactions <- eqnlist() |>
    addReaction("EGF + EGFR", "EGF_EGFR", "k_bind * EGF * EGFR")   |>
    addReaction("EGF_EGFR", "EGF + EGFR", "k_unbind * EGF_EGFR")   |>
    addReaction("MEK", "pMEK", "k_phos_MEK * EGF_EGFR * MEK")      |>
    addReaction("pMEK", "MEK", "k_dephos_MEK * pMEK")              |>
    addReaction("ERK", "pERK", "k_phos_ERK * pMEK * ERK")          |>
    addReaction("pERK", "ERK", "k_dephos_ERK * pERK")
  observables <- eqnvec(pMEK_obs = "scale_pMEK * pMEK",
                        pERK_obs = "scale_pERK * pERK")
  egf <- symdet2(reactions, observables, method = "observability",
                 reduceCQ = FALSE, reconstruct = TRUE)
  red <- redquiet(egf)
  # two clean unit scalings plus the receptor scaling merged with the rational
  # confounder into one curved block
  types <- vapply(red$blocks, `[[`, character(1), "type")
  expect_identical(sort(types), c("curved", "scaling", "scaling"))
  cb <- red$blocks[[which(types == "curved")]]
  expect_length(cb$labels, 2L)
  expect_length(cb$invariants, 4L)
  # the confounder invariant totalEGF*totalEGFR/EGF_EGFR^2 is among them
  expect_true(any(vapply(cb$invariants, function(iv) .symExprEqual(
    gsub("\\^", "**", iv),
    "(EGF + EGF_EGFR)*(EGFR + EGF_EGFR)/EGF_EGFR**2"), logical(1))))

  # the invariants obey sharp mutual bounds: a positive chart requires offsets
  expect_identical(cb$status, "reduced")
  expect_length(red$remaining, 0L)
  expect_true(any(grepl("carrier offset", cb$certificates)))

  # the chart stays positive for every positive outer value, extremes included
  vals <- c(0.7, 2.3, 1.9, 0.31, 1.7, 0.45, 1.1, 0.9)
  outer <- unique(c(red$coordinates, unlist(lapply(red$trafo, function(e)
    all.vars(parse(text = e))))))
  uAt <- function(scale) setNames(as.list(
    scale * rep_len(vals, length(outer))), outer)
  chart <- function(u) vapply(names(red$trafo), function(nm)
    eval(parse(text = red$trafo[[nm]]), u), numeric(1))
  for (scale in c(1e-3, 1, 1e3)) {
    z <- chart(uAt(scale))
    expect_true(all(is.finite(z)) && all(z > 0))
  }

  # survivorMeaning at the chart point returns the outer value
  u <- uAt(1)
  z <- as.list(chart(u))
  for (nm in names(cb$survivorMeaning))
    expect_equal(eval(parse(text = cb$survivorMeaning[[nm]]), z), u[[nm]],
                 tolerance = 1e-8)

  # coverage: the carried quantities are positive at arbitrary inner points
  for (z0 in list(list(EGF = 1.3, EGFR = 0.2, EGF_EGFR = 2.5, k_bind = 0.8,
                       k_phos_MEK = 0.3, k_unbind = 1.9),
                  list(EGF = 1e-4, EGFR = 1e-4, EGF_EGFR = 50, k_bind = 3,
                       k_phos_MEK = 2, k_unbind = 1e-3)))
    for (nm in names(cb$survivorMeaning))
      expect_gt(eval(parse(text = cb$survivorMeaning[[nm]]), z0), 0)
})


test_that("EGF cascade with fixed EGF: exp-stage argument unwraps to a rational
           invariant and the block reduces", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  # a rational invariant reachable only through the exp stage (mu = 0) is
  # reported unwrapped when sign-certified, keeping the section search available
  reactions <- eqnlist() |>
    addReaction("EGF + EGFR", "EGF_EGFR", "k_bind * EGF * EGFR")   |>
    addReaction("EGF_EGFR", "EGF + EGFR", "k_unbind * EGF_EGFR")   |>
    addReaction("MEK", "pMEK", "k_phos_MEK * EGF_EGFR * MEK")      |>
    addReaction("pMEK", "MEK", "k_dephos_MEK * pMEK")              |>
    addReaction("ERK", "pERK", "k_phos_ERK * pMEK * ERK")          |>
    addReaction("pERK", "ERK", "k_dephos_ERK * pERK")
  observables <- eqnvec(pMEK_obs = "scale_pMEK * pMEK",
                        pERK_obs = "scale_pERK * pERK")
  fx <- c("EGF", "scale_pMEK", "scale_pERK")
  det <- symdet2(reactions, observables, method = "observability",
                 reduceCQ = FALSE, reconstruct = TRUE, fixed = fx)
  red <- redquiet(det, fixed = fx)
  b <- red$blocks[[1]]
  expect_identical(b$status, "reduced")
  expect_length(red$remaining, 0L)
  expect_false(any(grepl("exp(", b$invariants, fixed = TRUE)))
  expect_true(any(vapply(b$invariants, function(iv) .symExprEqual(
    gsub("\\^", "**", iv),
    "(EGF*EGFR + EGF*EGF_EGFR + EGFR*EGF_EGFR)/EGF_EGFR**2"), logical(1))))
  # chart positive; carrier meanings round-trip (fixed coordinates stay outer)
  tr <- red$trafo[red$trafo != names(red$trafo)]
  unames <- unique(c(red$coordinates,
                     unlist(lapply(tr, function(e) all.vars(parse(text = e))))))
  u <- setNames(as.list(rep_len(c(0.7, 2.3, 1.9, 0.31, 1.7, 0.45),
                                length(unames))), unames)
  z <- vapply(names(red$trafo), function(nm)
    eval(parse(text = red$trafo[[nm]]), u), numeric(1))
  expect_true(all(is.finite(z)) && all(z > 0))
  zl <- c(as.list(z), u[setdiff(names(u), names(z))])
  for (nm in names(b$survivorMeaning))
    expect_equal(eval(parse(text = b$survivorMeaning[[nm]]), zl), u[[nm]],
                 tolerance = 1e-8)

  # with the coordinate substituted away instead of fixed, the polynomial stage
  # meets the target with a sign-indefinite invariant; the escalation continues
  # and a definite rational replacement wins the slot
  det2 <- symdet2(reactions, observables, method = "observability",
                  reduceCQ = FALSE, reconstruct = TRUE,
                  trafo = eqnvec(EGF = "1"))
  red2 <- redquiet(det2)
  b2 <- Filter(function(bb) identical(bb$type, "curved"), red2$blocks)[[1]]
  expect_identical(b2$status, "reduced")
  expect_length(b2$support, 5L)
  expect_true(all(vapply(b2$invariants, function(iv) {
    e <- parse(text = iv)
    v <- vapply(1:3, function(s) eval(e, setNames(as.list(
      s * c(0.3, 2.1, 0.7, 1.9, 0.4)), b2$support)), numeric(1))
    all(v > 0)
  }, logical(1))))
  u2names <- unique(c(red2$coordinates, unlist(lapply(red2$trafo, function(e)
    all.vars(parse(text = e))))))
  u2 <- setNames(as.list(rep_len(c(0.7, 2.3, 1.9, 0.31, 1.7, 0.45),
                                 length(u2names))), u2names)
  z2 <- vapply(names(red2$trafo), function(nm)
    eval(parse(text = red2$trafo[[nm]]), u2), numeric(1))
  expect_true(all(is.finite(z2)) && all(z2 > 0))
  # one fresh q_<k> per invariant, numbered gaplessly, and every carrier meaning
  # round-trips at the emitted point
  expect_identical(names(b2$survivorMeaning),
                   paste0("q_", seq_along(b2$invariants)))
  qs <- unique(unlist(lapply(red2$trafo, function(e) all.vars(parse(text = e)))))
  expect_setequal(qs[grepl("^q_[0-9]+$", qs)], names(b2$survivorMeaning))
  # the shifted carrier solves to its parameter outright
  expect_true(any(trimws(unlist(red2$trafo)) %in% names(b2$survivorMeaning)))
  zl2 <- c(as.list(z2), u2[setdiff(names(u2), names(z2))])
  for (nm in names(b2$survivorMeaning))
    expect_equal(eval(parse(text = b2$survivorMeaning[[nm]]), zl2), u2[[nm]],
                 tolerance = 1e-8)
})


test_that("EGF cascade with steady-state trafo: the rational stage finds the
           single-denominator invariant and the offset chart certifies", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  # the direction symmetryDetection(trafo = steadyStates(...), events = ...)
  # reports for the EGF/EGFR -> MEK/ERK cascade with an EGF dose event. Its 4th
  # invariant (EGFR^2*k_bind + EGFR*k_bind + EGFR*k_unbind + EGF_EGFR*k_unbind)/
  # EGFR is polynomial over the single coordinate EGFR -- the rational (Laurent)
  # stage's case, which no other stage reaches at the default caps. The chart
  # needs the balance section EGFR = 1 and a certified offset on the carrier of
  # the 4th invariant (q_4 rides I_4 - 2*sqrt(I_1), positive by AM-GM).
  gen <- list(EGFR = "2*EGFR*(EGFR + EGF_EGFR)",
              EGF_EGFR = "EGF_EGFR*(EGFR + EGF_EGFR)",
              k_bind = "-k_bind*(EGFR + EGF_EGFR)",
              k_phos_MEK = "-k_phos_MEK*(EGFR + EGF_EGFR)",
              k_unbind = "-EGFR**2*k_bind + EGFR*k_bind + EGF_EGFR*k_unbind")
  obj <- .mkobj(list(.mkdir(gen, "general")),
                c("EGFR", "EGF_EGFR", "pMEK", "pERK", "k_bind", "k_dephos_ERK",
                  "k_dephos_MEK", "k_phos_ERK", "k_phos_MEK", "k_unbind"))
  red <- redquiet(obj)                      # default caps: the point of the stage
  b <- red$blocks[[1]]
  expect_identical(b$status, "reduced")
  expect_identical(b$stage, "rational")
  expect_length(b$invariants, 4L)
  expect_true(any(grepl("carrier offset", b$certificates)))
  genR <- lapply(gen, function(x) gsub("\\*\\*", "^", x))
  pt <- c(EGFR = 0.7, EGF_EGFR = 2.3, k_bind = 1.9, k_phos_MEK = 0.31,
          k_unbind = 1.7)
  for (iv in b$invariants) expect_lt(abs(.lieAt(genR, iv, pt)), 1e-4)
  # the chart stays positive for every positive outer value, extremes included
  outer <- unique(c(red$coordinates, unlist(lapply(red$trafo, function(e)
    all.vars(parse(text = e))))))
  vals <- c(0.7, 2.3, 1.9, 0.31, 1.7, 0.45, 1.1, 0.9)
  for (scale in c(1e-3, 1, 1e3)) {
    u <- setNames(as.list(scale * rep_len(vals, length(outer))), outer)
    z <- vapply(names(red$trafo), function(nm)
      eval(parse(text = red$trafo[[nm]]), u), numeric(1))
    expect_true(all(is.finite(z)) && all(z > 0))
  }
  # survivorMeaning at the chart point returns the outer value
  u <- setNames(as.list(rep_len(vals, length(outer))), outer)
  z <- as.list(vapply(names(red$trafo), function(nm)
    eval(parse(text = red$trafo[[nm]]), u), numeric(1)))
  for (nm in names(b$survivorMeaning))
    expect_equal(eval(parse(text = b$survivorMeaning[[nm]]), z), u[[nm]],
                 tolerance = 1e-8)
})


test_that("exponential-factor invariant: found, certified, solved and printed", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  obj <- .mkobj(list(.mkdir(list(x = "x", y = "x + y"))), c("x", "y"))
  red <- redquiet(obj)
  b <- red$blocks[[1]]
  expect_identical(b$stage, "exp")
  expect_identical(b$status, "reduced")
  expect_length(b$invariants, 1L)
  # X(I) = 0 numerically (the verify layer already proved it symbolically)
  gen <- list(x = "x", y = "x + y")
  expect_lt(abs(.lieAt(gen, b$invariants[[1]], c(x = 1.7, y = 0.9))), 1e-4)
  # the emitted trafo reproduces the invariant value carried on the fresh
  # parameter: I(trafo(u)) is u or 1/u
  expect_identical(names(b$survivorMeaning), "q_1")
  z <- vapply(names(red$trafo), function(nm)
    eval(parse(text = red$trafo[[nm]]), list(x = 5, y = 3, q_1 = 3)),
    numeric(1))
  Iv <- eval(parse(text = b$invariants[[1]]), as.list(z))
  expect_true(abs(Iv - 3) < 1e-10 || abs(Iv - 1/3) < 1e-10)
  txt <- capture.output(summary(red, verbose = TRUE))
  expect_true(any(grepl("[exp]", txt, fixed = TRUE)))
  expect_true(any(grepl("gauge pin", txt)))                # wrap-safe fragment
  # the certificates stay reachable on the object, they are just not printed
  expect_true(any(grepl("exp stage", red$blocks[[1]]$certificates)))
  # dExp = 0 switches the stage off with a certificate
  red0 <- redquiet(obj, dExp = 0L)
  expect_identical(red0$blocks[[1]]$status, "unresolved")
  expect_true(any(grepl("exp stage skipped", red0$blocks[[1]]$certificates)))
})


test_that("root carrier: quadratic invariant solved with a fractional power", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  # the rotation a' = b, b' = -a has the single invariant a^2 + b^2. Pinning a = 1
  # would solve b = sqrt(I - 1), a chart valid only for I > 1 -- the orbit is the
  # circle of radius sqrt(I) and reaches a = 1 only when it is large enough. The
  # balance section a = b has no such wall: it meets every circle once, in the
  # positive quadrant, and puts both coordinates at sqrt(2*I)/2.
  obj <- .mkobj(list(.mkdir(list(a = "b", b = "-a"))), c("a", "b"))
  red <- redquiet(obj)
  b <- red$blocks[[1]]
  expect_identical(b$status, "reduced")
  expect_match(b$gaugeNote, "gauge section")
  expect_match(b$gaugeNote, "certified positive")
  expect_true(all(vapply(red$trafo, function(e)
    grepl("sqrt|\\^\\(1/2\\)", e), logical(1))))
  # numeric round-trip: the invariant value is carried by the fresh parameter
  z <- vapply(names(red$trafo), function(nm)
    eval(parse(text = red$trafo[[nm]]), list(a = 1, b = 1, q_1 = 9)),
    numeric(1))
  expect_equal(unname(z[["a"]]^2 + z[["b"]]^2), 9, tolerance = 1e-8)
  # and the emitted point is positive for every positive outer value, which the old
  # constant pin was not
  z2 <- vapply(names(red$trafo), function(nm)
    eval(parse(text = red$trafo[[nm]]), list(a = 1, b = 1, q_1 = 0.01)),
    numeric(1))
  expect_true(all(z2 > 0))
})


test_that("fixed removes curved directions that cannot avoid the coordinate", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  obj <- .mkobj(list(.mkdir(list(a = "a", b = "-b")),
                     .mkdir(list(b = "b", cc = "cc"))), c("a", "b", "cc"))
  red <- redquiet(obj, fixed = "a")
  st <- vapply(red$blocks, `[[`, character(1), "status")
  expect_true("fixed" %in% st)
  fb <- red$blocks[[which(st == "fixed")]]
  expect_identical(fb$labels, "X₁")
  kb <- red$blocks[[which(st != "fixed")]]
  expect_identical(kb$status, "reduced")
  expect_false("a" %in% kb$support)
  expect_setequal(red$removed, c("X₁", "X₂"))
  # a single direction moving the fixed coordinate disappears entirely
  obj2 <- .mkobj(list(.mkdir(list(a = "a^2", b = "a*b"))), c("a", "b"))
  red2 <- redquiet(obj2, fixed = "a")
  expect_identical(red2$blocks[[1]]$status, "fixed")
  expect_true(all(red2$trafo == names(red2$trafo)))
})


test_that("moved-only extactic basis: parameters no longer block the factor search", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  s <- "(k1 + k2 + k3 + k4)"
  obj <- .mkobj(list(.mkdir(list(x = paste0(s, "*x**2"), y = paste0(s, "*y**2")))),
                c("x", "y", "k1", "k2", "k3", "k4"))
  red <- redquiet(obj, dPoly = 0L, dDarboux = 1L, separable = FALSE)
  b <- red$blocks[[1]]
  # separable = FALSE and dPoly = 0: the generator decouples and its invariant is
  # Laurent, so the quadrature and rational stages would each answer before the
  # extactic basis is ever built -- this test is about the basis
  # 6 total variables used to hit "extactic skipped" (projected entry degree 13
  # over the all-variables basis); the moved-coordinate basis (2 coordinates)
  # passes the cap and finds the factor x - y
  expect_true(any(grepl("extactic complete", b$certificates)))
  expect_identical(b$stage, "darboux")
  expect_gt(length(b$invariants), 0L)
  gen <- list(x = paste0(s, "*x^2"), y = paste0(s, "*y^2"))
  pt <- c(x = 1.31, y = 0.77, k1 = 0.2, k2 = 0.3, k3 = 0.15, k4 = 0.4)
  expect_lt(abs(.lieAt(gen, b$invariants[1], pt)), 1e-4)
})


test_that("separable characteristics: quadrature before the factor stages", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  # Support 3, so the 2-coordinate integrating-factor stage cannot reach it, and
  # atan is outside the Darboux language: only the quadrature stage gets these.
  obj <- .mkobj(list(.mkdir(list(x = "1 + x**2", y = "1 + y**2", z = "1 + z**2"))),
                c("x", "y", "z"))
  b <- redquiet(obj)$blocks[[1]]
  expect_identical(b$stage, "separable")
  expect_length(b$invariants, 2L)          # target = 3 coordinates - rank 1
  gen <- list(x = "1 + x^2", y = "1 + y^2", z = "1 + z^2")
  for (iv in b$invariants)
    expect_lt(abs(.lieAt(gen, iv, c(x = 0.7, y = 1.3, z = 2.1))), 1e-4)
  # atan cannot be emitted as a trafo entry, so the block stays invariantOnly
  expect_identical(b$status, "invariantOnly")

  # a separable block whose quadrature IS rational reduces end-to-end
  obj2 <- .mkobj(list(.mkdir(list(x = "x**2", y = "y**2", z = "z**2"))),
                 c("x", "y", "z"))
  b2 <- redquiet(obj2)$blocks[[1]]
  expect_identical(b2$stage, "separable")
  expect_identical(b2$status, "reduced")
  expect_true(any(vapply(b2$invariants, function(iv) .symExprEqual(
    gsub("\\^", "**", iv), "(x - y)/(x*y)"), logical(1))))

  # a coupled generator is declined, not invented: xi_y involves x
  obj3 <- .mkobj(list(.mkdir(list(x = "x", y = "x + y"))), c("x", "y"))
  b3 <- redquiet(obj3)$blocks[[1]]
  expect_false(identical(b3$stage, "separable"))

  # the switch turns the stage off with a certificate (on the 2-coordinate block:
  # forcing the support-3 one through the factor stages costs a minute)
  obj4 <- .mkobj(list(.mkdir(list(x = "x**2", y = "y**2"))), c("x", "y"))
  b4 <- redquiet(obj4, separable = FALSE)$blocks[[1]]
  expect_false(identical(b4$stage, "separable"))
  expect_true(any(grepl("separable-characteristics stage skipped",
                        b4$certificates)))
})


test_that("integrating factor: Liouvillian integral by quadrature, verified", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  # damped oscillator: no rational or exp-factor invariant; M = 1/(x^2 + x*y + y^2)
  obj <- .mkobj(list(.mkdir(list(x = "y", y = "-x - y"))), c("x", "y"))
  red <- redquiet(obj)
  b <- red$blocks[[1]]
  expect_identical(b$stage, "intfactor")
  expect_identical(b$status, "invariantOnly")   # atan form: reported, not solved
  expect_length(b$invariants, 1L)
  expect_true(any(grepl("solved by quadrature", b$certificates)))
  expect_true(any(grepl("verified: X\\(I\\) = 0", b$certificates)))
  gen <- list(x = "y", y = "-x - y")
  pt <- c(x = 0.83, y = 0.41)
  expect_lt(abs(.lieAt(gen, b$invariants[[1]], pt)), 1e-4)
})


test_that("a carrier that takes both signs is declared real and gets a covering pin", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # Only the difference a - b enters, so a - b is the only invariant and it takes
  # both signs on the positive orthant. A constant gauge pin (b = 1) is reached only
  # by the orbits with a - b > -1, which leaves the rest of the positive orthant
  # outside the chart; the sum-of-squares pin clears every orbit at once.
  r <- symdet2(eqnvec(x = "-(a - b)*x"), eqnvec(y = "x"),
               method = "observability", reconstruct = TRUE)
  red <- redquiet(r, verbose = FALSE)
  b1 <- red$blocks[[1]]
  expect_equal(b1$status, "reduced")
  expect_equal(unname(b1$coverage), "total")
  expect_equal(unname(b1$carrierDomain[["q_1"]]), "real")

  # the chart stays in the positive orthant for EVERY real carrier value ...
  q <- c(-50, -3, -1, -0.4, 0, 0.4, 1, 3, 50)
  z <- lapply(red$trafo[c("a", "b")], function(e) eval(parse(text = e), list(q_1 = q)))
  expect_true(all(z$a > 0))
  expect_true(all(z$b > 0))
  # ... and lands on the orbit the carrier names
  expect_equal(eval(parse(text = b1$survivorMeaning[["q_1"]]), z), q)
})


test_that("a certified-positive carrier keeps its positive domain", {
  if (!.sympy_works()) skip("reticulate/sympy not available")

  # the moiety total A + B is positive on the positive orthant by inspection, so the
  # carrier stays log-fittable and no sum-of-squares pin is needed
  r <- symdet2(eqnvec(A = "-k1*A + k2*B", B = "k1*A - k2*B"), eqnvec(y = "A + B"),
               method = "observability", reconstruct = TRUE, reduceCQ = FALSE)
  red <- redquiet(r, verbose = FALSE)
  dom <- unlist(lapply(red$blocks, `[[`, "carrierDomain"))
  expect_true(length(dom) == 0L || all(dom == "positive"))
})


test_that("zero limits: unconditional ones reported, conditional ones stated", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  reactions <- eqnlist() |>
    addReaction("P",  "pP", "k_p*P",  "phosphorylation") |>
    addReaction("pP", "P",  "k_d*pP", "dephosphorylation")
  res <- symdet2(reactions, eqnvec(y = "s*pP"), method = "observability",
                 reconstruct = TRUE)
  red <- redquiet(res, fixed = "s", reportZeroCompatibility = TRUE)
  v <- red$zeroCompatibility
  rownames(v) <- v$coordinates

  # X2 is xi_P = P + pP, xi_k_d = k_p, xi_k_p = -k_p, so the flow is closed form:
  #   k_p(e) = k_p*exp(-e), k_d(e) = k_d + k_p*(1 - exp(-e)), P(e) = (P + pP)*exp(e) - pP
  # (pP unmoved). That is the ground truth every assertion below is checked against.
  kdZero <- function(z) -log1p(z[["k_d"]] / z[["k_p"]])   # eps where k_d hits 0
  flowP  <- function(e, z) (z[["P"]] + z[["pP"]]) * exp(e) - z[["pP"]]

  # k_p*(P + pP) is an invariant, strictly positive on the orthant: {k_p = 0} is
  # off every orbit, whatever the point -- and the flow agrees, k_p only decays
  expect_identical(v["k_p", "verdict"], "no")
  # the other two are reachable on complementary halves of the parameter space
  expect_identical(unname(v[c("P", "k_d"), "verdict"]), c("if", "if"))
  # neither face is invariant (v does not divide xi_v): finite eps, not a limit
  expect_false(any(v$limit[v$verdict == "if"]))
  expect_true(all(v$certain))

  t1 <- c(P = 1,   pP = 0.2, k_p = 0.3,  k_d = 0.05, s = 1)   # above steady state
  t2 <- c(P = 0.2, pP = 1,   k_p = 0.05, k_d = 0.3,  s = 1)   # below it
  # the flow decides it: at t1 the orbit reaches {k_d = 0} with P still positive,
  # at t2 P has left the orthant before k_d gets there
  expect_gt(flowP(kdZero(t1), t1), 0)
  expect_lt(flowP(kdZero(t2), t2), 0)
  # the reported conditions are R over the model's own names: one eval decides them,
  # and they agree with the flow at both points
  dec <- function(cond, z) eval(parse(text = cond), as.list(z))
  expect_true(dec(v["k_d", "condition"], t1))
  expect_false(dec(v["P", "condition"], t1))
  expect_false(dec(v["k_d", "condition"], t2))
  expect_true(dec(v["P", "condition"], t2))
  # the degenerate point that is reported IS where the flow lands
  pAt <- sub(".*\\bP = ([^,]+).*", "\\1", v["k_d", "at"])
  expect_equal(eval(parse(text = pAt), as.list(t1)), flowP(kdZero(t1), t1))

  # neither zero is announced as reachable; the condition is what is reported
  out <- capture.output(print(red))
  expect_true(any(grepl("k_d = 0 +where +P\\*k_p > k_d\\*pP", out)))
  expect_true(any(grepl("P += 0 +where +k_d\\*pP > P\\*k_p", out)))
  # a zero that cannot happen is not news: no k_p row (its trafo entry still prints)
  expect_false(any(grepl("^ +k_p += 0", out)))
  expect_false(any(grepl("everywhere", out)))
})


test_that("zero limits: unconditional case, and the switch", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  # A = k1 + u*k2 - kdeg*A with u fixed: the invariant is k1 + u*k2, so either
  # production term can be moved into the other from ANY positive point
  res <- symdet2(eqnvec(A = "k1 + u*k2 - kdeg*A"), eqnvec(y = "A"),
                 method = "observability", fixed = "u", reconstruct = TRUE)
  red <- redquiet(res, reportZeroCompatibility = TRUE)
  v <- red$zeroCompatibility
  expect_setequal(v$coordinates, c("k1", "k2"))
  expect_true(all(v$verdict == "yes"))
  expect_true(all(!nzchar(v$condition)))
  # the reported degenerate points carry the invariant unchanged
  pt <- c(k1 = 0.7, k2 = 1.3, u = 2, kdeg = 0.4, A = 1.1)
  inv <- function(z) z[["k1"]] + z[["u"]] * z[["k2"]]
  for (i in seq_len(nrow(v))) {
    z <- pt
    for (e in strsplit(v$at[i], ", ", fixed = TRUE)[[1]]) {
      kv <- strsplit(e, " = ", fixed = TRUE)[[1]]
      z[kv[1]] <- eval(parse(text = kv[2]), as.list(pt))
    }
    expect_equal(inv(z), inv(pt))
    expect_true(all(z[v$coordinates] >= 0))
  }
  expect_true(any(grepl("k1 = 0 +everywhere", capture.output(print(red)))))
  expect_null(redquiet(res)$zeroCompatibility)
})


test_that("zero limits: coordinates that can only vanish together", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  # a and b scale with the SAME weight: eps -> -Inf takes both to zero and moves
  # nothing else, so {a, b} is a zero limit and neither coordinate is one alone.
  # The invariant b/a is 0/0 there, which is why the set has to be found by growing
  # the face rather than by testing it.
  obj <- .mkobj(list(.mkdir(list(a = "a", b = "b"), "scaling",
                            list(a = "1", b = "1"))), c("a", "b", "other"))
  v <- redquiet(obj, reportZeroCompatibility = TRUE)$zeroCompatibility
  expect_identical(v$coordinates, "a, b")
  expect_identical(v$verdict, "yes")
  expect_true(v$limit)                       # a scaling never crosses {a = 0}
  # opposite weights: a -> 0 only with b -> Inf, which is no zero limit at all
  obj2 <- .mkobj(list(.mkdir(list(a = "a", b = "-b"), "scaling",
                             list(a = "1", b = "-1"))), c("a", "b", "other"))
  v2 <- redquiet(obj2, reportZeroCompatibility = TRUE)$zeroCompatibility
  expect_true(all(v2$verdict == "no"))
  # nothing can vanish, so the section does not appear at all
  expect_false(any(grepl("Zero limits",
                         capture.output(print(redquiet(obj2,
                           reportZeroCompatibility = TRUE))))))
})
