# Evaluation protocol behind the composition operators (R/fnProtocol.R):
# condition truth table, bundle alignment, and protocol-vs-ordinary equality.

## ---- Condition truth table ----------------------------------------------

test_that(".resolveConditions reproduces the leaf truth table", {
  rc <- dMod2:::.resolveConditions

  # unspecific leaf, unspecific request
  r <- rc(NULL, NULL)
  expect_null(r$conditions)
  expect_true(r$evaluate)
  expect_identical(r$slots, 1L)

  # unspecific leaf, several conditions: evaluate once, fill all slots
  r <- rc(c("C1", "C2"), NULL)
  expect_identical(r$conditions, c("C1", "C2"))
  expect_true(r$evaluate)
  expect_identical(r$slots, 1:2)

  # specific leaf, unspecific request
  r <- rc(NULL, "C1")
  expect_identical(r$conditions, "C1")
  expect_true(r$evaluate)
  expect_identical(r$slots, 1L)

  # specific leaf in a wider request: only its own slot
  r <- rc(c("C1", "C2"), "C1")
  expect_identical(r$conditions, c("C1", "C2"))
  expect_true(r$evaluate)
  expect_identical(r$slots, 1L)

  # disjoint: nothing evaluated
  r <- rc("C2", "C1")
  expect_identical(r$conditions, "C2")
  expect_false(r$evaluate)
  expect_length(r$slots, 0L)
})


test_that(".resolveConditions handles a leaf owning several conditions", {
  # Not reachable before the rebuild; the PEtab relabeller produces it.
  r <- dMod2:::.resolveConditions(c("C1", "C2", "C3"), c("C1", "C3"))
  expect_true(r$evaluate)
  expect_identical(r$slots, c(1L, 3L))
})


test_that(".emptySlots returns a NULL-holed named list", {
  s <- dMod2:::.emptySlots(c("C1", "C2"))
  expect_length(s, 2L)
  expect_identical(names(s), c("C1", "C2"))
  expect_true(all(vapply(s, is.null, logical(1))))

  # unnamed single slot
  s <- dMod2:::.emptySlots(NULL)
  expect_length(s, 1L)
  expect_null(names(s))
})


## ---- Bundles -------------------------------------------------------------

test_that("a bundle built from a public call broadcasts and stays aligned", {
  b <- dMod2:::.bundle_from_call(c("C1", "C2", "C3"), times = 0:3,
                                out = NULL, pars = c(a = 1), fixed = NULL)
  expect_true(b$shared)
  expect_identical(dMod2:::.bundle_n(b), 3L)
  expect_identical(dMod2:::.req_pars(b, 2L), c(a = 1))
  # shared times broadcasts
  expect_identical(dMod2:::.req_times(b, 3L), 0:3)
})


test_that(".bundle_subset keeps per-request alignment", {
  b <- dMod2:::.bundle(conds = c("C1", "C2", "C3"),
                      times = list(1, 2, 3),
                      pars  = list(c(a = 1), c(a = 2), c(a = 3)))
  s <- dMod2:::.bundle_subset(b, c("C3", "C1"))
  expect_identical(s$conds, c("C3", "C1"))
  expect_identical(dMod2:::.req_pars(s, 1L), c(a = 3))
  expect_identical(dMod2:::.req_times(s, 2L), 1)

  # unknown conditions are dropped, not NA
  s <- dMod2:::.bundle_subset(b, c("C2", "nope"))
  expect_identical(s$conds, "C2")
})


test_that(".splitParsFixed makes pars and fixed disjoint", {
  p <- as.parvec(c(a = 1, b = 2))
  r <- dMod2:::.splitParsFixed(p, c(b = 9))
  expect_identical(names(r$pars), "a")
  expect_identical(names(r$fixed), "b")

  # no fixed: pars unchanged
  r <- dMod2:::.splitParsFixed(p, NULL)
  expect_identical(r$pars, p)
  expect_null(r$fixed)
})


## ---- Driving a real chain through the protocol ---------------------------

test_that(".evalLegacy reproduces an ordinary multi-condition call", {
  skip_if_not_installed("cppDE")
  skip_on_cran()

  fx <- fx_decay_multicond_compiled()
  times <- seq(0, 5, by = 0.5)

  direct <- fx$prd(times, fx$outerpars, deriv = TRUE)

  b <- dMod2:::.bundle_from_call(fx$conditions, times = times, out = NULL,
                                pars = fx$outerpars, fixed = NULL)
  viaproto <- dMod2:::.evalLegacy(fx$prd, b, deriv = TRUE, deriv2 = FALSE,
                                 env = NULL)

  expect_identical(names(viaproto), fx$conditions)
  for (cn in fx$conditions) {
    expect_equal(unclass(viaproto[[cn]]), unclass(direct[[cn]]), tolerance = 0)
    expect_equal(attr(viaproto[[cn]], "deriv"), attr(direct[[cn]], "deriv"),
                 tolerance = 0)
  }
})


test_that(".evalLegacy respects a restricted condition set", {
  skip_if_not_installed("cppDE")
  skip_on_cran()

  fx <- fx_decay_multicond_compiled()
  times <- seq(0, 5, by = 0.5)

  b <- dMod2:::.bundle_from_call(c("C2", "C4"), times = times, out = NULL,
                                pars = fx$outerpars, fixed = NULL)
  got <- dMod2:::.evalLegacy(fx$prd, b, deriv = TRUE, deriv2 = FALSE, env = NULL)

  expect_identical(names(got), c("C2", "C4"))
  expect_false(any(vapply(got, is.null, logical(1))))

  direct <- fx$prd(times, fx$outerpars, deriv = TRUE,
                   conditions = c("C2", "C4"))
  expect_equal(unclass(got$C4), unclass(direct$C4), tolerance = 0)
})


test_that(".fnNode reports the node kind, and NULL for plain functions", {
  skip_if_not_installed("cppDE")
  skip_on_cran()

  fx <- fx_decay_multicond_compiled()
  expect_identical(dMod2:::.fnNode(fx$prd)$op, "*")     # gfn * xfn * pfn
  expect_identical(dMod2:::.fnNode(fx$pfn)$op, "+")     # P(branch(.)) sums
  expect_identical(dMod2:::.fnNode(fx$xfn)$op, "leaf")
  expect_null(dMod2:::.fnNode(function(x) x))
})


test_that("a composed chain flattens nested sums into one node", {
  skip_if_not_installed("cppDE")
  skip_on_cran()

  # P(branch(.)) builds its parfn with Reduce("+", .), so a+b+c must end up
  # as one node with four parts, not three nested ones.
  fx <- fx_decay_multicond_compiled()
  st <- dMod2:::.fnNode(fx$pfn)
  expect_length(st$parts, length(fx$conditions))
  expect_identical(names(st$owner), fx$conditions)
  expect_true(all(vapply(st$parts, function(p)
    identical(dMod2:::.fnNode(p)$op, "leaf"), logical(1))))
})


test_that("composed mappings carry the metadata their consumers read", {
  skip_if_not_installed("cppDE")
  skip_on_cran()

  fx <- fx_decay_multicond_compiled()
  m <- attr(fx$prd, "mappings")$C2

  # equations come from the prediction side, not the observables
  expect_s3_class(attr(m, "equations"), "eqnvec")
  # parameters are indexed by name, so C2 gets its own scale parameter
  expect_true("s_C2_log" %in% attr(m, "parameters"))
  expect_false("s_C1_log" %in% attr(m, "parameters"))

  eq <- getEquations(fx$prd)
  expect_length(eq, length(fx$conditions))
  expect_false(any(vapply(eq, is.null, logical(1))))
})


test_that("composition rejects class pairs it cannot handle", {
  skip_if_not_installed("cppDE")
  skip_on_cran()

  fx <- fx_decay_multicond_compiled()
  expect_error(fx$xfn * fx$xfn, "no composition defined")
  expect_error(fx$xfn + fx$gfn, "cannot add")
})


## ---- Objective composition ----------------------------------------------

test_that("objfn * parfn keeps parameters, modelname and the NLME handles", {
  skip_if_not_installed("cppDE")
  skip_on_cran()

  fx <- fx_decay_compiled()
  d  <- fx_decay_data()
  obj <- normL2(d, fx$prd_id)
  rep <- obj * fx$pfn_log

  expect_s3_class(rep, "objfn")
  expect_identical(attr(rep, "parameters"), attr(fx$pfn_log, "parameters"))
  expect_false(is.null(attr(rep, "modelname")))
  expect_identical(attr(rep, "data", exact = TRUE), d)

  # value in log coordinates equals the value at the transformed point
  # (solver tolerance only: the two paths seed different sensitivities)
  lp <- fx$outerpars_log
  expect_equal(rep(lp)$value,
               obj(setNames(exp(lp[c("A_log", "k_log")]), c("A", "k")))$value,
               tolerance = 1e-4)
})


test_that("deriv2 and cores survive a sum of objectives", {
  skip_if_not_installed("cppDE")
  skip_on_cran()

  fx  <- fx_decay_compiled()
  d   <- fx_decay_data()
  obj <- normL2(d, fx$prd_log)
  pri <- constraintL2(mu = fx$outerpars_log * 0, sigma = 10)
  p   <- fx$outerpars_log

  # Before the fix match.fnargs dropped deriv2 and the sum silently returned
  # the Gauss-Newton Hessian. Now it reaches the chain, which was built
  # without deriv2 and says so.
  expect_error((obj + pri)(p, deriv = TRUE, deriv2 = TRUE), "deriv2")
  expect_error(obj(p, deriv = TRUE, deriv2 = TRUE), "deriv2")

  # cores reaches the operands the same way
  expect_equal((obj + pri)(p, cores = 2L)$value, (obj + pri)(p)$value,
               tolerance = 0)
})


test_that(".fnWithConditions relabels a leaf so + can dispatch on it", {
  skip_if_not_installed("cppDE")
  skip_on_cran()

  fx <- fx_decay_compiled()
  relabelled <- dMod2:::.fnWithConditions(fx$gfn, c("A", "B"))

  expect_identical(attr(relabelled, "conditions"), c("A", "B"))
  expect_identical(names(attr(relabelled, "mappings")), c("A", "B"))
  expect_identical(dMod2:::.fnNode(relabelled)$condition, c("A", "B"))
  expect_identical(class(relabelled), class(fx$gfn))

  # and it now answers per condition rather than replicating one call
  pred <- fx$xfn(seq(0, 2, 0.5), c(A = 1, k = 0.5), deriv = FALSE)
  out  <- relabelled(pred[[1]], c(A = 1, k = 0.5), deriv = FALSE)
  expect_identical(names(out), c("A", "B"))
})


## ---- Repeated conditions in one request ---------------------------------

test_that(".predictMany sends the same condition several times", {
  skip_if_not_installed("cppDE")
  skip_on_cran()

  fx <- fx_decay_multicond_compiled()
  times <- seq(0, 5, by = 0.5)

  # three parameter sets, all routed to condition C2
  p1 <- fx$outerpars
  p2 <- fx$outerpars; p2["k_log"] <- p2["k_log"] + 0.1
  p3 <- fx$outerpars; p3["k_log"] <- p3["k_log"] - 0.1

  got <- dMod2:::.predictMany(fx$prd, times, list(p1, p2, p3),
                             conditions = rep("C2", 3), deriv = TRUE)
  expect_length(got, 3L)
  expect_identical(names(got), rep("C2", 3))

  for (j in 1:3) {
    ref <- fx$prd(times, list(p1, p2, p3)[[j]], deriv = TRUE,
                  conditions = "C2")$C2
    expect_equal(unclass(got[[j]]), unclass(ref), tolerance = 0)
    expect_equal(attr(got[[j]], "deriv"), attr(ref, "deriv"), tolerance = 0)
  }
  # the three differ, so they were not silently broadcast
  expect_false(isTRUE(all.equal(unclass(got[[1]]), unclass(got[[2]]))))
})


test_that(".predictMany mixes repeated and distinct conditions", {
  skip_if_not_installed("cppDE")
  skip_on_cran()

  fx <- fx_decay_multicond_compiled()
  times <- seq(0, 3, by = 0.5)
  pl <- list(fx$outerpars, fx$outerpars, fx$outerpars)
  conds <- c("C1", "C3", "C1")

  got <- dMod2:::.predictMany(fx$prd, times, pl, conditions = conds,
                             deriv = FALSE)
  expect_identical(names(got), conds)
  expect_equal(unclass(got[[1]]), unclass(got[[3]]), tolerance = 0)
  expect_false(isTRUE(all.equal(unclass(got[[1]]), unclass(got[[2]]))))
})


## ---- fixed handling across composition shapes ----------------------------

# `fixed` is handed from one operand to the next in five different ways across
# the six `*` branches. These pin what each shape actually does, so a future
# unification has to decide the differences deliberately rather than inherit
# them. Tolerances are solver-level: fixing a parameter shrinks the
# sensitivity system, so the integration is not bit-identical.

# Values only: prdframes also carry `parameters`, whose `fixed` marker and
# ordering legitimately differ between the two calls.
.vals <- function(z) { m <- unclass(z); attributes(m) <- list(dim = dim(m)); m }

test_that("with a trafo, fixed leaves the derivative basis", {
  skip_if_not_installed("cppDE")
  skip_on_cran()

  fx <- fx_decay_multicond_compiled()
  times <- seq(0, 4, by = 0.5)
  p <- fx$outerpars

  free <- fx$prd(times, p, deriv = TRUE)
  fixd <- fx$prd(times, p[setdiff(names(p), "k_log")], deriv = TRUE,
                 fixed = p["k_log"])

  for (cn in fx$conditions) {
    expect_equal(.vals(fixd[[cn]]), .vals(free[[cn]]), tolerance = 1e-4)
    keep <- setdiff(dimnames(attr(free[[cn]], "deriv"))[[3]], "k_log")
    expect_identical(dimnames(attr(fixd[[cn]], "deriv"))[[3]], keep)
    expect_equal(attr(fixd[[cn]], "deriv")[, , keep, drop = FALSE],
                 attr(free[[cn]], "deriv")[, , keep, drop = FALSE],
                 tolerance = 1e-4)
  }
})


test_that("without a trafo, fixed does NOT leave the derivative basis", {
  skip_if_not_installed("cppDE")
  skip_on_cran()

  # Xs seeds the identity over its full parameter set and hands the solver
  # fixed = NULL, so a run-time `fixed` only removes the value from `pars`.
  # With a trafo in the chain the parameter disappears through the trafo's
  # Jacobian instead. The two paths therefore disagree; pinned, not endorsed.
  fx <- fx_decay_compiled()
  gx <- fx$gfn * fx$xfn
  times <- seq(0, 4, by = 0.5)
  p <- c(A = 1.0, k = 0.5)

  free <- gx(times, p, deriv = TRUE)[[1]]
  fixd <- gx(times, p["A"], deriv = TRUE, fixed = p["k"])[[1]]

  expect_equal(.vals(fixd), .vals(free), tolerance = 1e-6)
  expect_identical(dimnames(attr(fixd, "deriv"))[[3]],
                   dimnames(attr(free, "deriv"))[[3]])
})


test_that("objfn * parfn reaches the same value with and without fixed", {
  skip_if_not_installed("cppDE")
  skip_on_cran()

  # classes.R passes fixed = NULL into the objective for this branch; the
  # value must still agree with the unfixed evaluation at the same point.
  fx  <- fx_decay_compiled()
  d   <- fx_decay_data()
  obj <- normL2(d, fx$prd_id)
  rep <- obj * fx$pfn_log
  lp  <- fx$outerpars_log

  expect_equal(rep(lp)$value,
               rep(lp[setdiff(names(lp), "k_log")], fixed = lp["k_log"])$value,
               tolerance = 1e-5)
})

test_that("normL2 honours a conditions= restriction", {
  fx <- fx_decay_multicond_compiled()
  d  <- fx_decay_data_multi(parslist = setNames(
    lapply(seq_along(fx$conditions), function(i) c(A = 1.0, k = 0.4 + 0.2 * i)),
    fx$conditions))
  obj <- normL2(d, fx$prd)
  sub <- fx$conditions[1]

  part <- obj(fx$outerpars, deriv = TRUE, conditions = sub)
  ref  <- normL2(d[sub], fx$prd)(fx$outerpars, deriv = TRUE)

  expect_equal(part$value, ref$value, tolerance = 0)
  expect_equal(part$gradient, ref$gradient[names(part$gradient)], tolerance = 0)
  expect_false(isTRUE(all.equal(part$value, obj(fx$outerpars, deriv = TRUE)$value)))
})

test_that("normL2 + constraintL2 restricts both terms to the same conditions", {
  fx <- fx_decay_multicond_compiled()
  d  <- fx_decay_data_multi(parslist = setNames(
    lapply(seq_along(fx$conditions), function(i) c(A = 1.0, k = 0.4 + 0.2 * i)),
    fx$conditions))
  cons <- constraintL2(mu = setNames(rep(0, length(fx$outerpars)),
                                     names(fx$outerpars)), sigma = 10)
  sub <- fx$conditions[1]

  tot <- (normL2(d, fx$prd) + cons)(fx$outerpars, deriv = TRUE, conditions = sub)
  ref <- normL2(d[sub], fx$prd)(fx$outerpars, deriv = TRUE)$value +
           cons(fx$outerpars, deriv = TRUE)$value
  expect_equal(tot$value, ref, tolerance = 1e-10)
})

test_that("the parvec C++ kernel reproduces the R subsetting and concatenation", {
  # `[.parvec` and `c.parvec` delegate to parvec_attach()/parvec_concat(); these
  # are the R bodies they replaced, kept here as the reference.
  sub_R <- function(x, i) {
    out <- .subset(x, i); nms <- names(out)
    deriv <- attr(x, "deriv")
    if (!is.null(deriv)) {
      rn <- rownames(deriv)
      if (!identical(rn, nms)) {
        av <- nms[match(nms, rn, 0L) > 0L]
        deriv <- if (length(av)) deriv[av, , drop = FALSE] else NULL
      }
    }
    d2 <- attr(x, "deriv2")
    if (!is.null(d2)) {
      rn2 <- dimnames(d2)[[1]]
      if (!identical(rn2, nms)) {
        av2 <- nms[match(nms, rn2, 0L) > 0L]
        d2 <- if (length(av2)) d2[av2, , , drop = FALSE] else NULL
      }
    }
    attr(out, "deriv") <- deriv; attr(out, "deriv2") <- d2
    attr(out, "fixed") <- if (!is.null(deriv) && nrow(deriv) < length(nms))
      setdiff(nms, rownames(deriv))
    class(out) <- c("parvec", "numeric"); out
  }
  mk <- function(nm, np, d2 = TRUE, rows = NULL, tag = "t") {
    v  <- setNames(seq_along(nm) + 0.5, nm)
    rr <- if (is.null(rows)) nm else rows
    J  <- matrix(seq_len(length(rr) * np) + 0.25, length(rr), np,
                 dimnames = list(rr, paste0(tag, seq_len(np))))
    H  <- if (d2) array(seq_len(length(rr) * np * np) + 0.125,
                        c(length(rr), np, np),
                        dimnames = list(rr, colnames(J), colnames(J))) else NULL
    as.parvec(v, deriv = J, deriv2 = H)
  }

  x  <- mk(letters[1:6], 4)
  xf <- mk(letters[1:6], 4, rows = letters[1:4])   # two names without a deriv row
  xn <- mk(letters[1:6], 4, d2 = FALSE)
  x0 <- as.parvec(setNames(c(1, 2, 3), c("p", "q", "r")))
  idx <- list(c("a", "c", "e"), letters[1:6], c("f", "a"), 1:3, integer(0),
              c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE), -1L)
  for (v in list(x, xf, xn)) for (i in idx)
    expect_equal(v[i], sub_R(v, i), tolerance = 0)
  expect_equal(x0[c("p", "r")], sub_R(x0, c("p", "r")), tolerance = 0)

  # concatenation, including a block that carries no deriv2 and one with a
  # deriv row missing
  parts <- list(
    list(mk(c("a", "b"), 4), mk(c("c", "d"), 4, tag = "t")),
    list(mk(c("a", "b"), 4), mk(c("c", "d"), 4, d2 = FALSE, tag = "t")),
    list(mk(c("a", "b"), 4), mk(c("c", "d"), 4, rows = "c", tag = "t")),
    list(mk(c("a", "b"), 4), x0))
  for (p in parts) {
    got <- do.call(c, p)
    expect_s3_class(got, "parvec")
    expect_identical(names(got), unlist(lapply(p, names), use.names = FALSE))
    expect_identical(as.numeric(got), unlist(lapply(p, unclass), use.names = FALSE))
    J <- do.call(rbind, Filter(Negate(is.null), lapply(p, attr, "deriv")))
    expect_equal(attr(got, "deriv"), J, tolerance = 0)
    expect_identical(attr(got, "fixed"),
                     if (nrow(J) < length(got)) setdiff(names(got), rownames(J)))
  }
})

test_that("an objective is bit-identical across thread counts", {
  # normL2_kernel scatters per-condition blocks serially in condition order, so
  # the summation order does not follow the thread count.
  fx <- fx_decay_multicond_compiled()
  d <- fx_decay_data_multi(parslist = setNames(
    lapply(seq_along(fx$conditions), function(i) c(A = 1.0, k = 0.4 + 0.2 * i)),
    fx$conditions))
  obj <- normL2(d, fx$prd)
  ref <- obj(fx$outerpars, deriv = TRUE, cores = 1L)
  for (k in c(2L, 3L)) {
    o <- obj(fx$outerpars, deriv = TRUE, cores = k)
    expect_identical(o$value, ref$value)
    expect_equal(o$gradient, ref$gradient, tolerance = 0)
    expect_equal(o$hessian,  ref$hessian,  tolerance = 0)
  }
})
