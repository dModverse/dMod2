## Cell-line specific parameters: the L_q scan of Hauber, Rosenblatt and Timmer (2023, PLoS
## Comput Biol 19, e1010867) against em() with Dirichlet-process random effects, dirichletProcess(), on simulated
## data of the Becker Epo receptor model with six lines. Truth: kon shared by all lines,
## ke {A, B, C} and {D, E, F}, kt {A, B}, {C, D} and {E, F}. Per dataset both methods report a
## partition of the lines per parameter; counted is the exact recovery and the time.
library(dMod2)
options(dMod2.outdir = file.path(tempdir(), "becker"))
dir.create(getOption("dMod2.outdir"), showWarnings = FALSE)

nsim     <- 10                              # simulated datasets
fitsFull <- 10                              # starts of the full model
fitsPath <- 4                               # starts per lambda and per refit
grid     <- 10^seq(-1, 3, length.out = 25)  # penalty strengths of the scan
q        <- 0.8                             # exponent of the L_q penalty
alpha    <- 0.05                            # level of the likelihood-ratio test

## Model ----------------------------------------------------------------------------------
reactions <- NULL
reactions <- addReaction(reactions, "Epo + EpoR", "Epo_EpoR", "Epo * EpoR * kon / init_Epo", "binding")
reactions <- addReaction(reactions, "Epo_EpoR", "Epo + EpoR", "Epo_EpoR * koff", "unbinding")
reactions <- addReaction(reactions, "", "EpoR", "4 * init_Epo * init_EpoR_rel * kt", "synthesis")
reactions <- addReaction(reactions, "EpoR", "", "EpoR * kt", "turnover")
reactions <- addReaction(reactions, "Epo_EpoR", "Epo_EpoR_i", "Epo_EpoR * ke", "internalisation")
reactions <- addReaction(reactions, "Epo_EpoR_i", "Epo + EpoR", "Epo_EpoR_i * kex", "recycling")
reactions <- addReaction(reactions, "Epo_EpoR_i", "dEpo_i", "Epo_EpoR_i * kdi", "intra degr.")
reactions <- addReaction(reactions, "Epo_EpoR_i", "dEpo_e", "Epo_EpoR_i * kde", "extra degr.")
x <- Xs(odemodel(reactions, modelname = "becker", compile = FALSE, deriv2 = TRUE))
g <- Y(eqnvec(y_ext = "log(Epo + dEpo_e + 1)", y_mem = "log(Epo_EpoR + 1)",
              y_int = "log(Epo_EpoR_i + dEpo_i + 1)"),
       x, modelname = "becker_obs", compile = FALSE, attachInput = FALSE, deriv2 = TRUE)
e <- Y(eqnvec(y_ext = "sigma", y_mem = "sigma", y_int = "sigma"), g, modelname = "becker_err",
       compile = FALSE, attachInput = FALSE, deriv2 = TRUE)

dose  <- 1347.49
lines <- c("A", "B", "C", "D", "E", "F")
cand  <- c("kon", "ke", "kt")
rates <- c("kon", "koff", "kt", "ke", "kex", "kdi", "kde", "init_EpoR_rel", "sigma")
tab   <- data.frame(line = lines, row.names = lines)
base  <- function() {
  tr <- branch(define(NULL, "x <- x", x = c(rates, "init_Epo", "Epo", "EpoR", "Epo_EpoR", "Epo_EpoR_i",
                                            "dEpo_i", "dEpo_e")), table = tab)
  tr <- insert(tr, "x <- d", x = c("init_Epo", "Epo"), d = dose)
  tr <- insert(tr, "EpoR <- 4 * d * init_EpoR_rel", d = dose)
  tr <- insert(tr, "x <- 0", x = c("Epo_EpoR", "Epo_EpoR_i", "dEpo_i", "dEpo_e"))
  insert(tr, "x <- exp(x)", x = rates)
}

# Hauber et al.: line A is the reference, every other line has a log fold change r_<par>_<line>
trH <- base()
for (l in lines[-1]) for (pp in cand)
  trH <- insert(trH, paste0(pp, " <- ", pp, " + r_", pp, "_", l), conditionMatch = paste0("^", l, "$"))
pH <- P(trH, method = "explicit", modelname = "becker_h", compile = FALSE, deriv2 = TRUE)

# EM: per line and candidate a deviation from the population value, shared by clusters of lines
pS <- P(insert(base(), x ~ 1 + (1 | gr(line, dist = dirichletProcess())), x = cand), method = "explicit",
        modelname = "becker_s", compile = FALSE, deriv2 = TRUE)
compile(g, x, e, pH, pS, output = "becker")

## Truth and the parameters held fixed ---------------------------------------------------
st <- c(kon = log(0.1513), koff = log(0.08069), kt = log(0.01601), ke = log(0.05546), kex = log(0.000578),
        kdi = log(0.001273), kde = log(0.01199), init_EpoR_rel = log(0.09205), sigma = log(0.03))
truth <- rbind(kon = c(A = 0, B = 0, C = 0, D = 0, E = 0, F = 0),
               ke  = c(A = 0, B = 0, C = 0, D = -1.6, E = -1.6, F = -1.6),
               kt  = c(A = 0, B = 0, C = -0.7, D = -0.7, E = -1.4, F = -1.4))
fixed <- st[c("koff", "kex", "kdi", "kde", "init_EpoR_rel")]
times <- c(3, 6, 12, 20, 30, 45, 60, 90, 120, 180, 240, 300)
r_true <- unlist(lapply(cand, function(pp) setNames(truth[pp, -1], paste0("r_", pp, "_", lines[-1]))))
blocks <- lapply(cand, function(pp) paste0("r_", pp, "_", lines[-1]))

# The partition of the lines of one parameter as labels in the order of first appearance.
partition <- function(v, tol = 1e-6) {
  lab <- integer(length(v)); k <- 0L
  for (i in seq_along(v)) if (!lab[i]) { k <- k + 1L; lab[abs(v - v[i]) < tol & !lab] <- k }
  setNames(lab, names(v))
}
target <- lapply(setNames(nm = cand), function(pp) partition(truth[pp, ]))

# The fold changes of line A (0) to F of a Hauber parameter vector, and the L_q penalty of all
# pairwise differences, line A included.
foldChanges <- function(par, pp) c(A = 0, setNames(par[paste0("r_", pp, "_", lines[-1])], lines[-1]))
penaltyLq <- function(par) sum(vapply(cand, function(pp) {
  d <- abs(outer(foldChanges(par, pp), foldChanges(par, pp), `-`))
  sum(d[upper.tri(d)]^q)
}, 0))

# The objective with the fold changes of a structure tied: one free value per cluster away
# from line A, the others 0; value, gradient and Hessian by the chain rule.
tied <- function(obj, structure, free) {
  map <- matrix(0, length(r_true), 0, dimnames = list(names(r_true), NULL))
  for (pp in cand) {
    lab <- structure[[pp]]
    for (k in setdiff(unique(lab), lab[["A"]])) {
      col <- setNames(numeric(length(r_true)), names(r_true))
      col[paste0("r_", pp, "_", names(lab)[lab == k])] <- 1
      map <- cbind(map, col)
    }
  }
  colnames(map) <- if (ncol(map)) paste0("c", seq_len(ncol(map))) else character(0)
  function(pars, fixed = NULL, ...) {
    full <- c(pars[free], drop(map %*% pars[colnames(map)]))
    o <- obj(full, fixed = fixed, ...)
    J <- rbind(cbind(diag(length(free)), matrix(0, length(free), ncol(map))),
               cbind(matrix(0, nrow(map), length(free)), map))
    dimnames(J) <- list(c(free, rownames(map)), c(free, colnames(map)))
    J <- J[names(o$gradient), , drop = FALSE]
    objlist(value = o$value, gradient = drop(crossprod(J, o$gradient)), hessian = crossprod(J, o$hessian %*% J))
  }
}

## Simulation and fits -------------------------------------------------------------------
results <- NULL
for (s in seq_len(nsim)) {
  set.seed(s)
  sim <- wide2long((g * x * pH)(times, c(st, r_true), deriv = FALSE))
  sim$value <- sim$value + rnorm(nrow(sim), 0, 0.03)
  sim$sigma <- NA_real_
  data <- as.datalist(sim)
  objH <- normL2(data, g * x * pH, errmodel = e)
  shared <- st[c("kon", "ke", "kt", "sigma")]

  # Hauber et al.: full fit, L_q path upward from the previous optimum and the sparse point,
  # L_q by reweighted L1 (three rounds), unpenalised refit of every structure on the path,
  # the structure of the largest lambda the likelihood-ratio test against the full fit keeps
  tH <- system.time({
    center <- c(shared, r_true * 0)
    full <- as.parframe(mstrust(objH, center, fits = fitsFull, sd = 0.5, fixed = fixed, startFromCenter = TRUE))
    best <- unlist(full[1, names(center)])
    path <- list(); prev <- best
    for (lam in grid) {
      fits <- list()
      for (start in list(prev, c(prev[names(shared)], r_true * 0))) {
        par <- start
        for (round in 1:3) {
          fuse <- lapply(blocks, function(b) {
            v <- c(par[b], 0)
            list(pars = b, lambda = lam, anchor = 0,
                 weights = q * (abs(outer(v, v, `-`)) + 0.01)^(q - 1))
          })
          fit <- trustL1(objH, par, fuse = fuse, rinit = 0.1, rmax = 10, iterlim = 200, fixed = fixed)
          par <- fit$argument
        }
        fits[[length(fits) + 1L]] <- fit
      }
      lq <- vapply(fits, function(f) objH(f$argument, fixed = fixed, deriv = FALSE)$value + lam * penaltyLq(f$argument), 0)
      fit <- fits[[which.min(lq)]]
      prev <- fit$argument
      path[[length(path) + 1L]] <- lapply(setNames(nm = cand), function(pp) partition(foldChanges(prev, pp)))
    }
    keys <- vapply(path, function(p) paste(unlist(p), collapse = ""), "")
    refits <- lapply(which(!duplicated(keys)), function(i) {
      structure <- path[[i]]
      f <- tied(objH, structure, names(shared))
      k <- sum(vapply(structure, function(l) length(setdiff(unique(l), l[["A"]])), 0))
      start <- c(best[names(shared)], setNames(rep(0, k), if (k) paste0("c", seq_len(k))))
      fit <- trust(f, start, rinit = 0.1, rmax = 10, iterlim = 200, fixed = fixed)
      for (j in seq_len(fitsPath - 1L)) {
        alt <- trust(f, start + rnorm(length(start), 0, 0.3), rinit = 0.1, rmax = 10, iterlim = 200, fixed = fixed)
        if (alt$value < fit$value) fit <- alt
      }
      list(lambda = max(grid[keys == keys[i]]), structure = structure, value = fit$value, k = k)
    })
    df <- length(r_true) - vapply(refits, `[[`, 0, "k")
    keep <- vapply(refits, `[[`, 0, "value") - full$value[1] <= stats::qchisq(1 - alpha, pmax(df, 1)) | df == 0
    hauber <- if (any(keep)) refits[keep][[which.max(vapply(refits[keep], `[[`, 0, "lambda"))]]$structure else
      lapply(setNames(nm = cand), function(pp) partition(foldChanges(best, pp)))
  })

  # EM with Dirichlet-process random effects, the concentration and the slab estimated
  objS <- normL2(data, g * x * pS, errmodel = e)
  tS <- system.time(fitS <- em(objS, randomEffects(pS), init = shared, fixed = fixed, verbose = FALSE))
  clustersOf <- function(fit) lapply(setNames(nm = cand), function(pp) {
    r <- fit$random[fit$random$eta == pp, ]
    partition(setNames(as.numeric(r$cluster), r$group)[lines])
  })

  for (m in list(list("Hauber L_q + LRT", hauber, tH), list("EM dirichletProcess()", clustersOf(fitS), tS))) {
    ok <- vapply(cand, function(pp) identical(unname(m[[2]][[pp]]), unname(target[[pp]])), TRUE)
    results <- rbind(results, data.frame(dataset = s, method = m[[1]], seconds = m[[3]][["elapsed"]],
                                         t(ok), exact = all(ok)))
  }
}

## Summary -------------------------------------------------------------------------------
results
aggregate(cbind(exact, kon, ke, kt) ~ method, results, sum)
aggregate(seconds ~ method, results, median)
