## The decay toy of Hauber, Rosenblatt and Timmer (2023, Fig 2): dx/dt = -p x, x(0) = 1, three
## cell types with log10 p = -1.5 (c1, known), -1.3 and -1.2. The L_q scan with the
## likelihood-ratio test against em() with Dirichlet-process random effects, dirichletProcess(); per dataset the
## partition of the types each method selects and the time.
library(dMod2)
options(dMod2.outdir = file.path(tempdir(), "toySel"))
dir.create(getOption("dMod2.outdir"), showWarnings = FALSE)

nsim     <- 50
grid     <- 10^seq(-1, 3, length.out = 25)
q        <- 0.8
alpha    <- 0.05
fitsPath <- 4

## Model ----------------------------------------------------------------------------------
x <- Xs(odemodel(eqnvec(x = "-p * x"), modelname = "toySel", compile = FALSE, derivMode = "forward-forward"))
g <- Y(eqnvec(y = "x"), x, modelname = "toySel_obs", compile = FALSE, attachInput = FALSE, derivMode = "forward-forward")
types <- c("c1", "c2", "c3")
tab <- data.frame(type = types, row.names = types)
base <- function() insert(branch(define(NULL, "x <- x", x = c("x", "p")), table = tab), "x <- 1")

# Hauber et al.: fold changes r_c2, r_c3 to the reference c1
trH <- insert(base(), "p <- 10^lp")
for (ct in types[-1]) trH <- insert(trH, paste0("lp <- lp + r_", ct), conditionMatch = paste0("^", ct, "$"))
pH <- P(trH, method = "explicit", modelname = "toySel_h", compile = FALSE, derivMode = "forward-forward")
# EM: a deviation of lp per type, shared by clusters of types
trS <- insert(base(), "p <- 10^lp")
pS <- P(insert(trS, lp ~ 1 + (1 | gr(type, dist = dirichletProcess()))), method = "explicit",
        modelname = "toySel_s", compile = FALSE, derivMode = "forward-forward")
compile(g, x, pH, pS, output = "toySel")

truth <- c(lp = -1.5, r_c2 = 0.2, r_c3 = 0.3)
sigma <- 10^-1.3
times <- seq(0, 100, by = 10)
block <- c("r_c2", "r_c3")
fixed <- truth["lp"]

# The partition of the types as labels in the order of first appearance, as one string.
partition <- function(v, tol = 1e-6) {
  lab <- integer(length(v)); k <- 0L
  for (i in seq_along(v)) if (!lab[i]) { k <- k + 1L; lab[abs(v - v[i]) < tol & !lab] <- k }
  paste(lab, collapse = "")
}
foldChanges <- function(par) c(c1 = 0, par[block])
penaltyLq <- function(par) { d <- abs(outer(foldChanges(par), foldChanges(par), `-`)); sum(d[upper.tri(d)]^q) }

## Simulation and fits -------------------------------------------------------------------
results <- NULL
for (s in seq_len(nsim)) {
  set.seed(s)
  sim <- wide2long((g * x * pH)(times, truth, deriv = FALSE))
  sim$value <- sim$value + rnorm(nrow(sim), 0, sigma)
  sim$sigma <- sigma
  data <- as.datalist(sim)
  objH <- normL2(data, g * x * pH)

  tH <- system.time({
    full <- trust(objH, c(r_c2 = 0, r_c3 = 0), rinit = 0.1, rmax = 10, iterlim = 200, fixed = fixed)
    path <- character(0); lams <- numeric(0); prev <- full$argument; refits <- list()
    for (lam in grid) {
      fits <- list()
      for (start in list(prev, c(r_c2 = 0, r_c3 = 0))) {
        par <- start
        for (round in 1:3) {
          v <- c(par[block], 0)
          fit <- trustL1(objH, par, fuse = list(list(pars = block, lambda = lam, anchor = 0,
                                                    weights = q * (abs(outer(v, v, `-`)) + 0.01)^(q - 1))),
                         rinit = 0.1, rmax = 10, iterlim = 200, fixed = fixed)
          par <- fit$argument
        }
        fits[[length(fits) + 1L]] <- fit
      }
      lq <- vapply(fits, function(f) objH(f$argument, fixed = fixed, deriv = FALSE)$value + lam * penaltyLq(f$argument), 0)
      prev <- fits[[which.min(lq)]]$argument
      path <- c(path, partition(foldChanges(prev))); lams <- c(lams, lam)
    }
    # refit per structure: one free value per cluster away from c1
    for (key in unique(path)) {
      lab <- as.integer(strsplit(key, "")[[1]])
      free <- setdiff(unique(lab), lab[1])
      f <- function(pars, fixed = NULL, ...) {
        r <- setNames(numeric(2), block)
        for (k in seq_along(free)) r[lab[-1] == free[k]] <- pars[[k]]
        o <- objH(r, fixed = fixed, ...)
        J <- vapply(seq_along(free), function(k) as.numeric(lab[-1] == free[k]), numeric(2))
        J <- matrix(J, 2, length(free), dimnames = list(block, paste0("c", seq_along(free))))
        objlist(value = o$value, gradient = drop(crossprod(J, o$gradient[block])),
                hessian = crossprod(J, o$hessian[block, block] %*% J))
      }
      value <- if (!length(free)) objH(c(r_c2 = 0, r_c3 = 0), fixed = fixed)$value else {
        start <- setNames(rep(0.2, length(free)), paste0("c", seq_along(free)))
        best <- trust(f, start, rinit = 0.1, rmax = 10, iterlim = 200, fixed = fixed)
        for (j in seq_len(fitsPath - 1L)) {
          alt <- trust(f, start + rnorm(length(start), 0, 0.3), rinit = 0.1, rmax = 10, iterlim = 200, fixed = fixed)
          if (alt$value < best$value) best <- alt
        }
        best$value
      }
      refits[[key]] <- c(lambda = max(lams[path == key]), value = value, df = 2 - length(free))
    }
    rf <- do.call(rbind, refits)
    keep <- rf[, "value"] - full$value <= stats::qchisq(1 - alpha, pmax(rf[, "df"], 1)) | rf[, "df"] == 0
    hauber <- if (any(keep)) rownames(rf)[keep][which.max(rf[keep, "lambda"])] else partition(foldChanges(full$argument))
  })

  objS <- normL2(data, g * x * pS)
  tS <- system.time(fitS <- em(objS, randomEffects(pS), fixed = fixed, verbose = FALSE))
  cl <- setNames(as.numeric(fitS$random$cluster), fitS$random$group)[types]
  results <- rbind(results,
                   data.frame(dataset = s, method = "Hauber L_q + LRT", structure = hauber, seconds = tH[["elapsed"]]),
                   data.frame(dataset = s, method = "EM dirichletProcess()", structure = partition(cl),
                              seconds = tS[["elapsed"]]))
}

## Summary -------------------------------------------------------------------------------
# "123": all types differ (the simulated values); "122": c2 and c3 share one value; "112" and
# "111" merge c2 or both with the reference
table(results$method, results$structure)
aggregate(seconds ~ method, results, median)
