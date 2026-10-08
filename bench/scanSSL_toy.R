## Spike-and-slab lasso path of scanL1() on the toy decay of Hauber, Rosenblatt
## and Timmer (2023, Fig 2): dx/dt = -p x, x(0) = 1, three cell types with
## log10 p = -1.5, -1.3, -1.2. Types 2 and 3 share a mutation that
## differs slightly; the reference value of p is known.
library(dMod2)
outdir <- tempdir()

x <- Xs(odemodel(eqnvec(x = "-p * x"), modelname = "toySSL", compile = FALSE,
                 outdir = outdir), compile = FALSE)
g <- Y(eqnvec(y = "x"), x, modelname = "toySSL_obs", compile = FALSE,
       attach.input = FALSE, outdir = outdir)
types <- c("c1", "c2", "c3")
trafo <- setNames(lapply(types, function(ct) eqnvec(
  x = "1", p = if (ct == "c1") "10^lp" else paste0("10^(lp + r_", ct, ")"))), types)
p <- P(trafo, modelname = "toySSL_p", compile = FALSE, outdir = outdir)
compile(g, x, p, output = "toySSL")
prd <- g * x * p

truth <- c(lp = -1.5, r_c2 = 0.2, r_c3 = 0.3)
sigma <- 10^-1.3
times <- seq(0, 100, by = 10)
block <- list(r = list(pars = c("r_c2", "r_c3"), anchor = 0))

## One dataset -----------------------------------------------------------------
set.seed(3)
sim <- wide2long(prd(times, truth, deriv = FALSE))
sim$value <- sim$value + rnorm(nrow(sim), 0, sigma)
sim$sigma <- sigma
obj <- normL2(as.datalist(sim), prd)
plot(prd(seq(0, 100, 1), truth), as.datalist(sim))

ssl <- scanL1(obj, c(r_c2 = 0, r_c3 = 0), groups = block, fixed = truth["lp"],
              lambda = 10^seq(0, 4, length.out = 25), fits = 10, pathFits = 1,
              ssl = list(lambda1 = 1), select = "plateau")
ssl
plot(ssl, type = "path")
plot(ssl, type = "clusters")
plot(ssl, type = "inclusion")

q08 <- scanL1(obj, c(r_c2 = 0, r_c3 = 0), groups = block, fixed = truth["lp"],
              lambda = 10^seq(0, 3, length.out = 13), fits = 10, q = 0.8)
q08
plot(q08, type = "path")

## Many datasets: how often each method lands on {c2, c3} ----------------------
res <- do.call(rbind, lapply(1:50, function(seed) {
  set.seed(seed)
  sim <- wide2long(prd(times, truth, deriv = FALSE))
  sim$value <- sim$value + rnorm(nrow(sim), 0, sigma)
  sim$sigma <- sigma
  obj <- normL2(as.datalist(sim), prd)
  s <- scanL1(obj, c(r_c2 = 0, r_c3 = 0), groups = block, fixed = truth["lp"],
              lambda = 10^seq(0, 4, length.out = 25), fits = 10, pathFits = 1,
              ssl = list(lambda1 = 1), select = "plateau")
  q <- scanL1(obj, c(r_c2 = 0, r_c3 = 0), groups = block, fixed = truth["lp"],
              lambda = 10^seq(0, 3, length.out = 13), fits = 10, q = 0.8)
  data.frame(seed = seed, ssl = s$selected, lq = q$selected,
             pSsl = s$refits$p[s$refits$key == s$selected])
}))
table(ssl = res$ssl, lq = res$lq)
summary(res$pSsl)
