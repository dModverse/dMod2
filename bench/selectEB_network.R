## selectEB on a reaction network: A, B, C, D with eight candidate reactions,
## three of them present (the chain A -> B -> C -> D), all species observed.
## Compared with the L1 path and the spike-and-slab path of scanL1.
library(dMod2)
outdir <- tempdir()

f <- eqnvec(A = "-k1*A - k4*A - k6*A + k7*C", B = "k1*A - k2*B - k5*B + k8*D",
            C = "k2*B + k4*A - k3*C - k7*C", D = "k3*C + k5*B + k6*A - k8*D")
x <- Xs(odemodel(f, modelname = "ebnet_x", compile = FALSE, outdir = outdir), compile = FALSE)
g <- Y(eqnvec(oA = "A", oB = "B", oC = "C", oD = "D"), f, modelname = "ebnet_g",
       compile = FALSE, attachInput = FALSE, outdir = outdir)
ks <- paste0("k", 1:8)
trafo <- gateL1(eqnvec(A = "1", B = "0", C = "0", D = "0",
                       setNames(paste0("10^log10_", ks), ks)), ks)
p <- P(trafo, condition = "C1", modelname = "ebnet_p", compile = FALSE, outdir = outdir)
compile(g, x, p, output = "ebnet_model")
prd <- g * x * p

truth <- c(setNames(rep(-1, 8), paste0("log10_", ks)), setNames(rep(0, 8), paste0("s_", ks)))
truth[c("log10_k1", "log10_k2", "log10_k3")] <- log10(c(0.5, 0.3, 0.2))
truth[c("s_k1", "s_k2", "s_k3")] <- 1
center <- setNames(rep(-0.5, 8), paste0("log10_", ks))
zero <- attr(trafo, "gates")
key <- "-s_k4 -s_k5 -s_k6 -s_k7 -s_k8"
sigma <- 0.1

## One dataset -----------------------------------------------------------------
set.seed(1)
sim <- wide2long(prd(seq(0.5, 15, 0.5), truth, deriv = FALSE))
sim$value <- sim$value + rnorm(nrow(sim), 0, sigma)
sim$sigma <- sigma
obj <- normL2(as.datalist(sim), prd)
plot(prd(seq(0, 15, 0.1), truth), as.datalist(sim))

eb <- selectEB(obj, center, zero = zero, fits = 10, control = list(hits = 3))
eb
eb$alternatives
plot(eb, type = "terms")
plot(eb, type = "waterfall")

ebA <- selectEB(obj, center, zero = zero, fits = 10, rule = "alpha")
ebA

l1 <- scanL1(obj, center, zero = zero, lambda = 10^seq(-1, 4, length.out = 11), fits = 4)
l1
plot(l1)

ssl <- scanL1(obj, center, zero = zero, lambda = 10^seq(0, 4, length.out = 13), fits = 4,
              pathFits = 1, ssl = list(lambda1 = 1), select = "plateau")
ssl

## Many datasets: how often each method selects the chain -----------------------
res <- do.call(rbind, lapply(1:20, function(seed) {
  set.seed(seed)
  sim <- wide2long(prd(seq(0.5, 15, 0.5), truth, deriv = FALSE))
  sim$value <- sim$value + rnorm(nrow(sim), 0, sigma)
  sim$sigma <- sigma
  obj <- normL2(as.datalist(sim), prd)
  t0 <- proc.time()[3]
  e <- selectEB(obj, center, zero = zero, fits = 6)
  t1 <- proc.time()[3]
  s <- scanL1(obj, center, zero = zero, lambda = 10^seq(0, 4, length.out = 13), fits = 4,
              pathFits = 1, ssl = list(lambda1 = 1), select = "plateau")
  t2 <- proc.time()[3]
  data.frame(seed = seed, eb = e$selected == key, ebBest = e$best$key == key,
             ssl = s$selected == key, g = e$best$g, tEB = t1 - t0, tSSL = t2 - t1)
}))
colMeans(res[, c("eb", "ebBest", "ssl")])
c(eb = min(res$tEB), ssl = min(res$tSSL))
