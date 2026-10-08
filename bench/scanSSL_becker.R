## Cluster identification across cell types with the spike-and-slab lasso path
## of scanL1(), on the simulated Becker Epo receptor study of Hauber,
## Rosenblatt and Timmer (2023, Fig 3; d2d Becker_Science2010/
## Setup_Regularization.m). Five cell lines, CL1 is the reference. True log10
## fold changes: init_Epo CL2, CL3 +0.4; kde CL2, CL4 -0.5; kon CL3, CL4 -0.5,
## CL5 +0.5; all other parameters shared.
## One multistart of the full model, then one warm-started chain over the spike
## strength; no multistart per lambda.
library(dMod2)
outdir <- tempdir()
set.seed(1)

## Model, log10 parameters as in d2d ------------------------------------------
reactions <- NULL
reactions <- addReaction(reactions, "Epo + EpoR", "Epo_EpoR",
                         "kon / init_Epo * Epo * EpoR", "binding")
reactions <- addReaction(reactions, "Epo_EpoR", "Epo + EpoR", "koff * Epo_EpoR", "unbinding")
reactions <- addReaction(reactions, "", "EpoR", "kt * 4 * init_Epo * init_EpoR_rel", "synthesis")
reactions <- addReaction(reactions, "EpoR", "", "kt * EpoR", "turnover")
reactions <- addReaction(reactions, "Epo_EpoR", "Epo_EpoR_i", "ke * Epo_EpoR", "internalisation")
reactions <- addReaction(reactions, "Epo_EpoR_i", "Epo + EpoR", "kex * Epo_EpoR_i", "recycling")
reactions <- addReaction(reactions, "Epo_EpoR_i", "dEpo_i", "kdi * Epo_EpoR_i", "intra degr.")
reactions <- addReaction(reactions, "Epo_EpoR_i", "dEpo_e", "kde * Epo_EpoR_i", "extra degr.")

## Binding curve: Michaelis-Menten in log10 free Epo, one observable per level.
ef  <- c(0.71, 1.76, 2.24, 2.59, 3.26, 3.6)
efN <- c(1, 3, 1, 3, 3, 1)
bind <- setNames(sprintf("log10(init_Epo * init_EpoR_rel * 10^%g / (koff / kon + 10^%g))",
                         ef, ef), paste0("y_bind", seq_along(ef)))
obs <- as.eqnvec(c(y_ext = "log10(offset + scale * (Epo + dEpo_e))",
                   y_mem = "log10(offset + scale * Epo_EpoR)",
                   y_int = "log10(offset + scale * (Epo_EpoR_i + dEpo_i))", bind))

x <- Xs(odemodel(reactions, modelname = "beckerSSL", compile = FALSE, outdir = outdir),
        compile = FALSE)
g <- Y(obs, x, modelname = "beckerSSL_obs", compile = FALSE, attach.input = FALSE,
       outdir = outdir)

lines <- paste0("CL", 1:5)
cand  <- c("init_Epo", "init_EpoR_rel", "kde", "kdi", "ke", "koff", "kon", "kt")
fc    <- function(par, l) if (l == "CL1") "" else paste0(" + r_", par, "_", l)
trafo <- setNames(lapply(lines, function(l) {
  tr <- setNames(paste0("10^(l", cand, vapply(cand, fc, "", l = l), ")"), cand)
  as.eqnvec(c(tr, kex = "10^lkex", offset = "10^loffset", scale = "10^lscale",
           Epo = tr[["init_Epo"]],
           EpoR = paste0("4 * ", tr[["init_Epo"]], " * ", tr[["init_EpoR_rel"]]),
           Epo_EpoR = "0", Epo_EpoR_i = "0", dEpo_i = "0", dEpo_e = "0"))
}), lines)
p <- P(trafo, modelname = "beckerSSL_p", compile = FALSE, outdir = outdir)
compile(g, x, p, output = "beckerSSL")
prd <- g * x * p

## Truth and simulated data ---------------------------------------------------
base <- c(linit_Epo = 3.129526, linit_EpoR_rel = -1.035969, lkde = -1.921247,
          lkdi = -2.895286, lke = -1.256020, lkex = -3.238069, lkoff = -1.093201,
          lkon = -0.8200363, lkt = -1.795540, loffset = -5, lscale = -0.00884944)
truth <- matrix(0, length(cand), 4, dimnames = list(cand, lines[-1]))
truth["init_Epo", c("CL2", "CL3")] <- 0.4
truth["kde", c("CL2", "CL4")]      <- -0.5
truth["kon", c("CL3", "CL4", "CL5")] <- c(-0.5, -0.5, 0.5)
r_true <- setNames(as.vector(t(truth)),
                   paste0("r_", rep(cand, each = 4), "_", lines[-1]))
sdlog <- c(y_ext = -2.074558, y_mem = -1.320326, y_int = -1.256882, y_bind = -1.400701)

times <- c(0.82, 5.82, 20.82, 60.82, 120.82, 180.82, 240.82, 300.82)
pred  <- prd(times, c(base, r_true), deriv = FALSE)
## Time courses in triplicate, the binding curve at the first time point
d  <- wide2long(pred)
tc <- d[d$name %in% c("y_ext", "y_mem", "y_int"), ]
tc <- tc[rep(seq_len(nrow(tc)), 3), ]
bd <- d[grepl("^y_bind", d$name) & d$time == times[1], ]
bd <- bd[rep(seq_len(nrow(bd)), efN[as.integer(sub("y_bind", "", bd$name))]), ]
d  <- rbind(tc, bd)
key <- ifelse(grepl("^y_bind", d$name), "y_bind", as.character(d$name))
d$sigma <- 10^sdlog[key]
d$value <- d$value + rnorm(nrow(d), 0, d$sigma)
data <- as.datalist(d)
plot(prd(seq(0, 310, 5), c(base, r_true)), data)

obj    <- normL2(data, prd)
fixed  <- base[c("lkex", "loffset")]
center <- c(base[setdiff(names(base), names(fixed))] * 0 - 1, r_true * 0)
blocks <- setNames(lapply(cand, function(pp)
  list(pars = paste0("r_", pp, "_", lines[-1]), anchor = 0)), cand)
grid <- 10^seq(0, 4, length.out = 25)

## Spike-and-slab path: one multistart of the full model, then one chain -----
t_ssl <- system.time(
  ssl <- scanL1(obj, center, groups = blocks, fixed = fixed, lambda = grid,
                fits = 100, pathFits = 1, cores = 10, sd = 1,
                ssl = list(lambda1 = 1), select = "plateau"))
t_ssl
ssl
plot(ssl, type = "waterfall")
plot(ssl, type = "clusters")
plot(ssl, type = "inclusion")
plot(ssl, type = "test")
ssl$structure[[ssl$selected]]$groups

## Truth as a structure key, for comparison.
truthKey <- dMod2:::.l1Structure(r_true, NULL, NULL, blocks)$key
identical(ssl$selected, truthKey)
ssl$path[, c("lambda", "key")]

## Reference: the L0.8 penalty of the paper with LRT, multistart per lambda.
t_q <- system.time(
  q08 <- scanL1(obj, center, groups = blocks, fixed = fixed, lambda = grid,
                fits = 100, pathFits = 10, cores = 10, sd = 1, q = 0.8))
t_q
q08
plot(q08, type = "clusters")
identical(q08$selected, truthKey)
