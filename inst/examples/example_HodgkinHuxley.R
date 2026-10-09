# Hodgkin-Huxley squid axon with piecewise() on time (current pulse) and on V
# (removable 0/0 gate rates); fits the three maximal conductances.
# piecewise() on a state is safe only where the right-hand side is continuous.

library(dMod2)

.outdir <- file.path(tempdir(), "hodgkinHuxley")
if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)


## Model -------------------------------------------------------------------
alpha_m <- "piecewise(1, V == -40, 0.1*(V + 40)/(1 - exp(-(V + 40)/10)))"
beta_m  <- "4*exp(-(V + 65)/18)"
alpha_h <- "0.07*exp(-(V + 65)/20)"
beta_h  <- "1/(1 + exp(-(V + 35)/10))"
alpha_n <- "piecewise(0.1, V == -55, 0.01*(V + 55)/(1 - exp(-(V + 55)/10)))"
beta_n  <- "0.125*exp(-(V + 65)/80)"

I_Na   <- "g_Na*m^3*h*(V - E_Na)"
I_K    <- "g_K*n^4*(V - E_K)"
I_L    <- "g_L*(V - E_L)"
I_stim <- "piecewise(I0, time >= t_on && time < t_on + dur, 0)"

f <- eqnvec(
  V = paste0("(", I_stim, " - ", I_Na, " - ", I_K, " - ", I_L, ")/C_m"),
  m = paste0(alpha_m, "*(1 - m) - ", beta_m, "*m"),
  h = paste0(alpha_h, "*(1 - h) - ", beta_h, "*h"),
  n = paste0(alpha_n, "*(1 - n) - ", beta_n, "*n")
)

model <- odemodel(f, modelname = "hh_ode", compile = FALSE, outdir = .outdir)
x <- Xs(model)


## Parameters --------------------------------------------------------------
# The gates start at their steady state at -65 mV
gateAtRest <- function(alpha, beta, V = -65) {
  env <- list(V = V, piecewise = function(value, condition, otherwise)
    if (condition) value else otherwise)
  a <- eval(parse(text = alpha), env); b <- eval(parse(text = beta), env)
  a/(a + b)
}

trafo <- eqnvec(
  V = "-65",
  m = gateAtRest(alpha_m, beta_m),
  h = gateAtRest(alpha_h, beta_h),
  n = gateAtRest(alpha_n, beta_n),
  g_Na = "10^lg_Na", g_K = "10^lg_K", g_L = "10^lg_L",
  E_Na = "50", E_K = "-77", E_L = "-54.387", C_m = "1",
  t_on = "2", dur = "1", I0 = "I0"
)

# A weak pulse stays below threshold, a strong one fires a spike
p <- P(list(weak = insert(trafo, "I0 ~ 3"), strong = insert(trafo, "I0 ~ 20")),
       modelname = "hh_trafo", compile = FALSE, outdir = .outdir)

compile(x, p, output = "hodgkinHuxley")

pars_true <- c(lg_Na = log10(120), lg_K = log10(36), lg_L = log10(0.3))


## Simulated data ----------------------------------------------------------
times <- seq(0, 20, by = 0.02)
plot((x*p)(times, pars_true))

set.seed(1)
data <- as.data.frame((x*p)(seq(0, 20, by = 0.25), pars_true, deriv = FALSE))
data <- subset(data, name == "V")
data$sigma <- 1
data$value <- data$value + rnorm(nrow(data), sd = data$sigma)
data <- as.datalist(data)


## Fit ---------------------------------------------------------------------
obj <- normL2(data, x*p)

starts <- msParframe(pars_true, n = 20, seed = 2, sd = 0.5)
fits <- mstrust(obj, starts, fits = 20, cores = 1, resultPath = .outdir)
fits <- as.parframe(fits)

plotValues(fits)
bestfit <- as.parvec(fits)
rbind(true = pars_true, fit = bestfit)

plotCombined((x*p)(times, bestfit), data, name == "V")


## Profiles ----------------------------------------------------------------
profiles <- profile(obj, bestfit, names(bestfit))
plotProfile(profiles)
