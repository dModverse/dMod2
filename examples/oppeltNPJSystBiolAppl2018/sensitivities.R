# Forward and reverse derivatives of the full model at its best fit: the gradient by both
# sweeps, their cost, and the exact Hessian by forward over reverse against Gauss-Newton. Saves
# results/sensitivities.rds, which the guide reads.

# model.R, unless it ran in this session
if (!exists("obj")) source("model.R")

## ---- reversechain
# the chain of model.R once more, every function object with the reverse sweep and second
# derivatives; Pimpl() passes the reverse sweep through the Jacobian of its root
dm <- c("forward", "reverse", "forward-reverse")
modelR <- odemodel(reactions, events = events, modelname = "oppelt_rev_ode", deriv2 = TRUE,
                   derivMode = dm, compile = FALSE)
xR <- Xs(modelR, options = list(atol = 1e-10, rtol = 1e-8))
gR <- Y(observables, xR, modelname = "oppelt_rev_obs", compile = FALSE, deriv2 = TRUE, derivMode = dm)
eR <- Y(errors, gR, modelname = "oppelt_rev_err", compile = FALSE, deriv2 = TRUE, derivMode = dm)
pSSR <- Pimpl(reactions, modelname = "oppelt_rev_ss", compile = FALSE, deriv2 = TRUE)
pR <- P(trafo, modelname = "oppelt_rev_par", compile = FALSE, deriv2 = TRUE, derivMode = dm)
compile(gR, xR, pSSR, pR, eR, output = "oppelt_rev", cores = getOption("dMod2.cores"))
objR <- normL2(data, gR * xR * pSSR * pR, eR) + constraintL2(mu, sigma = sigma)

## ---- gradients
# the data term, whose gradient the prior balances at the optimum, so that it does not vanish
best <- as.parvec(readRDS("results/fits.rds"))
dataTerm <- normL2(data, gR * xR * pSSR * pR, eR)
forward <- dataTerm(best, fixed = gauge, sweep = "forward")
reverse <- dataTerm(best, fixed = gauge, sweep = "reverse")
gradients <- data.frame(forward = forward$gradient,
                        reverse = reverse$gradient[names(forward$gradient)])
# the largest deviation relative to the largest entry of the gradient
max(abs(gradients$forward - gradients$reverse)) / max(abs(gradients$forward))

## ---- cost
# the shortest of ten evaluations, in seconds
cost <- function(...) {
  args <- c(list(best, fixed = gauge), list(...))
  min(replicate(10, system.time(do.call(dataTerm, args))[["elapsed"]]))
}
timing <- c(value = cost(deriv = FALSE), forward = cost(sweep = "forward"),
            reverse = cost(sweep = "reverse"),
            `forward-reverse` = cost(sweep = "reverse", deriv2 = TRUE))
timing

## ---- forwardforward
# one evaluation of the exact Hessian forward over forward, with second-order sensitivities in all
# pairs of parameters: its time and the peak memory of the process in GB
peak <- function() as.numeric(gsub("\\D", "", grep("^VmHWM", readLines("/proc/self/status"), value = TRUE))) / 2^20
before <- peak()
forwardForward <- system.time(ff <- dataTerm(best, fixed = gauge, sweep = "forward", deriv2 = TRUE))[["elapsed"]]
memory <- c(before = before, forwardForward = peak())
c(seconds = forwardForward, memory)

## ---- hessians
exact <- dataTerm(best, fixed = gauge, sweep = "reverse", deriv2 = TRUE)$hessian
# the two exact Hessians agree up to the error of the solver
max(abs(ff$hessian[rownames(exact), colnames(exact)] - exact)) / max(abs(exact))
gaussNewton <- forward$hessian[rownames(exact), colnames(exact)]
hessians <- data.frame(exact = sort(eigen(exact, symmetric = TRUE, only.values = TRUE)$values),
                       gaussNewton = sort(eigen(gaussNewton, symmetric = TRUE, only.values = TRUE)$values))

dir.create("results", showWarnings = FALSE)
saveRDS(list(gradients = gradients, timing = timing, forwardForward = forwardForward, memory = memory,
             hessians = hessians, exact = exact, forwardForwardHessian = ff$hessian,
             gaussNewton = gaussNewton), "results/sensitivities.rds")
