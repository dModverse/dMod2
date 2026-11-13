# Bile acid transport

Shows: `distributedComputing()`, `localNode()`, `profile()`, `datapointL2()`, `confint()`.

The model, the multistart and the first profiles are on the
[Get started](https://dmodverse.github.io/dMod2/articles/dMod2.html) page. This project runs the
multistart as tasks through `distributedComputing()`, removes the structural non-identifiability of
the profiles by a steady-state start and a fixed scale, and ends with a prediction band from
prediction profiles.

## Files

| file | content |
|---|---|
| `analysis.R` | the whole analysis, from the data to the prediction band; the parts below are taken from it |

## Run

With this folder as working directory: `source("analysis.R")`. The code below is shown without
being run; the fits, profiles and figures come from running the script.



## Multistart as tasks


``` r
# Multistart ----
# distributedComputing() saves the objects named in `input`, copies them with the
# code and the shared objects of the model to a job folder in the working
# directory and runs the code there in the background, as one task on the local
# machine. In the task, `cores` holds the cores of `on`, and mstrust() runs ten
# fits per core. It needs a Unix-like system; on Windows, call mstrust() itself
# with `cores = .cores`. ?distributedComputing shows how to send the same code to
# other machines with sshNodes() or to a SLURM cluster with slurmCluster().
job <- distributedComputing(
  mstrust(obj, pouter, sd = 4, fits = 10 * cores, cores = cores),
  on = localNode(.cores), jobname = "ba_fits", input = c("obj", "pouter")
)
job

# wait() reports the tasks done, running or failed until none runs. get() joins
# the parlists of the tasks; purge() deletes the job folders.
job$wait(interval = 5)
fits <- as.parframe(job$get())
job$purge(local = TRUE)
plotValues(fits)
```

## Steady-state start and fixed scale

The profiles of `S`, `TCA_CELL` and `TCA_CANA` are flat in the data. A steady-state start of the
canaliculi and a fixed scale remove the direction:


``` r
# Steady state ----
# The canaliculi start in balance with the cells, at the reflux of closed
# junctions, since they open only at time 0. The scale is fixed to 1, which
# removes the direction.
trafo <- eqnvec() |>
  define("x <- x", x = innerpars) |>
  define("TCA_buffer <- 0") |>
  define("TCA_cana <- TCA_cell*k_export_cana/k_reflux") |>
  insert("x <- 10^y", x = .currentSymbols, y = toupper(.currentSymbols)) |>
  branch(conditions = c("closed", "open")) |>
  define("k_reflux <- 10^K_REFLUX_OPEN", conditionMatch = "open") |>
  insert("S <- 0")
p <- P(trafo, modelname = "ba_trafo_ss", compile = TRUE)

outerpars <- getParameters(p)
pouter    <- structure(rep(-1, length(outerpars)), names = outerpars)
obj  <- normL2(data, g*x*p, attrName = "data") +
  constraintL2(pouter, sigma = 20, attrName = "prior")

# `wait = TRUE` waits for the tasks, returns their results and deletes the job.
fits <- distributedComputing(
  mstrust(obj, pouter, sd = 4, fits = 10 * cores, cores = cores),
  on = localNode(.cores), jobname = "ba_fits_ss", input = c("obj", "pouter"),
  wait = TRUE
)
plotValues(fits)
bestfit <- as.parvec(fits)

profiles <- profile(obj, bestfit, whichPar = names(bestfit), method = "integrate",
                    cores = .cores, limits = c(lower = -5, upper = 5),
                    stepControl = list(stop = "data"),
                    algoControl = list(reoptimize = TRUE))
plotProfile(profiles, mode %in% c("data", "prior"))
```


``` r
# Fixed reflux ----
# K_REFLUX_OPEN stays practically non-identifiable: the profile is open towards
# instantaneous reflux, so the reflux of open junctions is fixed to a fast rate.
# K_EXPORT_SINUS stays open towards slow export below the 68% level.
trafo <- eqnvec() |>
  define("x <- x", x = innerpars) |>
  define("TCA_buffer <- 0") |>
  define("TCA_cana <- TCA_cell*k_export_cana/k_reflux") |>
  insert("x <- 10^y", x = .currentSymbols, y = toupper(.currentSymbols)) |>
  branch(conditions = c("closed", "open")) |>
  define("k_reflux <- 10^6", conditionMatch = "open") |>
  insert("S <- 0")
p <- P(trafo, modelname = "ba_trafo_final", compile = TRUE)

obj     <- normL2(data, g*x*p, attrName = "data")
bestfit <- trust(obj, bestfit[getParameters(p)], rinit = 1)$argument
#+ fig.height = 3.5
plot((g*x*p)(times, bestfit, deriv = FALSE), data) + labs(color = "junctions")

# the parameters with their 95% confidence intervals, on log10
profiles <- profile(obj, bestfit, whichPar = names(bestfit), method = "integrate",
                    cores = .cores, limits = c(lower = -5, upper = 5),
                    algoControl = list(reoptimize = TRUE))
confint(profiles, level = 0.95)
```

## Prediction band


``` r
# Prediction profile ----
# The cellular bile acid at time t becomes a parameter v, tied to the model by a
# data point with a small sigma. Along the profile of v the data term is the
# prediction profile, and the profile interval the confidence interval of the prediction.
.predict <- function(t) {
  objv <- normL2(data, g*x*p, times = t, attrName = "data") +
    datapointL2(name = "cellular", time = t, parameter = "v", sigma = 0.1,
                attrName = "validation", condition = "closed")
  pred <- (g*x*p)(t, bestfit, deriv = FALSE)[["closed"]]
  v0   <- unname(pred[pred[, "time"] == t, "cellular"])
  fit  <- trust(objv, c(v = v0, bestfit), rinit = 1)
  profile(objv, fit$argument, "v", method = "integrate",
          stepControl = list(stop = "data"),
          algoControl = list(reoptimize = TRUE))
}

prof20 <- .predict(20)
plotProfile(prof20, mode %in% c("validation", "data")) + labs(x = "cellular bile acid at t = 20")
confint(prof20, val.column = "data")

# the interval at several times, drawn as a band around the prediction
.band <- do.call(rbind, lapply(c(0:10, seq(12.5, 50, by = 2.5)), function(t)
  data.frame(time = t, confint(.predict(t), val.column = "data")[c("lower", "upper")])))

.pred <- subset(as.data.frame((g*x*p)(times, bestfit, deriv = FALSE)),
                condition == "closed" & name == "cellular")
ggplot(.pred, aes(time, value)) +
  geom_ribbon(data = .band, aes(x = time, ymin = lower, ymax = upper),
              inherit.aes = FALSE, alpha = 0.3) +
  geom_line() +
  geom_point(data = .band, aes(y = lower), shape = 4) +
  geom_point(data = .band, aes(y = upper), shape = 4) +
  theme_dMod() +
  labs(x = "time", y = "cellular bile acid", title = "closed junctions: prediction band")
```
