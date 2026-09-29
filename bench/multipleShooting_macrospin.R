# -------------------------------------------------------------------------#
# Macrospin ringdown: single against multiple shooting from the same starts
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# The LLG macrospin of bench/symmetry_macrospin.R in both fields, fitted from
# 100 starts drawn uniformly from the box of prior knowledge, by single and by
# multiple shooting. In the tilted field the Kerr signal rings about eleven
# times in 1.5 ns; from a start with the wrong precession frequency single
# shooting ends in minima that fit the easy axis and average the ringing out.
# Multiple shooting lays one segment per half oscillation of the data, so no
# segment holds enough of the ringing for its phase to slip by a period.
# A start is a hit when the single-shooting objective at its end lies within 1
# of the best.
#
# On 29 Sep 2026, one core per fit, the least of three runs: single shooting
# hit from 52 of 100 starts at a median of 0.23 s, multiple shooting from 63 at
# 2.3 s; in Cartesian coordinates 52 and 83. The misses of single shooting fit
# the easy axis and average the ringing out.
#
# The spherical chart is singular at th = 0, which the easy-axis field relaxes
# to, and a segment that passes near the pole winds ph one turn more or less
# for an arbitrarily small change of its start. ph therefore lives on an angle
# chart, whose gaps count modulo 2 pi. The same physics in Cartesian
# coordinates, the second block, has no pole and needs no chart.
#
# [AUTHOR]
# Simon Beyer
#
# [Date]
# Tue 29 Sep 2026
#
# [Info]
# Units as in bench/symmetry_macrospin.R: time in ns, fields in mT, gamma in
# rad/(ns mT). The ODE models are built with includeTimeZero = FALSE, which
# multiple shooting needs.
# -------------------------------------------------------------------------#

library(dMod2)
library(ggplot2)

.modelname <- "msMacrospin"
.outdir    <- file.path(tempdir(), .modelname)
.cores     <- detectFreeCores()
.fits      <- as.integer(Sys.getenv("MS_FITS", "100"))

if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Model, spherical, as in bench/symmetry_macrospin.R
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
dEdth <- "(Hk*sin(th)*cos(th) + Hz*sin(th) - Hx*cos(th)*cos(ph))"
dEdph <- "(Hx*sin(th)*sin(ph))"

f <- eqnvec(
  th = sprintf("-gam/(1 + al^2)*(al*%s + %s/sin(th))", dEdth, dEdph),
  ph = sprintf("gam/(1 + al^2)*(%s/sin(th) - al*%s/sin(th)^2)", dEdth, dEdph))

opts <- list(atol = 1e-8, rtol = 1e-8, maxsteps = 1e5)
model <- odemodel(f, modelname = "msMacro_ode", compile = FALSE, outdir = .outdir,
                  method = "tsit5", includeTimeZero = FALSE)
x <- Xs(model, optionsOde = opts, optionsSens = opts)
g <- Y(eqnvec(y = "s*cos(th)"), f = x, modelname = "msMacro_obs", compile = FALSE,
       outdir = .outdir)

innerpars <- getParameters(g, x)
fields    <- data.frame(Hz = c(100, 30), Hx = c(0, 100), row.names = c("easy", "tilt"))

p <- Reduce(`+`, lapply(rownames(fields), function(cn) {
  eqnvec() |>
    define("x~x", x = innerpars) |>
    define("x~y", x = c("Hz", "Hx"), y = unlist(fields[cn, ])) |>
    define("ph~0") |>
    define(paste0("th~th0_", cn)) |>
    insert("x~10^x", x = c("gam", "al", "Hk", "s")) |>
    P(condition = cn, modelname = paste0("msMacro_p_", cn), compile = FALSE,
      outdir = .outdir)
}))

compile(x, g, p, output = .modelname, cores = .cores)
prd <- g*x*p


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Data, as in bench/symmetry_macrospin.R
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
truth <- c(gam = log10(0.176), al = log10(0.05), Hk = log10(300), s = 0,
           th0_easy = 0.6, th0_tilt = 0.6)
times <- seq(0, 1.5, by = 0.005)
.sigma <- 0.005
set.seed(1111)
sim <- as.data.frame(prd(seq(0, 1.5, by = 0.01), truth, deriv = FALSE))
sim$value <- sim$value + rnorm(nrow(sim), 0, .sigma)
sim$sigma <- .sigma
mydata <- as.datalist(sim[, c("name", "time", "value", "sigma", "condition")],
                      split.by = "condition")
plot(prd(times, truth, deriv = FALSE), mydata)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Single and multiple shooting from the same starts
#
# The layout multiple shooting starts from: ten segments on the easy axis,
# one per half oscillation in the tilted field.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
obj.ss <- normL2(mydata, prd, times = times)
obj.ms <- normL2(mydata, prd, times = times,
                 multipleShootingControl = list(charts = c(ph = "angle")))
lengths(attr(dMod2:::.shootObjOf(obj.ms), "spec")$tau)

lower <- c(gam = log10(0.088), al = -3, Hk = 1, s = -0.3, th0_easy = 0.1, th0_tilt = 0.1)
upper <- c(gam = log10(0.3), al = log10(0.3), Hk = 3, s = 0.3, th0_easy = 1.2,
           th0_tilt = 1.2)
set.seed(42)
starts <- t(replicate(.fits, runif(length(lower), lower, upper)))
colnames(starts) <- names(lower)

ss <- function(th) obj.ss(th, deriv = FALSE)$value
fitOne <- function(obj, s)
  suppressWarnings(trust(obj, s, rinit = 0.1, rmax = 10, iterlim = 1000, cores = 1L))

jobs <- expand.grid(method = c("single shooting", "multiple shooting"),
                    start = seq_len(.fits), stringsAsFactors = FALSE)
# one start per core, each fit on one
out <- do.call(rbind, parallel::mclapply(seq_len(nrow(jobs)), mc.cores = .cores,
                                         mc.preschedule = FALSE, function(j) {
  jb <- jobs[j, ]
  t0 <- Sys.time()
  r <- try(fitOne(if (jb$method == "single shooting") obj.ss else obj.ms,
                  starts[jb$start, ]), silent = TRUE)
  sec <- as.numeric(Sys.time() - t0, units = "secs")
  data.frame(jb, ss = if (inherits(r, "try-error")) NA_real_ else ss(r$argument),
             sec = sec, neval = if (inherits(r, "try-error")) NA_integer_ else r$neval)
}))
best <- min(out$ss, na.rm = TRUE)
out$hit <- !is.na(out$ss) & out$ss < best + 1
aggregate(cbind(hits = hit, sec, neval) ~ method, out,
          function(z) c(sum = sum(z), median = median(z)))

# the misses of single shooting: the easy axis fits, the ringing is averaged out
ggplot(out, aes(rank(ss, ties.method = "first"), pmin(ss - best, 1e4), colour = method)) +
  geom_point(size = 0.8) + facet_wrap(~method) +
  scale_y_continuous(trans = "log1p") +
  labs(x = "fit rank", y = "objective - best") + theme_dMod(base_size = 9)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# The same in Cartesian coordinates
#
# dm/dt = -g' m x H - g' al m x (m x H), g' = gam/(1 + al^2),
# H = (Hx, 0, Hz + Hk mz); no pole, every state on a linear chart.
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
mxH <- c(mx = "(my*(Hz + Hk*mz))", my = "(mz*Hx - mx*(Hz + Hk*mz))", mz = "(-my*Hx)")
mmxH <- c(mx = sprintf("(my*%s - mz*%s)", mxH[["mz"]], mxH[["my"]]),
          my = sprintf("(mz*%s - mx*%s)", mxH[["mx"]], mxH[["mz"]]),
          mz = sprintf("(mx*%s - my*%s)", mxH[["my"]], mxH[["mx"]]))
fc <- eqnvec(mx = sprintf("-gam/(1 + al^2)*(%s + al*%s)", mxH[["mx"]], mmxH[["mx"]]),
             my = sprintf("-gam/(1 + al^2)*(%s + al*%s)", mxH[["my"]], mmxH[["my"]]),
             mz = sprintf("-gam/(1 + al^2)*(%s + al*%s)", mxH[["mz"]], mmxH[["mz"]]))

modelc <- odemodel(fc, modelname = "msMacroC_ode", compile = FALSE, outdir = .outdir,
                   method = "tsit5", includeTimeZero = FALSE)
xc <- Xs(modelc, optionsOde = opts, optionsSens = opts)
gc <- Y(eqnvec(y = "s*mz"), f = xc, modelname = "msMacroC_obs", compile = FALSE,
        outdir = .outdir)
pc <- Reduce(`+`, lapply(rownames(fields), function(cn) {
  eqnvec() |>
    define("x~x", x = getParameters(gc, xc)) |>
    define("x~y", x = c("Hz", "Hx"), y = unlist(fields[cn, ])) |>
    define("mx~sin(th0)", th0 = paste0("th0_", cn)) |>
    define("my~0") |>
    define("mz~cos(th0)", th0 = paste0("th0_", cn)) |>
    insert("x~10^x", x = c("gam", "al", "Hk", "s")) |>
    P(condition = cn, modelname = paste0("msMacroC_p_", cn), compile = FALSE,
      outdir = .outdir)
}))
compile(xc, gc, pc, output = paste0(.modelname, "C"), cores = .cores)
prdc <- gc*xc*pc

objc.ss <- normL2(mydata, prdc, times = times)
objc.ms <- normL2(mydata, prdc, times = times, multipleShootingControl = TRUE)
c(spherical = ss(truth), cartesian = objc.ss(truth, deriv = FALSE)$value)

outc <- do.call(rbind, parallel::mclapply(seq_len(nrow(jobs)), mc.cores = .cores,
                                          mc.preschedule = FALSE, function(j) {
  jb <- jobs[j, ]
  t0 <- Sys.time()
  r <- try(fitOne(if (jb$method == "single shooting") objc.ss else objc.ms,
                  starts[jb$start, ]), silent = TRUE)
  sec <- as.numeric(Sys.time() - t0, units = "secs")
  data.frame(jb, ss = if (inherits(r, "try-error")) NA_real_
                      else objc.ss(r$argument, deriv = FALSE)$value,
             sec = sec, neval = if (inherits(r, "try-error")) NA_integer_ else r$neval)
}))
outc$hit <- !is.na(outc$ss) & outc$ss < min(outc$ss, na.rm = TRUE) + 1
aggregate(cbind(hits = hit, sec, neval) ~ method, outc,
          function(z) c(sum = sum(z), median = median(z)))
