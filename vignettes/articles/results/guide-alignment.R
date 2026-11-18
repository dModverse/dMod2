# The data, models and fits of guide chapter 9, guide-alignment.Rmd, which shows the labeled parts
# and reads the fits from guide-alignment.rds. Run in this folder.

library(dMod2)
options(dMod2.outdir = file.path(tempdir(), "guide-alignment"))
dir.create(getOption("dMod2.outdir"), showWarnings = FALSE)

## ---- data
set.seed(1)
times <- c(0, 10, 30, 60, 120)
truth <- c(0.3, 1, 4, 2, 0.5)
gels <- paste0("G", 1:6)
samples <- expand.grid(time = times, gel = gels, stringsAsFactors = FALSE)
samples$level <- truth[match(samples$time, times)] * rlnorm(nrow(samples), 0, 0.2)
scale <- setNames(rlnorm(6, 7, 0.5), gels)
background <- setNames(rlnorm(6, log(100), 0.3), gels)
lane <- rep(seq_len(nrow(samples)), each = 3)
d <- data.frame(samples[lane, c("time", "gel")], lane = 1:3, row.names = NULL)
d$value <- (scale[d$gel] * samples$level[lane] + background[d$gel]) * rlnorm(nrow(d), 0, 0.15)
head(d)

## ---- model
logs <- c("yi <- exp(yi)", "sj <- exp(sj)", "bj <- exp(bj)", "sigmaR <- exp(sigmaR)")
m <- alignmentModel(d, model = "yi * sj + bj", errmodel = "sigmaR * prediction",
                    biological = yi ~ time, scaling = sj ~ gel, error = bj + sigmaR ~ 1, trafo = logs,
                    technical = "lane")
m

## ---- fit
fit <- align(m, fits = 4)
fit

## ---- lanes
mLanes <- alignmentModel(d, model = "yi * sj + bj", errmodel = "sigmaR * prediction",
                         biological = yi ~ time, scaling = sj ~ gel, error = bj + sigmaR ~ 1,
                         trafo = logs)
fitLanes <- align(mLanes, fits = 4)
fitLanes

## ---- background
m2 <- alignmentModel(d, model = "yi * sj + bj", errmodel = "sigmaR * prediction",
                     biological = yi ~ time, scaling = list(sj ~ gel, bj ~ 1 + (1 | gel)),
                     error = sigmaR ~ 1, trafo = logs, technical = "lane",
                     prior = constraintL2(mu = c(omega_bj_bj = log(0.3)), sigma = 0.3))
fit2 <- align(m2, fits = 4)
fit2

## ---- store-replicates
result <- list(fit = fit, lanes = fitLanes, background = fit2)

## ---- qpcrdata
set.seed(5)
times <- c(0, 0.5, 1, 2, 4, 8, 24)
log2fold <- c(0, 0.3, 1.6, 3.6, 5.2, 5.6, 4.4)
runs <- list(E1 = times[1:6], E2 = times[c(1, 3:5, 7)], E3 = times[c(2, 3, 6, 7)], E4 = times[c(1, 2, 4:7)])
samples <- do.call(rbind, lapply(names(runs), function(e) data.frame(experiment = e, time = runs[[e]])))
offset <- setNames(rnorm(length(runs)), names(runs))
samples$ct <- 39.3 - log2fold[match(samples$time, times)] - offset[samples$experiment] +
  rnorm(nrow(samples), 0, 0.5)
well <- rep(seq_len(nrow(samples)), each = 2)
qpcr <- data.frame(samples[well, c("experiment", "time")], well = 1:2, row.names = NULL)
ct <- samples$ct[well] + rnorm(nrow(qpcr), 0, 0.25)
qpcr$value <- ifelse(ct < 40, 36 - ct, NA)
qpcr$lloq <- -4
head(qpcr)
sum(is.na(qpcr$value))

## ---- qpcrmodel
m <- alignmentModel(qpcr, model = "yi + sj", errmodel = "sigma",
                    biological = yi ~ time,
                    scaling    = sj ~ experiment,
                    error      = sigma ~ 1,
                    trafo = "sigma <- exp(sigma)", technical = "well")
m

## ---- qpcrfit
fit <- align(m, fits = 1)
fit

## ---- store-qpcr
result$qpcr <- fit

## ---- wbdata
set.seed(6)
short <- c(0, 10, 30, 60, 120, 240)
long <- c(0, 240, 480, 1440, 2880)
times <- sort(union(short, long))
conc <- c(1, 40, 85, 60, 25, 10, 6, 3, 2)
gels <- c(sprintf("S%d", 1:5), sprintf("L%d", 1:5))
course <- setNames(rep(c("short", "long"), each = 5), gels)
shift <- setNames(ifelse(course == "short", 6, 1) + rnorm(10, 0, 1.5), gels)
background <- setNames(1.5 + rnorm(10, 0, 0.3), gels)
wb <- do.call(rbind, lapply(gels, function(g)
  data.frame(method = "WB", gel = g, time = sort(sample(if (course[[g]] == "short") short else long, 4)))))
wb$value <- log2(conc[match(wb$time, times)] + 2^background[wb$gel]) + shift[wb$gel] +
  rnorm(nrow(wb), 0, 0.7)
ms <- data.frame(method = "MS", gel = NA, time = rep(c(0, 30, 240, 1440), each = 3))
ms$value <- log2(conc[match(ms$time, times)]) + rnorm(nrow(ms), 0, 0.5)
wbms <- rbind(cbind(wb, absolute = FALSE), cbind(ms, absolute = TRUE))
head(wbms)

## ---- wbmodel
m <- alignmentModel(wbms, model = "log2(2^yi + 2^(b + bj)) + sj", errmodel = "sigma",
                    biological = yi ~ time,
                    scaling    = list(sj ~ gel, b ~ 1, bj ~ 0 + (1 | gel)),
                    error      = sigma ~ method,
                    trafo = "sigma <- exp(sigma)",
                    prior = constraintL2(mu = c(omega_bj_bj = log(0.3), b = 0),
                                         sigma = c(omega_bj_bj = 0.3, b = 2)))
m

## ---- wbfit
fit <- align(m, fits = 1)
fit

## ---- save
result$wb <- fit
attr(result, "stamp") <- list(md5 = tools::md5sum("guide-alignment.R"),
                              dMod2 = as.character(packageVersion("dMod2")))
saveRDS(result, "guide-alignment.rds")
