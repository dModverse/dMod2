# Figures from the saved fits: the multistart, the parameters of the best fits and the best fit
# against the data.

# model.R, unless it ran in this session
if (!exists("obj")) source("model.R")

## ---- load
library(ggplot2)
fits <- readRDS("results/fits.rds")
best <- as.parvec(fits)

## ---- waterfall
# every fit on a pseudo-log scale; each plateau is a local optimum
plotValues(fits, pseudoLog = TRUE)

## ---- parameters
# the parameters of the fits on the best plateaus, one color per plateau
plotPars(fits, value < min(value) + 5) +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  labs(x = NULL, y = "log value", color = "plateau")

## ---- prediction
# the best fit with one standard deviation of the error model as a band, each observable up to
# the last time point of the data of the observable, at least 30 minutes
times <- seq(0, 360, length.out = 361)
last <- pmax(with(lbind(data), tapply(time, name, max)), 30)
plotCombined(prd(times, c(best, gauge), deriv = FALSE), data, time <= last[as.character(name)],
             errfn = e) +
  facet_wrap(~name, scales = "free", ncol = 4) +
  labs(x = "time [min]", y = "log signal [a.u.]", color = "Treatment", fill = "Treatment")

## ---- band
# the prediction band of nuclear IκBα in the stimulated cells, from prediction.R
if (!exists("prdReduced")) source("reduced.R")
controls(gReduced, name = "attachInput") <- TRUE
band <- readRDS("results/prediction.rds")$band
bestReduced <- as.parvec(readRDS("results/reduced.rds")$fits)
course <- data.frame(time = times,
                     value = prdReduced(times, c(bestReduced, gaugeReduced), deriv = FALSE)$TNF[, "nIkb"])
ggplot(course, aes(time, value)) +
  geom_ribbon(data = band, aes(x = time, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.3) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 360, 60)) +
  scale_y_log10() +
  labs(x = "time [min]", y = "nuclear IκBα [total NFκB]") +
  theme_dMod2()

## ---- dcfpaths
# the differences by DCF over the strength of the L1 penalty, from dcf.R, in color those kept at
# the selected strength
dcfL1 <- readRDS("results/dcfL1.rds")
paths <- data.frame(lambda = rep(as.numeric(rownames(dcfL1$path)), ncol(dcfL1$path)),
                    difference = rep(sub("^eta_(.*)_DCF$", "\\1", colnames(dcfL1$path)),
                                     each = nrow(dcfL1$path)),
                    value = as.vector(dcfL1$path))
paths <- paths[paths$lambda > 0, ]
kept <- sub("^eta_(.*)_DCF$", "\\1", dcfL1$selected)
ggplot(paths, aes(lambda, value, group = difference)) +
  geom_vline(xintercept = dcfL1$lambda, linetype = 2, color = "gray50") +
  geom_line(data = paths[!paths$difference %in% kept, ], color = "gray75") +
  geom_line(data = paths[paths$difference %in% kept, ], aes(color = difference)) +
  scale_x_log10(breaks = 10^(-2:4), labels = trimws(format(10^(-2:4), scientific = FALSE, drop0trailing = TRUE))) +
  labs(x = expression(lambda), y = expression(Delta ~ "[log]"), color = NULL) +
  theme_dMod2()
