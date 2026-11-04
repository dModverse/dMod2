# Figures from the saved fits: the multistart, the parameters of the best fits and the best fit
# against the data.

# model.R, unless it ran in this session
if (!exists("obj")) source("model.R")

## ---- load
library(ggplot2)
fits <- readRDS("results/fits.rds")
best <- as.parvec(fits)

## ---- waterfall
# every fit, on a pseudo-log scale above the best; each plateau is a local optimum
plotValues(fits, pseudoLog = TRUE)

## ---- parameters
# the parameters of the fits on the best plateaus, one colour per plateau
plotPars(fits, value < min(value) + 5) +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  labs(x = NULL, y = "log value", colour = "plateau")

## ---- prediction
# the best fit with one standard deviation of the error model as a band
times <- seq(0, 360, length.out = 361)
plotCombined(prd(times, c(best, gauge), deriv = FALSE), data, errfn = e) +
  facet_wrap(~name, scales = "free_y", ncol = 4) +
  scale_x_continuous(breaks = seq(0, 360, 120)) +
  labs(x = "time [min]", y = "log signal [a.u.]", colour = "Treatment", fill = "Treatment")

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
  theme_dMod()

## ---- changes
# the changes of the rate constants by the compounds, from compounds.R: colour the mode, opacity
# the evidence against the spike, -log10 p0 cut at 6; a circle where p0 < 0.05
changes <- readRDS("results/compounds.rds")$random
changes$eta <- factor(changes$eta, levels = unique(changes$eta))
changes$evidence <- -log10(pmax(changes$p0, 1e-6))
ggplot(changes, aes(group, eta)) +
  geom_tile(aes(fill = mode, alpha = evidence), colour = "white") +
  geom_point(data = changes[changes$p0 < 0.05, ], aes(shape = "p0 < 0.05"), size = 1.8, fill = "white") +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
  scale_alpha(range = c(0.1, 1), limits = c(0, 6), breaks = 1:3,
              labels = c("1, p0 = 0.1", "2, p0 = 0.01", "3, p0 = 0.001"),
              guide = guide_legend(override.aes = list(fill = "grey30"))) +
  scale_shape_manual(values = 21) +
  labs(x = NULL, y = NULL, fill = "change [log]", alpha = "evidence, -log10 p0", shape = NULL) +
  theme_dMod()
