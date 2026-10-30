# Aligning qPCR experiments and Western blots ----
#
# Replicates of relative measurements, each experiment on its own scale, aligned as in
# blotIt (Kemmer, Bang, Rosenblatt, Timmer, Kaschek 2022, PLOS ONE 17: e0264295,
# doi:10.1371/journal.pone.0264295): a free expression for the measured value in symbols
# whose predictors are formulas over the columns of the data. First an mRNA by qPCR with
# readings below the detection limit, aligned up to a common factor. Then a phosphoprotein
# on Western blots with a background per gel, random between gels, together with targeted
# mass spectrometry in nM, whose units the aligned values take. The blotIt of Kemmer et
# al. has fixed effects only; the background that varies between gels as a random effect
# goes beyond it.

library(dMod2)
library(ggplot2)
options(dMod2.outdir = file.path(tempdir(), "alignment"))
dir.create(getOption("dMod2.outdir"), showWarnings = FALSE)

# qPCR ----
#
# The mRNA of a target gene after stimulation at seven time points, measured in three
# experiments in technical duplicates. Each experiment has its own cDNA and plate and a
# subset of the time points; no time point is in every experiment. A well reports the
# cycle Ct at which its fluorescence crosses a threshold, a level set above the background
# while the signal still grows exponentially. The amount doubles with every cycle, so
# 36 - Ct is the log2 concentration in arbitrary units, shifted by an offset of the
# experiment. Untreated cells barely express the target: a well that has not crossed the
# threshold after 40 cycles reports no Ct, and its reading lies below 36 - 40 = -4.

## --- Data ---
#
# The column `lloq` holds the limit of every reading; a reading `NA` is censored there.
set.seed(5)
times <- c(0, 0.5, 1, 2, 4, 8, 24)
log2fold <- c(0, 0.3, 1.6, 3.6, 5.2, 5.6, 4.4)
runs <- list(E1 = times[1:6], E2 = times[c(1, 3:5, 7)], E3 = times[c(2, 6, 7)])
qpcr <- do.call(rbind, lapply(names(runs), function(e)
  data.frame(experiment = e, time = rep(runs[[e]], each = 2))))
offset <- setNames(rnorm(length(runs)), names(runs))
ct <- 39.3 - log2fold[match(qpcr$time, times)] - offset[qpcr$experiment] + rnorm(nrow(qpcr), 0, 0.8)
qpcr$value <- ifelse(ct < 40, 36 - ct, NA)
qpcr$lloq <- -4
head(qpcr)
sum(is.na(qpcr$value))

## --- Model ---
#
# The reading is `yi + sj`: `yi` the log2 concentration per time point, `sj` the offset of
# the experiment. The error is absolute.
m <- alignmentModel(qpcr, model = "yi + sj", errmodel = "sigma",
                    biological = yi ~ time,
                    scaling    = sj ~ experiment,
                    error      = sigma ~ 1,
                    trafo = "sigma <- exp(sigma)")
m

## --- Fit ---
#
# The mean of the values is held at 0, so they are log2 concentrations relative to their
# mean. Their differences to the untreated value are the log2 fold changes.
fit <- align(m, fits = 1)
fit
a <- fit$aligned
a$log2fc <- a$estimate - a$estimate[a$time == 0]
cbind(a[c("time", "estimate", "se", "log2fc")], simulated = log2fold[match(a$time, times)])

# The readings of every experiment before and after the alignment. A well without Ct has
# no reading; it enters the fit as censored and is not drawn.
s <- fit$scaled[!fit$scaled$bloq, ]
both <- rbind(data.frame(s[c("experiment", "time")], stage = "before", value = s$value),
              data.frame(s[c("experiment", "time")], stage = "after", value = s$scaled))
both$stage <- factor(both$stage, c("before", "after"))
#+ fig.height = 3.5
ggplot(both, aes(factor(time), value, colour = experiment)) +
  geom_point() + facet_wrap(~ stage) +
  labs(x = "time [h]", y = "log2(conc.) [a.u.]", colour = NULL) +
  scale_color_dMod() + theme_dMod()

#+ fig.height = 3.5
plotAligned(fit) + labs(x = "time [h]", y = "log2(conc.) [a.u.]")
#+ fig.height = 3
plotResiduals(fit)

## --- The aligned values as data for a model fit ---
as.datalist(fit)

# Western blot and targeted mass spectrometry ----
#
# A phosphoprotein after stimulation with a ligand, quantified on ten gels and by targeted
# mass spectrometry. Five gels hold a short time course up to 4 h, five a long one up to
# 48 h; the two courses share only the untreated lane and 4 h. Every gel holds four of the
# time points of its course and is imaged on a scale of its own, the long gels far lower
# on the whole. A gel adds its own background to the bands, the signal of the membrane
# around them, and its image scales both. The background shows on gels with weak lanes,
# the untreated one and the late ones, and hides under the bands of the others. With a
# background of its own on every gel, the untreated level and the typical background are
# told apart only by how the background varies between gels: a normal random effect per
# gel. The mass spectrometry measures the concentration in nM at four time points of both
# courses in triplicate. Its readings are absolute: marked in the column `absolute`, they
# are predicted by the expression `absolute`, by default the first biological symbol, and
# depend on no scaling parameter. Their gel is `NA`. The block of experiments they belong
# to is not normalized but comes out in nM.

## --- Data ---
#
# One row per reading, log2 of the band in arbitrary units or of the concentration in nM.
# The gels of the short course are S1 to S5, those of the long one L1 to L5.
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

## --- Model ---
#
# A band reads `log2(2^yi + 2^(b + bj)) + sj`: the log2 concentration `yi` in nM per time
# point plus the background, shifted by `sj` per gel, all on the log2 scale. The log2
# background is the fixed `b`, common to all gels, plus `bj`, a normal random effect around
# 0 per gel. A mass spectrometry reading is `yi` itself. The error is absolute, one per
# method. The data say little about the background: its spread between gels is small, and
# the bands hide most of it. A prior keeps the standard deviation of `bj` near 0.3, its log
# `omega_bj_bj` normal around log(0.3) with standard deviation 0.3, and a weak one puts `b`
# around 1 nM, normal around 0 with standard deviation 2.
m <- alignmentModel(wbms, model = "log2(2^yi + 2^(b + bj)) + sj", errmodel = "sigma",
                    biological = yi ~ time,
                    scaling    = list(sj ~ gel, b ~ 1, bj ~ 0 + (1 | gel)),
                    error      = sigma ~ method,
                    trafo = "sigma <- exp(sigma)",
                    prior = constraintL2(mu = c(omega_bj_bj = log(0.3), b = 0),
                                         sigma = c(omega_bj_bj = 0.3, b = 2)))
m

## --- Fit ---
#
# The values come out as log2 concentrations in nM, also at the time points the mass
# spectrometry did not measure:
fit <- align(m, fits = 1)
fit
a <- fit$aligned[order(fit$aligned$time), ]
cbind(a[c("time", "estimate", "se")], nM = round(2^a$estimate, 1), simulated = conc[match(a$time, times)])

# The background in nM and the standard deviations of the two methods on the log2 scale:
p <- fit$parameters
round(2^p[["b"]], 2)
c(MS = exp(p[["sigma"]]), WB = exp(p[["sigma"]] + p[["sigma_methodWB"]]))

# The readings of every gel and of the mass spectrometry before and after the alignment,
# one colour per gel. Before it, the bands are in arbitrary units on the scales of their
# gels and the mass spectrometry in nM; after it, all are in nM. The absolute readings stay
# as they are:
s <- fit$scaled
s$gel <- factor(ifelse(s$absolute, "MS", s$gel), c(gels, "MS"))
panel <- c("before: WB [a.u.]", "before: MS [nM]", "after [nM]")
both <- rbind(data.frame(s[c("gel", "method", "time")], reading = s$value,
                         panel = ifelse(s$absolute, panel[2], panel[1])),
              data.frame(s[c("gel", "method", "time")], reading = s$scaled, panel = panel[3]))
both$panel <- factor(both$panel, panel)
#+ fig.width = 11, fig.height = 4
ggplot(both[is.finite(both$reading), ], aes(factor(time), reading, colour = gel, shape = method)) +
  geom_point() + facet_wrap(~ panel, scales = "free_y") +
  guides(x = guide_axis(angle = 45)) +
  labs(x = "time [min]", y = "log2(conc.)", colour = NULL, shape = NULL) +
  scale_color_dMod(palette = "dMod") + theme_dMod()

#+ fig.height = 3.5
plotAligned(fit) + scale_x_sqrt(breaks = times) + guides(x = guide_axis(angle = 45)) +
  labs(x = "time [min]", y = "log2(conc.) [nM]")
#+ fig.height = 3
plotResiduals(fit)

## --- Random background ---
#
# The standard deviation of the log2 background between gels, and per gel the mode of its
# random effect with its conditional standard error. The standard deviation lies near the
# mean of its prior: its value comes mostly from the prior, not from the data.
fit$omega
fit$random

# The shrinkage, 1 - SD(modes) / omega, says how much the data know about the random
# effect of a single gel. Where the bands hide the background, the mode of a gel stays
# near 0, the population value. The aligned values, `b` and omega integrate over the
# random effects and do not depend on it. But the mode of a gel is no measurement of its
# background, and plots of the modes look too narrow.
fit$fit$shrinkage

# Draws from the conditional distributions of the random effects against the estimated
# normal distribution, which they follow whatever the shrinkage. At a high shrinkage the
# draws mostly repeat that normal distribution, so the plot could hardly reveal a
# background distributed otherwise.
#+ fig.height = 4.5
plotRandomEffects(fit$fit)

# The readings against the individual prediction (IPRED), with the random effects at
# their modes, and the population prediction (PRED), with every gel at the population
# background:
#+ fig.height = 3.5
plot(fit$fit) + theme(legend.position = "right") + labs(colour = NULL)
