# Figures from the saved fits: the multistart, the parameters of the best fits and the best fit
# against the data.

# model.R, unless it ran in this session
if (!exists("obj")) source("model.R")

## ---- load
library(ggplot2)
fits <- readRDS("results/fits.rds")
best <- as.parvec(fits)

## ---- waterfall
plotValues(fits)

## ---- parameters
# the dynamic parameters of the plateaus within 4 of the best value, one color per plateau
dynamic <- grep("^(scale|offset|sd)[0-9]*_", names(best), invert = TRUE, value = TRUE)
keep <- c(setdiff(names(fits), attr(fits, "parameters")), dynamic)
plotPars(fits[, keep], value < min(value) + 4) +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(breaks = seq(-4, 12, 2)) +
  guides(x = guide_axis()) +
  labs(x = NULL, y = "log10 value", color = "plateau")

## ---- helpers
# what each condition recorded; the error model `e` draws one standard deviation as a band
recorded <- unique(paste(bachmann$condition, bachmann$name))
times <- seq(0, 240, length.out = 500)
minutes <- scale_x_continuous(breaks = seq(0, 240, 60))
# panel titles without the unit, which the axis names
observable <- as_labeller(function(name) sub("_au", "", name))

## ---- reference
long <- prd(times, best, conditions = "long", deriv = FALSE)
plotCombined(long, data["long"], paste(condition, name) %in% recorded, errfn = e) +
  facet_wrap(~name, scales = "free", labeller = observable) +
  minutes +
  labs(x = "time [min]", y = "log10 signal [a.u.]") +
  theme(legend.position = "none")

## ---- perturbations
perturbed <- rownames(grid)[grid$experiment %in% c("actd", "cisoe", "socs3oe", "shp1oe")]
shown <- c("pJAK2_au", "pEpoR_au", "pSTAT5_au", "CIS_au")
pred <- prd(times, best, conditions = perturbed, deriv = FALSE)
perturbations <- c(actd = "actinomycin D", cisoe = "CIS overexpressed",
                   socs3oe = "SOCS3 overexpressed", shp1oe = "SHP1 overexpressed")
plotCombined(pred, data[perturbed], paste(condition, name) %in% recorded, name %in% shown,
             aesthetics = list(color = "perturbation", group = "condition"), errfn = e) +
  facet_grid(rows = vars(name), cols = vars(factor(experiment, names(perturbations), perturbations)),
             scales = "free_y", labeller = labeller(name = observable)) +
  scale_color_dMod2(limits = c("control", unname(perturbations))) +
  scale_fill_dMod2(limits = c("control", unname(perturbations))) +
  scale_x_continuous(breaks = c(0, 120, 240)) +
  labs(x = "time [min]", y = "log10 signal [a.u.]", color = NULL) +
  theme(legend.position = "bottom", panel.spacing.x = unit(1, "lines"))

## ---- dose
# the dose swept through one condition per experiment, with the scales and measurement time of the experiment
doseData <- subset(bachmann, grepl("^dr", experiment))
doses <- 10^seq(log10(min(doseData$epo_level)), log10(max(doseData$epo_level)), length.out = 300)
inner <- p(best, deriv = FALSE)
sweep <- function(ex) {
  cn   <- rownames(grid)[grid$experiment == ex][1]
  tt   <- unique(doseData$time[doseData$experiment == ex])
  seen <- unique(doseData$name[doseData$experiment == ex])
  pin  <- setNames(as.numeric(inner[[cn]]), names(inner[[cn]]))
  do.call(rbind, lapply(doses, function(d) {
    pin["epo_level"] <- d
    out <- as.data.frame((g * x)(sort(unique(c(0, tt))), pin, deriv = FALSE), errfn = e)
    cbind(out[out$time == tt & out$name %in% seen, c("name", "value", "sigma")], epo_level = d, experiment = ex)
  }))
}
dosePred <- do.call(rbind, lapply(unique(doseData$experiment), sweep))
# each dose response is measured at one time
at <- tapply(doseData$time, doseData$experiment, `[`, 1)
measuredAt <- function(ex) factor(paste(at[ex], "min"), paste(sort(at), "min"))
dosePred$measured <- measuredAt(dosePred$experiment)
doseData$measured <- measuredAt(doseData$experiment)
ggplot(dosePred, aes(epo_level, value, color = measured, fill = measured)) +
  geom_ribbon(aes(ymin = value - sigma, ymax = value + sigma), alpha = 0.25, linetype = 0,
              show.legend = FALSE) +
  geom_line() +
  geom_point(data = doseData) +
  facet_wrap(~name, scales = "free_y", labeller = observable) +
  scale_x_log10(breaks = 10^(-9:-5), labels = scales::label_log()) +
  scale_color_dMod2() +
  scale_fill_dMod2() +
  theme_dMod2() +
  labs(x = "Epo [units/cell]", y = "log10 signal [a.u.]", color = "measured at", fill = "measured at")
