## Five time points on six gels, each gel with its own scale and an additive background
set.seed(1)
d <- expand.grid(time = c("0", "10", "30", "60", "120"), gel = paste0("G", 1:6),
                 stringsAsFactors = FALSE)
truth <- c(0.02, 1, 4, 2, 0.5)[match(d$time, c("0", "10", "30", "60", "120"))]
d$value <- (rlnorm(6, 7, 0.5)[match(d$gel, unique(d$gel))] * truth + 50) * rlnorm(nrow(d), 0, 0.1)

\donttest{
## every symbol positive, on the log scale of its predictor; one value per time point,
## one scale per gel
logs <- c("yi <- exp(yi)", "sj <- exp(sj)", "bj <- exp(bj)", "sigmaR <- exp(sigmaR)")
m <- alignmentModel(d, model = "yi * sj + bj", errmodel = "sigmaR * prediction",
                    biological = yi ~ time, scaling = sj ~ gel, error = bj + sigmaR ~ 1,
                    trafo = logs)
## a background varying between gels: a list of formulas for one role
m2 <- alignmentModel(d, model = "yi * sj + bj", errmodel = "sigmaR * prediction",
                     biological = yi ~ time, scaling = list(sj ~ gel, bj ~ 1 + (1 | gel)),
                     error = sigmaR ~ 1, trafo = logs)
m
m2
fit <- align(m, fits = 4)
fit
fit$aligned
plotAligned(fit)
}
