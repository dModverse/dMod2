# Fits of the FitzHugh-Nagumo oscillator by single and multiple shooting from the same 40 starts,
# for the article on multiple shooting. Run in this directory; saves fitzhughNagumo.rds.

library(dMod2)

## ---- fhnmodel
fhnDir <- file.path(tempdir(), "fitzhughNagumo")
dir.create(fhnDir, showWarnings = FALSE)
fhn <- eqnvec(V = "c*(V - V^3/3 + R)", R = "-(V - a + b*R)/c")
fhnModel <- odemodel(fhn, modelname = "fhn_ode", includeTimeZero = FALSE, compile = FALSE,
                     outdir = fhnDir)
fx <- Xs(fhnModel, options = list(atol = 1e-10, rtol = 1e-10))
fg <- Y(c(y = "V"), f = fx, modelname = "fhn_obs", compile = FALSE, outdir = fhnDir)
fp <- P(eqnvec(V = "V", R = "R", a = "a", b = "b", c = "exp(lc)"), condition = "C",
        modelname = "fhn_p", compile = FALSE, outdir = fhnDir)
compile(fx, fg, fp, output = "fhn", cores = 6)
fhnPrd <- fg * fx * fp

fhnTruth <- c(V = -1, R = 1, a = 0.2, b = 0.2, lc = log(3))
tt <- seq(0, 20, by = 0.25)
set.seed(1)
yy <- fhnPrd(tt, fhnTruth, deriv = FALSE)[[1]][, "y"]
fhnData <- as.datalist(data.frame(name = "y", time = tt, value = yy + rnorm(length(tt), 0, 0.1),
                                  sigma = 0.1, condition = "C"))
fhnSS <- normL2(fhnData, fhnPrd)
fhnMS <- normL2(fhnData, fhnPrd, multipleShootingControl = TRUE)

## ---- fhnfits
set.seed(7)
fhnStarts <- parframe(do.call(rbind, lapply(1:40, function(i)
  fhnTruth + c(rnorm(2, 0, 1.5), rnorm(3, 0, 1)))))
fhnFit <- function(obj, ...)
  as.parframe(mstrust(obj, center = fhnStarts, cores = 4, iterlim = 300, ...))
fhnRuns <- list("single shooting" = fhnFit(fhnSS),
                "multiple shooting, filter" = fhnFit(fhnMS),
                "multiple shooting, merit" = fhnFit(fhnMS, stepControl = list(acceptance = "merit")))

## ---- fhnsave
# the md5 of the script, which the article checks on reading
attr(fhnRuns, "stamp") <- list(md5 = tools::md5sum("fitzhughNagumo.R"), dMod2 = packageVersion("dMod2"))
saveRDS(fhnRuns, "fitzhughNagumo.rds")
