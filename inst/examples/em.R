\donttest{
## Theophylline, 12 subjects: one-compartment model with first-order absorption in closed
## form, log parameters with a normal random intercept per subject, additive error
th <- datasets::Theoph
th$Subject <- as.character(th$Subject)
subjects <- unique(th$Subject)
data <- as.datalist(data.frame(name = "y", time = th$Time, value = th$conc, sigma = NA_real_,
                               condition = th$Subject))
## second derivatives forward over reverse, the default of em()
dm <- c("forward", "reverse", "forward-reverse")
x <- Xt()
g <- Y(c(y = "amt * Ka / (V * Ka - Cl) * (exp(-Cl / V * time) - exp(-Ka * time))"),
       states = "time", modelname = "theoph_obs", compile = FALSE, derivMode = dm,
       outdir = tempdir())
err <- Y(c(y = "sigma_add"), g, attachInput = FALSE, modelname = "theoph_err", compile = FALSE,
         derivMode = dm, outdir = tempdir())
## one condition per subject with the subject dose; the effect formula replaces each log parameter
## by an intercept plus a random effect per subject
doses <- tapply(th$Dose * th$Wt, th$Subject, `[`, 1)[subjects]
trafo <- branch(define(NULL, "x <- x", x = c("Ka", "V", "Cl", "amt", "sigma_add")),
                table = data.frame(ID = subjects, dose = doses, row.names = subjects))
trafo <- insert(trafo, "amt <- dose", dose = dose)
trafo <- insert(trafo, "x <- exp(x)", x = c("Ka", "V", "Cl", "sigma_add"))
p <- P(insert(trafo, Ka + V + Cl ~ 1 + (1 | ID)), method = "explicit", compile = FALSE,
       modelname = "theoph_p", derivMode = dm, outdir = tempdir())
compile(g, err, p, output = tempfile("theoph"))
obj <- normL2(data, g * x * p, errmodel = err)
init <- c(Ka = 0, V = 3, Cl = 0.7, sigma_add = log(0.5))

## FOCEI fit, the estimates and the random effects per subject
fit <- em(obj, randomEffects(p), init = init)
fit
summary(fit)
fit$omega
head(fit$random)

## quadrature started from the FOCEI estimates, and SAEM
fitQ <- em(obj, randomEffects(p), init = fit, method = "quadrature", verbose = FALSE)
fitS <- em(obj, randomEffects(p), init = init, method = "saem", verbose = FALSE)
c(focei = fit$value, quadrature = fitQ$value, saem = fitS$value)
}
\dontrun{
## Ka and Cl with correlated random effects, V with a Laplace one: a second transformation
p2 <- P(insert(trafo, Ka + Cl ~ 1 + (1 | p | ID), V ~ 1 + (1 | gr(ID, dist = laplace()))),
        method = "explicit", compile = FALSE, modelname = "theoph_p2", derivMode = dm,
        outdir = tempdir())
compile(p2, output = tempfile("theoph"))
fit2 <- em(normL2(data, g * x * p2, errmodel = err), randomEffects(p2), init = init)
fit2$omega
fit2$strengths
AIC(fit, fit2)
}
