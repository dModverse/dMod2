\donttest{
# The Brusselator rests at x = a, y = b / a, stable while b < 1 + a^2
f <- c(x = "a - (b + 1) * x + x^2 * y", y = "b * x - x^2 * y")
old <- setwd(tempdir())
pSS <- Pimpl(f, deriv2 = TRUE, modelname = "brusselator_ss")
p <- P(c(a = "exp(la)", b = "exp(lb)"), modelname = "brusselator_p")
compile(pSS, p, output = "brusselator")
stab <- constraintStability(pSS, p, onset = -0.5, weight = 10)
# far from the Hopf point nothing, near it the term and its gradient, beyond it the term on the
# unstable steady state
stab(c(la = 0, lb = log(0.5)))$value
stab(c(la = 0, lb = log(1.8)))
stab(c(la = 0, lb = log(2.5)))$value
setwd(old)
}
