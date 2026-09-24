# Models with states inside exponentials: symmetryDetection() and
# symmetryReduction(). exp(), exp10(), b^x, tanh(), cosh() and sinh() are
# supported by all engines; equilibrate = TRUE is not.

library(dMod2)

## 1. exp(A) drives the degradation of B ---------------------------------------
f <- eqnvec(A = "kA - kd*A", B = "kB - exp(A)*B")
g <- eqnvec(y = "s*B")
out <- symmetryDetection(f, g, reconstruct = TRUE)   # one scaling of B, kB, s
out <- symmetryDetection(f, g, method = "scaling")   # A cannot scale inside exp()

## 2. a prefactor turns A into a translation -----------------------------------
# A -> A + c, kA -> kA + kd*c, kE -> kE*exp(-c) leaves y unchanged
fE <- eqnvec(A = "kA - kd*A", B = "kB - kE*exp(A)*B")
obs <- symmetryDetection(fE, g, reconstruct = TRUE)
pol <- symmetryDetection(fE, g, method = "polynomial",
                         polynomial = polynomialControl(ansatz = "par", pMax = 1L))

# A and kA may take either sign; declared real, they give a chart that pins A
red <- symmetryReduction(obs, positive = c("B", "kB", "kd", "kE", "s"))
red
symmetryDetection(fE, g, trafo = red$trafo)$identifiable

# a known initial value removes the translation
out <- symmetryDetection(fE, g, trafo = eqnvec(A = "0"))

## 3. other bases and later events ---------------------------------------------
out <- symmetryDetection(eqnvec(A = "kA - kd*A", B = "kB - kE*exp10(A)*B"), g,
                         reconstruct = TRUE)          # kE moves with log(10)

# the analysis starts at the earliest event, so an event at t = 0 anchors it
ev <- eventlist() |>
  addEvent(var = "B", time = 0, value = "0", method = "add") |>
  addEvent(var = "A", time = 5, value = "dA", method = "add")
out <- symmetryDetection(fE, g, events = ev, reconstruct = TRUE)

## 4. Hodgkin-Huxley -----------------------------------------------------------
an <- "0.01*(V + 55)/(1 - exp(-(V + 55)/10))"
bn <- "0.125*exp(-(V + 65)/80)"
am <- "0.1*(V + 40)/(1 - exp(-(V + 40)/10))"
bm <- "4*exp(-(V + 65)/18)"
ah <- "0.07*exp(-(V + 65)/20)"
bh <- "1/(1 + exp(-(V + 35)/10))"
gate <- function(a, b, x) sprintf("(%s)*(1 - %s) - (%s)*%s", a, x, b, x)
hh <- eqnvec(
  V = "(I - gNa*m^3*h*(V - ENa) - gK*n^4*(V - EK) - gL*(V - EL))/C",
  m = gate(am, bm, "m"),
  h = gate(ah, bh, "h"),
  n = gate(an, bn, "n"))

# measuring V leaves the conductances and I relative to C, and only I + gL*EL
obs <- symmetryDetection(hh, eqnvec(y = "V"), reconstruct = TRUE)
out <- symmetryDetection(hh, eqnvec(y = "V"), method = "scaling")

# potentials and the current take either sign
red <- symmetryReduction(obs, positive = c("C", "gNa", "gK", "gL", "m", "h", "n"))
red
symmetryDetection(hh, eqnvec(y = "V"), trafo = red$trafo)$identifiable

## 5. Morris-Lecar: tanh() and cosh() are exponentials -------------------------
ml <- eqnvec(
  V = "(I - gL*(V - VL) - gCa*(1 + tanh((V - V1)/V2))/2*(V - VCa) - gK*w*(V - VK))/C",
  w = "phi*cosh((V - V3)/(2*V4))*((1 + tanh((V - V3)/V4))/2 - w)")
out <- symmetryDetection(ml, eqnvec(y = "V"), reconstruct = TRUE)

## 6. Boltzmann gate with parameters in the exponent ---------------------------
# a common shift of Vh and u, and a common scaling of Vh, u, k and a
bz <- eqnvec(x = "(1/(1 + exp((Vh - u)/k)) - x)/tau", u = "a - b*u")
out <- symmetryDetection(bz, eqnvec(y = "s*x"), reconstruct = TRUE)
red <- symmetryReduction(out)
red
symmetryDetection(bz, eqnvec(y = "s*x"), trafo = red$trafo)$identifiable
