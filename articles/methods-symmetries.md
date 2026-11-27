# Structural identifiability and Lie symmetries

[`symmetryDetection()`](https://dmodverse.github.io/dMod2/reference/symmetryDetection.md)
finds the Lie symmetries of the prediction with the Python package
[symident](https://dmodverse.github.io/symident/); a gauge or
[`symmetryReduction()`](https://dmodverse.github.io/dMod2/reference/symmetryReduction.md)
removes them. Guide:
[Symmetries](https://dmodverse.github.io/dMod2/articles/guide-symmetries.md).

## Structural identifiability

The model is

\\\dot x = f(x, u(t), \theta), \qquad x(t_0) = x_0(\theta), \qquad y =
g(x, \theta), \tag{1}\\

with states \\x\\, input \\u\\, observables \\y\\ and parameters
\\\theta \in \Theta \subseteq \mathbb{R}^d\\ of any kind, and prediction
\\y(t; \theta) = g(x(t), \theta)\\. It is *locally structurally
identifiable* at \\\theta\\ if no nearby \\\tilde\theta \ne \theta\\ has
the same prediction. The property depends on the equations alone, as if
every observable were measured exactly and continuously, and the answer
is the same for every \\\theta\\ outside a set of measure zero, the zero
set of a nonzero polynomial: for generic \\\theta\\. A structural
non-identifiability makes the likelihood exactly flat along a curve,
seen as a flat profile [\[1\]](#ref-raue2009structural).

## Lie symmetries

A *Lie symmetry* is a one-parameter group \\\Phi\_\varepsilon\\ on
\\\Theta\\, maps with \\\Phi_0\\ the identity and \\\Phi_a \circ \Phi_b
= \Phi\_{a + b}\\, with \\y(t; \Phi\_\varepsilon(\theta)) = y(t;
\theta)\\ for all \\t\\ and small \\\varepsilon\\. The orbit of
\\\theta\\, the curve \\\varepsilon \mapsto \Phi\_\varepsilon(\theta)\\,
is a continuous family of parameter values with the same prediction. The
orbits are the integral curves, with tangent \\\xi\\ at every point, of
the *generator*

\\\begin{aligned} X &= \sum\_{i=1}^d \xi_i(\theta)\\
\frac{\partial}{\partial \theta_i}, \\ \xi(\theta) &= \frac{\partial
\Phi\_\varepsilon(\theta)}{\partial \varepsilon}\Big\|\_{\varepsilon =
0}. \end{aligned} \tag{2}\\

The components \\\xi_i\\ are the rates of change of the coordinates
along an orbit. An *invariant*, a function \\I(\theta)\\ with \\X(I) =
0\\, is constant on every orbit, and the prediction depends on
\\\theta\\ only through invariants. With \\r\\ generators that span all
symmetries and are independent at generic \\\theta\\, the observables
determine \\d - r\\ parameter combinations.

*Scalings*, \\\xi_i(\theta) = w_i \theta_i\\ with integer weights
\\w_i\\, multiply coordinate \\i\\ by \\e^{w_i \varepsilon}\\. They are
the most common symmetries of reaction networks; on the log scale their
orbits are straight lines. Two kinds occur in almost every signaling
model:

- Observable `s * (A + B) + off`: concentrations of the moiety times
  \\\lambda\\; the scale `s` and the rate constant of every term of
  second order in the moiety, or of first order in it in the equation of
  another state, divided by \\\lambda\\. The generator moves the total,
  the scales and these rate constants together.
- Unobserved intermediate: level and production rate times \\\lambda\\,
  the rate it drives divided by \\\lambda\\.

General symmetries also occur, such as a Möbius transformation of an
input amplitude in an NFκB model [\[2\]](#ref-merkt2015higher); symident
finds generators with rational components.

## Detection

With \\z = (x, \theta, t)\\, \\\dot\theta = 0\\, \\\dot t = 1\\, the
model is the vector field \\F = f\\ \partial_x + \partial_t\\, started
at \\z_0(\theta) = (x_0(\theta), \theta, t_0)\\. The *Lie derivatives*
of the observables,

\\L_F^0 g = g, \qquad L_F^{k+1} g = F\bigl(L_F^k g\bigr), \tag{3}\\

at \\z_0\\ are the time derivatives of the prediction at \\t_0\\, so an
analytic prediction is

\\y(t; \theta) = \sum\_{k \ge 0} \frac{(t - t_0)^k}{k!}\\ L_F^k
g\bigl(z_0(\theta)\bigr). \tag{4}\\

Collected up to order \\N\\ they form \\\Psi_N(\theta)\\, with the
observability-identifiability matrix (OIM)

\\\mathrm{OIM}\_N(\theta) = \frac{\partial \Psi_N(\theta)}{\partial
\theta}. \tag{5}\\

With free initial values among the parameters, one matrix tests both
observability, the initial states follow from the observables, and
identifiability, the parameters do.

\\\Psi_N\\ is constant on every orbit,
\\\Psi_N(\Phi\_\varepsilon(\theta)) = \Psi_N(\theta)\\, and the
derivative at \\\varepsilon = 0\\ is, by the chain rule and (2),

\\\begin{aligned}
\frac{\mathrm{d}}{\mathrm{d}\varepsilon}\Psi_N\bigl(\Phi\_\varepsilon(\theta)\bigr)\Big\|\_{\varepsilon
= 0} &= \mathrm{OIM}\_N(\theta)\\ \xi(\theta) \\ &= 0. \end{aligned}
\tag{6}\\

At generic \\\theta\\ and large enough \\N\\ the kernel of
\\\mathrm{OIM}\_N\\ consists exactly of the generators, \\d -
\operatorname{rank} \mathrm{OIM}\_N\\ counts them, and full rank means
identifiable.

For scalings, (6) reads \\\mathrm{OIM}\_N(\theta)\\
\mathrm{diag}(\theta)\\ w = 0\\ for every \\\theta\\. With rational
right-hand sides, observables and initial values it holds where every
term of every equation has the weighted degree of the left-hand side,
the degree of a monomial with each variable counted by the weight of
that variable, with weights on the states as well: a linear system in
the weights, whose integer solution lattice symident computes exactly. A
scaling that holds only on an invariant set, such as a start at rest,
fails this test. It is among the remaining kernel vectors, the *general
directions*, generators that are no scaling, with rational closed forms
from interpolation ([symident
method](https://dmodverse.github.io/symident/method.html)).

symident computes (5) exactly over finite fields \\\mathbb{F}\_p\\ at a
fixed base point, the first primes 2, 3, 5, … as parameter values, and
takes the largest rank over several primes and a random point. The rank
modulo \\p\\ at a point can only fall below the generic value; the
checks make a lower one unlikely but do not exclude it. The Lie order
\\N\\, the highest order in \\\Psi_N\\, is sufficient when no higher
order raises the rank of \\\mathrm{OIM}\_N\\; symident proves this where
it can. The level of a direction is `"exact"` where the closed form is
proven over the rationals, `"verified"` or `"checked"` where it is
checked modulo a large prime at random points, with or without an error
bound, and `"open"` without closed form; by default
(`reconstruct = FALSE`) a general direction is returned as the support,
the coordinates with \\\xi_i \not\equiv 0\\ that the generator moves.

### Steady-state starts, inputs and conditions

| setting | effect |
|----|----|
| `equilibrate = TRUE` | initial state a root of \\f(x, 0, \theta) = 0\\, inputs at 0 |
| `reduceCQ = TRUE` | one species eliminated per conserved moiety; freedom reported on the totals |
| `forcings` | input states, 0 at rest |
| `conditions` | shared parameters; a direction is non-identifiable only if it is so in every condition |

symident writes an observable \\a \log g + c\\ through \\g\\. A
clock-time factor of \\g\\ such as \\e^{-\tau t}\\ is no function of the
states of a model at rest, so it is left out of the observables passed
to
[`symmetryDetection()`](https://dmodverse.github.io/dMod2/reference/symmetryDetection.md);
this is valid where the parameters of that factor take no part in a
symmetry.

## Gauge and reduction

Along every orbit the data term and the profiles of the data term are
flat; a prior such as
[`constraintL2()`](https://dmodverse.github.io/dMod2/reference/constraintL2.md)
is not invariant and tilts them. The Hessian is flat in the generator
directions at the optimum, in the Gauss-Newton approximation and for
scalings on the log scale. So symmetries are removed before the fit. A
*gauge* fixes one coordinate per scaling; each orbit meets the fixed
values in one point with the same prediction. With \\W\\ the weights of
the scalings on the fixed coordinates, one row per generator, the gauge
removes every scaling if

\\\det W \ne 0 . \tag{7}\\

[`symmetryReduction()`](https://dmodverse.github.io/dMod2/reference/symmetryReduction.md)
turns the generators into a parameter transformation: for scalings a
*section*, a set that meets every orbit in one point, here by fixing one
coordinate per direction, for general directions a chart on the set of
orbits, invariants \\q_k\\ as new parameters, returned only where valid
on the domain declared by `positive`
([reduction](https://dmodverse.github.io/symident/method/reduction.html)).

Reference:
[`symmetryDetection()`](https://dmodverse.github.io/dMod2/reference/symmetryDetection.md),
[`symmetryReduction()`](https://dmodverse.github.io/dMod2/reference/symmetryReduction.md)

## References

\[1\]

A. Raue *et al.*, “Structural and practical identifiability analysis of
partially observed dynamical models by exploiting the profile
likelihood,” *Bioinformatics*, vol. 25, no. 15, pp. 1923–1929, 2009,
doi:
[10.1093/bioinformatics/btp358](https://doi.org/10.1093/bioinformatics/btp358).

\[2\]

B. Merkt, J. Timmer, and D. Kaschek, “Higher-order Lie symmetries in
identifiability and predictability analysis of dynamic models,”
*Physical Review E*, vol. 92, no. 1, p. 012920, 2015, doi:
[10.1103/PhysRevE.92.012920](https://doi.org/10.1103/PhysRevE.92.012920).
