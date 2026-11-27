# Steady states

A steady state is a state \\x\\ at which the right-hand side \\f\\ of
the ODE \\\dot x = f(x, p)\\ of states \\x\\ and parameters \\p\\
vanishes, a root of \\f\\. A model at rest before a stimulus starts from
a steady state.
[`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md)
computes it numerically,
[`steadyStates()`](https://dmodverse.github.io/dMod2/reference/steadyStates.md)
in closed form. Guide: [Setting up
models](https://dmodverse.github.io/dMod2/articles/guide-model.html#steady-state-constraints).

## Conserved quantities

With compartment volumes \\V_i\\, the volume of the compartment that
holds species \\i\\, stoichiometric coefficients \\S\_{ri}\\, the net
number of molecules of species \\i\\ that reaction \\r\\ produces,
negative where it consumes them, the entries of the stoichiometry matrix
\\S\\, and rates \\v_r\\ in a frame of volume \\V\_{c(r)}\\, the
right-hand side is

\\\dot x_i = f_i(x, p) = \frac{1}{V_i} \sum_r S\_{ri}\\ V\_{c(r)}\\
v_r(x, p). \tag{1}\\

Every \\w\\ with \\S w = 0\\ conserves the amount

\\\sum_i w_i V_i x_i \tag{2}\\

by (1), a conserved quantity;
[`conservedQuantities()`](https://dmodverse.github.io/dMod2/reference/conservedQuantities.md)
lists these amounts. Then \\f(x, p) = 0\\ does not determine \\x\\:
every conservation class, the states with given totals, has separate
roots.

## Pseudo-transient continuation

For an `eqnlist`,
[`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md) adds
the conservation laws \\C(p)\\ x = T\\ of the stoichiometry, which set
the conserved quantities to totals \\T\\ that enter as parameters:

\\G(x, p) = \begin{pmatrix} f(x, p) \\ C(p)\\ x - T \end{pmatrix} = 0 .
\tag{3}\\

Equations given as an `eqnvec` or character vector have no stoichiometry
and are solved as they are, so they have to determine \\x\\ themselves.
With `flow = TRUE`, the default for an `eqnlist`, pseudo-transient
continuation solves (3) by implicit Euler steps of \\\dot x = f(x, p)\\
on the conservation manifold, the affine subspace \\C x = T\\ that a
trajectory never leaves, with a pseudo-time step that grows as the
residual falls. It follows the flow far from the root and becomes
Newton’s method near it. Newton steps may also reach an unstable root,
which the [stability test](#stability) catches. With `flow = FALSE`, the
default otherwise, the equations need not be the right-hand side of an
ODE: the solver takes damped Newton steps, without stability test, and
returns any regular root. A state no reaction changes keeps the given
value. The solver is cppDE’s `ptc()`, called from C++ ([step,
convergence and stability
tests](https://dmodverse.github.io/cppDE/articles/methods-ptc.html)).

Every root starts from the same guess, a function of parameters and
controls, never from an earlier root, so the objective is the same
function of the parameters in every call. Warm starts would change roots
in the last digits, which the ODE step-size control magnifies in the
objective.

## Derivatives of the root

From \\G(x(p), p) = 0\\ for all \\p\\, the implicit function theorem
gives

\\\frac{\partial G}{\partial x}\\ \frac{\partial x}{\partial p} =
-\frac{\partial G}{\partial p}, \tag{4}\\

solvable for \\\partial x / \partial p\\ at a regular root. A root is
regular where the Jacobian \\J = \partial f / \partial x\\ restricted to
the conservation manifold, \\N^\top J N\\, is nonsingular, with the
columns of \\N\\ a basis of the kernel \\\ker C\\, the vectors \\v\\
with \\C v = 0\\, the directions within the conservation class. A
regular root is isolated in the conservation class. Where the matrix is
singular,
[`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md) takes
the pseudo-inverse and warns. Differentiating (4) again gives the second
derivatives, `Pimpl(derivMode = "forward-forward")`, with the same
matrix. A reverse sweep applies the transpose of the solution \\\partial
x / \partial p\\ of (4) ([Sensitivities and
Hessians](https://dmodverse.github.io/dMod2/articles/methods-sensitivities.html#reverse-sweep)).

## Stability

A resting state is a stable steady state, the state of the model before
a stimulus. A root is a resting state only if every eigenvalue
\\\lambda_k\\ of \\J = \partial f / \partial x\\ on the conservation
manifold has negative real part. At a root whose leading eigenvalue has
real part above the larger of \\10^{-8} s\\, \\s\\ the largest
eigenvalue modulus on the conservation manifold and \\10^{-8} s\\ the
precision of the eigenvalues, and the stability tolerance
\\\lambda\_{tol}\\, `stabilityTol` with default 0, the solver continues
along the leading eigenvector
([cppDE](https://dmodverse.github.io/cppDE/articles/methods-ptc.html#stability)).
Where neither direction leads to a stable root, the next start is tried;
the evaluation fails only when no start reaches one.
`controlsPTC = list(stability = FALSE)` returns the unstable root.

[`constraintStability()`](https://dmodverse.github.io/dMod2/reference/constraintStability.md)
adds a term that rises on a ramp below \\\lambda\_{tol}\\, the sum of
squares of the residuals

\\\begin{aligned} \rho_k &= \sqrt{\kappa} \begin{cases} r_k^{3/2}, & r_k
\le 1, \\ 1 + \tfrac{3}{2}(r_k - 1), & r_k \> 1, \end{cases} \\ r_k &=
\frac{\max\bigl(0,\\ \mathrm{Re}\\\lambda_k -
\lambda_0\bigr)}{\lambda\_{tol} - \lambda_0}, \qquad \lambda_0 =
\lambda\_{tol} - \max(d, m s), \end{aligned} \tag{5}\\

with weight \\\kappa\\ and the width \\\max(d, m s)\\ of the ramp: at
least the absolute `width` \\d\\, default 0, and at least the `margin`
\\m\\ relative to the scale, default \\10^{-6}\\, the precision of the
eigenvalues. The term is zero on a root whose leading real part lies
below \\\lambda_0\\, however slowly a perturbation decays there. Each
eigenvalue whose real part reaches \\\lambda\_{tol}\\ adds \\\kappa\\,
and beyond, the term grows quadratically. Steady states that agree in
the parameters, resting states and start values they read count once.
For a simple eigenvalue with right and left eigenvectors \\u_k\\,
\\v_k\\ of \\J\\, \\v_k^{\mathsf H} u_k \ne 0\\, and (4),

\\ \frac{\mathrm{d}\lambda_k}{\mathrm{d}p} = \frac{v_k^{\mathsf H}
\bigl(\partial_p J + \partial_x J\\ \mathrm{d}x^\*/\mathrm{d}p\bigr)
u_k} {v_k^{\mathsf H} u_k}, \tag{6} \\

and \\s = \|\lambda_m\|\\, \\\lambda_m\\ the eigenvalue of largest
modulus, has the derivative

\\\frac{\mathrm{d}s}{\mathrm{d}p} = \frac{\mathrm{Re}(\bar\lambda_m\\
\mathrm{d}\lambda_m / \mathrm{d}p)}{s}\\

by (6), continuous except where the eigenvalue of largest modulus
changes. Eigenvalues that coincide within the error of the eigenvalues,
as a defective double eigenvalue of a chain of states with equal rates,
have no derivative each; each takes the mean \\\mathrm{tr}(P\\
\mathrm{d}J/\mathrm{d}p) / n_c\\ over the cluster of \\n_c\\ eigenvalues
with spectral projector \\P\\, which is (6) for a simple eigenvalue and
exact for a cluster that stays degenerate. Where two real eigenvalues on
the ramp meet and become a complex pair, the derivatives of both diverge
and the gradient of (5) has a kink. Both need
`Pimpl(derivMode = "forward-forward")`. The term is a model assumption,
not data: profiles with `stop = "data"` and likelihood ratios leave it
out.

The Gauss-Newton Hessian of (5),

\\\begin{aligned} &2 \sum\_{k \\:\\ \mathrm{Re}\\\lambda_k \> \lambda_0}
\frac{\mathrm{d} \rho_k}{\mathrm{d}p} \frac{\mathrm{d}
\rho_k}{\mathrm{d}p}^{\mathsf T}, \\ &\frac{\mathrm{d}
\rho_k}{\mathrm{d}p} = \frac{3}{2} \sqrt{\kappa}\\ \min(r_k, 1)^{1/2}\\
\frac{\mathrm{d} r_k}{\mathrm{d}p}, \end{aligned} \tag{7}\\

vanishes at the onset with the term and the gradient and is continuous
at \\r_k = 1\\. Outside the ramp the term and (7) vanish, so the model
of a trust-region step sees the boundary only from within the ramp; on a
ramp narrow against the change of \\\mathrm{Re}\\\lambda_k\\ over a
step, the steps cross it and are rejected, and the iterates creep along
the boundary. A `width` \\d\\ in the units of the rates, wide against
that change, avoids it and leaves a resting state untouched whose
slowest mode decays faster than \\d\\; the ramp then does not move with
\\s\\. A ramp \\m s\\ wide gives \\r_k\\ as accurate as the root
relative to that width. The eigenvalues above the onset, their
conditions, half-lives and the states with the largest entries of the
right eigenvectors are the attribute `"eigenvalues"` of a call of the
term.

## Analytical steady-state constraints

[`steadyStates()`](https://dmodverse.github.io/dMod2/reference/steadyStates.md)
solves \\f(x, p) = 0\\ in closed form
[\[1\]](#ref-rosenblatt2016customized), [\[2\]](#ref-kaschek2019dmod).
The balance of state \\i\\ is the equation \\V_i\\ \dot x_i = 0\\ in
amounts, a sum of flux terms \\S\_{ri} V\_{c(r)} v_r\\, so that a
quantity conserved across compartments is a dependent balance. First, a
state whose balance is linear in the state, with all terms of the state
on one side, is solved from that balance, influx over outflux rate,
unless the rate constants of the solution leave a cycle without a cut. A
cycle is a closed path of the graph of states with an edge from each
state to every state in the balance of that state; a cut is a state on a
cycle that stays a free parameter. The states of conserved quantities
are left to the cuts. Then the cycles are cut, each at a state of a
conserved quantity, which stays free, or at a state whose balance is
solved for a rate constant common to the flux terms of one side. The
cuts are ranked as by [\[1\]](#ref-rosenblatt2016customized): conserved
quantities first, then sides whose flux terms occur in no other balance,
fewer balances changed, fewer terms, more cycles through the state. Then
each unknown, a state or a kinetic parameter, is solved from a
combination of balances such that the solution is a ratio of sums of
positive terms; the unknown a balance is solved for is the pivot. The
solution is thus positive for positive parameters, also with
`solveQuadratic = TRUE`, which takes the positive root of a balance
quadratic in the state. A state whose balance is spent on a kinetic
parameter stays free, as the parameter of the steady-state value of the
state. States that are zero in every steady state are set to 0. A
conserved quantity keeps one state free. Volumes are never solved for.

Only if no solution of this form exists, a side whose flux terms have
separate rate constants \\k_0, \dots, k_m\\ is a cut as well. The fluxes
\\v_j\\ of the side get free ratios \\r\_{x,j} = v_j / v_0 \> 0\\, named
`r_<state>_<j>`; \\k_0\\ is solved, each \\k_j\\ follows from
\\r\_{x,j}\\ and \\k_0\\. The ratios are new parameters of the
transformation.

The result, one equation per solved unknown, replaces identity entries
of a parameter transformation, the map from outer to inner parameters
that [`P()`](https://dmodverse.github.io/dMod2/reference/P.md) builds:

``` r

steady <- steadyStates(reactions)
parameters <- getParameters(reactions)
trafo <- define(NULL, "x <- x", x = parameters)
trafo <- insert(trafo, "x <- y", x = names(steady), y = steady)
p <- P(trafo)
```

| argument | effect |
|----|----|
| `givenCQs` | conserved quantities, default from [`customTotals()`](https://dmodverse.github.io/dMod2/reference/customTotals.md) or derived |
| `neglect` | states and parameters never solved for |
| `priority` | states solved first, parameters preferred as pivots |
| `forcings` | forcings, states given as external inputs, held at zero |
| `solveQuadratic` | allow the positive root of a quadratic balance |
| `file` | save the result with [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) |

Without a solution of this form,
[`steadyStates()`](https://dmodverse.github.io/dMod2/reference/steadyStates.md)
returns `0`. The solver is a Python module (numpy, sympy, scipy) called
through reticulate.

[`P()`](https://dmodverse.github.io/dMod2/reference/P.md) or
[`Pexpl()`](https://dmodverse.github.io/dMod2/reference/Pexpl.md)
evaluate the closed form explicitly with symbolic derivatives: no root
search, no convergence failure, no solver tolerance in the objective,
one parameter less per solved equation. Rosenblatt et al. report more
reliable and faster optimization than with common numerical
implementations [\[1\]](#ref-rosenblatt2016customized).
[`steadyStates()`](https://dmodverse.github.io/dMod2/reference/steadyStates.md)
does not test stability;
[`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md) does
([Stability](#stability)).
[`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md) is
needed where no closed form of positive sums exists.

Reference:
[`conservedQuantities()`](https://dmodverse.github.io/dMod2/reference/conservedQuantities.md),
[`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md),
[`steadyStates()`](https://dmodverse.github.io/dMod2/reference/steadyStates.md),
[`constraintStability()`](https://dmodverse.github.io/dMod2/reference/constraintStability.md)

## References

\[1\]

M. Rosenblatt, J. Timmer, and D. Kaschek, “Customized steady-state
constraints for parameter estimation in non-linear ordinary differential
equation models,” *Frontiers in Cell and Developmental Biology*, vol. 4,
p. 41, 2016, doi:
[10.3389/fcell.2016.00041](https://doi.org/10.3389/fcell.2016.00041).

\[2\]

D. Kaschek, W. Mader, M. Fehling-Kaschek, M. Rosenblatt, and J. Timmer,
“Dynamic modeling, parameter estimation, and uncertainty analysis in R,”
*Journal of Statistical Software*, vol. 88, no. 10, pp. 1–32, 2019, doi:
[10.18637/jss.v088.i10](https://doi.org/10.18637/jss.v088.i10).
