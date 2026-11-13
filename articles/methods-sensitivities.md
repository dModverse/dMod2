# Sensitivities and Hessians

The derivatives of the objective come from the sensitivities of the ODE
solution and the chain rule through every function object, without
finite differences. Guide:
[Sensitivities](https://dmodverse.github.io/dMod2/articles/guide-sensitivities.md);
integrators: [cppDE](https://dmodverse.github.io/cppDE/).

| sweep                | returns                        |
|----------------------|--------------------------------|
| forward              | gradient, Gauss-Newton Hessian |
| forward over forward | exact Hessian                  |
| reverse              | gradient                       |
| forward over reverse | exact Hessian                  |

## Forward sensitivities

The sensitivities \\S = \partial x / \partial \theta\\ of \\\dot x =
f(x, \theta)\\, \\x(t_0) = x_0(\theta)\\, solve

\\\dot S = \frac{\partial f}{\partial x}\\ S + \frac{\partial
f}{\partial \theta}, \qquad S(t_0) = \frac{\partial x_0}{\partial
\theta}, \tag{1}\\

\\n_x n\_\theta\\ further equations. cppDE evaluates \\f\\ on dual
numbers, one tangent per parameter, and steps states and tangents
together. The step-size control includes the tangents (`sensErrCon`), so
the sensitivities are as accurate as the states ([forward
sensitivities](https://dmodverse.github.io/cppDE/articles/methods-forward.html)).
An event passes the tangents through the saltation matrix of the jump
([events](https://dmodverse.github.io/cppDE/articles/methods-events.html)).
Every further function object of `g*x*p` multiplies by the Jacobian of
that object, e.g. \\g_x S + g\_\theta\\ for \\y = g(x, \theta)\\.

## Gradient and Hessian of the objective

Each data point adds, up to a constant,

\\\ell_j = r_j^2 + \log \sigma_j^2, \qquad r_j = \frac{f_j -
y_j}{\sigma_j}, \tag{2}\\

the summand of [the data
term](https://dmodverse.github.io/dMod2/articles/methods-objective.html#data-term).
For fixed \\\sigma_j\\ and \\J = \partial r / \partial \theta\\, the
gradient is \\2 J^\top r\\ and the Hessian

\\\nabla^2 \sum_j r_j^2 = 2 J^\top J + 2 \sum_j r_j\\ \nabla^2 r_j .
\tag{3}\\

The Gauss-Newton term \\2 J^\top J\\ is positive semidefinite and needs
first-order sensitivities; the second term matters for large residuals
or strong nonlinearity.

### Gauss-Newton with an error model

For \\\sigma_j(\theta)\\, the Hessian of \\\ell_j\\ in \\(f_j,
\sigma_j)\\ is

\\L_j = \frac{2}{\sigma_j^2} \begin{pmatrix} 1 & -2 r_j \\ -2 r_j & 3
r_j^2 - 1 \end{pmatrix}, \qquad \det L_j = -\frac{4\\(1 +
r_j^2)}{\sigma_j^4} \< 0. \tag{4}\\

The Gauss-Newton Hessian of
[`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md) is
\\\sum_j J_j^\top L_j J_j\\, \\J_j\\ the Jacobian of \\(f_j,
\sigma_j)\\: exact curvature of \\\ell_j\\, no second derivatives of
\\f_j\\ and \\\sigma_j\\. The concave \\\log \sigma_j^2\\ gives \\L_j\\
one negative eigenvalue, so the matrix is in general indefinite along
directions that move a sigma where residuals are small. The exact
Hessian adds the second derivatives of \\f_j\\ and \\\sigma_j\\,
weighted by \\2 r_j / \sigma_j\\ and \\2 (1 - r_j^2) / \sigma_j\\.

## Forward over forward

Differentiating (1) again gives \\n_x n\_\theta (n\_\theta + 1) / 2\\
equations and the exact Hessian, for small models.

## Reverse sweep

The gradient of the objective alone follows from one backward sweep with
the multipliers \\\lambda(t)\\ and the data term \\\ell\\,

\\\dot\lambda = -\left(\frac{\partial f}{\partial x}\right)^{\\\top}
\lambda - \left(\frac{\partial \ell}{\partial x}\right)^{\\\top}, \qquad
\frac{\mathrm{d}(-2 \log L)}{\mathrm{d}\theta} = \int \lambda^\top
\frac{\partial f}{\partial \theta}\\ \mathrm{d}t + \dots, \tag{5}\\

with jumps of \\\lambda\\ at the measurement times. cppDE uses the
*discrete adjoint* of the integrator steps, the exact derivative of the
computed objective, with checkpoints of the forward trajectory ([reverse
mode](https://dmodverse.github.io/cppDE/articles/methods-reverse.html)).
The cost is a few times the value, independent of \\n\_\theta\\. The
value pass solves without sensitivities, and the sweep uses the step
grid of that pass.
[`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md) passes
the sweep through the Jacobian of the root ([Steady
states](https://dmodverse.github.io/dMod2/articles/methods-steady-states.html#derivatives-of-the-root)).

### Accuracy of a reverse gradient

By default the adjoint follows the step grid of the states without error
control.

| option of `optionsReverse` | effect |
|----|----|
| `refine = TRUE` | every step held to the CVODES backward error test under the solve tolerances; a failing step is swept again in substeps |
| `gradtol` | the share of each step in the gradient added to that test as absolute tolerance |

A trust region converges with a gradient accurate relative to the
gradient size [\[1\]](#ref-carter1991global),
[\[2\]](#ref-conn2000trust), so `gradtol` is a fraction of the gradient
scale ([checking the
sweep](https://dmodverse.github.io/cppDE/articles/methods-reverse-checks.html)).

A refined sweep differentiates the trajectory with substeps, not the
computed one. Under `refine` it therefore returns, per seed, the
first-order change of the seeded functional between the two: over the
steps, the adjoint at the step end times the difference of the end
states, plus the seeds times the difference of the observed values
([value
correction](https://dmodverse.github.io/cppDE/articles/methods-reverse-checks.html#value-correction)).
The [`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md)
seed gives the change of the objective, a second seed the change of
`chi2`. The reverse objective adds both, so value and gradient belong to
one function, and a value-only call also sweeps. The corrected value is
not necessarily closer to the exact one; accepted steps get no
correction.

## Forward over reverse

Forward tangents, one per parameter, ride on the backward sweep, which
returns one Hessian column per tangent. Cost grows like a forward solve;
nothing of order \\n\_\theta^2\\ is integrated or stored ([forward over
reverse](https://dmodverse.github.io/cppDE/articles/methods-forward-reverse.html)).

### Second derivatives along directions

FOCEI ([Nonlinear mixed-effects
models](https://dmodverse.github.io/dMod2/articles/guide-nlme.md)) needs
traces \\\operatorname{tr}(M\\ \partial H_i / \partial \psi)\\ for the
volume term and the joint Hessian along random-effect directions. Random
effects of different groups never share a condition, so the \\k\\-th
random effect of all groups shares one tangent \\v_k\\, and one forward
pass along \\K\\ tangents serves all groups. A backward sweep with seed
\\w\\ on the predictions and sigmas of the data rows, and \\\dot w_k\\
along tangent \\k\\, returns \\J^\top w\\ and

\\\sum_r w_r\\ \nabla^2 f_r\\ v_k + J^\top \dot w_k \tag{6}\\

for every \\k\\, \\f_r\\ prediction or sigma of row \\r\\, \\J\\ their
Jacobian. With the gradient of the row loss as seed and \\L_r J_r v_k\\
as tangent derivative, \\L_r\\ the Hessian (4), the sweep gives
\\\Phi_i^{\eta\eta}\\, \\\Phi_i^{\eta\psi}\\ and the mode shift
\\\mathrm{d}\hat\eta / \mathrm{d}\psi = -(\Phi_i^{\eta\eta})^{-1}
\Phi_i^{\eta\psi}\\. The trace gradients are linear in \\M\\: seed \\j\\
is \\2 L_r J\_{r,j}\\, with tangent derivative the derivative of \\L_r\\
contracted with \\J\_{r,j}\\ and \\J\_{r,k}\\. With \\U\_{jk}\\ slice
\\k\\ of the result for seed \\j\\,

\\\nabla \operatorname{tr}(M H_i) = \sum\_{j,k} M\_{jk}\\ U\_{jk}
\tag{7}\\

for every symmetric \\M\\: \\K\\ seeds for any number of matrices.
Hessian rows and trace seeds go back in one sweep, every cotangent with
\\K + 1\\ blocks, the seed and the derivative of the seed along each
tangent.

Reference:
[`odemodel()`](https://dmodverse.github.io/dMod2/reference/odemodel.md),
[`Xs()`](https://dmodverse.github.io/dMod2/reference/Xs.md),
[`Y()`](https://dmodverse.github.io/dMod2/reference/Y.md),
[`P()`](https://dmodverse.github.io/dMod2/reference/P.md),
[`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md),
[`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md)

## References

\[1\]

R. G. Carter, “On the global convergence of trust region algorithms
using inexact gradient information,” *SIAM Journal on Numerical
Analysis*, vol. 28, no. 1, pp. 251–265, 1991, doi:
[10.1137/0728014](https://doi.org/10.1137/0728014).

\[2\]

A. R. Conn, N. I. M. Gould, and P. L. Toint, *Trust-region methods*.
Philadelphia: SIAM, 2000. doi:
[10.1137/1.9780898719857](https://doi.org/10.1137/1.9780898719857).
