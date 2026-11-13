# Marginal likelihood of mixed-effects models

[`em()`](https://dmodverse.github.io/dMod2/reference/em.md) maximizes
the marginal likelihood of nonlinear mixed-effects models. Guide:
[Nonlinear mixed-effects
models](https://dmodverse.github.io/dMod2/articles/guide-nlme.md).

## Marginal likelihood

Group \\i\\ of \\G\\ has random effects \\\eta_i \in \mathbb{R}^K\\,
inner parameters \\\phi_i = h(\theta, \eta_i)\\ and the observation
model

\\\begin{aligned} y\_{ij} &= g\bigl(x(t\_{ij};\\ \phi_i)\bigr) +
\varepsilon\_{ij}, \\ \varepsilon\_{ij} &\sim \mathcal{N}\bigl(0,
\sigma^2(\phi_i, t\_{ij})\bigr), \qquad \eta_i \sim \mathcal{N}(0,
\Omega), \end{aligned} \tag{1}\\

with population parameters \\\psi = (\theta, \Omega, \sigma)\\, and

\\\ell(\psi) = \sum\_{i=1}^G \log L_i(\psi), \qquad L_i(\psi) = \int
p(y_i \mid \eta_i, \psi)\\ p(\eta_i \mid \Omega)\\ \mathrm{d}\eta_i .
\tag{2}\\

The integrand is `normL2(data, g * x * p, errmodel = e)`, with
sensitivities in \\\eta_i\\ by the chain rule, plus the prior \\\sum_i
\eta_i^\top \Omega^{-1} \eta_i + G \log \lvert \Omega \rvert\\ from
`randomEffects(p)`. Conditions without a group add their \\-2 \log L\\
unchanged. `normT(..., df = 4)` gives Student-t residuals.

## EM algorithm

With the posterior \\q_i(\eta_i) = p(\eta_i \mid y_i, \psi^{(t)})\\,
Jensen’s inequality gives [\[1\]](#ref-dempster1977maximum)

\\\log L_i(\psi) \ge \underbrace{\mathbb{E}\_{q_i}\bigl\[\log p(y_i,
\eta_i \mid \psi)\bigr\]}\_{Q_i(\psi \mid \psi^{(t)})} -
\mathbb{E}\_{q_i}\bigl\[\log q_i(\eta_i)\bigr\], \tag{3}\\

with equality at \\\psi = \psi^{(t)}\\. EM alternates forming \\Q =
\sum_i Q_i\\ and maximizing it, so \\\ell\\ never decreases; near the
optimum it converges linearly at the fraction of missing information.
ECM [\[2\]](#ref-meng1993maximum) splits the M-step: \\\Omega\\ in
closed form, structural parameters by optimization. [Sparse random
effects](https://dmodverse.github.io/dMod2/articles/methods-sparse-random-effects.html#hyperparameters)
minimizes \\-2Q\\ instead. The methods replace the posterior, which has
no closed form, by:

| method | posterior and expectation | M-step | reported \\-2 \log L\\ |
|----|----|----|----|
| `"focei"` | Gaussian at the mode | trust region, Laplace | Laplace |
| `"quadrature"` | adaptive sparse grid at the mode | trust region, frozen nodes | quadrature |
| `"laplaceEM"` | Laplace moments, nodes of degree 3 | ECM | quadrature, level \\K + 2\\ |
| `"saem"` | Metropolis draws, averaged | closed form, averaged Newton | quadrature, level \\K + 2\\ |

## FOCEI

### Laplace approximation

With \\\Phi_i(\eta) = -2 \log p(y_i \mid \eta, \psi) - 2 \log p(\eta
\mid \Omega)\\, mode \\\hat\eta_i\\ and curvature \\H_i = \nabla^2\_\eta
\Phi_i(\hat\eta_i)\\, Laplace’s approximation of (2) is \\-2 \log L_i
\approx \Phi_i(\hat\eta_i) + \log \lvert H_i/2 \rvert - K \log(2\pi)\\.
The prior omits \\K \log(2\pi)\\, so FOCEI minimizes

\\F(\psi) = \sum_i \Bigl\[ \Phi_i(\hat\eta_i) + \log \lvert H_i/2 \rvert
\Bigr\] \tag{4}\\

plus the full \\-2 \log L\\ of the conditions without a group. \\\log
\lvert H_i/2 \rvert\\ is the volume term. \\H_i\\ includes the
derivatives of \\\sigma\\ in \\\eta\\, the interaction of FOCEI
[\[3\]](#ref-wang2007derivation),
[\[4\]](#ref-pinheiro1995approximations). Conditional standard errors of
the modes are the square roots of the diagonal of \\2 H_i^{-1}\\. The
error has either sign and order \\K / n_i\\ for \\n_i\\ data points in
the group.

### Bilevel optimization

Each evaluation, in C++:

- **Modes.** One trust region per group, Gauss-Newton with SR1
  correction, started from the last modes moved to first order along the
  outer step, stopped at `innerControl$steptol`; eigenvalue floor
  \\10^{-10} \operatorname{tr}(H_i)/K\\; trial points of all groups in
  one batched prediction.
- **Gradient.** By the envelope theorem the mode shift drops out:
  \\\nabla\_\psi F = \sum_i \nabla\_\psi \Phi_i \big\|\_{\hat\eta_i} +
  \sum_i \mathrm{d} \log \lvert H_i \rvert / \mathrm{d}\psi\\, with a
  Newton correction for the gradient left by the inner search.
- **Model Hessian.** Gauss-Newton Schur complement \\\sum_i
  \[\Phi_i^{\psi\psi} - \Phi_i^{\psi\eta} (\Phi_i^{\eta\eta})^{-1}
  \Phi_i^{\eta\psi}\]\\, plus the volume-term Hessian with second-order
  sensitivities, plus SR1 from the last three accepted gradients.
- **Stop.** Model gain below `trustControl$mtol` in the region `rmax`,
  or five rejected trials flat within `ftol`.
- **Sparsity.** `sparse = TRUE`: sparse outer Hessians, Schur complement
  and volume term of a group local to the parameters of the group.

### Volume term

For Gauss-Newton \\H_i\\, from second-order sensitivities,

\\\frac{\mathrm{d} \log \lvert H_i \rvert}{\mathrm{d} \psi_k} =
\operatorname{tr}\Bigl(H_i^{-1} \frac{\partial H_i}{\partial
\psi_k}\Big\|\_\eta\Bigr) + \sum_l \operatorname{tr}\Bigl(H_i^{-1}
\frac{\partial H_i}{\partial \eta_l}\Bigr) \frac{\mathrm{d}
\hat\eta_l}{\mathrm{d} \psi_k}, \qquad
\frac{\mathrm{d}\hat\eta}{\mathrm{d}\psi} = -(\Phi_i^{\eta\eta})^{-1}
\Phi_i^{\eta\psi}, \tag{5}\\

| `control` | volume term |
|----|----|
| `derivMode = "forward-reverse"` (default) | traces from one forward pass and a backward sweep ([Sensitivities](https://dmodverse.github.io/dMod2/articles/guide-sensitivities.md)); exact gradient, curvature by SR1 |
| `derivMode = "forward-forward"` | `deriv2 = TRUE` sensitivities; exact Hessian in the outer model |
| `secondOrderCorrection = FALSE` | derivative dropped; biased estimates |

The reported Hessian is exact (forward over forward) where the model has
`deriv2`, as [`summary()`](https://rdrr.io/r/base/summary.html) states.
For \\\sigma\\ depending on \\\eta\\, quadrature gives better standard
errors of \\\Omega\\.

## Adaptive sparse-grid quadrature

Laplace biases \\\Omega\\ for skewed posteriors, as for oral doses
sampled near absorption. The cubature \\\int e^{-\Phi_i/2}\\
\mathrm{d}\eta \approx \sum_b W_b\\ e^{-\Phi_i(\eta_b)/2}\\ removes that
error level by level. It combines the nested Genz-Keister rules
[\[5\]](#ref-genz1996fully), 1, 3, 9, 19 and 35 points exact to degree
1, 5, 15, 29 and 51, by Smolyak’s sparse grid
[\[6\]](#ref-smolyak1963quadrature), [\[7\]](#ref-heiss2008likelihood),

\\A(L, K) = \sum\_{\lvert \mathbf{i} \rvert \le L} \Delta^{i_1} \otimes
\cdots \otimes \Delta^{i_K}, \tag{6}\\

\\\Delta^l\\ the difference of the one-dimensional rules of levels \\l\\
and \\l - 1\\. Level \\L\\ is exact for total degree \\2(L - K) + 1\\
against the Gaussian, with nodes polynomial in \\K\\ (for \\K = 3\\,
levels 3 to 7 have 1, 7, 19, 45, 117); level \\K\\ is the single Laplace
node. The grid sits on the FOCEI Gaussian \\\mathcal{N}(\hat\eta_i, 2
H_i^{-1})\\ and moves with \\\psi\\ [\[8\]](#ref-clausen2009maximum).
Negative weights make \\\log L_i\\ a signed log-sum-exp; a nonpositive
sum, an unresolved posterior, fails the evaluation, naming the group.
Gradient and Hessian are exact.

`method = "quadrature"` starts from the converged FOCEI fit and
alternates:

- **Rebuild.** Grids at the current \\\psi\\, at the lowest level (never
  falling) whose next level changes \\-2 \log L_i\\ by at most `tol`
  over the number of groups.
- **Step.** Trust region over all parameters on frozen nodes. For
  mu-referenced fixed effects, \\\phi_i = D_i \beta + \eta_i\\, nodes
  are frozen in \\\phi\\, need no sensitivities and hold over larger
  steps.
- **Stop.** A Newton step on the rebuilt nodes would gain less than
  `tolDecrement` or less than the nodes resolve.

At the estimate, importance sampling (\\t_5\\ proposal at the mode,
shifted lattice points) flags groups whose two estimates differ by more
than three standard errors plus the level error, or with Pareto shape
[\[9\]](#ref-vehtari2024pareto) above 0.7.

## Laplace EM

`method = "laplaceEM"` runs ECM: a trust-region step of the structural
parameters on frozen nodes of level \\K + 1\\, then every block of
\\\Omega\\ from \\S = G^{-1} \sum_i M_i\\ with the posterior second
moments \\M_i = \hat\eta_i \hat\eta_i^\top + 2 H_i^{-1}\\ (modes alone
would shrink \\\Omega\\). It stops when parameter change and
extrapolated remaining distance are below `epsEcm` (default `1e-5`), or
after `maxEcm` (200) iterations. Settings come from
`control$quadrature`.

## SAEM

`method = "saem"` uses posterior draws, with no Gaussian assumption and
equal cost for any \\K\\. Stochastic-approximation EM
[\[10\]](#ref-delyon1999convergence) averages the sufficient statistics
by a Robbins-Monro recursion [\[11\]](#ref-robbins1951stochastic),

\\s_k = s\_{k-1} + \gamma_k \bigl( S(y, \eta^{(k)}) - s\_{k-1} \bigr),
\qquad \eta^{(k)} \sim q(\cdot \mid \psi^{(k-1)}), \qquad \psi^{(k)} =
\hat\psi(s_k), \tag{7}\\

converging to a stationary point of \\\ell\\
[\[10\]](#ref-delyon1999convergence), [\[12\]](#ref-kuhn2004coupling).
\\\Omega\\ follows from \\\sum_i \eta_i \eta_i^\top\\, mu-referenced
fixed effects by generalized least squares, other structural parameters
from an average of quadratic models at the draws, whose Newton step has
the marginal score as fixed point.

| phase | rule |
|----|----|
| draws | few random-walk Metropolis steps and independence proposals from \\\mathcal{N}(\hat\eta_i, 1.5 \cdot 2 H_i^{-1})\\ [\[12\]](#ref-kuhn2004coupling); [`set.seed()`](https://rdrr.io/r/base/Random.html) reproduces |
| burn-in | \\\gamma_k = 1\\, variances shrink at most 5 % per iteration [\[13\]](#ref-lavielle2014mixed); ends when no parameter drifts over `window` |
| averaging | \\\gamma_k = k^{-0.7}\\, mean of iterates [\[14\]](#ref-polyak1992acceleration); ends when the expected Monte Carlo rise of \\-2 \log L\\ is below `mcTol` |

SAEM neither forms \\\ell\\ nor raises it in every iteration;
[`em()`](https://dmodverse.github.io/dMod2/reference/em.md) reports the
quadrature of level \\K + 2\\ with standard errors from the quadrature
Hessian. SAEM suits many random effects and skewed, heavy-tailed or
multimodal posteriors [\[15\]](#ref-kuhn2005maximum).

## Zero variance components

An unsupported variance runs to 0, minus infinity on the log scale. A
lone normal random effect in a block ending below standard deviation 0.1
is therefore refitted at 0. If \\-2 \log L\\ rises by at most
`control$boundary`, the refit is returned, the variance named in
`boundary`.

## References

\[1\]

A. P. Dempster, N. M. Laird, and D. B. Rubin, “Maximum likelihood from
incomplete data via the EM algorithm,” *Journal of the Royal Statistical
Society: Series B*, vol. 39, no. 1, pp. 1–38, 1977, doi:
[10.1111/j.2517-6161.1977.tb01600.x](https://doi.org/10.1111/j.2517-6161.1977.tb01600.x).

\[2\]

X.-L. Meng and D. B. Rubin, “Maximum likelihood estimation via the ECM
algorithm: A general framework,” *Biometrika*, vol. 80, no. 2, pp.
267–278, 1993, doi:
[10.1093/biomet/80.2.267](https://doi.org/10.1093/biomet/80.2.267).

\[3\]

Y. Wang, “Derivation of various NONMEM estimation methods,” *Journal of
Pharmacokinetics and Pharmacodynamics*, vol. 34, no. 5, pp. 575–593,
2007, doi:
[10.1007/s10928-007-9060-6](https://doi.org/10.1007/s10928-007-9060-6).

\[4\]

J. C. Pinheiro and D. M. Bates, “Approximations to the log-likelihood
function in the nonlinear mixed-effects model,” *Journal of
Computational and Graphical Statistics*, vol. 4, no. 1, pp. 12–35, 1995,
doi:
[10.1080/10618600.1995.10474663](https://doi.org/10.1080/10618600.1995.10474663).

\[5\]

A. Genz and B. D. Keister, “Fully symmetric interpolatory rules for
multiple integrals over infinite regions with Gaussian weight,” *Journal
of Computational and Applied Mathematics*, vol. 71, no. 2, pp. 299–309,
1996, doi:
[10.1016/0377-0427(95)00232-4](https://doi.org/10.1016/0377-0427(95)00232-4).

\[6\]

S. A. Smolyak, “Quadrature and interpolation formulas for tensor
products of certain classes of functions,” *Soviet Mathematics Doklady*,
vol. 4, pp. 240–243, 1963.

\[7\]

F. Heiss and V. Winschel, “Likelihood approximation by numerical
integration on sparse grids,” *Journal of Econometrics*, vol. 144, no.
1, pp. 62–80, 2008, doi:
[10.1016/j.jeconom.2007.12.004](https://doi.org/10.1016/j.jeconom.2007.12.004).

\[8\]

W. H. O. Clausen, B. B. Rønn, and I. M. Skovgaard, “Maximum likelihood
estimation in nonlinear mixed effect models: Adaptive Gaussian
quadrature by sparse grid sampling,” in *Population approach group in
europe (PAGE) 18*, 2009, pp. Abstract 1607.

\[9\]

A. Vehtari, D. Simpson, A. Gelman, Y. Yao, and J. Gabry, “Pareto
smoothed importance sampling,” *Journal of Machine Learning Research*,
vol. 25, no. 72, pp. 1–58, 2024.

\[10\]

B. Delyon, M. Lavielle, and E. Moulines, “Convergence of a stochastic
approximation version of the EM algorithm,” *The Annals of Statistics*,
vol. 27, no. 1, pp. 94–128, 1999, doi:
[10.1214/aos/1018031103](https://doi.org/10.1214/aos/1018031103).

\[11\]

H. Robbins and S. Monro, “A stochastic approximation method,” *The
Annals of Mathematical Statistics*, vol. 22, no. 3, pp. 400–407, 1951,
doi: [10.1214/aoms/1177729586](https://doi.org/10.1214/aoms/1177729586).

\[12\]

E. Kuhn and M. Lavielle, “Coupling a stochastic approximation version of
EM with an MCMC procedure,” *ESAIM: Probability and Statistics*, vol. 8,
pp. 115–131, 2004, doi:
[10.1051/ps:2004007](https://doi.org/10.1051/ps:2004007).

\[13\]

M. Lavielle, *Mixed effects models for the population approach: Models,
tasks, methods and tools*. Boca Raton: Chapman; Hall/CRC, 2014. doi:
[10.1201/b17203](https://doi.org/10.1201/b17203).

\[14\]

B. T. Polyak and A. B. Juditsky, “Acceleration of stochastic
approximation by averaging,” *SIAM Journal on Control and Optimization*,
vol. 30, no. 4, pp. 838–855, 1992, doi:
[10.1137/0330046](https://doi.org/10.1137/0330046).

\[15\]

E. Kuhn and M. Lavielle, “Maximum likelihood estimation in nonlinear
mixed effects models,” *Computational Statistics & Data Analysis*, vol.
49, no. 4, pp. 1020–1038, 2005, doi:
[10.1016/j.csda.2004.07.002](https://doi.org/10.1016/j.csda.2004.07.002).
