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
ECM [\[2\]](#ref-meng1993maximum) splits the M-step into conditional
steps, ECME [\[3\]](#ref-liu1994ecme) takes some of them on the
likelihood itself: \\\Omega\\ in closed form, structural parameters by
optimization. [Sparsity
priors](https://dmodverse.github.io/dMod2/articles/methods-sparsity-priors.html#hyperparameters)
minimizes \\-2Q\\ instead. The methods replace the posterior, which has
no closed form, by:

| method | posterior and expectation | M-step | reported \\-2 \log L\\ |
|----|----|----|----|
| `"focei"` | Gaussian at the mode | trust region, Laplace | Laplace |
| `"quadrature"` | adaptive sparse grid at the mode | trust region, frozen nodes | quadrature |
| `"laplaceEM"` | Laplace moments, nodes of degree 3 | ECME | quadrature, level \\K + 2\\ |
| `"saem"` | Metropolis draws, averaged | closed form, averaged Newton | quadrature, level \\K + 2\\ |

## FOCEI

### Laplace approximation

With \\\Phi_i(\eta) = -2 \log p(y_i \mid \eta, \psi) - 2 \log p(\eta
\mid \Omega) - K \log(2\pi)\\, the prior without the constant of the
normal density, mode \\\hat\eta_i\\ and the Gauss-Newton curvature
\\H_i\\ of \\\Phi_i\\ at \\\hat\eta_i\\, Laplace’s approximation of (2)
is \\-2 \log L_i \approx \Phi_i(\hat\eta_i) + \log \lvert H_i/2
\rvert\\, the constant canceling against that of the Gaussian integral.
FOCEI minimizes

\\F(\psi) = \sum_i \Bigl\[ \Phi_i(\hat\eta_i) + \log \lvert H_i/2 \rvert
\Bigr\] \tag{4}\\

plus the full \\-2 \log L\\ of the conditions without a group. \\\log
\lvert H_i/2 \rvert\\ is the volume term. \\H_i\\ includes the
derivatives of \\\sigma\\ in \\\eta\\, the interaction of FOCEI
[\[4\]](#ref-wang2007derivation),
[\[5\]](#ref-pinheiro1995approximations). Conditional standard errors of
the modes are the square roots of the diagonal of \\2 H_i^{-1}\\. The
error has either sign and order \\1 / n_i\\ for \\n_i\\ data points in
the group, with a factor that grows with \\K\\.

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
- **Hessian source.** With several groupings the outer loop is
  [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md), and
  `trustControl` takes the arguments of
  [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md),
  `phases` among them. A handover from Gauss-Newton to BFGS, after 20
  iterations or at a stalled radius ([Methods \> Trust
  region](https://dmodverse.github.io/dMod2/articles/methods-trust.html#hessian-phases)),
  reaches the optima of Gauss-Newton on the theophylline, warfarin and
  mavoglurant fits of `bench/nlme/` with more iterations and time: the
  evaluation builds the Gauss-Newton Hessian for the modes and the
  volume term anyway, so a BFGS step costs as much as a Gauss-Newton
  step.
- **Sparsity.** `sparse = TRUE`: the same outer loop on sparse outer
  Hessians, Schur complement and volume term of a group local to the
  parameters of the group, or, where the \\n_g\\ entries of \\\psi\\ a
  group couples satisfy \\n_g^2 \> k\\ n\_\psi\\ with \\n\_\psi\\ the
  number of entries of \\\psi\\, as a low-rank part of rank \\k\\: \\k =
  K\\ for the Schur complement, \\k = 2K + K^2\\ for the volume term.

### Group kernel

One C++ kernel evaluates all groups at once: the modes, the local model
of each group and the integral on it. At the mode, \\\Phi_i(\hat\eta_i +
\delta) \approx \Phi_i(\hat\eta_i) + g_i^\top \delta + \tfrac12
\delta^\top H_i \delta\\, \\g_i\\ the gradient the inner search leaves.
The block \\H\_{NN}\\ of the normal random effects is factored once;
eliminating it gives the outer model, \\\Phi_i^{\psi\psi} - \Phi_i^{\psi
N} H\_{NN}^{-1} \Phi_i^{N\psi}\\, and the curvature \\H\_{LL} - H\_{LN}
H\_{NN}^{-1} H\_{NL}\\ left on the other coordinates \\L\\. An
integrator maps the local model to the term \\-2 \log L_i -
\Phi_i(\hat\eta_i)\\ and to the derivative \\W_i\\ of the term in
\\H_i\\, which the volume derivative (5) contracts with \\\mathrm{d} H_i
/ \mathrm{d}\psi\\. FOCEI is Laplace on the normal block: the term
\\\log \lvert H\_{NN}/2 \rvert\\ and \\W_i = H\_{NN}^{-1}\\.

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
| `derivMode = "forward-forward"` | `deriv2 = TRUE` sensitivities; exact gradient, curvature without third-order terms |
| `secondOrderCorrection = FALSE` | derivative dropped; biased estimates |

Where the model has `deriv2`, the reported Hessian has the exact Schur
complement and the curvature of the volume term without third-order
sensitivities, as [`summary()`](https://rdrr.io/r/base/summary.html)
states. For \\\sigma\\ depending on \\\eta\\, quadrature gives better
standard errors of \\\Omega\\.

## Adaptive sparse-grid quadrature

Laplace biases \\\Omega\\ for skewed posteriors, as for oral doses
sampled near absorption. The cubature \\\int e^{-\Phi_i/2}\\
\mathrm{d}\eta \approx \sum_b W_b\\ e^{-\Phi_i(\eta_b)/2}\\ removes that
error level by level. It combines the nested Genz-Keister rules
[\[6\]](#ref-genz1996fully), 1, 3, 9, 19 and 35 points exact to degree
1, 5, 15, 29 and 51, by Smolyak’s sparse grid
[\[7\]](#ref-smolyak1963quadrature), [\[8\]](#ref-heiss2008likelihood),

\\A(L, K) = \sum\_{\lvert \mathbf{i} \rvert \le L} \Delta^{i_1} \otimes
\cdots \otimes \Delta^{i_K}, \tag{6}\\

\\\Delta^l\\ the difference of the one-dimensional rules of levels \\l\\
and \\l - 1\\. Level \\L\\ is exact for total degree \\2(L - K) + 1\\
against the Gaussian, with nodes polynomial in \\K\\ (for \\K = 3\\,
levels 3 to 7 have 1, 7, 19, 45, 117); level \\K\\ is the single Laplace
node. The grid sits on the FOCEI Gaussian \\\mathcal{N}(\hat\eta_i, 2
H_i^{-1})\\ and moves with \\\psi\\ [\[9\]](#ref-clausen2009maximum).
Negative weights make \\\log L_i\\ a signed log-sum-exp; a nonpositive
sum, an unresolved posterior, fails the evaluation, naming the group.
The gradient is the quadrature of the exact gradient of \\\log L_i\\ on
the current nodes (Fisher’s identity), the Hessian follows by Louis’
identity from data Hessians per node that are exact where the model has
`deriv2` and Gauss-Newton otherwise; the move of the grid with \\\psi\\
is not differentiated.

`method = "quadrature"` starts from the converged FOCEI fit and
alternates:

- **Rebuild.** Grids at the current \\\psi\\, at the lowest level (never
  falling) whose error is at most `tol` over the number of groups. The
  error of level \\l\\ is estimated by \\2\lvert \log L_i^{(l+1)} - \log
  L_i^{(l)} \rvert\\. A level of more than `maxNodes` nodes per group
  (default 1000) is never evaluated; there the estimate is the change
  from level \\l - 1\\, and a group above the tolerance is capped at
  \\l\\. For \\K = 20\\ levels \\K + 1\\ to \\K + 3\\ have 41, 801 and
  10041 nodes, so the default budget ends at \\K + 2\\.
- **Step.** Trust region over all parameters on frozen nodes. For
  mu-referenced fixed effects, \\\phi_i = D_i \beta + \eta_i\\, nodes
  are frozen in \\\phi\\, need no sensitivities and hold over larger
  steps.
- **Stop.** A Newton step on the rebuilt nodes would gain less than
  `tolDecrement`, or \\-2 \log L\\ rose after a step on unchanged
  levels, so the step lies below what the levels resolve. Either way the
  fit is stationary at the levels reached, capped ones included; it
  returns the iterate of lowest \\-2 \log L\\ on the highest levels, and
  the stage trace records the levels, errors and caps per group.

At the estimate, importance sampling (\\t_5\\ proposal at the mode,
shifted lattice points) flags groups whose two estimates differ by more
than three standard errors plus the level error, or with Pareto shape
[\[10\]](#ref-vehtari2024pareto) above 0.7.

## Laplace EM

`method = "laplaceEM"` runs ECME: a trust-region run of at most
`maxCm1Iter` (30) iterations of the structural parameters on the
likelihood over frozen nodes of level \\K + 1\\, then every block of
\\\Omega\\ from \\S = G^{-1} \sum_i M_i\\ with the posterior second
moments \\M_i = \hat\eta_i \hat\eta_i^\top + 2 H_i^{-1}\\ (modes alone
would shrink \\\Omega\\). It stops when parameter change and
extrapolated remaining distance are below `epsEcm` (default `1e-5`), or
after `maxEcm` (200) iterations. Settings come from
`control$quadrature`.

## SAEM

`method = "saem"` uses posterior draws, with no Gaussian assumption and
equal cost for any \\K\\. Stochastic-approximation EM
[\[11\]](#ref-delyon1999convergence) averages the sufficient statistics
by a Robbins-Monro recursion [\[12\]](#ref-robbins1951stochastic),

\\s_k = s\_{k-1} + \gamma_k \bigl( S(y, \eta^{(k)}) - s\_{k-1} \bigr),
\qquad \eta^{(k)} \sim q(\cdot \mid \psi^{(k-1)}), \qquad \psi^{(k)} =
\hat\psi(s_k), \tag{7}\\

converging to a stationary point of \\\ell\\
[\[11\]](#ref-delyon1999convergence), [\[13\]](#ref-kuhn2004coupling).
\\\Omega\\ follows from \\\sum_i \eta_i \eta_i^\top\\, mu-referenced
fixed effects by generalized least squares, other structural parameters
from an average of quadratic models at the draws, whose Newton step has
the marginal score as fixed point.

| phase | rule |
|----|----|
| draws | few random-walk Metropolis steps and independence proposals from \\\mathcal{N}(\hat\eta_i, 1.5 \cdot 2 H_i^{-1})\\ [\[13\]](#ref-kuhn2004coupling); [`set.seed()`](https://rdrr.io/r/base/Random.html) reproduces |
| burn-in | \\\gamma_k = 1\\, variances shrink at most 5 % per iteration; ends when no parameter drifts over `window` |
| averaging | \\\gamma_k = k^{-0.7}\\, \\k\\ counted from the end of burn-in, mean of iterates [\[14\]](#ref-kushner2003stochastic); ends when the expected Monte Carlo rise of \\-2 \log L\\ is below `mcTol` |

SAEM neither forms \\\ell\\ nor raises it in every iteration;
[`em()`](https://dmodverse.github.io/dMod2/reference/em.md) reports the
quadrature of level \\K + 2\\ with standard errors from the quadrature
Hessian. SAEM suits many random effects and skewed, heavy-tailed or
multimodal posteriors [\[13\]](#ref-kuhn2004coupling).

### Spike-and-slab random effects

A spike-and-slab random effect \\\eta\_{ij} = z\_{ij} v\_{ij}\\ has an
indicator \\z\_{ij}\\ with \\\Pr(z\_{ij} = 1) = p\\ and a slab value
\\v\_{ij} \sim \mathcal{N}(0, \tau^2)\\
[\[15\]](#ref-mitchell1988bayesian), [\[16\]](#ref-george1993variable).
Per iteration, after the normal random effects:

- **Local model.** Per group the quadratic model of the data term in the
  slab values, with the curvature clipped to be positive semidefinite,
  refreshed at the mean of the slab-only Gaussian.

- **Collapsed indicator moves.** Each indicator of a group once per
  sweep, in random order, all groups at once. The other slab values of
  the group move with \\v\_{ij}\\ along the regression of the local
  model, \\v_o = b_o + x\\u_o\\, so the state \\b\\ with \\z\_{ij} = 0\\
  and the line stay fixed. On this line the slab mass \\I = \int
  f(x)\\dx\\, \\f\\ the joint density of data and prior, is a
  Gauss-Hermite sum on the local Gaussian \\\mathcal{N}(m, s^2)\\ with
  the exact data term at `nodes` (3) points, which gives the conditional
  inclusion probability with the slab value integrated out
  [\[17\]](#ref-liu1994collapsed),

  \\\hat\pi\_{ij} = \frac{I}{I + f_0}, \qquad f_0 = f \text{ at }
  z\_{ij} = 0 . \tag{8}\\

  The proposal draws \\z\_{ij}\\ with probability \\0.01 +
  0.98\\\hat\pi\_{ij}\\ and \\x\\ from \\\mathcal{N}(m, 1.5\\s^2)\\,
  independently of the current state; an independence
  Metropolis-Hastings step [\[18\]](#ref-tierney1994markov) on
  \\(z\_{ij}, x)\\ against the model leaves the exact posterior
  invariant, whatever the error of the sum.

- **Slab values.** `nMcmc` joint random-walk Metropolis steps of the
  slab values of every group in the slab, scaled by the conditional
  standard deviations of the local model.

- **A-step and M-steps.** \\p\\ and \\\tau\\ follow from Robbins-Monro
  averages of \\\sum \hat\pi\_{ij}\\ and \\\sum
  \hat\pi\_{ij}\\\mathbb{E}(x^2 \mid z\_{ij} = 1)\\, the
  Rao-Blackwellized statistics [\[19\]](#ref-casella1996rao) of the
  sweep, in place of the drawn indicators. An effect the data do not
  inform has \\\hat\pi\_{ij} = p\\ and second moment \\\tau^2\\ and
  leaves both unchanged; counts of drawn indicators would instead let
  \\p\\ wander to 1 in burn-in, where \\\gamma_k = 1\\. In burn-in the
  odds \\p / (1 - p)\\ and \\\tau\\ change by a factor in \\\[0.95,
  1.25\]\\ per iteration.

\\p_0\\ is one minus the mean of \\\hat\pi\_{ij}\\ over the late sweeps,
an effect selected where it is below 1/2. An `em` fit as `init` starts
the chain at the selected values of that fit, the indicators set where
they are not 0, which checks a FOCEI optimum by SAEM. \\-2 \log L\\ at
the estimate is that of `method = "focei"`; with Dirichlet-process
random effects beside on more than 14 groups, or with
`integral = "importance"`, it comes from sequential importance sampling
over the partitions ([Sparsity
priors](https://dmodverse.github.io/dMod2/articles/methods-sparsity-priors.html#importance-sampling)),
each group summing over the subsets of the spike-and-slab effects that
reach the group data.

## Zero variance components

An unsupported variance runs to 0, minus infinity on the log scale. With
one grouping and several blocks, a normal random effect alone in a block
whose estimated standard deviation ends below 0.1 is therefore refitted
at 0. If \\-2 \log L\\ rises by at most `control$boundary`, the refit is
returned, the variance named in `boundary`.

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

C. Liu and D. B. Rubin, “The ECME algorithm: A simple extension of EM
and ECM with faster monotone convergence,” *Biometrika*, vol. 81, no. 4,
pp. 633–648, 1994, doi:
[10.1093/biomet/81.4.633](https://doi.org/10.1093/biomet/81.4.633).

\[4\]

Y. Wang, “Derivation of various NONMEM estimation methods,” *Journal of
Pharmacokinetics and Pharmacodynamics*, vol. 34, no. 5, pp. 575–593,
2007, doi:
[10.1007/s10928-007-9060-6](https://doi.org/10.1007/s10928-007-9060-6).

\[5\]

J. C. Pinheiro and D. M. Bates, “Approximations to the log-likelihood
function in the nonlinear mixed-effects model,” *Journal of
Computational and Graphical Statistics*, vol. 4, no. 1, pp. 12–35, 1995,
doi:
[10.1080/10618600.1995.10474663](https://doi.org/10.1080/10618600.1995.10474663).

\[6\]

A. Genz and B. D. Keister, “Fully symmetric interpolatory rules for
multiple integrals over infinite regions with Gaussian weight,” *Journal
of Computational and Applied Mathematics*, vol. 71, no. 2, pp. 299–309,
1996, doi:
[10.1016/0377-0427(95)00232-4](https://doi.org/10.1016/0377-0427(95)00232-4).

\[7\]

S. A. Smolyak, “Quadrature and Interpolation Formulas for Tensor
Products of Certain Classes of Functions,” *Doklady Akademii Nauk SSSR*,
vol. 148, no. 5, pp. 1042–1045, 1963, Available:
<https://www.mathnet.ru/eng/dan27586>

\[8\]

F. Heiss and V. Winschel, “Likelihood approximation by numerical
integration on sparse grids,” *Journal of Econometrics*, vol. 144, no.
1, pp. 62–80, 2008, doi:
[10.1016/j.jeconom.2007.12.004](https://doi.org/10.1016/j.jeconom.2007.12.004).

\[9\]

W. H. O. Clausen, B. B. Rønn, and I. M. Skovgaard, “Maximum likelihood
estimation in nonlinear mixed effect models: Adaptive Gaussian
quadrature by sparse grid sampling,” in *Population approach group in
europe (PAGE) 18*, 2009, pp. Abstract 1584. Available:
<https://www.page-meeting.org/?abstract=1584>

\[10\]

A. Vehtari, D. Simpson, A. Gelman, Y. Yao, and J. Gabry, “Pareto
smoothed importance sampling,” *Journal of Machine Learning Research*,
vol. 25, no. 72, pp. 1–58, 2024, Available:
<https://jmlr.org/papers/v25/19-556.html>

\[11\]

B. Delyon, M. Lavielle, and E. Moulines, “Convergence of a stochastic
approximation version of the EM algorithm,” *The Annals of Statistics*,
vol. 27, no. 1, pp. 94–128, 1999, doi:
[10.1214/aos/1018031103](https://doi.org/10.1214/aos/1018031103).

\[12\]

H. Robbins and S. Monro, “A stochastic approximation method,” *The
Annals of Mathematical Statistics*, vol. 22, no. 3, pp. 400–407, 1951,
doi: [10.1214/aoms/1177729586](https://doi.org/10.1214/aoms/1177729586).

\[13\]

E. Kuhn and M. Lavielle, “Coupling a stochastic approximation version of
EM with an MCMC procedure,” *ESAIM: Probability and Statistics*, vol. 8,
pp. 115–131, 2004, doi:
[10.1051/ps:2004007](https://doi.org/10.1051/ps:2004007).

\[14\]

H. J. Kushner and G. G. Yin, *Stochastic approximation and recursive
algorithms and applications*, 2nd ed., vol. 35. in Stochastic modelling
and applied probability, vol. 35. New York: Springer, 2003. doi:
[10.1007/b97441](https://doi.org/10.1007/b97441).

\[15\]

T. J. Mitchell and J. J. Beauchamp, “Bayesian variable selection in
linear regression,” *Journal of the American Statistical Association*,
vol. 83, no. 404, pp. 1023–1032, 1988, doi:
[10.1080/01621459.1988.10478694](https://doi.org/10.1080/01621459.1988.10478694).

\[16\]

E. I. George and R. E. McCulloch, “Variable selection via Gibbs
sampling,” *Journal of the American Statistical Association*, vol. 88,
no. 423, pp. 881–889, 1993, doi:
[10.1080/01621459.1993.10476353](https://doi.org/10.1080/01621459.1993.10476353).

\[17\]

J. S. Liu, “The collapsed Gibbs sampler in Bayesian computations with
applications to a gene regulation problem,” *Journal of the American
Statistical Association*, vol. 89, no. 427, pp. 958–966, 1994, doi:
[10.1080/01621459.1994.10476829](https://doi.org/10.1080/01621459.1994.10476829).

\[18\]

L. Tierney, “Markov chains for exploring posterior distributions,” *The
Annals of Statistics*, vol. 22, no. 4, pp. 1701–1728, 1994, doi:
[10.1214/aos/1176325750](https://doi.org/10.1214/aos/1176325750).

\[19\]

G. Casella and C. P. Robert, “Rao-Blackwellisation of sampling schemes,”
*Biometrika*, vol. 83, no. 1, pp. 81–94, 1996, doi:
[10.1093/biomet/83.1.81](https://doi.org/10.1093/biomet/83.1.81).
