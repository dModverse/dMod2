# Sparsity priors

[`em()`](https://dmodverse.github.io/dMod2/reference/em.md) integrates
random effects whose deviations are exactly zero unless the data need
them. `gr(group, dist = )` in an effect formula sets the distribution:

| distribution | deviation \\\eta\\ of a group | hyperparameters |
|----|----|----|
| [`spikeSlab()`](https://dmodverse.github.io/dMod2/reference/normal.md) | \\0\\ with probability \\1 - \rho\\, else \\\mathcal{N}(0, \tau^2)\\ | `prob_<grouping>`, `tau_<grouping>` |
| [`dirichletProcess()`](https://dmodverse.github.io/dMod2/reference/normal.md) | value of the cluster of the group: \\0\\ for the population cluster, else one \\\mathcal{N}(0, \tau^2)\\ draw per cluster | `concentration_<grouping>`, `tau_<grouping>` |

The spike and slab decides for every group alone whether it deviates.
The Dirichlet process also groups the groups: the groups of one cluster
share one deviation, so a fit tells which groups behave alike, not only
which differ from the population value.

Guide: [Sparse
deviations](https://dmodverse.github.io/dMod2/articles/guide-sparse-deviations.md);
normal random effects: [Marginal likelihood of mixed-effects
models](https://dmodverse.github.io/dMod2/articles/methods-nlme.md).

## Approximate marginal

\\D_i(\eta_i)\\ is the data \\-2 \log L\\ of group \\i\\, summed over
the conditions of the group. The random effects of a group split into
normal coordinates \\N\\ and sparse ones \\S\\.

### Joint mode

[`em()`](https://dmodverse.github.io/dMod2/reference/em.md) first finds
the joint mode \\\hat\eta\\, minimizing

\\\begin{aligned} \Phi(\eta) = {} & \sum_i \Bigl\[ D_i(\eta_i) +
\eta\_{iN}^\top \Omega^{-1} \eta\_{iN} \Bigr\] + \sum\_{\text{spike and
slab, Dirichlet}} \frac{\eta^2}{\tau^2} : \end{aligned} \tag{1}\\

data, normal prior, and the slab alone for spike-and-slab and
Dirichlet-process coordinates. Each group is minimized separately by
[`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md), on
Gauss-Newton plus a BFGS correction from the points that lowered
\\\Phi\\ (Gauss-Newton alone overestimates the curvature at large
residuals). The mode is smooth in the parameters, so the marginal does
not jump where the most probable structure changes. Random effects that
reach no data rest at 0 ([Reported structure](#reported-structure)).

| inner search | rule |
|----|----|
| start | modes of the last accepted point, moved to first order along the step ([Exact gradient](#exact-gradient)) |
| tolerance | `ftol` and `mtol` of the `tolControl` of the outer trust region where given, else \\10^{-10}\\; `gtol` where given, else that of [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md) |
| local minima | start from the lowest of the current mode and the minima at the last four distinct accepted points, since the lowest outweighs a higher one by \\e^{\Delta\Phi / 2}\\ |
| restarts | at the start and every 25 accepted steps also from 0 and from two fixed moves by the slab standard deviation; [`em()`](https://dmodverse.github.io/dMod2/reference/em.md) draws no random numbers |

### Normal coordinates

The data term is expanded to second order at the mode, with gradient
\\\nabla D_i\\ and Gauss-Newton curvature \\H_i\\. \\H\_{NN}\\ is the
normal block plus \\2 \Omega^{-1}\\, the curvature of [the Laplace
approximation](https://dmodverse.github.io/dMod2/articles/methods-nlme.html#laplace-approximation).
The normal coordinates integrate out exactly on the quadratic model,
leaving the volume \\\log \lvert H\_{NN}/2 \rvert\\ and for the sparse
coordinates the Schur complement

\\\tilde H_i = H\_{SS} - H\_{SN} H\_{NN}^{-1} H\_{NS}, \tag{2}\\

so the data term on the sparse coordinates is, up to a constant,
\\\exp(-\tfrac12 \eta\_{iS}^\top P_i\\ \eta\_{iS} + b_i^\top
\eta\_{iS})\\ with

\\P_i = \tilde H_i / 2, \qquad b_i = \tfrac12 \bigl( \tilde H_i\\
\hat\eta\_{iS} - \nabla_S D_i \bigr) \tag{3}\\

Unlike in FOCEI, \\P_i\\ holds no prior; the sparse prior enters the
integral \\I\\ exactly. \\P_i\\ may be singular where the data do not
see a coordinate; the integral needs no inverse.

### Value of the marginal

\\\begin{aligned} -2 \log L \approx {} & \sum_i \Bigl\[
D_i(\hat\eta_i) + \log \lvert H\_{NN}/2 \rvert + \tfrac12
\hat\eta\_{iS}^\top \tilde H_i\\ \hat\eta\_{iS} - \nabla_S D_i^\top
\hat\eta\_{iS} \\ & \qquad + \hat\eta\_{iN}^\top \Omega^{-1}
\hat\eta\_{iN} \Bigr\] + G \log \lvert \Omega \rvert - 2 \log I + D_0 ,
\end{aligned} \tag{4}\\

\\D_0\\ the \\-2 \log L\\ of conditions without random effects and of
other terms such as priors, and

\\I = \int \prod_i \exp\Bigl( -\tfrac12\\ \eta\_{iS}^\top P_i\\
\eta\_{iS} + b_i^\top \eta\_{iS} \Bigr)\\ p(\eta_S)\\ d\eta_S \tag{5}\\

\\p\\ the sparse prior. Without sparse coordinates (4) is the [FOCEI
value](https://dmodverse.github.io/dMod2/articles/methods-nlme.html#focei).
Every structure, a subset of spike coordinates or a partition into
clusters, is a Gaussian integral. Spike-and-slab priors factor over
groups, the Dirichlet process couples the groups of a symbol, and
off-diagonal entries of \\P_i\\ couple the coordinates of a group.

## Factor integrals

\\I\\ is evaluated in exactly integrated factors:

- every Dirichlet-process symbol, all groups jointly;
- the spike-and-slab symbols of a group jointly, within a budget of
  terms, beyond it each alone.

A factor is integrated against a tilt \\\exp(h^\top \eta - \tfrac12
\eta^\top T \eta)\\, \\T = P\_{ff}\\ uncoupled, diagonal under the
expectation consistent (EC) approximation; \\T = 0\\ needs no special
case. By Fisher’s identity,

\\\frac{\partial \log Z}{\partial h} = \mathbb{E}\\ \eta, \qquad
\frac{\partial \log Z}{\partial T} = -\tfrac12\\ \mathbb{E}\\ \eta
\eta^\top, \tag{6}\\

and a hyperparameter derivative is the posterior expectation of the
derivative of the log prior.

### Spike and slab

The prior [\[1\]](#ref-mitchell1988bayesian),
[\[2\]](#ref-george1993variable) is a point mass at 0 of weight \\1 -
\rho\\ and a normal slab of weight \\\rho\\. With \\V_F = (T\_{FF} +
\tau^{-2} \mathbb{1})^{-1}\\ for the subset \\F\\ of slab coordinates,

\\Z = \sum_F \rho^{\|F\|} (1 - \rho)^{K - \|F\|}\\ \frac{\lvert V_F
\rvert^{1/2}}{\tau^{\|F\|}}\\ \exp\bigl( \tfrac12 h_F^\top V_F\\ h_F
\bigr) , \tag{7}\\

\\2^K\\ Gaussian integrals, the factors \\2\pi\\ canceling. \\\Pr(F \mid
y)\\ is the term of \\F\\ over \\Z\\; given \\F\\, \\\eta_F \sim
\mathcal{N}(V_F h_F, V_F)\\ and the rest is 0. The mixture has

\\\begin{aligned} 1 - p\_{0k} &= \sum\_{F \ni k} \Pr(F \mid y), \\
\mathbb{E}\\ \eta &= \sum_F \Pr(F \mid y)\\ \mu_F, \\ \mathbb{E}\\ \eta
\eta^\top &= \sum_F \Pr(F \mid y) \bigl( V_F + \mu_F \mu_F^\top \bigr),
\end{aligned} \tag{8}\\

\\\mu_F\\, \\V_F\\ zero-filled outside \\F\\; the second moment is full
within the factor. Summed over coordinates sharing them,

\\\begin{aligned} \frac{\partial \log Z}{\partial \rho} &= \sum_k \Bigl(
\frac{1 - p\_{0k}}{\rho} - \frac{p\_{0k}}{1 - \rho} \Bigr), \\
\frac{\partial \log Z}{\partial \tau} &= \sum_k \Bigl( -\frac{1 -
p\_{0k}}{\tau} + \frac{\mathbb{E}\\ \eta_k^2}{\tau^3} \Bigr) .
\end{aligned} \tag{9}\\

The sum is a depth-first search, each coordinate in or out:

| element | rule |
|----|----|
| cost per term | \\O(K^2)\\: including a coordinate adds one row to \\W = L^{-1}\\, \\L\\ the Cholesky factor of \\V_F^{-1}\\ |
| second moment | \\V_F = W^\top W\\; each included coordinate adds the outer product of the new row times the weight of all terms below; means and rows enter by a symmetric rank-\\k\\ update |
| pruning | subtree skipped where a bound falls below \\e^{-40}\\ of the largest term, at most \\2^K e^{-40}\\ of \\Z\\; the bound uses the subset of slab and open coordinates and the conditional precisions as lower bounds on later diagonal entries of \\L\\ |
| order | branch with the larger bound first |
| single coordinate | needs only \\t \> -\tau^{-2}\\ |
| limit | \\2^{24}\\ nodes, beyond which the group goes to EC with one factor per coordinate |

### Clusters

A Dirichlet-process symbol integrates by a recursion over subsets of the
groups ([Partition sum](#partition-sum)).

## Expectation consistent approximation

Where no Dirichlet-process symbol joins two groups, the groups are
integrated one by one, a group exactly where the data do not couple the
factors of the group and the joint spike-and-slab sum stays within the
budget. Where a Dirichlet-process symbol joins two groups, all groups
are integrated together, exactly where this holds for every group.
Otherwise \\I\\ is the EC approximation
[\[3\]](#ref-opper2005expectation), of the group, or of all groups where
they are integrated together. With the coupling \\\mathcal{G}(\eta) =
\prod_i \exp(-\tfrac12 \eta_i^\top P_i \eta_i + b_i^\top \eta_i)\\ and
the statistics \\\phi(\eta) = (\eta_k, -\eta_k^2/2)\_k\\:

- \\q \propto p\\ e^{\lambda^\top \phi}\\, every factor exact against a
  diagonal cavity \\\lambda\\;
- \\r \propto \mathcal{G}\\ e^{\gamma^\top \phi}\\, the coupling exact
  under diagonal sites \\\gamma\\, a Gaussian of precision \\P_i +
  \operatorname{diag}(\gamma_t)\\ per group.

With \\A_q(\lambda) = \log \int p\\ e^{\lambda^\top \phi}\\,
\\A_r(\gamma) = \log \int \mathcal{G}\\ e^{\gamma^\top \phi}\\ and the
product of one-dimensional Gaussians \\A_s(\nu) = \log \int e^{\nu^\top
\phi}\\,

\\\log I \\\approx\\ A_q(\lambda) + A_r(\gamma) - A_s(\nu), \qquad \nu =
\lambda + \gamma, \tag{10}\\

where \\q\\, \\r\\ and \\e^{\nu^\top \phi}\\ agree in mean and variance
of every coordinate. The cavity is the self-consistent field of the
other factors. The result is exact without coupling and for Gaussian
factors, and no bound.

A double loop converges for any coupling
[\[4\]](#ref-heskes2005approximate). The outer step sets \\\nu_t = 1 /
\operatorname{Var}\\ \eta_k\\, \\\nu_h = \mathbb{E}\\ \eta_k /
\operatorname{Var}\\ \eta_k\\. The inner loop minimizes the convex

\\\Psi\_\nu(\lambda) = A_q(\lambda) + A_r(\nu - \lambda), \tag{11}\\

by Newton steps on \\\operatorname{Cov}\_q \phi + \operatorname{Cov}\_r
\phi\\ with backtracking, quadratically convergent.
\\\operatorname{Cov}\_q \phi\\ includes third and fourth factor moments,
and for a Dirichlet-process symbol of up to ten groups the moments
between groups, so a step solves all coordinates jointly. The cavity is
the variable, so no cavity precision is a difference \\1 /
\operatorname{Var}\_r \eta_k - \gamma\_{t,k}\\, and site precisions may
be negative. Steps stay in the domains: Gaussian part positive definite
(Cholesky of the Jacobi-scaled precision), \\\lambda_t \> -\tau^{-2}\\
for slab coordinates, \\\tau^{-2} + \sum_g \min(\lambda\_{t,g}, 0) \>
0\\ for a Dirichlet-process symbol. Evaluations start from the last
accepted sites.

The outer fixed-point step converges linearly, slowly for tight
coupling, so the outer loop takes Newton steps on \\\mathbb{E}\_q
\phi(\lambda^\*(\nu)) = \mathbb{E}\_s \phi(\nu)\\, \\s \propto
e^{\nu^\top \phi}\\. With \\\mathrm{d}\lambda^\* =
(\operatorname{Cov}\_q \phi + \operatorname{Cov}\_r \phi)^{-1}
\operatorname{Cov}\_r \phi\\ \mathrm{d}\nu\\ from the inner
stationarity, the Jacobian is

\\\operatorname{Cov}\_q \phi\\ (\operatorname{Cov}\_q \phi +
\operatorname{Cov}\_r \phi)^{-1} \operatorname{Cov}\_r \phi -
\operatorname{Cov}\_s \phi , \tag{12}\\

solved with row and column scaling, since pinned coordinates have tiny
second moments. The cavities move with \\\mathrm{d}\lambda^\*\\,
starting the next inner loop near the solution.

| loop | ends |
|----|----|
| inner | means within \\10^{-10}\\ standard deviations and variances within \\10^{-10}\\ relative; or within \\10^{-4}\\ when Newton no longer lowers \\\Psi\\ |
| outer | marginals of \\r\\ match the Gaussians of \\\nu\\ to \\10^{-10}\\; or within \\10^{-4}\\ when a fixed-point step no longer raises the free energy \\F = A_q + A_r - A_s\\ (rounding of an ill-conditioned Gaussian part) |
| outer Newton step | accepted where it raises \\F\\ by the Armijo condition, the step halved down to a sixteenth; else the fixed-point step, which always raises \\F\\ |
| failure | a failed evaluation; the trust region shortens the step |

`fit$coupling` is `"expectationConsistent"` where a group needed it,
`NULL` otherwise.

## Exact gradient

At the fixed point the value is stationary in sites and cavities, so
derivatives hold them fixed. For exact factors, Fisher’s identity gives

\\\frac{\partial \log I}{\partial b_i} = \mathbb{E}\\ \eta\_{iS}, \qquad
\frac{\partial \log I}{\partial P_i} = -\tfrac12 \bigl( \mathbb{E}\\
\eta\_{iS}\\ \mathbb{E}\\ \eta\_{iS}^\top + \Sigma_i \bigr), \tag{13}\\

\\\Sigma_i\\ the posterior covariance, block diagonal over factors.
Under EC, \\\mathbb{E}\\ is \\\mathbb{E}\_r\\, \\\Sigma_i\\ the
covariance of \\r\\, full within the group, and hyperparameter
derivatives are taken per factor at the cavity.

Through (3), with \\\delta_i = \mathbb{E}\\ \eta\_{iS} -
\hat\eta\_{iS}\\, the derivatives of \\-2 \log L\\ in \\\nabla_S D_i\\,
\\\tilde H_i\\ and the mode are

\\\delta_i, \qquad \tfrac12 \bigl( \delta_i \delta_i^\top + \Sigma_i
\bigr), \qquad -\tilde H_i\\ \delta_i - \nabla_S D_i , \tag{14}\\

and through (2) those in \\H_i\\ and \\\Omega\\. The chain then needs
\\\partial \nabla D_i / \partial \psi\\ from the exact data Hessian,
\\\partial H_i / \partial \psi\\ and \\\partial H_i / \partial \eta\\
from second derivatives as for [the FOCEI volume
term](https://dmodverse.github.io/dMod2/articles/methods-nlme.html#volume-term),
and the mode shift.

The mode moves in all coordinates, normal, spike-and-slab and
Dirichlet-process. The inner stationarity \\\nabla\_\eta \Phi = 0\\
gives

\\\frac{d \hat\eta}{d \psi} = -\bigl( \nabla^2\_\eta \Phi \bigr)^{-1}
\nabla^2\_{\eta \psi} \Phi , \tag{15}\\

group by group. The derivative of the slab gradient in \\\tau\\, \\-4
\eta / \tau^3\\, gives the mode shift in \\\tau\\; \\\rho\\ and
\\\alpha\\ do not move the mode. The shift also starts the inner search
at the next trial point. The gradient is exact for the approximate
marginal.

## Hyperparameters

### E-step

Every evaluation of the marginal is an E-step; the M-steps set the
hyperparameter curvature of the joint trust region. EM in general: [EM
algorithm](https://dmodverse.github.io/dMod2/articles/methods-nlme.html#em-algorithm).
The objective is \\-2Q\\, \\Q\\ the expected complete-data log-density.

With the indicator \\z\_{ik}\\ (1 slab, 0 spike), the complete-data \\-2
\log\\ prior over groups \\i\\ and symbols \\k\\ is, up to terms free of
\\\rho\\ and \\\tau\\,

\\\begin{aligned} -2 \log p(\eta, z \mid \rho, \tau) = \sum\_{i,k}
\Bigl\[ & -2 z\_{ik} \log \rho - 2 (1 - z\_{ik}) \log (1 - \rho) \\ & +
z\_{ik} \Bigl( 2 \log \tau + \frac{\eta\_{ik}^2}{\tau^2} \Bigr) \Bigr\]
. \end{aligned} \tag{16}\\

The posterior expectation is

\\\begin{aligned} -2Q(\rho, \tau) = {} & -2 \bar n \log \rho - 2 (n -
\bar n) \log (1 - \rho) \\ & + 2 \bar n \log \tau + \frac{\mathbb{E}
\sum v^2}{\tau^2} , \end{aligned} \tag{17}\\

\\n\\ the number of deviations of the term, and the statistics

\\\begin{aligned} \bar n &= \sum\_{i,k} \mathbb{E}\\ z\_{ik} =
\sum\_{i,k} (1 - p\_{0,ik}), \\ \mathbb{E} \sum v^2 &= \sum\_{i,k}
\mathbb{E}\\ z\_{ik} \eta\_{ik}^2 = \sum\_{i,k} \mathbb{E}\\
\eta\_{ik}^2 , \end{aligned} \tag{18}\\

since \\\eta\_{ik} = 0\\ wherever \\z\_{ik} = 0\\. The expectations are
under the posterior of the approximate marginal, the quadratic model at
the joint mode ([Normal coordinates](#normal-coordinates)) times the
exact prior, independent over groups:

\\p(\eta\_{iS} \mid y_i) \propto \exp\Bigl( -\tfrac12\\ \eta\_{iS}^\top
P_i\\ \eta\_{iS} + b_i^\top \eta\_{iS} \Bigr) \prod_k p(\eta\_{ik}) ,
\tag{19}\\

The spike-and-slab symbols of a group form one factor ([Spike and
slab](#spike-and-slab)), exact on \\(P_i, b_i)\\ without coupled sparse
coordinates. Under EC each coordinate is a factor on a one-dimensional
Gaussian of precision \\a\\ and center \\m\\, a sum of two terms, with
slab posterior

\\V = \frac{1}{a + \tau^{-2}}, \qquad \mu = V a\\ m . \tag{20}\\

deviation probability

\\1 - p_0 = \operatorname{logit}^{-1} \Bigl( \operatorname{logit} \rho +
\frac{\mu^2}{2 V} + \tfrac12 \log \frac{V}{\tau^2} \Bigr) : \tag{21}\\

prior log odds, exponent gain at the slab mean, and log width ratio. The
moments are

\\\begin{aligned} \mathbb{E}\\ \eta &= (1 - p_0)\\ \mu, \\ \mathbb{E}\\
\eta^2 &= (1 - p_0)\\ (\mu^2 + V), \\ \operatorname{Var} \eta &=
\mathbb{E}\\ \eta^2 - (\mathbb{E}\\ \eta)^2 . \end{aligned} \tag{22}\\

Here \\a = a^c\_{ik}\\, \\m = m^c\_{ik}\\ is the cavity
([EC](#expectation-consistent-approximation)); the moments are those of
\\q\\, equal to those of \\r\\ at the fixed point, and \\\Sigma_i =
(P_i + \operatorname{diag} \gamma\_{t,i})^{-1}\\ with the site
precisions \\\gamma\_{t,i}\\ of the group.

Given the sparse coordinates, on the expansion

\\\begin{aligned} \eta\_{iN} \mid \eta\_{iS} &\sim \mathcal{N}\bigl(
\hat\eta\_{iN} - R_i\\ (\eta\_{iS} - \hat\eta\_{iS}),\\ 2 H\_{NN}^{-1}
\bigr), \\ R_i &= H\_{NN}^{-1} H\_{NS} , \end{aligned} \tag{23}\\

so `fit$etaCov`, the posterior covariance of all random effects of the
group, has the blocks

\\\Sigma\_{NN} = 2 H\_{NN}^{-1} + R_i \Sigma\_{SS} R_i^\top, \qquad
\Sigma\_{NS} = -R_i \Sigma\_{SS} , \tag{24}\\

and `se` of `fit$random` is the square root of the diagonal. \\\Omega\\
has no M-step; it moves with the fixed parameters.

By Fisher’s identity the derivative of \\-2 \log I\\ at the current
point is that of \\-2Q\\,

\\\begin{aligned} \frac{\partial (-2Q)}{\partial \rho} &= -2 \Bigl(
\frac{\bar n}{\rho} - \frac{n - \bar n}{1 - \rho} \Bigr), \\
\frac{\partial (-2Q)}{\partial \tau} &= 2 \Bigl( \frac{\bar n}{\tau} -
\frac{\mathbb{E} \sum v^2}{\tau^3} \Bigr), \end{aligned} \tag{25}\\

the marginal gradient in \\\tau\\ adding the mode shift ([Exact
gradient](#exact-gradient)).

### M-steps

The closed-form M-step of a hyperparameter minimizes \\-2Q\\ at the
current posterior moments; the complete data are \\(\eta, z)\\ for the
spike and slab and partition plus cluster values for the Dirichlet
process. With \\n\\ the number of deviations of a term (groups times
symbols), \\\bar n\\ the expected number of slab deviations or of
clusters off the population value, and \\\mathbb{E} \sum v^2\\ the
expected sum of their squares,

\\\tau^2 = \frac{\mathbb{E} \sum v^2}{\bar n}, \qquad \rho = \frac{\bar
n}{n}, \qquad \bar n = \sum_k \sum\_{i=1}^{G} \frac{\alpha}{\alpha + i}
, \tag{26}\\

the last for the concentration, summed over the symbols \\k\\ of the
term; the sum increases in \\\alpha\\, and Newton’s method on \\\log
\alpha\\ solves it.

Near the optimum an M-step is a gradient step scaled by the
complete-data information \\\mathcal{I}\_c\\, while the data hold
\\\mathcal{I}\_{\text{obs}} = \mathcal{I}\_c -
\mathcal{I}\_{\text{mis}}\\. The error shrinks per iteration by the
fraction of missing information [\[5\]](#ref-dempster1977maximum),

\\\mathcal{I}\_c^{-1} \mathcal{I}\_{\text{mis}} , \tag{27}\\

close to one for a hyperparameter from a few groups. Newton on the
curvature of \\-2Q\\ is worse far from the optimum: on \\u = \log \tau\\
with M-step \\u^\* = \tfrac12 \log(\mathbb{E} \sum v^2 / \bar n)\\ it
steps \\\bigl(1 - e^{-2(u^\* - u)}\bigr) / 2\\, which is \\0.5\\ for
\\u^\*\\ five above and \\-11013\\ for five below.

### Joint trust region

[`em()`](https://dmodverse.github.io/dMod2/reference/em.md) moves fixed
parameters and hyperparameters, on \\\log \tau\\, \\\log \alpha\\,
\\\operatorname{logit} \rho\\, in one trust region with the exact
marginal gradient.

For the fixed parameters the model Hessian is the Schur complement over
all random effects, the sparse ones at the curvature of their posterior:
a factor with posterior covariance \\\Sigma\\ and data curvature
\\\tilde H\_{ff}\\ gets the prior curvature \\2 \Sigma^{-1} - \tilde
H\_{ff}\\, so the block of the factor after eliminating the normal
coordinates is \\2 \Sigma^{-1}\\; under EC the prior curvature is \\2
\operatorname{diag}(\gamma_t)\\, negative site precisions set to 0. The
data Hessian blocks are all Gauss-Newton (forward over reverse), making
the Schur complement positive semidefinite, or all exact (second-order
sensitivities); mixed blocks can be far from definite. Exact rows enter
gradient and mode shift. With normal random effects alone this is the
FOCEI Schur complement.

For every hyperparameter the model Hessian is the secant to the M-step,

\\\kappa_h = \frac{\partial_u (-2Q)}{u - u^\*}, \tag{28}\\

with

\\\partial_u (-2Q) = \begin{cases} 2 \bigl( \bar n - \mathbb{E} \sum v^2
/ \tau^2 \bigr) & \text{for } \tau, \\ -2 \\ (\bar n - n \rho) &
\text{for } \rho, \\ -2 \bigl( \bar n - \sum_k \sum\_{i=1}^{G} \alpha /
(\alpha + i) \bigr) & \text{for } \alpha . \end{cases} \tag{29}\\

\\-2Q\\ is convex in \\u\\, so \\\kappa_h \> 0\\; at \\u = u^\*\\,
\\\kappa_h\\ is the complete-data information

\\\partial_u^2 (-2Q) = \begin{cases} 4\\ \mathbb{E} \sum v^2 / \tau^2 &
\text{for } \tau, \\ 2 n \rho (1 - \rho) & \text{for } \rho, \\ 2 \sum_k
\sum\_{i=1}^{G} \alpha\\ i / (\alpha + i)^2 & \text{for } \alpha ,
\end{cases} \tag{30}\\

and \\\kappa_h \ge 10^{-8}\\. A Newton step on \\\kappa_h\\ alone moves
\\u\\ by

\\\Delta u = -\frac{1}{\kappa_h} \frac{\partial (-2 \log L)}{\partial u}
= u^\* - u - \frac{1}{\kappa_h}\\ \nabla\_{\hat\eta} (-2 \log L)^\top
\frac{d \hat\eta}{d u} , \tag{31}\\

with the mode shift \\d \hat\eta / d u\\ of [Exact
gradient](#exact-gradient): the exact M-step for \\\rho\\ and
\\\alpha\\, plus a term for \\\tau\\ that vanishes where the data term
is quadratic in the random effects.

A BFGS correction learns the missing information and the cross curvature
[\[6\]](#ref-lange1995quasi), [\[7\]](#ref-jamshidian1997acceleration)
from every step between two points that lowered \\-2 \log L\\. Pairs
with \\0 \< s^\top y \le 10^{-8} \lVert s \rVert \lVert y \rVert\\ are
skipped; a pair with \\s^\top y \le 0\\, or an update that would make
the model indefinite, drops the correction. Far from the optimum a step
is an EM step in the hyperparameters and a Gauss-Newton step in the
fixed parameters, at one marginal evaluation; near it the model
approaches the marginal Hessian and convergence is superlinear. Gradient
and bounded Hessian suffice for convergence [\[8\]](#ref-conn2000trust),
and no number of EM iterations is chosen.

The trust region stops at a predicted reduction below `mtol`; a run
ending short of it restarts with a fresh correction while runs lower
\\-2 \log L\\ by more than `mtol`.
[`msem()`](https://dmodverse.github.io/dMod2/reference/msem.md)
multistarts with concentrations and slab standard deviations perturbed
on the log scale, slab probabilities on the logit scale. No model
Hessian is the marginal Hessian, so `fit$hessian` is `NULL` and standard
errors come from the group scores ([Standard errors](#standard-errors)).

## Reported structure

The fit reports the most probable structure: the most probable subset
per spike-and-slab factor, the most probable partition per
Dirichlet-process symbol (the partition recursion with maximum for sum).
Spike deviations are exactly zero. A slab deviation is the slab mean
given that structure: in a joint factor \\\mu\_{F^\*} = V\_{F^\*}
b\_{F^\*}\\; for a coordinate alone or under EC, \\\mu = V a\\ m\\
([E-step](#e-step)) where \\1 - p_0 \> \tfrac12\\, at the cavity under
EC. Free clusters take the posterior mean given the partition.

| column of `fit$random` | content |
|----|----|
| `mode` | this value |
| `se` | posterior standard deviation over spike and slab |
| `cluster` | 0 for the population value, 1, 2, … for deviating groups or clusters in group order (spike and slab, Dirichlet process) |
| `p0` | posterior probability of a zero deviation (spike and slab, Dirichlet process) |
| `data` | whether the random effect reaches data of the group, where some do not |

Before the fit,
[`em()`](https://dmodverse.github.io/dMod2/reference/em.md) traces every
sparse random effect through parameter transformations, steady states,
right-hand sides, events, observables and error model to the measured
names of the group, a search in the symbol graph with the names that are
0 in the condition set to 0. A product with a zero factor vanishes; a
state that starts at 0, whose right-hand side vanishes with the other
zeros and that no event moves, stays 0 by uniqueness of the solution. A
steady-state root depends on an input only through the states whose
equations or conservation rows the input reaches, directly or through
moved states, where the remaining subsystem has full structural rank,
and through every root otherwise. No path means a sensitivity of 0 for
all parameter values.

A random effect without data rests at 0: no mode search, integral factor
1, no part in the hyperparameter updates, reported with the prior
standard deviation, for a Dirichlet-process one with the slab standard
deviation \\\tau\\, and, for the spike and slab, \\p_0 = 1 - \rho\\. A
Dirichlet-process group without data leaves the partition without
changing the others, since the Ewens prior of a subset is that of the
subset. This removes flat directions of mode and integral.

`fit$slab` lists the probabilities, concentrations and slab standard
deviations, each with the symbols, the estimate and whether it was
fixed.

## Selection of deviations

Choosing deviations is a choice among \\2^K\\ structures per group. A
likelihood-ratio test per coordinate depends on which others are free
and, for coupled coordinates, on the test order.

An L1 penalty \\\lambda \sum \lvert \eta \rvert\\
[\[9\]](#ref-tibshirani1996regression), as in
[`trustL1()`](https://dmodverse.github.io/dMod2/reference/trust.md),
sets small deviations to 0, shrinks every selected one by the same
amount, and compensates by selecting more small ones. Of coordinates the
data tie together it keeps one, unstably, since the lasso condition for
the right set fails for such correlations. It gives no probability of a
deviation.

The spike and slab separates decision and size: a deviation is 0 or from
the slab, which shrinks little, and (7) weighs every structure by the
marginal likelihood, with \\\rho\\ and \\\tau\\ estimated instead of
scanned. Coupled coordinates the data cannot separate share the
posterior, with intermediate \\p_0\\. Both statements hold within the
Laplace approximation of the data term.

Groups that deviate alike are a second structure. The Dirichlet process
([Dirichlet-process random effects](#dirichlet-process-random-effects))
weighs every partition of the groups into clusters with one common
deviation by the marginal likelihood, the population value one of the
clusters. The penalized counterpart is the clustered lasso
[\[10\]](#ref-hauber2023uncovering), a penalty on the pairwise
differences of the group parameters, for which Hauber et al. favor the
exponent 0.8 over the L1 norm; `fuse` of
[`trustL1()`](https://dmodverse.github.io/dMod2/reference/trust.md) is
the L1 form: groups with equal values move as one, and without an anchor
no group is the reference. Like the L1 penalty it is scanned over the
penalty strength and gives no probability of a cluster.

Where the data leave fixed-effect directions flat, the fixed-effect
prior decides between equally good structures, and the selection moves
with the width of that prior. Only a model reduced to what the data
determine resolves this.

## Standard errors

With spike-and-slab random effects, \\-2 \log L = -2 \sum_i \log L_i +
D_0\\. Mode, curvature and integral of a group depend on the data of
that group alone, and the mode shift is block diagonal, so the exact
gradient splits into the group gradients. With the scores \\s_i = \nabla
\log L_i\\ and their mean \\\bar s\\, the empirical information is

\\\mathcal{I} = \frac{G}{G - 1} \sum\_{i=1}^{G} (s_i - \bar s)(s_i -
\bar s)^\top + \tfrac12 H_0 , \tag{32}\\

\\H_0\\ the Hessian of \\D_0\\, priors included; `fit$information` holds
it, [`vcov()`](https://rdrr.io/r/stats/vcov.html) inverts it. At the
optimum the scores sum to \\\tfrac12 \nabla D_0\\, which the centering
removes: the first term has rank at most \\G - 1\\, and \\H_0\\ must
cover the remaining directions. Otherwise
[`vcov()`](https://rdrr.io/r/stats/vcov.html) is an error and standard
errors are `NA`. In the guide fit, five compounds give rank at most
four, and the conditions without a compound and the prior give the rest.
Under
[`dirichletProcess()`](https://dmodverse.github.io/dMod2/reference/normal.md)
the marginal is no sum over groups, so there is neither Hessian nor
empirical information; `confint(method = "profile")` gives intervals.

## Dirichlet-process random effects

`gr(group, dist = dirichletProcess())` lets groups share deviations:
groups and the population value, the anchor, form clusters; the anchor
cluster has deviation zero, every other one a common slab deviation
\\\mathcal{N}(0, \tau^2)\\. The reported partition is a grouping of the
groups by their behavior, found by the fit rather than given. \\\alpha\\
and \\\tau\\ are estimated by marginal maximum likelihood, without a
penalty scan. The reported structure is that of the spike and slab with
clusters for deviations. The guide compounds are fitted without it.

### Ewens prior

The partitions \\\mathcal{C}\\ of the \\G\\ groups and the anchor follow
the Ewens distribution [\[11\]](#ref-ewens1972sampling) of the Chinese
restaurant process, anchor seated first: a group joins a cluster with
probability proportional to the size, anchor counted, or opens one with
probability proportional to \\\alpha\\,

\\\Pr(\mathcal{C}) = \frac{\alpha^{\lvert \mathcal{C} \rvert - 1}
\prod\_{T \in \mathcal{C}} (\lvert T \rvert - 1)!} {\prod\_{i=1}^{G}
(\alpha + i)} , \tag{33}\\

\\\lvert \mathcal{C} \rvert\\ the number of clusters, \\\lvert T
\rvert\\ their sizes. The prior depends on the sizes only: it is
exchangeable, independent of the data, and consistent in \\G\\. \\\alpha
\to 0\\ puts every group at the population value, \\\alpha \to \infty\\
in separate clusters; the expected number of free clusters is
\\\sum\_{i=1}^{G} \alpha / (\alpha + i)\\.

### Partition sum

The data term is expanded at the mode under the slab alone. Per symbol,
on the group factors \\\exp(-\tfrac12 a_i \eta^2 + b_i \eta)\\ of the
natural form (5), a cluster \\T\\ with

\\A_T = \sum\_{i \in T} a_i, \qquad B_T = \sum\_{i \in T} b_i \tag{34}\\

has the weight

\\w(T) = (\lvert T \rvert - 1)! \tag{35}\\

with the anchor, at deviation 0, and without it, the cluster value
integrated over the slab,

\\w(T) = \alpha\\ (\lvert T \rvert - 1)!\\ \bigl(1 + \tau^2
A_T\bigr)^{-1/2} \exp\Bigl( \frac{B_T^2}{2 (A_T + \tau^{-2})} \Bigr) .
\tag{36}\\

With \\\mathcal{Z}(X)\\ the sum over the partitions of \\X\\ of the
products of cluster weights,

\\\mathcal{Z}(X) = \sum\_{T \subseteq X,\\ \min X \in T} w(T)\\
\mathcal{Z}(X \setminus T), \qquad I =
\frac{\mathcal{Z}(U)}{\prod\_{i=1}^{G} (\alpha + i)} , \tag{37}\\

\\3^{G+1}\\ terms, exact for up to 14 groups. \\T\\ is a cluster with
probability \\w(T)\\ \mathcal{Z}(U \setminus T) / \mathcal{Z}(U)\\, a
free cluster value then normal with mean \\B_T / (A_T + \tau^{-2})\\ and
variance \\1 / (A_T + \tau^{-2})\\; moments, \\p_0\\, M-step statistics
and derivatives follow. The most probable partition uses maximum for
sum.

### Gradient and M-steps

From the posterior moments,

\\\frac{\partial \log I}{\partial \alpha} = \frac{\bar n}{\alpha} -
\sum\_{i=1}^{G} \frac{1}{\alpha + i}, \qquad \frac{\partial \log
I}{\partial \tau} = -\frac{\bar n}{\tau} + \frac{\mathbb{E} \sum
v^2}{\tau^3} , \tag{38}\\

\\\bar n\\ the expected number of free clusters, \\\mathbb{E} \sum v^2\\
the expected sum of their squares; the chain through \\\nabla D_i\\,
\\H_i\\ and the mode is as above. The mode shifts with \\\tau\\, not
\\\alpha\\. The M-steps are those of [M-steps](#m-steps), summed over
symbols.

### SAEM

The exact partition sum is limited to 14 groups; `method = "saem"` draws
the partition instead, with normal random effects beside the
Dirichlet-process ones and no other sparse ones. Groups start in
separate clusters at the slab-only mode, so the sampler merges clusters.
Per iteration:

- **Normal random effects.** `nMcmc` random-walk Metropolis steps, step
  size adapted in burn-in toward acceptance 0.4.
- **Clusters.** One Gibbs sweep per symbol, Neal’s algorithm 8 with two
  auxiliary clusters [\[12\]](#ref-neal2000markov): each group, in
  random order, joins an existing cluster with probability proportional
  to size (anchor counted) times likelihood at the cluster value, or an
  auxiliary one with a slab draw with probability proportional to
  \\\alpha / 2\\ times likelihood. A group alone offers the current
  cluster value as one auxiliary. One prediction gives the likelihoods
  of all groups at one value.
- **Cluster values.** Three joint random-walk Metropolis steps, accepted
  per cluster on the data of the groups in the cluster and the slab.
- **A-step and M-steps.** Robbins-Monro averages of the free-cluster
  count and squared values give \\\tau\\ and \\\alpha\\; in burn-in they
  change by a factor in \\\[0.95, 1.25\]\\ per iteration and the
  variances of \\\Omega\\ shrink at most 5 %. The rest moves as in
  [SAEM](https://dmodverse.github.io/dMod2/articles/methods-nlme.html#saem).

| phase | rule |
|----|----|
| burn-in | ends when no parameter drifts over `window` iterations, or after `nBurnin`; `converged` reports which |
| averaging | `nEM` iterations, step sizes \\k^{-0.7}\\; estimate the mean after the first `window` |
| structure | most frequent partition after that window, cluster values averaged over the draws with that partition, \\p_0\\ the share of draws at the population value |

### Importance sampling

The \\-2 \log L\\ of an SAEM fit comes from sequential importance
sampling under the expansion [\[13\]](#ref-liu1996nonparametric). A
particle seats the groups in turn; group \\i\\ joins a cluster of every
symbol, existing or new, with probability proportional to the Ewens
prior times the predictive density given the seated groups, the cluster
values having a joint normal posterior. With \\\mathcal{C}\_{i-1}\\ the
partition after \\i - 1\\ groups and \\Z(\mathcal{C})\\ the integral of
the seated factors, the weight is

\\\zeta = \prod\_{i=1}^{G} \sum\_{T} \Pr(T \mid \mathcal{C}\_{i-1})\\
\frac{Z(\mathcal{C}\_{i-1} + T)}{Z(\mathcal{C}\_{i-1})} , \tag{39}\\

summed over the seats \\T\\ of group \\i\\. The mean \\\bar\zeta\\ over
\\M\\ particles is unbiased for \\I\\; by the delta method on \\\log
\bar\zeta\\,

\\\texttt{valueSE} = \frac{2\\ \operatorname{sd}(\zeta)}{\sqrt{M}\\
\bar\zeta} . \tag{40}\\

The symbols of a group are seated together on the full \\P_i\\.
`control$saem$particles` sets \\M\\ (500). An SAEM fit has no gradient,
Hessian or empirical information.

## References

\[1\]

T. J. Mitchell and J. J. Beauchamp, “Bayesian variable selection in
linear regression,” *Journal of the American Statistical Association*,
vol. 83, no. 404, pp. 1023–1032, 1988, doi:
[10.1080/01621459.1988.10478694](https://doi.org/10.1080/01621459.1988.10478694).

\[2\]

E. I. George and R. E. McCulloch, “Variable selection via Gibbs
sampling,” *Journal of the American Statistical Association*, vol. 88,
no. 423, pp. 881–889, 1993, doi:
[10.1080/01621459.1993.10476353](https://doi.org/10.1080/01621459.1993.10476353).

\[3\]

M. Opper and O. Winther, “Expectation consistent approximate inference,”
*Journal of Machine Learning Research*, vol. 6, pp. 2177–2204, 2005.

\[4\]

T. Heskes, M. Opper, W. Wiegerinck, O. Winther, and O. Zoeter,
“Approximate inference techniques with expectation constraints,”
*Journal of Statistical Mechanics: Theory and Experiment*, p. P11015,
2005, doi:
[10.1088/1742-5468/2005/11/P11015](https://doi.org/10.1088/1742-5468/2005/11/P11015).

\[5\]

A. P. Dempster, N. M. Laird, and D. B. Rubin, “Maximum likelihood from
incomplete data via the EM algorithm,” *Journal of the Royal Statistical
Society: Series B*, vol. 39, no. 1, pp. 1–38, 1977, doi:
[10.1111/j.2517-6161.1977.tb01600.x](https://doi.org/10.1111/j.2517-6161.1977.tb01600.x).

\[6\]

K. Lange, “A quasi-Newton acceleration of the EM algorithm,” *Statistica
Sinica*, vol. 5, no. 1, pp. 1–18, 1995.

\[7\]

M. Jamshidian and R. I. Jennrich, “Acceleration of the EM algorithm by
using quasi-Newton methods,” *Journal of the Royal Statistical Society:
Series B*, vol. 59, no. 3, pp. 569–587, 1997, doi:
[10.1111/1467-9868.00083](https://doi.org/10.1111/1467-9868.00083).

\[8\]

A. R. Conn, N. I. M. Gould, and P. L. Toint, *Trust-region methods*.
Philadelphia: SIAM, 2000. doi:
[10.1137/1.9780898719857](https://doi.org/10.1137/1.9780898719857).

\[9\]

R. Tibshirani, “Regression shrinkage and selection via the lasso,”
*Journal of the Royal Statistical Society: Series B*, vol. 58, no. 1,
pp. 267–288, 1996, doi:
[10.1111/j.2517-6161.1996.tb02080.x](https://doi.org/10.1111/j.2517-6161.1996.tb02080.x).

\[10\]

A. L. Hauber, M. Rosenblatt, and J. Timmer, “Uncovering specific
mechanisms across cell types in dynamical models,” *PLOS Computational
Biology*, vol. 19, no. 9, p. e1010867, 2023, doi:
[10.1371/journal.pcbi.1010867](https://doi.org/10.1371/journal.pcbi.1010867).

\[11\]

W. J. Ewens, “The sampling theory of selectively neutral alleles,”
*Theoretical Population Biology*, vol. 3, no. 1, pp. 87–112, 1972, doi:
[10.1016/0040-5809(72)90035-4](https://doi.org/10.1016/0040-5809(72)90035-4).

\[12\]

R. M. Neal, “Markov chain sampling methods for Dirichlet process mixture
models,” *Journal of Computational and Graphical Statistics*, vol. 9,
no. 2, pp. 249–265, 2000, doi:
[10.1080/10618600.2000.10474879](https://doi.org/10.1080/10618600.2000.10474879).

\[13\]

J. S. Liu, “Nonparametric hierarchical Bayes via sequential
imputations,” *The Annals of Statistics*, vol. 24, no. 3, pp. 911–930,
1996, doi:
[10.1214/aos/1032526949](https://doi.org/10.1214/aos/1032526949).
