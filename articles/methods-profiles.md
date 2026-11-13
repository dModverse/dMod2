# Profile likelihood

The profile likelihood gives confidence intervals of parameters and
predictions of nonlinear models without a quadratic approximation
[\[1\]](#ref-raue2009structural), [\[2\]](#ref-kreutz2012likelihood).
Guide: [Profile
likelihood](https://dmodverse.github.io/dMod2/articles/guide-profiles.md),
[Prediction
uncertainty](https://dmodverse.github.io/dMod2/articles/guide-prediction-uncertainty.md).

## Profile and threshold

The objective \\f(\theta) = -2 \log L(\theta) - 2 \log \pi(\theta)\\ is
the data term plus the prior \\\pi\\ ([Likelihood and
prior](https://dmodverse.github.io/dMod2/articles/methods-objective.md)).
Along the profile of \\\theta_i\\ the other parameters minimize \\f\\,

\\\hat\theta(c) = \operatorname\*{arg\\min}\_{\theta \\:\\ \theta_i = c}
f(\theta) , \tag{1}\\

and the profile likelihood is the data term along this path,

\\\mathrm{PL}\_i(c) = -2 \log L\bigl(\hat\theta(c)\bigr) . \tag{2}\\

Without a prior, (2) is the minimum of \\-2 \log L\\ with \\\theta_i =
c\\. By the likelihood ratio test, \\\mathrm{PL}\_i(c) - \min_c
\mathrm{PL}\_i(c)\\ below the \\1 - \alpha\\ quantile of \\\chi^2_1\\
(3.84 for 95 %) defines an asymptotic confidence interval, invariant
under reparametrization of \\\theta_i\\.
[`profileThreshold()`](https://dmodverse.github.io/dMod2/reference/profileThreshold.md)
computes the threshold, `method = "F"` the finite-sample one for
estimated standard deviations. A prior wide against the data curvature
keeps (1) finite where the data do not determine \\\theta_i\\ and
changes (2) little where they do.

[`profile()`](https://rdrr.io/r/stats/profile.html) stops on \\f\\ by
default; `stepControl = list(stop = "data")` ends a profile at the
threshold of (2), and `confint(val.column = "data")` reads the interval
off (2).

## Shapes

| shape of (2) | meaning |
|----|----|
| above the threshold on both sides | finite interval |
| below the threshold on one side | practically non-identifiable |
| exactly flat | structurally non-identifiable, a symmetry ([Symmetries](https://dmodverse.github.io/dMod2/articles/guide-symmetries.md)), removed by a gauge |

On the log scale, a one-sided profile towards \\-\infty\\ means the
reaction can be switched off, towards \\+\infty\\ that the process is in
equilibrium on the time scale of the data, replaceable by a large fixed
value or the quasi-steady state. Both limits give a smaller model with
the same fit [\[3\]](#ref-maiwald2016driving). Each profile tests one
limit with all other parameters free; several limits together need a
likelihood ratio test of the reduced model.

## Computing a profile

A profile is a sequence of constrained fits, each started from the last
point.

- `method = "optimize"` steps \\\theta_i\\ and refits the others with
  [`trust()`](https://dmodverse.github.io/dMod2/reference/trust.md), one
  converged fit per point.

- `method = "integrate"` follows the constrained optimum of the
  Lagrangian \\f(\theta) + \lambda\\(\theta_i - c)\\, where \\\nabla f +
  \lambda e_i = 0\\ and \\\theta_i = c\\. Differentiating in \\c\\ gives
  the path equation

  \\\begin{pmatrix} H & e_i \\ e_i^\top & 0 \end{pmatrix}
  \begin{pmatrix} \mathrm{d}\theta/\mathrm{d}c \\
  \mathrm{d}\lambda/\mathrm{d}c \end{pmatrix} = \begin{pmatrix} 0 \\ 1
  \end{pmatrix} \tag{3}\\

  with the Hessian \\H\\ of \\f\\: one evaluation and one linear solve
  per step. With a Gauss-Newton \\H\\ the path drifts;
  `algoControl = list(reoptimize = TRUE)` corrects each step by a short
  fit.

Step control and stops:

| rule | setting |
|----|----|
| step shrinks | mismatch of predicted and actual change above `atol` |
| step doubles | mismatch below `0.3 * atol`, or below `0.3 * rtol` relative |
| side ends | objective (or data term, `stop = "data"`) up by the threshold; `limits` relative to the optimum reached; `limit` steps |

A point below the start value is a better optimum, from which fit and
profiles are repeated. Where the curvature jumps, as at the onset of
[`constraintStability()`](https://dmodverse.github.io/dMod2/reference/constraintStability.md)
([Steady
states](https://dmodverse.github.io/dMod2/articles/methods-steady-states.html#stability)),
the step shrinks to `min` and stays on the profile.

## Prediction profiles

A prediction \\z = h(t, \theta)\\ has the profile likelihood

\\\mathrm{PPL}(z) = -2 \log L\bigl(\hat\theta(z)\bigr) , \quad
\hat\theta(z) = \operatorname\*{arg\\min}\_{\theta \\:\\ h(t, \theta) =
z} f(\theta) , \tag{4}\\

whose threshold gives a confidence interval of the prediction
[\[2\]](#ref-kreutz2012likelihood), the projection of the parameter
confidence region.

A validation data point makes \\z\\ a parameter \\v\\, tied to the
prediction by

\\\left( \frac{h(t, \theta) - v}{\sigma_v} \right)^2 . \tag{5}\\

The profile of \\v\\ is (4) for \\\sigma_v \to 0\\, and within
\\\sigma_v\\ of (4) for \\\sigma_v\\ small against the band.
[`datapointL2()`](https://dmodverse.github.io/dMod2/reference/datapointL2.md)
builds (5) as an `objfn` from the evaluation of
[`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md),
which therefore predicts at \\t\\ (`times = t`). With
`attrName = "data"`, (5) adds to the data term, so `stop = "data"` and
`confint(val.column = "data")` include it. Where the model cannot follow
\\v\\, (5) grows and ends the profile.

## References

\[1\]

A. Raue *et al.*, “Structural and practical identifiability analysis of
partially observed dynamical models by exploiting the profile
likelihood,” *Bioinformatics*, vol. 25, no. 15, pp. 1923–1929, 2009,
doi:
[10.1093/bioinformatics/btp358](https://doi.org/10.1093/bioinformatics/btp358).

\[2\]

C. Kreutz, A. Raue, and J. Timmer, “Likelihood based observability
analysis and confidence intervals for predictions of dynamic models,”
*BMC Systems Biology*, vol. 6, p. 120, 2012, doi:
[10.1186/1752-0509-6-120](https://doi.org/10.1186/1752-0509-6-120).

\[3\]

T. Maiwald *et al.*, “Driving the model to its limit: Profile likelihood
based model reduction,” *PLOS ONE*, vol. 11, no. 9, p. e0162366, 2016,
doi:
[10.1371/journal.pone.0162366](https://doi.org/10.1371/journal.pone.0162366).
