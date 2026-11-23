# Likelihood and prior

Without a prior, dMod2 estimates by maximum likelihood and minimizes
\\-2 \log L\\. With a prior, it minimizes \\-2 \log\\ of the posterior
up to a constant, and the minimum is the maximum a posteriori estimate
[\[1\]](#ref-kaschek2019dmod). The methods articles call the minimized
quantity the objective function \\\mathcal{J}(\theta)\\, \\-2 \log L\\
or \\-2 \log\\ of the posterior. Guide: [Setting up
models](https://dmodverse.github.io/dMod2/articles/guide-model.html#likelihood-and-prior);
derivatives: [Sensitivities and
Hessians](https://dmodverse.github.io/dMod2/articles/methods-sensitivities.md).

## Data term

[`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md)
computes twice the negative log-likelihood of independent Gaussian
errors, for data \\y_j\\, predictions \\g_j(\theta)\\ and standard
deviations \\\sigma_j(\theta)\\:

\\-2 \log L(\theta) = \sum_j \left( \frac{(y_j -
g_j(\theta))^2}{\sigma_j(\theta)^2} + \log \sigma_j(\theta)^2 \right) +
\text{const}. \tag{1}\\

\\\log \sigma_j^2\\ is constant for given \\\sigma_j\\ and bounds an
estimated \\\sigma_j\\: the first term of (1) falls as \\\sigma_j\\
grows, the second rises.

## Error model

\\\sigma_j\\ comes from the data column `sigma` where set, otherwise
from the error model, an `obsfn` \\\sigma = e(g, \theta)\\ of prediction
and parameters. A constant \\\sigma\\ on a log-scale observable \\y =
\log(s\\ x + o)\\ is a relative error: to first order, \\\sigma\\ in
\\\log z\\ is \\\sigma z\\ in \\z\\.

## Prior

[`constraintL2()`](https://dmodverse.github.io/dMod2/reference/constraintL2.md)
adds a Gaussian prior with means \\\mu_k\\ and standard deviations
\\s_k\\,

\\-2 \log \pi(\theta) = \sum_k \left( \left( \frac{\theta_k -
\mu_k}{s_k} \right)^2 + \log s_k^2 \right) + \text{const}, \tag{2}\\

added to the data term by `+`; \\\log s_k^2\\ enters only where \\s_k\\
is a parameter. It keeps an estimated \\s_k\\ from growing, but (2) then
has no minimum: \\\theta_k = \mu_k\\ with \\s_k \to 0\\ lowers it
without limit, so an estimated \\s_k\\ needs a lower bound. On the log
scale, a prior far wider than the data leaves determined parameters in
place and keeps undetermined ones finite, so the minimum exists without
bounds. The Hessian of (2) is \\2/s_k^2\\ on the diagonal.

## Further terms

- [`constraintStability()`](https://dmodverse.github.io/dMod2/reference/constraintStability.md):
  keeps the steady state of a
  [`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md)
  stable ([Steady
  states](https://dmodverse.github.io/dMod2/articles/methods-steady-states.html#stability)).
- [`datapointL2()`](https://dmodverse.github.io/dMod2/reference/datapointL2.md):
  a data point whose value is a parameter, for validation profiles
  [\[2\]](#ref-kreutz2012likelihood).

Each term is an `objfn`; a sum keeps the parts as attributes, so
profiles and likelihood ratios read the data term alone.

Reference:
[`normL2()`](https://dmodverse.github.io/dMod2/reference/normL2.md),
[`constraintL2()`](https://dmodverse.github.io/dMod2/reference/constraintL2.md),
[`datapointL2()`](https://dmodverse.github.io/dMod2/reference/datapointL2.md),
[`constraintStability()`](https://dmodverse.github.io/dMod2/reference/constraintStability.md),
[`Y()`](https://dmodverse.github.io/dMod2/reference/Y.md)

## References

\[1\]

D. Kaschek, W. Mader, M. Fehling-Kaschek, M. Rosenblatt, and J. Timmer,
“Dynamic modeling, parameter estimation, and uncertainty analysis in R,”
*Journal of Statistical Software*, vol. 88, no. 10, pp. 1–32, 2019, doi:
[10.18637/jss.v088.i10](https://doi.org/10.18637/jss.v088.i10).

\[2\]

C. Kreutz, A. Raue, and J. Timmer, “Likelihood based observability
analysis and confidence intervals for predictions of dynamic models,”
*BMC Systems Biology*, vol. 6, no. 1, p. 120, 2012, doi:
[10.1186/1752-0509-6-120](https://doi.org/10.1186/1752-0509-6-120).
