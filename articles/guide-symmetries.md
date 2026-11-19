# 2. Symmetries

If a continuous family of parameter values gives the same prediction for
every observable, no data of these observables can tell the members
apart.
[`symmetryDetection()`](https://dmodverse.github.io/dMod2/reference/symmetryDetection.md),
through the Python package
[symident](https://dmodverse.github.io/symident/), finds such families
in the model of [Setting up
models](https://dmodverse.github.io/dMod2/articles/guide-model.md); they
are removed before the fit. Code: `symmetries.R` and the end of
`model.R`.

## Lie symmetries and detection

Such a family is the orbit of a Lie symmetry with a generator, a vector
field on the parameters. Most symmetries of reaction networks are
scalings: the moved parameters are multiplied by powers of one factor,
and on log scale an orbit is a straight line. symident obtains the
generators as the kernel of the observability-identifiability matrix,
built from the Lie derivatives of the observables, and proves scalings
exactly; see [Methods \> Structural identifiability and Lie
symmetries](https://dmodverse.github.io/dMod2/articles/methods-symmetries.md).

## Detection settings

- `equilibrate = TRUE`: the cells start at the steady state without
  TNFα, so initial values are functions of the parameters.
- `reduceCQ = TRUE`: one species per conserved moiety is eliminated,
  symmetries act on the totals, the parameters of
  [`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md).
- `forcings = "TNF"`: TNFα is an input state, 0 at rest and set by the
  event at time 0.
- `conditions`: the two doses; a direction is non-identifiable only if
  it is so in both.

The imaging readouts bleach, \\y = \log(s x + o) - \tau t\\; \\e^{-\tau
t}\\ is no function of the states of a model at rest. In the control the
states are constant and the slope of an imaging readout is \\-\tau\\
alone, so the bleaching rates are in no symmetry, and `symmetries.R`
drops the bleaching terms.

## Symmetries of the full model

``` r

# at rest without TNFα, which an event adds at time 0 in the stimulated condition; the moieties
# of NFκB and IKK by their totals, the parameters of Pimpl()
doses <- data.frame(dose_TNF = c(0, 1), row.names = c("control", "TNF"))
# without the bleaching terms: in the resting control the states are constant, and the slope of
# an imaging readout is the bleaching rate alone
unbleached <- eqnvec(sub(" - tau_.*$", "", observables), names = names(observables))
symmetries <- symmetryDetection(reactions, unbleached, events = events, forcings = "TNF",
                                conditions = doses, equilibrate = TRUE, reduceCQ = TRUE, reconstruct = TRUE,
                                cores = getOption("dMod2.cores"))
```

``` r

symmetries <- readRDS("results/symmetries.rds")
symmetries
#> Result:  rank 71 / 79  |  8 non-identifiable directions
#> 
#> Generators  X = Σᵢ ξ(i) ∂ᵢ
#> 
#> Scalings:
#>   X₁  (exact)
#>       ξ(act_Ikk_by_TNF) :  + act_Ikk_by_TNF
#>       ξ(uptake)         :  - uptake
#> 
#>   X₂  (exact)
#>       ξ(act_Ikb_by_Ikk) :  + act_Ikb_by_Ikk
#>       ξ(act_Nfk_by_Ikk) :  + act_Nfk_by_Ikk
#>       ξ(s_pIKK)         :  + s_pIKK
#>       ξ(totalIkk)       :  - totalIkk
#> 
#>   X₃  (exact)
#>       ξ(Vnuc)                 :  + Vnuc
#>       ξ(eta_form_complex_nuc) :  + eta_form_complex_nuc
#>       ξ(form_complex)         :  - form_complex
#>       ξ(prod_Ikb)             :  + prod_Ikb
#>       ξ(s_IkBa)               :  - s_IkBa
#>       ξ(s_cytNFkB)            :  - s_cytNFkB
#>       ξ(s_ipIkBa_ibIkBa)      :  - s_ipIkBa_ibIkBa
#>       ξ(s_ipIkBa_ibNFkB)      :  - s_ipIkBa_ibNFkB
#>       ξ(s_ipIkBa_ibpNFkB)     :  - s_ipIkBa_ibpNFkB
#>       ξ(s_ipNFkB_ibNFkB)      :  - s_ipNFkB_ibNFkB
#>       ξ(s_ipNFkB_ibpIkBa)     :  - s_ipNFkB_ibpIkBa
#>       ξ(s_ipNFkB_ibpNFkB)     :  - s_ipNFkB_ibpNFkB
#>       ξ(s_pIkBa)              :  - s_pIkBa
#>       ξ(s_pNFkB)              :  - s_pNFkB
#>       ξ(totalNfk)             :  + totalNfk
#> 
#>   X₄  (exact)
#>       ξ(prod_Ikb)          :  + prod_Ikb
#>       ξ(prod_mIkb_by_nNfk) :  - prod_mIkb_by_nNfk
#> 
#>   X₅  (exact)
#>       ξ(Vnuc)                 :  + Vnuc
#>       ξ(build_RnaA20)         :  + build_RnaA20
#>       ξ(eta_form_complex_nuc) :  + eta_form_complex_nuc
#>       ξ(prod_mIkb_by_nNfk)    :  + prod_mIkb_by_nNfk
#>       ξ(s_nucNFkB)            :  + s_nucNFkB
#> 
#>   X₆  (exact)
#>       ξ(build_A20)    :  + build_A20
#>       ξ(build_RnaA20) :  - build_RnaA20
#> 
#>   X₇  (exact)
#>       ξ(build_A20) :  + build_A20
#>       ξ(inh_Ikk)   :  + inh_Ikk
#>       ξ(s_A20GFP)  :  - s_A20GFP
#> 
#>   X₈  (exact)
#>       ξ(eta_split_pNfkpIkb) :  + eta_split_pNfkpIkb
#>       ξ(split_NfkIkb)       :  - split_NfkIkb
#> 
#> exact = proven over Q.
#> 
#> 
#> Discrete symmetries: 0
```

The rank of the observability-identifiability matrix is 71 of 79
coordinates, proven; the kernel is spanned by eight scalings, each
proven over the rationals. The eight printed generators are one basis,
sums of them are symmetries too. Up to such sums there are four kinds.

- **Unobserved intermediates.** Active receptor `TNFR`, IκBα transcript
  `mIkb`, A20 transcript and A20 protein are unmeasured, A20 only up to
  a scale. Each acts through the product of level and driven rate
  constant: level and production rate times \\\lambda\\, driven rate
  divided by \\\lambda\\. \\X_1\\: receptor (`uptake`,
  `act_Ikk_by_TNF`); \\X_4\\: IκBα transcript (`prod_mIkb_by_nNfk`,
  `prod_Ikb`); \\X_6\\: A20 transcript (`build_RnaA20`, `build_A20`).
  \\X_7\\ moves A20 with the A20 production, the IKK inhibition constant
  `inh_Ikk` (entering as `A20 / inh_Ikk`) and the scale of the A20
  imaging readout.
- **Relative measurements.** Blots and imaging give no absolute amounts.
  \\X_2\\ divides the IKK total and every IKK form by \\\lambda\\ and
  multiplies the IKK-catalyzed phosphorylation rates and the pIKK blot
  scale by \\\lambda\\; the IKK cycle is linear in IKK. NFκB and IκBα
  scale together: every form times \\\lambda\\; complex formation
  (second order), the readout scales and the A20 transcription driven by
  nuclear NFκB divided by \\\lambda\\. The printed \\X_3\\ is this
  direction plus \\X_4\\ and \\X_5\\, hence it also moves nuclear volume
  and IκBα production.
- **Nuclear volume.** In \\X_5\\ `Vnuc` grows at fixed amounts, nuclear
  concentrations fall, and the rates reading them (transcription of IκBα
  and A20, nuclear complex formation) and the scale of nuclear NFκB grow
  to compensate.
- **A product.** `split_NfkIkb` appears only multiplied by
  `eta_split_pNfkpIkb` (\\X_8\\).

## Symmetry reduction

[`symmetryReduction()`](https://dmodverse.github.io/dMod2/reference/symmetryReduction.md)
returns in `$trafo` a parameter transformation that removes the
directions, for scalings by fixing one coordinate per direction
([Methods \> Gauge and
reduction](https://dmodverse.github.io/dMod2/articles/methods-symmetries.html#gauge-and-reduction)).

``` r

# a section of every orbit of the scalings
reduction <- symmetryReduction(symmetries, timeout = 120)
```

``` r

reduction <- readRDS("results/reduction.rds")
reduction
#> Reduced 8 of 8 directions.
#> 
#> Trafo (non-identity entries)
#>   act_Ikb_by_Ikk     = 1
#>   act_Ikk_by_TNF     = 1
#>   build_RnaA20       = 1
#>   eta_split_pNfkpIkb = 1
#>   form_complex       = 1
#>   inh_Ikk            = 1
#>   prod_Ikb           = 1
#>   s_nucNFkB          = 1
#> 
#> Sections:
#>   {X₁}                  act_Ikk_by_TNF = 1
#>   {X₂}                  act_Ikb_by_Ikk = 1
#>   {X₃, X₄, X₅, X₆, X₇}  form_complex = 1, s_nucNFkB = 1, inh_Ikk = 1, prod_Ikb = 1, build_RnaA20 = 1
#>   {X₈}                  eta_split_pNfkpIkb = 1
#>   summary(x, detailed = TRUE) explains each section.
#> 
#> Discrete symmetries: 0
```

## Gauge fixing

The objective is flat along every orbit, so the fit fixes one moved
coordinate per symmetry, the gauge. This chooses units and leaves the
fit unchanged. `model.R` fixes one parameter per scaling at 1, 0 on log
scale:

``` r

# one coordinate of each scaling (symmetries.R) fixed at 1, at 0 on the log scale
gauge <- c(Vnuc = 0, totalNfk = 0, totalIkk = 0, uptake = 0, prod_Ikb = 0, build_RnaA20 = 0,
           split_NfkIkb = 0, s_A20GFP = 0)
```

The eight coordinates remove all eight scalings if their weight matrix
is regular ([Methods \> Gauge and
reduction](https://dmodverse.github.io/dMod2/articles/methods-symmetries.html#gauge-and-reduction)):

``` r

W <- t(vapply(symmetries$symmetries, function(X)
  vapply(names(gauge), function(cn) as.numeric(X$weights[[cn]] %||% 0), 0), numeric(length(gauge))))
dimnames(W) <- list(paste0("X", seq_len(nrow(W))), names(gauge))
W
#>    Vnuc totalNfk totalIkk uptake prod_Ikb build_RnaA20 split_NfkIkb s_A20GFP
#> X1    0        0        0     -1        0            0            0        0
#> X2    0        0       -1      0        0            0            0        0
#> X3    1        1        0      0        1            0            0        0
#> X4    0        0        0      0        1            0            0        0
#> X5    1        0        0      0        0            1            0        0
#> X6    0        0        0      0        0           -1            0        0
#> X7    0        0        0      0        0            0            0       -1
#> X8    0        0        0      0        0            0           -1        0
det(W)
#> [1] -1
```

This is the gauge of Oppelt et al. [\[1\]](#ref-oppelt2018model): NFκB
and IKK totals and nuclear volume set the concentration units, so a
concentration is a fraction of the moiety total; receptor and A20
transcript are in units of their production rates, the IκBα transcript
in units of the translation it drives, A20 in units of the A20 imaging
scale. The reduction fixes other coordinates, a different choice of
units.
[Optimization](https://dmodverse.github.io/dMod2/articles/guide-optimization.md)
fits the remaining 63 parameters, [Profile
likelihood](https://dmodverse.github.io/dMod2/articles/guide-profiles.md)
checks which the data determine.

## Discrete symmetries

A locally identifiable model can still give the output of a parameter
value at other, isolated values: a discrete symmetry.
`symmetryDetection(discrete = TRUE)`, or
[`discreteSymmetries()`](https://dmodverse.github.io/dMod2/reference/discreteSymmetries.md)
on a result, counts these values in the domain and reconstructs the maps
between them where they are rational functions of the coordinates; the
search needs msolve and a model with one condition and no events. In the
two-compartment model of Godfrey [\[2\]](#ref-godfrey1983compartmental),
with the amounts \\x_1\\ and \\x_2\\, an input \\u\\ into the first
compartment and the second observed,

``` r

f <- eqnvec(x1 = "-(k21 + k01)*x1 + k12*x2 + u", x2 = "k21*x1 - k12*x2")
compartments <- symmetryDetection(f, eqnvec(y = "x2"), knownInputs = "u", discrete = TRUE)
```

``` r

compartments
#> Result:  structurally locally identifiable (rank 5 / 5)
#> 
#> Discrete symmetries: 1, 2 solutions over the algebraic closure
#>   k01 -> k12, x1 -> (k01*x2 - k12*x2 + k21*x1)/k21, k12 -> k01 [exact, in the domain where (k01*x2 - k12*x2 + k21*x1)/k21 > 0]
#>   fundamental domain: k01 < k12
```

the map exchanges the elimination rate \\k\_{01}\\ and the return rate
\\k\_{12}\\ and moves the initial value of \\x_1\\, so the model is
locally but not globally identifiable. `exact` means proven by a map of
the states. The fundamental domain \\k\_{01} \< k\_{12}\\ keeps one of
the two points; a fit restricted to it has a unique optimum. A model
with directions is reduced first, and the maps are stated in the
coordinates of the reduced model.

Reference:
[`symmetryDetection()`](https://dmodverse.github.io/dMod2/reference/symmetryDetection.md),
[`symmetryReduction()`](https://dmodverse.github.io/dMod2/reference/symmetryReduction.md),
[`discreteSymmetries()`](https://dmodverse.github.io/dMod2/reference/discreteSymmetries.md)

[← 1. Setting up
models](https://dmodverse.github.io/dMod2/articles/guide-model.md) \|
[3. Sensitivities
→](https://dmodverse.github.io/dMod2/articles/guide-sensitivities.md)

## References

\[1\]

A. Oppelt *et al.*, “Model-based identification of TNFα-induced
IKKβ-mediated and IκBα-mediated regulation of NFκB signal transduction
as a tool to quantify the impact of drug-induced liver injury
compounds,” *npj Systems Biology and Applications*, vol. 4, p. 23, 2018,
doi:
[10.1038/s41540-018-0058-z](https://doi.org/10.1038/s41540-018-0058-z).

\[2\]

K. Godfrey, *Compartmental models and their application*. London:
Academic Press, 1983.
