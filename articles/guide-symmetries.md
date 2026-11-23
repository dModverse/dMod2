# 2. Symmetries

[`symmetryDetection()`](https://dmodverse.github.io/dMod2/reference/symmetryDetection.md),
through the Python package
[symident](https://dmodverse.github.io/symident/), finds continuous
families of parameter values with identical observables in the model of
[Setting up
models](https://dmodverse.github.io/dMod2/articles/guide-model.md).
Code: `symmetries.R` and the end of `model.R`.

## Lie symmetries and detection

Each family is the orbit of a Lie symmetry. symident obtains the
generators as the kernel of the observability-identifiability matrix and
proves scalings exactly ([Methods \> Structural identifiability and Lie
symmetries](https://dmodverse.github.io/dMod2/articles/methods-symmetries.md)).

## Detection settings

- `equilibrate = TRUE`: initial values are the steady state without
  TNFα.
- `reduceCQ = TRUE`: symmetries act on the conserved totals, the
  parameters of
  [`Pimpl()`](https://dmodverse.github.io/dMod2/reference/Pimpl.md).
- `forcings = "TNF"`: TNFα is an input state, set by the event at time
  0.
- `conditions`: a direction is non-identifiable only if it is so in both
  doses.

The imaging readouts bleach, \\y = \log(s x + o) - \tau t\\. In the
control the states are constant, so the bleaching rates are in no
symmetry, and `symmetries.R` drops the bleaching terms.

## Symmetries of the full model

Computed by `symmetries.R`, loaded from `results/symmetries.rds`:

``` r

# at rest without TNFα, which an event adds at time 0 in the stimulated condition; the moieties
# of NFκB and IKK by their totals, the parameters of Pimpl()
doses <- data.frame(dose_TNF = c(0, 1), row.names = c("control", "TNF"))
# without the bleaching terms: in the resting control the states are constant, and the slope of
# an imaging readout is the bleaching rate alone
unbleached <- setNames(eqnvec(sub(" - tau_.*$", "", observables)), names(observables))
symmetries <- symmetryDetection(reactions, unbleached, events = events, forcings = "TNF",
                                conditions = doses, equilibrate = TRUE, reduceCQ = TRUE, reconstruct = TRUE,
                                cores = getOption("dMod2.cores"))
```

``` r

symmetries <- readRDS("results/symmetries.rds")
symmetries
#> Result:  rank 57 / 65  |  8 non-identifiable directions
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

The rank is 57 of 65 coordinates, proven; the kernel is spanned by eight
scalings, each proven over the rationals. The generators are one basis;
up to sums there are four kinds.

- **Unobserved intermediates**, level against driven rate: receptor
  \\X_1\\, IκBα transcript \\X_4\\, A20 transcript \\X_6\\, A20 \\X_7\\
  (with `inh_Ikk` and the A20 imaging scale).
- **Relative measurements**: IKK \\X_2\\; NFκB and IκBα together,
  printed as \\X_3\\ plus \\X_4\\ and \\X_5\\.
- **Nuclear volume** `Vnuc` against the rates reading nuclear
  concentrations, \\X_5\\.
- **A product**: `split_NfkIkb` appears only times `eta_split_pNfkpIkb`,
  \\X_8\\.

## Symmetry reduction

[`symmetryReduction()`](https://dmodverse.github.io/dMod2/reference/symmetryReduction.md)
returns in `$trafo` a transformation without the directions, for
scalings by fixing one coordinate per direction ([Methods \> Gauge and
reduction](https://dmodverse.github.io/dMod2/articles/methods-symmetries.html#gauge-and-reduction)).
Loaded from `results/reduction.rds`:

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

`model.R` fixes one moved parameter per scaling at 1, 0 on log scale,
the gauge:

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

This is the gauge of Oppelt et al. [\[1\]](#ref-oppelt2018model):
concentrations are fractions of the NFκB and IKK totals, receptor and
transcripts are in units of the rates they drive, A20 in units of the
A20 imaging scale.
[Optimization](https://dmodverse.github.io/dMod2/articles/guide-optimization.md)
fits the remaining 63 parameters, [Profile
likelihood](https://dmodverse.github.io/dMod2/articles/guide-profiles.md)
checks which the data determine.

## Discrete symmetries

A discrete symmetry maps a parameter value to isolated other values with
the same output. `symmetryDetection(mode = "global")` counts these
values and reconstructs rational maps between them (requires msolve).
`discrete.R` applies it to the reduced model of [Profile
likelihood](https://dmodverse.github.io/dMod2/articles/guide-profiles.html#reduced-model),
passed as `trafo`; loaded from `results/discrete.rds`:

``` r

# the settings of symmetries.R; the resting state and the fixed values of reduced.R enter as
# initial values and parameters of the detection
doses <- data.frame(dose_TNF = c(0, 1), row.names = c("control", "TNF"))
unbleachedReduced <- setNames(eqnvec(sub(" - tau_.*$", "", observablesReduced)),
                              names(observablesReduced))
known <- c(fixedReduced[names(fixedReduced) != "tau_A20GFP"], resting[names(resting) != "TNF"])
discrete <- symmetryDetection(reactionsReduced, unbleachedReduced,
                              trafo = setNames(as.character(known), names(known)), events = events,
                              forcings = "TNF", conditions = doses, mode = "global",
                              cores = getOption("dMod2.cores"))
```

``` r

saved <- readRDS("results/discrete.rds")
discrete <- saved$discrete
discrete$info$discrete[c("count", "countFinal", "fundamentalDomain")]
#> $count
#> [1] 2
#> 
#> $countFinal
#> [1] TRUE
#> 
#> $fundamentalDomain
#> [1] "degrad_A20 < degrad_RnaA20"
swap <- discrete$discreteSymmetries[[1]]
swap[c("map", "level", "stateMaps")]
#> $map
#> Idx         Inner <- Outer
#>   1    degrad_A20 <- degrad_RnaA20
#>   2 degrad_RnaA20 <- degrad_A20
#> 
#> $level
#> [1] "exact"
#> 
#> $stateMaps
#> $stateMaps[[1]]
#> $stateMaps[[1]]$RnaA20
#> [1] "-A20*degrad_A20 + A20*degrad_RnaA20 + RnaA20"
#> 
#> 
#> $stateMaps[[2]]
#> $stateMaps[[2]]$RnaA20
#> [1] "-A20*degrad_A20 + A20*degrad_RnaA20 + RnaA20"
```

There are two values with the output of a generic one; the map exchanges
the degradation rates of the A20 transcript and of A20, `exact` means
proven by the map of the states. The map moves no coordinate of the
scalings and holds in the original coordinates as well.

The publication does not name the swap and reports both rates with the
same 95 % interval, \\-5.97\\ to \\-3.32\\ on log scale
[\[1\]](#ref-oppelt2018model), Table S7, the union of the two images.
`symmetryReduction(branch = )` gives a chart of the fundamental domain
`degrad_A20 < degrad_RnaA20`:

``` r

# the side degrad_A20 < degrad_RnaA20 of the fixed set, in a chart of positive coordinates
fundamental <- symmetryReduction(discrete, branch = "degrad_A20 < degrad_RnaA20")
```

``` r

fundamental <- saved$fundamental
fundamental$trafo[c("degrad_A20", "degrad_RnaA20")]
#> Idx         Inner <- Outer
#>   1    degrad_A20 <- degrad_A20
#>   2 degrad_RnaA20 <- degrad_A20+delta
```

`delta` is the positive difference of the two rates; in this chart a fit
has one optimum instead of two ([Profile likelihood \> Profiles in the
fundamental
domain](https://dmodverse.github.io/dMod2/articles/guide-profiles.html#profiles-in-the-fundamental-domain)).

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
compounds,” *npj Systems Biology and Applications*, vol. 4, no. 1, p.
23, 2018, doi:
[10.1038/s41540-018-0058-z](https://doi.org/10.1038/s41540-018-0058-z).
