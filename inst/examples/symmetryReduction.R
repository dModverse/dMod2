\dontrun{
  library(dMod2)

  ## gene cascade observed through the protein: one scaling direction
  f <- eqnvec(m = "ktx - dm*m", p = "ktl*m - dp*p")
  g <- eqnvec(y = "p")
  res <- symmetryDetection(f, g, method = "observability", reconstruct = TRUE)
  red <- symmetryReduction(res)
  red                             # terse: trafo plus what the outer names carry
  summary(red)                    # invariants, stages, certificates, families
  red$trafo                       # ready for P() or symmetryDetection(trafo = )
  symmetryDetection(f, g, method = "observability", trafo = red$trafo)$identifiable

  ## pre-gauge with your own fixing; redundant fixings are reported
  symmetryReduction(res, fixed = "ktl")

  ## curved direction: the polynomial stage finds k1 + u*k2
  f2 <- eqnvec(A = "k1 + u*k2 - kdeg*A")
  res2 <- symmetryDetection(f2, eqnvec(y = "A"), method = "observability",
                            fixed = "u", reconstruct = TRUE)
  red2 <- symmetryReduction(res2)
  red2
  red2$trafo

  ## capped search: negative certificates instead of silence
  symmetryReduction(res2, dPoly = 0L, dDarboux = 0L, dExp = 0L)

  ## zero limits: k1 + u*k2 is the invariant, so either production term can be
  ## switched off from anywhere -- reported without a point
  symmetryReduction(res2, reportZeroCompatibility = TRUE)$zeroCompatibility

  ## P <-> pP through s*pP: neither zero holds everywhere, so nothing is claimed
  ## until a point decides it
  pp <- eqnlist() |>
    addReaction("P",  "pP", "k_p*P",  "phosphorylation") |>
    addReaction("pP", "P",  "k_d*pP", "dephosphorylation")
  resPP <- symmetryDetection(pp, eqnvec(y = "s*pP"), method = "observability",
                             reconstruct = TRUE)
  redPP <- symmetryReduction(resPP, fixed = "s", reportZeroCompatibility = TRUE)
  redPP                           # k_d and P each "where ...", k_p "nowhere"
  ## the conditions are R over the model's own names -- one eval decides them, and
  ## above the steady state the back reaction goes, below it the initial value
  cond <- setNames(redPP$zeroCompatibility$condition, redPP$zeroCompatibility$coordinates)
  sapply(cond[c("P", "k_d")], function(cc)
    eval(parse(text = cc), list(P = 1, pP = 0.2, k_p = 0.3, k_d = 0.05)))
  sapply(cond[c("P", "k_d")], function(cc)
    eval(parse(text = cc), list(P = 0.2, pP = 1, k_p = 0.05, k_d = 0.3)))

  ## scaling entangled with a general direction: one joint curved block
  f3 <- eqnvec(A = "-k1*A + k2*B", B = "k1*A - k2*B")
  g3 <- eqnvec(y = "alpha*A")
  res3 <- symmetryDetection(f3, g3, method = "observability", reconstruct = TRUE)
  red3 <- symmetryReduction(res3)
  symmetryDetection(f3, g3, method = "observability",
                    trafo = red3$trafo)$identifiable

  ## rational invariant via the separable quadrature stage
  f4 <- eqnvec(x = "-(b - a)/(a*b)*x")
  res4 <- symmetryDetection(f4, eqnvec(y = "x"), method = "observability",
                            reconstruct = TRUE)
  red4 <- symmetryReduction(res4)
  red4$blocks[[1]]$stage          # "separable"
  red4$blocks[[1]]$invariants     # "(a - b)/(a*b)"

  ## the factor stages reach it too: Darboux with the quadrature switched off
  symmetryReduction(res4, dDarboux = 1L, separable = FALSE)$blocks[[1]]$invariants

  ## both capped away: the exp stage takes over (log/exp entries suit P(),
  ## not symmetryDetection(trafo = ))
  symmetryReduction(res4, dDarboux = 0L, separable = FALSE)$blocks[[1]]$invariants

  ## EGF/EGFR -> MEK/ERK cascade: the curved block's invariants are carried on
  ## fresh q_<k> parameters with certified offsets
  reactions <- eqnlist() |>
    addReaction("EGF + EGFR", "EGF_EGFR", "k_bind * EGF * EGFR")   |>
    addReaction("EGF_EGFR", "EGF + EGFR", "k_unbind * EGF_EGFR")   |>
    addReaction("MEK", "pMEK", "k_phos_MEK * EGF_EGFR * MEK")      |>
    addReaction("pMEK", "MEK", "k_dephos_MEK * pMEK")              |>
    addReaction("ERK", "pERK", "k_phos_ERK * pMEK * ERK")          |>
    addReaction("pERK", "ERK", "k_dephos_ERK * pERK")
  
  reactions <- customTotals(reactions, list(
    totalEGF = "EGF + EGF_EGFR", totalEGFR = "EGFR + EGF_EGFR",
    totalMEK = "MEK + pMEK",
    totalERK = "ERK + pERK"))
  
  observables <- eqnvec(pMEK_obs = "pMEK",
                        pERK_obs = "pERK")
  
  mysteadies <- steadyStates(reactions)
  
  events <- eventlist() |> 
    addEvent(var = "EGF", time = 0, value = "1", method = "add")
  
  egf <- symmetryDetection(reactions, observables, method = "observability", events = events,
                           trafo = as.eqnvec(mysteadies), reduceCQ = FALSE, reconstruct = TRUE)
  
  redEgf <- symmetryReduction(egf, verbose = TRUE)
  summary(redEgf)                 # the certificates behind the positive chart
  redEgf$trafo                    # positive chart, ready for P()

}

