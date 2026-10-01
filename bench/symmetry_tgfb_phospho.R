# -------------------------------------------------------------------------#
# TGF-b signalling with non-canonical pathways: symmetries at scale
# -------------------------------------------------------------------------#
#
# [PURPOSE]
# A large signalling model for symmetryDetection() under equilibrate: the
# canonical Smad core of M017 (TGFbModelling, H1975) plus the MAPK (ERK, JNK,
# p38), PI3K/AKT/mTOR, Rho-GTPase and NF-kB branches, Smad linker crosstalk and
# a gene layer with the pathway feedbacks. Proteins and phosphosites are read by
# targeted proteomics and phosphoproteomics, in absolute units; mRNA carries one
# scale per gene. Kinase inhibitors, knockdowns and overexpressions make the
# conditions. Nothing is fitted: the script runs the detection only.
#
# [AUTHOR]
# Simon Beyer
#
# [Date]
# Wed 30 Sep 2026
#
# [Info]
# Pathways after Derynck & Zhang (2003) and Zhang (2017): ShcA-Grb2/SOS-Ras-Raf-
# MEK-ERK; TRAF6-TAK1 with MKK4-JNK, MKK3/6-p38 and IKK-NF-kB; p85 recruitment to
# TRAF6, AKT-mTORC1-S6K-S6, AKT sequestering Smad3; TbRII-Par6-Smurf1 degrading
# RhoA, RhoA-ROCK-MLC; ERK, JNK and p38 phosphorylating the Smad3 linker.
# Signalling proteins are conserved moieties; genes and receptors turn over.
# Needs dMod2 from devel-symmetry.
# -------------------------------------------------------------------------#

library(dMod2)

.modelname <- "symmetryTGFbPhospho"
.outdir    <- file.path(tempdir(), .modelname)
.cores     <- detectFreeCores()
if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)

# TGFB1 from the gene layer is secreted and binds the receptors
.autocrine <- TRUE
# transcription by its drivers: "linear", "mm" (Michaelis-Menten) or "hill" (exponent 2)
.transcription <- "hill"
# saturable phosphorylation and dephosphorylation in the non-canonical cycles
.mmCycles <- TRUE


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Canonical core (M017)
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
addRC <- function(eq, from, to, rate, ...) addReaction(eq, from, to, rate, compartment = "Cell", ...)
addRE <- function(eq, from, to, rate, ...) addReaction(eq, from, to, rate, compartment = "extraCell", ...)

.kmSmad7 <- "(Km_Smad7 + R1 + R1_int + R1_R2 + R1_R2_int + R1_R2_TGFB1 + R1_R2_TGFB1_int)"
.block   <- "(1 + k_inh_phos_PMEPA1 * PMEPA1)"
.surf    <- "(R1_R2 + R1_R2_TGFB1)"
.endo    <- "(R1_R2_int + R1_R2_TGFB1_int)"
.phos2   <- paste0("k_phos_Smad2_int * ", .endo, " * Smad2 / ((Km_phos_Smad2 + ", .endo, ") * ", .block, ")")
.phos3   <- paste0("(k_phos_Smad3_surface * ", .surf, " + k_phos_Smad3_int * ", .endo, ") * Smad3 / ", .block)
# the ligand-bound surface complex signals into the non-canonical branches
.lig     <- "R1_R2_TGFB1"
.CPLX <- list(C234 = c("pSmad2", "pSmad3", "Smad4"),
              C224 = c("pSmad2", "pSmad2", "Smad4"),
              C334 = c("pSmad3", "pSmad3", "Smad4"))
.pS2  <- "pSmad2 + C234 + 2 * C224"
.pS3  <- "pSmad3 + C234 + 2 * C334"
.cplx <- "C234 + C224 + C334"

reactions <- eqnlist() |>
  # inhibitors: dummy states, set by the events at t = -30
  addRC("", "bool_ActD",  "0") |>
  addRC("", "bool_CHX",   "0") |>
  addRC("", "bool_MG132", "0") |>
  addRC("", "bool_MEKi",  "0") |>
  addRC("", "bool_TAK1i", "0") |>
  addRC("", "bool_p38i",  "0") |>
  addRC("", "bool_IKKi",  "0") |>
  addRC("", "bool_PI3Ki", "0") |>
  addRC("", "bool_mTORi", "0") |>
  addRC("", "bool_ROCKi", "0") |>

  # ––– receptors –––
  addRC("", "R1mRNA", "k_transcr_R1mRNA * (1 + k_act_R1mRNA_FB3 * FB3) * (1 - bool_ActD)") |>
  addRC("", "R2mRNA", "k_transcr_R2mRNA * (1 - bool_ActD)") |>
  addRC("", "R1mRNA_OE", "k_transcr_R1mRNA_OE * plasmid_R1 * (1 - bool_ActD)") |>
  addRC("", "R2mRNA_OE", "k_transcr_R2mRNA_OE * plasmid_R2 * (1 - bool_ActD)") |>
  addRC("R1mRNA", "", "(k_dg_R1mRNA + k_si_R1mRNA * siRNA_R1) * R1mRNA") |>
  addRC("R2mRNA", "", "(k_dg_R2mRNA + k_si_R2mRNA * siRNA_R2 + k_dg_R2mRNA_FB4 * FB4) * R2mRNA") |>
  addRC("R1mRNA_OE", "", "k_dg_R1mRNA * R1mRNA_OE") |>
  addRC("R2mRNA_OE", "", "(k_dg_R2mRNA + k_dg_R2mRNA_FB4 * FB4) * R2mRNA_OE") |>
  addRC("", "R1", "(k_transl_R1 * R1mRNA + k_transl_R1_OE * R1mRNA_OE) * (1 - bool_CHX)") |>
  addRC("", "R2", "k_transl_R2 * R2mRNA * (1 - bool_CHX)") |>
  addRC("", "R2_ER", "k_transl_R2_OE * R2mRNA_OE * (1 - bool_CHX)") |>
  addRC("R2_ER", "R2", "k_exp_R2_ER * R2_ER") |>
  addRC("R2_ER", "", "k_dg_R2_ER * R2_ER") |>
  addRC("R2", "", "k_dg_R2_Transf * (plasmid_R1 + plasmid_R2) * R2") |>

  # ––– trafficking –––
  addRC("R1", "R1_int", "k_int_R * R1") |>
  addRC("R1_int", "R1", "k_rec_R * R1_int") |>
  addRC("R1_int", "", "k_dg_R_int * R1_int") |>
  addRC("R2", "R2_int", "k_int_R * R2") |>
  addRC("R2_int", "R2", "k_rec_R * R2_int") |>
  addRC("R2_int", "", "k_dg_R_int * R2_int") |>
  addRC("R1 + R2", "R1_R2", "k_form_R1_R2 * R1 * R2") |>
  addRC("R1_R2", "R1 + R2", "k_dec_R1_R2 * R1_R2") |>
  addRC("R1_R2", "R1_R2_int", "k_int_R * R1_R2") |>
  addRC("R1_R2_int", "R1_R2", "k_rec_R * R1_R2_int") |>
  addRC("R1_R2_int", "R1_int + R2_int", "k_dec_R1_R2_int * R1_R2_int") |>

  # ––– ligand capture and store –––
  addRC("R2 + TGFB1_ext", "R2_TGFB1", "k_form_R2_TGFB1 * R2 * TGFB1_ext", rateCompartment = "Cell") |>
  addRC("R2_TGFB1", "R2 + TGFB1_ext", "k_dec_R2_TGFB1 * R2_TGFB1", rateCompartment = "Cell") |>
  addRC("R1 + R2_TGFB1", "R1_R2_TGFB1", "k_form_R1_R2_TGFB1 * R1 * R2_TGFB1") |>
  addRC("R1_R2 + TGFB1_ext", "R1_R2_TGFB1", "k_form_R1_R2_TGFB1_R1_R2 * R1_R2 * TGFB1_ext", rateCompartment = "Cell") |>
  addRC("R1_R2_TGFB1", "R1_R2_TGFB1_int", "k_int_R * R1_R2_TGFB1") |>
  addRC("R1_R2_TGFB1_int", "R1_R2_TGFB1", "k_rec_R * R1_R2_TGFB1_int") |>
  addRC("R1_R2_TGFB1_int", "R1_int + R2_TGFB1_int", "k_dec_R1_R2_TGFB1_int * R1_R2_TGFB1_int") |>
  addRC("R2_TGFB1_int", "R2_TGFB1", "k_rec_R2_TGFB1_int * R2_TGFB1_int") |>
  addRC("R2_TGFB1_int", "R2_int", "k_dg_TGFB1_int * R2_TGFB1_int / (Km_dg_TGFB1_int + R2_TGFB1_int) * (1 - bool_MG132)") |>

  # ––– Smad7 removes TbRI from every pool –––
  addRC("R1 + Smad7", "", paste0("k_dg_R1_Smad7 * R1 * Smad7 / ", .kmSmad7, " * (1 - bool_MG132)")) |>
  addRC("R1_int + Smad7", "", paste0("k_dg_R1_Smad7 * R1_int * Smad7 / ", .kmSmad7, " * (1 - bool_MG132)")) |>
  addRC("R1_R2 + Smad7", "", paste0("k_dg_R1_Smad7 * R1_R2 * Smad7 / ", .kmSmad7, " * (1 - bool_MG132)")) |>
  addRC("R1_R2_int + Smad7", "", paste0("k_dg_R1_Smad7 * R1_R2_int * Smad7 / ", .kmSmad7, " * (1 - bool_MG132)")) |>
  addRC("R1_R2_TGFB1 + Smad7", "", paste0("k_dg_R1_Smad7 * R1_R2_TGFB1 * Smad7 / ", .kmSmad7, " * (1 - bool_MG132)")) |>
  addRC("R1_R2_TGFB1_int + Smad7", "", paste0("k_dg_R1_Smad7 * R1_R2_TGFB1_int * Smad7 / ", .kmSmad7, " * (1 - bool_MG132)")) |>

  # ––– R-Smad phosphorylation, PMEPA1 block on both –––
  addRC("Smad2", "pSmad2", .phos2) |>
  addRC("Smad3", "pSmad3", .phos3) |>
  addRC("pSmad2", "Smad2", "k_dephos_pSmad2 * pSmad2") |>
  addRC("pSmad3", "Smad3", "k_dephos_pSmad3 * pSmad3") |>

  # ––– FB3 / FB4, unmeasured, C234-driven –––
  addRC("", "FB3mRNA", "k_transcr_FB3mRNA_C234 * C234 * (1 - bool_ActD)") |>
  addRC("", "FB4mRNA", "k_transcr_FB4mRNA_C234 * C234 * (1 - bool_ActD)") |>
  addRC("FB3mRNA", "", "k_dg_FB3mRNA * FB3mRNA") |>
  addRC("FB4mRNA", "", "k_dg_FB4mRNA * FB4mRNA") |>
  addRC("", "FB3", "k_transl_FB3 * FB3mRNA * (1 - bool_CHX)") |>
  addRC("", "FB4", "k_transl_FB4 * FB4mRNA * (1 - bool_CHX)") |>
  addRC("FB3", "", "k_dg_FB3 * FB3") |>
  addRC("FB4", "", "k_dg_FB4 * FB4")

for (cn in names(.CPLX)) {
  su <- .CPLX[[cn]]
  reactions <- reactions |>
    addRC(paste(su, collapse = " + "), cn, paste(c(paste0("k_form_", cn), su), collapse = " * ")) |>
    addRC(cn, paste(su, collapse = " + "), paste0("k_dec_", cn, " * ", cn))
  for (ps in unique(grep("^pSmad", su, value = TRUE))) {
    n    <- sum(su == ps)
    rest <- su; rest[match(ps, rest)] <- sub("^p", "", ps)
    reactions <- addRC(reactions, cn, paste(rest, collapse = " + "),
                       paste0(if (n > 1L) paste0(n, " * "), "k_dephos_RSmad_cplx * ", cn))
  }
}


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Non-canonical branches
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# phosphorylation cycle: from -> to by kinase, back by a phosphatase (extra
# dephosphorylating species optional); an inhibitor multiplies the kinase term
cycle <- function(eq, from, to, kinase, inhibitor = NULL, phosphatase = NULL) {
  inh <- if (is.null(inhibitor)) "" else paste0(" * (1 - ", inhibitor, ")")
  back <- paste0("k_dephos_", to,
                 if (!is.null(phosphatase)) paste0(" + k_dephos_", to, "_", phosphatase, " * ", phosphatase))
  sat <- function(x, km) if (.mmCycles) paste0(x, "/(", km, " + ", x, ")") else x
  eq |>
    addRC(from, to, paste0("k_phos_", from, " * ", kinase, " * ", sat(from, paste0("Km_phos_", from)), inh)) |>
    addRC(to, from, paste0("(", back, ") * ", sat(to, paste0("Km_dephos_", to))))
}

reactions <- reactions |>
  # ––– ERK: ShcA pY317, Ras-GTP, Raf pS338 (ERK feedback), MEK pS217/221, ERK pT202/Y204 (DUSP6) –––
  cycle("ShcA", "pShcA", .lig) |>
  addRC("RasGDP", "RasGTP", "k_act_RasGDP * pShcA * RasGDP") |>
  addRC("RasGTP", "RasGDP", "k_inact_RasGTP * RasGTP") |>
  cycle("Raf", "pRaf", "RasGTP", phosphatase = "ppERK") |>
  cycle("MEK", "pMEK", "pRaf") |>
  cycle("ERK", "ppERK", "pMEK", inhibitor = "bool_MEKi", phosphatase = "DUSP6") |>

  # ––– TRAF6 K63 ubiquitination (A20 removes it), TAK1 pT187 –––
  addRC("TRAF6", "ubTRAF6", paste0("k_ub_TRAF6 * ", .lig, " * TRAF6")) |>
  addRC("ubTRAF6", "TRAF6", "(k_deub_ubTRAF6 + k_deub_ubTRAF6_A20 * A20) * ubTRAF6") |>
  cycle("TAK1", "pTAK1", "ubTRAF6", inhibitor = "bool_TAK1i") |>

  # ––– JNK and p38: MKK4 pS257, JNK pT183/Y185, MKK3/6 pS189/207, p38 pT180/Y182, c-Jun pS63 –––
  cycle("MKK4", "pMKK4", "pTAK1") |>
  cycle("JNK", "pJNK", "pMKK4") |>
  cycle("MKK36", "pMKK36", "pTAK1") |>
  cycle("p38", "pp38", "pMKK36") |>
  cycle("JUN", "pJUN", "pJNK") |>

  # ––– NF-kB: IKK pS176/180, IkBa phosphorylated and degraded, p65 pS536 –––
  cycle("IKK", "pIKK", "pTAK1") |>
  addRC("IkBa_p65", "p65", "k_dg_IkBa_p65 * pIKK * IkBa_p65 * (1 - bool_IKKi) * (1 - bool_MG132)") |>
  addRC("IkBa + p65", "IkBa_p65", "k_form_IkBa_p65 * IkBa * p65") |>
  cycle("p65", "pp65", "pIKK", inhibitor = "bool_IKKi") |>

  # ––– PI3K/AKT/mTOR: p85 on TRAF6, AKT pS473, mTOR pS2448, S6K pT389, S6 pS235/236 –––
  addRC("PI3K", "PI3Ka", "k_act_PI3K * ubTRAF6 * PI3K") |>
  addRC("PI3Ka", "PI3K", "k_inact_PI3Ka * PI3Ka") |>
  cycle("AKT", "pAKT", "PI3Ka", inhibitor = "bool_PI3Ki") |>
  cycle("mTOR", "pmTOR", "pAKT") |>
  cycle("S6K", "pS6K", "pmTOR", inhibitor = "bool_mTORi") |>
  cycle("S6", "pS6", "pS6K") |>
  # AKT holds Smad3 away from the receptor
  addRC("pAKT + Smad3", "pAKT_Smad3", "k_form_pAKT_Smad3 * pAKT * Smad3") |>
  addRC("pAKT_Smad3", "pAKT + Smad3", "k_dec_pAKT_Smad3 * pAKT_Smad3") |>

  # ––– Rho: Par6 pS345 by TbRII, Smurf1 on pPar6 degrades RhoA-GTP, ROCK on MLC pS19 –––
  cycle("Par6", "pPar6", .lig) |>
  addRC("", "RhoAGDP", "k_syn_RhoA * (1 - bool_CHX)") |>
  addRC("RhoAGDP", "RhoAGTP", paste0("k_act_RhoAGDP * ", .lig, " * RhoAGDP")) |>
  addRC("RhoAGTP", "RhoAGDP", "k_inact_RhoAGTP * RhoAGTP") |>
  addRC("RhoAGDP", "", "k_dg_RhoA * RhoAGDP") |>
  addRC("RhoAGTP", "", "(k_dg_RhoA + k_dg_RhoAGTP_pPar6 * pPar6 * (1 - bool_MG132)) * RhoAGTP") |>
  cycle("MLC", "pMLC", "RhoAGTP", inhibitor = "bool_ROCKi") |>

  # ––– Smad3 linker pS208/S213 by ERK, JNK and p38: out of the trimers –––
  addRC("pSmad3", "pSmad3L", paste0("(k_phos_pSmad3_ppERK * ppERK + k_phos_pSmad3_pJNK * pJNK + ",
                                    "k_phos_pSmad3_pp38 * pp38 * (1 - bool_p38i)) * pSmad3")) |>
  addRC("pSmad3L", "pSmad3", "k_dephos_pSmad3L * pSmad3L")


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Gene layer
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# drivers add to the basal transcription; ID1 is repressed by its own delayed
# chain. prot names the protein state, "" for an mRNA-only gene
.genes <- data.frame(
  gene = c("SMAD7", "TGFB1", "PMEPA1", "ID1", "JUNB", "SERPINE1", "FN1", "MMP2", "CDKN1A",
           "FOS", "EGR1", "DUSP6", "JUN", "NFKBIA", "TNFAIP3", "CCN2", "ACTA2", "SNAI1", "IL6", "ATF3",
           "CDH1", "CDH2", "VIM", "ZEB1", "COL1A1", "TGFBI", "SKIL", "BHLHE40", "LAMC2", "HMGA2"),
  drive = c("C334", "C334+pJUN", "C334+pJUN", "C334", "C334", "C234+pJUN", "C234+pJUN", "C224+pJUN", "C334",
            "ppERK", "ppERK", "ppERK", "pJUN", "pp65", "pp65", "C234+RhoAGTP", "C234+RhoAGTP", "C334+pp65",
            "pp65+pp38", "pp38+pJUN",
            "C334", "C234+pJUN", "C234+pp65", "C234+pp65", "C234+RhoAGTP", "C234", "C334",
            "C334+ppERK", "C234+ppERK", "C234+ppERK"),
  prot = c("Smad7", if (.autocrine) "TGFB1" else "", "PMEPA1", "", "JUNB", "SERPINE1", "FN1", "", "CDKN1A",
           "", "", "DUSP6", "JUN", "IkBa", "A20", "", "", "", "", "",
           "CDH1", "CDH2", "VIM", "", "", "", "", "", "", ""),
  stringsAsFactors = FALSE)
.PROTEASOMAL <- c("Smad7", "JUNB", "JUN", "IkBa")
.nREP <- 3L

for (i in seq_len(nrow(.genes))) {
  gn  <- .genes$gene[i]
  drv <- strsplit(.genes$drive[i], "+", fixed = TRUE)[[1]]
  K <- paste0("K_transcr_", gn, "mRNA_", drv)
  act <- switch(.transcription,
    linear = drv,
    mm     = paste0(drv, "/(", K, " + ", drv, ")"),
    hill   = paste0(drv, "^2/(", K, "^2 + ", drv, "^2)"))
  num <- paste(c(paste0("k_transcr_", gn, "mRNA_basal"),
                 paste0("k_transcr_", gn, "mRNA_", drv, " * ", act,
                        ifelse(drv == "pp38", " * (1 - bool_p38i)", ""))),
               collapse = " + ")
  den <- if (gn == "ID1") paste0(" / (1 + k_rep_ID1mRNA * ID1rep", .nREP, ")") else ""
  reactions <- reactions |>
    addRC("", paste0(gn, "mRNA"), paste0("(", num, ")", den, " * (1 - bool_ActD)")) |>
    addRC(paste0(gn, "mRNA"), "", paste0("k_dg_", gn, "mRNA * ", gn, "mRNA"))
  ps <- .genes$prot[i]
  if (!nzchar(ps)) next
  reactions <- addRC(reactions, "", ps, paste0("k_transl_", ps, " * ", gn, "mRNA * (1 - bool_CHX)"))
  mg <- if (ps %in% .PROTEASOMAL) " * (1 - bool_MG132)" else ""
  # c-Jun turns over in both phospho forms; secreted TGFB1 leaves by secretion only
  if (ps == "JUN")
    reactions <- reactions |>
      addRC("JUN", "", paste0("k_dg_JUN * JUN", mg)) |>
      addRC("pJUN", "", paste0("k_dg_JUN * pJUN", mg))
  else if (ps != "TGFB1")
    reactions <- addRC(reactions, ps, "", paste0("k_dg_", ps, " * ", ps, mg))
}

# ID1 repressor chain, driven by C334
.rp <- paste0("ID1rep", seq_len(.nREP))
reactions <- addRC(reactions, "", .rp[1], "k_transit_ID1rep * C334")
for (i in seq_len(.nREP))
  reactions <- addRC(reactions, .rp[i], if (i < .nREP) .rp[i + 1L] else "",
                     paste0("k_transit_ID1rep * ", .rp[i]))

if (.autocrine) reactions <- reactions |>
  addRE("TGFB1", "TGFB1_ext", "k_secr_TGFB1 * TGFB1", rateCompartment = "Cell")

reactions <- setCompartmentVolume(reactions, Cell = "volumeC", extraCell = "volumeEC")

# conserved signalling proteins
.pools <- list(
  tSmad2 = paste("Smad2 +", .pS2),
  tSmad3 = paste("Smad3 + pSmad3L + pAKT_Smad3 +", .pS3),
  tSmad4 = paste("Smad4 +", .cplx),
  tShcA = "ShcA + pShcA", tRas = "RasGDP + RasGTP", tRaf = "Raf + pRaf",
  tMEK = "MEK + pMEK", tERK = "ERK + ppERK", tTRAF6 = "TRAF6 + ubTRAF6",
  tTAK1 = "TAK1 + pTAK1", tMKK4 = "MKK4 + pMKK4", tJNK = "JNK + pJNK",
  tMKK36 = "MKK36 + pMKK36", tp38 = "p38 + pp38", tIKK = "IKK + pIKK",
  tp65 = "p65 + pp65 + IkBa_p65", tPI3K = "PI3K + PI3Ka", tAKT = "AKT + pAKT + pAKT_Smad3",
  tmTOR = "mTOR + pmTOR", tS6K = "S6K + pS6K", tS6 = "S6 + pS6", tPar6 = "Par6 + pPar6",
  tMLC = "MLC + pMLC")
reactions <- customTotals(reactions, .pools)
.states <- reactions$states

cat(sprintf("%d states, %d parameters\n", length(.states),
            length(setdiff(getParameters(reactions), .states))))


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Observables
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# targeted proteomics and phosphoproteomics in absolute units; mRNA with a scale
.R1tot <- "R1 + R1_int + R1_R2 + R1_R2_int + R1_R2_TGFB1 + R1_R2_TGFB1_int"
.R2tot <- "R2_ER + R2 + R2_int + R1_R2 + R1_R2_int + R2_TGFB1 + R2_TGFB1_int + R1_R2_TGFB1 + R1_R2_TGFB1_int"
.phospho <- c(
  pSmad2_S465_467 = .pS2, pSmad3_S423_425 = paste(.pS3, "+ pSmad3L"), pSmad3_S208_213 = "pSmad3L",
  pShcA_Y317 = "pShcA", pRaf_S338 = "pRaf", pMEK_S217_221 = "pMEK", pERK_T202_Y204 = "ppERK",
  pTAK1_T187 = "pTAK1", pMKK4_S257 = "pMKK4", pJNK_T183_Y185 = "pJNK", pMKK36_S189_207 = "pMKK36",
  pp38_T180_Y182 = "pp38", pJUN_S63 = "pJUN", pIKK_S176_180 = "pIKK", pp65_S536 = "pp65",
  pAKT_S473 = "pAKT + pAKT_Smad3", pmTOR_S2448 = "pmTOR", pS6K_T389 = "pS6K", pS6_S235_236 = "pS6",
  pPar6_S345 = "pPar6", pMLC_S19 = "pMLC")
.protein <- c(
  TGFBR1 = .R1tot, TGFBR2 = .R2tot, SMAD2 = paste("Smad2 +", .pS2),
  SMAD3 = paste("Smad3 + pSmad3L + pAKT_Smad3 +", .pS3), SMAD4 = paste("Smad4 +", .cplx),
  SMAD4_CoIP = .cplx, SMAD7 = "Smad7", PMEPA1 = "PMEPA1", JUNB = "JUNB", SERPINE1 = "SERPINE1",
  FN1 = "FN1", CDKN1A = "CDKN1A", DUSP6 = "DUSP6", JUN = "JUN + pJUN", NFKBIA = "IkBa + IkBa_p65",
  TNFAIP3 = "A20", RHOA = "RhoAGDP + RhoAGTP", CDH1 = "CDH1", CDH2 = "CDH2", VIM = "VIM",
  TGFB1_medium = "TGFB1_ext")
.mrna <- setNames(paste0("scale_", .genes$gene, "mRNA * ", .genes$gene, "mRNA"), .genes$gene)
.mrna <- c(.mrna, TGFBR1 = "scale_R1mRNA * (R1mRNA + R1mRNA_OE)",
           TGFBR2 = "scale_R2mRNA * (R2mRNA + R2mRNA_OE)")
observables <- as.eqnvec(c(
  setNames(.phospho, paste0(names(.phospho), "_obs")),
  setNames(.protein, paste0(names(.protein), "_prot_obs")),
  setNames(.mrna, paste0(names(.mrna), "_mRNA_obs"))))
cat(sprintf("%d observables\n", length(observables)))


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Conditions and events
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# inhibitors from t = -30, TGF-b1 (5 ng/ml, 25 kDa) at t = 0; one unstimulated control
.inhibitors <- c("ActD", "CHX", "MG132", "MEKi", "TAK1i", "p38i", "IKKi", "PI3Ki", "mTORi", "ROCKi")
.perts <- c("Ctrl0", "Ctrl", .inhibitors, "R1Knd", "R2Knd", "OE_R1", "OE_R2")
conditions <- data.frame(row.names = .perts,
                         plasmid_R1 = as.numeric(.perts == "OE_R1"),
                         plasmid_R2 = as.numeric(.perts == "OE_R2"),
                         siRNA_R1   = as.numeric(.perts == "R1Knd"),
                         siRNA_R2   = as.numeric(.perts == "R2Knd"),
                         var_TGFb_dose = ifelse(.perts == "Ctrl0", 0, 5 / 25))
for (inh in .inhibitors) conditions[[paste0("var_bool_", inh)]] <- as.numeric(.perts == inh)

.forcings <- paste0("bool_", .inhibitors)
events <- eventlist() |>
  addEvent(var = "TGFB1_ext", time = 0, value = "var_TGFb_dose", method = "add")
for (inh in .inhibitors)
  events <- addEvent(events, var = paste0("bool_", inh), time = -30,
                     value = paste0("var_bool_", inh), method = "replace")


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Symmetries
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
t0 <- Sys.time()
idResult <- symmetryDetection(reactions, observables, forcings = .forcings, events = events,
                              conditions = conditions, equilibrate = TRUE, reduceCQ = TRUE,
                              gaugePreference = NULL, reconstruct = TRUE,
                              cores = .cores, verbose = TRUE)
tDetect <- Sys.time() - t0
tDetect
summary(idResult)
