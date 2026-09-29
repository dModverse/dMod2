# numSteady bench model: minimal TGF-beta receptor model, H1975 data.
#
# Reduced from TGFbModelling M016_H1975TGFb: receptor transcripts with siRNA and
# plasmid-borne copies, receptor trafficking, ligand binding, Smad core with
# three conserved totals, Smad7 (negative feedback on all TbRI pools, shared
# Km) and FB3 (positive feedback on TbRI transcription).
#
# The steady state differs per arm (Ctrl, R1Knd, R2Knd, OE_R1, OE_R2) through
# the parameters siRNA_R1/R2 and plasmid_R1/R2. R1mRNA_OE and R2mRNA_OE are 0
# outside their own arm.
#
# numSteadyModel(dir) returns reactions, events, observables, data, condition
# grid, t0, forcings and the perturbation parameters. With ligand = FALSE the
# ligand-bound receptor states, which are 0 at steady state, are left out.

numSteadyModel <- function(dir = "bench/numSteady", ligand = TRUE) {
  MW_TGFB <- 25   # TGF-b1 dimer, kDa: ng/ml / kDa = pmol/ml

  data <- utils::read.csv(file.path(dir, "data.csv"), sep = ";")
  data$time <- as.numeric(data$time)
  data$lloq <- as.numeric(data$lloq)
  data$lloq[is.na(data$lloq)] <- -Inf
  datalist <- dMod2::as.datalist(data)

  grid <- attr(datalist, "condition.grid")
  grid$var_bool_ActD  <- as.numeric(grid$Pertubation == "ActD")
  grid$var_bool_CHX   <- as.numeric(grid$Pertubation == "CHX")
  grid$var_bool_MG132 <- as.numeric(grid$Pertubation == "MG132")
  grid$plasmid_R1 <- as.numeric(grid$Pertubation == "OE_R1_noninducible")
  grid$plasmid_R2 <- as.numeric(grid$Pertubation == "OE_R2_noninducible")
  grid$siRNA_R1   <- as.numeric(grid$Pertubation == "R1Knd")
  grid$siRNA_R2   <- as.numeric(grid$Pertubation == "R2Knd")
  grid$var_TGFb_dose <- as.numeric(grid$init_TGFb) / MW_TGFB
  # inhibitor arms start at their event (-30), all others at their steady state
  t0 <- setNames(ifelse(grid$var_bool_ActD + grid$var_bool_CHX + grid$var_bool_MG132 > 0, -30, 0),
                 rownames(grid))

  addRC <- function(eq, from, to, rate, ...)
    dMod2::addReaction(eq, from, to, rate, compartment = "Cell", ...)
  lg <- function(x) if (ligand) x else ""
  kmSmad7  <- paste0("(Km_Smad7 + R1 + R1_int", lg(" + R1_R2_TGFB1 + R1_R2_TGFB1_int"), ")")
  basal    <- "k_enc_R1_R2 * R1 * R2"
  basalInt <- "k_enc_R1_R2_int * R1_int * R2_int"
  active   <- paste0("(", basal, " + ", basalInt, lg(" + R1_R2_TGFB1 + R1_R2_TGFB1_int"), ")")
  phos <- function(k) paste0(
    "(k_phos_Smad", k, "_R1_R2 * (", basal, lg(" + R1_R2_TGFB1"), ") + k_phos_Smad", k,
    "_R1_R2_int * (", basalInt, lg(" + R1_R2_TGFB1_int"), ")) * Smad", k,
    " / (1 + ", active, " / Km_phos_Smad", k, ")")

  r <- dMod2::eqnlist() |>
    dMod2::assignCompartment(TGFB1_ext = "extraCell", volume = "volumeEC") |>
    addRC("", "bool_ActD",  "0") |>
    addRC("", "bool_CHX",   "0") |>
    addRC("", "bool_MG132", "0") |>
    # receptor transcripts: endogenous (TbRI under FB3) and plasmid-borne
    addRC("", "R1mRNA", "k_transcr_R1mRNA * (1 + k_act_R1mRNA_FB3 * FB3) * (1 - bool_ActD)") |>
    addRC("", "R2mRNA", "k_transcr_R2mRNA * (1 - bool_ActD)") |>
    addRC("", "R1mRNA_OE", "k_transcr_R1mRNA_OE * plasmid_R1 * (1 - bool_ActD)") |>
    addRC("", "R2mRNA_OE", "k_transcr_R2mRNA_OE * plasmid_R2 * (1 - bool_ActD)") |>
    addRC("R1mRNA", "", "(k_dg_R1mRNA + k_si_R1mRNA * siRNA_R1) * R1mRNA") |>
    addRC("R2mRNA", "", "(k_dg_R2mRNA + k_si_R2mRNA * siRNA_R2) * R2mRNA") |>
    addRC("R1mRNA_OE", "", "k_dg_R1mRNA * R1mRNA_OE") |>
    addRC("R2mRNA_OE", "", "k_dg_R2mRNA * R2mRNA_OE") |>
    addRC("", "R1", "(k_transl_R1 * R1mRNA + k_transl_R1_OE * R1mRNA_OE) * (1 - bool_CHX)") |>
    addRC("", "R2", "(k_transl_R2 * R2mRNA + k_transl_R2_OE * R2mRNA_OE) * (1 - bool_CHX)") |>
    addRC("R2", "", "k_dg_R2_Transf * (plasmid_R1 + plasmid_R2) * R2") |>
    # constitutive trafficking
    addRC("R1", "R1_int", "k_int_R1 * R1") |>
    addRC("R1_int", "R1", "k_rec_R1_int * R1_int") |>
    addRC("R1_int", "", "k_dg_R1_int * R1_int") |>
    addRC("R2", "R2_int", "k_int_R2 * R2") |>
    addRC("R2_int", "R2", "k_rec_R2_int * R2_int") |>
    addRC("R2_int", "", "k_dg_R2_int * R2_int") |>
    # Smad7 removes TbRI, one capacity Km_Smad7
    addRC("R1 + Smad7", "", paste0("k_dg_R1_Smad7 * R1 * Smad7 / ", kmSmad7, " * (1 - bool_MG132)")) |>
    addRC("R1_int + Smad7", "", paste0("k_dg_R1_int_Smad7 * R1_int * Smad7 / ", kmSmad7, " * (1 - bool_MG132)")) |>
    # R-Smad phosphorylation by basal and ligand-bound R1:R2, trimer formation
    addRC("Smad2", "pSmad2", phos(2)) |>
    addRC("Smad3", "pSmad3", phos(3)) |>
    addRC("pSmad2", "Smad2", "k_dephos_pSmad2 * pSmad2") |>
    addRC("pSmad3", "Smad3", "k_dephos_pSmad3 * pSmad3") |>
    addRC("pSmad2 + pSmad3 + Smad4", "C234", "k_form_C234 * pSmad2 * pSmad3 * Smad4") |>
    addRC("C234", "pSmad2 + pSmad3 + Smad4", "k_dec_C234 * C234") |>
    addRC("C234", "Smad2 + pSmad3 + Smad4", "k_dephos_RSmad_C234 * C234") |>
    addRC("C234", "pSmad2 + Smad3 + Smad4", "k_dephos_RSmad_C234 * C234") |>
    # feedbacks: Smad7 (negative) and FB3 (positive on TbRI transcription)
    addRC("", "SMAD7mRNA", "(k_transcr_SMAD7mRNA_basal + k_transcr_SMAD7mRNA_C234 * C234) * (1 - bool_ActD)") |>
    addRC("SMAD7mRNA", "", "k_dg_SMAD7mRNA * SMAD7mRNA") |>
    addRC("", "Smad7", "k_transl_Smad7 * SMAD7mRNA * (1 - bool_CHX)") |>
    addRC("Smad7", "", "k_dg_Smad7 * Smad7 * (1 - bool_MG132)") |>
    addRC("", "FB3mRNA", "k_transcr_FB3mRNA_C234 * C234 * (1 - bool_ActD)") |>
    addRC("FB3mRNA", "", "k_dg_FB3mRNA * FB3mRNA") |>
    addRC("", "FB3", "k_transl_FB3 * FB3mRNA * (1 - bool_CHX)") |>
    addRC("FB3", "", "k_dg_FB3 * FB3")
  if (ligand) r <- r |>
    # ligand capture and complex trafficking
    addRC("R2 + TGFB1_ext", "R2_TGFB1", "k_form_R2_TGFB1 * R2 * TGFB1_ext", rateCompartment = "Cell") |>
    addRC("R2_TGFB1", "R2 + TGFB1_ext", "k_dec_R2_TGFB1 * R2_TGFB1", rateCompartment = "Cell") |>
    addRC("R1 + R2_TGFB1", "R1_R2_TGFB1", "k_form_R1_R2_TGFB1 * R1 * R2_TGFB1") |>
    addRC("R1_R2_TGFB1", "R1_R2_TGFB1_int", "k_int_R1_R2_TGFB1 * R1_R2_TGFB1") |>
    addRC("R1_R2_TGFB1_int", "R1_R2_TGFB1", "k_rec_R1_R2_TGFB1_int * R1_R2_TGFB1_int") |>
    addRC("R1_R2_TGFB1_int", "R1_int + R2_TGFB1_int", "k_dec_R1_R2_TGFB1_int * R1_R2_TGFB1_int") |>
    addRC("R2_TGFB1_int", "R2_TGFB1", "k_rec_R2_TGFB1_int * R2_TGFB1_int") |>
    addRC("R2_TGFB1_int", "R2_int", "k_dg_TGFB1_int * R2_TGFB1_int * (1 - bool_MG132)") |>
    addRC("R1_R2_TGFB1 + Smad7", "", paste0("k_dg_R1_R2_TGFB1_Smad7 * R1_R2_TGFB1 * Smad7 / ", kmSmad7, " * (1 - bool_MG132)")) |>
    addRC("R1_R2_TGFB1_int + Smad7", "", paste0("k_dg_R1_R2_TGFB1_int_Smad7 * R1_R2_TGFB1_int * Smad7 / ", kmSmad7, " * (1 - bool_MG132)"))
  r <- dMod2::setCompartmentVolume(r, Cell = "volumeC")
  r <- dMod2::customTotals(r, list(tSmad2 = "Smad2 + pSmad2 + C234",
                                   tSmad3 = "Smad3 + pSmad3 + C234",
                                   tSmad4 = "Smad4 + C234"))

  events <- dMod2::eventlist() |>
    dMod2::addEvent(var = "TGFB1_ext",  time = 0,   value = "var_TGFb_dose",  method = "add") |>
    dMod2::addEvent(var = "bool_CHX",   time = -30, value = "var_bool_CHX",   method = "replace") |>
    dMod2::addEvent(var = "bool_MG132", time = -30, value = "var_bool_MG132", method = "replace") |>
    dMod2::addEvent(var = "bool_ActD",  time = -30, value = "var_bool_ActD",  method = "replace")

  R1tot <- "R1 + R1_int + R1_R2_TGFB1 + R1_R2_TGFB1_int"
  R2tot <- "R2 + R2_int + R2_TGFB1 + R2_TGFB1_int + R1_R2_TGFB1 + R1_R2_TGFB1_int"
  observables <- dMod2::eqnvec(
    pSmad2_obs = "log2(pSmad2 + C234 + 1e-4) + offset_pSmad2",
    pSmad3_obs = "log2(pSmad3 + C234 + 1e-4) + offset_pSmad3",
    TSmad2_obs = "log2(Smad2 + pSmad2 + C234 + 1e-4)",
    TSmad3_obs = "log2(Smad3 + pSmad3 + C234 + 1e-4)",
    TSmad4_obs = "log2(Smad4 + C234 + 1e-4)",
    TGFb_obs   = paste0("log10(", MW_TGFB, " * TGFB1_ext + 1e-3) + offset_TGFb"),
    TGFBR1_mRNA_obs      = "log2(R1mRNA + R1mRNA_OE + 1e-4) + offset_TGFBR1_mRNA",
    TGFBR2_mRNA_obs      = "log2(R2mRNA + R2mRNA_OE + 1e-4) + offset_TGFBR2_mRNA",
    TGFBR1_mRNA_obs_ActD = "log2(R1mRNA + R1mRNA_OE + 1e-4) + offset_TGFBR1_mRNA_ActD",
    TGFBR2_mRNA_obs_ActD = "log2(R2mRNA + R2mRNA_OE + 1e-4) + offset_TGFBR2_mRNA_ActD",
    TGFBR1_prot_obs_MS   = paste0("log2(", R1tot, " + 1e-4)"),
    TGFBR2_prot_obs_MS   = paste0("log2(", R2tot, " + 1e-4)"),
    SMAD7_mRNA_obs       = "log2(SMAD7mRNA + 1e-4) + offset_SMAD7_mRNA")

  list(reactions = r, events = events, observables = observables, data = datalist,
       grid = grid, t0 = t0, forcings = c("bool_ActD", "bool_CHX", "bool_MG132"),
       perturbation = c("plasmid_R1", "plasmid_R2", "siRNA_R1", "siRNA_R2"))
}
