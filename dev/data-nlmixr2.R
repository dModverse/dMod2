# Writes data/warfarin.tab and data/mavoglurant.tab from nlmixr2data (GPL >= 3): one row per
# reading, name, time, value, sigma, condition, dose and covariates as columns.
# Rscript dev/data-nlmixr2.R from the package root, with nlmixr2data installed.

w <- nlmixr2data::warfarin
w <- w[w$evid == 0 & !(w$dvid == "cp" & w$dv == 0), ]
dose <- tapply(nlmixr2data::warfarin$amt[nlmixr2data::warfarin$evid == 1],
               nlmixr2data::warfarin$id[nlmixr2data::warfarin$evid == 1], sum)
warfarin <- data.frame(name = as.character(w$dvid), time = w$time, value = w$dv, sigma = NA_real_,
                       condition = paste0("ID", w$id), id = w$id, dose = as.numeric(dose[as.character(w$id)]),
                       wt = w$wt, age = w$age, sex = as.character(w$sex))
utils::write.table(warfarin, "data/warfarin.tab", row.names = FALSE)

m <- nlmixr2data::mavoglurant
doses <- m[m$EVID == 1, ]
key <- paste(doses$ID, doses$OCC, sep = "_")
obs <- m[m$EVID == 0 & m$MDV == 0, ]
k <- match(paste(obs$ID, obs$OCC, sep = "_"), key)
mavoglurant <- data.frame(name = "cp", time = obs$TIME, value = obs$DV, sigma = NA_real_,
                          condition = key[k], ID = obs$ID, OCC = obs$OCC, dose = doses$AMT[k],
                          duration = doses$AMT[k] / doses$RATE[k], AGE = obs$AGE, SEX = obs$SEX,
                          WT = obs$WT, HT = obs$HT)
utils::write.table(mavoglurant, "data/mavoglurant.tab", row.names = FALSE)
