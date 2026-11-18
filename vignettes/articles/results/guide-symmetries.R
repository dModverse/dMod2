# The discrete symmetries of the two-compartment model of Godfrey, for guide chapter 2,
# guide-symmetries.Rmd, which shows the labeled part and reads the result from
# guide-symmetries.rds. The search needs msolve. Run in this folder.

library(dMod2)

## ---- godfrey
f <- eqnvec(x1 = "-(k21 + k01)*x1 + k12*x2 + u", x2 = "k21*x1 - k12*x2")
compartments <- symmetryDetection(f, eqnvec(y = "x2"), knownInputs = "u", discrete = TRUE)

## ---- save
# the md5 of the script, which the article checks on reading
attr(compartments, "stamp") <- list(md5 = tools::md5sum("guide-symmetries.R"),
                                    dMod2 = packageVersion("dMod2"))
saveRDS(compartments, "guide-symmetries.rds")
