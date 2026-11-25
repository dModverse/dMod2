el <- addReaction(eqnlist(), "A", "B", "k1*A", description = "conversion")
file <- tempfile(fileext = ".csv")
write.eqnlist(el, file = file)
read.csv(file)
