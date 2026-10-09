mydata1 <- data.frame(
  name = "A",
  time = 0:1,
  value = 1:2,
  sigma = .1,
  compound = c("DEM", "APAP"),
  dose = "0.1"
)

mydata2 <- data.frame(
  name = "A",
  time = 0:1,
  value = 3:4,
  sigma = .1,
  compound = c("APAP", "DCF"),
  dose = "0.1"
)

data1 <- as.datalist(mydata1, splitBy = c("compound", "dose"))
data2 <- as.datalist(mydata2, splitBy = c("compound", "dose"))

## Condition APAP_0.1 is taken from data2, with a warning
data <- data1 + data2
data
attr(data, "condition.grid")
