## Generate datalist from scratch
mydata1 <- data.frame(name = "A",
                      time = 0:5,
                      value = 0:5,
                      sigma = .1,
                      lloq = -0.5)

mydata2 <- data.frame(name = "A",
                      time = 0:5,
                      value = sin(0:5),
                      sigma = .1)

data <- datalist(C1 = mydata1, C2 = mydata2)
print(data)
plot(data)

## Split one data.frame into conditions, keeping a covariate per condition
t <- seq(0, 2*pi, length.out = 20)
long <- data.frame(name = "A", time = c(t, t), sigma = 0.1,
                   amplitude = rep(c(1, 1.5), each = 20),
                   phase = rep(c("early", "late"), each = 20),
                   stage = rep(rep(c("upper", "lower"), each = 10), 2))
long$value <- long$amplitude * sin(long$time)

data <- as.datalist(long, splitBy = c("stage", "phase"), keepCovariates = "amplitude")
print(data)
plot(data)
print(covariates(data))
