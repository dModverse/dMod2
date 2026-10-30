## A straight line with two gross outliers: the Student-t fit stays on the line
x <- Xd(data.frame(name = "y", time = 0:9, row.names = paste0("a", 0:9)), condition = "C1")
set.seed(1)
y <- 1 + 0.5 * (0:9) + rnorm(10, 0, 0.1)
y[c(3, 8)] <- y[c(3, 8)] + c(3, -4)
data <- as.datalist(data.frame(name = "y", time = 0:9, value = y, sigma = 0.1, condition = "C1"))
pars <- setNames(1 + 0.5 * (0:9), paste0("a", 0:9))
obj <- normT(data, x, df = 3)
obj(pars)$value
