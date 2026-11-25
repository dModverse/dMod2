## Two priors on overlapping parameters
obj1 <- constraintL2(mu = c(a = 0, b = 0), sigma = 1)
obj2 <- constraintL2(mu = c(b = 1, c = 1), sigma = 1)

## Gradients and Hessians are matched by parameter name
total <- obj1(c(a = 1, b = 2)) + obj2(c(b = 2, c = 3))
total$value
total$gradient
