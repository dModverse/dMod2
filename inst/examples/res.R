## Prediction with the derivative of the observable y in the parameter k
times <- 0:2
J <- array(-times * exp(-times), c(3, 1, 1), list(NULL, "y", "k"))
out <- prdframe(cbind(time = times, y = exp(-times)), deriv = J, parameters = c(k = 1))

data <- data.frame(time = times, name = "y", value = c(1.05, 0.35, 0.15),
                   sigma = 0.1, lloq = -Inf)
r <- res(data, out)
r
attr(r, "deriv")
