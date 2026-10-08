## Parameter vector without derivatives
p <- as.parvec(c(a = 1, b = 2, c = 3))
p

## Jacobian and second derivatives
J <- matrix(c(1, 0, 0,
              0, 2, 0,
              0, 0, 3),
            nrow = 3, byrow = TRUE,
            dimnames = list(c("a", "b", "c"), c("x", "y", "z")))

H <- array(0, dim = c(3, 3, 3),
           dimnames = list(c("a", "b", "c"),
                           c("x", "y", "z"),
                           c("x", "y", "z")))
H["a", "x", "x"] <- 1
H["b", "y", "y"] <- 2
H["c", "z", "z"] <- 3

p2 <- as.parvec(c(a = 10, b = 20, c = 30), deriv = J, deriv2 = H)
p2
getDerivs(p2)
getDerivs2(p2)

## Subsetting keeps the matching derivatives
p_sub <- p2[c("a", "b")]
getDerivs(p_sub)

## Concatenation
p3 <- as.parvec(c(d = 4, e = 5))
c(p2, p3)

## From numbers
parvec(a = 1, b = 2)
