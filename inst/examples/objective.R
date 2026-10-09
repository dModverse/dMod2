## Three objective functions
prior <- structure(rep(0, 5), names = letters[1:5])

obj1 <- constraintL2(mu = prior, attrName = "center")
obj2 <- constraintL2(mu = prior + 1, attrName = "right")
obj3 <- constraintL2(mu = prior - 1, attrName = "left")

## Evaluate the first objective function on a random vector
set.seed(1)
pouter <- prior + rnorm(length(prior))
obj1(pouter)

## Split into fixed and non-fixed part
fixed <- pouter[4:5]
pouter <- pouter[1:3]
obj1(pouter, fixed = fixed)

## Fit and profile
myfit <- trust(obj1, pouter, rinit = 1, rmax = 10, fixed = fixed)
myprof <- profile(obj1, myfit$argument, "a", fixed = fixed)
plotProfile(myprof)

## Sum of objective functions
obj <- obj1 + obj2 + obj3
obj(prior)
