# Discrete Symmetries

Searches the discrete symmetries of a
[`symmetryDetection()`](https://dmodverse.github.io/dMod2/reference/symmetryDetection.md)
result: other points of the domain with the output of a given one. The
Taylor coefficients of the observables up to `order` are equated with
their values at a point; msolve counts the solutions over the algebraic
closure modulo a prime, and over the reals and in the declared domain at
`points` random rational points, where the number in the domain may
differ between regions. Where the other points are rational functions
\\\sigma\\ of the coordinates, \\\sigma\\ is reconstructed, checked on
the output at random points (`"checked"`) and, where every state has a
free initial value, proven by a map of the states (`"exact"`). A model
with directions is reduced by
[`symmetryReduction()`](https://dmodverse.github.io/dMod2/reference/symmetryReduction.md)
first and the reduced model is analyzed.

## Usage

``` r
discreteSymmetries(
  object,
  order = NULL,
  points = 5L,
  timeout = 600,
  budget = 300,
  cores = 1
)
```

## Arguments

- object:

  A `symmetrydetection` result from
  [`symmetryDetection()`](https://dmodverse.github.io/dMod2/reference/symmetryDetection.md)
  of a model with one condition and no events, switches or unknown
  inputs, detected without `equilibrate` and `reduceCQ`.

- order:

  Highest Taylor coefficient; `NULL` (default) twice the number of
  coordinates.

- points:

  Number of random points of the count in the domain, or a list of
  points, each a named list of the coordinates. Defaults to 5.

- timeout:

  Seconds one call of msolve may take. Defaults to 600.

- budget:

  Seconds the reconstruction of the maps may take; past it, the counts
  are reported without closed forms. Defaults to 300.

- cores:

  Threads of msolve. Defaults to 1.

## Value

`object` with `discreteSymmetries`, per \\\sigma\\ the `map` of the
coordinates it changes (an
[eqnvec](https://dmodverse.github.io/dMod2/reference/eqnvec.md)),
`changes`, `level`, `preservesDomain` (`TRUE`, `FALSE`, or `NULL` with
the `domainConditions` under which \\\sigma(\theta)\\ lies in the
domain) and `fixedSet`, the factors on which \\\sigma\\ fixes
\\\theta\\; and `info$discrete` with `order`, `count` (the solutions
over the algebraic closure at a random point modulo a prime), `counts`
per point (`complex`, `real`, `inDomain`, `undecided`), `closedForms`,
`countFinal` (every solution is the image of a proven \\\sigma\\),
`fundamentalDomain` (inequalities that keep one point of each orbit),
`limits` and, for a model with directions, `reduced`.

## Details

The search needs msolve; see
[`symmetryDetection()`](https://dmodverse.github.io/dMod2/reference/symmetryDetection.md)
for its path.

## References

Godfrey, K. (1983). Compartmental Models and Their Application. Academic
Press, London.

## See also

[`symmetryDetection()`](https://dmodverse.github.io/dMod2/reference/symmetryDetection.md),
[`symmetryReduction()`](https://dmodverse.github.io/dMod2/reference/symmetryReduction.md)

## Examples

``` r
# two compartments, an input u into the first, the second observed
f <- eqnvec(x1 = "-(k21 + k01)*x1 + k12*x2 + u", x2 = "k21*x1 - k12*x2")
res <- symmetryDetection(f, eqnvec(y = "x2"), knownInputs = "u", discrete = TRUE)
#> Error: discreteSymmetries() needs msolve: confirm the msolve installation in an interactive session, or set SYMIDENT_MSOLVE to an msolve executable.
res
#> function (data, out, err = NULL) 
#> {
#>     obs <- as.character(data$name)
#>     n <- nrow(data)
#>     times <- sort(unique(data$time))
#>     names <- unique(obs)
#>     ti <- .matchNum(times, out[, 1])[.matchNum(data$time, times)]
#>     ni <- match(names, colnames(out))[match(obs, names)]
#>     if (anyNA(ni)) 
#>         stop("Observable not found: ", paste(setdiff(names, colnames(out)), 
#>             collapse = ", "))
#>     if (anyNA(ti)) 
#>         stop("Some data$time not found in out[,1]")
#>     pred <- out[cbind(ti, ni)]
#>     deriv <- NULL
#>     if (!is.null(d <- attr(out, "deriv"))) {
#>         oi <- match(obs, dimnames(d)[[2]])
#>         np <- dim(d)[3]
#>         deriv <- matrix(d[cbind(rep(ti, np), rep(oi, np), rep(seq_len(np), 
#>             each = n))], n, np, dimnames = list(NULL, dimnames(d)[[3]]))
#>     }
#>     deriv2 <- NULL
#>     if (!is.null(d2 <- attr(out, "deriv2"))) {
#>         oi2 <- match(obs, dimnames(d2)[[2]])
#>         np2 <- dim(d2)[3]
#>         idx <- cbind(rep(ti, np2 * np2), rep(oi2, np2 * np2), 
#>             rep(rep(seq_len(np2), each = n), np2), rep(seq_len(np2), 
#>                 each = n * np2))
#>         deriv2 <- array(d2[idx], c(n, np2, np2), dimnames = list(NULL, 
#>             dimnames(d2)[[3]], dimnames(d2)[[4]]))
#>     }
#>     sig <- data$sigma
#>     sNA <- is.na(sig)
#>     derr <- NULL
#>     derr2 <- NULL
#>     if (any(sNA)) {
#>         if (is.null(err)) 
#>             stop("NA sigmas but no errmodel")
#>         ti_e <- .matchNum(times, err[, 1])[.matchNum(data$time, 
#>             times)]
#>         ni_e <- match(names, colnames(err))[match(obs, names)]
#>         sig[sNA] <- err[cbind(ti_e, ni_e)][sNA]
#>         if (!is.null(de <- attr(err, "deriv"))) {
#>             oi <- match(obs, dimnames(de)[[2]])
#>             np <- dim(de)[3]
#>             ns <- sum(sNA)
#>             derr <- matrix(0, n, np, dimnames = list(NULL, dimnames(de)[[3]]))
#>             derr[sNA, ] <- matrix(de[cbind(rep(ti_e[sNA], np), 
#>                 rep(oi[sNA], np), rep(seq_len(np), each = ns))], 
#>                 ns, np)
#>         }
#>         if (!is.null(de2 <- attr(err, "deriv2"))) {
#>             oi <- match(obs, dimnames(de2)[[2]])
#>             np2 <- dim(de2)[3]
#>             ns <- sum(sNA)
#>             derr2 <- array(0, c(n, np2, np2), dimnames = list(NULL, 
#>                 dimnames(de2)[[3]], dimnames(de2)[[4]]))
#>             idx <- cbind(rep(ti_e[sNA], np2 * np2), rep(oi[sNA], 
#>                 np2 * np2), rep(rep(seq_len(np2), each = ns), 
#>                 np2), rep(seq_len(np2), each = ns * np2))
#>             derr2[sNA, , ] <- array(de2[idx], c(ns, np2, np2))
#>         }
#>     }
#>     val <- pmax(data$value, data$lloq)
#>     resi <- pred - val
#>     inv <- 1/sig
#>     objframe(data.table::data.table(time = data$time, name = obs, 
#>         value = val, prediction = pred, sigma = sig, residual = resi, 
#>         weighted.residual = resi * inv, bloq = val <= data$lloq, 
#>         weighted.0 = pred * inv), deriv = deriv, deriv.err = derr, 
#>         deriv2 = deriv2, deriv2.err = derr2)
#> }
#> <bytecode: 0x5594ad5312c8>
#> <environment: namespace:dMod2>
res$discreteSymmetries[[1]]$map
#> Error in res$discreteSymmetries: object of type 'closure' is not subsettable
res$info$discrete$fundamentalDomain
#> Error in res$info: object of type 'closure' is not subsettable
```
