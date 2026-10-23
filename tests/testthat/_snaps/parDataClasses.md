# summary of a parlist prints best and worst fit and the outcome statistics

    Code
      summary(.fitsFixture())
    Output
      Best and worst fit
            index value converged iterations stopReason neval iterations.first
      best      3     1      TRUE          3   gradient     4               NA
      worst     1     3      TRUE          5   gradient     7                2
            iterations.second elapsed a b  c
      best                 NA      NA 2 3 NA
      worst                 3     0.5 1 2 NA
      
      Statistics of fit outcome
      Fits aborted:       2
      Fits not converged: 1
      Fits converged:     2
      Fits total:         5
      
      
      Termination reason
      gradient            2
      iterlim             1

---

    Code
      summary(parlist(structure("boom", class = "try-error")))
    Output
      
      Statistics of fit outcome
      Fits aborted:       1
      Fits not converged: 0
      Fits converged:     0
      Fits total:         1

---

    Code
      summary(parlist())
    Output
      
      Statistics of fit outcome
      Fits aborted:       0
      Fits not converged: 0
      Fits converged:     0
      Fits total:         0

# as.parframe of a parlist

    Code
      .canonDput(as.parframe(fits))
    Output
      structure(list(index = c(3L, 4L, 1L), value = c(1, 2, 3), converged = c(TRUE, 
      FALSE, TRUE), iterations = c(3L, 10L, 5L), stopReason = c("gradient", 
      "iterlim", "gradient"), neval = c(4L, NA, 7L), iterations.first = c(NA, 
      NA, 2L), iterations.second = c(NA, NA, 3L), elapsed = c(NA, NA, 
      0.5), a = c(2, 5, 1), b = c(3, NA, 2), c = c(NA, 9, NA)), class = c("parframe", 
      "data.frame"), metanames = c("index", "value", "converged", "iterations", 
      "stopReason", "neval", "iterations.first", "iterations.second", 
      "elapsed"), parameters = c("a", "b", "c"), row.names = c(NA, 
      3L), run = "r1")

---

    Code
      .canonDput(as.parframe(fits, sortBy = "iterations"))
    Output
      structure(list(index = c(3L, 1L, 4L), value = c(1, 3, 2), converged = c(TRUE, 
      TRUE, FALSE), iterations = c(3L, 5L, 10L), stopReason = c("gradient", 
      "gradient", "iterlim"), neval = c(4L, 7L, NA), iterations.first = c(NA, 
      2L, NA), iterations.second = c(NA, 3L, NA), elapsed = c(NA, 0.5, 
      NA), a = c(2, 1, 5), b = c(3, 2, NA), c = c(NA, NA, 9)), class = c("parframe", 
      "data.frame"), metanames = c("index", "value", "converged", "iterations", 
      "stopReason", "neval", "iterations.first", "iterations.second", 
      "elapsed"), parameters = c("a", "b", "c"), row.names = c(NA, 
      3L), run = "r1")

---

    Code
      .canonDput(as.parframe(parlist(list(error = "x"))))
    Output
      structure(list(index = integer(0), value = numeric(0), converged = logical(0), 
          iterations = integer(0)), class = c("parframe", "data.frame"
      ), metanames = c("index", "value", "converged", "iterations"), parameters = character(0), row.names = integer(0))

---

    Code
      .canonDput(as.parframe(parlist()))
    Output
      structure(list(index = integer(0), value = numeric(0), converged = logical(0), 
          iterations = integer(0)), class = c("parframe", "data.frame"
      ), metanames = c("index", "value", "converged", "iterations"), parameters = character(0), row.names = integer(0))

---

    Code
      .canonDput(as.parframe(parlist(list(value = 1, argument = c(a = 1), iterations = 2))))
    Output
      structure(list(index = 1L, value = 1, converged = FALSE, iterations = 2L, 
          a = 1), class = c("parframe", "data.frame"), metanames = c("index", 
      "value", "converged", "iterations"), parameters = "a", row.names = 1L)

# plot of a parlist, initial against converged values and the paths

    Code
      .canonDput(as.data.frame(p$data))
    Output
      structure(list(parameter = c("a", "a", "b", "b"), x = c(0, 1, 
      0, 1), y = c(1, 2, 2, 3), run = c(1, 1, 1, 1)), class = "data.frame", row.names = c(NA, 
      4L))

---

    Code
      .canonDput(as.data.frame(q$data))
    Output
      structure(list(parameter = c("a", "a", "a", "b", "b", "b", "a", 
      "a", "b", "b"), iteration = c(1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 
      1L, 2L), path = c(0, 0.5, 1, 0, 1, 2, 1, 2, 1, 3), idx = c("1", 
      "1", "1", "1", "1", "1", "2", "2", "2", "2")), class = "data.frame", row.names = c(NA, 
      10L))

# plotPars and plotValues of a parframe

    Code
      .canonDput(attr(pp, "data"))
    Output
      structure(list(step = structure(c(1L, 1L, 2L, 1L, 1L, 2L, 1L, 
      1L, 2L), class = "factor", levels = c("1", "2")), value = c(1, 
      2, 3, 1, 2, 3, 1, 2, 3), name = structure(c(1L, 1L, 1L, 2L, 2L, 
      2L, 3L, 3L, 3L), class = "factor", levels = c("a", "b", "c")), 
          parvalue = c(2, 5, 1, 3, NA, 2, NA, 9, NA)), class = "data.frame", row.names = c(NA, 
      9L))

---

    Code
      .canonDput(attr(pv, "data"))
    Output
      structure(list(index = 1:3, value = c(1, 2, 3), converged = c(TRUE, 
      FALSE, TRUE), iterations = c(3L, 10L, 5L), stopReason = c("gradient", 
      "iterlim", "gradient"), neval = c(4L, NA, 7L), iterations.first = c(NA, 
      NA, 2L), iterations.second = c(NA, NA, 3L), elapsed = c(NA, NA, 
      0.5), a = c(2, 5, 1), b = c(3, NA, 2), c = c(NA, 9, NA), delta = c(0, 
      1, 2)), class = c("parframe", "data.frame"), metanames = c("index", 
      "value", "converged", "iterations", "stopReason", "neval", "iterations.first", 
      "iterations.second", "elapsed"), parameters = c("a", "b", "c"
      ), row.names = c(NA, 3L))

# the parlist overlay of a profile plot

    Code
      .canonDput(layer$data)
    Output
      structure(list(par = c(1, 2), name = c("a", "a"), delta = c(0.5, 
      3)), class = "data.frame", row.names = 1:2)

---

    Code
      .canonDput(dMod2:::.profilePoints(m, sets)$data)
    Output
      structure(list(par = c(1, 2), name = c("a", "a"), delta = c(0, 
      0)), class = "data.frame", row.names = 1:2)

# print of a parvec

    Code
      print(x)
    Output
      Parameter vector:
        a :       1
        b :  2.5e-07
        k :      -3
      
      Attributes:
        deriv  : 2 x 2 matrix
        deriv2 : 2 x 2 x 2 array
        fixed  : k

---

    Code
      print(as.parvec(c(a = 1)))
    Output
      Parameter vector:
        a :       1
      
      Attributes:
        deriv  : <none>
        deriv2 : <none>
        fixed  : <none>

---

    Code
      print(structure(numeric(0), class = c("parvec", "numeric")))
    Output
      Parameter vector:
      
      Attributes:
        deriv  : <none>
        deriv2 : <none>
        fixed  : <none>

# subsetting a parvec

    Code
      .canonDput(unclass(x[c("b", "c")]))
    Output
      structure(c(b = 2, c = 3), deriv = structure(c(0, 0, 2, 0), dim = c(2L, 
      2L), dimnames = list(c("b", "c"), c("u", "v"))))

---

    Code
      .canonDput(unclass(x["a", drop = TRUE]))
    Output
      structure(c(a = 1), deriv = structure(1, dim = c(1L, 1L), dimnames = list(
          "a", "u")))

# parframe() and is.parframe()

    Code
      .canonDput(parframe())
    Output
      structure(list(), class = c("parframe", "data.frame"), names = character(0), row.names = integer(0))

---

    Code
      .canonDput(parframe(df, metanames = "value"))
    Output
      structure(list(value = 1:2, a = 3:4), class = c("parframe", "data.frame"
      ), metanames = "value", parameters = "a", row.names = 1:2)

---

    Code
      .canonDput(parframe(cbind(value = 1, a = 2, b = 3), parameters = "a",
      metanames = "value", objAttributes = "b"))
    Output
      structure(list(value = 1, a = 2, b = 3), class = c("parframe", 
      "data.frame"), metanames = "value", obj.attributes = "b", parameters = "a", row.names = 1L)

---

    Code
      .canonDput(parframe(data.table::data.table(value = 1, a = 2)))
    Output
      structure(list(value = 1, a = 2), class = c("parframe", "data.frame"
      ), parameters = c("value", "a"), row.names = 1L)

# as.datalist of a list

    Code
      .canonDput(as.datalist(l))
    Output
      structure(list(C1 = structure(list(name = "A", time = 0, value = 1, 
          sigma = NA_real_, lloq = -Inf), class = "data.frame", row.names = 1L), 
          C2 = structure(list(name = "A", time = 1, value = 2, lloq = -Inf, 
              sigma = NA_real_), class = "data.frame", row.names = 1L)), class = c("datalist", 
      "list"), condition.grid = structure(list(condition = c("C1", 
      "C2")), class = "data.frame", row.names = c("C1", "C2")))

---

    Code
      .canonDput(as.datalist(unname(l), names = c("u", "v")))
    Output
      structure(list(u = structure(list(name = "A", time = 0, value = 1, 
          sigma = NA_real_, lloq = -Inf), class = "data.frame", row.names = 1L), 
          v = structure(list(name = "A", time = 1, value = 2, lloq = -Inf, 
              sigma = NA_real_), class = "data.frame", row.names = 1L)), class = c("datalist", 
      "list"), condition.grid = structure(list(condition = c("u", "v"
      )), class = "data.frame", row.names = c("u", "v")))

---

    Code
      .canonDput(as.datalist(list()))
    Output
      structure(list(), class = c("datalist", "list"), condition.grid = structure(list(), names = character(0), class = "data.frame", row.names = integer(0)))

# renaming a datalist renames its grid

    Code
      .canonDput(dl)
    Output
      structure(list(lo = structure(list(name = c("A", "B", "A", "B"
      ), time = c(0L, 0L, 1L, 1L), value = 1:4, sigma = c(0.5, 0.5, 
      0.5, 0.5), lloq = c(-Inf, -Inf, 3.5, -Inf)), class = "data.frame", row.names = c(NA, 
      4L)), hi = structure(list(name = c("A", "B", "A", "B"), time = c(0L, 
      0L, 1L, 1L), value = 5:8, sigma = c(0.5, 0.5, 0.5, 0.5), lloq = c(-Inf, 
      -Inf, 3.5, -Inf)), class = "data.frame", row.names = 5:8)), class = c("datalist", 
      "list"), condition.grid = structure(list(dose = c(1, 2)), row.names = c("lo", 
      "hi"), class = "data.frame"))

---

    Code
      .canonDput(attr(dl, "condition.grid"))
    Output
      structure(list(condition = "X"), class = "data.frame", row.names = "X")

# c() of datalists combines the grids

    Code
      .canonDput(c(d, e))
    Output
      structure(list(`1` = structure(list(name = c("A", "B", "A", "B"
      ), time = c(0L, 0L, 1L, 1L), value = 1:4, sigma = c(0.5, 0.5, 
      0.5, 0.5), lloq = c(-Inf, -Inf, 3.5, -Inf)), class = "data.frame", row.names = c(NA, 
      4L)), `2` = structure(list(name = c("A", "B", "A", "B"), time = c(0L, 
      0L, 1L, 1L), value = 5:8, sigma = c(0.5, 0.5, 0.5, 0.5), lloq = c(-Inf, 
      -Inf, 3.5, -Inf)), class = "data.frame", row.names = 5:8), C3 = structure(list(
          name = "A", time = 0, value = 9, sigma = 1, lloq = -Inf), class = "data.frame", row.names = 1L)), class = c("datalist", 
      "list"), condition.grid = structure(list(dose = c(1, 2, NA), 
          condition = c(NA, NA, "C3")), row.names = c("1", "2", "C3"
      ), class = "data.frame"))

---

    Code
      .canonDput(c(d[1], d[2]))
    Output
      structure(list(`1` = structure(list(name = c("A", "B", "A", "B"
      ), time = c(0L, 0L, 1L, 1L), value = 1:4, sigma = c(0.5, 0.5, 
      0.5, 0.5), lloq = c(-Inf, -Inf, 3.5, -Inf)), class = "data.frame", row.names = c(NA, 
      4L)), `2` = structure(list(name = c("A", "B", "A", "B"), time = c(0L, 
      0L, 1L, 1L), value = 5:8, sigma = c(0.5, 0.5, 0.5, 0.5), lloq = c(-Inf, 
      -Inf, 3.5, -Inf)), class = "data.frame", row.names = 5:8)), class = c("datalist", 
      "list"), condition.grid = structure(list(dose = c(1, 2)), row.names = c("1", 
      "2"), class = "data.frame"))

# plotData of a datalist

    Code
      .canonDput(attr(p, "data"))
    Output
      structure(list(condition = c("1", "1", "1", "1", "2", "2", "2", 
      "2"), name = c("A", "B", "A", "B", "A", "B", "A", "B"), time = c(0L, 
      0L, 1L, 1L, 0L, 0L, 1L, 1L), value = 1:8, sigma = c(0.5, 0.5, 
      0.5, 0.5, 0.5, 0.5, 0.5, 0.5), lloq = c(-Inf, -Inf, 3.5, -Inf, 
      -Inf, -Inf, 3.5, -Inf), dose = c(1, 1, 1, 1, 2, 2, 2, 2), bloq = c("no", 
      "no", "yes", "no", "no", "no", "no", "no")), class = "data.frame", row.names = c(NA, 
      8L))

