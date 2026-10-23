# mstrust prints the multi-start summary with stats = TRUE

    Code
      fits <- mstrust(quadObj, c(a = 0, b = 0), fits = 6, sd = 2, retry = FALSE,
      stats = TRUE)
    Output
      Multi start summary
      Outcome     : Occurrence
      Error       : 3
      Unconverged : 0
      Converged   : 3
                 -----------
      Total       : 6
      
      Termination reason
      gradient    : 3

---

    Code
      invisible(mstrust(quadObj, c(a = 0, b = 0), fits = 2, sd = 0.5, iterlim = 1,
      stats = TRUE))
    Condition
      Warning in `trust_impl()`:
      Maximum number of iterations exceeded. Fit is not converged.
      Warning in `trust_impl()`:
      Maximum number of iterations exceeded. Fit is not converged.
    Output
      Multi start summary
      Outcome     : Occurrence
      Error       : 0
      Unconverged : 2
      Converged   : 0
                 -----------
      Total       : 2
      
      Termination reason
      iterlim     : 2

# mstrust writes the run folder, log and parameter list with output = TRUE

    Code
      writeLines(readLines(log))
    Output
      Parameter assignment information
      mstrust     : name, center, fits, cores, optmethod, samplefun, resultPath, stats, output, cautiousMode, retry, nTries
      trust       : objfun, rinit, rmax
      rnorm       : n, sd
      
      -------
      Fit 1 completed
      --> iterations : 4
      -->  converged : TRUE
      --> obj. value : 0
      -------
      -------
      Fit 2 failed with error
      --> Error : parinit not feasible: objfun failed: outside
      -------
      -------
      Fit 3 failed with error
      --> Error : parinit not feasible: objfun failed: outside
      -------
      -------
      Fit 4 completed
      --> iterations : 6
      -->  converged : TRUE
      --> obj. value : 0
      -------
      Multi start summary
      Outcome     : Occurrence
      Error       : 2
      Unconverged : 0
      Converged   : 2
                 -----------
      Total       : 4
      
      Termination reason
      gradient    : 2
      

# reduceReplicates gives mean, standard error and the stable columns

    Code
      r <- reduceReplicates(d)
    Message
      Dropped columns: batch

---

    Code
      pinFrame(r)
    Output
      list(columns = list(time = c(0, 1, 2, 0, 1, 2), value = c(1.03333333333333, 
      2.2, 3, 5, 6.5, 8), sigma = c(0.0881917103688197, 0.2, NA, 0.288675134594813, 
      0.5, NA), n = c(3L, 2L, 1L, 3L, 2L, 1L), name = c("A", "A", "A", 
      "B", "B", "B"), condition = c("ctrl", "ctrl", "ctrl", "ctrl", 
      "ctrl", "ctrl"), dose = c(1, 1, 1, 2, 2, 2)), rows = c(1L, 4L, 
      6L, 7L, 10L, 12L), class = "data.frame")

---

    Code
      pinFrame(suppressMessages(reduceReplicates(d, select = c("condition", "dose"))))
    Output
      list(columns = list(time = c(0, 1, 2, 0, 1, 2), value = c(1.03333333333333, 
      2.2, 3, 5, 6.5, 8), sigma = c(0.0881917103688197, 0.2, NA, 0.288675134594813, 
      0.5, NA), n = c(3L, 2L, 1L, 3L, 2L, 1L), name = c("A", "A", "A", 
      "B", "B", "B"), condition = c("ctrl_1", "ctrl_1", "ctrl_1", "ctrl_2", 
      "ctrl_2", "ctrl_2")), rows = c(1L, 4L, 6L, 7L, 10L, 12L), class = "data.frame")

---

    Code
      pinFrame(suppressMessages(reduceReplicates(d, keep = c("batch", "nothere"))))
    Output
      list(columns = list(time = c(0, 1, 2, 0, 1, 2), value = c(1.03333333333333, 
      2.2, 3, 5, 6.5, 8), sigma = c(0.0881917103688197, 0.2, NA, 0.288675134594813, 
      0.5, NA), n = c(3L, 2L, 1L, 3L, 2L, 1L), name = c("A", "A", "A", 
      "B", "B", "B"), condition = c("ctrl", "ctrl", "ctrl", "ctrl", 
      "ctrl", "ctrl"), dose = c(1, 1, 1, 2, 2, 2), batch = c(1, 1, 
      1, 1, 1, 1)), rows = c(1L, 4L, 6L, 7L, 10L, 12L), class = "data.frame")

---

    Code
      pinFrame(suppressMessages(reduceReplicates(d, datatrans = "log(x)")))
    Output
      list(columns = list(time = c(0, 1, 2, 0, 1, 2), value = c(0.0256536803787094, 
      0.784307958956923, 1.09861228866811, 1.60608780048293, 1.86883480914168, 
      2.07944154167984), sigma = c(0.0840313972874974, 0.0911607783969773, 
      NA, 0.0579528527968493, 0.0770753399136291, NA), n = c(3L, 2L, 
      1L, 3L, 2L, 1L), name = c("A", "A", "A", "B", "B", "B"), condition = c("ctrl", 
      "ctrl", "ctrl", "ctrl", "ctrl", "ctrl"), dose = c(1, 1, 1, 2, 
      2, 2)), rows = c(1L, 4L, 6L, 7L, 10L, 12L), class = "data.frame")

---

    Code
      pinFrame(suppressMessages(reduceReplicates(d, weighted = TRUE)))
    Output
      list(columns = list(time = c(0, 1, 2, 0, 1, 2), value = c(0.977777777777778, 
      2.12307692307692, 3, 5.07894736842105, 6.69230769230769, 8), 
          sigma = c(0.0687184270936277, 0.2, NA, 0.298354070719928, 
          0.5, NA), n = c(3L, 2L, 1L, 3L, 2L, 1L), name = c("A", "A", 
          "A", "B", "B", "B"), condition = c("ctrl", "ctrl", "ctrl", 
          "ctrl", "ctrl", "ctrl"), dose = c(1, 1, 1, 2, 2, 2)), rows = c(1L, 
      4L, 6L, 7L, 10L, 12L), class = "data.frame")

---

    Code
      pinFrame(suppressMessages(reduceReplicates(d[, -7])))
    Output
      list(columns = list(time = c(0, 1, 2, 0, 1, 2), value = c(1.03333333333333, 
      2.2, 3, 5, 6.5, 8), sigma = c(0.0881917103688197, 0.2, NA, 0.288675134594813, 
      0.5, NA), n = c(3L, 2L, 1L, 3L, 2L, 1L), name = c("A", "A", "A", 
      "B", "B", "B"), condition = c("ctrl", "ctrl", "ctrl", "ctrl", 
      "ctrl", "ctrl"), dose = c(1, 1, 1, 2, 2, 2)), rows = c(1L, 4L, 
      6L, 7L, 10L, 12L), class = "data.frame")

# fitErrorModel fits per group and returns the standard error

    Code
      print(roundNumeric(out))
    Output
         time    value    sigma n name condition
      1     1   0.0538  0.05657 5    y      ctrl
      6     2   0.9984  0.08318 5    y      ctrl
      11    3   8.3500  0.51390 5    y      ctrl
      16    4  97.7200  5.97800 5    y      ctrl
      21    1   0.1846  0.01985 5    z      ctrl
      26    2   1.8230  0.13200 5    z      ctrl
      31    3  16.9500  1.22000 5    z      ctrl
      36    4 207.1000 14.91000 5    z      ctrl

---

    Code
      print(roundNumeric(full))
    Output
         time    value    sigma n name condition     s0   srel cbLower68 cbUpper68
      1     1   0.0538  0.05657 5    y      ctrl -4.139 -3.979 5.661e-03 2.640e-02
      6     2   0.9984  0.08318 5    y      ctrl -4.139 -3.979 1.224e-02 5.710e-02
      11    3   8.3500  0.51390 5    y      ctrl -4.139 -3.979 4.673e-01 2.179e+00
      16    4  97.7200  5.97800 5    y      ctrl -4.139 -3.979 6.324e+01 2.949e+02
      21    1   0.1846  0.01985 5    z      ctrl -6.824 -3.653 6.972e-04 3.252e-03
      26    2   1.8230  0.13200 5    z      ctrl -6.824 -3.653 3.084e-02 1.439e-01
      31    3  16.9500  1.22000 5    z      ctrl -6.824 -3.653 2.635e+00 1.229e+01
      36    4 207.1000 14.91000 5    z      ctrl -6.824 -3.653 3.932e+02 1.834e+03
         cbLower95 cbUpper95 condidnt  sigmaLS
      1  1.830e-03 4.556e-02        y  0.05032
      6  3.958e-03 9.853e-02        y  0.10030
      11 1.511e-01 3.761e+00        y  0.48350
      16 2.044e+01 5.089e+02        y  5.57400
      21 2.254e-04 5.611e-03        z  0.02000
      26 9.970e-03 2.482e-01        z  0.07783
      31 8.516e-01 2.120e+01        z  1.52400
      36 1.271e+02 3.164e+03        z 15.48000

# steadyStates refuses named givenCQs and returns lhs = rhs as a named vector

    Code
      ss
    Output
            k1        A        B 
      "B*k2/A"      "A"      "B" 

