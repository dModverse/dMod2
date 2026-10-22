# print.eqnvec lays out sorted, indexed, wrapped equations

    Code
      print(f)
    Output
      Idx     Inner <- Outer
        1         A <- -k1*A+k2*B
        2 longstate <- a+b
        3      zeta <- k1*A-k2*B

---

    Code
      print(long, width = 60)
    Output
      Idx Inner <- Outer
        1     x <- k1*x1+k2*x2+k3*x3+k4*x4+k5*x5+k6*x6+k7*x7+k8*x8+k9*x9+k10*x10+k11*x11+k12*x12+k13*x13+k14*x14+k15*x15+k16*x16+k17*x17+k18*x18+k19*x19+k20*x20+k21*x21+k22*x22+k23*x23+k24*x24+k25*x25+k26*x26+k27*x27+k28*x28+k29*x29+k30*x30+k31*x31+k32*x32+k33*x33+k34*x34+k35*x35+k36*x36+k37*x37+k38*x38+k39*x39+k40*x40
        2     y <- 1

---

    Code
      print(long, width = 5)
    Output
      Idx Inner <- Outer
        1     x <- k1*x1+k2*x2+k3*x3+k4*x4+k5*x5+k6*x6+k7*x7+k8*x8+k9*x9+k10*x10+k11*x11+k12*x12+k13*x13+k14*x14+k15*x15+k16*x16+k17*x17+k18*x18+k19*x19+k20*x20+k21*x21+k22*x22+k23*x23+k24*x24+k25*x25+k26*x26+k27*x27+k28*x28+k29*x29+k30*x30+k31*x31+k32*x32+k33*x33+k34*x34+k35*x35+k36*x36+k37*x37+k38*x38+k39*x39+k40*x40
        2     y <- 1

---

    Code
      print(as.eqnvec(eqnNet()))
    Output
      Idx Inner <- Outer
        1     A <- -1*(act_A*A*stimulus)+1*(deact_A*pA)
        2     B <- -1*(act_B*B*pA_pA)+1*(deact_B*pB)
        3    pA <- 1*(act_A*A*stimulus)-1*(deact_A*pA)-2*(form_complex_pA*pA^2)
        4 pA_pA <- 1*(form_complex_pA*pA^2)
        5    pB <- 1*(act_B*B*pA_pA)-1*(deact_B*pB)

---

    Code
      print(ws, width = 50)
    Output
      Idx Inner <- Outer
        1     a <- k*x+ k*x+ k*x+ k*x+ k*x+ k*x+ k*x+
                   k*x+ k*x+ k*x+ k*x+ k*x+
                   k*x+ k*x+ k*x+ k*x+ k*x+
                   k*x+ k*x+ k*x+ k*x+ k*x+
                   k*x+ k*x+ k*x+ k*x+ k*x+
                   k*x+ k*x+ k*x
        2     b <- a+ b
        3     c <- u -v

---

    Code
      print(ws, width = 0)
    Output
      Idx Inner <- Outer
        1     a <- k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x
        2     b <- a+
                   b
        3     c <- u
                   -v

---

    Code
      print(eqnvec(xy = paste(rep("k*x", 30), collapse = "+\n")), width = 30)
    Output
      Idx Inner <- Outer
        1    xy <- k*x+ k*x+ k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x+
                   k*x

---

    Code
      print(f, pander = TRUE)
    Output
      
      -----------------------------------
          &nbsp;      &nbsp;             
      --------------- -------------------
         **zeta**     k1 * A - k2 * B    
      
           **A**      - k1 * A + k2 * B  
      
       **longstate**  a + b              
      -----------------------------------
      

# print.eqnlist shows totals, compartments and reactions

    Code
      print(eq)
    Output
      Conserved quantities:
        total_1 = 0.5*A+0.5*pA+1*pA_pA
        total_2 = 1*B+1*pB
      
        Check Educt -> Product                   Rate             Description
      1           A ->      pA   act_A * A * stimulus    Phosphorylation of A
      2          pA ->       A           deact_A * pA Dephosphorylation of pA
      3        2*pA ->   pA_pA form_complex_pA * pA^2 Complex formation of pA
      4           B ->      pB      act_B * B * pA_pA    Phosphorylation of B
      5          pB ->       B           deact_B * pB Dephosphorylation of pB

---

    Code
      print(eq2)
    Output
      Conserved quantities:
        total_1 = 1*A+1*B
      
      Compartments:
        cyt (V=V1)
        nuc (V=V2)
      States by compartment:
        cyt: A
        nuc: B
      
        Check Educt -> Product Rate Description
      1           A ->       B  k*A           1

---

    Code
      print(eqnlist())
    Output
      
      NULL

# conservedQuantities returns the kernel as expressions

    Code
      print(cq)
    Output
        Conserved quantities: 2
      1    0.5*A+0.5*pA+1*pA_pA
      2                1*B+1*pB

