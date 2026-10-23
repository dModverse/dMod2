# print.objlist prints value, truncated gradient and Hessian, attributes

    Code
      print(o)
    Output
      value
      ==================
      3.5
      gradient[1:20] (full length = 25)
      ==================
       p1  p2  p3  p4  p5  p6  p7  p8  p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20 
      0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 
      
      hessian[1:6,1:6]
      ==================
         p1 p2 p3 p4 p5 p6
      p1  1  0  0  0  0  0
      p2  0  1  0  0  0  0
      p3  0  0  1  0  0  0
      p4  0  0  0  1  0  0
      p5  0  0  0  0  1  0
      p6  0  0  0  0  0  1
      
      
      attributes
      ==================
       chi2  2.25
       data  2.25

---

    Code
      print(o, n1 = 3, n2 = 2)
    Output
      value
      ==================
      3.5
      gradient[1:3] (full length = 25)
      ==================
       p1  p2  p3 
      0.1 0.2 0.3 
      
      hessian[1:2,1:2]
      ==================
         p1 p2
      p1  1  0
      p2  0  1
      
      
      attributes
      ==================
       chi2  2.25
       data  2.25

---

    Code
      print(objlist(1, c(a = 1), NULL))
    Output
      value
      ==================
      1
      gradient[1:1] (full length = 1)
      ==================
      a 
      1 
      
      
      attributes
      ==================

---

    Code
      print(structure(list(value = 4), class = "objlist"))
    Output
      value
      ==================
      4
      
      attributes
      ==================

