# constraintL2(Omega = ) -- the multivariate-normal path.
# Needs omega() from the NLME layer, so these live with it.

test_that("MVN value matches direct computation for diagonal Omega", {
  set.seed(2)
  subjects <- paste0("s", 1:5)
  om <- omega(eta = c("eta_a", "eta_b"), subjects = subjects)
  obj <- constraintL2(mu = 0, Omega = om)

  eta_vals <- runif(10, -0.5, 0.5)
  names(eta_vals) <- as.vector(om$subjectEtas)
  chol_vals <- c(omega_a_a = log(0.3), omega_b_b = log(0.4))
  pars <- c(eta_vals, chol_vals)

  L <- om$buildL(chol_vals)
  Omega_mat <- L %*% t(L)
  Omega_inv <- solve(Omega_mat)
  N <- length(subjects); K <- 2L
  eta_mat <- matrix(eta_vals, nrow = N, ncol = K)
  expected_quad <- sum(vapply(seq_len(N), function(i)
    drop(eta_mat[i, ] %*% Omega_inv %*% eta_mat[i, ]), 0.0))
  expected_logdet <- N * log(det(Omega_mat))
  expected_val <- expected_quad + expected_logdet

  expect_equal(obj(pars)$value, expected_val, tolerance = 1e-10)
})

test_that("MVN gradient (diagonal Omega) follows analytic closed form", {
  set.seed(3)

  subjects <- paste0("s", 1:4)
  om <- omega(eta = c("eta_a", "eta_b"), subjects = subjects)
  obj <- constraintL2(mu = 0, Omega = om)

  eta_vals <- rnorm(8, sd = 0.3)
  names(eta_vals) <- as.vector(om$subjectEtas)
  chol_vals <- c(omega_a_a = log(0.5), omega_b_b = log(0.4))
  pars <- c(eta_vals, chol_vals)

  # Diag Omega: omega_k = exp(chol_k).
  # value         = sum_i sum_k (eta_ik / omega_k)^2 + N * 2 * sum_k chol_k
  # dvalue/deta   = 2 * eta_ik / omega_k^2
  # dvalue/dchol  = -2 * sum_i (eta_ik / omega_k)^2 + 2 * N
  omega_vec <- exp(chol_vals)
  eta_mat <- matrix(eta_vals, nrow = length(subjects), ncol = 2,
                    dimnames = dimnames(om$subjectEtas))
  N <- nrow(eta_mat)
  g_eta <- 2 * eta_mat / matrix(omega_vec^2, N, 2, byrow = TRUE)
  g_chol <- c(omega_a_a = -2 * sum((eta_mat[, 1] / omega_vec[1])^2) + 2 * N,
              omega_b_b = -2 * sum((eta_mat[, 2] / omega_vec[2])^2) + 2 * N)
  g_ref <- c(as.vector(g_eta), g_chol)
  names(g_ref) <- c(as.vector(om$subjectEtas), names(chol_vals))

  analytic <- obj(pars)$gradient
  expect_equal(unname(analytic[names(g_ref)]), unname(g_ref),
               tolerance = 1e-10)
})

test_that("MVN gradient eta-block (full Omega) equals 2 * Omega^-1 * eta_i per subject", {
  set.seed(4)

  subjects <- paste0("s", 1:6)
  om <- omega(eta = c("eta_Cl", "eta_V", "eta_Ka"),
              structure = "full", subjects = subjects)
  obj <- constraintL2(mu = 0, Omega = om)

  eta_vals <- rnorm(18, sd = 0.4)
  names(eta_vals) <- as.vector(om$subjectEtas)
  chol_vals <- c(omega_Cl_Cl = log(0.3), omega_V_V = log(0.2), omega_Ka_Ka = log(0.5),
                 omega_V_Cl  = 0.05,     omega_Ka_Cl = -0.02,  omega_Ka_V  = 0.07)
  pars <- c(eta_vals, chol_vals)

  L <- om$buildL(chol_vals)
  Omega_mat <- L %*% t(L)
  Omega_inv <- solve(Omega_mat)
  eta_mat <- matrix(eta_vals, nrow = length(subjects), ncol = 3,
                    dimnames = dimnames(om$subjectEtas))
  g_eta_ref <- 2 * eta_mat %*% Omega_inv

  analytic <- obj(pars)$gradient
  eta_names <- as.vector(om$subjectEtas)
  expect_equal(unname(analytic[eta_names]),
               as.vector(g_eta_ref), tolerance = 1e-10)
})

test_that("MVN gradient eta-block (selective correlation) equals 2 * Omega^-1 * eta_i per subject", {
  set.seed(5)

  subjects <- paste0("s", 1:5)
  om <- omega(eta = c("eta_Cl", "eta_V", "eta_Ka"),
              correlate = list(c("eta_Cl", "eta_V")),
              subjects = subjects)
  obj <- constraintL2(mu = 0, Omega = om)

  eta_vals <- rnorm(15, sd = 0.3)
  names(eta_vals) <- as.vector(om$subjectEtas)
  chol_vals <- c(omega_Cl_Cl = log(0.3), omega_V_V = log(0.2), omega_Ka_Ka = log(0.4),
                 omega_V_Cl  = 0.04)
  pars <- c(eta_vals, chol_vals)

  L <- om$buildL(chol_vals)
  Omega_mat <- L %*% t(L)
  Omega_inv <- solve(Omega_mat)
  eta_mat <- matrix(eta_vals, nrow = length(subjects), ncol = 3,
                    dimnames = dimnames(om$subjectEtas))
  g_eta_ref <- 2 * eta_mat %*% Omega_inv

  analytic <- obj(pars)$gradient
  eta_names <- as.vector(om$subjectEtas)
  expect_equal(unname(analytic[eta_names]),
               as.vector(g_eta_ref), tolerance = 1e-10)
})

test_that("MVN Hessian eta-block at eta = 0 equals 2 * Omega^-1 block-diagonal per subject", {
  # Hessian of the quadratic term w.r.t. eta_i is 2 * Omega^-1 exactly. At
  # eta = 0 this matches the full Hessian because log|Omega| does not depend
  # on eta.
  set.seed(6)

  subjects <- paste0("s", 1:3)
  om <- omega(eta = c("eta_a", "eta_b"), structure = "full", subjects = subjects)
  obj <- constraintL2(mu = 0, Omega = om)

  eta_vals <- rep(0, 6)
  names(eta_vals) <- as.vector(om$subjectEtas)
  chol_vals <- c(omega_a_a = log(0.3), omega_b_b = log(0.4),
                 omega_b_a = 0.05)
  pars <- c(eta_vals, chol_vals)

  L <- om$buildL(chol_vals)
  Omega_inv <- solve(L %*% t(L))
  N <- length(subjects); K <- 2L
  eta_names <- as.vector(om$subjectEtas)
  H_ref <- matrix(0, length(eta_names), length(eta_names),
                  dimnames = list(eta_names, eta_names))
  for (i in seq_len(N)) {
    idx <- om$subjectEtas[i, ]
    H_ref[idx, idx] <- 2 * Omega_inv
  }
  analytic_h <- obj(pars)$hessian
  expect_equal(unname(analytic_h[eta_names, eta_names]),
               unname(H_ref), tolerance = 1e-10)
})

test_that("MVN value uses mu correctly when mu != 0", {
  set.seed(7)
  subjects <- paste0("s", 1:3)
  om <- omega(eta = c("eta_a", "eta_b"), subjects = subjects)
  obj <- constraintL2(mu = c(eta_a = 0.1, eta_b = -0.05), Omega = om)

  eta_vals <- c(eta_a_s1 = 0.1, eta_b_s1 = -0.05,
                eta_a_s2 = 0.1, eta_b_s2 = -0.05,
                eta_a_s3 = 0.1, eta_b_s3 = -0.05)
  chol_vals <- c(omega_a_a = log(0.3), omega_b_b = log(0.4))
  pars <- c(eta_vals, chol_vals)

  result <- obj(pars)
  expected_val <- 3 * log(det(om$buildL(chol_vals) %*% t(om$buildL(chol_vals))))
  expect_equal(result$value, expected_val, tolerance = 1e-10)
})

test_that("MVN errors helpfully when subject expansion is missing", {
  om <- omega(eta = c("eta_a", "eta_b"))   # no subjects
  expect_error(constraintL2(mu = 0, Omega = om), "subject expansion")
})

test_that("MVN summable with normL2-style objfn via +.objfn", {
  subjects <- paste0("s", 1:3)
  om <- omega(eta = c("eta_a"), subjects = subjects)

  obj_theta <- constraintL2(mu = c(theta_a = 0))
  obj_mvn   <- constraintL2(mu = 0, Omega = om)
  obj_sum   <- obj_theta + obj_mvn

  pars <- c(theta_a = 1.5,
            eta_a_s1 = 0.1, eta_a_s2 = -0.2, eta_a_s3 = 0.3,
            omega_a_a = log(0.4))

  vsum  <- obj_sum(pars)$value
  vsep  <- obj_theta(pars)$value + obj_mvn(pars)$value
  expect_equal(vsum, vsep, tolerance = 1e-10)
})

