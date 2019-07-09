context("Check simulate functions")

test_that("check sim_LDA_data", {
  N <- c(10, 22, 15, 31)
  alpha <- 1.2
  Beta <- matrix(c(0.1, 0.1, 0.8, 0.2, 0.6, 0.2), 2, 3, byrow = TRUE)
  expect_is(sim_LDA_data(N, alpha, Beta), "matrix")
  expect_equal(dim(sim_LDA_data(N, alpha, Beta)), c(4,3))
  expect_error(sim_LDA_data("ok", alpha, Beta))
  expect_error(sim_LDA_data(N + 1.1, alpha, Beta))
  expect_error(sim_LDA_data(N, rep(alpha, 2), Beta))
  expect_error(sim_LDA_data(N, "ok", Beta))
  expect_error(sim_LDA_data(N, -1, Beta))
  expect_error(sim_LDA_data(N, alpha, "ok"))
  expect_error(sim_LDA_data(N, alpha, alpha))
  Beta2 <- matrix(c(0.1, 0.2, 0.8, 0.2, 0.6, 0.2), 2, 3, byrow = TRUE)
  Beta3 <- matrix(c(0.3, -0.1, 0.8, 0.2, 0.6, 0.2), 2, 3, byrow = TRUE)
  expect_error(sim_LDA_data(N, alpha, Beta2))
  expect_error(sim_LDA_data(N, alpha, Beta3))
  Beta4 <- matrix(c(0.1, 0.1, 0.1, 0.1, 0.4, 0.2), 1, 6, byrow = TRUE)
  expect_is(sim_LDA_data(N, alpha, Beta4), "matrix")
  expect_equal(dim(sim_LDA_data(N, alpha, Beta4)), c(4,6))
})

test_that("check sim_TS_data", {
   tD <- c(1, 3, 4, 6)
   rho <- 3
   X <- matrix(c(1,1,0,0,1,2,0,0,0,0,1,3,0,0,1,4), 
          nrow = length(tD), ncol = 4, byrow = TRUE)
   Eta <- matrix(c(0.5, 1.2, 0.3, 1.1, 0.9, 0.1, 0.5, 0.5), 
           nrow = ncol(X), ncol = 2, byrow = TRUE)
   expect_is(sim_TS_data(X, Eta, rho, tD), "matrix")
   expect_equal(dim(sim_TS_data(X, Eta, rho, tD)), c(4,2))
   expect_equal(round(sim_TS_data(X, Eta, rho, tD), 2)[1,1], 0.18)

   expect_error(sim_TS_data(X, Eta, rho, "ok"))
   expect_error(sim_TS_data(X, Eta, rho, matrix(1, 2, 2)))
   expect_error(sim_TS_data(X, Eta, rho, tD + 0.1))
   expect_error(sim_TS_data("ok", Eta, rho, tD))
   expect_error(sim_TS_data(array(1, dim=c(2,2,2)), Eta, rho, tD))
   expect_error(sim_TS_data(X, "ok", rho, tD))
   expect_error(sim_TS_data(X, array(1, dim=c(2,2,2)), rho, tD))
   expect_error(sim_TS_data(X, Eta, "ok", tD))
   expect_error(sim_TS_data(X, Eta, 1.1, tD))
   expect_error(sim_TS_data(X, Eta, matrix(1, 2, 1), tD))
   expect_is(sim_TS_data(X, Eta, NULL, tD), "matrix")

   expect_error(sim_TS_data(X, Eta, rho, tD, -1))
   expect_error(sim_TS_data(X, Eta, rho, tD, c(1,2)))
   expect_error(sim_TS_data(X, Eta, rho, tD, "ok"))
})