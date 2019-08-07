   tD <- c(1, 3, 4, 6)
   rho <- 3
   X <- matrix(c(1,1,1,2,1,3,1,4), 
          nrow = length(tD), ncol = 2, byrow = TRUE)
   Eta <- matrix(c(0.5, 1.2, 0.3, 1.1, 0.9, 0.1, 0.5, 0.5), 
           nrow = ncol(X), ncol = 2, byrow = TRUE)
   err = 0
   seed = 1


# consider adding a check for Eta and X conformity to the top end of the
sim functions

devtools::document()
devtools::test()
devtools::check()
devtools::check(run_dont_test=T)