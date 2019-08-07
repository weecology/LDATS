   tD <- c(1, 3, 4, 6)
   rho <- 3
   X <- matrix(c(1,1,1,2,1,3,1,4), 
          nrow = length(tD), ncol = 2, byrow = TRUE)
   Eta <- matrix(c(0.5, 1.2, 0.3, 1.1, 0.9, 0.1, 0.5, 0.5), 
           nrow = ncol(X), ncol = 2, byrow = TRUE)
   err = 0
   seed = 1

I think if we set up in2 better it will mean X can be simple and not sparse

ACK
:(C * S) should have been :(C * s)
this didnt actually influence anything since how i had it set up

# consider adding a check for Eta and X conformity to the top end of the
sim functions