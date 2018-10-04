x <- cbind(1, dct)
eta <- matrix(rep(c(0,0), 3), 
         nrow = ncol(x), ncol = ncol(gamma))
eta_p <- eta[,-1] # eta_p is what we're optimizing 
y <- gamma 

eta_p2 <- matrix(c(1, 0.5, 1, 0.5), 2, 2)
k <- ncol(eta)
fn(eta_p2, x, y, k)
gn(eta_p2, x, y, k)

optim(eta_p, fn, x = x, y = y, k = k)
optim(eta_p, fn, gn, x = x, y = y, k = k, method = "BFGS")

#####
fn <- function(eta_p, x, y, k){

  M <- nrow(x)
  eta_p <- matrix(eta_p, nrow = ncol(x), ncol = k - 1)
  eta <- cbind(0, eta_p)
  predicted_prob <- matrix(0, nrow = nrow(x), ncol = ncol(eta))
  err <- matrix(0, nrow = nrow(x), ncol = ncol(eta))

  for(d in 1:M){
    exp_x_eta <- numeric(k)
    for(i in 1:k){
      x_eta <- as.matrix(x[d, ]) %*% as.matrix(eta[ , i])
      exp_x_eta[i] <- exp(x_eta)
    }
    sum_exp_x_eta <- sum(exp_x_eta)
    for(i in 1:k){
      predicted_prob[d, i] <- exp_x_eta[i]/sum_exp_x_eta
      err[d, i] <- -y[d, i] * log(predicted_prob[d, i])
    }
  }
  sum(err)
}

# trying to write the gradient function to check i get it, but the code
#isnt right yet.

gn <- function(eta_p, x, y, k){

  M <- nrow(x)
  eta_p <- matrix(eta_p, nrow = ncol(x), ncol = k - 1)
  eta <- cbind(0, eta_p)
  predicted_prob <- matrix(0, nrow = nrow(x), ncol = ncol(eta))
  y_minus_pred <- matrix(0, nrow = nrow(x), ncol = ncol(eta))

  for(d in 1:M){
    exp_x_eta <- numeric(k)
    for(i in 1:k){
      x_eta <- as.matrix(x[d, ]) %*% as.matrix(eta[ , i])
      exp_x_eta[i] <- exp(x_eta)
    }
    sum_exp_x_eta <- sum(exp_x_eta)
    for(i in 1:k){
      predicted_prob[d, i] <- exp_x_eta[i]/sum_exp_x_eta
      y_minus_pred[d, i] <- (y[d, i] - predicted_prob[d, i])
    }
  }
  grad_mat <-  t(x) %*% y_minus_pred[,-1]
  as.double(grad_mat)
}