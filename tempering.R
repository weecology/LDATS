a work in progress toy example of tempering

i think the pseudo priors should be inverse temps?
basically how theyre set now
messing around, have the basic idea, i think
verify the math for the temps and such tho
library(mcmc)

set.seed(1)
xx <- c(rnorm(2e5, 0, 0.5), rnorm(1e5, 4, 0.7))
plot(density(xx))
plot(density(xx)$x, log(density(xx)$y), type = "l")

prixx <- c(rnorm(2e4, 0, 1), rnorm(1e4, 5, 2.5))
plot(density(prixx))
plot(density(prixx)$x, log(density(prixx)$y), type = "l")


temp_fun <- function(chain){
  c(1, 2, 5)[chain]
}

lik_fun <- function(params){
  out <- approxfun(density(xx))(params)
  out[is.na(out)] <- 1e-100
  out
}

prior_fun <- function(params){
  out <- approxfun(density(prixx))(params)
  out[is.na(out)] <- 1e-100
  out
}

pseudoprior_fun <- function(chain){
  c(1, 0.5, 0.2)[chain]
}

eval_fun <- function(inputs, lik_fun, prior_fun, temp_fun, pseudoprior_fun){
  chain <- inputs[1]
  params <- inputs[-1]
  temp <- temp_fun(chain)
  lik <- lik_fun(params)
  prior <- prior_fun(params)
  pseudoprior <- pseudoprior_fun(chain)
  1/temp * log(lik) + log(prior) + log(pseudoprior)
}


xv<-rep(NA,512)
yv1<-rep(NA,512)
yv2<-rep(NA,512)
yv3<-rep(NA,512)

for(i in 1:512){
  xv[i] <- density(xx)$x[i]
  yv1[i] <- eval_fun(c(1, xv[i]), lik_fun, prior_fun, temp_fun,
                     pseudoprior_fun)
  yv2[i] <- eval_fun(c(2, xv[i]), lik_fun, prior_fun, temp_fun,
                     pseudoprior_fun)
  yv3[i] <- eval_fun(c(3, xv[i]), lik_fun, prior_fun, temp_fun,
                     pseudoprior_fun)
}

yvs <- c(yv1, yv2, yv3)
plot(xv, yv1, type = "l", 
     ylim = c(min(yvs, na.rm = TRUE), max(yvs, na.rm = TRUE)))
points(xv, yv2, type = "l", lty = 2)
points(xv, yv3, type = "l", lty = 3)



neighbors <- matrix(FALSE, 3, 3)
neighbors[row(neighbors) == col(neighbors) + 1] <- TRUE
neighbors[row(neighbors) == col(neighbors) - 1] <- TRUE


#serial tempering
initial <- c(3, 0)
tout <- temper(eval_fun, initial = initial, neighbors = neighbors,
               nbatch = 1000, blen = 1, 
               temp_fun = temp_fun, lik_fun = lik_fun, prior_fun = prior_fun,
               pseudoprior_fun = pseudoprior_fun)

#parallel tempering
initial <- matrix(c(0,0,0), 3, 1)
ptout <- temper(eval_fun, initial = initial, neighbors = neighbors,
               nbatch = 1000, blen = 1, parallel = TRUE,
               temp_fun = temp_fun, lik_fun = lik_fun, prior_fun = prior_fun,
               pseudoprior_fun = pseudoprior_fun)




# notes on temper

so, in looking through the mcmc package's code, the main calculator function
logh takes 3 arguments: the objective function (returns log unnormalized
density), the state, the R environment at the time of the C call within 
the temper function

basically, the temper function creates a function that is an evaluation
of the objective function at the state variable, but it doesn't actually do 
the evaluation, it just declares the relationship. this way, the function
(via the relationship) can be used on any state variable value, most
especially the proposal values. 
it then creates an object that is the environment (including things passed 
into temper via ...). this allows both more global environment variables and 
generalized inputs via ... to be passed down to where the execution of the
function actually happens 

this is because mcmc uses R to call C to call R


