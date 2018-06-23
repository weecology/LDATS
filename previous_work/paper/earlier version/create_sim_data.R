

#' create simulated data for demonstrating LDA analysis
#'    2 topics
#'     
#' @param nspecies = number of species in all topic groups
#' @param tsteps = number of [monthly] time steps
#' 
#' @return 
#'    beta = matrix of species composition of the groups
#'    gamma = matrix of topic composition over time
#'            3 simulations of gamma: uniform, slow transition, and fast transition
#' @export 

create_sim_data_2topic = function(nspecies=24,tsteps=400) {

  topics = 2
  
  # beta: species composition of topics -- uniform distribution, nonoverlapping species composition
  beta = matrix(rep(0,topics*nspecies),nrow=topics,ncol=nspecies)
  beta[1,] = c(rep(1/(nspecies/2),nspecies/2),rep(0,nspecies/2))
  beta[2,] = c(rep(0,nspecies/2),rep(1/(nspecies/2),nspecies/2))
  
  # gamma for a constant topic prevalence through time: topic1 at 90% and topic2 at 10%
  gamma_constant = matrix(rep(0,tsteps*topics),nrow=tsteps,ncol=topics)
  gamma_constant[,1] = rep(1,tsteps)
  gamma_constant[,2] = rep(0,tsteps)
  
  # gamma for a fast transition from topic1 to topic2 (one year/12 time steps)
  gamma_fast = matrix(rep(0,tsteps*topics),nrow=tsteps,ncol=topics)
  # proportions are constant for first 200 time steps
  gamma_fast[1:200,1] = rep(1)
  gamma_fast[1:200,2] = rep(0)
  # fast transition from tstep 201-212
  gamma_fast[201:212,1] = seq(12)*(-1/12)+1
  gamma_fast[201:212,2] = seq(12)*(1/12)+0
  # proportions are constant for rest of time series
  gamma_fast[213:tsteps,1] = rep(0)
  gamma_fast[213:tsteps,2] = rep(1) 
  
  # gamma for a slow transition from topic1 to topic2
  gamma_slow = matrix(rep(0,tsteps*topics),nrow=tsteps,ncol=topics)
  # brief period of constant values at beginning and end of series
  gamma_slow[1:50,1] = rep(1)
  gamma_slow[1:50,2] = rep(0)
  gamma_slow[351:400,1] = rep(0)
  gamma_slow[351:400,2] = rep(1)
  gamma_slow[51:350,1] = seq(300)*(-1/(tsteps-100))+1
  gamma_slow[51:350,2] = seq(300)*(1/(tsteps-100))+0
  
  return(list(beta,gamma_constant,gamma_fast,gamma_slow))
}

#' variation: create simulated data for demonstrating LDA analysis
#' 2 topics, nonuniform distribution of species in two community-types
#'     
#' @param tsteps = number of [monthly] time steps
#' 
#' @return 
#'    beta = matrix of species composition of the groups
#'    gamma = matrix of topic composition over time
#'            3 simulations of gamma: uniform, slow transition, and fast transition
#' @export 

create_sim_data_2topic_nonuniform = function(tsteps=400) {
  
  topics = 2
  nspecies = 12
  
  # beta: species composition of topics
  # I calculated this distribution by taking the average of each Portal sampling sp distribution (periods 1:436)
  distribution = c(27,13,7, 5, 3, 2, 1, 1, 1, 0, 0, 0)
  # simple permutation of the first distribution
  distribution2 = c(3,1, 0, 1, 0, 13,2, 0, 1,27, 5, 7)
  
  beta = matrix(rep(0,topics*nspecies),nrow=topics,ncol=nspecies)
  beta[1,] = distribution/sum(distribution)
  beta[2,] = distribution2/sum(distribution2)
  
  # gamma for a constant topic prevalence through time: topic1 at 90% and topic2 at 10%
  gamma_constant = matrix(rep(0,tsteps*topics),nrow=tsteps,ncol=topics)
  gamma_constant[,1] = rep(1,tsteps)
  gamma_constant[,2] = rep(0,tsteps)
  
  # gamma for a fast transition from topic1 to topic2 (one year/12 time steps)
  gamma_fast = matrix(rep(0,tsteps*topics),nrow=tsteps,ncol=topics)
  # proportions are constant for first 200 time steps
  gamma_fast[1:200,1] = rep(1)
  gamma_fast[1:200,2] = rep(0)
  # fast transition from tstep 201-212
  gamma_fast[201:212,1] = seq(12)*(-1/12)+1
  gamma_fast[201:212,2] = seq(12)*(1/12)+0
  # proportions are constant for rest of time series
  gamma_fast[213:tsteps,1] = rep(0)
  gamma_fast[213:tsteps,2] = rep(1) 
  
  # gamma for a slow transition from topic1 to topic2
  gamma_slow = matrix(rep(0,tsteps*topics),nrow=tsteps,ncol=topics)
  # brief period of constant values at beginning and end of series
  gamma_slow[1:50,1] = rep(1)
  gamma_slow[1:50,2] = rep(0)
  gamma_slow[351:400,1] = rep(0)
  gamma_slow[351:400,2] = rep(1)
  gamma_slow[51:350,1] = seq(300)*(-1/(tsteps-100))+1
  gamma_slow[51:350,2] = seq(300)*(1/(tsteps-100))+0
  
  return(list(beta,gamma_constant,gamma_fast,gamma_slow))
}

#' create simulated data for demonstrating LDA analysis
#'    3 topics
#'     
#' @param nspecies = number of species in all topic groups
#' @param tsteps = number of [monthly] time steps
#' 
#' @export 

create_sim_data_3topic = function(nspecies=24,tsteps=400) {

  topics = 3
  
  beta = matrix(rep(0,topics*nspecies),nrow=topics,ncol=nspecies)
  evencomp = 1/(nspecies/3)
  beta[1,] = c(rep(evencomp,nspecies/3),rep(0,nspecies/3),rep(0,nspecies/3))
  beta[2,] = c(rep(0,nspecies/3),rep(evencomp,nspecies/3),rep(0,nspecies/3))
  beta[3,] = c(rep(0,nspecies/3),rep(0,nspecies/3),rep(evencomp,nspecies/3))
  
  # gamma for a constant topic prevalence through time
  gamma_constant = matrix(rep(0,tsteps*topics),nrow=tsteps,ncol=topics)
  gamma_constant[,1] = rep(.7,tsteps)
  gamma_constant[,2] = rep(.2,tsteps)
  gamma_constant[,3] = rep(.1,tsteps)
  
  # gamma for a fast transition from topic1 to topic2 (one year/12 time steps)
  gamma_fast = matrix(rep(0,tsteps*topics),nrow=tsteps,ncol=topics)
  # proportions are constant for first 1/3 of the series
  gamma_fast[1:150,1] = rep(1)
  # fast transition from tstep 151-163
  gamma_fast[151:162,1] = seq(12)*(-1/12)+1
  gamma_fast[151:162,2] = seq(12)*(1/12)
  # topic 2 prevails for middle
  gamma_fast[163:250,2] = rep(1)
  # fast transition from 251-263
  gamma_fast[251:262,2] = seq(12)*(-1/12)+1
  gamma_fast[251:262,3] = seq(12)*(1/12)
  # proportions are constant for rest of time series
  gamma_fast[263:400,3] = rep(1)
  
  # gamma for a slow transition from topic1 to topic2
  gamma_slow = matrix(rep(0,tsteps*topics),nrow=tsteps,ncol=topics)
  gamma_slow[,1] = c(seq(tsteps/2)*(-1/(tsteps/2))+1,rep(0,tsteps/2))
  gamma_slow[,2] = c(seq(tsteps/2)*(1/(tsteps/2)),seq(tsteps/2)*(-1/(tsteps/2))+1)
  gamma_slow[,3] = c(rep(0,tsteps/2),seq(tsteps/2)*(1/(tsteps/2)))
  
  return(list(beta,gamma_constant,gamma_fast,gamma_slow))
  
}