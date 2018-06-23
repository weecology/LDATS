
#' @title create series of simulated data
#' 
#' @param total_len total length of time series (number of steps)
#' @param change_len length of gradiated change (number of steps)
#' @param ntopics number of "topics" (columns in output frame)
#' 
#' @return matrix with columns representing topics and rows representing time steps
#' 
create_sim_series = function(total_len,change_len,ntopics) {
  
  start_change = floor(total_len/2 - change_len/2)
  end_change = start_change + change_len
  
  sim_series = matrix(rep(0,total_len*ntopics),nrow=total_len,ncol=ntopics)
  # proportions are constant before change
  sim_series[1:start_change,1] = rep(1)
  sim_series[1:start_change,2] = rep(0)
  # duration of change
  sim_series[(start_change+1):end_change,1] = seq(change_len)*(-1/change_len)+1
  sim_series[(start_change+1):end_change,2] = seq(change_len)*(1/change_len)+0
  # proportions are constant for rest of time series
  sim_series[(end_change+1):total_len,1] = rep(0)
  sim_series[(end_change+1):total_len,2] = rep(1)
  
  return(sim_series)
}

#' create simulated data for demonstrating LDA analysis
#'    2 topics
#'     
#' @param nspecies = number of species in all topic groups
#' @param tsteps = number of [monthly] time steps
#' @param N = total number of individuals of all species
#' 
#' @return 
#'    beta = matrix of species composition of the groups
#'    gamma = matrix of topic composition over time
#'            3 simulations of gamma: uniform, slow transition, and fast transition
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
  gamma_fast = create_sim_series(total_len=tsteps,change_len=12,topics)
  
  # gamma for a slow transition from topic1 to topic2
  gamma_slow = create_sim_series(total_len=tsteps,change_len=300,topics)
  
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
  gamma_1yr = create_sim_series(total_len=tsteps,change_len=12,topics)
  
  # gamma for a fast transition from topic1 to topic2 (two year/24 time steps)
  gamma_2yr = create_sim_series(total_len=tsteps,change_len=24,topics)
  
  # gamma for a slow transition from topic1 to topic2: 5 years, 60 time steps
  gamma_5yr = create_sim_series(total_len=tsteps,change_len=60,topics)
  
  # gamma for a slow transition from topic1 to topic2
  gamma_25yr = create_sim_series(total_len=tsteps,change_len=300,topics)
  
  return(list(beta,gamma_constant,gamma_1yr,gamma_2yr,gamma_5yr,gamma_25yr))
}



create_sim_data_3topic = function(nspecies=24,tsteps=400) {
  # create beta and gammas; 3 topics
  #  Inputs:
  #    nspecies = number of species in all topic groups
  #    tsteps = number of [monthly] time steps
  #    N = total number of individuals of all species
  #  Outputs:
  #    beta = matrix of species composition of the groups
  #    gamma = matrix of topic composition over time
  #            3 simulations of gamma: uniform, slow transition, and fast transition

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