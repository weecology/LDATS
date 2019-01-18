# Functions for using AIC for model seclection on LDA models with different numbers of topics

library(topicmodels)
library(dplyr)


#' AIC model selection for LDA using VEM method
#' 
#' Runs LDA using different numbers of topics and VEM method, calculates AIC values for comparison
#' 
#' @param dat Table of integer data (species counts by period)
#' @param SEED set seed to keep LDA model runs consistent (default 2010)
#' @param topic_min lowest number of topics; must be >=2
#' @param topic_max highest number of topics
#' 
#' @return data frame of number of topics (k) and aic value (aic)
#' 
#' @example aic_values = aic_model(dat,2010,2,10)


aic_model = function(dat,SEED,topic_min,topic_max) {
  aic_values = data.frame(seed = integer(), k = integer(), aic = numeric())
  for (k in seq(topic_min,topic_max)) {
    #run LDA
    VEM=LDA(dat,k,control=list(seed=SEED),method='VEM')
  
    #get parameter estimates
    z=posterior(VEM)
    commun.plot=z$topics
    commun.spp=z$terms
  
    #calculate AIC
    max.logl=sum(VEM@loglikelihood) #extract estimate of maximum loglikelihood
    nparam=(nrow(commun.plot)-1)*(ncol(commun.plot))+nrow(commun.spp)*(ncol(commun.spp)-1) #number of parameters
    aic=2*nparam-2*max.logl   #aic calculation
    aic_values = rbind(aic_values,data.frame(SEED, k,aic))
  }
  return(aic_values)
}





#' AIC model selection for LDA using Gibbs sampler method
#' 
#' Runs LDA using different numbers of topics and Gibbs method, calculates AIC values for comparison
#' 
#' @param dat Table of integer data (species counts by period)
#' @param ngibbs number of iterations for gibbs sampler -- must be greater than 200 (default 1000)
#' @param topic_min lowest number of topics; must be >=2
#' @param topic_max highest number of topics
#' @param save_runs T/F whether to save each run of LDA for later retrieval (these take a long time to run)
#' 
#' @return data frame containing column for number of topics (k) and aic values (aic)
#' 
#' @example aic_values = aic_model_gibbs(dat,500,2,3,T)
#' 
#' 

aic_model_gibbs = function(dat,ngibbs=1000,topic_min,topic_max,save_runs=T) {
  source('gibbs_functions.R')
  
  nspp = ncol(dat)   # number of species
  tsteps = nrow(dat) # number of time steps

  aic_values = data.frame()
  # run LDA using gibbs
  for (k in seq(topic_min,topic_max)) {
    results=gibbs.samp(dat.agg=dat,ngibbs=ngibbs,ncommun=k,a.betas=1,a.theta=1)
    if (save_runs==T) {save(results,file=paste('gibbs_results_',k,'topics'))}
    max.logl=max(results$logL) #extract estimate of maximum loglikelihood
    nparam=(tsteps-1)*(k)+nspp*(k-1) #number of parameters
    aic=2*nparam-2*max.logl   #aic calculation
    aic_values = rbind(aic_values,c(k,aic))
  }
  return(aic_values)
}



#' repeat VEM a bunch of times with different seeds and calculate AICs to find distribution of "best" ntopics
#' 
#' 
#' @param dat
#' @param seeds vector of seeds to use for analysis
#' @param topic_min lowest number of topics; must be >=2
#' @param topic_max highest number of topics
#' 
#' @example best_ntopic = repeat_VEM(dat,1:500,2,6)

repeat_VEM = function(dat,seeds,topic_min,topic_max) {
  purrr::map_df(seeds, 
         ~ aic_model(dat,SEED=.x,topic_min,topic_max) %>% 
           filter(aic == min(aic))) %>% 
    return()
}
