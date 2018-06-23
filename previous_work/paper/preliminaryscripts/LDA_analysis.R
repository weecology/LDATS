
#  1. prepare data
#      - specific script for each data set
#  2a. run LDA with different number of topics and use AIC to select best model
#  2b. run LDA with best number of topics (determined by 2a)
#  3. run changepoint model
#  4. produce figures



# =======================================================================================

#' Run through whole pipeline of LDA analysis, VEM method
#'
#'
#'
#' @param dat Table of integer data (species counts by time step)
#' @param SEED set seed to keep LDA model runs consistent (default 2010)
#' @param test_topics numbers of topics to test, of the form (topic_min,topic_max)
#' @param dates vector of dates that correspond to time steps from 'dat' -- for plotting purposes
#' @param n_chpoints number of change points the changepoint model should look for
#' @param maxit max iterations for changepoint model (more=slower)
#' 
#' @example LDA_analysis(dat,2010,c(2,3),dates,2)

source('AIC_model_selection.R')
source('gibbs_functions.R')
source('changepointmodel.r')


LDA_analysis_VEM = function(dat,SEED,test_topics) {
  
  # choose number of topics -- model selection using AIC
  aic_values = aic_model(dat,SEED,test_topics[1],test_topics[2])
  
  # run LDA model
  ntopics = filter(aic_values,aic==min(aic)) %>% select(k) %>% as.numeric()
  ldamodel = LDA(dat,ntopics, control = list(seed = SEED),method='VEM')

  return(ldamodel)
}

# =======================================================================================

#' Run through whole pipeline of LDA analysis, Gibbs method --- WIP
#'
#'
#'
#' @param dat Table of integer data (species counts by time step)
#' @param ngibbs number of iterations for gibbs sampler -- must be greater than 200 (default 1000)
#' @param test_topics numbers of topics to test, of the form (topic_min,topic_max)
#' @param dates vector of dates that correspond to time steps from 'dat' -- for plotting purposes
#' 
#' @example LDA_analysis_gibbs(dat,200,c(2,3),dates)

source('AIC_model_selection.R')
source('gibbs_functions.R')
source('changepointmodel.r')
source('LDA_figure_scripts.R')



LDA_analysis_gibbs = function(dat,ngibbs,test_topics,dates) {
  
  # choose number of topics -- model selection using AIC
  aic_values = aic_model_gibbs(dat,ngibbs,test_topics[1],test_topics[2],F)
  
  # run LDA model
  ntopics = filter(aic_values,aic==min(aic)) %>% select(k) %>% as.numeric()
  ldamodel = gibbs.samp(dat.agg=dat,ngibbs=ngibbs,ncommun=ntopics,a.betas=1,a.theta=1)

  
}
