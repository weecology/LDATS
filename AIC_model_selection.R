
# AIC and model selection

aic_model = function(dat) {
  
  #run LDA
  SEED=2010
  # run model for number of groups ranging from 2 to 10
  aic_values = c()
  for (k in seq(2,19)) {
    VEM=LDA(dat,k, control = list(seed = SEED),method='VEM')
  
    #get parameter estimates
    z=posterior(VEM)
    commun.plot=z$topics
    commun.spp=z$terms
  
    #calculate AIC
    max.logl=sum(VEM@loglikelihood) #extract estimate of maximum loglikelihood
    nparam=(nrow(commun.plot)-1)*(ncol(commun.plot))+nrow(commun.spp)*(ncol(commun.spp)-1) #number of parameters
    aic=2*nparam-2*max.logl   #aic calculation
    aic_values = rbind(aic_values,c(k,aic))
  }
  return(aic_values)
}


# AIC and model selection for LDA using gibbs sampler

aic_model_gibbs = function(dat,nspp,tsteps) {
  source('gibbs_functions.R')
  
  ngibbs=1000 #has to be greater than 200
  aic_values = c()
  # run LDA using gibbs
  for (k in seq(2,7)) {
    results=gibbs.samp(dat.agg=dat,ngibbs=ngibbs,ncommun=k,a.betas=1,a.theta=1)
    save(results,file=paste('gibbs_results_',k,'topics'))
    max.logl=max(results$logL) #extract estimate of maximum loglikelihood SUM or MAX?
    nparam=(tsteps-1)*(k)+nspp*(k-1) #number of parameters
    aic=2*nparam-2*max.logl   #aic calculation
    aic_values = rbind(aic_values,c(k,aic))
  }
  return(aic_values)
}


# WAIC -- doesn't work yet

waic_model_gibbs = function(lda) {
  ll = purrr::map_dbl(lda@fitted, function(x) x@loglikelihood)
  # "Effective sample size:" should probably be at least 1000.
  # Should also run multiple MCMC chains and evaluate coda::gelman.diag
  coda::effectiveSize(ll)
  # Numerically stable way to compute the mean likelihood (as opposed to mean
  # _log_-likelihood)
  log_mean_exp = function(x) {
    m = max(x)
    x = x - m
    log(mean(exp(x))) + m
  }
  # WAIC (specifically, WAIC2 as defined by Equations 11 & 13 in
  # http://www.stat.columbia.edu/~gelman/research/published/waic_understand3.pdf)
  lppd = log_mean_exp(ll)
  p_waic2 = var(ll)
  -2 * (lppd - p_waic2)
  
}
