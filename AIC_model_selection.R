
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
