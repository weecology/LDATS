# Use simulated data to investigate the results of the LDA model with known dynamics:
#    1. fast transition between 2 topics with steady-state before and after
#    2. slow transition between 2 topics with little to no steady-state before/after
#    3. fixed proportion of 2 topics through time (no transition)


# Started EMC 3/2017

# ====================================================================================
# Create simulated data:
#  beta = matrix of species composition of the topics
#  gamma = matrix of prevalence of topics through time
#  assume even species composition and no overlap of species between topics

setwd('C:/Users/EC/Desktop/git/Extreme-events-LDA')
source('LDA_figure_scripts.R')
source('AIC_model_selection.R')
source('changepointmodel.r')

source('gibbs_functions.r')


# ====================================================
# create beta and gammas; 2 topics

nspecies = 24
topics = 2
tsteps = 400 #I think of these as monthly time steps
N = 200      # total number of "animals" -- LDA model depends on integer counts as data


beta = matrix(rep(0,topics*nspecies),nrow=topics,ncol=nspecies)
beta[1,] = c(rep(.1,nspecies/2),rep(0,nspecies/2))
beta[2,] = c(rep(0,nspecies/2),rep(.1,nspecies/2))

# gamma for a constant topic prevalence through time: topic1 at 90% and topic2 at 10%
gamma_constant = matrix(rep(0,tsteps*topics),nrow=tsteps,ncol=topics)
gamma_constant[,1] = rep(1,tsteps)

# gamma for a fast transition from topic1 to topic2 (one year/12 time steps)
gamma_fast = matrix(rep(0,tsteps*topics),nrow=tsteps,ncol=topics)
# proportions are constant for first 200 time steps
gamma_fast[1:200,1] = rep(.9)
gamma_fast[1:200,2] = rep(.1)
# fast transition from tstep 201-212
gamma_fast[201:212,1] = seq(12)*(-.8/12)+.9
gamma_fast[201:212,2] = seq(12)*(.8/12)+.1
# proportions are constant for rest of time series
gamma_fast[213:tsteps,1] = rep(.1)
gamma_fast[213:tsteps,2] = rep(.9) 

# gamma for a slow transition from topic1 to topic2 that takes the entire time series
gamma_slow = matrix(rep(0,tsteps*topics),nrow=tsteps,ncol=topics)
gamma_slow[,1] = seq(tsteps)*(-.8/tsteps)+.9
gamma_slow[,2] = seq(tsteps)*(.8/tsteps)+.1

# ====================================================
# create beta and gammas; 3 topics

nspecies = 24
topics = 3
tsteps = 400 #I think of these as monthly time steps
N = 200      # total number of "animals" -- LDA model depends on integer counts as data


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

# gamma for a slow transition from topic1 to topic2 that takes the entire time series
gamma_slow = matrix(rep(0,tsteps*topics),nrow=tsteps,ncol=topics)
gamma_slow[,1] = c(seq(tsteps/2)*(-1/(tsteps/2))+1,rep(0,tsteps/2))
gamma_slow[,2] = c(seq(tsteps/2)*(1/(tsteps/2)),seq(tsteps/2)*(-1/(tsteps/2))+1)
gamma_slow[,3] = c(rep(0,tsteps/2),seq(tsteps/2)*(1/(tsteps/2)))


# ==================================================================================
# plot beta and gammas
plot_community_composition(beta,c(0,.2))

plot_gammas = function(gamma_constant,gamma_fast,gamma_slow) {
  par(mfrow=c(1,3))
  topics = dim(gamma_constant)[2]
  plot(gamma_constant[,1],type='l',xlab='time',ylab='topic',ylim=c(0,1),main='constant')
  for (j in seq(2,topics)) {lines(gamma_constant[,j],col=j)}
  plot(gamma_fast[,1],type='l',xlab='time',ylab='topic',ylim=c(0,1),main='fast')
  for (j in seq(2,topics)) {lines(gamma_fast[,j],col=j)}
  plot(gamma_slow[,1],type='l',xlab='time',ylab='topic',ylim=c(0,1),main='slow')
  for (j in seq(2,topics)) {lines(gamma_slow[,j],col=j)}
}

plot_gammas(gamma_constant,gamma_fast,gamma_slow)


# =================================================================================
# create data set from beta and gamma; data must be in integer form
dataset1 = round(as.data.frame(gamma_constant %*% beta) *N,digits=0)
dataset2 = ceiling(as.data.frame(gamma_fast %*% beta) *N)
dataset3 = round(as.data.frame(gamma_slow %*% beta) *N,digits=0)

# ================================================================================
# option to add noise to datasets
gamma_fast_noise = gamma_fast + rnorm(n=length(gamma_fast),mean=0,sd=.05)



# =================================================================================
# run LDA model -- VEM
nstart = 20 # For the final analysis, maybe do 1000
ldamodel = LDA(dataset1,2,control=list(estimate.alpha=F,alpha=.5, nstart = nstart),method="VEM")

plot_component_communities(ldamodel,3,seq(400))

# AIC model selection for number of topics
aic_values1 = aic_model(dataset1)
aic_values2 = aic_model(dataset2)
aic_values3 = aic_model(dataset3)


# =================================================================================
# run LDA model -- Gibbs

ngibbs=1000 #has to be greater than 200
ncommun=3
results=gibbs.samp(dat.agg=dataset1,ngibbs=ngibbs,ncommun=ncommun,a.betas=1,a.theta=1)

# plots
beta1=matrix(apply(results$beta,2,mean),ncommun,nspecies)
plot_community_composition(beta1,c(0,1))

plot_component_communities_gibbs(results,ncommun,seq(400))

# AIC
aic_values1 = aic_model_gibbs(dataset1,nspecies,tsteps)
aic_values2 = aic_model_gibbs(dataset2,nspecies,tsteps)
aic_values3 = aic_model_gibbs(dataset3,nspecies,tsteps)

# =================================================================================
# changepoint model  -- doesn't work yet
#year_continuous = 1970+seq(400)/12
#x = data.frame(
#  year_continuous = year_continuous,
#  sin_year = sin(year_continuous * 2 * pi),
#  cos_year = cos(year_continuous * 2 * pi)
#)
#results = changepoint_model(ldamodel, x, 2)
