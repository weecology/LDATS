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
source('create_sim_data.R')

source('gibbs_functions.r')


# ==================================================================================
# create data
N = 200   # total number of individuals
output = create_sim_data_2topic()

beta = as.matrix(as.data.frame(output[1]))
gamma_constant = as.matrix(as.data.frame(output[2]))
gamma_fast = as.matrix(as.data.frame(output[3]))
gamma_slow = as.matrix(as.data.frame(output[4]))

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
#gamma_fast_noise = gamma_fast + rnorm(n=length(gamma_fast),mean=0,sd=.05)



# =================================================================================
# run LDA model -- VEM
nstart = 20 # For the final analysis, maybe do 1000
ldamodel = LDA(dataset2,2,control=list(estimate.alpha=F,alpha=.1, nstart = nstart),method="VEM")

plot_component_communities(ldamodel,2,seq(400))

# # AIC model selection for number of topics
# aic_values1 = aic_model(dataset1)
# aic_values2 = aic_model(dataset2)
# aic_values3 = aic_model(dataset3)

# =================================================================================
# run LDA model -- Gibbs
# nstart = 20 # For the final analysis, maybe do 1000
# ldamodel = LDA(dataset2,2,control=list(iter=10000,delta=1),method="Gibbs")
# ldamodel = LDA(dataset3,3,control=list(iter=10000,alpha=.1),method="Gibbs")
# 
# plot_component_communities(ldamodel,3,seq(400))
# 
# plot_community_composition(community_composition(ldamodel2),c(0,.2))
# 
# 
# # AIC model selection for number of topics
# aic_values1 = aic_model(dataset1)
# aic_values2 = aic_model(dataset2)
# aic_values3 = aic_model(dataset3)
# 
# # =================================================================================
# # run LDA model -- Denis Valle's version of Gibbs
# 
# ngibbs=1000 #has to be greater than 200
# ncommun=2
# results=gibbs.samp(dat.agg=dataset3,ngibbs=ngibbs,ncommun=ncommun,a.betas=1,a.theta=25)
# 
# # plots
# beta1=matrix(apply(results$beta,2,mean),ncommun,nspecies)
# plot_community_composition(beta1,c(0,.2))
# 
# plot_component_communities_gibbs(results,ncommun,seq(400))
# 
# # AIC
# aic_values1 = aic_model_gibbs(dataset1,nspecies,tsteps)
# aic_values2 = aic_model_gibbs(dataset2,nspecies,tsteps)
# aic_values3 = aic_model_gibbs(dataset3,nspecies,tsteps)

# =================================================================================
# changepoint model 
year_continuous = seq(400)/12
x = data.frame(
  year_continuous=year_continuous,
  sin_year = sin(year_continuous * 2 * pi),
  cos_year = cos(year_continuous * 2 * pi))
dat = dataset2

cp_results = changepoint_model(ldamodel, x, 2)

save(cp_results,file='C:/Users/EC/Desktop/git/Extreme-events-LDA/changepoint results/chpoint_2topics_fastgamma_VEM')
# changepoint visualizations
par(mfrow=c(1,1))
annual_hist(cp_results,year_continuous)
