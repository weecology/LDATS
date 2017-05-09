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

library(gridExtra)


source('LDA_figure_scripts.R')
source('AIC_model_selection.R')
source('changepointmodel.r')
source('create_sim_data.R')
source('LDA_analysis.R')

#source('gibbs_functions.r')

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
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
plot_community_composition(beta)

sim_dates = seq.Date(from=as.Date('1977-01-01'),by=30,length.out = 400) 


const = data.frame(date = rep(sim_dates,dim(gamma_constant)[2]),
                   relabund = as.vector(gamma_constant),
                   community = as.factor(c(rep(1,dim(gamma_constant)[1]),rep(2,dim(gamma_constant)[1]),rep(3,dim(gamma_constant)[1]))))
slow = data.frame(date = rep(sim_dates,dim(gamma_slow)[2]),
                   relabund = as.vector(gamma_slow),
                   community = as.factor(c(rep(1,dim(gamma_slow)[1]),rep(2,dim(gamma_slow)[1]),rep(3,dim(gamma_slow)[1]))))
fast = data.frame(date = rep(sim_dates,dim(gamma_fast)[2]),
                   relabund = as.vector(gamma_fast),
                   community = as.factor(c(rep(1,dim(gamma_fast)[1]),rep(2,dim(gamma_fast)[1]),rep(3,dim(gamma_fast)[1]))))

g_1 = plot_gamma(const)
g_2 = plot_gamma(fast)
g_3 = plot_gamma(slow)
grid.arrange(g_1,g_2,g_3,nrow=1)
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
SEED  = 1
topic_min = 2
topic_max = 6

#nstart = 20 # For the final analysis, maybe do 1000
ldamodel1 = LDA_analysis_VEM(dataset1,SEED,c(topic_min,topic_max))
ldamodel2 = LDA_analysis_VEM(dataset2,SEED,c(topic_min,topic_max))
ldamodel3 = LDA_analysis_VEM(dataset3,SEED,c(topic_min,topic_max))


g1 = plot_component_communities(ldamodel1,2,sim_dates)
g2 = plot_component_communities(ldamodel2,3,sim_dates)
g3 = plot_component_communities(ldamodel3,3,sim_dates)
grid.arrange(g1,g2,g3,nrow=1)

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
#year_continuous = sim_dates
x = data.frame(
  year_continuous=year_continuous,
  sin_year = sin(year_continuous * 2 * pi),
  cos_year = cos(year_continuous * 2 * pi))
#dat = dataset2

cp_results1 = changepoint_model(ldamodel1, x, 1)
cp_results2 = changepoint_model(ldamodel2, x, 1)
cp_results3 = changepoint_model(ldamodel3, x, 1)

save(cp_results1,file='C:/Users/EC/Desktop/git/Extreme-events-LDA/changepoint results/chpoint_2topics_constgamma_VEM')
save(cp_results2,file='C:/Users/EC/Desktop/git/Extreme-events-LDA/changepoint results/chpoint_2topics_fastgamma_VEM')
save(cp_results3,file='C:/Users/EC/Desktop/git/Extreme-events-LDA/changepoint results/chpoint_2topics_slowgamma_VEM')
# changepoint visualizations
par(mfrow=c(1,1))
annual_hist(cp_results3,year_continuous)
