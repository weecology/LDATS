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
                   community = as.factor(c(rep(1,dim(gamma_constant)[1]),rep(2,dim(gamma_constant)[1]))))#,rep(3,dim(gamma_constant)[1]))))
slow = data.frame(date = rep(sim_dates,dim(gamma_slow)[2]),
                   relabund = as.vector(gamma_slow),
                   community = as.factor(c(rep(1,dim(gamma_slow)[1]),rep(2,dim(gamma_slow)[1]))))#,rep(3,dim(gamma_slow)[1]))))
fast = data.frame(date = rep(sim_dates,dim(gamma_fast)[2]),
                   relabund = as.vector(gamma_fast),
                   community = as.factor(c(rep(1,dim(gamma_fast)[1]),rep(2,dim(gamma_fast)[1]))))#,rep(3,dim(gamma_fast)[1]))))

g_1 = plot_gamma(const,2,ylab='Model Input')
g_2 = plot_gamma(fast,2)
g_3 = plot_gamma(slow,2)
grid.arrange(g_1,g_2,g_3,nrow=1)
# =================================================================================
# create data set from beta and gamma; data must be in integer form
dataset1 = round(as.data.frame(gamma_constant %*% beta) *N,digits=0)
dataset2 = ceiling(as.data.frame(gamma_fast %*% beta) *N)
dataset3 = round(as.data.frame(gamma_slow %*% beta) *N,digits=0)

# =================================================================================
# run LDA model -- VEM
SEED  = 1
topic_min = 2
topic_max = 6

#nstart = 20 # For the final analysis, maybe do 1000
ldamodel1 = LDA_analysis_VEM(dataset1,SEED,c(topic_min,topic_max))
ldamodel2 = LDA_analysis_VEM(dataset2,SEED,c(topic_min,topic_max))
ldamodel3 = LDA_analysis_VEM(dataset3,SEED,c(topic_min,topic_max))


g1 = plot_component_communities(ldamodel1,2,sim_dates,ylab='LDA model output')
g2 = plot_component_communities(ldamodel2,2,sim_dates)
g3 = plot_component_communities(ldamodel3,2,sim_dates)
grid.arrange(g1,g2,g3,nrow=1)

# =======================
# panel of figures: simulation data

c = capture_base_plot(plot_community_composition(beta))
(figure <- multi_panel_figure(
  width = c(40,40,40,40,40,40),
  height = c(40,50,50),
  panel_label_type = "upper-roman"))
figure %<>% fill_panel(
  c,
  row = 1, column = 2:5)
figure %<>% fill_panel(
  g_1,
  row = 2, column = 1:2)
figure %<>% fill_panel(
  g_2,
  row = 2, column = 3:4)
figure %<>% fill_panel(
  g_3,
  row = 2, column = 5:6)
figure %<>% fill_panel(
  g1,
  row = 3, column = 1:2)
figure %<>% fill_panel(
  g2,
  row = 3, column = 3:4)
figure %<>% fill_panel(
  g3,
  row = 3, column = 5:6)

figure


#grid.arrange(g_1,g_2,g_3,g1,g2,g3,nrow=2)

 
# =================================================================================
# changepoint model 
year_continuous = (seq(400)/12) +1977
x = data.frame(
  year_continuous=year_continuous,
  sin_year = sin(year_continuous * 2 * pi),
  cos_year = cos(year_continuous * 2 * pi))

cp_results1 = changepoint_model(ldamodel1, x, 1, weights = rep(1,length(year_continuous)))
cp_results2 = changepoint_model(ldamodel2, x, 1, weights = rep(1,length(year_continuous)))
cp_results3 = changepoint_model(ldamodel3, x, 1, weights = rep(1,length(year_continuous)))

# changepoint visualizations
par(mfrow=c(1,3))
hist(year_continuous[cp_results1$saved[,1,]],breaks=seq(1977,2016),xlab='',main='Changepoint Estimate',ylim=c(0,800))
hist(year_continuous[cp_results2$saved[,1,]],breaks=seq(1977,2016),xlab='',main='Changepoint Estimate',ylim=c(0,800))
hist(year_continuous[cp_results3$saved[,1,]],breaks=seq(1977,2016),xlab='',main='Changepoint Estimate',ylim=c(0,800))
par(mfrow=c(1,1))
