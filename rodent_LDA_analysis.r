# LDA analysis (with figures) for the extreme events/ LDA project
library(topicmodels)
library(ggplot2)
library(lubridate)

source('LDA_figure_scripts.R')

# =========================================================================================
# LDA analyses
set.seed(20)

setwd('C:/Users/EC/Desktop/git/Extreme-events-LDA')

# load csv of abundance of each sp per plot
dat = read.csv('Rodent_table_dat.csv',na.strings = '',as.is=T)

# dates of trapping periods
perdat = read.csv('period_dates_single.csv')

perdat$date = as.Date(perdat$date,format='%m/%d/%Y')



# LDA models: groups from 2 to 7
nstart = 200 # For the final analysis, maybe do 1000
ldamodel2 = LDA(dat,2,control=list(estimate.alpha=F,alpha=.5, nstart = nstart),method="VEM")
ldamodel3 = LDA(dat,3,control=list(estimate.alpha=F,alpha=.5, nstart = nstart),method="VEM")
ldamodel4 = LDA(dat,4,control=list(estimate.alpha=F,alpha=.5, nstart = nstart),method="VEM")
ldamodel5 = LDA(dat,5,control=list(estimate.alpha=F,alpha=.5, nstart = nstart),method="VEM")
ldamodel5 = LDA(dat,5,control=list(estimate.alpha=F,alpha=.5, nstart = nstart),method="VEM")
ldamodel6 = LDA(dat,6,control=list(estimate.alpha=F,alpha=.5, nstart = nstart),method="VEM")
ldamodel7 = LDA(dat,7,control=list(estimate.alpha=F,alpha=.5, nstart = nstart),method="VEM")


# model selection
aic_values = aic_model(dat)

# gibbs method and aic model selection
source('AIC_model_selection.R')
nspecies = ncol(dat)
tsteps = nrow(dat)
aic_values = aic_model_gibbs(dat,nspecies,tsteps)



# ==========================================================================================
# figures

# look at sp comp of topics
comp_matrix = community_composition(ldamodel2)
comp_matrix

# make figure of sp comp of topics
plot_community_composition(comp_matrix,c(0,.6))

# time series plot of topics
dates = as.Date(perdat$date[1:length(dat[,1])])

plot_component_communities(ldamodel2,2,dates)
plot_component_communities(ldamodel3,3,dates)
plot_component_communities(ldamodel4,4,dates)
plot_component_communities(ldamodel5,5,dates)
plot_component_communities(ldamodel6,6,dates)
plot_component_communities(ldamodel7,7,dates)

# gibbs method
ncommun=5
plot_component_communities_gibbs(results,ncommun,dates)

# time series plot of topics, smoothed
plot_component_communities_smooth(ldamodel2,2,dates)

# =========================================================================================
# running changepoint model

source('changepointmodel.r')

d = ymd(as.Date(perdat$date[1:length(dat[,1])]))
year_continuous = 1970 + as.integer(julian(d)) / 365.25
x = data.frame(
  year_continuous = year_continuous,
  sin_year = sin(year_continuous * 2 * pi),
  cos_year = cos(year_continuous * 2 * pi)
)

results = changepoint_model(ldamodel3, x, 2, maxit = 1000)
 