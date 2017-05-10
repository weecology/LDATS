# LDA and changepoint analysis pipeline on rodent data -- uses VEM method
#
#  1. prepare data
#      - specific script for each data set
#  2a. run LDA with different number of topics and use AIC to select best model
#  2b. run LDA with best number of topics (determined by 2a)
#  3. run changepoint model
#  4. produce figures


library(topicmodels)
library(RCurl)


source('rodent_data_for_LDA.R')
source('AIC_model_selection.R')
source('LDA_figure_scripts.R')


# ===================================================================
# prepare rodent data
dat = create_rodent_table(period_first = 1,
                          period_last = 436,
                          selected_plots = c(2,4,8,11,12,14,17,22),
                          selected_species = c('BA','DM','DO','DS','NA','OL','OT','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO','SF','SH','SO'))

# dates to go with count data
moondat = read.csv(text=getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/moon_dates.csv"),stringsAsFactors = F)
moondat$date = as.Date(moondat$CensusDate)

period_dates = filter(moondat,Period %in% rownames(dat)) %>% select(Period,date)
dates = period_dates$date



# ===================================
# model parameters:
SEED = 113052032
topic_min = 2
topic_max = 9
nspp=length(dat)
# ==================================================================
# select number of topics

# choose number of topics -- model selection using AIC
aic_values = aic_model(dat,SEED,topic_min,topic_max)
ntopics = filter(aic_values,aic==min(aic)) %>% select(k) %>% as.numeric()

# repeat aic calculation with a bunch of different seeds to test robustness of the analysis
best_ntopic = repeat_VEM(dat,
                         ntimes=200,
                         topic_min=2,
                         topic_max=9)
#hist(best_ntopic[,1],breaks=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5),xlab='best # of topics')

# ==================================================================
# run LDA model

ldamodel = LDA_analysis_VEM(dat,SEED,c(topic_min,topic_max))

# ==================================================================
# change point model 

 year_continuous = 1970 + as.integer(julian(dates)) / 365.25
 x = data.frame(
   year_continuous = year_continuous,
   sin_year = sin(year_continuous * 2 * pi),
   cos_year = cos(year_continuous * 2 * pi)
 )
par(mfrow=c(1,1))
cp_results_rodent = changepoint_model(ldamodel, x, 1, weights = rep(1,length(year_continuous)))
cp_results_rodent2 = changepoint_model(ldamodel, x, 2, weights = rep(1,length(year_continuous)))
cp_results_rodent3 = changepoint_model(ldamodel, x, 3, weights = rep(1,length(year_continuous)))
cp_results_rodent4 = changepoint_model(ldamodel, x, 4, weights = rep(1,length(year_continuous)))
hist(year_continuous[cp_results_rodent$saved[,1,]],breaks = seq(1977,2016),xlab='',main='Changepoint Estimate')

# =================================================================
# figures

beta1 = community_composition(ldamodel)
# put columns in order of largest species to smallest
composition = beta1[,c('NA','DS','SH','SF','SO','DO','DM','PB','PH','OL','OT','PL','PM','PE','PP','PI','RF','RM','RO','BA','PF')]
plot_community_composition(composition,c(3,4,2,1))

plot_component_communities(ldamodel,ntopics,dates)

# ===================================================================
# appendix: LDA with 3 and 5 topics

# 3 topics
ldamodel3topic = LDA(dat,3, control = list(seed = SEED),method='VEM')
plot_component_communities(ldamodel3topic,3,dates)
beta13topic = community_composition(ldamodel3topic)
composition3 = beta13topic[,c('NA','DS','SH','SF','SO','DO','DM','PB','PH','OL','OT','PL','PM','PE','PP','PI','RF','RM','RO','BA','PF')]
plot_community_composition(composition3,c(3,2,1))

# 5 topics
ldamodel5topic = LDA(dat,5, control = list(seed = SEED),method='VEM')
plot_component_communities(ldamodel5topic,5,dates) # 800 x 300
beta15topic = community_composition(ldamodel5topic)
composition5 = beta15topic[,c('NA','DS','SH','SF','SO','DO','DM','PB','PH','OL','OT','PL','PM','PE','PP','PI','RF','RM','RO','BA','PF')]
plot_community_composition(composition5,c(3,4,1,2,5))

