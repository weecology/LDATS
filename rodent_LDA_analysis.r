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
library(multipanelfigure)


source('rodent_data_for_LDA.R')
source('AIC_model_selection.R')
source('LDA_figure_scripts.R')
source('changepointmodel.r')
source('LDA_analysis.R')

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
topic_min = 2
topic_max = 9
nspp=length(dat)
# ==================================================================
# select number of topics

# Fit a bunch of LDA models with different seeds
# Only use every other seed because consecutive seeds give identical results (!?)
seeds = 2*seq(200)

# repeat aic calculation with a bunch of different seeds to test robustness of the analysis
best_ntopic = repeat_VEM(dat,
                         seeds,
                         topic_min=2,
                         topic_max=6)

# plot histogram of how many seeds chose how many topics
hist(best_ntopic[,1],breaks=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5),xlab='best # of topics', main='')

# ==================================================================
# how different is species composition of 4 community-types when LDA is run with different seeds?

# get list of 100 seeds where 4 topics was the best LDA model
seeds_4topics = data.frame(best_ntopic) %>% filter(X1 == 4) %>% select(X2) %>% head(100) %>% unlist() %>% as.numeric()

best_seed = calculate_LDA_distance(dat,seeds_4topics)


# ==================================================================
# run LDA model

ntopics = 4
#SEED = 113052032
SEED = best_seed
ldamodel = LDA(dat,ntopics, control = list(seed = SEED),method='VEM')


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
hist(year_continuous[cp_results_rodent3$saved[,1,]],breaks = seq(1977,2016,.25),xlab='',main='Changepoint Estimate')
annual_hist(cp_results_rodent4,year_continuous)

# =================================================================
# figures

beta1 = community_composition(ldamodel)

# put columns in order of largest species to smallest
composition = beta1[,c('NA','DS','SH','SF','SO','DO','DM','PB','PH','OL','OT','PL','PM','PE','PP','PI','RF','RM','RO','BA','PF')]
plot_community_composition(composition,c(3,4,1,2))


# with grassland communities highlighted
P = plot_community_composition_gg(composition,c(3,4,1,2))
(figure <- multi_panel_figure(
  width = c(70,70,70,70),
  height = c(70,10),
  panel_label_type = "none",
  column_spacing = 0))
figure %<>% fill_panel(
  P[[1]],
  row = 1, column = 1)
figure %<>% fill_panel(
  P[[2]],
  row = 1, column = 2)
figure %<>% fill_panel(
  P[[3]],
  row = 1, column = 3)
figure %<>% fill_panel(
  P[[4]],
  row = 1, column = 4)
figure

cc = plot_component_communities(ldamodel,ntopics,dates)

D = capture_base_plot(plot_community_composition(composition,c(3,4,2,1)))

H = capture_base_plot(annual_hist(cp_results_rodent4,year_continuous))


# Hot mess -- fix this!! (changepoint histogram)
df = as.data.frame(t(cp_results_rodent3$saved[,1,]))
df1 = data.frame(value = year_continuous[df$V1])
df2 = data.frame(value = year_continuous[df$V2])
df3 = data.frame(value = year_continuous[df$V3])
H = ggplot(data = df, aes(x=value)) +
  geom_histogram(data=df1,aes(y=..count../sum(..count..)),binwidth = .25,fill='gray1',alpha=.3) +
  geom_histogram(data=df2,aes(y=..count../sum(..count..)),binwidth = .25,fill='gray2',alpha=.5) +
  geom_histogram(data=df3,aes(y=..count../sum(..count..)),binwidth = .25,fill='gray3',alpha=1) +
  labs(x='',y='') +
  xlim(range(year_continuous)) +
  theme(axis.text=element_text(size=12),
      panel.border=element_rect(colour='black',fill=NA))
  
H

# 95% confidence intervals on changepoints:
quantile(df1$value,probs=c(.025,.975))
quantile(df2$value,probs=c(.025,.975))
quantile(df3$value,probs=c(.025,.975))
format(date_decimal(1983.933), "%d-%m-%Y")




# changepoint histogram w 4 cpts
df_4 = as.data.frame(t(cp_results_rodent4$saved[,1,]))
df1_4 = data.frame(value = year_continuous[df_4$V1])
df2_4 = data.frame(value = year_continuous[df_4$V2])
df3_4 = data.frame(value = year_continuous[df_4$V3])
df4_4 = data.frame(value = year_continuous[df_4$V4])
H_4 = ggplot(data = df_4, aes(x=value)) +
  geom_histogram(data=df1_4,aes(y=..count../sum(..count..)),binwidth = .25,fill='gray1',alpha=.3) +
  geom_histogram(data=df2_4,aes(y=..count../sum(..count..)),binwidth = .25,fill='gray2',alpha=.5) +
  geom_histogram(data=df3_4,aes(y=..count../sum(..count..)),binwidth = .25,fill='gray3',alpha=.7) +
  geom_histogram(data=df4_4,aes(y=..count../sum(..count..)),binwidth = .25,fill='gray',alpha=.9) +
  labs(x='',y='') +
  xlim(range(year_continuous)) +
  theme(axis.text=element_text(size=12),
        panel.border=element_rect(colour='black',fill=NA))

H_4


# changepoint model plot
cpts = find_changepoint_location(cp_results_rodent4)

get_ll_non_memoized(ldamodel,x,cpts,make_plot=T,weights=rep(1,length(year_continuous)))


# Figure 3 -- LDA model, changepoint histogram, changepoint timeseries
(figure <- multi_panel_figure(
  width = c(120,120),
  height = c(80,80),
  panel_label_type = "upper-roman"))
figure %<>% fill_panel(
  cc,
  row = 1, column = 1:2)
figure %<>% fill_panel(
  H_4,
  row = 2, column = 1:2)
figure

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

