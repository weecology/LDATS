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


source('rodent_data_for_LDA.r')
source('AIC_model_selection.R')
source('LDA_figure_scripts.R')
source('changepointmodel.r')
source('LDA_analysis.R')
source('LDA-distance.R')

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
cp_results_rodent5 = changepoint_model(ldamodel, x, 5, weights = rep(1,length(year_continuous)))

hist(year_continuous[cp_results_rodent3$saved[,1,]],breaks = seq(1977,2016,.25),xlab='',main='Changepoint Estimate')
annual_hist(cp_results_rodent4,year_continuous)

# change point model selection

# mean deviance ( -2 * log likelihood) + 2*(#parameters)
mean(cp_results_rodent$saved_lls * -2) + 2*(3*(n_topics-1)*(1+1)+(1))
mean(cp_results_rodent2$saved_lls * -2)+ 2*(3*(n_topics-1)*(2+1)+(2))
mean(cp_results_rodent3$saved_lls * -2)+ 2*(3*(n_topics-1)*(3+1)+(3))
mean(cp_results_rodent4$saved_lls * -2)+ 2*(3*(n_topics-1)*(4+1)+(4))
mean(cp_results_rodent5$saved_lls * -2)+ 2*(3*(n_topics-1)*(5+1)+(5))

# =================================================================
# figures

beta1 = community_composition(ldamodel)

# put columns in order of largest species to smallest
composition = beta1[,c('NA','DS','SH','SF','SO','DO','DM','PB','PH','OL','OT','PL','PM','PE','PP','PI','RF','RM','RO','BA','PF')]
plot_community_composition(composition,c(3,4,1,2))


# with grassland communities highlighted
P = plot_community_composition_gg(composition,c(3,4,2,1))

(figure_spcomp <- multi_panel_figure(
  width = c(70,70,70,70),
  height = c(70,10),
  panel_label_type = "none",
  column_spacing = 0))
figure_spcomp %<>% fill_panel(
  P[[1]],
  row = 1, column = 1)
figure_spcomp %<>% fill_panel(
  P[[2]],
  row = 1, column = 2)
figure_spcomp %<>% fill_panel(
  P[[3]],
  row = 1, column = 3)
figure_spcomp %<>% fill_panel(
  P[[4]],
  row = 1, column = 4)
figure_spcomp

# plot of component communities over time
cc = plot_component_communities(ldamodel,ntopics,dates)



# changepoint histogram w 4 cpts
df_4 = as.data.frame(t(cp_results_rodent4$saved[,1,]))
df1_4 = data.frame(value = year_continuous[df_4$V1])
df2_4 = data.frame(value = year_continuous[df_4$V2])
df3_4 = data.frame(value = year_continuous[df_4$V3])
df4_4 = data.frame(value = year_continuous[df_4$V4])
H_4 = ggplot(data = df_4, aes(x=value)) +
  geom_histogram(data=df1_4,aes(y=..count../sum(..count..)),binwidth = .5,fill='gray1',alpha=.3) +
  geom_histogram(data=df2_4,aes(y=..count../sum(..count..)),binwidth = .5,fill='gray2',alpha=.5) +
  geom_histogram(data=df3_4,aes(y=..count../sum(..count..)),binwidth = .5,fill='gray3',alpha=.7) +
  geom_histogram(data=df4_4,aes(y=..count../sum(..count..)),binwidth = .5,fill='gray',alpha=.9) +
  labs(x='',y='') +
  xlim(range(year_continuous)) +
  theme(axis.text=element_text(size=12),
        panel.border=element_rect(colour='black',fill=NA)) +
  scale_y_continuous(labels=c('0.00','0.20','0.40','0.06','0.08'),breaks = c(0,.2,.4,.6,.8))
H_4


# find 95% confidence intervals on changepoints:
quantile(df1_4$value,probs=c(.025,.975))
quantile(df2_4$value,probs=c(.025,.975))
quantile(df3_4$value,probs=c(.025,.975))
quantile(df4_4$value,probs=c(.025,.975))
format(date_decimal(1996.07), "%d-%m-%Y")


# changepoint model plot
cpts = find_changepoint_location(cp_results_rodent4)

cpt_plot = get_ll_non_memoized_plot(ldamodel,x,cpts,make_plot=T,weights=rep(1,length(year_continuous)))


# Figure 3 -- community composition, LDA model, changepoint histogram, changepoint timeseries
(figure <- multi_panel_figure(
  width = c(70,70,70,70),
  height = c(60,60,60,60),
  column_spacing = 0))
figure %<>% fill_panel(
  figure_spcomp,
  row = 1, column = 1:4)
figure %<>% fill_panel(
  cc,
  row = 2, column = 1:4)
figure %<>% fill_panel(
  H_4,
  row = 3, column = 1:4)
figure %<>% fill_panel(
  cpt_plot,
  row = 4, column = 1:4)
figure


# ===================================================================
# appendix: LDA with 3 and 5 topics

# 3 topics
ldamodel3topic = LDA(dat,3, control = list(seed = 46),method='VEM')
cc3 = plot_component_communities(ldamodel3topic,3,dates)
beta13topic = community_composition(ldamodel3topic)
composition3 = beta13topic[,c('NA','DS','SH','SF','SO','DO','DM','PB','PH','OL','OT','PL','PM','PE','PP','PI','RF','RM','RO','BA','PF')]
P3topic = plot_community_composition_gg(composition3,c(3,2,1))

(figure_spcomp3 <- multi_panel_figure(
  width = c(70,70,70),
  height = c(70,10),
  panel_label_type = "none",
  column_spacing = 0))
figure_spcomp3 %<>% fill_panel(
  P3topic[[1]],
  row = 1, column = 1)
figure_spcomp3 %<>% fill_panel(
  P3topic[[2]],
  row = 1, column = 2)
figure_spcomp3 %<>% fill_panel(
  P3topic[[3]],
  row = 1, column = 3)
figure_spcomp3

(figure_s2 <- multi_panel_figure(
  width = c(70,70,70,70),
  height = c(60,60),
  column_spacing = 0))
figure_s2 %<>% fill_panel(
  figure_spcomp3,
  row = 1, column = 1:4)
figure_s2 %<>% fill_panel(
  cc3,
  row = 2, column = 1:4)
figure_s2

# 4 topics
(figure_s3 <- multi_panel_figure(
  width = c(70,70,70,70),
  height = c(60,60),
  column_spacing = 0))
figure_s3 %<>% fill_panel(
  figure_spcomp,
  row = 1, column = 1:4)
figure_s3 %<>% fill_panel(
  cc,
  row = 2, column = 1:4)
figure_s3

# 5 topics
ldamodel5topic = LDA(dat,5, control = list(seed = 110),method='VEM')
cc5 = plot_component_communities(ldamodel5topic,5,dates,'',c(1,5,3,4,2)) # 800 x 300
beta15topic = community_composition(ldamodel5topic)
composition5 = beta15topic[c(1,5,3,4,2),c('NA','DS','SH','SF','SO','DO','DM','PB','PH','OL','OT','PL','PM','PE','PP','PI','RF','RM','RO','BA','PF')]
P5topic = plot_community_composition_gg(composition5,c(3,4,5,2,1))

(figure_spcomp5 <- multi_panel_figure(
  width = c(60,60,60,60,60),
  height = c(60,10),
  panel_label_type = "none",
  column_spacing = 0))
figure_spcomp5 %<>% fill_panel(
  P5topic[[1]],
  row = 1, column = 1)
figure_spcomp5 %<>% fill_panel(
  P5topic[[2]],
  row = 1, column = 2)
figure_spcomp5 %<>% fill_panel(
  P5topic[[3]],
  row = 1, column = 3)
figure_spcomp5 %<>% fill_panel(
  P5topic[[4]],
  row = 1, column = 4)
figure_spcomp5 %<>% fill_panel(
  P5topic[[5]],
  row = 1, column = 5)
figure_spcomp5

(figure_s4 <- multi_panel_figure(
  width = c(70,70,70,70),
  height = c(60,60),
  column_spacing = 0))
figure_s4 %<>% fill_panel(
  figure_spcomp5,
  row = 1, column = 1:4)
figure_s4 %<>% fill_panel(
  cc5,
  row = 2, column = 1:4)
figure_s4


# three changepoints
# changepoint model plot
cpts_3 = find_changepoint_location(cp_results_rodent3)

cpt_plot_3pts = get_ll_non_memoized_plot(ldamodel,x,cpts_3,make_plot=T,weights=rep(1,length(year_continuous)))

# changepoint histogram w 4 cpts
df_3 = as.data.frame(t(cp_results_rodent3$saved[,1,]))
df1_3 = data.frame(value = year_continuous[df_3$V1])
df2_3 = data.frame(value = year_continuous[df_3$V2])
df3_3 = data.frame(value = year_continuous[df_3$V3])
H_3 = ggplot(data = df_3, aes(x=value)) +
  geom_histogram(data=df1_3,aes(y=..count../sum(..count..)),binwidth = .5,fill='gray1',alpha=.3) +
  geom_histogram(data=df2_3,aes(y=..count../sum(..count..)),binwidth = .5,fill='gray2',alpha=.6) +
  geom_histogram(data=df3_3,aes(y=..count../sum(..count..)),binwidth = .5,fill='gray3',alpha=.9) +
  labs(x='',y='') +
  xlim(range(year_continuous)) +
  theme(axis.text=element_text(size=12),
        panel.border=element_rect(colour='black',fill=NA)) +
  scale_y_continuous(labels=c('0.00','0.20','0.40','0.06','0.08'),breaks = c(0,.2,.4,.6,.8))
H_3

# Figure 3 -- community composition, LDA model, changepoint histogram, changepoint timeseries
(figure_s6 <- multi_panel_figure(
  width = c(70,70,70,70),
  height = c(60,60,60,60),
  column_spacing = 0))
figure_s6 %<>% fill_panel(
  figure_spcomp,
  row = 1, column = 1:4)
figure_s6 %<>% fill_panel(
  cc,
  row = 2, column = 1:4)
figure_s6 %<>% fill_panel(
  H_3,
  row = 3, column = 1:4)
figure_s6 %<>% fill_panel(
  cpt_plot_3pts,
  row = 4, column = 1:4)
figure_s6

