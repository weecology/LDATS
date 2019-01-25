# LDA and changepoint analysis pipeline on rodent data -- uses VEM method
#
#  1. prepare data
#  2a. run LDA with different number of topics and use AIC to select best model
#  2b. how much does seed choice affect species composition of component communities?
#  3. run LDA with best number of topics (determined by 2a)
#  4. run changepoint model
#  5. produce figures
#  6. appendix figures


library(topicmodels)
library(RCurl)
library(multipanelfigure)
library(reshape2)


source('rodent_data_for_LDA.r')
source('AIC_model_selection.R')
source('LDA_figure_scripts.R')
source('changepointmodel.r')
source('LDA-distance.R')

# ===================================================================
# 1. prepare rodent data
# ===================================================================
dat = create_rodent_table(period_first = 1,
                          period_last = 436,
                          selected_plots = c(2,4,8,11,12,14,17,22),
                          selected_species = c('BA','DM','DO','DS','NA','OL','OT','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO','SF','SH','SO'))

# dates to go with count data
moondat = read.csv(text=getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/moon_dates.csv"),stringsAsFactors = F)
moondat$date = as.Date(moondat$censusdate)

period_dates = filter(moondat,period %in% rownames(dat)) %>% select(period,date)
dates = period_dates$date


# ==================================================================
# 2a. select number of topics
# ==================================================================

# Fit a bunch of LDA models with different seeds
# Only use even numbers for seeds because consecutive seeds give identical results
seeds = 2*seq(200)

# repeat LDA model fit and AIC calculation with a bunch of different seeds to test robustness of the analysis
best_ntopic = repeat_VEM(dat,
                         seeds,
                         topic_min=2,
                         topic_max=6)

# histogram of how many seeds chose how many topics
hist(best_ntopic$k,breaks=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5),xlab='best # of topics', main='')

# ==================================================================
# 2b. how different is species composition of 4 community-types when LDA is run with different seeds?
# ==================================================================
# get the best 100 seeds where 4 topics was the best LDA model
seeds_4topics = best_ntopic %>% 
  filter(k == 4) %>% 
  arrange(aic) %>% 
  head(100) %>% 
  pull(SEED)

# choose seed with highest log likelihood for all following analyses
#    (also produces plot of community composition for 'best' run compared to 'worst')
best_seed = calculate_LDA_distance(dat,seeds_4topics)
mean_dist = unlist(best_seed)[2]
max_dist = unlist(best_seed)[3]

# ==================================================================
# 3. run LDA model
# ==================================================================
ntopics = 4
SEED = unlist(best_seed)[1]  # For the paper, I use seed 206
ldamodel = LDA(dat,ntopics, control = list(seed = SEED),method='VEM')


# ==================================================================
# 4. change point model 
# ==================================================================

# set up parameters for model
year_continuous = 1970 + as.integer(julian(dates)) / 365.25
x = data.frame(
  year_continuous = year_continuous,
  sin_year = sin(year_continuous * 2 * pi),
  cos_year = cos(year_continuous * 2 * pi)
)

# run models with 1, 2, 3, 4, 5 changepoints
cp_results_rodent = changepoint_model(ldamodel, x, 1, weights = rep(1,length(year_continuous)))
cp_results_rodent2 = changepoint_model(ldamodel, x, 2, weights = rep(1,length(year_continuous)))
cp_results_rodent3 = changepoint_model(ldamodel, x, 3, weights = rep(1,length(year_continuous)))
cp_results_rodent4 = changepoint_model(ldamodel, x, 4, weights = rep(1,length(year_continuous)))
cp_results_rodent5 = changepoint_model(ldamodel, x, 5, weights = rep(1,length(year_continuous)))

# some quick histograms of changepoint model results
hist(year_continuous[cp_results_rodent4$saved[,1,]],breaks = seq(1977,2016,.25),xlab='',main='Changepoint Estimate')
annual_hist(cp_results_rodent4,year_continuous)

# turn changepoint results into data frame
df_4 = as.data.frame(t(cp_results_rodent4$saved[,1,])) %>% melt()
df_4$value = year_continuous[df_4$value]

# find 95% confidence intervals on each changepoint:
quantile(df_4[df_4$variable=='V1','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
quantile(df_4[df_4$variable=='V2','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
quantile(df_4[df_4$variable=='V3','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
quantile(df_4[df_4$variable=='V4','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')

# change point model selection
# mean deviance ( -2 * log likelihood) + 2*(#parameters)
mean(cp_results_rodent$saved_lls * -2) + 2*(3*(ntopics-1)*(1+1)+(1))
mean(cp_results_rodent2$saved_lls * -2)+ 2*(3*(ntopics-1)*(2+1)+(2))
mean(cp_results_rodent3$saved_lls * -2)+ 2*(3*(ntopics-1)*(3+1)+(3))
mean(cp_results_rodent4$saved_lls * -2)+ 2*(3*(ntopics-1)*(4+1)+(4))
mean(cp_results_rodent5$saved_lls * -2)+ 2*(3*(ntopics-1)*(5+1)+(5))

# =================================================================
# 5. figures
# =================================================================

# plot community compositions
beta1 = community_composition(ldamodel)
# change column names to modern taxonomy of species
colnames(beta1)[colnames(beta1)=='PB'] <- 'CB'
colnames(beta1)[colnames(beta1)=='PH'] <- 'CH'
colnames(beta1)[colnames(beta1)=='PP'] <- 'CP'
colnames(beta1)[colnames(beta1)=='PI'] <- 'CI'
# put columns in order of largest species to smallest
composition = beta1[,c('NA','DS','SH','SF','SO','DO','DM','CB','CH','OL','OT','PL','PM','PE','CP','CI','RF','RM','RO','BA','PF')]
plot_community_composition(composition,c(3,4,1,2))


# community composition with grassland communities highlighted
P = plot_community_composition_gg(composition,c(3,4,1,2),ylim=c(0,.8))

(figure_spcomp <- multi_panel_figure(
  width = c(38,38,38,38),
  height = c(40,5),
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
H_4 = ggplot(data = df_4, aes(x=value)) +
  geom_histogram(data=subset(df_4,variable=='V1'),aes(y=..count../sum(..count..)),binwidth = .5,fill='black',alpha=.2) +
  geom_histogram(data=subset(df_4,variable=='V2'),aes(y=..count../sum(..count..)),binwidth = .5,fill='black',alpha=.4) +
  geom_histogram(data=subset(df_4,variable=='V3'),aes(y=..count../sum(..count..)),binwidth = .5,fill='black',alpha=.6) +
  geom_histogram(data=subset(df_4,variable=='V4'),aes(y=..count../sum(..count..)),binwidth = .5,fill='black',alpha=1) +
  labs(x='',y='') +
  xlim(range(year_continuous)) +
  scale_y_continuous(labels=c('0.00','0.20','0.40','0.60','0.80'),breaks = c(0,.2,.4,.6,.8)) +
  theme(axis.text=element_text(size=12),
        panel.border=element_rect(colour='black',fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour='grey90'),
        panel.grid.minor = element_line(colour='grey90')) 
  #theme_bw()
H_4

 
# changepoint model plot
cpts = find_changepoint_location(cp_results_rodent4)
cpt_plot = get_ll_non_memoized_plot(ldamodel,x,cpts,make_plot=T,weights=rep(1,length(year_continuous)))


# Figure 1 -- community composition, LDA model, changepoint histogram, changepoint timeseries
(figure <- multi_panel_figure(
  width = c(3,38,38,38,38),
  height = c(5,40,42,42,42),
  column_spacing = 0,
  row_spacing = 0,
  panel_label_type = "lower-alpha"))
figure %<>% fill_panel(
  figure_spcomp,
  row = 2, column = 2:5)
figure %<>% fill_panel(
  cc,
  row = 3, column = 2:5)
figure %<>% fill_panel(
  H_4,
  row = 4, column = 2:5)
figure %<>% fill_panel(
  cpt_plot,
  row = 5, column = 2:5)
figure

save_multi_panel_figure(figure,'Figure1.tiff',dpi=600,compression='lzw')

# ===================================================================
# 6. appendix: LDA with 3 and 5 topics
# ===================================================================

# 3 topics  size = 1100x650
ldamodel3topic = LDA(dat,3, control = list(seed = 46),method='VEM')
cc3 = plot_component_communities(ldamodel3topic,3,dates)
beta13topic = community_composition(ldamodel3topic)
# change column names to modern taxonomy of species
colnames(beta13topic)[colnames(beta13topic)=='PB'] <- 'CB'
colnames(beta13topic)[colnames(beta13topic)=='PH'] <- 'CH'
colnames(beta13topic)[colnames(beta13topic)=='PP'] <- 'CP'
colnames(beta13topic)[colnames(beta13topic)=='PI'] <- 'CI'
composition3 = beta13topic[,c('NA','DS','SH','SF','SO','DO','DM','CB','CH','OL','OT','PL','PM','PE','CP','CI','RF','RM','RO','BA','PF')]
P3topic = plot_community_composition_gg(composition3,c(3,2,1),c(0,.8))

(figure_spcomp3 <- multi_panel_figure(
  width = c(60,60,60),
  height = c(60,10),
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
  width = c(5,70,70,70,70),
  height = c(60,60),
  column_spacing = 0,
  panel_label_type = "upper-alpha"))
figure_s2 %<>% fill_panel(
  figure_spcomp3,
  row = 1, column = 2:5)
figure_s2 %<>% fill_panel(
  cc3,
  row = 2, column = 2:5)
figure_s2
save_multi_panel_figure(figure_s2,'FigureB-3.tiff',dpi=600,compression='lzw')

# 4 topics
(figure_spcomp4 <- multi_panel_figure(
  width = c(65,65,65,65),
  height = c(60,10),
  panel_label_type = "none",
  column_spacing = 0))
figure_spcomp4 %<>% fill_panel(
  P[[1]],
  row = 1, column = 1)
figure_spcomp4 %<>% fill_panel(
  P[[2]],
  row = 1, column = 2)
figure_spcomp4 %<>% fill_panel(
  P[[3]],
  row = 1, column = 3)
figure_spcomp4 %<>% fill_panel(
  P[[4]],
  row = 1, column = 4)
figure_spcomp4

(figure_s3 <- multi_panel_figure(
  width = c(5,65,65,65,65),
  height = c(60,60),
  column_spacing = 0,
  panel_label_type = "upper-alpha"))
figure_s3 %<>% fill_panel(
  figure_spcomp4,
  row = 1, column = 2:5)
figure_s3 %<>% fill_panel(
  cc,
  row = 2, column = 2:5)
figure_s3
save_multi_panel_figure(figure_s3,'FigureB-4.tiff',dpi=600,compression='lzw')

# 5 topics
ldamodel5topic = LDA(dat,5, control = list(seed = 110),method='VEM')
cc5 = plot_component_communities(ldamodel5topic,5,dates,'',c(1,5,3,4,2))
beta15topic = community_composition(ldamodel5topic)
# change column names to modern taxonomy of species
colnames(beta15topic)[colnames(beta15topic)=='PB'] <- 'CB'
colnames(beta15topic)[colnames(beta15topic)=='PH'] <- 'CH'
colnames(beta15topic)[colnames(beta15topic)=='PP'] <- 'CP'
colnames(beta15topic)[colnames(beta15topic)=='PI'] <- 'CI'
composition5 = beta15topic[c(1,5,3,4,2),c('NA','DS','SH','SF','SO','DO','DM','CB','CH','OL','OT','PL','PM','PE','CP','CI','RF','RM','RO','BA','PF')]
P5topic = plot_community_composition_gg(composition5,c(3,4,5,2,1),c(0,.8))

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
  width = c(70,70,70,70,5),
  height = c(60,60),
  panel_label_type = "upper-alpha"))
figure_s4 %<>% fill_panel(
  figure_spcomp5,
  row = 1, column = 1:4)
figure_s4 %<>% fill_panel(
  cc5,
  row = 2, column = 1:4)
figure_s4
save_multi_panel_figure(figure_s4,'FigureB-5.tiff',dpi=600,compression='lzw')

# =======================================================
# two, three, and five changepoints
cols = viridis_pal()(5)

# create dataframes from model outputs
df_2 = as.data.frame(t(cp_results_rodent2$saved[,1,])) %>% melt()
df_2$value = year_continuous[df_2$value]
df_3 = as.data.frame(t(cp_results_rodent3$saved[,1,])) %>% melt()
df_3$value = year_continuous[df_3$value]
df_5 = as.data.frame(t(cp_results_rodent5$saved[,1,])) %>% melt()
df_5$value = year_continuous[df_5$value]

# changepoint histogram
H_2 = ggplot(data = df_2, aes(x=value)) +
  geom_histogram(data=subset(df_2,variable=='V1'),aes(y=..count../sum(..count..)),binwidth = .5,fill=cols[1],alpha=.5,color='black') +
  geom_histogram(data=subset(df_2,variable=='V2'),aes(y=..count../sum(..count..)),binwidth = .5,fill=cols[2],alpha=.5,color='black') +
  labs(x='',y='') +
  scale_y_continuous(labels=c('0.00','0.20','0.40','0.60','0.80'),breaks = c(0,.2,.4,.6,.8)) +
  theme(axis.text=element_text(size=12),
        panel.border=element_rect(colour='black',fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour='grey90'),
        panel.grid.minor = element_line(colour='grey90')) +
  xlim(range(year_continuous))
H_2

H_3 = ggplot(data = df_3, aes(x=value)) +
  geom_histogram(data=subset(df_3,variable=='V1'),aes(y=..count../sum(..count..)),binwidth = .5,fill=cols[1],alpha=.5,color='black') +
  geom_histogram(data=subset(df_3,variable=='V2'),aes(y=..count../sum(..count..)),binwidth = .5,fill=cols[2],alpha=.5,color='black') +
  geom_histogram(data=subset(df_3,variable=='V3'),aes(y=..count../sum(..count..)),binwidth = .5,fill=cols[3],alpha=.5,color='black') +
  labs(x='',y='') +
  scale_y_continuous(labels=c('0.00','0.20','0.40','0.60','0.80'),breaks = c(0,.2,.4,.6,.8)) +
  theme(axis.text=element_text(size=12),
        panel.border=element_rect(colour='black',fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour='grey90'),
        panel.grid.minor = element_line(colour='grey90')) +
  xlim(range(year_continuous))
H_3

H_4b = ggplot(data = df_4, aes(x=value)) +
  geom_histogram(data=subset(df_4,variable=='V1'),aes(y=..count../sum(..count..)),binwidth = .5,fill=cols[1],alpha=.5,color='black') +
  geom_histogram(data=subset(df_4,variable=='V3'),aes(y=..count../sum(..count..)),binwidth = .5,fill=cols[2],alpha=.5,color='black') +
  geom_histogram(data=subset(df_4,variable=='V4'),aes(y=..count../sum(..count..)),binwidth = .5,fill=cols[3],alpha=.5,color='black') +
  geom_histogram(data=subset(df_4,variable=='V2'),aes(y=..count../sum(..count..)),binwidth = .5,fill=cols[4],alpha=.5,color='black') +
  labs(x='',y='') +
  xlim(range(year_continuous)) +
  scale_y_continuous(labels=c('0.00','0.20','0.40','0.60','0.80'),breaks = c(0,.2,.4,.6,.8)) +
  theme(axis.text=element_text(size=12),
        panel.border=element_rect(colour='black',fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour='grey90'),
        panel.grid.minor = element_line(colour='grey90')) 
H_4b

H_5 = ggplot(data = df_5, aes(x=value)) +
  geom_histogram(data=subset(df_5,variable=='V1'),aes(y=..count../sum(..count..)),binwidth = .5,fill=cols[1],alpha=.5,color='black') +
  geom_histogram(data=subset(df_5,variable=='V4'),aes(y=..count../sum(..count..)),binwidth = .5,fill=cols[2],alpha=.5,color='black') +
  geom_histogram(data=subset(df_5,variable=='V5'),aes(y=..count../sum(..count..)),binwidth = .5,fill=cols[3],alpha=.5,color='black') +
  geom_histogram(data=subset(df_5,variable=='V2'),aes(y=..count../sum(..count..)),binwidth = .5,fill=cols[4],alpha=.5,color='black') +
  geom_histogram(data=subset(df_5,variable=='V3'),aes(y=..count../sum(..count..)),binwidth = .5,fill=cols[5],alpha=.5,color='black') +
  labs(x='',y='') +
  scale_y_continuous(labels=c('0.00','0.20','0.40','0.60','0.80'),breaks = c(0,.2,.4,.6,.8)) +
  theme(axis.text=element_text(size=12),
        panel.border=element_rect(colour='black',fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour='grey90'),
        panel.grid.minor = element_line(colour='grey90')) +
  xlim(range(year_continuous))
H_5



# ======================
# stacked histograms
# H_2 = ggplot(data = df_2, aes(x=value)) +
#   geom_histogram(data=df_2,aes(y=..count../sum(..count..),fill=variable),binwidth = .5,color='black') +
#   labs(x='',y='') +
#   scale_y_continuous(labels=c('0.00','0.20','0.40','0.60','0.80'),breaks = c(0,.2,.4,.6,.8)) +
#   theme(axis.text=element_text(size=12),
#         panel.border=element_rect(colour='black',fill=NA),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(colour='grey90'),
#         panel.grid.minor = element_line(colour='grey90'),
#         legend.position = 'none') +
#   xlim(range(year_continuous))
# H_2
# H_3 = ggplot(data = df_3, aes(x=value)) +
#   geom_histogram(data=df_3,aes(y=..count../sum(..count..),fill=variable),binwidth = .5,color='black') +
#   labs(x='',y='') +
#   scale_y_continuous(labels=c('0.00','0.20','0.40','0.60','0.80'),breaks = c(0,.2,.4,.6,.8)) +
#   theme(axis.text=element_text(size=12),
#         panel.border=element_rect(colour='black',fill=NA),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(colour='grey90'),
#         panel.grid.minor = element_line(colour='grey90'),
#         legend.position = 'none') +
#   xlim(range(year_continuous)) 
# H_3
# H_4b = ggplot(data = df_4, aes(x=value)) +
#   geom_histogram(data=df_4,aes(y=..count../sum(..count..),fill=variable),binwidth = .5,color='black') +
#   labs(x='',y='') +
#   scale_y_continuous(labels=c('0.00','0.20','0.40','0.60','0.80'),breaks = c(0,.2,.4,.6,.8)) +
#   theme(axis.text=element_text(size=12),
#         panel.border=element_rect(colour='black',fill=NA),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(colour='grey90'),
#         panel.grid.minor = element_line(colour='grey90'),
#         legend.position = 'none') +
#   xlim(range(year_continuous))
# H_4b
# H_5 = ggplot(data = df_5, aes(x=value)) +
#   geom_histogram(data=df_5,aes(y=..count../sum(..count..),fill=variable),binwidth = .5,color='black') +
#   labs(x='',y='') +
#   scale_y_continuous(labels=c('0.00','0.20','0.40','0.60','0.80','1.00','1.20'),breaks = c(0,.2,.4,.6,.8,1,1.2)) +
#   theme(axis.text=element_text(size=12),
#         panel.border=element_rect(colour='black',fill=NA),
#         panel.background = element_blank(),
#         panel.grid.major = element_line(colour='grey90'),
#         panel.grid.minor = element_line(colour='grey90'),
#         legend.position = 'none') +
#   xlim(range(year_continuous))
# H_5



(figure_s6 <- multi_panel_figure(
  width = c(5,60,60,60,60),
  height = c(60,60,60,60),
  column_spacing = 0,
  panel_label_type = "upper-alpha"))
figure_s6 %<>% fill_panel(
  H_2,
  row = 1, column = 2:5)
figure_s6 %<>% fill_panel(
  H_3,
  row = 2, column = 2:5)
figure_s6 %<>% fill_panel(
  H_4b,
  row = 3, column = 2:5)
figure_s6 %<>% fill_panel(
  H_5,
  row = 4, column = 2:5)
figure_s6
save_multi_panel_figure(figure_s6,'FigureB-2.tiff',dpi=600,compression='lzw')

# ============================================================
# figures not in manuscript
# changepoint model plot
cpts_3 = find_changepoint_location(cp_results_rodent3)
cpt_plot_3pts = get_ll_non_memoized_plot(ldamodel,x,cpts_3,make_plot=T,weights=rep(1,length(year_continuous)))

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
