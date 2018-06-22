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
library(portalr)
library(dplyr)
library(tidyr)
#source('rodent_data_for_LDA.r')
source('RodentAbundancesAdjustable.R')
# try this
# source('/Users/renatadiaz/Documents/GitHub/portalr/R/RodentAbundancesAdjustable.R') # not sure how to source the thing I added to portalr
source('AIC_model_selection.R')
source('LDA_figure_scripts.R')
source('changepointmodel.r')
source('LDA-distance.R')


# to skip running the changepoint and just load results
#source('readResults.R')

# ===================================================================
# 1. prepare rodent data -- needs to be adjusted to work with new changes to portalr
# ===================================================================
# 
# dat.full <- abundance.adjustable(path = 'repo', level="treatment.adj",type="Rodents",
#                                              length="longterm",unknowns=F,incomplete=T,
#                                              shape="list",time="period", dates = T)
# 
# 
# 
# period_first = 1
# period_last = 436
# selected_treatment = 'exclosure'
# 
# dat = dat.full %>% 
#   filter(period >= period_first, period <= period_last, 
#                 treatment == selected_treatment) %>% 
#   mutate(usual.n = (ceiling(mean(n)))) %>% 
#   mutate(abundance.adj = round(abundance.perplot * usual.n)) %>% 
#   #select(species, abundance.adj, period, treatment, n, usual.n, censusdate) %>% 
#   select(species, period, abundance.adj, censusdate) %>% 
#   spread(species, abundance.adj) 
# 
# 
# dates = select(dat, censusdate)
# dates = as.Date(dates$censusdate)
# 
# dates.periods = select(dat, censusdate, period)
# dates.periods$censusdate = as.Date(dates.periods$censusdate)
# 
# 
# dat = dat[,3:23]
# write.csv(dat, 'exclosures_dat.csv', row.names = FALSE)
# 
# 
# write.csv(dates.periods, 'exclosures_dates.csv', row.names = FALSE)
# 

dat = read.csv('exclosures_dat.csv')
dates = read.csv('exclosures_dates.csv')
dates$censusdate = as.Date(dates$censusdate)
# ==================================================================
# 2a. select number of topics
# ==================================================================

# Fit a bunch of LDA models with different seeds
# Only use even numbers for seeds because consecutive seeds give identical results
seeds = 2*seq(20)

# repeat LDA model fit and AIC calculation with a bunch of different seeds to test robustness of the analysis
best_ntopic = repeat_VEM(dat,
                         seeds,
                         topic_min=2,
                         topic_max=6)

# histogram of how many seeds chose how many topics
hist(best_ntopic$k,breaks=c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5),xlab='best # of topics', main='')

# with exclosures, 3 topics wins by far. 


# ==================================================================
# 2b. how different is species composition of ##3## community-types when LDA is run with different seeds?
# ==================================================================
# get the best 100 seeds where 4 topics was the best LDA model
seeds_3topics = best_ntopic %>% 
  filter(k == 3) %>% 
  arrange(aic) %>% 
  head(100) %>% 
  pull(SEED)

# choose seed with highest log likelihood for all following analyses
#    (also produces plot of community composition for 'best' run compared to 'worst')
best_seed = calculate_LDA_distance(dat,seeds_3topics,3)
mean_dist = unlist(best_seed)[2]
max_dist = unlist(best_seed)[3]

# ==================================================================
# 3. run LDA model
# ==================================================================
ntopics = 3
SEED = unlist(best_seed)[1]
ldamodel = LDA(dat,ntopics, control = list(seed = SEED),method='VEM')


# ==================================================================
# 4. change point model 
# ==================================================================

# set up parameters for model
year_continuous = 1970 + as.integer(julian(dates$censusdate)) / 365.25
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
hist(year_continuous[cp_results_rodent3$saved[,1,]],breaks = seq(1977,2016,.25),xlab='',main='Changepoint Estimate')
annual_hist(cp_results_rodent3,year_continuous)

# turn changepoint results into data frame

df_1 = as.data.frame(t(cp_results_rodent$saved[,1,])) %>% melt()
df_1$value = year_continuous[df_1$value]

df_2 = as.data.frame(t(cp_results_rodent2$saved[,1,])) %>% melt()
df_2$value = year_continuous[df_2$value]

df_3 = as.data.frame(t(cp_results_rodent3$saved[,1,])) %>% melt()
df_3$value = year_continuous[df_3$value]

df_4 = as.data.frame(t(cp_results_rodent4$saved[,1,])) %>% melt()
df_4$value = year_continuous[df_4$value]

df_5 = as.data.frame(t(cp_results_rodent5$saved[,1,])) %>% melt()
df_5$value = year_continuous[df_5$value]

# saving inputs and results because this took a couple of days to run 
saveThese = list(dat, seeds, best_ntopic, seeds_3topics, best_seed, mean_dist, max_dist, ntopics, SEED, 
                 ldamodel, year_continuous, x, cp_results_rodent, cp_results_rodent2, cp_results_rodent3,
                 cp_results_rodent4, cp_results_rodent5, df_1, df_2, df_3, df_4, df_5, dates)
saveRDS(saveThese, 'exclosuresResults.rds')

source('readResults.R')

# find 95% confidence intervals on each changepoint:
quantile(df_3[df_3$variable=='V1','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
quantile(df_3[df_3$variable=='V2','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
quantile(df_3[df_3$variable=='V3','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
# 
# > quantile(df_3[df_3$variable=='V1','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
# 2.5%        97.5% 
#   "22-05-1990" "03-02-1995" 
# > quantile(df_3[df_3$variable=='V2','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
# 2.5%        97.5% 
#   "08-02-1997" "22-11-1997" 
# > quantile(df_3[df_3$variable=='V3','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
# 2.5%        97.5% 
#   "19-09-2009" "14-05-2010" 


# do this before finding the confidence intervals?
# change point model selection
# mean deviance ( -2 * log likelihood) + 2*(#parameters)
mean(cp_results_rodent$saved_lls * -2) + 2*(3*(ntopics-1)*(1+1)+(1))
mean(cp_results_rodent2$saved_lls * -2)+ 2*(3*(ntopics-1)*(2+1)+(2))
mean(cp_results_rodent3$saved_lls * -2)+ 2*(3*(ntopics-1)*(3+1)+(3))
mean(cp_results_rodent4$saved_lls * -2)+ 2*(3*(ntopics-1)*(4+1)+(4))
mean(cp_results_rodent5$saved_lls * -2)+ 2*(3*(ntopics-1)*(5+1)+(5))

# lowest deviance w/3 change points
# > mean(cp_results_rodent$saved_lls * -2) + 2*(3*(ntopics-1)*(1+1)+(1))
# [1] 701.7584
# > mean(cp_results_rodent2$saved_lls * -2)+ 2*(3*(ntopics-1)*(2+1)+(2))
# [1] 671.2621
# > mean(cp_results_rodent3$saved_lls * -2)+ 2*(3*(ntopics-1)*(3+1)+(3))
# [1] 662.1949
# > mean(cp_results_rodent4$saved_lls * -2)+ 2*(3*(ntopics-1)*(4+1)+(4))
# [1] 664.2978
# > mean(cp_results_rodent5$saved_lls * -2)+ 2*(3*(ntopics-1)*(5+1)+(5))
# [1] 666.7015

# =================================================================
# 5. figures
# =================================================================
library(cowplot)
# plot community compositions
beta1 = community_composition(ldamodel)
# put columns in order of largest species to smallest
composition = beta1[,c('NA','DS','SH','SF','SO','DO','DM','PB','PH','OL','OT','PL','PM','PE','PP','PI','RF','RM','RO','BA','PF')]
comp_plots = plot_community_composition_gg(composition,c(3,1,2), ylim=c(0,1), grassland = FALSE)


# community composition with grassland communities highlighted
comp_plots_grassland = plot_community_composition_gg(composition,c(3,1,2),ylim=c(0,1), grassland = TRUE)


for (i in 1:3) {
setEPS()
postscript(paste0("exclosure_community", i, ".eps"))
print(comp_plots[[i]])
dev.off()

setEPS()
postscript(paste0("exclosure_community", i, "_grassland.eps"))
print(comp_plots_grassland[[i]])
dev.off()

}

together <- plot_grid(plotlist = comp_plots, align = "v", axis = 'l', nrow = 3, ncol = 1)
setEPS()
postscript("exclosure_community_all.eps")
print(together)
dev.off()



(spcomp3= multi_panel_figure(
  width = c(70,70,70),
  height = c(70,10),
  panel_label_type = "none",
  column_spacing = 0))
spcomp3 %<>% fill_panel(
  comp_plots[[1]],
  row = 1, column = 1)
spcomp3 %<>% fill_panel(
  comp_plots[[2]],
  row = 1, column = 2)
spcomp3 %<>% fill_panel(
  comp_plots[[3]],
  row = 1, column = 3)
spcomp3

# plot of component communities over time
cc = plot_component_communities(ldamodel,ntopics,dates)
cc

setEPS()
postscript("component_communities_time.eps")
print(cc)
dev.off()

# changepoint histogram w 3 cpts
H_3 = ggplot(data = df_3, aes(x=value)) +
  geom_histogram(data=subset(df_3,variable=='V1'),aes(y=..count../sum(..count..)),binwidth = .5,fill='grey') + #,alpha=.2) +
  geom_histogram(data=subset(df_3,variable=='V2'),aes(y=..count../sum(..count..)),binwidth = .5,fill='darkgrey') + #,alpha=.4) +
  geom_histogram(data=subset(df_3,variable=='V3'),aes(y=..count../sum(..count..)),binwidth = .5,fill='black') + #,alpha=.6) +
  labs(x='',y='') +
  xlim(range(year_continuous)) +
  scale_y_continuous(labels=c('0.00','0.20','0.40','0.60','0.80'),breaks = c(0,.2,.4,.6,.8)) +
  theme(axis.text=element_text(size=12),
        panel.border=element_rect(colour='black',fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour='grey90'),
        panel.grid.minor = element_line(colour='grey90'))+  
  theme_bw()
H_3

setEPS()
postscript("changepoint_histogram.eps")
print(H_3)
dev.off()

 
# changepoint model plot
cpts = find_changepoint_location(cp_results_rodent3)
cpt_plot = get_ll_non_memoized_plot(ldamodel,x,cpts,make_plot=T,weights=rep(1,length(year_continuous)))

setEPS()
postscript('changepoint_plot.eps')
print(cpt_plot)
dev.off()

cpt_plot
# Figure 3 -- community composition, LDA model, changepoint histogram, changepoint timeseries
(figure <- multi_panel_figure(
  width = c(70,70,70,70),
  height = c(60,60,60,60),
  column_spacing = 0))
figure %<>% fill_panel(
  spcomp3,
  row = 1, column = 1:4)
figure %<>% fill_panel(
  cc,
  row = 2, column = 1:4)
figure %<>% fill_panel(
  H_3,
  row = 3, column = 1:4)
figure %<>% fill_panel(
  cpt_plot,
  row = 4, column = 1:4)
figure

setEPS()
postscript('figure3_exclosures.eps', width = 12, height = 12)
print(figure)
dev.off()


## saving pdfs

for (i in 1:3) {
  pdf(paste0("exclosure_community", i, ".pdf"))
  print(comp_plots[[i]])
  dev.off()
  
  pdf(paste0("exclosure_community", i, "_grassland.pdf"))
  print(comp_plots_grassland[[i]])
  dev.off()
  
}

together <- plot_grid(plotlist = comp_plots, align = "v", axis = 'l', nrow = 3, ncol = 1)
pdf("exclosure_community_all.pdf")
print(together)
dev.off()



(spcomp3= multi_panel_figure(
  width = c(70,70,70),
  height = c(70,10),
  panel_label_type = "none",
  column_spacing = 0))
spcomp3 %<>% fill_panel(
  comp_plots[[1]],
  row = 1, column = 1)
spcomp3 %<>% fill_panel(
  comp_plots[[2]],
  row = 1, column = 2)
spcomp3 %<>% fill_panel(
  comp_plots[[3]],
  row = 1, column = 3)
spcomp3

# plot of component communities over time
cc = plot_component_communities(ldamodel,ntopics,dates)
cc


pdf("component_communities_time.pdf")
print(cc)
dev.off()

# changepoint histogram w 3 cpts
H_3 = ggplot(data = df_3, aes(x=value)) +
  geom_histogram(data=subset(df_3,variable=='V1'),aes(y=..count../sum(..count..)),binwidth = .5,fill='grey') + #,alpha=.2) +
  geom_histogram(data=subset(df_3,variable=='V2'),aes(y=..count../sum(..count..)),binwidth = .5,fill='darkgrey') + #,alpha=.4) +
  geom_histogram(data=subset(df_3,variable=='V3'),aes(y=..count../sum(..count..)),binwidth = .5,fill='black') + #,alpha=.6) +
  labs(x='',y='') +
  xlim(range(year_continuous)) +
  scale_y_continuous(labels=c('0.00','0.20','0.40','0.60','0.80'),breaks = c(0,.2,.4,.6,.8)) +
  theme(axis.text=element_text(size=12),
        panel.border=element_rect(colour='black',fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour='grey90'),
        panel.grid.minor = element_line(colour='grey90'))+  
  theme_bw()
H_3


pdf("changepoint_histogram.pdf")
print(H_3)
dev.off()


# changepoint model plot
cpts = find_changepoint_location(cp_results_rodent3)
cpt_plot = get_ll_non_memoized_plot(ldamodel,x,cpts,make_plot=T,weights=rep(1,length(year_continuous)))


pdf('changepoint_plot.pdf')
print(cpt_plot)
dev.off()

cpt_plot
# Figure 3 -- community composition, LDA model, changepoint histogram, changepoint timeseries
(figure <- multi_panel_figure(
  width = c(70,70,70,70),
  height = c(60,60,60,60),
  column_spacing = 0))
figure %<>% fill_panel(
  spcomp3,
  row = 1, column = 1:4)
figure %<>% fill_panel(
  cc,
  row = 2, column = 1:4)
figure %<>% fill_panel(
  H_3,
  row = 3, column = 1:4)
figure %<>% fill_panel(
  cpt_plot,
  row = 4, column = 1:4)
figure


pdf('figure3_exclosures.pdf', width = 12, height = 12)
print(figure)
dev.off()


# ===================================================================
# 6. appendix: LDA with 3 and 5 topics
# ===================================================================

# 3 topics
ldamodel3topic = LDA(dat,3, control = list(seed = 46),method='VEM')
cc3 = plot_component_communities(ldamodel3topic,3,dates)
beta13topic = community_composition(ldamodel3topic)
composition3 = beta13topic[,c('NA','DS','SH','SF','SO','DO','DM','PB','PH','OL','OT','PL','PM','PE','PP','PI','RF','RM','RO','BA','PF')]
P3topic = plot_community_composition_gg(composition3,c(3,2,1),c(0,.8))

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
cc5 = plot_component_communities(ldamodel5topic,5,dates,'',c(1,5,3,4,2))
beta15topic = community_composition(ldamodel5topic)
composition5 = beta15topic[c(1,5,3,4,2),c('NA','DS','SH','SF','SO','DO','DM','PB','PH','OL','OT','PL','PM','PE','PP','PI','RF','RM','RO','BA','PF')]
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

# =======================================================
# two, three, and five changepoints
cols = viridis_pal()(5)

# create dataframes from model outputs
df_2 = as.data.frame(t(cp_results_rodent2$saved[,1,])) %>% melt()
df_2$value = year_continuous[df_2$value]
df_4 = as.data.frame(t(cp_results_rodent4$saved[,1,])) %>% melt()
df_4$value = year_continuous[df_4$value]
df_5 = as.data.frame(t(cp_results_rodent5$saved[,1,])) %>% melt()
df_5$value = year_continuous[df_5$value]

# changepoint histogram
H_2 = ggplot(data = df_2, aes(x=value)) +
  geom_histogram(data=subset(df_2,variable=='V1'),aes(y=..count../sum(..count..)),binwidth = .5,fill=cols[1],alpha=.5) +
  geom_histogram(data=subset(df_2,variable=='V2'),aes(y=..count../sum(..count..)),binwidth = .5,fill=cols[2],alpha=.5) +
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
  geom_histogram(data=subset(df_3,variable=='V1'),aes(y=..count../sum(..count..)),binwidth = .5,fill=cols[1],alpha=.5) +
  geom_histogram(data=subset(df_3,variable=='V2'),aes(y=..count../sum(..count..)),binwidth = .5,fill=cols[2],alpha=.5) +
  geom_histogram(data=subset(df_3,variable=='V3'),aes(y=..count../sum(..count..)),binwidth = .5,fill=cols[3],alpha=.5) +
  labs(x='',y='') +
  scale_y_continuous(labels=c('0.00','0.20','0.40','0.60','0.80'),breaks = c(0,.2,.4,.6,.8)) +
  theme(axis.text=element_text(size=12),
        panel.border=element_rect(colour='black',fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour='grey90'),
        panel.grid.minor = element_line(colour='grey90')) +
  xlim(range(year_continuous))
H_3

H_4 = ggplot(data = df_4, aes(x=value)) +
  geom_histogram(data=subset(df_4,variable=='V1'),aes(y=..count../sum(..count..)),binwidth = .5,alpha=.5) +
  geom_histogram(data=subset(df_4,variable=='V3'),aes(y=..count../sum(..count..)),binwidth = .5,alpha=.5) +
  geom_histogram(data=subset(df_4,variable=='V4'),aes(y=..count../sum(..count..)),binwidth = .5,alpha=.5) +
  geom_histogram(data=subset(df_4,variable=='V2'),aes(y=..count../sum(..count..)),binwidth = .5,alpha=.5) +
  labs(x='',y='') +
  xlim(range(year_continuous)) +
  scale_y_continuous(labels=c('0.00','0.20','0.40','0.60','0.80'),breaks = c(0,.2,.4,.6,.8)) +
  theme(axis.text=element_text(size=12),
        panel.border=element_rect(colour='black',fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour='grey90'),
        panel.grid.minor = element_line(colour='grey90')) 
H_4

H_5 = ggplot(data = df_5, aes(x=value)) +
  geom_histogram(data=subset(df_5,variable=='V1'),aes(y=..count../sum(..count..)),binwidth = .5,fill=cols[1],alpha=.5) +
  geom_histogram(data=subset(df_5,variable=='V4'),aes(y=..count../sum(..count..)),binwidth = .5,fill=cols[2],alpha=.5) +
  geom_histogram(data=subset(df_5,variable=='V5'),aes(y=..count../sum(..count..)),binwidth = .5,fill=cols[3],alpha=.5) +
  geom_histogram(data=subset(df_5,variable=='V2'),aes(y=..count../sum(..count..)),binwidth = .5,fill=cols[4],alpha=.5) +
  geom_histogram(data=subset(df_5,variable=='V3'),aes(y=..count../sum(..count..)),binwidth = .5,fill=cols[5],alpha=.5) +
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
H_2 = ggplot(data = df_2, aes(x=value)) +
  geom_histogram(data=df_2,aes(y=..count../sum(..count..),fill=variable),binwidth = .5,color='black') +
  labs(x='',y='') +
  scale_y_continuous(labels=c('0.00','0.20','0.40','0.60','0.80'),breaks = c(0,.2,.4,.6,.8)) +
  theme(axis.text=element_text(size=12),
        panel.border=element_rect(colour='black',fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour='grey90'),
        panel.grid.minor = element_line(colour='grey90'),
        legend.position = 'none') +
  xlim(range(year_continuous))
H_2
H_3 = ggplot(data = df_3, aes(x=value)) +
  geom_histogram(data=df_3,aes(y=..count../sum(..count..),fill=variable),binwidth = .5,color='black') +
  labs(x='',y='') +
  scale_y_continuous(labels=c('0.00','0.20','0.40','0.60','0.80'),breaks = c(0,.2,.4,.6,.8)) +
  theme(axis.text=element_text(size=12),
        panel.border=element_rect(colour='black',fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour='grey90'),
        panel.grid.minor = element_line(colour='grey90'),
        legend.position = 'none') +
  xlim(range(year_continuous))
H_3
H_4b = ggplot(data = df_4, aes(x=value)) +
  geom_histogram(data=df_4,aes(y=..count../sum(..count..),fill=variable),binwidth = .5,color='black') +
  labs(x='',y='') +
  scale_y_continuous(labels=c('0.00','0.20','0.40','0.60','0.80'),breaks = c(0,.2,.4,.6,.8)) +
  theme(axis.text=element_text(size=12),
        panel.border=element_rect(colour='black',fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour='grey90'),
        panel.grid.minor = element_line(colour='grey90'),
        legend.position = 'none') +
  xlim(range(year_continuous))
H_4b
H_5 = ggplot(data = df_5, aes(x=value)) +
  geom_histogram(data=df_5,aes(y=..count../sum(..count..),fill=variable),binwidth = .5,color='black') +
  labs(x='',y='') +
  scale_y_continuous(labels=c('0.00','0.20','0.40','0.60','0.80','1.00','1.20'),breaks = c(0,.2,.4,.6,.8,1,1.2)) +
  theme(axis.text=element_text(size=12),
        panel.border=element_rect(colour='black',fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour='grey90'),
        panel.grid.minor = element_line(colour='grey90'),
        legend.position = 'none') +
  xlim(range(year_continuous))
H_5



(figure_s6 <- multi_panel_figure(
  width = c(60,60,60,60),
  height = c(60,60,60,60),
  column_spacing = 0))
figure_s6 %<>% fill_panel(
  H_2,
  row = 1, column = 1:4)
figure_s6 %<>% fill_panel(
  H_3,
  row = 2, column = 1:4)
figure_s6 %<>% fill_panel(
  H_4b,
  row = 3, column = 1:4)
figure_s6 %<>% fill_panel(
  H_5,
  row = 4, column = 1:4)
figure_s6


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
