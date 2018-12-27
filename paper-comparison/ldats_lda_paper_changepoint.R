library(LDATS)
#### Data #### 

data(rodents)

ldats_dat <- rodents[[1]]
ldats_dates <- rodents[[2]]
head(ldats_dat)
source('paper-comparison/christensen-ecology-files/rodent_data_for_LDA.r')
dat = create_rodent_table(period_first = 1,
                          period_last = 436,
                          selected_plots = c(2,4,8,11,12,14,17,22),
                          selected_species = c('BA','DM','DO','DS','NA','OL','OT','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO','SF','SH','SO'), diagnose = F)

# dates to go with count data
moondat = read.csv(text=getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/moon_dates.csv"),stringsAsFactors = F)
moondat$date = as.Date(moondat$censusdate)

period_dates = filter(moondat,period %in% rownames(dat)) %>% select(period,date)
dates = period_dates$date

ch_dat <- dat

ch_dates <- dates

rm(list = c('dat', 'dates', 'period_dates', 'moondat',
            'create_rodent_table'))

ldats_dat <- as.matrix(ch_dat)
ldats_dat <- apply(ldats_dat, c(1,2), FUN = as.integer)
rodents$document_term_table <- ldats_dat
rm(ldats_dat)
rm(ldats_dates)
rm(ch_dat)

#### Load LDA ####
load('/Users/renatadiaz/Documents/model-stash/ldats_lda.Rds')

#### Run changepoint ####
source('paper-comparison/christensen-ecology-files/changepointmodel.r')
# set up parameters for model
year_continuous = 1970 + as.integer(julian(ch_dates)) / 365.25
x = data.frame(
  year_continuous = year_continuous,
  sin_year = sin(year_continuous * 2 * pi),
  cos_year = cos(year_continuous * 2 * pi)
)

# # run models with 1, 2, 3, 4, 5 changepoints
# cp_results_rodent = changepoint_model(ldats_lda_selected[[1]], x, 1, weights = rep(1,length(year_continuous)))
# save(cp_results_rodent, file = '/Users/renatadiaz/Documents/model-stash/paper_cpt_ldats_lda1.Rds')
# rm(cp_results_rodent)
# 
# cp_results_rodent2 = changepoint_model(ldats_lda_selected[[1]], x, 2, weights = rep(1,length(year_continuous)))
# save(cp_results_rodent2, file = '/Users/renatadiaz/Documents/model-stash/paper_cpt_ldats_lda2.Rds')
# rm(cp_results_rodent2)
# 
# cp_results_rodent3 = changepoint_model(ldats_lda_selected[[1]], x, 3, weights = rep(1,length(year_continuous)))
# save(cp_results_rodent3, file = '/Users/renatadiaz/Documents/model-stash/paper_cpt_ldats_lda3.Rds')
# rm(cp_results_rodent3)
# 
# cp_results_rodent4 = changepoint_model(ldats_lda_selected[[1]], x, 4, weights = rep(1,length(year_continuous)))
# save(cp_results_rodent4, file = '/Users/renatadiaz/Documents/model-stash/paper_cpt_ldats_lda4.Rds')
# rm(cp_results_rodent4)
# 
# cp_results_rodent5 = changepoint_model(ldats_lda_selected[[1]], x, 5, weights = rep(1,length(year_continuous)))
# save(cp_results_rodent5, file = '/Users/renatadiaz/Documents/model-stash/paper_cpt_ldats_lda5.Rds')
# rm(cp_results_rodent5)
# 
# cp_results_rodent6 = changepoint_model(ldats_lda_selected[[1]], x, 6, weights = rep(1,length(year_continuous)))
# save(cp_results_rodent6, file = '/Users/renatadiaz/Documents/model-stash/paper_cpt_ldats_lda6.Rds')
# rm(cp_results_rodent6)
# 

#### Changepoint model selection ####
load('/Users/renatadiaz/Documents/model-stash/paper_cpt_ldats_lda1.Rds')
load('/Users/renatadiaz/Documents/model-stash/paper_cpt_ldats_lda2.Rds')
load('/Users/renatadiaz/Documents/model-stash/paper_cpt_ldats_lda3.Rds')
load('/Users/renatadiaz/Documents/model-stash/paper_cpt_ldats_lda4.Rds')
load('/Users/renatadiaz/Documents/model-stash/paper_cpt_ldats_lda5.Rds')
load('/Users/renatadiaz/Documents/model-stash/paper_cpt_ldats_lda6.Rds')

ntopics = ldats_lda_selected[[1]]@k

# change point model selection
# mean deviance ( -2 * log likelihood) + 2*(#parameters)
mean(cp_results_rodent$saved_lls * -2) + 2*(3*(ntopics-1)*(1+1)+(1))
mean(cp_results_rodent2$saved_lls * -2)+ 2*(3*(ntopics-1)*(2+1)+(2))
mean(cp_results_rodent3$saved_lls * -2)+ 2*(3*(ntopics-1)*(3+1)+(3))
mean(cp_results_rodent4$saved_lls * -2)+ 2*(3*(ntopics-1)*(4+1)+(4))
mean(cp_results_rodent5$saved_lls * -2)+ 2*(3*(ntopics-1)*(5+1)+(5))
mean(cp_results_rodent6$saved_lls * -2)+ 2*(3*(ntopics-1)*(5+1)+(5))

# lowest deviance is for 4 (or 6) changepoints

rm(list = c('cp_results_rodent', 'cp_results_rodent2', 
            'cp_results_rodent3', 'cp_results_rodent5'))
hist(year_continuous[cp_results_rodent4$saved[,1,]],breaks = seq(1977,2016,.25),xlab='',main='Changepoint Estimate')
hist(year_continuous[cp_results_rodent6$saved[,1,]],breaks = seq(1977,2016,.25),xlab='',main='Changepoint Estimate')

annual_hist(cp_results_rodent4,year_continuous)

# turn changepoint results into data frame
df_4 = as.data.frame(t(cp_results_rodent4$saved[,1,])) %>% melt()
df_4$value = year_continuous[df_4$value]

# find 95% confidence intervals on each changepoint:
quantile(df_4[df_4$variable=='V1','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
quantile(df_4[df_4$variable=='V2','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
quantile(df_4[df_4$variable=='V3','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
quantile(df_4[df_4$variable=='V4','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')

## 6 cpts 
annual_hist(cp_results_rodent6,year_continuous)
# turn changepoint results into data frame
df_6 = as.data.frame(t(cp_results_rodent6$saved[,1,])) %>% melt()
df_6$value = year_continuous[df_6$value]

# find 95% confidence intervals on each changepoint:
quantile(df_6[df_6$variable=='V1','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
quantile(df_6[df_6$variable=='V2','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
quantile(df_6[df_6$variable=='V3','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
quantile(df_6[df_6$variable=='V4','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
quantile(df_6[df_6$variable=='V5','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
quantile(df_6[df_6$variable=='V6','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
