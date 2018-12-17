library(LDATS)
data(rodents)

#### swap raw data for adjusted data (to match paper) ####
you can skip this

source('/Users/renatadiaz/Documents/GitHub/weecology/LDATS/paper-comparison/christensen-ecology-files/rodent_data_for_LDA.r')
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

#### load paper LDA ####

load('/Users/renatadiaz/Documents/GitHub/weecology/LDATS/paper-comparison/paper_lda.Rds')


#### try changepoint ####

ldats_cpt_paper_lda <- LDATS::TS_on_LDA(LDA_models = ldamodel, rodents[[2]],
                                        formulas =c( ~ sin_year + cos_year), nchangepoints = c(2, 3, 4, 5, 6))

save.image(file = 'ldats_cpt_on_paper_lda.Rdata')

# 'Error in LDA_models[[1]] : this S4 class is not subsettable'

ldamodel_list <- list(ldamodel)

ldats_cpt_paper_lda <- LDATS::TS_on_LDA(LDA_models = ldamodel_list, rodents[[2]],
                                        formulas =c( ~ sin_year + cos_year), nchangepoints = c(2, 3, 4, 5, 6))

# 'Error in check_LDA_models(LDA_models) : 
#   LDA_models is not an LDA object or LDA_set object'
