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

period_dates = filter(moondat,period %in% rownames(dat)) %>% select(newmoonnumber, period,date)
dates = period_dates$date

ch_dat <- dat

ch_dates <- dates
ch_newmoons <- period_dates$newmoonnumber

rm(list = c('dat', 'dates', 'period_dates', 'moondat',
            'create_rodent_table'))

ldats_dat <- as.matrix(ch_dat)
ldats_dat <- apply(ldats_dat, c(1,2), FUN = as.integer)
rodents$document_term_table <- ldats_dat
rm(ldats_dat)
rm(ldats_dates)
rm(ch_dat)

#### Load LDA ####
load('/~/Dropbox/ldats-models/ldats_lda.Rds')
year_continuous = 1970 + as.integer(julian(ch_dates)) / 365.25
x = data.frame(
  year_continuous = year_continuous,
  sin_year = sin(year_continuous * 2 * pi),
  cos_year = cos(year_continuous * 2 * pi)
)

x$newmoon <- ch_newmoons
x <- select(x, newmoon, sin_year, cos_year)

#### Run LDATS changepoint ####

ldats_ldats_cpt <- TS_on_LDA(LDA_models = ldats_lda_selected, 
                             document_covariate_table = x,
                             formulas = ~ sin_year + cos_year,
                             nchangepoints = 1:6,
                             weights = NULL)


save(ldats_ldats_cpt, file = '/~/Dropbox/ldats-models/ldats_lda_ldats_cpt.Rds')
rm(ldats_ldats_cpt)
rm(ldats_lda_selected)

#### Load paper LDA ####
load('/~/Dropbox/ldats-models/paper_lda.Rds')
paper_ldats_cpt <- TS_on_LDA(LDA_models = ldamodel, 
                             document_covariate_table = x,
                             formulas = ~ sin_year + cos_year,
                             nchangepoints = 1:6,
                             weights = NULL)
save(paper_ldats_cpt, file = '/~/Dropbox/ldats-models/paper_lda_ldats_cpt.Rds')
