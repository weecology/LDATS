Paper changepoint
================
Renata Diaz
12/26/2018

``` r
library(LDATS)
#### Data #### 

data(rodents)

ldats_dat <- rodents[[1]]
ldats_dates <- rodents[[2]]
head(ldats_dat)
```

    ##   BA DM DO DS NA. OL OT PB PE PF PH PI PL PM PP RF RM RO SF SH SO
    ## 1  0 13  0  2   2  0  0  0  1  1  0  0  0  0  3  0  0  0  0  0  0
    ## 2  0 20  1  3   2  0  0  0  0  4  0  0  0  0  2  0  0  0  0  0  0
    ## 3  0 21  0  8   4  0  0  0  1  2  0  0  0  0  1  0  0  0  0  0  0
    ## 4  0 21  3 12   4  2  3  0  1  1  0  0  0  0  0  0  0  0  0  0  0
    ## 5  0 16  1  9   5  2  1  0  0  2  0  0  0  0  0  0  1  0  0  0  0
    ## 6  0 17  1 13   5  1  5  0  0  3  0  0  0  0  0  0  0  0  0  0  0

``` r
source('christensen-ecology-files/rodent_data_for_LDA.r')
```

    ## Loading required package: bitops

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
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
source('christensen-ecology-files/changepointmodel.r')
```

    ## 
    ## Attaching package: 'lubridate'

    ## The following object is masked from 'package:base':
    ## 
    ##     date

    ## Loading required package: viridisLite

``` r
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
```

    ## [1] 1342.239

``` r
mean(cp_results_rodent2$saved_lls * -2)+ 2*(3*(ntopics-1)*(2+1)+(2))
```

    ## [1] 1306.275

``` r
mean(cp_results_rodent3$saved_lls * -2)+ 2*(3*(ntopics-1)*(3+1)+(3))
```

    ## [1] 1283.039

``` r
mean(cp_results_rodent4$saved_lls * -2)+ 2*(3*(ntopics-1)*(4+1)+(4))
```

    ## [1] 1274.954

``` r
mean(cp_results_rodent5$saved_lls * -2)+ 2*(3*(ntopics-1)*(5+1)+(5))
```

    ## [1] 1287.404

``` r
mean(cp_results_rodent6$saved_lls * -2)+ 2*(3*(ntopics-1)*(5+1)+(5))
```

    ## [1] 1270.919

``` r
# lowest deviance is for 4 (or 6) changepoints

rm(list = c('cp_results_rodent', 'cp_results_rodent2', 
            'cp_results_rodent3', 'cp_results_rodent5'))
hist(year_continuous[cp_results_rodent4$saved[,1,]],breaks = seq(1977,2016,.25),xlab='',main='Changepoint Estimate')
```

![](paper_changepoint_files/figure-markdown_github/ldats%20lda-1.png)

``` r
hist(year_continuous[cp_results_rodent6$saved[,1,]],breaks = seq(1977,2016,.25),xlab='',main='Changepoint Estimate')
```

![](paper_changepoint_files/figure-markdown_github/ldats%20lda-2.png)

``` r
annual_hist(cp_results_rodent4,year_continuous)
```

![](paper_changepoint_files/figure-markdown_github/ldats%20lda-3.png)

``` r
# turn changepoint results into data frame
df_4 = as.data.frame(t(cp_results_rodent4$saved[,1,])) %>% melt()
```

    ## No id variables; using all as measure variables

``` r
df_4$value = year_continuous[df_4$value]

# find means on each changepoint:
mean(df_4[df_4$variable=='V1','value']) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ## [1] "12-01-1984"

``` r
mean(df_4[df_4$variable=='V2','value']) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ## [1] "02-06-1991"

``` r
mean(df_4[df_4$variable=='V3','value']) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ## [1] "04-01-1999"

``` r
mean(df_4[df_4$variable=='V4','value']) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ## [1] "10-11-2009"

``` r
# find 95% confidence intervals on each changepoint:
quantile(df_4[df_4$variable=='V1','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ##         2.5%        97.5% 
    ## "13-05-1983" "30-07-1984"

``` r
quantile(df_4[df_4$variable=='V2','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ##         2.5%        97.5% 
    ## "14-12-1990" "09-10-1991"

``` r
quantile(df_4[df_4$variable=='V3','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ##         2.5%        97.5% 
    ## "22-11-1997" "08-10-1999"

``` r
quantile(df_4[df_4$variable=='V4','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ##         2.5%        97.5% 
    ## "23-05-2009" "14-05-2010"

``` r
## 6 cpts 
annual_hist(cp_results_rodent6,year_continuous)
```

![](paper_changepoint_files/figure-markdown_github/ldats%20lda-4.png)

``` r
# turn changepoint results into data frame
df_6 = as.data.frame(t(cp_results_rodent6$saved[,1,])) %>% melt()
```

    ## No id variables; using all as measure variables

``` r
df_6$value = year_continuous[df_6$value]

# find 95% confidence intervals on each changepoint:
quantile(df_6[df_6$variable=='V1','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ##         2.5%        97.5% 
    ## "25-02-1983" "02-07-1984"

``` r
quantile(df_6[df_6$variable=='V2','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ##         2.5%        97.5% 
    ## "30-01-1987" "06-08-1991"

``` r
quantile(df_6[df_6$variable=='V3','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ##         2.5%        97.5% 
    ## "12-03-1991" "21-12-1998"

``` r
quantile(df_6[df_6$variable=='V4','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ##         2.5%        97.5% 
    ## "14-07-1993" "05-03-2003"

``` r
quantile(df_6[df_6$variable=='V5','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ##         2.5%        97.5% 
    ## "17-07-1998" "12-12-2009"

``` r
quantile(df_6[df_6$variable=='V6','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ##         2.5%        97.5% 
    ## "20-06-2009" "25-04-2014"

``` r
rm(list = ls())
```

``` r
library(LDATS)
#### Data #### 

data(rodents)

ldats_dat <- rodents[[1]]
ldats_dates <- rodents[[2]]
head(ldats_dat)
```

    ##   BA DM DO DS NA. OL OT PB PE PF PH PI PL PM PP RF RM RO SF SH SO
    ## 1  0 13  0  2   2  0  0  0  1  1  0  0  0  0  3  0  0  0  0  0  0
    ## 2  0 20  1  3   2  0  0  0  0  4  0  0  0  0  2  0  0  0  0  0  0
    ## 3  0 21  0  8   4  0  0  0  1  2  0  0  0  0  1  0  0  0  0  0  0
    ## 4  0 21  3 12   4  2  3  0  1  1  0  0  0  0  0  0  0  0  0  0  0
    ## 5  0 16  1  9   5  2  1  0  0  2  0  0  0  0  0  0  1  0  0  0  0
    ## 6  0 17  1 13   5  1  5  0  0  3  0  0  0  0  0  0  0  0  0  0  0

``` r
source('christensen-ecology-files/rodent_data_for_LDA.r')
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
load('/Users/renatadiaz/Documents/model-stash/paper_lda.Rds')

#### Run changepoint ####
source('christensen-ecology-files/changepointmodel.r')
# set up parameters for model
year_continuous = 1970 + as.integer(julian(ch_dates)) / 365.25
x = data.frame(
  year_continuous = year_continuous,
  sin_year = sin(year_continuous * 2 * pi),
  cos_year = cos(year_continuous * 2 * pi)
)

# # run models with 1, 2, 3, 4, 5 changepoints
# cp_results_rodent = changepoint_model(ldamodel, x, 1, weights = rep(1,length(year_continuous)))
# save(cp_results_rodent, file = '/Users/renatadiaz/Documents/model-stash/paper_cpt_paper_lda1.Rds')
# rm(cp_results_rodent)
# 
# 
# cp_results_rodent2 = changepoint_model(ldamodel, x, 2, weights = rep(1,length(year_continuous)))
# save(cp_results_rodent2, file = '/Users/renatadiaz/Documents/model-stash/paper_cpt_paper_lda2.Rds')
# rm(cp_results_rodent2)
# 
# cp_results_rodent3 = changepoint_model(ldamodel, x, 3, weights = rep(1,length(year_continuous)))
# save(cp_results_rodent3, file = '/Users/renatadiaz/Documents/model-stash/paper_cpt_paper_lda3.Rds')
# rm(cp_results_rodent3)
# 
# cp_results_rodent4 = changepoint_model(ldamodel, x, 4, weights = rep(1,length(year_continuous)))
# save(cp_results_rodent4, file = '/Users/renatadiaz/Documents/model-stash/paper_cpt_paper_lda4.Rds')
# rm(cp_results_rodent4)
# 
# cp_results_rodent5 = changepoint_model(ldamodel, x, 5, weights = rep(1,length(year_continuous)))
# save(cp_results_rodent5, file = '/Users/renatadiaz/Documents/model-stash/paper_cpt_paper_lda5.Rds')
# rm(cp_results_rodent5)
# 
# cp_results_rodent6 = changepoint_model(ldamodel, x, 6, weights = rep(1,length(year_continuous)))
# save(cp_results_rodent6, file = '/Users/renatadiaz/Documents/model-stash/paper_cpt_paper_lda6.Rds')
# rm(cp_results_rodent6)

## ----paper paper selection-----------------------------------------------
ntopics = ldamodel@k
# Load models and compare deviance to select the best one
load('/Users/renatadiaz/Documents/model-stash/paper_cpt_paper_lda1.Rds')
load('/Users/renatadiaz/Documents/model-stash/paper_cpt_paper_lda2.Rds')
load('/Users/renatadiaz/Documents/model-stash/paper_cpt_paper_lda3.Rds')
load('/Users/renatadiaz/Documents/model-stash/paper_cpt_paper_lda4.Rds')
load('/Users/renatadiaz/Documents/model-stash/paper_cpt_paper_lda5.Rds')
load('/Users/renatadiaz/Documents/model-stash/paper_cpt_paper_lda6.Rds')

# change point model selection
# mean deviance ( -2 * log likelihood) + 2*(#parameters)
mean(cp_results_rodent$saved_lls * -2) + 2*(3*(ntopics-1)*(1+1)+(1))
```

    ## [1] 898.5209

``` r
mean(cp_results_rodent2$saved_lls * -2)+ 2*(3*(ntopics-1)*(2+1)+(2))
```

    ## [1] 834.8776

``` r
mean(cp_results_rodent3$saved_lls * -2)+ 2*(3*(ntopics-1)*(3+1)+(3))
```

    ## [1] 818.1249

``` r
mean(cp_results_rodent4$saved_lls * -2)+ 2*(3*(ntopics-1)*(4+1)+(4))
```

    ## [1] 813.6164

``` r
mean(cp_results_rodent5$saved_lls * -2)+ 2*(3*(ntopics-1)*(5+1)+(5))
```

    ## [1] 817.8308

``` r
mean(cp_results_rodent6$saved_lls * -2)+ 2*(3*(ntopics-1)*(5+1)+(5))
```

    ## [1] 806.4913

``` r
# lowest deviance is again for 4 (or 6) changepoints

rm(list = c('cp_results_rodent', 'cp_results_rodent2', 
            'cp_results_rodent3', 'cp_results_rodent5'))
hist(year_continuous[cp_results_rodent4$saved[,1,]],breaks = seq(1977,2016,.25),xlab='',main='Changepoint Estimate')
```

![](paper_changepoint_files/figure-markdown_github/paper%20lda-1.png)

``` r
hist(year_continuous[cp_results_rodent6$saved[,1,]],breaks = seq(1977,2016,.25),xlab='',main='Changepoint Estimate')
```

![](paper_changepoint_files/figure-markdown_github/paper%20lda-2.png)

``` r
annual_hist(cp_results_rodent4,year_continuous)
```

![](paper_changepoint_files/figure-markdown_github/paper%20lda-3.png)

``` r
# turn changepoint results into data frame
df_4 = as.data.frame(t(cp_results_rodent4$saved[,1,])) %>% melt()
```

    ## No id variables; using all as measure variables

``` r
df_4$value = year_continuous[df_4$value]

# find means on each changepoint:
mean(df_4[df_4$variable=='V1','value']) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ## [1] "14-04-1984"

``` r
mean(df_4[df_4$variable=='V2','value']) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ## [1] "24-11-1992"

``` r
mean(df_4[df_4$variable=='V3','value']) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ## [1] "29-05-1999"

``` r
mean(df_4[df_4$variable=='V4','value']) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ## [1] "01-01-2010"

``` r
# find 95% confidence intervals on each changepoint:
quantile(df_4[df_4$variable=='V1','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ##         2.5%        97.5% 
    ## "11-11-1983" "24-08-1984"

``` r
quantile(df_4[df_4$variable=='V2','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ##         2.5%        97.5% 
    ## "08-08-1988" "26-01-1996"

``` r
quantile(df_4[df_4$variable=='V3','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ##         2.5%        97.5% 
    ## "21-08-1998" "05-11-1999"

``` r
quantile(df_4[df_4$variable=='V4','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ##         2.5%        97.5% 
    ## "20-06-2009" "05-11-2010"

``` r
## 6 cpts 
annual_hist(cp_results_rodent6,year_continuous)
```

![](paper_changepoint_files/figure-markdown_github/paper%20lda-4.png)

``` r
# turn changepoint results into data frame
df_6 = as.data.frame(t(cp_results_rodent6$saved[,1,])) %>% melt()
```

    ## No id variables; using all as measure variables

``` r
df_6$value = year_continuous[df_6$value]

# find 95% confidence intervals on each changepoint:
quantile(df_6[df_6$variable=='V1','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ##         2.5%        97.5% 
    ## "17-02-1978" "30-07-1984"

``` r
quantile(df_6[df_6$variable=='V2','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ##         2.5%        97.5% 
    ## "03-02-1984" "03-03-1995"

``` r
quantile(df_6[df_6$variable=='V3','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ##         2.5%        97.5% 
    ## "13-05-1988" "10-09-1999"

``` r
quantile(df_6[df_6$variable=='V4','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ##         2.5%        97.5% 
    ## "16-08-1994" "20-01-2007"

``` r
quantile(df_6[df_6$variable=='V5','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ##         2.5%        97.5% 
    ## "19-02-1999" "17-10-2009"

``` r
quantile(df_6[df_6$variable=='V6','value'],probs=c(.025,.975)) %>% date_decimal() %>% format('%d-%m-%Y')
```

    ##         2.5%        97.5% 
    ## "22-08-2009" "25-05-2011"
