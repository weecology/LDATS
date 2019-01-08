LDATS changepoint
================
Renata Diaz
1/8/2019

LDATS LDA
---------

### Setup

To run changepoint models:

``` r
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

period_dates = filter(moondat,period %in% rownames(dat)) %>% select(newmoonnumber, period,date)
dates = period_dates$date

ch_dat <- dat

ch_dates <- dates
ch_newmoons <- period_dates$newmoonnumber

rm(list = c('dat', 'dates',
            'create_rodent_table'))

ldats_dat <- as.matrix(ch_dat)
ldats_dat <- apply(ldats_dat, c(1,2), FUN = as.integer)
rodents$document_term_table <- ldats_dat
rm(ldats_dat)
rm(ldats_dates)
rm(ch_dat)

#### Load LDA ####
load('/Users/renatadiaz/Documents/model-stash/ldats_lda.Rds')
year_continuous = 1970 + as.integer(julian(ch_dates)) / 365.25
x = data.frame(
  year_continuous = year_continuous,
  sin_year = sin(year_continuous * 2 * pi),
  cos_year = cos(year_continuous * 2 * pi)
)

x$newmoon <- ch_newmoons
x <- select(x, newmoon, sin_year, cos_year)
```

### Run the TS models

``` r
#### Run LDATS changepoint ####

ldats_ldats_cpt <- TS_on_LDA(LDA_models = ldats_lda_selected, 
                             document_covariate_table = x,
                             formulas = ~ sin_year + cos_year,
                             nchangepoints = 1:6,
                             weights = NULL)


save(ldats_ldats_cpt, file = '/Users/renatadiaz/Documents/model-stash/ldats_lda_ldats_cpt.Rds')
rm(ldats_ldats_cpt)
rm(ldats_lda_selected)
```

### Select model

``` r
ldats_ldats_cpt_selected <- select_TS(ldats_ldats_cpt)
```

### Show results

``` r
plot(ldats_ldats_cpt_selected)
```

![](ldats_changepoint_files/figure-markdown_github/show%20ldats%20lda%20+%20ldats%20cpt-1.png)

``` r
ldats_ldats_cpt_selected$formula
```

    ## gamma ~ sin_year + cos_year
    ## <environment: 0x7fc317e62980>

``` r
ldats_ldats_cpt_selected$nchangepoints
```

    ## [1] 4

``` r
cpt_means <- ldats_ldats_cpt_selected$rho_summary$Mean

cpt_means <- floor(cpt_means)
cpt_means <- as.data.frame(cpt_means)
colnames(cpt_means) <- 'newmoonnumber'
cpt_means$cpt <- 1:ldats_ldats_cpt_selected$nchangepoints

cpt_dates <- as.data.frame(moondat) %>%
  dplyr::select(newmoonnumber, newmoondate) %>%
  dplyr::left_join(cpt_means, by = 'newmoonnumber') %>%
  dplyr::filter(!is.na(cpt))

cpt_dates
```

    ##   newmoonnumber newmoondate cpt
    ## 1            81  1984-01-02   1
    ## 2           172  1991-05-13   2
    ## 3           266  1998-12-18   3
    ## 4           401  2009-11-16   4

``` r
write.csv(cpt_dates, 'ldats_ldats_dates.csv')
```

Paper LDA
---------

### Load LDA

``` r
#### Load paper LDA ####
load('/Users/renatadiaz/Documents/model-stash/paper_lda.Rds')
```

### Run TS models

``` r
paper_ldats_cpt <- TS_on_LDA(LDA_models = ldamodel, document_covariate_table = x, formulas = ~ sin_year + cos_year,nchangepoints = 1:6, weights = NULL)
save(paper_ldats_cpt, file = '/Users/renatadiaz/Documents/model-stash/paper_lda_ldats_cpt.Rds')
```

``` r
paper_ldats_cpt_selected <- select_TS(paper_ldats_cpt)
```

### Show results

``` r
plot(paper_ldats_cpt_selected)
```

![](ldats_changepoint_files/figure-markdown_github/show%20paper%20lda%20+%20ldats%20cpt-1.png)

``` r
paper_ldats_cpt_selected$formula
```

    ## gamma ~ sin_year + cos_year
    ## <environment: 0x7fc310be3660>

``` r
paper_ldats_cpt_selected$nchangepoints
```

    ## [1] 4

``` r
cpt_means <- paper_ldats_cpt_selected$rho_summary$Mean


cpt_means <- floor(cpt_means)
cpt_means <- as.data.frame(cpt_means)
colnames(cpt_means) <- 'newmoonnumber'
cpt_means$cpt <- 1:paper_ldats_cpt_selected$nchangepoints

cpt_dates <- as.data.frame(moondat) %>%
  dplyr::select(newmoonnumber, newmoondate) %>%
  dplyr::left_join(cpt_means, by = 'newmoonnumber') %>%
  dplyr::filter(!is.na(cpt))

cpt_dates
```

    ##   newmoonnumber newmoondate cpt
    ## 1            85  1984-04-30   1
    ## 2           190  1992-10-25   2
    ## 3           273  1999-07-12   3
    ## 4           403  2010-01-15   4

``` r
write.csv(cpt_dates, 'paper_ldats_dates.csv')
```
