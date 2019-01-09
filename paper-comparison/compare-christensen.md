Comparison with Christensen paper
================
Renata Diaz
12/13/2018

Side-by-side comparison of LDATS results (and Portal rodents vignette) with analysis from Christensen et al 2018.

LDATS Installation
------------------

To obtain the most recent version of **LDATS**, install the most recent version from GitHub:

``` r
install.packages("devtools")
devtools::install_github("weecology/LDATS")
```

``` r
library(LDATS)
```

Christensen 2018 analysis files
-------------------------------

Download Christensen 2018 analysis scripts & data files from [Extreme-events-LDA repo:](https://github.com/emchristensen/Extreme-events-LDA)

Main Analysis Scripts:

-   rodent\_LDA\_analysis.R main script for analyzing rodent community change using LDA

-   rodent\_data\_for\_LDA.R contains a function that creates the rodent data table used in analyses

-   AIC\_model\_selection.R contains functions for calculating AIC for different candidate LDA models

-   changepointmodel.r contains change-point model code

-   LDA-distance.R function for computing Hellinger distance analyses

Data:

-   Rodent\_table\_dat.csv table of rodent data, created by rodent\_data\_for\_LDA.R

Figure scripts:

-   LDA\_figure\_scripts.R contains functions for making main plots in manuscript (Fig 1). Called from rodent\_LDA\_analysis.R

``` r

files_to_download <- c('rodent_LDA_analysis.r', 'rodent_data_for_LDA.r', 'AIC_model_selection.R', 'changepointmodel.r', 'LDA-distance.R', 'Rodent_table_dat.csv', 'LDA_figure_scripts.R')

for(i in 1:length(files_to_download)) {
  download.file(url = paste0("https://raw.githubusercontent.com/emchristensen/Extreme-events-LDA/master/", files_to_download[i]),
                destfile = paste0('christensen-ecology-files/', files_to_download[i]))
}

rm(files_to_download)
rm(i)
```

Data
----

The Portal rodents control data is included in the LDATS package:

``` r

data(rodents)

ldats_dat <- rodents[[1]]
ldats_dates <- rodents[[2]]
head(ldats_dat)
#>   BA DM DO DS NA. OL OT PB PE PF PH PI PL PM PP RF RM RO SF SH SO
#> 1  0 13  0  2   2  0  0  0  1  1  0  0  0  0  3  0  0  0  0  0  0
#> 2  0 20  1  3   2  0  0  0  0  4  0  0  0  0  2  0  0  0  0  0  0
#> 3  0 21  0  8   4  0  0  0  1  2  0  0  0  0  1  0  0  0  0  0  0
#> 4  0 21  3 12   4  2  3  0  1  1  0  0  0  0  0  0  0  0  0  0  0
#> 5  0 16  1  9   5  2  1  0  0  2  0  0  0  0  0  0  1  0  0  0  0
#> 6  0 17  1 13   5  1  5  0  0  3  0  0  0  0  0  0  0  0  0  0  0
```

Load the data used in Christensen et al:

``` r
source('christensen-ecology-files/rodent_data_for_LDA.r')
#> Loading required package: bitops
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
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
```

Compare paper data to LDATS data:

``` r

compare <- ldats_dat == ch_dat

compare_rows <- vector(length = nrow(compare)) 

for(i in 1:nrow(compare)) {
  if(sum(compare[i, ]) == 21) {
    compare_rows[i] <- TRUE
  } else {
    compare_rows[i] <- FALSE
  }
}

length(which(compare_rows == F))
ldats_dat[ which(compare_rows == F), ]
ch_dat[ which(compare_rows == F), ]
```

There are 16 rows where the data included in LDATS differs from the paper data. This looks like it is because the LDATS data is not adjusted to account for trapping effort, but the paper data has divided all census counts by the actual number of plots trapped and multiplied by 8 to account for incompletely-trapped censuses.

Double-check:

``` r

source('christensen-ecology-files/rodent_data_for_LDA.r')
diagnose_dat = create_rodent_table(period_first = 1,
                          period_last = 436,
                          selected_plots = c(2,4,8,11,12,14,17,22),
                          selected_species = c('BA','DM','DO','DS','NA','OL','OT','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO','SF','SH','SO'), diagnose = T)

adjusted_table <- diagnose_dat[[1]]
raw_table <- diagnose_dat[[2]]
nplots <- diagnose_dat[[3]]
nplots <- nplots[ which(nplots$period <= 436), ]

compare_raw <- raw_table == ldats_dat
length(which(compare_raw == F))

which(compare_rows == F) 
which(nplots$x < 8)
```

I added an argument 'diagnose' to the create\_rodent\_table() function. If T, returns the raw table & \# of plots tables as well as the adjusted table.

The lines that are different are the lines that were adjusted; if you don't adjust the paper data you get the same table as the LDATS data.

For now I will use the *adjusted rodent table* because this will avoid artificially low abundances for incompletely trapped periods.

``` r

ldats_dat <- as.matrix(ch_dat)
ldats_dat <- apply(ldats_dat, c(1,2), FUN = as.integer)
rodents$document_term_table <- ldats_dat
rm(ldats_dat)
rm(ldats_dates)
rm(ch_dat)
```

LDA
---

While LDATS can run start-to-finish with `LDATS::LDA_TS`, I'm first going to step through function-by-function to isolate differences with the paper.

``` r

ldats_ldas <- LDATS::LDA_set(document_term_table = rodents$document_term_table, topics = c(2:6), nseeds = 100)
ldats_lda_selected <- LDATS::select_LDA(LDA_models = ldats_ldas)
```

Paper LDAS:

``` r
source('christensen-ecology-files/AIC_model_selection.R')
source('christensen-ecology-files/LDA-distance.R')

# Fit a bunch of LDA models with different seeds
# Only use even numbers for seeds because consecutive seeds give identical results
seeds = 2*seq(200)

# repeat LDA model fit and AIC calculation with a bunch of different seeds to test robustness of the analysis
best_ntopic = repeat_VEM(rodents[[1]],
                         seeds,
                         topic_min=2,
                         topic_max=6)
hist(best_ntopic$k,breaks=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5),xlab='best # of topics', main='')

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
best_seed = calculate_LDA_distance(rodents[[1]],seeds_4topics)
mean_dist = unlist(best_seed)[2]
max_dist = unlist(best_seed)[3]

# ==================================================================
# 3. run LDA model
# ==================================================================
ntopics = 4
SEED = unlist(best_seed)[1]  # For the paper, I use seed 206
ldamodel = LDA(rodents[[1]],ntopics, control = list(seed = SEED),method='VEM')
```

### Plot LDAS

``` r
# Paper
plot(ldamodel, cols = NULL, option = "D")
```

![](compare-christensen_files/figure-markdown_github/plot%20LDAs-1.png)

``` r

# LDATS
plot(ldats_lda_selected[[1]], cols = NULL, option = "D")
```

![](compare-christensen_files/figure-markdown_github/plot%20LDAs-2.png)

Changepoint models
------------------

I am going to compare four combinations of LDA + changepoint models:

-   LDATS LDA + LDATS changepoint
-   LDATS LDA + paper changepoint
-   Paper LDA + LDATS changepoint
-   Paper LDA + paper changepoint

There is the additional wrinkle of `document_term_weights`. The paper weighted all sample periods equally, wheras LDATS can weight sample periods according to how many individuals were captured. We now believe it is more appropriate to weight periods proportional to captures. However, for the purposes of comparison, I will continue to set all weights = 1 for both changepoint models. For an example of LDATS run with proportional weights, see the rodents vignette. \[?\]

For narratives running all these models, see `ldats_changepoint.md` and `paper_changepoint.md`.

All the models found 4 changepoints.

### Compare model outcomes

``` r
ldats_ldats <- read.csv('ldats_ldats_dates.csv')
ldats_paper <- read.csv('ldats_paper_dates.csv')
paper_ldats <- read.csv('paper_ldats_dates.csv')
paper_paper <- read.csv('paper_paper_dates.csv')


estimates <- c(1:4)
estimates <- as.data.frame(estimates)

# Naming convention: ldamethod_cptmethod
# so ldats_paper is the ldats lda and paper cpt, etc
estimates$ldats_ldats <- ldats_ldats$newmoondate
estimates$ldats_paper <- ldats_paper$cpt_dates
estimates$paper_ldats <- paper_ldats$newmoondate
estimates$paper_paper <- paper_paper$cpt_dates

estimates
#>   estimates ldats_ldats ldats_paper paper_ldats paper_paper
#> 1         1  1984-01-02  12-01-1984  1984-04-30  14-04-1984
#> 2         2  1991-05-13  02-06-1991  1992-10-25  24-11-1992
#> 3         3  1998-12-18  04-01-1999  1999-07-12  29-05-1999
#> 4         4  2009-11-16  10-11-2009  2010-01-15  01-01-2010
```

All of the outcomes are reasonably close.

The choice of changepoint model seems to matter less than the choice of LDA model.
