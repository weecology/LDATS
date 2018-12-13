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

Christensen 2018 analysis
-------------------------

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
#> [1] 16
ldats_dat[ which(compare_rows == F), ]
#>     BA DM DO DS NA. OL OT PB PE PF PH PI PL PM PP RF RM RO SF SH SO
#> 241  0 31  0  0   1  0  2  3  0  0  0  0  0  0  3  0  0  0  0  0  0
#> 267  0  9  2  0   1  0  0 18  0  0  0  0  0  0 11  0  0  0  0  1  0
#> 277  0 15  2  0   3  0  0  6  0  0  0  0  0  0 11  0  0  0  0  1  0
#> 278  0  9  5  0   0  0  3 11  0  2  0  0  0  0 18  0  0  0  0  0  0
#> 283  0  9  5  0   3  0  7  7  2  0  0  0  0  0  2  0  0  0  0  0  0
#> 284  0  9  6  0   1  0  7  6  3  0  0  0  0  0  0  0  0  0  0  0  0
#> 311  0 13  4  0   1  0  4  2  0  0  0  0  0  0 10  0  0  0  0  0  0
#> 318  0  2  2  0   0  0  3 18  2  0  0  0  0  0  6  0  0  0  0  0  0
#> 321  0 13 13  0   0  0  1  8  0  0  0  0  0  0 11  0  0  0  0  0  0
#> 337  0  6  6  0   0  0  0  2  0  0  0  0  0  0 29  0  0  0  0  0  0
#> 339  0  5  3  0   0  0  1  2  0  0  0  0  0  0  9  0  0  0  0  0  0
#> 344  0  5  6  0   0  4  3 17  3  0  0  0  0  0 13  0  1  0  0  3  0
#> 351  1  5  7  0   3  1  5 36  0  0  0  0  0  1 35  0  2  0  0  1  0
#> 400  0 10  0  0   0  1  1  1  5  0  0  0  0  0  2  0  0  0  0  0  0
#> 405  0 28  3  0   1  1  0  3  3  2  0  0  0  0 47  0  1  0  0  0  0
#> 433  4 10  5  0   1  4  4  3  6  0  0  0  0  0  9  0  1  1  1  0  0
ch_dat[ which(compare_rows == F), ]
#>     BA DM DO DS NA OL OT PB PE PF PH PI PL PM PP RF RM RO SF SH SO
#> 241  0 62  0  0  2  0  4  6  0  0  0  0  0  0  6  0  0  0  0  0  0
#> 267  0 14  3  0  2  0  0 29  0  0  0  0  0  0 18  0  0  0  0  2  0
#> 277  0 30  4  0  6  0  0 12  0  0  0  0  0  0 22  0  0  0  0  2  0
#> 278  0 18 10  0  0  0  6 22  0  4  0  0  0  0 36  0  0  0  0  0  0
#> 283  0 18 10  0  6  0 14 14  4  0  0  0  0  0  4  0  0  0  0  0  0
#> 284  0 18 12  0  2  0 14 12  6  0  0  0  0  0  0  0  0  0  0  0  0
#> 311  0 15  5  0  1  0  5  2  0  0  0  0  0  0 11  0  0  0  0  0  0
#> 318  0  4  4  0  0  0  6 36  4  0  0  0  0  0 12  0  0  0  0  0  0
#> 321  0 15 15  0  0  0  1  9  0  0  0  0  0  0 13  0  0  0  0  0  0
#> 337  0 12 12  0  0  0  0  4  0  0  0  0  0  0 58  0  0  0  0  0  0
#> 339  0 10  6  0  0  0  2  4  0  0  0  0  0  0 18  0  0  0  0  0  0
#> 344  0 10 12  0  0  8  6 34  6  0  0  0  0  0 26  0  2  0  0  6  0
#> 351  2 10 14  0  6  2 10 72  0  0  0  0  0  2 70  0  4  0  0  2  0
#> 400  0 20  0  0  0  2  2  2 10  0  0  0  0  0  4  0  0  0  0  0  0
#> 405  0 37  4  0  1  1  0  4  4  3  0  0  0  0 63  0  1  0  0  0  0
#> 433  8 20 10  0  2  8  8  6 12  0  0  0  0  0 18  0  2  2  2  0  0
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

ldats_dat <- ch_dat
rodents$document_term_table <- ch_dat
```
