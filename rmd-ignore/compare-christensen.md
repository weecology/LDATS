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
                destfile = paste0('christensen-ecology/', files_to_download[i]))
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
