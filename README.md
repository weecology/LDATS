# Two-Stage LDA - Time Series Analyses

## Overview

The `LDATS` package provides functionality for analyzing time series of 
multivariate data using a two-stage approach comprised of Latent Dirichlet
Allocation (LDA) and Bayesian time series (change point and continuous change)
analyses.

## Status: In Development

The package is currently in development by the [Weecology 
Team](https://www.weecology.org). As such, any output should be considered
***extremely provisional***. 

Previous versions of code and implementation (Christensen et al. manuscript
and downsampling analyes) are located within `/previous_work`.

## Contributing

Folks interested in contributing to development, should see 
[issues](https://github.com/weecology/LDATS/issues) for specific pre-package 
development tasks.

### Style
For all code, please follow the [Weecology Code Style 
Guide](https://github.com/weecology/lab-wiki/wiki/Code-style-guide)

## Installation

Install the `devtools` package and then run:

```
devtools::install_github("weecology/LDATS")
```

## Current Usage

Here is an example of an LDA and an LDATS using the Portal rodent data:

```
data(rodents)
lda_data <- dplyr::select(rodents, -c(newmoon, date, plots, traps))
ts_data <- data.frame("time" = data.frame(rodents)[ , "newmoon"])

r_LDA <- LDATS::LDA(data = lda_data, ntopics = 2:5, nseeds = 2, ncores = 4)
r_LDATS <- LDATS::LDA_TS(lda_data, ts_data, formula = c(~1, ~time),
                         ntopics = 2:5, nseeds = 2, ncores = 4, nit = 100)
```

## More Information 

Based on inital work using [LDA to analyze time-series data at Portal by Erica
Christensen, Dave Harris, and Morgan 
Ernest](https://github.com/emchristensen/Extreme-events-LDA).
