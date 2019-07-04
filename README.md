# Latent Dirichlet Allocation coupled with Bayesian Time Series analyses

[![Build Status](https://travis-ci.org/weecology/LDATS.svg?branch=master)](https://travis-ci.org/weecology/LDATS)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://raw.githubusercontent.com/weecology/LDATS/master/LICENSE)
[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Codecov test coverage](https://img.shields.io/codecov/c/github/weecology/LDATS/master.svg)](https://codecov.io/github/weecology/LDATS/branch/master)

## Overview

The `LDATS` package provides functionality for analyzing time series of 
high-dimensional data using a two-stage approach comprised of Latent 
Dirichlet Allocation (LDA) and Bayesian time series (TS) analyses.

For a full description of the math underlying the `LDATS` package, see the
[*draft* technical manuscript](https://github.com/weecology/LDATS/blob/master/manuscript/simonis_et_al.pdf).

## Status: Stable Version Available, In Development

The package is currently in active development by the 
[Weecology Team](https://www.weecology.org), although a stable version (0.1.0)
has been submitted to CRAN. The API is well defined at this point and should 
not change substantially.

## Installation

Install the `devtools` package and then run:

```r
devtools::install_github("weecology/LDATS")
```

## Usage

Here is an example of a full LDA-TS analysis using the 
[Portal rodent data](https://github.com/weecology/PortalData):

```r
library(LDATS)
data(rodents)
dtt <- rodents$document_term_table
dct <- rodents$document_covariate_table
weights <- document_weights(dtt)
r_LDATS <- LDA_TS(dtt, dct, topics = 2:5, nseeds = 2, 
                  formulas = c(~1), nchangepoints = 0:1,
                  weights = weights, timename = "newmoon")
```
Which conducts two replicates (`nseeds`) for each of two to five topics in an
LDA model using the document term table, selects the best (AIC) of those, 
then conducts six time series models on it (each of an intercept only and
newmoon-based regression under 0, 1, and 2 changepoints), then selects the 
best (AIC) of the time series, and packages all the models together. This uses
the document term table to weight the samples by their sizes (number of words)
and instructs the function to use the column named `"newmoon"` in the
document covariates table as the time variable.

The resulting object is of class `LDA_TS`, which has a few basic routines 
available:

```r
print(r_LDATS)
```
prints the selected LDA and TS models and 
```r
plot(r_LDATS, timename = "newmoon")
```
produces a 4-panel figure of them a la Figure 1 from
[Christensen et al. 2018](https://doi.org/10.1002/ecy.2373).

## More Information 

Based on inital work using [LDA to analyze time-series data at Portal by Erica
M. Christensen, David J. Harris, and S. K. Morgan 
Ernest](https://github.com/emchristensen/Extreme-events-LDA), which has been
[published in *Ecology*](https://doi.org/10.1002/ecy.2373)

## Acknowledgements 

The motivating study—the Portal Project—has been funded nearly continuously 
since 1977 by the [National Science Foundation](http://nsf.gov/), 
most recently by 
[DEB-1622425](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1622425) 
to S. K. M. Ernest, which also supported (in part) E. Christensen’s
time. Much of the computational work (including time of J. Simonis, D. Harris,
and H. Ye) was supported by the [Gordon and Betty Moore Foundation’s 
Data-Driven Discovery 
Initiative](http://www.moore.org/programs/science/data-driven-discovery) 
through [Grant GBMF4563](http://www.moore.org/grants/list/GBMF4563) to E. P. 
White. R. Diaz was supported in part by a [National Science 
Foundation Graduate Research Fellowship](https://www.nsfgrfp.org/) 
(No. [DGE-1315138](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1315138) 
and [DGE-1842473](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1842473)).

## Author Contributions

**J. L. Simonis** provided insight on LDA applications and feedback on 
technical writing during development of the first version of the LDATS model
and application, led the coding and mathematical development of the model into
an R package, and led writing on the technical manuscript. 
**E. M. Christensen** led the project during development of the first version
of the LDATS model and its application to the Portal data, specifically 
conceiving the project, coding the pipeline wrappers of the analysis, and 
writing and editing the first description of the model and its application 
([Christensen *et al.* 2018](https://doi.org/10.1002/ecy.2373)). **D. J.
Harris** was involved in developing and applying the first version of the 
LDATS model, specifically suggesting the LDA and change point approaches, 
coding the first version of the change point model, and writing and editing 
the first description of the model 
([Christensen *et al.* 2018](https://doi.org/10.1002/ecy.2373)). **R. Diaz**
contributed code to the LDATS package, wrote vignettes, provided insight into
model development, and conducted end-user code application testing. **H. Ye**
contributed code to the LDATS package and insight into data structures and
LDA algorithms. **E. P. White** helped design, troubleshoot, and supervise 
initial methods development and provided big-picture feedback on 
development of the R package. **S. K. Morgan Ernest** provided managerial
oversight and feedback on the project in both the initial and second stages
of LDATS development, tested applications of the code to data sets, and
assisted with writing and editing of the first description of the
model and its application 
([Christensen *et al.* 2018](https://doi.org/10.1002/ecy.2373)).

