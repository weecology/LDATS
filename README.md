# Latent Dirichlet Allocation coupled with Bayesian Time Series analyses

[![Build Status](https://travis-ci.org/weecology/LDATS.svg?branch=master)](https://travis-ci.org/weecology/LDATS)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://raw.githubusercontent.com/weecology/LDATS/master/LICENSE)
[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)

## Overview

The `LDATS` package provides functionality for analyzing time series of 
high-dimensional data using a two-stage approach comprised of Latent 
Dirichlet Allocation (LDA) and Bayesian time series (TS) analyses.

## Status: In Development

The package is currently in active development by the 
[Weecology Team](https://www.weecology.org), in advance of submission to 
CRAN. The API is well defined at this point and should not change 
substantially, but we are presently engaged in testing and documentation,
and further edits may occur before submission to CRAN. Therefore, any 
output should still be considered ***provisional***. 

Previous versions of code and implementation can be found in the 
[repository for Christensen *et al.* 2018](https://github.com/emchristensen/Extreme-events-LDA)


## Mathematical background

For a full description of the math underlying the `LDATS` package, see the
*draft* [technical manuscript](https://github.com/weecology/LDATS/blob/master/manuscript/simonis_et_al.pdf).

## Contributing

Folks interested in contributing to development, should see 
[issues](https://github.com/weecology/LDATS/issues) for specific pre-package 
development tasks.

## Installation

Install the `devtools` package and then run:

```
devtools::install_github("weecology/LDATS")
```

## Current Usage

Here is an example of an LDA and an LDATS using the Portal rodent data:

```
data(rodents)
dtt <- rodents$document_term_table
dct <- rodents$document_covariate_table

r_LDA <- LDATS::LDA_set(dtt, topics = 2:5, nseeds = 2)
r_LDATS <- LDATS::LDA_TS(dtt, dct, topics = 2:5, nseeds = 2, 
                         formulas = ~ 1, nchangepoints = 1)

```

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
contributed code to the LDATS package, provided insight into model 
development, and conducted end-user code application testing. **H. Ye**
contributed code to the LDATS package and insight into data structures and
LDA algorithms. **E. P. White** helped design, troubleshoot, and supervise 
initial methods development and provided big-picture feedback on 
development of the R package. **S. K. Morgan Ernest** provided managerial
oversight and feedback on the project in both the initial and second stages
of LDATS development, tested applications of the code to data sets, and
assisted with writing and editing of the first description of the
model and its application 
([Christensen *et al.* 2018](https://doi.org/10.1002/ecy.2373)).

