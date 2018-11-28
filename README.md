# Two-Stage LDA - Time Series Analyses

[![Build Status](https://travis-ci.org/weecology/LDATS.svg?branch=master)](https://travis-ci.org/weecology/LDATS)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://raw.githubusercontent.com/weecology/LDATS/master/LICENSE)

## Overview

The `LDATS` package provides functionality for analyzing time series of 
multivariate data using a two-stage approach comprised of Latent Dirichlet
Allocation (LDA) and Bayesian time series (TS) analyses.

## Status: In Development

The package is currently *in development* by the [Weecology 
Team](https://www.weecology.org). As such, any output should be considered
***provisional***. 

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
lda_data <- rodents$document_term_table
ts_data <- rodents$document_covariate_table

r_LDA <- LDATS::LDA(data = lda_data, ntopics = 2:5, nseeds = 2, ncores = 4)
```

## More Information 

Based on inital work using [LDA to analyze time-series data at Portal by Erica
Christensen, Dave Harris, and Morgan 
Ernest](https://github.com/emchristensen/Extreme-events-LDA).

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
and DGE-1842473).

## Author Contributions

J. L. Simonis provided insight on LDA applications and feedback on technical
writing during development of the first version of the LDATS model and
application, led the coding and mathematical development of the model into an
R package, and led writing on the technical manuscript. E. M. Christensen led
the project during development of the first version of the LDATS model and its
application to the Portal data, specifically conceiving the project, coding
the pipeline wrappers of the analysis, and writing and editing the first
description of the model and its application (Christensen et al. 2018). D. J.
Harris was involved in developing and applying the first version of the LDATS
model, specifically suggesting the LDA and change point approaches, coding the
first version of the change point model, and writing and editing the first
description of the model (Christensen et al. 2018). R. Diaz contributed code
to the LDATS package, provided insight into model development, and conducted
end-user code application testing. H. Ye contributed code to the LDATS package 
and insight into data structures and LDA algorithms. E. P. White helped design,
troubleshoot, and supervise initial methods development and provided 
big-picture feedback on development of the R package. S. K. Morgan Ernest 
provided managerial oversight and feedback on the project in both the initial 
and second stages of LDATS development, tested applications of the code to data
sets, and assisted with writing and editing of the first description of the
model and its application (Christensen et al. 2018).

