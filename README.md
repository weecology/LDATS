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

The current development code work flow includes two scripts in the highest 
level of the repo: `package_dev.R` and `working_space.R`. `package_dev.R` 
faciliates building an in-development local version of the package. Any 
dependency packages should be included within the `pkg_depends` vector, 
as this populates the `DESCRIPTION` file. `working_space.R` is a sandbox 
script that walks through example implementations of the code for testing
and development purposes. Once functions are complete and documented, they
should be located in a script within `/R`. 

### Style
For all code, please follow the [Weecology Code Style 
Guide](https://github.com/weecology/lab-wiki/wiki/Code-style-guide)

## Installation

Install the `devtools` package and then run:

```
devtools::install_github("weecology/portalr")
```

## Current Usage

To be added shortly. 

## More Information 

Based on inital work using [LDA to analyze time-series data at Portal by Erica
Christensen, Dave Harris, and Morgan 
Ernest](https://github.com/emchristensen/Extreme-events-LDA).
