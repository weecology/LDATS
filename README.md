# Latent Dirichlet Allocation coupled with Bayesian Time Series analyses 

[![R-CMD-check](https://github.com/weecology/LDATS/actions/workflows/r-cmd-check.yaml/badge.svg)](https://github.com/weecology/LDATS/actions/workflows/r-cmd-check.yaml) 
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://raw.githubusercontent.com/weecology/LDATS/main/LICENSE)
[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html#maturing)
[![Codecov test coverage](https://img.shields.io/codecov/c/github/weecology/LDATS/main.svg)](https://app.codecov.io/github/weecology/LDATS/branch/main)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/LDATS)](https://CRAN.R-project.org/package=LDATS)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3286617.svg)](https://zenodo.org/record/3715386)

## Overview

The **`LDATS`** package provides functionality for analyzing time series of high-dimensional data using a two-stage approach comprised of Latent Dirichlet Allocation (LDA) and Bayesian time series (TS) analyses.

For a full description of the math underlying the **`LDATS`** package, see the [technical document](https://github.com/weecology/LDATS/blob/main/LDATS_model.pdf).

## Status: Stable Version Available, Continuing Development

A stable version of LDATS is available on [CRAN](https://CRAN.R-project.org/package=LDATS), but the package is actively being developed by the [Weecology Team](https://www.weecology.org).
The API is well defined at this point and should not change substantially.

## Installation

You can install the stable version of **`LDATS`** from CRAN with:

To obtain the current development version of **`LDATS`** from GitHub, install the **`devtools`** package and then use it to install **`LDATS`**:

```r
install.packages("devtools")
devtools::install_github("weecology/LDATS")
```

## Usage

Here is an example of a full LDA-TS analysis using the 
[Portal rodent data](https://github.com/weecology/PortalData):

```r
library(LDATS)
data(rodents)
r_LDATS <- LDA_TS(rodents, topics = 2:5, nseeds = 2, formulas = ~1,  
                  nchangepoints = 0:1, timename = "newmoon")
```
Which conducts two replicates (`nseeds`) for each of two to five topics in an LDA model using the document term table, selects the best (AIC) of those, then conducts two time series models on it (an intercept-only model under 0 and 1 changepoints), then selects the best (AIC) of the time series, and packages all the models together. 
This uses the document term table to weight the samples by their sizes (number of words) and instructs the function to use the column named `"newmoon"` in the document covariates table as the time variable.

The resulting object is of class `LDA_TS`, which has a few basic routines 
available:

```r
print(r_LDATS)
```
prints the selected LDA and TS models and 
```r
plot(r_LDATS)
```
produces a 4-panel figure of them a la Figure 1 from
[Christensen et al. 2018](https://pubmed.ncbi.nlm.nih.gov/29718539/).

## More Information 

Based on initial work using [LDA to analyze time-series data at Portal by Erica M. Christensen, David J. Harris, and S. K. Morgan Ernest](https://github.com/emchristensen/Extreme-events-LDA), which has been [published in *Ecology*](https://pubmed.ncbi.nlm.nih.gov/29718539/)

## Acknowledgements 

The motivating study—the Portal Project—has been funded nearly continuously since 1977 by the [National Science Foundation](https://www.nsf.gov/), most recently by [DEB-1622425](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1622425) to S. K. M. Ernest, which also supported (in part) E. Christensen’s time. 
Much of the computational work (including time of J. Simonis, D. Harris, and H. Ye) was supported by the [Gordon and Betty Moore Foundation’s Data-Driven Discovery Initiative](https://www.moore.org/initiative-strategy-detail?initiativeId=data-driven-discovery) through [Grant GBMF4563](https://www.moore.org/grant-detail?grantId=GBMF4563) to E. P. White. 
R. Diaz was supported in part by a [National Science Foundation Graduate Research Fellowship](http://www.nsfgrfp.org/) (No. [DGE-1315138](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1315138) and [DGE-1842473](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1842473)).

## Author Contributions

**J. L. Simonis** provided insight on LDA applications and feedback on technical writing during development of the first version of the LDATS model and application, led the coding and mathematical development of the model into an R package, and led writing on the technical model document. 
**E. M. Christensen** led the project during development of the first version of the LDATS model and its application to the Portal data, specifically conceiving the project, coding the pipeline wrappers of the analysis, and writing and editing the first description of the model and its application ([Christensen *et al.* 2018](https://pubmed.ncbi.nlm.nih.gov/29718539/)). 
**D. J. Harris** was involved in developing and applying the first version of the LDATS model, specifically suggesting the LDA and change point approaches, coding the first version of the change point model, and writing and editing the first description of the model ([Christensen *et al.* 2018](https://pubmed.ncbi.nlm.nih.gov/29718539/)). 
**R. Diaz** contributed code to the LDATS package, wrote vignettes, provided insight into model development, and conducted extensive end-user code application testing. 
**H. Ye** contributed code to the LDATS package, insight into data structures and LDA algorithms, and significant feedback on vignettes. 
**E. P. White** helped design, troubleshoot, and supervise initial methods development; provided big-picture feedback on development of the R package; contributed end-user application testing; and gave substantial editing feedback on the technical document. 
**S. K. Morgan Ernest** provided managerial oversight and feedback on the project in both the initial and second stages of LDATS development, tested applications of the code to data sets, and assisted with writing and editing of the first description of the model and its application ([Christensen *et al.* 2018](https://pubmed.ncbi.nlm.nih.gov/29718539/)) as well as the technical model document.

