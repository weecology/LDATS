# LDATS (development version)

# LDATS 0.1.0

## Code structure
* Creation of a [standard API and code 
pipeline](https://weecology.github.io/LDATS/articles/LDATS_codebase.html) for 
all components of the LDATS analysis.
* Substantial refactor of the underlying code from hardcoded to generalized functions.
* Development of checking functions used to run the basic structural checks on the function inputs.
* Inclusion of control options lists for the LDA stage, TS stage, and overall to reduce the length of input lists.

## Full inclusion of functions
* All functions used in the code base are now exported, documented, and tested.

## LDA model AIC calculation
* `AIC.LDA_VEM()` now uses the number of parameters as reported from `logLik` to calculate AIC.
* Previous by-hand calculations of AIC included variational parameters that are integrated out of the model in the total parameter count.

## Regressor estimates
* Time series models allow for flexible covariate set for regression via formula inputs to the top-level functions.
* The time series model code now also includes estimation of the parameters defining the between-change point regressions (*i.e.*, the regressors).
* Regressor estimates come as marginal posterior distributions, and are calculated by unconditioning the estimates generated under known change points. 

## Document weighting
* `document_weights()` function is provided to allow for appropriate weighting of documents by their sizes (number of words) so that an average-length document is 1.
* Document weighting is done automatically by default, which is easily undone by using `weights = NULL`. 

## ptMCMC functionality
* The ptMCMC code has been refactored into functions, many of which are generalized to use in other contexts.
* The temperature schema is fully controllable via arguments to the TS controls list
* Burn-in removal and thinning of final chains is controllable via the TS controls list

## Optional memoisation
* Memoisation of `multinom_TS()` and `multinom_TS_chunk()` now is optional via `memoise_fun()` and is controlled through the TS controls list.

## Plotting functions
* `LDA_set()`, `LDA_TS()`, and `TS()` now all have default plotting options on their outputs.
* `plot.TS()` provides MCMC diagnostic plots and summary plots.
* `plot.LDA_TS()` plots produce the combination of plots.

## Rodents data set
* Portal rodent data from [Christensen *et al.* (2018)](https://doi.org/10.1002/ecy.2373) are now provided in a pre-formatted and ready-to-roll data object.
* The data in Christensen *et al.* 2018 are scaled according to trapping effort. The data included in LDATS are not, to allow for appropriate weighting.
* See `data(rodents)`

## Comparison with [Christensen *et al.* (2018)](https://doi.org/10.1002/ecy.2373)
* The paper-comparison vignette provides a step-by-step comparison of the LDATS pipeline to the analysis in Christensen *et al.* 2018. The differences are as follows:

      * The `document_term_table` in Christensen *et al.* 2018 was adjusted to account for variable trapping effort. The data included in LDATS are not adjusted, so that sampling periods can be weighted appropriately.
      * The LDA model selection criterion has changed (see LDA model AIC calculation, above), so that LDATS now identifies 6 topics compared to the 4 topics found in the paper.
      * LDATS will by default weight sampling periods according to the number of terms (see Document weighting, above). 
      * Despite these changes, the updated LDATS pipeline gives qualitatively similar results to the analysis in Christensen *et al.* 2018. 

# [LDATS 0.0.1](https://github.com/weecology/LDATS/commit/326506b9d7fb3e0223948d0245381963f83a2b37) 2017-11-16

* Beginning initial development of package from [original
code]((https://github.com/emchristensen/Extreme-events-LDA)) used in 
[Christensen *et al.* (2018)](https://doi.org/10.1002/ecy.2373).
