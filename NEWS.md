# LDATS (development version)

# LDATS 0.1.0

## Code structure
* Creation of a [standard API and code 
pipeline](https://weecology.github.io/LDATS/articles/LDATS_codebase.html) for 
all components of the LDATS analysis.
* Substantial refactor of the underlying code from hardcoded to generalized functions.
* Development of checking functions used to run the basic structural checks on the function inputs.
* Inclusion of control options lists for the LDA stage, TS stage, and overall to reduce the length of input lists.

## LDA Model AIC calculation
* `AIC.LDA_VEM` now uses the number of parameters as reported from `logLik` to calculate AIC.
* Previous by-hand calculations of AIC included variational parameters that are integrated out of the model in the total parameter count.

## Regressor estimates
* The time series model code now also includes estimation of the parameters defining the between-change point regressions (*i.e.*, the regressors).
* Regressor estimates come as marginal posterior distributions, and are calculated by unconditioning the estimates generated under known change points. 

## Plotting functions
* `LDA_set`, `LDA_TS`, and `TS` now all have default plotting options on their outputs.
* `TS` provides MCMC diagnostic plots and summary plots.
* `LDA_TS` plots produce the combination of plots.

# rodents data set
* Portal rodent data from [Christensen *et al.* (2018)](https://doi.org/10.1002/ecy.2373) are now provided in a pre-formatted and ready-to-roll data object.

# [LDATS 0.0.1](https://github.com/weecology/LDATS/commit/326506b9d7fb3e0223948d0245381963f83a2b37) 2017-11-16

* Beginning initial development of package from [original
code]((https://github.com/emchristensen/Extreme-events-LDA)) used in 
[Christensen *et al.* (2018)](https://doi.org/10.1002/ecy.2373).